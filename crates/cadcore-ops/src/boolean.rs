//! Solid Boolean half-space cut for cadcore B-Rep solids.
//!
//! # Strategy
//!
//! Each *cylinder face* in the [`BRep`] is classified relative to the cutting
//! plane and rebuilt independently.  This works because each cylinder face in
//! the existing sweep pipeline corresponds to one straight segment of the
//! polyline вЂ” solids built segment-by-segment can be trimmed face-by-face.
//!
//! Three cases per cylinder:
//!
//! | axis вЉҐ plane.normal? | dist(axis, plane) | result                              |
//! |----------------------|-------------------|-------------------------------------|
//! | yes (LATERAL)        | > radius          | kept                                |
//! | yes (LATERAL)        | < в€’radius         | dropped                             |
//! | yes (LATERAL)        | |dist| < radius   | partial cylinder + flat chord       |
//! | no  (AXIAL)          | both ends kept    | kept                                |
//! | no  (AXIAL)          | both ends discard | dropped                             |
//! | no  (AXIAL)          | crosses plane     | truncated + flat disk cap           |
//!
//! The same cylinder may be cut by multiple planes; planes are processed
//! sequentially.

use cadcore_geom::{CylSurf, Plane3};
use cadcore_math::{Frame3, Point3, UnitVec3};
use cadcore_topo::{
    BRep, Face, FaceBoundary, FaceExtent, FaceGeom, FaceId, FaceNormal, Shell, Solid, SolidId,
};

use crate::sweep::ClipPlane;

// в”Ђв”Ђ Public API в”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђ

/// Apply a half-space cut to every solid in `brep`.
///
/// Iterates over every solid and over every face inside that solid:
/// each cylinder face is classified relative to `plane` and replaced with the
/// surviving partial geometry (partial cylinder + chord face for lateral cuts,
/// truncated cylinder + flat disk cap for axial cuts).
///
/// Solids whose ALL cylinder faces are entirely on the discarded side are
/// removed; solids whose ALL cylinder faces are entirely on the kept side are
/// left unchanged; mixed solids are rebuilt face-by-face.
///
/// Returns the number of solids remaining after the cut.
pub fn half_space_cut_brep(brep: &mut BRep, plane: &ClipPlane) -> usize {
    let solid_ids: Vec<SolidId> = brep.solids.keys().collect();
    let mut to_drop: Vec<SolidId> = Vec::new();
    let mut to_add: Vec<NewSolidParts> = Vec::new();

    for solid_id in solid_ids {
        match process_solid(brep, solid_id, plane) {
            SolidOutcome::Unchanged => {}
            SolidOutcome::Drop => to_drop.push(solid_id),
            SolidOutcome::Replace(parts) => {
                to_drop.push(solid_id);
                to_add.push(parts);
            }
        }
    }

    for id in to_drop {
        brep.solids.remove(id);
    }
    for parts in to_add {
        materialise_solid(brep, parts);
    }

    brep.solids.len()
}

// в”Ђв”Ђ Implementation в”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђ

enum SolidOutcome {
    Unchanged,
    Drop,
    Replace(NewSolidParts),
}

/// Rebuilt geometry to be inserted as a new solid.
struct NewSolidParts {
    faces: Vec<FaceTemplate>,
    /// Non-cylinder faces from the original solid that lie on the kept side of
    /// the cut plane (sphere caps, existing flat disks, torus fillets, вЂ¦).
    /// These are copied verbatim into the replacement solid so the mesh stays closed.
    copied_face_ids: Vec<FaceId>,
    name: Option<String>,
}

/// A face waiting to be materialised in the BRep.
enum FaceTemplate {
    /// Full cylinder (kept as-is, copied from the original solid).
    FullCylinder {
        cyl: CylSurf,
        length: f64,
        start: FaceBoundary,
        end: FaceBoundary,
    },
    /// Cylinder with the axis parallel to the cut plane вЂ” chord cut.
    PartialCylinder {
        axis_start: Point3,
        axis_end: Point3,
        radius: f64,
        arc_half_angle: f64,
        up: UnitVec3, // cut plane normal (into kept half-space)
    },
    /// Cylinder with the axis perpendicular to the cut plane вЂ” truncated.
    /// The cylinder runs from `axis_start` to `axis_end`; one of those points is
    /// the truncation point (exactly on the plane), the other is the original,
    /// untouched endpoint.
    ///
    /// `start`/`end` carry the boundary curves for each end.  The **cut** end
    /// gets a fresh plain `Circle` (it pairs with the new flat disk cap).  The
    /// **uncut** end keeps the *original* boundary (which may be a miter
    /// ellipse) so the shared junction edge with a kept neighbour face is
    /// preserved вЂ” without this the shell would open up at every truncated
    /// filament that still joins a surviving connector.
    AxialTruncated {
        axis_start: Point3,
        axis_end: Point3,
        radius: f64,
        start: FaceBoundary,
        end: FaceBoundary,
    },
    /// Polygonal flat face (e.g. chord face of a PartialCylinder).
    Polygon { plane: Plane3, points: Vec<Point3> },
    /// Flat disk (full circle, used for AxialTruncated end caps).
    Disk { plane: Plane3, radius: f64 },
    /// Partial disk (used for PartialCylinder end caps).
    PartialDisk {
        plane: Plane3,
        radius: f64,
        start_angle: f64,
        end_angle: f64,
    },
}

/// Process one solid: classify each face and build the replacement parts.
fn process_solid(brep: &BRep, solid_id: SolidId, plane: &ClipPlane) -> SolidOutcome {
    let solid = match brep.solids.get(solid_id) {
        Some(s) => s,
        None => return SolidOutcome::Unchanged,
    };

    // Collect face references.
    let mut cyl_faces: Vec<(CylSurf, f64, FaceBoundary, FaceBoundary)> = Vec::new();
    // Non-cylinder faces (sphere caps, existing flat disks, torus fillets) on the kept side.
    let mut kept_other_face_ids: Vec<FaceId> = Vec::new();
    // Template faces that genuinely cross the plane but have no cut form.
    let mut straddling = 0usize;
    for &shell_id in &solid.shells {
        let shell = match brep.shells.get(shell_id) {
            Some(s) => s,
            None => continue,
        };
        for &face_id in &shell.faces {
            let face = match brep.faces.get(face_id) {
                Some(f) => f,
                None => continue,
            };
            if let (FaceGeom::Cylinder(cyl), FaceExtent::Cylinder { length, start, end }) =
                (&face.geom, &face.extent)
            {
                cyl_faces.push((*cyl, *length, start.clone(), end.clone()));
            } else if let Some((d_min, d_max)) = face_signed_range(face, plane) {
                // Face template with a closed-form region: classify EXACTLY.
                if d_min >= -COINCIDENT_TOL {
                    kept_other_face_ids.push(face_id); // wholly on the kept side
                } else if d_max <= COINCIDENT_TOL {
                    // wholly on the discarded side в†’ drop
                } else {
                    // Genuinely straddles the plane and no template can express
                    // the cut region.  Keeping it whole leaves a little material
                    // past the plane; DROPPING it would open the shell, which is
                    // far worse вЂ” so keep and say so loudly.
                    straddling += 1;
                    kept_other_face_ids.push(face_id);
                }
            } else {
                // No closed-form region (explicitly trimmed loops / templates we
                // do not model): fall back to the carrier's representative point.
                // NB this is coarse вЂ” it is deliberately biased toward KEEPING,
                // because a face we cannot classify must never be silently
                // deleted; that is exactly how a shell opens up.
                let rep = match &face.geom {
                    FaceGeom::Sphere(s) => s.centre,
                    FaceGeom::Plane(p) => p.frame.origin,
                    FaceGeom::Torus(t) => t.frame.origin,
                    FaceGeom::Cylinder(c) => c.frame.origin,
                };
                if plane.normal.dot_vec(rep - plane.origin) >= -COINCIDENT_TOL {
                    kept_other_face_ids.push(face_id);
                }
            }
        }
    }

    if cyl_faces.is_empty() {
        return SolidOutcome::Unchanged;
    }

    // Classify each cylinder face.
    let mut new_faces: Vec<FaceTemplate> = Vec::new();
    let mut any_kept = false;
    let mut any_cut = false;
    let mut all_dropped = true;

    for (cyl, length, start, end) in cyl_faces {
        let outcome = classify_cylinder(&cyl, length, plane);
        match outcome {
            CylinderOutcome::EntirelyKept => {
                new_faces.push(FaceTemplate::FullCylinder {
                    cyl,
                    length,
                    start,
                    end,
                });
                any_kept = true;
                all_dropped = false;
            }
            CylinderOutcome::EntirelyDiscarded => {
                // Drop this cylinder.
            }
            CylinderOutcome::LateralCut {
                arc_half_angle,
                up,
                axis_start,
                axis_end,
                radius,
            } => {
                new_faces.push(FaceTemplate::PartialCylinder {
                    axis_start,
                    axis_end,
                    radius,
                    arc_half_angle,
                    up,
                });
                // Add chord face + 2 partial disk caps.
                add_lateral_cut_caps(
                    &mut new_faces,
                    axis_start,
                    axis_end,
                    radius,
                    arc_half_angle,
                    up,
                );
                any_cut = true;
                all_dropped = false;
            }
            CylinderOutcome::AxialCut {
                new_start,
                new_end,
                kept_end_at_plane,
                radius,
                axis_dir,
            } => {
                // The cut end gets a fresh plain circle (it pairs with the new
                // flat disk cap below).  The uncut end keeps its *original*
                // boundary curve so the shared junction edge with any surviving
                // neighbour face (e.g. a kept connector's miter ellipse) is
                // preserved вЂ” otherwise the shell opens up exactly at the
                // truncated filament ends.
                let cut_circle = |centre: Point3| {
                    FaceBoundary::Circle(cadcore_geom::Circle3::new(centre, axis_dir, radius))
                };
                let (start_bound, end_bound) = if kept_end_at_plane {
                    // axis_end is the cut end; axis_start is original/untouched.
                    (start.clone(), cut_circle(new_end))
                } else {
                    // axis_start is the cut end; axis_end is original/untouched.
                    (cut_circle(new_start), end.clone())
                };
                new_faces.push(FaceTemplate::AxialTruncated {
                    axis_start: new_start,
                    axis_end: new_end,
                    radius,
                    start: start_bound,
                    end: end_bound,
                });
                // Add flat disk cap at the cut end (replacing the original hemisphere).
                // Cap normal points OUTWARD = away from the kept solid:
                //   = +axis_dir if the cut was at the END,
                //   = -axis_dir if the cut was at the START.
                let cap_normal = if kept_end_at_plane {
                    // axis_end is the cut end в†’ outward = +axis_dir
                    axis_dir
                } else {
                    // axis_start is the cut end в†’ outward = -axis_dir
                    -axis_dir
                };
                let cap_centre = if kept_end_at_plane {
                    new_end
                } else {
                    new_start
                };
                new_faces.push(FaceTemplate::Disk {
                    plane: Plane3::from_origin_normal(cap_centre, cap_normal),
                    radius,
                });
                any_cut = true;
                all_dropped = false;
            }
        }
    }

    if all_dropped {
        return SolidOutcome::Drop;
    }
    if any_kept && !any_cut {
        return SolidOutcome::Unchanged;
    }

    if straddling > 0 {
        eprintln!(
            "[cadcore::boolean] {straddling} template face(s) cross the cut plane with no cut \
             form (kept whole; material extends past the plane) in solid {:?}",
            solid.name
        );
    }

    SolidOutcome::Replace(NewSolidParts {
        faces: new_faces,
        copied_face_ids: kept_other_face_ids,
        name: solid.name.clone(),
    })
}

// в”Ђв”Ђ Cylinder classification в”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђ

enum CylinderOutcome {
    EntirelyKept,
    EntirelyDiscarded,
    /// Cylinder axis is approximately parallel to the cut plane; centre-line
    /// runs along the plane, and the cylinder envelope partially intersects it.
    LateralCut {
        axis_start: Point3,
        axis_end: Point3,
        radius: f64,
        arc_half_angle: f64,
        up: UnitVec3,
    },
    /// Cylinder axis is approximately perpendicular to the cut plane; the axis
    /// crosses the plane.  The surviving cylinder runs from `new_start` to
    /// `new_end`; one of those points is on the plane (the truncation point).
    AxialCut {
        new_start: Point3,
        new_end: Point3,
        /// `true` when `new_end` is the truncation point (original `end` was discarded).
        /// `false` when `new_start` is the truncation point.
        kept_end_at_plane: bool,
        radius: f64,
        axis_dir: UnitVec3,
    },
}

const PARALLEL_TOL: f64 = 0.05; // ~3В° from parallel
const COINCIDENT_TOL: f64 = 1.0e-6;

fn classify_cylinder(cyl: &CylSurf, length: f64, plane: &ClipPlane) -> CylinderOutcome {
    let axis_dir = cyl.frame.z;
    let axis_start = cyl.frame.origin;
    let axis_end = axis_start + axis_dir.as_vec() * length;
    let radius = cyl.radius;

    let axis_dot_n = plane.normal.dot_vec(axis_dir.as_vec());

    // в”Ђв”Ђ LATERAL case (axis вЉҐ plane.normal, axis в€Ґ plane) в”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђ
    if axis_dot_n.abs() < PARALLEL_TOL {
        // Distance from any axis point to the plane (axis is parallel to plane).
        let dist = plane.normal.dot_vec(axis_start - plane.origin);

        if dist >= radius - COINCIDENT_TOL {
            return CylinderOutcome::EntirelyKept;
        }
        if dist <= -radius + COINCIDENT_TOL {
            return CylinderOutcome::EntirelyDiscarded;
        }

        // arc_half_angle = the half-span of the surviving arc, measured from
        // the "up" direction (= plane.normal, pointing into kept half-space).
        //   dist =  radius в†’ arc_half_angle = 0     (tangent from outside, nothing kept)
        //   dist =  0      в†’ arc_half_angle = ПЂ/2   (half-cylinder)
        //   dist = -radius в†’ arc_half_angle = ПЂ     (whole cylinder, but EntirelyDiscarded above)
        let cos_a = (-dist / radius).clamp(-1.0, 1.0);
        let arc_half_angle = cos_a.acos();

        return CylinderOutcome::LateralCut {
            axis_start,
            axis_end,
            radius,
            arc_half_angle,
            up: plane.normal,
        };
    }

    // в”Ђв”Ђ AXIAL case (axis в€Ґ plane.normal) в”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђ
    let d_start = plane.normal.dot_vec(axis_start - plane.origin);
    let d_end = plane.normal.dot_vec(axis_end - plane.origin);

    let s_kept = d_start >= -COINCIDENT_TOL;
    let e_kept = d_end >= -COINCIDENT_TOL;

    if s_kept && e_kept {
        return CylinderOutcome::EntirelyKept;
    }
    if !s_kept && !e_kept {
        return CylinderOutcome::EntirelyDiscarded;
    }

    // Crossing: compute the intersection.
    let denom = d_end - d_start;
    if denom.abs() < 1.0e-12 {
        return CylinderOutcome::EntirelyKept; // degenerate; bail
    }
    let t = -d_start / denom;
    let intersect = axis_start + (axis_end - axis_start) * t;

    if s_kept {
        // start kept, end discarded в†’ truncate at end
        CylinderOutcome::AxialCut {
            new_start: axis_start,
            new_end: intersect,
            kept_end_at_plane: true,
            radius,
            axis_dir,
        }
    } else {
        // start discarded, end kept в†’ truncate at start
        CylinderOutcome::AxialCut {
            new_start: intersect,
            new_end: axis_end,
            kept_end_at_plane: false,
            radius,
            axis_dir,
        }
    }
}

// в”Ђв”Ђ Exact face-region classification в”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђ
//
// A face that is not a full cylinder band still has to be classified against
// the plane, and a *representative point* is not good enough: a U-turn torus
// fillet's centre sits far from the fillet itself, and a chord face's plane
// origin is one of its corners.  These helpers give the EXACT signed-distance
// range of each template's material region, so the keep/drop decision is right
// and only a genuine straddle is reported as such.

/// Min/max of `a·cos О± + b·sin О±` over `О± в€€ [lo, hi]`.
///
/// `f(О±) = R·cos(О± в€’ П†)` with `R = hypot(a, b)`, `П† = atan2(b, a)`: the extrema
/// are the interval endpoints plus, when they fall inside the interval, the
/// stationary points `О± в‰Ў П†` (maximum `+R`) and `О± в‰Ў П† + ПЂ` (minimum `в€’R`).
fn harmonic_range(a: f64, b: f64, lo: f64, hi: f64) -> (f64, f64) {
    use std::f64::consts::{PI, TAU};
    let (lo, hi) = if lo <= hi { (lo, hi) } else { (hi, lo) };
    let r = a.hypot(b);
    if r < 1.0e-15 {
        return (0.0, 0.0);
    }
    if hi - lo >= TAU {
        return (-r, r);
    }
    let f = |ang: f64| a * ang.cos() + b * ang.sin();
    let (f_lo, f_hi) = (f(lo), f(hi));
    let mut min = f_lo.min(f_hi);
    let mut max = f_lo.max(f_hi);
    let phi = b.atan2(a);
    for (base, is_max) in [(phi, true), (phi + PI, false)] {
        // Smallest representative of `base` (mod 2ПЂ) that is в‰Ґ lo.
        let ang = base + TAU * ((lo - base) / TAU).ceil();
        if ang <= hi {
            if is_max {
                max = max.max(r);
            } else {
                min = min.min(-r);
            }
        }
    }
    (min, max)
}

/// Signed-distance range of a cylinder band, optionally limited to the arc
/// `[a0, a1]` measured from `ref_dir` about the axis (the chord-cut /
/// quarter-cylinder templates).
fn cylinder_band_range(
    cyl: &CylSurf,
    length: f64,
    arc: Option<(UnitVec3, f64, f64)>,
    plane: &ClipPlane,
) -> (f64, f64) {
    let n = plane.normal;
    let d0 = n.dot_vec(cyl.frame.origin - plane.origin);
    let na = n.dot_vec(cyl.axis().as_vec()) * length;
    let (ax_lo, ax_hi) = if na >= 0.0 { (0.0, na) } else { (na, 0.0) };

    let (r_lo, r_hi) = match arc {
        // p(О±) = axis_point + r·(cos О±·ref_dir + sin О±·(axis Г— ref_dir)) вЂ” the
        // basis `add_lateral_cut_caps` and the STEP writer both use.
        Some((ref_dir, a0, a1)) => {
            let right = cyl.axis().cross(ref_dir);
            harmonic_range(
                cyl.radius * n.dot_vec(ref_dir.as_vec()),
                cyl.radius * n.dot_vec(right),
                a0,
                a1,
            )
        }
        // Full cross-section: В± the radius projected into the plane вЉҐ axis.
        None => {
            let na_unit = n.dot_vec(cyl.axis().as_vec());
            let r = cyl.radius * (1.0 - na_unit * na_unit).max(0.0).sqrt();
            (-r, r)
        }
    };

    (d0 + ax_lo + r_lo, d0 + ax_hi + r_hi)
}

/// Signed-distance range `(min, max)` of the material region of a *template*
/// face against `plane`.
///
/// `None` means the extent carries no closed-form region (explicitly trimmed
/// loops, or a template we do not model here) вЂ” the caller falls back to a
/// representative point.
fn face_signed_range(face: &Face, plane: &ClipPlane) -> Option<(f64, f64)> {
    use std::f64::consts::TAU;
    let n = plane.normal;
    let dist = |p: Point3| n.dot_vec(p - plane.origin);

    // Range of a flat circular region (disk / circular boundary), optionally
    // limited to an arc.
    let planar_circle = |frame: &Frame3, radius: f64, arc: Option<(f64, f64)>| {
        let c = dist(frame.origin);
        let (a, b) = (
            radius * n.dot_vec(frame.x.as_vec()),
            radius * n.dot_vec(frame.y.as_vec()),
        );
        let (lo, hi) = match arc {
            Some((a0, a1)) => harmonic_range(a, b, a0, a1),
            None => harmonic_range(a, b, 0.0, TAU),
        };
        (c + lo, c + hi)
    };

    match (&face.geom, &face.extent) {
        (FaceGeom::Cylinder(cyl), FaceExtent::Cylinder { length, .. }) => {
            Some(cylinder_band_range(cyl, *length, None, plane))
        }
        (
            FaceGeom::Cylinder(cyl),
            FaceExtent::PartialCylinder {
                length,
                arc_start_angle,
                arc_end_angle,
                arc_ref_dir,
            }
            | FaceExtent::CylinderArcFace {
                length,
                arc_start_angle,
                arc_end_angle,
                arc_ref_dir,
            },
        ) => Some(cylinder_band_range(
            cyl,
            *length,
            Some((*arc_ref_dir, *arc_start_angle, *arc_end_angle)),
            plane,
        )),

        (FaceGeom::Plane(_), FaceExtent::Polygon { points }) if !points.is_empty() => {
            let mut lo = f64::MAX;
            let mut hi = f64::MIN;
            for p in points {
                let t = dist(*p);
                lo = lo.min(t);
                hi = hi.max(t);
            }
            Some((lo, hi))
        }
        (FaceGeom::Plane(pl), FaceExtent::Disk { radius }) => {
            Some(planar_circle(&pl.frame, *radius, None))
        }
        // A partial disk is the region between the arc and its chord вЂ” the
        // convex hull of the arc, so the arc's own range is exact.
        (
            FaceGeom::Plane(pl),
            FaceExtent::PartialDisk {
                radius,
                start_angle,
                end_angle,
            },
        ) => Some(planar_circle(
            &pl.frame,
            *radius,
            Some((*start_angle, *end_angle)),
        )),
        (
            FaceGeom::Plane(_),
            FaceExtent::PlanarBoundary {
                boundary: FaceBoundary::Circle(c),
            },
        ) => Some(planar_circle(&c.frame, c.radius, None)),

        // Torus fillet: every surface point is within `minor_radius` of the
        // centre-line circle, so the centre-line's range padded by the minor
        // radius is exact.  (The torus CENTRE вЂ” what this used to test вЂ” is not
        // even on the solid.)
        (FaceGeom::Torus(t), FaceExtent::TorusFillet { start_circle, end_circle }) => {
            let theta_of = |p: Point3| {
                let w = p - t.frame.origin;
                t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w))
            };
            let (lo, hi) = if (start_circle.frame.origin - end_circle.frame.origin).length() < 1e-9
            {
                (0.0, TAU) // full donut
            } else {
                let lo = theta_of(start_circle.frame.origin);
                let mut hi = theta_of(end_circle.frame.origin);
                hi += TAU * ((lo - hi) / TAU).round();
                (lo, hi)
            };
            let c = dist(t.frame.origin);
            let (a, b) = (
                t.major_radius * n.dot_vec(t.frame.x.as_vec()),
                t.major_radius * n.dot_vec(t.frame.y.as_vec()),
            );
            let (r_lo, r_hi) = harmonic_range(a, b, lo, hi);
            Some((c + r_lo - t.minor_radius, c + r_hi + t.minor_radius))
        }

        // A sphere cap is bounded by its whole sphere.
        (FaceGeom::Sphere(s), _) => {
            let c = dist(s.centre);
            Some((c - s.radius, c + s.radius))
        }

        _ => None,
    }
}

// в”Ђв”Ђ Chord-face / partial-disk cap construction в”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђ

fn add_lateral_cut_caps(
    out: &mut Vec<FaceTemplate>,
    axis_start: Point3,
    axis_end: Point3,
    radius: f64,
    arc_half_angle: f64,
    up: UnitVec3,
) {
    let axis_vec = axis_end - axis_start;
    let axis_dir = match UnitVec3::try_from_vec(axis_vec) {
        Some(u) => u,
        None => return,
    };

    let right = match UnitVec3::try_from_vec(axis_dir.cross(up)) {
        Some(u) => u,
        None => return,
    };

    let cos_a = arc_half_angle.cos();
    let sin_a = arc_half_angle.sin();

    // Arc endpoints in cross-section at axis_start:
    //   p = axis_start + radius*(cos(О±)*up + В±sin(О±)*right)
    let chord_offset = up.as_vec() * (radius * cos_a);
    let chord_p_neg_start = axis_start + chord_offset - right.as_vec() * (radius * sin_a);
    let chord_p_pos_start = axis_start + chord_offset + right.as_vec() * (radius * sin_a);
    let chord_p_neg_end = axis_end + chord_offset - right.as_vec() * (radius * sin_a);
    let chord_p_pos_end = axis_end + chord_offset + right.as_vec() * (radius * sin_a);

    // в”Ђв”Ђ Chord rectangle face (flat) в”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђ
    // Outward normal points into the DISCARDED half-space = в€’up.
    let chord_normal = -up;
    let chord_plane = Plane3::from_origin_normal(chord_p_neg_start, chord_normal);
    out.push(FaceTemplate::Polygon {
        plane: chord_plane,
        // CCW order viewed from -up direction:
        points: vec![
            chord_p_neg_start,
            chord_p_neg_end,
            chord_p_pos_end,
            chord_p_pos_start,
        ],
    });

    // в”Ђв”Ђ Two partial-disk end caps в”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђ
    let cap_start_frame = Frame3 {
        origin: axis_start,
        x: up,
        y: right, // must match facet_partial_cylinder's +right basis
        z: -axis_dir,
    };
    let cap_start_plane = Plane3 {
        frame: cap_start_frame,
    };

    let cap_end_frame = Frame3 {
        origin: axis_end,
        x: up,
        y: right,
        z: axis_dir,
    };
    let cap_end_plane = Plane3 {
        frame: cap_end_frame,
    };

    out.push(FaceTemplate::PartialDisk {
        plane: cap_start_plane,
        radius,
        start_angle: -arc_half_angle,
        end_angle: arc_half_angle,
    });
    out.push(FaceTemplate::PartialDisk {
        plane: cap_end_plane,
        radius,
        start_angle: -arc_half_angle,
        end_angle: arc_half_angle,
    });
}

// в”Ђв”Ђ Materialisation в”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђв”Ђ

fn materialise_solid(brep: &mut BRep, parts: NewSolidParts) {
    let mut face_ids = Vec::new();

    for tpl in parts.faces {
        let fid = match tpl {
            FaceTemplate::FullCylinder {
                cyl,
                length,
                start,
                end,
            } => {
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom: FaceGeom::Cylinder(cyl),
                    normal: FaceNormal::Same,
                    outer_loop: loop_id,
                    inner_loops: vec![],
                    shell: cadcore_topo::ShellId::default(),
                    extent: FaceExtent::Cylinder { length, start, end },
                })
            }
            FaceTemplate::PartialCylinder {
                axis_start,
                axis_end,
                radius,
                arc_half_angle,
                up,
            } => {
                let axis_vec = axis_end - axis_start;
                let length = axis_vec.length();
                let axis_dir = match UnitVec3::try_from_vec(axis_vec) {
                    Some(u) => u,
                    None => continue,
                };
                let cyl = CylSurf::new(axis_start, axis_dir, radius);
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom: FaceGeom::Cylinder(cyl),
                    normal: FaceNormal::Same,
                    outer_loop: loop_id,
                    inner_loops: vec![],
                    shell: cadcore_topo::ShellId::default(),
                    extent: FaceExtent::PartialCylinder {
                        length,
                        arc_start_angle: -arc_half_angle,
                        arc_end_angle: arc_half_angle,
                        arc_ref_dir: up,
                    },
                })
            }
            FaceTemplate::AxialTruncated {
                axis_start,
                axis_end,
                radius,
                start,
                end,
            } => {
                let axis_vec = axis_end - axis_start;
                let length = axis_vec.length();
                let axis_dir = match UnitVec3::try_from_vec(axis_vec) {
                    Some(u) => u,
                    None => continue,
                };
                let cyl = CylSurf::new(axis_start, axis_dir, radius);
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom: FaceGeom::Cylinder(cyl),
                    normal: FaceNormal::Same,
                    outer_loop: loop_id,
                    inner_loops: vec![],
                    shell: cadcore_topo::ShellId::default(),
                    extent: FaceExtent::Cylinder { length, start, end },
                })
            }
            FaceTemplate::Polygon { plane, points } => {
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom: FaceGeom::Plane(plane),
                    normal: FaceNormal::Same,
                    outer_loop: loop_id,
                    inner_loops: vec![],
                    shell: cadcore_topo::ShellId::default(),
                    extent: FaceExtent::Polygon { points },
                })
            }
            FaceTemplate::Disk { plane, radius } => {
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom: FaceGeom::Plane(plane),
                    normal: FaceNormal::Same,
                    outer_loop: loop_id,
                    inner_loops: vec![],
                    shell: cadcore_topo::ShellId::default(),
                    extent: FaceExtent::Disk { radius },
                })
            }
            FaceTemplate::PartialDisk {
                plane,
                radius,
                start_angle,
                end_angle,
            } => {
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom: FaceGeom::Plane(plane),
                    normal: FaceNormal::Same,
                    outer_loop: loop_id,
                    inner_loops: vec![],
                    shell: cadcore_topo::ShellId::default(),
                    extent: FaceExtent::PartialDisk {
                        radius,
                        start_angle,
                        end_angle,
                    },
                })
            }
        };
        face_ids.push(fid);
    }

    // Copy preserved non-cylinder faces (sphere caps, existing disks, fillets).
    for src_id in parts.copied_face_ids {
        if let Some(src) = brep.faces.get(src_id).cloned() {
            let loop_id = build_placeholder_loop(brep);
            let fid = brep.add_face(Face {
                geom: src.geom,
                normal: src.normal,
                outer_loop: loop_id,
                inner_loops: vec![],
                shell: cadcore_topo::ShellId::default(),
                extent: src.extent,
            });
            face_ids.push(fid);
        }
    }

    if face_ids.is_empty() {
        return;
    }

    let shell_id = brep.add_shell(Shell {
        faces: face_ids.clone(),
        is_outer: true,
        solid: cadcore_topo::SolidId::default(),
    });
    for &fid in &face_ids {
        if let Some(f) = brep.faces.get_mut(fid) {
            f.shell = shell_id;
        }
    }
    let solid_id = brep.add_solid(Solid {
        shells: vec![shell_id],
        name: parts.name,
    });
    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.solid = solid_id;
    }
}

fn build_placeholder_loop(brep: &mut BRep) -> cadcore_topo::LoopId {
    brep.add_loop(cadcore_topo::Loop {
        start: cadcore_topo::CoEdgeId::default(),
        face: cadcore_topo::FaceId::default(),
    })
}

// Faceted mesh Boolean union was removed from cadcore-ops.
// CADCore intentionally stays pure Rust and keeps exact analytic sweep surfaces here.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sweep::{
        sweep_circle_along_polyline, sweep_circle_along_rounded_polyline, SweepOptions,
    };

    /// Faces reachable from a LIVE solid.  `half_space_cut_brep` drops solids
    /// but leaves their faces in the arena, so counting `brep.faces` directly
    /// counts orphans too (and the STEP writer only walks solids).
    fn live_faces(brep: &BRep) -> Vec<&Face> {
        let mut out = Vec::new();
        for solid in brep.solids.values() {
            for &sh in &solid.shells {
                if let Some(shell) = brep.shells.get(sh) {
                    for &fid in &shell.faces {
                        if let Some(f) = brep.faces.get(fid) {
                            out.push(f);
                        }
                    }
                }
            }
        }
        out
    }

    fn count_live<F: Fn(&Face) -> bool>(brep: &BRep, pred: F) -> usize {
        live_faces(brep).into_iter().filter(|f| pred(f)).count()
    }

    /// X-direction cylinder centred at y=0.20 (= radius), cut by plane y=0.20.
    /// Axis runs through the plane в†’ LATERAL case, half-cylinder survives.
    #[test]
    fn x_cylinder_tangent_to_plane_survives_as_halfcyl() {
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(0.0, 0.20, 0.0), Point3::new(20.0, 0.20, 0.0)],
            0.20,
            &SweepOptions::default(),
        )
        .unwrap();
        let plane = ClipPlane {
            origin: Point3::new(0.0, 0.20, 0.0),
            normal: UnitVec3::Y, // keep y >= 0.20
        };
        let n = half_space_cut_brep(&mut brep, &plane);
        assert_eq!(n, 1, "cylinder must survive as a partial cylinder");

        let total_faces = brep.faces.len();
        assert!(
            total_faces >= 3,
            "partial-cyl solid needs >= 3 faces (cyl + chord + 2 caps), got {total_faces}"
        );

        // Must contain a PartialCylinder extent.
        let has_partial_cyl = brep
            .faces
            .values()
            .any(|f| matches!(f.extent, FaceExtent::PartialCylinder { .. }));
        assert!(has_partial_cyl, "no PartialCylinder face produced");

        // Must contain at least one Polygon (the chord face).
        let has_polygon = brep
            .faces
            .values()
            .any(|f| matches!(f.extent, FaceExtent::Polygon { .. }));
        assert!(has_polygon, "no chord polygon face produced");
    }

    /// Y-direction cylinder running from y=0 to y=20, cut by plane y=0.40.
    /// Axis в€Ґ plane.normal в†’ AXIAL case, truncated cylinder + flat disk cap.
    #[test]
    fn y_cylinder_perpendicular_to_plane_truncates() {
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(5.0, 0.0, 0.0), Point3::new(5.0, 20.0, 0.0)],
            0.20,
            &SweepOptions::default(),
        )
        .unwrap();
        let plane = ClipPlane {
            origin: Point3::new(0.0, 0.40, 0.0),
            normal: UnitVec3::Y, // keep y >= 0.40
        };
        let n = half_space_cut_brep(&mut brep, &plane);
        assert_eq!(n, 1, "cylinder must survive as truncated");

        // Flat disk cap should be present.
        let has_disk = brep
            .faces
            .values()
            .any(|f| matches!(f.extent, FaceExtent::Disk { .. }));
        assert!(has_disk, "no flat disk cap produced after axial cut");
    }

    /// Cylinder fully below the plane в†’ dropped.
    #[test]
    fn cylinder_below_plane_dropped() {
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(0.0, -5.0, 0.0), Point3::new(20.0, -5.0, 0.0)],
            0.20,
            &SweepOptions::default(),
        )
        .unwrap();
        let plane = ClipPlane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: UnitVec3::Y,
        };
        let n = half_space_cut_brep(&mut brep, &plane);
        assert_eq!(n, 0, "cylinder fully below plane must be dropped");
    }

    /// Two cut planes applied in sequence (the default `trim_side = "both"`).
    ///
    /// The chord (`PartialCylinder`) face produced by the FIRST plane used to be
    /// silently dropped by the second one вЂ” it is a `FaceGeom::Cylinder` face
    /// with a non-`Cylinder` extent, so it fell through both collection arms вЂ”
    /// while its chord polygon and partial-disk caps survived.  Result: a hole
    /// where the half-tube was, i.e. an open shell in the exported STEP.
    #[test]
    fn sequential_cuts_keep_earlier_chord_faces() {
        let mut brep = BRep::new();
        // U-shape: X-leg at y=1, Y-connector, X-leg at y=19.
        sweep_circle_along_polyline(
            &mut brep,
            &[
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(20.0, 1.0, 0.0),
                Point3::new(20.0, 19.0, 0.0),
                Point3::new(0.0, 19.0, 0.0),
            ],
            0.2,
            &SweepOptions::default(),
        )
        .unwrap();

        let is_partial_cyl = |f: &Face| matches!(f.extent, FaceExtent::PartialCylinder { .. });
        let is_polygon = |f: &Face| matches!(f.extent, FaceExtent::Polygon { .. });
        let is_partial_disk = |f: &Face| matches!(f.extent, FaceExtent::PartialDisk { .. });

        // Plane 1 chord-cuts the y = 1 leg (its axis lies IN the plane).
        half_space_cut_brep(
            &mut brep,
            &ClipPlane {
                origin: Point3::new(0.0, 1.0, 0.0),
                normal: UnitVec3::Y,
            },
        );
        assert_eq!(count_live(&brep, is_partial_cyl), 1, "plane 1 chord face");

        // Plane 2 chord-cuts the y = 19 leg.  BOTH chord faces must survive.
        half_space_cut_brep(
            &mut brep,
            &ClipPlane {
                origin: Point3::new(0.0, 19.0, 0.0),
                normal: -UnitVec3::Y,
            },
        );
        assert_eq!(
            count_live(&brep, is_polygon),
            2,
            "one chord polygon per plane"
        );
        assert_eq!(
            count_live(&brep, is_partial_disk),
            4,
            "two partial-disk caps per plane"
        );
        assert_eq!(
            count_live(&brep, is_partial_cyl),
            2,
            "plane 2 dropped the chord cylinder made by plane 1 \
             (its flat caps survived вЂ” the shell is open)"
        );
    }

    /// A U-turn torus fillet lying entirely on the discarded side must be
    /// dropped.  It used to be kept because the test point was the torus
    /// CENTRE, which for a U-turn sits on the kept side вЂ” leaving the fillet
    /// floating past the cut plane with an unshared junction edge.
    #[test]
    fn torus_fillet_past_the_plane_is_dropped() {
        let mut brep = BRep::new();
        // Hairpin in the XY plane: up to y = 20, U-turn, back down.
        sweep_circle_along_rounded_polyline(
            &mut brep,
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 20.0, 0.0),
                Point3::new(2.0, 20.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            ],
            0.2,
            1.0, // corner radius в†’ two torus fillets around y = 19..20
            &SweepOptions::default(),
        )
        .unwrap();
        let torus_count = |b: &BRep| count_live(b, |f| matches!(f.geom, FaceGeom::Torus(_)));
        assert!(torus_count(&brep) >= 2, "rounded corners must be tori");

        // Keep y <= 15: the whole U-turn is well past the plane.
        half_space_cut_brep(
            &mut brep,
            &ClipPlane {
                origin: Point3::new(0.0, 15.0, 0.0),
                normal: -UnitVec3::Y,
            },
        );
        assert_eq!(
            torus_count(&brep),
            0,
            "fillets entirely on the discarded side must be dropped"
        );
    }

    /// The mirror case: a fillet entirely on the kept side stays.
    #[test]
    fn torus_fillet_before_the_plane_is_kept() {
        let mut brep = BRep::new();
        sweep_circle_along_rounded_polyline(
            &mut brep,
            &[
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 20.0, 0.0),
                Point3::new(2.0, 20.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            ],
            0.2,
            1.0,
            &SweepOptions::default(),
        )
        .unwrap();
        let torus_count = |b: &BRep| count_live(b, |f| matches!(f.geom, FaceGeom::Torus(_)));
        let before = torus_count(&brep);

        // Keep y >= 5: the legs are truncated, the U-turn is untouched.
        half_space_cut_brep(
            &mut brep,
            &ClipPlane {
                origin: Point3::new(0.0, 5.0, 0.0),
                normal: UnitVec3::Y,
            },
        );
        assert_eq!(torus_count(&brep), before, "kept-side fillets must survive");
    }

    /// Cylinder fully above the plane в†’ kept unchanged.
    #[test]
    fn cylinder_above_plane_kept() {
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(0.0, 5.0, 0.0), Point3::new(20.0, 5.0, 0.0)],
            0.20,
            &SweepOptions::default(),
        )
        .unwrap();
        let face_count_before = brep.faces.len();
        let plane = ClipPlane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: UnitVec3::Y,
        };
        let n = half_space_cut_brep(&mut brep, &plane);
        assert_eq!(n, 1);
        assert_eq!(
            brep.faces.len(),
            face_count_before,
            "fully-above cylinder must not gain or lose faces"
        );
    }
}
