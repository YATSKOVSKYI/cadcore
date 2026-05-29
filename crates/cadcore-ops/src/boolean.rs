//! Solid Boolean half-space cut for cadcore B-Rep solids.
//!
//! # Strategy
//!
//! Each *cylinder face* in the [`BRep`] is classified relative to the cutting
//! plane and rebuilt independently.  This works because each cylinder face in
//! the existing sweep pipeline corresponds to one straight segment of the
//! polyline — solids built segment-by-segment can be trimmed face-by-face.
//!
//! Three cases per cylinder:
//!
//! | axis ⊥ plane.normal? | dist(axis, plane) | result                              |
//! |----------------------|-------------------|-------------------------------------|
//! | yes (LATERAL)        | > radius          | kept                                |
//! | yes (LATERAL)        | < −radius         | dropped                             |
//! | yes (LATERAL)        | |dist| < radius   | partial cylinder + flat chord       |
//! | no  (AXIAL)          | both ends kept    | kept                                |
//! | no  (AXIAL)          | both ends discard | dropped                             |
//! | no  (AXIAL)          | crosses plane     | truncated + flat disk cap           |
//!
//! The same cylinder may be cut by multiple planes; planes are processed
//! sequentially.

use cadcore_geom::{CylSurf, Plane3};
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{
    BRep, Face, FaceBoundary, FaceExtent, FaceGeom, FaceNormal, Shell, Solid, SolidId,
};

use crate::sweep::ClipPlane;

// ── Public API ────────────────────────────────────────────────────────────────

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

// ── Implementation ────────────────────────────────────────────────────────────

enum SolidOutcome {
    Unchanged,
    Drop,
    Replace(NewSolidParts),
}

/// Rebuilt geometry to be inserted as a new solid.
struct NewSolidParts {
    faces: Vec<FaceTemplate>,
    name:  Option<String>,
}

/// A face waiting to be materialised in the BRep.
enum FaceTemplate {
    /// Full cylinder (kept as-is, copied from the original solid).
    FullCylinder {
        cyl:    CylSurf,
        length: f64,
        start:  FaceBoundary,
        end:    FaceBoundary,
    },
    /// Cylinder with the axis parallel to the cut plane — chord cut.
    PartialCylinder {
        axis_start:     Point3,
        axis_end:       Point3,
        radius:         f64,
        arc_half_angle: f64,
        up:             UnitVec3,   // cut plane normal (into kept half-space)
    },
    /// Cylinder with the axis perpendicular to the cut plane — truncated.
    /// The cylinder runs from `axis_start` to `axis_end`; `axis_end` was
    /// originally past the plane and is now exactly on it.
    AxialTruncated {
        axis_start: Point3,
        axis_end:   Point3,
        radius:     f64,
    },
    /// Polygonal flat face (e.g. chord face of a PartialCylinder).
    Polygon {
        plane:  Plane3,
        points: Vec<Point3>,
    },
    /// Flat disk (full circle, used for AxialTruncated end caps).
    Disk {
        plane:  Plane3,
        radius: f64,
    },
    /// Partial disk (used for PartialCylinder end caps).
    PartialDisk {
        plane:       Plane3,
        radius:      f64,
        start_angle: f64,
        end_angle:   f64,
    },
}

/// Process one solid: classify each face and build the replacement parts.
fn process_solid(brep: &BRep, solid_id: SolidId, plane: &ClipPlane) -> SolidOutcome {
    let solid = match brep.solids.get(solid_id) {
        Some(s) => s,
        None    => return SolidOutcome::Unchanged,
    };

    // Collect face references.
    let mut cyl_faces: Vec<(CylSurf, f64, FaceBoundary, FaceBoundary)> = Vec::new();
    for &shell_id in &solid.shells {
        let shell = match brep.shells.get(shell_id) { Some(s) => s, None => continue };
        for &face_id in &shell.faces {
            let face = match brep.faces.get(face_id) { Some(f) => f, None => continue };
            if let (FaceGeom::Cylinder(cyl), FaceExtent::Cylinder { length, start, end }) =
                (&face.geom, &face.extent)
            {
                cyl_faces.push((*cyl, *length, start.clone(), end.clone()));
            }
        }
    }

    if cyl_faces.is_empty() {
        return SolidOutcome::Unchanged;
    }

    // Classify each cylinder face.
    let mut new_faces: Vec<FaceTemplate> = Vec::new();
    let mut any_kept   = false;
    let mut any_cut    = false;
    let mut all_dropped = true;

    for (cyl, length, start, end) in cyl_faces {
        let outcome = classify_cylinder(&cyl, length, plane);
        match outcome {
            CylinderOutcome::EntirelyKept => {
                new_faces.push(FaceTemplate::FullCylinder { cyl, length, start, end });
                any_kept    = true;
                all_dropped = false;
            }
            CylinderOutcome::EntirelyDiscarded => {
                // Drop this cylinder.
            }
            CylinderOutcome::LateralCut { arc_half_angle, up, axis_start, axis_end, radius } => {
                new_faces.push(FaceTemplate::PartialCylinder {
                    axis_start, axis_end, radius, arc_half_angle, up,
                });
                // Add chord face + 2 partial disk caps.
                add_lateral_cut_caps(&mut new_faces, axis_start, axis_end, radius, arc_half_angle, up);
                any_cut     = true;
                all_dropped = false;
            }
            CylinderOutcome::AxialCut { new_start, new_end, kept_end_at_plane, radius, axis_dir } => {
                new_faces.push(FaceTemplate::AxialTruncated {
                    axis_start: new_start,
                    axis_end:   new_end,
                    radius,
                });
                // Add flat disk cap at the cut end (replacing the original hemisphere).
                // Cap normal points OUTWARD = away from the kept solid:
                //   = +axis_dir if the cut was at the END,
                //   = -axis_dir if the cut was at the START.
                let cap_normal = if kept_end_at_plane {
                    // axis_end is the cut end → outward = +axis_dir
                    axis_dir
                } else {
                    // axis_start is the cut end → outward = -axis_dir
                    -axis_dir
                };
                let cap_centre = if kept_end_at_plane { new_end } else { new_start };
                new_faces.push(FaceTemplate::Disk {
                    plane:  Plane3::from_origin_normal(cap_centre, cap_normal),
                    radius,
                });
                any_cut     = true;
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

    SolidOutcome::Replace(NewSolidParts {
        faces: new_faces,
        name:  solid.name.clone(),
    })
}

// ── Cylinder classification ───────────────────────────────────────────────────

enum CylinderOutcome {
    EntirelyKept,
    EntirelyDiscarded,
    /// Cylinder axis is approximately parallel to the cut plane; centre-line
    /// runs along the plane, and the cylinder envelope partially intersects it.
    LateralCut {
        axis_start:     Point3,
        axis_end:       Point3,
        radius:         f64,
        arc_half_angle: f64,
        up:             UnitVec3,
    },
    /// Cylinder axis is approximately perpendicular to the cut plane; the axis
    /// crosses the plane.  The surviving cylinder runs from `new_start` to
    /// `new_end`; one of those points is on the plane (the truncation point).
    AxialCut {
        new_start:           Point3,
        new_end:             Point3,
        /// `true` when `new_end` is the truncation point (original `end` was discarded).
        /// `false` when `new_start` is the truncation point.
        kept_end_at_plane:   bool,
        radius:              f64,
        axis_dir:            UnitVec3,
    },
}

const PARALLEL_TOL: f64 = 0.05;        // ~3° from parallel
const COINCIDENT_TOL: f64 = 1.0e-6;

fn classify_cylinder(cyl: &CylSurf, length: f64, plane: &ClipPlane) -> CylinderOutcome {
    let axis_dir   = cyl.frame.z;
    let axis_start = cyl.frame.origin;
    let axis_end   = axis_start + axis_dir.as_vec() * length;
    let radius     = cyl.radius;

    let axis_dot_n = plane.normal.dot_vec(axis_dir.as_vec());

    // ── LATERAL case (axis ⊥ plane.normal, axis ∥ plane) ──────────────────────
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
        //   dist =  radius → arc_half_angle = 0     (tangent from outside, nothing kept)
        //   dist =  0      → arc_half_angle = π/2   (half-cylinder)
        //   dist = -radius → arc_half_angle = π     (whole cylinder, but EntirelyDiscarded above)
        let cos_a = (-dist / radius).clamp(-1.0, 1.0);
        let arc_half_angle = cos_a.acos();

        return CylinderOutcome::LateralCut {
            axis_start, axis_end, radius,
            arc_half_angle,
            up: plane.normal,
        };
    }

    // ── AXIAL case (axis ∥ plane.normal) ───────────────────────────────────────
    let d_start = plane.normal.dot_vec(axis_start - plane.origin);
    let d_end   = plane.normal.dot_vec(axis_end   - plane.origin);

    let s_kept = d_start >= -COINCIDENT_TOL;
    let e_kept = d_end   >= -COINCIDENT_TOL;

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
        // start kept, end discarded → truncate at end
        CylinderOutcome::AxialCut {
            new_start: axis_start,
            new_end:   intersect,
            kept_end_at_plane: true,
            radius,
            axis_dir,
        }
    } else {
        // start discarded, end kept → truncate at start
        CylinderOutcome::AxialCut {
            new_start: intersect,
            new_end:   axis_end,
            kept_end_at_plane: false,
            radius,
            axis_dir,
        }
    }
}

// ── Chord-face / partial-disk cap construction ────────────────────────────────

fn add_lateral_cut_caps(
    out:            &mut Vec<FaceTemplate>,
    axis_start:     Point3,
    axis_end:       Point3,
    radius:         f64,
    arc_half_angle: f64,
    up:             UnitVec3,
) {
    let axis_vec = axis_end - axis_start;
    let axis_dir = match UnitVec3::try_from_vec(axis_vec) {
        Some(u) => u,
        None    => return,
    };

    let right = match UnitVec3::try_from_vec(axis_dir.cross(up)) {
        Some(u) => u,
        None    => return,
    };

    let cos_a = arc_half_angle.cos();
    let sin_a = arc_half_angle.sin();

    // Arc endpoints in cross-section at axis_start:
    //   p = axis_start + radius*(cos(α)*up + ±sin(α)*right)
    let chord_offset = up.as_vec() * (radius * cos_a);
    let chord_p_neg_start = axis_start + chord_offset - right.as_vec() * (radius * sin_a);
    let chord_p_pos_start = axis_start + chord_offset + right.as_vec() * (radius * sin_a);
    let chord_p_neg_end   = axis_end   + chord_offset - right.as_vec() * (radius * sin_a);
    let chord_p_pos_end   = axis_end   + chord_offset + right.as_vec() * (radius * sin_a);

    // ── Chord rectangle face (flat) ──────────────────────────────────────────
    // Outward normal points into the DISCARDED half-space = −up.
    let chord_normal = -up;
    let chord_plane  = Plane3::from_origin_normal(chord_p_neg_start, chord_normal);
    out.push(FaceTemplate::Polygon {
        plane:  chord_plane,
        // CCW order viewed from -up direction:
        points: vec![chord_p_neg_start, chord_p_neg_end, chord_p_pos_end, chord_p_pos_start],
    });

    // ── Two partial-disk end caps ────────────────────────────────────────────
    use std::f64::consts::PI;
    let cap_start_plane = Plane3::from_origin_normal(axis_start, -axis_dir);
    let cap_end_plane   = Plane3::from_origin_normal(axis_end,    axis_dir);
    // Arc angles are measured from "up" (plane normal); for the disk cap we
    // need angles measured from the plane's frame x-axis.  As an approximation
    // we use 0..arc_full_angle and rely on the writer's frame orientation.
    let arc_full = 2.0 * arc_half_angle;
    out.push(FaceTemplate::PartialDisk {
        plane:       cap_start_plane,
        radius,
        start_angle: PI * 0.5 - arc_half_angle,
        end_angle:   PI * 0.5 + arc_half_angle,
    });
    out.push(FaceTemplate::PartialDisk {
        plane:       cap_end_plane,
        radius,
        start_angle: PI * 0.5 - arc_half_angle,
        end_angle:   PI * 0.5 + arc_half_angle,
    });
    let _ = arc_full;
}

// ── Materialisation ───────────────────────────────────────────────────────────

fn materialise_solid(brep: &mut BRep, parts: NewSolidParts) {
    let mut face_ids = Vec::new();

    for tpl in parts.faces {
        let fid = match tpl {
            FaceTemplate::FullCylinder { cyl, length, start, end } => {
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom:        FaceGeom::Cylinder(cyl),
                    normal:      FaceNormal::Same,
                    outer_loop:  loop_id,
                    inner_loops: vec![],
                    shell:       cadcore_topo::ShellId::default(),
                    extent:      FaceExtent::Cylinder { length, start, end },
                })
            }
            FaceTemplate::PartialCylinder { axis_start, axis_end, radius, arc_half_angle, up } => {
                use std::f64::consts::PI;
                let axis_vec = axis_end - axis_start;
                let length = axis_vec.length();
                let axis_dir = match UnitVec3::try_from_vec(axis_vec) {
                    Some(u) => u,
                    None    => continue,
                };
                let cyl = CylSurf::new(axis_start, axis_dir, radius);
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom:        FaceGeom::Cylinder(cyl),
                    normal:      FaceNormal::Same,
                    outer_loop:  loop_id,
                    inner_loops: vec![],
                    shell:       cadcore_topo::ShellId::default(),
                    extent:      FaceExtent::PartialCylinder {
                        length,
                        arc_start_angle: PI * 0.5 - arc_half_angle,
                        arc_end_angle:   PI * 0.5 + arc_half_angle,
                        arc_ref_dir:     up,
                    },
                })
            }
            FaceTemplate::AxialTruncated { axis_start, axis_end, radius } => {
                let axis_vec = axis_end - axis_start;
                let length = axis_vec.length();
                let axis_dir = match UnitVec3::try_from_vec(axis_vec) {
                    Some(u) => u,
                    None    => continue,
                };
                let cyl = CylSurf::new(axis_start, axis_dir, radius);
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom:        FaceGeom::Cylinder(cyl),
                    normal:      FaceNormal::Same,
                    outer_loop:  loop_id,
                    inner_loops: vec![],
                    shell:       cadcore_topo::ShellId::default(),
                    extent:      FaceExtent::Cylinder {
                        length,
                        start: FaceBoundary::Circle(cadcore_geom::Circle3::new(axis_start, axis_dir, radius)),
                        end:   FaceBoundary::Circle(cadcore_geom::Circle3::new(axis_end,   axis_dir, radius)),
                    },
                })
            }
            FaceTemplate::Polygon { plane, points } => {
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom:        FaceGeom::Plane(plane),
                    normal:      FaceNormal::Same,
                    outer_loop:  loop_id,
                    inner_loops: vec![],
                    shell:       cadcore_topo::ShellId::default(),
                    extent:      FaceExtent::Polygon { points },
                })
            }
            FaceTemplate::Disk { plane, radius } => {
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom:        FaceGeom::Plane(plane),
                    normal:      FaceNormal::Same,
                    outer_loop:  loop_id,
                    inner_loops: vec![],
                    shell:       cadcore_topo::ShellId::default(),
                    extent:      FaceExtent::Disk { radius },
                })
            }
            FaceTemplate::PartialDisk { plane, radius, start_angle, end_angle } => {
                let loop_id = build_placeholder_loop(brep);
                brep.add_face(Face {
                    geom:        FaceGeom::Plane(plane),
                    normal:      FaceNormal::Same,
                    outer_loop:  loop_id,
                    inner_loops: vec![],
                    shell:       cadcore_topo::ShellId::default(),
                    extent:      FaceExtent::PartialDisk { radius, start_angle, end_angle },
                })
            }
        };
        face_ids.push(fid);
    }

    if face_ids.is_empty() {
        return;
    }

    let shell_id = brep.add_shell(Shell {
        faces:    face_ids.clone(),
        is_outer: true,
        solid:    cadcore_topo::SolidId::default(),
    });
    for &fid in &face_ids {
        if let Some(f) = brep.faces.get_mut(fid) {
            f.shell = shell_id;
        }
    }
    let solid_id = brep.add_solid(Solid {
        shells: vec![shell_id],
        name:   parts.name,
    });
    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.solid = solid_id;
    }
}

fn build_placeholder_loop(brep: &mut BRep) -> cadcore_topo::LoopId {
    brep.add_loop(cadcore_topo::Loop {
        start: cadcore_topo::CoEdgeId::default(),
        face:  cadcore_topo::FaceId::default(),
    })
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sweep::{sweep_circle_along_polyline, SweepOptions};

    /// X-direction cylinder centred at y=0.20 (= radius), cut by plane y=0.20.
    /// Axis runs through the plane → LATERAL case, half-cylinder survives.
    #[test]
    fn x_cylinder_tangent_to_plane_survives_as_halfcyl() {
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(0.0, 0.20, 0.0), Point3::new(20.0, 0.20, 0.0)],
            0.20,
            &SweepOptions::default(),
        ).unwrap();
        let plane = ClipPlane {
            origin: Point3::new(0.0, 0.20, 0.0),
            normal: UnitVec3::Y,  // keep y >= 0.20
        };
        let n = half_space_cut_brep(&mut brep, &plane);
        assert_eq!(n, 1, "cylinder must survive as a partial cylinder");

        let total_faces = brep.faces.len();
        assert!(total_faces >= 3,
            "partial-cyl solid needs >= 3 faces (cyl + chord + 2 caps), got {total_faces}");

        // Must contain a PartialCylinder extent.
        let has_partial_cyl = brep.faces.values()
            .any(|f| matches!(f.extent, FaceExtent::PartialCylinder { .. }));
        assert!(has_partial_cyl, "no PartialCylinder face produced");

        // Must contain at least one Polygon (the chord face).
        let has_polygon = brep.faces.values()
            .any(|f| matches!(f.extent, FaceExtent::Polygon { .. }));
        assert!(has_polygon, "no chord polygon face produced");
    }

    /// Y-direction cylinder running from y=0 to y=20, cut by plane y=0.40.
    /// Axis ∥ plane.normal → AXIAL case, truncated cylinder + flat disk cap.
    #[test]
    fn y_cylinder_perpendicular_to_plane_truncates() {
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(5.0, 0.0, 0.0), Point3::new(5.0, 20.0, 0.0)],
            0.20,
            &SweepOptions::default(),
        ).unwrap();
        let plane = ClipPlane {
            origin: Point3::new(0.0, 0.40, 0.0),
            normal: UnitVec3::Y,  // keep y >= 0.40
        };
        let n = half_space_cut_brep(&mut brep, &plane);
        assert_eq!(n, 1, "cylinder must survive as truncated");

        // Flat disk cap should be present.
        let has_disk = brep.faces.values()
            .any(|f| matches!(f.extent, FaceExtent::Disk { .. }));
        assert!(has_disk, "no flat disk cap produced after axial cut");
    }

    /// Cylinder fully below the plane → dropped.
    #[test]
    fn cylinder_below_plane_dropped() {
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(0.0, -5.0, 0.0), Point3::new(20.0, -5.0, 0.0)],
            0.20,
            &SweepOptions::default(),
        ).unwrap();
        let plane = ClipPlane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: UnitVec3::Y,
        };
        let n = half_space_cut_brep(&mut brep, &plane);
        assert_eq!(n, 0, "cylinder fully below plane must be dropped");
    }

    /// Cylinder fully above the plane → kept unchanged.
    #[test]
    fn cylinder_above_plane_kept() {
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(0.0, 5.0, 0.0), Point3::new(20.0, 5.0, 0.0)],
            0.20,
            &SweepOptions::default(),
        ).unwrap();
        let face_count_before = brep.faces.len();
        let plane = ClipPlane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: UnitVec3::Y,
        };
        let n = half_space_cut_brep(&mut brep, &plane);
        assert_eq!(n, 1);
        assert_eq!(brep.faces.len(), face_count_before,
            "fully-above cylinder must not gain or lose faces");
    }
}
