//! Solid Boolean half-space cut for cadcore B-Rep solids.
//!
//! # Algorithm
//!
//! For each solid in the [`BRep`]:
//!
//! 1.  **Determine spatial relationship** between the solid's bounding cylinder
//!     and the cut plane by inspecting the solid's faces.
//! 2.  **Entirely on kept side** → solid is unchanged.
//! 3.  **Entirely on discarded side** → solid is removed.
//! 4.  **Crossing case: axis perpendicular to plane normal**
//!     (cylinder axis ⊥ plane normal, i.e. the axis is parallel to the plane)
//!     → the cylinder is cut laterally — a "partial cylinder" (half-pipe) is
//!     created with a flat chord face at the cut plane.
//! 5.  **Crossing case: axis parallel to plane normal**
//!     (cylinder axis ∥ plane normal, i.e. the axis crosses the plane)
//!     → the cylinder is truncated at the plane; the end cap becomes a flat circle.
//!     This case is already handled correctly by `sweep_circle_along_path_with_caps`
//!     in the polyline-clip pipeline — it is kept unchanged here.
//!
//! # Why this is needed
//!
//! The existing `clip_polyline_with_radius` clips polyline *centre-lines* before
//! sweeping.  A filament whose centre-line is **parallel** to the cut plane and
//! whose centre is within one radius of the plane cannot be represented by
//! any clipped centre-line: the centre-line is either fully kept or fully dropped.
//!
//! `half_space_cut_brep` operates on the already-swept *solid* and correctly
//! handles the lateral case by building a new partial-cylinder solid.

use std::f64::consts::PI;

use cadcore_geom::{CylSurf, Plane3};
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{
    BRep, Face, FaceExtent, FaceGeom, FaceNormal, Shell, Solid, SolidId,
};

use crate::sweep::ClipPlane;

// ── Public API ────────────────────────────────────────────────────────────────

/// Apply a half-space cut to every solid in `brep`.
///
/// For each solid:
/// - If entirely on the **kept** side (`plane.normal · (centre − plane.origin) ≥ radius`):
///   the solid is left unchanged.
/// - If entirely on the **discarded** side: the solid (and all its topology) is dropped.
/// - If **crossing** (the plane intersects the solid's bounding cylinder laterally):
///   the solid is replaced with a partial-cylinder solid whose flat face lies on the plane.
///
/// Returns the number of solids remaining after the cut.
///
/// # Notes
///
/// This function works on the solid metadata encoded in [`FaceExtent`] variants.
/// Solids whose faces do not carry recognisable extent information are left as-is.
pub fn half_space_cut_brep(brep: &mut BRep, plane: &ClipPlane) -> usize {
    // Collect solid ids to avoid borrow conflicts.
    let solid_ids: Vec<SolidId> = brep.solids.keys().collect();
    let mut to_drop: Vec<SolidId>   = Vec::new();
    let mut to_add:  Vec<CutResult> = Vec::new();

    for solid_id in solid_ids {
        match classify_solid(brep, solid_id, plane) {
            SolidRelation::EntirelyKept => {
                // Keep unchanged.
            }
            SolidRelation::EntirelyDiscarded => {
                to_drop.push(solid_id);
            }
            SolidRelation::LateralCrossing { axis_start, axis_end, radius, dist } => {
                to_drop.push(solid_id);
                if let Some(result) = build_partial_cylinder(
                    axis_start, axis_end, radius, plane, dist,
                ) {
                    to_add.push(result);
                }
            }
            SolidRelation::Unknown => {
                // Leave as-is; don't know how to cut.
            }
        }
    }

    // Drop marked solids (shells + faces are orphaned — acceptable for STEP export).
    for id in to_drop {
        brep.solids.remove(id);
    }

    // Insert new partial-cylinder solids.
    for result in to_add {
        add_partial_cylinder_solid(brep, result);
    }

    brep.solids.len()
}

// ── Solid classification ──────────────────────────────────────────────────────

enum SolidRelation {
    EntirelyKept,
    EntirelyDiscarded,
    /// The solid is a cylinder whose axis is approximately PARALLEL to the cut
    /// plane (axis · plane_normal ≈ 0).  The cylinder centre is within the
    /// lateral intersection zone: |dist| < radius.
    LateralCrossing {
        axis_start: Point3,
        axis_end:   Point3,
        radius:     f64,
        /// Signed distance from cylinder centre to plane.
        /// Positive = centre is on the kept side.
        dist: f64,
    },
    Unknown,
}

/// Inspect the faces of `solid_id` to classify its spatial relationship with
/// `plane`.  Only cylindrical solids produced by `sweep_circle_along_path_with_caps`
/// are handled; others return [`SolidRelation::Unknown`].
fn classify_solid(brep: &BRep, solid_id: SolidId, plane: &ClipPlane) -> SolidRelation {
    let solid = match brep.solids.get(solid_id) {
        Some(s) => s,
        None    => return SolidRelation::Unknown,
    };

    // Gather cylinder faces to extract axis and radius.
    let mut cyl_axis_dir: Option<UnitVec3> = None;
    let mut cyl_origin:   Option<Point3>   = None;
    let mut cyl_length:   Option<f64>      = None;
    let mut cyl_radius:   Option<f64>      = None;

    for &shell_id in &solid.shells {
        let shell = match brep.shells.get(shell_id) {
            Some(s) => s,
            None    => continue,
        };
        for &face_id in &shell.faces {
            let face = match brep.faces.get(face_id) {
                Some(f) => f,
                None    => continue,
            };
            if let FaceGeom::Cylinder(cyl) = &face.geom {
                if cyl_axis_dir.is_none() {
                    cyl_axis_dir = Some(cyl.frame.z);
                    cyl_origin   = Some(cyl.frame.origin);
                    cyl_radius   = Some(cyl.radius);
                    if let FaceExtent::Cylinder { length, .. } = &face.extent {
                        cyl_length = Some(*length);
                    }
                }
            }
        }
    }

    let (axis_dir, origin, radius, length) = match (cyl_axis_dir, cyl_origin, cyl_radius, cyl_length) {
        (Some(d), Some(o), Some(r), Some(l)) => (d, o, r, l),
        _ => return SolidRelation::Unknown,
    };

    // axis_end is origin + axis_dir * length (CylSurf convention).
    let axis_start = origin;
    let axis_end   = origin + axis_dir.as_vec() * length;

    // Distance from cylinder centre-line to the cutting plane.
    // For a straight cylinder the centre-line is the axis; use its midpoint.
    let mid = axis_start + (axis_end - axis_start) * 0.5;
    let dist = plane.normal.dot_vec(mid - plane.origin);

    // Check if axis is parallel to the plane (axis · plane_normal ≈ 0).
    let axis_dot_normal = plane.normal.dot_vec(axis_dir.as_vec()).abs();
    let is_lateral      = axis_dot_normal < 0.05; // within ~3° of parallel

    if is_lateral {
        // Lateral (parallel) case — classify by cylinder envelope.
        if dist >= radius - 1.0e-7 {
            SolidRelation::EntirelyKept
        } else if dist <= -radius + 1.0e-7 {
            SolidRelation::EntirelyDiscarded
        } else {
            SolidRelation::LateralCrossing { axis_start, axis_end, radius, dist }
        }
    } else {
        // Axial (perpendicular) case.
        // Bounding AABB along the axis direction.
        let d_start = plane.normal.dot_vec(axis_start - plane.origin);
        let d_end   = plane.normal.dot_vec(axis_end   - plane.origin);
        // Physical extent in the plane-normal direction (add ±radius slop).
        let d_min = d_start.min(d_end) - radius;
        let d_max = d_start.max(d_end) + radius;

        if d_min >= -1.0e-7 {
            SolidRelation::EntirelyKept
        } else if d_max <= 1.0e-7 {
            SolidRelation::EntirelyDiscarded
        } else {
            // Let the existing polyline-clip pipeline handle this case
            // (it already produces the correct flat-cap solid).
            SolidRelation::Unknown
        }
    }
}

// ── Partial cylinder solid construction ───────────────────────────────────────

/// Intermediate representation of a partial-cylinder solid to be inserted.
struct CutResult {
    axis_start:      Point3,
    axis_end:        Point3,
    radius:          f64,
    arc_half_angle:  f64,   // π/2 = semicircle, π = full (shouldn't happen)
    plane_normal:    UnitVec3,
    /// "Up" direction in the cross-section: the direction toward the kept half.
    up:              UnitVec3,
}

/// Build a `CutResult` for a cylinder that is laterally crossed by `plane`.
///
/// `dist` = signed distance from cylinder centre to plane (positive = kept side).
fn build_partial_cylinder(
    axis_start: Point3,
    axis_end:   Point3,
    radius:     f64,
    plane:      &ClipPlane,
    dist:       f64,
) -> Option<CutResult> {
    // Arc half-angle of the surviving portion.
    // The kept arc spans where the cross-section is on the kept side of the plane.
    // If dist = 0: half of the circle → arc_half_angle = π/2.
    // If dist > 0: more than half → arc_half_angle > π/2.
    // If dist < 0: less than half → arc_half_angle < π/2.
    let cos_half = -(dist / radius).clamp(-1.0, 1.0); // negate: cos(π - θ) trick
    let arc_half_angle = cos_half.acos(); // in [0, π]

    if arc_half_angle < 1.0e-6 {
        // Tangent case: no surviving material.
        return None;
    }

    // "Up" direction is the plane normal pointing INTO the kept half-space.
    let up = plane.normal;

    Some(CutResult {
        axis_start,
        axis_end,
        radius,
        arc_half_angle,
        plane_normal: plane.normal,
        up,
    })
}

/// Materialise a `CutResult` as new faces/shells/solid in `brep`.
///
/// The solid consists of:
/// 1. **Partial cylinder** face — the arc of the CylSurf above the cut plane.
/// 2. **Chord face** — flat rectangle at the cut plane boundary.
/// 3. **Start arc cap** — partial disk at the axis start.
/// 4. **End arc cap** — partial disk at the axis end.
fn add_partial_cylinder_solid(brep: &mut BRep, r: CutResult) {
    let axis_vec = r.axis_end - r.axis_start;
    let length = axis_vec.length();
    if length < 1.0e-9 { return; }
    let axis_dir = match UnitVec3::try_from_vec(axis_vec) {
        Some(u) => u,
        None    => return,
    };

    let mut face_ids = Vec::new();

    // ── 1. Partial cylinder face ─────────────────────────────────────────────
    {
        let cyl = CylSurf::new(r.axis_start, axis_dir, r.radius);
        let loop_id = build_placeholder_loop(brep);
        let face_id = brep.add_face(Face {
            geom:        FaceGeom::Cylinder(cyl),
            normal:      FaceNormal::Same,
            outer_loop:  loop_id,
            inner_loops: vec![],
            shell:       cadcore_topo::ShellId::default(),
            extent:      FaceExtent::PartialCylinder {
                length,
                arc_start_angle: PI * 0.5 - r.arc_half_angle,
                arc_end_angle:   PI * 0.5 + r.arc_half_angle,
                arc_ref_dir:     r.up,
            },
        });
        face_ids.push(face_id);
    }

    // ── 2. Chord face (flat rectangle at cut plane) ──────────────────────────
    //
    // The chord rectangle is built from the arc endpoints at each end of the
    // cylinder.  The arc spans angle = ±arc_half_angle from the "up" direction
    // (= plane normal, pointing into the kept half-space) in the cylinder
    // cross-section plane.
    //
    // Local frame in the cross-section:
    //   up    = r.up           (toward kept half-space)
    //   right = axis_dir × up  (right-handed)
    {
        let right = match UnitVec3::try_from_vec(axis_dir.cross(r.up)) {
            Some(u) => u,
            None    => return,
        };
        let alpha = r.arc_half_angle;
        // The arc endpoints at axis_start are:
        //   p_left  = axis_start + up * cos(α) * r + right * sin(α) * r ... no
        // Wait: angle measured from UP (not from right). At angle=0 → up direction.
        // At angle=+α → rotated CCW by α from up in the (up, right) plane:
        //   p = axis_start + r*(sin(α)*right + cos(α)*up)   ... no, that's angle from up toward right
        // The chord connects the two arc endpoints:
        //   p1 = axis_start + r * (sin(-α) * right + cos(-α) * up)  (one side)
        //      = axis_start + r * (-sin(α) * right + cos(α) * up)
        //   p2 = axis_start + r * (sin(+α) * right + cos(α) * up)
        //
        // But cos(α) * r is the distance from the chord to the cylinder centre in the "up" direction.
        // And the chord is at: axis_start + r*cos(α)*up + (±r*sin(α))*right
        // The chord face (flat rectangle) is at y = axis_start.y + r*cos(α)  [for y-plane]
        // which matches: the plane is at this position when we built arc_half_angle.
        let cos_a = alpha.cos();
        let sin_a = alpha.sin();
        let chord_offset = r.up.as_vec() * (r.radius * cos_a); // shift from axis toward plane

        let chord_p1s = r.axis_start + chord_offset - right.as_vec() * (r.radius * sin_a);
        let chord_p2s = r.axis_start + chord_offset + right.as_vec() * (r.radius * sin_a);
        let chord_p1e = r.axis_end   + chord_offset - right.as_vec() * (r.radius * sin_a);
        let chord_p2e = r.axis_end   + chord_offset + right.as_vec() * (r.radius * sin_a);

        // Chord face normal points in the -up direction (into the discarded side):
        // But we want outward normal from the solid = -up (solid is above the plane).
        let chord_normal = -r.plane_normal;
        let chord_plane = Plane3::from_origin_normal(chord_p1s, chord_normal);
        let loop_id = build_placeholder_loop(brep);
        let face_id = brep.add_face(Face {
            geom:        FaceGeom::Plane(chord_plane),
            normal:      FaceNormal::Same,
            outer_loop:  loop_id,
            inner_loops: vec![],
            shell:       cadcore_topo::ShellId::default(),
            // CCW order when viewed from outside (from discarded side, i.e. -up direction):
            extent:      FaceExtent::Polygon {
                points: vec![chord_p1s, chord_p1e, chord_p2e, chord_p2s],
            },
        });
        face_ids.push(face_id);
    }

    // ── 3. Start arc end cap ─────────────────────────────────────────────────
    {
        let cap_normal = -axis_dir; // outward = toward negative axis
        let cap_plane  = Plane3::from_origin_normal(r.axis_start, cap_normal);
        let loop_id    = build_placeholder_loop(brep);
        let face_id    = brep.add_face(Face {
            geom:        FaceGeom::Plane(cap_plane),
            normal:      FaceNormal::Same,
            outer_loop:  loop_id,
            inner_loops: vec![],
            shell:       cadcore_topo::ShellId::default(),
            extent:      FaceExtent::PartialDisk {
                radius:      r.radius,
                start_angle: PI * 0.5 - r.arc_half_angle,
                end_angle:   PI * 0.5 + r.arc_half_angle,
            },
        });
        face_ids.push(face_id);
    }

    // ── 4. End arc end cap ───────────────────────────────────────────────────
    {
        let cap_normal = axis_dir;
        let cap_plane  = Plane3::from_origin_normal(r.axis_end, cap_normal);
        let loop_id    = build_placeholder_loop(brep);
        let face_id    = brep.add_face(Face {
            geom:        FaceGeom::Plane(cap_plane),
            normal:      FaceNormal::Same,
            outer_loop:  loop_id,
            inner_loops: vec![],
            shell:       cadcore_topo::ShellId::default(),
            extent:      FaceExtent::PartialDisk {
                radius:      r.radius,
                start_angle: PI * 0.5 - r.arc_half_angle,
                end_angle:   PI * 0.5 + r.arc_half_angle,
            },
        });
        face_ids.push(face_id);
    }

    // ── Assemble solid ────────────────────────────────────────────────────────
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
        name:   Some("partial_cylinder_cut".to_string()),
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
