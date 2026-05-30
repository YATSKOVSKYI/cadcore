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

use std::collections::HashMap;

use rayon::prelude::*;

use cadcore_geom::{Circle3, CylSurf, Plane3, SphereSurf, TorusSurf};
use cadcore_math::{Frame3, Point3, UnitVec3, Vec3, PI, TAU};
use cadcore_topo::{
    BRep, Face, FaceBoundary, FaceExtent, FaceGeom, FaceId, FaceNormal, Shell, ShellId, Solid,
    SolidId,
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
    let cap_start_frame = Frame3 {
        origin: axis_start,
        x:      up,
        y:      -right,
        z:      -axis_dir,
    };
    let cap_start_plane = Plane3 { frame: cap_start_frame };

    let cap_end_frame = Frame3 {
        origin: axis_end,
        x:      up,
        y:      right,
        z:      axis_dir,
    };
    let cap_end_plane = Plane3 { frame: cap_end_frame };

    out.push(FaceTemplate::PartialDisk {
        plane:       cap_start_plane,
        radius,
        start_angle: -arc_half_angle,
        end_angle:   arc_half_angle,
    });
    out.push(FaceTemplate::PartialDisk {
        plane:       cap_end_plane,
        radius,
        start_angle: -arc_half_angle,
        end_angle:   arc_half_angle,
    });
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
                        arc_start_angle: -arc_half_angle,
                        arc_end_angle:   arc_half_angle,
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

// ════════════════════════════════════════════════════════════════════════════
//  Boolean UNION (faceted / watertight)
// ════════════════════════════════════════════════════════════════════════════
//
// A general analytic Boolean union is not expressible in this kernel: the STEP
// writer derives face bounds from `FaceExtent` (not the half-edge topology), and
// the edge vocabulary is limited to Line / Circle / Ellipse — which cannot carry
// the quartic intersection curve of two unequal cylinders.  We therefore take the
// robust faceted route:
//
//   1. Tessellate each input solid into a closed, outward-oriented triangle mesh.
//   2. Combine the two meshes with a BSP CSG union (the classic `csg.js`
//      algorithm) — splits faces along the partner solid and keeps the outside.
//   3. Weld coincident vertices (tol = 1e-7) and heal T-junctions so the result
//      is a true closed 2-manifold (every edge shared by exactly two triangles
//      with opposite orientation).
//   4. Materialise the result as `FaceExtent::Polygon` faces, which the existing
//      STEP writer emits as a watertight `MANIFOLD_SOLID_BREP`.
//
// The result is faceted, not analytic, but it is genuinely watertight for the
// general case (including unequal-radius perpendicular cylinders).

/// Vertex-coincidence tolerance (mm-scale geometry: 10 nm).
const WELD_TOL: f64 = 1.0e-5;
/// BSP plane-classification epsilon (matched to weld tolerance).
const BSP_EPS: f64 = 1.0e-5;
/// Parametric edge-interior threshold for T-junction detection.
const EDGE_EPS: f64 = 1.0e-4;
/// Quantisation scale for the welding spatial hash.
const WELD_SCALE: f64 = 1.0 / WELD_TOL;
/// AABB expansion used to decide whether two solids may overlap (and thus need
/// a real CSG union).  Clearly-disjoint solids are never sent to the BSP.
const AABB_EPS: f64 = 1.0e-6;

/// Axis-aligned bounding box used for spatial culling of solid pairs.
#[derive(Clone, Copy)]
struct Aabb {
    min: [f64; 3],
    max: [f64; 3],
}

impl Aabb {
    fn empty() -> Self {
        Self { min: [f64::INFINITY; 3], max: [f64::NEG_INFINITY; 3] }
    }

    fn add(&mut self, p: Point3) {
        let c = [p.x, p.y, p.z];
        for k in 0..3 {
            if c[k] < self.min[k] {
                self.min[k] = c[k];
            }
            if c[k] > self.max[k] {
                self.max[k] = c[k];
            }
        }
    }

    fn from_polys(polys: &[CsgPoly]) -> Self {
        let mut bb = Aabb::empty();
        for p in polys {
            for &v in &p.verts {
                bb.add(v);
            }
        }
        bb
    }

    /// True if the two boxes overlap when expanded by `eps` on each side.
    fn overlaps(&self, o: &Aabb, eps: f64) -> bool {
        for k in 0..3 {
            if self.min[k] - eps > o.max[k] || self.max[k] + eps < o.min[k] {
                return false;
            }
        }
        true
    }
}

/// Options controlling the faceted Boolean union.
#[derive(Clone, Copy, Debug)]
pub struct UnionOptions {
    /// Number of facets around a full circle (cylinders, disks, …).
    pub facets: usize,
}

impl Default for UnionOptions {
    fn default() -> Self {
        Self { facets: 32 }
    }
}

/// Errors returned by the faceted Boolean operations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BooleanError {
    /// An input solid was missing or tessellated to nothing.
    EmptyInput,
    /// A face used a `FaceExtent`/`FaceGeom` combination the faceter cannot handle.
    UnsupportedFace,
    /// Tessellation produced no usable triangles.
    FacetingFailed,
}

/// Watertightness / manifold report for a faceted solid (task §3.4).
#[derive(Clone, Copy, Debug)]
pub struct WatertightReport {
    /// Number of triangles in the solid.
    pub triangles: usize,
    /// Number of welded (unique) vertices.
    pub vertices: usize,
    /// Number of undirected edges.
    pub edges: usize,
    /// Edges whose two directed uses do not cancel (open / naked edges).
    pub open_edges: usize,
    /// Edges shared by a number of faces other than two (non-manifold).
    pub non_manifold_edges: usize,
}

impl WatertightReport {
    /// `true` when the solid is a closed 2-manifold: no naked edges and every
    /// edge shared by exactly two oppositely-oriented faces.
    pub fn is_watertight(&self) -> bool {
        self.open_edges == 0 && self.non_manifold_edges == 0
    }
}

/// Boolean union of two solids in `brep`.
///
/// Both input solids are removed and replaced by a single watertight faceted
/// solid; its id is returned.  The result inherits `a`'s name.
pub fn union_solids(
    brep: &mut BRep,
    a: SolidId,
    b: SolidId,
    opts: &UnionOptions,
) -> Result<SolidId, BooleanError> {
    let pa = solid_to_polys(brep, a, opts.facets)?;
    let pb = solid_to_polys(brep, b, opts.facets)?;
    if pa.is_empty() || pb.is_empty() {
        return Err(BooleanError::EmptyInput);
    }
    let name = brep.solids.get(a).and_then(|s| s.name.clone());
    let merged = polys_union(pa, pb);
    remove_solid(brep, a);
    remove_solid(brep, b);
    finalize(brep, merged, name)
}

/// Boolean union of *every* solid in `brep` into one watertight faceted solid.
///
/// This is the honest replacement for the old container-merge `fuse_solids`:
/// overlapping interior faces are trimmed away instead of being left as
/// non-manifold internal partitions.  Returns `Ok(None)` if there are no solids.
pub fn union_all(brep: &mut BRep, opts: &UnionOptions) -> Result<Option<SolidId>, BooleanError> {
    let ids: Vec<SolidId> = brep.solids.keys().collect();
    if ids.is_empty() {
        return Ok(None);
    }
    let name = brep.solids.get(ids[0]).and_then(|s| s.name.clone());
    let facets = opts.facets.max(8);

    // 1. Tessellate every solid in parallel (one outward-oriented mesh + AABB each).
    let meshes: Vec<(Vec<CsgPoly>, Aabb)> = ids
        .par_iter()
        .map(|&id| {
            let polys = solid_to_polys(brep, id, facets)?;
            let aabb = Aabb::from_polys(&polys);
            Ok::<_, BooleanError>((polys, aabb))
        })
        .collect::<Result<Vec<_>, _>>()?;

    if meshes.iter().all(|(p, _)| p.is_empty()) {
        return Err(BooleanError::EmptyInput);
    }

    // 2. Spatial culling: group solids into connected components by AABB overlap.
    //    Disjoint solids never enter a CSG union — they are merely collected.
    let aabbs: Vec<Aabb> = meshes.iter().map(|(_, a)| *a).collect();
    let components = connected_components(&aabbs, AABB_EPS);

    // 3. Union each component independently (components in parallel; inside a
    //    component a parallel tree-reduction replaces the O(n²) linear fold).
    let merged: Vec<CsgPoly> = components
        .par_iter()
        .map(|comp| {
            let layers: Vec<Vec<CsgPoly>> =
                comp.iter().map(|&i| meshes[i].0.clone()).collect();
            union_reduce(layers)
        })
        .reduce(Vec::new, |mut a, b| {
            a.extend(b);
            a
        });

    if merged.is_empty() {
        return Err(BooleanError::FacetingFailed);
    }

    for &id in &ids {
        remove_solid(brep, id);
    }
    let sid = finalize(brep, merged, name)?;
    Ok(Some(sid))
}

/// Group solids into connected components: two solids are connected when their
/// AABBs overlap (the only pairs that can geometrically intersect).  Union-find
/// over the overlap graph.
fn connected_components(aabbs: &[Aabb], eps: f64) -> Vec<Vec<usize>> {
    let n = aabbs.len();
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]]; // path halving
            x = parent[x];
        }
        x
    }

    for i in 0..n {
        for j in (i + 1)..n {
            if aabbs[i].overlaps(&aabbs[j], eps) {
                let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }

    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    groups.into_values().collect()
}

/// Union a component's solids via a parallel balanced tree-reduction.
///
/// Linear folding (`acc = acc ∪ next`) is O(n) deep and repeatedly re-clips the
/// growing accumulator → ~O(n · P).  A balanced reduction is O(log n) deep,
/// pairs at each level run on separate cores (rayon), and total work drops to
/// ~O(P · log n).
fn union_reduce(mut layers: Vec<Vec<CsgPoly>>) -> Vec<CsgPoly> {
    layers.retain(|l| !l.is_empty());
    if layers.is_empty() {
        return Vec::new();
    }
    while layers.len() > 1 {
        let half = layers.len() / 2;
        let mut next: Vec<Vec<CsgPoly>> = (0..half)
            .into_par_iter()
            .map(|i| polys_union(layers[2 * i].clone(), layers[2 * i + 1].clone()))
            .collect();
        if layers.len() % 2 == 1 {
            next.push(layers.last().unwrap().clone());
        }
        layers = next;
    }
    layers.pop().unwrap_or_default()
}

/// Validate that a faceted solid (made of `FaceExtent::Polygon` faces) is a
/// closed, oriented 2-manifold.  See [`WatertightReport`].
pub fn validate_watertight(brep: &BRep, solid: SolidId) -> WatertightReport {
    let mut tris: Vec<[Point3; 3]> = Vec::new();
    if let Some(s) = brep.solids.get(solid) {
        for &sh in &s.shells {
            let shell = match brep.shells.get(sh) {
                Some(x) => x,
                None => continue,
            };
            for &fid in &shell.faces {
                if let Some(f) = brep.faces.get(fid) {
                    if let FaceExtent::Polygon { points } = &f.extent {
                        if points.len() >= 3 {
                            for i in 1..points.len() - 1 {
                                tris.push([points[0], points[i], points[i + 1]]);
                            }
                        }
                    }
                }
            }
        }
    }
    let (verts, idx) = weld(&tris);
    report_from(&idx, verts.len())
}

// ── BSP CSG (port of the csg.js union algorithm) ────────────────────────────────

/// A convex polygon used by the BSP CSG, carrying its support plane.
#[derive(Clone)]
struct CsgPoly {
    verts: Vec<Point3>,
    normal: Vec3,
    w: f64,
}

impl CsgPoly {
    /// Build from a vertex ring, deriving the plane via Newell's method.
    fn new(verts: Vec<Point3>) -> Option<Self> {
        let m = verts.len();
        if m < 3 {
            return None;
        }
        let mut n = Vec3::ZERO;
        for i in 0..m {
            let a = verts[i];
            let b = verts[(i + 1) % m];
            n.x += (a.y - b.y) * (a.z + b.z);
            n.y += (a.z - b.z) * (a.x + b.x);
            n.z += (a.x - b.x) * (a.y + b.y);
        }
        let nn = n.try_normalize()?;
        let w = nn.dot(verts[0].to_vec());
        Some(Self { verts, normal: nn, w })
    }

    fn flip(&mut self) {
        self.verts.reverse();
        self.normal = -self.normal;
        self.w = -self.w;
    }
}

const C_COPLANAR: u8 = 0;
const C_FRONT: u8 = 1;
const C_BACK: u8 = 2;
const C_SPANNING: u8 = 3;

/// Split `poly` against the plane `(pn, pw)`; append the pieces to the four
/// output buckets (coplanar-front, coplanar-back, front, back).
fn split_poly(
    pn: Vec3,
    pw: f64,
    poly: &CsgPoly,
    coplanar_front: &mut Vec<CsgPoly>,
    coplanar_back: &mut Vec<CsgPoly>,
    front: &mut Vec<CsgPoly>,
    back: &mut Vec<CsgPoly>,
) {
    let mut poly_type: u8 = 0;
    let mut types: Vec<u8> = Vec::with_capacity(poly.verts.len());
    for v in &poly.verts {
        let t = pn.dot(v.to_vec()) - pw;
        let ty = if t < -BSP_EPS {
            C_BACK
        } else if t > BSP_EPS {
            C_FRONT
        } else {
            C_COPLANAR
        };
        poly_type |= ty;
        types.push(ty);
    }

    match poly_type {
        C_COPLANAR => {
            if pn.dot(poly.normal) > 0.0 {
                coplanar_front.push(poly.clone());
            } else {
                coplanar_back.push(poly.clone());
            }
        }
        C_FRONT => front.push(poly.clone()),
        C_BACK => back.push(poly.clone()),
        _ => {
            // Spanning — split along the plane.
            let mut fv: Vec<Point3> = Vec::new();
            let mut bv: Vec<Point3> = Vec::new();
            let n = poly.verts.len();
            for i in 0..n {
                let j = (i + 1) % n;
                let ti = types[i];
                let tj = types[j];
                let vi = poly.verts[i];
                let vj = poly.verts[j];
                if ti != C_BACK {
                    fv.push(vi);
                }
                if ti != C_FRONT {
                    bv.push(vi);
                }
                if (ti | tj) == C_SPANNING {
                    let denom = pn.dot(vj - vi);
                    if denom.abs() > 1.0e-18 {
                        let tt = (pw - pn.dot(vi.to_vec())) / denom;
                        let v = vi + (vj - vi) * tt;
                        fv.push(v);
                        bv.push(v);
                    }
                }
            }
            if let Some(p) = CsgPoly::new(fv) {
                front.push(p);
            }
            if let Some(p) = CsgPoly::new(bv) {
                back.push(p);
            }
        }
    }
}

/// A node in the BSP tree.
struct Node {
    plane: Option<(Vec3, f64)>,
    front: Option<Box<Node>>,
    back: Option<Box<Node>>,
    polys: Vec<CsgPoly>,
}

impl Node {
    fn new() -> Self {
        Self { plane: None, front: None, back: None, polys: Vec::new() }
    }

    fn from_polys(polys: Vec<CsgPoly>) -> Self {
        let mut n = Node::new();
        n.build(polys);
        n
    }

    fn build(&mut self, polys: Vec<CsgPoly>) {
        if polys.is_empty() {
            return;
        }
        if self.plane.is_none() {
            self.plane = Some((polys[0].normal, polys[0].w));
        }
        let (pn, pw) = self.plane.unwrap();
        let mut cf = Vec::new();
        let mut cb = Vec::new();
        let mut f = Vec::new();
        let mut b = Vec::new();
        for p in &polys {
            split_poly(pn, pw, p, &mut cf, &mut cb, &mut f, &mut b);
        }
        self.polys.extend(cf);
        self.polys.extend(cb);
        if !f.is_empty() {
            self.front.get_or_insert_with(|| Box::new(Node::new())).build(f);
        }
        if !b.is_empty() {
            self.back.get_or_insert_with(|| Box::new(Node::new())).build(b);
        }
    }

    fn clip_polygons(&self, polys: Vec<CsgPoly>) -> Vec<CsgPoly> {
        let (pn, pw) = match self.plane {
            Some(p) => p,
            None => return polys,
        };
        let mut cf = Vec::new();
        let mut cb = Vec::new();
        let mut f = Vec::new();
        let mut b = Vec::new();
        for p in &polys {
            split_poly(pn, pw, p, &mut cf, &mut cb, &mut f, &mut b);
        }
        f.extend(cf);
        b.extend(cb);
        let mut f = match &self.front {
            Some(node) => node.clip_polygons(f),
            None => f,
        };
        let b = match &self.back {
            Some(node) => node.clip_polygons(b),
            None => Vec::new(),
        };
        f.extend(b);
        f
    }

    fn clip_to(&mut self, bsp: &Node) {
        self.polys = bsp.clip_polygons(std::mem::take(&mut self.polys));
        if let Some(n) = self.front.as_mut() {
            n.clip_to(bsp);
        }
        if let Some(n) = self.back.as_mut() {
            n.clip_to(bsp);
        }
    }

    fn invert(&mut self) {
        for p in self.polys.iter_mut() {
            p.flip();
        }
        if let Some((n, w)) = self.plane {
            self.plane = Some((-n, -w));
        }
        if let Some(n) = self.front.as_mut() {
            n.invert();
        }
        if let Some(n) = self.back.as_mut() {
            n.invert();
        }
        std::mem::swap(&mut self.front, &mut self.back);
    }

    fn all_polygons(&self) -> Vec<CsgPoly> {
        let mut out = self.polys.clone();
        if let Some(n) = &self.front {
            out.extend(n.all_polygons());
        }
        if let Some(n) = &self.back {
            out.extend(n.all_polygons());
        }
        out
    }
}

/// CSG union of two polygon sets (`a ∪ b`).
fn polys_union(a: Vec<CsgPoly>, b: Vec<CsgPoly>) -> Vec<CsgPoly> {
    let mut a = Node::from_polys(a);
    let mut b = Node::from_polys(b);
    a.clip_to(&b);
    b.clip_to(&a);
    b.invert();
    b.clip_to(&a);
    b.invert();
    a.build(b.all_polygons());
    a.all_polygons()
}

// ── Tessellation (FaceExtent → outward-oriented triangles) ──────────────────────

/// Flip an axis to a canonical sign so that a ring is identical regardless of the
/// axis orientation it was generated from (lets cap and cylinder-end rings weld).
fn canon_axis(a: UnitVec3) -> UnitVec3 {
    let v = a.as_vec();
    let s = if v.x.abs() > 1.0e-12 {
        v.x
    } else if v.y.abs() > 1.0e-12 {
        v.y
    } else {
        v.z
    };
    if s < 0.0 {
        -a
    } else {
        a
    }
}

/// A canonical ring of `n` points on a circle (centre, axis, radius).
fn ring(center: Point3, axis: UnitVec3, r: f64, n: usize) -> Vec<Point3> {
    let a = canon_axis(axis);
    let (u, v) = a.perp_basis();
    (0..n)
        .map(|i| {
            let th = TAU * (i as f64) / (n as f64);
            center + u * (r * th.cos()) + v * (r * th.sin())
        })
        .collect()
}

/// Outward radial direction at `p` on a cylinder of axis `(c0, axis)`.
fn radial(c0: Point3, axis: UnitVec3, p: Point3, flip: bool) -> Vec3 {
    let foot = c0 + axis * axis.dot_vec(p - c0);
    let d = p - foot;
    let u = d.try_normalize().unwrap_or_else(|| axis.as_vec());
    if flip {
        -u
    } else {
        u
    }
}

/// Push a triangle oriented so its geometric normal agrees with `outward`.
fn push_tri(out: &mut Vec<CsgPoly>, a: Point3, b: Point3, c: Point3, outward: Vec3) {
    let g = (b - a).cross(c - a);
    if g.length_sq() < 1.0e-24 {
        return;
    }
    let (b, c) = if g.dot(outward) < 0.0 { (c, b) } else { (b, c) };
    if let Some(p) = CsgPoly::new(vec![a, b, c]) {
        out.push(p);
    }
}

/// Tessellate every face of a solid into outward-oriented triangles.
fn solid_to_polys(brep: &BRep, solid: SolidId, n: usize) -> Result<Vec<CsgPoly>, BooleanError> {
    let n = n.max(8);
    let solid = brep.solids.get(solid).ok_or(BooleanError::EmptyInput)?;
    let mut out = Vec::new();
    for &sh in &solid.shells {
        let shell = match brep.shells.get(sh) {
            Some(x) => x,
            None => continue,
        };
        for &fid in &shell.faces {
            let face = match brep.faces.get(fid) {
                Some(x) => x,
                None => continue,
            };
            facet_face(face, n, &mut out)?;
        }
    }
    Ok(out)
}

fn facet_face(face: &Face, n: usize, out: &mut Vec<CsgPoly>) -> Result<(), BooleanError> {
    let flip = matches!(face.normal, FaceNormal::Reversed);
    match (&face.geom, &face.extent) {
        (FaceGeom::Cylinder(cyl), FaceExtent::Cylinder { length, .. }) => {
            facet_cylinder(*cyl, *length, n, flip, out);
        }
        (FaceGeom::Cylinder(cyl), FaceExtent::PartialCylinder {
            length,
            arc_start_angle,
            arc_end_angle,
            arc_ref_dir,
        }) => {
            facet_partial_cylinder(
                *cyl, *length, *arc_start_angle, *arc_end_angle, *arc_ref_dir, n, flip, out,
            );
        }
        (FaceGeom::Plane(p), FaceExtent::Disk { radius }) => {
            facet_disk(*p, *radius, n, flip, out);
        }
        (FaceGeom::Plane(p), FaceExtent::PartialDisk { radius, start_angle, end_angle }) => {
            facet_partial_disk(*p, *radius, *start_angle, *end_angle, n, flip, out);
        }
        (FaceGeom::Plane(p), FaceExtent::Polygon { points }) => {
            facet_polygon(*p, points, flip, out);
        }
        (FaceGeom::Plane(p), FaceExtent::PlanarBoundary { boundary }) => {
            facet_planar_boundary(*p, boundary, n, flip, out);
        }
        (FaceGeom::Torus(t), FaceExtent::TorusFillet { start_circle, end_circle }) => {
            facet_torus_fillet(*t, *start_circle, *end_circle, n, flip, out);
        }
        (FaceGeom::Sphere(s), _) => {
            facet_sphere(*s, n, flip, out);
        }
        _ => return Err(BooleanError::UnsupportedFace),
    }
    Ok(())
}

fn facet_cylinder(cyl: CylSurf, length: f64, n: usize, flip: bool, out: &mut Vec<CsgPoly>) {
    let axis = cyl.frame.z;
    let c0 = cyl.frame.origin;
    let c1 = c0 + axis * length;
    let r = cyl.radius;
    let r0 = ring(c0, axis, r, n);
    let r1 = ring(c1, axis, r, n);
    for i in 0..n {
        let j = (i + 1) % n;
        let (a, b, c, d) = (r0[i], r0[j], r1[j], r1[i]);
        push_tri(out, a, b, c, radial(c0, axis, b, flip));
        push_tri(out, a, c, d, radial(c0, axis, a, flip));
    }
}

fn facet_partial_cylinder(
    cyl: CylSurf,
    length: f64,
    start: f64,
    end: f64,
    arc_ref_dir: UnitVec3,
    n: usize,
    flip: bool,
    out: &mut Vec<CsgPoly>,
) {
    let axis = cyl.frame.z;
    let c0 = cyl.frame.origin;
    let r = cyl.radius;
    let right = UnitVec3::try_from_vec(axis.cross(arc_ref_dir)).unwrap_or(cyl.frame.x);
    let range = end - start;
    let mut a_ring = Vec::with_capacity(n + 1);
    let mut b_ring = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let th = start + range * (i as f64 / n as f64);
        let local = arc_ref_dir.as_vec() * (r * th.cos()) + right.as_vec() * (r * th.sin());
        a_ring.push(c0 + local);
        b_ring.push(c0 + axis * length + local);
    }
    for i in 0..n {
        let (a, b, c, d) = (a_ring[i], a_ring[i + 1], b_ring[i + 1], b_ring[i]);
        push_tri(out, a, b, c, radial(c0, axis, b, flip));
        push_tri(out, a, c, d, radial(c0, axis, a, flip));
    }
}

fn facet_disk(plane: Plane3, r: f64, n: usize, flip: bool, out: &mut Vec<CsgPoly>) {
    let c = plane.frame.origin;
    let nrm = if flip { -plane.frame.z } else { plane.frame.z };
    let pts = ring(c, plane.frame.z, r, n);
    for i in 0..n {
        let j = (i + 1) % n;
        push_tri(out, c, pts[i], pts[j], nrm.as_vec());
    }
}

fn facet_partial_disk(
    plane: Plane3,
    r: f64,
    start: f64,
    end: f64,
    n: usize,
    flip: bool,
    out: &mut Vec<CsgPoly>,
) {
    let c = plane.frame.origin;
    let nrm = if flip { -plane.frame.z } else { plane.frame.z };
    let x = plane.frame.x;
    let y = plane.frame.y;
    let range = end - start;
    let mut prev: Option<Point3> = None;
    for i in 0..=n {
        let th = start + range * (i as f64 / n as f64);
        let p = c + x * (r * th.cos()) + y * (r * th.sin());
        if let Some(pp) = prev {
            push_tri(out, c, pp, p, nrm.as_vec());
        }
        prev = Some(p);
    }
}

fn facet_polygon(plane: Plane3, points: &[Point3], flip: bool, out: &mut Vec<CsgPoly>) {
    if points.len() < 3 {
        return;
    }
    let nrm = if flip { -plane.frame.z } else { plane.frame.z };
    for i in 1..points.len() - 1 {
        push_tri(out, points[0], points[i], points[i + 1], nrm.as_vec());
    }
}

fn facet_planar_boundary(
    plane: Plane3,
    boundary: &FaceBoundary,
    n: usize,
    flip: bool,
    out: &mut Vec<CsgPoly>,
) {
    let c = plane.frame.origin;
    let nrm = if flip { -plane.frame.z } else { plane.frame.z };
    match boundary {
        FaceBoundary::Circle(circ) => {
            let pts = ring(c, plane.frame.z, circ.radius, n);
            for i in 0..n {
                let j = (i + 1) % n;
                push_tri(out, c, pts[i], pts[j], nrm.as_vec());
            }
        }
        FaceBoundary::Ellipse(e) => {
            let mut prev = e.point_at(0.0);
            for i in 1..=n {
                let th = TAU * (i as f64) / (n as f64);
                let p = e.point_at(th);
                push_tri(out, e.frame.origin, prev, p, nrm.as_vec());
                prev = p;
            }
        }
    }
}

fn torus_theta(t: &TorusSurf, p: Point3) -> f64 {
    let d = p - t.frame.origin;
    let x = t.frame.x.dot_vec(d);
    let y = t.frame.y.dot_vec(d);
    y.atan2(x)
}

fn facet_torus_fillet(
    t: TorusSurf,
    start_circle: Circle3,
    end_circle: Circle3,
    n: usize,
    flip: bool,
    out: &mut Vec<CsgPoly>,
) {
    let ts = torus_theta(&t, start_circle.frame.origin);
    let te = torus_theta(&t, end_circle.frame.origin);
    let mut dte = te - ts;
    while dte > PI {
        dte -= TAU;
    }
    while dte <= -PI {
        dte += TAU;
    }
    let tsteps = n.max(8);
    let psteps = n.max(8);
    for i in 0..tsteps {
        let th0 = ts + dte * (i as f64 / tsteps as f64);
        let th1 = ts + dte * ((i + 1) as f64 / tsteps as f64);
        for k in 0..psteps {
            let ph0 = TAU * (k as f64) / (psteps as f64);
            let ph1 = TAU * ((k + 1) as f64) / (psteps as f64);
            let p00 = t.point_at(th0, ph0);
            let p10 = t.point_at(th1, ph0);
            let p11 = t.point_at(th1, ph1);
            let p01 = t.point_at(th0, ph1);
            let nm = t.normal_at((th0 + th1) * 0.5, (ph0 + ph1) * 0.5);
            let nm = if flip { -nm } else { nm };
            push_tri(out, p00, p10, p11, nm.as_vec());
            push_tri(out, p00, p11, p01, nm.as_vec());
        }
    }
}

fn facet_sphere(s: SphereSurf, n: usize, flip: bool, out: &mut Vec<CsgPoly>) {
    let lat = n.max(4);
    let lon = n.max(6);
    for i in 0..lat {
        let ph0 = -PI / 2.0 + PI * (i as f64 / lat as f64);
        let ph1 = -PI / 2.0 + PI * ((i + 1) as f64 / lat as f64);
        for k in 0..lon {
            let th0 = TAU * (k as f64) / (lon as f64);
            let th1 = TAU * ((k + 1) as f64) / (lon as f64);
            let p00 = s.point_at(th0, ph0);
            let p10 = s.point_at(th1, ph0);
            let p11 = s.point_at(th1, ph1);
            let p01 = s.point_at(th0, ph1);
            let outward = |p: Point3| -> Vec3 {
                let d = (p - s.centre).try_normalize().unwrap_or(Vec3::Z);
                if flip {
                    -d
                } else {
                    d
                }
            };
            push_tri(out, p00, p10, p11, outward(p10));
            push_tri(out, p00, p11, p01, outward(p01));
        }
    }
}

// ── Welding, T-junction healing, validation, materialisation ────────────────────

/// Index a point into a shared vertex pool, welding within `WELD_TOL`.
fn weld_index(
    verts: &mut Vec<Point3>,
    map: &mut HashMap<(i64, i64, i64), Vec<usize>>,
    p: Point3,
) -> usize {
    let key = (
        (p.x * WELD_SCALE).round() as i64,
        (p.y * WELD_SCALE).round() as i64,
        (p.z * WELD_SCALE).round() as i64,
    );
    for dx in -1..=1 {
        for dy in -1..=1 {
            for dz in -1..=1 {
                let k = (key.0 + dx, key.1 + dy, key.2 + dz);
                if let Some(list) = map.get(&k) {
                    for &i in list {
                        if (verts[i] - p).length() < WELD_TOL {
                            return i;
                        }
                    }
                }
            }
        }
    }
    let i = verts.len();
    verts.push(p);
    map.entry(key).or_default().push(i);
    i
}

/// Weld a triangle soup into a shared vertex pool + index triples.
fn weld(tris: &[[Point3; 3]]) -> (Vec<Point3>, Vec<[usize; 3]>) {
    let mut verts: Vec<Point3> = Vec::new();
    let mut map: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    let mut idx: Vec<[usize; 3]> = Vec::with_capacity(tris.len());
    for t in tris {
        let a = weld_index(&mut verts, &mut map, t[0]);
        let b = weld_index(&mut verts, &mut map, t[1]);
        let c = weld_index(&mut verts, &mut map, t[2]);
        if a == b || b == c || c == a {
            continue;
        }
        idx.push([a, b, c]);
    }
    (verts, idx)
}

/// Heal T-junctions: where a welded vertex lies on the interior of a triangle
/// edge, re-triangulate that triangle (fan from its centroid through the split
/// boundary) so every edge is shared by exactly two triangles.
fn heal(verts: &mut Vec<Point3>, tris: &[[usize; 3]]) -> Vec<[usize; 3]> {
    let nv = verts.len();
    if nv == 0 {
        return tris.to_vec();
    }

    // Spatial grid over the original vertices so each edge only tests nearby
    // candidates.  Without this, heal is O(triangles × vertices) and becomes the
    // dominant cost on large meshes.  Cell size ≈ mean edge length → ~O(1) per cell.
    let h = {
        let mut sum = 0.0f64;
        let mut cnt = 0usize;
        for t in tris {
            for e in 0..3 {
                sum += (verts[t[(e + 1) % 3]] - verts[t[e]]).length();
                cnt += 1;
            }
        }
        if cnt == 0 { 1.0 } else { (sum / cnt as f64).max(1.0e-6) }
    };
    let inv = 1.0 / h;
    let cell = |p: Point3| -> (i64, i64, i64) {
        ((p.x * inv).floor() as i64, (p.y * inv).floor() as i64, (p.z * inv).floor() as i64)
    };
    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();
    for i in 0..nv {
        grid.entry(cell(verts[i])).or_default().push(i);
    }

    let mut out: Vec<[usize; 3]> = Vec::with_capacity(tris.len());
    let mut near: Vec<usize> = Vec::new();
    for t in tris {
        let corners = [t[0], t[1], t[2]];
        let mut loop_pts: Vec<usize> = Vec::new();
        let mut any = false;
        for e in 0..3 {
            let i = corners[e];
            let j = corners[(e + 1) % 3];
            loop_pts.push(i);
            let a = verts[i];
            let b = verts[j];
            let ab = b - a;
            let len2 = ab.length_sq();
            if len2 > 1.0e-24 {
                // Gather candidate vertices from the grid cells spanning this edge.
                near.clear();
                let (c0, c1) = (cell(a), cell(b));
                for cx in (c0.0.min(c1.0) - 1)..=(c0.0.max(c1.0) + 1) {
                    for cy in (c0.1.min(c1.1) - 1)..=(c0.1.max(c1.1) + 1) {
                        for cz in (c0.2.min(c1.2) - 1)..=(c0.2.max(c1.2) + 1) {
                            if let Some(list) = grid.get(&(cx, cy, cz)) {
                                near.extend_from_slice(list);
                            }
                        }
                    }
                }
                let mut on: Vec<(f64, usize)> = Vec::new();
                for &k in &near {
                    if k == i || k == j {
                        continue;
                    }
                    let p = verts[k];
                    let tp = (p - a).dot(ab) / len2;
                    if tp > EDGE_EPS && tp < 1.0 - EDGE_EPS {
                        let proj = a + ab * tp;
                        if (p - proj).length() < WELD_TOL {
                            on.push((tp, k));
                        }
                    }
                }
                if !on.is_empty() {
                    any = true;
                    on.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
                    on.dedup_by_key(|x| x.1);
                    for (_, k) in on {
                        loop_pts.push(k);
                    }
                }
            }
        }
        if !any {
            out.push(*t);
            continue;
        }
        // Fan from the centroid of the (now subdivided) boundary loop.
        let m = loop_pts.len() as f64;
        let (mut cx, mut cy, mut cz) = (0.0, 0.0, 0.0);
        for &p in &loop_pts {
            let v = verts[p];
            cx += v.x;
            cy += v.y;
            cz += v.z;
        }
        let cen = Point3::new(cx / m, cy / m, cz / m);
        let ci = verts.len();
        verts.push(cen);
        let nl = loop_pts.len();
        for e in 0..nl {
            let p0 = loop_pts[e];
            let p1 = loop_pts[(e + 1) % nl];
            if p0 != p1 {
                out.push([ci, p0, p1]);
            }
        }
    }
    out
}

/// Build a [`WatertightReport`] from welded index triples.
fn report_from(idx: &[[usize; 3]], nverts: usize) -> WatertightReport {
    let mut dir: HashMap<(usize, usize), i32> = HashMap::new();
    for t in idx {
        for k in 0..3 {
            let a = t[k];
            let b = t[(k + 1) % 3];
            if a == b {
                continue;
            }
            *dir.entry((a, b)).or_insert(0) += 1;
        }
    }
    let mut und: HashMap<(usize, usize), i32> = HashMap::new();
    for (&(a, b), &c) in &dir {
        *und.entry((a.min(b), a.max(b))).or_insert(0) += c;
    }
    let mut open = 0;
    for (&(a, b), &c) in &dir {
        let rev = *dir.get(&(b, a)).unwrap_or(&0);
        if c != rev {
            open += 1;
        }
    }
    let mut nm = 0;
    for (_, &c) in &und {
        if c != 2 {
            nm += 1;
        }
    }
    WatertightReport {
        triangles: idx.len(),
        vertices: nverts,
        edges: und.len(),
        open_edges: open,
        non_manifold_edges: nm,
    }
}

/// Triangulate, weld, heal, and materialise CSG polygons into a new solid built
/// from `FaceExtent::Polygon` faces.
fn finalize(
    brep: &mut BRep,
    polys: Vec<CsgPoly>,
    name: Option<String>,
) -> Result<SolidId, BooleanError> {
    let mut tris: Vec<[Point3; 3]> = Vec::new();
    for p in &polys {
        if p.verts.len() < 3 {
            continue;
        }
        for i in 1..p.verts.len() - 1 {
            tris.push([p.verts[0], p.verts[i], p.verts[i + 1]]);
        }
    }
    if tris.is_empty() {
        return Err(BooleanError::FacetingFailed);
    }
    let (mut verts, idx) = weld(&tris);
    let idx = heal(&mut verts, &idx);

    let mut face_ids: Vec<FaceId> = Vec::with_capacity(idx.len());
    for t in &idx {
        let a = verts[t[0]];
        let b = verts[t[1]];
        let c = verts[t[2]];
        let g = (b - a).cross(c - a);
        let nrm = match g.try_normalize() {
            Some(u) => u,
            None => continue,
        };
        let normal = UnitVec3::new_unchecked(nrm);
        let plane = Plane3::from_origin_normal(a, normal);
        let loop_id = brep.add_loop(cadcore_topo::Loop {
            start: cadcore_topo::CoEdgeId::default(),
            face: FaceId::default(),
        });
        let fid = brep.add_face(Face {
            geom: FaceGeom::Plane(plane),
            normal: FaceNormal::Same,
            outer_loop: loop_id,
            inner_loops: vec![],
            shell: ShellId::default(),
            extent: FaceExtent::Polygon { points: vec![a, b, c] },
        });
        face_ids.push(fid);
    }
    if face_ids.is_empty() {
        return Err(BooleanError::FacetingFailed);
    }
    let shell_id = brep.add_shell(Shell {
        faces: face_ids.clone(),
        is_outer: true,
        solid: SolidId::default(),
    });
    for &fid in &face_ids {
        if let Some(f) = brep.faces.get_mut(fid) {
            f.shell = shell_id;
        }
    }
    let solid_id = brep.add_solid(Solid { shells: vec![shell_id], name });
    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.solid = solid_id;
    }
    Ok(solid_id)
}

/// Remove a solid and its shells/faces from the B-Rep.
fn remove_solid(brep: &mut BRep, id: SolidId) {
    if let Some(solid) = brep.solids.get(id).cloned() {
        for sh in &solid.shells {
            if let Some(shell) = brep.shells.get(*sh).cloned() {
                for f in &shell.faces {
                    brep.faces.remove(*f);
                }
                brep.shells.remove(*sh);
            }
        }
        brep.solids.remove(id);
    }
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

    // ── Boolean UNION (faceted) ─────────────────────────────────────────────

    /// A single swept cylinder (lateral surface + two disk caps) tessellates to
    /// a closed manifold mesh.  This is the precondition for a correct CSG union.
    #[test]
    fn single_cylinder_tessellates_watertight() {
        let mut brep = BRep::new();
        let cyl = sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(-2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            0.5,
            &SweepOptions::default(),
        )
        .unwrap();
        let polys = solid_to_polys(&brep, cyl, 32).unwrap();
        let mut tris = Vec::new();
        for p in &polys {
            for i in 1..p.verts.len() - 1 {
                tris.push([p.verts[0], p.verts[i], p.verts[i + 1]]);
            }
        }
        let (verts, idx) = weld(&tris);
        let rep = report_from(&idx, verts.len());
        assert!(
            rep.is_watertight(),
            "tessellated cylinder must be closed: {rep:?}"
        );
    }

    /// Two equal-radius cylinders whose axes cross perpendicularly through the
    /// origin must fuse into a SINGLE watertight (closed, oriented 2-manifold)
    /// solid — the classic Steinmetz-cross union.
    #[test]
    fn perpendicular_cylinders_fuse_watertight() {
        let mut brep = BRep::new();
        let r = 0.5;
        // Cylinder A along X, cylinder B along Z, both through the origin.
        let a = sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(-2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            r,
            &SweepOptions::default(),
        )
        .unwrap();
        let b = sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(0.0, 0.0, -2.0), Point3::new(0.0, 0.0, 2.0)],
            r,
            &SweepOptions::default(),
        )
        .unwrap();
        assert_eq!(brep.solids.len(), 2);

        let fused = union_solids(&mut brep, a, b, &UnionOptions { facets: 32 }).unwrap();

        // Exactly one solid remains.
        assert_eq!(brep.solids.len(), 1, "union must leave a single solid");

        // The fused solid is a closed, oriented 2-manifold: no naked edges,
        // every edge shared by exactly two oppositely-oriented faces (task §3.4).
        let rep = validate_watertight(&brep, fused);
        assert!(
            rep.is_watertight(),
            "fused perpendicular cylinders must be watertight: {rep:?}"
        );
        assert!(rep.triangles > 0 && rep.vertices > 0);

        // Euler characteristic of a genus-0 closed surface: V - E + F = 2.
        let euler = rep.vertices as i64 - rep.edges as i64 + rep.triangles as i64;
        assert_eq!(euler, 2, "fused cross must be a topological sphere: {rep:?}");
    }

    /// Two disjoint cylinders fuse to two shells worth of geometry, still a valid
    /// closed surface set (each component watertight).
    #[test]
    fn perpendicular_cylinders_offset_still_watertight() {
        let mut brep = BRep::new();
        let r = 0.4;
        let a = sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(-2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
            r,
            &SweepOptions::default(),
        )
        .unwrap();
        // Unequal radius crossing — exercises the general (non-elliptic) seam.
        let b = sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(0.0, 0.0, -2.0), Point3::new(0.0, 0.0, 2.0)],
            r * 1.7,
            &SweepOptions::default(),
        )
        .unwrap();
        let fused = union_solids(&mut brep, a, b, &UnionOptions { facets: 28 }).unwrap();
        let rep = validate_watertight(&brep, fused);
        assert!(
            rep.is_watertight(),
            "unequal-radius perpendicular cylinders must fuse watertight: {rep:?}"
        );
    }
}
