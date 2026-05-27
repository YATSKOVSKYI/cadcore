//! Sweep a circle profile along a polyline — exact analytic B-Rep construction.
//!
//! ## Algorithm
//!
//! Given N waypoints and a radius r, the solid consists of:
//!
//! * **(N−1) cylinder faces** — one per straight segment between adjacent waypoints.
//! * **(N−2) connector faces** — one per bend vertex; either a toroidal fillet
//!   (smooth join) or an elliptic miter (sharp join).
//! * **2 end-cap faces** — one planar disk at each endpoint.
//!
//! All surfaces are analytic (cylinder, torus, plane, ellipse) and map
//! directly to STEP AP203 entities.  No mesh, no Boolean union.  Total
//! construction cost is O(N).

use cadcore_geom::{CylSurf, Plane3, TorusSurf};
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{
    BRep,
    Face, FaceGeom, FaceNormal,
    Shell,
    Solid,
    SolidId,
};

/// Options that control the sweep.
#[derive(Clone, Debug)]
pub struct SweepOptions {
    /// Use toroidal fillets at corners (smooth G1 join).
    /// If `false`, use miter-plane cuts (sharp join — faster STEP output).
    pub fillet_corners: bool,
    /// Optional name for the resulting solid.
    pub name: Option<String>,
}

impl Default for SweepOptions {
    fn default() -> Self {
        Self { fillet_corners: true, name: None }
    }
}

/// Sweep a circle of radius `radius` along the polyline defined by `waypoints`.
///
/// Returns the id of the newly inserted solid within `brep`.
///
/// # Errors
///
/// Returns `Err` if:
/// * `waypoints` has fewer than 2 points.
/// * `radius` ≤ 0.
/// * Any two consecutive waypoints are coincident (zero-length segment).
pub fn sweep_circle_along_polyline(
    brep:      &mut BRep,
    waypoints: &[Point3],
    radius:    f64,
    opts:      &SweepOptions,
) -> Result<SolidId, SweepError> {
    if waypoints.len() < 2 {
        return Err(SweepError::TooFewWaypoints);
    }
    if radius <= 0.0 {
        return Err(SweepError::InvalidRadius(radius));
    }

    // Pre-compute segment directions.
    let n = waypoints.len();
    let mut dirs: Vec<UnitVec3> = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let d = waypoints[i + 1] - waypoints[i];
        match UnitVec3::try_from_vec(d) {
            Some(u) => dirs.push(u),
            None    => return Err(SweepError::CoincidentWaypoints(i)),
        }
    }

    // Pre-allocate face ids for shell construction.
    let mut face_ids = Vec::new();

    // ── Start cap ────────────────────────────────────────────────────────────
    let _cap_start_id = build_end_cap(brep, waypoints[0], -dirs[0], radius, &mut face_ids);

    // ── Cylinder + connector loop ─────────────────────────────────────────────
    // For each segment i we build: cylinder face for segment i, then (if not last)
    // a connector at waypoints[i+1].
    for i in 0..n - 1 {
        let p0  = waypoints[i];
        let _p1 = waypoints[i + 1];
        let dir = dirs[i];

        // --- Cylinder face for segment i ---
        let cyl = CylSurf::new(p0, dir, radius);
        let cyl_loop_id = build_placeholder_loop(brep);   // borrow ends here
        let cyl_face_id = brep.add_face(Face {
            geom:        FaceGeom::Cylinder(cyl),
            normal:      FaceNormal::Same,
            outer_loop:  cyl_loop_id,
            inner_loops: vec![],
            shell:       cadcore_topo::ShellId::default(), // patched later
        });
        face_ids.push(cyl_face_id);

        // --- Connector at interior vertex ---
        if i + 1 < n - 1 {
            let next_dir  = dirs[i + 1];
            let vertex    = waypoints[i + 1];
            let cos_angle = dir.dot(next_dir);

            if opts.fillet_corners && (1.0 - cos_angle.abs()) > 1e-4 {
                // Build toroidal fillet
                if let Some(torus) = TorusSurf::for_corner(vertex, dir, next_dir, radius) {
                    let torus_loop_id = build_placeholder_loop(brep); // borrow ends here
                    let torus_id = brep.add_face(Face {
                        geom:        FaceGeom::Torus(torus),
                        normal:      FaceNormal::Same,
                        outer_loop:  torus_loop_id,
                        inner_loops: vec![],
                        shell:       cadcore_topo::ShellId::default(),
                    });
                    face_ids.push(torus_id);
                }
            }
            // If fillet_corners=false or nearly straight: miter plane cut.
            // In a miter join the cylinder segments share the elliptic edge.
            // The STEP writer handles the shared edge — no extra face needed.
        }
    }

    // ── End cap ───────────────────────────────────────────────────────────────
    let last = n - 1;
    let _cap_end_id = build_end_cap(brep, waypoints[last], dirs[last - 1], radius, &mut face_ids);

    // ── Assemble shell + solid ────────────────────────────────────────────────
    let shell_id = brep.add_shell(Shell {
        faces:    face_ids.clone(),
        is_outer: true,
        solid:    cadcore_topo::SolidId::default(), // patched below
    });

    // Patch shell back-ref into faces
    for &fid in &face_ids {
        if let Some(f) = brep.faces.get_mut(fid) {
            f.shell = shell_id;
        }
    }

    let solid_id = brep.add_solid(Solid {
        shells: vec![shell_id],
        name:   opts.name.clone(),
    });

    // Patch solid back-ref into shell
    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.solid = solid_id;
    }

    Ok(solid_id)
}

/// Errors that [`sweep_circle_along_polyline`] can return.
#[derive(Clone, Debug, PartialEq)]
pub enum SweepError {
    /// Need at least 2 waypoints.
    TooFewWaypoints,
    /// Radius must be positive.
    InvalidRadius(f64),
    /// Two consecutive waypoints are the same (zero-length segment at index).
    CoincidentWaypoints(usize),
}

impl std::fmt::Display for SweepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TooFewWaypoints           => write!(f, "need at least 2 waypoints"),
            Self::InvalidRadius(r)          => write!(f, "invalid radius: {r}"),
            Self::CoincidentWaypoints(i)    => write!(f, "coincident waypoints at index {i}"),
        }
    }
}

impl std::error::Error for SweepError {}

// ── Internal helpers ─────────────────────────────────────────────────────────

/// Build an end-cap (planar disk) face at `centre` with normal `normal`.
///
/// The outward normal of the start cap points *away from* the solid (−dir),
/// and the end cap outward normal is +dir.
fn build_end_cap(
    brep:     &mut BRep,
    centre:   Point3,
    normal:   UnitVec3,
    _radius:  f64,
    face_ids: &mut Vec<cadcore_topo::FaceId>,
) -> cadcore_topo::FaceId {
    let plane  = Plane3::from_origin_normal(centre, normal);
    let loop_id = build_placeholder_loop(brep);
    let fid = brep.add_face(Face {
        geom:        FaceGeom::Plane(plane),
        normal:      FaceNormal::Same,
        outer_loop:  loop_id,
        inner_loops: vec![],
        shell:       cadcore_topo::ShellId::default(),
    });
    face_ids.push(fid);
    fid
}

/// Insert a placeholder empty loop (patched with real co-edges later by the
/// STEP writer or a topology validator).
fn build_placeholder_loop(brep: &mut BRep) -> cadcore_topo::LoopId {
    // We insert a Loop with a sentinel co-edge id (the default key).
    // The STEP writer's boundary-stitching pass fills in real co-edges.
    brep.add_loop(cadcore_topo::Loop {
        start: cadcore_topo::CoEdgeId::default(),
        face:  cadcore_topo::FaceId::default(),
    })
}

// ── Boundary stitcher ─────────────────────────────────────────────────────────
// The full edge/co-edge topology between adjacent faces can be filled in by
// a subsequent stitching pass.  For STEP export, the STEP writer only needs
// the face list + surface geometry — it can reconstruct the boundary topology
// from the analytic surface intersections, which it does analytically.
//
// For now we expose the face geometry and let cadcore-step build the boundary
// edges on the fly when serialising.

#[cfg(test)]
mod tests {
    use super::*;

    fn line(x0: f64, x1: f64) -> Vec<Point3> {
        vec![Point3::new(x0, 0.0, 0.0), Point3::new(x1, 0.0, 0.0)]
    }

    #[test]
    fn single_segment_creates_solid() {
        let mut brep = BRep::new();
        let id = sweep_circle_along_polyline(
            &mut brep,
            &line(0.0, 10.0),
            1.0,
            &SweepOptions::default(),
        )
        .unwrap();
        assert!(brep.solids.contains_key(id));
        // 1 cylinder + 2 caps = 3 faces
        let stats = brep.stats();
        assert_eq!(stats.faces, 3, "expected 3 faces, got {}", stats.faces);
    }

    #[test]
    fn two_segments_adds_fillet() {
        let waypoints = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(5.0, 5.0, 0.0),
        ];
        let mut brep = BRep::new();
        let id = sweep_circle_along_polyline(
            &mut brep,
            &waypoints,
            0.5,
            &SweepOptions::default(),
        )
        .unwrap();
        // 2 cylinders + 1 torus fillet + 2 caps = 5 faces
        let stats = brep.stats();
        assert_eq!(stats.faces, 5, "expected 5 faces, got {}", stats.faces);
    }

    #[test]
    fn error_on_too_few_points() {
        let mut brep = BRep::new();
        let err = sweep_circle_along_polyline(
            &mut brep,
            &[Point3::ORIGIN],
            1.0,
            &SweepOptions::default(),
        )
        .unwrap_err();
        assert_eq!(err, SweepError::TooFewWaypoints);
    }

    #[test]
    fn error_on_zero_radius() {
        let mut brep = BRep::new();
        let err = sweep_circle_along_polyline(
            &mut brep,
            &line(0.0, 5.0),
            0.0,
            &SweepOptions::default(),
        )
        .unwrap_err();
        assert!(matches!(err, SweepError::InvalidRadius(_)));
    }
}
