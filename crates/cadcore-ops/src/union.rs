//! Analytic-surface B-Rep Boolean **UNION** for DIW woodpile scaffolds.
//!
//! # Goal
//!
//! Fuse the independent, mutually-overlapping analytic sweep solids produced by
//! the filament sweep pipeline into a **single watertight manifold B-Rep solid**,
//! so downstream CAD tools (ANSYS SpaceClaim) never need to run "Combine".
//!
//! This is **not** a generic OCCT-style boolean kernel.  It targets exactly the
//! scaffold geometry: axis-aligned, **equal-radius cylinders** with **planar end
//! caps**.  Sphere caps and torus corner fillets are out of v1 scope and are
//! preserved untouched (their filaments are merged into the result shell but not
//! geometrically trimmed).
//!
//! **Electrodes are never unioned** — any solid named `electrode*` is left as an
//! independent solid.
//!
//! # Method
//!
//! For each cylindrical face, every crossing cylinder removes a *window* (the
//! region of this cylinder that lies inside the other).  The window boundary is
//! the exact cyl∩cyl intersection loop ([`cadcore_geom::cyl_cyl_intersection`]),
//! added as an inner hole loop.  The same loop is shared (opposite sense) by the
//! two crossing faces, so the result is watertight.  The cylinder's original end
//! boundaries (circle / miter ellipse) are reproduced so they keep sharing with
//! the untouched caps / neighbouring segments.

use std::f64::consts::PI;

use cadcore_geom::{
    cyl_cyl_intersection, cyl_plane_intersection, orient_loop_as_hole,
    surface_surface_intersection, swept_swept_intersection, CenterlineSeg, Circle3, CylPlaneCurve,
    CylSurf, Ellipse3, IntersectionPolyline, Line3, ParamSurface, Plane3, SsiOptions,
    SweptTubeSurface,
    TorusSurf,
};
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{
    BRep, CoEdge, CoEdgeId, CoEdgeSense, Edge, EdgeGeom, EdgeId, Face, FaceBoundary, FaceExtent,
    FaceGeom, FaceId, FaceNormal, Loop, LoopId, Shell, ShellId, Solid, SolidId, Vertex, VertexId,
};

const SAMPLES: usize = 64;
const OPEN_REJECT_LOG_LIMIT: usize = 16;

/// Configurable options for the analytic Union operation.
#[derive(Clone, Copy, Debug)]
pub struct UnionOptions {
    pub linear_tolerance: f64,
    pub angular_tolerance: f64,
    pub min_edge_length: f64,
    pub min_face_area: f64,
    pub preserve_analytic_edges: bool,
    pub enable_tolerant_vertex_fusion: bool,
    pub enable_tolerant_edge_fusion: bool,
}

impl Default for UnionOptions {
    fn default() -> Self {
        Self {
            linear_tolerance: 1e-4,       // 0.1 micron
            angular_tolerance: 1e-6,      // rad
            min_edge_length: 1e-4,        // mm
            min_face_area: 1e-8,          // mm²
            preserve_analytic_edges: true,
            enable_tolerant_vertex_fusion: true,
            enable_tolerant_edge_fusion: true,
        }
    }
}

/// Diagnostics from one analytic union run.
#[derive(Clone, Copy, Debug, Default)]
pub struct UnionReport {
    /// Number of solids after the union step.
    pub solids: usize,
    /// Broad-phase elbow-leg candidate pairs.
    pub elbow_broad_phase: usize,
    /// Elbow-leg pairs whose sampled surfaces really overlap.
    pub elbow_real_overlaps: usize,
    /// Real-overlap pairs for which SSI returned at least one bounded curve.
    pub elbow_overlap_found: usize,
    /// SSI curves returned for elbow-leg pairs.
    pub elbow_curves: usize,
    /// Open-to-boundary SSI curves.
    pub elbow_open_curves: usize,
    /// Open-to-boundary curve groups converted into trim topology.
    pub elbow_open_trims: usize,
    /// Closed elbow-leg trim loops converted into topology.
    pub elbow_closed_trims: usize,
    /// Straight-cylinder crossing trim loops converted into topology.
    pub cylinder_cross_trims: usize,
    /// SSI tracer fallback invocations.
    pub marching_fallbacks: usize,
}

/// One cylindrical face slated for trimming.
struct CylFace {
    face_id: FaceId,
    surf: CylSurf,
    length: f64,
    start: FaceBoundary,
    end: FaceBoundary,
    /// Owning filament solid name (`"path_N"`), to skip self-crossings.
    solid_name: Option<String>,
}

/// A torus *fillet* face (a curved U-turn elbow) slated for trimming.  Carries
/// the analytic torus + the arc range + its two minor end circles.
struct TorusFace {
    face_id: FaceId,
    surf: TorusSurf,
    theta_lo: f64,
    theta_hi: f64,
    start_circle: Circle3,
    end_circle: Circle3,
}

/// A hole (window) carved into a cylinder by a crossing cylinder, shared with
/// the partner face.  Built once per crossing pair; `reversed` makes the two
/// faces traverse the same two edges in opposite directions.
#[derive(Clone, Copy)]
struct HoleRef {
    e0: EdgeId,
    e1: EdgeId,
    reversed: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum BoundarySide {
    Start,
    End,
}

#[derive(Clone)]
struct OpenNotchRef {
    /// Trim-curve edges of this notch in chain order (stored directions chain
    /// head-to-tail; a single edge in the one-sided cases).
    edges: Vec<EdgeId>,
    side: BoundarySide,
    reversed: bool,
    /// Canonical junction circle to sample the gap arc on, so BOTH faces sharing
    /// this junction produce the identical (point-shared) gap polyline.
    gap_circle: Option<Circle3>,
    /// `true` ⇒ the bitten interval runs CCW from the notch edge's v_start
    /// angle to its v_end angle on the junction circle.  Computed ONCE in
    /// Pass 2c by sampling the crossing tube, so the elbow's notched wire and
    /// the connector strips agree on which side is clear (identical gap arcs).
    bite_ccw_from_start: bool,
}

/// Which side of an open-boundary split to keep.
///
/// `StartToEnd` keeps the region bounded by the outer-loop path from the
/// cutter edge start vertex to its end vertex plus the cutter edge traversed
/// backward. `EndToStart` keeps the complementary region.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OpenBoundaryKeep {
    StartToEnd,
    EndToStart,
}

/// Topology created by one open-boundary face split.
#[derive(Clone, Copy, Debug)]
pub struct OpenBoundarySplit {
    pub face_id: FaceId,
    pub loop_id: LoopId,
    pub open_coedge: CoEdgeId,
    pub boundary_coedges: usize,
    pub dropped_boundary_coedges: usize,
}

/// Fuse every scaffold filament solid in `brep` into a single analytic union
/// solid.  Electrode solids (and any solid with no cylindrical face) are left
/// untouched.  Returns the number of solids remaining.
pub fn union_solids(brep: &mut BRep) -> usize {
    union_solids_with_report(brep).solids
}

/// Same as [`union_solids`], with diagnostics for export/UI reporting.
pub fn union_solids_with_report(brep: &mut BRep) -> UnionReport {
    union_solids_with_centerlines(brep, &[])
}

/// Continuous-filament union.  `filaments` carries each swept filament's FULL
/// centre-line keyed by its solid name (`"path_N"`), so a crossing leg can be
/// intersected with the WHOLE continuous swept tube of a filament instead of one
/// segment at a time — turning the artificial open curves at segment junctions
/// into proper closed bite loops.  Pass `&[]` to disable (per-segment behaviour).
pub fn union_solids_with_centerlines(
    brep: &mut BRep,
    filaments: &[(String, Vec<crate::sweep::SweepPathSegment>)],
) -> UnionReport {
    union_solids_with_centerlines_and_options(brep, filaments, &UnionOptions::default())
}

/// Same as [`union_solids`], with custom configuration options.
pub fn union_solids_with_options(
    brep: &mut BRep,
    options: &UnionOptions,
) -> usize {
    union_solids_with_centerlines_and_options(brep, &[], options).solids
}

/// Same as [`union_solids_with_centerlines`], with custom configuration options.
pub fn union_solids_with_centerlines_and_options(
    brep: &mut BRep,
    filaments: &[(String, Vec<crate::sweep::SweepPathSegment>)],
    options: &UnionOptions,
) -> UnionReport {
    let tolerance = options.linear_tolerance;
    let mut vertex_map = VertexMap::new(brep, tolerance);
    let mut edge_map = EdgeMap::new(brep, tolerance);

    // ── Pass 1: classify solids and collect filament faces ────────────────────
    let mut cyl_faces: Vec<CylFace> = Vec::new();
    let mut torus_faces: Vec<TorusFace> = Vec::new();
    let mut keep_faces: Vec<FaceId> = Vec::new(); // caps / etc. to re-home
    let mut old_shells: Vec<ShellId> = Vec::new();
    let mut old_solids: Vec<SolidId> = Vec::new();

    for (sid, solid) in brep.solids.iter() {
        if is_electrode(solid) {
            continue; // electrodes never participate in the union
        }
        // A filament solid has at least one cylindrical face.
        let is_filament = solid.shells.iter().any(|sh| {
            brep.shells.get(*sh).map_or(false, |s| {
                s.faces.iter().any(|f| {
                    matches!(
                        brep.faces.get(*f).map(|x| &x.geom),
                        Some(FaceGeom::Cylinder(_))
                    )
                })
            })
        });
        if !is_filament {
            continue;
        }
        old_solids.push(sid);
        let solid_name = solid.name.clone();
        for &sh in &solid.shells {
            old_shells.push(sh);
            let Some(shell) = brep.shells.get(sh) else {
                continue;
            };
            for &fid in &shell.faces {
                let Some(face) = brep.faces.get(fid) else {
                    continue;
                };
                match (&face.geom, &face.extent) {
                    (FaceGeom::Cylinder(c), FaceExtent::Cylinder { length, start, end }) => {
                        cyl_faces.push(CylFace {
                            face_id: fid,
                            surf: *c,
                            length: *length,
                            start: start.clone(),
                            end: end.clone(),
                            solid_name: solid_name.clone(),
                        });
                    }
                    (
                        FaceGeom::Torus(t),
                        FaceExtent::TorusFillet {
                            start_circle,
                            end_circle,
                        },
                    ) => {
                        let theta = |c: &Circle3| -> f64 {
                            let d = c.frame.origin - t.frame.origin;
                            t.frame.y.dot_vec(d).atan2(t.frame.x.dot_vec(d))
                        };
                        let lo = theta(start_circle);
                        let mut d = theta(end_circle) - lo;
                        let pi = std::f64::consts::PI;
                        while d > pi {
                            d -= 2.0 * pi;
                        }
                        while d < -pi {
                            d += 2.0 * pi;
                        }
                        torus_faces.push(TorusFace {
                            face_id: fid,
                            surf: *t,
                            theta_lo: lo,
                            theta_hi: lo + d,
                            start_circle: *start_circle,
                            end_circle: *end_circle,
                        });
                    }
                    _ => keep_faces.push(fid),
                }
            }
        }
    }

    if cyl_faces.len() < 2 {
        return UnionReport {
            solids: brep.solids.len(),
            ..UnionReport::default()
        }; // nothing crosses → leave the model untouched
    }

    // ── Continuous-filament tubes (Step 3b): one analytic swept tube per
    //    filament, so a crossing leg meets the WHOLE filament (closed loops)
    //    rather than a single segment (artificial open curves at junctions). ──
    let filament_radius = params_radius_or(&cyl_faces);
    // (name, tube, has_arc) — has_arc marks elbow-bearing filaments.
    let filament_tubes: Vec<(String, SweptTubeSurface, bool)> = filaments
        .iter()
        .filter_map(|(name, segs)| {
            let has_arc = segs
                .iter()
                .any(|s| matches!(s, crate::sweep::SweepPathSegment::Arc { .. }));
            swept_tube_from_segs(segs, filament_radius).map(|t| (name.clone(), t, has_arc))
        })
        .collect();
    if std::env::var("CADCORE_DUMP_CONT").is_ok() && !filament_tubes.is_empty() {
        // Sanity: every elbow's start-circle point lies ON exactly one tube.
        let mut matched = 0usize;
        for tf in &torus_faces {
            let p = tf.start_circle.point_at(0.0);
            if filament_tubes
                .iter()
                .any(|(_, t, _)| t.signed_distance(p).map_or(false, |d| d.abs() < 1e-3))
            {
                matched += 1;
            }
        }
        // Count CLOSED bite loops the continuous approach finds on elbow
        // filaments × every leg (the per-segment pass got these OPEN → 0 trims).
        let ssi_opts = SsiOptions {
            step: (filament_radius * 0.25).max(0.01),
            ..Default::default()
        };
        let mut cont_closed = 0usize;
        for (_, tube, has_arc) in &filament_tubes {
            if !has_arc {
                continue;
            }
            for cf in &cyl_faces {
                let leg = swept_from_cyl(cf);
                for c in surface_surface_intersection(&leg, tube, &ssi_opts) {
                    if c.closed && c.points.len() >= 6 {
                        cont_closed += 1;
                    }
                }
            }
        }
        eprintln!(
            "[union][cont] tubes={} (r={:.4}) elbows_on_a_tube={}/{} continuous_closed_loops={}",
            filament_tubes.len(),
            filament_radius,
            matched,
            torus_faces.len(),
            cont_closed
        );
    }

    // ── Pass 2: pairwise cyl∩cyl windows (built once, shared by both faces) ────
    let mut holes: Vec<Vec<HoleRef>> = vec![Vec::new(); cyl_faces.len()];
    let mut torus_holes: Vec<Vec<HoleRef>> = vec![Vec::new(); torus_faces.len()];
    let mut open_notches: Vec<Vec<OpenNotchRef>> = vec![Vec::new(); cyl_faces.len()];
    let mut torus_open_notches: Vec<Vec<OpenNotchRef>> = vec![Vec::new(); torus_faces.len()];
    let mut multi_holes: Vec<Vec<MultiHoleRef>> = vec![Vec::new(); cyl_faces.len()];
    let mut torus_multi_holes: Vec<Vec<MultiHoleRef>> = vec![Vec::new(); torus_faces.len()];
    let mut torus_through_pieces: Vec<Vec<ThroughPiece>> = vec![Vec::new(); torus_faces.len()];
    // Butt-end (capped) trims: collected in Pass 2c, processed in Pass 2e.
    let mut cap_bites: Vec<(usize, BoundarySide, usize, Vec<Point3>)> = Vec::new();
    // Deterministic iteration order: a randomized-`HashMap` order here makes the
    // near-tangent vertex/edge merge (`get_or_create_vertex`) resolve marginal
    // dense-fill junctions differently every run, flipping otherwise-identical
    // layers between manifold and open.  `BTreeMap` iterates by key so the union
    // is reproducible (a prerequisite for regression-testing the tangent fixes).
    let mut cap_notches: std::collections::BTreeMap<FaceId, Vec<OpenNotchRef>> =
        std::collections::BTreeMap::new();
    // (leg, side) -> planar Disk cap face at that free end.
    let cap_map: std::collections::BTreeMap<(usize, BoundarySide), FaceId> = {
        let mut m = std::collections::BTreeMap::new();
        for &fid in &keep_faces {
            let Some(face) = brep.faces.get(fid) else { continue };
            let (FaceGeom::Plane(p), FaceExtent::Disk { radius }) = (&face.geom, &face.extent)
            else {
                continue;
            };
            for (ci, cf) in cyl_faces.iter().enumerate() {
                for (side, b) in [(BoundarySide::Start, &cf.start), (BoundarySide::End, &cf.end)]
                {
                    if let FaceBoundary::Circle(c) = b {
                        if (c.frame.origin - p.frame.origin).length() < 1e-4
                            && (c.radius - radius).abs() < 1e-4
                        {
                            m.insert((ci, side), fid);
                        }
                    }
                }
            }
        }
        if std::env::var("CADCORE_DUMP_CAP").is_ok() {
            let n_disk = keep_faces
                .iter()
                .filter(|&&fid| {
                    brep.faces.get(fid).map_or(false, |f| {
                        matches!(f.extent, FaceExtent::Disk { .. })
                    })
                })
                .count();
            let mut keys: Vec<(usize, BoundarySide)> = m.keys().copied().collect();
            keys.sort_by_key(|k| k.0);
            eprintln!("[cap-map] entries={} disk_caps={} keep_faces={} keys={keys:?}", m.len(), n_disk, keep_faces.len());
        }
        m
    };
    let mut collars: Vec<Vec<CollarRef>> = vec![Vec::new(); cyl_faces.len()];
    let mut through_pieces: Vec<Vec<ThroughPiece>> = vec![Vec::new(); cyl_faces.len()];
    let mut bite_counter = 0usize;
    // Faces per filament (for the N-face spanning-loop split): legs by solid
    // name; elbows attached transitively via shared junction circles.
    let legs_of: std::collections::HashMap<String, Vec<usize>> = {
        let mut m: std::collections::HashMap<String, Vec<usize>> = std::collections::HashMap::new();
        for (ci, cf) in cyl_faces.iter().enumerate() {
            if let Some(nm) = &cf.solid_name {
                m.entry(nm.clone()).or_default().push(ci);
            }
        }
        m
    };
    let elbows_for_legs = |legs: &[usize]| -> Vec<usize> {
        torus_faces
            .iter()
            .enumerate()
            .filter(|(_, tf)| {
                legs.iter().any(|&ci| {
                    let cf = &cyl_faces[ci];
                    boundary_circle_matches(&cf.start, &tf.start_circle)
                        || boundary_circle_matches(&cf.start, &tf.end_circle)
                        || boundary_circle_matches(&cf.end, &tf.start_circle)
                        || boundary_circle_matches(&cf.end, &tf.end_circle)
                })
            })
            .map(|(i, _)| i)
            .collect()
    };
    let mut n_trimmed = 0usize;
    let mut n_dropped = 0usize;
    let mut n_should_cross = 0usize;
    for i in 0..cyl_faces.len() {
        for j in (i + 1)..cyl_faces.len() {
            let a = &cyl_faces[i];
            let b = &cyl_faces[j];
            // Diagnostic: perpendicular, equal-radius, axis distance < 2r.
            {
                let da = a.surf.axis();
                let db = b.surf.axis();
                let perp = da.dot(db).abs() < 0.999; // not near-parallel
                let eqr = (a.surf.radius - b.surf.radius).abs() < 1e-9;
                let w = b.surf.frame.origin - a.surf.frame.origin;
                // distance between the two (perpendicular) axes
                let along_a = da.dot_vec(w);
                let along_b = db.dot_vec(w);
                let perp_vec = w - da.as_vec() * along_a - db.as_vec() * along_b;
                let axis_dist = perp_vec.length();
                if perp && eqr && axis_dist < 2.0 * a.surf.radius {
                    n_should_cross += 1;
                }
            }
            for loop_pts in cyl_cyl_intersection(&a.surf, &b.surf, SAMPLES) {
                if !loop_pts.closed || loop_pts.points.len() < 6 {
                    continue;
                }
                // Both segments must contain the loop (interior crossing).
                if !within_axial(&a.surf, a.length, &loop_pts.points)
                    || !within_axial(&b.surf, b.length, &loop_pts.points)
                {
                    n_dropped += 1;
                    continue;
                }
                n_trimmed += 1;
                // Wind the loop as a proper hole (CW in cylinder i's (θ,v) param,
                // i.e. opposite the outer bound) so it REMOVES the patch of i
                // inside j.  Cylinder j (reversed) then gets the correct opposite
                // winding for free, which also satisfies the shared-edge rule.
                let oriented = orient_hole_for_cylinder(&a.surf, &loop_pts.points);
                let (e0, e1) = build_hole_edges(brep, &oriented, &mut vertex_map, &mut edge_map, tolerance, "pass2_cylcyl");
                holes[i].push(HoleRef {
                    e0,
                    e1,
                    reversed: false,
                });
                holes[j].push(HoleRef {
                    e0,
                    e1,
                    reversed: true,
                });
            }
        }
    }

    // ── Pass 2b: elbow(curved swept) ∩ leg(cylinder) via the GENERIC SSI ──────
    // The curved U-turn elbows at the perimeter cross adjacent-layer legs.  Both
    // are swept tubes; the generic numeric surface∩surface engine handles the
    // straight×curved case (the analytic cyl∩cyl above is just the fast path for
    // straight×straight).  Bite loop → hole on the leg + the SAME loop as a hole
    // on the elbow (shared edges).  Each face's hole is wound CW on its own
    // surface via `orient_loop_as_hole`.
    let ssi_opts = SsiOptions {
        step: (params_radius_or(&cyl_faces) * 0.25).max(0.01),
        ..Default::default()
    };
    // DIAGNOSTIC: arc-angle distribution of the elbows.
    if !torus_faces.is_empty() {
        let mut amin = f64::MAX;
        let mut amax = f64::MIN;
        let mut asum = 0.0;
        for tf in &torus_faces {
            let a = (tf.theta_hi - tf.theta_lo).abs().to_degrees();
            amin = amin.min(a);
            amax = amax.max(a);
            asum += a;
        }
        eprintln!(
            "[union][diag] elbow arc° min={:.1} max={:.1} avg={:.1}",
            amin,
            amax,
            asum / torus_faces.len() as f64
        );
    }
    // DIAGNOSTIC: minimum REAL elbow↔leg surface gap over broad-phase pairs
    // (negative = surfaces actually overlap = a true crossing exists).
    let mut min_real_gap = f64::MAX;

    let mut n_torus_trim = 0usize;
    let mut n_bp = 0usize; // broad-phase candidates
    let mut n_curves = 0usize; // SSI returned a bounded curve
    let mut n_open_boundary = 0usize;
    let mut n_open_converted = 0usize;
    let mut n_marching_fallback = 0usize;
    let mut n_axial_drop = 0usize;
    let mut n_overlap = 0usize; // pairs that truly overlap (probe)
    let mut n_overlap_found = 0usize; // ... and SSI found a curve
    let mut open_diag_logged = 0usize;
    let mut n_open_reject_not_same_boundary = 0usize;
    let mut pair_seq = 0usize;
    for ti in 0..torus_faces.len() {
        let elbow = swept_from_torus(&torus_faces[ti]);
        for ci in 0..cyl_faces.len() {
            pair_seq += 1;
            let cf = &cyl_faces[ci];
            let tf = &torus_faces[ti];
            // Broad-phase: leg axis must pass within R+2r of the elbow centre.
            let reach = tf.surf.major_radius + 2.0 * tf.surf.minor_radius;
            let w = tf.surf.frame.origin - cf.surf.frame.origin;
            let along = cf.surf.axis().dot_vec(w);
            let perp = (w - cf.surf.axis().as_vec() * along).length();
            if perp > reach || along < -reach || along > cf.length + reach {
                continue;
            }
            n_bp += 1;
            // Skip G1-ADJACENT pairs (this elbow CONTINUES this leg — they
            // share a junction circle): their "intersection" is the
            // degenerate tangent contact along the junction, which the SSI
            // tracer turns into a phantom collar-like loop wrapped around the
            // leg start (vertices at the tube's top/bottom on the junction
            // plane).  Wired as a hole it winds in (u,v) → importers read it
            // as an island and fragment both faces.
            {
                let tfx = &torus_faces[ti];
                if boundary_circle_matches(&cf.start, &tfx.start_circle)
                    || boundary_circle_matches(&cf.start, &tfx.end_circle)
                    || boundary_circle_matches(&cf.end, &tfx.start_circle)
                    || boundary_circle_matches(&cf.end, &tfx.end_circle)
                {
                    continue;
                }
            }
            // Per-pair real-gap probe: does the elbow actually overlap the leg?
            let mut pair_gap = f64::MAX;
            {
                let r_leg = cf.surf.radius;
                let lax = cf.surf.axis();
                let lo = cf.surf.frame.origin;
                let (u0, u1) = elbow.u_domain();
                for iu in 0..=16 {
                    for iv in 0..16 {
                        let uu = u0 + (u1 - u0) * iu as f64 / 16.0;
                        let vv = std::f64::consts::PI * 2.0 * iv as f64 / 16.0;
                        let p = elbow.point(uu, vv);
                        let w = p - lo;
                        let ax = lax.dot_vec(w);
                        if ax < 0.0 || ax > cf.length {
                            continue;
                        }
                        let radial = (w - lax.as_vec() * ax).length();
                        pair_gap = pair_gap.min(radial - r_leg);
                    }
                }
            }
            min_real_gap = min_real_gap.min(pair_gap);
            let truly_overlaps = pair_gap < -0.02;
            if truly_overlaps {
                n_overlap += 1;
                if n_overlap == 1 {
                    let t = &torus_faces[ti];
                    eprintln!(
                        "[union][pair0] elbow c={:?} z={:?} x={:?} R={:.4} r={:.4} th=[{:.4},{:.4}] | leg o={:?} ax={:?} len={:.4} r={:.4} gap={:.4}",
                        t.surf.frame.origin, t.surf.frame.z.as_vec(), t.surf.frame.x.as_vec(),
                        t.surf.major_radius, t.surf.minor_radius, t.theta_lo, t.theta_hi,
                        cf.surf.frame.origin, cf.surf.axis().as_vec(), cf.length, cf.surf.radius,
                        pair_gap
                    );
                }
            }
            let leg = swept_from_cyl(cf);
            let mut found_here = false;
            let traced = surface_surface_intersection(&elbow, &leg, &ssi_opts);
            let tracer_failed = traced.is_empty()
                || traced
                    .iter()
                    .any(|c| c.points.len() >= ssi_opts.max_pts || c.points.len() < 6 || !c.closed);
            let curves = if tracer_failed {
                n_marching_fallback += 1;
                swept_swept_intersection(&elbow, &leg, 96, 48)
            } else {
                traced
            };
            let mut open_curves: Vec<IntersectionPolyline> = Vec::new();
            for curve in curves {
                if curve.points.len() < 6 {
                    continue;
                }
                n_curves += 1;
                found_here = true;
                if !curve.closed {
                    n_open_boundary += 1;
                    if open_diag_logged < OPEN_REJECT_LOG_LIMIT {
                        log_open_curve_diag(
                            "queued_open",
                            pair_seq,
                            ti,
                            ci,
                            tf.face_id,
                            cf.face_id,
                            &elbow,
                            &leg,
                            &curve,
                        );
                        open_diag_logged += 1;
                    }
                    // Pass 2b open-notch conversion DISABLED: the per-segment
                    // SSI yields artificially-open curves at tangent
                    // elbow↔leg junctions; wiring them as boundary notches
                    // (without a canonical junction circle) produced
                    // self-intersecting leg wires.  The continuous-filament
                    // junction split (Pass 2c, with the on-J endpoint guard)
                    // is the correct mechanism for these crossings.
                    if false {
                        if let Some((edge, elbow_side, leg_side)) =
                            build_open_notch_if_same_boundary(&curve, &elbow, &leg, brep, &mut vertex_map, &mut edge_map, tolerance)
                        {
                            open_notches[ci].push(OpenNotchRef {
                                edges: vec![edge],
                                side: leg_side,
                                reversed: false,
                                gap_circle: None,
                                bite_ccw_from_start: false,
                            });
                            torus_open_notches[ti].push(OpenNotchRef {
                                edges: vec![edge],
                                side: elbow_side,
                                reversed: true,
                                gap_circle: None,
                                bite_ccw_from_start: false,
                            });
                            n_open_converted += 1;
                            n_torus_trim += 1;
                            continue;
                        }
                    }
                    if open_diag_logged < OPEN_REJECT_LOG_LIMIT {
                        n_open_reject_not_same_boundary += 1;
                        eprintln!(
                            "[union][open:reject] pair={} ti={} ci={} reason=not_same_boundary_on_both_faces_or_snap_failed",
                            pair_seq, ti, ci
                        );
                        open_diag_logged += 1;
                    } else {
                        n_open_reject_not_same_boundary += 1;
                    }
                    open_curves.push(curve);
                    continue;
                }
                if !within_axial(&cf.surf, cf.length, &curve.points) {
                    n_axial_drop += 1;
                    continue;
                }
                // A REAL bite loop lies on BOTH surfaces.  The numeric tracer
                // can degenerate near the parallel/tangent zone (U-connector
                // over a crossing leg) and emit a zero-area "ridge" loop that
                // sits on the leg but nowhere near the elbow — wiring it as a
                // hole fragments the leg face in importers.
                if !loop_on_tube_surface(&leg, &curve.points, 1e-3)
                    || !loop_on_tube_surface(&elbow, &curve.points, 1e-3)
                {
                    n_axial_drop += 1;
                    continue;
                }
                // Wind CW on the leg; the elbow takes the same edges reversed
                // (opposite outward normals at the bite ⇒ that is CW on the elbow).
                let oriented = orient_loop_as_hole(&leg, &curve.points);
                let (e0, e1) = build_hole_edges(brep, &oriented, &mut vertex_map, &mut edge_map, tolerance, "pass2b_closed");
                holes[ci].push(HoleRef {
                    e0,
                    e1,
                    reversed: false,
                });
                torus_holes[ti].push(HoleRef {
                    e0,
                    e1,
                    reversed: true,
                });
                n_torus_trim += 1;
            }
            for loop_pts in open_boundary_loops(&torus_faces[ti], &elbow, open_curves) {
                if !within_axial(&cf.surf, cf.length, &loop_pts) {
                    n_axial_drop += 1;
                    continue;
                }
                if !loop_on_tube_surface(&leg, &loop_pts, 1e-4)
                    || !loop_on_tube_surface(&elbow, &loop_pts, 1e-3)
                {
                    continue;
                }
                let oriented = orient_loop_as_hole(&leg, &loop_pts);
                let (e0, e1) = build_hole_edges(brep, &oriented, &mut vertex_map, &mut edge_map, tolerance, "pass2b_open");
                holes[ci].push(HoleRef {
                    e0,
                    e1,
                    reversed: false,
                });
                torus_holes[ti].push(HoleRef {
                    e0,
                    e1,
                    reversed: true,
                });
                n_torus_trim += 1;
                n_open_converted += 1;
            }
            if truly_overlaps && found_here {
                n_overlap_found += 1;
            }
        }
    }

    // Discard the OLD Pass 2b open-notch attempts (spurious tangent junctions =
    // the 166-violation case).  Only the continuous-SSI junction split (Pass 2c)
    // below feeds correct notches.
    for v in open_notches.iter_mut() {
        v.clear();
    }
    for v in torus_open_notches.iter_mut() {
        v.clear();
    }

    // ── Pass 2c: CONTINUOUS-filament elbow trim ───────────────────────────────
    // Intersect each crossing leg with the WHOLE filament tube (closed loops),
    // route each loop to the elbow face it lies on, and add it as a shared hole.
    // For junction-spanning loops the crossing leg chooses `HoleRef::reversed`
    // from its own UV winding; the elbow / adjacent leg notches then take the
    // opposite sense to keep the shared edges manifold-clean.
    let cont_ssi_opts = SsiOptions {
        step: (filament_radius * 0.25).max(0.01),
        ..Default::default()
    };
    let mut n_cont_trim = 0usize;
    let mut n_cont_closed = 0usize;
    let mut n_cont_on_elbow = 0usize;
    let mut n_cont_spans = 0usize; // touches an elbow but not entirely on one
                                   // DEFAULT ON since the full chain went externally valid (FreeCAD/OCC:
                                   // fused 15×15×12 = 1 closed solid, 0 invalid faces of 712): N-face
                                   // junction split + connector strips + derived senses + wire-chained
                                   // pcurves + G1-adjacent-pair skip.  Set CADCORE_ELBOW_TRIM=0 to
                                   // disable for debugging.
    let elbow_trim_enabled = std::env::var("CADCORE_ELBOW_TRIM").map_or(true, |v| v != "0");
    for (fname, tube, has_arc) in &filament_tubes {
        if !has_arc || !elbow_trim_enabled {
            continue;
        }
        for ci in 0..cyl_faces.len() {
            // Skip legs that belong to THIS filament (tangent self-junctions).
            if cyl_faces[ci].solid_name.as_deref() == Some(fname.as_str()) {
                continue;
            }
            let leg = swept_from_cyl(&cyl_faces[ci]);
            for loop_pts in surface_surface_intersection(&leg, tube, &cont_ssi_opts) {
                if !loop_pts.closed || loop_pts.points.len() < 6 {
                    // OPEN curve fragments are stitched in Pass 2e (the
                    // tracer splits cap-plane bites into pieces).
                    if !loop_pts.closed && loop_pts.points.len() >= 2 {
                        if let Some(fidx) =
                            filament_tubes.iter().position(|(t, _, _)| t == fname)
                        {
                            cap_bites.push((ci, BoundarySide::Start, fidx, loop_pts.points.clone()));
                        }
                    }
                    continue;
                }
                if !within_axial(&cyl_faces[ci].surf, cyl_faces[ci].length, &loop_pts.points) {
                    continue;
                }
                // PHANTOM-LOOP guard: where a straight section of F runs
                // PARALLEL to this leg (U-connector over the crossing leg at
                // the apex, axis gap < r_L + r_F), the surfaces are
                // near-tangent along the leg's "ridge" line facing F.  The
                // numeric SSI tracer degenerates there and stitches fake
                // loops across the ridge (the real intersection is the pair
                // of analytic parallel lines).  A run of chain points hugging
                // a ridge ⇒ reject the loop (leave untrimmed).
                let phantom = {
                    let lsurf = &cyl_faces[ci].surf;
                    let lax = lsurf.axis();
                    let mut ridges: Vec<Point3> = Vec::new(); // ridge line origins
                    if let Some((_, segs)) = filaments.iter().find(|(n, _)| n == fname) {
                        for s in segs {
                            if let crate::sweep::SweepPathSegment::Line { start, end } = s {
                                let d = *end - *start;
                                let len = d.length();
                                if len < 1e-9 {
                                    continue;
                                }
                                let dir = d * (1.0 / len);
                                if dir.dot(lax.as_vec()).abs() < 0.99 {
                                    continue; // not parallel to the leg
                                }
                                let w = *start - lsurf.frame.origin;
                                let perp = w - lax.as_vec() * lax.dot_vec(w);
                                let gap = perp.length();
                                if gap < 1e-9 || gap > 2.0 * lsurf.radius + filament_radius {
                                    continue;
                                }
                                // ridge = leg-surface line nearest F's segment
                                ridges.push(
                                    lsurf.frame.origin + perp * (lsurf.radius / gap),
                                );
                            }
                        }
                    }
                    if ridges.is_empty() {
                        false
                    } else {
                        let ridge_d = |p: Point3| -> f64 {
                            ridges
                                .iter()
                                .map(|&ro| {
                                    let w = p - ro;
                                    (w - lax.as_vec() * lax.dot_vec(w)).length()
                                })
                                .fold(f64::MAX, f64::min)
                        };
                        let mut run = 0usize;
                        let mut worst = 0usize;
                        let mut dmin = f64::MAX;
                        for &p in &loop_pts.points {
                            let d = ridge_d(p);
                            dmin = dmin.min(d);
                            if d < 0.03 {
                                run += 1;
                                worst = worst.max(run);
                            } else {
                                run = 0;
                            }
                        }
                        if std::env::var("CADCORE_DUMP_CONT").is_ok() && dmin < 0.2 {
                            eprintln!(
                                "[phantom-stat] ci={ci} ridges={} dmin={dmin:.4} worst_run={worst} pts={}",
                                ridges.len(),
                                loop_pts.points.len()
                            );
                        }
                        // NOTE: dmin-based rejection was tried and rejected:
                        // EVERY U-apex bite approaches the tangent zone
                        // (dmin≈0.02), so a distance threshold kills all 220
                        // legitimate trims.  Only a sustained RUN along the
                        // ridge marks a truly degenerate trace.
                        worst >= 3
                    }
                };
                if phantom {
                    if std::env::var("CADCORE_DUMP_CONT").is_ok() {
                        eprintln!("[phantom-reject] ci={ci} f={fname}");
                    }
                    continue;
                }
                n_cont_closed += 1;
                // Which elbow face does this loop sit on?  (Loops on the straight
                // parts are already handled by the cyl∩cyl Pass 2.)
                let on_elbow = torus_faces
                    .iter()
                    .position(|tf| loop_on_torus_face(tf, &loop_pts.points, 1e-2));
                if let Some(ti) = on_elbow {
                    // Entirely on one elbow → simple shared hole (leg ↔ elbow).
                    n_cont_on_elbow += 1;
                    let oriented = orient_loop_as_hole(&leg, &loop_pts.points);
                    let (e0, e1) = build_hole_edges(brep, &oriented, &mut vertex_map, &mut edge_map, tolerance, "pass2c_onelbow");
                    holes[ci].push(HoleRef {
                        e0,
                        e1,
                        reversed: false,
                    });
                    torus_holes[ti].push(HoleRef {
                        e0,
                        e1,
                        reversed: true,
                    });
                    n_cont_trim += 1;
                } else if let Some(pieces) = {
                    let legs = legs_of.get(fname.as_str()).cloned().unwrap_or_default();
                    let elbows = elbows_for_legs(&legs);
                    split_multi_span(&loop_pts.points, &elbows, &legs, &torus_faces, &cyl_faces)
                } {
                    // The loop spans several faces of the bitten filament
                    // (elbow / short connector / elbow at a U-turn apex).
                    // Distribute: crossing leg ← one N-edge hole; elbow and
                    // same-junction leg pieces ← boundary notches; fully
                    // crossed connectors ← a through-cut wire.  Commit only
                    // when EVERY part is representable.
                    n_cont_spans += 1;
                    let same_j = |p: &SpanPiece| {
                        (p.j_in.frame.origin - p.j_out.frame.origin).length() < 1e-3
                    };
                    // tangent / sliver guard
                    let mut ok = pieces.iter().all(|p| {
                        (p.pts[0] - *p.pts.last().unwrap()).length() > 0.01 || !same_j(p)
                    });
                    // elbow pieces must stay within one junction
                    ok &= pieces.iter().all(|p| match p.owner {
                        PieceOwner::Elbow(_) => same_j(p),
                        PieceOwner::Leg(_) => true,
                    });
                    // group through-cut pieces per leg face
                    // BTreeMap: deterministic iteration (see cap_notches note).
                    let mut tc_legs: std::collections::BTreeMap<usize, Vec<usize>> =
                        std::collections::BTreeMap::new();
                    for (i, p) in pieces.iter().enumerate() {
                        if let PieceOwner::Leg(lci) = p.owner {
                            if !same_j(p) {
                                tc_legs.entry(lci).or_default().push(i);
                            }
                        }
                    }
                    for (&lci, list) in &tc_legs {
                        // A through-cut connector may host MANY bites (one per
                        // crossing layer), but no other trim kinds.
                        ok &= list.len() == 2
                            && holes[lci].is_empty()
                            && open_notches[lci].is_empty();
                        if ok {
                            let a = &pieces[list[0]];
                            let b = &pieces[list[1]];
                            ok &= (a.j_in.frame.origin - b.j_out.frame.origin).length() < 1e-3
                                && (a.j_out.frame.origin - b.j_in.frame.origin).length() < 1e-3;
                            // Both connector boundaries must be plain circles.
                            ok &= matches!(cyl_faces[lci].start, FaceBoundary::Circle(_))
                                && matches!(cyl_faces[lci].end, FaceBoundary::Circle(_));
                        }
                    }
                    if !ok {
                        if std::env::var("CADCORE_DUMP_CONT").is_ok() {
                            let sliver = pieces
                                .iter()
                                .map(|p| {
                                    ((p.pts[0] - *p.pts.last().unwrap()).length() * 1000.0).round()
                                        / 1000.0
                                })
                                .collect::<Vec<_>>();
                            let ej = pieces
                                .iter()
                                .map(|p| same_j(p))
                                .collect::<Vec<_>>();
                            let tcinfo = tc_legs
                                .iter()
                                .map(|(lci, l)| {
                                    (
                                        *lci,
                                        l.len(),
                                        through_pieces[*lci].len(),
                                        holes[*lci].len(),
                                        multi_holes[*lci].len(),
                                        open_notches[*lci].len(),
                                    )
                                })
                                .collect::<Vec<_>>();
                            eprintln!(
                                "[union][cont:multi-skip] pieces={} owners={:?} endsep={:?} same_j={:?} tc={:?}",
                                pieces.len(),
                                pieces.iter().map(|p| p.owner).collect::<Vec<_>>(),
                                sliver,
                                ej,
                                tcinfo
                            );
                        }
                        continue;
                    }
                    // Build the shared piece edges.
                    let mut piece_edges: Vec<EdgeId> = Vec::with_capacity(pieces.len());
                    for p in &pieces {
                        piece_edges.push(build_single_edge(
                            brep, &p.pts, &mut vertex_map, &mut edge_map, tolerance,
                        ));
                    }
                    // Crossing-leg hole orientation: must be CW in its (θ,v).
                    let mut chain: Vec<Point3> = Vec::new();
                    for p in &pieces {
                        chain.extend(p.pts.iter().copied());
                    }
                    // Hole orientation by MATERIAL-ON-LEFT (not signed area —
                    // the folded apex chains wrap most of the leg's θ and the
                    // area sign becomes unreliable): walking the hole with the
                    // leg's kept material (outside the biting tube F) on the
                    // LEFT is the correct direction; `hole_rev` = stored chain
                    // order fails that test.
                    let material_left_ok = {
                        let s = &cyl_faces[ci].surf;
                        let mut ok = 0usize;
                        let mut n = 0usize;
                        let step = (chain.len() / 9).max(1);
                        for k in (0..chain.len().saturating_sub(1)).step_by(step) {
                            let p = chain[k];
                            let t = chain[k + 1] - p;
                            if t.length() < 1e-9 {
                                continue;
                            }
                            let w = p - s.frame.origin;
                            let ax = s.axis().dot_vec(w);
                            let nrm = (w - s.axis().as_vec() * ax).normalize();
                            let left = nrm.cross(t).normalize();
                            let probe = p + left * 0.05;
                            n += 1;
                            if tube.signed_distance(probe).map_or(true, |d| d >= 0.0) {
                                ok += 1; // material (outside F) on the left
                            }
                        }
                        n == 0 || ok * 2 >= n
                    };
                    let wu = chain_winding_theta(&cyl_faces[ci].surf, &chain);
                    let hole_rev = if wu.abs() > PI {
                        // WINDING collar (the U wraps the leg's circumference):
                        // not a hole — it will split the leg into bands.  The
                        // leg-side traversal: material-on-the-left ⇒ the KEPT
                        // axial side decides (+θ when kept above).  Sample
                        // both sides against the bitten filament's tube.
                        let s = &cyl_faces[ci].surf;
                        let mean_v = chain
                            .iter()
                            .map(|p| s.axis().dot_vec(*p - s.frame.origin))
                            .sum::<f64>()
                            / chain.len() as f64;
                        let probe = |dv: f64| -> usize {
                            (0..8)
                                .filter(|k| {
                                    let th = *k as f64 * PI / 4.0;
                                    let p = s.point_at(th, mean_v + dv);
                                    tube.signed_distance(p).map_or(false, |d| d < 0.0)
                                })
                                .count()
                        };
                        let r = s.radius;
                        let kept_above = probe(1.2 * r) <= probe(-1.2 * r);
                        let stored_ccw = wu > 0.0;
                        let leg_reversed = stored_ccw != kept_above;
                        collars[ci].push(CollarRef {
                            edges: piece_edges.clone(),
                            stored_ccw,
                            mean_v,
                            pts: chain.clone(),
                        });
                        leg_reversed
                    } else {
                        let rev = !material_left_ok;
                        multi_holes[ci].push(MultiHoleRef {
                            edges: piece_edges.clone(),
                            reversed: rev,
                        });
                        rev
                    };
                    // Bitten interval on a junction circle, sampled against
                    // the CROSSING tube (single source of truth for both
                    // sides of every junction): does the CCW arc from the
                    // piece's first point to its last contain the bite?
                    let bite_ccw_of = |p: &SpanPiece| -> bool {
                        let j = &p.j_in;
                        let a1 = circle_angle(j, p.pts[0]);
                        let a2 = circle_angle(j, *p.pts.last().unwrap());
                        let mid = a1 + (a2 - a1).rem_euclid(2.0 * PI) * 0.5;
                        let q = j.point_at(mid);
                        let s = &cyl_faces[ci].surf;
                        let w = q - s.frame.origin;
                        let ax = s.axis().dot_vec(w);
                        (w - s.axis().as_vec() * ax).length() < s.radius
                    };
                    // Notches (opposite-of-hole traversal).
                    for (pi, p) in pieces.iter().enumerate() {
                        match p.owner {
                            PieceOwner::Elbow(ti2) => {
                                let tf2 = &torus_faces[ti2];
                                let side = if (p.j_in.frame.origin
                                    - tf2.start_circle.frame.origin)
                                    .length()
                                    < 1e-3
                                {
                                    BoundarySide::Start
                                } else {
                                    BoundarySide::End
                                };
                                torus_open_notches[ti2].push(OpenNotchRef {
                                    edges: vec![piece_edges[pi]],
                                    side,
                                    reversed: !hole_rev,
                                    gap_circle: Some(p.j_in),
                                    bite_ccw_from_start: bite_ccw_of(p),
                                });
                            }
                            PieceOwner::Leg(lci) if same_j(p) => {
                                let cf2 = &cyl_faces[lci];
                                let side = if boundary_circle_matches(&cf2.start, &p.j_in) {
                                    BoundarySide::Start
                                } else {
                                    BoundarySide::End
                                };
                                if std::env::var("CADCORE_DUMP_CONT").is_ok() {
                                    let vr: Vec<f64> = p
                                        .pts
                                        .iter()
                                        .map(|&q| cf2.surf.axis().dot_vec(q - cf2.surf.frame.origin))
                                        .collect();
                                    let vmin = vr.iter().fold(f64::MAX, |a, &b| a.min(b));
                                    let vmax = vr.iter().fold(f64::MIN, |a, &b| a.max(b));
                                    let matches_end = boundary_circle_matches(&cf2.end, &p.j_in);
                                    eprintln!(
                                        "[notch-push] lci={lci} side={side:?} m_start={} m_end={} v=[{vmin:.3},{vmax:.3}] len={:.3} j_o={:?}",
                                        boundary_circle_matches(&cf2.start, &p.j_in),
                                        matches_end,
                                        cf2.length,
                                        p.j_in.frame.origin
                                    );
                                }
                                open_notches[lci].push(OpenNotchRef {
                                    edges: vec![piece_edges[pi]],
                                    side,
                                    reversed: !hole_rev,
                                    gap_circle: Some(p.j_in),
                                    bite_ccw_from_start: bite_ccw_of(p),
                                });
                            }
                            PieceOwner::Leg(_) => {} // through-cut, below
                        }
                    }
                    for (&lci, list) in &tc_legs {
                        let start_origin = match &cyl_faces[lci].start {
                            FaceBoundary::Circle(c) => c.frame.origin,
                            FaceBoundary::Ellipse(e) => e.frame.origin,
                        };
                        for &pi in list {
                            let p = &pieces[pi];
                            let first = p.pts[0];
                            let last = *p.pts.last().unwrap();
                            // Endpoint on the connector's START circle = p_j1.
                            let first_on_start = (p.j_in.frame.origin - start_origin).length()
                                < (p.j_out.frame.origin - start_origin).length();
                            let (p_j1, p_j2) = if first_on_start {
                                (first, last)
                            } else {
                                (last, first)
                            };
                            let crossing_tube = cyl_faces[ci]
                                .solid_name
                                .as_deref()
                                .and_then(|nm| {
                                    filament_tubes.iter().position(|(t, _, _)| t == nm)
                                })
                                .unwrap_or(usize::MAX);
                            through_pieces[lci].push(ThroughPiece {
                                edges: vec![piece_edges[pi]],
                                bite: bite_counter,
                                crossing: crossing_tube,
                                p_j1,
                                p_j2,
                                strip_same: hole_rev,
                            });
                        }
                    }
                    bite_counter += 1;
                    n_cont_trim += 1;
                }
            }
        }
    }
    if std::env::var("CADCORE_DUMP_CONT").is_ok() {
        eprintln!(
            "[union][cont] closed_loops={} on_elbow={} spans_junction={} trims_applied={}",
            n_cont_closed, n_cont_on_elbow, n_cont_spans, n_cont_trim
        );
    }
    let _ = (n_cont_closed, n_cont_on_elbow, n_cont_spans);

    // ── Pass 2d: ELBOW as the CROSSING side × foreign U-tube (corner U-on-U) ──
    // At the structure corners two perpendicular U-turns stack: the upper U's
    // ELBOW bites the lower U's elbow/connector.  Pass 2c only uses straight
    // cylinders as the crossing side; this pass mirrors it with a torus face
    // crossing a foreign filament tube, reusing the same multi-span
    // distribution (the bitten side) and an N-edge hole on the elbow.
    let corner_trim_enabled = std::env::var("CADCORE_CORNER_TRIM").map_or(true, |v| v != "0");
    let mut n_corner_trim = 0usize;
    let mut d_bp = 0usize;
    let mut d_loops = 0usize;
    let mut d_not_on_elbow = 0usize;
    let mut d_one_leg = 0usize;
    let mut d_pair = 0usize;
    let mut d_split_none = 0usize;
    let mut d_guard = 0usize;
    if elbow_trim_enabled && corner_trim_enabled {
        // elbow → owning filament tube (junction-circle adjacency).
        let elbow_tube: Vec<Option<usize>> = torus_faces
            .iter()
            .map(|tf| {
                filament_tubes.iter().position(|(nm, _, _)| {
                    legs_of.get(nm).map_or(false, |legs| {
                        legs.iter().any(|&ci| {
                            let c = &cyl_faces[ci];
                            boundary_circle_matches(&c.start, &tf.start_circle)
                                || boundary_circle_matches(&c.start, &tf.end_circle)
                                || boundary_circle_matches(&c.end, &tf.start_circle)
                                || boundary_circle_matches(&c.end, &tf.end_circle)
                        })
                    })
                })
            })
            .collect();
        // Side classification on a filament's faces.
        enum SideKind {
            OneLeg,
            OneElbow(usize),
            Multi,
        }
        let on_leg_all = |lci: usize, pts: &[Point3]| -> bool {
            let c = &cyl_faces[lci];
            pts.iter().all(|&p| {
                let w = p - c.surf.frame.origin;
                let ax = c.surf.axis().dot_vec(w);
                let rad = (w - c.surf.axis().as_vec() * ax).length();
                (rad - c.surf.radius).abs() < 1e-2 && ax > -1e-3 && ax < c.length + 1e-3
            })
        };
        let classify_side =
            |legs: &[usize], elbows: &[usize], pts: &[Point3]| -> SideKind {
                if legs.iter().any(|&lci| on_leg_all(lci, pts)) {
                    return SideKind::OneLeg;
                }
                for &ti in elbows {
                    if loop_on_torus_face(&torus_faces[ti], pts, 1e-2) {
                        return SideKind::OneElbow(ti);
                    }
                }
                SideKind::Multi
            };
        // Project a point onto the nearest ANALYTIC face surface of a filament
        // (leg cylinders / elbow tori).  The corner SSI marches on the SAMPLED
        // swept tubes; the writer emits exact analytic surfaces and SolidWorks
        // knits at ~1 µm — trim curves must lie on the EMITTED geometry (the
        // raw SSI points sit 1–5 µm off it, which SW punishes by dropping the
        // whole elbow face).
        let proj_filament = |p: Point3, legs: &[usize], elbows: &[usize]| -> Point3 {
            let mut best: Option<(f64, Point3)> = None;
            for &lci in legs {
                let c = &cyl_faces[lci];
                let w = p - c.surf.frame.origin;
                let ax = c.surf.axis().dot_vec(w);
                if ax < -0.1 || ax > c.length + 0.1 {
                    continue;
                }
                let rad = w - c.surf.axis().as_vec() * ax;
                let rl = rad.length();
                if rl < 1e-9 {
                    continue;
                }
                let q =
                    c.surf.frame.origin + c.surf.axis().as_vec() * ax + rad * (c.surf.radius / rl);
                let d = (q - p).length();
                if best.map_or(true, |(bd, _)| d < bd) {
                    best = Some((d, q));
                }
            }
            for &ti in elbows {
                let t = &torus_faces[ti].surf;
                let w = p - t.frame.origin;
                let h = t.frame.z.dot_vec(w);
                let q2 = w - t.frame.z.as_vec() * h;
                let a = q2.length();
                if a < 1e-9 {
                    continue;
                }
                let c = t.frame.origin + q2 * (t.major_radius / a);
                let v = p - c;
                let vl = v.length();
                if vl < 1e-9 {
                    continue;
                }
                let q = c + v * (t.minor_radius / vl);
                let d = (q - p).length();
                if best.map_or(true, |(bd, _)| d < bd) {
                    best = Some((d, q));
                }
            }
            best.map_or(p, |(_, q)| q)
        };
        let n_t = filament_tubes.len();
        for gi in 0..n_t {
            if !filament_tubes[gi].2 {
                continue;
            }
            let legs_g = legs_of
                .get(filament_tubes[gi].0.as_str())
                .cloned()
                .unwrap_or_default();
            let elbows_g = elbows_for_legs(&legs_g);
            for fi2 in (gi + 1)..n_t {
                if !filament_tubes[fi2].2 {
                    continue;
                }
                let legs_f2 = legs_of
                    .get(filament_tubes[fi2].0.as_str())
                    .cloned()
                    .unwrap_or_default();
                let elbows_f2 = elbows_for_legs(&legs_f2);
                d_bp += 1;
                // CADCORE_CORNER_TRACER=1 routes corner-loop discovery through
                // the new engine: predictor–corrector tracing on COMPOSITE
                // analytic filaments (cadcore-union).  Loops are exact on the
                // emitted geometry by construction, junction-spanning by
                // design, and found from broad-phase seeds — no sampled-sweep
                // SSI involved.  Default stays on the legacy marcher until
                // the matrix comparison is signed off.
                let corner_loops: Vec<IntersectionPolyline> = if std::env::var(
                    "CADCORE_CORNER_TRACER",
                )
                .map_or(false, |v| v != "0")
                {
                    use cadcore_union::geom::composite::{closed_loops_between, CompositeTube};
                    use cadcore_union::geom::intersect::TraceOptions;
                    let build = |legs: &[usize], elbows: &[usize]| {
                        let mut t = CompositeTube::new();
                        for &lci in legs {
                            t = t.with_leg(cyl_faces[lci].surf, cyl_faces[lci].length);
                        }
                        for &ti2 in elbows {
                            let tf = &torus_faces[ti2];
                            t = t.with_elbow(tf.surf, tf.theta_lo, tf.theta_hi);
                        }
                        t
                    };
                    let ca = build(&legs_g, &elbows_g);
                    let cb = build(&legs_f2, &elbows_f2);
                    let opts = TraceOptions {
                        step: (filament_tubes[gi].1.radius() * 0.2).max(0.02),
                        ..Default::default()
                    };
                    let traced = closed_loops_between(&ca, &cb, &opts);
                    // HYBRID: thick filaments with tight turns self-overlap at
                    // the elbow inner side; composite projection legitimately
                    // refuses there.  Backfill loops the tracer missed from
                    // the legacy marcher — downstream refinement makes them
                    // exact anyway.
                    let mut out: Vec<IntersectionPolyline> = traced
                        .iter()
                        .map(|c| IntersectionPolyline {
                            points: c.points.clone(),
                            closed: true,
                        })
                        .collect();
                    for l in surface_surface_intersection(
                        &filament_tubes[gi].1,
                        &filament_tubes[fi2].1,
                        &cont_ssi_opts,
                    ) {
                        if !l.closed || l.points.len() < 6 {
                            continue;
                        }
                        let probe = l.points[l.points.len() / 2];
                        let known = traced.iter().any(|c| {
                            c.points
                                .iter()
                                .any(|&q| (q - probe).length() < opts.step * 3.0)
                        });
                        if !known {
                            out.push(l);
                        }
                    }
                    out
                } else {
                    surface_surface_intersection(
                        &filament_tubes[gi].1,
                        &filament_tubes[fi2].1,
                        &cont_ssi_opts,
                    )
                };
                for mut loop_pts in corner_loops {
                    if !loop_pts.closed || loop_pts.points.len() < 6 {
                        continue;
                    }
                    // Refine every loop point onto the analytic surfaces of
                    // BOTH filaments (cyclic projection to the curve of
                    // intersection of the emitted geometry).
                    for p in loop_pts.points.iter_mut() {
                        let mut q = *p;
                        for _ in 0..30 {
                            let q0 = q;
                            q = proj_filament(q, &legs_g, &elbows_g);
                            q = proj_filament(q, &legs_f2, &elbows_f2);
                            if (q - q0).length() < 1e-12 {
                                break;
                            }
                        }
                        *p = q;
                    }
                    d_loops += 1;
                    let sg = classify_side(&legs_g, &elbows_g, &loop_pts.points);
                    let sf = classify_side(&legs_f2, &elbows_f2, &loop_pts.points);
                    // Straight-leg sides are Pass 2b/2c territory.
                    if matches!(sg, SideKind::OneLeg) || matches!(sf, SideKind::OneLeg) {
                        d_one_leg += 1;
                        continue;
                    }
                    // Dispatch: which side is the single crossing ELBOW.
                    let (ti, fi, leg_list, elbow_list) = match (&sg, &sf) {
                        (SideKind::OneElbow(a), SideKind::OneElbow(b)) => {
                            // elbow × elbow — simple shared 2-edge hole.
                            let e_surf = swept_from_torus(&torus_faces[*a]);
                            let oriented = orient_loop_as_hole(&e_surf, &loop_pts.points);
                            let (e0, e1) = build_hole_edges(
                                brep, &oriented, &mut vertex_map, &mut edge_map, tolerance,
                                "pass2d_elbow_elbow",
                            );
                            torus_holes[*a].push(HoleRef { e0, e1, reversed: false });
                            torus_holes[*b].push(HoleRef { e0, e1, reversed: true });
                            n_corner_trim += 1;
                            d_pair += 1;
                            continue;
                        }
                        (SideKind::OneElbow(a), SideKind::Multi) => {
                            (*a, fi2, legs_f2.clone(), elbows_f2.clone())
                        }
                        (SideKind::Multi, SideKind::OneElbow(b)) => {
                            (*b, gi, legs_g.clone(), elbows_g.clone())
                        }
                        (SideKind::Multi, SideKind::Multi) => {
                            // TWO-SIDED corner bite: both U's span several
                            // faces.  Split the loop by BOTH filaments' faces,
                            // merge the two cut sets into shared ATOMIC edges,
                            // and give every face its multi-edge run.
                            let (Some(pg), Some(pf)) = (
                                split_multi_span(
                                    &loop_pts.points,
                                    &elbows_g,
                                    &legs_g,
                                    &torus_faces,
                                    &cyl_faces,
                                ),
                                split_multi_span(
                                    &loop_pts.points,
                                    &elbows_f2,
                                    &legs_f2,
                                    &torus_faces,
                                    &cyl_faces,
                                ),
                            ) else {
                                d_split_none += 1;
                                continue;
                            };
                            let same_j = |p: &SpanPiece| {
                                (p.j_in.frame.origin - p.j_out.frame.origin).length() < 1e-3
                            };
                            // Validity guards on BOTH sides.
                            let side_ok = |ps: &[SpanPiece]| -> bool {
                                let mut ok = ps.iter().all(|p| {
                                    (p.pts[0] - *p.pts.last().unwrap()).length() > 0.01
                                        || !same_j(p)
                                });
                                // Through-pieces (j_in != j_out) must come in
                                // pairs per face - legs AND elbows alike.
                                // BTreeMap: deterministic iteration (see cap_notches note).
                                let mut tcl: std::collections::BTreeMap<usize, Vec<usize>> =
                                    std::collections::BTreeMap::new();
                                let mut tce: std::collections::BTreeMap<usize, Vec<usize>> =
                                    std::collections::BTreeMap::new();
                                for (i, p) in ps.iter().enumerate() {
                                    if same_j(p) {
                                        continue;
                                    }
                                    match p.owner {
                                        PieceOwner::Leg(lci) => {
                                            tcl.entry(lci).or_default().push(i);
                                        }
                                        PieceOwner::Elbow(ti2) => {
                                            tce.entry(ti2).or_default().push(i);
                                        }
                                    }
                                }
                                for (&lci, list) in &tcl {
                                    ok &= list.len() == 2
                                        && holes[lci].is_empty()
                                        && open_notches[lci].is_empty()
                                        && matches!(cyl_faces[lci].start, FaceBoundary::Circle(_))
                                        && matches!(cyl_faces[lci].end, FaceBoundary::Circle(_));
                                }
                                for (&ti2, list) in &tce {
                                    ok &= list.len() == 2
                                        && torus_holes[ti2].is_empty()
                                        && torus_open_notches[ti2].is_empty();
                                    // an elbow must not mix notches (same loop)
                                    ok &= !ps.iter().any(|p| {
                                        matches!(p.owner, PieceOwner::Elbow(t) if t == ti2)
                                            && same_j(p)
                                    });
                                }
                                ok
                            };
                            if !side_ok(&pg) || !side_ok(&pf) {
                                d_guard += 1;
                                if std::env::var("CADCORE_DUMP_CONT").is_ok() {
                                    eprintln!(
                                        "[corner:guard-side] g_ok={} f_ok={} g={:?} f={:?}",
                                        side_ok(&pg),
                                        side_ok(&pf),
                                        pg.iter().map(|p| (p.owner, same_j(p), ((p.pts[0] - *p.pts.last().unwrap()).length()*1000.0).round()/1000.0)).collect::<Vec<_>>(),
                                        pf.iter().map(|p| (p.owner, same_j(p), ((p.pts[0] - *p.pts.last().unwrap()).length()*1000.0).round()/1000.0)).collect::<Vec<_>>()
                                    );
                                }
                                continue;
                            }
                            // Merged cut set: precise circular position of each
                            // side's joints along the original loop.
                            let lp = &loop_pts.points;
                            let nl = lp.len();
                            let cut_pos = |i1: usize, s: Point3| -> f64 {
                                let p0 = lp[i1];
                                let p1 = lp[(i1 + 1) % nl];
                                let d = p1 - p0;
                                let l2 = d.dot(d);
                                let t = if l2 < 1e-18 {
                                    0.5
                                } else {
                                    ((s - p0).dot(d) / l2).clamp(0.0, 1.0)
                                };
                                i1 as f64 + t
                            };
                            // Cut points are interpolated junction crossings —
                            // refine each onto its junction circle ∩ the OTHER
                            // filament's analytic surface, so the atomic edge
                            // endpoints lie exactly on all incident faces.
                            let refine_cut = |q0: Point3,
                                              jc: &Circle3,
                                              ol: &[usize],
                                              oe: &[usize]|
                             -> Point3 {
                                let mut q = q0;
                                for _ in 0..30 {
                                    let prev = q;
                                    q = jc.point_at(circle_angle(jc, q));
                                    q = proj_filament(q, ol, oe);
                                    if (q - prev).length() < 1e-12 {
                                        break;
                                    }
                                }
                                // final word: exactly on the junction circle
                                jc.point_at(circle_angle(jc, q))
                            };
                            let mut cuts: Vec<(f64, Point3)> = Vec::new();
                            let mut g_cut_pos: Vec<f64> = Vec::new();
                            let mut f_cut_pos: Vec<f64> = Vec::new();
                            let mut g_cut_pts: Vec<Point3> = Vec::new();
                            let mut f_cut_pts: Vec<Point3> = Vec::new();
                            for p in &pg {
                                let s =
                                    refine_cut(*p.pts.last().unwrap(), &p.j_out, &legs_f2, &elbows_f2);
                                let pos = cut_pos(p.i1, s);
                                g_cut_pos.push(pos);
                                g_cut_pts.push(s);
                                cuts.push((pos, s));
                            }
                            for p in &pf {
                                let s =
                                    refine_cut(*p.pts.last().unwrap(), &p.j_out, &legs_g, &elbows_g);
                                let pos = cut_pos(p.i1, s);
                                f_cut_pos.push(pos);
                                f_cut_pts.push(s);
                                cuts.push((pos, s));
                            }
                            cuts.sort_by(|a, b| a.0.total_cmp(&b.0));
                            cuts.dedup_by(|a, b| (a.1 - b.1).length() < 1e-6);
                            let nc = cuts.len();
                            if nc < 2 {
                                d_guard += 1;
                                continue;
                            }
                            // Atomic edges between consecutive cuts.
                            let mut atomic: Vec<EdgeId> = Vec::with_capacity(nc);
                            let mut atomic_ok = true;
                            for c in 0..nc {
                                let (p0, s0) = cuts[c];
                                let (p1m, s1) = cuts[(c + 1) % nc];
                                let p1 = if c + 1 == nc { p1m + nl as f64 } else { p1m };
                                let mut pts: Vec<Point3> = vec![s0];
                                let mut i = p0.floor() as usize + 1;
                                while (i as f64) < p1 {
                                    pts.push(lp[i % nl]);
                                    i += 1;
                                }
                                pts.push(s1);
                                pts.dedup_by(|a, b| (*a - *b).length() < 1e-9);
                                if pts.len() < 2 {
                                    atomic_ok = false;
                                    break;
                                }
                                atomic.push(build_single_edge(
                                    brep, &pts, &mut vertex_map, &mut edge_map, tolerance,
                                ));
                            }
                            if !atomic_ok {
                                d_guard += 1;
                                continue;
                            }
                            // Run of atomic edges covering (cut_a .. cut_b).
                            let merged_idx = |pos: f64| -> usize {
                                cuts.iter()
                                    .enumerate()
                                    .min_by(|(_, a), (_, b)| {
                                        (a.0 - pos).abs().total_cmp(&(b.0 - pos).abs())
                                    })
                                    .map(|(i, _)| i)
                                    .unwrap_or(0)
                            };
                            let run_of = |from_pos: f64, to_pos: f64| -> Vec<EdgeId> {
                                let a = merged_idx(from_pos);
                                let b = merged_idx(to_pos);
                                let mut v = Vec::new();
                                let mut c = a;
                                loop {
                                    v.push(atomic[c]);
                                    c = (c + 1) % nc;
                                    if c == b {
                                        break;
                                    }
                                    if v.len() > nc {
                                        break;
                                    }
                                }
                                v
                            };
                            // Side traversal: material (outside the OTHER
                            // tube) on the LEFT; owner-face normal per point.
                            let owner_normal = |ps: &[SpanPiece], q: Point3| -> Option<cadcore_math::Vec3> {
                                let mut best: Option<(f64, PieceOwner)> = None;
                                for p in ps {
                                    for &pp in &p.pts {
                                        let d = (pp - q).length();
                                        if best.map_or(true, |(bd, _)| d < bd) {
                                            best = Some((d, p.owner));
                                        }
                                    }
                                }
                                let (_, o) = best?;
                                let geom = match o {
                                    PieceOwner::Elbow(ti2) => FaceGeom::Torus(torus_faces[ti2].surf),
                                    PieceOwner::Leg(ci2) => FaceGeom::Cylinder(cyl_faces[ci2].surf),
                                };
                                Some(face_outward_normal(&geom, q))
                            };
                            let side_same = |ps: &[SpanPiece], other: &SweptTubeSurface| -> bool {
                                let mut okc = 0usize;
                                let mut n = 0usize;
                                let step = (lp.len() / 9).max(1);
                                for k in (0..lp.len().saturating_sub(1)).step_by(step) {
                                    let p = lp[k];
                                    let t = lp[k + 1] - p;
                                    if t.length() < 1e-9 {
                                        continue;
                                    }
                                    let Some(nrm) = owner_normal(ps, p) else { continue };
                                    let left = nrm.cross(t).normalize();
                                    let probe = p + left * 0.05;
                                    n += 1;
                                    if other.signed_distance(probe).map_or(true, |d| d >= 0.0) {
                                        okc += 1;
                                    }
                                }
                                n == 0 || okc * 2 >= n
                            };
                            let g_same = side_same(&pg, &filament_tubes[fi2].1);
                            let f_same = side_same(&pf, &filament_tubes[gi].1);
                            if g_same == f_same {
                                d_guard += 1; // sides must traverse oppositely
                                if std::env::var("CADCORE_DUMP_CONT").is_ok() {
                                    eprintln!("[corner:guard-dir] g_same={g_same} f_same={f_same}");
                                }
                                continue;
                            }
                            // Distribute each side's pieces as multi-edge runs.
                            for (ps, cut_pos_v, cut_pts_v, x_same, other_idx) in [
                                (&pg, &g_cut_pos, &g_cut_pts, g_same, fi2),
                                (&pf, &f_cut_pos, &f_cut_pts, f_same, gi),
                            ] {
                                let other = &filament_tubes[other_idx].1;
                                let m = ps.len();
                                for (k, p) in ps.iter().enumerate() {
                                    let from_pos = cut_pos_v[(k + m - 1) % m];
                                    let to_pos = cut_pos_v[k];
                                    // refined endpoints: entry = previous cut,
                                    // exit = own cut
                                    let first_ref = cut_pts_v[(k + m - 1) % m];
                                    let last_ref = cut_pts_v[k];
                                    let run = run_of(from_pos, to_pos);
                                    if run.is_empty() {
                                        continue;
                                    }
                                    let bite_ccw = {
                                        let j = &p.j_in;
                                        let a1 = circle_angle(j, first_ref);
                                        let a2 = circle_angle(j, last_ref);
                                        let mid = a1 + (a2 - a1).rem_euclid(2.0 * PI) * 0.5;
                                        other
                                            .signed_distance(j.point_at(mid))
                                            .map_or(false, |d| d < 0.0)
                                    };
                                    match p.owner {
                                        PieceOwner::Elbow(ti2) if same_j(p) => {
                                            let tf2 = &torus_faces[ti2];
                                            let side = if (p.j_in.frame.origin
                                                - tf2.start_circle.frame.origin)
                                                .length()
                                                < 1e-3
                                            {
                                                BoundarySide::Start
                                            } else {
                                                BoundarySide::End
                                            };
                                            torus_open_notches[ti2].push(OpenNotchRef {
                                                edges: run,
                                                side,
                                                reversed: !x_same,
                                                gap_circle: Some(p.j_in),
                                                bite_ccw_from_start: bite_ccw,
                                            });
                                        }
                                        PieceOwner::Elbow(ti2) => {
                                            // THROUGH-elbow piece: this elbow
                                            // becomes strip bands like a
                                            // connector.
                                            let tf2 = &torus_faces[ti2];
                                            let start_origin =
                                                tf2.start_circle.frame.origin;
                                            let first = first_ref;
                                            let last = last_ref;
                                            let first_on_start = (p.j_in.frame.origin
                                                - start_origin)
                                                .length()
                                                < (p.j_out.frame.origin - start_origin)
                                                    .length();
                                            let (p_j1, p_j2) = if first_on_start {
                                                (first, last)
                                            } else {
                                                (last, first)
                                            };
                                            torus_through_pieces[ti2].push(ThroughPiece {
                                                edges: run,
                                                bite: bite_counter,
                                                crossing: other_idx,
                                                p_j1,
                                                p_j2,
                                                strip_same: x_same,
                                            });
                                        }
                                        PieceOwner::Leg(lci) if same_j(p) => {
                                            let cf2 = &cyl_faces[lci];
                                            let side =
                                                if boundary_circle_matches(&cf2.start, &p.j_in) {
                                                    BoundarySide::Start
                                                } else {
                                                    BoundarySide::End
                                                };
                                            open_notches[lci].push(OpenNotchRef {
                                                edges: run,
                                                side,
                                                reversed: !x_same,
                                                gap_circle: Some(p.j_in),
                                                bite_ccw_from_start: bite_ccw,
                                            });
                                        }
                                        PieceOwner::Leg(lci) => {
                                            let start_origin = match &cyl_faces[lci].start {
                                                FaceBoundary::Circle(c) => c.frame.origin,
                                                FaceBoundary::Ellipse(e) => e.frame.origin,
                                            };
                                            let first = first_ref;
                                            let last = last_ref;
                                            let first_on_start = (p.j_in.frame.origin
                                                - start_origin)
                                                .length()
                                                < (p.j_out.frame.origin - start_origin).length();
                                            let (p_j1, p_j2) = if first_on_start {
                                                (first, last)
                                            } else {
                                                (last, first)
                                            };
                                            through_pieces[lci].push(ThroughPiece {
                                                edges: run,
                                                bite: bite_counter,
                                                crossing: other_idx,
                                                p_j1,
                                                p_j2,
                                                strip_same: x_same,
                                            });
                                        }
                                    }
                                }
                            }
                            bite_counter += 1;
                            n_corner_trim += 1;
                            continue;
                        }
                        _ => continue,
                    };
                    let tf = &torus_faces[ti];
                    let tube = &filament_tubes[fi].1;
                    let legs = leg_list;
                    let elbows_f = elbow_list;
                    // Bitten side spans several foreign faces → multi-span.
                    let Some(pieces) = split_multi_span(
                        &loop_pts.points,
                        &elbows_f,
                        &legs,
                        &torus_faces,
                        &cyl_faces,
                    ) else {
                        d_split_none += 1;
                        continue;
                    };
                    let same_j = |p: &SpanPiece| {
                        (p.j_in.frame.origin - p.j_out.frame.origin).length() < 1e-3
                    };
                    let mut ok = pieces.iter().all(|p| {
                        (p.pts[0] - *p.pts.last().unwrap()).length() > 0.01 || !same_j(p)
                    });
                    ok &= pieces.iter().all(|p| match p.owner {
                        PieceOwner::Elbow(_) => same_j(p),
                        PieceOwner::Leg(_) => true,
                    });
                    // BTreeMap: deterministic iteration (see cap_notches note).
                    let mut tc_legs: std::collections::BTreeMap<usize, Vec<usize>> =
                        std::collections::BTreeMap::new();
                    for (i, p) in pieces.iter().enumerate() {
                        if let PieceOwner::Leg(lci) = p.owner {
                            if !same_j(p) {
                                tc_legs.entry(lci).or_default().push(i);
                            }
                        }
                    }
                    for (&lci, list) in &tc_legs {
                        ok &= list.len() == 2
                            && holes[lci].is_empty()
                            && open_notches[lci].is_empty();
                        if ok {
                            let a = &pieces[list[0]];
                            let b = &pieces[list[1]];
                            ok &= (a.j_in.frame.origin - b.j_out.frame.origin).length() < 1e-3
                                && (a.j_out.frame.origin - b.j_in.frame.origin).length() < 1e-3
                                && matches!(cyl_faces[lci].start, FaceBoundary::Circle(_))
                                && matches!(cyl_faces[lci].end, FaceBoundary::Circle(_));
                        }
                    }
                    if !ok {
                        d_guard += 1;
                        continue;
                    }
                    let mut piece_edges: Vec<EdgeId> = Vec::with_capacity(pieces.len());
                    for p in &pieces {
                        piece_edges.push(build_single_edge(
                            brep, &p.pts, &mut vertex_map, &mut edge_map, tolerance,
                        ));
                    }
                    // Material-on-the-left in the ELBOW's surface decides the
                    // crossing-side hole direction.
                    let mut chain: Vec<Point3> = Vec::new();
                    for p in &pieces {
                        chain.extend(p.pts.iter().copied());
                    }
                    let geom_e = FaceGeom::Torus(tf.surf);
                    let material_left_ok = {
                        let mut okc = 0usize;
                        let mut n = 0usize;
                        let step = (chain.len() / 9).max(1);
                        for k in (0..chain.len().saturating_sub(1)).step_by(step) {
                            let p = chain[k];
                            let t = chain[k + 1] - p;
                            if t.length() < 1e-9 {
                                continue;
                            }
                            let nrm = face_outward_normal(&geom_e, p);
                            let left = nrm.cross(t).normalize();
                            let probe = p + left * 0.05;
                            n += 1;
                            if tube.signed_distance(probe).map_or(true, |d| d >= 0.0) {
                                okc += 1;
                            }
                        }
                        n == 0 || okc * 2 >= n
                    };
                    let hole_rev = !material_left_ok;
                    torus_multi_holes[ti].push(MultiHoleRef {
                        edges: piece_edges.clone(),
                        reversed: hole_rev,
                    });
                    // Bitten interval flag sampled against the CROSSING
                    // elbow's filament tube.
                    let cross_tube = elbow_tube[ti].map(|g| &filament_tubes[g].1);
                    let bite_ccw_of = |p: &SpanPiece| -> bool {
                        let j = &p.j_in;
                        let a1 = circle_angle(j, p.pts[0]);
                        let a2 = circle_angle(j, *p.pts.last().unwrap());
                        let mid = a1 + (a2 - a1).rem_euclid(2.0 * PI) * 0.5;
                        let q = j.point_at(mid);
                        cross_tube
                            .and_then(|t| t.signed_distance(q))
                            .map_or(false, |d| d < 0.0)
                    };
                    for (pi, p) in pieces.iter().enumerate() {
                        match p.owner {
                            PieceOwner::Elbow(ti2) => {
                                let tf2 = &torus_faces[ti2];
                                let side = if (p.j_in.frame.origin
                                    - tf2.start_circle.frame.origin)
                                    .length()
                                    < 1e-3
                                {
                                    BoundarySide::Start
                                } else {
                                    BoundarySide::End
                                };
                                torus_open_notches[ti2].push(OpenNotchRef {
                                    edges: vec![piece_edges[pi]],
                                    side,
                                    reversed: !hole_rev,
                                    gap_circle: Some(p.j_in),
                                    bite_ccw_from_start: bite_ccw_of(p),
                                });
                            }
                            PieceOwner::Leg(lci) if same_j(p) => {
                                let cf2 = &cyl_faces[lci];
                                let side = if boundary_circle_matches(&cf2.start, &p.j_in) {
                                    BoundarySide::Start
                                } else {
                                    BoundarySide::End
                                };
                                open_notches[lci].push(OpenNotchRef {
                                    edges: vec![piece_edges[pi]],
                                    side,
                                    reversed: !hole_rev,
                                    gap_circle: Some(p.j_in),
                                    bite_ccw_from_start: bite_ccw_of(p),
                                });
                            }
                            PieceOwner::Leg(_) => {}
                        }
                    }
                    for (&lci, list) in &tc_legs {
                        let start_origin = match &cyl_faces[lci].start {
                            FaceBoundary::Circle(c) => c.frame.origin,
                            FaceBoundary::Ellipse(e) => e.frame.origin,
                        };
                        for &pi in list {
                            let p = &pieces[pi];
                            let first = p.pts[0];
                            let last = *p.pts.last().unwrap();
                            let first_on_start = (p.j_in.frame.origin - start_origin).length()
                                < (p.j_out.frame.origin - start_origin).length();
                            let (p_j1, p_j2) = if first_on_start {
                                (first, last)
                            } else {
                                (last, first)
                            };
                            through_pieces[lci].push(ThroughPiece {
                                edges: vec![piece_edges[pi]],
                                bite: bite_counter,
                                crossing: elbow_tube[ti].unwrap_or(usize::MAX),
                                p_j1,
                                p_j2,
                                strip_same: hole_rev,
                            });
                        }
                    }
                    bite_counter += 1;
                    n_corner_trim += 1;
                }
            }
        }
    }
    if std::env::var("CADCORE_DUMP_CONT").is_ok() {
        eprintln!("[union][corner] trims={n_corner_trim} bp={d_bp} loops={d_loops} not_on_elbow={d_not_on_elbow} one_leg={d_one_leg} pair={d_pair} split_none={d_split_none} guard={d_guard}");
    }
    let _ = (n_corner_trim, d_bp, d_loops, d_not_on_elbow, d_one_leg, d_pair, d_split_none, d_guard);

    eprintln!(
        "[union] cyl={} tor={} cross_trim={} | elbow: bp={} overlap={} overlap_found={} curves={} open={} open_trim={} fallback={} trim={} min_gap={:.4}",
        cyl_faces.len(),
        torus_faces.len(),
        n_trimmed,
        n_bp,
        n_overlap,
        n_overlap_found,
        n_curves,
        n_open_boundary,
        n_open_converted,
        n_marching_fallback,
        n_torus_trim,
        min_real_gap
    );
    if n_open_boundary > n_open_converted {
        eprintln!(
            "[union][open:summary] converted={} rejected={} not_same_boundary_or_snap={} other_reasons={}",
            n_open_converted,
            n_open_boundary.saturating_sub(n_open_converted),
            n_open_reject_not_same_boundary,
            n_open_boundary.saturating_sub(n_open_converted).saturating_sub(n_open_reject_not_same_boundary),
        );
    }
    let _ = (n_should_cross, n_dropped, n_axial_drop);

    // ── Pass 2e: BUTT-END joints (filament end cap inside a foreign tube) ─────
    // A capped free end overlaps the perpendicular legs of the layers above /
    // below.  Both surfaces are ANALYTIC here (cylinder×cylinder + plane), so
    // the trim is built in closed form — the numeric SSI tracer runs away on
    // these near-tangent end regions (runaway 4001-point chains).
    //   * LATERAL: the cyl∩cyl bite loop clipped at the cap plane → an
    //     end-boundary notch (e_lat);
    //   * CAP: trimmed plane — disk minus the foreign cylinder's conic
    //     (e_cap, between the same two snap points);
    //   * FOREIGN leg: 2-edge hole [e_lat, e_cap].
    let cap_trim_enabled = std::env::var("CADCORE_CAP_TRIM").map_or(true, |v| v != "0");
    let mut n_cap_trim = 0usize;
    let mut c_bp = 0usize;
    let mut c_loops = 0usize;
    let mut c_cross = 0usize;
    let mut c_norun = 0usize;
    let mut c_sliv = 0usize;
    let mut c_lines = 0usize;
    for (&(ci, sd), _cap_fid) in cap_map.clone().iter().filter(|_| cap_trim_enabled) {
        let cf = &cyl_faces[ci];
        let bnd = match sd {
            BoundarySide::Start => &cf.start,
            BoundarySide::End => &cf.end,
        };
        let FaceBoundary::Circle(jc) = bnd else { continue };
        let jc = *jc;
        let cap_fid = cap_map[&(ci, sd)];
        let c_ax = axial(&cf.surf, jc.frame.origin);
        // Inward = direction along the leg axis from the cap INTO material.
        let inward: f64 = if c_ax < cf.length * 0.5 { 1.0 } else { -1.0 };
        for fci in 0..cyl_faces.len() {
            if fci == ci
                || cyl_faces[fci].solid_name == cf.solid_name
            {
                continue;
            }
            let fs = &cyl_faces[fci].surf;
            // Broad-phase: foreign axis near the cap centre.
            let w = jc.frame.origin - fs.frame.origin;
            let f_ax = fs.axis().dot_vec(w);
            let radial = (w - fs.axis().as_vec() * f_ax).length();
            if radial > cf.surf.radius + fs.radius
                || f_ax < -cf.surf.radius
                || f_ax > cyl_faces[fci].length + cf.surf.radius
            {
                continue;
            }
            c_bp += 1;
            for lp in cyl_cyl_intersection(&cf.surf, fs, SAMPLES) {
                if !lp.closed || lp.points.len() < 8 {
                    continue;
                }
                c_loops += 1;
                // Signed side of the cap plane (positive = inside material).
                let sideof = |p: Point3| (axial(&cf.surf, p) - c_ax) * inward;
                let n_in = lp.points.iter().filter(|&&p| sideof(p) > 0.0).count();
                if n_in == 0 || n_in == lp.points.len() {
                    continue; // fully outside / fully interior (Pass 2's case)
                }
                c_cross += 1;
                // Extract the INSIDE run with exact plane crossings.
                let n = lp.points.len();
                let start = match (0..n)
                    .find(|&i| sideof(lp.points[i]) <= 0.0 && sideof(lp.points[(i + 1) % n]) > 0.0)
                {
                    Some(s) => s,
                    None => continue,
                };
                let lerp_cross = |a: Point3, b: Point3| -> Point3 {
                    let sa = sideof(a);
                    let sb = sideof(b);
                    let t = (sa / (sa - sb)).clamp(0.0, 1.0);
                    a + (b - a) * t
                };
                let mut pts: Vec<Point3> = Vec::new();
                pts.push(lerp_cross(lp.points[start], lp.points[(start + 1) % n]));
                let mut i = (start + 1) % n;
                let mut crossed_back = false;
                let mut guard = 0;
                while sideof(lp.points[i]) > 0.0 {
                    pts.push(lp.points[i]);
                    let nx = (i + 1) % n;
                    if sideof(lp.points[nx]) <= 0.0 {
                        pts.push(lerp_cross(lp.points[i], lp.points[nx]));
                        crossed_back = true;
                        break;
                    }
                    i = nx;
                    guard += 1;
                    if guard > n {
                        break;
                    }
                }
                if !crossed_back || pts.len() < 4 {
                    c_norun += 1;
                    continue;
                }
                // Snap the crossings onto the end circle.
                let mut p1 = snap_to_circle(&jc, pts[0]);
                let mut p2 = snap_to_circle(&jc, *pts.last().unwrap());
                if (p1 - p2).length() < 0.02 {
                    c_sliv += 1;
                    continue; // tangent sliver
                }
                pts[0] = p1;
                *pts.last_mut().unwrap() = p2;
                // Conic on the cap plane (foreign cylinder ∩ plane), the
                // branch inside the disk, p2 → p1.
                let plane = Plane3::from_origin_normal(jc.frame.origin, jc.frame.z);
                let arc_pts: Vec<Point3> = match cyl_plane_intersection(fs, &plane) {
                    CylPlaneCurve::Ellipse(e) => {
                        let ang = |p: Point3| -> f64 {
                            let d = p - e.frame.origin;
                            (e.frame.y.dot_vec(d) / e.semi_minor)
                                .atan2(e.frame.x.dot_vec(d) / e.semi_major)
                        };
                        let t2 = ang(p2);
                        let t1 = ang(p1);
                        let sample = |dirn: f64| -> Vec<Point3> {
                            let sweep = if dirn > 0.0 {
                                (t1 - t2).rem_euclid(2.0 * PI)
                            } else {
                                -((t2 - t1).rem_euclid(2.0 * PI))
                            };
                            let nseg = 24usize;
                            (0..=nseg)
                                .map(|k| e.point_at(t2 + sweep * k as f64 / nseg as f64))
                                .collect()
                        };
                        let inside = |v: &Vec<Point3>| -> bool {
                            let mid = v[v.len() / 2];
                            (mid - jc.frame.origin).length() < jc.radius - 1e-4
                        };
                        let s1 = sample(1.0);
                        if inside(&s1) {
                            s1
                        } else {
                            let s2 = sample(-1.0);
                            if inside(&s2) {
                                s2
                            } else {
                                continue;
                            }
                        }
                    }
                    CylPlaneCurve::Circle(c) => {
                        let t2 = circle_angle(&c, p2);
                        let t1 = circle_angle(&c, p1);
                        let sweep = (t1 - t2).rem_euclid(2.0 * PI);
                        let nseg = 24usize;
                        let v: Vec<Point3> = (0..=nseg)
                            .map(|k| c.point_at(t2 + sweep * k as f64 / nseg as f64))
                            .collect();
                        if (v[v.len() / 2] - jc.frame.origin).length() < jc.radius - 1e-4 {
                            v
                        } else {
                            let sweep2 = sweep - 2.0 * PI;
                            (0..=nseg)
                                .map(|k| c.point_at(t2 + sweep2 * k as f64 / nseg as f64))
                                .collect()
                        }
                    }
                    CylPlaneCurve::Lines(ls) => {
                        // Perpendicular foreign leg: the cap plane is parallel
                        // to its axis → the conic degenerates to ruling LINES
                        // and the cap boundary is a straight CHORD.  Pick the
                        // nearest line and recompute the endpoints EXACTLY as
                        // line ∩ end circle (the snapped crossings carry
                        // sampling noise).
                        let line_d = |l: &cadcore_geom::Line3, p: Point3| {
                            let d = p - l.origin;
                            (d - l.direction.as_vec() * l.direction.dot_vec(d)).length()
                        };
                        let Some(l) = ls
                            .iter()
                            .min_by(|a, b| {
                                (line_d(a, p1) + line_d(a, p2))
                                    .total_cmp(&(line_d(b, p1) + line_d(b, p2)))
                            })
                        else {
                            c_lines += 1;
                            continue;
                        };
                        if line_d(l, p1).max(line_d(l, p2)) > 0.08 {
                            c_lines += 1;
                            continue; // crossings on different lines (deep strip)
                        }
                        // line ∩ circle: |o + t·d − c|² = R² (in the cap plane)
                        let oc = l.origin - jc.frame.origin;
                        let b_half = l.direction.dot_vec(oc);
                        let c0 = oc.dot(oc) - jc.radius * jc.radius;
                        let disc = b_half * b_half - c0;
                        if disc <= 0.0 {
                            c_lines += 1;
                            continue;
                        }
                        let sq = disc.sqrt();
                        let q1 = l.origin + l.direction * (-b_half - sq);
                        let q2 = l.origin + l.direction * (-b_half + sq);
                        let (qa, qb) = if (q1 - p1).length() < (q2 - p1).length() {
                            (q1, q2)
                        } else {
                            (q2, q1)
                        };
                        p1 = qa;
                        p2 = qb;
                        pts[0] = p1;
                        *pts.last_mut().unwrap() = p2;
                        let nseg = 8usize;
                        (0..=nseg)
                            .map(|k| p2 + (p1 - p2) * (k as f64 / nseg as f64))
                            .collect()
                    }
                    _ => continue,
                };
                let mut cap_pts = arc_pts;
                cap_pts[0] = p2;
                *cap_pts.last_mut().unwrap() = p1;
                // OVERLAP GUARD: the end-notch must not cross any existing
                // trim on this lateral (a perpendicular bite near the end) —
                // overlapping wires are invalid B-Rep and SolidWorks drops
                // the whole face (OCC silently heals).  Compare in (θ,v).
                let uv_of = |p: Point3| -> (f64, f64) {
                    let w2 = p - cf.surf.frame.origin;
                    let ax2 = cf.surf.axis().dot_vec(w2);
                    let th = cf
                        .surf
                        .frame
                        .y
                        .dot_vec(w2)
                        .atan2(cf.surf.frame.x.dot_vec(w2));
                    (th, ax2)
                };
                let chain_uv = |pts3: &[Point3]| -> Vec<(f64, f64)> {
                    let mut out: Vec<(f64, f64)> = Vec::with_capacity(pts3.len());
                    for &p in pts3 {
                        let (mut u, v) = uv_of(p);
                        if let Some(&(pu, _)) = out.last() {
                            while u - pu > PI {
                                u -= 2.0 * PI;
                            }
                            while u - pu < -PI {
                                u += 2.0 * PI;
                            }
                        }
                        out.push((u, v));
                    }
                    let mu = out.iter().map(|x| x.0).sum::<f64>() / out.len() as f64;
                    let k = (mu / (2.0 * PI)).round();
                    out.iter().map(|&(u, v)| (u - k * 2.0 * PI, v)).collect()
                };
                let seg_int = |a1: (f64, f64),
                               a2: (f64, f64),
                               b1: (f64, f64),
                               b2: (f64, f64)|
                 -> bool {
                    let d1 = (a2.0 - a1.0, a2.1 - a1.1);
                    let d2 = (b2.0 - b1.0, b2.1 - b1.1);
                    let den = d1.0 * d2.1 - d1.1 * d2.0;
                    if den.abs() < 1e-15 {
                        return false;
                    }
                    let t = ((b1.0 - a1.0) * d2.1 - (b1.1 - a1.1) * d2.0) / den;
                    let u = ((b1.0 - a1.0) * d1.1 - (b1.1 - a1.1) * d1.0) / den;
                    t > 1e-9 && t < 1.0 - 1e-9 && u > 1e-9 && u < 1.0 - 1e-9
                };
                let chains_cross = |a: &[(f64, f64)], b: &[(f64, f64)]| -> bool {
                    for da in [-2.0 * PI, 0.0, 2.0 * PI] {
                        for ia in 0..a.len().saturating_sub(1) {
                            for ib in 0..b.len().saturating_sub(1) {
                                if seg_int(
                                    (a[ia].0 + da, a[ia].1),
                                    (a[ia + 1].0 + da, a[ia + 1].1),
                                    b[ib],
                                    b[ib + 1],
                                ) {
                                    return true;
                                }
                            }
                        }
                    }
                    false
                };
                let notch_uv = chain_uv(&pts);
                let edge_pts = |brep: &BRep, e: EdgeId| -> Vec<Point3> {
                    brep.edges
                        .get(e)
                        .map(|edge| match &edge.geom {
                            EdgeGeom::Polyline(p) => p.clone(),
                            _ => Vec::new(),
                        })
                        .unwrap_or_default()
                };
                let mut conflict = false;
                {
                    let mut check_edges: Vec<EdgeId> = Vec::new();
                    for h in &holes[ci] {
                        check_edges.push(h.e0);
                        check_edges.push(h.e1);
                    }
                    for mh in &multi_holes[ci] {
                        check_edges.extend(mh.edges.iter().copied());
                    }
                    for nt in &open_notches[ci] {
                        check_edges.extend(nt.edges.iter().copied());
                    }
                    for e in check_edges {
                        let p3 = edge_pts(brep, e);
                        if p3.len() >= 2 && chains_cross(&notch_uv, &chain_uv(&p3)) {
                            conflict = true;
                            break;
                        }
                    }
                }
                // MUTUAL-BITE guard: when the foreign leg is ALSO capped
                // and ends at the same corner, the hole would run past ITS
                // end — that needs the full multi-face split; skip for now.
                {
                    let lenf = cyl_faces[fci].length;
                    let mut bad = false;
                    for &p in pts.iter().chain(cap_pts.iter()) {
                        let ax2 = fs.axis().dot_vec(p - fs.frame.origin);
                        for (sdx, lim) in [(BoundarySide::Start, 0.0), (BoundarySide::End, lenf)]
                        {
                            if cap_map.contains_key(&(fci, sdx)) && (ax2 - lim).abs() < 0.08 {
                                bad = true;
                            }
                        }
                        if ax2 < -0.02 || ax2 > lenf + 0.02 {
                            bad = true;
                        }
                    }
                    if bad {
                        c_lines += 1;
                        continue;
                    }
                }
                // SYMMETRIC check on the FOREIGN leg: the new hole
                // [e_lat, e_cap] must not cross fci's existing trims either.
                if !conflict {
                    let uv_f = |p: Point3| -> (f64, f64) {
                        let w2 = p - fs.frame.origin;
                        let ax2 = fs.axis().dot_vec(w2);
                        let th = fs.frame.y.dot_vec(w2).atan2(fs.frame.x.dot_vec(w2));
                        (th, ax2)
                    };
                    let chain_uv_f = |pts3: &[Point3]| -> Vec<(f64, f64)> {
                        let mut out: Vec<(f64, f64)> = Vec::with_capacity(pts3.len());
                        for &p in pts3 {
                            let (mut u, v) = uv_f(p);
                            if let Some(&(pu, _)) = out.last() {
                                while u - pu > PI {
                                    u -= 2.0 * PI;
                                }
                                while u - pu < -PI {
                                    u += 2.0 * PI;
                                }
                            }
                            out.push((u, v));
                        }
                        let mu = out.iter().map(|x| x.0).sum::<f64>() / out.len() as f64;
                        let k = (mu / (2.0 * PI)).round();
                        out.iter().map(|&(u, v)| (u - k * 2.0 * PI, v)).collect()
                    };
                    let mut hole_chain = pts.clone();
                    hole_chain.extend(cap_pts.iter().copied());
                    let hole_uv = chain_uv_f(&hole_chain);
                    let mut check_edges: Vec<EdgeId> = Vec::new();
                    for h in &holes[fci] {
                        check_edges.push(h.e0);
                        check_edges.push(h.e1);
                    }
                    for mh in &multi_holes[fci] {
                        check_edges.extend(mh.edges.iter().copied());
                    }
                    for nt in &open_notches[fci] {
                        check_edges.extend(nt.edges.iter().copied());
                    }
                    for e in check_edges {
                        let p3 = edge_pts(brep, e);
                        if p3.len() >= 2 && chains_cross(&hole_uv, &chain_uv_f(&p3)) {
                            conflict = true;
                            break;
                        }
                    }
                }
                if conflict {
                    c_lines += 1; // counted as skipped (overlap)
                    continue;
                }
                // Shared edges: lateral p1→p2, cap conic p2→p1 (chained).
                let e_lat =
                    build_single_edge(brep, &pts, &mut vertex_map, &mut edge_map, tolerance);
                let e_cap =
                    build_single_edge(brep, &cap_pts, &mut vertex_map, &mut edge_map, tolerance);
                // Hole orientation on the FOREIGN leg: the compact bite
                // loop's signed area in F's (θ,v) — the material-left probe
                // is unreliable here (the chord lies ON the cap plane, so
                // probes exit the capped solid regardless of direction).
                let hole_rev = {
                    let mut chain: Vec<Point3> = pts.clone();
                    chain.extend(cap_pts.iter().copied());
                    signed_area_theta_v(fs, &chain) > 0.0
                };
                multi_holes[fci].push(MultiHoleRef {
                    edges: vec![e_lat, e_cap],
                    reversed: hole_rev,
                });
                // Bitten interval on the end circle, sampled against the
                // foreign CYLINDER (radial test — robust and local).
                let bite_of = |a1: f64, a2: f64| -> bool {
                    let mid = a1 + (a2 - a1).rem_euclid(2.0 * PI) * 0.5;
                    let q = jc.point_at(mid);
                    let w2 = q - fs.frame.origin;
                    let ax2 = fs.axis().dot_vec(w2);
                    (w2 - fs.axis().as_vec() * ax2).length() < fs.radius
                };
                let a_p1 = circle_angle(&jc, p1);
                let a_p2 = circle_angle(&jc, p2);
                open_notches[ci].push(OpenNotchRef {
                    edges: vec![e_lat],
                    side: sd,
                    reversed: !hole_rev,
                    gap_circle: Some(jc),
                    bite_ccw_from_start: bite_of(a_p1, a_p2),
                });
                cap_notches.entry(cap_fid).or_default().push(OpenNotchRef {
                    edges: vec![e_cap],
                    side: sd,
                    reversed: !hole_rev,
                    gap_circle: Some(jc),
                    bite_ccw_from_start: bite_of(a_p2, a_p1),
                });
                n_cap_trim += 1;
            }
        }
    }
    let _ = std::mem::take(&mut cap_bites);
    if std::env::var("CADCORE_DUMP_CONT").is_ok() {
        eprintln!("[union][cap] trims={n_cap_trim} bp={c_bp} loops={c_loops} cross={c_cross} norun={c_norun} sliv={c_sliv} lines={c_lines}");
    }
    let _ = (n_cap_trim, c_bp, c_loops, c_cross, c_norun, c_sliv, c_lines);

    // ── Pass 2f: MUTUAL butt-end bite (two capped legs crossing at a corner) ──
    // Two capped legs L1 ⊥ L2 ending at one corner bite EACH OTHER.  The seam
    // is chained from fully ANALYTIC pieces, each shared by exactly two faces:
    //   ll   = lat1∩lat2 (cyl∩cyl, clipped to BOTH slabs)
    //   lc12 = lat1∩cap2 (rulings of L1 in cap2's plane, clipped to disk2+slab1)
    //   lc21 = lat2∩cap1 (symmetric)
    //   cc   = cap1∩cap2 (plane∩plane line, clipped to both disks)
    // Every joint lies on one of the two end circles, so the laterals get
    // multi-edge END notches and the caps trimmed-plane notches — the existing
    // machinery — with gap arcs shared lat↔cap on each junction circle.
    let mut n_mut_trim = 0usize;
    let mut m_pairs = 0usize;
    let mut m_chain_fail = 0usize;
    let mut m_conflict = 0usize;
    {
        let cap_keys: Vec<(usize, BoundarySide)> = cap_map.keys().copied().collect();
        for a in 0..cap_keys.len() {
            for b in (a + 1)..cap_keys.len() {
                let (ci, sdi) = cap_keys[a];
                let (cj, sdj) = cap_keys[b];
                if cyl_faces[ci].solid_name == cyl_faces[cj].solid_name {
                    continue;
                }
                let (FaceBoundary::Circle(jc1), FaceBoundary::Circle(jc2)) = (
                    match sdi {
                        BoundarySide::Start => &cyl_faces[ci].start,
                        BoundarySide::End => &cyl_faces[ci].end,
                    },
                    match sdj {
                        BoundarySide::Start => &cyl_faces[cj].start,
                        BoundarySide::End => &cyl_faces[cj].end,
                    },
                ) else {
                    continue;
                };
                let (jc1, jc2) = (*jc1, *jc2);
                let s1 = cyl_faces[ci].surf;
                let s2 = cyl_faces[cj].surf;
                let (len1, len2) = (cyl_faces[ci].length, cyl_faces[cj].length);
                // Broad-phase: cap centres near each other.
                if (jc1.frame.origin - jc2.frame.origin).length() > 2.0 * (s1.radius + s2.radius)
                {
                    continue;
                }
                // Diagnostic bisection: CADCORE_MUT_ZBAND="lo,hi" commits only
                // pairs with both cap centres inside [lo,hi] (z); with
                // CADCORE_MUT_ZBAND_EXCL=1 the band is skipped instead.
                if let Ok(band) = std::env::var("CADCORE_MUT_ZBAND") {
                    if let Some((lo, hi)) = band.split_once(',').and_then(|(a, b)| {
                        Some((a.trim().parse::<f64>().ok()?, b.trim().parse::<f64>().ok()?))
                    }) {
                        let z1 = jc1.frame.origin.z;
                        let z2 = jc2.frame.origin.z;
                        let inside = z1 >= lo && z1 <= hi && z2 >= lo && z2 <= hi;
                        let excl = std::env::var("CADCORE_MUT_ZBAND_EXCL").is_ok();
                        if inside == excl {
                            continue;
                        }
                    }
                }
                m_pairs += 1;
                let ax1 = |p: Point3| s1.axis().dot_vec(p - s1.frame.origin);
                let ax2 = |p: Point3| s2.axis().dot_vec(p - s2.frame.origin);
                let in1 = |p: Point3| ax1(p) > -1e-9 && ax1(p) < len1 + 1e-9;
                let in2 = |p: Point3| ax2(p) > -1e-9 && ax2(p) < len2 + 1e-9;
                let inside_solid1 = |p: Point3| -> bool {
                    let w = p - s1.frame.origin;
                    let t = s1.axis().dot_vec(w);
                    (w - s1.axis().as_vec() * t).length() < s1.radius && t > 0.0 && t < len1
                };
                let inside_solid2 = |p: Point3| -> bool {
                    let w = p - s2.frame.origin;
                    let t = s2.axis().dot_vec(w);
                    (w - s2.axis().as_vec() * t).length() < s2.radius && t > 0.0 && t < len2
                };
                // piece = (pts, face_a, face_b) where faces: 0=lat1 1=cap1 2=lat2 3=cap2
                let mut pieces: Vec<(Vec<Point3>, usize, usize)> = Vec::new();
                // ll: cyl∩cyl loops clipped to both slabs.
                for lp in cyl_cyl_intersection(&s1, &s2, SAMPLES) {
                    if !lp.closed || lp.points.len() < 8 {
                        continue;
                    }
                    let pn = &lp.points;
                    let n = pn.len();
                    let ok = |p: Point3| in1(p) && in2(p);
                    if pn.iter().all(|&p| ok(p)) {
                        continue; // fully interior crossing — Pass 2 territory
                    }
                    if !pn.iter().any(|&p| ok(p)) {
                        continue;
                    }
                    // maximal inside-runs with interpolated crossings
                    let cross = |pa: Point3, pb: Point3| -> Point3 {
                        // first constraint that flips decides t
                        let mut t_best = 1.0f64;
                        for (f, lim) in [(0usize, 0.0f64), (0, len1), (1, 0.0), (1, len2)] {
                            let (va, vb) = if f == 0 {
                                (ax1(pa) - lim, ax1(pb) - lim)
                            } else {
                                (ax2(pa) - lim, ax2(pb) - lim)
                            };
                            if (va < 0.0) != (vb < 0.0) {
                                let t = va / (va - vb);
                                if t > 1e-9 && t < t_best {
                                    t_best = t;
                                }
                            }
                        }
                        pa + (pb - pa) * t_best.clamp(0.0, 1.0)
                    };
                    let start = match (0..n).find(|&i| !ok(pn[i]) && ok(pn[(i + 1) % n])) {
                        Some(s) => s,
                        None => continue,
                    };
                    let mut i = start;
                    let mut guard = 0;
                    loop {
                        // enter
                        let mut run: Vec<Point3> = vec![cross(pn[i % n], pn[(i + 1) % n])];
                        let mut j = (i + 1) % n;
                        while ok(pn[j]) {
                            run.push(pn[j]);
                            let nx = (j + 1) % n;
                            if !ok(pn[nx]) {
                                run.push(cross(pn[j], pn[nx]));
                                break;
                            }
                            j = nx;
                        }
                        if run.len() >= 3 {
                            pieces.push((run, 0, 2)); // lat1 + lat2
                        }
                        // advance to next entry
                        i = j;
                        loop {
                            i = (i + 1) % n;
                            guard += 1;
                            if guard > 2 * n {
                                break;
                            }
                            if !ok(pn[i % n]) && ok(pn[(i + 1) % n]) {
                                break;
                            }
                        }
                        if guard > 2 * n || i == start {
                            break;
                        }
                    }
                }
                // lc12: rulings of L1 in cap2's plane, clipped to disk2 ∧ slab1.
                let plane2 = Plane3::from_origin_normal(jc2.frame.origin, jc2.frame.z);
                let plane1 = Plane3::from_origin_normal(jc1.frame.origin, jc1.frame.z);
                // cyl_plane_intersection treats the plane as parallel to the
                // axis only below cos 1e-6; path data carries ~1e-5 tilt noise,
                // which flips the result to a degenerate giant ellipse.  Up to
                // cos 1e-3 (~0.06°) project the axis into the plane and build
                // the ruling lines directly — error well under a micron.
                let near_parallel_lines = |s: &CylSurf, plane: &Plane3| -> CylPlaneCurve {
                    let n = plane.normal();
                    let d = s.axis().as_vec();
                    let cos = n.dot_vec(d).abs();
                    if cos >= 1e-3 {
                        return cyl_plane_intersection(s, plane);
                    }
                    let dist = plane.signed_distance(s.frame.origin);
                    if dist.abs() > s.radius + 1e-12 {
                        return CylPlaneCurve::None;
                    }
                    let dp = d - n.as_vec() * n.dot_vec(d);
                    let Some(dpn) = UnitVec3::try_from_vec(dp) else {
                        return CylPlaneCurve::None;
                    };
                    let Some(perp) = UnitVec3::try_from_vec(n.as_vec().cross(dpn.as_vec()))
                    else {
                        return CylPlaneCurve::None;
                    };
                    let foot = s.frame.origin - n.as_vec() * dist;
                    let half = (s.radius * s.radius - dist * dist).max(0.0).sqrt();
                    let mut lines = Vec::new();
                    if half < 1e-12 {
                        lines.push(Line3::new(foot, dpn));
                    } else {
                        lines.push(Line3::new(foot + perp.as_vec() * half, dpn));
                        lines.push(Line3::new(foot - perp.as_vec() * half, dpn));
                    }
                    CylPlaneCurve::Lines(lines)
                };
                let mut add_line_pieces =
                    |conic: CylPlaneCurve, disk_c: Point3, disk_r: f64, axf: &dyn Fn(Point3) -> f64, len: f64, fa: usize, fb: usize, pieces: &mut Vec<(Vec<Point3>, usize, usize)>| {
                        let dbg = std::env::var("CADCORE_DUMP_MUT").is_ok();
                        let CylPlaneCurve::Lines(ls) = conic else {
                            if dbg {
                                eprintln!("[mut-lc] fa={fa} fb={fb} conic_not_lines");
                            }
                            return;
                        };
                        for l in ls {
                            // |o + t d − c|² ≤ R²  ∧  0 ≤ ax(o + t d) ≤ len
                            let oc = l.origin - disk_c;
                            let bh = l.direction.dot_vec(oc);
                            let c0 = oc.dot(oc) - disk_r * disk_r;
                            let disc = bh * bh - c0;
                            if disc <= 1e-12 {
                                if dbg {
                                    eprintln!("[mut-lc] fa={fa} fb={fb} disc={disc:.3e} miss_disk");
                                }
                                continue;
                            }
                            let sq = disc.sqrt();
                            let (mut t0, mut t1) = (-bh - sq, -bh + sq);
                            // slab on ax: ax is affine in t
                            let a_at = |t: f64| axf(l.origin + l.direction * t);
                            let (a0, a1) = (a_at(t0), a_at(t1));
                            let da = a1 - a0;
                            if da.abs() > 1e-12 {
                                let t_for = |aval: f64| t0 + (t1 - t0) * (aval - a0) / da;
                                let (mut lo, mut hi) = (t_for(0.0), t_for(len));
                                if lo > hi {
                                    std::mem::swap(&mut lo, &mut hi);
                                }
                                t0 = t0.max(lo);
                                t1 = t1.min(hi);
                            } else if a0 < 0.0 || a0 > len {
                                if dbg {
                                    eprintln!("[mut-lc] fa={fa} fb={fb} a0={a0:.4} len={len:.3} off_slab");
                                }
                                continue;
                            }
                            if t1 - t0 < 5e-3 {
                                if dbg {
                                    eprintln!(
                                        "[mut-lc] fa={fa} fb={fb} span={:.4e} a0={a0:.4} a1={a1:.4} len={len:.3} too_short",
                                        t1 - t0
                                    );
                                }
                                continue;
                            }
                            let nseg = 8usize;
                            let run: Vec<Point3> = (0..=nseg)
                                .map(|k| l.origin + l.direction * (t0 + (t1 - t0) * k as f64 / nseg as f64))
                                .collect();
                            pieces.push((run, fa, fb));
                        }
                    };
                if std::env::var("CADCORE_DUMP_MUT").is_ok() {
                    eprintln!(
                        "[mut-cos] pair=({ci},{sdi:?})x({cj},{sdj:?}) cos12={:.3e} cos21={:.3e}",
                        plane2.normal().dot(s1.axis()).abs(),
                        plane1.normal().dot(s2.axis()).abs()
                    );
                }
                add_line_pieces(
                    near_parallel_lines(&s1, &plane2),
                    jc2.frame.origin,
                    jc2.radius,
                    &ax1,
                    len1,
                    0,
                    3,
                    &mut pieces,
                );
                add_line_pieces(
                    near_parallel_lines(&s2, &plane1),
                    jc1.frame.origin,
                    jc1.radius,
                    &ax2,
                    len2,
                    2,
                    1,
                    &mut pieces,
                );
                // cc: plane∩plane line clipped to both disks.
                {
                    let n1 = jc1.frame.z.as_vec();
                    let n2 = jc2.frame.z.as_vec();
                    let dir = n1.cross(n2);
                    if dir.length() > 1e-6 {
                        let dirn = dir.normalize();
                        // a point on both planes: solve via projection
                        let d1 = n1.dot(jc1.frame.origin - Point3::new(0.0, 0.0, 0.0));
                        let d2 = n2.dot(jc2.frame.origin - Point3::new(0.0, 0.0, 0.0));
                        let n1n2 = n1.dot(n2);
                        let det = 1.0 - n1n2 * n1n2;
                        let c1 = (d1 - d2 * n1n2) / det;
                        let c2 = (d2 - d1 * n1n2) / det;
                        let p0 = Point3::new(
                            n1.x * c1 + n2.x * c2,
                            n1.y * c1 + n2.y * c2,
                            n1.z * c1 + n2.z * c2,
                        );
                        let clip = |c: Point3, r: f64, t0: f64, t1: f64| -> Option<(f64, f64)> {
                            let oc = p0 + dirn * 0.0 - c;
                            let bh = dirn.dot(oc);
                            let c0 = oc.dot(oc) - r * r;
                            let disc = bh * bh - c0;
                            if disc <= 1e-12 {
                                return None;
                            }
                            let sq = disc.sqrt();
                            Some(((-bh - sq).max(t0), (-bh + sq).min(t1)))
                        };
                        if let Some((t0, t1)) =
                            clip(jc1.frame.origin, jc1.radius, -1e9, 1e9).and_then(|(a0, a1)| {
                                clip(jc2.frame.origin, jc2.radius, a0, a1)
                            })
                        {
                            if t1 - t0 > 5e-3 {
                                let nseg = 6usize;
                                let run: Vec<Point3> = (0..=nseg)
                                    .map(|k| p0 + dirn * (t0 + (t1 - t0) * k as f64 / nseg as f64))
                                    .collect();
                                pieces.push((run, 1, 3));
                            }
                        }
                    }
                }
                if pieces.len() < 2 {
                    continue;
                }
                // Chain pieces into ONE closed loop by endpoint proximity.
                let tolj = 2e-2;
                let mut order: Vec<(usize, bool)> = vec![(0, true)];
                let mut used = vec![false; pieces.len()];
                used[0] = true;
                let mut endp = *pieces[0].0.last().unwrap();
                let startp = pieces[0].0[0];
                let mut closed_chain = false;
                for _ in 0..pieces.len() {
                    if (endp - startp).length() < tolj && order.len() == pieces.len() {
                        closed_chain = true;
                        break;
                    }
                    let mut found = None;
                    for (k, p) in pieces.iter().enumerate() {
                        if used[k] {
                            continue;
                        }
                        if (p.0[0] - endp).length() < tolj {
                            found = Some((k, true));
                            break;
                        }
                        if (*p.0.last().unwrap() - endp).length() < tolj {
                            found = Some((k, false));
                            break;
                        }
                    }
                    let Some((k, fwd)) = found else { break };
                    used[k] = true;
                    endp = if fwd {
                        *pieces[k].0.last().unwrap()
                    } else {
                        pieces[k].0[0]
                    };
                    order.push((k, fwd));
                }
                if !closed_chain {
                    closed_chain =
                        order.len() == pieces.len() && (endp - startp).length() < tolj;
                }
                if !closed_chain {
                    if std::env::var("CADCORE_DUMP_MUT").is_ok() {
                        let dump: Vec<String> = pieces
                            .iter()
                            .map(|(run, fa, fb)| {
                                let a = run[0];
                                let b = *run.last().unwrap();
                                format!(
                                    "({fa},{fb}) n={} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                                    run.len(),
                                    a.x, a.y, a.z, b.x, b.y, b.z
                                )
                            })
                            .collect();
                        eprintln!(
                            "[mut-chainfail] pair=({ci},{sdi:?})x({cj},{sdj:?}) pieces={} ordered={} endp=({:.3},{:.3},{:.3}) startp=({:.3},{:.3},{:.3})\n  {}",
                            pieces.len(),
                            order.len(),
                            endp.x, endp.y, endp.z,
                            startp.x, startp.y, startp.z,
                            dump.join("\n  ")
                        );
                    }
                    m_chain_fail += 1;
                    continue;
                }
                // Conflict guard: no piece may cross existing trims on its faces.
                let face_ci = |f: usize| -> Option<usize> {
                    match f {
                        0 => Some(ci),
                        2 => Some(cj),
                        _ => None,
                    }
                };
                let mut conflict = false;
                for &(k, _) in &order {
                    let (run, fa, fb) = &pieces[k];
                    for f in [*fa, *fb] {
                        let Some(ck) = face_ci(f) else { continue };
                        let sck = &cyl_faces[ck].surf;
                        let uvk = |p: Point3| -> (f64, f64) {
                            let w = p - sck.frame.origin;
                            let t = sck.axis().dot_vec(w);
                            (sck.frame.y.dot_vec(w).atan2(sck.frame.x.dot_vec(w)), t)
                        };
                        let chain_uvk = |pp: &[Point3]| -> Vec<(f64, f64)> {
                            let mut out: Vec<(f64, f64)> = Vec::new();
                            for &p in pp {
                                let (mut u, v) = uvk(p);
                                if let Some(&(pu, _)) = out.last() {
                                    while u - pu > PI {
                                        u -= 2.0 * PI;
                                    }
                                    while u - pu < -PI {
                                        u += 2.0 * PI;
                                    }
                                }
                                out.push((u, v));
                            }
                            out
                        };
                        let seg_int2 = |a1: (f64, f64), a2: (f64, f64), b1: (f64, f64), b2: (f64, f64)| -> bool {
                            let d1 = (a2.0 - a1.0, a2.1 - a1.1);
                            let d2 = (b2.0 - b1.0, b2.1 - b1.1);
                            let den = d1.0 * d2.1 - d1.1 * d2.0;
                            if den.abs() < 1e-15 {
                                return false;
                            }
                            let t = ((b1.0 - a1.0) * d2.1 - (b1.1 - a1.1) * d2.0) / den;
                            let u = ((b1.0 - a1.0) * d1.1 - (b1.1 - a1.1) * d1.0) / den;
                            t > 1e-9 && t < 1.0 - 1e-9 && u > 1e-9 && u < 1.0 - 1e-9
                        };
                        let me = chain_uvk(run);
                        let mut check_edges: Vec<EdgeId> = Vec::new();
                        for h in &holes[ck] {
                            check_edges.push(h.e0);
                            check_edges.push(h.e1);
                        }
                        for mh in &multi_holes[ck] {
                            check_edges.extend(mh.edges.iter().copied());
                        }
                        for nt in &open_notches[ck] {
                            check_edges.extend(nt.edges.iter().copied());
                        }
                        for e in check_edges {
                            let p3 = brep
                                .edges
                                .get(e)
                                .map(|edge| match &edge.geom {
                                    EdgeGeom::Polyline(p) => p.clone(),
                                    _ => Vec::new(),
                                })
                                .unwrap_or_default();
                            if p3.len() < 2 {
                                continue;
                            }
                            let other = chain_uvk(&p3);
                            'outer: for da in [-2.0 * PI, 0.0, 2.0 * PI] {
                                for ia in 0..me.len() - 1 {
                                    for ib in 0..other.len() - 1 {
                                        if seg_int2(
                                            (me[ia].0 + da, me[ia].1),
                                            (me[ia + 1].0 + da, me[ia + 1].1),
                                            other[ib],
                                            other[ib + 1],
                                        ) {
                                            conflict = true;
                                            break 'outer;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if conflict {
                    m_conflict += 1;
                    continue;
                }
                // Weld chain joints + refine every chain point EXACTLY onto its
                // incident surfaces.  The raw pieces are only approximate: the
                // cyl×cyl polyline is exact on one cylinder and ~2 µm off the
                // other, and chain joints can disagree by several µm.  OCC
                // forgives that; SolidWorks knits at ~1 µm and silently drops
                // the affected faces (filaments disappear).  Cyclic projection
                // onto the 2 (piece interior) or 3 (joint) incident surfaces
                // converges to machine precision for these near-orthogonal
                // surface sets.
                let proj_face = |p: Point3, f: usize| -> Point3 {
                    match f {
                        0 | 2 => {
                            let s = if f == 0 { &s1 } else { &s2 };
                            let w = p - s.frame.origin;
                            let ax = s.axis().dot_vec(w);
                            let rad = w - s.axis().as_vec() * ax;
                            let rl = rad.length();
                            if rl < 1e-12 {
                                p
                            } else {
                                s.frame.origin + s.axis().as_vec() * ax + rad * (s.radius / rl)
                            }
                        }
                        1 | 3 => {
                            let pl = if f == 1 { &plane1 } else { &plane2 };
                            p - pl.normal().as_vec() * pl.signed_distance(p)
                        }
                        _ => p,
                    }
                };
                let refine = |mut p: Point3, fs: &[usize]| -> Point3 {
                    for _ in 0..40 {
                        let q0 = p;
                        for &f in fs {
                            p = proj_face(p, f);
                        }
                        if (p - q0).length() < 1e-12 {
                            break;
                        }
                    }
                    p
                };
                let mut pps: Vec<(usize, usize, Vec<Point3>)> = Vec::new();
                for &(k, fwd) in &order {
                    let (run, fa, fb) = &pieces[k];
                    let mut pp: Vec<Point3> = if fwd {
                        run.clone()
                    } else {
                        run.iter().rev().copied().collect()
                    };
                    for p in pp.iter_mut() {
                        *p = refine(*p, &[*fa, *fb]);
                    }
                    pps.push((*fa, *fb, pp));
                }
                let nm = pps.len();
                for i in 0..nm {
                    let j = (i + 1) % nm;
                    let a = *pps[i].2.last().unwrap();
                    let b = pps[j].2[0];
                    // joint lies on the union of both pieces' incident faces
                    let mut fs: Vec<usize> = vec![pps[i].0, pps[i].1];
                    for f in [pps[j].0, pps[j].1] {
                        if !fs.contains(&f) {
                            fs.push(f);
                        }
                    }
                    let w = refine(a + (b - a) * 0.5, &fs);
                    *pps[i].2.last_mut().unwrap() = w;
                    pps[j].2[0] = w;
                }
                // Build edges in CHAIN order (stored direction = chain dir).
                let mut chain_edges: Vec<(EdgeId, usize, usize, Vec<Point3>)> = Vec::new();
                for (fa, fb, pp) in pps {
                    let e = build_single_edge(brep, &pp, &mut vertex_map, &mut edge_map, tolerance);
                    chain_edges.push((e, fa, fb, pp));
                }
                // Per-face runs (consecutive chain pieces touching the face).
                let touches = |idx: usize, f: usize| -> bool {
                    chain_edges[idx].1 == f || chain_edges[idx].2 == f
                };
                let m = chain_edges.len();
                for f in 0..4usize {
                    // find maximal circular runs of pieces on face f
                    let any: Vec<usize> = (0..m).filter(|&i| touches(i, f)).collect();
                    if any.is_empty() {
                        continue;
                    }
                    let allf = any.len() == m;
                    let mut starts: Vec<usize> = Vec::new();
                    if allf {
                        starts.push(0);
                    } else {
                        for i in 0..m {
                            if touches(i, f) && !touches((i + m - 1) % m, f) {
                                starts.push(i);
                            }
                        }
                    }
                    for &st in &starts {
                        let mut run_edges: Vec<EdgeId> = Vec::new();
                        let mut run_pts: Vec<Point3> = Vec::new();
                        let mut i = st;
                        loop {
                            run_edges.push(chain_edges[i].0);
                            run_pts.extend(chain_edges[i].3.iter().copied());
                            let nx = (i + 1) % m;
                            if !touches(nx, f) || (allf && nx == st) {
                                break;
                            }
                            i = nx;
                            if run_edges.len() > m {
                                break;
                            }
                        }
                        // Traversal: material (outside the OTHER solid) on left.
                        let other_inside: &dyn Fn(Point3) -> bool = if f <= 1 {
                            &inside_solid2
                        } else {
                            &inside_solid1
                        };
                        let out1 = s1.axis().as_vec()
                            * if ax1(jc1.frame.origin) < len1 * 0.5 { -1.0 } else { 1.0 };
                        let out2 = s2.axis().as_vec()
                            * if ax2(jc2.frame.origin) < len2 * 0.5 { -1.0 } else { 1.0 };
                        let mut okc = 0usize;
                        let mut nn = 0usize;
                        let step = (run_pts.len() / 9).max(1);
                        for k2 in (0..run_pts.len().saturating_sub(1)).step_by(step) {
                            let p = run_pts[k2];
                            let t = run_pts[k2 + 1] - p;
                            if t.length() < 1e-9 {
                                continue;
                            }
                            let nrm = match f {
                                0 => face_outward_normal(&FaceGeom::Cylinder(s1), p),
                                2 => face_outward_normal(&FaceGeom::Cylinder(s2), p),
                                1 => out1,
                                _ => out2,
                            };
                            let left = nrm.cross(t).normalize();
                            nn += 1;
                            if !other_inside(p + left * 0.04) {
                                okc += 1;
                            }
                        }
                        let same = nn == 0 || okc * 2 >= nn;
                        // jc of this face's notch:
                        let jcf = if f <= 1 { jc1 } else { jc2 };
                        let other_solid_mid = |a1: f64, a2: f64| -> bool {
                            let mid = a1 + (a2 - a1).rem_euclid(2.0 * PI) * 0.5;
                            other_inside(jcf.point_at(mid))
                        };
                        let pa = run_pts[0];
                        let pb = *run_pts.last().unwrap();
                        let aa = circle_angle(&jcf, pa);
                        let ab = circle_angle(&jcf, pb);
                        let nref = OpenNotchRef {
                            edges: run_edges,
                            side: if f == 0 {
                                sdi
                            } else if f == 2 {
                                sdj
                            } else {
                                sdi // unused for caps
                            },
                            reversed: !same,
                            gap_circle: Some(jcf),
                            bite_ccw_from_start: other_solid_mid(aa, ab),
                        };
                        if std::env::var("CADCORE_DUMP_MUT").is_ok() {
                            let a = run_pts[0];
                            let b = *run_pts.last().unwrap();
                            eprintln!(
                                "[mut-push] f={f} run_edges={} same={same} st={st} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                                nref.edges.len(),
                                a.x, a.y, a.z, b.x, b.y, b.z
                            );
                        }
                        match f {
                            0 => open_notches[ci].push(nref),
                            2 => open_notches[cj].push(nref),
                            1 => cap_notches.entry(cap_map[&(ci, sdi)]).or_default().push(nref),
                            _ => cap_notches.entry(cap_map[&(cj, sdj)]).or_default().push(nref),
                        }
                    }
                }
                if std::env::var("CADCORE_DUMP_MUT").is_ok() {
                    let kinds: Vec<String> = chain_edges
                        .iter()
                        .map(|(_, fa, fb, pp)| {
                            let a = pp[0];
                            let b = *pp.last().unwrap();
                            format!(
                                "({fa},{fb}) ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                                a.x, a.y, a.z, b.x, b.y, b.z
                            )
                        })
                        .collect();
                    eprintln!(
                        "[mut] pair=({ci},{sdi:?})x({cj},{sdj:?}) pieces={} chain: {}",
                        pieces.len(),
                        kinds.join(" | ")
                    );
                }
                n_mut_trim += 1;
            }
        }
    }
    if std::env::var("CADCORE_DUMP_CONT").is_ok() {
        eprintln!(
            "[union][mutual] trims={n_mut_trim} pairs={m_pairs} chain_fail={m_chain_fail} conflict={m_conflict}"
        );
    }
    let _ = (n_mut_trim, m_pairs, m_chain_fail, m_conflict);

    // ── Pass 3: build the fused shell ─────────────────────────────────────────
    let solid_id = brep.add_solid(Solid {
        shells: vec![],
        name: Some("scaffold_union".to_string()),
    });
    let shell_id = brep.add_shell(Shell {
        faces: vec![],
        is_outer: true,
        solid: solid_id,
    });

    let mut new_faces: Vec<FaceId> = Vec::new();

    // Cylinders: trimmed (if they have windows) or re-homed unchanged.
    for (idx, cf) in cyl_faces.iter().enumerate() {
        if !collars[idx].is_empty() {
            // Crossing leg wrapped by U-turn collars: split into kept bands.
            let bands = build_collared_leg(
                brep,
                cf,
                &std::mem::take(&mut collars[idx]),
                &std::mem::take(&mut holes[idx]),
                &std::mem::take(&mut multi_holes[idx]),
                &std::mem::take(&mut open_notches[idx]),
                &filament_tubes,
                shell_id,
                &mut vertex_map,
                &mut edge_map,
                tolerance,
            );
            if bands.is_empty() {
                if let Some(face) = brep.faces.get_mut(cf.face_id) {
                    face.shell = shell_id;
                }
                new_faces.push(cf.face_id);
            } else {
                new_faces.extend(bands);
                brep.faces.remove(cf.face_id);
            }
        } else if !through_pieces[idx].is_empty() {
            // Fully-crossed connector: replaced by its kept strips.
            let strips = match (&cf.start, &cf.end) {
                (FaceBoundary::Circle(a), FaceBoundary::Circle(b)) => build_through_strips(
                    brep,
                    *a,
                    *b,
                    FaceGeom::Cylinder(cf.surf),
                    true,
                    &std::mem::take(&mut through_pieces[idx]),
                    &filament_tubes,
                    shell_id,
                    &mut vertex_map,
                    &mut edge_map,
                    tolerance,
                ),
                _ => Vec::new(),
            };
            if strips.is_empty() {
                // Could not assemble — keep the untouched template (safe).
                if let Some(face) = brep.faces.get_mut(cf.face_id) {
                    face.shell = shell_id;
                }
                new_faces.push(cf.face_id);
            } else {
                new_faces.extend(strips);
                brep.faces.remove(cf.face_id);
            }
        } else if holes[idx].is_empty()
            && open_notches[idx].is_empty()
            && multi_holes[idx].is_empty()
        {
            // Untouched cylinder — keep the original template face, re-home it.
            if let Some(face) = brep.faces.get_mut(cf.face_id) {
                face.shell = shell_id;
            }
            new_faces.push(cf.face_id);
        } else {
            let new_id = build_trimmed_swept(
                brep,
                FaceGeom::Cylinder(cf.surf),
                &cf.start,
                &cf.end,
                // Face interior: the axis midpoint (material-side reference
                // for the boundary-notch walk direction).
                cf.surf.frame.origin + cf.surf.axis() * (cf.length * 0.5),
                &holes[idx].clone(),
                &open_notches[idx].clone(),
                &std::mem::take(&mut multi_holes[idx]),
                shell_id,
                &mut vertex_map,
                &mut edge_map,
                tolerance,
            );
            new_faces.push(new_id);
            brep.faces.remove(cf.face_id); // drop the old full-cylinder template
        }
    }
    // Torus elbows: trimmed (if a leg bit them) or re-homed unchanged.
    for (idx, tf) in torus_faces.iter().enumerate() {
        if !torus_through_pieces[idx].is_empty() {
            // Elbow fully crossed by corner bites: replaced by strip bands.
            let strips = build_through_strips(
                brep,
                tf.start_circle,
                tf.end_circle,
                FaceGeom::Torus(tf.surf),
                false,
                &std::mem::take(&mut torus_through_pieces[idx]),
                &filament_tubes,
                shell_id,
                &mut vertex_map,
                &mut edge_map,
                tolerance,
            );
            if strips.is_empty() {
                if let Some(face) = brep.faces.get_mut(tf.face_id) {
                    face.shell = shell_id;
                }
                new_faces.push(tf.face_id);
            } else {
                new_faces.extend(strips);
                brep.faces.remove(tf.face_id);
            }
            continue;
        }
        if torus_holes[idx].is_empty()
            && torus_open_notches[idx].is_empty()
            && torus_multi_holes[idx].is_empty()
        {
            if let Some(face) = brep.faces.get_mut(tf.face_id) {
                face.shell = shell_id;
            }
            new_faces.push(tf.face_id);
        } else {
            let theta_mid = (tf.theta_lo + tf.theta_hi) * 0.5;
            let arc_mid = tf.surf.frame.origin
                + (tf.surf.frame.x.as_vec() * theta_mid.cos()
                    + tf.surf.frame.y.as_vec() * theta_mid.sin())
                    * tf.surf.major_radius;
            let new_id = build_trimmed_swept(
                brep,
                FaceGeom::Torus(tf.surf),
                &FaceBoundary::Circle(tf.start_circle),
                &FaceBoundary::Circle(tf.end_circle),
                // Face interior: the elbow-arc midpoint on the centerline
                // (unambiguous even for 180° U-turns).
                arc_mid,
                &torus_holes[idx].clone(),
                &torus_open_notches[idx].clone(),
                &std::mem::take(&mut torus_multi_holes[idx]),
                shell_id,
                &mut vertex_map,
                &mut edge_map,
                tolerance,
            );
            new_faces.push(new_id);
            brep.faces.remove(tf.face_id);
        }
    }
    // Caps / other faces: re-home unchanged (or rebuild bitten caps).
    for &fid in &keep_faces {
        if let Some(notches) = cap_notches.remove(&fid) {
            let Some(face) = brep.faces.get(fid).cloned() else { continue };
            let jc = notches[0].gap_circle.unwrap();
            let outer = boundary_multi_notch_loop(
                brep,
                &face.geom,
                &FaceBoundary::Circle(jc),
                jc.frame.origin,
                &notches,
                &mut vertex_map,
                &mut edge_map,
                tolerance,
            );
            if outer.is_none() && std::env::var("CADCORE_DUMP_CONT").is_ok() {
                eprintln!("[cap-build-fail] fid={fid:?} notches={}", notches.len());
            }
            if let Some(outer) = outer {
                let new_id = brep.add_face(Face {
                    geom: face.geom.clone(),
                    normal: face.normal,
                    outer_loop: outer,
                    inner_loops: vec![],
                    shell: shell_id,
                    extent: FaceExtent::Trimmed,
                });
                if let Some(lp) = brep.loops.get_mut(outer) {
                    lp.face = new_id;
                }
                new_faces.push(new_id);
                brep.faces.remove(fid);
                continue;
            }
        }
        if let Some(face) = brep.faces.get_mut(fid) {
            face.shell = shell_id;
        }
        new_faces.push(fid);
    }
    if !cap_notches.is_empty() && std::env::var("CADCORE_DUMP_CONT").is_ok() {
        eprintln!(
            "[cap-unconsumed] {} cap faces never matched keep_faces!",
            cap_notches.len()
        );
    }

    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.faces = new_faces;
    }
    if let Some(sd) = brep.solids.get_mut(solid_id) {
        sd.shells = vec![shell_id];
    }

    // Remove the now-empty original filament shells & solids.
    for sh in old_shells {
        brep.shells.remove(sh);
    }
    for sd in old_solids {
        brep.solids.remove(sd);
    }

    UnionReport {
        solids: brep.solids.len(),
        elbow_broad_phase: n_bp,
        elbow_real_overlaps: n_overlap,
        elbow_overlap_found: n_overlap_found,
        elbow_curves: n_curves,
        elbow_open_curves: n_open_boundary,
        elbow_open_trims: n_open_converted,
        elbow_closed_trims: n_torus_trim.saturating_sub(n_open_converted),
        cylinder_cross_trims: n_trimmed,
        marching_fallbacks: n_marching_fallback,
    }
}

// ── helpers ───────────────────────────────────────────────────────────────────

type GridKey = [i64; 3];

struct VertexMap {
    tolerance: f64,
    grid: std::collections::HashMap<GridKey, Vec<VertexId>>,
}

impl VertexMap {
    fn new(brep: &BRep, tolerance: f64) -> Self {
        let mut map = VertexMap {
            tolerance,
            grid: std::collections::HashMap::new(),
        };
        for (vid, v) in brep.vertices.iter() {
            let kx = (v.point.x / tolerance).round() as i64;
            let ky = (v.point.y / tolerance).round() as i64;
            let kz = (v.point.z / tolerance).round() as i64;
            map.grid.entry([kx, ky, kz]).or_default().push(vid);
        }
        map
    }
}

struct EdgeMap {
    tolerance: f64,
    edges: std::collections::HashMap<(VertexId, VertexId), Vec<EdgeId>>,
}

impl EdgeMap {
    fn new(brep: &BRep, tolerance: f64) -> Self {
        let mut map = EdgeMap {
            tolerance,
            edges: std::collections::HashMap::new(),
        };
        for (eid, edge) in brep.edges.iter() {
            let v0 = edge.v_start;
            let v1 = edge.v_end;
            let key = if v0 <= v1 { (v0, v1) } else { (v1, v0) };
            map.edges.entry(key).or_default().push(eid);
        }
        map
    }
}

fn edges_are_tolerant_equal(
    geom1: &EdgeGeom,
    geom2: &EdgeGeom,
    tolerance: f64,
) -> bool {
    match (geom1, geom2) {
        (EdgeGeom::Line(l1), EdgeGeom::Line(l2)) => {
            let dir_dot = l1.direction.dot(l2.direction).abs();
            let collinear = dir_dot >= 0.9999;
            let dist = l1.dist_sq(l2.origin).sqrt();
            collinear && dist <= tolerance
        }
        (EdgeGeom::Circle(c1), EdgeGeom::Circle(c2)) => {
            let center_dist = (c1.frame.origin - c2.frame.origin).length();
            let radius_diff = (c1.radius - c2.radius).abs();
            let normal_dot = c1.frame.z.dot(c2.frame.z).abs();
            center_dist <= tolerance
                && radius_diff <= tolerance
                && normal_dot >= 0.9999
        }
        (EdgeGeom::Ellipse(e1), EdgeGeom::Ellipse(e2)) => {
            let center_dist = (e1.frame.origin - e2.frame.origin).length();
            let major_diff = (e1.semi_major - e2.semi_major).abs();
            let minor_diff = (e1.semi_minor - e2.semi_minor).abs();
            let normal_dot = e1.frame.z.dot(e2.frame.z).abs();
            center_dist <= tolerance
                && major_diff <= tolerance
                && minor_diff <= tolerance
                && normal_dot >= 0.9999
        }
        (EdgeGeom::Polyline(p1), EdgeGeom::Polyline(p2)) => {
            if p1.len() != p2.len() {
                return false;
            }
            let n = p1.len();
            if n == 0 {
                return true;
            }
            // Check forward match
            let mut forward = true;
            for i in 0..n {
                if (p1[i] - p2[i]).length() > tolerance {
                    forward = false;
                    break;
                }
            }
            if forward {
                return true;
            }
            // Check backward match
            let mut backward = true;
            for i in 0..n {
                if (p1[i] - p2[n - 1 - i]).length() > tolerance {
                    backward = false;
                    break;
                }
            }
            backward
        }
        _ => false,
    }
}

fn get_or_create_vertex(
    brep: &mut BRep,
    point: Point3,
    vertex_map: &mut VertexMap,
    tolerance: f64,
) -> VertexId {
    let kx = (point.x / tolerance).round() as i64;
    let ky = (point.y / tolerance).round() as i64;
    let kz = (point.z / tolerance).round() as i64;

    let mut best_vtx: Option<(VertexId, f64)> = None;
    for dx in -1..=1 {
        for dy in -1..=1 {
            for dz in -1..=1 {
                let neighbor_key = [kx + dx, ky + dy, kz + dz];
                if let Some(v_ids) = vertex_map.grid.get(&neighbor_key) {
                    for &v_id in v_ids {
                        if let Some(v) = brep.vertices.get(v_id) {
                            let dist = (v.point - point).length();
                            if dist <= tolerance {
                                if let Some((_, best_d)) = best_vtx {
                                    if dist < best_d {
                                        best_vtx = Some((v_id, dist));
                                    }
                                } else {
                                    best_vtx = Some((v_id, dist));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if let Some((v_id, _)) = best_vtx {
        v_id
    } else {
        let v_id = brep.add_vertex(Vertex { point });
        vertex_map.grid.entry([kx, ky, kz]).or_default().push(v_id);
        v_id
    }
}

fn ellipse_angle(e: &Ellipse3, p: Point3) -> f64 {
    let d = p - e.frame.origin;
    let x = e.frame.x.dot_vec(d) / e.semi_major;
    let y = e.frame.y.dot_vec(d) / e.semi_minor;
    y.atan2(x)
}

fn get_or_create_edge(
    brep: &mut BRep,
    curve: EdgeGeom,
    v0: VertexId,
    v1: VertexId,
    edge_map: &mut EdgeMap,
    tolerance: f64,
) -> EdgeId {
    let key = if v0 <= v1 { (v0, v1) } else { (v1, v0) };
    if let Some(edge_ids) = edge_map.edges.get(&key) {
        for &eid in edge_ids {
            if let Some(existing_edge) = brep.edges.get(eid) {
                if edges_are_tolerant_equal(&existing_edge.geom, &curve, tolerance) {
                    return eid;
                }
            }
        }
    }

    let (t_start, t_end) = match &curve {
        EdgeGeom::Line(l) => {
            let p0 = brep.vertices.get(v0).unwrap().point;
            let p1 = brep.vertices.get(v1).unwrap().point;
            (l.project(p0), l.project(p1))
        }
        EdgeGeom::Circle(c) => {
            let p0 = brep.vertices.get(v0).unwrap().point;
            let p1 = brep.vertices.get(v1).unwrap().point;
            (circle_angle(c, p0), circle_angle(c, p1))
        }
        EdgeGeom::Ellipse(e) => {
            let p0 = brep.vertices.get(v0).unwrap().point;
            let p1 = brep.vertices.get(v1).unwrap().point;
            (ellipse_angle(e, p0), ellipse_angle(e, p1))
        }
        EdgeGeom::Polyline(_) => (0.0, 1.0),
    };

    let eid = brep.add_edge(Edge {
        geom: curve,
        v_start: v0,
        v_end: v1,
        t_start,
        t_end,
        partner: None,
    });
    edge_map.edges.entry(key).or_default().push(eid);
    eid
}

fn is_electrode(solid: &Solid) -> bool {
    solid
        .name
        .as_deref()
        .map_or(false, |n| n.starts_with("electrode"))
}

/// Representative filament radius (all equal in a scaffold) for SSI step sizing.
fn params_radius_or(cyl_faces: &[CylFace]) -> f64 {
    cyl_faces.first().map(|c| c.surf.radius).unwrap_or(0.275)
}

/// The leg as a (straight) swept tube — the generic-SSI view of a cylinder.
fn swept_from_cyl(cf: &CylFace) -> SweptTubeSurface {
    SweptTubeSurface::new(
        vec![CenterlineSeg::Line {
            p0: cf.surf.frame.origin,
            dir: cf.surf.axis(),
            length: cf.length,
        }],
        cf.surf.radius,
    )
}

/// Build the CONTINUOUS swept tube for a whole filament from its analytic
/// centre-line segments (`Line`/`Arc`).  Intersecting a crossing leg with this
/// (rather than one segment) yields CLOSED bite loops even when the crossing is
/// near a segment junction.
fn swept_tube_from_segs(
    segs: &[crate::sweep::SweepPathSegment],
    minor_radius: f64,
) -> Option<SweptTubeSurface> {
    use crate::sweep::SweepPathSegment as S;
    let mut cl: Vec<CenterlineSeg> = Vec::new();
    for s in segs {
        match s {
            S::Line { start, end } => {
                let d = *end - *start;
                let len = d.length();
                if len < 1e-9 {
                    continue;
                }
                cl.push(CenterlineSeg::Line {
                    p0: *start,
                    dir: UnitVec3::try_from_vec(d)?,
                    length: len,
                });
            }
            S::Arc {
                start,
                end,
                center,
                normal,
            } => {
                let xref = UnitVec3::try_from_vec(*start - *center)?;
                let radius = (*start - *center).length();
                let yref = normal.as_vec().cross(xref.as_vec());
                let ve = *end - *center;
                // CCW angle from xref (=start dir) to end dir about `normal`.
                let mut ang1 = ve.dot(yref).atan2(ve.dot(xref.as_vec()));
                if ang1 <= 1e-9 {
                    ang1 += 2.0 * PI; // keep a forward (CCW) sweep
                }
                cl.push(CenterlineSeg::Arc {
                    center: *center,
                    axis: *normal,
                    xref,
                    radius,
                    ang0: 0.0,
                    ang1,
                });
            }
        }
    }
    if cl.is_empty() {
        None
    } else {
        Some(SweptTubeSurface::new(cl, minor_radius))
    }
}

/// The elbow as an (arc) swept tube — the generic-SSI view of a torus fillet.
fn swept_from_torus(tf: &TorusFace) -> SweptTubeSurface {
    SweptTubeSurface::new(
        vec![CenterlineSeg::Arc {
            center: tf.surf.frame.origin,
            axis: tf.surf.frame.z,
            xref: tf.surf.frame.x,
            radius: tf.surf.major_radius,
            ang0: tf.theta_lo,
            ang1: tf.theta_hi,
        }],
        tf.surf.minor_radius,
    )
}

/// Axial parameter of `p` along the cylinder axis (0 at `frame.origin`).
fn axial(surf: &CylSurf, p: Point3) -> f64 {
    surf.axis().dot_vec(p - surf.frame.origin)
}

/// Signed area of a loop projected into the cylinder's `(θ, v)` parameter
/// domain (θ around the axis — CCW gives the outward normal; v along the axis).
/// Positive ⇒ CCW, negative ⇒ CW.
fn signed_area_theta_v(surf: &CylSurf, pts: &[Point3]) -> f64 {
    let o = surf.frame.origin;
    let ax = surf.axis();
    let fx = surf.frame.x;
    let fy = surf.frame.y;
    let mut uv: Vec<(f64, f64)> = pts
        .iter()
        .map(|p| {
            let d = *p - o;
            (fy.dot_vec(d).atan2(fx.dot_vec(d)), ax.dot_vec(d))
        })
        .collect();
    // Unwrap θ so the localized loop is continuous (no 2π jumps).
    for k in 1..uv.len() {
        while uv[k].0 - uv[k - 1].0 > PI {
            uv[k].0 -= 2.0 * PI;
        }
        while uv[k].0 - uv[k - 1].0 < -PI {
            uv[k].0 += 2.0 * PI;
        }
    }
    let mut a = 0.0;
    for k in 0..uv.len() {
        let (x0, y0) = uv[k];
        let (x1, y1) = uv[(k + 1) % uv.len()];
        a += x0 * y1 - x1 * y0;
    }
    a * 0.5
}

/// `HoleRef::reversed` value that makes a two-edge spanning junction loop run
/// CCW in the crossing cylinder's own `(theta, axial)` coordinates.
///
/// The forward loop is `edge_ee` followed by `edge_ea`.  `HoleRef::reversed`
/// traverses the same two shared edges in the opposite order/sense, so it flips
/// the signed area while preserving the topological edge identity.
fn hole_ref_reversed_for_cylinder_ccw(
    surf: &CylSurf,
    edge_ee_pts: &[Point3],
    edge_ea_pts: &[Point3],
) -> bool {
    let mut pts = Vec::with_capacity(edge_ee_pts.len() + edge_ea_pts.len());
    pts.extend(edge_ee_pts.iter().copied());
    pts.extend(edge_ea_pts.iter().skip(1).copied());
    signed_area_theta_v(surf, &pts) < 0.0
}

/// If every point of `pts` lies on the elbow `tf` (within `tol` of its torus
/// surface AND inside its `[theta_lo, theta_hi]` arc band), return `true`.  Used
/// to route a continuous-filament bite loop to the elbow face it sits on.
fn loop_on_torus_face(tf: &TorusFace, pts: &[Point3], tol: f64) -> bool {
    let t = &tf.surf;
    let (lo, hi) = if tf.theta_lo <= tf.theta_hi {
        (tf.theta_lo, tf.theta_hi)
    } else {
        (tf.theta_hi, tf.theta_lo)
    };
    let band_tol = 1e-3;
    pts.iter().all(|&p| {
        let l = t.frame.to_local_point(p);
        let theta = l.y.atan2(l.x);
        // distance from p to the torus surface
        let (ct, st) = (theta.cos(), theta.sin());
        let spine = cadcore_math::Vec3::new(ct * t.major_radius, st * t.major_radius, 0.0);
        let d = (cadcore_math::Vec3::new(l.x, l.y, l.z) - spine).length();
        let surf_dist = (d - t.minor_radius).abs();
        // θ within the arc band (allow a small margin; handle wrap)
        let mut in_band = theta >= lo - band_tol && theta <= hi + band_tol;
        if !in_band {
            // try the ±2π alias
            let ta = theta + 2.0 * PI;
            let tb = theta - 2.0 * PI;
            in_band = (ta >= lo - band_tol && ta <= hi + band_tol)
                || (tb >= lo - band_tol && tb <= hi + band_tol);
        }
        surf_dist < tol && in_band
    })
}

/// Return the loop points wound **CW** in cylinder `surf`'s `(θ, v)` domain so
/// they form a proper hole (interior removed).
fn orient_hole_for_cylinder(surf: &CylSurf, pts: &[Point3]) -> Vec<Point3> {
    if signed_area_theta_v(surf, pts) > 0.0 {
        pts.iter().rev().copied().collect()
    } else {
        pts.to_vec()
    }
}

/// The loop is relevant to this cylinder segment: its points lie within the
/// `[0, length]` axial extent **expanded by one radius**.  The swept solid
/// extends ~one radius past each centreline endpoint (cap / miter / adjoining
/// elbow), so a crossing near a leg end still produces a valid bite there.
/// Without this tolerance, perimeter crossings were dropped → untrimmed
/// self-intersecting overlaps that break SpaceClaim's stitching (model-
/// dependent, which is why some scaffolds combined and others did not).
fn within_axial(surf: &CylSurf, length: f64, pts: &[Point3]) -> bool {
    let tol = 1e-6;
    pts.iter().all(|&p| {
        let t = axial(surf, p);
        t >= -tol && t <= length + tol
    })
}

fn open_boundary_loops(
    tf: &TorusFace,
    elbow: &SweptTubeSurface,
    curves: Vec<IntersectionPolyline>,
) -> Vec<Vec<Point3>> {
    let mut usable = Vec::new();
    for curve in curves {
        let Some(oriented) = orient_open_curve_u0_to_u1(elbow, &curve.points) else {
            continue;
        };
        let mid = oriented[oriented.len() / 2];
        let (_, v) = elbow.project_point(mid).unwrap_or((0.0, 0.0));
        usable.push((v, oriented));
    }
    usable.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut loops = Vec::new();
    let mut i = 0;
    while i + 1 < usable.len() {
        let a = &usable[i].1;
        let b = &usable[i + 1].1;
        if let Some(loop_pts) = close_open_pair_with_elbow_boundaries(tf, a, b) {
            loops.push(loop_pts);
        }
        i += 2;
    }
    loops
}

fn orient_open_curve_u0_to_u1(elbow: &SweptTubeSurface, pts: &[Point3]) -> Option<Vec<Point3>> {
    if pts.len() < 6 {
        return None;
    }
    let len = elbow.length();
    let tol = len.max(1.0) * 1e-5;
    let u_first = elbow.project_point(*pts.first()?).map(|p| p.0)?;
    let u_last = elbow.project_point(*pts.last()?).map(|p| p.0)?;
    let first_lo = u_first <= tol;
    let first_hi = u_first >= len - tol;
    let last_lo = u_last <= tol;
    let last_hi = u_last >= len - tol;
    if first_lo && last_hi {
        Some(pts.to_vec())
    } else if first_hi && last_lo {
        let mut rev = pts.to_vec();
        rev.reverse();
        Some(rev)
    } else {
        None
    }
}

fn close_open_pair_with_elbow_boundaries(
    tf: &TorusFace,
    a: &[Point3],
    b: &[Point3],
) -> Option<Vec<Point3>> {
    let a0 = *a.first()?;
    let a1 = *a.last()?;
    let b0 = *b.first()?;
    let b1 = *b.last()?;

    let mut pts = Vec::with_capacity(a.len() + b.len() + 32);
    pts.extend_from_slice(a);
    append_circle_arc(&mut pts, &tf.end_circle, a1, b1, 12);
    for p in b.iter().rev() {
        push_distinct(&mut pts, *p);
    }
    append_circle_arc(&mut pts, &tf.start_circle, b0, a0, 12);
    if (pts[0] - *pts.last()?).length() > 1e-7 {
        pts.push(pts[0]);
    }
    if polyline_len(&pts) <= tf.surf.minor_radius * 0.1 {
        return None;
    }
    Some(pts)
}

fn append_circle_arc(out: &mut Vec<Point3>, c: &Circle3, from: Point3, to: Point3, samples: usize) {
    let a0 = circle_angle(c, from);
    let a1 = circle_angle(c, to);
    let mut d = a1 - a0;
    while d > PI {
        d -= 2.0 * PI;
    }
    while d < -PI {
        d += 2.0 * PI;
    }
    let n = samples.max(2);
    for i in 1..=n {
        let t = i as f64 / n as f64;
        push_distinct(out, c.point_at(a0 + d * t));
    }
}

fn circle_angle(c: &Circle3, p: Point3) -> f64 {
    let d = p - c.frame.origin;
    c.frame.y.dot_vec(d).atan2(c.frame.x.dot_vec(d))
}

fn push_distinct(out: &mut Vec<Point3>, p: Point3) {
    if out.last().map_or(true, |q| (p - *q).length() > 1e-8) {
        out.push(p);
    }
}

fn polyline_len(pts: &[Point3]) -> f64 {
    pts.windows(2).map(|w| (w[1] - w[0]).length()).sum()
}

fn loop_on_tube_surface(tube: &SweptTubeSurface, pts: &[Point3], tol: f64) -> bool {
    pts.iter()
        .all(|p| tube.signed_distance(*p).map_or(false, |d| d.abs() <= tol))
}

fn same_boundary_side(tube: &SweptTubeSurface, pts: &[Point3]) -> Option<BoundarySide> {
    let len = tube.length();
    let tol = (len.max(1.0) * 1e-5).max(1e-5);
    let a = tube.project_point(*pts.first()?)?.0;
    let b = tube.project_point(*pts.last()?)?.0;
    if a <= tol && b <= tol {
        Some(BoundarySide::Start)
    } else if a >= len - tol && b >= len - tol {
        Some(BoundarySide::End)
    } else {
        None
    }
}
fn build_open_notch_if_same_boundary(
    curve: &IntersectionPolyline,
    elbow: &SweptTubeSurface,
    leg: &SweptTubeSurface,
    brep: &mut BRep,
    vertex_map: &mut VertexMap,
    edge_map: &mut EdgeMap,
    tolerance: f64,
) -> Option<(EdgeId, BoundarySide, BoundarySide)> {
    if curve.closed || curve.points.len() < 4 {
        return None;
    }
    let elbow_side = same_boundary_side(elbow, &curve.points)?;
    let leg_side = same_boundary_side(leg, &curve.points)?;
    if !loop_on_tube_surface(elbow, &curve.points, 1e-4)
        || !loop_on_tube_surface(leg, &curve.points, 1e-4)
    {
        return None;
    }

    let start = *curve.points.first()?;
    let end = *curve.points.last()?;
    if (end - start).length() <= 1e-7 {
        return None;
    }
    let v0 = get_or_create_vertex(brep, start, vertex_map, tolerance);
    let v1 = get_or_create_vertex(brep, end, vertex_map, tolerance);
    let edge = get_or_create_edge(
        brep,
        EdgeGeom::Polyline(curve.points.clone()),
        v0,
        v1,
        edge_map,
        tolerance,
    );
    Some((edge, elbow_side, leg_side))
}
fn endpoint_domain_label(tube: &SweptTubeSurface, p: Point3) -> String {
    let len = tube.length();
    let tol = (len.max(1.0) * 1e-5).max(1e-5);
    let Some(pr) = tube.project_point_diagnostics(p) else {
        return "project=fail".to_string();
    };
    let u0 = pr.u;
    let v = pr.v.rem_euclid(2.0 * PI);
    let u_hit = if u0 <= tol {
        "u_min"
    } else if u0 >= len - tol {
        "u_max"
    } else {
        "u_inside"
    };
    let seam_tol = 1e-4;
    let v_hit = if v <= seam_tol || (2.0 * PI - v) <= seam_tol {
        "v_seam"
    } else {
        "v_inside"
    };
    format!(
        "{u_hit}/{v_hit} u={:.6}/{:.6} v={:.6} sd={:.3e} clamped={} conv={} it={} center_d={:.3e}",
        pr.u,
        len,
        v,
        pr.surface_distance,
        pr.clamped,
        pr.converged,
        pr.iterations,
        pr.centerline_distance,
    )
}

fn endpoint_boundary_distance(tube: &SweptTubeSurface, p: Point3) -> f64 {
    let Some((u, v)) = tube.project_point(p) else {
        return f64::INFINITY;
    };
    let len = tube.length();
    let u_boundary = u.min((len - u).abs());
    let arc_dist = u_boundary.abs();
    let surface_dist = tube.signed_distance(p).map_or(f64::INFINITY, f64::abs);
    let _ = v;
    (arc_dist * arc_dist + surface_dist * surface_dist).sqrt()
}

fn log_open_curve_diag(
    stage: &str,
    pair_seq: usize,
    ti: usize,
    ci: usize,
    torus_face: FaceId,
    cyl_face: FaceId,
    elbow: &SweptTubeSurface,
    leg: &SweptTubeSurface,
    curve: &IntersectionPolyline,
) {
    let Some(start) = curve.points.first().copied() else {
        return;
    };
    let Some(end) = curve.points.last().copied() else {
        return;
    };
    eprintln!(
        "[union][open:{stage}] pair={} ti={} ci={} tor_face={:?} cyl_face={:?} pts={} closed={} len={:.6}",
        pair_seq,
        ti,
        ci,
        torus_face,
        cyl_face,
        curve.points.len(),
        curve.closed,
        polyline_len(&curve.points)
    );
    eprintln!(
        "[union][open:{stage}]   start={:?} elbow={} leg={} d_boundary(elbow)={:.3e} d_boundary(leg)={:.3e}",
        start,
        endpoint_domain_label(elbow, start),
        endpoint_domain_label(leg, start),
        endpoint_boundary_distance(elbow, start),
        endpoint_boundary_distance(leg, start),
    );
    eprintln!(
        "[union][open:{stage}]   end  ={:?} elbow={} leg={} d_boundary(elbow)={:.3e} d_boundary(leg)={:.3e}",
        end,
        endpoint_domain_label(elbow, end),
        endpoint_domain_label(leg, end),
        endpoint_boundary_distance(elbow, end),
        endpoint_boundary_distance(leg, end),
    );
}

/// Build the two shared polyline edges for a closed intersection loop, split at
/// two distinct seam vertices so the two crossing faces can traverse them in
/// opposite directions (the single-closed-edge sense trick does not work for
/// polylines).  Returns `(e0: seam0→seam1, e1: seam1→seam0)`.
/// Result of splitting a continuous-SSI loop that spans an elbow↔leg junction.
struct JunctionSplit {
    /// Elbow (torus) face index the loop partly lies on.
    ti: usize,
    /// Adjacent leg (cyl) face index sharing the junction circle, if found.
    aj: Option<usize>,
    /// `true` ⇒ junction is the elbow's start circle, else its end circle.
    j_is_start: bool,
    /// Junction minor circle.
    j: Circle3,
    /// Elbow-part curve P1→P2 (snapped endpoints on `j`).
    e_e: Vec<Point3>,
    /// Adjacent-leg-part curve P2→P1 (snapped endpoints on `j`).
    e_a: Vec<Point3>,
}

fn snap_to_circle(c: &Circle3, p: Point3) -> Point3 {
    c.point_at(circle_angle(c, p))
}

/// Rotate a trimmed cylinder face's parameterization about its axis so the
/// θ=0 **seam** lands in the middle of the largest angular gap between the
/// face's trim features (`avoid` = points of its polyline trim edges).
///
/// Why: pcurves are unwrapped per edge and anchored independently; when the
/// seam runs through a hole/strip, edges of one wire land on different ±2π
/// branches and OCC sees the wire as open / self-intersecting in (u,v).
/// Pointing the seam into removed (or feature-free) material keeps every
/// wire of the face on a single branch — no explicit SEAM_EDGE needed.
fn seam_rotated_cyl(surf: &CylSurf, avoid: &[Point3]) -> CylSurf {
    if avoid.is_empty() {
        return *surf;
    }
    let two_pi = 2.0 * PI;
    let mut angs: Vec<f64> = avoid
        .iter()
        .map(|&p| {
            let d = p - surf.frame.origin;
            surf.frame
                .y
                .dot_vec(d)
                .atan2(surf.frame.x.dot_vec(d))
                .rem_euclid(two_pi)
        })
        .collect();
    angs.sort_by(f64::total_cmp);
    let n = angs.len();
    let (mut best_gap, mut best_mid) = (-1.0_f64, 0.0_f64);
    for i in 0..n {
        let a = angs[i];
        let b = if i + 1 < n { angs[i + 1] } else { angs[0] + two_pi };
        if b - a > best_gap {
            best_gap = b - a;
            best_mid = (a + b) * 0.5;
        }
    }
    let (s, c) = best_mid.sin_cos();
    let new_x = surf.frame.x.as_vec() * c + surf.frame.y.as_vec() * s;
    let Some(x) = cadcore_math::UnitVec3::try_from_vec(new_x) else {
        return *surf;
    };
    let Some(y) = cadcore_math::UnitVec3::try_from_vec(surf.frame.z.cross(x)) else {
        return *surf;
    };
    CylSurf {
        frame: cadcore_math::Frame3 {
            origin: surf.frame.origin,
            x,
            y,
            z: surf.frame.z,
        },
        radius: surf.radius,
    }
}

// ── N-face spanning-loop distribution ─────────────────────────────────────────
// A bite loop at a U-turn apex covers THREE faces of the bitten filament
// (elbow + short connector + elbow).  The loop is split at EVERY junction
// crossing; each face gets its pieces (elbow → boundary notch, connector →
// through-cut wire) and the crossing leg gets one N-edge hole.

/// Which face of the bitten filament owns a run of loop points.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum PieceOwner {
    Elbow(usize),
    Leg(usize),
}

/// One face-run of a multi-face spanning loop (endpoints snapped onto the
/// junction circles crossed at its start/end).
struct SpanPiece {
    owner: PieceOwner,
    pts: Vec<Point3>,
    j_in: Circle3,
    j_out: Circle3,
    /// Original-loop index range of the run (inclusive, circular): the run
    /// covers loop points `i0 ..= i1`.  Used to merge the cut sets of the TWO
    /// sides of a corner bite into shared atomic edges.
    i0: usize,
    i1: usize,
}

/// One full-length crossing piece over a short connector face.  A connector
/// crossed by `k` bites carries `2k` of these; the face is rebuilt as the `k`
/// kept STRIPS between bites (piece + clear-arc(J2) + piece + clear-arc(J1)).
#[derive(Clone)]
struct ThroughPiece {
    /// Crossing-piece edges in chain order.
    edges: Vec<EdgeId>,
    /// Pieces of one bite share this id (the strip BETWEEN them is removed).
    bite: usize,
    /// Index into `filament_tubes` of the CROSSING filament (kept/removed
    /// sampling — works for straight-leg and elbow crossings alike).
    crossing: usize,
    p_j1: Point3,
    p_j2: Point3,
    /// `true` ⇒ the kept strip traverses the piece in its stored direction
    /// (the crossing-leg hole takes the opposite).
    strip_same: bool,
}

/// A WINDING bite chain on the crossing leg: the U-turn wraps the leg's full
/// circumference, so the chain is not a hole — it is a "collar" that splits
/// the lateral surface into bands.  Bands alternate inside/outside the U.
#[derive(Clone)]
struct CollarRef {
    /// Piece edges in chain order (stored directions chain head-to-tail).
    edges: Vec<EdgeId>,
    /// Net θ of the stored chain is positive (CCW) on this leg.
    stored_ccw: bool,
    /// Mean axial position (band ordering).
    mean_v: f64,
    /// All chain points (band kept/removed sampling).
    pts: Vec<Point3>,
}

/// Net unwrapped θ travelled by a closed chain on a cylinder (≈ ±2π for a
/// winding collar, ≈ 0 for an ordinary hole).
fn chain_winding_theta(surf: &CylSurf, pts: &[Point3]) -> f64 {
    let o = surf.frame.origin;
    let (fx, fy) = (surf.frame.x, surf.frame.y);
    let mut prev: Option<f64> = None;
    let mut acc = 0.0;
    let mut first = 0.0;
    for p in pts {
        let d = *p - o;
        let mut th = fy.dot_vec(d).atan2(fx.dot_vec(d));
        if let Some(pr) = prev {
            let mut dd = th - pr;
            while dd > PI {
                dd -= 2.0 * PI;
            }
            while dd < -PI {
                dd += 2.0 * PI;
            }
            th = pr + dd;
            acc += dd;
        } else {
            first = th;
        }
        prev = Some(th);
    }
    // closing step
    if let Some(pr) = prev {
        let mut dd = first - pr;
        while dd > PI {
            dd -= 2.0 * PI;
        }
        while dd < -PI {
            dd += 2.0 * PI;
        }
        acc += dd;
    }
    acc
}

/// N-edge hole loop on the crossing leg (generalizes the 2-edge `hole_loop`).
#[derive(Clone)]
struct MultiHoleRef {
    /// Piece edges in chain order (stored directions chain head-to-tail).
    edges: Vec<EdgeId>,
    /// `true` ⇒ traverse the chain backwards (each edge Opposite).
    reversed: bool,
}

fn multi_hole_loop(brep: &mut BRep, h: &MultiHoleRef) -> Option<LoopId> {
    let mut coedges = Vec::with_capacity(h.edges.len());
    let ordered: Vec<EdgeId> = if h.reversed {
        h.edges.iter().rev().copied().collect()
    } else {
        h.edges.clone()
    };
    // Derive each sense from the actual chaining (edge stored direction vs the
    // wire's running endpoint).
    let mut cur: Option<VertexId> = None;
    for &e in &ordered {
        let edge = brep.edges.get(e)?.clone();
        let sense = match cur {
            None => {
                if h.reversed {
                    CoEdgeSense::Opposite
                } else {
                    CoEdgeSense::Same
                }
            }
            Some(v) if v == edge.v_start => CoEdgeSense::Same,
            Some(v) if v == edge.v_end => CoEdgeSense::Opposite,
            Some(_) => return None, // chain broken
        };
        cur = Some(match sense {
            CoEdgeSense::Same => edge.v_end,
            CoEdgeSense::Opposite => edge.v_start,
        });
        coedges.push(brep.add_coedge(CoEdge {
            edge: e,
            sense,
            next: CoEdgeId::default(),
            prev: CoEdgeId::default(),
            loop_id: LoopId::default(),
        }));
    }
    brep.patch_coedge_links(&coedges);
    let loop_id = brep.add_loop(Loop {
        start: coedges[0],
        face: FaceId::default(),
    });
    for &ce in &coedges {
        if let Some(c) = brep.coedges.get_mut(ce) {
            c.loop_id = loop_id;
        }
    }
    Some(loop_id)
}

/// Rebuild a fully-crossed connector face as its kept STRIPS.  Pieces are
/// sorted by angle on the start circle; the region between the two pieces of
/// one bite is removed material, every other gap is a kept strip face:
/// `piece_a + clear-arc(J2) + piece_b(rev) + clear-arc(J1)`.
fn build_through_strips(
    brep: &mut BRep,
    j1: Circle3,
    j2: Circle3,
    geom_base: FaceGeom,
    seam_rotate: bool,
    pieces: &[ThroughPiece],
    tubes: &[(String, SweptTubeSurface, bool)],
    shell_id: ShellId,
    vertex_map: &mut VertexMap,
    edge_map: &mut EdgeMap,
    tolerance: f64,
) -> Vec<FaceId> {
    let mut ps: Vec<&ThroughPiece> = pieces.iter().collect();
    ps.sort_by(|a, b| {
        circle_angle(&j1, a.p_j1)
            .rem_euclid(2.0 * PI)
            .total_cmp(&circle_angle(&j1, b.p_j1).rem_euclid(2.0 * PI))
    });
    let n = ps.len();
    let mut faces = Vec::new();
    if n < 2 || n % 2 != 0 {
        return faces;
    }
    // Inside-any-crossing-tube test (kept/removed gap classification).
    let crossing_idx: Vec<usize> = {
        let mut idx: Vec<usize> = pieces.iter().map(|p| p.crossing).collect();
        idx.sort_unstable();
        idx.dedup();
        idx
    };
    let inside_crossing = |p: Point3| -> bool {
        crossing_idx.iter().any(|&i| {
            tubes
                .get(i)
                .and_then(|(_, t, _)| t.signed_distance(p))
                .map_or(false, |d| d < 0.0)
        })
    };
    let all_j1: Vec<Point3> = ps.iter().map(|p| p.p_j1).collect();
    let all_j2: Vec<Point3> = ps.iter().map(|p| p.p_j2).collect();
    // Gap arc direction from X to Y on a junction circle: prefer the rotation
    // whose interval contains no other piece point AND whose midpoint (ON the
    // junction circle — same sample the notch flags use) is OUTSIDE the
    // crossing tubes (decides the k=1 case, where no blocking points exist).
    let pick_dir = |j: &Circle3, from: Point3, to: Point3, all: &[Point3]| -> f64 {
        let af = circle_angle(j, from);
        let at = circle_angle(j, to);
        for dir in [1.0_f64, -1.0] {
            let blocked = all.iter().any(|&q| {
                let a = circle_angle(j, q);
                (q - from).length() > 1e-6
                    && (q - to).length() > 1e-6
                    && arc_contains(af, at, a, dir)
            });
            if blocked {
                continue;
            }
            let sweep = if dir > 0.0 {
                (at - af).rem_euclid(2.0 * PI)
            } else {
                -((af - at).rem_euclid(2.0 * PI))
            };
            if !inside_crossing(j.point_at(af + sweep * 0.5)) {
                return dir;
            }
        }
        1.0
    };

    for i in 0..n {
        let a = ps[i];
        let b = ps[(i + 1) % n];
        // Kept or removed gap?  Sample the gap middle at the face's axial
        // middle: inside a crossing tube ⇒ removed material (works for k=1,
        // where both adjacencies join the same bite's pieces).
        {
            let af = circle_angle(&j1, a.p_j1);
            let at = circle_angle(&j1, b.p_j1);
            let sweep = (at - af).rem_euclid(2.0 * PI);
            if inside_crossing(j1.point_at(af + sweep * 0.5)) {
                continue; // removed material
            }
        }
        // Run endpoints (edges chain head→tail in stored order).
        let run_ends = |edges: &Vec<EdgeId>| -> Option<(VertexId, VertexId)> {
            let f = brep.edges.get(*edges.first()?)?;
            let l = brep.edges.get(*edges.last()?)?;
            Some((f.v_start, l.v_end))
        };
        let Some((ea_s, ea_e)) = run_ends(&a.edges) else { continue };
        let Some((eb_s, eb_e)) = run_ends(&b.edges) else { continue };
        // Traverse `a` in its required strip direction; note which junction
        // its traversal ENDS on.
        let a_start_is_vstart = a.strip_same;
        let (a_from, a_to, sa) = if a_start_is_vstart {
            (ea_s, ea_e, CoEdgeSense::Same)
        } else {
            (ea_e, ea_s, CoEdgeSense::Opposite)
        };
        let a_to_p = match brep.vertices.get(a_to) {
            Some(v) => v.point,
            None => continue,
        };
        let a_from_p = match brep.vertices.get(a_from) {
            Some(v) => v.point,
            None => continue,
        };
        // Which circle does a's traversal end on?
        let a_ends_on_j2 = (a_to_p - a.p_j2).length() < (a_to_p - a.p_j1).length();
        // b must be traversed so it STARTS on the same circle a ended on.
        let (b_start_p, b_from, b_to, sb) = {
            let bs = match brep.vertices.get(eb_s) {
                Some(v) => v.point,
                None => continue,
            };
            let b_vstart_on_j2 = (bs - b.p_j2).length() < (bs - b.p_j1).length();
            if b_vstart_on_j2 == a_ends_on_j2 {
                (bs, eb_s, eb_e, CoEdgeSense::Same)
            } else {
                let bep = match brep.vertices.get(eb_e) {
                    Some(v) => v.point,
                    None => continue,
                };
                (bep, eb_e, eb_s, CoEdgeSense::Opposite)
            }
        };
        // Consistency with b's required direction: if it disagrees, the strip
        // cannot satisfy both shared-edge constraints — skip (diagnosed).
        let b_required_same = b.strip_same;
        let b_actual_same = sb == CoEdgeSense::Same;
        if b_required_same != b_actual_same {
            if std::env::var("CADCORE_DUMP_CONT").is_ok() {
                eprintln!("[union][strip-skip] sense conflict on connector");
            }
            continue;
        }
        let b_to_p = match brep.vertices.get(b_to) {
            Some(v) => v.point,
            None => continue,
        };
        let (jx, allx) = if a_ends_on_j2 { (j2, &all_j2) } else { (j1, &all_j1) };
        let (jy, ally) = if a_ends_on_j2 { (j1, &all_j1) } else { (j2, &all_j2) };

        let mut coedges = Vec::new();
        let push_run = |brep: &mut BRep, coedges: &mut Vec<CoEdgeId>, edges: &Vec<EdgeId>, sense: CoEdgeSense| {
            let iter: Vec<EdgeId> = match sense {
                CoEdgeSense::Same => edges.clone(),
                CoEdgeSense::Opposite => edges.iter().rev().copied().collect(),
            };
            for e in iter {
                coedges.push(brep.add_coedge(CoEdge {
                    edge: e,
                    sense,
                    next: CoEdgeId::default(),
                    prev: CoEdgeId::default(),
                    loop_id: LoopId::default(),
                }));
            }
        };
        push_run(brep, &mut coedges, &a.edges, sa);
        {
            let dir = pick_dir(&jx, a_to_p, b_start_p, allx);
            let af = circle_angle(&jx, a_to_p);
            let at = circle_angle(&jx, b_start_p);
            add_gap_arc_directed(
                brep, &jx, a_to, af, b_from, at, dir, &mut coedges, vertex_map, edge_map,
                tolerance,
            );
        }
        push_run(brep, &mut coedges, &b.edges, sb);
        {
            let dir = pick_dir(&jy, b_to_p, a_from_p, ally);
            let af = circle_angle(&jy, b_to_p);
            let at = circle_angle(&jy, a_from_p);
            add_gap_arc_directed(
                brep, &jy, b_to, af, a_from, at, dir, &mut coedges, vertex_map, edge_map,
                tolerance,
            );
        }
        if std::env::var("CADCORE_DUMP_STRIP").is_ok() {
            eprintln!("[strip] face j1 o={:?}", j1.frame.origin);
            for &ce in &coedges {
                let c = &brep.coedges[ce];
                let e = &brep.edges[c.edge];
                let (vf, vt) = match c.sense {
                    CoEdgeSense::Same => (e.v_start, e.v_end),
                    CoEdgeSense::Opposite => (e.v_end, e.v_start),
                };
                let pf = brep.vertices[vf].point;
                let pt = brep.vertices[vt].point;
                eprintln!(
                    "   {:?} {:?} ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                    c.edge, c.sense, pf.x, pf.y, pf.z, pt.x, pt.y, pt.z
                );
            }
        }
        brep.patch_coedge_links(&coedges);
        let loop_id = brep.add_loop(Loop {
            start: coedges[0],
            face: FaceId::default(),
        });
        for &ce in &coedges {
            if let Some(c) = brep.coedges.get_mut(ce) {
                c.loop_id = loop_id;
            }
        }
        // Per-strip parameterization: the seam must land in REMOVED material.
        // The strip's wire covers the pieces AND the whole kept gap (its long
        // arcs) — the only wire-free interval is the adjacent removed gap
        // (b → next sorted piece), so point the seam at ITS middle.
        let seam_surf = {
            let nxt = ps[(i + 2) % n];
            let ab = circle_angle(&j1, b.p_j1);
            let an = circle_angle(&j1, nxt.p_j1);
            let gap = if (b.p_j1 - nxt.p_j1).length() < 1e-9 {
                PI // k=1 degenerate: only 2 pieces; removed gap = (b → a)
            } else {
                (an - ab).rem_euclid(2.0 * PI)
            };
            let seam_mid = if n == 2 {
                // k=1: removed gap is (b → a) ascending.
                let aa = circle_angle(&j1, a.p_j1);
                ab + (aa - ab).rem_euclid(2.0 * PI) * 0.5
            } else {
                ab + gap * 0.5
            };
            // Rotate so θ=0 sits at seam_mid: reuse helper with a single avoid
            // point OPPOSITE the desired seam (largest gap centres on seam).
            let p_opposite = j1.point_at(seam_mid + PI);
            match &geom_base {
                FaceGeom::Cylinder(c) => seam_rotated_cyl(c, &[p_opposite]),
                _ => CylSurf::new(j1.frame.origin, j1.frame.z, j1.radius), // unused
            }
        };
        let strip_geom = match (&geom_base, seam_rotate) {
            (FaceGeom::Cylinder(_), true) => FaceGeom::Cylinder(seam_surf),
            _ => geom_base.clone(),
        };
        let face_id = brep.add_face(Face {
            geom: strip_geom,
            normal: FaceNormal::Same,
            outer_loop: loop_id,
            inner_loops: vec![],
            shell: shell_id,
            extent: FaceExtent::Trimmed,
        });
        if let Some(lp) = brep.loops.get_mut(loop_id) {
            lp.face = face_id;
        }
        faces.push(face_id);
    }
    faces
}

/// Loop from a collar chain traversed with net `+θ` (`want_ccw`) or `−θ`.
fn collar_loop(brep: &mut BRep, c: &CollarRef, want_ccw: bool) -> Option<LoopId> {
    multi_hole_loop(
        brep,
        &MultiHoleRef {
            edges: c.edges.clone(),
            reversed: c.stored_ccw != want_ccw,
        },
    )
}

/// Rebuild a crossing leg wrapped by winding U-turn collars as its kept
/// BANDS.  Collars sorted axially split the lateral surface; bands alternate
/// inside/outside the U-tubes — each kept band gets the lower boundary
/// traversed `+θ` (material above) and the upper traversed `−θ`.  Ordinary
/// (non-winding) holes are assigned to their band by axial position.
fn build_collared_leg(
    brep: &mut BRep,
    cf: &CylFace,
    collars: &[CollarRef],
    holes: &[HoleRef],
    multi_holes: &[MultiHoleRef],
    notches: &[OpenNotchRef],
    tubes: &[(String, SweptTubeSurface, bool)],
    shell_id: ShellId,
    vertex_map: &mut VertexMap,
    edge_map: &mut EdgeMap,
    tolerance: f64,
) -> Vec<FaceId> {
    let s = &cf.surf;
    let axial_of = |p: Point3| s.axis().dot_vec(p - s.frame.origin);
    let mut cs: Vec<&CollarRef> = collars.iter().collect();
    cs.sort_by(|a, b| a.mean_v.total_cmp(&b.mean_v));

    // v(θ) of a collar: nearest chain point by wrapped angular distance.
    let theta_of = |p: Point3| {
        let d = p - s.frame.origin;
        s.frame.y.dot_vec(d).atan2(s.frame.x.dot_vec(d))
    };
    let v_at = |c: &CollarRef, th: f64| -> f64 {
        let mut best = (f64::MAX, 0.0);
        for &p in &c.pts {
            let mut dt = (theta_of(p) - th).abs();
            if dt > PI {
                dt = 2.0 * PI - dt;
            }
            if dt < best.0 {
                best = (dt, axial_of(p));
            }
        }
        best.1
    };
    enum Bnd<'a> {
        Start,
        Collar(&'a CollarRef),
        End,
    }
    let mut bounds: Vec<Bnd> = vec![Bnd::Start];
    bounds.extend(cs.iter().map(|c| Bnd::Collar(c)));
    bounds.push(Bnd::End);
    let bnd_v = |b: &Bnd, th: f64| -> f64 {
        match b {
            Bnd::Start => 0.0,
            Bnd::End => cf.length,
            Bnd::Collar(c) => v_at(c, th),
        }
    };
    let inside_other = |p: Point3| -> bool {
        tubes.iter().any(|(name, t, _)| {
            cf.solid_name.as_deref() != Some(name.as_str())
                && t.signed_distance(p).map_or(false, |d| d < -1e-4)
        })
    };
    // Axial mid of a hole's first edge (band assignment).
    let edge_mid_v = |brep: &BRep, e: EdgeId| -> f64 {
        brep.edges
            .get(e)
            .map(|edge| match &edge.geom {
                EdgeGeom::Polyline(pts) => {
                    pts.iter().map(|&p| axial_of(p)).sum::<f64>() / pts.len() as f64
                }
                _ => 0.0,
            })
            .unwrap_or(0.0)
    };

    let mut faces = Vec::new();
    for w in 0..bounds.len() - 1 {
        let lo = &bounds[w];
        let hi = &bounds[w + 1];
        // Kept band?  Majority of mid samples outside every other tube.
        let inside_count = (0..8)
            .filter(|k| {
                let th = *k as f64 * PI / 4.0;
                let vm = (bnd_v(lo, th) + bnd_v(hi, th)) * 0.5;
                inside_other(s.point_at(th, vm))
            })
            .count();
        if inside_count > 4 {
            continue; // removed band
        }
        // Lower boundary loop: +θ (material above).
        let outer = match lo {
            Bnd::Start => {
                let start_notch: Vec<OpenNotchRef> = notches
                    .iter()
                    .filter(|n| n.side == BoundarySide::Start)
                    .cloned()
                    .collect();
                if start_notch.is_empty() {
                    boundary_loop(brep, &cf.start, vertex_map, edge_map, tolerance)
                } else {
                    boundary_multi_notch_loop(
                        brep,
                        &FaceGeom::Cylinder(*s),
                        &cf.start,
                        s.frame.origin + s.axis() * (cf.length * 0.5),
                        &start_notch,
                        vertex_map,
                        edge_map,
                        tolerance,
                    )
                    .unwrap_or_else(|| boundary_loop(brep, &cf.start, vertex_map, edge_map, tolerance))
                }
            }
            Bnd::Collar(c) => match collar_loop(brep, c, true) {
                Some(l) => l,
                None => continue,
            },
            Bnd::End => continue,
        };
        // Upper boundary loop: −θ (material below).
        let inner0 = match hi {
            Bnd::End => {
                let end_notch: Vec<OpenNotchRef> = notches
                    .iter()
                    .filter(|n| n.side == BoundarySide::End)
                    .cloned()
                    .collect();
                if end_notch.is_empty() {
                    boundary_loop(brep, &cf.end, vertex_map, edge_map, tolerance)
                } else {
                    boundary_multi_notch_loop(
                        brep,
                        &FaceGeom::Cylinder(*s),
                        &cf.end,
                        s.frame.origin + s.axis() * (cf.length * 0.5),
                        &end_notch,
                        vertex_map,
                        edge_map,
                        tolerance,
                    )
                    .unwrap_or_else(|| boundary_loop(brep, &cf.end, vertex_map, edge_map, tolerance))
                }
            }
            Bnd::Collar(c) => match collar_loop(brep, c, false) {
                Some(l) => l,
                None => continue,
            },
            Bnd::Start => continue,
        };
        let mut inner = vec![inner0];
        let mut avoid: Vec<Point3> = Vec::new();
        // Exclusive band interval (collar means partition [0, L]).
        let lo_v = match lo {
            Bnd::Start => 0.0,
            Bnd::Collar(c) => c.mean_v,
            Bnd::End => cf.length,
        };
        let hi_v = match hi {
            Bnd::End => cf.length + 1e-6,
            Bnd::Collar(c) => c.mean_v,
            Bnd::Start => 0.0,
        };
        let in_band = |v: f64| -> bool { v >= lo_v && v < hi_v };
        for h in holes {
            if in_band(edge_mid_v(brep, h.e0)) {
                inner.push(hole_loop(brep, h));
                if let Some(EdgeGeom::Polyline(pts)) = brep.edges.get(h.e0).map(|e| e.geom.clone())
                {
                    avoid.extend(pts);
                }
            }
        }
        for mh in multi_holes {
            if in_band(edge_mid_v(brep, mh.edges[0])) {
                if let Some(l) = multi_hole_loop(brep, mh) {
                    inner.push(l);
                }
            }
        }
        for c in collars {
            avoid.extend(c.pts.iter().copied());
        }
        let face_id = brep.add_face(Face {
            geom: FaceGeom::Cylinder(seam_rotated_cyl(s, &avoid)),
            normal: FaceNormal::Same,
            outer_loop: outer,
            inner_loops: inner.clone(),
            shell: shell_id,
            extent: FaceExtent::Trimmed,
        });
        for lp in std::iter::once(outer).chain(inner) {
            if let Some(l) = brep.loops.get_mut(lp) {
                l.face = face_id;
            }
        }
        faces.push(face_id);
    }
    faces
}

/// Split a closed bite loop into face-runs across ALL junctions it crosses.
/// `elbows` / `legs` are the face indices of the BITTEN filament (the loop
/// lies on its tube).  Returns `None` unless every point classifies to a face
/// and every transition maps to a shared junction circle.
fn split_multi_span(
    pts: &[Point3],
    elbows: &[usize],
    legs: &[usize],
    torus_faces: &[TorusFace],
    cyl_faces: &[CylFace],
) -> Option<Vec<SpanPiece>> {
    let n = pts.len();
    if n < 6 {
        return None;
    }
    let on_elbow = |ti: usize, p: Point3| loop_on_torus_face(&torus_faces[ti], std::slice::from_ref(&p), 1e-2);
    let on_leg = |ci: usize, p: Point3| -> bool {
        let cf = &cyl_faces[ci];
        let w = p - cf.surf.frame.origin;
        let ax = cf.surf.axis().dot_vec(w);
        let rad = (w - cf.surf.axis().as_vec() * ax).length();
        (rad - cf.surf.radius).abs() < 1e-2 && ax > -1e-3 && ax < cf.length + 1e-3
    };
    let matches = |o: PieceOwner, p: Point3| -> bool {
        match o {
            PieceOwner::Elbow(ti) => on_elbow(ti, p),
            PieceOwner::Leg(ci) => on_leg(ci, p),
        }
    };
    let classify = |p: Point3, prev: Option<PieceOwner>| -> Option<PieceOwner> {
        // Hysteresis: keep the previous owner while it still matches, so
        // points near a junction (on both faces) don't flip-flop.
        if let Some(o) = prev {
            if matches(o, p) {
                return Some(o);
            }
        }
        for &ti in elbows {
            if on_elbow(ti, p) {
                return Some(PieceOwner::Elbow(ti));
            }
        }
        for &ci in legs {
            if on_leg(ci, p) {
                return Some(PieceOwner::Leg(ci));
            }
        }
        None
    };

    let mut owners: Vec<PieceOwner> = Vec::with_capacity(n);
    let mut prev = None;
    for &p in pts {
        let o = classify(p, prev)?;
        owners.push(o);
        prev = Some(o);
    }
    // Rotate so index 0 starts a fresh run.
    let start = (0..n).find(|&i| owners[i] != owners[(i + n - 1) % n])?;
    // Collect circular runs.
    let mut runs: Vec<(PieceOwner, Vec<Point3>, usize, usize)> = Vec::new();
    for k in 0..n {
        let i = (start + k) % n;
        if runs.last().map_or(true, |(o, _, _, _)| *o != owners[i]) {
            runs.push((owners[i], Vec::new(), i, i));
        }
        let last = runs.last_mut().unwrap();
        last.1.push(pts[i]);
        last.3 = i;
    }
    if runs.len() < 2 {
        return None;
    }
    // Junction circle shared between two consecutive owners.
    let circles_of = |o: PieceOwner| -> Vec<Circle3> {
        match o {
            PieceOwner::Elbow(ti) => {
                vec![torus_faces[ti].start_circle, torus_faces[ti].end_circle]
            }
            PieceOwner::Leg(ci) => {
                let mut v = Vec::new();
                if let FaceBoundary::Circle(c) = cyl_faces[ci].start {
                    v.push(c);
                }
                if let FaceBoundary::Circle(c) = cyl_faces[ci].end {
                    v.push(c);
                }
                v
            }
        }
    };
    let junction_between = |a: PieceOwner, b: PieceOwner, near: Point3| -> Option<Circle3> {
        let ca = circles_of(a);
        let cb = circles_of(b);
        let mut best: Option<(Circle3, f64)> = None;
        for x in &ca {
            for y in &cb {
                if (x.frame.origin - y.frame.origin).length() < 1e-3
                    && (x.radius - y.radius).abs() < 1e-3
                {
                    let d = (x.frame.origin - near).length();
                    if best.map_or(true, |(_, bd)| d < bd) {
                        // CANONICAL instance: prefer the LEG's circle so that
                        // gap arcs built by the elbow's notched wire and by the
                        // connector strips sample identical waypoints (the 90°
                        // grid lives in the circle's own frame).
                        let canonical = if matches!(a, PieceOwner::Leg(_)) {
                            *x
                        } else if matches!(b, PieceOwner::Leg(_)) {
                            *y
                        } else {
                            *x
                        };
                        best = Some((canonical, d));
                    }
                }
            }
        }
        best.map(|(c, _)| c)
    };

    let m = runs.len();
    let mut joints: Vec<(Circle3, Point3)> = Vec::with_capacity(m); // joint AFTER run k
    for k in 0..m {
        let (oa, ra, _, _) = &runs[k];
        let (ob, rb, _, _) = &runs[(k + 1) % m];
        let pa = *ra.last().unwrap();
        let pb = rb[0];
        let mid = Point3::new(
            (pa.x + pb.x) * 0.5,
            (pa.y + pb.y) * 0.5,
            (pa.z + pb.z) * 0.5,
        );
        let j = junction_between(*oa, *ob, mid)?;
        joints.push((j, snap_to_circle(&j, mid)));
    }
    let mut pieces = Vec::with_capacity(m);
    for k in 0..m {
        let (o, run, i0, i1) = &runs[k];
        let (j_in, p_in) = joints[(k + m - 1) % m];
        let (j_out, p_out) = joints[k];
        let mut p = Vec::with_capacity(run.len() + 2);
        p.push(p_in);
        p.extend(run.iter().copied());
        p.push(p_out);
        // On-face guard for the interior.
        if !run.iter().all(|&q| matches(*o, q)) {
            return None;
        }
        pieces.push(SpanPiece {
            owner: *o,
            pts: p,
            j_in,
            j_out,
            i0: *i0,
            i1: *i1,
        });
    }
    Some(pieces)
}

fn boundary_circle_matches(b: &FaceBoundary, j: &Circle3) -> bool {
    match b {
        FaceBoundary::Circle(c) => {
            (c.frame.origin - j.frame.origin).length() < 1e-3 && (c.radius - j.radius).abs() < 1e-3
        }
        FaceBoundary::Ellipse(_) => false,
    }
}
/// Build a single polyline EDGE from `pts` (start→end vertices).
fn build_single_edge(
    brep: &mut BRep,
    pts: &[Point3],
    vertex_map: &mut VertexMap,
    edge_map: &mut EdgeMap,
    tolerance: f64,
) -> EdgeId {
    if std::env::var("CADCORE_PROV").is_ok() {
        let a = pts[0];
        let b = *pts.last().unwrap();
        eprintln!(
            "[prov:single] ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}) n={}",
            a.x, a.y, a.z, b.x, b.y, b.z, pts.len()
        );
    }
    let v0 = get_or_create_vertex(brep, pts[0], vertex_map, tolerance);
    let v1 = get_or_create_vertex(brep, *pts.last().unwrap(), vertex_map, tolerance);
    get_or_create_edge(
        brep,
        EdgeGeom::Polyline(pts.to_vec()),
        v0,
        v1,
        edge_map,
        tolerance,
    )
}
/// Split a closed loop that spans an elbow↔leg junction into its elbow-part and
/// adjacent-leg-part, with the two crossing points snapped onto the junction
/// minor circle.  Returns `None` for loops that are not a clean single-junction
/// span (≠ 2 elbow/leg transitions).
fn split_spanning_loop(
    pts: &[Point3],
    torus_faces: &[TorusFace],
    cyl_faces: &[CylFace],
) -> Option<JunctionSplit> {
    let n = pts.len();
    if n < 6 {
        return None;
    }
    let on = |tf: &TorusFace, p: Point3| loop_on_torus_face(tf, std::slice::from_ref(&p), 1e-2);
    // Elbow with the most (but not all) points on it.
    let mut ti = 0usize;
    let mut best = 0usize;
    for (i, tf) in torus_faces.iter().enumerate() {
        let c = pts.iter().filter(|&&p| on(tf, p)).count();
        if c > best && c < n {
            best = c;
            ti = i;
        }
    }
    if best == 0 {
        return None;
    }
    let tf = &torus_faces[ti];
    let on_e: Vec<bool> = pts.iter().map(|&p| on(tf, p)).collect();
    let trans: Vec<usize> = (0..n).filter(|&i| on_e[i] != on_e[(i + 1) % n]).collect();
    if trans.len() != 2 {
        return None;
    }
    let (t0, t1) = (trans[0], trans[1]);
    let collect_run = |from: usize, to: usize| -> Vec<Point3> {
        let mut v = Vec::new();
        let mut i = (from + 1) % n;
        loop {
            v.push(pts[i]);
            if i == to {
                break;
            }
            i = (i + 1) % n;
        }
        v
    };
    let seg1 = collect_run(t0, t1);
    let seg2 = collect_run(t1, t0);
    let (e_arc, a_arc) = if on_e[(t0 + 1) % n] {
        (seg1, seg2)
    } else {
        (seg2, seg1)
    };
    if e_arc.is_empty() || a_arc.is_empty() {
        return None;
    }
    // Junction circle: whichever elbow end the e-arc endpoints sit near.
    let theta_of = |p: Point3| {
        let l = tf.surf.frame.to_local_point(p);
        l.y.atan2(l.x)
    };
    let ths = theta_of(*e_arc.first().unwrap());
    let the = theta_of(*e_arc.last().unwrap());
    let d_lo = (ths - tf.theta_lo).abs().min((the - tf.theta_lo).abs());
    let d_hi = (ths - tf.theta_hi).abs().min((the - tf.theta_hi).abs());
    let j_is_start = d_lo <= d_hi;
    let j = if j_is_start {
        tf.start_circle
    } else {
        tf.end_circle
    };
    let p1 = snap_to_circle(&j, *e_arc.first().unwrap());
    let p2 = snap_to_circle(&j, *e_arc.last().unwrap());
    if (p1 - p2).length() < 1e-6 {
        return None;
    }
    let aj = cyl_faces.iter().position(|cf| {
        boundary_circle_matches(&cf.start, &j) || boundary_circle_matches(&cf.end, &j)
    });
    let mut e_e = vec![p1];
    e_e.extend(e_arc.iter().copied());
    e_e.push(p2);
    let mut e_a = vec![p2];
    e_a.extend(a_arc.iter().copied());
    e_a.push(p1);
    Some(JunctionSplit {
        ti,
        aj,
        j_is_start,
        j,
        e_e,
        e_a,
    })
}
fn build_hole_edges(
    brep: &mut BRep,
    loop_pts: &[Point3],
    vertex_map: &mut VertexMap,
    edge_map: &mut EdgeMap,
    tolerance: f64,
    site: &'static str,
) -> (EdgeId, EdgeId) {
    if std::env::var("CADCORE_PROV").is_ok() {
        let a = loop_pts[0];
        let m = loop_pts[loop_pts.len() / 2];
        eprintln!(
            "[prov:hole2:{site}] seam0=({:.4},{:.4},{:.4}) seam1=({:.4},{:.4},{:.4}) n={}",
            a.x, a.y, a.z, m.x, m.y, m.z, loop_pts.len()
        );
    }
    // Drop a duplicated closing point if present.
    let mut pts: Vec<Point3> = loop_pts.to_vec();
    if pts.len() >= 2 && (pts[0] - pts[pts.len() - 1]).length() < 1e-9 {
        pts.pop();
    }
    // (Seam-vertex rotation was tried here and REVERTED: moving the split
    // points off the loop extremes worsened sewing — the fragments come from
    // mixed pcurve/no-pcurve wires, not vertex placement.)
    let n = pts.len();
    let mid = n / 2;

    let seam0 = pts[0];
    let seam1 = pts[mid];
    let v0 = get_or_create_vertex(brep, seam0, vertex_map, tolerance);
    let v1 = get_or_create_vertex(brep, seam1, vertex_map, tolerance);

    // e0: seam0 → seam1 (points 0..=mid)
    let mut p0: Vec<Point3> = pts[0..=mid].to_vec();
    if (p0[0] - seam0).length() > 1e-12 {
        p0.insert(0, seam0);
    }
    let e0 = get_or_create_edge(
        brep,
        EdgeGeom::Polyline(p0),
        v0,
        v1,
        edge_map,
        tolerance,
    );

    // e1: seam1 → seam0 (points mid..n, then wrap to seam0)
    let mut p1: Vec<Point3> = pts[mid..n].to_vec();
    p1.push(seam0);
    let e1 = get_or_create_edge(
        brep,
        EdgeGeom::Polyline(p1),
        v1,
        v0,
        edge_map,
        tolerance,
    );

    (e0, e1)
}
fn coedge_vertices(brep: &BRep, ce_id: CoEdgeId) -> Option<(VertexId, VertexId)> {
    let ce = brep.coedges.get(ce_id)?;
    let edge = brep.edges.get(ce.edge)?;
    Some(match ce.sense {
        CoEdgeSense::Same => (edge.v_start, edge.v_end),
        CoEdgeSense::Opposite => (edge.v_end, edge.v_start),
    })
}

fn coedge_points(brep: &BRep, ce_id: CoEdgeId) -> Option<(VertexId, VertexId, Point3, Point3)> {
    let (from, to) = coedge_vertices(brep, ce_id)?;
    let p0 = brep.vertices.get(from)?.point;
    let p1 = brep.vertices.get(to)?.point;
    Some((from, to, p0, p1))
}

#[derive(Clone, Copy)]
struct BoundaryLoc {
    edge_idx: usize,
    alpha: f64,
    vertex: VertexId,
}

#[derive(Clone, Copy)]
enum BoundaryUse {
    Existing(CoEdgeId),
    SplitLine(EdgeId),
}

fn point_on_segment_alpha(p0: Point3, p1: Point3, p: Point3, tol: f64) -> Option<f64> {
    let d = p1 - p0;
    let len2 = d.length_sq();
    if len2 <= tol * tol {
        return None;
    }
    let alpha = d.dot(p - p0) / len2;
    if alpha < -tol || alpha > 1.0 + tol {
        return None;
    }
    let q = p0 + d * alpha.clamp(0.0, 1.0);
    if (q - p).length() <= tol {
        Some(alpha.clamp(0.0, 1.0))
    } else {
        None
    }
}

fn locate_boundary_location(
    brep: &BRep,
    ordered: &[CoEdgeId],
    vertex: VertexId,
    tol: f64,
) -> Option<BoundaryLoc> {
    let p = brep.vertices.get(vertex)?.point;
    for (idx, &ce) in ordered.iter().enumerate() {
        let (from, _, p0, _) = coedge_points(brep, ce)?;
        if from == vertex || (p0 - p).length() <= tol {
            return Some(BoundaryLoc {
                edge_idx: idx,
                alpha: 0.0,
                vertex,
            });
        }
    }
    for (idx, &ce) in ordered.iter().enumerate() {
        let edge = brep.edges.get(brep.coedges.get(ce)?.edge)?;
        if !matches!(edge.geom, EdgeGeom::Line(_)) {
            continue;
        }
        let (_, _, p0, p1) = coedge_points(brep, ce)?;
        if let Some(alpha) = point_on_segment_alpha(p0, p1, p, tol) {
            if alpha <= tol || alpha >= 1.0 - tol {
                continue;
            }
            return Some(BoundaryLoc {
                edge_idx: idx,
                alpha,
                vertex,
            });
        }
    }
    None
}

fn add_split_line_edge(brep: &mut BRep, from: VertexId, to: VertexId) -> Option<EdgeId> {
    let p0 = brep.vertices.get(from)?.point;
    let p1 = brep.vertices.get(to)?.point;
    let dir = UnitVec3::try_from_vec(p1 - p0)?;
    Some(brep.add_edge(Edge {
        geom: EdgeGeom::Line(cadcore_geom::Line3::new(p0, dir)),
        v_start: from,
        v_end: to,
        t_start: 0.0,
        t_end: (p1 - p0).length(),
        partner: None,
    }))
}

fn push_boundary_piece(
    brep: &mut BRep,
    out: &mut Vec<BoundaryUse>,
    ordered: &[CoEdgeId],
    edge_idx: usize,
    alpha0: f64,
    from: VertexId,
    alpha1: f64,
    to: VertexId,
) -> Option<()> {
    if from == to {
        return Some(());
    }
    if alpha0 <= 1e-9 && alpha1 >= 1.0 - 1e-9 {
        out.push(BoundaryUse::Existing(ordered[edge_idx]));
    } else {
        let edge = add_split_line_edge(brep, from, to)?;
        out.push(BoundaryUse::SplitLine(edge));
    }
    Some(())
}

fn boundary_path_with_splits(
    brep: &mut BRep,
    ordered: &[CoEdgeId],
    start: BoundaryLoc,
    end: BoundaryLoc,
) -> Option<Vec<BoundaryUse>> {
    let n = ordered.len();
    let mut out = Vec::new();
    if start.edge_idx == end.edge_idx && start.alpha < end.alpha {
        push_boundary_piece(
            brep,
            &mut out,
            ordered,
            start.edge_idx,
            start.alpha,
            start.vertex,
            end.alpha,
            end.vertex,
        )?;
        return Some(out);
    }

    let (_, start_edge_to, _, _) = coedge_points(brep, ordered[start.edge_idx])?;
    push_boundary_piece(
        brep,
        &mut out,
        ordered,
        start.edge_idx,
        start.alpha,
        start.vertex,
        1.0,
        start_edge_to,
    )?;

    let mut idx = (start.edge_idx + 1) % n;
    while idx != end.edge_idx {
        out.push(BoundaryUse::Existing(ordered[idx]));
        idx = (idx + 1) % n;
    }

    let (end_edge_from, _, _, _) = coedge_points(brep, ordered[end.edge_idx])?;
    push_boundary_piece(
        brep,
        &mut out,
        ordered,
        end.edge_idx,
        0.0,
        end_edge_from,
        end.alpha,
        end.vertex,
    )?;
    Some(out)
}

fn add_split_region_loop(
    brep: &mut BRep,
    boundary_path: &[BoundaryUse],
    open_edge: EdgeId,
    open_sense: CoEdgeSense,
) -> Option<(LoopId, CoEdgeId)> {
    if boundary_path.is_empty() {
        return None;
    }
    let mut coedges = Vec::with_capacity(boundary_path.len() + 1);
    for &piece in boundary_path {
        let (edge, sense) = match piece {
            BoundaryUse::Existing(old_ce) => {
                let ce = brep.coedges.get(old_ce)?;
                (ce.edge, ce.sense)
            }
            BoundaryUse::SplitLine(edge) => (edge, CoEdgeSense::Same),
        };
        let new_ce = brep.add_coedge(CoEdge {
            edge,
            sense,
            next: CoEdgeId::default(),
            prev: CoEdgeId::default(),
            loop_id: LoopId::default(),
        });
        coedges.push(new_ce);
    }
    let open_coedge = brep.add_coedge(CoEdge {
        edge: open_edge,
        sense: open_sense,
        next: CoEdgeId::default(),
        prev: CoEdgeId::default(),
        loop_id: LoopId::default(),
    });
    coedges.push(open_coedge);
    brep.patch_coedge_links(&coedges);
    let loop_id = brep.add_loop(Loop {
        start: coedges[0],
        face: FaceId::default(),
    });
    for &ce in &coedges {
        if let Some(c) = brep.coedges.get_mut(ce) {
            c.loop_id = loop_id;
        }
    }
    Some((loop_id, open_coedge))
}

/// Split one trimmed face by an open cutter edge whose endpoints are on the
/// face's outer loop, then replace the original face with the selected region.
/// Endpoints may be existing boundary vertices or interior points on linear
/// boundary edges; the latter are promoted to real split boundary segments.
pub fn split_face_by_open_edge_at_existing_vertices(
    brep: &mut BRep,
    face_id: FaceId,
    open_edge: EdgeId,
    keep: OpenBoundaryKeep,
) -> Option<OpenBoundarySplit> {
    let face = brep.faces.get(face_id)?.clone();
    if !face.inner_loops.is_empty() {
        return None;
    }
    let open = brep.edges.get(open_edge)?;
    if open.v_start == open.v_end {
        return None;
    }
    let ordered = brep.loop_coedges(face.outer_loop)?;
    let n = ordered.len();
    if n < 2 {
        return None;
    }

    for i in 0..n {
        let (_, to) = coedge_vertices(brep, ordered[i])?;
        let (next_from, _) = coedge_vertices(brep, ordered[(i + 1) % n])?;
        if to != next_from {
            return None;
        }
    }

    let tol = 1e-7;
    let start = locate_boundary_location(brep, &ordered, open.v_start, tol)?;
    let end = locate_boundary_location(brep, &ordered, open.v_end, tol)?;
    if start.edge_idx == end.edge_idx && (start.alpha - end.alpha).abs() <= tol {
        return None;
    }

    let path_start_to_end = boundary_path_with_splits(brep, &ordered, start, end)?;
    let path_end_to_start = boundary_path_with_splits(brep, &ordered, end, start)?;
    let start_to_end_len = path_start_to_end.len();
    let end_to_start_len = path_end_to_start.len();
    let (boundary_path, dropped_boundary_coedges, open_sense) = match keep {
        OpenBoundaryKeep::StartToEnd => {
            (path_start_to_end, end_to_start_len, CoEdgeSense::Opposite)
        }
        OpenBoundaryKeep::EndToStart => (path_end_to_start, start_to_end_len, CoEdgeSense::Same),
    };
    let boundary_coedges = boundary_path.len();
    let (loop_id, open_coedge) =
        add_split_region_loop(brep, &boundary_path, open_edge, open_sense)?;

    let new_face_id = brep.add_face(Face {
        geom: face.geom,
        normal: face.normal,
        outer_loop: loop_id,
        inner_loops: vec![],
        shell: face.shell,
        extent: FaceExtent::Trimmed,
    });
    if let Some(lp) = brep.loops.get_mut(loop_id) {
        lp.face = new_face_id;
    }
    for ce_id in brep.loop_coedges(loop_id)? {
        if let Some(ce) = brep.coedges.get_mut(ce_id) {
            ce.loop_id = loop_id;
        }
    }
    if let Some(shell) = brep.shells.get_mut(face.shell) {
        if let Some(slot) = shell.faces.iter_mut().find(|fid| **fid == face_id) {
            *slot = new_face_id;
        } else {
            shell.faces.push(new_face_id);
        }
    }
    brep.faces.remove(face_id);

    Some(OpenBoundarySplit {
        face_id: new_face_id,
        loop_id,
        open_coedge,
        boundary_coedges,
        dropped_boundary_coedges,
    })
}


fn add_closed_analytic_loop(
    brep: &mut BRep,
    geom: EdgeGeom,
    seam: Point3,
    vertex_map: &mut VertexMap,
    edge_map: &mut EdgeMap,
    tolerance: f64,
) -> LoopId {
    let v = get_or_create_vertex(brep, seam, vertex_map, tolerance);
    let e = get_or_create_edge(brep, geom, v, v, edge_map, tolerance);
    let ce = brep.add_coedge(CoEdge {
        edge: e,
        sense: CoEdgeSense::Same,
        next: CoEdgeId::default(),
        prev: CoEdgeId::default(),
        loop_id: LoopId::default(),
    });
    brep.patch_coedge_links(&[ce]);
    brep.add_loop(Loop {
        start: ce,
        face: FaceId::default(),
    })
}

/// Loop for an end boundary (circle or miter ellipse) of a cylinder segment.
fn boundary_loop(
    brep: &mut BRep,
    b: &FaceBoundary,
    vertex_map: &mut VertexMap,
    edge_map: &mut EdgeMap,
    tolerance: f64,
) -> LoopId {
    match b {
        FaceBoundary::Circle(c) => {
            let seam = c.point_at(0.0);
            add_closed_analytic_loop(brep, EdgeGeom::Circle(*c), seam, vertex_map, edge_map, tolerance)
        }
        FaceBoundary::Ellipse(e) => {
            let seam = e.point_at(0.0);
            add_closed_analytic_loop(brep, EdgeGeom::Ellipse(*e), seam, vertex_map, edge_map, tolerance)
        }
    }
}

/// OUTWARD surface normal of an analytic face at a surface point `p`.
fn face_outward_normal(geom: &FaceGeom, p: Point3) -> cadcore_math::Vec3 {
    match geom {
        FaceGeom::Cylinder(c) => {
            let d = p - c.frame.origin;
            let axial = c.frame.z.dot_vec(d);
            (d - c.frame.z.as_vec() * axial).normalize()
        }
        FaceGeom::Torus(t) => {
            let d = p - t.frame.origin;
            let axial = t.frame.z.dot_vec(d);
            let e = (d - t.frame.z.as_vec() * axial).normalize();
            let centerline = t.frame.origin + e * t.major_radius;
            (p - centerline).normalize()
        }
        _ => cadcore_math::Vec3::new(0.0, 0.0, 1.0),
    }
}

/// Does travelling from `from` in rotational direction `dir` (+1 CCW / −1 CW)
/// reach `mid` before `to`?
fn arc_contains(from: f64, to: f64, mid: f64, dir: f64) -> bool {
    let two_pi = 2.0 * PI;
    if dir > 0.0 {
        (mid - from).rem_euclid(two_pi) <= (to - from).rem_euclid(two_pi)
    } else {
        (from - mid).rem_euclid(two_pi) <= (from - to).rem_euclid(two_pi)
    }
}

/// Build the notched boundary wire of a face whose end circle is bitten by one
/// or more junction-split trim curves.
///
/// All directions and co-edge senses are **derived from geometry**:
/// 1. each notch is oriented so its angular interval (entry→exit along the
///    walk direction) covers the bitten section (the notch curve's midpoint);
/// 2. the clear gaps between consecutive notches are walked in the same
///    rotational direction, split at FIXED absolute angles (multiples of 90° on
///    the canonical junction circle) so the partner face — which walks the same
///    circle in the OPPOSITE direction — produces the identical sub-arc edges
///    and shares them with opposite senses;
/// 3. the walk direction (CW/CCW on the canonical circle) is chosen by the
///    material-on-the-left rule: walking direction `t` must satisfy
///    `(n_outward × t) · (interior − p) > 0` — the elbow and the adjacent leg
///    have their interiors on opposite sides of the junction, so they derive
///    opposite walks.
fn boundary_multi_notch_loop(
    brep: &mut BRep,
    geom: &FaceGeom,
    b: &FaceBoundary,
    interior_pt: Point3,
    notches: &[OpenNotchRef],
    vertex_map: &mut VertexMap,
    edge_map: &mut EdgeMap,
    tolerance: f64,
) -> Option<LoopId> {
    // The CANONICAL junction circle (shared by both faces) so all gap-arc
    // geometry is identical on the elbow and the adjacent leg.
    let circle = match notches.iter().find_map(|n| n.gap_circle) {
        Some(c) => c,
        None => match b {
            FaceBoundary::Circle(c) => *c,
            FaceBoundary::Ellipse(_) => return None,
        },
    };
    if notches.is_empty() {
        return None;
    }

    struct RawNotch {
        edges: Vec<EdgeId>, // chain order (stored directions chain head→tail)
        va: VertexId,
        vb: VertexId,
        aa: f64,
        ab: f64,
        mid_a: f64,
        pts: Vec<Point3>, // stored order = va → vb
    }
    let mut raw = Vec::with_capacity(notches.len());
    for notch in notches {
        let first = brep.edges.get(*notch.edges.first()?)?.clone();
        let last = brep.edges.get(*notch.edges.last()?)?.clone();
        let va = first.v_start;
        let vb = last.v_end;
        let pa = brep.vertices.get(va)?.point;
        let pb = brep.vertices.get(vb)?.point;
        let mut pts: Vec<Point3> = Vec::new();
        for &e in &notch.edges {
            let edge = brep.edges.get(e)?.clone();
            match &edge.geom {
                EdgeGeom::Polyline(p) => pts.extend(p.iter().copied()),
                _ => {
                    pts.push(brep.vertices.get(edge.v_start)?.point);
                    pts.push(brep.vertices.get(edge.v_end)?.point);
                }
            }
        }
        let mid_p = pts[pts.len() / 2];
        raw.push(RawNotch {
            edges: notch.edges.clone(),
            va,
            vb,
            aa: circle_angle(&circle, pa),
            ab: circle_angle(&circle, pb),
            mid_a: circle_angle(&circle, mid_p),
            pts,
        });
    }

    // Each notch edge is SHARED with the crossing tube's hole loop, whose
    // traversal (CW in the crossing tube's parameters) is already fixed —
    // so this wire MUST traverse every notch opposite to the hole.  That
    // direction is carried in `notch.reversed` (set in Pass 2c as the
    // complement of the hole orientation).  Only the GAP PATH between notches
    // is derived geometrically: from a notch's exit, the gap leaves through
    // the CLEAR side of the junction circle (away from that notch's bitten
    // interval) and runs to the nearest notch entry.
    struct Oriented {
        ri: usize,
        entry_v: VertexId,
        exit_v: VertexId,
        entry_a: f64,
        exit_a: f64,
        sense: CoEdgeSense,
        /// rotational direction (+1 CCW / −1 CW) that LEAVES the bite at exit
        clear_dir: f64,
    }
    let mut oriented = Vec::with_capacity(raw.len());
    for (i, (n, notch)) in raw.iter().zip(notches.iter()).enumerate() {
        let (entry_v, exit_v, entry_a, exit_a, sense) = if notch.reversed {
            (n.vb, n.va, n.ab, n.aa, CoEdgeSense::Opposite)
        } else {
            (n.va, n.vb, n.aa, n.ab, CoEdgeSense::Same)
        };
        // Clear side from the SHARED Pass 2c flag (crossing-tube sampled):
        // bite ⊂ CCW(v_start→v_end) iff `bite_ccw_from_start`.  Moving CCW
        // from EXIT walks the complement of (entry→exit):
        //  - not reversed (exit = v_end): CCW from v_end = complement of
        //    CCW(v_start→v_end) → bite there iff !flag → clear CCW iff flag;
        //  - reversed (exit = v_start): CCW from v_start = CCW(v_start→v_end)
        //    → bite there iff flag → clear CCW iff !flag.
        let clear_dir = if notch.reversed == notch.bite_ccw_from_start {
            -1.0
        } else {
            1.0
        };
        let _ = n.mid_a;
        oriented.push(Oriented { ri: i, entry_v, exit_v, entry_a, exit_a, sense, clear_dir });
    }

    // Chain notches: from each exit, walk the clear side to the NEAREST other
    // endpoint, which must be some notch's entry; repeat until the wire closes
    // over all notches.  Anything inconsistent → None (caller falls back).
    let two_pi = 2.0 * PI;
    let eps_ang = 1e-9;
    let n = oriented.len();
    let mut order: Vec<usize> = Vec::with_capacity(n);
    let mut visited = vec![false; n];
    let mut cur = 0usize;
    for _ in 0..n {
        visited[cur] = true;
        order.push(cur);
        let d = oriented[cur].clear_dir;
        let from = oriented[cur].exit_a;
        // nearest endpoint (entry of any unvisited notch, or the start notch's
        // entry to close) in direction d
        let mut best: Option<(usize, f64)> = None;
        for (j, o) in oriented.iter().enumerate() {
            let candidates = [(j, o.entry_a)];
            for &(jj, a) in &candidates {
                let delta = (d * (a - from)).rem_euclid(two_pi);
                let delta = if delta < eps_ang { two_pi } else { delta };
                let eligible = if order.len() == n { jj == order[0] } else { !visited[jj] };
                if eligible && best.map_or(true, |(_, bd)| delta < bd) {
                    best = Some((jj, delta));
                }
            }
        }
        let Some((next, _)) = best else {
            if std::env::var("CADCORE_DUMP_CONT").is_ok() {
                eprintln!("[notch-fail] no next entry (n={n}, walked={})", order.len());
            }
            return None;
        };
        if order.len() == n {
            if next != order[0] {
                if std::env::var("CADCORE_DUMP_CONT").is_ok() {
                    eprintln!("[notch-fail] chain does not close (n={n})");
                }
                return None;
            }
            break;
        }
        cur = next;
    }
    if order.len() != n {
        if std::env::var("CADCORE_DUMP_CONT").is_ok() {
            eprintln!("[notch-fail] visited {} of {n}", order.len());
        }
        return None;
    }

    let mut coedges = Vec::with_capacity(n * 3);
    for k in 0..n {
        let o = &oriented[order[k]];
        let nx = &oriented[order[(k + 1) % n]];
        // Emit the notch RUN: stored chain order when sense=Same, reversed
        // order with Opposite senses otherwise.
        match o.sense {
            CoEdgeSense::Same => {
                for &e in &raw[o.ri].edges {
                    coedges.push(brep.add_coedge(CoEdge {
                        edge: e,
                        sense: CoEdgeSense::Same,
                        next: CoEdgeId::default(),
                        prev: CoEdgeId::default(),
                        loop_id: LoopId::default(),
                    }));
                }
            }
            CoEdgeSense::Opposite => {
                for &e in raw[o.ri].edges.iter().rev() {
                    coedges.push(brep.add_coedge(CoEdge {
                        edge: e,
                        sense: CoEdgeSense::Opposite,
                        next: CoEdgeId::default(),
                        prev: CoEdgeId::default(),
                        loop_id: LoopId::default(),
                    }));
                }
            }
        }
        add_gap_arc_directed(
            brep,
            &circle,
            o.exit_v,
            o.exit_a,
            nx.entry_v,
            nx.entry_a,
            o.clear_dir,
            &mut coedges,
            vertex_map,
            edge_map,
            tolerance,
        );
    }
    let _ = interior_pt;

    brep.patch_coedge_links(&coedges);
    let loop_id = brep.add_loop(Loop {
        start: coedges[0],
        face: FaceId::default(),
    });
    for &ce in &coedges {
        if let Some(c) = brep.coedges.get_mut(ce) {
            c.loop_id = loop_id;
        }
    }

    Some(loop_id)
}

/// Interior waypoint angles for a gap arc from `af` to `at` along `dir`:
/// the multiples of 90° strictly inside the swept interval, in traversal
/// order, plus the end angle.  Absolute (direction-independent) split points
/// guarantee the two faces walking the same gap in opposite directions create
/// the IDENTICAL sub-arc edges.
fn gap_waypoints(af: f64, at: f64, dir: f64) -> Vec<f64> {
    let two_pi = 2.0 * PI;
    let step = PI / 2.0;
    let eps = 1e-9;
    let sweep = if dir > 0.0 {
        (at - af).rem_euclid(two_pi)
    } else {
        -((af - at).rem_euclid(two_pi))
    };
    let mut out = Vec::new();
    if sweep.abs() < eps {
        return out;
    }
    if dir > 0.0 {
        let mut k = (af / step).floor() * step + step;
        while k < af + sweep - eps {
            if k > af + eps {
                out.push(k);
            }
            k += step;
        }
    } else {
        let mut k = (af / step).ceil() * step - step;
        while k > af + sweep + eps {
            if k < af - eps {
                out.push(k);
            }
            k -= step;
        }
    }
    out.push(af + sweep);
    out
}

/// Append the clear-gap circle-arc co-edges from `v_from` (angle `af`) to
/// `v_to` (angle `at`) walking direction `dir` on `circle`.  Sub-arc senses
/// are DERIVED from each shared edge's stored direction (the partner face may
/// have created it first, traversing the other way).
fn add_gap_arc_directed(
    brep: &mut BRep,
    circle: &Circle3,
    v_from: VertexId,
    af: f64,
    v_to: VertexId,
    at: f64,
    dir: f64,
    coedges: &mut Vec<CoEdgeId>,
    vertex_map: &mut VertexMap,
    edge_map: &mut EdgeMap,
    tolerance: f64,
) {
    let angles = gap_waypoints(af, at, dir);
    if angles.is_empty() {
        return; // zero-length gap (shared vertex)
    }
    let last = angles.len() - 1;
    let mut prev_v = v_from;
    for (k, &a) in angles.iter().enumerate() {
        let next_v = if k == last {
            v_to
        } else {
            get_or_create_vertex(brep, circle.point_at(a), vertex_map, tolerance)
        };
        if next_v == prev_v {
            continue;
        }
        let e = get_or_create_edge(
            brep,
            EdgeGeom::Circle(*circle),
            prev_v,
            next_v,
            edge_map,
            tolerance,
        );
        let sense = match brep.edges.get(e) {
            Some(edge) if edge.v_start == prev_v => CoEdgeSense::Same,
            _ => CoEdgeSense::Opposite,
        };
        coedges.push(brep.add_coedge(CoEdge {
            edge: e,
            sense,
            next: CoEdgeId::default(),
            prev: CoEdgeId::default(),
            loop_id: LoopId::default(),
        }));
        prev_v = next_v;
    }
}
/// Hole loop referencing the two shared edges, in the orientation for one of the
/// two faces (`reversed` flips traversal so the partner gets opposite senses).
fn hole_loop(brep: &mut BRep, h: &HoleRef) -> LoopId {
    let (a, sa, b, sb) = if !h.reversed {
        // seam0 → seam1 → seam0
        (h.e0, CoEdgeSense::Same, h.e1, CoEdgeSense::Same)
    } else {
        // opposite traversal: seam0 →(e1 rev)→ seam1 →(e0 rev)→ seam0
        (h.e1, CoEdgeSense::Opposite, h.e0, CoEdgeSense::Opposite)
    };
    let c0 = brep.add_coedge(CoEdge {
        edge: a,
        sense: sa,
        next: CoEdgeId::default(),
        prev: CoEdgeId::default(),
        loop_id: LoopId::default(),
    });
    let c1 = brep.add_coedge(CoEdge {
        edge: b,
        sense: sb,
        next: CoEdgeId::default(),
        prev: CoEdgeId::default(),
        loop_id: LoopId::default(),
    });
    brep.patch_coedge_links(&[c0, c1]);
    brep.add_loop(Loop {
        start: c0,
        face: FaceId::default(),
    })
}

/// Build a trimmed cylindrical face: original end boundaries + window holes.
/// Build a trimmed swept face for ANY carrier (cylinder leg OR torus elbow):
/// its two end boundaries (shared with adjoining segments/caps) + window holes.
/// The carrier surface stays analytic (`CYLINDRICAL_SURFACE` / `TOROIDAL_SURFACE`);
/// only the trim loops are explicit topology.
fn build_trimmed_swept(
    brep: &mut BRep,
    geom: FaceGeom,
    start: &FaceBoundary,
    end: &FaceBoundary,
    interior_pt: Point3,
    holes: &[HoleRef],
    open_notches: &[OpenNotchRef],
    multi_holes: &[MultiHoleRef],
    shell_id: ShellId,
    vertex_map: &mut VertexMap,
    edge_map: &mut EdgeMap,
    tolerance: f64,
) -> FaceId {
    // Trimmed cylinder: rotate the parameterization so the θ=0 seam threads
    // the largest gap between trim features (all wires stay on one pcurve
    // branch — see `seam_rotated_cyl`).
    let geom = if let FaceGeom::Cylinder(c) = &geom {
        let mut avoid: Vec<Point3> = Vec::new();
        let mut collect = |e: EdgeId| {
            if let Some(edge) = brep.edges.get(e) {
                if let EdgeGeom::Polyline(pts) = &edge.geom {
                    avoid.extend(pts.iter().copied());
                }
            }
        };
        for h in holes {
            collect(h.e0);
            collect(h.e1);
        }
        for mh in multi_holes {
            for &e in &mh.edges {
                collect(e);
            }
        }
        for n in open_notches {
            for &e in &n.edges {
                collect(e);
            }
        }
        FaceGeom::Cylinder(seam_rotated_cyl(c, &avoid))
    } else {
        geom
    };

    // Notches from the continuous-SSI junction split (Pass 2c).
    let start_notches: Vec<OpenNotchRef> = open_notches
        .iter()
        .filter(|n| n.side == BoundarySide::Start)
        .cloned()
        .collect();
    let end_notches: Vec<OpenNotchRef> = open_notches
        .iter()
        .filter(|n| n.side == BoundarySide::End)
        .cloned()
        .collect();
    let outer = if start_notches.is_empty() {
        boundary_loop(brep, start, vertex_map, edge_map, tolerance)
    } else {
        boundary_multi_notch_loop(
            brep,
            &geom,
            start,
            interior_pt,
            &start_notches,
            vertex_map,
            edge_map,
            tolerance,
        )
        .unwrap_or_else(|| boundary_loop(brep, start, vertex_map, edge_map, tolerance))
    };
    let first_inner = if end_notches.is_empty() {
        boundary_loop(brep, end, vertex_map, edge_map, tolerance)
    } else {
        boundary_multi_notch_loop(
            brep,
            &geom,
            end,
            interior_pt,
            &end_notches,
            vertex_map,
            edge_map,
            tolerance,
        )
        .unwrap_or_else(|| boundary_loop(brep, end, vertex_map, edge_map, tolerance))
    };
    let mut inner = vec![first_inner];
    for h in holes {
        inner.push(hole_loop(brep, h));
    }
    for mh in multi_holes {
        if let Some(l) = multi_hole_loop(brep, mh) {
            inner.push(l);
        }
    }
    let face_id = brep.add_face(Face {
        geom,
        normal: FaceNormal::Same,
        outer_loop: outer,
        inner_loops: inner.clone(),
        shell: shell_id,
        extent: FaceExtent::Trimmed,
    });
    // Back-link loops to the face (keep topology self-consistent).
    for lp in std::iter::once(outer).chain(inner) {
        if let Some(l) = brep.loops.get_mut(lp) {
            l.face = face_id;
        }
    }
    face_id
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spanning_hole_reversed_makes_crossing_cylinder_loop_ccw() {
        let cyl = CylSurf::new(Point3::ORIGIN, UnitVec3::Z, 1.0);
        let p00 = cyl.point_at(0.0, 0.0);
        let p10 = cyl.point_at(0.4, 0.0);
        let p11 = cyl.point_at(0.4, 1.0);
        let p01 = cyl.point_at(0.0, 1.0);

        let ccw_edge_ee = [p00, p10, p11];
        let ccw_edge_ea = [p11, p01, p00];
        assert!(!hole_ref_reversed_for_cylinder_ccw(
            &cyl,
            &ccw_edge_ee,
            &ccw_edge_ea
        ));

        let cw_edge_ee = [p00, p01, p11];
        let cw_edge_ea = [p11, p10, p00];
        assert!(hole_ref_reversed_for_cylinder_ccw(
            &cyl,
            &cw_edge_ee,
            &cw_edge_ea
        ));
    }
}
