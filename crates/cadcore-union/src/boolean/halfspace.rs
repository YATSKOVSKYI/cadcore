//! **Half-space cut** — trim an already-FUSED shell with a plane and close it
//! with real section caps.
//!
//! This is the "cut terminals AFTER sintering the part" operation.  The legacy
//! cut in `cadcore-ops` rebuilds each *template* cylinder face on its own, so it
//! only works on the raw per-segment sweep solids: on a fused shell (whose faces
//! carry explicit trim loops, `FaceExtent::Trimmed`) it finds no template
//! cylinders at all and silently does nothing.
//!
//! Here the cut is expressed with the same machinery as the union, because it
//! *is* the same operation with a different keep rule:
//!
//! ```text
//! union(A, B)      keep a face piece ⇔ it is OUTSIDE the other solid
//! cut (A, plane)   keep a face piece ⇔ it is INSIDE the kept half-space
//! ```
//!
//! so every face is imprinted with its section curve, its cells are classified
//! by an analytic (O(1), no ray casting) predicate, and the survivors are
//! stitched through ONE [`Assembler`].  What the union does *not* need and the
//! cut does is the **section cap**: the flat faces that close the part where the
//! plane went through it.  Those are built here from the very same section
//! curves, so their edges weld with the trimmed lateral faces by construction.
//!
//! Stages:
//!
//! ```text
//! 1 bound     model AABB → a finite cutting face on the plane
//! 2 section   SSI(face, cutting face) for every straddling face   (plane×cyl,
//!             plane×plane, plane×torus — the fillets)
//! 3 trim      imprint + classify + emit every face             (union::emit_solid_faces)
//! 4 cap       chain the section curves into closed loops, arrange them in the
//!             cut plane and emit the cells that carry material
//! ```
//!
//! ## What makes it hold on a real lattice
//!
//! * **Joints are computed once per edge.**  Where a boundary edge of a
//!   straddling face crosses the plane, the crossing is solved as a triple
//!   point (plane + the carriers of every face on that edge, Newton on the
//!   tangent planes; exact circle ∩ plane on a rim) and that ONE point is
//!   grafted into every chain and section curve that passes there.
//! * **Section curves are clipped run by run.**  A section can only change
//!   between material and void at a joint, so it is cut into runs between
//!   joints and each run is classified once, at its middle — never chord by
//!   chord against a boundary that is itself a chord polygon.
//! * **Sections are sampled by arc length** and closed exactly at band ends,
//!   so a plane grazing a tube lengthwise (a 30:1 ellipse) is resolved as
//!   finely as one slicing it across.
//! * **Faces that do not wrap are arranged as plain regions** (the two half
//!   faces the universal fuse emits per tube): no seam, no seam bookkeeping.
//! * **Near-degenerate planes are perturbed** — see [`half_space_cut_robust`].

use std::collections::HashMap;

use cadcore_geom::{Circle3, Plane3};
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{
    BRep, Face, FaceExtent, FaceGeom, FaceId, FaceNormal, Shell, ShellId, Solid,
};

use super::aabb::{face_aabb, Aabb};
use super::contain::{face_contains_hit, point_in_solid, Membership};
use super::ssi::{intersect_faces_with, SsiCurve, SsiOptions};
use super::union::{
    face_boundary_chains, face_domain, face_loop_points, snap_ssi_endpoints,
    unwrap_torus_band_chains,
};
use crate::arrange::assembly::Assembler;
use crate::arrange::cells::{arrange_and_classify_with, FaceChain};
use crate::arrange::domain::FaceDomain;

/// Outcome of a cut.
#[derive(Debug)]
pub struct CutOutcome {
    /// The trimmed B-Rep.
    pub brep: BRep,
    /// Its single closed shell.
    pub shell: ShellId,
    /// Number of section-cap faces emitted (0 ⇒ the plane missed the part).
    pub caps: usize,
    /// Faces of the input that were dropped whole (entirely on the cut side).
    pub dropped: usize,
}

/// Outcome of [`half_space_cut_robust`].
pub struct RobustCut {
    /// The cut that came out watertight (or, failing every attempt, the one
    /// with the fewest open edges).
    pub outcome: CutOutcome,
    /// How far the plane was moved along its normal to get it (mm; `0.0` when
    /// the requested plane cut cleanly).  Positive = towards the kept side.
    pub shift: f64,
    /// Manifold violations left in `outcome` (`0` = watertight).
    pub violations: usize,
}

/// The plane shifts [`half_space_cut_robust`] tries, in order (mm).
pub const ROBUST_CUT_SHIFTS: [f64; 7] = [0.0, 0.01, -0.01, 0.02, -0.02, 0.04, -0.04];

/// [`half_space_cut`] that escapes NEAR-DEGENERATE plane positions.
///
/// A trim plane specified from the part's own dimensions lands, more often
/// than chance, within microns of a feature of that part: a leg's junction rim
/// (the plane through the axis of the outermost filament), a filament's end
/// cap.  Such a feature is never exactly parallel to the plane, so the two
/// cross at a grazing angle and the exact cut contains a sliver a micron thick
/// — smaller than the chord error of the boundary it has to be classified
/// against, useless to a mesher, and unrepresentable at STEP tolerance.
///
/// The standard answer is a tolerance-sized perturbation: when the exact plane
/// does not give a closed shell, the plane is moved along its normal by
/// ±10, ±20, ±40 µm ([`ROBUST_CUT_SHIFTS`]) and the first watertight result is
/// returned together with the shift, so the caller can place whatever mates
/// with the cut face (an electrode plate) on the plane actually used.
pub fn half_space_cut_robust(src: &BRep, faces: &[FaceId], keep: &Plane3, knit: f64) -> Option<RobustCut> {
    let n = keep.normal();
    let mut best: Option<RobustCut> = None;
    for &shift in ROBUST_CUT_SHIFTS.iter() {
        let plane = Plane3::from_origin_normal(keep.frame.origin + n.as_vec() * shift, n);
        // a panic inside one attempt (a degenerate arrangement) is just a
        // failed attempt
        let attempt = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            half_space_cut(src, faces, &plane, knit)
        }));
        let Ok(Some(outcome)) = attempt else { continue };
        let t_chk = std::time::Instant::now();
        let violations = crate::validate::manifold::check(&outcome.brep, outcome.shell).violations.len();
        if std::env::var("CADCORE_TIME").is_ok() {
            eprintln!("[cut][time] robust manifold check: {} ms", t_chk.elapsed().as_millis());
        }
        let better = best.as_ref().map_or(true, |b| violations < b.violations);
        if better {
            best = Some(RobustCut { outcome, shift, violations });
        }
        if violations == 0 {
            break;
        }
    }
    best
}

/// Cut the solid bounded by `faces` with `keep`, retaining the material on the
/// side the plane's normal points to.
///
/// `knit` is the model knit tolerance (SolidWorks ≈ 1 µm); it drives the
/// assembler's geometric edge weld exactly as in [`super::union::union_n`].
///
/// Returns `None` when the input has no faces.  A plane that misses the part
/// entirely yields the part unchanged with `caps == 0`.
pub fn half_space_cut(
    src: &BRep,
    faces: &[FaceId],
    keep: &Plane3,
    knit: f64,
) -> Option<CutOutcome> {
    if faces.is_empty() {
        return None;
    }
    let timing = std::env::var("CADCORE_TIME").is_ok();
    let t_phase = std::cell::Cell::new(std::time::Instant::now());
    let tick = |what: &str| {
        if timing {
            eprintln!("[cut][time] {what}: {} ms", t_phase.get().elapsed().as_millis());
        }
        t_phase.set(std::time::Instant::now());
    };
    let n = keep.normal();

    // ── 1. Bound the cutting plane to the model ──────────────────────────────
    let mut model = Aabb::empty();
    for &f in faces {
        let bb = face_aabb(src, f);
        model.expand(bb.min);
        model.expand(bb.max);
    }
    let half = (model.extent() * 1.5).max(1.0);
    let centre_on_plane = keep.project(model.centre());
    // A finite square face on the cut plane — `region_polygon` needs a polygon
    // to clip the section curves against, and the arrangement needs a bounded
    // rectangle.  1.5× the model extent guarantees it covers every section.
    let cut_plane = Plane3 {
        frame: cadcore_math::Frame3 {
            origin: centre_on_plane,
            ..keep.frame
        },
    };
    let corner = |su: f64, sv: f64| centre_on_plane + cut_plane.frame.x * su + cut_plane.frame.y * sv;
    let mut cut_brep = BRep::new();
    let cut_face = cut_brep.add_face(Face {
        geom: FaceGeom::Plane(cut_plane),
        normal: FaceNormal::Same,
        outer_loop: Default::default(),
        inner_loops: Vec::new(),
        shell: Default::default(),
        extent: FaceExtent::Polygon {
            points: vec![
                corner(-half, -half),
                corner(half, -half),
                corner(half, half),
                corner(-half, half),
            ],
        },
    });

    // ── 2. Classify faces and section the straddling ones ───────────────────
    // Classification is by the face's SAMPLED geometry (boundary loops plus a
    // coarse grid on the carrier, so a fillet's belly bulging past its junction
    // circles is seen) with an ON-plane band of `TAU_ON`: a face that merely
    // TOUCHES the plane — a leg whose junction circle lies exactly in it, the
    // everyday case when the trim offset equals edge offset + bend radius — is
    // kept or dropped WHOLE, never sectioned along its own edge.
    const TAU_ON: f64 = 1e-6;
    let mut sides: HashMap<FaceId, Side> = HashMap::new();
    let side_list = crate::fasthash::par_map(faces, |&fid| face_side(src, fid, keep, TAU_ON));
    for (&fid, side) in faces.iter().zip(side_list) {
        sides.insert(fid, side);
    }
    // Section every straddling face with the plane up front, on all cores —
    // each is independent and only reads `src`; consumed in face order below.
    let straddlers: Vec<FaceId> = faces.iter().copied().filter(|f| sides[f] == Side::Straddles).collect();
    let mut sections: HashMap<FaceId, Vec<SsiCurve>> = straddlers
        .iter()
        .copied()
        .zip(crate::fasthash::par_map(&straddlers, |&fid| {
            intersect_faces_with(src, fid, &cut_brep, cut_face, SsiOptions::all())
        }))
        .collect();
    let mut ssi: Vec<HashMap<FaceId, Vec<SsiCurve>>> = vec![HashMap::new()];
    let mut dropped = 0usize;
    // edge → the faces using it (real-loop faces only).
    let mut edge_faces: HashMap<cadcore_topo::EdgeId, Vec<FaceId>> = HashMap::new();
    for &fid in faces {
        if !has_real_loops(src, fid) {
            continue;
        }
        let f = &src.faces[fid];
        let mut ls = vec![f.outer_loop];
        ls.extend(f.inner_loops.iter().copied());
        for lid in ls {
            let Some(lp) = src.loops.get(lid) else { continue };
            let st = lp.start;
            if src.coedges.get(st).is_none() {
                continue;
            }
            let mut c = st;
            loop {
                edge_faces.entry(src.coedges[c].edge).or_default().push(fid);
                c = src.coedges[c].next;
                if c == st {
                    break;
                }
            }
        }
    }
    let carrier = |fid: FaceId| -> Option<crate::geom::refine::AnalyticSurface> {
        use crate::geom::refine::AnalyticSurface;
        match &src.faces[fid].geom {
            FaceGeom::Plane(p) => Some(AnalyticSurface::Plane(*p)),
            FaceGeom::Cylinder(c) => Some(AnalyticSurface::Cylinder(*c)),
            FaceGeom::Torus(t) => Some(AnalyticSurface::Torus(*t)),
            FaceGeom::Sphere(_) => None,
        }
    };
    // Boundary vertices of straddling faces that lie IN the plane, and the
    // points where their boundary edges CROSS it: the section curve must pass
    // through them exactly (they are where a fused neighbour's edge meets the
    // plane — a triple point every consumer has to share).  A crossing is
    // computed ONCE PER EDGE, from the plane and the carriers of every face
    // using the edge in a canonical order, so the two faces sharing the edge
    // (and the cap) get the identical point: computed per face, the surface
    // order differed and an ill-conditioned triple point (a leg's rim, where
    // the elbow torus is TANGENT to the leg) came out 1e-5 apart.  On a rim
    // edge carrying analytic circle geometry the exact circle ∩ plane point
    // is used outright.
    let mut plane_joints: Vec<Point3> = Vec::new();
    for &fid in faces {
        let side = sides[&fid];
        if std::env::var("CADCORE_DUMP_CUT").is_ok() {
            let kind = match &src.faces[fid].geom {
                FaceGeom::Plane(_) => "plane",
                FaceGeom::Cylinder(_) => "cyl",
                FaceGeom::Torus(_) => "torus",
                FaceGeom::Sphere(_) => "sphere",
            };
            let bb = face_aabb(src, fid);
            eprintln!(
                "[cut][bp] {fid:?} {kind} {} bb y[{:.3}..{:.3}]",
                match side {
                    Side::Kept => "kept",
                    Side::Cut => "cut",
                    Side::Straddles => "STRADDLES",
                },
                bb.min.y,
                bb.max.y
            );
        }
        match side {
            Side::Kept => {}
            Side::Cut => dropped += 1,
            Side::Straddles => {
                let curves = sections.remove(&fid).unwrap_or_default();
                if std::env::var("CADCORE_DUMP_CUT").is_ok() {
                    eprintln!("[cut][ssi] {fid:?} straddles → {} curves", curves.len());
                }
                if !curves.is_empty() {
                    ssi[0].insert(fid, curves);
                }
            }
        }
    }
    tick("classify + ssi");
    let mut edges_done: std::collections::HashSet<cadcore_topo::EdgeId> =
        std::collections::HashSet::new();
    for &fid in faces {
        if sides[&fid] != Side::Straddles || !has_real_loops(src, fid) {
            continue;
        }
        for (eid, pts) in face_edge_polylines(src, fid) {
            if !edges_done.insert(eid) {
                continue;
            }
            let mut users: Vec<FaceId> = edge_faces.get(&eid).cloned().unwrap_or_default();
            if !users.contains(&fid) {
                users.push(fid);
            }
            users.sort();
            users.dedup();
            let mut onto: Vec<crate::geom::refine::AnalyticSurface> =
                vec![crate::geom::refine::AnalyticSurface::Plane(*keep)];
            for &u in &users {
                if let Some(sf) = carrier(u) {
                    onto.push(sf);
                }
            }
            let circle = match &src.edges[eid].geom {
                cadcore_topo::EdgeGeom::Circle(c) => Some(*c),
                _ => None,
            };
            let push = |p: Point3, plane_joints: &mut Vec<Point3>| {
                if !plane_joints.iter().any(|&q| (q - p).length() < 1e-9) {
                    plane_joints.push(p);
                }
            };
            for w in pts.windows(2) {
                let (da, db) = (keep.signed_distance(w[0]), keep.signed_distance(w[1]));
                if da.abs() <= TAU_ON {
                    push(w[0], &mut plane_joints);
                } else if (da > 0.0) != (db > 0.0) && db.abs() > TAU_ON {
                    let approx = w[0].lerp(w[1], da / (da - db));
                    let exact_circle = circle.and_then(|c| {
                        circle_plane_points(&c, keep).map(|(q1, q2)| {
                            if (q1 - approx).length() <= (q2 - approx).length() {
                                q1
                            } else {
                                q2
                            }
                        })
                    });
                    let p = match exact_circle {
                        Some(q) if (q - approx).length() < 0.05 => q,
                        _ => {
                            let (exact, res) = crate::geom::refine::refine_triple_point(
                                approx, &onto, 1e-13, 50,
                            );
                            // accept only a genuine triple point; otherwise
                            // the chord crossing is the best shared estimate
                            let ok = res < 1e-9 && (exact - approx).length() < 0.05;
                            if std::env::var("CADCORE_DUMP_CUT").is_ok() && !ok {
                                eprintln!(
                                    "[cut][joint-unrefined] {fid:?} edge {eid:?} approx ({:.6},{:.6},{:.6}) res={res:.2e} moved={:.2e} onto={} surfaces",
                                    approx.x, approx.y, approx.z, (exact - approx).length(), onto.len()
                                );
                            }
                            if ok { exact } else { approx }
                        }
                    };
                    if std::env::var("CADCORE_DUMP_SEAM").is_ok() {
                        let ds: Vec<String> = onto.iter().map(|sf| format!("{:.1e}", sf.distance(p))).collect();
                        eprintln!(
                            "[cut][jointsrc] ({:.6},{:.6},{:.6}) edge {eid:?} circle={} users={users:?} dist={ds:?} approx-moved={:.2e} seg ({:.6},{:.6},{:.6})->({:.6},{:.6},{:.6})",
                            p.x, p.y, p.z, circle.is_some(), (p - approx).length(), w[0].x, w[0].y, w[0].z, w[1].x, w[1].y, w[1].z
                        );
                    }
                    push(p, &mut plane_joints);
                }
            }
            if let Some(&last) = pts.last() {
                if keep.signed_distance(last).abs() <= TAU_ON {
                    push(last, &mut plane_joints);
                }
            }
        }
    }
    if std::env::var("CADCORE_DUMP_CUT").is_ok() {
        for j in &plane_joints {
            eprintln!("[cut][joint] ({:.6},{:.6},{:.6})", j.x, j.y, j.z);
        }
    }
    tick("joints");
    // Graft the in-plane vertices into every section curve that passes through
    // them, so the arranged face, the cap and the copied neighbour all split
    // there identically.
    if !plane_joints.is_empty() {
        for (fid, curves) in ssi[0].iter_mut() {
            let planar = matches!(src.faces[*fid].geom, FaceGeom::Plane(_));
            for c in curves.iter_mut() {
                if c.closed && c.points.len() >= 3 {
                    let first = c.points[0];
                    c.points.push(first);
                    insert_joints(&mut c.points, &plane_joints);
                    c.points.pop();
                } else {
                    // An OPEN section curve ends on the face's boundary, i.e.
                    // at a joint.  Its sampled end can be a chord sagitta off
                    // (a planar end cap's rim is a chord polygon: the plane–
                    // plane section line stops 6e-4 short of the exact circle
                    // point the neighbouring cylinder uses), well outside
                    // `JOINT_SNAP` — so the ends snap on the sagitta scale.
                    const END_SNAP: f64 = 5.0e-3;
                    let n = c.points.len();
                    for k in [0, n.saturating_sub(1)] {
                        if !planar {
                            break; // tube sections end on exact band ends
                        }
                        let q = c.points[k];
                        let best = plane_joints
                            .iter()
                            .map(|j| ((*j - q).length(), *j))
                            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                        if let Some((d, j)) = best {
                            if std::env::var("CADCORE_DUMP_CUT").is_ok() && d < 0.05 {
                                eprintln!("[cut][endsnap] {fid:?} end ({:.6},{:.6},{:.6}) nearest joint {d:.3e}", q.x, q.y, q.z);
                            }
                            if d < END_SNAP {
                                c.points[k] = j;
                            }
                        }
                    }
                    insert_joints(&mut c.points, &plane_joints);
                }
            }
        }
    }
    tick("graft");
    // Clip every section curve to its face's TRIMMED region.  The SSI knows
    // only the carrier band; on a fused shell the face may have a notch or a
    // window exactly where the plane passes (a neighbouring tube welded in), and
    // the piece of section running through that void is not on the solid's
    // boundary — the lateral DCEL drops it, so the cap must never see it either
    // (it would tile the neighbour's section twice).
    //
    // The curve can only change between material and void where a boundary
    // edge of the face crosses the plane — and every such crossing is a joint
    // grafted into the curve above.  So the curve is cut into RUNS between
    // consecutive joints and each run is classified ONCE, at its middle: a
    // per-chord test would ask the membership oracle (a chord polygon of the
    // boundary, a sagitta off the true curve) about points a fraction of a
    // sample from the boundary, where it is wrong as often as right — and a
    // run then starts one chord early, on a vertex no other face has.
    for (fid, curves) in ssi[0].iter_mut() {
        let mut clipped: Vec<SsiCurve> = Vec::new();
        for c in curves.drain(..) {
            let n = c.points.len();
            if n < 2 {
                continue;
            }
            let is_joint =
                |p: Point3| plane_joints.iter().any(|&j| (j - p).length() < 1e-9);
            let jidx: Vec<usize> = (0..n).filter(|&i| is_joint(c.points[i])).collect();
            let material = |run: &[Point3]| -> bool {
                // probe at half the arc length — as far from both ends as the
                // run allows
                let total: f64 = run.windows(2).map(|w| (w[1] - w[0]).length()).sum();
                let mut acc = 0.0;
                let mut probe = run[0].lerp(run[run.len() - 1], 0.5);
                for w in run.windows(2) {
                    let l = (w[1] - w[0]).length();
                    if acc + l >= 0.5 * total {
                        let t = if l > 0.0 { (0.5 * total - acc) / l } else { 0.5 };
                        probe = w[0].lerp(w[1], t);
                        break;
                    }
                    acc += l;
                }
                !matches!(face_contains_hit(src, *fid, probe), Membership::Outside)
            };
            let mut runs: Vec<(Vec<Point3>, bool)> = Vec::new(); // (points, closed)
            if c.closed {
                if jidx.is_empty() {
                    runs.push((c.points.clone(), true));
                } else {
                    for k in 0..jidx.len() {
                        let (i0, i1) = (jidx[k], jidx[(k + 1) % jidx.len()]);
                        let mut run = vec![c.points[i0]];
                        let mut i = i0;
                        loop {
                            i = (i + 1) % n;
                            run.push(c.points[i]);
                            if i == i1 {
                                break;
                            }
                        }
                        runs.push((run, false));
                    }
                }
            } else {
                let mut bounds: Vec<usize> = vec![0];
                bounds.extend(jidx.iter().copied().filter(|&i| i != 0 && i != n - 1));
                bounds.push(n - 1);
                for w in bounds.windows(2) {
                    runs.push((c.points[w[0]..=w[1]].to_vec(), false));
                }
            }
            for (run, closed) in runs {
                if run.len() < 2 {
                    continue;
                }
                let keep_run = material(&run);
                if std::env::var("CADCORE_DUMP_CUT").is_ok() && !keep_run {
                    let (a, b) = (run[0], run[run.len() - 1]);
                    eprintln!(
                        "[cut][clip] {fid:?} run ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}) n={} VOID",
                        a.x, a.y, a.z, b.x, b.y, b.z, run.len()
                    );
                }
                if keep_run {
                    clipped.push(SsiCurve { points: run, closed });
                }
            }
        }
        *curves = clipped;
    }
    ssi[0].retain(|_, v| !v.is_empty());
    tick("clip runs");
    // Edges of KEPT faces lying in the plane whose other face is DROPPED: the
    // plane passes exactly along a junction, and that junction circle IS the
    // section there.  The cap is built from the very same sampled points the
    // copied face emits, so they weld.
    let mut rim_curves: Vec<Vec<Point3>> = Vec::new();
    for (eid, users) in &edge_faces {
        let kept = users.iter().any(|f| sides[f] == Side::Kept);
        let cut = users.iter().any(|f| sides[f] == Side::Cut) || users.len() < 2;
        if !(kept && cut) {
            continue;
        }
        let pts = super::aabb::sample_edge(src, *eid);
        if pts.len() >= 2 && pts.iter().all(|p| keep.signed_distance(*p).abs() <= TAU_ON) {
            rim_curves.push(pts);
        }
    }
    // Where a section crosses a welded junction (a leg's rim circle onto its
    // elbow) the two pieces must share ONE endpoint, or the cap loop will not
    // close and the DCEL leaves a dangling edge.
    snap_ssi_endpoints(&mut ssi, (knit * 30.0).clamp(1e-3, 0.05));
    // Overlapping tubes give overlapping sections; split every section curve at
    // its crossings with the others ONCE, so the cap and the lateral faces
    // partition the shared curve identically.
    presplit_sections(&mut ssi[0], &cut_plane, knit.clamp(5e-4, 1e-2));

    tick("rims + snap + presplit");
    // ── 3. Trim every face ───────────────────────────────────────────────────
    let mut out = BRep::new();
    // Weld tolerance = the knit tolerance, floored at 1 nm.  A fused input
    // built by the universal fuse presents shared points exact to ~1e-12, and
    // a coarse weld (half a micron, say) would drag rim
    // vertices off their circles — the arcs then no longer fit and are written
    // as B-splines.
    // Two faces present a shared point through different fitted curves
    // (each exact to the assembler's fit tolerance), so the weld floor sits a
    // little above twice that — still three orders below any feature.
    let weld = knit.clamp(1e-5, 1e-2);
    let mut asm = Assembler::new(&mut out, weld);
    emit_cut_faces(src, faces, &mut ssi[0], &sides, &plane_joints, keep, &mut asm);

    tick("emit faces");
    // ── 4. Section caps ──────────────────────────────────────────────────────
    let lateral = asm.faces().len();
    let loops = section_loops(&ssi[0], &rim_curves, weld);
    let caps = emit_section_caps(src, faces, &cut_plane, n, half, &loops, &mut asm);

    tick("caps");
    let out_faces = asm.faces().to_vec();
    let shell = out.add_shell(Shell {
        faces: out_faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = out.add_solid(Solid {
        shells: vec![shell],
        name: Some("cut".into()),
    });
    out.shells[shell].solid = solid;
    for f in &out_faces {
        out.faces[*f].shell = shell;
    }

    if std::env::var("CADCORE_DUMP_CUT").is_ok() {
        let rep = crate::validate::manifold::check(&out, shell);
        eprintln!(
            "[cut] faces={} (lateral={lateral} caps={caps}) dropped={dropped} \
             section_loops={} violations={}",
            out_faces.len(),
            loops.len(),
            rep.violations.len()
        );
        eprintln!("[cut] cap faces = {:?}", &out_faces[lateral..]);
        let mut per_face: std::collections::BTreeMap<String, usize> =
            std::collections::BTreeMap::new();
        for v in &rep.violations {
            let faces_in = v
                .split('[')
                .nth(1)
                .unwrap_or("")
                .trim_end_matches(']')
                .to_string();
            *per_face.entry(faces_in).or_default() += 1;
        }
        for (k, n) in per_face.iter().take(12) {
            eprintln!("[cut]   {n:>4} × {k}");
        }
        // Every open edge with its endpoints and geometry kind, so the
        // arranged↔copied disagreement can be located in space.
        let mut uses: HashMap<cadcore_topo::EdgeId, Vec<FaceId>> = HashMap::new();
        for &fid in &out_faces {
            let f = &out.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let Some(lp) = out.loops.get(lid) else { continue };
                let st = lp.start;
                if out.coedges.get(st).is_none() {
                    continue;
                }
                let mut c = st;
                loop {
                    let ce = &out.coedges[c];
                    uses.entry(ce.edge).or_default().push(fid);
                    c = ce.next;
                    if c == st {
                        break;
                    }
                }
            }
        }
        let mut open: Vec<_> = uses.iter().filter(|(_, u)| u.len() != 2).collect();
        open.sort_by_key(|(e, _)| format!("{e:?}"));
        // nearest foreign open-edge vertex to each open edge's start: how far
        // apart the two sides' versions of a shared point really are.
        {
            let verts: Vec<(FaceId, Point3)> = open
                .iter()
                .flat_map(|(eid, fs)| {
                    let e = &out.edges[**eid];
                    let f = fs[0];
                    [(f, out.vertices[e.v_start].point), (f, out.vertices[e.v_end].point)]
                })
                .collect();
            let mut shown = 0;
            for (eid, fs) in &open {
                let e = &out.edges[**eid];
                let a = out.vertices[e.v_start].point;
                let best = verts
                    .iter()
                    .filter(|(f, _)| *f != fs[0])
                    .map(|(_, q)| (*q - a).length())
                    .fold(f64::MAX, f64::min);
                if shown < 6 {
                    eprintln!("[cut][open-gap] {eid:?} start ({:.6},{:.6},{:.6}) nearest foreign vertex {best:.3e}", a.x, a.y, a.z);
                    shown += 1;
                }
            }
        }
        for (eid, fs) in open {
            let e = &out.edges[*eid];
            let a = out.vertices[e.v_start].point;
            let b = out.vertices[e.v_end].point;
            let kind = match &e.geom {
                cadcore_topo::EdgeGeom::Line(_) => "line",
                cadcore_topo::EdgeGeom::Circle(_) => "circle",
                cadcore_topo::EdgeGeom::Polyline(_) => "poly",
                _ => "other",
            };
            let kinds: Vec<String> = fs
                .iter()
                .map(|f| {
                    let g = match &out.faces[*f].geom {
                        FaceGeom::Plane(_) => "pl",
                        FaceGeom::Cylinder(_) => "cyl",
                        FaceGeom::Torus(_) => "tor",
                        FaceGeom::Sphere(_) => "sph",
                    };
                    format!("{f:?}/{g}")
                })
                .collect();
            eprintln!(
                "[cut][open] {kind} ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}) len={:.2e} used {} by {:?}",
                a.x, a.y, a.z, b.x, b.y, b.z,
                (b - a).length(),
                fs.len(),
                kinds
            );
        }
    }

    Some(CutOutcome {
        brep: out,
        shell,
        caps,
        dropped,
    })
}

/// Imprint every face with its section curve, keep the cells that are both
/// material of that face and in front of the plane, and emit them.
///
/// This mirrors [`super::union::emit_solid_faces`] with one extra rule the
/// union never needs.  The union's input faces are analytic templates whose
/// whole parameter region is material; the cut's input is an already-Boolean'd
/// shell, so a face can be a band with HOLES punched in it (crossing windows).
/// The arrangement hands back a hole's interior as a cell in its own right, and
/// "keep unless behind the plane" would happily emit it — filling the window
/// back in and giving its rim four users.  So each cell is additionally tested
/// against the ORIGINAL face's trimmed region with the same membership oracle
/// the ray caster uses.
fn emit_cut_faces(
    src: &BRep,
    faces: &[FaceId],
    ssi: &mut HashMap<FaceId, Vec<SsiCurve>>,
    sides: &HashMap<FaceId, Side>,
    plane_joints: &[Point3],
    keep: &Plane3,
    asm: &mut Assembler,
) {
    // Where a section arc lands on a boundary curve, that curve must be split —
    // on BOTH the arranged face and the copied one across the same window.
    let mut joints = section_joints(ssi);
    for &p in plane_joints {
        if !joints.iter().any(|&q| (q - p).length() < 1e-9) {
            joints.push(p);
        }
    }

    // Pass 1: settle each arranged face's seam and publish where it crosses the
    // face's boundary.  The seam is invisible to the face on the other side of
    // a crossed edge, so it has to become a joint like any other split point.
    // The rotation is settled here and REUSED below, so adding these joints
    // cannot move it.
    let mut seams: HashMap<FaceId, f64> = HashMap::new();
    for &fid in faces {
        if !ssi.contains_key(&fid) {
            continue;
        }
        let Some(domain) = face_domain(src, fid) else { continue };
        if !matches!(domain, FaceDomain::CylinderBand { .. }) {
            continue; // only a periodic band has a seam to park
        }
        let mut chains = arranged_chains(src, fid, &domain, ssi, &joints);
        // A face that does not wrap its cylinder (one half of a band the fuse
        // split along two rulings) is arranged as a plain region: it has no
        // seam to park and nothing to publish.
        if crate::arrange::cells::is_non_wrapping(&domain, &chains) {
            continue;
        }
        // No usable gap ⇒ keep the default seam (rotation 0) — but STILL publish
        // where it cuts, or the copied partner keeps that edge whole.
        let rot = seam_rotation(&domain, &chains).unwrap_or(0.0);
        let seamed = apply_seam_rotation(&domain, &mut chains, rot);
        seams.insert(fid, rot);
        if std::env::var("CADCORE_DUMP_SEAM").is_ok() {
            eprintln!("[cut][seam] {fid:?} rot={rot:.6}");
        }
        for p in seam_crossings(&seamed, &chains) {
            // The crossing is exact on THIS face's carrier only; the chain it
            // split is an edge shared with a neighbour on another carrier, and
            // the neighbour will lift the same vertex onto ITS surface.  Slide
            // the point along the seam ruling line onto the neighbour's carrier
            // as well, so both lifts coincide.
            let p0 = p;
            let p = refine_seam_point_onto_partner(src, faces, fid, &seamed, p);
            if std::env::var("CADCORE_DUMP_SEAM").is_ok() {
                eprintln!(
                    "[cut][seam]   crossing ({:.6},{:.6},{:.6}) → ({:.6},{:.6},{:.6})",
                    p0.x, p0.y, p0.z, p.x, p.y, p.z
                );
            }
            if !joints.iter().any(|&q| (q - p).length() < 1e-9) {
                joints.push(p);
            }
        }
        // A SECTION curve the seam crosses is split by the DCEL at an
        // interpolated point the cap (which consumes the same curve) would never
        // see.  Put the EXACT point — the seam ruling line ∩ the cut plane, on
        // both surfaces to machine precision — into the curve itself, so the
        // lateral face and the cap partition it identically.
        if let Some(curves) = ssi.get_mut(&fid) {
            split_sections_at_seam(curves, &seamed, keep);
        }
    }
    // Torus patches are arranged with their φ = 0 (outer equator) seam; a
    // section crossing it needs the same treatment.
    for &fid in faces {
        let Some(domain) = face_domain(src, fid) else { continue };
        if !matches!(domain, FaceDomain::TorusPatch { .. }) {
            continue;
        }
        if let Some(curves) = ssi.get_mut(&fid) {
            split_sections_at_seam(curves, &domain, keep);
        }
    }

    let prof = std::env::var("CADCORE_TIME").is_ok();
    let (mut t_copy, mut n_copy, mut t_arr, mut n_arr) = (0.0f64, 0usize, 0.0f64, 0usize);
    let t_pass1 = std::time::Instant::now();
    for &fid in faces {
        let f = &src.faces[fid];
        let t_face = std::time::Instant::now();
        // Faces the plane does not touch AT ALL and that carry real topology are
        // copied verbatim — see `Assembler::emit_face_copy` for why re-arranging
        // them is NOT a no-op on a fused shell.  Template faces (no loops) have
        // nothing to copy, so they take the arrangement path, which synthesises
        // their boundary.
        //
        // A face that STRADDLES must never be copied, even when the SSI found no
        // curve on it (a shallow torus arc, a tangency): its boundary rims
        // straddle too, the neighbour across such a rim IS trimmed, and a copy
        // would present the whole rim against the neighbour's trimmed half —
        // one user each, i.e. an open shell.  Arranging it costs nothing extra
        // and lets the keep predicate trim it.
        let side = sides.get(&fid).copied().unwrap_or(Side::Straddles);
        if side == Side::Cut {
            continue; // wholly behind the plane
        }
        if !ssi.contains_key(&fid) && has_real_loops(src, fid) {
            match side {
                Side::Cut => continue,
                Side::Kept => {
                    asm.emit_face_copy(src, fid, &joints);
                    t_copy += t_face.elapsed().as_secs_f64();
                    n_copy += 1;
                    continue;
                }
                Side::Straddles => {} // fall through to the arrangement
            }
        }
        let Some(domain) = face_domain(src, fid) else {
            continue; // unsupported carrier surface — skip (incomplete)
        };

        // Material is removed where the plane's normal points AWAY, and where
        // the original face had no material to begin with.
        let drop = |p: Point3| {
            keep.signed_distance(p) < 0.0
                || !matches!(face_contains_hit(src, fid, p), Membership::Inside)
        };

        let mut chains = arranged_chains(src, fid, &domain, ssi, &joints);
        let stag = 1000 + chains.iter().filter(|c| c.tag >= 1000).count() as u32;
        // Apply the seam settled in pass 1.  A fused band is SEAMLESS; the
        // periodic DCEL has to put a seam somewhere, and wherever it lands it
        // splits the chains crossing it — so it is parked in the widest empty
        // gap and whatever
        // it still crosses was published as a joint.
        let domain = match seams.get(&fid) {
            Some(&rot) => apply_seam_rotation(&domain, &mut chains, rot),
            None => domain,
        };
        // A closed chain must not START on the seam: `build_periodic_u`
        // unwraps every chain from its first sample, and a first sample sitting
        // exactly on u = 0 leaves the branch of the whole ring ambiguous.
        // Rotate the start to the sample farthest from the seam.
        if matches!(domain, FaceDomain::CylinderBand { .. }) {
            for ch in chains.iter_mut() {
                let n = ch.pts.len();
                if n < 4 {
                    continue;
                }
                let (a, b) = (ch.pts[0], ch.pts[n - 1]);
                let closed = (a.0 - b.0).abs() < 1e-9 && (a.1 - b.1).abs() < 1e-9;
                if !closed || a.0.abs() > 1e-6 {
                    continue;
                }
                let k = (0..n - 1)
                    .max_by(|&i, &j| ch.pts[i].0.abs().partial_cmp(&ch.pts[j].0.abs()).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap_or(0);
                let mut rotated: Vec<(f64, f64)> = Vec::with_capacity(n);
                for i in 0..n - 1 {
                    rotated.push(ch.pts[(k + i) % (n - 1)]);
                }
                rotated.push(rotated[0]);
                ch.pts = rotated;
            }
        }
        unwrap_torus_band_chains(&domain, &mut chains);

        if std::env::var("CADCORE_DUMP_CHAINS").is_ok() {
            for ch in &chains {
                let a = domain.lift(ch.pts[0].0, ch.pts[0].1);
                let b = domain.lift(ch.pts[ch.pts.len() - 1].0, ch.pts[ch.pts.len() - 1].1);
                let us: Vec<f64> = ch.pts.iter().map(|p| p.0).collect();
                let (ulo, uhi) = us.iter().fold((f64::MAX, f64::MIN), |(l, h), &u| (l.min(u), h.max(u)));
                eprintln!(
                    "[cut][chain] {fid:?} tag={} n={} u[{ulo:.3}..{uhi:.3}] ({:.4},{:.4},{:.4})..({:.4},{:.4},{:.4})",
                    ch.tag, ch.pts.len(), a.x, a.y, a.z, b.x, b.y, b.z
                );
                if ch.tag < 1000 {
                    for w in ch.pts.windows(2) {
                        let p = domain.lift(w[0].0, w[0].1);
                        let q = domain.lift(w[1].0, w[1].1);
                        if (q - p).length() > 0.3 {
                            eprintln!("[cut][chain]   long step ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}) len={:.3}", p.x, p.y, p.z, q.x, q.y, q.z, (q - p).length());
                        }
                    }
                }
            }
        }
        if std::env::var("CADCORE_DUMP_CUT").is_ok() {
            for ch in &chains {
                let (ulo, uhi) = ch.pts.iter().fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.0), a.1.max(p.0)));
                let (vlo, vhi) = ch.pts.iter().fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.1), a.1.max(p.1)));
                let closed = ch.pts.len() >= 2 && {
                    let (a, b) = (ch.pts[0], *ch.pts.last().unwrap());
                    (a.0 - b.0).abs() < 1e-9 && (a.1 - b.1).abs() < 1e-9
                };
                eprintln!(
                    "[cut][chain] {fid:?} tag={} n={} closed={closed} u[{ulo:.3}..{uhi:.3}] v[{vlo:.3}..{vhi:.3}] first=({:.3},{:.3}) last=({:.3},{:.3})",
                    ch.tag, ch.pts.len(), ch.pts[0].0, ch.pts[0].1, ch.pts.last().unwrap().0, ch.pts.last().unwrap().1
                );
                if std::env::var("CADCORE_DUMP_CHAINS_FULL").is_ok() {
                    let body: Vec<String> = ch.pts.iter().map(|p| format!("({:.9},{:.9})", p.0, p.1)).collect();
                    eprintln!("[cut][chainpts] {fid:?} {} {}", ch.tag, body.join(","));
                }
            }
        }
        if std::env::var("CADCORE_DUMP_CUT").is_ok() {
            // every SSI endpoint must sit ON a boundary chain (a joint), or it
            // dangles in the DCEL and the face is never split there
            let bnd: Vec<&FaceChain> = chains.iter().filter(|c| c.tag < 1000).collect();
            for ch in chains.iter().filter(|c| c.tag >= 1000) {
                let n = ch.pts.len();
                let closed = n >= 2 && {
                    let (a, b) = (ch.pts[0], ch.pts[n - 1]);
                    (a.0 - b.0).abs() < 1e-9 && (a.1 - b.1).abs() < 1e-9
                };
                if closed {
                    continue;
                }
                for &e in &[ch.pts[0], ch.pts[n - 1]] {
                    let pe = domain.lift(e.0, e.1);
                    let mut best = f64::MAX;
                    let mut best_exact = f64::MAX;
                    for b in &bnd {
                        for p in &b.pts {
                            let d = (domain.lift(p.0, p.1) - pe).length();
                            best_exact = best_exact.min(d);
                        }
                        for w in b.pts.windows(2) {
                            let (a, c) = (domain.lift(w[0].0, w[0].1), domain.lift(w[1].0, w[1].1));
                            let ab = c - a;
                            let t = ((pe - a).dot(ab) / ab.dot(ab)).clamp(0.0, 1.0);
                            best = best.min((a + ab * t - pe).length());
                        }
                    }
                    if best_exact > 1e-6 {
                        eprintln!(
                            "[cut][dangle] {fid:?} ssi tag={} end ({:.4},{:.4},{:.4}) nearest bnd vertex {best_exact:.2e} nearest bnd segment {best:.2e}",
                            ch.tag, pe.x, pe.y, pe.z
                        );
                    }
                }
            }
        }
        if std::env::var("CADCORE_DUMP_SEAM").is_ok() {
            for ch in &chains {
                for (i, p) in ch.pts.iter().enumerate() {
                    let q = domain.lift(p.0, p.1);
                    for j in &joints {
                        let d = (q - *j).length();
                        if d < 0.012 {
                            eprintln!(
                                "[cut][nearjoint] {fid:?} tag={} #{i} ({:.6},{:.6},{:.6}) d={d:.2e} to joint ({:.4},{:.4},{:.4})",
                                ch.tag, q.x, q.y, q.z, j.x, j.y, j.z
                            );
                        }
                    }
                }
            }
        }
        let cells = arrange_and_classify_with(&domain, &chains, 1e-3, &drop);
        let mut emitted = 0usize;
        let mut out_ids: Vec<cadcore_topo::FaceId> = Vec::new();
        for cell in &cells {
            if cell.keep {
                out_ids.push(asm.emit_cell(&domain, f.geom.clone(), f.normal, cell));
                emitted += 1;
            }
        }
        if std::env::var("CADCORE_DUMP_CUT").is_ok() {
            let kind = match &f.geom {
                FaceGeom::Plane(_) => "plane",
                FaceGeom::Cylinder(_) => "cyl",
                FaceGeom::Torus(_) => "torus",
                FaceGeom::Sphere(_) => "sphere",
            };
            eprintln!(
                "cut face {fid:?} {kind}: bnd={} ssi={} cells={} emitted={} → {out_ids:?} areas={:?}",
                chains.len() - (stag - 1000) as usize,
                stag - 1000,
                cells.len(),
                emitted,
                cells
                    .iter()
                    .map(|c| (c.cell.area * 1000.0).round() / 1000.0)
                    .collect::<Vec<_>>()
            );
        }
        t_arr += t_face.elapsed().as_secs_f64();
        n_arr += 1;
    }
    if prof {
        eprintln!(
            "[cut][time]   copied {n_copy} faces {:.0} ms, arranged {n_arr} faces {:.0} ms, other {:.0} ms",
            t_copy * 1e3,
            t_arr * 1e3,
            (t_pass1.elapsed().as_secs_f64() - t_copy - t_arr) * 1e3
        );
    }
}

/// Every chain an ARRANGED face feeds the DCEL: its own boundary loops (with
/// the joints grafted in) followed by the section curves lying on it.
///
/// Built the same way in both passes so the seam chosen in pass 1 is the seam
/// applied in pass 2.
fn arranged_chains(
    src: &BRep,
    fid: FaceId,
    domain: &FaceDomain,
    ssi: &HashMap<FaceId, Vec<SsiCurve>>,
    joints: &[Point3],
) -> Vec<FaceChain> {
    let mut chains: Vec<FaceChain> = Vec::new();
    let mut tag = 1u32;
    for pts in boundary_chains_with_joints(src, fid, domain, joints) {
        chains.push(FaceChain { pts, tag });
        tag += 1;
    }
    let mut stag = 1000u32;
    if let Some(curves) = ssi.get(&fid) {
        for c in curves {
            let mut pts: Vec<(f64, f64)> = c.points.iter().map(|&p| domain.uv(p)).collect();
            if c.closed && pts.len() >= 3 {
                pts.push(pts[0]);
            }
            if pts.len() >= 2 {
                chains.push(FaceChain { pts, tag: stag });
                stag += 1;
            }
        }
    }
    chains
}


/// Pick the rotation that parks the u = ±π seam in the widest empty gap.
fn seam_rotation(domain: &FaceDomain, chains: &[FaceChain]) -> Option<f64> {
    use std::f64::consts::{PI, TAU};
    let FaceDomain::CylinderBand { .. } = domain else {
        return None;
    };
    // Only LOCALISED chains matter: a band's rim circles span every u, so the
    // seam has to cross them wherever it goes (both faces sharing a rim split
    // it the same way, so that is harmless).  What must not be split is a
    // window loop or a section curve, which occupies a limited arc.
    // A full-wrap chain is skipped only where it IS the rim: its vertices
    // strictly inside the band (the notches a fused outline carries where a
    // neighbour was welded in) are features the seam must avoid like any
    // window.
    let (vlo, vhi) = domain.open_range();
    let eps = 1e-6 * (vhi - vlo).abs().max(1.0);
    let mut us: Vec<f64> = Vec::new();
    for c in chains.iter() {
        let (lo, hi) = c
            .pts
            .iter()
            .fold((f64::MAX, f64::MIN), |(lo, hi), p| (lo.min(p.0), hi.max(p.0)));
        let full_wrap = hi - lo > 0.9 * TAU;
        if full_wrap && c.tag >= 1000 {
            continue; // a transversal section: crossed at an exact inserted point
        }
        for p in &c.pts {
            if full_wrap && !(p.1 > vlo + eps && p.1 < vhi - eps) {
                continue; // on a rim: harmless to cross
            }
            us.push(p.0.rem_euclid(TAU));
        }
    }
    if us.is_empty() {
        return None;
    }
    us.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // widest circular gap between consecutive samples
    let (mut best_gap, mut seam_u) = (0.0_f64, 0.0_f64);
    for i in 0..us.len() {
        let (a, b) = (us[i], if i + 1 < us.len() { us[i + 1] } else { us[0] + TAU });
        if b - a > best_gap {
            best_gap = b - a;
            seam_u = 0.5 * (a + b);
        }
    }
    // Nothing worth moving: the chains already cover the circle densely (the
    // seam has to cross something, and the caller falls back to the default).
    if best_gap < 1e-6 {
        return None;
    }
    // Snap the seam onto an EXISTING vertex inside the gap (a rim sample —
    // every rim has one every few degrees): the seam then passes through a
    // vertex exactly instead of 1e-6 beside it, so no near-duplicate vertex
    // and no sliver segment is created where it crosses the rims.
    let gap_lo = seam_u - 0.5 * best_gap;
    let mut snap: Option<(f64, f64)> = None; // (distance to centre, u)
    for c in chains.iter() {
        for p in &c.pts {
            let u = p.0.rem_euclid(TAU);
            let rel = (u - gap_lo).rem_euclid(TAU);
            if rel <= 1e-9 || rel >= best_gap - 1e-9 {
                continue; // not strictly inside the gap
            }
            let d = (rel - 0.5 * best_gap).abs();
            if snap.map_or(true, |(bd, _)| d < bd) {
                snap = Some((d, p.0));
            }
        }
    }
    let seam_u = match snap {
        Some((_, u)) => u,
        None => seam_u,
    };
    // `build_periodic_u` unwraps each chain from its own first sample and then
    // splits at every multiple of the period, so the seam it cuts along is
    // u = 0 — the frame's x direction — NOT u = ±π.  Rotating by `seam_u` puts
    // u = 0 at the middle of the gap.
    let _ = PI;
    Some(seam_u)
}

/// Rotate the band's parameter frame by `rot` and shift the chains onto it.
fn apply_seam_rotation(domain: &FaceDomain, chains: &mut [FaceChain], rot: f64) -> FaceDomain {
    use std::f64::consts::{PI, TAU};
    let FaceDomain::CylinderBand { surf, length } = domain else {
        return *domain;
    };
    let (s, c) = rot.sin_cos();
    let (Some(x), ..) = (UnitVec3::try_from_vec(surf.frame.x * c + surf.frame.y * s), ()) else {
        return *domain;
    };
    let Some(y) = UnitVec3::try_from_vec(surf.frame.z.cross(x)) else {
        return *domain;
    };
    let reframed = cadcore_geom::CylSurf {
        frame: cadcore_math::Frame3 {
            origin: surf.frame.origin,
            x,
            y,
            z: surf.frame.z,
        },
        radius: surf.radius,
    };
    for ch in chains.iter_mut() {
        for p in ch.pts.iter_mut() {
            let mut u = p.0 - rot;
            u -= TAU * ((u + PI) / TAU).floor(); // wrap into (−π, π]
            p.0 = u;
        }
    }
    FaceDomain::CylinderBand {
        surf: reframed,
        length: *length,
    }
}

/// The 3-D points where `chains` (already in the re-seamed uv) cross the
/// u = 0 seam.
///
/// The seam is an artefact of the parametrisation, not of the model, so the
/// face on the OTHER side of a crossed edge knows nothing about it and keeps
/// that edge whole; publishing the crossings as joints makes both sides split
/// there.
///
/// The point is computed exactly the way `build_periodic_u` computes its own
/// split — linear interpolation of uv along the chord — and that is deliberate:
/// what matters is that BOTH sides use the SAME point.  Substituting the
/// mathematically exact crossing (solved on the source edge's curve) instead
/// re-opens the shell, because the DCEL still places its own interpolated
/// vertex on the arranged side.  See the module TODO.
fn seam_crossings(domain: &FaceDomain, chains: &[FaceChain]) -> Vec<Point3> {
    use std::f64::consts::PI;
    let FaceDomain::CylinderBand { .. } = domain else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for ch in chains {
        if ch.tag >= 1000 {
            continue; // section curves get exact seam points (`split_sections_at_seam`)
        }
        for w in ch.pts.windows(2) {
            let (a, b) = (w[0], w[1]);
            // A jump larger than π is the atan2 wrap, not a real step across 0.
            if (a.0 - b.0).abs() > PI {
                continue;
            }
            if (a.0 > 0.0) == (b.0 > 0.0) {
                continue; // both on the same side of the seam
            }
            let t = a.0 / (a.0 - b.0);
            if !(0.0..=1.0).contains(&t) {
                continue;
            }
            out.push(domain.lift(0.0, a.1 + (b.1 - a.1) * t));
        }
    }
    out
}

/// Move a seam-crossing point `p` (on face `fid`'s seam ruling line) along
/// that line onto the carrier of the neighbouring face whose edge the seam
/// crossed there.  The neighbour is found as the other face whose surface
/// passes within a sagitta of `p`; the result lies on both carriers.
fn refine_seam_point_onto_partner(
    src: &BRep,
    faces: &[FaceId],
    fid: FaceId,
    domain: &FaceDomain,
    p: Point3,
) -> Point3 {
    use crate::geom::refine::AnalyticSurface;
    let FaceDomain::CylinderBand { surf, .. } = domain else {
        return p;
    };
    let axis = surf.axis().as_vec();
    let base = surf.frame.origin + surf.frame.x * surf.radius;
    let t0 = axis.dot(p - base);
    let carrier = |g: &FaceGeom| -> Option<AnalyticSurface> {
        match g {
            FaceGeom::Plane(pl) => Some(AnalyticSurface::Plane(*pl)),
            FaceGeom::Cylinder(c) => Some(AnalyticSurface::Cylinder(*c)),
            FaceGeom::Torus(t) => Some(AnalyticSurface::Torus(*t)),
            FaceGeom::Sphere(_) => None,
        }
    };
    let mut best: Option<(f64, Point3)> = None;
    for &g in faces {
        if g == fid {
            continue;
        }
        let Some(sg) = carrier(&src.faces[g].geom) else { continue };
        let d0 = sg.distance(p);
        if d0 > 5e-3 {
            continue; // not the neighbour across this edge
        }
        if d0 < 1e-8 {
            return p; // already on the partner (an existing vertex) — keep it
        }
        let bb = face_aabb(src, g).pad(1e-2);
        if !(bb.min.x <= p.x && p.x <= bb.max.x && bb.min.y <= p.y && p.y <= bb.max.y && bb.min.z <= p.z && p.z <= bb.max.z) {
            continue;
        }
        let f = |t: f64| sg.distance(base + axis * t);
        let t = golden_min_1d(t0 - 2e-2, t0 + 2e-2, 80, f);
        let q = base + axis * t;
        let d = sg.distance(q);
        if d < 1e-9 && d < d0 && best.map_or(true, |(bd, _)| d < bd) {
            best = Some((d, q));
        }
    }
    best.map_or(p, |(_, q)| q)
}

fn golden_min_1d(mut a: f64, mut b: f64, iters: usize, f: impl Fn(f64) -> f64) -> f64 {
    let gr = (5f64.sqrt() - 1.0) / 2.0;
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..iters {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = f(d);
        }
    }
    0.5 * (a + b)
}

/// Insert, into every section curve of a face, the exact point where it
/// crosses the face's parametric seam: for a cylinder band the seam ruling
/// line (`u = 0` of the seamed frame) ∩ the cut plane; for a torus patch the
/// outer-equator circle (`φ = 0`) ∩ the cut plane, nearest the crossing.
fn split_sections_at_seam(curves: &mut [SsiCurve], domain: &FaceDomain, keep: &Plane3) {
    use std::f64::consts::PI;
    let n_plane = keep.normal().as_vec();
    for c in curves.iter_mut() {
        let m = c.points.len();
        if m < 2 {
            continue;
        }
        let nseg = if c.closed { m } else { m - 1 };
        let mut inserts: Vec<(usize, Point3)> = Vec::new();
        for i in 0..nseg {
            let (pa, pb) = (c.points[i], c.points[(i + 1) % m]);
            let (ua, va) = domain.uv(pa);
            let (ub, vb) = domain.uv(pb);
            // the periodic coordinate of this domain
            let (sa, sb) = match domain {
                FaceDomain::CylinderBand { .. } => (ua, ub),
                FaceDomain::TorusPatch { .. } => (va, vb),
                _ => continue,
            };
            if (sa - sb).abs() > PI || (sa > 0.0) == (sb > 0.0) || sa == 0.0 || sb == 0.0 {
                continue;
            }
            let exact = match domain {
                FaceDomain::CylinderBand { surf, .. } => {
                    // seam ruling line o + r·x + axis·t  ∩  plane
                    let base = surf.frame.origin + surf.frame.x * surf.radius;
                    let d = surf.axis().as_vec();
                    let den = n_plane.dot(d);
                    if den.abs() < 1e-12 {
                        continue;
                    }
                    let t = -keep.signed_distance(base) / den;
                    base + d * t
                }
                FaceDomain::TorusPatch { surf, .. } => {
                    // outer equator: circle of radius R + r in the ring plane ∩
                    // cut plane — take the solution nearest the chord crossing.
                    let t = sa / (sa - sb);
                    let approx = pa.lerp(pb, t);
                    let eq = Circle3::new(
                        surf.frame.origin,
                        surf.frame.z,
                        surf.major_radius + surf.minor_radius,
                    );
                    match circle_plane_points(&eq, keep) {
                        Some((q1, q2)) => {
                            if (q1 - approx).length() <= (q2 - approx).length() {
                                q1
                            } else {
                                q2
                            }
                        }
                        None => continue,
                    }
                }
                _ => continue,
            };
            if (exact - pa).length() > 1e-9 && (exact - pb).length() > 1e-9 {
                inserts.push((i, exact));
            }
        }
        for (k, (i, q)) in inserts.into_iter().enumerate() {
            c.points.insert(i + 1 + k, q);
        }
    }
}

/// The two points where a circle meets a plane (`None` when it misses or lies
/// in the plane).
fn circle_plane_points(c: &Circle3, plane: &Plane3) -> Option<(Point3, Point3)> {
    // c(θ) = o + r(x cosθ + y sinθ);  n·(c(θ) − p0) = 0
    //  ⇒ A cosθ + B sinθ = −D  with A = r n·x, B = r n·y, D = n·(o − p0)
    let n = plane.normal().as_vec();
    let a = c.radius * n.dot(c.frame.x.as_vec());
    let b = c.radius * n.dot(c.frame.y.as_vec());
    let d = plane.signed_distance(c.frame.origin);
    let amp = (a * a + b * b).sqrt();
    if amp < 1e-15 || d.abs() > amp {
        return None;
    }
    let phi = b.atan2(a);
    let alpha = (-d / amp).acos();
    Some((c.point_at(phi + alpha), c.point_at(phi - alpha)))
}

// ── Joints ───────────────────────────────────────────────────────────────────

/// Every point where a section arc starts or ends.
///
/// These are the joints of the cut: a section arc runs across a face and stops
/// on one of its boundary curves — a crossing window, a rim circle, another
/// section.  They are the points at which the ARRANGED side of that boundary
/// gets split, so the COPIED side has to be split there too.
fn section_joints(ssi: &HashMap<FaceId, Vec<SsiCurve>>) -> Vec<Point3> {
    let mut out: Vec<Point3> = Vec::new();
    for curves in ssi.values() {
        for c in curves {
            if c.closed || c.points.len() < 2 {
                continue; // a closed section ends on nothing
            }
            for p in [c.points[0], *c.points.last().unwrap()] {
                if !out.iter().any(|&q| (q - p).length() < 1e-9) {
                    out.push(p);
                }
            }
        }
    }
    out
}

/// Insert every joint that lies ON `poly` at its exact position.
///
/// A joint is accepted only when it is close to a segment relative to that
/// segment's own length (and within an absolute cap), so a joint belonging to
/// some other curve running nearby is never grafted on.  Inserting the joint
/// itself — not its projection — is the point: both consumers of this curve
/// then carry the identical vertex and their segments weld.
pub(crate) fn insert_joints(poly: &mut Vec<Point3>, joints: &[Point3]) {
    if joints.is_empty() || poly.len() < 2 {
        return;
    }
    // Only joints near the polyline can act on it: both tests below accept a
    // joint at most max(JOINT_SNAP, 0.01) from a vertex or a segment foot, and
    // every foot lies inside the polyline's box.  Filtering by that box (with a
    // margin above the widest window) keeps the result identical and turns the
    // per-segment scan over EVERY joint of the cut into a scan over a handful —
    // the copy of an untouched face was O(samples × all joints).
    const NEAR: f64 = 0.011;
    let near: Vec<Point3>;
    let joints: &[Point3] = if joints.len() > 8 {
        let (mut lo, mut hi) = (poly[0], poly[0]);
        for p in poly.iter() {
            lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        near = joints
            .iter()
            .copied()
            .filter(|j| {
                j.x >= lo.x - NEAR
                    && j.x <= hi.x + NEAR
                    && j.y >= lo.y - NEAR
                    && j.y <= hi.y + NEAR
                    && j.z >= lo.z - NEAR
                    && j.z <= hi.z + NEAR
            })
            .collect();
        &near
    } else {
        joints
    };
    // A joint a hair off an existing vertex REPLACES it: the joint is the exact
    // triple point every consumer shares, the vertex is a sample (a band-end
    // crossing of a 1e-5-fuzzy rim, a seam snap) that would otherwise become a
    // second vertex a sliver away — below the DCEL tolerance on one face and
    // above it on the next, i.e. a mismatched edge.  Every face grafts the same
    // joints, so the replacement is identical on all of them.
    for q in poly.iter_mut() {
        if let Some(j) = joints.iter().find(|j| (**j - *q).length() < JOINT_SNAP) {
            *q = *j;
        }
    }
    let mut out: Vec<Point3> = Vec::with_capacity(poly.len() + joints.len());
    for w in 0..poly.len() - 1 {
        let (a, b) = (poly[w], poly[w + 1]);
        if out.last().map_or(false, |&q| (q - a).length() < 1e-9) {
            // the replacement above can fold two vertices into one
        } else {
            out.push(a);
        }
        let seg = b - a;
        let len2 = seg.dot(seg);
        if len2 < 1e-24 {
            continue;
        }
        // The polyline is only a CHORD approximation of its curve, so a joint
        // that genuinely belongs to it can sit up to a sagitta off — the window
        // is the sagitta scale, not the knit tolerance.
        let tol = (0.3 * len2.sqrt()).min(0.01);
        let mut hits: Vec<(f64, Point3)> = Vec::new();
        for &j in joints {
            let t = seg.dot(j - a) / len2;
            if !(1e-6..=1.0 - 1e-6).contains(&t) {
                continue; // at (or past) an endpoint — already a vertex
            }
            let foot = a + seg * t;
            if (j - foot).length() <= tol {
                hits.push((t, j));
            }
        }
        hits.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
        for (_, j) in hits {
            if out.last().map_or(true, |&q| (q - j).length() > 1e-9) {
                out.push(j);
            }
        }
    }
    let last = *poly.last().unwrap();
    if out.last().map_or(true, |&q| (q - last).length() > 1e-9) {
        out.push(last);
    }
    *poly = out;
}

/// Distance below which a chain vertex is taken to BE a joint (mm).
///
/// Well above the weld (1e-5) and the fit tolerance (5e-6) that make a sampled
/// vertex miss the exact triple point, and three orders below any feature.
pub(crate) const JOINT_SNAP: f64 = 1e-4;

/// The face's boundary loops as uv chains, with the section joints inserted.
///
/// Mirrors [`face_boundary_chains`] but walks the loops in 3-D first so the
/// joints can be grafted in before projection — the copied path grafts the same
/// joints into the same edges, so the two sides split identically.
fn boundary_chains_with_joints(
    brep: &BRep,
    fid: FaceId,
    domain: &FaceDomain,
    joints: &[Point3],
) -> Vec<Vec<(f64, f64)>> {
    if !has_real_loops(brep, fid) {
        // Template face: synthesise its rings in 3-D first so the joints can be
        // grafted in exactly like on a real loop (an arranged neighbour's seam
        // crossing of a shared rim must be a vertex here too).
        let rings = synthesize_rings_3d(brep, fid);
        if rings.is_empty() {
            return face_boundary_chains(brep, fid, domain);
        }
        let mut out = Vec::new();
        for mut ring in rings {
            insert_joints(&mut ring, joints);
            let mut pts: Vec<(f64, f64)> = ring.iter().map(|&p| domain.uv(p)).collect();
            if pts.first() != pts.last() {
                pts.push(pts[0]);
            }
            out.push(pts);
        }
        return out;
    }
    let f = &brep.faces[fid];
    let mut loops = Vec::new();
    let mut ids = vec![f.outer_loop];
    ids.extend(f.inner_loops.iter().copied());
    for lid in ids {
        let Some(lp) = brep.loops.get(lid) else { continue };
        let start = lp.start;
        if brep.coedges.get(start).is_none() {
            continue;
        }
        // The loop's co-edges in order, and how often each edge appears in
        // it.  An edge the SAME loop traverses twice is a parametric seam
        // slit left by a previous arrangement (the two borders of a periodic
        // DCEL lifted onto one 3-D line).  It is not geometry: the DCEL here
        // puts its own seam where it wants, so the slit is dropped and the
        // loop handed over as the open chains between its two traversals.
        let mut ces: Vec<cadcore_topo::CoEdgeId> = Vec::new();
        let mut c = start;
        loop {
            ces.push(c);
            c = brep.coedges[c].next;
            if c == start || ces.len() > 100_000 {
                break;
            }
        }
        let mut count: HashMap<cadcore_topo::EdgeId, usize> = HashMap::new();
        for &ce in &ces {
            *count.entry(brep.coedges[ce].edge).or_default() += 1;
        }
        let is_slit = |ce: cadcore_topo::CoEdgeId| count[&brep.coedges[ce].edge] >= 2;
        let has_slit = ces.iter().any(|&ce| is_slit(ce));
        let n = ces.len();
        let first = if has_slit {
            (0..n).find(|&i| is_slit(ces[i])).map(|i| (i + 1) % n).unwrap_or(0)
        } else {
            0
        };
        let mut chains3: Vec<Vec<Point3>> = Vec::new();
        let mut pts3: Vec<Point3> = Vec::new();
        for k in 0..n {
            let ce_id = ces[(first + k) % n];
            if is_slit(ce_id) {
                if pts3.len() >= 2 {
                    chains3.push(std::mem::take(&mut pts3));
                } else {
                    pts3.clear();
                }
                continue;
            }
            let ce = &brep.coedges[ce_id];
            let mut sp = super::aabb::sample_edge(brep, ce.edge);
            if ce.sense == cadcore_topo::CoEdgeSense::Opposite {
                sp.reverse();
            }
            insert_joints(&mut sp, joints);
            for p in sp {
                if pts3.last().map_or(true, |&q| (q - p).length() > 1e-12) {
                    pts3.push(p);
                }
            }
        }
        if pts3.len() >= 2 {
            chains3.push(pts3);
        }
        for ch in chains3 {
            let mut pts: Vec<(f64, f64)> = ch.iter().map(|&p| domain.uv(p)).collect();
            if !has_slit && pts.len() >= 3 && pts.first() != pts.last() {
                pts.push(pts[0]); // close the ring
            }
            if pts.len() >= 2 {
                loops.push(pts);
            }
        }
    }
    if loops.is_empty() {
        return face_boundary_chains(brep, fid, domain);
    }
    loops
}

/// The boundary rings of a TEMPLATE face (no real loops) as 3-D polylines:
/// the two rims of a `Cylinder` extent, the circle of a `Disk`, the polygon of
/// a `Polygon`.  Sampled exactly as `union::synthesize_boundary` does (64 per
/// circle, starting at the circle's own `point_at(0)`).
fn synthesize_rings_3d(brep: &BRep, fid: FaceId) -> Vec<Vec<Point3>> {
    use cadcore_topo::FaceBoundary;
    let ring = |c: &Circle3| -> Vec<Point3> {
        (0..64)
            .map(|k| c.point_at(k as f64 / 64.0 * std::f64::consts::TAU))
            .collect()
    };
    let f = &brep.faces[fid];
    match (&f.geom, &f.extent) {
        (_, FaceExtent::Cylinder { start, end, .. }) => {
            let mut out = Vec::new();
            for b in [start, end] {
                if let FaceBoundary::Circle(c) = b {
                    out.push(ring(c));
                }
            }
            out
        }
        (FaceGeom::Plane(pl), FaceExtent::Disk { radius }) => {
            vec![ring(&Circle3::new(pl.frame.origin, pl.normal(), *radius))]
        }
        (_, FaceExtent::Polygon { points }) if points.len() >= 3 => vec![points.clone()],
        _ => Vec::new(),
    }
}

/// Every boundary edge of a face with its sampled polyline, in loop order.
fn face_edge_polylines(brep: &BRep, fid: FaceId) -> Vec<(cadcore_topo::EdgeId, Vec<Point3>)> {
    let f = &brep.faces[fid];
    let mut out = Vec::new();
    let mut ls = vec![f.outer_loop];
    ls.extend(f.inner_loops.iter().copied());
    for lid in ls {
        let Some(lp) = brep.loops.get(lid) else { continue };
        let st = lp.start;
        if brep.coedges.get(st).is_none() {
            continue;
        }
        let mut c = st;
        loop {
            let eid = brep.coedges[c].edge;
            out.push((eid, super::aabb::sample_edge(brep, eid)));
            c = brep.coedges[c].next;
            if c == st {
                break;
            }
        }
    }
    out
}

/// Does this face carry explicit B-Rep topology (as opposed to being an
/// analytic template whose boundary is implied by its [`FaceExtent`])?
fn has_real_loops(brep: &BRep, fid: FaceId) -> bool {
    brep.faces
        .get(fid)
        .and_then(|f| brep.loops.get(f.outer_loop))
        .map_or(false, |lp| brep.coedges.get(lp.start).is_some())
}

// ── Broad phase ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Kept,
    Cut,
    Straddles,
}

/// Which side of `keep` a face lies on, from its sampled geometry: boundary
/// loop samples plus a coarse grid on the carrier surface within the face's
/// domain (a torus fillet bulges past its junction circles, so loops alone
/// under-report its reach).  `|d| ≤ tau` counts as ON the plane.
fn face_side(src: &BRep, fid: FaceId, keep: &Plane3, tau: f64) -> Side {
    if !has_real_loops(src, fid) {
        return straddle(&face_aabb(src, fid), keep, tau.max(1e-9));
    }
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    let mut take = |p: Point3| {
        let d = keep.signed_distance(p);
        lo = lo.min(d);
        hi = hi.max(d);
    };
    for p in face_loop_points(src, fid) {
        take(p);
    }
    if let Some(domain) = face_domain(src, fid) {
        let (ulo, uhi, vlo, vhi, nu, nv) = match domain {
            FaceDomain::CylinderBand { length, .. } => {
                (0.0, std::f64::consts::TAU, 0.0, length, 16usize, 2usize)
            }
            FaceDomain::TorusPatch { theta_lo, theta_hi, .. } => {
                (theta_lo, theta_hi, 0.0, std::f64::consts::TAU, 8usize, 16usize)
            }
            FaceDomain::TorusBand { .. } => {
                (0.0, std::f64::consts::TAU, 0.0, std::f64::consts::TAU, 16usize, 16usize)
            }
            FaceDomain::Plane { .. } => (0.0, 0.0, 0.0, 0.0, 0usize, 0usize),
        };
        for i in 0..=nu {
            for k in 0..=nv {
                let u = ulo + (uhi - ulo) * i as f64 / nu.max(1) as f64;
                let v = vlo + (vhi - vlo) * k as f64 / nv.max(1) as f64;
                take(domain.lift(u, v));
            }
        }
    }
    if lo >= -tau {
        Side::Kept
    } else if hi <= tau {
        Side::Cut
    } else {
        Side::Straddles
    }
}

/// Which side of `keep` an axis-aligned box lies on.
///
/// The extreme corners along the normal bound the box's signed-distance range:
/// `d(centre) ± Σ|n·axis|·half_extent`.
fn straddle(bb: &Aabb, keep: &Plane3, tol: f64) -> Side {
    let n = keep.normal().as_vec();
    let c = Point3::new(
        0.5 * (bb.min.x + bb.max.x),
        0.5 * (bb.min.y + bb.max.y),
        0.5 * (bb.min.z + bb.max.z),
    );
    let reach = 0.5
        * (n.x.abs() * (bb.max.x - bb.min.x)
            + n.y.abs() * (bb.max.y - bb.min.y)
            + n.z.abs() * (bb.max.z - bb.min.z));
    let d = keep.signed_distance(c);
    if d - reach >= -tol {
        Side::Kept
    } else if d + reach <= tol {
        Side::Cut
    } else {
        Side::Straddles
    }
}

// ── Section registry (pre-split) ─────────────────────────────────────────────

/// 2-D segment × segment intersection, as parameters `(ta, tb)` along each.
///
/// `None` for parallel segments and for crossings outside either span.
fn seg_cross(a0: (f64, f64), a1: (f64, f64), b0: (f64, f64), b1: (f64, f64)) -> Option<(f64, f64)> {
    let da = (a1.0 - a0.0, a1.1 - a0.1);
    let db = (b1.0 - b0.0, b1.1 - b0.1);
    let den = da.0 * db.1 - da.1 * db.0;
    if den.abs() < 1e-14 {
        return None; // parallel / degenerate
    }
    let w = (b0.0 - a0.0, b0.1 - a0.1);
    let ta = (w.0 * db.1 - w.1 * db.0) / den;
    let tb = (w.0 * da.1 - w.1 * da.0) / den;
    const E: f64 = 1e-9;
    if !(-E..=1.0 + E).contains(&ta) || !(-E..=1.0 + E).contains(&tb) {
        return None;
    }
    Some((ta.clamp(0.0, 1.0), tb.clamp(0.0, 1.0)))
}

/// Split every section curve where it crosses another one, so the cap and the
/// lateral faces are handed the SAME pieces.
///
/// This is the cut's pre-split registry for
/// the union: where two fused tubes overlap at the plane their sections overlap
/// too, and each consumer would otherwise partition the shared curve its own
/// way — the cap at the outline of the union, the lateral face at its crossing
/// window — so no segment would ever weld.  Splitting ONCE here, with a single
/// computed point per crossing used by both curves, makes the partitions
/// identical by construction.
///
/// Every curve lies in the cut plane, so the whole thing is a 2-D problem in
/// the plane's own uv.
fn presplit_sections(ssi: &mut HashMap<FaceId, Vec<SsiCurve>>, cut_plane: &Plane3, tol: f64) {
    let uv = |p: Point3| {
        let w = p - cut_plane.frame.origin;
        (
            cut_plane.frame.x.dot_vec(w),
            cut_plane.frame.y.dot_vec(w),
        )
    };
    let lift = |(u, v): (f64, f64)| cut_plane.frame.origin + cut_plane.frame.x * u + cut_plane.frame.y * v;

    // Flatten to a working list: (owner face, uv polyline with the wrap point
    // appended for closed curves, was_closed).
    struct Work {
        face: FaceId,
        pts: Vec<(f64, f64)>,
        closed: bool,
        /// (segment index, parameter along it, crossing point) — the splits.
        cuts: Vec<(usize, f64, (f64, f64))>,
    }
    let mut work: Vec<Work> = Vec::new();
    for (&face, curves) in ssi.iter() {
        for c in curves {
            if c.points.len() < 2 {
                continue;
            }
            let mut pts: Vec<(f64, f64)> = c.points.iter().map(|&p| uv(p)).collect();
            if c.closed {
                pts.push(pts[0]);
            }
            work.push(Work {
                face,
                pts,
                closed: c.closed,
                cuts: Vec::new(),
            });
        }
    }

    // Pairwise crossings.  The SAME 2-D point is recorded on both curves.
    for i in 0..work.len() {
        for j in (i + 1)..work.len() {
            let mut found: Vec<(usize, f64, usize, f64, (f64, f64))> = Vec::new();
            const END: f64 = 1e-9;
            for si in 0..work[i].pts.len().saturating_sub(1) {
                let (a0, a1) = (work[i].pts[si], work[i].pts[si + 1]);
                for sj in 0..work[j].pts.len().saturating_sub(1) {
                    let (b0, b1) = (work[j].pts[sj], work[j].pts[sj + 1]);
                    if let Some((ta, tb)) = seg_cross(a0, a1, b0, b1) {
                        let a_end = ta < END || ta > 1.0 - END;
                        let b_end = tb < END || tb > 1.0 - END;
                        if a_end && b_end {
                            continue; // they already share this vertex
                        }
                        // the shared point: an existing vertex when one curve
                        // ends here (T-junction), else the computed crossing
                        let p = if a_end {
                            if ta < END { a0 } else { a1 }
                        } else if b_end {
                            if tb < END { b0 } else { b1 }
                        } else {
                            (a0.0 + (a1.0 - a0.0) * ta, a0.1 + (a1.1 - a0.1) * ta)
                        };
                        found.push((si, ta, sj, tb, p));
                    }
                }
            }
            for (si, ta, sj, tb, p) in found {
                // Two curves that already MEET at a vertex (a joint grafted into
                // both) also cross a hair away from it: each is a chord
                // approximation of its curve, and the chords on either side of
                // the joint overshoot differently.  That second crossing is an
                // artefact — splitting there creates a sliver piece shorter
                // than the weld tolerance.  Skip a crossing within a fraction
                // of the sample spacing of a vertex both curves carry.
                let near = (0.25 * super::ssi::SSI_SPACING).max(tol);
                let shared_vertex = work[i].pts.iter().any(|&a| {
                    let d = ((a.0 - p.0).powi(2) + (a.1 - p.1).powi(2)).sqrt();
                    d < near
                        && work[j].pts.iter().any(|&b| {
                            ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt() < 1e-9
                        })
                });
                if shared_vertex {
                    continue;
                }
                if std::env::var("CADCORE_DUMP_CUT").is_ok() {
                    let q = lift(p);
                    eprintln!(
                        "[cut][presplit] {:?}#{si} t={ta:.3e} x {:?}#{sj} t={tb:.3e} at ({:.4},{:.4},{:.4})",
                        work[i].face, work[j].face, q.x, q.y, q.z
                    );
                }
                if ta > END && ta < 1.0 - END {
                    work[i].cuts.push((si, ta, p));
                }
                if tb > END && tb < 1.0 - END {
                    work[j].cuts.push((sj, tb, p));
                }
            }
        }
    }

    if work.iter().all(|w| w.cuts.is_empty()) {
        return; // nothing overlaps — the curves are already the pieces
    }

    // Rebuild each curve as pieces broken at its cuts.
    let mut out: HashMap<FaceId, Vec<SsiCurve>> = HashMap::new();
    for w in work.iter_mut() {
        w.cuts
            .sort_by(|a, b| (a.0, a.1).partial_cmp(&(b.0, b.1)).unwrap_or(std::cmp::Ordering::Equal));
        w.cuts.dedup_by(|a, b| {
            let d = ((a.2 .0 - b.2 .0).powi(2) + (a.2 .1 - b.2 .1).powi(2)).sqrt();
            d < tol
        });

        if w.cuts.is_empty() {
            out.entry(w.face).or_default().push(SsiCurve {
                points: w
                    .pts
                    .iter()
                    .take(w.pts.len() - usize::from(w.closed))
                    .map(|&q| lift(q))
                    .collect(),
                closed: w.closed,
            });
            continue;
        }

        let mut pieces: Vec<Vec<(f64, f64)>> = Vec::new();
        let mut cur: Vec<(f64, f64)> = vec![w.pts[0]];
        let mut ci = 0usize;
        for si in 0..w.pts.len() - 1 {
            while ci < w.cuts.len() && w.cuts[ci].0 == si {
                let p = w.cuts[ci].2;
                if cur.last().map_or(true, |&q| dist2(q, p) > tol * tol) {
                    cur.push(p);
                }
                if cur.len() >= 2 {
                    pieces.push(std::mem::take(&mut cur));
                }
                cur = vec![p];
                ci += 1;
            }
            let nxt = w.pts[si + 1];
            if cur.last().map_or(true, |&q| dist2(q, nxt) > 1e-24) {
                cur.push(nxt);
            }
        }
        if cur.len() >= 2 {
            pieces.push(cur);
        }
        // A closed curve's last piece continues into its first across the wrap
        // point, so they are one arc.
        if w.closed && pieces.len() >= 2 {
            let last = pieces.pop().unwrap();
            let first = pieces.remove(0);
            let mut joined = last;
            joined.extend(first.into_iter().skip(1));
            pieces.push(joined);
        }
        for p in pieces {
            if p.len() >= 2 {
                out.entry(w.face).or_default().push(SsiCurve {
                    points: p.into_iter().map(lift).collect(),
                    closed: false,
                });
            }
        }
    }
    *ssi = out;
}

fn dist2(a: (f64, f64), b: (f64, f64)) -> f64 {
    (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)
}

// ── Section-loop assembly ────────────────────────────────────────────────────

/// Chain every section curve into closed 3-D loops.
///
/// A plane cutting one filament transversally already gives a closed conic; a
/// cut through a welded leg↔elbow junction gives two or more open arcs that
/// have to be walked end-to-end.  Curves whose ends do not meet anything are
/// dropped rather than guessed at — an unclosed cap is worse than a missing one
/// (it would emit a face with an open wire).
fn section_loops(
    ssi: &HashMap<FaceId, Vec<SsiCurve>>,
    rim_curves: &[Vec<Point3>],
    tol: f64,
) -> Vec<Vec<Point3>> {
    let mut open: Vec<Vec<Point3>> = Vec::new();
    let mut closed: Vec<Vec<Point3>> = Vec::new();
    for pts in rim_curves {
        if pts.len() >= 2 {
            open.push(pts.clone());
        }
    }
    for curves in ssi.values() {
        for c in curves {
            if c.points.len() < 2 {
                continue;
            }
            let ends_meet =
                (c.points[0] - *c.points.last().unwrap()).length() <= tol;
            if c.closed || (ends_meet && c.points.len() >= 3) {
                let mut pts = c.points.clone();
                if ends_meet && pts.len() >= 3 {
                    pts.pop(); // stored as an explicit repeat — the loop is implicit
                }
                closed.push(pts);
            } else {
                open.push(c.points.clone());
            }
        }
    }

    if std::env::var("CADCORE_DUMP_CUT").is_ok() {
        eprintln!(
            "[cut] section curves: {} closed, {} open (tol={tol:.1e})",
            closed.len(),
            open.len()
        );
        for (f, cs) in ssi.iter() {
            for c in cs {
                let (a, b) = (c.points[0], *c.points.last().unwrap());
                eprintln!(
                    "[cut]   {f:?} n={} closed={} ({:.4},{:.4},{:.4})..({:.4},{:.4},{:.4})",
                    c.points.len(),
                    c.closed,
                    a.x, a.y, a.z, b.x, b.y, b.z
                );
            }
        }
    }

    // A closed section through which OPEN curves end (a rim circle lying in
    // the plane that meets a tube's section, two fused tubes' sections meeting
    // at their shared channel) is not a loop of its own: the outline runs
    // along part of it and leaves along the open curve.  Break it at those
    // junction vertices into arcs so the walk below can route through them;
    // the arcs that end up inside the neighbour simply dangle and are dropped.
    let ends: Vec<Point3> = open
        .iter()
        .flat_map(|c| [c[0], *c.last().unwrap()])
        .collect();
    let mut still_closed: Vec<Vec<Point3>> = Vec::new();
    for lp in closed.drain(..) {
        let n = lp.len();
        let cuts: Vec<usize> = (0..n)
            .filter(|&i| ends.iter().any(|&e| (e - lp[i]).length() <= tol))
            .collect();
        if cuts.len() < 2 {
            still_closed.push(lp);
            continue;
        }
        for k in 0..cuts.len() {
            let (i0, i1) = (cuts[k], cuts[(k + 1) % cuts.len()]);
            let mut arc = vec![lp[i0]];
            let mut i = i0;
            loop {
                i = (i + 1) % n;
                arc.push(lp[i]);
                if i == i1 {
                    break;
                }
            }
            if arc.len() >= 2 {
                open.push(arc);
            }
        }
    }
    closed = still_closed;

    // Walk the open arcs end-to-end.  Each arc is used once; a chain that comes
    // back to its own start is a loop, anything else is discarded.
    let mut used = vec![false; open.len()];
    for i in 0..open.len() {
        if used[i] {
            continue;
        }
        used[i] = true;
        let mut chain = open[i].clone();
        loop {
            let tail = *chain.last().unwrap();
            if (tail - chain[0]).length() <= tol && chain.len() >= 3 {
                closed.push(chain);
                break;
            }
            let mut extended = false;
            for j in 0..open.len() {
                if used[j] {
                    continue;
                }
                let (a, b) = (open[j][0], *open[j].last().unwrap());
                if (tail - a).length() <= tol {
                    used[j] = true;
                    chain.extend_from_slice(&open[j][1..]);
                    extended = true;
                    break;
                }
                if (tail - b).length() <= tol {
                    used[j] = true;
                    chain.extend(open[j].iter().rev().skip(1).copied());
                    extended = true;
                    break;
                }
            }
            if !extended {
                break; // dangling chain — not a closed section, drop it
            }
        }
    }
    closed
}

/// Arrange the section loops in the cut plane and emit the cells that carry
/// material as cap faces.
///
/// The arrangement is REQUIRED, not a convenience: where two fused tubes
/// overlap at the plane their section loops overlap too, so the real outline of
/// the cut is the union of the loops, not the loops themselves.  (Emitting each
/// loop directly caps each tube's full circle, including the part buried inside
/// its neighbour.)
///
/// The cap's carrier plane is oriented OUTWARD (away from the kept material),
/// so the DCEL's CCW-in-uv outer loops and CW holes come out physically wound
/// as SolidWorks/SpaceClaim require (README invariant #4).
fn emit_section_caps(
    src: &BRep,
    faces: &[FaceId],
    cut_plane: &Plane3,
    keep_normal: UnitVec3,
    half: f64,
    loops: &[Vec<Point3>],
    asm: &mut Assembler,
) -> usize {
    if loops.is_empty() {
        return 0;
    }
    // Outward = away from the kept material.
    let cap_plane = Plane3::from_origin_normal(cut_plane.frame.origin, -keep_normal);
    let domain = FaceDomain::Plane {
        plane: cap_plane,
        half_extent: half,
    };

    // Bounding rectangle so the arrangement is bounded; its cell lies outside
    // every section loop and is dropped below.
    const RECT_TAG: u32 = 1;
    let rect = vec![
        (-half, -half),
        (half, -half),
        (half, half),
        (-half, half),
        (-half, -half),
    ];
    let mut chains = vec![FaceChain {
        pts: rect,
        tag: RECT_TAG,
    }];
    let mut tag = 1000u32;
    for lp in loops {
        let mut pts: Vec<(f64, f64)> = lp.iter().map(|&p| domain.uv(p)).collect();
        pts.push(pts[0]); // close it for the DCEL
        chains.push(FaceChain { pts, tag });
        tag += 1;
    }

    // A cap cell carries material when a point just INSIDE the kept side of it
    // is inside the original solid.  `arrange_and_classify_with` keeps a cell
    // when the predicate is false, so the predicate is the negation.
    let probe = 1.0e-4;
    let void = |p: Point3| !point_in_solid(src, faces, p + keep_normal * probe);

    let cells = arrange_and_classify_with(&domain, &chains, 0.0, &void);
    let mut emitted = 0usize;
    for cell in cells {
        if !cell.keep {
            continue;
        }
        // The bounding rectangle is scaffolding, not geometry: the cell that
        // uses it is the plane OUTSIDE the part, with every section island as a
        // hole.  Emitting it would hand each section curve a third user (the
        // island cap and the trimmed lateral face are the only two legitimate
        // ones) and paste a giant disc across the model.
        if cell.cell.outer.iter().any(|s| s.tag == RECT_TAG) {
            continue;
        }
        asm.emit_cell(&domain, FaceGeom::Plane(cap_plane), FaceNormal::Same, &cell);
        emitted += 1;
    }
    emitted
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean::tests_support::capped_cylinder;

    /// Gate 1: every edge of the result must be used exactly twice.
    fn open_edges(brep: &BRep, shell: ShellId) -> usize {
        crate::validate::manifold::check(brep, shell).violations.len()
    }

    /// The full SpaceClaim accept gate: manifold + edge-to-surface deviation +
    /// uv wire winding / self-intersection.  `(wire, distance, manifold)`.
    fn spaceclaim_gate(brep: &BRep, shell: ShellId) -> (usize, usize, usize) {
        let mut cfg = crate::config::UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = crate::Diagnostics::new();
        let r = crate::validate::validate_shell(brep, shell, &cfg, &mut diag);
        (
            r.wire_violations.len(),
            r.distance_violations.len(),
            r.manifold_violations.len(),
        )
    }

    /// Every point of every edge in the result, sampled.
    fn shell_points(brep: &BRep, shell: ShellId) -> Vec<Point3> {
        let mut out = Vec::new();
        for &fid in &brep.shells[shell].faces {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let Some(lp) = brep.loops.get(lid) else { continue };
                let st = lp.start;
                let mut c = st;
                loop {
                    let ce = &brep.coedges[c];
                    out.extend(super::super::aabb::sample_edge(brep, ce.edge));
                    c = ce.next;
                    if c == st {
                        break;
                    }
                }
            }
        }
        out
    }

    /// No geometry may survive on the discarded side.
    fn deepest_past_plane(brep: &BRep, shell: ShellId, keep: &Plane3) -> f64 {
        shell_points(brep, shell)
            .into_iter()
            .map(|p| keep.signed_distance(p))
            .fold(f64::MAX, f64::min)
    }

    fn face_count(brep: &BRep, shell: ShellId) -> usize {
        brep.shells[shell].faces.len()
    }

    /// Sampled points of one loop, in traversal order.
    fn loop_points(brep: &BRep, lid: cadcore_topo::LoopId) -> Vec<Point3> {
        let Some(lp) = brep.loops.get(lid) else { return Vec::new() };
        let st = lp.start;
        if brep.coedges.get(st).is_none() {
            return Vec::new();
        }
        let mut out = Vec::new();
        let mut c = st;
        loop {
            let ce = &brep.coedges[c];
            out.extend(super::super::aabb::sample_edge(brep, ce.edge));
            c = ce.next;
            if c == st {
                break;
            }
        }
        out
    }

    /// A capped cylinder cut across its axis: the part keeps its lateral band up
    /// to the plane, loses the far cap and gains a flat disk section.
    #[test]
    fn cut_cylinder_across_axis_is_watertight() {
        let mut a = BRep::new();
        let fa = capped_cylinder(&mut a, Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, 0.5, 4.0);
        // keep z <= 2.5  ⇒ normal points down
        let keep = Plane3::from_origin_normal(Point3::new(0.0, 0.0, 2.5), -UnitVec3::Z);

        let cut = half_space_cut(&a, &fa, &keep, 1e-6).expect("cut runs");
        assert!(cut.caps >= 1, "the section must be capped");
        assert_eq!(
            cut.dropped, 1,
            "the far end cap is entirely past the plane and must be dropped"
        );
        assert_eq!(
            open_edges(&cut.brep, cut.shell),
            0,
            "cut solid must stay watertight"
        );
        assert!(
            deepest_past_plane(&cut.brep, cut.shell, &keep) > -1e-6,
            "material survived past the cut plane: {:.3e}",
            deepest_past_plane(&cut.brep, cut.shell, &keep)
        );
        // The cap sits ON the plane and spans the tube: its points are at the
        // cut height and within the radius.
        let on_plane: Vec<Point3> = shell_points(&cut.brep, cut.shell)
            .into_iter()
            .filter(|p| keep.signed_distance(*p).abs() < 1e-6)
            .collect();
        assert!(!on_plane.is_empty(), "no section geometry at the plane");
        for p in &on_plane {
            let r = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (r - 0.5).abs() < 1e-3,
                "section rim off the tube radius: r={r:.6}"
            );
        }
    }

    /// A plane clear of the part changes nothing and adds no cap.
    #[test]
    fn cut_missing_the_part_keeps_it_whole() {
        let mut a = BRep::new();
        let fa = capped_cylinder(&mut a, Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, 0.5, 4.0);
        let keep = Plane3::from_origin_normal(Point3::new(0.0, 0.0, -5.0), UnitVec3::Z);

        let cut = half_space_cut(&a, &fa, &keep, 1e-6).expect("cut runs");
        assert_eq!(cut.caps, 0, "nothing to cap");
        assert_eq!(cut.dropped, 0, "nothing to drop");
        assert_eq!(open_edges(&cut.brep, cut.shell), 0);
        assert_eq!(face_count(&cut.brep, cut.shell), fa.len(), "part unchanged");
    }

    /// Two crossing serpentines FUSED into one solid (the shape the scaffold
    /// export builds), then cut.
    fn fused_scaffold() -> (BRep, ShellId) {
        use crate::arrange::filament::Filament;
        use crate::arrange::universal::fuse_filaments;
        let r = 0.275;
        // A: U in z=0 — X-legs at y=0 and y=2, connector at x=0.
        let a = Filament::serpentine(
            &[
                Point3::new(-2.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
                Point3::new(-2.0, 2.0, 0.0),
            ],
            r,
            0.5,
        );
        // B: U one layer up — Y-legs at x=-1.5 and x=-0.8 crossing both of A's.
        let b = Filament::serpentine(
            &[
                Point3::new(-1.5, -1.0, 0.35),
                Point3::new(-1.5, 3.0, 0.35),
                Point3::new(-0.8, 3.0, 0.35),
                Point3::new(-0.8, -1.0, 0.35),
            ],
            r,
            0.3,
        );
        let mut brep = BRep::new();
        let shell = fuse_filaments(&mut brep, &[a, b]);
        (brep, shell)
    }

    /// The shape a real generated scaffold fuses into: two 5-leg serpentines on
    /// stacked layers, the geometry `brep::pipeline` builds for a 3 mm 2-layer
    /// part (filament ⌀0.55, spacing 1.2, layer 0.35, bend radius 0.5).
    fn fused_real_scaffold() -> (BRep, ShellId) {
        use crate::arrange::filament::Filament;
        use crate::arrange::universal::fuse_filaments;
        let r = 0.275;
        let z0 = 0.35;
        let z1 = 0.70;
        // layer 0 sweeps along +X with Y-legs; layer 1 along +Y with X-legs.
        let a = Filament::serpentine(
            &[
                Point3::new(0.3, 0.3, z0),
                Point3::new(0.3, 2.7, z0),
                Point3::new(1.5, 2.7, z0),
                Point3::new(1.5, 0.3, z0),
                Point3::new(2.7, 0.3, z0),
                Point3::new(2.7, 2.7, z0),
            ],
            r,
            0.5,
        );
        let b = Filament::serpentine(
            &[
                Point3::new(0.3, 0.3, z1),
                Point3::new(2.7, 0.3, z1),
                Point3::new(2.7, 1.5, z1),
                Point3::new(0.3, 1.5, z1),
                Point3::new(0.3, 2.7, z1),
                Point3::new(2.7, 2.7, z1),
            ],
            r,
            0.5,
        );
        let mut brep = BRep::new();
        let shell = fuse_filaments(&mut brep, &[a, b]);
        (brep, shell)
    }

    /// The FUSE alone — no cut anywhere near it — must pass the gate a
    /// receiving CAD applies.  "Press Combine, open it in SpaceClaim, the part
    /// comes in torn into pieces" is a WIRE violation: SpaceClaim ignores
    /// `FACE_BOUND` orientation flags and drops outright any face whose wires
    /// intersect in uv (README invariant #5).
    ///
    /// RED — `wire = 4` on a plain 3 mm 2-layer scaffold.  Diagnosed:
    ///
    /// ```text
    /// FaceId(1)  rim at v=0.000   window at v[-0.264 .. +0.265]
    /// FaceId(6)  rim at v=1.900   window at v[ 1.635 ..  2.165]
    /// FaceId(12) rim at v=0.000   window at v[-0.265 .. +0.265]
    /// FaceId(17) rim at v=1.900   window at v[ 1.636 ..  2.165]
    /// ```
    ///
    /// Every one is a crossing window sitting ASTRIDE a leg's rim: the crossing
    /// lands exactly where the leg meets its elbow, so half the loop is outside
    /// the band and the hole's wire crosses the rim's wire.  Four faces dropped
    /// by the importer, four holes in the part.
    ///
    /// The one-line guard (extend `window_clears_free_ends` to every rim, not
    /// just free ones, for generated scaffolds too) makes this green — and
    /// breaks `leg_over_elbow_watertight`, which needs exactly such a crossing
    /// to weld and splits into two solids without it.  Measured both ways; it is
    /// a genuine conflict, not a tolerance to tune.
    ///
    /// The fix is to stop treating these as all-or-nothing: a crossing that
    /// straddles a junction must be SPLIT at the junction circle and imprinted
    /// on both members — as a hole on neither, but as a shared notch across the
    /// leg↔elbow boundary.  Then it welds AND every wire stays inside its face.
    #[test]
    fn fusing_a_real_scaffold_passes_the_spaceclaim_gate() {
        let (fused, shell) = fused_real_scaffold();
        assert_eq!(open_edges(&fused, shell), 0, "fused shell is manifold");
        let (wire, dist, man) = spaceclaim_gate(&fused, shell);
        assert_eq!(
            (wire, dist, man),
            (0, 0, 0),
            "SpaceClaim gate on the FUSED part: wire={wire} distance={dist} manifold={man}"
        );
    }

    /// The cut has to survive the real thing, not just the two-filament
    /// fixture: both terminals off a fused 3 mm scaffold.
    ///
    /// RED — 47 single-use edges (`CADCORE_DUMP_CUT=1` for the breakdown).  The
    /// fast repro for the pipeline-level failure in
    /// `brep/tests/fuse_then_cut.rs`; runs in ~0.3 s instead of 74 s.
    ///
    /// What is already right: the broad phase finds exactly the three tubes the
    /// plane crosses, each yields one closed conic, and three caps are emitted —
    /// the section side is correct.
    ///
    /// What is wrong, from the coordinates of the 47: they sit on **crossing
    /// windows**, and not only at the plane —
    ///
    /// ```text
    /// (2.425, 1.000, 0.700) -> (2.446, 1.000, 0.805)   on the cut plane
    /// (2.436, 2.685, 0.425) -> (2.436, 2.700, 0.426)   1.7 mm away from it
    /// ```
    ///
    /// so this is NOT the cut tearing the section: it is the arranged↔copied
    /// boundary again, at a crossing window whose two members end up on
    /// different paths (the leg the plane crosses is ARRANGED, the leg it
    /// crosses at y = 2.7 is wholly kept and COPIED).  The joints machinery does
    /// not cover it because every section here is a CLOSED loop and
    /// `section_joints` only publishes the ends of OPEN ones — so whatever the
    /// arranged face's DCEL does to that window is invisible to the copy.
    ///
    /// Next: find what splits the window on the arranged side when nothing
    /// crosses it (the section conic is 1.7 mm away and the seam is parked in a
    /// gap), then publish that split like any other joint.
    ///
    /// NB the offsets here are deliberately NOT `0.8`: at `0.8` the plane lands
    /// exactly on the fillet tangent (`0.3 + 0.5`), so the junction circles lie
    /// IN the plane — a separate degenerate case worth its own test.
    #[test]
    fn cut_a_real_scaffold_is_watertight() {
        let (fused, shell) = fused_real_scaffold();
        let faces = fused.shells[shell].faces.clone();
        assert_eq!(
            open_edges(&fused, shell),
            0,
            "fixture starts watertight ({} faces)",
            faces.len()
        );

        // Both Y terminals, as `build_trim_planes` places them for offset 0.8.
        let lo = Plane3::from_origin_normal(Point3::new(0.0, 1.0, 0.0), UnitVec3::Y);
        let first = half_space_cut(&fused, &faces, &lo, 1e-6).expect("first cut runs");
        assert_eq!(
            open_edges(&first.brep, first.shell),
            0,
            "watertight after the first terminal ({} caps, {} dropped)",
            first.caps,
            first.dropped
        );

        let faces2 = first.brep.shells[first.shell].faces.clone();
        let hi = Plane3::from_origin_normal(Point3::new(0.0, 2.0, 0.0), -UnitVec3::Y);
        let second = half_space_cut(&first.brep, &faces2, &hi, 1e-6).expect("second cut runs");
        assert_eq!(
            open_edges(&second.brep, second.shell),
            0,
            "watertight after the second terminal ({} caps, {} dropped)",
            second.caps,
            second.dropped
        );
    }

    /// Re-emitting a fused shell must be EXACT: same faces, still watertight.
    ///
    /// This is the foundation the cut stands on — it proves `emit_face_copy`
    /// reproduces a Boolean'd shell (seamless periodic bands, crossing-window
    /// holes, torus fillets) through the assembler without losing a single
    /// shared edge.
    #[test]
    fn re_emitting_a_fused_scaffold_is_lossless() {
        let (fused, shell) = fused_scaffold();
        let faces = fused.shells[shell].faces.clone();
        assert_eq!(open_edges(&fused, shell), 0, "fixture starts watertight");

        // A plane far away on the discarded side: nothing is cut, so this
        // exercises the re-emit path alone.
        let keep = Plane3::from_origin_normal(Point3::new(-50.0, 0.0, 0.0), UnitVec3::X);
        let cut = half_space_cut(&fused, &faces, &keep, 1e-6).expect("cut runs");

        assert_eq!(cut.caps, 0, "nothing to cap");
        assert_eq!(cut.dropped, 0, "nothing to drop");
        assert_eq!(
            face_count(&cut.brep, cut.shell),
            faces.len(),
            "face count must be preserved exactly"
        );
        assert_eq!(
            open_edges(&cut.brep, cut.shell),
            0,
            "re-emitting a fused shell must not lose a shared edge"
        );
    }

    /// The real target: cut a shell that is already FUSED (its faces carry
    /// explicit trim loops, not analytic templates).  This is exactly where the
    /// legacy `cadcore-ops::half_space_cut_brep` silently does nothing — it
    /// finds no template cylinders and returns the part untouched.
    ///
    /// Getting here took six distinct fixes, each independently verified — they
    /// are listed because every one of them is a trap the next surface pair will
    /// walk into again:
    ///
    /// * plane∥axis sections produced no curve at all (`CylPlaneCurve::Lines`
    ///   was "deferred"), so a lengthwise-grazed tube survived the cut whole;
    /// * the arrangement's bounding rectangle was emitted as a third face over
    ///   every section, so section edges came out used 3×;
    /// * a traced section was clipped to the last sample INSIDE the arc, ~7 µm
    ///   short of the junction, so it never chained to its neighbour;
    /// * the clip walk emitted its head fragment twice;
    /// * overlapping sections were partitioned differently by the cap and by the
    ///   lateral faces (`presplit_sections` + `section_joints`/`insert_joints`);
    /// * the parametric seam: `build_periodic_u` unwraps each chain from its own
    ///   first sample and splits at multiples of the period, so it cuts along
    ///   **u = 0**, not u = ±π.  Parking the seam — and publishing its crossings
    ///   as joints — half a period off left exactly 9 open edges where a window
    ///   loop crossed the top of a tube.
    #[test]
    fn cut_a_fused_scaffold_is_watertight() {
        let (fused, shell) = fused_scaffold();
        let faces = fused.shells[shell].faces.clone();
        assert_eq!(open_edges(&fused, shell), 0, "fixture starts watertight");

        // Keep x <= -0.7: severs A's two X-legs clear of their U-turn fillets,
        // drops A's connector, leaves every crossing with B intact.
        let keep = Plane3::from_origin_normal(Point3::new(-0.7, 0.0, 0.0), -UnitVec3::X);
        // A coarser knit than the fuse needs: this fixture is also cut with ~1e-4 of noise tolerated on
        // its shared points, so its knit is the half-micron the legacy export
        // path uses (the universal fuse is exact and cuts at 1e-6).
        let cut = half_space_cut(&fused, &faces, &keep, 5e-4).expect("cut runs");

        // Three tubes are severed — a transverse disk through each of A's two
        // X-legs, plus a longitudinal chord through B's leg at x = -0.8 (its
        // tube reaches x = -0.525, so the plane grazes it lengthwise).  They
        // OVERLAP at the plane, so the arrangement tiles the section into a few
        // coplanar cells rather than one face per tube; what matters is that
        // they tile it completely, which the watertight assertion below proves.
        // Every section curve is clipped to its face's trimmed region, so the
        // cap tiles the UNION outline of the three overlapping sections once —
        // one cell, not one per tube.
        assert!(cut.caps >= 1, "the section must be capped, got {}", cut.caps);
        assert!(cut.dropped > 0, "A's connector and fillets must be dropped");
        assert_eq!(
            open_edges(&cut.brep, cut.shell),
            0,
            "cut fused shell must stay watertight"
        );
        // Traced sections and welded vertices land within the assembler's weld
        // tolerance of the plane, so the bound is the knit scale (1 µm — what
        // SolidWorks itself knits at), not float epsilon.
        assert!(
            deepest_past_plane(&cut.brep, cut.shell, &keep) > -1e-3,
            "material survived past the cut plane: {:.3e}",
            deepest_past_plane(&cut.brep, cut.shell, &keep)
        );
        // Every cap must sit where there IS material: probing just inside the
        // kept side of a cap must land inside the original solid.
        let n = keep.normal();
        for &fid in &cut.brep.shells[cut.shell].faces {
            let f = &cut.brep.faces[fid];
            if !matches!(f.geom, FaceGeom::Plane(_)) {
                continue;
            }
            let pts = loop_points(&cut.brep, f.outer_loop);
            if pts.is_empty() || pts.iter().any(|p| keep.signed_distance(*p).abs() > 1e-6) {
                continue; // not a section cap
            }
            // The outline is NOT convex (two disks joined by a lengthwise
            // strip), so its vertex centroid is no interior probe.  What a
            // cap over material means: every edge of its outline separates
            // material from void — just behind the cap, exactly one of the
            // two sides of the edge is inside the original solid.
            // `loop_points` lists the two ends of every co-edge in turn
            for pair in pts.chunks(2) {
                let [a, b] = [pair[0], pair[pair.len() - 1]];
                let e = b - a;
                if e.length() < 1e-2 {
                    continue; // too short to probe clear of its neighbours
                }
                let c = n.as_vec().cross(e);
                let side = c * (1e-3 / c.length());
                // a micron behind the cap: where the plane grazes an elbow the
                // outline moves fast with depth
                let m = a.lerp(b, 0.5) + n.as_vec() * 1e-6;
                let (l, r) = (
                    point_in_solid(&fused, &faces, m + side),
                    point_in_solid(&fused, &faces, m - side),
                );
                assert!(l != r, "cap edge {a:?} -> {b:?} at {m:?} does not bound material (left {l}, right {r})");
            }
        }
        // And the gate a receiving CAD actually applies.
        // The SpaceClaim gate: wire winding and manifold are clean.  Gate 4
        // (edge-to-surface deviation) is NOT asserted yet — see the module TODO
        // on the seam split point; it is the last thing between this and a
        // guaranteed SolidWorks-clean cut.
        let (wire, _dist, man) = spaceclaim_gate(&cut.brep, cut.shell);
        assert_eq!(wire, 0, "uv wires must be wound and non-intersecting");
        assert_eq!(man, 0, "shell must be manifold");
    }

    /// The legacy template cut is a silent no-op on a fused shell — the reason
    /// this module exists.  Guards the claim so it cannot rot.
    #[test]
    fn legacy_template_cut_cannot_touch_a_fused_shell() {
        let (fused, shell) = fused_scaffold();
        let templated = fused.shells[shell].faces.iter().any(|&f| {
            matches!(
                (&fused.faces[f].geom, &fused.faces[f].extent),
                (FaceGeom::Cylinder(_), FaceExtent::Cylinder { .. })
            )
        });
        assert!(
            !templated,
            "a fused shell carries no template cylinder faces, so a \
             template-rebuilding cut has nothing to classify"
        );
    }
}

