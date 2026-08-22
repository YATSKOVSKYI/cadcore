//! The scaffold orchestrator — fuse an arbitrary set of serpentine filaments
//! (each: straight legs + torus-fillet elbows) that cross each other into ONE
//! watertight, seamless solid.
//!
//! This is the universal entry point the app pipeline targets: every leg of
//! every filament is trimmed by every leg of every OTHER filament it
//! overlaps, junctions weld legs to elbows, free leg ends are capped, and the
//! whole thing is assembled through one [`Assembler`] (shared edges → one
//! solid).  All cylinders are seamless (analytic circle rims), all elbows are
//! seamless torus fillets.

use cadcore_geom::{Circle3, CylSurf, Plane3, TorusSurf};
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{BRep, Shell, ShellId, Solid};

use super::assembly::Assembler;
use super::filament::{ElbowSeg, LegSeg};
use crate::boolean::Aabb;
use crate::geom::composite::{closed_loops_between, CompositeTube};
use crate::geom::intersect::TraceOptions;

/// One serpentine filament path: `legs[i]` joined to `legs[i+1]` by
/// `elbows[i]`.  The two free ends (first leg incoming, last leg outgoing)
/// are capped.  Each leg is built with `surf.frame.origin` at its incoming
/// end and `surf.axis()` pointing toward the outgoing end.
#[derive(Debug, Clone)]
pub struct Filament {
    /// Straight legs in path order.
    pub legs: Vec<LegSeg>,
    /// Elbows (`legs.len() - 1` of them); `elbows[i]` joins legs i and i+1.
    pub elbows: Vec<ElbowSeg>,
}

impl Filament {
    /// Build a planar serpentine filament from a polyline path: straight legs
    /// joined by constant-radius circular fillets (one per interior vertex).
    /// `r` = filament radius, `fillet` = elbow bend radius.  Each leg runs
    /// `surf.frame.origin` → `+axis`, matching [`fuse_scaffold`]'s convention.
    pub fn serpentine(verts: &[Point3], r: f64, fillet: f64) -> Filament {
        let n = verts.len();
        assert!(n >= 2, "a filament needs at least two path vertices");
        let dir = |a: Point3, b: Point3| UnitVec3::try_from_vec(b - a).unwrap();
        let mut tin = vec![Point3::new(0.0, 0.0, 0.0); n];
        let mut tout = vec![Point3::new(0.0, 0.0, 0.0); n];
        let mut elbows = Vec::new();
        for k in 1..n.saturating_sub(1) {
            let din = dir(verts[k - 1], verts[k]);
            let dout = dir(verts[k], verts[k + 1]);
            let turn = din.dot_vec(dout.as_vec()).clamp(-1.0, 1.0).acos();
            let t = fillet / (turn / 2.0).tan();
            let pin = verts[k] - din.as_vec() * t;
            let pout = verts[k] + dout.as_vec() * t;
            let axis = UnitVec3::try_from_vec(din.as_vec().cross(dout.as_vec())).unwrap();
            let inward = UnitVec3::try_from_vec(axis.as_vec().cross(din.as_vec())).unwrap();
            let centre = pin + inward.as_vec() * fillet;
            let torus = TorusSurf::new(centre, axis, fillet, r);
            let theta = |p: Point3| {
                let w = p - torus.frame.origin;
                let rd = w - torus.frame.z.as_vec() * torus.frame.z.dot_vec(w);
                torus.frame.y.dot_vec(rd).atan2(torus.frame.x.dot_vec(rd))
            };
            let lo = theta(pin);
            let mut d = theta(pout) - lo;
            while d > std::f64::consts::PI { d -= std::f64::consts::TAU; }
            while d < -std::f64::consts::PI { d += std::f64::consts::TAU; }
            elbows.push(ElbowSeg { surf: torus, theta_lo: lo, theta_hi: lo + d });
            tin[k] = pin;
            tout[k] = pout;
        }
        let mut legs = Vec::new();
        for k in 0..n - 1 {
            let start = if k == 0 { verts[0] } else { tout[k] };
            let end = if k + 1 == n - 1 { verts[n - 1] } else { tin[k + 1] };
            legs.push(LegSeg {
                surf: CylSurf::new(start, dir(start, end), r),
                length: (end - start).length(),
            });
        }
        Filament { legs, elbows }
    }

    /// The composite tube (legs + elbows) for crossing / containment queries.
    fn tube(&self) -> CompositeTube {
        let mut t = CompositeTube::new();
        for l in &self.legs {
            t = t.with_leg(l.surf, l.length);
        }
        for e in &self.elbows {
            t = t.with_elbow(e.surf, e.theta_lo, e.theta_hi);
        }
        t
    }
}

/// Subsample a traced crossing loop to at most `max` points (keeping the first,
/// so the closed loop stays anchored).  The SSI tracer emits hundreds of points
/// per bite loop; a smooth crossing curve needs only a few dozen, and the
/// unfiltered points otherwise explode the STEP (each becomes a B-spline control
/// point + a dense pcurve sample → 100 KB/face).  The SAME decimated loop is
/// imprinted on both crossing members, so they stay welded.
fn decimate_loop(pts: &[Point3], max: usize) -> Vec<Point3> {
    let n = pts.len();
    if n <= max || max < 2 {
        return pts.to_vec();
    }
    let step = (n + max - 1) / max; // ceil(n / max)
    let mut out: Vec<Point3> = pts.iter().step_by(step).copied().collect();
    // Keep the true last point so the loop closes on the same span it traced.
    if let (Some(&last), Some(&outl)) = (pts.last(), out.last()) {
        if (last - outl).length() > 1e-9 {
            out.push(last);
        }
    }
    out
}

/// Centroid of a window loop (used to dedup coincident imprints).
fn win_centroid(p: &[Point3]) -> Point3 {
    let n = p.len().max(1) as f64;
    let mut acc = cadcore_math::Vec3::new(0.0, 0.0, 0.0);
    let o = Point3::new(0.0, 0.0, 0.0);
    for q in p {
        acc = acc + (*q - o);
    }
    o + acc * (1.0 / n)
}

/// Bounding radius of a window loop (max distance from its centroid).
fn win_radius(w: &[Point3], c: Point3) -> f64 {
    w.iter().map(|p| (*p - c).length()).fold(0.0, f64::max)
}

/// A window must clear any FREE (disk-capped) leg end by its bounding radius:
/// a crossing right at a short split-piece's capped end would imprint a hole
/// overlapping the rim circle → intersecting inner wires → SpaceClaim drops the
/// face.  `free_start`/`free_end` mark which ends carry a disk cap (a filament's
/// first-leg start and last-leg end); INTERIOR leg ends are elbow junctions
/// (welded, no cap), so a window there is fine and must NOT be dropped — else
/// tight scaffolds stop fusing into one solid.  A dropped crossing simply
/// interpenetrates (still watertight).
fn window_clears_free_ends(
    leg: &LegSeg,
    cen: Point3,
    rad: f64,
    free_start: bool,
    free_end: bool,
) -> bool {
    let ax = leg.surf.axis().dot_vec(cen - leg.surf.frame.origin);
    (!free_start || ax > rad) && (!free_end || ax < leg.length - rad)
}

/// True if the new window (centroid `c`, bounding radius `rad`) would overlap
/// any window already on the face: their bounding spheres meet, scaled by
/// `factor`.  Two overlapping window HOLES on one face make its inner wires
/// intersect in uv → SolidWorks/SpaceClaim drop the whole face (README #5).  A
/// dropped (too-close) crossing simply interpenetrates — both tubes stay closed
/// manifolds, so the result is still watertight, just not pre-welded there.
/// Sphere overlap (vs a fixed centroid distance) adapts to the elongated windows
/// of shallow crossings and the tightly-clustered crossings of dense solid fill.
fn overlaps_existing(existing: &[(u32, Vec<Point3>)], c: Point3, rad: f64, factor: f64) -> bool {
    existing.iter().any(|(_, w)| {
        let wc = win_centroid(w);
        (wc - c).length() < (win_radius(w, wc) + rad) * factor
    })
}

/// Every point of `pts` lies within the elbow's arc span (torus θ range) plus a
/// small margin — guards against a traced loop that wandered past the fillet's
/// actual arc onto phantom torus surface (which would imprint a hole where
/// there is no material).
fn points_within_arc(e: &ElbowSeg, pts: &[Point3]) -> bool {
    let t = &e.surf;
    let sweep = e.theta_hi - e.theta_lo;
    if sweep.abs() < 1e-9 {
        return false;
    }
    let margin = 0.15_f64; // radians
    pts.iter().all(|&p| {
        let w = p - t.frame.origin;
        let ang = t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w));
        let rel = (ang - e.theta_lo) * sweep.signum();
        rel > -margin && rel < sweep.abs() + margin
    })
}

/// Fuse a set of crossing serpentine filaments into one closed, seamless,
/// watertight solid.  Returns the shell.
///
/// Generated scaffolds fuse maximally (crossings kept everywhere) — their
/// woodpile crossings sit mid-leg and never overlap a rim, so this stays
/// SpaceClaim-clean AND one solid.
pub fn fuse_scaffold(brep: &mut BRep, filaments: &[Filament]) -> ShellId {
    fuse_scaffold_impl(brep, filaments, false)
}

/// Like [`fuse_scaffold`] but SpaceClaim-safe for IMPORTED parts: a crossing is
/// dropped (the tubes interpenetrate there) whenever its window would overlap a
/// leg-end rim circle — dense solid-fill caps are built from many short split
/// pieces whose crossings land right on the rims, and an overlapping hole makes
/// SpaceClaim drop the face ("torn caps").  The trade-off (a few interpenetrating
/// tubes vs a single fused solid) is invisible to a downstream Boolean/mesh.
pub fn fuse_scaffold_rimsafe(brep: &mut BRep, filaments: &[Filament]) -> ShellId {
    fuse_scaffold_impl(brep, filaments, true)
}

fn fuse_scaffold_impl(brep: &mut BRep, filaments: &[Filament], strict_rims: bool) -> ShellId {
    let nf = filaments.len();
    // Per-leg tube (one straight leg) for crossing traces, and the whole
    // filament tube for keep/drop containment.
    let fil_tubes: Vec<CompositeTube> = filaments.iter().map(|f| f.tube()).collect();
    let leg_tubes: Vec<Vec<CompositeTube>> = filaments
        .iter()
        .map(|f| {
            f.legs
                .iter()
                .map(|l| CompositeTube::new().with_leg(l.surf, l.length))
                .collect()
        })
        .collect();

    // Per-leg axis-aligned bounding box (segment endpoints, padded by the tube
    // radius) — the broad phase: two legs can only intersect if their boxes
    // overlap, so we skip the (expensive, fine-step) SSI trace for the vast
    // majority of non-crossing leg pairs.  Without this cull `fuse_scaffold`
    // is O(legs²) traces and does not scale to a real multi-layer scaffold
    // (hundreds of legs → tens of thousands of traces).
    let leg_box = |l: &LegSeg| -> Aabb {
        let a = l.surf.frame.origin;
        let b = a + l.surf.axis().as_vec() * l.length;
        let mut bb = Aabb::empty();
        bb.expand(a);
        bb.expand(b);
        bb.pad(l.surf.radius)
    };
    let leg_boxes: Vec<Vec<Aabb>> = filaments
        .iter()
        .map(|f| f.legs.iter().map(leg_box).collect())
        .collect();

    // Per-elbow single-member composite tubes + their arc bounding boxes, for
    // the elbow↔leg crossing pass (a U-turn fillet overlapping a neighbouring
    // leg — the dense-fill / hairpin case leg↔leg tracing alone cannot see).
    let elbow_tubes: Vec<Vec<CompositeTube>> = filaments
        .iter()
        .map(|f| {
            f.elbows
                .iter()
                .map(|e| CompositeTube::new().with_elbow(e.surf, e.theta_lo, e.theta_hi))
                .collect()
        })
        .collect();
    let elbow_box = |e: &ElbowSeg| -> Aabb {
        let mut bb = Aabb::empty();
        let steps = 8;
        for k in 0..=steps {
            let th = e.theta_lo + (e.theta_hi - e.theta_lo) * k as f64 / steps as f64;
            let dir = e.surf.frame.x * th.cos() + e.surf.frame.y * th.sin();
            bb.expand(e.surf.frame.origin + dir * e.surf.major_radius);
        }
        bb.pad(e.surf.minor_radius)
    };
    let elbow_boxes: Vec<Vec<Aabb>> = filaments
        .iter()
        .map(|f| f.elbows.iter().map(elbow_box).collect())
        .collect();

    // windows[fi][li] = crossing loops on leg (fi,li); elbow_windows[fi][ei] =
    // crossing loops on elbow (fi,ei).  Every crossing is traced ONCE and the
    // SAME 3-D loop is imprinted on BOTH crossing members — the single-source-
    // of-truth invariant that makes the shared edges weld with opposite senses
    // (README #2).  Imprinting symmetrically (both sides or neither) is what
    // keeps the shell manifold even when a marginal crossing is missed.
    // Each window carries a unique crossing id: the SAME id is imprinted on the
    // two crossing members so their pieces weld (single source of truth), while
    // two DISTINCT crossings that merely coincide in space (a self-overlapping
    // serpentine → a leg crossed by two coincident legs) get different ids and
    // stay separate twice-used edges instead of one four-used edge.
    let mut windows: Vec<Vec<Vec<(u32, Vec<Point3>)>>> = filaments
        .iter()
        .map(|f| vec![Vec::new(); f.legs.len()])
        .collect();
    let mut elbow_windows: Vec<Vec<Vec<(u32, Vec<Point3>)>>> = filaments
        .iter()
        .map(|f| vec![Vec::new(); f.elbows.len()])
        .collect();
    let opts = TraceOptions::default();
    let mut win_id: u32 = 0;

    let r_fil = filaments
        .iter()
        .flat_map(|f| f.legs.iter())
        .map(|l| l.surf.radius)
        .next()
        .unwrap_or(0.5);
    let min_pts = 12usize;
    // Max control points kept per crossing window (decimated from the dense SSI
    // trace) — caps the STEP size / write time.
    let win_max_pts: usize = std::env::var("CADCORE_WIN_MAX")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(40);
    // Coincident-window spacing: a face must never receive two OVERLAPPING holes
    // — SolidWorks drops a face whose inner wires intersect (README #5) — so a
    // window is imprinted only when well away from every window already on BOTH
    // target members.  Applied to every pass (a self-overlapping toolpath yields
    // coincident crossings even leg↔leg).
    // Bounding-sphere overlap factor for the window dedup (see `overlaps_existing`):
    // a new window is dropped if its bounding sphere meets an existing one on the
    // same face, scaled by this factor.  1.0 = spheres just touching; a little
    // above 1 leaves a safety gap so SpaceClaim never sees intersecting holes.
    let dedup_factor = std::env::var("CADCORE_DEDUP_R")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(1.15);
    let _ = r_fil;
    // A/B switch: set CADCORE_NO_XPASS to disable the self- and elbow↔leg
    // crossing passes (measure the leg↔leg-only baseline).
    let extra_xpass = std::env::var("CADCORE_NO_XPASS").is_err();

    // Two legs whose axes are within ~10° of parallel never form a transverse
    // bite loop: they either interpenetrate collinearly or merely touch, both of
    // which are already manifold (each tube is a closed solid).  Tracing them is
    // the dense-fill perf trap — adjacent raster lines (parallel, one filament
    // radius apart) would each spawn dense seeds and divergent traces along the
    // WHOLE overlap.  Reject them cheaply before the SSI trace.
    const PARALLEL_COS: f64 = 0.985;
    let leg_axis = |fi: usize, li: usize| filaments[fi].legs[li].surf.axis().as_vec();

    let timing = std::env::var("CADCORE_TIME").is_ok();
    let (mut n_bp1, mut n_trace1, mut n_bp3, mut n_trace3) = (0u64, 0u64, 0u64, 0u64);
    let t_start = std::time::Instant::now();

    // (1) leg↔leg between DIFFERENT filaments — transverse grid crossings.
    for fi in 0..nf {
        for fj in (fi + 1)..nf {
            for li in 0..filaments[fi].legs.len() {
                for lj in 0..filaments[fj].legs.len() {
                    if !leg_boxes[fi][li].overlaps(&leg_boxes[fj][lj]) {
                        continue; // broad-phase reject: boxes disjoint → no crossing
                    }
                    // Perf reject only for imported dense fill (many legs); a
                    // generated scaffold traces every overlapping pair so no tight
                    // crossing is missed.
                    if strict_rims && leg_axis(fi, li).dot(leg_axis(fj, lj)).abs() > PARALLEL_COS {
                        continue; // near-parallel → no transverse loop
                    }
                    n_bp1 += 1;
                    n_trace1 += 1;
                    for c in closed_loops_between(&leg_tubes[fi][li], &leg_tubes[fj][lj], &opts) {
                        // Scaffolds keep even small transverse loops (tight HL
                        // grids); only imported dense fill needs the larger floor.
                        let floor = if strict_rims { min_pts } else { 6 };
                        if c.points.len() < floor {
                            continue;
                        }
                        let cen = win_centroid(&c.points);
                        let rad = win_radius(&c.points, cen);
                        if strict_rims
                            && (!window_clears_free_ends(&filaments[fi].legs[li], cen, rad, true, true)
                                || !window_clears_free_ends(&filaments[fj].legs[lj], cen, rad, true, true)
                                || overlaps_existing(&windows[fi][li], cen, rad, dedup_factor)
                                || overlaps_existing(&windows[fj][lj], cen, rad, dedup_factor))
                        {
                            continue;
                        }
                        let id = win_id;
                        win_id += 1;
                        let dp = decimate_loop(&c.points, win_max_pts);
                        windows[fi][li].push((id, dp.clone()));
                        windows[fj][lj].push((id, dp));
                    }
                }
            }
        }
    }

    if timing {
        eprintln!(
            "[scaffold][time] pass1 done: {} legpair-traces, {} ms",
            n_trace1,
            t_start.elapsed().as_millis()
        );
    }

    // (2) leg↔leg within the SAME filament, non-adjacent legs — self-crossings
    // (a continuous serpentine that weaves back over itself, e.g. an imported
    // infill grid printed as ONE path).  Adjacent legs (|li-lj| == 1) share an
    // elbow junction: their contact is the tangent junction circle, not a
    // crossing, so it is skipped.
    for fi in 0..nf {
        if !extra_xpass {
            break;
        }
        let nl = filaments[fi].legs.len();
        for li in 0..nl {
            for lj in (li + 2)..nl {
                if !leg_boxes[fi][li].overlaps(&leg_boxes[fi][lj]) {
                    continue;
                }
                if leg_axis(fi, li).dot(leg_axis(fi, lj)).abs() > PARALLEL_COS {
                    continue; // near-parallel → no transverse loop
                }
                for c in closed_loops_between(&leg_tubes[fi][li], &leg_tubes[fi][lj], &opts) {
                    if c.points.len() < min_pts {
                        continue;
                    }
                    let cen = win_centroid(&c.points);
                    let rad = win_radius(&c.points, cen);
                    if strict_rims
                        && (!window_clears_free_ends(&filaments[fi].legs[li], cen, rad, true, true)
                            || !window_clears_free_ends(&filaments[fi].legs[lj], cen, rad, true, true)
                            || overlaps_existing(&windows[fi][li], cen, rad, dedup_factor)
                            || overlaps_existing(&windows[fi][lj], cen, rad, dedup_factor))
                    {
                        continue;
                    }
                    let id = win_id;
                    win_id += 1;
                    let dp = decimate_loop(&c.points, win_max_pts);
                    windows[fi][li].push((id, dp.clone()));
                    windows[fi][lj].push((id, dp));
                }
            }
        }
    }

    // (3) elbow↔leg — a U-turn fillet overlapping a neighbouring leg (any
    // filament).  Skip the elbow's OWN two junction legs (ei and ei+1 of the
    // same filament): those meet the elbow tangentially at the junction
    // circles, already welded seamlessly.
    for fi in 0..nf {
        if !extra_xpass {
            break;
        }
        for ei in 0..filaments[fi].elbows.len() {
            for fj in 0..nf {
                for lj in 0..filaments[fj].legs.len() {
                    if fi == fj && (lj == ei || lj == ei + 1) {
                        continue;
                    }
                    if !elbow_boxes[fi][ei].overlaps(&leg_boxes[fj][lj]) {
                        continue;
                    }
                    // NB: NO parallel-reject here.  A leg parallel to one of the
                    // elbow's arms can still cross the fillet transversally (e.g.
                    // a neighbouring layer's leg passing over the bend), so only
                    // the trace itself can decide — reject by axis would drop
                    // real elbow crossings (leg-over-elbow fusion).
                    n_bp3 += 1;
                    n_trace3 += 1;
                    for c in closed_loops_between(&elbow_tubes[fi][ei], &leg_tubes[fj][lj], &opts) {
                        let cen = win_centroid(&c.points);
                        let rad = win_radius(&c.points, cen);
                        let in_arc = points_within_arc(&filaments[fi].elbows[ei], &c.points);
                        let rim_ok = !strict_rims
                            || (window_clears_free_ends(&filaments[fj].legs[lj], cen, rad, true, true)
                                && !overlaps_existing(&elbow_windows[fi][ei], cen, rad, dedup_factor)
                                && !overlaps_existing(&windows[fj][lj], cen, rad, dedup_factor));
                        if std::env::var("CADCORE_DUMP_FIL").is_ok() {
                            eprintln!("[p3] elbow({fi},{ei}) x leg({fj},{lj}): pts={} rad={:.3} in_arc={in_arc} rim_ok={rim_ok}", c.points.len(), rad);
                        }
                        if c.points.len() < min_pts || !in_arc {
                            continue;
                        }
                        if !rim_ok {
                            continue;
                        }
                        let id = win_id;
                        win_id += 1;
                        let dp = decimate_loop(&c.points, win_max_pts);
                        elbow_windows[fi][ei].push((id, dp.clone()));
                        windows[fj][lj].push((id, dp));
                    }
                }
            }
        }
    }

    if timing {
        let nlegs: usize = filaments.iter().map(|f| f.legs.len()).sum();
        eprintln!(
            "[scaffold][time] {} filaments, {} legs | pass1 {} traces (bp {}), pass3 {} traces (bp {}) | crossings pre-emit {} ms",
            nf, nlegs, n_trace1, n_bp1, n_trace3, n_bp3, t_start.elapsed().as_millis()
        );
    }
    let t_emit = std::time::Instant::now();

    let mut asm = Assembler::new(brep, 1e-7);
    let circ_at = |l: &LegSeg, v: f64| -> Circle3 {
        Circle3::new(
            l.surf.frame.origin + l.surf.axis().as_vec() * v,
            l.surf.axis(),
            l.surf.radius,
        )
    };

    // Global leg index base per filament, so every leg in the whole model gets
    // a unique index.  Each leg L owns two junction weld-groups: `2·L` (its
    // incoming rim) and `2·L+1` (its outgoing rim).  A leg and its mating
    // elbow/cap pass the SAME group for their shared circle, so they weld; two
    // DISTINCT junctions that merely coincide in space (a self-touching
    // serpentine, or two tangent filaments) carry different groups and stay
    // separate twice-used edges — the manifold invariant (README #2).
    let mut leg_base = vec![0u32; nf];
    {
        let mut acc = 0u32;
        for fi in 0..nf {
            leg_base[fi] = acc;
            acc += filaments[fi].legs.len() as u32;
        }
    }
    let g_in = |gbase: u32, li: usize| 2 * (gbase + li as u32);
    let g_out = |gbase: u32, li: usize| 2 * (gbase + li as u32) + 1;

    let dbg_fil = std::env::var("CADCORE_DUMP_FIL").is_ok();
    for fi in 0..nf {
        let gbase = leg_base[fi];
        let before = asm.faces().len();
        let f = &filaments[fi];
        let n = f.legs.len();
        let j_in: Vec<Circle3> = f.legs.iter().map(|l| circ_at(l, 0.0)).collect();
        let j_out: Vec<Circle3> = f.legs.iter().map(|l| circ_at(l, l.length)).collect();

        for li in 0..n {
            asm.emit_cylinder_face(
                f.legs[li].surf,
                f.legs[li].length,
                j_in[li],
                j_out[li],
                &windows[fi][li],
                g_in(gbase, li),
                g_out(gbase, li),
            );
            if li == 0 {
                let nrm = UnitVec3::try_from_vec(f.legs[li].surf.axis().as_vec() * -1.0).unwrap();
                asm.emit_disk_cap(
                    Plane3::from_origin_normal(f.legs[li].surf.frame.origin, nrm),
                    j_in[li],
                    g_in(gbase, li),
                );
            }
            if li + 1 == n {
                let c = f.legs[li].surf.frame.origin
                    + f.legs[li].surf.axis().as_vec() * f.legs[li].length;
                asm.emit_disk_cap(
                    Plane3::from_origin_normal(c, f.legs[li].surf.axis()),
                    j_out[li],
                    g_out(gbase, li),
                );
            }
        }
        for (ei, e) in f.elbows.iter().enumerate() {
            // elbow ei joins leg ei (its outgoing rim) and leg ei+1 (its incoming rim).
            asm.emit_torus_fillet_windows(
                e.surf,
                j_out[ei],
                j_in[ei + 1],
                &elbow_windows[fi][ei],
                g_out(gbase, ei),
                g_in(gbase, ei + 1),
            );
        }
        if dbg_fil {
            let ids: Vec<String> = asm.faces()[before..].iter().map(|f| format!("{:?}", f)).collect();
            eprintln!("[fil {fi}] legs={} elbows={} faces=[{}]", f.legs.len(), f.elbows.len(), ids.join(","));
        }
    }
    let _ = &fil_tubes; // (retained for future elbow↔elbow crossing support)

    let faces = asm.faces().to_vec();
    let shell = brep.add_shell(Shell {
        faces: faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = brep.add_solid(Solid {
        shells: vec![shell],
        name: Some("scaffold".into()),
    });
    brep.shells[shell].solid = solid;
    for f in &faces {
        brep.faces[*f].shell = shell;
    }
    if timing {
        eprintln!(
            "[scaffold][time] emit {} faces: {} ms | total {} ms",
            faces.len(),
            t_emit.elapsed().as_millis(),
            t_start.elapsed().as_millis()
        );
    }
    shell
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::UnionConfig;
    use crate::validate::validate_shell;
    use crate::Diagnostics;
    use std::collections::HashMap;

    fn open_edges(brep: &BRep, shell: ShellId) -> usize {
        let mut uses: HashMap<cadcore_topo::EdgeId, usize> = HashMap::new();
        for &fid in &brep.shells[shell].faces {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let st = brep.loops[lid].start;
                let mut c = st;
                loop {
                    let ce = &brep.coedges[c];
                    *uses.entry(ce.edge).or_insert(0) += 1;
                    c = ce.next;
                    if c == st { break; }
                }
            }
        }
        uses.values().filter(|&&n| n == 1).count()
    }

    /// Number of connected components of the shell, faces joined by shared
    /// edges.  A genuine fusion of crossing filaments is ONE component; two
    /// merely-coexisting (interpenetrating) solids are two.
    fn components(brep: &BRep, shell: ShellId) -> usize {
        let faces = &brep.shells[shell].faces;
        // edge -> faces touching it
        let mut by_edge: HashMap<cadcore_topo::EdgeId, Vec<usize>> = HashMap::new();
        for (i, &fid) in faces.iter().enumerate() {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let st = brep.loops[lid].start;
                let mut c = st;
                loop {
                    by_edge.entry(brep.coedges[c].edge).or_default().push(i);
                    c = brep.coedges[c].next;
                    if c == st { break; }
                }
            }
        }
        let mut parent: Vec<usize> = (0..faces.len()).collect();
        fn find(p: &mut [usize], x: usize) -> usize {
            let mut r = x;
            while p[r] != r { r = p[r]; }
            let mut c = x;
            while p[c] != c { let n = p[c]; p[c] = r; c = n; }
            r
        }
        for fs in by_edge.values() {
            for w in fs.windows(2) {
                let (a, b) = (find(&mut parent, w[0]), find(&mut parent, w[1]));
                parent[a] = b;
            }
        }
        (0..faces.len()).filter(|&i| find(&mut parent, i) == i).count()
    }

    /// Two serpentine filaments on stacked layers that cross transversally
    /// (X-running vs Y-running) — fused into ONE watertight, manifold solid.
    #[test]
    fn two_crossing_serpentines_watertight() {
        let r = 0.275;
        // filament A: U in z=0, legs run along X at y=0 and y=2 (connector at x=0)
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
        // filament B: U one layer up (z=0.35), legs run along Y at x=-1.5 and
        // x=-0.8 so each crosses BOTH of A's X-legs transversally; its elbows
        // sit at y=3, clear of A.
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
        let shell = fuse_scaffold(&mut brep, &[a, b]);
        assert_eq!(open_edges(&brep, shell), 0, "scaffold watertight");
        assert_eq!(components(&brep, shell), 1, "filaments fused into one solid");

        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let report = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(report.manifold_violations.is_empty(), "manifold");
        assert!(report.distance_violations.is_empty(),
            "distance max {:.2}um", report.max_edge_deviation * 1e3);
    }

    /// A leg of one filament passing directly OVER the U-turn elbow of another
    /// (stacked layers) must fuse watertight — the elbow↔leg crossing the old
    /// leg↔leg-only pass could not see (the dense-fill / hairpin-over-neighbour
    /// case).  Regression for `emit_torus_fillet_windows` + the elbow pass.
    #[test]
    fn leg_over_elbow_watertight() {
        let r = 0.275;
        // A: L-corner in z=0 (+X then 90° elbow then +Y).  Its elbow torus is
        // centred at (-0.5, 0.5, 0); the fillet tube arcs through ~(-0.15,0.15).
        let a = Filament::serpentine(
            &[
                Point3::new(-2.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
            ],
            r,
            0.5,
        );
        // B: a straight leg one layer up (z=0.35) running along +X across the
        // elbow apex — its tube (z∈[0.075,0.625]) overlaps the elbow tube
        // (z∈[-0.275,0.275]) and crosses it transversally.
        let b = Filament::serpentine(
            &[
                Point3::new(-2.0, 0.15, 0.35),
                Point3::new(2.0, 0.15, 0.35),
            ],
            r,
            0.5,
        );
        let mut brep = BRep::new();
        let shell = fuse_scaffold(&mut brep, &[a, b]);
        assert_eq!(open_edges(&brep, shell), 0, "leg-over-elbow watertight");
        assert_eq!(components(&brep, shell), 1, "fused into one solid");

        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let report = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(report.manifold_violations.is_empty(), "manifold");
        assert!(
            report.distance_violations.is_empty(),
            "distance max {:.2}um",
            report.max_edge_deviation * 1e3
        );
    }

    /// A single continuous serpentine that weaves back over ITSELF (non-adjacent
    /// legs cross) must fuse watertight — the self-crossing pass.  Two stacked
    /// weaves of the same filament reproduce the imported infill-grid-as-one-path
    /// case.
    ///
    /// NOTE: this fixture was RECONSTRUCTED after the working copy of this file
    /// was lost; the geometry follows the description in the comment below, but
    /// it is not byte-identical to the original.
    #[test]
    fn self_crossing_weave_watertight() {
        let r = 0.275;
        // One filament: X-legs at y=0 and y=1 (z=0), then it climbs and lays
        // Y-legs at x=0 and x=1 (z=0.35) — the later Y-legs cross the earlier
        // X-legs of the SAME filament (non-adjacent legs).
        let f = Filament::serpentine(
            &[
                Point3::new(-1.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(2.0, 1.0, 0.0),
                Point3::new(-1.0, 1.0, 0.0),
                Point3::new(-1.0, 1.0, 0.35),
                Point3::new(0.0, 1.0, 0.35),
                Point3::new(0.0, -1.0, 0.35),
                Point3::new(1.0, -1.0, 0.35),
                Point3::new(1.0, 2.0, 0.35),
            ],
            r,
            0.3,
        );
        let mut brep = BRep::new();
        let shell = fuse_scaffold(&mut brep, &[f]);
        assert_eq!(open_edges(&brep, shell), 0, "self-crossing weave watertight");
        assert_eq!(components(&brep, shell), 1, "one solid");

        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let report = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(report.manifold_violations.is_empty(), "manifold");
    }

    /// Four stacked layers (the smallest model with crossings on three
    /// interfaces) must fuse into ONE watertight solid, and the broad-phase cull
    /// must keep it fast enough to be usable on a real part.
    ///
    /// NOTE: this fixture was RECONSTRUCTED after the working copy of this file
    /// was lost; it exercises the same path (multi-layer alternating weave) but
    /// is not byte-identical to the original.
    #[test]
    fn four_layer_weave_fuses_fast_and_watertight() {
        let r = 0.275;
        let hl = 0.35;
        let mut fils = Vec::new();
        for layer in 0..4 {
            let z = layer as f64 * hl;
            // alternate the running direction per layer
            let verts: Vec<Point3> = if layer % 2 == 0 {
                vec![
                    Point3::new(0.0, 0.0, z),
                    Point3::new(3.0, 0.0, z),
                    Point3::new(3.0, 1.2, z),
                    Point3::new(0.0, 1.2, z),
                    Point3::new(0.0, 2.4, z),
                    Point3::new(3.0, 2.4, z),
                ]
            } else {
                vec![
                    Point3::new(0.0, 0.0, z),
                    Point3::new(0.0, 3.0, z),
                    Point3::new(1.2, 3.0, z),
                    Point3::new(1.2, 0.0, z),
                    Point3::new(2.4, 0.0, z),
                    Point3::new(2.4, 3.0, z),
                ]
            };
            fils.push(Filament::serpentine(&verts, r, 0.5));
        }
        let mut brep = BRep::new();
        let t = std::time::Instant::now();
        let shell = fuse_scaffold(&mut brep, &fils);
        let ms = t.elapsed().as_millis();

        assert_eq!(open_edges(&brep, shell), 0, "four-layer weave watertight");
        assert_eq!(components(&brep, shell), 1, "fused into one solid");
        assert!(ms < 60_000, "fuse took {ms} ms — the broad phase regressed");
    }
}
