//! The general union orchestrator: `union(A, B)` over two arbitrary B-Reps.
//!
//! Composes the solid-agnostic primitives:
//!   broad phase ([`super::aabb`]) → SSI ([`super::ssi`]) → imprint+classify
//!   ([`crate::arrange::cells`]) → stitch ([`crate::arrange::assembly`]).
//!
//! For each face we build its parameter domain, feed the DCEL its own boundary
//! plus every SSI curve lying on it, and keep the cells whose material is
//! OUTSIDE the other solid (the union rule).  Because each SSI curve is handed
//! — as identical 3-D geometry — to both faces it came from, the assembler's
//! geometric edge cache welds the two kept pieces into one shared edge, so the
//! result is watertight with no special-casing.

use std::collections::HashMap;

use cadcore_math::Point3;
use cadcore_topo::{BRep, FaceBoundary, FaceExtent, FaceGeom, FaceId, Shell, ShellId, Solid};

use super::aabb::{candidate_pairs, face_aabb};
use super::contain::point_in_solid;
use super::ssi::{intersect_faces, SsiCurve};
use crate::arrange::assembly::Assembler;
use crate::arrange::cells::{arrange_and_classify_with, FaceChain};
use crate::arrange::domain::FaceDomain;
use crate::geom::refine::AnalyticSurface;

/// Boolean **union** of solid `a` (faces `a_faces`) and solid `b` (faces
/// `b_faces`) into a fresh watertight B-Rep.  Returns the output B-Rep and the
/// id of its single closed shell.
///
/// Only the planar-face case (e.g. boxes / prisms) is wired today; faces whose
/// carrier surface or SSI is not yet supported are emitted untrimmed when they
/// don't meet the other solid, and skipped when they do (incomplete, never
/// wrong-but-silent — the caller validates the result).
pub fn union(
    a: &BRep,
    a_faces: &[FaceId],
    b: &BRep,
    b_faces: &[FaceId],
    knit: f64,
) -> (BRep, ShellId) {
    union_n(&[(a, a_faces), (b, b_faces)], knit).expect("two solids")
}

/// Boolean union of **N** solids in ONE direct pass — the scaffold path.
///
/// Every cross-solid face pair is intersected ONCE, its shared curve pre-split
/// at all incident cylinder seams, and stored on BOTH faces.  Then each
/// *original* face is imprinted and classified against ALL other solids
/// (material kept ⇔ outside every other solid) and emitted through one
/// assembler.  Because no intermediate result is ever re-processed, existing
/// welds are never re-split — unlike a pairwise fold, this stays watertight on
/// periodic (cylinder) faces.
///
/// `None` for fewer than two solids.
pub fn union_n(solids: &[(&BRep, &[FaceId])], knit: f64) -> Option<(BRep, ShellId)> {
    let n = solids.len();
    if n < 2 {
        return None;
    }

    // Each solid's curved lateral surface(s) — for splitting intersection curves
    // at TRIPLE points where a third body's surface crosses them.
    let solid_surfs: Vec<Vec<AnalyticSurface>> = solids
        .iter()
        .map(|(b, faces)| {
            faces
                .iter()
                .filter_map(|&f| match &b.faces[f].geom {
                    FaceGeom::Cylinder(s) => Some(AnalyticSurface::Cylinder(*s)),
                    FaceGeom::Torus(s) => Some(AnalyticSurface::Torus(*s)),
                    _ => None,
                })
                .collect()
        })
        .collect();

    // SSI per (solid index, face), each cross-solid pair traced ONCE and
    // pre-split identically for both faces it lies on.
    let mut ssi: Vec<HashMap<FaceId, Vec<SsiCurve>>> = vec![HashMap::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            let (bi, fi) = solids[i];
            let (bj, fj) = solids[j];
            // every OTHER body's lateral surface: where this (i,j) curve crosses
            // one of them, all three meet → a triple point that must be a shared
            // vertex on every curve through it.
            let crossers: Vec<AnalyticSurface> = (0..n)
                .filter(|&k| k != i && k != j)
                .flat_map(|k| solid_surfs[k].iter().copied())
                .collect();
            for (fa, fb) in candidate_pairs(bi, fi, bj, fj, knit) {
                let curves = intersect_faces(bi, fa, bj, fb);
                if curves.is_empty() {
                    continue;
                }
                // periodic carriers the curve lies on (their seams must be
                // pre-split so both faces subdivide it identically).
                let mut seams: Vec<AnalyticSurface> = Vec::new();
                for g in [&bi.faces[fa].geom, &bj.faces[fb].geom] {
                    match g {
                        FaceGeom::Cylinder(s) => seams.push(AnalyticSurface::Cylinder(*s)),
                        FaceGeom::Torus(s) => seams.push(AnalyticSurface::Torus(*s)),
                        _ => {}
                    }
                }
                for c in curves {
                    for arc in split_curve(&c.points, c.closed, &seams, &crossers) {
                        let piece = SsiCurve { points: arc, closed: false };
                        ssi[i].entry(fa).or_default().push(piece.clone());
                        ssi[j].entry(fb).or_default().push(piece);
                    }
                }
            }
        }
    }

    // Triple-point registry: weld the ENDPOINTS of all SSI curves so that
    // several curves meeting at a triple point (two surfaces ∩ a cap plane —
    // e.g. a ruling line and an end-cap arc on parallel tubes) share the exact
    // same vertex.  Otherwise their endpoints differ by the sampling jitter,
    // the DCEL leaves dangling edges, and the face never splits.
    snap_ssi_endpoints(&mut ssi, 0.03);

    // imprint + classify (outside EVERY other solid) + stitch, one assembler.
    let mut out = BRep::new();
    // Weld tolerance: merges the residual seam-interpolation noise that
    // `build_periodic_u` leaves at a curve's high-curvature turns (up to a few
    // 1e-4 where two periodic carriers meet, e.g. torus×cylinder), while
    // staying far below the curves' real sample spacing (~1e-2).
    let weld = knit.clamp(5e-4, 1e-2);
    let mut asm = Assembler::new(&mut out, weld);
    for i in 0..n {
        let (bi, fi) = solids[i];
        let inside_other = |p: Point3| {
            (0..n).any(|j| j != i && point_in_solid(solids[j].0, solids[j].1, p))
        };
        emit_solid_faces(bi, fi, &ssi[i], &inside_other, &mut asm);
    }

    // finalise the shell.
    let faces = asm.faces().to_vec();
    let shell = out.add_shell(Shell {
        faces: faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = out.add_solid(Solid {
        shells: vec![shell],
        name: Some("union".into()),
    });
    out.shells[shell].solid = solid;
    for f in &faces {
        out.faces[*f].shell = shell;
    }
    Some((out, shell))
}

/// Boolean union of **N** solids.
///
/// Planar solids may also be folded pairwise (each `Trimmed` result feeds the
/// next), but [`union_n`]'s single direct pass is the general path and the only
/// one that stays watertight when periodic (cylinder) faces meet, so this
/// simply delegates to it.
pub fn union_many(solids: &[(&BRep, &[FaceId])], knit: f64) -> Option<(BRep, ShellId)> {
    union_n(solids, knit)
}

/// Imprint + classify + emit every face of one solid.
fn emit_solid_faces(
    brep: &BRep,
    faces: &[FaceId],
    ssi: &HashMap<FaceId, Vec<SsiCurve>>,
    inside_other: &dyn Fn(Point3) -> bool,
    asm: &mut Assembler,
) {
    for &fid in faces {
        let f = &brep.faces[fid];
        let Some(domain) = face_domain(brep, fid) else {
            continue; // unsupported carrier surface — skip (incomplete)
        };

        // boundary chains: the face's own loops, projected into uv & closed.
        let mut chains: Vec<FaceChain> = Vec::new();
        let mut tag = 1u32;
        for chain in face_boundary_chains(brep, fid, &domain) {
            chains.push(FaceChain { pts: chain, tag });
            tag += 1;
        }
        // SSI chains lying on this face.
        let mut stag = 1000u32;
        if let Some(curves) = ssi.get(&fid) {
            for c in curves {
                let mut pts: Vec<(f64, f64)> = c.points.iter().map(|&p| domain.uv(p)).collect();
                if c.closed && pts.len() >= 3 {
                    pts.push(pts[0]); // close the loop for the DCEL
                }
                if pts.len() >= 2 {
                    chains.push(FaceChain { pts, tag: stag });
                    stag += 1;
                }
            }
        }

        // A full-donut TorusBand is arranged non-periodically in θ, so every
        // chain's θ must be continuous (no atan2 ±π jump): unwrap each chain's
        // u onto the branch nearest its previous point.
        if matches!(domain, FaceDomain::TorusBand { .. }) {
            let tau = std::f64::consts::TAU;
            for ch in &mut chains {
                for i in 1..ch.pts.len() {
                    let du = ch.pts[i].0 - ch.pts[i - 1].0;
                    ch.pts[i].0 -= tau * (du / tau).round();
                    let dv = ch.pts[i].1 - ch.pts[i - 1].1;
                    ch.pts[i].1 -= tau * (dv / tau).round();
                }
                // An SSI loop that wraps the full ring (every θ, with φ = v(θ) —
                // coaxial stack: constant v; offset stack: a varying band) is a
                // chord spanning the [0,2π]² rectangle left→right.  Whichever way
                // it was traced, rebuild it as v(u) sampled cleanly on u∈[0,2π]
                // so it touches BOTH θ-seams (welds) and carries the true φ.
                if ch.tag >= 1000 {
                    let umin = ch.pts.iter().cloned().fold(f64::MAX, |m, p| m.min(p.0));
                    let umax = ch.pts.iter().cloned().fold(f64::MIN, |m, p| m.max(p.0));
                    if umax - umin >= tau - 0.25 {
                        ch.pts = resample_fullwrap_chord(&ch.pts, 96);
                    }
                }
            }
        }

        if std::env::var("CADCORE_DUMP_CHAINS").is_ok() && matches!(domain, FaceDomain::TorusBand { .. }) {
            for ch in &chains {
                let us: Vec<f64> = ch.pts.iter().map(|p| p.0).collect();
                let vs: Vec<f64> = ch.pts.iter().map(|p| p.1).collect();
                let closed = ch.pts.len() >= 2 && (ch.pts[0].0 - ch.pts[ch.pts.len()-1].0).abs() < 1e-6 && (ch.pts[0].1 - ch.pts[ch.pts.len()-1].1).abs() < 1e-6;
                eprintln!("  chain tag={} n={} closed={} u[{:.2}..{:.2}] v[{:.2}..{:.2}]", ch.tag, ch.pts.len(), closed,
                    us.iter().cloned().fold(f64::MAX, f64::min), us.iter().cloned().fold(f64::MIN, f64::max),
                    vs.iter().cloned().fold(f64::MAX, f64::min), vs.iter().cloned().fold(f64::MIN, f64::max));
            }
        }

        let cells = arrange_and_classify_with(&domain, &chains, 1e-3, inside_other);
        if std::env::var("CADCORE_DUMP_UNION").is_ok() {
            let k = match &f.geom { FaceGeom::Plane(_) => "plane", FaceGeom::Cylinder(_) => "cyl", FaceGeom::Torus(_) => "TORUS", FaceGeom::Sphere(_) => "sph" };
            eprintln!("face {:?} {}: bchains={} ssi={} cells={} kept={}", fid, k,
                chains.len() - ssi.get(&fid).map(|v| v.len()).unwrap_or(0),
                ssi.get(&fid).map(|v| v.len()).unwrap_or(0), cells.len(), cells.iter().filter(|c| c.keep).count());
        }
        for cell in cells {
            if cell.keep {
                asm.emit_cell(&domain, f.geom.clone(), f.normal, &cell);
            }
        }
    }
}

/// Resample a θ-wrapping SSI loop (uv points, u already unwrapped) into a clean
/// chord on `u ∈ [0, 2π]` carrying the true `v = φ(u)`.  The loop wraps θ once
/// (possibly traced backwards or in a shifted ±2π window); we fold u back into a
/// single period, sort, and linearly interpolate `v` at a uniform u-grid with a
/// periodic wrap so both endpoints land exactly on the θ-seams (u=0 and u=2π)
/// with the same v — making the chord weld to both seam edges of the rectangle.
fn resample_fullwrap_chord(pts: &[(f64, f64)], n: usize) -> Vec<(f64, f64)> {
    let tau = std::f64::consts::TAU;
    let mut s: Vec<(f64, f64)> = pts.iter().map(|&(u, v)| (u.rem_euclid(tau), v)).collect();
    s.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    s.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-9);
    if s.len() < 2 {
        let v = pts.iter().map(|p| p.1).sum::<f64>() / pts.len().max(1) as f64;
        return (0..=n).map(|k| (k as f64 / n as f64 * tau, v)).collect();
    }
    // Extend by one wrapped sample on each side so u=0 and u=2π are bracketed.
    let first = *s.first().unwrap();
    let last = *s.last().unwrap();
    let mut ext = Vec::with_capacity(s.len() + 2);
    ext.push((last.0 - tau, last.1));
    ext.extend_from_slice(&s);
    ext.push((first.0 + tau, first.1));
    (0..=n)
        .map(|k| {
            let u = k as f64 / n as f64 * tau;
            let mut v = ext[0].1;
            for w in ext.windows(2) {
                if u >= w[0].0 - 1e-12 && u <= w[1].0 + 1e-12 {
                    let d = w[1].0 - w[0].0;
                    let t = if d.abs() < 1e-12 { 0.0 } else { (u - w[0].0) / d };
                    v = w[0].1 + t * (w[1].1 - w[0].1);
                    break;
                }
            }
            (u, v)
        })
        .collect()
}

/// Weld the endpoints of every SSI curve (across ALL solids' faces) onto a
/// shared set of representative points within `tol`.  Curves that meet at a
/// triple point then carry the IDENTICAL endpoint, so the per-face DCEL
/// connects them instead of leaving dangling edges.
fn snap_ssi_endpoints(ssi: &mut [HashMap<FaceId, Vec<SsiCurve>>], tol: f64) {
    // 1) gather every endpoint, tagged exact (from a 2-point curve — a straight
    //    ruling line whose ends ARE the true triple points) or approximate
    //    (from a sampled arc).
    let mut pts: Vec<(Point3, bool)> = Vec::new();
    for map in ssi.iter() {
        for curves in map.values() {
            for c in curves.iter() {
                if c.points.len() >= 2 {
                    let exact = c.points.len() == 2;
                    pts.push((c.points[0], exact));
                    pts.push((*c.points.last().unwrap(), exact));
                }
            }
        }
    }
    if pts.is_empty() {
        return;
    }
    // 2) build representatives DETERMINISTICALLY, preferring EXACT endpoints so
    //    sampled arc ends snap onto the true triple points (not vice versa).
    pts.sort_by(|a, b| {
        (!a.1, a.0.x, a.0.y, a.0.z)
            .partial_cmp(&(!b.1, b.0.x, b.0.y, b.0.z))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut reps: Vec<Point3> = Vec::new();
    for &(p, _) in &pts {
        if !reps.iter().any(|&r| (r - p).length() < tol) {
            reps.push(p);
        }
    }
    // 3) snap each curve endpoint to its nearest representative within tol.
    let nearest = |p: Point3| -> Point3 {
        reps.iter()
            .copied()
            .filter(|&r| (r - p).length() < tol)
            .min_by(|&u, &v| {
                (u - p)
                    .length()
                    .partial_cmp(&(v - p).length())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(p)
    };
    for map in ssi.iter_mut() {
        for curves in map.values_mut() {
            for c in curves.iter_mut() {
                let m = c.points.len();
                if m >= 2 {
                    c.points[0] = nearest(c.points[0]);
                    c.points[m - 1] = nearest(c.points[m - 1]);
                }
            }
        }
    }
}

/// Split a closed loop into open arcs at every periodic carrier's seam
/// crossing — a cylinder's `u=0` half-plane (`fy=0 & fx>0`) and a torus's
/// `φ=0` outer equator (`h=0 & a>R`).  Split points are refined onto ALL the
/// carriers so they are exact shared vertices; both faces consume the SAME
/// arcs ⇒ identical subdivision.
/// Signed distance of `p` to a crosser surface (0 ⇔ on it; sign = in/out).
fn signed_dist(s: &AnalyticSurface, p: Point3) -> f64 {
    match s {
        AnalyticSurface::Cylinder(c) => {
            let w = p - c.frame.origin;
            let axial = c.axis().dot_vec(w);
            let radial = (w - c.axis().as_vec() * axial).length();
            radial - c.radius
        }
        AnalyticSurface::Torus(t) => {
            let w = p - t.frame.origin;
            let h = t.frame.z.dot_vec(w);
            let a = (w - t.frame.z.as_vec() * h).length();
            ((a - t.major_radius).powi(2) + h * h).sqrt() - t.minor_radius
        }
        AnalyticSurface::Plane(_) => f64::MAX,
    }
}

/// Append cuts where the polyline TRANSVERSALLY crosses a third body's surface
/// (the signed distance changes sign) — a triple point.  Each cut is refined
/// onto the curve's carriers PLUS that surface, so it is an exact shared vertex
/// every curve through the triple point converges to.  `nseg` segments are
/// scanned (`n` for a closed loop, `n-1` for an open arc).
fn add_crosser_cuts(
    pts: &[Point3],
    n: usize,
    nseg: usize,
    carriers: &[AnalyticSurface],
    crossers: &[AnalyticSurface],
    cuts: &mut Vec<(usize, f64, Point3)>,
) {
    use crate::geom::refine::refine_point;
    for s in crossers {
        let onto: Vec<AnalyticSurface> =
            carriers.iter().chain(std::iter::once(s)).copied().collect();
        for i in 0..nseg {
            let p0 = pts[i];
            let p1 = pts[(i + 1) % n];
            let v0 = signed_dist(s, p0);
            let v1 = signed_dist(s, p1);
            if (v0 > 0.0) != (v1 > 0.0) && (v0 - v1).abs() > 1e-15 {
                let t = v0 / (v0 - v1);
                let approx = p0 + (p1 - p0) * t;
                let (exact, _) = refine_point(approx, &onto, 1e-12, 40);
                cuts.push((i, t, exact));
            }
        }
    }
}

/// Split an SSI curve at all carrier seams (closed loops) and all third-body
/// crossings (triple points), returning shared arcs.
fn split_curve(
    points: &[Point3],
    closed: bool,
    carriers: &[AnalyticSurface],
    crossers: &[AnalyticSurface],
) -> Vec<Vec<Point3>> {
    if closed && !carriers.is_empty() {
        presplit_loop(points, carriers, crossers)
    } else {
        presplit_open(points, carriers, crossers)
    }
}

/// Split an OPEN arc at third-body crossings (triple points), linear walk.
fn presplit_open(
    points: &[Point3],
    carriers: &[AnalyticSurface],
    crossers: &[AnalyticSurface],
) -> Vec<Vec<Point3>> {
    let n = points.len();
    if n < 2 {
        return vec![points.to_vec()];
    }
    let mut cuts: Vec<(usize, f64, Point3)> = Vec::new();
    add_crosser_cuts(points, n, n - 1, carriers, crossers, &mut cuts);
    if cuts.is_empty() {
        return vec![points.to_vec()];
    }
    cuts.sort_by(|a, b| (a.0, a.1).partial_cmp(&(b.0, b.1)).unwrap());
    cuts.dedup_by(|b, a| (a.2 - b.2).length() < 1e-3);
    let mut arcs: Vec<Vec<Point3>> = Vec::new();
    let mut cur: Vec<Point3> = vec![points[0]];
    let mut ci = 0usize;
    for i in 0..(n - 1) {
        while ci < cuts.len() && cuts[ci].0 == i {
            cur.push(cuts[ci].2);
            arcs.push(std::mem::take(&mut cur));
            cur = vec![cuts[ci].2];
            ci += 1;
        }
        let endp = points[i + 1];
        if cur.last().map_or(true, |&q| (q - endp).length() > 1e-12) {
            cur.push(endp);
        }
    }
    if cur.len() >= 2 {
        arcs.push(cur);
    }
    arcs
}

fn presplit_loop(
    loop_pts: &[Point3],
    surfs: &[AnalyticSurface],
    crossers: &[AnalyticSurface],
) -> Vec<Vec<Point3>> {
    use crate::geom::refine::refine_point;
    let n = loop_pts.len();
    if n < 3 {
        return vec![loop_pts.to_vec()];
    }
    // (segment index, fraction along it) of each seam crossing → refined point.
    let mut cuts: Vec<(usize, f64, Point3)> = Vec::new();
    add_crosser_cuts(loop_pts, n, n, surfs, crossers, &mut cuts);
    for s in surfs {
        // f(p) = signed coord whose zero (with the side guard) marks the seam.
        let seam_value = |p: Point3| -> (f64, bool) {
            match s {
                AnalyticSurface::Cylinder(c) => {
                    let w = p - c.frame.origin;
                    (c.frame.y.dot_vec(w), c.frame.x.dot_vec(w) > 0.0) // u=0: fy=0, fx>0
                }
                AnalyticSurface::Torus(t) => {
                    let w = p - t.frame.origin;
                    let h = t.frame.z.dot_vec(w);
                    let a = (w - t.frame.z.as_vec() * h).length();
                    (h, a > t.major_radius) // φ=0: h=0, a>R (outer equator)
                }
                AnalyticSurface::Plane(_) => (1.0, false), // no seam
            }
        };
        for i in 0..n {
            let p0 = loop_pts[i];
            let p1 = loop_pts[(i + 1) % n];
            let (v0, _) = seam_value(p0);
            let (v1, _) = seam_value(p1);
            if (v0 > 0.0) != (v1 > 0.0) && (v0 - v1).abs() > 1e-15 {
                let t = v0 / (v0 - v1);
                let approx = p0 + (p1 - p0) * t;
                if seam_value(approx).1 {
                    let (exact, _) = refine_point(approx, surfs, 1e-12, 40);
                    cuts.push((i, t, exact));
                }
            }
        }
    }
    if cuts.is_empty() {
        return vec![loop_pts.to_vec()];
    }
    cuts.sort_by(|a, b| (a.0, a.1).partial_cmp(&(b.0, b.1)).unwrap());
    // Merge cuts that nearly coincide in 3-D: when two carriers' seams cross
    // the loop at almost the same point (e.g. a torus φ-seam and a cylinder
    // u-seam meeting at a loop extreme), keeping both would leave a
    // sub-resolution arc that the two faces' arrangements split differently.
    cuts.dedup_by(|b, a| (a.2 - b.2).length() < 1e-3);
    if (cuts[0].2 - cuts.last().unwrap().2).length() < 1e-3 && cuts.len() > 1 {
        cuts.pop();
    }
    // walk the loop, starting at the first cut, accumulating arcs between cuts
    let mut arcs: Vec<Vec<Point3>> = Vec::new();
    let mut cur: Vec<Point3> = vec![cuts[0].2];
    let mut ci = 1usize; // next cut to consume
    let first_seg = cuts[0].0;
    // iterate segments cyclically starting just after the first cut
    for step in 0..n {
        let seg = (first_seg + step) % n;
        // any cuts on this segment beyond the one we started at
        while ci < cuts.len() && cuts[ci].0 == seg {
            // push the segment's far endpoint up to this cut, close the arc
            cur.push(cuts[ci].2);
            arcs.push(std::mem::take(&mut cur));
            cur = vec![cuts[ci].2];
            ci += 1;
        }
        // add the segment's end vertex (start of next segment)
        let end_v = loop_pts[(seg + 1) % n];
        // skip if essentially the last pushed point
        if cur.last().map_or(true, |&q| (q - end_v).length() > 1e-12) {
            cur.push(end_v);
        }
    }
    // close the final arc back to the first cut
    cur.push(cuts[0].2);
    if cur.len() >= 2 {
        arcs.push(cur);
    }
    arcs
}

/// Build a [`FaceDomain`] for a face (planar + cylinder-band supported).
fn face_domain(brep: &BRep, fid: FaceId) -> Option<FaceDomain> {
    let f = &brep.faces[fid];
    match (&f.geom, &f.extent) {
        (FaceGeom::Plane(pl), _) => {
            let bb = face_aabb(brep, fid);
            let half = (bb.extent() * 0.75).max(1.0);
            Some(FaceDomain::Plane {
                plane: *pl,
                half_extent: half,
            })
        }
        (FaceGeom::Cylinder(surf), FaceExtent::Cylinder { length, .. }) => {
            Some(FaceDomain::CylinderBand { surf: *surf, length: *length })
        }
        (FaceGeom::Cylinder(surf), _) => {
            // Trimmed cylinder (prior-union output): recover the band length.
            let axis = surf.axis();
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            for p in face_loop_points(brep, fid) {
                let t = axis.dot_vec(p - surf.frame.origin);
                lo = lo.min(t);
                hi = hi.max(t);
            }
            if hi <= lo {
                return None;
            }
            let shifted = cadcore_geom::CylSurf::new(surf.frame.origin + axis.as_vec() * lo, axis, surf.radius);
            Some(FaceDomain::CylinderBand { surf: shifted, length: hi - lo })
        }
        (FaceGeom::Torus(surf), ext) => {
            // FULL torus (a donut): start_circle == end_circle ⇒ θ wraps the
            // whole 2π — use the θ-periodic, φ-open band.
            if let FaceExtent::TorusFillet { start_circle, end_circle } = ext {
                if (start_circle.frame.origin - end_circle.frame.origin).length() < 1e-6 {
                    return Some(FaceDomain::TorusBand { surf: *surf });
                }
            }
            let theta_of = |p: Point3| {
                let w = p - surf.frame.origin;
                surf.frame.y.dot_vec(w).atan2(surf.frame.x.dot_vec(w))
            };
            let (lo, hi) = if let FaceExtent::TorusFillet { start_circle, end_circle } = ext {
                (theta_of(start_circle.frame.origin), theta_of(end_circle.frame.origin))
            } else {
                // Trimmed torus: recover θ-range from the loop vertices.
                let mut lo = f64::MAX;
                let mut hi = f64::MIN;
                let base = face_loop_points(brep, fid)
                    .first()
                    .map(|&p| theta_of(p))
                    .unwrap_or(0.0);
                for p in face_loop_points(brep, fid) {
                    // unwrap onto the branch nearest the first sample
                    let raw = theta_of(p);
                    let t = raw + std::f64::consts::TAU * ((base - raw) / std::f64::consts::TAU).round();
                    lo = lo.min(t);
                    hi = hi.max(t);
                }
                if hi <= lo {
                    return None;
                }
                (lo, hi)
            };
            // put theta_hi on the branch within π of theta_lo (elbow arc < π)
            let hi = hi + std::f64::consts::TAU * ((lo - hi) / std::f64::consts::TAU).round();
            Some(FaceDomain::TorusPatch { surf: *surf, theta_lo: lo, theta_hi: hi })
        }
        _ => None,
    }
}

/// Every loop vertex of a face (outer + inner), sampled from its real edges.
pub(crate) fn face_loop_points(brep: &BRep, fid: FaceId) -> Vec<Point3> {
    let f = &brep.faces[fid];
    let mut out = Vec::new();
    let mut ls = vec![f.outer_loop];
    ls.extend(f.inner_loops.iter().copied());
    for lid in ls {
        let Some(lp) = brep.loops.get(lid) else { continue };
        let start = lp.start;
        if brep.coedges.get(start).is_none() {
            continue;
        }
        let mut c = start;
        loop {
            let ce = &brep.coedges[c];
            out.extend(super::aabb::sample_edge(brep, ce.edge));
            c = ce.next;
            if c == start {
                break;
            }
        }
    }
    out
}

/// The face's boundary loops as closed uv polylines (outer + inners).
fn face_boundary_chains(brep: &BRep, fid: FaceId, domain: &FaceDomain) -> Vec<Vec<(f64, f64)>> {
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
        let mut pts: Vec<(f64, f64)> = Vec::new();
        let mut c = start;
        loop {
            let ce = &brep.coedges[c];
            // sample this edge's 3-D points in traversal order, project to uv
            let mut sp = super::aabb::sample_edge(brep, ce.edge);
            if ce.sense == cadcore_topo::CoEdgeSense::Opposite {
                sp.reverse();
            }
            for p in sp {
                pts.push(domain.uv(p));
            }
            c = ce.next;
            if c == start {
                break;
            }
        }
        if pts.len() >= 3 {
            // close the ring
            if pts.first() != pts.last() {
                pts.push(pts[0]);
            }
            loops.push(pts);
        }
    }
    if !loops.is_empty() {
        return loops;
    }
    // No real loops (analytic template faces) — synthesise the boundary from
    // the FaceExtent so the face still gets a domain rectangle to imprint.
    synthesize_boundary(brep, fid, domain)
}

/// Boundary chains for a face that has no explicit loops, derived from its
/// [`FaceExtent`] template (cylinder rims, disk circle).
fn synthesize_boundary(brep: &BRep, fid: FaceId, domain: &FaceDomain) -> Vec<Vec<(f64, f64)>> {
    let ring = |c: &cadcore_geom::Circle3| -> Vec<(f64, f64)> {
        let mut pts: Vec<(f64, f64)> = (0..64)
            .map(|k| domain.uv(c.point_at(k as f64 / 64.0 * std::f64::consts::TAU)))
            .collect();
        pts.push(pts[0]);
        pts
    };
    match &brep.faces[fid].extent {
        FaceExtent::Cylinder { start, end, .. } => {
            let mut out = Vec::new();
            for b in [start, end] {
                if let FaceBoundary::Circle(c) = b {
                    out.push(ring(c));
                }
            }
            out
        }
        FaceExtent::Disk { radius } => {
            if let FaceGeom::Plane(pl) = &brep.faces[fid].geom {
                let c = cadcore_geom::Circle3::new(pl.frame.origin, pl.normal(), *radius);
                vec![ring(&c)]
            } else {
                Vec::new()
            }
        }
        // FULL torus (TorusBand, non-periodic [0,2π]²): cut at BOTH seams.
        // φ-seam = outer-equator circle at v=0 and v=2π; θ-seam = the minor
        // circle at u=0 and u=2π.  Each seam's two copies are the SAME 3-D
        // circle, so they weld and close the donut in both θ and φ.
        FaceExtent::TorusFillet { start_circle, end_circle }
            if (start_circle.frame.origin - end_circle.frame.origin).length() < 1e-6 =>
        {
            // Built DIRECTLY in (θ,φ) so all four rectangle corners are the
            // single 3-D seam point (θ=0,φ=0) and align exactly.  φ-seam = the
            // φ=0/2π circles (constant v); θ-seam = the θ=0/2π circles (const u).
            let tau = std::f64::consts::TAU;
            let thetas: Vec<f64> = (0..=96).map(|k| k as f64 / 96.0 * tau).collect();
            let phis: Vec<f64> = (0..=64).map(|k| k as f64 / 64.0 * tau).collect();
            vec![
                thetas.iter().map(|&th| (th, 0.0)).collect(),
                thetas.iter().map(|&th| (th, tau)).collect(),
                phis.iter().map(|&ph| (0.0, ph)).collect(),
                phis.iter().map(|&ph| (tau, ph)).collect(),
            ]
        }
        _ => Vec::new(),
    }
}
