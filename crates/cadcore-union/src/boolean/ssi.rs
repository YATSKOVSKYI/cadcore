//! Surface×surface intersection (SSI) between two faces of different solids.
//!
//! Each candidate face pair (from the broad phase) is intersected here.  The
//! result is a set of 3-D curves *clipped to the region common to both faces*
//! — so the SAME curve is handed to both faces' arrangements and welds into
//! one shared edge at stitch time (README invariant: joints are a single
//! source of truth).
//!
//! Dispatch is on the carrier-surface pair.  Plane×Plane (the box-union case)
//! is exact (a clipped line).  Curved pairs route through the analytic
//! intersectors / the predictor–corrector tracer; they are added incrementally.

use cadcore_geom::{CylPlaneCurve, CylSurf, Plane3, TorusSurf};
use cadcore_math::{Point3, UnitVec3, Vec3};
use cadcore_topo::{BRep, FaceExtent, FaceGeom, FaceId};

/// One intersection curve, as a 3-D polyline (≥2 points).  A straight segment
/// is just its two endpoints; a `closed` curve (e.g. a circle) repeats no
/// point but is understood to loop back to the start.
#[derive(Debug, Clone)]
pub struct SsiCurve {
    /// Polyline points, start → end.
    pub points: Vec<Point3>,
    /// Whether the curve is a closed loop.
    pub closed: bool,
}

/// Intersect face `fa` of solid `a` with face `fb` of solid `b`, returning the
/// shared intersection curve(s) clipped to both faces' regions.
pub fn intersect_faces(a: &BRep, fa: FaceId, b: &BRep, fb: FaceId) -> Vec<SsiCurve> {
    let ga = &a.faces[fa].geom;
    let gb = &b.faces[fb].geom;
    match (ga, gb) {
        (FaceGeom::Plane(pa), FaceGeom::Plane(pb)) => {
            plane_plane(pa, region_polygon(a, fa), pb, region_polygon(b, fb))
        }
        (FaceGeom::Plane(pl), FaceGeom::Cylinder(_)) => {
            let (cy, len) = match cyl_band(b, fb) { Some(x) => x, None => return Vec::new() };
            plane_cylinder(pl, region_polygon(a, fa), &cy, Some(len))
        }
        (FaceGeom::Cylinder(_), FaceGeom::Plane(pl)) => {
            let (cy, len) = match cyl_band(a, fa) { Some(x) => x, None => return Vec::new() };
            plane_cylinder(pl, region_polygon(b, fb), &cy, Some(len))
        }
        (FaceGeom::Cylinder(_), FaceGeom::Cylinder(_)) => {
            let (Some((ca, la)), Some((cb, lb))) = (cyl_band(a, fa), cyl_band(b, fb)) else {
                return Vec::new();
            };
            cyl_cyl(ca, Some(la), cb, Some(lb))
        }
        (FaceGeom::Torus(_), FaceGeom::Cylinder(_)) => {
            let (Some((tor, lo, hi)), Some((cy, _))) = (torus_arc(a, fa), cyl_band(b, fb)) else {
                return Vec::new();
            };
            torus_cyl(&tor, lo, hi, &cy)
        }
        (FaceGeom::Cylinder(_), FaceGeom::Torus(_)) => {
            let (Some((cy, _)), Some((tor, lo, hi))) = (cyl_band(a, fa), torus_arc(b, fb)) else {
                return Vec::new();
            };
            torus_cyl(&tor, lo, hi, &cy)
        }
        (FaceGeom::Torus(_), FaceGeom::Torus(_)) => {
            let (Some((ta, alo, ahi)), Some((tb, blo, bhi))) = (torus_arc(a, fa), torus_arc(b, fb))
            else {
                return Vec::new();
            };
            torus_torus(&ta, alo, ahi, &tb, blo, bhi)
        }
        (FaceGeom::Plane(pl), FaceGeom::Torus(_)) if std::env::var("CADCORE_PT").is_ok() => {
            let Some((tor, lo, hi)) = torus_arc(b, fb) else { return Vec::new() };
            plane_torus(pl, region_polygon(a, fa), &tor, lo, hi)
        }
        (FaceGeom::Torus(_), FaceGeom::Plane(pl)) if std::env::var("CADCORE_PT").is_ok() => {
            let Some((tor, lo, hi)) = torus_arc(a, fa) else { return Vec::new() };
            plane_torus(pl, region_polygon(b, fb), &tor, lo, hi)
        }
        // Remaining curved pair: none left for the scaffold primitives.
        _ => Vec::new(),
    }
}

/// Torus×Torus intersection via the tracer, seeded by sampling torus A's arc
/// surface and refining onto both tori.  Exact, smooth closed loops shared by
/// both faces; de-duplicated by centroid.
fn torus_torus(
    ta: &TorusSurf,
    alo: f64,
    ahi: f64,
    tb: &TorusSurf,
    blo: f64,
    bhi: f64,
) -> Vec<SsiCurve> {
    use crate::geom::intersect::{intersection_point_near, trace_closed_loop, TraceOptions};
    use crate::geom::refine::AnalyticSurface;
    let sa = AnalyticSurface::Torus(*ta);
    let sb = AnalyticSurface::Torus(*tb);
    // Two doubly-curved tori make the cyclic-projection corrector less stable
    // than torus×cylinder; a smaller step keeps the predictor near the curve.
    let opts = TraceOptions { step: 0.02, ..TraceOptions::default() };
    let mut out: Vec<SsiCurve> = Vec::new();
    // de-dup signature: (centroid, mean radius from centroid).  Two COAXIAL
    // intersection circles share a centroid (on the axis) but differ in radius,
    // so centroid alone would wrongly merge them.
    let mut sigs: Vec<(Point3, f64)> = Vec::new();
    // sample torus A's surface over its arc × full minor circle for seeds.
    let m = 48;
    let n = 32;
    for i in 0..=m {
        let theta = alo + (ahi - alo) * (i as f64 / m as f64);
        for j in 0..n {
            let phi = std::f64::consts::TAU * (j as f64 / n as f64);
            let p = ta.point_at(theta, phi);
            if sb.distance(p) > 0.5 * ta.minor_radius {
                continue; // far from torus B — not near an intersection
            }
            let Some(seed) = intersection_point_near(&sa, &sb, p, 1e-12) else { continue };
            let Some(curve) = trace_closed_loop(&sa, &sb, seed, &opts) else { continue };
            if curve.points.len() < 3 {
                continue;
            }
            let c = curve.points.iter().fold(Point3::new(0.0, 0.0, 0.0), |a, &q| {
                a + (q - Point3::new(0.0, 0.0, 0.0)) * (1.0 / curve.points.len() as f64)
            });
            let mr = curve.points.iter().map(|&q| (q - c).length()).sum::<f64>() / curve.points.len() as f64;
            let tol = 0.5 * ta.minor_radius;
            if sigs.iter().any(|&(q, r)| (q - c).length() < tol && (r - mr).abs() < tol) {
                continue;
            }
            sigs.push((c, mr));
            // Clip the traced loop to BOTH tori's θ-arcs: a crossing keeps the
            // whole closed loop, while coaxial/stacked tori meet in a circle
            // that wraps the full major angle — each θ-arc face keeps only an
            // open arc of it.
            out.extend(clip_loop_to_arcs(&curve.points, ta, alo, ahi, tb, blo, bhi));
        }
    }
    out
}

/// Torus×Cylinder intersection.  Marching squares locates each loop and seeds
/// the predictor–corrector tracer, which returns an EXACT, smooth loop (≤1e-9
/// on both surfaces) — the raw marching contour is too jagged in uv and
/// shatters the arrangement.  Loops are de-duplicated by centroid.
fn torus_cyl(tor: &TorusSurf, theta_lo: f64, theta_hi: f64, cy: &CylSurf) -> Vec<SsiCurve> {
    use crate::geom::intersect::{intersection_point_near, trace_closed_loop, TraceOptions};
    use crate::geom::refine::AnalyticSurface;
    let sa = AnalyticSurface::Torus(*tor);
    let sb = AnalyticSurface::Cylinder(*cy);
    let opts = TraceOptions::default();

    let mut out: Vec<SsiCurve> = Vec::new();
    let mut centres: Vec<Point3> = Vec::new();
    for rough in cadcore_geom::torus_cyl_intersection(tor, cy, theta_lo, theta_hi, 96) {
        if rough.points.len() < 3 {
            continue;
        }
        // seed from a rough point, refined onto both surfaces, then trace.
        let Some(seed) = rough
            .points
            .iter()
            .find_map(|&p| intersection_point_near(&sa, &sb, p, 1e-12))
        else {
            continue;
        };
        let Some(curve) = trace_closed_loop(&sa, &sb, seed, &opts) else {
            continue;
        };
        if curve.points.len() < 3 {
            continue;
        }
        // de-dup: marching squares yields several rough loops that seed the
        // SAME exact loop (traced with slightly different sampling), so compare
        // centroids with a tolerance well below the spacing of the distinct
        // entry/exit loops (~a tube diameter) but above the sampling jitter.
        let c = curve.points.iter().fold(Point3::new(0.0, 0.0, 0.0), |a, &p| {
            a + (p - Point3::new(0.0, 0.0, 0.0)) * (1.0 / curve.points.len() as f64)
        });
        if centres.iter().any(|&q| (q - c).length() < 0.5 * tor.minor_radius) {
            continue;
        }
        centres.push(c);
        out.push(SsiCurve { points: curve.points, closed: true });
    }
    out
}

/// Plane×Torus intersection (an end cap cutting a neighbouring torus tube):
/// trace the planar section of the torus, clipped to the torus's θ-arc AND the
/// cap's disk polygon.  Open arcs where it exits either.
fn plane_torus(plane: &Plane3, poly: Option<Vec<Point3>>, tor: &TorusSurf, lo: f64, hi: f64) -> Vec<SsiCurve> {
    use crate::geom::intersect::{intersection_point_near, trace_closed_loop, TraceOptions};
    use crate::geom::refine::AnalyticSurface;
    let Some(poly) = poly else { return Vec::new() };
    let sa = AnalyticSurface::Plane(*plane);
    let sb = AnalyticSurface::Torus(*tor);
    let opts = TraceOptions { step: 0.02, ..TraceOptions::default() };
    let to_uv = |q: Point3| {
        let w = q - plane.frame.origin;
        (plane.frame.x.dot_vec(w), plane.frame.y.dot_vec(w))
    };
    let poly2: Vec<(f64, f64)> = poly.iter().map(|&q| to_uv(q)).collect();
    let theta_t = |p: Point3| {
        let w = p - tor.frame.origin;
        tor.frame.y.dot_vec(w).atan2(tor.frame.x.dot_vec(w))
    };
    let in_arc = |th: f64| (th - lo).rem_euclid(std::f64::consts::TAU) <= (hi - lo) + 1e-9;

    let mut out: Vec<SsiCurve> = Vec::new();
    let mut sigs: Vec<(Point3, f64)> = Vec::new();
    let m = 64;
    let nn = 48;
    for i in 0..=m {
        let theta = lo + (hi - lo) * (i as f64 / m as f64);
        for j in 0..nn {
            let phi = std::f64::consts::TAU * (j as f64 / nn as f64);
            let p = tor.point_at(theta, phi);
            if sa.distance(p) > 0.3 * tor.minor_radius || !point_in_polygon_2d(&poly2, to_uv(p)) {
                continue;
            }
            let Some(seed) = intersection_point_near(&sa, &sb, p, 1e-12) else { continue };
            let Some(curve) = trace_closed_loop(&sa, &sb, seed, &opts) else { continue };
            if curve.points.len() < 3 {
                continue;
            }
            let c = curve.points.iter().fold(Point3::new(0.0, 0.0, 0.0), |a, &q| a + (q - Point3::new(0.0, 0.0, 0.0)) * (1.0 / curve.points.len() as f64));
            let mr = curve.points.iter().map(|&q| (q - c).length()).sum::<f64>() / curve.points.len() as f64;
            let tol = 0.5 * tor.minor_radius;
            if sigs.iter().any(|&(q, r)| (q - c).length() < tol && (r - mr).abs() < tol) {
                continue;
            }
            sigs.push((c, mr));
            // clip the traced section to the torus θ-arc AND the cap disk.
            let n = curve.points.len();
            let inside: Vec<bool> = curve.points.iter().map(|&p| in_arc(theta_t(p)) && point_in_polygon_2d(&poly2, to_uv(p))).collect();
            if inside.iter().all(|&b| b) {
                out.push(SsiCurve { points: curve.points, closed: true });
                continue;
            }
            if inside.iter().all(|&b| !b) {
                continue;
            }
            let Some(start) = (0..n).find(|&k| inside[k] && !inside[(k + n - 1) % n]) else { continue };
            let mut cur: Vec<Point3> = Vec::new();
            for k in 0..=n {
                let idx = (start + k) % n;
                if inside[idx] {
                    cur.push(curve.points[idx]);
                } else if cur.len() >= 2 {
                    out.push(SsiCurve { points: std::mem::take(&mut cur), closed: false });
                } else {
                    cur.clear();
                }
            }
        }
    }
    out
}

/// Keep the part of a closed loop where both tori's major angle θ lies in their
/// arc.  Fully-inside ⇒ one closed curve; otherwise the inside runs as open
/// arcs.
fn clip_loop_to_arcs(
    pts: &[Point3],
    ta: &TorusSurf,
    alo: f64,
    ahi: f64,
    tb: &TorusSurf,
    blo: f64,
    bhi: f64,
) -> Vec<SsiCurve> {
    let theta = |t: &TorusSurf, p: Point3| {
        let w = p - t.frame.origin;
        t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w))
    };
    let in_arc = |th: f64, lo: f64, hi: f64| (th - lo).rem_euclid(std::f64::consts::TAU) <= (hi - lo) + 1e-9;
    let inside_p = |p: Point3| in_arc(theta(ta, p), alo, ahi) && in_arc(theta(tb, p), blo, bhi);
    let n = pts.len();
    let inside: Vec<bool> = pts.iter().map(|&p| inside_p(p)).collect();
    if inside.iter().all(|&b| b) {
        return vec![SsiCurve { points: pts.to_vec(), closed: true }];
    }
    if inside.iter().all(|&b| !b) {
        return Vec::new();
    }
    let start = (0..n).find(|&i| inside[i] && !inside[(i + n - 1) % n]);
    let Some(start) = start else { return Vec::new() };
    let mut arcs = Vec::new();
    let mut cur: Vec<Point3> = Vec::new();
    for k in 0..=n {
        let i = (start + k) % n;
        if inside[i] {
            cur.push(pts[i]);
        } else if !cur.is_empty() {
            if cur.len() >= 2 {
                arcs.push(SsiCurve { points: std::mem::take(&mut cur), closed: false });
            } else {
                cur.clear();
            }
        }
    }
    arcs
}

/// The torus arc of a face: `(TorusSurf, theta_lo, theta_hi)` from the
/// `TorusFillet` template or — for a `Trimmed` torus — its loop vertices.
fn torus_arc(brep: &BRep, fid: FaceId) -> Option<(TorusSurf, f64, f64)> {
    let f = &brep.faces[fid];
    let FaceGeom::Torus(surf) = &f.geom else { return None };
    let theta_of = |p: Point3| {
        let w = p - surf.frame.origin;
        surf.frame.y.dot_vec(w).atan2(surf.frame.x.dot_vec(w))
    };
    if let FaceExtent::TorusFillet { start_circle, end_circle } = &f.extent {
        // full torus (start == end circle) ⇒ θ ∈ [0, 2π].
        if (start_circle.frame.origin - end_circle.frame.origin).length() < 1e-6 {
            return Some((*surf, 0.0, std::f64::consts::TAU));
        }
        let lo = theta_of(start_circle.frame.origin);
        let mut hi = theta_of(end_circle.frame.origin);
        hi += std::f64::consts::TAU * ((lo - hi) / std::f64::consts::TAU).round();
        return Some((*surf, lo, hi));
    }
    // Trimmed torus: θ-span from loop vertices (unwrapped to the first sample).
    let pts = super::union::face_loop_points(brep, fid);
    let base = theta_of(*pts.first()?);
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    for p in &pts {
        let raw = theta_of(*p);
        let t = raw + std::f64::consts::TAU * ((base - raw) / std::f64::consts::TAU).round();
        lo = lo.min(t);
        hi = hi.max(t);
    }
    if hi <= lo {
        return None;
    }
    Some((*surf, lo, hi))
}

/// Cylinder×Cylinder intersection via the predictor–corrector tracer (exact to
/// ~1e-9, clipped to both finite bands).  Each closed loop lies on BOTH
/// cylinders and is handed to both — welding into one shared edge ring.
fn cyl_cyl(ca: CylSurf, len_a: Option<f64>, cb: CylSurf, len_b: Option<f64>) -> Vec<SsiCurve> {
    let (Some(len_a), Some(len_b)) = (len_a, len_b) else {
        return Vec::new();
    };
    // PARALLEL tubes (filaments lying on each other): the surfaces meet in two
    // straight RULING lines, not a closed loop — the tracer would find nothing.
    if ca.axis().dot(cb.axis()).abs() > 1.0 - 1e-7 {
        return parallel_cyl_cyl(&ca, len_a, &cb, len_b);
    }
    use crate::geom::composite::{closed_loops_between, CompositeTube};
    use crate::geom::intersect::TraceOptions;
    let ta = CompositeTube::new().with_leg(ca, len_a);
    let tb = CompositeTube::new().with_leg(cb, len_b);
    closed_loops_between(&ta, &tb, &TraceOptions::default())
        .into_iter()
        .filter(|c| c.points.len() >= 3)
        .map(|c| SsiCurve { points: c.points, closed: true })
        .collect()
}

/// Two parallel overlapping cylinders → the two ruling lines, each clipped to
/// the axial overlap of both bands.  Open segments (run rim-to-rim when the
/// bands coincide), shared by both faces.
fn parallel_cyl_cyl(ca: &CylSurf, len_a: f64, cb: &CylSurf, len_b: f64) -> Vec<SsiCurve> {
    let axis = ca.axis();
    let mut out = Vec::new();
    for line in cadcore_geom::cyl_cyl_parallel_lines(ca, cb) {
        // axial coord on each cylinder of the line at t=0, slope is 1 along axis.
        // axial coord on each cylinder is (offset + t); valid t keeps both in
        // their band: t ∈ [max(−a0,−b0), min(len_a−a0, len_b−b0)].
        let a0 = axis.dot_vec(line.origin - ca.frame.origin);
        let b0 = ca.axis().dot_vec(line.origin - cb.frame.origin);
        let lo = (-a0).max(-b0);
        let hi = (len_a - a0).min(len_b - b0);
        if hi - lo > 1e-9 {
            out.push(SsiCurve {
                points: vec![line.origin + axis.as_vec() * lo, line.origin + axis.as_vec() * hi],
                closed: false,
            });
        }
    }
    out
}

/// The cylinder band of a face: a [`CylSurf`] whose `frame.origin` sits at the
/// band's start plus its axial `length`.  Recovered from the `Cylinder` extent
/// template, or — for a `Trimmed` cylinder face (output of a prior union) —
/// from the axial span of its real loop vertices, so the union algebra closes
/// on crossing cylinders too.
fn cyl_band(brep: &BRep, fid: FaceId) -> Option<(CylSurf, f64)> {
    let f = &brep.faces[fid];
    let FaceGeom::Cylinder(surf) = &f.geom else { return None };
    if let FaceExtent::Cylinder { length, .. } = f.extent {
        return Some((*surf, length));
    }
    // walk all loops, project vertices onto the axis to find [min, max] axial.
    let axis = surf.axis();
    let mut lo = f64::MAX;
    let mut hi = f64::MIN;
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
            for p in super::aabb::sample_edge(brep, ce.edge) {
                let t = axis.dot_vec(p - surf.frame.origin);
                lo = lo.min(t);
                hi = hi.max(t);
            }
            c = ce.next;
            if c == start {
                break;
            }
        }
    }
    if hi <= lo {
        return None;
    }
    let shifted = CylSurf::new(surf.frame.origin + axis.as_vec() * lo, axis, surf.radius);
    Some((shifted, hi - lo))
}

/// The planar region of a face as a polygon of 3-D points, if it has one.
///
/// Prefers the explicit `Polygon` template; otherwise walks the face's real
/// outer loop (so the **output of a previous union** — `FaceExtent::Trimmed`
/// planar faces — can feed the next union: the algebra closes).
fn region_polygon(brep: &BRep, fid: FaceId) -> Option<Vec<Point3>> {
    let f = &brep.faces[fid];
    match (&f.geom, &f.extent) {
        (_, FaceExtent::Polygon { points }) => return Some(points.clone()),
        // A disk cap (no real loops): sample its boundary circle as a polygon
        // so it can clip a plane×cylinder conic (e.g. a stacked tube's end cap
        // cutting the neighbouring tube).
        (FaceGeom::Plane(pl), FaceExtent::Disk { radius }) => {
            return Some(
                (0..48)
                    .map(|k| {
                        let a = std::f64::consts::TAU * k as f64 / 48.0;
                        pl.frame.origin + (pl.frame.x * a.cos() + pl.frame.y * a.sin()) * *radius
                    })
                    .collect(),
            );
        }
        _ => {}
    }
    // walk the outer loop into a 3-D polygon
    let f = &brep.faces[fid];
    let lp = brep.loops.get(f.outer_loop)?;
    let start = lp.start;
    brep.coedges.get(start)?;
    let mut pts: Vec<Point3> = Vec::new();
    let mut c = start;
    loop {
        let ce = &brep.coedges[c];
        let mut sp = super::aabb::sample_edge(brep, ce.edge);
        if ce.sense == cadcore_topo::CoEdgeSense::Opposite {
            sp.reverse();
        }
        // append without duplicating the shared vertex between consecutive edges
        for p in sp {
            if pts.last().map_or(true, |&q| (q - p).length() > 1e-9) {
                pts.push(p);
            }
        }
        c = ce.next;
        if c == start {
            break;
        }
    }
    // drop trailing duplicate of the first point
    if pts.len() > 1 && (pts[0] - *pts.last().unwrap()).length() < 1e-9 {
        pts.pop();
    }
    if pts.len() >= 3 {
        Some(pts)
    } else {
        None
    }
}

/// Plane×Plane intersection clipped to both faces' polygons → segment(s).
fn plane_plane(
    pa: &Plane3,
    poly_a: Option<Vec<Point3>>,
    pb: &Plane3,
    poly_b: Option<Vec<Point3>>,
) -> Vec<SsiCurve> {
    let (Some(poly_a), Some(poly_b)) = (poly_a, poly_b) else {
        return Vec::new();
    };
    let na = pa.normal().as_vec();
    let nb = pb.normal().as_vec();
    let dir = na.cross(nb);
    let dlen = dir.length();
    if dlen < 1e-12 {
        return Vec::new(); // parallel / coincident — no transversal line
    }
    let dir = UnitVec3::try_from_vec(dir).unwrap();

    // A point on both planes: solve the 3×3 [na; nb; dir]·x = [da; db; 0].
    let da = na.dot(pa.frame.origin - Point3::new(0.0, 0.0, 0.0));
    let db = nb.dot(pb.frame.origin - Point3::new(0.0, 0.0, 0.0));
    let Some(origin) = three_plane_point(na, da, nb, db, dir.as_vec(), 0.0) else {
        return Vec::new();
    };

    // intervals (in line parameter t) where the line lies inside each polygon
    let ia = line_in_polygon(origin, dir, pa, &poly_a);
    let ib = line_in_polygon(origin, dir, pb, &poly_b);
    let mut out = Vec::new();
    for &(a0, a1) in &ia {
        for &(b0, b1) in &ib {
            let lo = a0.max(b0);
            let hi = a1.min(b1);
            if hi - lo > 1e-9 {
                out.push(SsiCurve {
                    points: vec![
                        origin + dir.as_vec() * lo,
                        origin + dir.as_vec() * hi,
                    ],
                    closed: false,
                });
            }
        }
    }
    out
}

/// Plane×Cylinder intersection (the conic) clipped to the planar face's
/// polygon and the cylinder's axial band.
///
/// Currently handles the **fully-contained closed conic** (circle / ellipse
/// lying entirely inside the plane polygon and within `[0, length]`): the case
/// of a round peg passing through a plate face.  Partial arcs (conic clipped
/// by a polygon edge) and the plane-parallel `Lines` case are deferred —
/// returning empty there is incomplete, never wrong.
fn plane_cylinder(
    plane: &Plane3,
    poly: Option<Vec<Point3>>,
    cyl: &CylSurf,
    length: Option<f64>,
) -> Vec<SsiCurve> {
    let (Some(poly), Some(length)) = (poly, length) else {
        return Vec::new();
    };
    let conic = cadcore_geom::cyl_plane_intersection(cyl, plane);
    let pts: Vec<Point3> = match &conic {
        CylPlaneCurve::Circle(c) => (0..96).map(|i| c.point_at(i as f64 / 96.0 * std::f64::consts::TAU)).collect(),
        CylPlaneCurve::Ellipse(e) => (0..96).map(|i| e.point_at(i as f64 / 96.0 * std::f64::consts::TAU)).collect(),
        // Lines (plane ∥ axis) — deferred.
        _ => return Vec::new(),
    };
    let to_uv = |q: Point3| {
        let w = q - plane.frame.origin;
        (plane.frame.x.dot_vec(w), plane.frame.y.dot_vec(w))
    };
    let poly2: Vec<(f64, f64)> = poly.iter().map(|&q| to_uv(q)).collect();
    let axis = cyl.axis();
    // a sample is "in" when it lies inside the plane polygon AND the cyl band.
    let inside: Vec<bool> = pts
        .iter()
        .map(|&p| {
            let v = axis.dot_vec(p - cyl.frame.origin);
            v >= -1e-9 && v <= length + 1e-9 && point_in_polygon_2d(&poly2, to_uv(p))
        })
        .collect();
    if inside.iter().all(|&b| b) {
        return vec![SsiCurve { points: pts, closed: true }]; // fully contained loop
    }
    if inside.iter().all(|&b| !b) {
        return Vec::new(); // no overlap
    }
    // PARTIAL: keep the maximal runs of in-samples as open arcs (the conic
    // clipped to the polygon / band — e.g. a stacked tube's end cap cutting its
    // neighbour).  The shared sample points weld on both faces.
    let n = pts.len();
    let start = (0..n).find(|&i| inside[i] && !inside[(i + n - 1) % n]);
    let Some(start) = start else { return Vec::new() };
    let mut arcs = Vec::new();
    let mut cur: Vec<Point3> = Vec::new();
    for k in 0..=n {
        let i = (start + k) % n;
        if inside[i] {
            cur.push(pts[i]);
        } else if !cur.is_empty() {
            if cur.len() >= 2 {
                arcs.push(SsiCurve { points: std::mem::take(&mut cur), closed: false });
            } else {
                cur.clear();
            }
        }
    }
    arcs
}

/// Solve the linear system of three plane equations `ni·x = di` for x.
fn three_plane_point(n0: Vec3, d0: f64, n1: Vec3, d1: f64, n2: Vec3, d2: f64) -> Option<Point3> {
    // Cramer's rule on the 3×3 [n0;n1;n2].
    let det = n0.dot(n1.cross(n2));
    if det.abs() < 1e-12 {
        return None;
    }
    let dvec = Vec3::new(d0, d1, d2);
    // x = (d0 (n1×n2) + d1 (n2×n0) + d2 (n0×n1)) / det
    let x = (n1.cross(n2) * dvec.x + n2.cross(n0) * dvec.y + n0.cross(n1) * dvec.z) * (1.0 / det);
    Some(Point3::new(0.0, 0.0, 0.0) + x)
}

/// Parameter intervals `[t_lo, t_hi]` along the line `origin + t·dir` where it
/// lies inside the planar polygon (projected into the plane's uv).
fn line_in_polygon(origin: Point3, dir: UnitVec3, plane: &Plane3, poly3: &[Point3]) -> Vec<(f64, f64)> {
    let to_uv = |q: Point3| {
        let w = q - plane.frame.origin;
        (plane.frame.x.dot_vec(w), plane.frame.y.dot_vec(w))
    };
    let poly: Vec<(f64, f64)> = poly3.iter().map(|&q| to_uv(q)).collect();
    let o2 = to_uv(origin);
    // dir projected into the plane's uv (line lies in the plane, so exact)
    let d2 = (plane.frame.x.dot_vec(dir.as_vec()), plane.frame.y.dot_vec(dir.as_vec()));

    // crossings of the line with each polygon edge → line parameters t
    let mut ts: Vec<f64> = Vec::new();
    let n = poly.len();
    for i in 0..n {
        let (ax, ay) = poly[i];
        let (bx, by) = poly[(i + 1) % n];
        let ex = bx - ax;
        let ey = by - ay;
        // O + t D = A + s E  →  solve for t, s
        let det = d2.0 * (-ey) - (-ex) * d2.1;
        if det.abs() < 1e-12 {
            continue; // line parallel to this edge
        }
        let rx = ax - o2.0;
        let ry = ay - o2.1;
        let t = (rx * (-ey) - (-ex) * ry) / det;
        let s = (d2.0 * ry - d2.1 * rx) / det;
        if s >= -1e-9 && s <= 1.0 + 1e-9 {
            ts.push(t);
        }
    }
    ts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ts.dedup_by(|a, b| (*a - *b).abs() < 1e-9);

    let mut intervals = Vec::new();
    for w in ts.windows(2) {
        let mid = 0.5 * (w[0] + w[1]);
        let mp = (o2.0 + d2.0 * mid, o2.1 + d2.1 * mid);
        if point_in_polygon_2d(&poly, mp) {
            intervals.push((w[0], w[1]));
        }
    }
    intervals
}

/// Crossing-number point-in-polygon in 2-D.
fn point_in_polygon_2d(poly: &[(f64, f64)], (px, py): (f64, f64)) -> bool {
    let n = poly.len();
    let mut inside = false;
    for i in 0..n {
        let (x0, y0) = poly[i];
        let (x1, y1) = poly[(i + 1) % n];
        if (y0 > py) != (y1 > py) {
            let xint = x0 + (py - y0) / (y1 - y0) * (x1 - x0);
            if px < xint {
                inside = !inside;
            }
        }
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean::tests_support::axis_box;

    /// Two overlapping axis boxes: the SSI of a face pair that genuinely
    /// crosses must be a segment; non-crossing pairs give nothing.
    #[test]
    fn box_box_plane_plane_segment() {
        let mut a = BRep::new();
        let fa = axis_box(&mut a, Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 2.0));
        let mut b = BRep::new();
        // shifted so it pokes into A from +x/+y/+z corner
        let fb = axis_box(&mut b, Point3::new(1.0, 1.0, 1.0), Point3::new(3.0, 3.0, 3.0));

        // collect all SSI segments across all candidate pairs
        let pairs = crate::boolean::candidate_pairs(&a, &fa, &b, &fb, 1e-6);
        assert!(!pairs.is_empty());
        let mut total = 0usize;
        let mut seg_len_sum = 0.0;
        for (pa, pb) in &pairs {
            for c in intersect_faces(&a, *pa, &b, *pb) {
                total += 1;
                seg_len_sum += (c.points[1] - c.points[0]).length();
                // every SSI point lies on BOTH carrier planes
                for &p in &c.points {
                    if let (FaceGeom::Plane(pla), FaceGeom::Plane(plb)) =
                        (&a.faces[*pa].geom, &b.faces[*pb].geom)
                    {
                        assert!(pla.signed_distance(p).abs() < 1e-9);
                        assert!(plb.signed_distance(p).abs() < 1e-9);
                    }
                }
            }
        }
        assert!(total > 0, "overlapping boxes produce SSI segments");
        assert!(seg_len_sum > 0.0);
    }

    /// Faces of separated boxes never intersect.
    #[test]
    fn disjoint_boxes_no_ssi() {
        let mut a = BRep::new();
        let fa = axis_box(&mut a, Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
        let mut b = BRep::new();
        let fb = axis_box(&mut b, Point3::new(5.0, 5.0, 5.0), Point3::new(6.0, 6.0, 6.0));
        for &pa in &fa {
            for &pb in &fb {
                assert!(intersect_faces(&a, pa, &b, pb).is_empty());
            }
        }
    }
}
