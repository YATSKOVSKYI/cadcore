//! Point-in-solid containment for an *arbitrary* B-Rep, by ray casting.
//!
//! This is the classification heart of the general union: to decide whether a
//! face-piece of A survives, we ask whether a point just inside it lies inside
//! solid B (swallowed → dropped) or outside (kept).  Unlike the scaffold
//! engine's `CompositeTube::contains`, this makes no assumption about the
//! solid's shape — it works on any closed shell of analytic faces.
//!
//! Method: shoot a ray from the query point in a pseudo-random direction and
//! count how many times it crosses the shell's boundary.  Odd ⇒ inside.  Each
//! face contributes a crossing when the ray meets its carrier surface *within
//! the face's trimmed region*.  Rays that graze a surface tangentially or pass
//! suspiciously close to a face boundary are rejected and re-cast in a fresh
//! direction (the standard robustness trick).

use cadcore_math::{Point3, UnitVec3, Vec3};
use cadcore_topo::{BRep, FaceExtent, FaceGeom, FaceId};

/// Is `p` strictly inside the solid bounded by `faces` (a closed shell)?
///
/// Points exactly on the boundary are reported per the ray parity and should
/// be avoided by probing slightly inside the material — callers classifying a
/// cell lift an interior point, never a boundary point.
pub fn point_in_solid(brep: &BRep, faces: &[FaceId], p: Point3) -> bool {
    // A few quasi-random, well-separated directions; the first that yields a
    // clean (no grazing / no boundary-grazing) count wins.
    const DIRS: [[f64; 3]; 6] = [
        [0.412_3, 0.713_7, 0.566_1],
        [-0.701_2, 0.219_9, 0.678_5],
        [0.337_1, -0.811_4, 0.477_3],
        [0.901_2, 0.121_7, -0.415_9],
        [-0.258_8, -0.532_1, 0.806_0],
        [0.557_7, -0.301_2, -0.773_9],
    ];
    for d in DIRS {
        let dir = UnitVec3::try_from_vec(Vec3::new(d[0], d[1], d[2])).unwrap();
        if let Some(inside) = cast_once(brep, faces, p, dir) {
            return inside;
        }
    }
    // Extremely unlucky: every direction grazed something.  Fall back to the
    // last cast's parity without the cleanliness guard.
    let dir = UnitVec3::try_from_vec(Vec3::new(DIRS[0][0], DIRS[0][1], DIRS[0][2])).unwrap();
    cast_once_unchecked(brep, faces, p, dir)
}

/// One ray cast.  Returns `Some(inside)` if the count was clean, `None` if the
/// ray grazed a surface or skimmed a face boundary (caller re-casts).
fn cast_once(brep: &BRep, faces: &[FaceId], p: Point3, dir: UnitVec3) -> Option<bool> {
    let mut crossings = 0usize;
    for &fid in faces {
        let f = &brep.faces[fid];
        for (t, n) in surface_hits(&f.geom, p, dir) {
            if t <= 1e-9 {
                continue;
            }
            // grazing: ray nearly tangent to the surface here → unreliable.
            if n.dot(dir.as_vec()).abs() < 1e-6 {
                return None;
            }
            let hit = p + dir.as_vec() * t;
            match face_contains_hit(brep, fid, hit) {
                Membership::Inside => crossings += 1,
                Membership::Outside => {}
                Membership::OnBoundary => return None, // skimmed the trim edge
            }
        }
    }
    Some(crossings % 2 == 1)
}

/// Parity-only cast with no cleanliness guard (last-resort fallback).
fn cast_once_unchecked(brep: &BRep, faces: &[FaceId], p: Point3, dir: UnitVec3) -> bool {
    let mut crossings = 0usize;
    for &fid in faces {
        let f = &brep.faces[fid];
        for (t, _n) in surface_hits(&f.geom, p, dir) {
            if t <= 1e-9 {
                continue;
            }
            let hit = p + dir.as_vec() * t;
            if matches!(face_contains_hit(brep, fid, hit), Membership::Inside) {
                crossings += 1;
            }
        }
    }
    crossings % 2 == 1
}

/// Ray × carrier-surface hits: `(t, outward_normal_at_hit)` for every positive
/// root.  The normal is the surface's natural normal (orientation-agnostic;
/// only its angle to the ray matters for the grazing test).
fn surface_hits(geom: &FaceGeom, p: Point3, dir: UnitVec3) -> Vec<(f64, Vec3)> {
    match geom {
        FaceGeom::Plane(pl) => {
            let n = pl.normal().as_vec();
            let denom = n.dot(dir.as_vec());
            if denom.abs() < 1e-12 {
                return vec![];
            }
            let t = (pl.frame.origin - p).dot(n) / denom;
            vec![(t, n)]
        }
        FaceGeom::Cylinder(s) => {
            // distance-to-axis = r, in the plane perpendicular to the axis.
            let ax = s.axis().as_vec();
            let perp = |v: Vec3| v - ax * ax.dot(v);
            let d = perp(dir.as_vec());
            let oc = perp(p - s.frame.origin);
            let a = d.dot(d);
            if a < 1e-18 {
                return vec![]; // ray parallel to axis
            }
            let b = 2.0 * d.dot(oc);
            let c = oc.dot(oc) - s.radius * s.radius;
            roots(a, b, c)
                .into_iter()
                .map(|t| {
                    let hit = p + dir.as_vec() * t;
                    let radial = perp(hit - s.frame.origin);
                    (t, radial.try_normalize().unwrap_or(Vec3::ZERO))
                })
                .collect()
        }
        FaceGeom::Sphere(s) => {
            let oc = p - s.centre;
            let a = dir.as_vec().dot(dir.as_vec());
            let b = 2.0 * oc.dot(dir.as_vec());
            let c = oc.dot(oc) - s.radius * s.radius;
            roots(a, b, c)
                .into_iter()
                .map(|t| {
                    let hit = p + dir.as_vec() * t;
                    (t, (hit - s.centre).try_normalize().unwrap_or(Vec3::ZERO))
                })
                .collect()
        }
        FaceGeom::Torus(tor) => torus_hits(tor, p, dir),
    }
}

/// Real roots of `a t² + b t + c = 0`, ascending.
fn roots(a: f64, b: f64, c: f64) -> Vec<f64> {
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return vec![];
    }
    let s = disc.sqrt();
    let (t0, t1) = ((-b - s) / (2.0 * a), (-b + s) / (2.0 * a));
    if t0 <= t1 {
        vec![t0, t1]
    } else {
        vec![t1, t0]
    }
}

/// Ray × torus by scanning the implicit function for sign changes, then
/// bisecting — robust without a quartic solver.
fn torus_hits(tor: &cadcore_geom::TorusSurf, p: Point3, dir: UnitVec3) -> Vec<(f64, Vec3)> {
    // implicit f(x) in the torus frame: (|q|² + R² − r²)² − 4R²(qx²+qy²), with
    // q = x − centre expressed in (frame.x, frame.y, frame.z).
    let big_r = tor.major_radius;
    let small_r = tor.minor_radius;
    let f = |x: Point3| {
        let w = x - tor.frame.origin;
        let qx = tor.frame.x.dot_vec(w);
        let qy = tor.frame.y.dot_vec(w);
        let qz = tor.frame.z.dot_vec(w);
        let ss = qx * qx + qy * qy + qz * qz + big_r * big_r - small_r * small_r;
        ss * ss - 4.0 * big_r * big_r * (qx * qx + qy * qy)
    };
    // scan t over [0, span]; span = a few tube diameters past the far side.
    let reach = 2.0 * (big_r + small_r)
        + (p - tor.frame.origin).length()
        + 4.0 * small_r;
    let steps = 256;
    let mut hits = Vec::new();
    let mut t_prev = 1e-9;
    let mut f_prev = f(p + dir.as_vec() * t_prev);
    for k in 1..=steps {
        let t = 1e-9 + reach * k as f64 / steps as f64;
        let fc = f(p + dir.as_vec() * t);
        if f_prev == 0.0 || (f_prev < 0.0) != (fc < 0.0) {
            // bracket [t_prev, t]; bisect.
            let (mut lo, mut hi) = (t_prev, t);
            let (mut flo, _fhi) = (f_prev, fc);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                let fm = f(p + dir.as_vec() * mid);
                if (flo < 0.0) != (fm < 0.0) {
                    hi = mid;
                } else {
                    lo = mid;
                    flo = fm;
                }
            }
            let troot = 0.5 * (lo + hi);
            let hit = p + dir.as_vec() * troot;
            // outward normal via project (degenerate loci skipped).
            let an = crate::geom::refine::AnalyticSurface::Torus(*tor).normal_at(hit);
            if let Some(nv) = an {
                hits.push((troot, nv));
            }
        }
        t_prev = t;
        f_prev = fc;
    }
    hits
}

/// Where a ray hit sits relative to a face's trimmed region.
pub(crate) enum Membership {
    Inside,
    Outside,
    OnBoundary,
}

/// Does `hit` (already on the carrier surface) lie within the face's trimmed
/// region?  Uses the analytic [`FaceExtent`] for primitive input faces; falls
/// back to the real loop polygon when the extent is `Trimmed`/`None`.
pub(crate) fn face_contains_hit(brep: &BRep, fid: FaceId, hit: Point3) -> Membership {
    let f = &brep.faces[fid];
    const ON: f64 = 1e-7;
    let inside_outside = |d_in: f64| {
        if d_in.abs() < ON {
            Membership::OnBoundary
        } else if d_in > 0.0 {
            Membership::Inside
        } else {
            Membership::Outside
        }
    };
    match (&f.geom, &f.extent) {
        (FaceGeom::Plane(pl), FaceExtent::Polygon { points }) => {
            point_in_plane_polygon(pl, points, hit)
        }
        (FaceGeom::Plane(pl), FaceExtent::Disk { radius }) => {
            let r = (hit - pl.frame.origin).length();
            inside_outside(radius - r)
        }
        (FaceGeom::Plane(pl), FaceExtent::PlanarBoundary { boundary }) => {
            use cadcore_topo::FaceBoundary;
            match boundary {
                FaceBoundary::Circle(c) => inside_outside(c.radius - (hit - pl.frame.origin).length()),
                FaceBoundary::Ellipse(_) => point_in_loops(brep, fid, hit),
            }
        }
        (FaceGeom::Cylinder(s), FaceExtent::Cylinder { length, .. }) => {
            let ax = s.axis().dot_vec(hit - s.frame.origin);
            // inside the band ⇔ 0 < ax < length
            let d_in = ax.min(*length - ax);
            inside_outside(d_in)
        }
        // Trimmed cylinder band (output of a prior union): a periodic band is
        // NOT a uv polygon — test the axial range recovered from its loops and
        // exclude any LOCALISED window holes.  (Falling back to point_in_loops
        // would treat a rim circle as a containing polygon, which is wrong.)
        (FaceGeom::Cylinder(s), _) => cylinder_band_membership(brep, fid, s, hit),
        (FaceGeom::Sphere(_), _) => Membership::Inside, // a full sphere face
        (FaceGeom::Torus(t), _) => torus_patch_membership(brep, fid, t, hit),
        // Trimmed / unknown templates: defer to the explicit loops.
        _ => point_in_loops(brep, fid, hit),
    }
}

/// Membership for a trimmed cylinder band: axial range from its loops, minus
/// localised window holes.  Distinguishes a window (inner loop spanning a
/// limited u-arc) from a far-rim boundary (inner loop wrapping the full 2π).
fn cylinder_band_membership(
    brep: &BRep,
    fid: FaceId,
    surf: &cadcore_geom::CylSurf,
    hit: Point3,
) -> Membership {
    use std::f64::consts::TAU;
    const ON: f64 = 1e-7;
    let f = &brep.faces[fid];
    let axis = surf.axis();
    let to_uv = |p: Point3| {
        let w = p - surf.frame.origin;
        let u = surf.frame.y.dot_vec(w).atan2(surf.frame.x.dot_vec(w));
        (u, axis.dot_vec(w))
    };
    // axial span of the whole band (all loops).
    let mut lo = f64::MAX;
    let mut hi = f64::MIN;
    let mut loops = Vec::new();
    let mut ids = vec![f.outer_loop];
    ids.extend(f.inner_loops.iter().copied());
    for lid in ids {
        let Some(lp) = brep.loops.get(lid) else { continue };
        if brep.coedges.get(lp.start).is_none() {
            continue;
        }
        let mut uv: Vec<(f64, f64)> = Vec::new();
        let mut c = lp.start;
        // u is unwrapped CONTINUOUSLY along the loop.  Raw atan2 values make a
        // window that happens to sit on the u = ±π seam look like it spans the
        // whole period, so the full-wrap test below would mistake it for a rim
        // and stop excluding it — the hole then reads as solid material.
        let mut prev_u: Option<f64> = None;
        loop {
            let ce = &brep.coedges[c];
            for p in super::aabb::sample_edge(brep, ce.edge) {
                let (raw_u, v) = to_uv(p);
                let u = match prev_u {
                    Some(pu) => raw_u - TAU * ((raw_u - pu) / TAU).round(),
                    None => raw_u,
                };
                prev_u = Some(u);
                lo = lo.min(v);
                hi = hi.max(v);
                uv.push((u, v));
            }
            c = ce.next;
            if c == lp.start {
                break;
            }
        }
        loops.push(uv);
    }
    let (_, qv) = to_uv(hit);
    let d_in = (qv - lo).min(hi - qv);
    if d_in < -ON {
        return Membership::Outside;
    }
    if d_in < ON {
        return Membership::OnBoundary;
    }
    // exclude localised window holes (u-span < ~full period).
    let (qu, qv) = to_uv(hit);
    for uv in &loops {
        if uv.len() < 3 {
            continue;
        }
        let umin = uv.iter().map(|p| p.0).fold(f64::MAX, f64::min);
        let umax = uv.iter().map(|p| p.0).fold(f64::MIN, f64::max);
        if umax - umin > 0.9 * TAU {
            continue; // a full-wrap rim, not a window
        }
        // Shift the window onto the branch nearest the query as ONE rigid
        // polygon (a per-vertex shift would tear a window that spans the seam).
        let umean = uv.iter().map(|p| p.0).sum::<f64>() / uv.len() as f64;
        let k = ((qu - umean) / TAU).round();
        let poly: Vec<(f64, f64)> = uv.iter().map(|&(u, v)| (u + TAU * k, v)).collect();
        if let Membership::Inside = pip(&poly, qu, qv) {
            return Membership::Outside; // inside a window ⇒ off the band
        }
    }
    Membership::Inside
}

/// Membership for a torus patch (fillet): θ within the arc `[lo, hi]`, minus
/// localised window holes.  Boundaries are minor circles (full-φ at constant
/// θ); a window is an inner loop localised in φ.
fn torus_patch_membership(
    brep: &BRep,
    fid: FaceId,
    t: &cadcore_geom::TorusSurf,
    hit: Point3,
) -> Membership {
    use std::f64::consts::TAU;
    const ON: f64 = 1e-7;
    let f = &brep.faces[fid];
    let theta = |p: Point3| {
        let w = p - t.frame.origin;
        t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w))
    };
    let phi = |p: Point3| {
        let w = p - t.frame.origin;
        let h = t.frame.z.dot_vec(w);
        let a = (w - t.frame.z.as_vec() * h).length();
        h.atan2(a - t.major_radius)
    };
    // θ-arc from all loop vertices (unwrapped to the first sample).
    let mut ids = vec![f.outer_loop];
    ids.extend(f.inner_loops.iter().copied());
    let mut loops: Vec<Vec<(f64, f64)>> = Vec::new(); // (θ,φ) per loop
    let mut base = None;
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    for lid in ids {
        let Some(lp) = brep.loops.get(lid) else { continue };
        if brep.coedges.get(lp.start).is_none() {
            continue;
        }
        let mut uv = Vec::new();
        let mut c = lp.start;
        loop {
            let ce = &brep.coedges[c];
            for p in super::aabb::sample_edge(brep, ce.edge) {
                let raw = theta(p);
                let b = *base.get_or_insert(raw);
                let th = raw + TAU * ((b - raw) / TAU).round();
                lo = lo.min(th);
                hi = hi.max(th);
                uv.push((th, phi(p)));
            }
            c = ce.next;
            if c == lp.start {
                break;
            }
        }
        loops.push(uv);
    }
    if hi <= lo {
        return Membership::Inside;
    }
    let b = base.unwrap_or(0.0);
    let raw = theta(hit);
    let qth = raw + TAU * ((b - raw) / TAU).round();
    let qphi = phi(hit);
    let d_in = (qth - lo).min(hi - qth);
    if d_in < -ON {
        return Membership::Outside;
    }
    if d_in < ON {
        return Membership::OnBoundary;
    }
    // exclude localised windows (inner loops not spanning the full φ circle).
    for uv in &loops {
        if uv.len() < 3 {
            continue;
        }
        let pmin = uv.iter().map(|p| p.1).fold(f64::MAX, f64::min);
        let pmax = uv.iter().map(|p| p.1).fold(f64::MIN, f64::max);
        if pmax - pmin > 0.9 * TAU {
            continue; // full-φ minor circle (a boundary), not a window
        }
        let poly: Vec<(f64, f64)> = uv
            .iter()
            .map(|&(th, ph)| (th, ph + TAU * ((qphi - ph) / TAU).round()))
            .collect();
        if let Membership::Inside = pip(&poly, qth, qphi) {
            return Membership::Outside;
        }
    }
    Membership::Inside
}

/// Point-in-polygon for a planar polygon face (project to the plane's uv).
fn point_in_plane_polygon(pl: &cadcore_geom::Plane3, points: &[Point3], hit: Point3) -> Membership {
    let to_uv = |q: Point3| {
        let w = q - pl.frame.origin;
        (pl.frame.x.dot_vec(w), pl.frame.y.dot_vec(w))
    };
    let poly: Vec<(f64, f64)> = points.iter().map(|&q| to_uv(q)).collect();
    let (px, py) = to_uv(hit);
    pip(&poly, px, py)
}

/// Crossing-number point-in-polygon with an on-edge guard.
fn pip(poly: &[(f64, f64)], px: f64, py: f64) -> Membership {
    const ON: f64 = 1e-7;
    let n = poly.len();
    if n < 3 {
        return Membership::Outside;
    }
    let mut inside = false;
    for i in 0..n {
        let (x0, y0) = poly[i];
        let (x1, y1) = poly[(i + 1) % n];
        // distance to the segment for the on-edge guard
        let dx = x1 - x0;
        let dy = y1 - y0;
        let len2 = dx * dx + dy * dy;
        if len2 > 0.0 {
            let t = (((px - x0) * dx + (py - y0) * dy) / len2).clamp(0.0, 1.0);
            let cx = x0 + dx * t;
            let cy = y0 + dy * t;
            if ((px - cx).powi(2) + (py - cy).powi(2)).sqrt() < ON {
                return Membership::OnBoundary;
            }
        }
        if (y0 > py) != (y1 > py) {
            let xint = x0 + (py - y0) / (y1 - y0) * (x1 - x0);
            if px < xint {
                inside = !inside;
            }
        }
    }
    if inside {
        Membership::Inside
    } else {
        Membership::Outside
    }
}

/// Fallback membership via the face's real loops sampled into a uv polygon.
fn point_in_loops(brep: &BRep, fid: FaceId, hit: Point3) -> Membership {
    let f = &brep.faces[fid];
    let dom = FaceUv::from_face(&f.geom);
    let (qu, qv) = dom.uv(hit);
    let sample_loop = |lid: cadcore_topo::LoopId| -> Vec<(f64, f64)> {
        let mut out = Vec::new();
        let Some(lp) = brep.loops.get(lid) else { return out };
        let start = lp.start;
        if brep.coedges.get(start).is_none() {
            return out;
        }
        let mut c = start;
        loop {
            let ce = &brep.coedges[c];
            for p3 in super::aabb::sample_edge(brep, ce.edge) {
                out.push(dom.uv(p3));
            }
            c = ce.next;
            if c == start {
                break;
            }
        }
        dom.unwrap_periodic(&mut out, qu, qv);
        out
    };
    let outer = sample_loop(f.outer_loop);
    if outer.len() < 3 {
        return Membership::Inside; // no usable boundary → treat as covering
    }
    match pip(&outer, qu, qv) {
        Membership::Outside => return Membership::Outside,
        Membership::OnBoundary => return Membership::OnBoundary,
        Membership::Inside => {}
    }
    for &lid in &f.inner_loops {
        let hole = sample_loop(lid);
        if hole.len() >= 3 {
            if let Membership::Inside = pip(&hole, qu, qv) {
                return Membership::Outside; // inside a hole ⇒ outside the face
            }
        }
    }
    Membership::Inside
}

/// Minimal uv mapper for the loop-fallback path (mirrors `arrange::FaceDomain`
/// but without an extent, and with a periodic-unwrap helper).
struct FaceUv {
    geom: FaceGeom,
}

impl FaceUv {
    fn from_face(geom: &FaceGeom) -> Self {
        FaceUv { geom: geom.clone() }
    }
    fn uv(&self, p: Point3) -> (f64, f64) {
        match &self.geom {
            FaceGeom::Plane(pl) => {
                let w = p - pl.frame.origin;
                (pl.frame.x.dot_vec(w), pl.frame.y.dot_vec(w))
            }
            FaceGeom::Cylinder(s) => {
                let w = p - s.frame.origin;
                let u = s.frame.y.dot_vec(w).atan2(s.frame.x.dot_vec(w));
                (u, s.axis().dot_vec(w))
            }
            FaceGeom::Sphere(s) => {
                let w = p - s.centre;
                (w.y.atan2(w.x), (w.z / s.radius).clamp(-1.0, 1.0).asin())
            }
            FaceGeom::Torus(t) => {
                let w = p - t.frame.origin;
                let theta = t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w));
                let h = t.frame.z.dot_vec(w);
                let q = w - t.frame.z.as_vec() * h;
                let phi = h.atan2(q.length() - t.major_radius);
                (theta, phi)
            }
        }
    }
    /// Unwrap a periodic-axis polygon onto the branch nearest the query so
    /// crossing-number tests stay continuous.
    fn unwrap_periodic(&self, poly: &mut [(f64, f64)], qu: f64, qv: f64) {
        use std::f64::consts::TAU;
        let wrap = |val: &mut f64, around: f64| {
            *val += TAU * ((around - *val) / TAU).round();
        };
        match &self.geom {
            FaceGeom::Cylinder(_) | FaceGeom::Sphere(_) => {
                for (u, _) in poly.iter_mut() {
                    wrap(u, qu);
                }
            }
            FaceGeom::Torus(_) => {
                for (u, v) in poly.iter_mut() {
                    wrap(u, qu);
                    wrap(v, qv);
                }
            }
            FaceGeom::Plane(_) => {}
        }
    }
}
