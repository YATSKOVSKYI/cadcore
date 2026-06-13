//! Generic numeric **surface–surface intersection** (SSI) for any pair of
//! [`ParamSurface`]s (see `docs/swept-tube-ssi-pipeline.md`).
//!
//! Pipeline: seed search → predictor–corrector tracing.
//!
//! * **Corrector** — given a guess `(u₁,v₁,u₂,v₂)`, drive the 3-D residual
//!   `F = S₁(u₁,v₁) − S₂(u₂,v₂)` to zero with a Levenberg–Marquardt step on the
//!   4 parameters (`JᵀJ` is 4×4, rank ≤ 3 along the curve; `λI` regularises the
//!   free tangent direction so the corrector projects *onto* the curve).
//! * **Predictor** — at a converged point the 3-D curve tangent is
//!   `d = n₁ × n₂`; advance the parameters along `d` projected into each
//!   surface's tangent plane, then re-correct.
//!
//! Output curves are sampled polylines, each marked closed/open.

use cadcore_math::{Point3, Vec3};

use crate::swept::{ParamSurface, SweptTubeSurface};
use crate::IntersectionPolyline;

const TAU: f64 = std::f64::consts::PI * 2.0;

/// **Robust** swept-tube ∩ swept-tube intersection by marching squares.
///
/// The predictor–corrector [`surface_surface_intersection`] is elegant but can
/// wander on shallow / curved / large-coordinate crossings (it runs away
/// without closing).  This variant grids surface `a`'s `(u, v)` domain, samples
/// the signed distance to `b`'s solid (`g = dist(a(u,v), b_centerline) − r_b`),
/// and extracts the `g = 0` contour by marching squares — no tracing, no
/// runaway.  `v` is periodic.  Returns closed loops in `a`'s parameter domain,
/// lifted to 3-D (points lie on `a` exactly, on `b` to grid tolerance).
pub fn swept_swept_intersection(
    a: &SweptTubeSurface,
    b: &SweptTubeSurface,
    n_u: usize,
    n_v: usize,
) -> Vec<IntersectionPolyline> {
    let (u0, u1) = (0.0, a.length());
    let nu = n_u.max(8);
    let nv = n_v.max(8);

    // Signed distance from a(u,v) to b's solid (negative = inside b).
    let g = |u: f64, v: f64| -> f64 {
        b.signed_distance(a.point(u.clamp(u0, u1), v.rem_euclid(TAU)))
            .unwrap_or(1.0)
    };

    // Sample g on a (nu+1) × nv grid (v wraps).
    let uu = |i: usize| u0 + (u1 - u0) * i as f64 / nu as f64;
    let vv = |j: usize| TAU * j as f64 / nv as f64;
    let mut grid = vec![0.0f64; (nu + 1) * nv];
    for i in 0..=nu {
        for j in 0..nv {
            grid[i * nv + j] = g(uu(i), vv(j));
        }
    }
    let at = |i: usize, j: usize| grid[i * nv + (j % nv)];

    // Marching squares -> undirected (u,v) segments.
    // Linear-interpolate the crossing on each cell edge.
    let mut segs: Vec<((f64, f64), (f64, f64))> = Vec::new();
    let interp = |ua: f64, va: f64, ga: f64, ub: f64, vb: f64, gb: f64| -> (f64, f64) {
        let denom = ga - gb;
        let t = if denom.abs() < 1e-18 {
            0.5
        } else {
            (ga / denom).clamp(0.0, 1.0)
        };
        (ua + (ub - ua) * t, va + (vb - va) * t)
    };
    let crosses = |ga: f64, gb: f64| -> bool {
        (ga <= 0.0 && gb >= 0.0 || ga >= 0.0 && gb <= 0.0) && (ga - gb).abs() > 1e-14
    };
    for i in 0..nu {
        for j in 0..nv {
            let (jp, jn) = (j, (j + 1) % nv);
            let g00 = at(i, jp);
            let g10 = at(i + 1, jp);
            let g11 = at(i + 1, jn);
            let g01 = at(i, jn);
            let (u_lo, u_hi) = (uu(i), uu(i + 1));
            let (v_lo, v_hi) = (vv(j), vv(j) + TAU / nv as f64);
            let mut pts = Vec::new();
            if crosses(g00, g10) {
                pts.push(interp(u_lo, v_lo, g00, u_hi, v_lo, g10));
            }
            if crosses(g10, g11) {
                pts.push(interp(u_hi, v_lo, g10, u_hi, v_hi, g11));
            }
            if crosses(g11, g01) {
                pts.push(interp(u_hi, v_hi, g11, u_lo, v_hi, g01));
            }
            if crosses(g01, g00) {
                pts.push(interp(u_lo, v_hi, g01, u_lo, v_lo, g00));
            }
            if pts.len() == 2 {
                segs.push((pts[0], pts[1]));
            } else if pts.len() == 4 {
                segs.push((pts[0], pts[1]));
                segs.push((pts[2], pts[3]));
            }
        }
    }
    if segs.is_empty() {
        return vec![];
    }

    // Link segments into contour chains by endpoint matching (quantised key).
    let key = |p: (f64, f64)| -> (i64, i64) {
        let vv = if (p.1.rem_euclid(TAU) - TAU).abs() < 1e-9 {
            0.0
        } else {
            p.1.rem_euclid(TAU)
        };
        ((p.0 * 1e7).round() as i64, (vv * 1e7).round() as i64)
    };
    use std::collections::HashMap;
    let mut adj: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (s, seg) in segs.iter().enumerate() {
        adj.entry(key(seg.0)).or_default().push(s);
        adj.entry(key(seg.1)).or_default().push(s);
    }
    let mut used = vec![false; segs.len()];
    let mut loops = Vec::new();
    for start in 0..segs.len() {
        if used[start] {
            continue;
        }
        let mut chain: Vec<(f64, f64)> = vec![segs[start].0, segs[start].1];
        used[start] = true;
        // Grow tail.
        loop {
            let tail = *chain.last().unwrap();
            let mut advanced = false;
            if let Some(cands) = adj.get(&key(tail)) {
                for &s in cands {
                    if used[s] {
                        continue;
                    }
                    let (p0, p1) = segs[s];
                    let next = if key(p0) == key(tail) {
                        p1
                    } else if key(p1) == key(tail) {
                        p0
                    } else {
                        continue;
                    };
                    chain.push(next);
                    used[s] = true;
                    advanced = true;
                    break;
                }
            }
            if !advanced {
                break;
            }
        }
        // Grow head.
        loop {
            let head = chain[0];
            let mut advanced = false;
            if let Some(cands) = adj.get(&key(head)) {
                for &s in cands {
                    if used[s] {
                        continue;
                    }
                    let (p0, p1) = segs[s];
                    let next = if key(p0) == key(head) {
                        p1
                    } else if key(p1) == key(head) {
                        p0
                    } else {
                        continue;
                    };
                    chain.insert(0, next);
                    used[s] = true;
                    advanced = true;
                    break;
                }
            }
            if !advanced {
                break;
            }
        }

        if chain.len() >= 4 {
            let refined: Vec<((f64, f64), Point3)> = chain
                .iter()
                .map(|&(u, v)| refine_uv_on_zero(&g, u, v, u0, u1))
                .filter_map(|(u, v)| refine_point_to_both_surfaces(a, b, u, v))
                .collect();
            if refined.len() < 4 {
                continue;
            }
            let refined_uv: Vec<(f64, f64)> = refined.iter().map(|(uv, _)| *uv).collect();
            let closed = uv_closed(refined_uv[0], *refined_uv.last().unwrap(), u1 - u0)
                || key(refined_uv[0]) == key(*refined_uv.last().unwrap());
            let open_to_boundary = !closed
                && is_domain_boundary(refined_uv[0], u0, u1)
                && is_domain_boundary(*refined_uv.last().unwrap(), u0, u1);
            if !closed && !open_to_boundary {
                continue;
            }

            let pts: Vec<Point3> = refined.iter().map(|(_, p)| *p).collect();
            if polyline_length(&pts, closed) > a.radius() * 0.05 {
                loops.push(IntersectionPolyline {
                    points: pts,
                    closed,
                });
            }
        }
    }
    loops
}

fn refine_point_to_both_surfaces(
    a: &SweptTubeSurface,
    b: &SweptTubeSurface,
    u: f64,
    v: f64,
) -> Option<((f64, f64), Point3)> {
    let p = a.point(u, v);
    let (bu, bv) = b.project_point(p)?;
    let seed = correct(a, b, (u, v, bu, bv), 1e-8)?;
    let pa = a.point(seed.u1, seed.v1);
    let pb = b.point(seed.u2, seed.v2);
    if (pa - pb).length() <= 1e-6 {
        Some(((seed.u1, seed.v1), seed.point))
    } else {
        None
    }
}

fn refine_uv_on_zero(
    g: &dyn Fn(f64, f64) -> f64,
    mut u: f64,
    mut v: f64,
    u0: f64,
    u1: f64,
) -> (f64, f64) {
    let hu = ((u1 - u0) * 1e-5).max(1e-8);
    let hv = 1e-5;
    for _ in 0..20 {
        let gv = g(u, v);
        if gv.abs() < 1e-9 {
            break;
        }
        let dgu = if u <= u0 + hu {
            (g((u + hu).min(u1), v) - gv) / hu
        } else if u >= u1 - hu {
            (gv - g((u - hu).max(u0), v)) / hu
        } else {
            (g(u + hu, v) - g(u - hu, v)) / (2.0 * hu)
        };
        let dgv = (g(u, v + hv) - g(u, v - hv)) / (2.0 * hv);
        let n2 = dgu * dgu + dgv * dgv;
        if n2 < 1e-18 {
            break;
        }

        let du = -gv * dgu / n2;
        let dv = -gv * dgv / n2;
        let max_du = ((u1 - u0) * 0.05).max(1e-8);
        let max_dv = 0.25;
        let scale = (max_du / du.abs().max(1e-18))
            .min(max_dv / dv.abs().max(1e-18))
            .min(1.0);
        let mut accepted = false;
        let mut alpha = scale;
        for _ in 0..8 {
            let nu = (u + du * alpha).clamp(u0, u1);
            let nv = (v + dv * alpha).rem_euclid(TAU);
            if g(nu, nv).abs() < gv.abs() {
                u = nu;
                v = nv;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    (u.clamp(u0, u1), v.rem_euclid(TAU))
}

fn uv_closed(a: (f64, f64), b: (f64, f64), u_span: f64) -> bool {
    let du = (a.0 - b.0).abs();
    let dv = (a.1 - b.1).abs().min(TAU - (a.1 - b.1).abs());
    du < u_span.max(1.0) * 1e-6 && dv < 1e-5
}

fn is_domain_boundary(p: (f64, f64), u0: f64, u1: f64) -> bool {
    let tol = (u1 - u0).max(1.0) * 1e-6;
    p.0 <= u0 + tol || p.0 >= u1 - tol
}

fn polyline_length(pts: &[Point3], closed: bool) -> f64 {
    if pts.len() < 2 {
        return 0.0;
    }
    let mut len = pts.windows(2).map(|w| (w[1] - w[0]).length()).sum::<f64>();
    if closed {
        len += (pts[0] - pts[pts.len() - 1]).length();
    }
    len
}

/// Tunables for the tracer.
#[derive(Clone, Copy, Debug)]
pub struct SsiOptions {
    /// 3-D step length between traced points.
    pub step: f64,
    /// Convergence tolerance on `|S₁ − S₂|`.
    pub tol: f64,
    /// Max points per traced curve.
    pub max_pts: usize,
    /// Seed-grid resolution in `u` (per surface).
    pub seed_u: usize,
    /// Seed-grid resolution in `v` (per surface).
    pub seed_v: usize,
}

impl Default for SsiOptions {
    fn default() -> Self {
        Self {
            step: 0.05,
            tol: 1e-9,
            max_pts: 4000,
            seed_u: 64,
            seed_v: 24,
        }
    }
}

/// Intersect two parametric surfaces.  Returns the intersection curve(s).
pub fn surface_surface_intersection(
    a: &dyn ParamSurface,
    b: &dyn ParamSurface,
    opts: &SsiOptions,
) -> Vec<IntersectionPolyline> {
    let seeds = find_seeds(a, b, opts);
    let mut curves: Vec<IntersectionPolyline> = Vec::new();

    for seed in seeds {
        // Skip seeds already covered by a previously-traced curve.
        if curves.iter().any(|c| {
            c.points
                .iter()
                .any(|p| (*p - seed.point).length() < opts.step * 1.5)
        }) {
            continue;
        }
        if let Some(curve) = trace(a, b, seed, opts) {
            if curve.points.len() >= 4 {
                curves.push(curve);
            }
        }
    }
    curves
}

// ── seeds ────────────────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct Seed {
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
    point: Point3,
}

fn find_seeds(a: &dyn ParamSurface, b: &dyn ParamSurface, opts: &SsiOptions) -> Vec<Seed> {
    let (u0, u1) = a.u_domain();
    let (v0, v1) = a.v_domain();
    let mut seeds = Vec::new();
    for i in 0..=opts.seed_u {
        for j in 0..opts.seed_v {
            let u = u0 + (u1 - u0) * i as f64 / opts.seed_u as f64;
            let v = v0 + (v1 - v0) * j as f64 / opts.seed_v as f64;
            let p = a.point(u, v);
            let Some((bu, bv)) = b.project_point(p) else {
                continue;
            };
            let pb = b.point(bu, bv);
            if (p - pb).length() < opts.step {
                // Refine to the curve before storing.
                if let Some(s) = correct(a, b, (u, v, bu, bv), opts.tol) {
                    seeds.push(s);
                }
            }
        }
    }
    seeds
}

// ── corrector ────────────────────────────────────────────────────────────────

/// Levenberg–Marquardt projection of a 4-param guess onto the intersection.
fn correct(
    a: &dyn ParamSurface,
    b: &dyn ParamSurface,
    init: (f64, f64, f64, f64),
    tol: f64,
) -> Option<Seed> {
    let (mut u1, mut v1, mut u2, mut v2) = init;
    let (au0, au1) = a.u_domain();
    let (bu0, bu1) = b.u_domain();
    let lambda = 1e-8;

    for _ in 0..40 {
        let s1 = a.point(u1, v1);
        let s2 = b.point(u2, v2);
        let f = s1 - s2;
        if f.length() < tol {
            let mid = Point3::new(
                0.5 * (s1.x + s2.x),
                0.5 * (s1.y + s2.y),
                0.5 * (s1.z + s2.z),
            );
            return Some(Seed {
                u1,
                v1,
                u2,
                v2,
                point: mid,
            });
        }
        let da = a.derivatives(u1, v1);
        let db = b.derivatives(u2, v2);
        // Columns of J (∂F/∂param): F = S1 - S2.
        let cols = [da.du, da.dv, db.du * -1.0, db.dv * -1.0];
        // JᵀJ + λI  and  −Jᵀf
        let mut m = [[0.0f64; 4]; 4];
        let mut rhs = [0.0f64; 4];
        for i in 0..4 {
            for j in 0..4 {
                m[i][j] = cols[i].dot(cols[j]);
            }
            m[i][i] += lambda;
            rhs[i] = -cols[i].dot(f);
        }
        let Some(d) = solve4(m, rhs) else { return None };
        // Damp overly large steps.
        let dmax = d.iter().fold(0.0f64, |a, &x| a.max(x.abs()));
        let scale = if dmax > 0.5 { 0.5 / dmax } else { 1.0 };
        u1 = (u1 + d[0] * scale).clamp(au0, au1);
        v1 = wrap_v(v1 + d[1] * scale);
        u2 = (u2 + d[2] * scale).clamp(bu0, bu1);
        v2 = wrap_v(v2 + d[3] * scale);
    }
    None
}

// ── tracer ───────────────────────────────────────────────────────────────────

fn trace(
    a: &dyn ParamSurface,
    b: &dyn ParamSurface,
    seed: Seed,
    opts: &SsiOptions,
) -> Option<IntersectionPolyline> {
    let start = seed;
    let mut pts: Vec<Point3> = vec![start.point];
    let (mut u1, mut v1, mut u2, mut v2) = (start.u1, start.v1, start.u2, start.v2);
    let mut prev_dir: Option<Vec3> = None;

    for step_i in 0..opts.max_pts {
        let n1 = a.normal(u1, v1);
        let n2 = b.normal(u2, v2);
        let mut d = n1.cross(n2);
        if d.length() < 1e-9 {
            break; // tangent / degenerate
        }
        d = d.normalize();
        if let Some(pd) = prev_dir {
            if d.dot(pd) < 0.0 {
                d = d * -1.0;
            }
        }

        // Predict in parameter space: project the 3-D step into each tangent plane.
        let da = a.derivatives(u1, v1);
        let db = b.derivatives(u2, v2);
        let step_vec = d * opts.step;
        let (p1u, p1v) = project_into_plane(da.du, da.dv, step_vec);
        let (p2u, p2v) = project_into_plane(db.du, db.dv, step_vec);

        let guess = (u1 + p1u, v1 + p1v, u2 + p2u, v2 + p2v);
        let Some(next) = correct(a, b, guess, opts.tol) else {
            break;
        };

        // Domain exit (open curve end).
        let (au0, au1) = a.u_domain();
        let (bu0, bu1) = b.u_domain();
        let eps = 1e-9;
        let on_edge = next.u1 <= au0 + eps
            || next.u1 >= au1 - eps
            || next.u2 <= bu0 + eps
            || next.u2 >= bu1 - eps;

        pts.push(next.point);
        prev_dir = Some(d);
        u1 = next.u1;
        v1 = next.v1;
        u2 = next.u2;
        v2 = next.v2;

        // Closure: returned near the start after a few steps.
        if step_i > 2 && (next.point - start.point).length() < opts.step * 0.9 {
            return Some(IntersectionPolyline {
                points: pts,
                closed: true,
            });
        }
        if on_edge {
            return Some(IntersectionPolyline {
                points: pts,
                closed: false,
            });
        }
    }
    Some(IntersectionPolyline {
        points: pts,
        closed: false,
    })
}

// ── linear algebra helpers ────────────────────────────────────────────────────

/// Least-squares parameter step that best realises 3-D vector `w` in the tangent
/// plane spanned by `su, sv`: solve `[su sv]·(a,b) ≈ w` (2×2 normal equations).
fn project_into_plane(su: Vec3, sv: Vec3, w: Vec3) -> (f64, f64) {
    let a = su.dot(su);
    let b = su.dot(sv);
    let c = sv.dot(sv);
    let d = su.dot(w);
    let e = sv.dot(w);
    let det = a * c - b * b;
    if det.abs() < 1e-18 {
        return (0.0, 0.0);
    }
    ((d * c - b * e) / det, (a * e - b * d) / det)
}

/// Solve a 4×4 linear system by Gaussian elimination with partial pivoting.
fn solve4(mut m: [[f64; 4]; 4], mut b: [f64; 4]) -> Option<[f64; 4]> {
    for col in 0..4 {
        // pivot
        let mut piv = col;
        for r in (col + 1)..4 {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if m[piv][col].abs() < 1e-14 {
            return None;
        }
        m.swap(col, piv);
        b.swap(col, piv);
        // eliminate
        for r in (col + 1)..4 {
            let f = m[r][col] / m[col][col];
            for c in col..4 {
                m[r][c] -= f * m[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = [0.0f64; 4];
    for i in (0..4).rev() {
        let mut s = b[i];
        for j in (i + 1)..4 {
            s -= m[i][j] * x[j];
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

#[inline]
fn wrap_v(v: f64) -> f64 {
    v.rem_euclid(TAU)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::swept::{CenterlineSeg, SweptTubeSurface};
    use crate::{cyl_cyl_intersection, cyl_surface_distance, CylSurf};
    use cadcore_math::{Point3, UnitVec3, Vec3};

    fn straight_tube(p0: Point3, p1: Point3, r: f64) -> SweptTubeSurface {
        let d = p1 - p0;
        SweptTubeSurface::new(
            vec![CenterlineSeg::Line {
                p0,
                dir: UnitVec3::try_from_vec(d).unwrap(),
                length: d.length(),
            }],
            r,
        )
    }

    /// Generic SSI on two straight tubes must reproduce the analytic
    /// cylinder∩cylinder loop (the woodpile crossing).
    #[test]
    fn ssi_two_straight_tubes_matches_analytic() {
        let r = 0.275;
        let a = straight_tube(Point3::new(-2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0), r);
        let b = straight_tube(Point3::new(0.0, -2.0, 0.35), Point3::new(0.0, 2.0, 0.35), r);

        let opts = SsiOptions {
            step: 0.03,
            ..Default::default()
        };
        let curves = surface_surface_intersection(&a, &b, &opts);
        assert!(!curves.is_empty(), "SSI found no curve");
        // Expect one closed loop near the crossing.
        let loop_ = curves.iter().max_by_key(|c| c.points.len()).unwrap();
        assert!(loop_.closed, "expected a closed loop");

        // All points lie on both analytic cylinders.
        let ca = CylSurf::new(Point3::new(-2.0, 0.0, 0.0), UnitVec3::X, r);
        let cb = CylSurf::new(Point3::new(0.0, -2.0, 0.35), UnitVec3::Y, r);
        for p in &loop_.points {
            assert!(cyl_surface_distance(&ca, *p).abs() < 1e-4, "on A");
            assert!(cyl_surface_distance(&cb, *p).abs() < 1e-4, "on B");
        }

        // Same crossing the analytic way → compare bounding boxes (same curve).
        let analytic = cyl_cyl_intersection(&ca, &cb, 96);
        assert_eq!(analytic.len(), 1);
        let bbox = |pts: &[Point3]| {
            let mut lo = Point3::new(f64::MAX, f64::MAX, f64::MAX);
            let mut hi = Point3::new(f64::MIN, f64::MIN, f64::MIN);
            for p in pts {
                lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
                hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
            }
            (lo, hi)
        };
        let (la, ha) = bbox(&loop_.points);
        let (lb, hb) = bbox(&analytic[0].points);
        for (x, y) in [
            (la.x, lb.x),
            (la.y, lb.y),
            (la.z, lb.z),
            (ha.x, hb.x),
            (ha.y, hb.y),
            (ha.z, hb.z),
        ] {
            assert!((x - y).abs() < 0.02, "bbox mismatch {x} vs {y}");
        }
    }

    /// Generic SSI on a straight tube ∩ a curved (arc) tube — the perimeter
    /// U-turn case that the analytic cyl-only path could not handle.
    #[test]
    fn ssi_straight_tube_meets_arc_tube() {
        let r = 0.275;
        // Arc elbow in the z=0 plane.
        let elbow = SweptTubeSurface::new(
            vec![CenterlineSeg::Arc {
                center: Point3::ORIGIN,
                axis: UnitVec3::Z,
                xref: UnitVec3::X,
                radius: 0.5,
                ang0: 0.0,
                ang1: std::f64::consts::FRAC_PI_2,
            }],
            r,
        );
        // A straight leg just above the elbow, crossing it near θ≈45°.
        let mid = elbow.point(elbow.length() * 0.5, 0.0); // a point on the elbow tube
        let leg = straight_tube(
            Point3::new(mid.x - 1.0, mid.y, 0.35),
            Point3::new(mid.x + 1.0, mid.y, 0.35),
            r,
        );

        let opts = SsiOptions {
            step: 0.02,
            ..Default::default()
        };
        let curves = surface_surface_intersection(&elbow, &leg, &opts);
        assert!(!curves.is_empty(), "no intersection found for straight∩arc");
        for c in &curves {
            let _ = c;
        }
        // Same crossing but a LONG leg (20 mm) crossing near one end — the
        // 20×20-scaffold case where coverage was being lost.
        let long_leg = straight_tube(
            Point3::new(mid.x - 0.5, mid.y, 0.35),
            Point3::new(mid.x + 19.5, mid.y, 0.35),
            r,
        );
        let curves_long = surface_surface_intersection(&elbow, &long_leg, &opts);
        assert!(
            !curves_long.is_empty(),
            "LONG-leg straight∩arc found no intersection (coverage bug)"
        );
        // Same long leg but with the union's actual step (= radius*0.25).
        let opts_big = SsiOptions {
            step: r * 0.25,
            ..Default::default()
        };
        let curves_big = surface_surface_intersection(&elbow, &long_leg, &opts_big);
        assert!(
            !curves_big.is_empty(),
            "LONG-leg with step={:.3} found nothing (THIS is the union coverage bug)",
            r * 0.25
        );

        // Elbow with an OFFSET, NEGATIVE-direction arc (ang0 > ang1) — exactly
        // how `swept_from_torus` parametrises real scaffold elbows (theta_lo is
        // an arbitrary atan2 value and d can be negative).
        let ang0 = 2.5_f64;
        let ang1 = 2.5 - std::f64::consts::FRAC_PI_2; // negative sweep
        let elbow2 = SweptTubeSurface::new(
            vec![CenterlineSeg::Arc {
                center: Point3::new(3.0, 4.0, 0.0),
                axis: UnitVec3::Z,
                xref: UnitVec3::X,
                radius: 0.5,
                ang0,
                ang1,
            }],
            r,
        );
        let mid2 = elbow2.point(elbow2.length() * 0.5, 0.0);
        let (_, t2) = elbow2.centerline(elbow2.length() * 0.5);
        // Leg perpendicular-ish to the elbow tangent, just above it, long.
        let perp = Vec3::new(-t2.as_vec().y, t2.as_vec().x, 0.0).normalize();
        let leg2 = straight_tube(
            Point3::new(mid2.x - perp.x * 0.5, mid2.y - perp.y * 0.5, mid2.z + 0.35),
            Point3::new(
                mid2.x + perp.x * 10.0,
                mid2.y + perp.y * 10.0,
                mid2.z + 0.35,
            ),
            r,
        );
        let curves2 = surface_surface_intersection(&elbow2, &leg2, &opts_big);
        assert!(
            !curves2.is_empty(),
            "OFFSET/NEGATIVE-arc elbow ∩ leg found nothing (THE union coverage bug)"
        );
        for c in &curves {
            for p in &c.points {
                // On the leg cylinder.
                let lc = CylSurf::new(Point3::new(mid.x - 1.0, mid.y, 0.35), UnitVec3::X, r);
                assert!(cyl_surface_distance(&lc, *p).abs() < 1e-3, "on leg");
                // On the elbow tube (distance to nearest centerline pt ≈ r).
                let (u, _) = elbow.project_point(*p).unwrap();
                let (cc, _) = elbow.centerline(u);
                assert!(((*p - cc).length() - r).abs() < 5e-3, "on elbow tube");
            }
        }
    }

    /// The EXACT 20×20-scaffold elbow×leg pair the union missed (dumped from
    /// Pass 2b): large coordinates (~19.6), crossing near both ends.
    #[test]
    fn ssi_real_20x20_pair() {
        let r = 0.275;
        let elbow = SweptTubeSurface::new(
            vec![CenterlineSeg::Arc {
                center: Point3::new(0.9, 19.225, 0.35),
                axis: UnitVec3::try_from_vec(Vec3::new(0.0, 0.0, -1.0)).unwrap(),
                xref: UnitVec3::try_from_vec(Vec3::new(0.0, -1.0, 0.0)).unwrap(),
                radius: 0.5,
                ang0: 1.5708,
                ang1: 3.1416,
            }],
            r,
        );
        let leg = straight_tube(
            Point3::new(19.725, 19.6, 0.7),
            Point3::new(0.775, 19.6, 0.7),
            r,
        );
        // Robust marching-squares SSI (grid surface a = the elbow).
        let curves = swept_swept_intersection(&elbow, &leg, 96, 48);
        eprintln!(
            "[test] real pair (marching) → {} curves: {:?}",
            curves.len(),
            curves
                .iter()
                .map(|c| (c.closed, c.points.len()))
                .collect::<Vec<_>>()
        );
        assert!(
            !curves.is_empty(),
            "real 20×20 pair: marching SSI found nothing"
        );
        // The key Task-2 fix: NO runaway (the predictor-corrector traced 4001
        // points here).  Marching squares yields small bounded loops.
        for c in &curves {
            assert!(c.points.len() < 200, "runaway loop: {} pts", c.points.len());
        }
        // Every accepted point has been corrected onto both parametric surfaces.
        let big = curves.iter().max_by_key(|c| c.points.len()).unwrap();
        let lc = CylSurf::new(Point3::new(19.725, 19.6, 0.7), UnitVec3::X, r);
        let max_leg_err = big
            .points
            .iter()
            .map(|p| cyl_surface_distance(&lc, *p).abs())
            .fold(0.0f64, f64::max);
        let max_elbow_err = big
            .points
            .iter()
            .map(|p| elbow.signed_distance(*p).unwrap().abs())
            .fold(0.0f64, f64::max);
        assert!(
            max_leg_err < 1e-5,
            "leg surface error too large: {max_leg_err}"
        );
        assert!(
            max_elbow_err < 1e-5,
            "elbow surface error too large: {max_elbow_err}"
        );
    }
}
