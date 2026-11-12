//! Global **cubic B-spline interpolation** for tolerance-sampled intersection
//! curves (Piegl & Tiller, *The NURBS Book*, Algorithm A9.1).
//!
//! A surface–surface intersection is sampled as a polyline of points lying
//! exactly on both carriers.  Emitting it as a degree-1 chord curve leaves a
//! deviation of the chord sagitta (`≈ h²/8R`).  A degree-3 curve that
//! *interpolates* the same on-surface points hugs the true smooth curve to
//! `≈ h⁴` — orders of magnitude tighter at the same sample spacing — which is
//! how analytic kernels (OCC, Parasolid) carry intersection edges.
//!
//! The 3-D curve and every pcurve share ONE parameterisation (`ū`) and knot
//! vector (`U`), so a STEP `SURFACE_CURVE` pairs them consistently.

use cadcore_math::Point3;

use crate::entities::{Ctx, StepError};
use crate::nurbs::bspline_basis;

/// Chord-length parameters `ū` and the clamped degree-3 knot vector `U` (flat,
/// with repeats) for interpolating `pts`.  Needs `pts.len() >= 4`.
pub(crate) fn cubic_params_knots(pts: &[Point3]) -> (Vec<f64>, Vec<f64>) {
    let np = pts.len();
    let n = np - 1;
    // chord-length parameters in [0,1]
    let mut seg = vec![0.0f64; np];
    let mut total = 0.0;
    for i in 1..np {
        seg[i] = (pts[i] - pts[i - 1]).length();
        total += seg[i];
    }
    let mut t = vec![0.0f64; np];
    if total < 1e-30 {
        for (i, ti) in t.iter_mut().enumerate() {
            *ti = i as f64 / n as f64;
        }
    } else {
        let mut acc = 0.0;
        for i in 1..np {
            acc += seg[i];
            t[i] = acc / total;
        }
        t[n] = 1.0;
    }
    // clamped knots: U0..U3 = 0, U_{n+1}..U_{n+4} = 1, interior = averaged params
    let m = n + 4; // last knot index
    let mut u = vec![0.0f64; m + 1];
    for k in (n + 1)..=(n + 4) {
        u[k] = 1.0;
    }
    for j in 1..=(n.saturating_sub(3)) {
        u[j + 3] = (t[j] + t[j + 1] + t[j + 2]) / 3.0;
    }
    (t, u)
}

/// `(n+1)×(n+1)` basis matrix `N[i][k] = N_{k,3}(ū_i)`.
fn basis_matrix(t: &[f64], u: &[f64]) -> Vec<Vec<f64>> {
    let np = t.len();
    let mut a = vec![vec![0.0f64; np]; np];
    for i in 0..np {
        for k in 0..np {
            a[i][k] = bspline_basis(k, 3, t[i], u);
        }
    }
    a
}

/// Solve `A·x = b` by Gaussian elimination with partial pivoting (A is `m×m`).
fn solve(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let m = b.len();
    let mut a: Vec<Vec<f64>> = a.iter().map(|r| r.clone()).collect();
    let mut b = b.to_vec();
    for c in 0..m {
        // pivot
        let mut piv = c;
        for r in (c + 1)..m {
            if a[r][c].abs() > a[piv][c].abs() {
                piv = r;
            }
        }
        a.swap(c, piv);
        b.swap(c, piv);
        let d = a[c][c];
        if d.abs() < 1e-14 {
            continue;
        }
        for r in (c + 1)..m {
            let f = a[r][c] / d;
            if f == 0.0 {
                continue;
            }
            for k in c..m {
                a[r][k] -= f * a[c][k];
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = vec![0.0f64; m];
    for r in (0..m).rev() {
        let mut s = b[r];
        for k in (r + 1)..m {
            s -= a[r][k] * x[k];
        }
        x[r] = if a[r][r].abs() < 1e-14 { 0.0 } else { s / a[r][r] };
    }
    x
}

/// Solve the interpolation system for one coordinate (values sampled per point).
fn ctrl_coord(a: &[Vec<f64>], vals: &[f64]) -> Vec<f64> {
    solve(a, vals)
}

/// Flatten knots → (distinct knots, multiplicities) for STEP.
fn distinct_knots(u: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut kv: Vec<f64> = Vec::new();
    let mut mult: Vec<usize> = Vec::new();
    for &k in u {
        if let Some(last) = kv.last() {
            if (k - *last).abs() < 1e-12 {
                *mult.last_mut().unwrap() += 1;
                continue;
            }
        }
        kv.push(k);
        mult.push(1);
    }
    (kv, mult)
}

/// Emit a degree-3 3-D `B_SPLINE_CURVE_WITH_KNOTS` interpolating `pts`, returning
/// `(entity_id, ū, U)` so a matching pcurve can reuse the parameterisation.
/// Returns `None` for fewer than 4 points (caller falls back to a polyline).
pub(crate) fn emit_cubic_curve3d(
    ctx: &mut Ctx,
    pts: &[Point3],
) -> Result<Option<(usize, Vec<f64>, Vec<f64>)>, StepError> {
    if pts.len() < 4 {
        return Ok(None);
    }
    let (t, u) = cubic_params_knots(pts);
    let a = basis_matrix(&t, &u);
    let cx = ctrl_coord(&a, &pts.iter().map(|p| p.x).collect::<Vec<_>>());
    let cy = ctrl_coord(&a, &pts.iter().map(|p| p.y).collect::<Vec<_>>());
    let cz = ctrl_coord(&a, &pts.iter().map(|p| p.z).collect::<Vec<_>>());
    let cp_ids: Vec<usize> = (0..pts.len())
        .map(|i| crate::entities::emit_point(ctx, Point3::new(cx[i], cy[i], cz[i]), "bsp"))
        .collect::<Result<_, _>>()?;
    let id = emit_bspline_with_knots(ctx, &cp_ids, &u)?;
    Ok(Some((id, t, u)))
}

/// Emit a degree-3 2-D pcurve interpolating `uv` at the SAME parameters `t` and
/// knots `u` as its 3-D curve (so the SURFACE_CURVE stays consistent).
pub(crate) fn emit_cubic_curve2d(
    ctx: &mut Ctx,
    uv: &[(f64, f64)],
    t: &[f64],
    u: &[f64],
) -> Result<usize, StepError> {
    let a = basis_matrix(t, u);
    let cu = ctrl_coord(&a, &uv.iter().map(|p| p.0).collect::<Vec<_>>());
    let cv = ctrl_coord(&a, &uv.iter().map(|p| p.1).collect::<Vec<_>>());
    let cp_ids: Vec<usize> = (0..uv.len())
        .map(|i| crate::entities::emit_point2d(ctx, cu[i], cv[i]))
        .collect::<Result<_, _>>()?;
    emit_bspline_with_knots(ctx, &cp_ids, u)
}

/// Emit the `B_SPLINE_CURVE_WITH_KNOTS` entity (degree 3) from control-point ids
/// and a flat knot vector.
fn emit_bspline_with_knots(ctx: &mut Ctx, cp_ids: &[usize], u: &[f64]) -> Result<usize, StepError> {
    let (kv, mult) = distinct_knots(u);
    let refs: String = cp_ids
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let mults: String = mult
        .iter()
        .map(|m| m.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let knots: String = kv
        .iter()
        .map(|k| format!("{k:.6}"))
        .collect::<Vec<_>>()
        .join(",");
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        &format!(
            "B_SPLINE_CURVE_WITH_KNOTS('',3,({refs}),.UNSPECIFIED.,.F.,.F.,({mults}),({knots}),.UNSPECIFIED.)"
        ),
    )?;
    Ok(id)
}

// ── Least-squares cubic approximation sharing one knot vector ────────────────
//
// Interpolating every sample (above) keeps as many control points as samples:
// a sampled intersection edge of a fused lattice then carries ~40 control
// points in 3-D and again in every pcurve, and CAD importers (SpaceClaim) spend
// minutes re-fitting them.  A smooth intersection curve needs far fewer: fit
// the samples in the least-squares sense (Piegl & Tiller A9.7, ends clamped
// onto the end samples) and add control points until EVERY sample of the 3-D
// curve and of every pcurve is within tolerance.  All curves of one edge use
// the same parameters and knots, so a SURFACE_CURVE stays consistent.

/// The fitted edge: flat clamped knot vector, 3-D control points and one set of
/// 2-D control points per pcurve (same order as the input `uvs`).
pub(crate) struct SharedCubicFit {
    pub u: Vec<f64>,
    pub ctrl3: Vec<Point3>,
    pub ctrl2: Vec<Vec<(f64, f64)>>,
}

/// Chord-length parameters in `[0,1]`.
fn chord_params(pts: &[Point3]) -> Vec<f64> {
    let np = pts.len();
    let mut t = vec![0.0f64; np];
    let mut total = 0.0;
    for i in 1..np {
        total += (pts[i] - pts[i - 1]).length();
        t[i] = total;
    }
    if total < 1e-30 {
        for (i, ti) in t.iter_mut().enumerate() {
            *ti = i as f64 / (np - 1) as f64;
        }
    } else {
        for ti in t.iter_mut() {
            *ti /= total;
        }
        t[np - 1] = 1.0;
    }
    t
}

/// Knot span index of `x` in the clamped cubic knot vector with `nc` control points.
fn find_span(nc: usize, x: f64, u: &[f64]) -> usize {
    let n = nc - 1;
    if x >= u[n + 1] {
        return n;
    }
    let (mut lo, mut hi) = (3usize, n + 1);
    let mut mid = (lo + hi) / 2;
    while x < u[mid] || x >= u[mid + 1] {
        if x < u[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
    }
    mid
}

/// The four non-zero cubic basis values at `x` in span `s` (Piegl & Tiller A2.2).
fn basis4(s: usize, x: f64, u: &[f64]) -> [f64; 4] {
    let mut n = [0.0f64; 4];
    let mut left = [0.0f64; 4];
    let mut right = [0.0f64; 4];
    n[0] = 1.0;
    for j in 1..=3 {
        left[j] = x - u[s + 1 - j];
        right[j] = u[s + j] - x;
        let mut saved = 0.0;
        for r in 0..j {
            let den = right[r + 1] + left[j - r];
            let tmp = if den.abs() < 1e-300 { 0.0 } else { n[r] / den };
            n[r] = saved + right[r + 1] * tmp;
            saved = left[j - r] * tmp;
        }
        n[j] = saved;
    }
    n
}

/// Least-squares fit of `vals` (`dim` coordinates per sample) with `nc` control
/// points on knots `u`, ends clamped onto the first/last sample.  Returns the
/// control points and the max sample error (Euclidean).
fn lsq_fit(vals: &[Vec<f64>], nc: usize, rows: &[(usize, [f64; 4])]) -> (Vec<Vec<f64>>, f64) {
    let m = vals.len() - 1;
    let n = nc - 1;
    let dim = vals[0].len();
    let unk = n - 1; // interior control points P1..P(n-1)
    let mut ctrl = vec![vec![0.0f64; dim]; nc];
    ctrl[0] = vals[0].clone();
    ctrl[n] = vals[m].clone();
    if unk > 0 {
        // normal equations (NᵀN) P = Nᵀ R
        let mut a = vec![vec![0.0f64; unk]; unk];
        let mut b = vec![vec![0.0f64; unk]; dim];
        for k in 1..m {
            let (s, nb) = rows[k];
            // R_k = Q_k - N0 Q0 - Nn Qm
            let mut r = vals[k].clone();
            for (j, &w) in nb.iter().enumerate() {
                let i = s - 3 + j;
                if i == 0 {
                    for d in 0..dim {
                        r[d] -= w * vals[0][d];
                    }
                } else if i == n {
                    for d in 0..dim {
                        r[d] -= w * vals[m][d];
                    }
                }
            }
            for (j1, &w1) in nb.iter().enumerate() {
                let i1 = s - 3 + j1;
                if i1 == 0 || i1 == n {
                    continue;
                }
                for d in 0..dim {
                    b[d][i1 - 1] += w1 * r[d];
                }
                for (j2, &w2) in nb.iter().enumerate() {
                    let i2 = s - 3 + j2;
                    if i2 == 0 || i2 == n {
                        continue;
                    }
                    a[i1 - 1][i2 - 1] += w1 * w2;
                }
            }
        }
        for d in 0..dim {
            let x = solve(&a, &b[d]);
            for i in 0..unk {
                ctrl[i + 1][d] = x[i];
            }
        }
    }
    let mut err = 0.0f64;
    for k in 0..=m {
        let (s, nb) = rows[k];
        let mut e2 = 0.0;
        for d in 0..dim {
            let mut c = 0.0;
            for (j, &w) in nb.iter().enumerate() {
                c += w * ctrl[s - 3 + j][d];
            }
            e2 += (c - vals[k][d]).powi(2);
        }
        err = err.max(e2.sqrt());
    }
    (ctrl, err)
}

/// Fit `pts` and every pcurve sample list in `uvs` (same length and order as
/// `pts`) with ONE cubic knot vector, adding control points until each sample
/// is within `tol` (mm in 3-D, parameter units on the pcurves).  `None` when
/// no fit with fewer control points than samples meets the tolerance — the
/// caller then interpolates every sample as before.
pub(crate) fn approx_cubic_shared(pts: &[Point3], uvs: &[Vec<(f64, f64)>], tol: f64) -> Option<SharedCubicFit> {
    let np = pts.len();
    if np < 8 || uvs.iter().any(|w| w.len() != np) {
        return None;
    }
    let t = chord_params(pts);
    let m = np - 1;
    let v3: Vec<Vec<f64>> = pts.iter().map(|p| vec![p.x, p.y, p.z]).collect();
    let v2: Vec<Vec<Vec<f64>>> = uvs.iter().map(|w| w.iter().map(|q| vec![q.0, q.1]).collect()).collect();
    let mut nc = 4usize.max(np / 8);
    while nc < np {
        // knots (Piegl & Tiller eq. 9.69): every span holds samples
        let n = nc - 1;
        let mut u = vec![0.0f64; nc + 4];
        for k in (n + 1)..(nc + 4) {
            u[k] = 1.0;
        }
        let d = (m + 1) as f64 / (n - 2) as f64;
        for j in 1..=(n - 3) {
            let jd = j as f64 * d;
            let i = jd.floor() as usize;
            let a = jd - i as f64;
            u[3 + j] = (1.0 - a) * t[i - 1] + a * t[i];
        }
        let rows: Vec<(usize, [f64; 4])> = t
            .iter()
            .map(|&x| {
                let s = find_span(nc, x, &u);
                (s, basis4(s, x, &u))
            })
            .collect();
        let (c3, e3) = lsq_fit(&v3, nc, &rows);
        let mut ok = e3 <= tol;
        let mut c2s: Vec<Vec<(f64, f64)>> = Vec::with_capacity(uvs.len());
        if ok {
            for vals in &v2 {
                let (c2, e2) = lsq_fit(vals, nc, &rows);
                if e2 > tol {
                    ok = false;
                    break;
                }
                c2s.push(c2.into_iter().map(|c| (c[0], c[1])).collect());
            }
        }
        if ok {
            return Some(SharedCubicFit {
                u,
                ctrl3: c3.into_iter().map(|c| Point3::new(c[0], c[1], c[2])).collect(),
                ctrl2: c2s,
            });
        }
        let next = (nc * 3 / 2).max(nc + 2);
        if next >= np {
            return None;
        }
        nc = next;
    }
    None
}

/// Emit a fitted 3-D cubic from [`approx_cubic_shared`].
pub(crate) fn emit_fit_curve3d(ctx: &mut Ctx, fit: &SharedCubicFit) -> Result<usize, StepError> {
    let ids: Vec<usize> = fit
        .ctrl3
        .iter()
        .map(|&p| crate::entities::emit_point(ctx, p, "bsp"))
        .collect::<Result<_, _>>()?;
    emit_bspline_with_knots(ctx, &ids, &fit.u)
}

/// Emit the `k`-th fitted pcurve from [`approx_cubic_shared`].
pub(crate) fn emit_fit_curve2d(ctx: &mut Ctx, fit: &SharedCubicFit, k: usize) -> Result<usize, StepError> {
    let ids: Vec<usize> = fit.ctrl2[k]
        .iter()
        .map(|&(a, b)| crate::entities::emit_point2d(ctx, a, b))
        .collect::<Result<_, _>>()?;
    emit_bspline_with_knots(ctx, &ids, &fit.u)
}

#[cfg(test)]
mod fit_tests {
    use super::*;

    fn eval(fit_u: &[f64], ctrl: &[Point3], x: f64) -> Point3 {
        let nc = ctrl.len();
        let s = find_span(nc, x, fit_u);
        let nb = basis4(s, x, fit_u);
        let mut p = (0.0, 0.0, 0.0);
        for (j, &w) in nb.iter().enumerate() {
            let c = ctrl[s - 3 + j];
            p.0 += w * c.x;
            p.1 += w * c.y;
            p.2 += w * c.z;
        }
        Point3::new(p.0, p.1, p.2)
    }

    /// A densely sampled helix (smooth, like an intersection curve) is carried
    /// by far fewer control points, every sample and every mid-sample point
    /// stays on the true curve, and the ends are exact.
    #[test]
    fn helix_is_reduced_and_stays_on_curve() {
        let r = 0.275;
        let curve = |s: f64| Point3::new(r * (3.0 * s).cos(), r * (3.0 * s).sin(), 0.4 * s);
        let n = 60;
        let pts: Vec<Point3> = (0..n).map(|i| curve(i as f64 / (n - 1) as f64)).collect();
        let uv: Vec<(f64, f64)> = (0..n).map(|i| (3.0 * i as f64 / (n - 1) as f64, 0.4 * i as f64 / (n - 1) as f64)).collect();
        let fit = approx_cubic_shared(&pts, &[uv], 1e-6).expect("fit");
        assert!(fit.ctrl3.len() < n / 2, "control points {}", fit.ctrl3.len());
        assert!((fit.ctrl3[0] - pts[0]).length() < 1e-15);
        assert!((*fit.ctrl3.last().unwrap() - pts[n - 1]).length() < 1e-15);
        // true-curve deviation everywhere (sampled finely along the fit)
        let t = chord_params(&pts);
        for k in 0..n - 1 {
            for f in [0.0, 0.25, 0.5, 0.75] {
                let x = t[k] + f * (t[k + 1] - t[k]);
                let p = eval(&fit.u, &fit.ctrl3, x);
                let s = (k as f64 + f) / (n - 1) as f64;
                assert!((p - curve(s)).length() < 2e-6, "dev {} at {k}+{f}", (p - curve(s)).length());
            }
        }
        assert_eq!(fit.ctrl2.len(), 1);
        assert_eq!(fit.ctrl2[0].len(), fit.ctrl3.len());
    }

    /// A curve with a kink cannot be reduced at a tight tolerance → `None`
    /// (the caller keeps interpolating every sample).
    #[test]
    fn short_runs_are_left_alone() {
        let pts: Vec<Point3> = (0..5).map(|i| Point3::new(i as f64, 0.0, 0.0)).collect();
        assert!(approx_cubic_shared(&pts, &[], 1e-6).is_none());
    }
}
