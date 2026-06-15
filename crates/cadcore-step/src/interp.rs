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
