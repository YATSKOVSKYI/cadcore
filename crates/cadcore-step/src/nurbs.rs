//! Exact **rational B-spline (NURBS)** representation of a *trimmed torus elbow*
//! patch, for STEP export.
//!
//! # Why
//!
//! A `TOROIDAL_SURFACE` is periodic in the major angle θ.  A face on it bounded
//! only by its two minor end-circles is **ambiguous**: the same two circles
//! bound both the short (≤180°) and the long (≥180°) θ-band, so importers
//! (SpaceClaim / ANSYS) routinely pick the wrong band and the U-turn elbow
//! "turns the wrong way" / fails to stitch to the adjoining straight tube.
//!
//! Representing the elbow as a NURBS surface that **only spans** the actual arc
//! `θ ∈ [θ_lo, θ_hi]` removes the ambiguity: the surface simply does not exist
//! outside the band.  Using the *rational* (exact) torus NURBS — rather than a
//! polynomial approximation — keeps the boundary minor circles lying exactly on
//! the surface, so the watertight edge-sharing is preserved.
//!
//! # Construction
//!
//! The patch is `make_revolved_surface` (Piegl & Tiller, *The NURBS Book*,
//! Algorithm A8.1) specialised to a torus: revolve the full minor circle
//! (a 9-control-point degree-2 rational full circle, the *v* / φ direction)
//! about the torus axis through the arc `[θ_lo, θ_hi]` (the *u* / θ direction,
//! 1–4 degree-2 rational arc segments of ≤90° each).

use cadcore_geom::TorusSurf;
use cadcore_math::Point3;

const C45: f64 = std::f64::consts::FRAC_1_SQRT_2; // cos(45°) = 1/√2

/// A tensor-product rational B-spline surface, ready for STEP emission.
///
/// `ctrl[i][j]` / `weights[i][j]` are indexed `[u_index][v_index]` — `u` is the
/// revolution (major-angle θ) direction, `v` is the profile (minor-angle φ)
/// direction.  Knot vectors are stored as distinct values + multiplicities, as
/// STEP's `B_SPLINE_SURFACE_WITH_KNOTS` expects.
#[derive(Clone, Debug)]
pub(crate) struct RationalBSplineSurface {
    pub u_degree: usize,
    pub v_degree: usize,
    pub u_knots: Vec<f64>,
    pub u_mult: Vec<usize>,
    pub v_knots: Vec<f64>,
    pub v_mult: Vec<usize>,
    pub u_closed: bool,
    pub v_closed: bool,
    pub ctrl: Vec<Vec<Point3>>,
    pub weights: Vec<Vec<f64>>,
}

impl RationalBSplineSurface {
    /// Number of control points in the u (revolution) direction.
    pub fn n_u(&self) -> usize {
        self.ctrl.len()
    }
    /// Number of control points in the v (profile) direction.
    pub fn n_v(&self) -> usize {
        self.ctrl.first().map_or(0, |row| row.len())
    }

    /// `(u,v)` parameter domain bounds `((u_lo,u_hi),(v_lo,v_hi))`.
    pub fn domain(&self) -> ((f64, f64), (f64, f64)) {
        (
            (
                *self.u_knots.first().unwrap_or(&0.0),
                *self.u_knots.last().unwrap_or(&1.0),
            ),
            (
                *self.v_knots.first().unwrap_or(&0.0),
                *self.v_knots.last().unwrap_or(&4.0),
            ),
        )
    }

    /// Invert the surface: find `(u,v)` with `S(u,v) ≈ p` (nearest point).
    /// Coarse grid seed + Gauss-Newton refinement using finite-difference
    /// surface derivatives.  Used to build pcurves for elbow (B-spline) faces.
    pub fn project_point(&self, p: Point3) -> (f64, f64) {
        let ((ulo, uhi), (vlo, vhi)) = self.domain();
        // Coarse seed.
        let (mut bu, mut bv, mut bd) = (ulo, vlo, f64::MAX);
        let (gu, gv) = (8usize, 24usize);
        for i in 0..=gu {
            for j in 0..=gv {
                let u = ulo + (uhi - ulo) * i as f64 / gu as f64;
                let v = vlo + (vhi - vlo) * j as f64 / gv as f64;
                let d = (self.eval(u, v) - p).length_sq();
                if d < bd {
                    bd = d;
                    bu = u;
                    bv = v;
                }
            }
        }
        // Gauss-Newton refine.  Derivative sampling stays strictly inside the
        // domain: `eval` outside the knot range has zero basis sum (0/0 = NaN),
        // and trim-curve endpoints land exactly on the patch boundary.
        let (mut u, mut v) = (bu, bv);
        let eps = 1e-6;
        for _ in 0..40 {
            let r = self.eval(u, v) - p; // residual Vec3
            let uc = u.clamp(ulo + eps, uhi - eps);
            let vc = v.clamp(vlo + eps, vhi - eps);
            let su = (self.eval(uc + eps, v) - self.eval(uc - eps, v)) * (0.5 / eps);
            let sv = (self.eval(u, vc + eps) - self.eval(u, vc - eps)) * (0.5 / eps);
            let (a, b, c) = (su.dot(su), su.dot(sv), sv.dot(sv));
            let (g1, g2) = (su.dot(r), sv.dot(r));
            let det = a * c - b * b;
            if det.abs() < 1e-14 {
                break;
            }
            let du = -(c * g1 - b * g2) / det;
            let dv = -(-b * g1 + a * g2) / det;
            if !(du.is_finite() && dv.is_finite()) {
                break;
            }
            u = (u + du).clamp(ulo, uhi);
            v = (v + dv).clamp(vlo, vhi);
            if du.abs() + dv.abs() < 1e-12 {
                break;
            }
        }
        if u.is_finite() && v.is_finite() {
            (u, v)
        } else {
            (bu, bv)
        }
    }

    /// Flatten a (distinct-knot, multiplicity) pair into the full knot vector.
    fn flat_knots(knots: &[f64], mult: &[usize]) -> Vec<f64> {
        let mut v = Vec::new();
        for (&k, &m) in knots.iter().zip(mult.iter()) {
            for _ in 0..m {
                v.push(k);
            }
        }
        v
    }

    /// Evaluate the rational surface point at parameters `(u, v)`.
    ///
    /// Used by tests to verify the patch reproduces the analytic torus exactly.
    pub fn eval(&self, u: f64, v: f64) -> Point3 {
        let uk = Self::flat_knots(&self.u_knots, &self.u_mult);
        let vk = Self::flat_knots(&self.v_knots, &self.v_mult);
        let nu = self.n_u();
        let nv = self.n_v();

        let mut num = cadcore_math::Vec3::new(0.0, 0.0, 0.0);
        let mut den = 0.0;
        for i in 0..nu {
            let bu = bspline_basis(i, self.u_degree, u, &uk);
            if bu == 0.0 {
                continue;
            }
            for j in 0..nv {
                let bv = bspline_basis(j, self.v_degree, v, &vk);
                if bv == 0.0 {
                    continue;
                }
                let w = self.weights[i][j] * bu * bv;
                let p = self.ctrl[i][j];
                num = num + cadcore_math::Vec3::new(p.x, p.y, p.z) * w;
                den += w;
            }
        }
        Point3::new(num.x / den, num.y / den, num.z / den)
    }
}

/// Cox-de Boor basis function `N_{i,p}(u)` over the flat knot vector `knots`.
///
/// Includes the standard half-open-interval convention plus the right-endpoint
/// fix so `N` of the last function is 1 at `u == knots.last()`.
fn bspline_basis(i: usize, p: usize, u: f64, knots: &[f64]) -> f64 {
    let m = knots.len() - 1;
    // Right-endpoint special case.
    if i == m - p - 1 && (u - knots[m]).abs() < 1e-12 {
        return 1.0;
    }
    if p == 0 {
        return if u >= knots[i] && u < knots[i + 1] {
            1.0
        } else {
            0.0
        };
    }
    let mut left = 0.0;
    let d1 = knots[i + p] - knots[i];
    if d1 > 0.0 {
        left = (u - knots[i]) / d1 * bspline_basis(i, p - 1, u, knots);
    }
    let mut right = 0.0;
    let d2 = knots[i + p + 1] - knots[i + 1];
    if d2 > 0.0 {
        right = (knots[i + p + 1] - u) / d2 * bspline_basis(i + 1, p - 1, u, knots);
    }
    left + right
}

/// Build the exact rational-NURBS patch for the torus arc `θ ∈ [theta_lo,
/// theta_hi]` (signed; the arc follows the short signed sweep between them) over
/// the full minor circle `φ ∈ [0, 2π]`.
pub(crate) fn torus_patch_nurbs(
    t: &TorusSurf,
    theta_lo: f64,
    theta_hi: f64,
) -> RationalBSplineSurface {
    let r_major = t.major_radius;
    let r_minor = t.minor_radius;

    // ── v direction: full minor circle, 9 control points (a = radial distance
    //    from the torus axis, h = height along the axis) + weights. ───────────
    // On-circle points have weight 1; the 4 corner points have weight cos45 and
    // sit at the tangent intersections (distance r_minor/cos45 from the minor
    // centre, i.e. offset ±r_minor in both a and h).
    let profile: [(f64, f64, f64); 9] = [
        (r_major + r_minor, 0.0, 1.0),      // φ=0
        (r_major + r_minor, r_minor, C45),  // φ=45 (corner)
        (r_major, r_minor, 1.0),            // φ=90
        (r_major - r_minor, r_minor, C45),  // φ=135 (corner)
        (r_major - r_minor, 0.0, 1.0),      // φ=180
        (r_major - r_minor, -r_minor, C45), // φ=225 (corner)
        (r_major, -r_minor, 1.0),           // φ=270
        (r_major + r_minor, -r_minor, C45), // φ=315 (corner)
        (r_major + r_minor, 0.0, 1.0),      // φ=360 == φ=0
    ];

    // ── u direction: revolve through [theta_lo, theta_hi]. ────────────────────
    // A single degree-2 *rational* Bézier arc is exact for any sweep < 180°; we
    // cap each segment at 120° for good numerical conditioning (corner weight
    // cos60° = 0.5).  This keeps a genuine 90° elbow as ONE segment (float noise
    // pushing it a hair past 90° must not force a needless split) while a 180°
    // U-turn correctly becomes two segments.
    let total = theta_hi - theta_lo;
    let seg_max = 2.0 * std::f64::consts::PI / 3.0; // 120°
    let narcs = ((total.abs() / (seg_max - 1e-9)).ceil() as usize).clamp(1, 6);
    let dtheta = total / narcs as f64;
    let wm = (dtheta * 0.5).cos().abs(); // corner weight per segment

    let n_u = 2 * narcs + 1;
    let mut ctrl = vec![Vec::with_capacity(9); n_u];
    let mut weights = vec![Vec::with_capacity(9); n_u];

    let radial =
        |theta: f64| -> cadcore_math::Vec3 { t.frame.x * theta.cos() + t.frame.y * theta.sin() };
    let zaxis = t.frame.z;

    for iu in 0..n_u {
        // Even index = on-arc boundary (weight factor 1, radius factor 1).
        // Odd index  = corner between two boundaries (weight wm, radius 1/wm).
        let (theta_i, radius_factor, u_weight) = if iu % 2 == 0 {
            let seg = iu / 2;
            (theta_lo + dtheta * seg as f64, 1.0, 1.0)
        } else {
            let seg = iu / 2;
            (theta_lo + dtheta * (seg as f64 + 0.5), 1.0 / wm, wm)
        };
        let rad = radial(theta_i);
        for &(a, h, vw) in &profile {
            let p = t.frame.origin + rad * (a * radius_factor) + zaxis * h;
            ctrl[iu].push(p);
            weights[iu].push(u_weight * vw);
        }
    }

    // Knot vectors.  Degree 2 in both directions.
    // u: `narcs` Bézier segments → distinct knots 0..narcs, interior mult 2.
    let u_knots: Vec<f64> = (0..=narcs).map(|k| k as f64).collect();
    let mut u_mult = vec![2usize; narcs + 1];
    u_mult[0] = 3;
    u_mult[narcs] = 3;
    // v: full circle → 4 Bézier segments, distinct knots 0..4, interior mult 2.
    let v_knots: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0];
    let v_mult = vec![3usize, 2, 2, 2, 3];

    RationalBSplineSurface {
        u_degree: 2,
        v_degree: 2,
        u_knots,
        u_mult,
        v_knots,
        v_mult,
        u_closed: false,
        v_closed: true,
        ctrl,
        weights,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_math::{Point3, UnitVec3};

    fn approx(a: Point3, b: Point3, tol: f64) -> bool {
        (a - b).length() < tol
    }

    /// The patch must reproduce the analytic torus exactly at the interpolated
    /// boundary parameters (segment corners are on-surface) for a 90° elbow.
    #[test]
    fn torus_patch_90_matches_analytic() {
        let t = TorusSurf::new(Point3::new(1.0, 2.0, 3.0), UnitVec3::Z, 0.5, 0.275);
        let lo = 0.3;
        let hi = lo + std::f64::consts::FRAC_PI_2; // 90°
        let s = torus_patch_nurbs(&t, lo, hi);
        assert_eq!(s.n_u(), 3);
        assert_eq!(s.n_v(), 9);

        // u knots {0,1}; v knots {0,1,2,3,4}. On-knot params interpolate the
        // on-circle control points = exact analytic torus points.
        // u=0 → θ_lo, u=1 → θ_hi. v=k → φ = k·90°.
        let cases = [
            (0.0, 0.0, lo, 0.0),
            (0.0, 1.0, lo, std::f64::consts::FRAC_PI_2),
            (0.0, 2.0, lo, std::f64::consts::PI),
            (1.0, 0.0, hi, 0.0),
            (1.0, 3.0, hi, 3.0 * std::f64::consts::FRAC_PI_2),
        ];
        for (u, v, theta, phi) in cases {
            let got = s.eval(u, v);
            let want = t.point_at(theta, phi);
            assert!(
                approx(got, want, 1e-9),
                "u={u} v={v}: got {got:?} want {want:?}"
            );
        }
    }

    /// Interior parameters (not on knots) must ALSO lie exactly on the torus —
    /// this is what distinguishes the exact rational patch from an approximation.
    #[test]
    fn torus_patch_interior_on_surface() {
        let t = TorusSurf::new(Point3::ORIGIN, UnitVec3::Z, 0.5, 0.275);
        let lo = 0.0;
        let hi = std::f64::consts::FRAC_PI_2;
        let s = torus_patch_nurbs(&t, lo, hi);
        // Sample a grid of interior params; every point must be at distance
        // r_minor from the spine circle of radius r_major.
        for iu in 0..=8 {
            for iv in 0..=16 {
                let u = iu as f64 / 8.0;
                let v = iv as f64 * 4.0 / 16.0; // v domain is [0,4]
                let p = s.eval(u, v);
                // Distance from torus axis (z through origin):
                let a = (p.x * p.x + p.y * p.y).sqrt();
                let h = p.z;
                // Must satisfy (a - R)² + h² = r²  (point on the tube).
                let d = ((a - t.major_radius).powi(2) + h * h).sqrt();
                assert!(
                    (d - t.minor_radius).abs() < 1e-9,
                    "u={u} v={v}: tube dist {d} != r {}",
                    t.minor_radius
                );
            }
        }
    }

    /// `project_point` recovers the `(u,v)` of a known on-surface point.
    #[test]
    fn torus_patch_project_roundtrip() {
        let t = TorusSurf::new(Point3::new(1.0, -2.0, 0.5), UnitVec3::Z, 0.5, 0.275);
        let s = torus_patch_nurbs(&t, 0.2, 0.2 + std::f64::consts::FRAC_PI_2);
        for &(u0, v0) in &[(0.3_f64, 0.7_f64), (0.5, 2.3), (0.85, 3.6), (0.1, 1.5)] {
            let p = s.eval(u0, v0);
            let (u, v) = s.project_point(p);
            // The recovered params must map back to the same 3D point.
            assert!(
                (s.eval(u, v) - p).length() < 1e-7,
                "(u0,v0)=({u0},{v0}) → ({u},{v}); pt err {}",
                (s.eval(u, v) - p).length()
            );
        }
    }

    /// Points exactly ON the patch boundary (trim-curve endpoints land on the
    /// junction circle = u-domain edge) must project to FINITE (u,v) — the
    /// Gauss-Newton FD sampling used to step outside the knot range, where the
    /// rational eval is 0/0 = NaN, poisoning the emitted pcurve.
    #[test]
    fn torus_patch_project_boundary_is_finite() {
        let t = TorusSurf::new(Point3::new(1.0, -2.0, 0.5), UnitVec3::Z, 0.5, 0.275);
        let s = torus_patch_nurbs(&t, 0.2, 0.2 + std::f64::consts::FRAC_PI_2);
        let ((ulo, uhi), (vlo, vhi)) = s.domain();
        for &(u0, v0) in &[(ulo, vlo), (uhi, vlo), (ulo, vhi), (uhi, vhi), (uhi, 1.7), (ulo, 3.1)] {
            let p = s.eval(u0, v0);
            let (u, v) = s.project_point(p);
            assert!(u.is_finite() && v.is_finite(), "(u0,v0)=({u0},{v0}) → NaN");
            assert!(
                (s.eval(u, v) - p).length() < 1e-7,
                "(u0,v0)=({u0},{v0}) → ({u},{v}); pt err {}",
                (s.eval(u, v) - p).length()
            );
        }
    }

    /// A 180° U-turn must split into 2 segments (5 u control points) and stay
    /// exact.
    #[test]
    fn torus_patch_180_splits_and_matches() {
        let t = TorusSurf::new(Point3::ORIGIN, UnitVec3::Z, 0.6, 0.3);
        let lo = 0.0;
        let hi = std::f64::consts::PI; // 180°
        let s = torus_patch_nurbs(&t, lo, hi);
        assert_eq!(s.n_u(), 5, "180° → 2 segments → 5 u control points");
        // Midpoint of the arc (u at the internal knot = 1.0) is θ = 90°.
        let got = s.eval(1.0, 0.0);
        let want = t.point_at(std::f64::consts::FRAC_PI_2, 0.0);
        assert!(approx(got, want, 1e-9), "got {got:?} want {want:?}");
    }
}
