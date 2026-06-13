//! Phase 3 of the swept-tube union: **classification** (is a point inside a
//! filament solid?) and **p-curve** helpers (project an intersection loop into a
//! surface's `(u,v)` domain, measure its winding, orient it as a hole).
//!
//! These generalise the cylinder-specific helpers in `cadcore-ops/union.rs` to
//! any [`SweptTubeSurface`] / [`ParamSurface`], so face splitting and Boolean
//! classification work for straight *and* curved filament regions uniformly.

use cadcore_math::Point3;

use crate::swept::{ParamSurface, SweptTubeSurface};

const TAU: f64 = std::f64::consts::PI * 2.0;

/// Is `p` inside the solid swept tube?  A point is inside iff it lies within the
/// tube radius of the (closest point on the) centerline:
/// `dist(p, C(t_min)) < r − tol`.
///
/// Note: the centerline-distance test models the tube ends as hemispheres, which
/// is exact for interior crossings (the union case that matters); the very tips
/// are off by at most the cap fillet — irrelevant for trimming classification.
pub fn point_inside_tube(tube: &SweptTubeSurface, p: Point3, tol: f64) -> bool {
    let Some((u, _)) = tube.project_point(p) else {
        return false;
    };
    let (c, _) = tube.centerline(u);
    (p - c).length() < tube.radius() - tol
}

/// Project a 3-D loop onto a surface's `(u, v)` parameter domain (drops points
/// that fail to project).
pub fn project_loop_to_uv(surf: &dyn ParamSurface, pts: &[Point3]) -> Vec<(f64, f64)> {
    pts.iter().filter_map(|p| surf.project_point(*p)).collect()
}

/// Signed area of a `(u, v)` loop, with `v` (period `v_period`) unwrapped so a
/// localized loop is continuous.  Positive = CCW, negative = CW.
pub fn loop_signed_area_uv(uv: &[(f64, f64)], v_period: f64) -> f64 {
    if uv.len() < 3 {
        return 0.0;
    }
    let mut w = uv.to_vec();
    for k in 1..w.len() {
        while w[k].1 - w[k - 1].1 > v_period * 0.5 {
            w[k].1 -= v_period;
        }
        while w[k].1 - w[k - 1].1 < -v_period * 0.5 {
            w[k].1 += v_period;
        }
    }
    let mut a = 0.0;
    for k in 0..w.len() {
        let (x0, y0) = w[k];
        let (x1, y1) = w[(k + 1) % w.len()];
        a += x0 * y1 - x1 * y0;
    }
    a * 0.5
}

/// Reorder a 3-D loop so its projection onto `surf`'s `(u, v)` is wound **CW**
/// (a proper hole that removes the enclosed region — opposite the outer bound).
pub fn orient_loop_as_hole(surf: &dyn ParamSurface, pts: &[Point3]) -> Vec<Point3> {
    let uv = project_loop_to_uv(surf, pts);
    if loop_signed_area_uv(&uv, TAU) > 0.0 {
        pts.iter().rev().copied().collect()
    } else {
        pts.to_vec()
    }
}

/// A point on `surf`'s surface at the centroid of a loop's `(u, v)` projection —
/// i.e. a sample *inside* the loop region.  Used to classify which side of the
/// loop is interior to the other solid.
pub fn loop_interior_sample(surf: &dyn ParamSurface, pts: &[Point3]) -> Option<Point3> {
    let uv = project_loop_to_uv(surf, pts);
    if uv.len() < 3 {
        return None;
    }
    // Unwrap v, average, wrap back.
    let mut w = uv.clone();
    for k in 1..w.len() {
        while w[k].1 - w[k - 1].1 > TAU * 0.5 {
            w[k].1 -= TAU;
        }
        while w[k].1 - w[k - 1].1 < -TAU * 0.5 {
            w[k].1 += TAU;
        }
    }
    let n = w.len() as f64;
    let uc = w.iter().map(|p| p.0).sum::<f64>() / n;
    let vc = w.iter().map(|p| p.1).sum::<f64>() / n;
    Some(surf.point(uc, vc.rem_euclid(TAU)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssi::{surface_surface_intersection, SsiOptions};
    use crate::swept::{CenterlineSeg, SweptTubeSurface};
    use cadcore_math::UnitVec3;

    fn straight(p0: Point3, p1: Point3, r: f64) -> SweptTubeSurface {
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

    #[test]
    fn inside_outside_straight_tube() {
        let r = 0.275;
        let t = straight(Point3::new(-2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0), r);
        assert!(point_inside_tube(&t, Point3::new(0.0, 0.0, 0.0), 1e-9)); // on axis
        assert!(point_inside_tube(&t, Point3::new(1.0, 0.1, 0.0), 1e-9)); // radial 0.1 < r
        assert!(!point_inside_tube(&t, Point3::new(0.0, 0.5, 0.0), 1e-9)); // radial 0.5 > r
        assert!(!point_inside_tube(&t, Point3::new(10.0, 0.0, 0.0), 1e-9)); // beyond end
    }

    #[test]
    fn inside_outside_arc_tube() {
        let r = 0.275;
        let arc = SweptTubeSurface::new(
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
        // Point on the spine at θ=45° is inside.
        let a = std::f64::consts::FRAC_PI_4;
        let spine = Point3::new(0.5 * a.cos(), 0.5 * a.sin(), 0.0);
        assert!(point_inside_tube(&arc, spine, 1e-9));
        // Far above is outside.
        assert!(!point_inside_tube(
            &arc,
            Point3::new(spine.x, spine.y, 1.0),
            1e-9
        ));
    }

    #[test]
    fn signed_area_sign_and_orientation() {
        // CCW unit square in (u,v).
        let ccw = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        assert!(loop_signed_area_uv(&ccw, TAU) > 0.0);
        let cw: Vec<_> = ccw.iter().rev().copied().collect();
        assert!(loop_signed_area_uv(&cw, TAU) < 0.0);
    }

    /// Integration: at a crossing, the bite loop on tube A encloses a region of
    /// A that is *inside* tube B (so the union removes it), and
    /// `orient_loop_as_hole` makes that loop CW on A.
    #[test]
    fn bite_loop_interior_is_inside_other() {
        let r = 0.275;
        let a = straight(Point3::new(-2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0), r);
        let b = straight(Point3::new(0.0, -2.0, 0.35), Point3::new(0.0, 2.0, 0.35), r);
        let opts = SsiOptions {
            step: 0.03,
            ..Default::default()
        };
        let curves = surface_surface_intersection(&a, &b, &opts);
        let loop_ = curves.iter().max_by_key(|c| c.points.len()).unwrap();
        assert!(loop_.closed);

        // Interior sample of the loop on A must be inside B.
        let sample = loop_interior_sample(&a, &loop_.points).unwrap();
        assert!(
            point_inside_tube(&b, sample, 1e-6),
            "loop interior on A should be inside B"
        );

        // Oriented as a hole → CW (negative area) in A's (u,v).
        let oriented = orient_loop_as_hole(&a, &loop_.points);
        let uv = project_loop_to_uv(&a, &oriented);
        assert!(loop_signed_area_uv(&uv, TAU) <= 0.0, "hole must be CW on A");
    }
}
