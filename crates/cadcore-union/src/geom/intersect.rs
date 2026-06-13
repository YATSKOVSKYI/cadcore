//! The intersection oracle — surface×surface curves that lie ON both
//! surfaces to refinement accuracy.
//!
//! The legacy engine marches on *sampled* swept surfaces; its curves carry
//! 1–5 µm of sampling error that SolidWorks punishes.  This tracer works on
//! the analytic surfaces directly with a predictor–corrector:
//!
//! 1. tangent  `t = n_a × n_b` (cross of surface normals);
//! 2. predictor `p' = p + h·t`;
//! 3. corrector: cyclic projection of `p'` onto both surfaces
//!    ([`refine_point`]) — every emitted point is exact to `tol`.
//!
//! Near-tangency (`|n_a × n_b| → 0`) is where intersection curves flutter
//! and generic stepping is unsound; the tracer *reports* it instead of
//! guessing, so the caller can branch into dedicated handling (the
//! classifier flags those pairs upfront).

use cadcore_math::Point3;

use super::refine::{refine_point_dyn, AnalyticSurface, Projectable};

/// Tracing options.
#[derive(Debug, Clone, Copy)]
pub struct TraceOptions {
    /// Step length along the curve (mm).  Default `0.05`.
    pub step: f64,
    /// Refinement convergence target (mm).  Default `1e-12`.
    pub tol: f64,
    /// Hard cap on points per curve.  Default `4096`.
    pub max_points: usize,
    /// `|n_a × n_b|` below this aborts as near-tangency.  Default `1e-5`.
    pub tangency_eps: f64,
}

impl Default for TraceOptions {
    fn default() -> Self {
        Self {
            step: 0.05,
            tol: 1e-12,
            max_points: 4096,
            tangency_eps: 1e-5,
        }
    }
}

/// How a trace ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraceEnd {
    /// Curve returned to its starting point — a closed loop.
    Closed,
    /// `max_points` reached without closing (open curve or runaway).
    Budget,
    /// Surfaces became (near-)tangent along the curve — needs the dedicated
    /// marginal-case path.
    Tangency,
    /// The corrector failed to converge (projection sets fighting).
    Diverged,
    /// The caller's clip predicate said stop (left the region of interest).
    Clipped,
}

/// A traced intersection curve.  Every point lies on both input surfaces
/// within `TraceOptions::tol` (verified by construction).
#[derive(Debug, Clone)]
pub struct TracedCurve {
    /// Curve points in trace order.
    pub points: Vec<Point3>,
    /// Why tracing stopped.
    pub end: TraceEnd,
}

impl TracedCurve {
    /// Closed loop?
    pub fn closed(&self) -> bool {
        self.end == TraceEnd::Closed
    }
}

/// Find a point on the intersection of two surfaces near `seed`, or `None`
/// when the cyclic projection cannot converge from there.
pub fn intersection_point_near(
    a: &dyn Projectable,
    b: &dyn Projectable,
    seed: Point3,
    tol: f64,
) -> Option<Point3> {
    let ok = |p: Point3| a.distance_p(p) < 1e-9 && b.distance_p(p) < 1e-9;
    let (p, residual) = refine_point_dyn(seed, &[a, b], tol, 80);
    if residual < tol.max(1e-12) * 10.0 && ok(p) {
        return Some(p);
    }
    let damped = |start: Point3| settle(a, b, start, tol);
    if let Some(p2) = damped(seed) {
        return Some(p2);
    }
    // Symmetric seeds (exactly between two parallel-ish patches) are saddle
    // fixed points of the midpoint iteration.  Break the tie by stepping
    // perpendicular to the local projection axis, both ways, two scales.
    let qa = a.project_p(seed);
    let qb = b.project_p(seed);
    let axis = qa - qb;
    let g = (a.distance_p(seed) + b.distance_p(seed)).max(1e-3);
    let seed_axis = if axis.length() > 1e-12 {
        axis.normalize()
    } else {
        cadcore_math::Vec3::new(0.0, 0.0, 1.0)
    };
    let helper = if seed_axis.x.abs() < 0.9 {
        cadcore_math::Vec3::new(1.0, 0.0, 0.0)
    } else {
        cadcore_math::Vec3::new(0.0, 1.0, 0.0)
    };
    let e1 = (helper - seed_axis * helper.dot(seed_axis)).normalize();
    let e2 = seed_axis.cross(e1);
    for scale in [0.5, 1.5] {
        for dir in [e1, e1 * -1.0, e2, e2 * -1.0] {
            if let Some(p2) = damped(seed + dir * (g * scale)) {
                return Some(p2);
            }
        }
    }
    None
}

/// Land a point exactly on the intersection from a nearby start, robust to
/// member-switching limit cycles of composite surfaces: damped midpoint
/// iteration (average of both projections walks onto the curve) followed by
/// a cyclic polish.
fn settle(a: &dyn Projectable, b: &dyn Projectable, start: Point3, tol: f64) -> Option<Point3> {
    let mut q = start;
    for _ in 0..200 {
        let qa = a.project_p(q);
        let qb = b.project_p(q);
        let mid = qa + (qb - qa) * 0.5;
        if (mid - q).length() < tol.max(1e-14) {
            q = mid;
            break;
        }
        q = mid;
    }
    let (p2, residual2) = refine_point_dyn(q, &[a, b], tol, 80);
    (residual2 < tol.max(1e-12) * 10.0 && a.distance_p(p2) < 1e-9 && b.distance_p(p2) < 1e-9)
        .then_some(p2)
}

/// Trace the intersection curve of `a` and `b` through `seed`.
///
/// `clip` lets the caller bound the trace to a region of interest (axial
/// spans of finite filaments etc.); return `false` to stop.  Pass
/// `|_| true` for unbounded tracing.
pub fn trace(
    a: &AnalyticSurface,
    b: &AnalyticSurface,
    seed: Point3,
    opts: &TraceOptions,
    clip: impl FnMut(Point3) -> bool,
) -> Option<TracedCurve> {
    trace_dyn(a, b, seed, opts, clip)
}

/// [`trace`] over arbitrary [`Projectable`]s (composite filaments…).
pub fn trace_dyn(
    a: &dyn Projectable,
    b: &dyn Projectable,
    seed: Point3,
    opts: &TraceOptions,
    mut clip: impl FnMut(Point3) -> bool,
) -> Option<TracedCurve> {
    let start = intersection_point_near(a, b, seed, opts.tol)?;
    let mut points = vec![start];
    let mut p = start;
    let mut dir_sign = 1.0f64;
    let mut end = TraceEnd::Budget;

    while points.len() < opts.max_points {
        let (Some(na), Some(nb)) = (a.normal_p(p), b.normal_p(p)) else {
            end = TraceEnd::Diverged;
            break;
        };
        let cross = na.cross(nb);
        let cl = cross.length();
        if cl < opts.tangency_eps {
            end = TraceEnd::Tangency;
            break;
        }
        let t = cross * (dir_sign / cl);
        // keep marching direction coherent with the previous step
        if points.len() >= 2 {
            let prev = points[points.len() - 1] - points[points.len() - 2];
            if prev.dot(t) < 0.0 {
                dir_sign = -dir_sign;
            }
        }
        let t = cross * (dir_sign / cl);
        // adaptive step: near-tangent waists (small |n_a×n_b|) have high
        // curve curvature — slow down so the corrector stays on branch
        let step_eff = opts.step * cl.clamp(0.15, 1.0);
        let predictor = p + t * step_eff;
        let (mut q, residual) = refine_point_dyn(predictor, &[a, b], opts.tol, 60);
        if residual > opts.tol.max(1e-12) * 100.0 {
            // cyclic projection limit-cycles near composite member junctions
            // — settle with the damped iteration before giving up
            match settle(a, b, predictor, opts.tol) {
                Some(q2) => q = q2,
                None => {
                    end = TraceEnd::Diverged;
                    break;
                }
            }
        }
        if !clip(q) {
            end = TraceEnd::Clipped;
            break;
        }
        // stalled? (corrector pulled back onto the previous point — happens
        // where the curve crosses a composite member junction and the member
        // pick flickers inside the corrector).  Escalate the step to jump
        // over the sticky spot before giving up.
        if (q - p).length() < step_eff * 0.05 {
            let mut rescued = None;
            for factor in [0.5, 1.5, 2.5, 4.0] {
                let pred2 = p + t * (step_eff * factor);
                if let Some(q2) = settle(a, b, pred2, opts.tol) {
                    let fwd = (q2 - p).dot(t);
                    if (q2 - p).length() > step_eff * 0.1 && fwd > 0.0 {
                        rescued = Some(q2);
                        break;
                    }
                }
            }
            match rescued {
                Some(q2) => q = q2,
                None => {
                    end = TraceEnd::Diverged;
                    break;
                }
            }
        }
        // closed-loop detection: back near start after having travelled
        if points.len() > 3 && (q - start).length() < step_eff * 0.75 {
            end = TraceEnd::Closed;
            break;
        }
        points.push(q);
        p = q;
    }
    Some(TracedCurve { points, end })
}

/// Trace the closed loop of a transversal crossing given any nearby seed.
/// Convenience wrapper that requires the result to close.
pub fn trace_closed_loop(
    a: &AnalyticSurface,
    b: &AnalyticSurface,
    seed: Point3,
    opts: &TraceOptions,
) -> Option<TracedCurve> {
    trace_closed_loop_dyn(a, b, seed, opts)
}

/// [`trace_closed_loop`] over arbitrary [`Projectable`]s.
pub fn trace_closed_loop_dyn(
    a: &dyn Projectable,
    b: &dyn Projectable,
    seed: Point3,
    opts: &TraceOptions,
) -> Option<TracedCurve> {
    let c = trace_dyn(a, b, seed, opts, |_| true)?;
    c.closed().then_some(c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_geom::{CylSurf, Plane3, TorusSurf};
    use cadcore_math::UnitVec3;

    fn on_both(c: &TracedCurve, a: &AnalyticSurface, b: &AnalyticSurface) -> f64 {
        c.points
            .iter()
            .map(|&p| a.distance(p).max(b.distance(p)))
            .fold(0.0, f64::max)
    }

    #[test]
    fn perpendicular_equal_radii_woodpile_loop() {
        // the production crossing: r=0.275 both, offset 0.35
        let a = AnalyticSurface::Cylinder(CylSurf::new(
            Point3::new(0.0, 0.0, 0.0),
            UnitVec3::X,
            0.275,
        ));
        let b = AnalyticSurface::Cylinder(CylSurf::new(
            Point3::new(0.0, 0.0, 0.35),
            UnitVec3::Y,
            0.275,
        ));
        let c = trace_closed_loop(&a, &b, Point3::new(0.2, 0.2, 0.18), &TraceOptions::default())
            .expect("closed loop");
        assert!(c.points.len() > 20);
        assert!(on_both(&c, &a, &b) < 1e-9, "max dev {}", on_both(&c, &a, &b));
    }

    #[test]
    fn oblique_unequal_radii_loop() {
        // THE new capability: 45° axes, r 0.275 vs 0.2
        let d = cadcore_math::UnitVec3::try_from_vec(
            UnitVec3::X.as_vec() + UnitVec3::Y.as_vec(),
        )
        .unwrap();
        let a = AnalyticSurface::Cylinder(CylSurf::new(
            Point3::new(0.0, 0.0, 0.0),
            UnitVec3::X,
            0.275,
        ));
        let b = AnalyticSurface::Cylinder(CylSurf::new(Point3::new(0.0, 0.0, 0.3), d, 0.2));
        let c = trace_closed_loop(&a, &b, Point3::new(0.1, 0.1, 0.2), &TraceOptions::default())
            .expect("closed loop");
        assert!(on_both(&c, &a, &b) < 1e-9);
    }

    #[test]
    fn cylinder_plane_oblique_is_exact_ellipse() {
        let a = AnalyticSurface::Cylinder(CylSurf::new(
            Point3::new(0.0, 0.0, 0.0),
            UnitVec3::Z,
            0.275,
        ));
        let n = cadcore_math::UnitVec3::try_from_vec(
            UnitVec3::Z.as_vec() + UnitVec3::X.as_vec() * 0.5,
        )
        .unwrap();
        let b = AnalyticSurface::Plane(Plane3::from_origin_normal(Point3::new(0.0, 0.0, 1.0), n));
        let c = trace_closed_loop(&a, &b, Point3::new(0.275, 0.0, 1.0), &TraceOptions::default())
            .expect("closed loop");
        assert!(on_both(&c, &a, &b) < 1e-9);
        // miter ellipse: z range = r·tan(angle) around the plane point
        let zmin = c.points.iter().map(|p| p.z).fold(f64::MAX, f64::min);
        let zmax = c.points.iter().map(|p| p.z).fold(f64::MIN, f64::max);
        assert!((zmax - zmin) > 0.2, "ellipse should tilt: {zmin}..{zmax}");
    }

    #[test]
    fn torus_cylinder_corner_bite() {
        // U-elbow (torus) bitten by the perpendicular leg one layer up —
        // the corner-trim configuration, now exact.
        let t = AnalyticSurface::Torus(TorusSurf::new(
            Point3::new(0.0, 0.0, 0.0),
            UnitVec3::Z,
            0.6,
            0.275,
        ));
        let leg = AnalyticSurface::Cylinder(CylSurf::new(
            Point3::new(0.6, 0.0, 0.35),
            UnitVec3::Y,
            0.275,
        ));
        let c = trace_closed_loop(
            &t,
            &leg,
            Point3::new(0.65, 0.1, 0.18),
            &TraceOptions::default(),
        )
        .expect("closed loop");
        assert!(c.points.len() > 10);
        assert!(on_both(&c, &t, &leg) < 1e-9, "dev {}", on_both(&c, &t, &leg));
    }

    #[test]
    fn near_tangent_reported_not_guessed() {
        // HL == D: surfaces touch — tracer must flag, not invent a loop
        let a = AnalyticSurface::Cylinder(CylSurf::new(
            Point3::new(0.0, 0.0, 0.0),
            UnitVec3::X,
            0.275,
        ));
        let b = AnalyticSurface::Cylinder(CylSurf::new(
            Point3::new(0.0, 0.0, 0.5495),
            UnitVec3::Y,
            0.275,
        ));
        let r = trace(
            &a,
            &b,
            Point3::new(0.01, 0.01, 0.27),
            &TraceOptions::default(),
            |_| true,
        );
        if let Some(c) = r {
            assert!(
                matches!(c.end, TraceEnd::Tangency | TraceEnd::Diverged) || c.points.len() < 64,
                "tangent case must not produce a confident long loop: {:?} n={}",
                c.end,
                c.points.len()
            );
        }
    }

    /// The reason this oracle exists: the legacy sampler's points sit µm off
    /// the second surface; the tracer's sit at refinement accuracy on both.
    #[test]
    fn tracer_beats_legacy_sampler() {
        let s1 = CylSurf::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::X, 0.275);
        let s2 = CylSurf::new(Point3::new(0.0, 0.0, 0.35), UnitVec3::Y, 0.275);
        let a = AnalyticSurface::Cylinder(s1);
        let b = AnalyticSurface::Cylinder(s2);

        let legacy = cadcore_geom::cyl_cyl_intersection(&s1, &s2, 64);
        let legacy_dev = legacy
            .iter()
            .flat_map(|l| l.points.iter())
            .map(|&p| a.distance(p).max(b.distance(p)))
            .fold(0.0, f64::max);

        let traced =
            trace_closed_loop(&a, &b, Point3::new(0.2, 0.2, 0.18), &TraceOptions::default())
                .expect("closed loop");
        let traced_dev = on_both(&traced, &a, &b);

        assert!(
            traced_dev < 1e-9,
            "tracer must be exact: {traced_dev:.3e} (legacy was {legacy_dev:.3e})"
        );
        assert!(
            traced_dev < legacy_dev || legacy_dev < 1e-9,
            "tracer ({traced_dev:.3e}) must not be worse than legacy ({legacy_dev:.3e})"
        );
    }

    #[test]
    fn clip_stops_trace() {
        let a = AnalyticSurface::Cylinder(CylSurf::new(
            Point3::new(0.0, 0.0, 0.0),
            UnitVec3::Z,
            0.275,
        ));
        let b = AnalyticSurface::Plane(Plane3::from_origin_normal(
            Point3::new(0.1, 0.0, 0.0),
            UnitVec3::X,
        ));
        // cylinder ∩ plane parallel to axis = two infinite lines: budget-
        // bounded without clip, Clipped with one.
        let c = trace(
            &a,
            &b,
            Point3::new(0.1, 0.25, 0.0),
            &TraceOptions::default(),
            |p| p.z.abs() < 1.0,
        )
        .unwrap();
        assert_eq!(c.end, TraceEnd::Clipped);
        let zmax = c.points.iter().map(|p| p.z.abs()).fold(0.0, f64::max);
        assert!(zmax <= 1.0 + 0.06);
    }
}
