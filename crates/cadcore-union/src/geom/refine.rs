//! The projection oracle — refinement of trim geometry onto emitted surfaces.
//!
//! **The single most important lesson of this kernel** (see README invariant
//! #1): intersection curves computed by marching on sampled swept surfaces sit
//! 1–5 µm off the analytic surfaces that the writer emits, and SolidWorks
//! (knitting at ≈1 µm) responds by silently dropping the affected faces.
//!
//! Every point that the engine commits to the B-Rep must therefore take a
//! final trip through this module:
//!
//! * curve interior points → [`refine_point`] onto the 2 incident surfaces;
//! * joints between curve pieces → [`refine_point`] onto **all** incident
//!   surfaces (a triple point), and the result becomes the single source of
//!   truth for every face that references the joint;
//! * cut points on shared boundary circles → [`refine_on_circle`], so both
//!   faces partition the circle at *identical* parameters.
//!
//! Cyclic projection converges fast for the near-orthogonal surface sets that
//! occur at filament junctions; `max_iter` is a hard stop for the pathological
//! (near-tangent) configurations, which the caller must detect via
//! [`crate::geom::classify`] and handle explicitly.

use cadcore_geom::{Circle3, CylSurf, Plane3, TorusSurf};
use cadcore_math::{Point3, Vec3};

/// Anything the refinement and tracing machinery can project onto.
///
/// [`AnalyticSurface`] implements this for single surfaces; composite
/// geometry (a whole filament: leg cylinders + elbow tori) implements it by
/// delegating to the nearest member, which is exactly what trim curves that
/// span junctions need.
pub trait Projectable {
    /// Closest point on the surface (degenerate loci return the input).
    fn project_p(&self, p: Point3) -> Point3;
    /// Outward surface normal at (the projection of) `p`.
    fn normal_p(&self, p: Point3) -> Option<Vec3>;
    /// Unsigned distance from `p` to the surface.
    fn distance_p(&self, p: Point3) -> f64 {
        (self.project_p(p) - p).length()
    }
}

impl Projectable for AnalyticSurface {
    fn project_p(&self, p: Point3) -> Point3 {
        self.project(p)
    }
    fn normal_p(&self, p: Point3) -> Option<Vec3> {
        self.normal_at(p)
    }
}

/// An analytic surface the engine can project onto exactly.
///
/// This enum mirrors what the STEP writer can emit *exactly* (plane, cylinder,
/// and the rational-NURBS torus patch which is pointwise identical to the
/// analytic torus).  Anything not represented here must be refined against
/// its emitted approximation instead.
#[derive(Debug, Clone, Copy)]
pub enum AnalyticSurface {
    /// An infinite plane.
    Plane(Plane3),
    /// An infinite cylinder.
    Cylinder(CylSurf),
    /// A full torus.
    Torus(TorusSurf),
}

impl AnalyticSurface {
    /// Closest point on the surface.
    ///
    /// Degenerate inputs (a point on a cylinder axis, on the torus axis, or
    /// at the torus centre-circle) have no unique projection; the point is
    /// returned unchanged — callers at those loci are already in a
    /// configuration the classifier must have flagged.
    pub fn project(&self, p: Point3) -> Point3 {
        match self {
            AnalyticSurface::Plane(pl) => p - pl.normal().as_vec() * pl.signed_distance(p),
            AnalyticSurface::Cylinder(s) => {
                let w = p - s.frame.origin;
                let ax = s.axis().dot_vec(w);
                let rad = w - s.axis().as_vec() * ax;
                let rl = rad.length();
                if rl < 1e-12 {
                    return p;
                }
                s.frame.origin + s.axis().as_vec() * ax + rad * (s.radius / rl)
            }
            AnalyticSurface::Torus(t) => {
                let w = p - t.frame.origin;
                let h = t.frame.z.dot_vec(w);
                let q = w - t.frame.z.as_vec() * h;
                let a = q.length();
                if a < 1e-12 {
                    return p;
                }
                let ring = t.frame.origin + q * (t.major_radius / a);
                let v = p - ring;
                let vl = v.length();
                if vl < 1e-12 {
                    return p;
                }
                ring + v * (t.minor_radius / vl)
            }
        }
    }

    /// Unsigned distance from `p` to the surface.
    pub fn distance(&self, p: Point3) -> f64 {
        (self.project(p) - p).length()
    }

    /// Surface normal at (the projection of) `p`, pointing away from the
    /// material-defining axis (outward for cylinder/torus tubes, the plane's
    /// stored normal for planes).  Returns `None` at degenerate loci.
    pub fn normal_at(&self, p: Point3) -> Option<cadcore_math::Vec3> {
        match self {
            AnalyticSurface::Plane(pl) => Some(pl.normal().as_vec()),
            AnalyticSurface::Cylinder(s) => {
                let w = p - s.frame.origin;
                let ax = s.axis().dot_vec(w);
                let rad = w - s.axis().as_vec() * ax;
                let rl = rad.length();
                if rl < 1e-12 {
                    return None;
                }
                Some(rad * (1.0 / rl))
            }
            AnalyticSurface::Torus(t) => {
                let w = p - t.frame.origin;
                let h = t.frame.z.dot_vec(w);
                let q = w - t.frame.z.as_vec() * h;
                let a = q.length();
                if a < 1e-12 {
                    return None;
                }
                let ring = t.frame.origin + q * (t.major_radius / a);
                let v = p - ring;
                let vl = v.length();
                if vl < 1e-12 {
                    return None;
                }
                Some(v * (1.0 / vl))
            }
        }
    }
}

/// Refine `p` onto the common intersection of `surfaces` by cyclic projection.
///
/// Returns the refined point and the residual displacement of the final
/// sweep (0 ⇒ converged to `tol`).  With 2 surfaces this converges to the
/// intersection *curve*; with 3 to the triple *point*.
pub fn refine_point(
    p: Point3,
    surfaces: &[AnalyticSurface],
    tol: f64,
    max_iter: usize,
) -> (Point3, f64) {
    let dyns: Vec<&dyn Projectable> = surfaces.iter().map(|s| s as &dyn Projectable).collect();
    refine_point_dyn(p, &dyns, tol, max_iter)
}

/// [`refine_point`] over arbitrary [`Projectable`]s (composite filaments…).
pub fn refine_point_dyn(
    mut p: Point3,
    surfaces: &[&dyn Projectable],
    tol: f64,
    max_iter: usize,
) -> (Point3, f64) {
    let mut residual = f64::INFINITY;
    for _ in 0..max_iter {
        let before = p;
        for s in surfaces {
            p = s.project_p(p);
        }
        residual = (p - before).length();
        if residual < tol {
            break;
        }
    }
    (p, residual)
}

/// Refine a cut point that must live ON a shared boundary circle while also
/// lying on the surfaces of `others`.
///
/// The final word always belongs to the circle: both faces that share the
/// circle take their partition parameters from the returned point, which is
/// exactly on it (README invariant #2).
pub fn refine_on_circle(
    p: Point3,
    circle: &Circle3,
    others: &[AnalyticSurface],
    tol: f64,
    max_iter: usize,
) -> Point3 {
    let on_circle = |q: Point3| -> Point3 {
        let w = q - circle.frame.origin;
        let x = circle.frame.x.dot_vec(w);
        let y = circle.frame.y.dot_vec(w);
        let ang = y.atan2(x);
        circle.point_at(ang)
    };
    let mut q = p;
    for _ in 0..max_iter {
        let before = q;
        q = on_circle(q);
        for s in others {
            q = s.project(q);
        }
        if (q - before).length() < tol {
            break;
        }
    }
    on_circle(q)
}

/// Refine every interior point of a polyline curve onto its two incident
/// surfaces.  Endpoints are left alone — they are joints and must be refined
/// by the caller with the full incident-surface set, then written back.
pub fn refine_curve_interior(
    pts: &mut [Point3],
    surfaces: &[AnalyticSurface],
    tol: f64,
    max_iter: usize,
) {
    let n = pts.len();
    if n <= 2 {
        return;
    }
    for p in &mut pts[1..n - 1] {
        let (q, _) = refine_point(*p, surfaces, tol, max_iter);
        *p = q;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_math::UnitVec3;

    fn cyl_x(y: f64, z: f64, r: f64) -> CylSurf {
        CylSurf::new(Point3::new(0.0, y, z), UnitVec3::X, r)
    }

    fn cyl_y(x: f64, z: f64, r: f64) -> CylSurf {
        CylSurf::new(Point3::new(x, 0.0, z), UnitVec3::Y, r)
    }

    #[test]
    fn project_onto_cylinder() {
        let s = AnalyticSurface::Cylinder(cyl_x(0.0, 0.0, 0.275));
        let q = s.project(Point3::new(1.0, 0.3, 0.4));
        assert!((s.distance(q)) < 1e-12);
        // axial coordinate preserved
        assert!((q.x - 1.0).abs() < 1e-12);
    }

    #[test]
    fn project_onto_torus() {
        let t = TorusSurf::new(Point3::new(1.0, 2.0, 3.0), UnitVec3::Z, 0.5, 0.275);
        let s = AnalyticSurface::Torus(t);
        let q = s.project(Point3::new(1.7, 2.1, 3.4));
        assert!(s.distance(q) < 1e-12);
    }

    #[test]
    fn two_perpendicular_cylinders_converges_to_curve() {
        // classic woodpile: equal radii, perpendicular axes, offset 0.35
        let a = AnalyticSurface::Cylinder(cyl_x(0.0, 0.0, 0.275));
        let b = AnalyticSurface::Cylinder(cyl_y(0.0, 0.35, 0.275));
        // a noisy point near the intersection curve (µm-level, like SSI output)
        let p = Point3::new(0.11, 0.13, 0.21);
        let (q, res) = refine_point(p, &[a, b], 1e-12, 60);
        assert!(res < 1e-12, "did not converge: {res:.3e}");
        assert!(a.distance(q) < 1e-10);
        assert!(b.distance(q) < 1e-10);
    }

    #[test]
    fn triple_point_cyl_cyl_plane() {
        let a = AnalyticSurface::Cylinder(cyl_x(0.0, 0.0, 0.275));
        let b = AnalyticSurface::Cylinder(cyl_y(0.0, 0.35, 0.275));
        let c = AnalyticSurface::Plane(Plane3::from_origin_normal(
            Point3::new(0.1, 0.0, 0.0),
            UnitVec3::X,
        ));
        let (q, res) = refine_point(Point3::new(0.12, 0.1, 0.2), &[a, b, c], 1e-12, 80);
        assert!(res < 1e-11, "residual {res:.3e}");
        for s in [&a, &b, &c] {
            assert!(s.distance(q) < 1e-9, "off surface by {:.2e}", s.distance(q));
        }
    }

    #[test]
    fn circle_keeps_final_word() {
        let circle = Circle3::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::X, 0.275);
        let other = AnalyticSurface::Cylinder(cyl_y(0.0, 0.35, 0.275));
        let q = refine_on_circle(
            Point3::new(0.01, 0.2, 0.21),
            &circle,
            &[other],
            1e-12,
            60,
        );
        // exactly on the circle
        let w = q - circle.frame.origin;
        assert!(w.dot(circle.frame.z.as_vec()).abs() < 1e-12);
        assert!((w.length() - 0.275).abs() < 1e-12);
    }

    #[test]
    fn curve_interior_refined_endpoints_untouched() {
        let a = AnalyticSurface::Cylinder(cyl_x(0.0, 0.0, 0.275));
        let b = AnalyticSurface::Cylinder(cyl_y(0.0, 0.35, 0.275));
        let start = Point3::new(0.2, 0.1, 0.15);
        let end = Point3::new(-0.2, 0.1, 0.15);
        let mut pts = vec![start, Point3::new(0.1, 0.14, 0.2), end];
        refine_curve_interior(&mut pts, &[a, b], 1e-12, 60);
        assert_eq!(pts[0], start);
        assert_eq!(pts[2], end);
        assert!(a.distance(pts[1]) < 1e-10);
        assert!(b.distance(pts[1]) < 1e-10);
    }
}
