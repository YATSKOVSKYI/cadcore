//! Analytic curve types.

use cadcore_math::{Frame3, Point3, UnitVec3, Vec3, TAU, EPS};

// ── Line3 ────────────────────────────────────────────────────────────────────

/// An infinite directed line in 3-D space.
///
/// Parameterised as `P(t) = origin + t · direction`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Line3 {
    /// A point on the line.
    pub origin: Point3,
    /// Direction (unit vector).
    pub direction: UnitVec3,
}

impl Line3 {
    /// Construct from origin and direction.
    pub fn new(origin: Point3, direction: UnitVec3) -> Self { Self { origin, direction } }

    /// Point at parameter `t`.
    #[inline]
    pub fn point_at(self, t: f64) -> Point3 { self.origin + self.direction * t }

    /// Closest parameter value to world point `p`.
    #[inline]
    pub fn project(self, p: Point3) -> f64 {
        self.direction.dot_vec(p - self.origin)
    }

    /// Closest point on the line to `p`.
    #[inline]
    pub fn closest_point(self, p: Point3) -> Point3 {
        self.point_at(self.project(p))
    }

    /// Squared distance from `p` to the line.
    #[inline]
    pub fn dist_sq(self, p: Point3) -> f64 {
        (p - self.closest_point(p)).length_sq()
    }
}

// ── Circle3 ──────────────────────────────────────────────────────────────────

/// A planar circle.
///
/// Parameterised as `P(θ) = frame.origin + frame.x·cos(θ)·r + frame.y·sin(θ)·r`.
/// `frame.z` is the outward normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Circle3 {
    /// Placement — origin is the centre; Z is the outward normal; X is the
    /// 0-angle reference direction.
    pub frame: Frame3,
    /// Radius (mm).
    pub radius: f64,
}

impl Circle3 {
    /// Construct from centre, normal, and radius.
    ///
    /// The X axis of the frame is chosen arbitrarily.
    pub fn new(centre: Point3, normal: UnitVec3, radius: f64) -> Self {
        Self { frame: Frame3::from_origin_z(centre, normal), radius }
    }

    /// Point at angle `theta` (radians).
    #[inline]
    pub fn point_at(self, theta: f64) -> Point3 {
        let (s, c) = theta.sin_cos();
        self.frame.origin
            + self.frame.x * (c * self.radius)
            + self.frame.y * (s * self.radius)
    }

    /// Tangent direction at angle `theta` (unit vector, counter-clockwise).
    #[inline]
    pub fn tangent_at(self, theta: f64) -> UnitVec3 {
        let (s, c) = theta.sin_cos();
        UnitVec3::new_unchecked(-self.frame.x * s + self.frame.y * c)
    }

    /// Circumference.
    #[inline]
    pub fn circumference(self) -> f64 { TAU * self.radius }

    /// True if `p` lies on the circle (within tolerance `tol`).
    pub fn contains(self, p: Point3, tol: f64) -> bool {
        let local = self.frame.to_local_point(p);
        let r2 = local.x * local.x + local.y * local.y;
        local.z.abs() < tol && (r2.sqrt() - self.radius).abs() < tol
    }
}

// ── Ellipse3 ─────────────────────────────────────────────────────────────────

/// A planar ellipse.
///
/// Parameterised as `P(θ) = centre + x_axis·a·cos(θ) + y_axis·b·sin(θ)`.
/// In STEP terms: `ELLIPSE` entity with semi_axis_1 = `a` and semi_axis_2 = `b`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ellipse3 {
    /// Placement.  X = semi-major direction, Y = semi-minor direction.
    pub frame: Frame3,
    /// Semi-major axis length (mm).
    pub semi_major: f64,
    /// Semi-minor axis length (mm).
    pub semi_minor: f64,
}

impl Ellipse3 {
    /// Construct an ellipse at `centre` in a plane with the given normal.
    ///
    /// The major axis is `x_dir` (must be perpendicular to `normal`).
    pub fn new(
        centre:    Point3,
        normal:    UnitVec3,
        x_dir:     UnitVec3,
        semi_major: f64,
        semi_minor: f64,
    ) -> Self {
        let y_dir = UnitVec3::try_from_vec(normal.as_vec().cross(x_dir.as_vec()))
            .unwrap_or_else(|| normal.perp_basis().1);
        Self {
            frame: Frame3 { origin: centre, x: x_dir, y: y_dir, z: normal },
            semi_major,
            semi_minor,
        }
    }

    /// The miter ellipse for a cylinder of radius `r` cut by a plane whose
    /// normal makes an angle `half_angle` with the cylinder axis.
    ///
    /// This is the exact profile curve at the end of a cylinder segment when
    /// it meets the bisecting miter plane at a bend in a filament path.
    pub fn miter(radius: f64, half_angle: f64) -> (f64, f64) {
        // semi-minor = r, semi-major = r / cos(half_angle)
        let semi_major = if half_angle.cos().abs() > EPS {
            radius / half_angle.cos()
        } else {
            radius  // degenerate (180° turn)
        };
        (semi_major, radius)
    }

    /// Point at parameter `theta` (radians).
    #[inline]
    pub fn point_at(self, theta: f64) -> Point3 {
        let (s, c) = theta.sin_cos();
        self.frame.origin
            + self.frame.x * (c * self.semi_major)
            + self.frame.y * (s * self.semi_minor)
    }
}

// ── BezierCubic ──────────────────────────────────────────────────────────────

/// A degree-3 (cubic) Bézier curve in 3-D.
///
/// Used for approximation in contexts where NURBS export is preferred.
/// Not used in the primary STEP output path (which is all analytic).
#[derive(Clone, Copy, Debug)]
pub struct BezierCubic {
    /// Four control points: P0 (start), P1 (start tangent), P2 (end tangent), P3 (end).
    pub ctrl: [Point3; 4],
}

impl BezierCubic {
    /// Evaluate the curve at `t ∈ [0, 1]`.
    pub fn point_at(self, t: f64) -> Point3 {
        let _u = 1.0 - t;
        // De-Casteljau formula (numerically stable)
        let p01 = self.ctrl[0].lerp(self.ctrl[1], t);
        let p12 = self.ctrl[1].lerp(self.ctrl[2], t);
        let p23 = self.ctrl[2].lerp(self.ctrl[3], t);
        let p012 = p01.lerp(p12, t);
        let p123 = p12.lerp(p23, t);
        p012.lerp(p123, t)
    }

    /// Tangent direction at `t`.
    pub fn tangent_at(self, t: f64) -> Vec3 {
        let u = 1.0 - t;
        // Derivative = 3 * ((1-t)² * (P1-P0) + 2(1-t)t * (P2-P1) + t² * (P3-P2))
        let d0 = (self.ctrl[1] - self.ctrl[0]) * (3.0 * u * u);
        let d1 = (self.ctrl[2] - self.ctrl[1]) * (6.0 * u * t);
        let d2 = (self.ctrl[3] - self.ctrl[2]) * (3.0 * t * t);
        d0 + d1 + d2
    }

    /// Chord-length approximation of arc length using `n` subdivisions.
    pub fn approx_length(self, n: usize) -> f64 {
        let n = n.max(1);
        let mut len = 0.0;
        let mut prev = self.ctrl[0];
        for i in 1..=n {
            let p = self.point_at(i as f64 / n as f64);
            len += (p - prev).length();
            prev = p;
        }
        len
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circle_roundtrip() {
        let c = Circle3::new(Point3::ORIGIN, UnitVec3::Z, 2.0);
        let p = c.point_at(0.0);
        // The point must lie on the circle: radius-2 from centre, in the Z=0 plane.
        let dist = (p - Point3::ORIGIN).length();
        assert!((dist - 2.0).abs() < 1e-10, "dist={dist}");
        assert!(p.z.abs() < 1e-10, "z={}", p.z);
    }

    #[test]
    fn miter_ellipse_semi_axes() {
        let (a, b) = Ellipse3::miter(1.0, std::f64::consts::FRAC_PI_4);
        // half_angle = 45° => semi_major = 1/cos(45°) ≈ 1.4142
        assert!((a - std::f64::consts::SQRT_2).abs() < 1e-10);
        assert!((b - 1.0).abs() < 1e-10);
    }
}
