//! Analytic surface types.
//!
//! Each surface has a canonical `(u, v)` parameterisation and maps directly
//! to a STEP AP203 surface entity.

use cadcore_math::{Frame3, Point3, UnitVec3};

// ── Plane3 ───────────────────────────────────────────────────────────────────

/// An infinite plane.
///
/// STEP entity: `PLANE`.
/// Parameterised as `P(u,v) = origin + u·x + v·y`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane3 {
    /// Frame: origin on plane, Z = outward normal, X/Y span the plane.
    pub frame: Frame3,
}

impl Plane3 {
    /// Construct a plane from origin and normal.
    pub fn from_origin_normal(origin: Point3, normal: UnitVec3) -> Self {
        Self { frame: Frame3::from_origin_z(origin, normal) }
    }

    /// Plane's outward normal.
    #[inline] pub fn normal(self) -> UnitVec3 { self.frame.z }

    /// Signed distance from point `p` to the plane (positive = normal side).
    #[inline]
    pub fn signed_distance(self, p: Point3) -> f64 {
        self.normal().dot_vec(p - self.frame.origin)
    }

    /// Project `p` onto the plane.
    #[inline]
    pub fn project(self, p: Point3) -> Point3 {
        p - self.normal() * self.signed_distance(p)
    }

    /// Evaluate the parameterised surface at `(u, v)`.
    #[inline]
    pub fn point_at(self, u: f64, v: f64) -> Point3 {
        self.frame.origin + self.frame.x * u + self.frame.y * v
    }
}

// ── CylSurf ──────────────────────────────────────────────────────────────────

/// Right circular cylinder (infinite in the axis direction).
///
/// STEP entity: `CYLINDRICAL_SURFACE`.
/// Parameterised as:
/// ```text
/// P(θ, z) = axis_origin + x·r·cos(θ) + y·r·sin(θ) + axis·z
/// ```
/// `θ ∈ [0, 2π)`, `z` unbounded.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CylSurf {
    /// Frame: origin on axis, Z = axis direction, X = 0° reference.
    pub frame: Frame3,
    /// Cylinder radius (mm).
    pub radius: f64,
}

impl CylSurf {
    /// Construct from an axis origin, axis direction, and radius.
    pub fn new(axis_origin: Point3, axis_dir: UnitVec3, radius: f64) -> Self {
        Self { frame: Frame3::from_origin_z(axis_origin, axis_dir), radius }
    }

    /// Axis direction (unit vector).
    #[inline] pub fn axis(self) -> UnitVec3 { self.frame.z }

    /// Evaluate at `(theta, z)`.
    #[inline]
    pub fn point_at(self, theta: f64, z: f64) -> Point3 {
        let (s, c) = theta.sin_cos();
        self.frame.origin
            + self.frame.x * (c * self.radius)
            + self.frame.y * (s * self.radius)
            + self.frame.z * z
    }

    /// Outward normal at `(theta, z)`.
    #[inline]
    pub fn normal_at(self, theta: f64, _z: f64) -> UnitVec3 {
        let (s, c) = theta.sin_cos();
        UnitVec3::new_unchecked(self.frame.x * c + self.frame.y * s)
    }

    /// The axis line.
    #[inline]
    pub fn axis_line(self) -> cadcore_math::Frame3 { self.frame }
}

// ── ConeSurf ─────────────────────────────────────────────────────────────────

/// Right circular cone.
///
/// STEP entity: `CONICAL_SURFACE`.
/// Apex at `frame.origin`, axis along `frame.z`, half-angle `semi_angle`.
/// Parameterised as:
/// ```text
/// P(θ, z) = apex + (x·cos(θ) + y·sin(θ))·z·tan(semi_angle) + axis·z
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConeSurf {
    /// Frame: origin = apex, Z = axis, X = 0° reference.
    pub frame: Frame3,
    /// Half-angle at the apex (radians).
    pub semi_angle: f64,
}

impl ConeSurf {
    /// Construct from apex, axis direction, and semi-angle.
    pub fn new(apex: Point3, axis_dir: UnitVec3, semi_angle: f64) -> Self {
        Self { frame: Frame3::from_origin_z(apex, axis_dir), semi_angle }
    }

    /// Radius at axial distance `z`.
    #[inline] pub fn radius_at_z(self, z: f64) -> f64 { z * self.semi_angle.tan() }

    /// Point at `(theta, z)`.
    #[inline]
    pub fn point_at(self, theta: f64, z: f64) -> Point3 {
        let r = self.radius_at_z(z);
        let (s, c) = theta.sin_cos();
        self.frame.origin
            + self.frame.x * (c * r)
            + self.frame.y * (s * r)
            + self.frame.z * z
    }
}

// ── SphereSurf ───────────────────────────────────────────────────────────────

/// Sphere.
///
/// STEP entity: `SPHERICAL_SURFACE`.
/// Parameterised as:
/// ```text
/// P(θ, φ) = centre + r·(cos(φ)·cos(θ)·x + cos(φ)·sin(θ)·y + sin(φ)·z)
/// ```
/// `θ ∈ [0, 2π)`, `φ ∈ [-π/2, π/2]` (latitude).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SphereSurf {
    /// Centre.
    pub centre: Point3,
    /// Radius (mm).
    pub radius: f64,
    /// Orientation frame.
    pub frame: Frame3,
}

impl SphereSurf {
    /// Construct a sphere at `centre` with `radius`.
    pub fn new(centre: Point3, radius: f64) -> Self {
        Self {
            centre,
            radius,
            frame: Frame3::from_origin_z(centre, UnitVec3::Z),
        }
    }

    /// Point at longitude `theta` (0-based) and latitude `phi`.
    #[inline]
    pub fn point_at(self, theta: f64, phi: f64) -> Point3 {
        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();
        self.centre
            + self.frame.x * (cp * ct * self.radius)
            + self.frame.y * (cp * st * self.radius)
            + self.frame.z * (sp * self.radius)
    }
}

// ── TorusSurf ────────────────────────────────────────────────────────────────

/// Ring torus (donut shape).
///
/// STEP entity: `TOROIDAL_SURFACE`.
/// Parameterised as:
/// ```text
/// P(θ, φ) = (R + r·cos(φ))·(cos(θ)·x + sin(θ)·y) + r·sin(φ)·z
/// ```
/// `θ` = rotation around torus axis, `φ` = angle around the tube.
/// * `major_radius R` — distance from torus centre to tube centre
/// * `minor_radius r` — tube radius (= filament radius for corners)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TorusSurf {
    /// Frame: origin = torus centre, Z = axis, X = 0° reference.
    pub frame: Frame3,
    /// Distance from the torus centre-line to the tube centre (mm).
    pub major_radius: f64,
    /// Tube cross-section radius (= filament radius, mm).
    pub minor_radius: f64,
}

impl TorusSurf {
    /// Construct with centre, axis, and both radii.
    pub fn new(centre: Point3, axis: UnitVec3, major_radius: f64, minor_radius: f64) -> Self {
        Self {
            frame: Frame3::from_origin_z(centre, axis),
            major_radius,
            minor_radius,
        }
    }

    /// Build the torus used at a filament path corner.
    ///
    /// * `vertex`         — the bend vertex in the path
    /// * `incoming_dir`   — normalised direction *arriving at* the vertex
    /// * `outgoing_dir`   — normalised direction *leaving* the vertex
    /// * `filament_radius`
    ///
    /// For a 90° corner with radius `r`, the major radius is `r` and the
    /// minor radius is the filament cross-section radius.
    pub fn for_corner(
        vertex:         Point3,
        incoming_dir:   UnitVec3,
        outgoing_dir:   UnitVec3,
        filament_radius: f64,
    ) -> Option<Self> {
        // Bisector of the corner
        let bisect = UnitVec3::try_from_vec(
            incoming_dir.as_vec() + outgoing_dir.as_vec()
        )?;

        // Torus axis is perpendicular to both dirs (the bend plane normal)
        let axis = UnitVec3::try_from_vec(
            incoming_dir.cross(outgoing_dir)
        )?;

        // Half-bend angle (β = angle between bisector and either incoming/outgoing dir)
        let cos_beta = bisect.dot(outgoing_dir);
        let major_radius = if cos_beta.abs() > 1e-6 {
            filament_radius / cos_beta  // ensures the torus arc meets the cylinders tangentially
        } else {
            filament_radius
        };

        // Arc centre is along the bisector
        let centre = vertex + bisect * major_radius;

        Some(Self::new(centre, axis, major_radius, filament_radius))
    }

    /// Point at `(theta, phi)`.
    #[inline]
    pub fn point_at(self, theta: f64, phi: f64) -> Point3 {
        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();
        let spine_point = self.frame.origin
            + self.frame.x * (ct * self.major_radius)
            + self.frame.y * (st * self.major_radius);
        // Radial direction at angle theta
        let radial = self.frame.x * ct + self.frame.y * st;
        spine_point + radial * (cp * self.minor_radius) + self.frame.z * (sp * self.minor_radius)
    }

    /// Outward surface normal at `(theta, phi)`.
    #[inline]
    pub fn normal_at(self, theta: f64, phi: f64) -> UnitVec3 {
        let (st, ct) = theta.sin_cos();
        let (sp, cp) = phi.sin_cos();
        let radial = self.frame.x * ct + self.frame.y * st;
        UnitVec3::new_unchecked(radial * cp + self.frame.z * sp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_math::PI;

    #[test]
    fn cylinder_point_on_surface() {
        let cyl = CylSurf::new(Point3::ORIGIN, UnitVec3::Z, 2.0);
        let p = cyl.point_at(0.0, 5.0);
        // The point must be at radius 2 from the Z-axis and at z=5.
        let dist_from_axis = (p.x * p.x + p.y * p.y).sqrt();
        assert!((dist_from_axis - 2.0).abs() < 1e-10, "radial dist={dist_from_axis}");
        assert!((p.z - 5.0).abs() < 1e-10, "z={}", p.z);
    }

    #[test]
    fn torus_minor_circle() {
        let torus = TorusSurf::new(Point3::ORIGIN, UnitVec3::Z, 3.0, 1.0);
        // At theta=0 the spine point is at distance R from the Z-axis.
        // phi=0 → outermost point: distance R+r from axis.
        // phi=π → innermost point: distance R-r from axis.
        let p0 = torus.point_at(0.0, 0.0);
        let p1 = torus.point_at(0.0, PI);
        let d0 = (p0.x * p0.x + p0.y * p0.y).sqrt();
        let d1 = (p1.x * p1.x + p1.y * p1.y).sqrt();
        assert!((d0 - 4.0).abs() < 1e-10, "outer dist={d0}");
        assert!((d1 - 2.0).abs() < 1e-10, "inner dist={d1}");
    }
}
