//! Right-handed coordinate frame (`Frame3`).
//!
//! A `Frame3` is an orthonormal basis + origin — a local coordinate system
//! embedded in world space.  It is the fundamental "placement" primitive used
//! throughout the CAD kernel (surface parameterisation, joint construction, etc.)

use crate::{Mat3, Point3, UnitVec3, Vec3};

/// Right-handed orthonormal frame: origin + three unit axes.
///
/// `x`, `y`, `z` form a right-handed basis: `x.cross(y) == z`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame3 {
    /// Frame origin in world space.
    pub origin: Point3,
    /// Local X axis.
    pub x: UnitVec3,
    /// Local Y axis.
    pub y: UnitVec3,
    /// Local Z axis (= x × y).
    pub z: UnitVec3,
}

impl Frame3 {
    /// World frame: origin at (0,0,0), axes aligned with world axes.
    pub const WORLD: Self = Self {
        origin: Point3::ORIGIN,
        x: UnitVec3::X,
        y: UnitVec3::Y,
        z: UnitVec3::Z,
    };

    /// Build a frame from an origin and a forward (+Z) direction.
    ///
    /// The X axis is chosen to be perpendicular to `forward` (arbitrary but
    /// consistent).  Y is derived as `z × x` (right-handed).
    pub fn from_origin_z(origin: Point3, forward: UnitVec3) -> Self {
        let (x, y) = forward.perp_basis();
        Self { origin, x, y, z: forward }
    }

    /// Build a frame from origin, forward (+Z), and a suggested up vector.
    ///
    /// `up_hint` need not be perpendicular to `forward`; it is projected to
    /// produce the Y axis.  Falls back to the arbitrary-perp method if `up_hint`
    /// is nearly parallel to `forward`.
    pub fn from_origin_z_up(origin: Point3, forward: UnitVec3, up_hint: Vec3) -> Self {
        let up_proj = up_hint - forward.as_vec() * forward.dot_vec(up_hint);
        let y = match UnitVec3::try_from_vec(up_proj) {
            Some(u) => u,
            None    => forward.perp_basis().1,  // fallback
        };
        let x = UnitVec3::try_from_vec(y.as_vec().cross(forward.as_vec()))
            .unwrap_or_else(|| forward.perp_basis().0);
        Self { origin, x, y, z: forward }
    }

    /// Express a world-space point in local coordinates.
    #[inline]
    pub fn to_local_point(self, p: Point3) -> Point3 {
        let d = p - self.origin;
        Point3::new(
            self.x.dot_vec(d),
            self.y.dot_vec(d),
            self.z.dot_vec(d),
        )
    }

    /// Express a local-space point in world coordinates.
    #[inline]
    pub fn to_world_point(self, p: Point3) -> Point3 {
        self.origin + self.x * p.x + self.y * p.y + self.z * p.z
    }

    /// Express a world-space vector in local coordinates.
    #[inline]
    pub fn to_local_vec(self, v: Vec3) -> Vec3 {
        Vec3::new(self.x.dot_vec(v), self.y.dot_vec(v), self.z.dot_vec(v))
    }

    /// Express a local-space vector in world coordinates.
    #[inline]
    pub fn to_world_vec(self, v: Vec3) -> Vec3 {
        self.x * v.x + self.y * v.y + self.z * v.z
    }

    /// Rotation matrix that maps world axes to local axes.
    #[inline]
    pub fn rotation(self) -> Mat3 {
        Mat3::from_axes(self.x, self.y, self.z)
    }

    /// Move the frame origin by `delta` (in world space).
    #[inline]
    pub fn translate(self, delta: Vec3) -> Self {
        Self { origin: self.origin + delta, ..self }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_point() {
        let frame = Frame3::from_origin_z(
            Point3::new(1.0, 2.0, 3.0),
            UnitVec3::try_from_vec(Vec3::new(1.0, 1.0, 1.0).normalize()).unwrap(),
        );
        let world_pt = Point3::new(5.0, 7.0, -2.0);
        let local_pt = frame.to_local_point(world_pt);
        let back     = frame.to_world_point(local_pt);
        assert!((world_pt - back).length() < 1e-10);
    }
}
