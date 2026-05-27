//! Unit vector (`UnitVec3`) — a `Vec3` guaranteed to have length 1.

use std::fmt;
use std::ops::{Mul, Neg};

use crate::{Vec3, EPS};

/// A `Vec3` that is normalised to unit length at construction time.
///
/// Stored as a plain `Vec3` internally — the invariant is only enforced
/// at the construction boundary.
#[derive(Clone, Copy, PartialEq)]
pub struct UnitVec3(Vec3);

impl UnitVec3 {
    /// +X axis.
    pub const X: Self = Self(Vec3::X);
    /// +Y axis.
    pub const Y: Self = Self(Vec3::Y);
    /// +Z axis.
    pub const Z: Self = Self(Vec3::Z);

    /// Construct from a `Vec3`, normalising it.  Returns `None` for zero vectors.
    #[inline]
    pub fn try_from_vec(v: Vec3) -> Option<Self> {
        let len = v.length();
        if len < EPS { None } else { Some(Self(v / len)) }
    }

    /// Construct without checking — caller asserts `|v| ≈ 1`.
    ///
    /// # Safety
    /// Violating the unit-length invariant leads to incorrect geometric results,
    /// not unsound memory, so this is not `unsafe` in the Rust sense.
    #[inline]
    pub fn new_unchecked(v: Vec3) -> Self { Self(v) }

    /// Access the underlying `Vec3`.
    #[inline]
    pub fn as_vec(self) -> Vec3 { self.0 }

    /// Dot product with another unit vector.
    #[inline]
    pub fn dot(self, rhs: Self) -> f64 { self.0.dot(rhs.0) }

    /// Dot product with a free `Vec3`.
    #[inline]
    pub fn dot_vec(self, v: Vec3) -> f64 { self.0.dot(v) }

    /// Cross product — result is a plain `Vec3` (may not be unit length).
    #[inline]
    pub fn cross(self, rhs: Self) -> Vec3 { self.0.cross(rhs.0) }

    /// Return the opposite direction.
    #[inline]
    pub fn flip(self) -> Self { Self(-self.0) }

    /// Angle between two unit vectors (radians, always in `[0, π]`).
    #[inline]
    pub fn angle_to(self, other: Self) -> f64 {
        self.dot(other).clamp(-1.0, 1.0).acos()
    }

    /// Build an orthonormal basis `(self, u, v)` where `u` and `v` span the
    /// plane perpendicular to `self`.  `u` is arbitrary but consistent.
    pub fn perp_basis(self) -> (Self, Self) {
        let u = self.0.any_perp().normalize();
        let v_vec = self.0.cross(u);
        (Self(u), Self(v_vec))
    }
}

impl Neg for UnitVec3 {
    type Output = Self;
    #[inline] fn neg(self) -> Self { self.flip() }
}

impl Mul<f64> for UnitVec3 {
    type Output = Vec3;
    #[inline] fn mul(self, s: f64) -> Vec3 { self.0 * s }
}
impl Mul<UnitVec3> for f64 {
    type Output = Vec3;
    #[inline] fn mul(self, u: UnitVec3) -> Vec3 { u.0 * self }
}

impl From<UnitVec3> for Vec3 {
    fn from(u: UnitVec3) -> Vec3 { u.0 }
}

impl fmt::Debug for UnitVec3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "UnitVec3({:.6}, {:.6}, {:.6})", self.0.x, self.0.y, self.0.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_length() {
        let v = Vec3::new(1.0, 2.0, 3.0);
        let u = UnitVec3::try_from_vec(v).unwrap();
        assert!((u.as_vec().length() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn zero_returns_none() {
        assert!(UnitVec3::try_from_vec(Vec3::ZERO).is_none());
    }

    #[test]
    fn perp_basis_orthogonal() {
        let u = UnitVec3::Z;
        let (a, b) = u.perp_basis();
        assert!(a.dot(u).abs() < 1e-12);
        assert!(b.dot(u).abs() < 1e-12);
        assert!(a.dot(b).abs() < 1e-12);
    }
}
