//! Rigid-body transformation (`Transform3`): rotation + translation.
//!
//! Represents the map  `p ↦ R·p + t`  where `R` is an orthogonal 3×3 matrix
//! and `t` is a translation vector.

use std::ops::Mul;
use crate::{Mat3, Point3, Vec3};

/// A rigid (isometric) transformation in 3-D space.
///
/// Composes as: `(B ∘ A).apply(p) == B.apply(A.apply(p))`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Transform3 {
    /// Linear (rotation / reflection) part.
    pub rotation:    Mat3,
    /// Translation part (applied *after* the rotation).
    pub translation: Vec3,
}

impl Transform3 {
    /// Identity transform.
    pub const IDENTITY: Self = Self {
        rotation:    Mat3::IDENTITY,
        translation: Vec3::ZERO,
    };

    /// Pure translation.
    #[inline]
    pub fn from_translation(t: Vec3) -> Self {
        Self { rotation: Mat3::IDENTITY, translation: t }
    }

    /// Pure rotation around `axis` by `angle` radians.
    #[inline]
    pub fn from_rotation(axis: crate::UnitVec3, angle: f64) -> Self {
        Self { rotation: Mat3::rotation(axis, angle), translation: Vec3::ZERO }
    }

    /// Apply the transform to a point.
    #[inline]
    pub fn apply_point(self, p: Point3) -> Point3 {
        Point3::from_vec(self.rotation * p.to_vec()) + self.translation
    }

    /// Apply only the rotation part to a vector (no translation).
    #[inline]
    pub fn apply_vec(self, v: Vec3) -> Vec3 {
        self.rotation * v
    }

    /// Inverse transform.
    pub fn inverse(self) -> Self {
        let r_inv = self.rotation.transpose(); // orthogonal ⟹ R⁻¹ = Rᵀ
        let t_inv = -(r_inv * self.translation);
        Self { rotation: r_inv, translation: t_inv }
    }

    /// Compose two transforms: `self ∘ rhs` (apply `rhs` first).
    #[inline]
    pub fn compose(self, rhs: Self) -> Self {
        Self {
            rotation:    self.rotation * rhs.rotation,
            translation: self.rotation * rhs.translation + self.translation,
        }
    }
}

/// Composition operator: `lhs * rhs` applies `rhs` first.
impl Mul for Transform3 {
    type Output = Self;
    #[inline] fn mul(self, rhs: Self) -> Self { self.compose(rhs) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::UnitVec3;

    #[test]
    fn inverse_round_trip() {
        let t = Transform3 {
            rotation:    Mat3::rotation(UnitVec3::Z, 0.7),
            translation: Vec3::new(1.0, 2.0, 3.0),
        };
        let p = Point3::new(4.0, 5.0, 6.0);
        let q = t.inverse().apply_point(t.apply_point(p));
        assert!((p - q).length() < 1e-10);
    }
}
