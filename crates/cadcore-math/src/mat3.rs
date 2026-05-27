//! 3×3 matrix (`Mat3`) for rotations and linear maps.

use std::fmt;
use std::ops::{Mul, MulAssign};

use crate::{Vec3, UnitVec3, EPS};

/// Column-major 3×3 matrix.
///
/// Columns are stored as `col[0..2]`, so `mat * v` is `sum_i col[i]*v[i]`.
#[derive(Clone, Copy, PartialEq)]
pub struct Mat3 {
    /// Columns: `col[j]` is the j-th column.
    pub col: [Vec3; 3],
}

impl Mat3 {
    /// Identity matrix.
    pub const IDENTITY: Self = Self {
        col: [Vec3::X, Vec3::Y, Vec3::Z],
    };

    /// Zero matrix.
    pub const ZERO: Self = Self {
        col: [Vec3::ZERO, Vec3::ZERO, Vec3::ZERO],
    };

    /// Construct from columns.
    #[inline]
    pub const fn from_cols(c0: Vec3, c1: Vec3, c2: Vec3) -> Self {
        Self { col: [c0, c1, c2] }
    }

    /// Construct from rows.
    #[inline]
    pub fn from_rows(r0: Vec3, r1: Vec3, r2: Vec3) -> Self {
        Self::from_cols(
            Vec3::new(r0.x, r1.x, r2.x),
            Vec3::new(r0.y, r1.y, r2.y),
            Vec3::new(r0.z, r1.z, r2.z),
        )
    }

    /// Rotation matrix that maps the +Z axis to `dir` and +X axis to `right`,
    /// forming a right-handed orthonormal frame.
    ///
    /// `right` must be perpendicular to `dir`.  Use [`UnitVec3::perp_basis`]
    /// to construct a consistent `right` if you don't have one.
    pub fn from_axes(right: UnitVec3, up: UnitVec3, forward: UnitVec3) -> Self {
        Self::from_cols(right.as_vec(), up.as_vec(), forward.as_vec())
    }

    /// Rotation by `angle` radians around `axis` (must be unit length).
    pub fn rotation(axis: UnitVec3, angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;
        let Vec3 { x, y, z } = axis.as_vec();
        Self::from_rows(
            Vec3::new(t*x*x + c,   t*x*y - s*z, t*x*z + s*y),
            Vec3::new(t*x*y + s*z, t*y*y + c,   t*y*z - s*x),
            Vec3::new(t*x*z - s*y, t*y*z + s*x, t*z*z + c  ),
        )
    }

    /// Transpose.
    #[inline]
    pub fn transpose(self) -> Self {
        Self::from_rows(self.col[0], self.col[1], self.col[2])
    }

    /// Determinant.
    pub fn det(self) -> f64 {
        let [c0, c1, c2] = self.col;
        c0.dot(c1.cross(c2))
    }

    /// Inverse (returns `None` if singular).
    ///
    /// Uses the cofactor/adjugate method: `M⁻¹ = adj(M)ᵀ / det`.
    /// The three `Vec3`s below are the *columns* of the inverse.
    pub fn try_inverse(self) -> Option<Self> {
        let d = self.det();
        if d.abs() < EPS { return None; }
        let inv_d = 1.0 / d;
        let [c0, c1, c2] = self.col;

        // Each variable is a *column* of the inverse (= a column of adj(M) / det).
        // Column 0 of adj(M) = cofactors of row 0: C_{00}, C_{10}, C_{20}
        let ic0 = Vec3::new(
            ( c1.y * c2.z - c1.z * c2.y) * inv_d,  //  C_{00} = +(e·i - f·h)
            -(c0.y * c2.z - c0.z * c2.y) * inv_d,  //  C_{10} = -(b·i - c·h) wait...
            ( c0.y * c1.z - c0.z * c1.y) * inv_d,  //  C_{20} = +(b·f - c·e)
        );
        // Column 1 of adj(M)
        let ic1 = Vec3::new(
            -(c1.x * c2.z - c1.z * c2.x) * inv_d,  //  C_{01} = -(d·i - f·g)
            ( c0.x * c2.z - c0.z * c2.x) * inv_d,  //  C_{11} = +(a·i - c·g)
            -(c0.x * c1.z - c0.z * c1.x) * inv_d,  //  C_{21} = -(a·f - c·d)
        );
        // Column 2 of adj(M)
        let ic2 = Vec3::new(
            ( c1.x * c2.y - c1.y * c2.x) * inv_d,  //  C_{02} = +(d·h - e·g)
            -(c0.x * c2.y - c0.y * c2.x) * inv_d,  //  C_{12} = -(a·h - b·g)
            ( c0.x * c1.y - c0.y * c1.x) * inv_d,  //  C_{22} = +(a·e - b·d)
        );
        Some(Self::from_cols(ic0, ic1, ic2))
    }

    /// Apply the matrix to a vector: `M * v`.
    #[inline]
    pub fn apply(self, v: Vec3) -> Vec3 {
        self.col[0] * v.x + self.col[1] * v.y + self.col[2] * v.z
    }
}

impl Mul for Mat3 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self::from_cols(
            self.apply(rhs.col[0]),
            self.apply(rhs.col[1]),
            self.apply(rhs.col[2]),
        )
    }
}

impl MulAssign for Mat3 {
    fn mul_assign(&mut self, rhs: Self) { *self = *self * rhs; }
}

impl Mul<Vec3> for Mat3 {
    type Output = Vec3;
    #[inline] fn mul(self, v: Vec3) -> Vec3 { self.apply(v) }
}

impl fmt::Debug for Mat3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let [c0, c1, c2] = self.col;
        write!(f,
            "Mat3[\n  [{:.4} {:.4} {:.4}]\n  [{:.4} {:.4} {:.4}]\n  [{:.4} {:.4} {:.4}]\n]",
            c0.x, c1.x, c2.x,
            c0.y, c1.y, c2.y,
            c0.z, c1.z, c2.z,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PI;

    #[test]
    fn rotation_90_deg_around_z() {
        let m = Mat3::rotation(UnitVec3::Z, PI / 2.0);
        let v = m * Vec3::X;
        assert!((v - Vec3::Y).length() < 1e-12);
    }

    #[test]
    fn inverse_times_original_is_identity() {
        let m = Mat3::rotation(
            UnitVec3::try_from_vec(Vec3::new(1.0, 2.0, 3.0).normalize()).unwrap(),
            1.23,
        );
        let inv = m.try_inverse().unwrap();
        let prod = m * inv;
        for i in 0..3 {
            let diff = prod.col[i] - Mat3::IDENTITY.col[i];
            assert!(diff.length() < 1e-10, "column {i}: {:?}", diff);
        }
    }
}
