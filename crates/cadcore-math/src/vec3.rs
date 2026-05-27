//! 3-D vector (`Vec3`).

use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::EPS;

/// A free 3-D vector with f64 components.
///
/// Distinct from [`crate::Point3`]: a `Vec3` represents a *direction + magnitude*,
/// while a `Point3` represents a *location*.  The two types are deliberately
/// separated so the compiler catches confusion between them.
#[derive(Clone, Copy, PartialEq)]
#[repr(C)]
pub struct Vec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl Vec3 {
    /// Zero vector.
    pub const ZERO: Self = Self { x: 0.0, y: 0.0, z: 0.0 };
    /// Unit vector along +X.
    pub const X: Self = Self { x: 1.0, y: 0.0, z: 0.0 };
    /// Unit vector along +Y.
    pub const Y: Self = Self { x: 0.0, y: 1.0, z: 0.0 };
    /// Unit vector along +Z.
    pub const Z: Self = Self { x: 0.0, y: 0.0, z: 1.0 };

    /// Construct from components.
    #[inline] pub const fn new(x: f64, y: f64, z: f64) -> Self { Self { x, y, z } }

    /// Squared Euclidean length.
    #[inline] pub fn length_sq(self) -> f64 { self.dot(self) }

    /// Euclidean length.
    #[inline] pub fn length(self) -> f64 { self.length_sq().sqrt() }

    /// Dot product.
    #[inline]
    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    /// Cross product: `self × rhs`.
    #[inline]
    pub fn cross(self, rhs: Self) -> Self {
        Self {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }

    /// Return the normalised version of this vector, or `None` if near-zero.
    #[inline]
    pub fn try_normalize(self) -> Option<Self> {
        let len = self.length();
        if len < EPS { None } else { Some(self / len) }
    }

    /// Normalise without a safety check (panics on zero in debug builds).
    #[inline]
    pub fn normalize(self) -> Self {
        self.try_normalize().expect("Vec3::normalize called on zero vector")
    }

    /// Component-wise absolute value.
    #[inline]
    pub fn abs(self) -> Self { Self::new(self.x.abs(), self.y.abs(), self.z.abs()) }

    /// Reflect `self` about a unit normal `n`.
    #[inline]
    pub fn reflect(self, n: Self) -> Self { self - n * (2.0 * self.dot(n)) }

    /// Return a vector perpendicular to `self` (arbitrary, consistent).
    pub fn any_perp(self) -> Self {
        // Choose the axis with the smallest absolute component to maximise
        // numerical stability of the cross product.
        if self.x.abs() <= self.y.abs() && self.x.abs() <= self.z.abs() {
            self.cross(Self::X)
        } else if self.y.abs() <= self.z.abs() {
            self.cross(Self::Y)
        } else {
            self.cross(Self::Z)
        }
    }

    /// Linearly interpolate: `(1-t)*self + t*rhs`.
    #[inline]
    pub fn lerp(self, rhs: Self, t: f64) -> Self { self + (rhs - self) * t }
}

// ── Arithmetic operators ─────────────────────────────────────────────────────

impl Add for Vec3 {
    type Output = Self;
    #[inline] fn add(self, r: Self) -> Self { Self::new(self.x+r.x, self.y+r.y, self.z+r.z) }
}
impl Sub for Vec3 {
    type Output = Self;
    #[inline] fn sub(self, r: Self) -> Self { Self::new(self.x-r.x, self.y-r.y, self.z-r.z) }
}
impl Neg for Vec3 {
    type Output = Self;
    #[inline] fn neg(self) -> Self { Self::new(-self.x, -self.y, -self.z) }
}
impl Mul<f64> for Vec3 {
    type Output = Self;
    #[inline] fn mul(self, s: f64) -> Self { Self::new(self.x*s, self.y*s, self.z*s) }
}
impl Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline] fn mul(self, v: Vec3) -> Vec3 { v * self }
}
impl Div<f64> for Vec3 {
    type Output = Self;
    #[inline] fn div(self, s: f64) -> Self { self * (1.0 / s) }
}

impl AddAssign for Vec3 { #[inline] fn add_assign(&mut self, r: Self) { *self = *self + r; } }
impl SubAssign for Vec3 { #[inline] fn sub_assign(&mut self, r: Self) { *self = *self - r; } }
impl MulAssign<f64> for Vec3 { #[inline] fn mul_assign(&mut self, s: f64) { *self = *self * s; } }
impl DivAssign<f64> for Vec3 { #[inline] fn div_assign(&mut self, s: f64) { *self = *self / s; } }

impl fmt::Debug for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Vec3({:.6}, {:.6}, {:.6})", self.x, self.y, self.z)
    }
}
impl fmt::Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.4}, {:.4}, {:.4})", self.x, self.y, self.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cross_product() {
        let x = Vec3::X;
        let y = Vec3::Y;
        let z = x.cross(y);
        assert!((z - Vec3::Z).length() < 1e-12);
    }

    #[test]
    fn normalize_unit_length() {
        let v = Vec3::new(3.0, 4.0, 0.0).normalize();
        assert!((v.length() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn reflect() {
        let v = Vec3::new(1.0, -1.0, 0.0);
        let n = Vec3::Y;
        let r = v.reflect(n);
        assert!((r - Vec3::new(1.0, 1.0, 0.0)).length() < 1e-12);
    }
}
