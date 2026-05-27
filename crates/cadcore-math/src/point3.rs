//! 3-D point (`Point3`).

use std::fmt;
use std::ops::{Add, AddAssign, Sub};

use crate::Vec3;

/// A location in 3-D space.
///
/// Separate from [`Vec3`] so the type system enforces correct usage
/// (e.g. `Point3 - Point3 → Vec3`, `Point3 + Vec3 → Point3`).
#[derive(Clone, Copy, PartialEq)]
#[repr(C)]
pub struct Point3 {
    /// X coordinate (mm).
    pub x: f64,
    /// Y coordinate (mm).
    pub y: f64,
    /// Z coordinate (mm).
    pub z: f64,
}

impl Point3 {
    /// Origin (0, 0, 0).
    pub const ORIGIN: Self = Self { x: 0.0, y: 0.0, z: 0.0 };

    /// Construct from coordinates.
    #[inline] pub const fn new(x: f64, y: f64, z: f64) -> Self { Self { x, y, z } }

    /// Displacement vector from `other` to `self`: `self - other`.
    #[inline]
    pub fn offset_from(self, other: Self) -> Vec3 {
        Vec3::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    /// Euclidean distance to `other`.
    #[inline]
    pub fn distance_to(self, other: Self) -> f64 { self.offset_from(other).length() }

    /// Linearly interpolate between two points: `(1-t)*self + t*other`.
    #[inline]
    pub fn lerp(self, other: Self, t: f64) -> Self {
        Self::new(
            self.x + t * (other.x - self.x),
            self.y + t * (other.y - self.y),
            self.z + t * (other.z - self.z),
        )
    }

    /// Convert to a vector from the origin.
    #[inline]
    pub fn to_vec(self) -> Vec3 { Vec3::new(self.x, self.y, self.z) }

    /// Construct from a vector (treating it as an offset from origin).
    #[inline]
    pub fn from_vec(v: Vec3) -> Self { Self::new(v.x, v.y, v.z) }
}

// ── Arithmetic ───────────────────────────────────────────────────────────────

/// `Point3 + Vec3 → Point3`
impl Add<Vec3> for Point3 {
    type Output = Self;
    #[inline]
    fn add(self, v: Vec3) -> Self {
        Self::new(self.x + v.x, self.y + v.y, self.z + v.z)
    }
}
impl AddAssign<Vec3> for Point3 {
    #[inline] fn add_assign(&mut self, v: Vec3) { *self = *self + v; }
}

/// `Point3 - Point3 → Vec3`
impl Sub for Point3 {
    type Output = Vec3;
    #[inline]
    fn sub(self, other: Self) -> Vec3 { self.offset_from(other) }
}

/// `Point3 - Vec3 → Point3`
impl Sub<Vec3> for Point3 {
    type Output = Self;
    #[inline]
    fn sub(self, v: Vec3) -> Self { Self::new(self.x - v.x, self.y - v.y, self.z - v.z) }
}

impl fmt::Debug for Point3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Point3({:.6}, {:.6}, {:.6})", self.x, self.y, self.z)
    }
}
impl fmt::Display for Point3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.4}, {:.4}, {:.4})", self.x, self.y, self.z)
    }
}

// ── Conversions ──────────────────────────────────────────────────────────────

impl From<[f64; 3]> for Point3 {
    fn from(a: [f64; 3]) -> Self { Self::new(a[0], a[1], a[2]) }
}
impl From<Point3> for [f64; 3] {
    fn from(p: Point3) -> [f64; 3] { [p.x, p.y, p.z] }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subtraction_gives_vec() {
        let a = Point3::new(3.0, 0.0, 0.0);
        let b = Point3::ORIGIN;
        let v = a - b;
        assert_eq!(v, Vec3::new(3.0, 0.0, 0.0));
    }

    #[test]
    fn add_vec_moves_point() {
        let p = Point3::ORIGIN + Vec3::Z;
        assert_eq!(p, Point3::new(0.0, 0.0, 1.0));
    }
}
