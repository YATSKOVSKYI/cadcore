//! Closed interval `[lo, hi]` on the real line.

use std::fmt;
use crate::EPS;

/// A closed real interval `[lo, hi]`.
#[derive(Clone, Copy, PartialEq)]
pub struct Interval {
    /// Lower bound.
    pub lo: f64,
    /// Upper bound.
    pub hi: f64,
}

impl Interval {
    /// The empty interval (lo > hi).
    pub const EMPTY: Self = Self { lo: f64::INFINITY, hi: f64::NEG_INFINITY };

    /// The whole real line.
    pub const UNIVERSE: Self = Self { lo: f64::NEG_INFINITY, hi: f64::INFINITY };

    /// Construct `[lo, hi]`.  Panics in debug if `lo > hi`.
    #[inline]
    pub fn new(lo: f64, hi: f64) -> Self {
        debug_assert!(lo <= hi, "Interval::new: lo ({lo}) > hi ({hi})");
        Self { lo, hi }
    }

    /// Construct from any two values, sorting them.
    #[inline]
    pub fn from_unordered(a: f64, b: f64) -> Self {
        Self { lo: a.min(b), hi: a.max(b) }
    }

    /// Smallest interval containing both `self` and `other`.
    #[inline]
    pub fn union(self, other: Self) -> Self {
        Self { lo: self.lo.min(other.lo), hi: self.hi.max(other.hi) }
    }

    /// Intersection (may be empty).
    #[inline]
    pub fn intersect(self, other: Self) -> Self {
        Self {
            lo: self.lo.max(other.lo),
            hi: self.hi.min(other.hi),
        }
    }

    /// Length of the interval.
    #[inline]
    pub fn length(self) -> f64 { (self.hi - self.lo).max(0.0) }

    /// `true` if the interval contains no points.
    #[inline]
    pub fn is_empty(self) -> bool { self.lo > self.hi + EPS }

    /// `true` if `v` is inside (inclusive with tolerance).
    #[inline]
    pub fn contains(self, v: f64) -> bool { v >= self.lo - EPS && v <= self.hi + EPS }

    /// Clamp `v` to the interval.
    #[inline]
    pub fn clamp(self, v: f64) -> f64 { v.max(self.lo).min(self.hi) }

    /// Midpoint.
    #[inline]
    pub fn midpoint(self) -> f64 { 0.5 * (self.lo + self.hi) }

    /// Expand both ends by `margin`.
    #[inline]
    pub fn expand(self, margin: f64) -> Self { Self::new(self.lo - margin, self.hi + margin) }
}

impl fmt::Debug for Interval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{:.6}, {:.6}]", self.lo, self.hi)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn union_span() {
        let a = Interval::new(0.0, 3.0);
        let b = Interval::new(2.0, 5.0);
        let u = a.union(b);
        assert_eq!(u.lo, 0.0);
        assert_eq!(u.hi, 5.0);
    }

    #[test]
    fn intersection_empty() {
        let a = Interval::new(0.0, 1.0);
        let b = Interval::new(2.0, 3.0);
        assert!(a.intersect(b).is_empty());
    }
}
