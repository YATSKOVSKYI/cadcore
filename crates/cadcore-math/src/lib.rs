//! # cadcore-math
//!
//! Fundamental numeric primitives for the cadcore CAD kernel.
//!
//! **Design rules**
//! * Zero external dependencies — everything is `std` + const-generic arithmetic.
//! * All lengths are in **millimetres** unless noted otherwise.
//! * Angles are in **radians** unless noted otherwise.
//! * No global state; every function is pure or takes explicit parameters.

#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

mod interval;
mod mat3;
mod point3;
mod transform;
mod unitvec;
mod vec3;
mod frame;

pub use interval::Interval;
pub use mat3::Mat3;
pub use point3::Point3;
pub use transform::Transform3;
pub use unitvec::UnitVec3;
pub use vec3::Vec3;
pub use frame::Frame3;

/// Machine epsilon for `f64` comparisons.
pub const EPS: f64 = 1e-10;
/// π
pub const PI: f64 = std::f64::consts::PI;
/// 2π
pub const TAU: f64 = std::f64::consts::TAU;

/// Return `true` when two `f64` values are within [`EPS`] of each other.
#[inline]
pub fn approx_eq(a: f64, b: f64) -> bool { (a - b).abs() < EPS }

/// Clamp `v` into `[lo, hi]`.
#[inline]
pub fn clamp(v: f64, lo: f64, hi: f64) -> f64 { v.max(lo).min(hi) }

/// Linearly interpolate: `(1-t)*a + t*b`.
#[inline]
pub fn lerp(a: f64, b: f64, t: f64) -> f64 { a + t * (b - a) }
