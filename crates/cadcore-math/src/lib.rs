//! # cadcore-math
//!
//! Fundamental numeric primitives for the **cadcore** CAD kernel.
//!
//! This crate has **zero external dependencies** — everything is built on
//! `std` and const-generic arithmetic.
//!
//! ## Design rules
//!
//! * All lengths are in **millimetres** unless noted otherwise.
//! * All angles are in **radians** unless noted otherwise.
//! * [`Point3`] and [`Vec3`] are distinct types — the compiler will not let
//!   you confuse a location with a direction.
//! * Every function is pure; no global mutable state.
//!
//! ## Types at a glance
//!
//! | Type | Description |
//! |---|---|
//! | [`Point3`] | A location in 3-D space |
//! | [`Vec3`] | A free vector (direction + magnitude) |
//! | [`UnitVec3`] | A `Vec3` guaranteed to have length 1 |
//! | [`Mat3`] | Column-major 3×3 matrix (rotations, linear maps) |
//! | [`Frame3`] | Right-handed orthonormal frame (origin + 3 axes) |
//! | [`Transform3`] | Rigid-body transform: rotation + translation |
//! | [`Interval`] | Closed real interval `[lo, hi]` |
//!
//! ## Quick example
//!
//! ```rust
//! use cadcore_math::{Point3, Vec3, UnitVec3, Frame3};
//!
//! let origin = Point3::new(1.0, 2.0, 3.0);
//! let dir    = UnitVec3::try_from_vec(Vec3::new(0.0, 0.0, 1.0)).unwrap();
//! let frame  = Frame3::from_origin_z(origin, dir);
//!
//! let world_pt = Point3::new(4.0, 5.0, 6.0);
//! let local_pt = frame.to_local_point(world_pt);
//! let back     = frame.to_world_point(local_pt);
//! assert!((world_pt - back).length() < 1e-10);
//! ```

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
