//! # cadcore-ops
//!
//! Geometric operations that build B-Rep topology from high-level inputs.
//!
//! The main entry point is [`sweep::sweep_circle_along_polyline`], which takes
//! a sequence of 3-D waypoints and a radius and returns an exact B-Rep solid —
//! no Boolean operations, O(N) construction.

#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod sweep;

pub use sweep::{sweep_circle_along_polyline, SweepOptions};
