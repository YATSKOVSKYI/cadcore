//! # cadcore-geom
//!
//! Analytic geometry primitives for the cadcore CAD kernel.
//!
//! ## Curves
//! * [`Line3`]     — infinite line through a point with a direction
//! * [`Circle3`]   — planar circle defined by a frame and radius
//! * [`Ellipse3`]  — planar ellipse defined by a frame and two radii
//! * [`BezierCubic`] — degree-3 Bézier curve (for approximation use-cases)
//!
//! ## Surfaces
//! * [`Plane3`]    — infinite plane (normal + offset)
//! * [`CylSurf`]   — right circular cylinder
//! * [`ConeSurf`]  — right circular cone
//! * [`SphereSurf`]— full sphere
//! * [`TorusSurf`] — ring torus
//!
//! All analytic surfaces have `point_at(u, v)` methods and can write exact
//! STEP AP203 entities via `cadcore-step`.

#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod curves;
pub mod surfaces;

pub use curves::{BezierCubic, Circle3, Ellipse3, Line3};
pub use surfaces::{ConeSurf, CylSurf, Plane3, SphereSurf, TorusSurf};
