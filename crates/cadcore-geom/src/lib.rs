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

pub mod arrangement;
pub mod curves;
pub mod intersect;
pub mod ssi;
pub mod surfaces;
pub mod swept;
pub mod trim;

pub use arrangement::{Arrangement, ArrangementOptions, Cell, ChainInput, LoopStep};
pub use curves::{BezierCubic, Circle3, Ellipse3, Line3};
pub use intersect::{
    cyl_cyl_intersection, cyl_cyl_parallel_lines, cyl_plane_intersection, cyl_surface_distance,
    sample_cyl_plane, torus_cyl_intersection, CylPlaneCurve, IntersectionPolyline,
};
pub use ssi::{surface_surface_intersection, swept_swept_intersection, SsiOptions};
pub use surfaces::{ConeSurf, CylSurf, Plane3, SphereSurf, TorusSurf};
pub use swept::{
    CenterlineSeg, ParamSurface, SurfaceDerivatives, SweptProjection, SweptTubeSurface,
};
pub use trim::{
    loop_interior_sample, loop_signed_area_uv, orient_loop_as_hole, point_inside_tube,
    project_loop_to_uv,
};
