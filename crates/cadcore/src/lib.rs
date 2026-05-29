//! # cadcore
//!
//! Pure-Rust CAD geometry kernel.
//!
//! This is the **facade crate** — it re-exports the public API of every
//! sub-crate so downstream users only need a single `cadcore` dependency.
//!
//! ## Quick start
//!
//! ```rust,no_run
//! use cadcore::{
//!     math::Point3,
//!     topo::BRep,
//!     ops::{sweep_circle_along_polyline, SweepOptions},
//!     step::brep_to_step,
//! };
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let waypoints = vec![
//!         Point3::new(0.0, 0.0, 0.0),
//!         Point3::new(0.0, 0.0, 10.0),
//!     ];
//!     let mut brep = BRep::new();
//!     sweep_circle_along_polyline(&mut brep, &waypoints, 0.5, &SweepOptions::default())?;
//!     let step_text = brep_to_step(&brep)?;
//!     std::fs::write("rod.step", &step_text)?;
//!     Ok(())
//! }
//! ```
//!
//! ## Crate layout
//!
//! | Module         | Crate                | Contents                                  |
//! |----------------|----------------------|-------------------------------------------|
//! | [`math`]       | `cadcore-math`       | Point3, Vec3, UnitVec3, Mat3, Frame3, …   |
//! | [`geom`]       | `cadcore-geom`       | Plane3, CylSurf, TorusSurf, Circle3, …    |
//! | [`topo`]       | `cadcore-topo`       | BRep, Solid, Shell, Face, Edge, Vertex    |
//! | [`ops`]        | `cadcore-ops`        | sweep_circle_along_polyline, …            |
//! | [`step`]       | `cadcore-step`       | brep_to_step, StepWriter                  |

#![warn(missing_docs)]
#![doc(html_root_url = "https://docs.rs/cadcore/0.1.11")]

/// Re-export of `cadcore-math`.
pub mod math {
    pub use cadcore_math::*;
}

/// Re-export of `cadcore-geom`.
pub mod geom {
    pub use cadcore_geom::*;
}

/// Re-export of `cadcore-topo`.
pub mod topo {
    pub use cadcore_topo::*;
}

/// Re-export of `cadcore-ops`.
pub mod ops {
    pub use cadcore_ops::*;
}

/// Re-export of `cadcore-step`.
pub mod step {
    pub use cadcore_step::*;
}

// Flat convenience re-exports for the most-used types.

pub use cadcore_math::{Frame3, Point3, UnitVec3, Vec3};
pub use cadcore_ops::{
    analytic_path_from_polyline_samples, rounded_path_from_polyline,
    sharp_path_from_polyline_samples, sweep_circle_along_path, sweep_circle_along_polyline,
    sweep_circle_along_rounded_polyline, PathApproxOptions, SweepOptions, SweepPathSegment,
    ClipPlane, clip_polyline, clip_polyline_with_radius, sweep_circle_along_path_with_caps,
    sweep_circle_along_polyline_with_caps, build_solid_box,
};
pub use cadcore_step::brep_to_step;
pub use cadcore_topo::BRep;
