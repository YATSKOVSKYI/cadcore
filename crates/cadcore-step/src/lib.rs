//! # cadcore-step
//!
//! Pure-Rust STEP AP203 writer.
//!
//! Converts a [`cadcore_topo::BRep`] into a compliant ISO 10303-21
//! exchange file (`.step` / `.stp`).
//!
//! ## STEP entity mapping
//!
//! | cadcore type       | STEP entity                   |
//! |--------------------|-------------------------------|
//! | `Plane3`           | `PLANE`                       |
//! | `CylSurf`          | `CYLINDRICAL_SURFACE`         |
//! | `SphereSurf`       | `SPHERICAL_SURFACE`           |
//! | `TorusSurf`        | `TOROIDAL_SURFACE`            |
//! | `Line3`            | `LINE`                        |
//! | `Circle3`          | `CIRCLE`                      |
//! | `Ellipse3`         | `ELLIPSE`                     |
//! | `Point3`           | `CARTESIAN_POINT`             |
//! | `UnitVec3`/`Vec3`  | `DIRECTION` / `VECTOR`        |
//! | `Frame3`           | `AXIS2_PLACEMENT_3D`          |

#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

mod writer;
mod entities;

pub use writer::{brep_to_step, StepWriter};
pub use entities::StepError;
