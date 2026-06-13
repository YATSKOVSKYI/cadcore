//! # cadcore-union — the Boolean-union engine
//!
//! Industrial-grade Boolean union for cadcore: a staged pipeline over
//! `cadcore-topo` B-Reps with an analytic intersection oracle, refinement of
//! every trim curve onto the *emitted* geometry, and in-process validation
//! gates that reproduce SolidWorks' accept/reject behaviour.
//!
//! See `README.md` for the architecture and the hard-won invariants.
//!
//! ## Status
//!
//! P0: tolerance model, structured diagnostics, the refinement oracle
//! (`geom::refine`) and validation gates 1 (manifold) and 4 (edge-to-surface
//! distance) are live.  The legacy union in `cadcore-ops` remains the
//! production path; use [`validate::validate_shell`] on its output today.

pub mod arrange;
pub mod boolean;
pub mod config;
pub mod diag;
pub mod geom;
pub mod pipeline;
pub mod tolerance;
pub mod validate;

pub use config::UnionConfig;
pub use diag::{Diagnostics, Event, Severity};
pub use pipeline::{UnionOutcome, UnionPipeline};
pub use tolerance::Tolerances;
pub use validate::{validate_shell, ValidationReport};
