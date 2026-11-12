//! # cadcore-union — the Boolean-union engine
//!
//! Industrial-grade Boolean union for cadcore: a staged pipeline over
//! `cadcore-topo` B-Reps with an analytic intersection oracle, refinement of
//! every trim curve onto the *emitted* geometry, and in-process validation
//! gates that reproduce SolidWorks' accept/reject behaviour.
//!
//! See `README.md` for the architecture and the hard-won invariants.
//!
//! ## The API in one screen
//!
//! | Call | What it does |
//! |---|---|
//! | [`fuse_and_trim`] | filaments → one fused solid → planar trims → merged edges |
//! | [`fuse_filaments_with`] | universal fuse of swept filaments into one closed shell |
//! | [`trim_shell`] | robust planar trims of an already fused shell + edge merge |
//! | [`half_space_cut_robust`] / [`half_space_cut`] | one planar cut with section caps |
//! | [`merge_shell_edges`] | fold per-step edges back into whole curves |
//! | [`validate_shell`] | the manifold / wire / distance gates |
//!
//! Input is a list of [`Filament`]s ([`LegSeg`] cylinders joined by
//! [`ElbowSeg`] tori; [`Filament::serpentine`] builds one from a polyline).

pub mod arrange;
pub mod boolean;
pub mod config;
pub mod diag;
pub mod geom;
pub mod lattice;
pub mod pipeline;
pub mod tolerance;
pub mod validate;

pub use arrange::filament::{ElbowSeg, LegSeg};
pub use arrange::merge::merge_shell_edges;
pub use arrange::scaffold::Filament;
pub use arrange::universal::{fuse_filaments, fuse_filaments_with, FuseOptions, FuseReport};
pub use boolean::{half_space_cut, half_space_cut_robust, CutOutcome, RobustCut, ROBUST_CUT_SHIFTS};
pub use config::UnionConfig;
pub use lattice::{fuse_and_trim, trim_shell, LatticeOutcome, TrimInfo, TrimOutcome};
pub use diag::{Diagnostics, Event, Severity};
pub use pipeline::{UnionOutcome, UnionPipeline};
pub use tolerance::Tolerances;
pub use validate::{validate_shell, ValidationReport};
