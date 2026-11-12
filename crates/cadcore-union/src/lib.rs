//! # cadcore-union — Booleans for the cadcore kernel
//!
//! Fuse swept filaments into one watertight solid, trim it with planes, and
//! ask the same solid whether a point is material — in pure Rust, on exact
//! analytic surfaces, with in-process validation gates that reproduce what a
//! CAD importer accepts.
//!
//! See `README.md` for how it works, what it was verified against and its
//! limits.
//!
//! ## The API in one screen
//!
//! | Call | What it does |
//! |---|---|
//! | [`fuse_and_trim`] | filaments → one fused solid → planar trims → merged edges |
//! | [`fuse_filaments_with`] | universal fuse of swept filaments into one closed shell |
//! | [`trim_shell`] | robust planar trims of an already fused shell + edge merge |
//! | [`LatticeModel`] | exact point membership of the same solid, for volume meshers |
//! | [`half_space_cut_robust`] / [`half_space_cut`] | one planar cut with section caps |
//! | [`merge_shell_edges`] | fold per-step edges back into whole curves |
//! | [`boolean::union_n`] | union of analytic solids (boxes, cylinders, filaments, coaxial tori) |
//! | [`validate_shell`] | the manifold / wire / distance gates |
//!
//! Input is a list of [`Filament`]s ([`LegSeg`] cylinders joined by
//! [`ElbowSeg`] tori; [`Filament::serpentine`] builds one from a polyline).

pub mod arrange;
pub mod boolean;
pub mod config;
pub mod diag;
pub(crate) mod fasthash;
pub mod geom;
pub mod lattice;
pub mod tolerance;
pub mod validate;

pub use arrange::filament::{ElbowSeg, Filament, LegSeg};
pub use arrange::merge::merge_shell_edges;
pub use arrange::universal::{fuse_filaments, fuse_filaments_with, FuseOptions, FuseReport};
pub use boolean::{half_space_cut, half_space_cut_robust, CutOutcome, RobustCut, ROBUST_CUT_SHIFTS};
pub use config::UnionConfig;
pub use lattice::{fuse_and_trim, trim_shell, LatticeModel, LatticeOutcome, TrimInfo, TrimOutcome};
pub use diag::{Diagnostics, Event, Severity};
pub use tolerance::Tolerances;
pub use validate::{validate_shell, ValidationReport};
