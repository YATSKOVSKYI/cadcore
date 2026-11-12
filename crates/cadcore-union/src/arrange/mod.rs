//! The arrangement layer — from intersection curves to authoritative face
//! loops, and the lattice fuse built on it.
//!
//! ```text
//! filament  — the filament model (legs + elbows) and a single-filament solid
//! universal — the universal fuse: any number of filaments → one closed shell
//! merge     — fold per-step edges back into whole curves
//! domain    — per-face uv parametrisations (cylinder band, torus patch,
//!             plane) and the lift back to 3-D
//! cells     — per-face DCEL arrangement of boundary wires + curve pieces,
//!             keep/drop classification of every cell
//! assembly  — kept cells → topology (welded vertices, shared edges, loops)
//! ```

pub mod assembly;
pub mod cells;
pub mod filament;
pub mod fixtures;
pub mod merge;
pub mod domain;
pub mod universal;
