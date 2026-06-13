//! The arrangement layer — from intersection curves to authoritative face
//! loops.
//!
//! This replaces the legacy per-case notch/hole/strip machinery with one
//! general algorithm:
//!
//! ```text
//! registry  — split every loop ONCE at all member-junction crossings of
//!             both solids; refined joints are the single source of truth
//!             for every consumer                                  [LIVE]
//! domain    — per-face uv parametrisations (cylinder band, torus patch,
//!             plane), lift back to 3D                             [LIVE]
//! cells     — per-face DCEL arrangement of boundary wires + curve pieces
//!             (cadcore-geom DCEL), keep/drop classification by lifting
//!             cell interior points against the other solid        [LIVE]
//! assembly  — kept cells → topo loops/coedges with derived senses [LIVE]
//! ```

pub mod assembly;
pub mod cells;
pub mod filament;
pub mod domain;
pub mod registry;
pub mod scaffold;
pub mod solid;
