//! # cadcore-topo
//!
//! Arena-based B-Rep (Boundary Representation) topology.
//!
//! ## Data model
//!
//! ```text
//! Solid
//!  └─ Shell (one outer + optional inner voids)
//!      └─ Face (a bounded region of a surface)
//!          └─ Loop (a closed chain of co-edges — outer loop + optional holes)
//!              └─ CoEdge (directed use of an Edge by a Face)
//!                  └─ Edge (geometric curve + two vertex endpoints)
//!                      └─ Vertex (a 3-D point)
//! ```
//!
//! All entities are stored in typed arenas (`slotmap::SlotMap`).
//! Cross-references are typed IDs, not raw pointers, so the borrow
//! checker cannot cause use-after-free and IDs remain stable across
//! insertions.

#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

mod brep;
mod entities;
mod ids;

pub use brep::BRep;
pub use entities::{
    CoEdge, CoEdgeSense, Edge, EdgeGeom, Face, FaceBoundary, FaceExtent, FaceGeom, FaceNormal,
    Loop, Shell, Solid, Vertex,
};
pub use ids::{CoEdgeId, EdgeId, FaceId, LoopId, ShellId, SolidId, VertexId};
