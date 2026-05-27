//! Typed stable identifiers for B-Rep entities.
//!
//! Using `slotmap::new_key_type!` gives us zero-cost distinct types that
//! the compiler treats as separate — you can't accidentally pass a `FaceId`
//! where a `VertexId` is expected.

use slotmap::new_key_type;

new_key_type! {
    /// Unique identifier for a [`crate::Vertex`].
    pub struct VertexId;
}

new_key_type! {
    /// Unique identifier for an [`crate::Edge`].
    pub struct EdgeId;
}

new_key_type! {
    /// Unique identifier for a [`crate::CoEdge`].
    pub struct CoEdgeId;
}

new_key_type! {
    /// Unique identifier for a [`crate::Loop`].
    pub struct LoopId;
}

new_key_type! {
    /// Unique identifier for a [`crate::Face`].
    pub struct FaceId;
}

new_key_type! {
    /// Unique identifier for a [`crate::Shell`].
    pub struct ShellId;
}

new_key_type! {
    /// Unique identifier for a [`crate::Solid`].
    pub struct SolidId;
}
