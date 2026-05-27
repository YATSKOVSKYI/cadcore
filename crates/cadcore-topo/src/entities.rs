//! B-Rep entity structs.

use cadcore_geom::{Circle3, CylSurf, Ellipse3, Line3, Plane3, SphereSurf, TorusSurf};
use cadcore_math::Point3;

use crate::ids::*;

// ── Vertex ───────────────────────────────────────────────────────────────────

/// A point in 3-D space — the lowest-level topological entity.
#[derive(Clone, Debug)]
pub struct Vertex {
    /// World-space position (mm).
    pub point: Point3,
}

// ── Edge geometry ─────────────────────────────────────────────────────────────

/// The geometric curve underlying an [`Edge`].
#[derive(Clone, Debug)]
pub enum EdgeGeom {
    /// Straight line segment.
    Line(Line3),
    /// A circular arc.
    Circle(Circle3),
    /// An elliptic arc (e.g. miter intersection of cylinder and plane).
    Ellipse(Ellipse3),
}

/// An undirected topological edge: a bounded segment of a curve between two
/// vertices.
///
/// The *sense* of traversal (which vertex is start vs. end) is defined per
/// face-use in [`CoEdge`].
#[derive(Clone, Debug)]
pub struct Edge {
    /// Curve the edge lies on.
    pub geom: EdgeGeom,
    /// Start vertex (at parameter `t_start`).
    pub v_start: VertexId,
    /// End vertex (at parameter `t_end`).
    pub v_end: VertexId,
    /// Parameter value at the start vertex.
    pub t_start: f64,
    /// Parameter value at the end vertex.
    pub t_end: f64,
    /// Back-reference to paired `CoEdge` on the neighbouring face (may be
    /// `None` for open boundary / naked edge).
    pub partner: Option<CoEdgeId>,
}

// ── CoEdge ───────────────────────────────────────────────────────────────────

/// Whether a co-edge traverses its underlying [`Edge`] in the same or
/// opposite direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoEdgeSense {
    /// Same direction as the edge (v_start → v_end).
    Same,
    /// Opposite direction (v_end → v_start).
    Opposite,
}

/// A *directed* use of an [`Edge`] by a single [`Loop`].
///
/// Two co-edges with the same underlying edge but opposite sense are the
/// shared boundary between two adjacent faces.
#[derive(Clone, Debug)]
pub struct CoEdge {
    /// The underlying geometric edge.
    pub edge: EdgeId,
    /// Traversal sense.
    pub sense: CoEdgeSense,
    /// Next co-edge in the loop (forms a circular linked list).
    pub next: CoEdgeId,
    /// Previous co-edge in the loop.
    pub prev: CoEdgeId,
    /// The loop this co-edge belongs to.
    pub loop_id: LoopId,
}

// ── Loop ─────────────────────────────────────────────────────────────────────

/// A closed chain of co-edges bounding (part of) a face.
///
/// A face has exactly one *outer* loop and zero or more *inner* loops (holes).
#[derive(Clone, Debug)]
pub struct Loop {
    /// Arbitrary start co-edge (walk `next` to traverse).
    pub start: CoEdgeId,
    /// The face this loop belongs to.
    pub face: FaceId,
}

// ── Face geometry ─────────────────────────────────────────────────────────────

/// The carrier surface of a [`Face`].
#[derive(Clone, Debug)]
pub enum FaceGeom {
    /// Flat planar face.
    Plane(Plane3),
    /// Portion of a cylinder.
    Cylinder(CylSurf),
    /// Portion of a sphere.
    Sphere(SphereSurf),
    /// Portion of a torus.
    Torus(TorusSurf),
}

/// Whether the face's outward normal agrees with the surface's natural normal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FaceNormal {
    /// Outward normal = surface normal.
    Same,
    /// Outward normal = −surface normal.
    Reversed,
}

/// A bounded region of a surface — the main building block of a B-Rep.
#[derive(Clone, Debug)]
pub struct Face {
    /// Carrier surface.
    pub geom: FaceGeom,
    /// Normal orientation relative to the surface.
    pub normal: FaceNormal,
    /// Outer boundary loop.
    pub outer_loop: LoopId,
    /// Inner loops (holes/voids within the face).
    pub inner_loops: Vec<LoopId>,
    /// The shell this face belongs to.
    pub shell: ShellId,
}

// ── Shell ─────────────────────────────────────────────────────────────────────

/// A connected, orientable set of faces forming a closed or open surface.
///
/// A *closed* shell (no naked edges) encloses a volume.
#[derive(Clone, Debug)]
pub struct Shell {
    /// All faces in the shell.
    pub faces: Vec<FaceId>,
    /// Whether this shell is the outer boundary (`true`) or an inner void
    /// (`false`) of its solid.
    pub is_outer: bool,
    /// The solid this shell belongs to.
    pub solid: SolidId,
}

// ── Solid ─────────────────────────────────────────────────────────────────────

/// A single connected solid body defined by its boundary shells.
///
/// A valid manifold solid has exactly one outer shell and zero or more
/// inner void shells.
#[derive(Clone, Debug)]
pub struct Solid {
    /// All shells (outer first, then voids).
    pub shells: Vec<ShellId>,
    /// Optional human-readable name (used in STEP output).
    pub name: Option<String>,
}
