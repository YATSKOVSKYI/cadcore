//! B-Rep entity structs.
//!
//! Extended with `PartialCylinder`, `PartialDisk`, and `Arc` boundary
//! for solid Boolean half-space cut support.

use cadcore_geom::{Circle3, CylSurf, Ellipse3, Line3, Plane3, SphereSurf, TorusSurf};
use cadcore_math::{Point3, UnitVec3};

use crate::ids::*;

// ── FaceExtent ────────────────────────────────────────────────────────────────

/// Geometric extent of a face — data *not* captured by the infinite carrier surface.
#[derive(Clone, Debug)]
pub enum FaceExtent {
    /// No extent info available (placeholder — face will have empty bounds).
    None,
    /// Finite cylinder extending from z = 0 to z = `length` along the carrier
    /// cylinder's axis.  Full circular cross-section.
    Cylinder {
        /// Length along the cylinder axis (mm).
        length: f64,
        /// Boundary at z = 0.
        start: FaceBoundary,
        /// Boundary at z = `length`.
        end: FaceBoundary,
    },
    /// **NEW** — Finite cylinder with an arc cross-section (chord cut by a plane
    /// parallel to the cylinder axis).
    ///
    /// The surviving arc spans from `arc_start_angle` to `arc_end_angle`
    /// (CCW, in radians) measured in the cylinder cross-section plane.
    /// Angle 0 points in the `arc_ref_dir` direction.
    PartialCylinder {
        /// Length along the cylinder axis (mm).
        length: f64,
        /// Arc start angle (radians, CCW from `arc_ref_dir`).
        arc_start_angle: f64,
        /// Arc end angle (radians, CCW from `arc_ref_dir`).
        arc_end_angle: f64,
        /// Reference direction in the cross-section (usually the cut plane normal).
        arc_ref_dir: UnitVec3,
    },
    /// Planar disk at the carrier plane's origin with the given radius.
    Disk {
        /// Disk radius (mm).
        radius: f64,
    },
    /// **NEW** — Planar arc sector (partial disk): the region bounded by an arc
    /// from `start_angle` to `end_angle` and the chord connecting the endpoints.
    ///
    /// Used for the end caps of a `PartialCylinder` solid.
    PartialDisk {
        /// Outer arc radius (mm).
        radius: f64,
        /// Arc start angle (radians, CCW from the plane's x-axis).
        start_angle: f64,
        /// Arc end angle (radians).
        end_angle: f64,
    },
    /// Torus fillet arc: the boundary consists of two minor circles.
    TorusFillet {
        /// Minor circle bounding the arc at the incoming-cylinder junction.
        start_circle: Circle3,
        /// Minor circle bounding the arc at the outgoing-cylinder junction.
        end_circle: Circle3,
    },
    /// Planar boundary described by a FaceBoundary (Circle or Ellipse).
    PlanarBoundary {
        /// The boundary curve.
        boundary: FaceBoundary,
    },
    /// Flat polygonal face defined by its vertices (useful for box/plate solids).
    Polygon {
        /// The polygon vertices in counter-clockwise order.
        points: Vec<Point3>,
    },
}

/// Closed curve used as a face boundary.
#[derive(Clone, Debug)]
pub enum FaceBoundary {
    /// Full circle boundary.
    Circle(Circle3),
    /// Full ellipse boundary, typically a cylinder/miter-plane intersection.
    Ellipse(Ellipse3),
}

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

/// An undirected topological edge.
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
    /// Back-reference to paired `CoEdge`.
    pub partner: Option<CoEdgeId>,
}

// ── CoEdge ───────────────────────────────────────────────────────────────────

/// Whether a co-edge traverses its underlying [`Edge`] in the same or opposite direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoEdgeSense {
    /// Same direction as the edge.
    Same,
    /// Opposite direction.
    Opposite,
}

/// A *directed* use of an [`Edge`] by a single [`Loop`].
#[derive(Clone, Debug)]
pub struct CoEdge {
    /// The underlying geometric edge.
    pub edge:    EdgeId,
    /// Traversal sense.
    pub sense:   CoEdgeSense,
    /// Next co-edge in the loop.
    pub next:    CoEdgeId,
    /// Previous co-edge in the loop.
    pub prev:    CoEdgeId,
    /// The loop this co-edge belongs to.
    pub loop_id: LoopId,
}

// ── Loop ─────────────────────────────────────────────────────────────────────

/// A closed chain of co-edges bounding (part of) a face.
#[derive(Clone, Debug)]
pub struct Loop {
    /// Arbitrary start co-edge.
    pub start: CoEdgeId,
    /// The face this loop belongs to.
    pub face:  FaceId,
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
    pub geom:        FaceGeom,
    /// Normal orientation relative to the surface.
    pub normal:      FaceNormal,
    /// Outer boundary loop.
    pub outer_loop:  LoopId,
    /// Inner loops (holes / voids within the face).
    pub inner_loops: Vec<LoopId>,
    /// The shell this face belongs to.
    pub shell:       ShellId,
    /// Geometric extent — used by the STEP writer for face bounds.
    pub extent:      FaceExtent,
}

// ── Shell ─────────────────────────────────────────────────────────────────────

/// A connected, orientable set of faces forming a closed or open surface.
#[derive(Clone, Debug)]
pub struct Shell {
    /// All faces in the shell.
    pub faces:    Vec<FaceId>,
    /// Outer boundary (`true`) or inner void (`false`) of its solid.
    pub is_outer: bool,
    /// The solid this shell belongs to.
    pub solid:    SolidId,
}

// ── Solid ─────────────────────────────────────────────────────────────────────

/// A single connected solid body defined by its boundary shells.
#[derive(Clone, Debug)]
pub struct Solid {
    /// All shells (outer first, then voids).
    pub shells: Vec<ShellId>,
    /// Optional human-readable name (used in STEP output).
    pub name:   Option<String>,
}
