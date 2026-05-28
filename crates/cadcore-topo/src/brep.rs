//! The [`BRep`] struct — the main container holding all B-Rep entities.

use slotmap::SlotMap;

use crate::{
    ids::*,
    entities::*,
};

/// The complete B-Rep model.
///
/// All entities are stored in typed arena maps (`SlotMap`).  Ids remain
/// valid as long as the corresponding entity has not been explicitly
/// removed.
///
/// # Typical usage
///
/// ```rust
/// use cadcore_topo::BRep;
/// let mut brep = BRep::new();
/// // ... use cadcore_ops to populate it ...
/// ```
#[derive(Debug, Default)]
pub struct BRep {
    /// Vertex arena.
    pub vertices:  SlotMap<VertexId,  Vertex>,
    /// Edge arena.
    pub edges:     SlotMap<EdgeId,    Edge>,
    /// CoEdge arena.
    pub coedges:   SlotMap<CoEdgeId,  CoEdge>,
    /// Loop arena.
    pub loops:     SlotMap<LoopId,    Loop>,
    /// Face arena.
    pub faces:     SlotMap<FaceId,    Face>,
    /// Shell arena.
    pub shells:    SlotMap<ShellId,   Shell>,
    /// Solid arena.
    pub solids:    SlotMap<SolidId,   Solid>,
}

impl BRep {
    /// Create an empty B-Rep.
    pub fn new() -> Self { Self::default() }

    // ── Insertion helpers ────────────────────────────────────────────────────

    /// Add a vertex and return its id.
    #[inline]
    pub fn add_vertex(&mut self, v: Vertex) -> VertexId { self.vertices.insert(v) }

    /// Add an edge and return its id.
    #[inline]
    pub fn add_edge(&mut self, e: Edge) -> EdgeId { self.edges.insert(e) }

    /// Add a co-edge and return its id.
    ///
    /// **Note:** the caller is responsible for wiring up `next`/`prev` after
    /// all co-edges in the loop are inserted.  Use a placeholder id and then
    /// [`BRep::patch_coedge_links`] when the loop is complete.
    #[inline]
    pub fn add_coedge(&mut self, ce: CoEdge) -> CoEdgeId { self.coedges.insert(ce) }

    /// Add a loop and return its id.
    #[inline]
    pub fn add_loop(&mut self, lp: Loop) -> LoopId { self.loops.insert(lp) }

    /// Add a face and return its id.
    #[inline]
    pub fn add_face(&mut self, f: Face) -> FaceId { self.faces.insert(f) }

    /// Add a shell and return its id.
    #[inline]
    pub fn add_shell(&mut self, s: Shell) -> ShellId { self.shells.insert(s) }

    /// Add a solid and return its id.
    #[inline]
    pub fn add_solid(&mut self, s: Solid) -> SolidId { self.solids.insert(s) }

    // ── Link patching ────────────────────────────────────────────────────────

    /// Fix up the `next`/`prev` links for a complete loop given the co-edges
    /// in traversal order.
    ///
    /// This is called after all co-edges have been inserted, once their ids
    /// are known.
    pub fn patch_coedge_links(&mut self, ordered: &[CoEdgeId]) {
        let n = ordered.len();
        for i in 0..n {
            let prev_id = ordered[(i + n - 1) % n];
            let next_id = ordered[(i + 1) % n];
            if let Some(ce) = self.coedges.get_mut(ordered[i]) {
                ce.prev = prev_id;
                ce.next = next_id;
            }
        }
    }

    // ── Topology queries ─────────────────────────────────────────────────────

    /// Collect the co-edges belonging to a loop in traversal order.
    ///
    /// Returns `None` if the loop id is invalid or the loop contains a broken
    /// `next` link.
    pub fn loop_coedges(&self, loop_id: LoopId) -> Option<Vec<CoEdgeId>> {
        let lp = self.loops.get(loop_id)?;
        let start = lp.start;
        let mut ids = Vec::new();
        let mut cur  = start;
        loop {
            ids.push(cur);
            let ce = self.coedges.get(cur)?;
            cur = ce.next;
            if cur == start || ids.len() > self.coedges.len() {
                break;
            }
        }
        Some(ids)
    }

    /// Count of topological entities.
    pub fn stats(&self) -> BRepStats {
        BRepStats {
            vertices: self.vertices.len(),
            edges:    self.edges.len(),
            coedges:  self.coedges.len(),
            loops:    self.loops.len(),
            faces:    self.faces.len(),
            shells:   self.shells.len(),
            solids:   self.solids.len(),
        }
    }

    /// Merge all solids in the B-Rep into a single solid.
    pub fn fuse_solids(&mut self) {
        if self.solids.is_empty() {
            return;
        }
        let mut all_shells = Vec::new();
        for (_, solid) in self.solids.drain() {
            all_shells.extend(solid.shells);
        }
        let merged_solid = Solid {
            shells: all_shells,
            name: Some("fused_solid".to_string()),
        };
        let solid_id = self.add_solid(merged_solid);
        let shells = self.solids.get(solid_id).unwrap().shells.clone();
        for shell_id in shells {
            if let Some(shell) = self.shells.get_mut(shell_id) {
                shell.solid = solid_id;
            }
        }
    }
}

/// Summary statistics about a [`BRep`].
#[derive(Debug, Clone, Copy)]
pub struct BRepStats {
    /// Number of vertices.
    pub vertices: usize,
    /// Number of edges.
    pub edges:    usize,
    /// Number of co-edges.
    pub coedges:  usize,
    /// Number of loops.
    pub loops:    usize,
    /// Number of faces.
    pub faces:    usize,
    /// Number of shells.
    pub shells:   usize,
    /// Number of solids.
    pub solids:   usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_math::Point3;

    #[test]
    fn empty_brep_has_zero_stats() {
        let b = BRep::new();
        let s = b.stats();
        assert_eq!(s.vertices, 0);
        assert_eq!(s.faces, 0);
    }

    #[test]
    fn add_vertex_increases_count() {
        let mut b = BRep::new();
        b.add_vertex(Vertex { point: Point3::ORIGIN });
        assert_eq!(b.stats().vertices, 1);
    }
}
