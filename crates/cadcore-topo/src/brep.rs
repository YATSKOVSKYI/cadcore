//! The [`BRep`] struct — the main container holding all B-Rep entities.

use slotmap::SlotMap;

use crate::{entities::*, ids::*};

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
    pub vertices: SlotMap<VertexId, Vertex>,
    /// Edge arena.
    pub edges: SlotMap<EdgeId, Edge>,
    /// CoEdge arena.
    pub coedges: SlotMap<CoEdgeId, CoEdge>,
    /// Loop arena.
    pub loops: SlotMap<LoopId, Loop>,
    /// Face arena.
    pub faces: SlotMap<FaceId, Face>,
    /// Shell arena.
    pub shells: SlotMap<ShellId, Shell>,
    /// Solid arena.
    pub solids: SlotMap<SolidId, Solid>,
}

impl BRep {
    /// Create an empty B-Rep.
    pub fn new() -> Self {
        Self::default()
    }

    // ── Insertion helpers ────────────────────────────────────────────────────

    /// Add a vertex and return its id.
    #[inline]
    pub fn add_vertex(&mut self, v: Vertex) -> VertexId {
        self.vertices.insert(v)
    }

    /// Add an edge and return its id.
    #[inline]
    pub fn add_edge(&mut self, e: Edge) -> EdgeId {
        self.edges.insert(e)
    }

    /// Add a co-edge and return its id.
    ///
    /// **Note:** the caller is responsible for wiring up `next`/`prev` after
    /// all co-edges in the loop are inserted.  Use a placeholder id and then
    /// [`BRep::patch_coedge_links`] when the loop is complete.
    #[inline]
    pub fn add_coedge(&mut self, ce: CoEdge) -> CoEdgeId {
        self.coedges.insert(ce)
    }

    /// Add a loop and return its id.
    #[inline]
    pub fn add_loop(&mut self, lp: Loop) -> LoopId {
        self.loops.insert(lp)
    }

    /// Add a face and return its id.
    #[inline]
    pub fn add_face(&mut self, f: Face) -> FaceId {
        self.faces.insert(f)
    }

    /// Add a shell and return its id.
    #[inline]
    pub fn add_shell(&mut self, s: Shell) -> ShellId {
        self.shells.insert(s)
    }

    /// Add a solid and return its id.
    #[inline]
    pub fn add_solid(&mut self, s: Solid) -> SolidId {
        self.solids.insert(s)
    }

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
        let mut cur = start;
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
            edges: self.edges.len(),
            coedges: self.coedges.len(),
            loops: self.loops.len(),
            faces: self.faces.len(),
            shells: self.shells.len(),
            solids: self.solids.len(),
        }
    }

    /// Deep-copy a single solid from `src` into `self`, remapping every entity
    /// id, and return the new solid's id (`None` if `solid_id` is invalid).
    ///
    /// Used to graft an independent body (e.g. a named electrode) into the
    /// result of a Boolean union that was computed in its own fresh `BRep`.
    /// The whole reachable sub-graph — vertices → edges → co-edges → loops →
    /// faces → shells → solid — is cloned, then every internal id reference
    /// (edge endpoints, co-edge `next`/`prev`/`edge`/`loop_id`, loop
    /// `start`/`face`, face loops/shell, shell faces/solid) is rewritten to the
    /// freshly-minted ids.  Geometry is value-copied, so the two BReps stay
    /// independent.
    pub fn append_solid(&mut self, src: &BRep, solid_id: SolidId) -> Option<SolidId> {
        use std::collections::HashMap;
        let solid = src.solids.get(solid_id)?;

        // ── Gather the reachable sub-graph ────────────────────────────────────
        let mut edge_set: Vec<EdgeId> = Vec::new();
        let mut coedge_set: Vec<CoEdgeId> = Vec::new();
        let mut loop_set: Vec<LoopId> = Vec::new();
        let mut face_set: Vec<FaceId> = Vec::new();
        let mut vert_set: Vec<VertexId> = Vec::new();

        macro_rules! push_unique {
            ($v:expr, $id:expr) => {{
                let id = $id;
                if !$v.contains(&id) {
                    $v.push(id);
                }
            }};
        }

        for &sh in &solid.shells {
            let Some(shell) = src.shells.get(sh) else { continue };
            for &fid in &shell.faces {
                push_unique!(face_set, fid);
                let Some(face) = src.faces.get(fid) else { continue };
                let mut loops = vec![face.outer_loop];
                loops.extend(face.inner_loops.iter().copied());
                for lp in loops {
                    if src.loops.get(lp).is_none() {
                        continue;
                    }
                    push_unique!(loop_set, lp);
                    if let Some(coedges) = src.loop_coedges(lp) {
                        for ce_id in coedges {
                            push_unique!(coedge_set, ce_id);
                            if let Some(ce) = src.coedges.get(ce_id) {
                                push_unique!(edge_set, ce.edge);
                                if let Some(edge) = src.edges.get(ce.edge) {
                                    push_unique!(vert_set, edge.v_start);
                                    push_unique!(vert_set, edge.v_end);
                                }
                            }
                        }
                    }
                }
            }
        }

        // ── Clone every entity, recording old → new id ────────────────────────
        let mut vmap: HashMap<VertexId, VertexId> = HashMap::new();
        for &v in &vert_set {
            if let Some(vert) = src.vertices.get(v) {
                vmap.insert(v, self.add_vertex(vert.clone()));
            }
        }
        let mut emap: HashMap<EdgeId, EdgeId> = HashMap::new();
        for &e in &edge_set {
            if let Some(edge) = src.edges.get(e) {
                emap.insert(e, self.add_edge(edge.clone()));
            }
        }
        let mut cmap: HashMap<CoEdgeId, CoEdgeId> = HashMap::new();
        for &c in &coedge_set {
            if let Some(ce) = src.coedges.get(c) {
                cmap.insert(c, self.add_coedge(ce.clone()));
            }
        }
        let mut lmap: HashMap<LoopId, LoopId> = HashMap::new();
        for &l in &loop_set {
            if let Some(lp) = src.loops.get(l) {
                lmap.insert(l, self.add_loop(lp.clone()));
            }
        }
        let mut fmap: HashMap<FaceId, FaceId> = HashMap::new();
        for &f in &face_set {
            if let Some(face) = src.faces.get(f) {
                fmap.insert(f, self.add_face(face.clone()));
            }
        }
        let mut shmap: HashMap<ShellId, ShellId> = HashMap::new();
        for &sh in &solid.shells {
            if let Some(shell) = src.shells.get(sh) {
                shmap.insert(sh, self.add_shell(shell.clone()));
            }
        }
        let new_solid = self.add_solid(solid.clone());

        // ── Rewrite every internal id reference to the new ids ────────────────
        let remap_v = |id| vmap.get(&id).copied().unwrap_or(id);
        let remap_e = |id| emap.get(&id).copied().unwrap_or(id);
        let remap_c = |id| cmap.get(&id).copied().unwrap_or(id);
        let remap_l = |id| lmap.get(&id).copied().unwrap_or(id);
        let remap_f = |id| fmap.get(&id).copied().unwrap_or(id);
        let remap_sh = |id| shmap.get(&id).copied().unwrap_or(id);

        for &new in emap.values() {
            if let Some(edge) = self.edges.get_mut(new) {
                edge.v_start = remap_v(edge.v_start);
                edge.v_end = remap_v(edge.v_end);
                edge.partner = edge.partner.and_then(|p| cmap.get(&p).copied());
            }
        }
        for &new in cmap.values() {
            if let Some(ce) = self.coedges.get_mut(new) {
                ce.edge = remap_e(ce.edge);
                ce.next = remap_c(ce.next);
                ce.prev = remap_c(ce.prev);
                ce.loop_id = remap_l(ce.loop_id);
            }
        }
        for &new in lmap.values() {
            if let Some(lp) = self.loops.get_mut(new) {
                lp.start = remap_c(lp.start);
                lp.face = remap_f(lp.face);
            }
        }
        for &new in fmap.values() {
            if let Some(face) = self.faces.get_mut(new) {
                face.outer_loop = remap_l(face.outer_loop);
                for il in &mut face.inner_loops {
                    *il = remap_l(*il);
                }
                face.shell = remap_sh(face.shell);
            }
        }
        for &new in shmap.values() {
            if let Some(shell) = self.shells.get_mut(new) {
                for f in &mut shell.faces {
                    *f = remap_f(*f);
                }
                shell.solid = new_solid;
            }
        }
        if let Some(s) = self.solids.get_mut(new_solid) {
            for sh in &mut s.shells {
                *sh = remap_sh(*sh);
            }
        }

        Some(new_solid)
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
    pub edges: usize,
    /// Number of co-edges.
    pub coedges: usize,
    /// Number of loops.
    pub loops: usize,
    /// Number of faces.
    pub faces: usize,
    /// Number of shells.
    pub shells: usize,
    /// Number of solids.
    pub solids: usize,
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
        b.add_vertex(Vertex {
            point: Point3::ORIGIN,
        });
        assert_eq!(b.stats().vertices, 1);
    }

    /// Build a degenerate-but-well-linked one-face solid (triangle loop) for
    /// exercising topology copies.
    fn triangle_solid() -> BRep {
        use cadcore_geom::{Line3, Plane3};
        use cadcore_math::{UnitVec3, Vec3};
        let mut b = BRep::new();
        let pts = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let v: Vec<_> = pts
            .iter()
            .map(|&p| b.add_vertex(Vertex { point: p }))
            .collect();
        let face = b.add_face(Face {
            geom: FaceGeom::Plane(Plane3::from_origin_normal(
                Point3::ORIGIN,
                UnitVec3::new_unchecked(Vec3::Z),
            )),
            normal: FaceNormal::Same,
            outer_loop: LoopId::default(),
            inner_loops: Vec::new(),
            shell: ShellId::default(),
            extent: FaceExtent::None,
        });
        let lp = b.add_loop(Loop {
            start: CoEdgeId::default(),
            face,
        });
        let mut coedges = Vec::new();
        for i in 0..3 {
            let a = v[i];
            let c = v[(i + 1) % 3];
            let e = b.add_edge(Edge {
                geom: EdgeGeom::Line(Line3::new(pts[i], UnitVec3::new_unchecked(Vec3::X))),
                v_start: a,
                v_end: c,
                t_start: 0.0,
                t_end: 1.0,
                partner: None,
            });
            coedges.push(b.add_coedge(CoEdge {
                edge: e,
                sense: CoEdgeSense::Same,
                next: CoEdgeId::default(),
                prev: CoEdgeId::default(),
                loop_id: lp,
            }));
        }
        b.patch_coedge_links(&coedges);
        b.loops[lp].start = coedges[0];
        b.faces[face].outer_loop = lp;
        let shell = b.add_shell(Shell {
            faces: vec![face],
            is_outer: true,
            solid: SolidId::default(),
        });
        b.faces[face].shell = shell;
        let solid = b.add_solid(Solid {
            shells: vec![shell],
            name: Some("electrode_test".into()),
        });
        b.shells[shell].solid = solid;
        b
    }

    #[test]
    fn append_solid_deep_copies_and_remaps() {
        let src = triangle_solid();
        let mut dst = BRep::new();
        // pre-populate dst so the new ids cannot accidentally equal src ids.
        dst.add_vertex(Vertex { point: Point3::ORIGIN });

        let src_solid = src.solids.keys().next().unwrap();
        let new_solid = dst.append_solid(&src, src_solid).expect("append");

        // Counts: dst gained the whole sub-graph.
        assert_eq!(dst.stats().faces, 1);
        assert_eq!(dst.stats().edges, 3);
        assert_eq!(dst.stats().coedges, 3);
        assert_eq!(dst.solids.len(), 1);

        // The copied solid's loop walks 3 co-edges, each referencing a copied
        // edge that lives in dst (not src).
        let solid = dst.solids.get(new_solid).unwrap();
        let shell = dst.shells.get(solid.shells[0]).unwrap();
        assert_eq!(shell.solid, new_solid);
        let face = dst.faces.get(shell.faces[0]).unwrap();
        let coedges = dst.loop_coedges(face.outer_loop).expect("walks");
        assert_eq!(coedges.len(), 3);
        for ce_id in coedges {
            let ce = dst.coedges.get(ce_id).unwrap();
            assert!(dst.edges.contains_key(ce.edge), "edge remapped into dst");
            assert_eq!(ce.loop_id, face.outer_loop);
        }
        assert_eq!(solid.name.as_deref(), Some("electrode_test"));
    }
}
