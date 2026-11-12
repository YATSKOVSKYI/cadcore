//! **Edge merging** — collapse chains of consecutive edges that separate the
//! same two faces into single edges.
//!
//! The arrangement emits one edge per 2-point step of every chain, so a fused
//! or cut shell carries thousands of hair-short edges: a junction circle is 48
//! arcs, a crossing window 60 polyline segments.  Every one of them becomes an
//! `EDGE_CURVE` with its own pcurves in STEP — a multi-megabyte file that
//! receiving CADs chew on for seconds and, worse, a wire of dozens of tiny
//! B-splines where one smooth curve belongs.
//!
//! Merging is purely topological: two consecutive edges of a loop are fused
//! when their shared vertex belongs to exactly those two edges, the partner
//! face sees the same two edges consecutively (in the opposite order, as a
//! consistently oriented shell has it) and their geometry continues — two
//! polylines, two arcs of one circle, two collinear lines.  Vertices where a
//! third edge arrives (joints, seam touches, triple points) are kept, so the
//! shell's topology stays exactly what the Boolean computed.

use std::collections::HashMap;

use cadcore_geom::Line3;
use cadcore_math::Point3;
use cadcore_topo::{BRep, CoEdgeId, CoEdgeSense, Edge, EdgeGeom, EdgeId, LoopId, ShellId, VertexId};

/// Merge every mergeable edge pair of `shell`.  Returns how many edges were
/// removed.
pub fn merge_shell_edges(brep: &mut BRep, shell: ShellId) -> usize {
    let Some(sh) = brep.shells.get(shell) else { return 0 };
    let mut loops: Vec<LoopId> = Vec::new();
    for &fid in &sh.faces {
        if let Some(f) = brep.faces.get(fid) {
            loops.push(f.outer_loop);
            loops.extend(f.inner_loops.iter().copied());
        }
    }
    let mut removed = 0usize;
    loop {
        let (edge_coedges, vertex_edges) = incidence(brep, &loops);
        let mut did = 0usize;
        for &lid in &loops {
            // Walk the ring; after a merge the ring changed, so restart it.
            'ring: loop {
                let Some(lp) = brep.loops.get(lid) else { break };
                let start = lp.start;
                if brep.coedges.get(start).is_none() {
                    break;
                }
                let mut c = start;
                let mut guard = 0usize;
                loop {
                    let n = brep.coedges[c].next;
                    if n == c || n == start && guard == 0 && brep.coedges[start].next == start {
                        break 'ring;
                    }
                    if try_merge(brep, lid, c, n, &edge_coedges, &vertex_edges) {
                        did += 1;
                        // incidence is stale now — rebuild before continuing
                        break 'ring;
                    }
                    c = n;
                    guard += 1;
                    if c == start || guard > 200_000 {
                        break 'ring;
                    }
                }
            }
        }
        if did == 0 {
            break;
        }
        removed += did;
    }
    removed
}

/// `edge → its co-edges` and `vertex → edges meeting there` over `loops`.
fn incidence(
    brep: &BRep,
    loops: &[LoopId],
) -> (HashMap<EdgeId, Vec<CoEdgeId>>, HashMap<VertexId, Vec<EdgeId>>) {
    let mut edge_coedges: HashMap<EdgeId, Vec<CoEdgeId>> = HashMap::new();
    for &lid in loops {
        let Some(lp) = brep.loops.get(lid) else { continue };
        let start = lp.start;
        if brep.coedges.get(start).is_none() {
            continue;
        }
        let mut c = start;
        let mut guard = 0usize;
        loop {
            edge_coedges.entry(brep.coedges[c].edge).or_default().push(c);
            c = brep.coedges[c].next;
            guard += 1;
            if c == start || guard > 200_000 {
                break;
            }
        }
    }
    let mut vertex_edges: HashMap<VertexId, Vec<EdgeId>> = HashMap::new();
    for &eid in edge_coedges.keys() {
        let e = &brep.edges[eid];
        for v in [e.v_start, e.v_end] {
            let list = vertex_edges.entry(v).or_default();
            if !list.contains(&eid) {
                list.push(eid);
            }
        }
    }
    (edge_coedges, vertex_edges)
}

fn from_to(brep: &BRep, c: CoEdgeId) -> (VertexId, VertexId) {
    let ce = &brep.coedges[c];
    let e = &brep.edges[ce.edge];
    match ce.sense {
        CoEdgeSense::Same => (e.v_start, e.v_end),
        CoEdgeSense::Opposite => (e.v_end, e.v_start),
    }
}

/// The edge's geometry as a polyline in the co-edge's traversal direction.
fn traversal_points(brep: &BRep, c: CoEdgeId) -> Vec<Point3> {
    let ce = &brep.coedges[c];
    let e = &brep.edges[ce.edge];
    let (vf, vt) = from_to(brep, c);
    let (pf, pt) = (brep.vertices[vf].point, brep.vertices[vt].point);
    match &e.geom {
        EdgeGeom::Polyline(p) => {
            let mut pts = p.clone();
            if ce.sense == CoEdgeSense::Opposite {
                pts.reverse();
            }
            // the stored polyline runs v_start → v_end; make the ends exact
            if let Some(first) = pts.first_mut() {
                *first = pf;
            }
            if let Some(last) = pts.last_mut() {
                *last = pt;
            }
            pts
        }
        _ => vec![pf, pt],
    }
}

/// Signed parameter sweep of an edge in the co-edge's direction.
fn sweep(brep: &BRep, c: CoEdgeId) -> f64 {
    let ce = &brep.coedges[c];
    let e = &brep.edges[ce.edge];
    let d = e.t_end - e.t_start;
    if ce.sense == CoEdgeSense::Same {
        d
    } else {
        -d
    }
}

/// Try to fuse the edge of `c1` with the edge of `c2` (= `c1.next`).
fn try_merge(
    brep: &mut BRep,
    lid: LoopId,
    c1: CoEdgeId,
    c2: CoEdgeId,
    edge_coedges: &HashMap<EdgeId, Vec<CoEdgeId>>,
    vertex_edges: &HashMap<VertexId, Vec<EdgeId>>,
) -> bool {
    let e1 = brep.coedges[c1].edge;
    let e2 = brep.coedges[c2].edge;
    if e1 == e2 {
        return false; // a seam slit traversed twice in a row — not a chain
    }
    let (a, v) = from_to(brep, c1);
    let (v2, b) = from_to(brep, c2);
    if v != v2 || a == b {
        return false; // broken ring, or the merge would close the edge on itself
    }
    // exactly these two edges meet at v
    match vertex_edges.get(&v) {
        Some(list) if list.len() == 2 && list.contains(&e1) && list.contains(&e2) => {}
        _ => return false,
    }
    // partner co-edges, consecutive in the opposite order on the other face
    let (Some(u1), Some(u2)) = (edge_coedges.get(&e1), edge_coedges.get(&e2)) else {
        return false;
    };
    if u1.len() != 2 || u2.len() != 2 {
        return false;
    }
    let p1 = if u1[0] == c1 { u1[1] } else { u1[0] };
    let p2 = if u2[0] == c2 { u2[1] } else { u2[0] };
    if p1 == c1 || p2 == c2 || p1 == c2 || p2 == c1 {
        return false;
    }
    if brep.coedges[p2].next != p1 {
        return false; // partner does not traverse e2 then e1 → different faces or inconsistent
    }
    let plid = brep.coedges[p1].loop_id;
    if brep.coedges[p2].loop_id != plid {
        return false;
    }
    // same two faces on both sides (the partner loop may be this very loop —
    // a seam slit — which is fine)
    let _ = lid;

    // ── geometry ────────────────────────────────────────────────────────────
    let geom_a = brep.edges[e1].geom.clone();
    let geom_b = brep.edges[e2].geom.clone();
    let pa = brep.vertices[a].point;
    let pv = brep.vertices[v].point;
    let pb = brep.vertices[b].point;
    let (geom, t_start, t_end) = match (&geom_a, &geom_b) {
        (EdgeGeom::Polyline(_), EdgeGeom::Polyline(_))
        | (EdgeGeom::Polyline(_), EdgeGeom::Line(_))
        | (EdgeGeom::Line(_), EdgeGeom::Polyline(_)) => {
            let mut pts = traversal_points(brep, c1);
            let more = traversal_points(brep, c2);
            for p in more.into_iter().skip(1) {
                if pts.last().map_or(true, |&q| (q - p).length() > 1e-12) {
                    pts.push(p);
                }
            }
            if pts.len() < 2 {
                return false;
            }
            (EdgeGeom::Polyline(pts), 0.0, 1.0)
        }
        (EdgeGeom::Line(_), EdgeGeom::Line(_)) => {
            let d = pb - pa;
            let len = d.length();
            let Some(dir) = d.try_normalize() else { return false };
            // collinear: the middle vertex lies on the chord
            let w = pv - pa;
            let off = (w - dir * w.dot(dir)).length();
            if off > 1e-9 {
                return false;
            }
            (
                EdgeGeom::Line(Line3::new(pa, cadcore_math::UnitVec3::new_unchecked(dir))),
                0.0,
                len,
            )
        }
        (EdgeGeom::Circle(ca), EdgeGeom::Circle(cb)) => {
            let same_axis = ca.frame.z.dot(cb.frame.z).abs() > 1.0 - 1e-9;
            if (ca.frame.origin - cb.frame.origin).length() > 1e-9
                || (ca.radius - cb.radius).abs() > 1e-9
                || !same_axis
            {
                return false;
            }
            // parameterise on circle `ca`, running a → v → b
            let ang = |p: Point3| {
                let w = p - ca.frame.origin;
                ca.frame.y.dot_vec(w).atan2(ca.frame.x.dot_vec(w))
            };
            let t0 = ang(pa);
            let s1 = sweep(brep, c1);
            let s2 = sweep(brep, c2);
            // c2's sweep is measured on cb; if cb's axis is flipped relative to
            // ca its angle runs the other way
            let s2 = if ca.frame.z.dot(cb.frame.z) < 0.0 { -s2 } else { s2 };
            let total = s1 + s2;
            // An arc is written between its two end vertices; past 180° the
            // short way round is the WRONG way, so arcs stay under it.
            if total.abs() > 0.9 * std::f64::consts::PI {
                return false;
            }
            (EdgeGeom::Circle(*ca), t0, t0 + total)
        }
        _ => return false,
    };

    // ── apply ───────────────────────────────────────────────────────────────
    let merged = brep.add_edge(Edge {
        geom,
        v_start: a,
        v_end: b,
        t_start,
        t_end,
        partner: None,
    });
    // this side: c1 takes the merged edge, c2 leaves the ring
    let c2_next = brep.coedges[c2].next;
    brep.coedges[c1].edge = merged;
    brep.coedges[c1].sense = CoEdgeSense::Same;
    brep.coedges[c1].next = c2_next;
    brep.coedges[c2_next].prev = c1;
    if brep.loops[lid].start == c2 {
        brep.loops[lid].start = c1;
    }
    // partner side: p2 (first in its order) takes the merged edge reversed,
    // p1 leaves the ring
    let p1_next = brep.coedges[p1].next;
    brep.coedges[p2].edge = merged;
    brep.coedges[p2].sense = CoEdgeSense::Opposite;
    brep.coedges[p2].next = p1_next;
    brep.coedges[p1_next].prev = p2;
    if brep.loops[plid].start == p1 {
        brep.loops[plid].start = p2;
    }
    brep.coedges.remove(c2);
    brep.coedges.remove(p1);
    brep.edges.remove(e1);
    brep.edges.remove(e2);
    brep.vertices.remove(v);
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrange::filament::Filament;
    use crate::arrange::universal::fuse_filaments;
    use crate::config::UnionConfig;
    use crate::validate::validate_shell;
    use crate::Diagnostics;

    #[test]
    fn merging_keeps_the_fused_pair_watertight_and_shrinks_it() {
        let r = 0.275;
        let a = Filament::serpentine(
            &[
                Point3::new(-2.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
                Point3::new(-2.0, 2.0, 0.0),
            ],
            r,
            0.5,
        );
        let b = Filament::serpentine(
            &[
                Point3::new(-1.5, -1.0, 0.35),
                Point3::new(-1.5, 3.0, 0.35),
                Point3::new(-0.8, 3.0, 0.35),
                Point3::new(-0.8, -1.0, 0.35),
            ],
            r,
            0.3,
        );
        let mut brep = BRep::new();
        let shell = fuse_filaments(&mut brep, &[a, b]);
        let before = brep.edges.len();
        let removed = merge_shell_edges(&mut brep, shell);
        let after = brep.edges.len();
        eprintln!("edges {before} → {after} ({removed} merged)");
        assert!(after * 4 < before, "merging must remove most step edges ({before} → {after})");
        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let rep = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(rep.manifold_violations.is_empty(), "manifold after merge: {:?}", &rep.manifold_violations[..rep.manifold_violations.len().min(5)]);
        assert!(rep.wire_violations.is_empty(), "wires after merge: {:?}", &rep.wire_violations[..rep.wire_violations.len().min(5)]);
        assert!(rep.distance_violations.is_empty(), "distance after merge");
    }
}
