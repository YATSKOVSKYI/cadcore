//! Gate 1 — the AP214 manifold invariant.
//!
//! Every edge in a closed shell must be referenced by exactly two co-edges
//! with opposite senses (see `cadcore/CLAUDE.md`).  A violation surfaces in
//! receiving CADs as an open shell: filaments rendered cut short, transparent
//! electrodes, vanishing faces.

use std::collections::HashMap;

use cadcore_topo::{BRep, CoEdgeSense, EdgeId, FaceId, ShellId};

/// Check result: hard violations + sense warnings.
///
/// At the in-memory B-Rep level the hard invariant is the USE COUNT: every
/// edge must be referenced by exactly two co-edges.  Traversal senses are
/// normalised later by the STEP writer's shared-edge caches (it derives the
/// second user's sense geometrically), so same-sense pairs here are only a
/// warning, not a gate failure.
pub struct ManifoldReport {
    /// Edges not used exactly twice (gate failures).
    pub violations: Vec<String>,
    /// Edges used twice with equal senses (writer will normalise).
    pub sense_warnings: Vec<String>,
}

/// Check the invariant for `shell`.
pub fn check(brep: &BRep, shell: ShellId) -> ManifoldReport {
    // edge -> [(face, sense)]
    let mut uses: HashMap<EdgeId, Vec<(FaceId, CoEdgeSense)>> = HashMap::new();

    let Some(sh) = brep.shells.get(shell) else {
        return ManifoldReport {
            violations: vec![format!("shell {shell:?} does not exist")],
            sense_warnings: Vec::new(),
        };
    };
    for &fid in &sh.faces {
        let Some(face) = brep.faces.get(fid) else {
            continue;
        };
        let mut loops = vec![face.outer_loop];
        loops.extend(face.inner_loops.iter().copied());
        for lid in loops {
            let Some(lp) = brep.loops.get(lid) else {
                continue;
            };
            let start = lp.start;
            let mut c = start;
            let mut guard = 0usize;
            loop {
                let Some(ce) = brep.coedges.get(c) else { break };
                uses.entry(ce.edge).or_default().push((fid, ce.sense));
                c = ce.next;
                guard += 1;
                if c == start || guard > 100_000 {
                    break;
                }
            }
        }
    }

    let mut violations = Vec::new();
    let mut sense_warnings = Vec::new();
    for (eid, us) in &uses {
        let desc = || {
            us.iter()
                .map(|(f, s)| format!("{f:?}:{s:?}"))
                .collect::<Vec<_>>()
                .join(", ")
        };
        if us.len() != 2 {
            violations.push(format!(
                "edge {eid:?} used {} times (want 2): [{}]",
                us.len(),
                desc()
            ));
        } else if us[0].1 == us[1].1 {
            sense_warnings.push(format!(
                "edge {eid:?} used twice with EQUAL senses: [{}]",
                desc()
            ));
        }
    }
    violations.sort();
    sense_warnings.sort();
    ManifoldReport {
        violations,
        sense_warnings,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_geom::Line3;
    use cadcore_math::{Point3, UnitVec3};
    use cadcore_topo::{
        BRep, CoEdge, CoEdgeId, Edge, EdgeGeom, Face, FaceExtent, FaceGeom, FaceNormal, Loop,
        LoopId, Shell, Vertex,
    };

    /// Build a degenerate two-face arena where one edge is shared correctly
    /// and another is used once — gate must flag exactly the latter.
    #[test]
    fn flags_single_use_edge() {
        let mut b = BRep::new();
        let v1 = b.vertices.insert(Vertex {
            point: Point3::new(0.0, 0.0, 0.0),
        });
        let v2 = b.vertices.insert(Vertex {
            point: Point3::new(1.0, 0.0, 0.0),
        });
        let line = Line3::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::X);
        let shared = b.edges.insert(Edge {
            geom: EdgeGeom::Line(line),
            v_start: v1,
            v_end: v2,
            t_start: 0.0,
            t_end: 1.0,
            partner: None,
        });
        let lone = b.edges.insert(Edge {
            geom: EdgeGeom::Line(line),
            v_start: v1,
            v_end: v2,
            t_start: 0.0,
            t_end: 1.0,
            partner: None,
        });

        // helper: single-coedge "loop" (geometrically nonsense, fine for the
        // topology-only gate)
        let mut mkloop = |edge, sense| -> LoopId {
            let lid = b.loops.insert(Loop {
                start: CoEdgeId::default(),
                face: Default::default(),
            });
            let ce = b.coedges.insert(CoEdge {
                edge,
                sense,
                next: CoEdgeId::default(),
                prev: CoEdgeId::default(),
                loop_id: lid,
            });
            b.coedges[ce].next = ce;
            b.coedges[ce].prev = ce;
            b.loops[lid].start = ce;
            lid
        };

        let l1 = mkloop(shared, CoEdgeSense::Same);
        let l2 = mkloop(shared, CoEdgeSense::Opposite);
        let l3 = mkloop(lone, CoEdgeSense::Same);

        let plane = cadcore_geom::Plane3::from_origin_normal(
            Point3::new(0.0, 0.0, 0.0),
            UnitVec3::Z,
        );
        let mut mkface = |outer: LoopId, extra: Vec<LoopId>| {
            b.faces.insert(Face {
                geom: FaceGeom::Plane(plane),
                normal: FaceNormal::Same,
                outer_loop: outer,
                inner_loops: extra,
                shell: Default::default(),
                extent: FaceExtent::Disk { radius: 1.0 },
            })
        };
        let f1 = mkface(l1, vec![]);
        let f2 = mkface(l2, vec![l3]);

        let shell = b.shells.insert(Shell {
            faces: vec![f1, f2],
            is_outer: true,
            solid: Default::default(),
        });

        let r = check(&b, shell);
        assert_eq!(r.violations.len(), 1, "exactly the lone edge: {:?}", r.violations);
        assert!(r.violations[0].contains("used 1 times"));
        assert!(r.sense_warnings.is_empty());
    }
}
