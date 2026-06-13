//! Gates 2+3 — uv wire sanity per face.
//!
//! SolidWorks rejects faces whose trimming wires misbehave in parameter
//! space, even when OpenCASCADE accepts them (README invariants #4, #5):
//!
//! * every wire must CLOSE in uv;
//! * wires must be *physically* wound — outer counter-clockwise, holes
//!   clockwise, as seen along the outward face normal (orientation flags are
//!   ignored by SolidWorks' importer);
//! * no wire may intersect another wire or itself.
//!
//! Port of `out/strict_audit.py` (cylinders) and `out/audit_planes.py`
//! (planes), operating on the in-memory B-Rep.  Faces whose surfaces have no
//! natural global uv here (sphere, torus) are skipped for now — torus faces
//! get their own treatment in a later phase.

use std::f64::consts::{PI, TAU};

use cadcore_topo::{BRep, CoEdgeSense, EdgeGeom, FaceGeom, FaceId, FaceNormal, LoopId, ShellId};
use cadcore_math::Point3;

/// One wire's uv polyline plus bookkeeping.
struct WireUv {
    lp: LoopId,
    /// Chained uv points (u branch-continuous on periodic surfaces).
    pts: Vec<(f64, f64)>,
    /// Total u advance around the loop (≈ ±2π for boundary wires on a
    /// cylinder, ≈ 0 for compact wires).
    winding: f64,
    /// Signed area (shoelace) in uv.
    area: f64,
    /// Is this the face's outer loop?
    is_outer: bool,
}

/// Sample an edge's 3-D geometry (oriented `v_start → v_end`).
fn edge_points(brep: &BRep, eid: cadcore_topo::EdgeId) -> Vec<Point3> {
    let Some(e) = brep.edges.get(eid) else {
        return Vec::new();
    };
    match &e.geom {
        EdgeGeom::Line(_) => {
            let a = brep.vertices.get(e.v_start).map(|v| v.point);
            let b = brep.vertices.get(e.v_end).map(|v| v.point);
            [a, b].into_iter().flatten().collect()
        }
        EdgeGeom::Circle(c) => {
            let n = 24;
            (0..=n)
                .map(|k| c.point_at(e.t_start + (e.t_end - e.t_start) * k as f64 / n as f64))
                .collect()
        }
        EdgeGeom::Ellipse(el) => {
            let n = 24;
            (0..=n)
                .map(|k| el.point_at(e.t_start + (e.t_end - e.t_start) * k as f64 / n as f64))
                .collect()
        }
        EdgeGeom::Polyline(pts) => pts.clone(),
    }
}

/// uv projection for the supported surfaces.  `periodic_u` says whether u is
/// an angle (cylinder) needing branch chaining.
enum UvMap {
    Cylinder(cadcore_geom::CylSurf),
    Plane(cadcore_geom::Plane3),
}

impl UvMap {
    fn of(geom: &FaceGeom) -> Option<Self> {
        match geom {
            FaceGeom::Cylinder(c) => Some(UvMap::Cylinder(*c)),
            FaceGeom::Plane(p) => Some(UvMap::Plane(*p)),
            _ => None,
        }
    }

    fn periodic_u(&self) -> bool {
        matches!(self, UvMap::Cylinder(_))
    }

    fn uv(&self, p: Point3) -> (f64, f64) {
        match self {
            UvMap::Cylinder(s) => {
                let w = p - s.frame.origin;
                let u = s.frame.y.dot_vec(w).atan2(s.frame.x.dot_vec(w));
                let v = s.axis().dot_vec(w);
                (u, v)
            }
            UvMap::Plane(pl) => {
                let w = p - pl.frame.origin;
                (pl.frame.x.dot_vec(w), pl.frame.y.dot_vec(w))
            }
        }
    }
}

fn shoelace(pts: &[(f64, f64)]) -> f64 {
    let n = pts.len();
    let mut s = 0.0;
    for i in 0..n {
        let (x0, y0) = pts[i];
        let (x1, y1) = pts[(i + 1) % n];
        s += x0 * y1 - x1 * y0;
    }
    s / 2.0
}

fn seg_intersect(p1: (f64, f64), p2: (f64, f64), q1: (f64, f64), q2: (f64, f64)) -> bool {
    let d1 = (p2.0 - p1.0, p2.1 - p1.1);
    let d2 = (q2.0 - q1.0, q2.1 - q1.1);
    let den = d1.0 * d2.1 - d1.1 * d2.0;
    if den.abs() < 1e-15 {
        return false;
    }
    let t = ((q1.0 - p1.0) * d2.1 - (q1.1 - p1.1) * d2.0) / den;
    let u = ((q1.0 - p1.0) * d1.1 - (q1.1 - p1.1) * d1.0) / den;
    t > 1e-9 && t < 1.0 - 1e-9 && u > 1e-9 && u < 1.0 - 1e-9
}

/// Build the chained uv polyline of one loop.  Returns `None` when the loop
/// is empty; `Err(gap)` style issues are reported via the out-param.
fn wire_uv(
    brep: &BRep,
    map: &UvMap,
    lp: LoopId,
    issues: &mut Vec<String>,
) -> Option<Vec<(f64, f64)>> {
    let loop_rec = brep.loops.get(lp)?;
    let mut chain: Vec<(f64, f64)> = Vec::new();
    let mut chain3: Vec<Point3> = Vec::new();
    let start = loop_rec.start;
    let mut c = start;
    let mut guard = 0usize;
    loop {
        let ce = brep.coedges.get(c)?;
        let mut pts3 = edge_points(brep, ce.edge);
        if ce.sense == CoEdgeSense::Opposite {
            pts3.reverse();
        }
        // The writer normalises traversal senses geometrically (shared-edge
        // caches), so CoEdge.sense is a HINT here, not truth: re-orient by
        // 3-D continuity with the chain tail, like the external auditors do.
        if let Some(&last3) = chain3.last() {
            if !pts3.is_empty() {
                let d_fwd = (pts3[0] - last3).length();
                let d_rev = (*pts3.last().unwrap() - last3).length();
                if d_rev < d_fwd {
                    pts3.reverse();
                }
            }
        }
        chain3.extend(pts3.iter().copied());
        let mut seg: Vec<(f64, f64)> = pts3.iter().map(|&p| map.uv(p)).collect();
        if map.periodic_u() {
            // unwrap u within the segment
            for i in 1..seg.len() {
                let mut du = seg[i].0 - seg[i - 1].0;
                while du > PI {
                    du -= TAU;
                }
                while du < -PI {
                    du += TAU;
                }
                seg[i].0 = seg[i - 1].0 + du;
            }
            // shift the whole segment onto the chain's branch
            if let Some(&(lu, _)) = chain.last() {
                let mut du = seg[0].0 - lu;
                let mut shift = 0.0;
                while du > PI {
                    du -= TAU;
                    shift -= TAU;
                }
                while du < -PI {
                    du += TAU;
                    shift += TAU;
                }
                for p in &mut seg {
                    p.0 += shift;
                }
            }
        }
        if let Some(&(lu, lv)) = chain.last() {
            let gap = ((seg[0].0 - lu).powi(2) + (seg[0].1 - lv).powi(2)).sqrt();
            if gap > 0.05 {
                issues.push(format!("loop {lp:?}: uv gap {gap:.4} between edges"));
            }
        }
        chain.extend(seg);
        c = ce.next;
        guard += 1;
        if c == start || guard > 100_000 {
            break;
        }
    }
    Some(chain)
}

/// Run gates 2+3 on every supported face of `shell`.
pub fn check(brep: &BRep, shell: ShellId) -> Vec<String> {
    let mut issues = Vec::new();
    let Some(sh) = brep.shells.get(shell) else {
        return vec![format!("shell {shell:?} does not exist")];
    };
    for &fid in &sh.faces {
        check_face(brep, fid, &mut issues);
    }
    issues
}

/// Run gates 2+3 on one face.
pub fn check_face(brep: &BRep, fid: FaceId, issues: &mut Vec<String>) {
    let Some(face) = brep.faces.get(fid) else {
        return;
    };
    let Some(map) = UvMap::of(&face.geom) else {
        return; // torus/sphere handled in a later phase
    };
    // outward normal sign: physical winding is CCW about the OUTWARD normal
    let outward_sign = match face.normal {
        FaceNormal::Same => 1.0,
        FaceNormal::Reversed => -1.0,
    };

    let mut wires: Vec<WireUv> = Vec::new();
    let mut loops = vec![(face.outer_loop, true)];
    loops.extend(face.inner_loops.iter().map(|&l| (l, false)));
    for (lp, is_outer) in loops {
        let Some(pts) = wire_uv(brep, &map, lp, issues) else {
            continue;
        };
        if pts.len() < 2 {
            continue;
        }
        let first = pts[0];
        let last = *pts.last().unwrap();
        let mut wind = last.0 - first.0;
        let mut close_du = first.0 - last.0;
        if map.periodic_u() {
            while close_du > PI {
                close_du -= TAU;
            }
            while close_du < -PI {
                close_du += TAU;
            }
            wind += close_du;
        }
        let gap = (close_du.powi(2) + (first.1 - last.1).powi(2)).sqrt();
        if gap > 0.05 {
            issues.push(format!(
                "face {fid:?} loop {lp:?}: wire not closed in uv (gap {gap:.4})"
            ));
        }
        let area = shoelace(&pts);
        wires.push(WireUv {
            lp,
            pts,
            winding: wind,
            area,
            is_outer,
        });
    }

    // physical winding of compact wires
    for w in &wires {
        if map.periodic_u() && w.winding.abs() > PI {
            continue; // boundary wire on a periodic surface — no area role
        }
        let eff = w.area * outward_sign;
        if w.is_outer && eff < -1e-9 {
            issues.push(format!(
                "face {fid:?} loop {:?}: OUTER wire physically clockwise (area {:.4})",
                w.lp, w.area
            ));
        }
        if !w.is_outer && eff > 1e-9 {
            issues.push(format!(
                "face {fid:?} loop {:?}: HOLE wire physically counter-clockwise (area {:.4})",
                w.lp, w.area
            ));
        }
    }

    // wire × wire intersections (with ±2π branch shifts on periodic surfaces)
    let shifts: &[f64] = if map.periodic_u() {
        &[-TAU, 0.0, TAU]
    } else {
        &[0.0]
    };
    for i in 0..wires.len() {
        for j in (i + 1)..wires.len() {
            'shift: for &sh in shifts {
                let a = &wires[i].pts;
                let b = &wires[j].pts;
                for ai in 0..a.len().saturating_sub(1) {
                    for bi in 0..b.len().saturating_sub(1) {
                        let q1 = (b[bi].0 + sh, b[bi].1);
                        let q2 = (b[bi + 1].0 + sh, b[bi + 1].1);
                        if seg_intersect(a[ai], a[ai + 1], q1, q2) {
                            issues.push(format!(
                                "face {fid:?}: loops {:?} and {:?} INTERSECT in uv",
                                wires[i].lp, wires[j].lp
                            ));
                            break 'shift;
                        }
                    }
                }
            }
        }
    }
    // self-intersection
    for w in &wires {
        let a = &w.pts;
        let n = a.len();
        'outer: for ai in 0..n.saturating_sub(1) {
            for bi in (ai + 2)..n.saturating_sub(1) {
                if ai == 0 && bi == n - 2 {
                    continue;
                }
                if seg_intersect(a[ai], a[ai + 1], a[bi], a[bi + 1]) {
                    issues.push(format!(
                        "face {fid:?} loop {:?}: wire SELF-INTERSECTS in uv",
                        w.lp
                    ));
                    break 'outer;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_geom::{Circle3, Plane3};
    use cadcore_math::UnitVec3;
    use cadcore_topo::{
        BRep, CoEdge, CoEdgeId, Edge, Face, FaceExtent, Loop, Shell, Vertex,
    };

    /// Helper: single-circle loop on a plane disk.
    fn disk_face(b: &mut BRep, radius: f64, hole: Option<(f64, bool)>) -> FaceId {
        let plane = Plane3::from_origin_normal(Point3::new(0.0, 0.0, 0.0), UnitVec3::Z);
        let mut mk_loop = |r: f64, ccw: bool| -> LoopId {
            // circle about +Z; reversing the coedge makes it physically CW
            let circ = Circle3::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, r);
            let v = b.vertices.insert(Vertex {
                point: circ.point_at(0.0),
            });
            let e = b.edges.insert(Edge {
                geom: EdgeGeom::Circle(circ),
                v_start: v,
                v_end: v,
                t_start: 0.0,
                t_end: TAU,
                partner: None,
            });
            let lid = b.loops.insert(Loop {
                start: CoEdgeId::default(),
                face: Default::default(),
            });
            let ce = b.coedges.insert(CoEdge {
                edge: e,
                sense: if ccw {
                    CoEdgeSense::Same
                } else {
                    CoEdgeSense::Opposite
                },
                next: CoEdgeId::default(),
                prev: CoEdgeId::default(),
                loop_id: lid,
            });
            b.coedges[ce].next = ce;
            b.coedges[ce].prev = ce;
            b.loops[lid].start = ce;
            lid
        };
        let outer = mk_loop(radius, true);
        let inner: Vec<LoopId> = hole.map(|(r, ccw)| mk_loop(r, ccw)).into_iter().collect();
        b.faces.insert(Face {
            geom: FaceGeom::Plane(plane),
            normal: FaceNormal::Same,
            outer_loop: outer,
            inner_loops: inner,
            shell: Default::default(),
            extent: FaceExtent::None,
        })
    }

    #[test]
    fn good_disk_with_cw_hole_passes() {
        let mut b = BRep::new();
        let f = disk_face(&mut b, 1.0, Some((0.3, false)));
        let mut issues = Vec::new();
        check_face(&b, f, &mut issues);
        assert!(issues.is_empty(), "{issues:?}");
    }

    #[test]
    fn ccw_hole_flagged() {
        let mut b = BRep::new();
        let f = disk_face(&mut b, 1.0, Some((0.3, true)));
        let mut issues = Vec::new();
        check_face(&b, f, &mut issues);
        assert_eq!(issues.len(), 1, "{issues:?}");
        assert!(issues[0].contains("counter-clockwise"));
    }

    #[test]
    fn overlapping_wires_flagged() {
        // hole radius larger than... make two intersecting holes by giving
        // the hole the same radius as the outer wire shifted — simplest:
        // outer r=1 and "hole" r=1.2 (CW so winding role passes) → wires
        // intersect?  Concentric circles never intersect; instead build a
        // hole centred off-origin crossing the outer circle.
        let mut b = BRep::new();
        let plane = Plane3::from_origin_normal(Point3::new(0.0, 0.0, 0.0), UnitVec3::Z);
        let mut mk = |centre: Point3, r: f64, ccw: bool| -> LoopId {
            let circ = Circle3::new(centre, UnitVec3::Z, r);
            let v = b.vertices.insert(Vertex {
                point: circ.point_at(0.0),
            });
            let e = b.edges.insert(Edge {
                geom: EdgeGeom::Circle(circ),
                v_start: v,
                v_end: v,
                t_start: 0.0,
                t_end: TAU,
                partner: None,
            });
            let lid = b.loops.insert(Loop {
                start: CoEdgeId::default(),
                face: Default::default(),
            });
            let ce = b.coedges.insert(CoEdge {
                edge: e,
                sense: if ccw {
                    CoEdgeSense::Same
                } else {
                    CoEdgeSense::Opposite
                },
                next: CoEdgeId::default(),
                prev: CoEdgeId::default(),
                loop_id: lid,
            });
            b.coedges[ce].next = ce;
            b.coedges[ce].prev = ce;
            b.loops[lid].start = ce;
            lid
        };
        let outer = mk(Point3::new(0.0, 0.0, 0.0), 1.0, true);
        let hole = mk(Point3::new(0.9, 0.0, 0.0), 0.3, false); // sticks out
        let f = b.faces.insert(Face {
            geom: FaceGeom::Plane(plane),
            normal: FaceNormal::Same,
            outer_loop: outer,
            inner_loops: vec![hole],
            shell: Default::default(),
            extent: FaceExtent::None,
        });
        let mut issues = Vec::new();
        check_face(&b, f, &mut issues);
        assert!(
            issues.iter().any(|s| s.contains("INTERSECT")),
            "{issues:?}"
        );
    }

    #[test]
    fn shell_runner_walks_all_faces(){
        let mut b = BRep::new();
        let good = disk_face(&mut b, 1.0, None);
        let bad = disk_face(&mut b, 1.0, Some((0.4, true)));
        let shell = b.shells.insert(Shell {
            faces: vec![good, bad],
            is_outer: true,
            solid: Default::default(),
        });
        let issues = check(&b, shell);
        assert_eq!(issues.len(), 1, "{issues:?}");
    }
}
