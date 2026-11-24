//! Gate 4 — edge-to-surface deviation.
//!
//! Every edge must lie within the knit tolerance of **every** face it bounds.
//! This is the gate that finally explained the "vanishing filaments": OCC
//! accepts µm-level drift between trim curves and surfaces, SolidWorks does
//! not.  Port of `out/dist_to_spline.py` (the analytic half; spline faces are
//! emitted as exact rational torus patches, so checking against the analytic
//! torus is equivalent).

use std::collections::HashMap;

use crate::geom::refine::AnalyticSurface;
use cadcore_topo::{BRep, EdgeGeom, EdgeId, FaceGeom, FaceId, ShellId};
use cadcore_math::Point3;

/// Sample an edge's geometry as 3-D points (coarse but sufficient: trim
/// polylines carry their defining points, analytic curves are exact).
fn sample_edge(brep: &BRep, eid: EdgeId) -> Vec<Point3> {
    let Some(e) = brep.edges.get(eid) else {
        return Vec::new();
    };
    match &e.geom {
        EdgeGeom::Line(_) => {
            let a = brep.vertices.get(e.v_start).map(|v| v.point);
            let b = brep.vertices.get(e.v_end).map(|v| v.point);
            [a, b].into_iter().flatten().collect()
        }
        EdgeGeom::Circle(c) => (0..16)
            .map(|k| c.point_at(e.t_start + (e.t_end - e.t_start) * k as f64 / 15.0))
            .collect(),
        EdgeGeom::Ellipse(el) => (0..16)
            .map(|k| el.point_at(e.t_start + (e.t_end - e.t_start) * k as f64 / 15.0))
            .collect(),
        EdgeGeom::Polyline(pts) => pts.clone(),
    }
}

/// Analytic surface of a face, when the writer emits one (or an exact
/// rational equivalent, as it does for torus patches).
fn face_surface(geom: &FaceGeom) -> Option<AnalyticSurface> {
    match geom {
        FaceGeom::Plane(p) => Some(AnalyticSurface::Plane(*p)),
        FaceGeom::Cylinder(c) => Some(AnalyticSurface::Cylinder(*c)),
        FaceGeom::Torus(t) => Some(AnalyticSurface::Torus(*t)),
        FaceGeom::Sphere(_) => None,
    }
}

/// Largest distance of the edge's two vertices from the curve it is drawn on.
fn vertex_off_curve(brep: &BRep, eid: EdgeId) -> Option<f64> {
    let e = brep.edges.get(eid)?;
    let ends = [brep.vertices.get(e.v_start)?.point, brep.vertices.get(e.v_end)?.point];
    let dist = |p: Point3| -> f64 {
        match &e.geom {
            EdgeGeom::Circle(c) => {
                let w = p - c.frame.origin;
                let planar = c.frame.z.dot_vec(w);
                let radial = (w - c.frame.z.as_vec() * planar).length() - c.radius;
                planar.hypot(radial)
            }
            EdgeGeom::Line(l) => {
                let w = p - l.origin;
                (w - l.direction.as_vec() * l.direction.dot_vec(w)).length()
            }
            EdgeGeom::Polyline(pts) => pts
                .first()
                .zip(pts.last())
                .map_or(0.0, |(&a, &b)| (p - a).length().min((p - b).length())),
            EdgeGeom::Ellipse(_) => 0.0,
        }
    };
    Some(dist(ends[0]).max(dist(ends[1])))
}

/// Check the gate.  Returns `(violations, max_deviation_mm)`.
pub fn check(brep: &BRep, shell: ShellId, knit_tol: f64) -> (Vec<String>, f64) {
    // edge -> faces using it
    let mut edge_faces: HashMap<EdgeId, Vec<FaceId>> = HashMap::new();
    let Some(sh) = brep.shells.get(shell) else {
        return (vec![format!("shell {shell:?} does not exist")], 0.0);
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
                edge_faces.entry(ce.edge).or_default().push(fid);
                c = ce.next;
                guard += 1;
                if c == start || guard > 100_000 {
                    break;
                }
            }
        }
    }

    let mut violations = Vec::new();
    let mut max_dev = 0.0f64;
    for (eid, faces) in &edge_faces {
        // An edge runs between its two VERTICES along its curve, so both
        // vertices have to lie on that curve.  A vertex grafted onto a circle
        // it merely passes (8 µm off, in a trimmed ±30° diamond) leaves the
        // curve and the surfaces untouched — the samples below are taken on
        // the curve — yet no reader can build the edge, and the faces around
        // it come out invalid.
        if let Some(d) = vertex_off_curve(brep, *eid) {
            max_dev = max_dev.max(d);
            if d > knit_tol {
                violations.push(format!(
                    "edge {eid:?}: a vertex is {:.2} µm off the edge's curve (limit {:.2} µm)",
                    d * 1e3,
                    knit_tol * 1e3
                ));
            }
        }
        let pts = sample_edge(brep, *eid);
        if pts.is_empty() {
            continue;
        }
        for fid in faces {
            let Some(face) = brep.faces.get(*fid) else {
                continue;
            };
            let Some(surf) = face_surface(&face.geom) else {
                continue;
            };
            let d = pts
                .iter()
                .map(|&p| surf.distance(p))
                .fold(0.0f64, f64::max);
            max_dev = max_dev.max(d);
            if d > knit_tol {
                violations.push(format!(
                    "edge {eid:?} deviates {:.2} µm from face {fid:?} (limit {:.2} µm)",
                    d * 1e3,
                    knit_tol * 1e3
                ));
            }
        }
    }
    violations.sort();
    (violations, max_dev)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_geom::{Circle3, CylSurf};
    use cadcore_math::UnitVec3;
    use cadcore_topo::{
        BRep, CoEdge, CoEdgeId, CoEdgeSense, Edge, Face, FaceExtent, FaceNormal, Loop, Shell,
        Vertex,
    };

    /// A circle edge exactly on its cylinder passes; one displaced by 5 µm in
    /// radius fails.
    #[test]
    fn detects_offset_edge() {
        let mut b = BRep::new();
        let cyl = CylSurf::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, 0.275);

        let mut add = |radius: f64| {
            let circ = Circle3::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, radius);
            let v = b.vertices.insert(Vertex {
                point: circ.point_at(0.0),
            });
            let e = b.edges.insert(Edge {
                geom: cadcore_topo::EdgeGeom::Circle(circ),
                v_start: v,
                v_end: v,
                t_start: 0.0,
                t_end: std::f64::consts::TAU,
                partner: None,
            });
            let lid = b.loops.insert(Loop {
                start: CoEdgeId::default(),
                face: Default::default(),
            });
            let ce = b.coedges.insert(CoEdge {
                edge: e,
                sense: CoEdgeSense::Same,
                next: CoEdgeId::default(),
                prev: CoEdgeId::default(),
                loop_id: lid,
            });
            b.coedges[ce].next = ce;
            b.coedges[ce].prev = ce;
            b.loops[lid].start = ce;
            b.faces.insert(Face {
                geom: cadcore_topo::FaceGeom::Cylinder(cyl),
                normal: FaceNormal::Same,
                outer_loop: lid,
                inner_loops: vec![],
                shell: Default::default(),
                extent: FaceExtent::None,
            })
        };

        let good = add(0.275);
        let bad = add(0.275 + 5e-3); // 5 µm out

        let s_good = b.shells.insert(Shell {
            faces: vec![good],
            is_outer: true,
            solid: Default::default(),
        });
        let s_bad = b.shells.insert(Shell {
            faces: vec![bad],
            is_outer: true,
            solid: Default::default(),
        });

        let (v, dev) = check(&b, s_good, 1e-3);
        assert!(v.is_empty(), "{v:?}");
        assert!(dev < 1e-9);

        let (v, dev) = check(&b, s_bad, 1e-3);
        assert_eq!(v.len(), 1);
        assert!((dev - 5e-3).abs() < 1e-9, "dev={dev}");
    }

    /// A vertex that is not on its edge's curve fails, even though the curve
    /// itself lies exactly on the face.
    #[test]
    fn detects_a_vertex_off_its_edge_curve() {
        let mut b = BRep::new();
        let cyl = CylSurf::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, 0.275);
        let circ = Circle3::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, 0.275);
        let mut add = |off: f64| {
            // an arc of the rim whose END vertex was moved `off` along the axis
            let v0 = b.vertices.insert(Vertex { point: circ.point_at(0.0) });
            let v1 = b.vertices.insert(Vertex { point: circ.point_at(1.0) + UnitVec3::Z * off });
            let e = b.edges.insert(Edge {
                geom: cadcore_topo::EdgeGeom::Circle(circ),
                v_start: v0,
                v_end: v1,
                t_start: 0.0,
                t_end: 1.0,
                partner: None,
            });
            let lid = b.loops.insert(Loop { start: CoEdgeId::default(), face: Default::default() });
            let ce = b.coedges.insert(CoEdge {
                edge: e,
                sense: CoEdgeSense::Same,
                next: CoEdgeId::default(),
                prev: CoEdgeId::default(),
                loop_id: lid,
            });
            b.coedges[ce].next = ce;
            b.coedges[ce].prev = ce;
            b.loops[lid].start = ce;
            let f = b.faces.insert(Face {
                geom: cadcore_topo::FaceGeom::Cylinder(cyl),
                normal: FaceNormal::Same,
                outer_loop: lid,
                inner_loops: vec![],
                shell: Default::default(),
                extent: FaceExtent::None,
            });
            b.shells.insert(Shell { faces: vec![f], is_outer: true, solid: Default::default() })
        };
        let good = add(0.0);
        let bad = add(8.0e-3); // 8 µm, the trimmed ±30° diamond
        assert!(check(&b, good, 1e-3).0.is_empty());
        let (v, dev) = check(&b, bad, 1e-3);
        assert_eq!(v.len(), 1, "{v:?}");
        assert!(v[0].contains("off the edge's curve"), "{v:?}");
        assert!((dev - 8.0e-3).abs() < 1e-9, "dev={dev}");
    }
}
