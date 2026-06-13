//! Assembly — kept cells become cadcore-topo faces, loops and co-edges.
//!
//! This is the final arrangement stage: it materialises the 2-D result into
//! topology that the STEP writer can serialise.  The watertight invariant is
//! enforced by a **shared-edge cache** keyed on the 3-D geometry of each edge
//! (`EdgeKey`): the first user creates the [`Edge`], the second references it
//! with the opposite [`CoEdgeSense`].  Because the registry already split
//! every curve once and refined joints to a single source of truth, the two
//! faces meeting at an edge present byte-identical endpoints, so the cache
//! collapses them into exactly one twice-used edge — README invariants #2,#3.

use std::collections::HashMap;

use cadcore_geom::arrangement::LoopStep;
use cadcore_geom::{Circle3, CylSurf, Plane3, TorusSurf};
use cadcore_math::Point3;
use cadcore_topo::{
    BRep, CoEdge, CoEdgeId, CoEdgeSense, Edge, EdgeGeom, EdgeId, Face, FaceExtent, FaceGeom,
    FaceNormal, Loop, LoopId, VertexId,
};

use super::cells::ClassifiedCell;
use super::domain::FaceDomain;

/// Quantised 3-D point used as a vertex/edge dictionary key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct PtKey(i64, i64, i64);

/// If `poly` is a closed, planar, near-constant-radius loop, fit and return the
/// circle it traces.  Used to weld rim / junction circles by geometry (centre,
/// sign-normalised axis, radius) regardless of where each face's seam split it.
///
/// The axis is sign-normalised (largest component positive) so two faces that
/// wind the circle in opposite directions fit the IDENTICAL `Circle3` and thus
/// share one analytic edge with correct opposite senses.
fn fit_circle(poly: &[Point3]) -> Option<Circle3> {
    let mut pts = poly.to_vec();
    if pts.len() > 1 && (pts[0] - *pts.last().unwrap()).length() < 1e-9 {
        pts.pop();
    }
    let n = pts.len();
    if n < 6 {
        return None;
    }
    // centroid
    let mut c = cadcore_math::Vec3::new(0.0, 0.0, 0.0);
    for p in &pts {
        c = c + (*p - Point3::new(0.0, 0.0, 0.0));
    }
    let centre = Point3::new(0.0, 0.0, 0.0) + c * (1.0 / n as f64);
    // radius mean + variance check
    let radii: Vec<f64> = pts.iter().map(|p| (*p - centre).length()).collect();
    let r = radii.iter().sum::<f64>() / n as f64;
    if r < 1e-9 {
        return None;
    }
    if radii.iter().any(|d| (d - r).abs() > 1e-3 * r) {
        return None; // not a constant-radius loop (e.g. a window / miter)
    }
    // plane normal via the summed cross products of consecutive spokes
    let mut nrm = cadcore_math::Vec3::new(0.0, 0.0, 0.0);
    for i in 0..n {
        let a = pts[i] - centre;
        let b = pts[(i + 1) % n] - centre;
        nrm = nrm + a.cross(b);
    }
    let axis = nrm.try_normalize()?;
    // planarity: every point close to the plane through `centre`
    if pts.iter().any(|p| axis.dot(*p - centre).abs() > 1e-3 * r) {
        return None;
    }
    // sign-normalise the axis (largest |component| positive)
    let (ax, ay, az) = (axis.x, axis.y, axis.z);
    let flip = if ax.abs() >= ay.abs() && ax.abs() >= az.abs() {
        ax < 0.0
    } else if ay.abs() >= az.abs() {
        ay < 0.0
    } else {
        az < 0.0
    };
    let axis = if flip { axis * -1.0 } else { axis };
    Some(Circle3::new(centre, cadcore_math::UnitVec3::try_from_vec(axis)?, r))
}

/// Midpoint of a polyline by arc length (independent of traversal sense).
fn arc_midpoint(poly: &[Point3]) -> Point3 {
    let total: f64 = poly.windows(2).map(|w| (w[1] - w[0]).length()).sum();
    if total < 1e-30 {
        return poly[0];
    }
    let half = total * 0.5;
    let mut acc = 0.0;
    for w in poly.windows(2) {
        let seg = (w[1] - w[0]).length();
        if acc + seg >= half {
            let t = (half - acc) / seg.max(1e-30);
            return w[0] + (w[1] - w[0]) * t;
        }
        acc += seg;
    }
    poly[poly.len() / 2]
}

fn pt_key(p: Point3, tol: f64) -> PtKey {
    let q = 1.0 / tol;
    PtKey(
        (p.x * q).round() as i64,
        (p.y * q).round() as i64,
        (p.z * q).round() as i64,
    )
}

/// Canonical, direction-agnostic key for a shared edge: the two endpoint
/// keys sorted, plus a quantised midpoint so two different curves between the
/// same endpoints (e.g. the two halves of a split circle) stay distinct.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct EdgeKey {
    lo: PtKey,
    hi: PtKey,
    mid: PtKey,
}

/// Shared topology builder across all faces of one union solid.
///
/// Hold ONE of these for the whole shell so edges are shared between
/// neighbouring faces; dropping it per-face would make every edge single-use.
pub struct Assembler<'b> {
    brep: &'b mut BRep,
    tol: f64,
    verts: HashMap<PtKey, VertexId>,
    /// edge key → (edge id, first user's start vertex)
    edges: HashMap<EdgeKey, (EdgeId, VertexId)>,
    /// circle key → (edge id, first user's start vertex) — analytic rim/cap
    /// circles, shared SEAMLESSLY between a cylinder rim and its cap / elbow.
    circles: HashMap<CircleKey, (EdgeId, VertexId)>,
    faces: Vec<cadcore_topo::FaceId>,
}

/// Direction-agnostic canonical key for a full circle (centre, axis-line,
/// radius).  The axis sign is normalised so a rim and its cap (opposite
/// normals) share one edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct CircleKey {
    centre: PtKey,
    axis: (i64, i64, i64),
    radius: i64,
}

impl<'b> Assembler<'b> {
    /// New assembler writing into `brep`.
    pub fn new(brep: &'b mut BRep, tol: f64) -> Self {
        Self {
            brep,
            tol,
            verts: HashMap::new(),
            edges: HashMap::new(),
            circles: HashMap::new(),
            faces: Vec::new(),
        }
    }

    /// Faces produced so far.
    pub fn faces(&self) -> &[cadcore_topo::FaceId] {
        &self.faces
    }

    fn vertex(&mut self, p: Point3) -> VertexId {
        let k = pt_key(p, self.tol);
        if let Some(&v) = self.verts.get(&k) {
            return v;
        }
        let v = self.brep.add_vertex(cadcore_topo::Vertex { point: p });
        self.verts.insert(k, v);
        v
    }

    /// Get-or-create the shared edge for a 3-D polyline; returns the co-edge
    /// sense this user must traverse it with (`Same` for the first user,
    /// `Opposite` for the second).
    fn shared_edge(&mut self, poly3: &[Point3]) -> (EdgeId, CoEdgeSense) {
        let a = poly3[0];
        let b = *poly3.last().unwrap();
        // geometric arc midpoint by length — direction-agnostic, so the two
        // halves of a split circle stay distinct while a curve and its
        // reverse map to the same key
        let mid = arc_midpoint(poly3);
        let (ka, kb) = (pt_key(a, self.tol), pt_key(b, self.tol));
        let (lo, hi) = if (ka.0, ka.1, ka.2) <= (kb.0, kb.1, kb.2) {
            (ka, kb)
        } else {
            (kb, ka)
        };
        let key = EdgeKey {
            lo,
            hi,
            mid: pt_key(mid, self.tol),
        };
        if let Some(&(eid, first_start)) = self.edges.get(&key) {
            // second user: opposite sense iff it starts at the same vertex
            let my_start = self.vertex(a);
            let sense = if my_start == first_start {
                CoEdgeSense::Same
            } else {
                CoEdgeSense::Opposite
            };
            return (eid, sense);
        }
        let v_start = self.vertex(a);
        let v_end = self.vertex(b);
        let eid = self.brep.add_edge(Edge {
            geom: EdgeGeom::Polyline(poly3.to_vec()),
            v_start,
            v_end,
            t_start: 0.0,
            t_end: 1.0,
            partner: None,
        });
        self.edges.insert(key, (eid, v_start));
        (eid, CoEdgeSense::Same)
    }

    /// Materialise one loop (a sequence of [`LoopStep`]s in uv) into a
    /// cadcore-topo [`Loop`] of co-edges, lifting each step to 3-D via the
    /// face `domain`.  Returns the new loop id.
    fn build_loop(&mut self, domain: &FaceDomain, steps: &[LoopStep], face: cadcore_topo::FaceId) -> LoopId {
        let lp = self.brep.add_loop(Loop {
            start: CoEdgeId::default(),
            face,
        });
        let mut coedges: Vec<CoEdgeId> = Vec::with_capacity(steps.len());
        for step in steps {
            // lift this step's uv polyline to 3-D in traversal order
            let poly3: Vec<Point3> = step.pts.iter().map(|&(u, v)| domain.lift(u, v)).collect();
            if poly3.len() < 2 {
                continue;
            }
            // A full closed circle (a rim / junction circle that survived as one
            // step) is welded by its centre/axis/radius via the analytic circle
            // cache — rotation-invariant, so a cylinder rim and a torus minor
            // circle split at different seams still share ONE edge.
            let (edge, sense) = if let Some(c) = fit_circle(&poly3) {
                self.shared_circle(c)
            } else {
                self.shared_edge(&poly3)
            };
            let ce = self.brep.add_coedge(CoEdge {
                edge,
                sense,
                next: CoEdgeId::default(),
                prev: CoEdgeId::default(),
                loop_id: lp,
            });
            coedges.push(ce);
        }
        // link the ring
        let n = coedges.len();
        for i in 0..n {
            let cur = coedges[i];
            let nxt = coedges[(i + 1) % n];
            let prv = coedges[(i + n - 1) % n];
            self.brep.coedges[cur].next = nxt;
            self.brep.coedges[cur].prev = prv;
        }
        if let Some(&first) = coedges.first() {
            self.brep.loops[lp].start = first;
        }
        lp
    }

    /// Get-or-create a shared analytic full-circle edge; returns the co-edge
    /// sense this user must traverse it with.  Welds a cylinder rim to its
    /// cap / elbow with NO seam (the periodic surface closes implicitly).
    fn shared_circle(&mut self, c: Circle3) -> (EdgeId, CoEdgeSense) {
        let q = 1.0 / self.tol;
        let ax = c.frame.z.as_vec();
        // normalise axis sign (largest component positive) so ±normal match
        let flip = {
            let (x, y, z) = (ax.x, ax.y, ax.z);
            if x.abs() >= y.abs() && x.abs() >= z.abs() {
                x < 0.0
            } else if y.abs() >= z.abs() {
                y < 0.0
            } else {
                z < 0.0
            }
        };
        let a = if flip { ax * -1.0 } else { ax };
        let key = CircleKey {
            centre: pt_key(c.frame.origin, self.tol),
            axis: ((a.x * q).round() as i64, (a.y * q).round() as i64, (a.z * q).round() as i64),
            radius: (c.radius * q).round() as i64,
        };
        let start_pt = c.point_at(0.0);
        if let Some(&(eid, first_start)) = self.circles.get(&key) {
            let my_start = self.vertex(start_pt);
            let sense = if my_start == first_start {
                CoEdgeSense::Same
            } else {
                CoEdgeSense::Opposite
            };
            return (eid, sense);
        }
        let v = self.vertex(start_pt);
        let eid = self.brep.add_edge(Edge {
            geom: EdgeGeom::Circle(c),
            v_start: v,
            v_end: v,
            t_start: 0.0,
            t_end: std::f64::consts::TAU,
            partner: None,
        });
        self.circles.insert(key, (eid, v));
        (eid, CoEdgeSense::Same)
    }

    /// Build a loop holding exactly ONE closed analytic circle edge.
    fn circle_loop(&mut self, c: Circle3, face: cadcore_topo::FaceId) -> LoopId {
        let lp = self.brep.add_loop(Loop { start: CoEdgeId::default(), face });
        let (edge, sense) = self.shared_circle(c);
        let ce = self.brep.add_coedge(CoEdge {
            edge,
            sense,
            next: CoEdgeId::default(),
            prev: CoEdgeId::default(),
            loop_id: lp,
        });
        self.brep.coedges[ce].next = ce;
        self.brep.coedges[ce].prev = ce;
        self.brep.loops[lp].start = ce;
        lp
    }

    /// Emit a SEAMLESS periodic cylinder face: the two v-boundary rims as
    /// analytic circle bounds (shared with caps/elbows) plus the crossing
    /// `windows` as interior holes.  No parametric seam — the cylinder closes
    /// in u implicitly.  Each window is a closed 3-D loop shared with the
    /// crossing tube's matching window (single source of truth → welded).
    pub fn emit_cylinder_face(
        &mut self,
        surf: CylSurf,
        _length: f64,
        rim_lo: Circle3,
        rim_hi: Circle3,
        windows: &[Vec<Point3>],
    ) -> cadcore_topo::FaceId {
        let face = self.brep.add_face(Face {
            geom: FaceGeom::Cylinder(surf),
            normal: FaceNormal::Same,
            outer_loop: LoopId::default(),
            inner_loops: Vec::new(),
            shell: Default::default(),
            extent: FaceExtent::Trimmed,
        });
        // outer = lo rim circle; inner = hi rim circle + window holes
        let outer = self.circle_loop(rim_lo, face);
        let mut inner = vec![self.circle_loop(rim_hi, face)];
        for w in windows {
            inner.push(self.window_hole_loop(w, &surf, face));
        }
        self.brep.faces[face].outer_loop = outer;
        self.brep.faces[face].inner_loops = inner;
        self.faces.push(face);
        face
    }

    /// Build a window hole loop on a cylinder.  The loop is wound CW about the
    /// cylinder's OUTWARD normal (the STEP hole convention) and split into a
    /// FEW shared edges — a single closed edge cannot carry an opposite sense
    /// between the two faces sharing it (start == end vertex), so we split so
    /// the two legs' pieces weld with opposite senses → a manifold edge.
    fn window_hole_loop(&mut self, pts3: &[Point3], surf: &CylSurf, face: cadcore_topo::FaceId) -> LoopId {
        // drop a trailing duplicate of the start
        let mut poly: Vec<Point3> = pts3.to_vec();
        if poly.len() > 1 && (poly[0] - *poly.last().unwrap()).length() < 1e-9 {
            poly.pop();
        }
        let m = poly.len();
        // orient CW about the outward normal: area vector A = 0.5 Σ p_i×p_{i+1}
        let mut area = cadcore_math::Vec3::new(0.0, 0.0, 0.0);
        let o = Point3::new(0.0, 0.0, 0.0);
        for i in 0..m {
            let a = poly[i] - o;
            let b = poly[(i + 1) % m] - o;
            area = area + a.cross(b);
        }
        // outward normal at the window centroid (radial from the cyl axis)
        let cen = {
            let mut c = Point3::new(0.0, 0.0, 0.0);
            for p in &poly { c = c + (*p - o) * (1.0 / m as f64); }
            c
        };
        let w = cen - surf.frame.origin;
        let ax = surf.axis().dot_vec(w);
        let out = w - surf.axis().as_vec() * ax;
        // CW about outward ⇒ traverse the canonical loop BACKWARD when the
        // forward (traced) loop is CCW about the outward normal.
        let backward = area.dot(out) > 0.0;
        // Split the CANONICAL poly (same for both legs → identical 3-D pieces)
        // into fixed pieces; only the TRAVERSAL direction differs per leg, so
        // the two legs weld each shared piece with opposite senses.
        let parts = 4usize;
        let mut pieces: Vec<Vec<Point3>> = Vec::with_capacity(parts);
        for k in 0..parts {
            let i0 = m * k / parts;
            let i1 = m * (k + 1) / parts;
            let mut seg: Vec<Point3> = Vec::new();
            let mut i = i0;
            loop {
                seg.push(poly[i % m]);
                if i == i1 { break; }
                i += 1;
            }
            pieces.push(seg);
        }
        let lp = self.brep.add_loop(Loop { start: CoEdgeId::default(), face });
        let mut coedges: Vec<CoEdgeId> = Vec::new();
        let order: Vec<usize> = if backward {
            (0..parts).rev().collect()
        } else {
            (0..parts).collect()
        };
        for &k in &order {
            let mut seg = pieces[k].clone();
            if backward { seg.reverse(); }
            if seg.len() < 2 { continue; }
            let (edge, sense) = self.shared_edge(&seg);
            let ce = self.brep.add_coedge(CoEdge {
                edge, sense, next: CoEdgeId::default(), prev: CoEdgeId::default(), loop_id: lp,
            });
            coedges.push(ce);
        }
        let n = coedges.len();
        for i in 0..n {
            self.brep.coedges[coedges[i]].next = coedges[(i + 1) % n];
            self.brep.coedges[coedges[i]].prev = coedges[(i + n - 1) % n];
        }
        if let Some(&first) = coedges.first() { self.brep.loops[lp].start = first; }
        lp
    }

    /// Emit a SEAMLESS elbow torus-fillet face bounded by its two junction
    /// minor-circles (shared with the mating legs).  No φ-seam — the writer's
    /// `FaceExtent::TorusFillet` template emits the two circles directly, and
    /// the BRep loops (via `shared_circle`) weld them to the legs.
    pub fn emit_torus_fillet(
        &mut self,
        surf: TorusSurf,
        junction_lo: Circle3,
        junction_hi: Circle3,
    ) -> cadcore_topo::FaceId {
        let face = self.brep.add_face(Face {
            geom: FaceGeom::Torus(surf),
            normal: FaceNormal::Same,
            outer_loop: LoopId::default(),
            inner_loops: Vec::new(),
            shell: Default::default(),
            extent: FaceExtent::TorusFillet {
                start_circle: junction_lo,
                end_circle: junction_hi,
            },
        });
        let outer = self.circle_loop(junction_lo, face);
        let inner = self.circle_loop(junction_hi, face);
        self.brep.faces[face].outer_loop = outer;
        self.brep.faces[face].inner_loops = vec![inner];
        self.faces.push(face);
        face
    }

    /// Emit a flat end-cap disk bounded by a single shared circle (welds to
    /// the cylinder rim with no seam).
    pub fn emit_disk_cap(&mut self, plane: Plane3, rim: Circle3) -> cadcore_topo::FaceId {
        let face = self.brep.add_face(Face {
            geom: FaceGeom::Plane(plane),
            normal: FaceNormal::Same,
            outer_loop: LoopId::default(),
            inner_loops: Vec::new(),
            shell: Default::default(),
            extent: FaceExtent::Disk { radius: rim.radius },
        });
        let outer = self.circle_loop(rim, face);
        self.brep.faces[face].outer_loop = outer;
        self.faces.push(face);
        face
    }

    /// Emit one kept cell as a face on `geom`/`normal`, with its outer loop
    /// and holes.  Boundary chains tagged with the reserved seam id are
    /// lifted like any other (they carry true 3-D geometry at the seam).
    pub fn emit_cell(
        &mut self,
        domain: &FaceDomain,
        geom: FaceGeom,
        normal: FaceNormal,
        cell: &ClassifiedCell,
    ) -> cadcore_topo::FaceId {
        let face = self.brep.add_face(Face {
            geom,
            normal,
            outer_loop: LoopId::default(),
            inner_loops: Vec::new(),
            shell: Default::default(),
            extent: FaceExtent::Trimmed,
        });
        let outer = self.build_loop(domain, &cell.cell.outer, face);
        let holes: Vec<LoopId> = cell
            .cell
            .holes
            .iter()
            .map(|h| self.build_loop(domain, h, face))
            .collect();
        self.brep.faces[face].outer_loop = outer;
        self.brep.faces[face].inner_loops = holes;
        self.faces.push(face);
        face
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrange::cells::{arrange_and_classify, FaceChain};
    use crate::geom::composite::{closed_loops_between, CompositeTube};
    use crate::geom::intersect::TraceOptions;
    use cadcore_geom::CylSurf;
    use cadcore_math::UnitVec3;
    use std::collections::HashMap;
    use std::f64::consts::{PI, TAU};

    /// One leg crossed by a perpendicular leg: assemble its lateral band as
    /// faces and assert the window edges are SHARED (used exactly twice)
    /// between the kept band and the dropped window — the watertight unit.
    #[test]
    fn crossing_band_assembles_shared_edges() {
        let leg = CylSurf::new(Point3::new(-2.0, 0.0, 0.0), UnitVec3::X, 0.275);
        let a = CompositeTube::new().with_leg(leg, 4.0);
        let b = CompositeTube::new().with_leg(
            CylSurf::new(Point3::new(0.0, -2.0, 0.35), UnitVec3::Y, 0.275),
            4.0,
        );
        let domain = FaceDomain::CylinderBand { surf: leg, length: 4.0 };

        let rim = |v: f64, tag: u32| FaceChain {
            pts: (0..=64).map(|k| (-PI + TAU * k as f64 / 64.0, v)).collect(),
            tag,
        };
        let loops = closed_loops_between(&a, &b, &TraceOptions::default());
        let mut window: Vec<(f64, f64)> = loops[0].points.iter().map(|&p| domain.uv(p)).collect();
        window.push(window[0]);
        let chains = vec![rim(0.0, 1), rim(4.0, 2), FaceChain { pts: window, tag: 10 }];

        // KEEP both sides so shared edges have two users (the union keeps the
        // band; the window is the neighbour's material — here we assemble the
        // whole partition to verify edge sharing topology).
        let cells = arrange_and_classify(&domain, &chains, &b, 1e-3);
        assert!(cells.len() >= 2);

        let mut brep = BRep::new();
        let mut asm = Assembler::new(&mut brep, 1e-7);
        for c in &cells {
            asm.emit_cell(&domain, FaceGeom::Cylinder(leg), FaceNormal::Same, c);
        }
        let faces = asm.faces().to_vec();
        assert_eq!(faces.len(), cells.len());

        // count edge uses across all emitted faces
        let mut uses: HashMap<EdgeId, usize> = HashMap::new();
        for &fid in &faces {
            let f = &brep.faces[fid];
            let mut loops = vec![f.outer_loop];
            loops.extend(f.inner_loops.iter().copied());
            for lid in loops {
                let start = brep.loops[lid].start;
                let mut c = start;
                loop {
                    let ce = &brep.coedges[c];
                    *uses.entry(ce.edge).or_insert(0) += 1;
                    c = ce.next;
                    if c == start {
                        break;
                    }
                }
            }
        }
        // The window boundary is shared by the band cell and the window
        // cell → those edges are used exactly twice.  The outer rims (v=0,
        // v=4) have no neighbouring face in this isolated band, so they are
        // single-use here (in a full solid the adjacent cylinder section
        // shares them).  The manifold invariant is: NOTHING exceeds 2.
        let twice = uses.values().filter(|&&n| n == 2).count();
        assert!(twice > 0, "window boundary edges must be shared: {uses:?}");
        assert!(
            uses.values().all(|&n| n <= 2),
            "no edge over-shared (manifold): {uses:?}"
        );
        // Watertight test for an isolated band: every SINGLE-USE edge must
        // lie on an open rim (v≈0 or v≈4) — the only true boundary.  Window
        // and seam edges are all shared (used twice).  An interior edge used
        // once would mean an open shell.
        for (&eid, &n) in &uses {
            if n != 1 {
                continue;
            }
            let e = &brep.edges[eid];
            let v0 = domain.uv(brep.vertices[e.v_start].point).1;
            let v1 = domain.uv(brep.vertices[e.v_end].point).1;
            let on_rim = |v: f64| v.abs() < 1e-6 || (v - 4.0).abs() < 1e-6;
            assert!(
                on_rim(v0) && on_rim(v1),
                "single-use edge off the open rim → open shell: v=({v0:.3},{v1:.3})"
            );
        }
    }

    /// Two faces sharing one straight edge: the second user gets the opposite
    /// sense, so the edge is manifold.
    #[test]
    fn shared_edge_opposite_senses() {
        let mut brep = BRep::new();
        let mut asm = Assembler::new(&mut brep, 1e-7);
        let poly = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let (e1, s1) = asm.shared_edge(&poly);
        let rev: Vec<Point3> = poly.iter().rev().copied().collect();
        let (e2, s2) = asm.shared_edge(&rev);
        assert_eq!(e1, e2, "same geometry → one edge");
        assert_eq!(s1, CoEdgeSense::Same);
        assert_eq!(s2, CoEdgeSense::Opposite, "second user flips");
    }
}
