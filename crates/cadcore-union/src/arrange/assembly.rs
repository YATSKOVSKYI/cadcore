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

/// Outward (radial) direction of a cylinder surface at `p` — from the axis to
/// `p`, in the plane through `p` perpendicular to the axis.  Not normalised;
/// only its sign relative to a loop's area vector is used.
fn cyl_outward(surf: &CylSurf, p: Point3) -> cadcore_math::Vec3 {
    let w = p - surf.frame.origin;
    let ax = surf.axis().dot_vec(w);
    w - surf.axis().as_vec() * ax
}

/// Outward direction of a torus (tube) surface at `p` — from the nearest point
/// on the tube centre-line (the major circle) to `p`.  Not normalised.
fn torus_outward(surf: &TorusSurf, p: Point3) -> cadcore_math::Vec3 {
    let w = p - surf.frame.origin;
    let h = surf.frame.z.dot_vec(w);
    let planar = w - surf.frame.z.as_vec() * h;
    let a = planar.length();
    if a < 1e-12 {
        return surf.frame.z.as_vec();
    }
    let centre_pt = surf.frame.origin + planar * (surf.major_radius / a);
    p - centre_pt
}

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

/// An analytic shape recovered from a lifted 3-D polyline edge.
///
/// SSI curves and trimmed boundary arcs arrive as tolerance-sampled polylines.
/// Emitting them verbatim makes the STEP writer produce degree-1 chord
/// B-splines whose deviation from the true surface (the chord sagitta, up to
/// ~0.2 mm here) inflates the OCC edge tolerance — which strict importers
/// (Parasolid / SolidWorks) read as self-intersections.  Recovering the exact
/// `Line` / `Circle` collapses that tolerance to ~1e-7.
enum EdgeShape {
    /// Collinear points → exact line segment.
    Line,
    /// Co-circular points → exact circular arc on this circle.
    Arc(Circle3),
    /// Genuinely free-form (transcendental SSI curve) — keep as a polyline.
    Free,
}

#[inline]
fn dot3(a: cadcore_math::Vec3, b: cadcore_math::Vec3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

/// Circumcircle of three 3-D points: `(centre, axis, radius)`.  `None` when the
/// points are (near-)collinear.
fn circumcircle(a: Point3, b: Point3, c: Point3) -> Option<(Point3, cadcore_math::UnitVec3, f64)> {
    let ab = b - a;
    let ac = c - a;
    let n = ab.cross(ac);
    let n2 = dot3(n, n);
    if n2 < 1e-24 {
        return None; // collinear — no finite circle
    }
    // Vector from `a` to the circumcentre (Ericson, Real-Time Collision Detection).
    let to_centre =
        (n.cross(ab) * dot3(ac, ac) + ac.cross(n) * dot3(ab, ab)) * (1.0 / (2.0 * n2));
    let centre = a + to_centre;
    let axis = cadcore_math::UnitVec3::try_from_vec(n)?;
    let r = to_centre.length();
    Some((centre, axis, r))
}

/// Classify a lifted 3-D polyline as an exact `Line` or `Arc` when its points
/// lie on one within `ftol`; otherwise `Free`.  `ftol` is tight so genuine
/// transcendental SSI curves (which deviate macroscopically from any
/// line/circle) are never mis-recovered.
fn classify_edge(poly: &[Point3], ftol: f64) -> EdgeShape {
    let n = poly.len();
    if n < 2 {
        return EdgeShape::Free;
    }
    let p0 = poly[0];
    let pe = poly[n - 1];
    let chord = pe - p0;
    let chord_len = chord.length();

    // ── Line: every point within ftol of the chord p0→pe (open edge only;
    //    a closed loop has p0≈pe and cannot be a line). ──
    if chord_len > 1e-9 {
        if let Some(dir) = chord.try_normalize() {
            let max_perp = poly
                .iter()
                .map(|&p| {
                    let w = p - p0;
                    let along = dot3(dir, w);
                    (w - dir * along).length()
                })
                .fold(0.0_f64, f64::max);
            if max_perp < ftol {
                return EdgeShape::Line;
            }
        }
    }

    // ── Circle / arc: fit through three spread points, verify all. ──
    if n >= 3 {
        let (a, b, c) = (poly[0], poly[n / 2], poly[n - 1]);
        if let Some((centre, axis, r)) = circumcircle(a, b, c) {
            // Reject implausibly flat fits (huge radius from near-collinear
            // points) — those are better served as a polyline / spline.
            if r.is_finite() && r > 1e-9 && r < 1e6 {
                let ok = poly.iter().all(|&p| {
                    let radial = ((p - centre).length() - r).abs();
                    let planar = dot3(axis.as_vec(), p - centre).abs();
                    radial < ftol && planar < ftol
                });
                if ok {
                    return EdgeShape::Arc(Circle3::new(centre, axis, r));
                }
            }
        }
    }

    EdgeShape::Free
}

/// The analytic curve a run of same-`tag` segments was recovered to.  Every
/// segment of the run is then emitted as a trimmed piece of this curve, keeping
/// the original per-segment topology (and thus the watertight invariant) while
/// upgrading degree-1 chord B-splines to exact LINE / CIRCLE geometry.
#[derive(Clone, Copy)]
enum RunCurve {
    /// Run lies on this circle (full circle or arc) — segments become arcs.
    Circle(Circle3),
    /// Run is straight — segments become line segments.
    Line,
    /// Run is genuinely free-form — segments stay polylines.
    Free,
}

/// Angle of `p` in `c`'s frame (radians).
#[inline]
fn circle_angle(c: &Circle3, p: Point3) -> f64 {
    let w = p - c.frame.origin;
    dot3(c.frame.y.as_vec(), w).atan2(dot3(c.frame.x.as_vec(), w))
}

/// For each step, the analytic curve of the same-`tag` run it belongs to.
/// On a torus carrier, a run at (near-)constant φ is an arc of the MAJOR circle
/// at that φ (the φ-seam / equator of the band); a run at constant θ is a minor
/// circle.  Both are built EXACTLY from the torus geometry instead of fitting —
/// which is the only way to recover the analytic arc from a 2-point seam segment
/// (a point fit needs ≥3 spread points).  Returns `None` when the run varies in
/// both θ and φ (a genuine transcendental imprint curve).
fn torus_run_circle(domain: &FaceDomain, run_pts: &[Point3], ftol: f64) -> Option<Circle3> {
    use std::f64::consts::{PI, TAU};
    let surf = match domain {
        FaceDomain::TorusPatch { surf, .. } | FaceDomain::TorusBand { surf } => surf,
        _ => return None,
    };
    // The constructed circle must actually pass through EVERY run point — a short
    // transcendental imprint segment can look locally constant-φ/θ yet not lie on
    // the major/minor circle.  Verifying rejects those (they fall through to the
    // fine-sampled polyline path).
    let verify = |c: &Circle3| -> bool {
        run_pts.iter().all(|&p| {
            let w = p - c.frame.origin;
            let planar = c.frame.z.dot_vec(w).abs();
            let radial = (w - c.frame.z.as_vec() * c.frame.z.dot_vec(w)).length() - c.radius;
            planar < ftol && radial.abs() < ftol
        })
    };
    let uvs: Vec<(f64, f64)> = run_pts.iter().map(|&p| domain.uv(p)).collect();
    let circ_mean = |angs: &[f64]| -> f64 {
        let (mut sx, mut sy) = (0.0, 0.0);
        for &a in angs {
            sx += a.cos();
            sy += a.sin();
        }
        sy.atan2(sx)
    };
    let spread = |angs: &[f64], mean: f64| -> f64 {
        angs.iter()
            .map(|&a| {
                let mut d = a - mean;
                while d > PI {
                    d -= TAU;
                }
                while d < -PI {
                    d += TAU;
                }
                d.abs()
            })
            .fold(0.0_f64, f64::max)
    };
    let thetas: Vec<f64> = uvs.iter().map(|p| p.0).collect();
    let phis: Vec<f64> = uvs.iter().map(|p| p.1).collect();
    let th_mean = circ_mean(&thetas);
    let ph_mean = circ_mean(&phis);
    let th_spread = spread(&thetas, th_mean);
    let ph_spread = spread(&phis, ph_mean);
    let flat = 1e-4;
    if ph_spread < flat && th_spread > flat {
        // constant φ → major circle (radius R + r·cosφ, centred on the axis at
        // height r·sinφ, normal = torus axis).
        let radius = surf.major_radius + surf.minor_radius * ph_mean.cos();
        if radius <= 1e-9 {
            return None;
        }
        let centre = surf.frame.origin + surf.frame.z * (surf.minor_radius * ph_mean.sin());
        let c = Circle3::new(centre, surf.frame.z, radius);
        return verify(&c).then_some(c);
    }
    if th_spread < flat && ph_spread > flat {
        // constant θ → minor circle (radius r, centred on the ring, normal =
        // ring tangent = ring_dir × axis).
        let ring_dir = surf.frame.x * th_mean.cos() + surf.frame.y * th_mean.sin();
        let centre = surf.frame.origin + ring_dir * surf.major_radius;
        let axis = cadcore_math::UnitVec3::try_from_vec(ring_dir.cross(surf.frame.z.as_vec()))?;
        let c = Circle3::new(centre, axis, surf.minor_radius);
        return verify(&c).then_some(c);
    }
    None
}

/// Order of step indices for run-grouping: rotated so the loop starts at a tag
/// boundary, ensuring a run that wraps past the loop's start index is grouped
/// as ONE run.  Returns `0..n` when every step shares one tag (a single closed
/// curve filling the whole loop).
fn tag_run_order(steps: &[LoopStep]) -> Vec<usize> {
    let n = steps.len();
    if n == 0 {
        return Vec::new();
    }
    match (0..n).find(|&i| steps[i].tag != steps[(i + n - 1) % n].tag) {
        Some(s) => (0..n).map(|k| (s + k) % n).collect(),
        None => (0..n).collect(),
    }
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
    /// Topological weld group (crossing identity for window edges; 0 for the
    /// globally-welded cell/boundary edges).  Two coincident polyline pieces
    /// weld only within the same group — a leg crossed by two coincident legs
    /// (a self-overlapping serpentine) then gets two DISTINCT twice-used window
    /// edges instead of one four-used edge.
    group: u32,
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
    /// (weld group, circle key) → (edge id, first user's start vertex) — analytic
    /// rim/cap circles, shared SEAMLESSLY between a cylinder rim and its cap /
    /// elbow WITHIN one weld group.
    ///
    /// A full junction/cap circle is a *topological* junction of ONE filament
    /// (leg↔elbow, leg↔cap), shared by exactly two faces.  Two DISTINCT
    /// filaments that merely touch tangentially (e.g. equal-radius tubes whose
    /// centre-lines sit 2r apart, or two U-turns that coincide) have
    /// geometrically-identical junction circles — welding those by geometry
    /// alone collapses them into ONE edge used by four faces → non-manifold.
    /// Keying by `weld_group` (set per filament, see [`Self::set_weld_group`])
    /// keeps each filament's junctions separate, so tangent filaments produce
    /// two independent twice-used edges the way a manifold shell requires.
    /// Crossing *windows* (which MUST weld across filaments) go through the
    /// polyline `edges` cache instead and stay global.
    circles: HashMap<(u32, CircleKey), (EdgeId, VertexId)>,
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
        // Neighbour-checking weld: a point within `tol` of an existing vertex
        // can round into an ADJACENT grid cell, so scan the 27 cells around `k`
        // and reuse any vertex within `tol`.  Without this, two faces' versions
        // of a shared point (≈ a few 1e-4 apart at a curve's seam) straddle a
        // cell boundary and fail to merge → unshared slivers.
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    if let Some(&v) = self.verts.get(&PtKey(k.0 + dx, k.1 + dy, k.2 + dz)) {
                        if (self.brep.vertices[v].point - p).length() < self.tol {
                            return v;
                        }
                    }
                }
            }
        }
        let v = self.brep.add_vertex(cadcore_topo::Vertex { point: p });
        self.verts.insert(k, v);
        v
    }

    /// Snap a point to its welded vertex's canonical position, so two faces
    /// that compute a shared point slightly differently emit IDENTICAL
    /// polylines (and therefore identical edge keys → one shared edge).
    fn snap(&mut self, p: Point3) -> Point3 {
        let v = self.vertex(p);
        self.brep.vertices[v].point
    }

    /// Get-or-create the shared edge for a 3-D polyline; returns the co-edge
    /// sense this user must traverse it with (`Same` for the first user,
    /// `Opposite` for the second).
    fn shared_edge(&mut self, poly3: &[Point3], group: u32) -> (EdgeId, CoEdgeSense) {
        let geom = EdgeGeom::Polyline(poly3.to_vec());
        self.shared_edge_geom(poly3, geom, 0.0, 1.0, group)
    }

    /// Emit one segment edge, trimmed onto its run's recovered analytic curve.
    /// `poly3` is the segment's own (usually 2) lifted points; `rc` is the
    /// curve fitted from the whole same-`tag` run it belongs to.
    fn emit_step_edge(&mut self, poly3: &[Point3], rc: RunCurve) -> (EdgeId, CoEdgeSense) {
        let from = poly3[0];
        let to = *poly3.last().unwrap();
        match rc {
            RunCurve::Circle(c) => {
                let t0 = circle_angle(&c, from);
                let mut t1 = circle_angle(&c, to);
                // a single segment subtends < π — unwrap the endpoint onto the
                // short arc so `t_start..t_end` describes the real sweep.
                let pi = std::f64::consts::PI;
                while t1 - t0 > pi {
                    t1 -= 2.0 * pi;
                }
                while t1 - t0 < -pi {
                    t1 += 2.0 * pi;
                }
                self.shared_edge_geom(poly3, EdgeGeom::Circle(c), t0, t1, 0)
            }
            RunCurve::Line => match (to - from).try_normalize() {
                Some(dir) => {
                    let len = (to - from).length();
                    let udir = cadcore_math::UnitVec3::new_unchecked(dir);
                    self.shared_edge_geom(poly3, EdgeGeom::Line(cadcore_geom::Line3::new(from, udir)), 0.0, len, 0)
                }
                None => self.shared_edge(poly3, 0),
            },
            RunCurve::Free => self.shared_edge(poly3, 0),
        }
    }

    /// Like [`shared_edge`] but stores an explicit `geom` and parameter range —
    /// used to weld an analytic arc / line recovered from the lifted polyline.
    /// Keying stays on the 3-D endpoints + arc midpoint, so the analytic edge
    /// shares ONE `Edge` with a neighbour that presents the same byte-identical
    /// points, regardless of which `EdgeGeom` form each chose.
    fn shared_edge_geom(
        &mut self,
        poly3: &[Point3],
        geom: EdgeGeom,
        t_start: f64,
        t_end: f64,
        group: u32,
    ) -> (EdgeId, CoEdgeSense) {
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
            group,
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
            geom,
            v_start,
            v_end,
            t_start,
            t_end,
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
        // Lift every step's uv polyline to 3-D, snapping each point to its
        // welded vertex so a shared edge is byte-identical on both faces, then
        // drop consecutive duplicates the snap creates.  One entry per step,
        // aligned with `steps`.
        let mut step_poly: Vec<Vec<Point3>> = Vec::with_capacity(steps.len());
        for step in steps {
            let mut poly3: Vec<Point3> = Vec::with_capacity(step.pts.len());
            for &(u, v) in &step.pts {
                let p = self.snap(domain.lift(u, v));
                if poly3.last().map_or(true, |&q| (q - p).length() > 1e-12) {
                    poly3.push(p);
                }
            }
            step_poly.push(poly3);
        }

        // The arrangement shatters every input chain into one `LoopStep` per
        // 2-point segment, so a single source curve (a minor circle, a straight
        // boundary, an SSI arc) arrives as a RUN of consecutive steps sharing
        // one `tag`.  We fit each run as a whole to recover the exact line /
        // circle it came from, then keep the ORIGINAL per-segment topology
        // (one edge per step) — only the stored geometry changes.  This is what
        // preserves the watertight invariant: the segment endpoints, the
        // shared-edge keying and the co-edge structure are byte-identical to the
        // proven polyline path; the neighbouring face fits the same run to the
        // same circle (the writer reconstructs the arc from centre/axis/radius +
        // the two vertices, independent of the frame's arbitrary X), so both
        // emit one shared analytic EDGE_CURVE.  A 2-point arc segment is exact
        // because its circle is known from the parent run, not inferred from the
        // two points.
        let ftol = (50.0 * self.tol).max(1e-5);
        let order = tag_run_order(steps);
        let mut coedges: Vec<CoEdgeId> = Vec::with_capacity(steps.len());
        let mut push_coedge = |brep: &mut BRep, coedges: &mut Vec<CoEdgeId>, edge, sense| {
            let ce = brep.add_coedge(CoEdge {
                edge,
                sense,
                next: CoEdgeId::default(),
                prev: CoEdgeId::default(),
                loop_id: lp,
            });
            coedges.push(ce);
        };
        let mut gi = 0;
        while gi < order.len() {
            let tag = steps[order[gi]].tag;
            let seg_lo = gi;
            let mut run_pts: Vec<Point3> = Vec::new();
            while gi < order.len() && steps[order[gi]].tag == tag {
                for &p in &step_poly[order[gi]] {
                    if run_pts.last().map_or(true, |&q| (q - p).length() > 1e-12) {
                        run_pts.push(p);
                    }
                }
                gi += 1;
            }
            let seg_hi = gi;
            if run_pts.len() < 2 {
                continue;
            }
            let rc = if let Some(c) = torus_run_circle(domain, &run_pts, ftol) {
                RunCurve::Circle(c)
            } else if let Some(c) = fit_circle(&run_pts) {
                RunCurve::Circle(c)
            } else {
                match classify_edge(&run_pts, ftol) {
                    EdgeShape::Arc(c) => RunCurve::Circle(c),
                    EdgeShape::Line => RunCurve::Line,
                    EdgeShape::Free => RunCurve::Free,
                }
            };
            for k in seg_lo..seg_hi {
                let poly3 = &step_poly[order[k]];
                if poly3.len() < 2 {
                    continue;
                }
                let (edge, sense) = self.emit_step_edge(poly3, rc);
                push_coedge(self.brep, &mut coedges, edge, sense);
            }
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

    /// Re-emit an EXISTING face unchanged through this assembler.
    ///
    /// The half-space cut uses this for every face the plane does not touch.
    /// Running such a face back through the arrangement is not neutral: a fused
    /// periodic band is *seamless* (its boundary is two rim circles, no
    /// parametric seam), and the periodic DCEL re-introduces a seam — splitting
    /// each crossing-window loop at a seam the partner face does not share, so
    /// the two halves never weld.  Copying keeps the source topology.
    ///
    /// Every edge still goes through the shared-edge cache, at the SAME
    /// per-segment granularity [`Self::build_loop`] uses, so a copied face and a
    /// re-arranged one weld along the boundary between them.
    ///
    /// `joints` are extra points that must become vertices wherever they lie on
    /// one of this face's edges: they are where the neighbouring ARRANGED face
    /// was split, so the copy has to be split there too or the two disagree.
    pub fn emit_face_copy(
        &mut self,
        src: &BRep,
        fid: cadcore_topo::FaceId,
        joints: &[Point3],
    ) -> Option<cadcore_topo::FaceId> {
        let sf = src.faces.get(fid)?;
        let (geom, normal, extent) = (sf.geom.clone(), sf.normal, sf.extent.clone());
        let inner_ids: Vec<LoopId> = sf.inner_loops.clone();
        let outer_id = sf.outer_loop;
        let face = self.brep.add_face(Face {
            geom,
            normal,
            outer_loop: LoopId::default(),
            inner_loops: Vec::new(),
            shell: Default::default(),
            extent,
        });
        let outer = self.copy_loop(src, outer_id, face, joints)?;
        let inners: Vec<LoopId> = inner_ids
            .into_iter()
            .filter_map(|l| self.copy_loop(src, l, face, joints))
            .collect();
        self.brep.faces[face].outer_loop = outer;
        self.brep.faces[face].inner_loops = inners;
        self.faces.push(face);
        Some(face)
    }

    /// Copy one loop of `src` into this assembler's B-Rep, welding every edge.
    fn copy_loop(
        &mut self,
        src: &BRep,
        lid: LoopId,
        face: cadcore_topo::FaceId,
        joints: &[Point3],
    ) -> Option<LoopId> {
        let slp = src.loops.get(lid)?;
        let start = slp.start;
        src.coedges.get(start)?;
        let lp = self.brep.add_loop(Loop {
            start: CoEdgeId::default(),
            face,
        });
        let mut coedges: Vec<CoEdgeId> = Vec::new();
        let mut c = start;
        let mut guard = 0usize;
        loop {
            let ce = &src.coedges[c];
            let (edge_id, sense_src, next) = (ce.edge, ce.sense, ce.next);
            let rc = match &src.edges[edge_id].geom {
                EdgeGeom::Circle(c3) => RunCurve::Circle(*c3),
                EdgeGeom::Line(_) => RunCurve::Line,
                _ => RunCurve::Free,
            };
            let mut pts = crate::boolean::aabb::sample_edge(src, edge_id);
            if sense_src == CoEdgeSense::Opposite {
                pts.reverse();
            }
            crate::boolean::halfspace::insert_joints(&mut pts, joints);
            let mut poly: Vec<Point3> = Vec::with_capacity(pts.len());
            for p in pts {
                let q = self.snap(p);
                if poly.last().map_or(true, |&r| (r - q).length() > 1e-12) {
                    poly.push(q);
                }
            }
            for w in poly.windows(2) {
                let (edge, sense) = self.emit_step_edge(w, rc);
                let ceid = self.brep.add_coedge(CoEdge {
                    edge,
                    sense,
                    next: CoEdgeId::default(),
                    prev: CoEdgeId::default(),
                    loop_id: lp,
                });
                coedges.push(ceid);
            }
            c = next;
            guard += 1;
            if c == start || guard > 100_000 {
                break;
            }
        }
        if coedges.is_empty() {
            return None;
        }
        let n = coedges.len();
        for i in 0..n {
            let cur = coedges[i];
            self.brep.coedges[cur].next = coedges[(i + 1) % n];
            self.brep.coedges[cur].prev = coedges[(i + n - 1) % n];
        }
        self.brep.loops[lp].start = coedges[0];
        Some(lp)
    }

    /// Get-or-create a shared analytic full-circle edge; returns the co-edge
    /// sense this user must traverse it with.  Welds a cylinder rim to its
    /// cap / elbow with NO seam (the periodic surface closes implicitly).
    ///
    /// `group` is the TOPOLOGICAL junction identity: two coincident circles
    /// weld into one edge only if they carry the same group.  The mating faces
    /// of ONE junction (leg rim ↔ its elbow / cap) pass the same group; two
    /// DISTINCT junctions that merely coincide in space (a self-touching
    /// serpentine, or two tangent filaments) pass different groups and stay
    /// separate twice-used edges — the manifold invariant.  Callers that don't
    /// care (single features, no coincidence) pass `0`.
    fn shared_circle(&mut self, c: Circle3, group: u32) -> (EdgeId, CoEdgeSense) {
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
        let gkey = (group, key);
        let start_pt = c.point_at(0.0);
        if let Some(&(eid, first_start)) = self.circles.get(&gkey) {
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
        self.circles.insert(gkey, (eid, v));
        (eid, CoEdgeSense::Same)
    }

    /// Build a loop holding exactly ONE closed analytic circle edge, welded
    /// within the topological junction `group` (see [`Self::shared_circle`]).
    ///
    /// `outward`: when `Some(n)`, the loop is wound CCW about `n` (flip the
    /// shared co-edge if the circle's own axis opposes it).  A planar DISK cap
    /// MUST be CCW about its outward normal — SolidWorks/SpaceClaim IGNORE the
    /// `FACE_OUTER_BOUND` orientation flag and read the physical winding, so a
    /// clockwise outer circle is taken as "keep the outside" and the cap face is
    /// dropped (torn caps).  Periodic cylinder/torus rims pass `None` (their
    /// boundary wire winds around the surface period, so it has no area role).
    fn circle_loop(
        &mut self,
        c: Circle3,
        face: cadcore_topo::FaceId,
        group: u32,
        outward: Option<cadcore_math::UnitVec3>,
    ) -> LoopId {
        let lp = self.brep.add_loop(Loop { start: CoEdgeId::default(), face });
        let (edge, mut sense) = self.shared_circle(c, group);
        if let Some(n) = outward {
            if c.frame.z.dot(n) < 0.0 {
                sense = match sense {
                    CoEdgeSense::Same => CoEdgeSense::Opposite,
                    CoEdgeSense::Opposite => CoEdgeSense::Same,
                };
            }
        }
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
        windows: &[(u32, Vec<Point3>)],
        g_lo: u32,
        g_hi: u32,
    ) -> cadcore_topo::FaceId {
        let face = self.brep.add_face(Face {
            geom: FaceGeom::Cylinder(surf),
            normal: FaceNormal::Same,
            outer_loop: LoopId::default(),
            inner_loops: Vec::new(),
            shell: Default::default(),
            extent: FaceExtent::Trimmed,
        });
        // outer = lo rim circle; inner = hi rim circle + window holes.  `g_lo`,
        // `g_hi` = the topological junction groups of the two rims; each window
        // carries its crossing id so coincident windows stay distinct edges.
        let outer = self.circle_loop(rim_lo, face, g_lo, None);
        let mut inner = vec![self.circle_loop(rim_hi, face, g_hi, None)];
        for (wg, w) in windows {
            inner.push(self.window_hole_loop(w, |p| cyl_outward(&surf, p), face, *wg));
        }
        self.brep.faces[face].outer_loop = outer;
        self.brep.faces[face].inner_loops = inner;
        self.faces.push(face);
        face
    }

    /// Build a window hole loop on a carrier surface.  The loop is wound CW
    /// about the surface's OUTWARD normal (the STEP hole convention) and split
    /// into a FEW shared edges — a single closed edge cannot carry an opposite
    /// sense between the two faces sharing it (start == end vertex), so we split
    /// so the two members' pieces weld with opposite senses → a manifold edge.
    ///
    /// `outward_at` returns the surface's outward normal at a point (radial from
    /// the axis for a cylinder, radial from the tube centre-line for a torus);
    /// it decides the traversal direction so BOTH faces sharing this window wind
    /// it CW about their own outward normal and thus weld each shared edge with
    /// opposite senses regardless of carrier type (cylinder↔cylinder,
    /// cylinder↔torus).
    fn window_hole_loop(
        &mut self,
        pts3: &[Point3],
        outward_at: impl Fn(Point3) -> cadcore_math::Vec3,
        face: cadcore_topo::FaceId,
        group: u32,
    ) -> LoopId {
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
        let out = outward_at(cen);
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
            let (edge, sense) = self.shared_edge(&seg, group);
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
        g_lo: u32,
        g_hi: u32,
    ) -> cadcore_topo::FaceId {
        self.emit_torus_fillet_windows(surf, junction_lo, junction_hi, &[], g_lo, g_hi)
    }

    /// Emit an elbow torus-fillet face bounded by its two junction minor-circles
    /// PLUS the crossing `windows` as interior holes — the torus counterpart of
    /// [`Self::emit_cylinder_face`].  Each window is a closed 3-D loop shared
    /// with the crossing tube's matching window (single source of truth → the
    /// two faces weld it with opposite senses → manifold).  Used where an elbow
    /// (a U-turn fillet) overlaps a neighbouring leg — the dense-fill case the
    /// leg↔leg pass alone cannot resolve.
    ///
    /// The extent stays `TorusFillet`: `face_has_real_loops` is true (the two
    /// minor circles are real), so the writer emits every inner loop (end circle
    /// + windows) as `FACE_BOUND`s, and the θ-band is read from the extent
    /// circles — the window polylines never pollute it (they are not
    /// minor-radius circles).
    pub fn emit_torus_fillet_windows(
        &mut self,
        surf: TorusSurf,
        junction_lo: Circle3,
        junction_hi: Circle3,
        windows: &[(u32, Vec<Point3>)],
        g_lo: u32,
        g_hi: u32,
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
        let outer = self.circle_loop(junction_lo, face, g_lo, None);
        let mut inner = vec![self.circle_loop(junction_hi, face, g_hi, None)];
        for (wg, w) in windows {
            inner.push(self.window_hole_loop(w, |p| torus_outward(&surf, p), face, *wg));
        }
        self.brep.faces[face].outer_loop = outer;
        self.brep.faces[face].inner_loops = inner;
        self.faces.push(face);
        face
    }

    /// Emit a flat end-cap disk bounded by a single shared circle (welds to
    /// the cylinder rim with no seam).
    pub fn emit_disk_cap(&mut self, plane: Plane3, rim: Circle3, group: u32) -> cadcore_topo::FaceId {
        let face = self.brep.add_face(Face {
            geom: FaceGeom::Plane(plane),
            normal: FaceNormal::Same,
            outer_loop: LoopId::default(),
            inner_loops: Vec::new(),
            shell: Default::default(),
            extent: FaceExtent::Disk { radius: rim.radius },
        });
        // A planar disk cap MUST wind CCW about its outward (plane) normal or
        // SpaceClaim/SolidWorks read it as "keep the outside" and drop the face
        // (torn caps).  Pass the plane normal so `circle_loop` flips the shared
        // co-edge when the rim circle's own axis opposes the outward normal.
        let outer = self.circle_loop(rim, face, group, Some(plane.normal()));
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
        let (e1, s1) = asm.shared_edge(&poly, 0);
        let rev: Vec<Point3> = poly.iter().rev().copied().collect();
        let (e2, s2) = asm.shared_edge(&rev, 0);
        assert_eq!(e1, e2, "same geometry → one edge");
        assert_eq!(s1, CoEdgeSense::Same);
        assert_eq!(s2, CoEdgeSense::Opposite, "second user flips");
    }
}
