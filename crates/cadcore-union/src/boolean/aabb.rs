//! Axis-aligned bounding boxes and the broad phase.
//!
//! The general union starts by finding which faces of A *could* intersect
//! which faces of B.  A cheap AABB overlap test culls the O(n·m) face pairs
//! down to the handful that actually meet, so the (expensive) surface×surface
//! intersector only runs where it must.

use cadcore_math::Point3;
use cadcore_topo::{BRep, FaceBoundary, FaceExtent, FaceGeom, FaceId};

/// An axis-aligned bounding box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// Minimum corner.
    pub min: Point3,
    /// Maximum corner.
    pub max: Point3,
}

impl Aabb {
    /// An empty box (inverted bounds) that absorbs points via [`Aabb::expand`].
    pub fn empty() -> Self {
        Aabb {
            min: Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            max: Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        }
    }

    /// Grow the box to include `p`.
    pub fn expand(&mut self, p: Point3) {
        self.min.x = self.min.x.min(p.x);
        self.min.y = self.min.y.min(p.y);
        self.min.z = self.min.z.min(p.z);
        self.max.x = self.max.x.max(p.x);
        self.max.y = self.max.y.max(p.y);
        self.max.z = self.max.z.max(p.z);
    }

    /// Grow the box uniformly by `pad` on every side.
    pub fn pad(mut self, pad: f64) -> Self {
        self.min = self.min - cadcore_math::Vec3::new(pad, pad, pad);
        self.max = self.max + cadcore_math::Vec3::new(pad, pad, pad);
        self
    }

    /// Whether two boxes overlap (closed test; touching counts).
    pub fn overlaps(&self, o: &Aabb) -> bool {
        self.min.x <= o.max.x
            && self.max.x >= o.min.x
            && self.min.y <= o.max.y
            && self.max.y >= o.min.y
            && self.min.z <= o.max.z
            && self.max.z >= o.min.z
    }

    /// Centre of the box.
    pub fn centre(&self) -> Point3 {
        Point3::new(
            0.5 * (self.min.x + self.max.x),
            0.5 * (self.min.y + self.max.y),
            0.5 * (self.min.z + self.max.z),
        )
    }

    /// Longest side length (used to size probe rays / fallback extents).
    pub fn extent(&self) -> f64 {
        (self.max.x - self.min.x)
            .max(self.max.y - self.min.y)
            .max(self.max.z - self.min.z)
    }
}

/// Bounding box of one face.
///
/// Uses the face's real loops when present (the authoritative trimmed
/// boundary); otherwise falls back to the [`FaceExtent`] template (legacy
/// faces whose loops are vestigial) so we never under-bound a face.
pub fn face_aabb(brep: &BRep, face: FaceId) -> Aabb {
    let mut bb = Aabb::empty();
    let f = &brep.faces[face];

    // A torus arc's loops are only its END circles — they do NOT bound the
    // SWEEP of the tube between them, so a crosser meeting the arc mid-span
    // would be missed by the broad phase.  Sample the whole arc surface.
    if let FaceGeom::Torus(t) = &f.geom {
        if let Some((lo, hi)) = torus_theta_range(brep, face, t) {
            for i in 0..=24 {
                let theta = lo + (hi - lo) * (i as f64 / 24.0);
                for j in 0..16 {
                    let phi = std::f64::consts::TAU * (j as f64 / 16.0);
                    bb.expand(t.point_at(theta, phi));
                }
            }
            return bb;
        }
    }

    // 1) real loop vertices, if any.
    let mut had_loop = false;
    let mut loops = Vec::new();
    if brep.loops.get(f.outer_loop).is_some() {
        loops.push(f.outer_loop);
    }
    loops.extend(f.inner_loops.iter().copied());
    for lid in loops {
        let Some(lp) = brep.loops.get(lid) else { continue };
        let start = lp.start;
        if brep.coedges.get(start).is_none() {
            continue;
        }
        let mut c = start;
        loop {
            let ce = &brep.coedges[c];
            for p in sample_edge(brep, ce.edge) {
                bb.expand(p);
                had_loop = true;
            }
            c = ce.next;
            if c == start {
                break;
            }
        }
    }
    if had_loop {
        return bb;
    }

    // 2) fall back to the analytic extent template.
    extent_aabb(&f.geom, &f.extent, &mut bb);
    bb
}

/// The θ-arc range of a torus face, from its `TorusFillet` template circles or
/// from its loop vertices (a Trimmed torus).
fn torus_theta_range(brep: &BRep, face: FaceId, t: &cadcore_geom::TorusSurf) -> Option<(f64, f64)> {
    let theta_of = |p: Point3| {
        let w = p - t.frame.origin;
        t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w))
    };
    let f = &brep.faces[face];
    if let FaceExtent::TorusFillet { start_circle, end_circle } = &f.extent {
        if (start_circle.frame.origin - end_circle.frame.origin).length() < 1e-6 {
            return Some((0.0, std::f64::consts::TAU)); // full torus
        }
        let lo = theta_of(start_circle.frame.origin);
        let mut hi = theta_of(end_circle.frame.origin);
        hi += std::f64::consts::TAU * ((lo - hi) / std::f64::consts::TAU).round();
        return Some(if lo <= hi { (lo, hi) } else { (hi, lo) });
    }
    // Trimmed torus: unwrap loop-vertex θ onto the first sample's branch.
    let mut ids = vec![f.outer_loop];
    ids.extend(f.inner_loops.iter().copied());
    let mut base = None;
    let (mut lo, mut hi) = (f64::MAX, f64::MIN);
    for lid in ids {
        let Some(lp) = brep.loops.get(lid) else { continue };
        if brep.coedges.get(lp.start).is_none() {
            continue;
        }
        let mut c = lp.start;
        loop {
            let ce = &brep.coedges[c];
            for p in sample_edge(brep, ce.edge) {
                let raw = theta_of(p);
                let b = *base.get_or_insert(raw);
                let th = raw + std::f64::consts::TAU * ((b - raw) / std::f64::consts::TAU).round();
                lo = lo.min(th);
                hi = hi.max(th);
            }
            c = ce.next;
            if c == lp.start {
                break;
            }
        }
    }
    (hi > lo).then_some((lo, hi))
}

/// Sample an edge's curve into a few world-space points (endpoints + interior
/// for curved edges).  Shared with the containment loop-fallback path.
pub(crate) fn sample_edge(brep: &BRep, edge: cadcore_topo::EdgeId) -> Vec<Point3> {
    use cadcore_topo::EdgeGeom;
    let e = &brep.edges[edge];
    match &e.geom {
        EdgeGeom::Line(l) => vec![l.point_at(e.t_start), l.point_at(e.t_end)],
        EdgeGeom::Circle(c) => sample_param(e.t_start, e.t_end, 16, |t| c.point_at(t)),
        EdgeGeom::Ellipse(el) => sample_param(e.t_start, e.t_end, 16, |t| el.point_at(t)),
        EdgeGeom::Polyline(pts) => pts.clone(),
    }
}

fn sample_param(t0: f64, t1: f64, n: usize, f: impl Fn(f64) -> Point3) -> Vec<Point3> {
    (0..=n).map(|k| f(t0 + (t1 - t0) * k as f64 / n as f64)).collect()
}

/// Bound a face from its [`FaceExtent`] template (no real loops).
fn extent_aabb(geom: &FaceGeom, extent: &FaceExtent, bb: &mut Aabb) {
    let bound_pts = |b: &FaceBoundary, bb: &mut Aabb| match b {
        FaceBoundary::Circle(c) => {
            for k in 0..24 {
                bb.expand(c.point_at(std::f64::consts::TAU * k as f64 / 24.0));
            }
        }
        FaceBoundary::Ellipse(e) => {
            for k in 0..24 {
                bb.expand(e.point_at(std::f64::consts::TAU * k as f64 / 24.0));
            }
        }
    };
    match extent {
        FaceExtent::Cylinder { start, end, .. } => {
            bound_pts(start, bb);
            bound_pts(end, bb);
        }
        FaceExtent::Disk { radius } | FaceExtent::PartialDisk { radius, .. } => {
            if let FaceGeom::Plane(pl) = geom {
                for k in 0..24 {
                    let a = std::f64::consts::TAU * k as f64 / 24.0;
                    bb.expand(pl.frame.origin + (pl.frame.x * a.cos() + pl.frame.y * a.sin()) * *radius);
                }
            }
        }
        FaceExtent::TorusFillet { start_circle, end_circle } => {
            for c in [start_circle, end_circle] {
                for k in 0..24 {
                    bb.expand(c.point_at(std::f64::consts::TAU * k as f64 / 24.0));
                }
            }
        }
        FaceExtent::PlanarBoundary { boundary } => bound_pts(boundary, bb),
        FaceExtent::Polygon { points } => {
            for &p in points {
                bb.expand(p);
            }
        }
        _ => {
            // None / templates we don't special-case: bound by the carrier
            // surface's reference point at least, so the box is never NaN.
            if let FaceGeom::Plane(pl) = geom {
                bb.expand(pl.frame.origin);
            }
        }
    }
}

/// Candidate intersecting face pairs `(face_a, face_b)` after AABB culling.
///
/// `pad` widens each box by the knit tolerance so grazing contacts aren't
/// missed.
pub fn candidate_pairs(
    a: &BRep,
    a_faces: &[FaceId],
    b: &BRep,
    b_faces: &[FaceId],
    pad: f64,
) -> Vec<(FaceId, FaceId)> {
    let boxes_a: Vec<Aabb> = a_faces.iter().map(|&f| face_aabb(a, f).pad(pad)).collect();
    let boxes_b: Vec<Aabb> = b_faces.iter().map(|&f| face_aabb(b, f).pad(pad)).collect();
    let mut out = Vec::new();
    for (i, fa) in a_faces.iter().enumerate() {
        for (j, fb) in b_faces.iter().enumerate() {
            if boxes_a[i].overlaps(&boxes_b[j]) {
                out.push((*fa, *fb));
            }
        }
    }
    out
}
