//! 2D planar **arrangement** (DCEL) on a face's `(u,v)` parameter domain.
//!
//! This is the core of the face-splitting Boolean: a face's existing boundary
//! loops plus the imprinted surface∩surface intersection curves are inserted
//! as polyline segments into one planar subdivision; the subdivision's cells
//! are then classified in/out against the other solid and the kept cells are
//! stitched back into B-Rep faces.
//!
//! Key property: each cell's outer loop comes out **CCW** and its holes
//! **CW** directly from the half-edge traversal — loop winding and shared
//! edge senses are *derived*, never hand-tuned (see
//! `docs/face-arrangement-boolean.md`).
//!
//! Stage S0: plain rectangular domain, no seam periodicity (S1 adds the
//! periodic `u` glue).

use std::collections::HashMap;

/// One input polyline chain in `(u,v)` space.
///
/// `tag` is an opaque caller id used to map output edges back to their source
/// (an original boundary edge, or an SSI intersection curve).
#[derive(Debug, Clone)]
pub struct ChainInput {
    /// Polyline vertices in `(u,v)`; at least 2 points.
    pub pts: Vec<(f64, f64)>,
    /// Opaque source id (original boundary edge / imprint curve).
    pub tag: u32,
}

/// A directed reference to an arrangement edge inside a cell loop.
#[derive(Debug, Clone, PartialEq)]
pub struct LoopStep {
    /// Source chain tag this edge came from.
    pub tag: u32,
    /// Edge geometry in traversal order (first point = loop direction).
    pub pts: Vec<(f64, f64)>,
    /// `true` when traversal direction matches the input chain direction.
    pub forward: bool,
}

/// One cell (2D face) of the arrangement.
#[derive(Debug, Clone)]
pub struct Cell {
    /// Outer boundary, CCW (positive signed area).
    pub outer: Vec<LoopStep>,
    /// Holes, each CW (negative signed area).
    pub holes: Vec<Vec<LoopStep>>,
    /// Signed area of the outer loop minus holes.
    pub area: f64,
}

impl Cell {
    /// A point strictly inside the cell (centroid-biased; falls back to edge
    /// midpoint offset).  Used for in/out classification against the other
    /// solid.
    #[must_use]
    pub fn interior_point(&self) -> Option<(f64, f64)> {
        let poly = loop_polygon(&self.outer);
        interior_point_of(&poly, &self.holes.iter().map(|h| loop_polygon(h)).collect::<Vec<_>>())
    }
}

/// Numeric options for arrangement construction.
#[derive(Debug, Clone, Copy)]
pub struct ArrangementOptions {
    /// Distance below which two `(u,v)` points are the same vertex.
    pub tol: f64,
}

impl Default for ArrangementOptions {
    fn default() -> Self {
        Self { tol: 1e-9 }
    }
}

// ── internal primitive types ──────────────────────────────────────────────────

type P = (f64, f64);

fn sub(a: P, b: P) -> P {
    (a.0 - b.0, a.1 - b.1)
}
fn cross(a: P, b: P) -> f64 {
    a.0 * b.1 - a.1 * b.0
}
fn dist2(a: P, b: P) -> f64 {
    let d = sub(a, b);
    d.0 * d.0 + d.1 * d.1
}

/// One atomic (already split) segment.
#[derive(Debug, Clone)]
struct Seg {
    a: usize, // vertex index
    b: usize,
    tag: u32,
    forward: bool, // does a→b match the input chain direction
}

/// Segment-segment intersection: returns parameters `(t, s)` on a/b when the
/// open interiors properly intersect or touch within `tol`.
fn seg_seg_params(p1: P, p2: P, q1: P, q2: P, tol: f64) -> Option<(f64, f64)> {
    let r = sub(p2, p1);
    let s = sub(q2, q1);
    let denom = cross(r, s);
    let qp = sub(q1, p1);
    if denom.abs() < 1e-18 {
        return None; // parallel/collinear handled by endpoint snapping
    }
    let t = cross(qp, s) / denom;
    let u = cross(qp, r) / denom;
    let len_r = (r.0 * r.0 + r.1 * r.1).sqrt().max(1e-300);
    let len_s = (s.0 * s.0 + s.1 * s.1).sqrt().max(1e-300);
    let tt = tol / len_r;
    let ts = tol / len_s;
    if t > -tt && t < 1.0 + tt && u > -ts && u < 1.0 + ts {
        Some((t.clamp(0.0, 1.0), u.clamp(0.0, 1.0)))
    } else {
        None
    }
}

/// Vertex pool with tolerance snapping.
struct VertexPool {
    pts: Vec<P>,
    grid: HashMap<(i64, i64), Vec<usize>>,
    tol: f64,
}

impl VertexPool {
    fn new(tol: f64) -> Self {
        Self { pts: Vec::new(), grid: HashMap::new(), tol }
    }
    fn key(&self, p: P) -> (i64, i64) {
        let s = 1.0 / (self.tol.max(1e-12) * 4.0);
        ((p.0 * s).round() as i64, (p.1 * s).round() as i64)
    }
    fn insert(&mut self, p: P) -> usize {
        let k = self.key(p);
        let t2 = self.tol * self.tol;
        for dk in -1..=1i64 {
            for dl in -1..=1i64 {
                if let Some(ids) = self.grid.get(&(k.0 + dk, k.1 + dl)) {
                    for &i in ids {
                        if dist2(self.pts[i], p) <= t2 {
                            return i;
                        }
                    }
                }
            }
        }
        let id = self.pts.len();
        self.pts.push(p);
        self.grid.entry(k).or_default().push(id);
        id
    }
}

/// Planar arrangement of polyline chains; computes the cells of the induced
/// subdivision.
pub struct Arrangement {
    verts: Vec<P>,
    segs: Vec<Seg>,
}

impl Arrangement {
    /// Build the arrangement: snap vertices, split chains at all mutual
    /// intersections, drop zero-length pieces.
    #[must_use]
    pub fn build(chains: &[ChainInput], opt: &ArrangementOptions) -> Self {
        let tol = opt.tol;
        // 1. Flatten chains into raw segments (endpoint coords).
        struct Raw {
            p: P,
            q: P,
            tag: u32,
        }
        let mut raw: Vec<Raw> = Vec::new();
        for ch in chains {
            for w in ch.pts.windows(2) {
                if dist2(w[0], w[1]) > tol * tol {
                    raw.push(Raw { p: w[0], q: w[1], tag: ch.tag });
                }
            }
        }
        // 2. Collect split parameters per raw segment.
        let mut cuts: Vec<Vec<f64>> = vec![Vec::new(); raw.len()];
        for i in 0..raw.len() {
            for j in (i + 1)..raw.len() {
                if let Some((t, s)) =
                    seg_seg_params(raw[i].p, raw[i].q, raw[j].p, raw[j].q, tol)
                {
                    if t > 0.0 && t < 1.0 {
                        cuts[i].push(t);
                    }
                    if s > 0.0 && s < 1.0 {
                        cuts[j].push(s);
                    }
                }
            }
        }
        // 3. Emit atomic segments through the vertex pool.
        let mut pool = VertexPool::new(tol);
        let mut segs: Vec<Seg> = Vec::new();
        for (i, r) in raw.iter().enumerate() {
            let mut ts = std::mem::take(&mut cuts[i]);
            ts.push(0.0);
            ts.push(1.0);
            ts.sort_by(f64::total_cmp);
            ts.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
            let lerp = |t: f64| -> P {
                (r.p.0 + (r.q.0 - r.p.0) * t, r.p.1 + (r.q.1 - r.p.1) * t)
            };
            for w in ts.windows(2) {
                let a = pool.insert(lerp(w[0]));
                let b = pool.insert(lerp(w[1]));
                if a != b {
                    segs.push(Seg { a, b, tag: r.tag, forward: true });
                }
            }
        }
        // 4. Drop exact duplicates (same vertex pair) keeping the first tag.
        let mut seen: HashMap<(usize, usize), ()> = HashMap::new();
        segs.retain(|s| {
            let k = (s.a.min(s.b), s.a.max(s.b));
            seen.insert(k, ()).is_none()
        });
        Self { verts: pool.pts, segs }
    }

    /// Extract the bounded cells of the subdivision.
    ///
    /// Each returned [`Cell`] has a CCW outer loop and CW holes.  The single
    /// unbounded outer cycle is discarded.
    #[must_use]
    pub fn cells(&self) -> Vec<Cell> {
        let nv = self.verts.len();
        // Half-edges: 2 per segment. he = 2*seg (a→b), 2*seg+1 (b→a).
        let nh = self.segs.len() * 2;
        let he_from = |h: usize| -> usize {
            let s = &self.segs[h / 2];
            if h % 2 == 0 { s.a } else { s.b }
        };
        let he_to = |h: usize| -> usize {
            let s = &self.segs[h / 2];
            if h % 2 == 0 { s.b } else { s.a }
        };
        // Outgoing half-edges per vertex, sorted CCW by angle.
        let mut out: Vec<Vec<usize>> = vec![Vec::new(); nv];
        for h in 0..nh {
            out[he_from(h)].push(h);
        }
        for (v, list) in out.iter_mut().enumerate() {
            let vp = self.verts[v];
            list.sort_by(|&h1, &h2| {
                let a1 = {
                    let q = self.verts[he_to(h1)];
                    (q.1 - vp.1).atan2(q.0 - vp.0)
                };
                let a2 = {
                    let q = self.verts[he_to(h2)];
                    (q.1 - vp.1).atan2(q.0 - vp.0)
                };
                a1.total_cmp(&a2)
            });
        }
        // index of each half-edge within its origin's sorted ring
        let mut ring_pos: Vec<usize> = vec![0; nh];
        for list in &out {
            for (i, &h) in list.iter().enumerate() {
                ring_pos[h] = i;
            }
        }
        // next(h): at w = to(h), take twin(h) (outgoing at w) and step to the
        // ring predecessor (CW neighbour) → traverses faces with interior on
        // the LEFT (CCW outer loops).
        let twin = |h: usize| -> usize { h ^ 1 };
        let next = |h: usize| -> usize {
            let w = he_to(h);
            let r = twin(h);
            let ring = &out[w];
            let i = ring_pos[r];
            ring[(i + ring.len() - 1) % ring.len()]
        };
        // Walk cycles.
        let mut visited = vec![false; nh];
        struct Cycle {
            hes: Vec<usize>,
            area: f64,
        }
        let mut cycles: Vec<Cycle> = Vec::new();
        for h0 in 0..nh {
            if visited[h0] {
                continue;
            }
            let mut hes = Vec::new();
            let mut area = 0.0;
            let mut h = h0;
            loop {
                visited[h] = true;
                hes.push(h);
                let p = self.verts[he_from(h)];
                let q = self.verts[he_to(h)];
                area += p.0 * q.1 - q.0 * p.1;
                h = next(h);
                if h == h0 {
                    break;
                }
            }
            cycles.push(Cycle { hes, area: area * 0.5 });
        }
        // Positive cycles = cell outer loops; negative = holes or the outer
        // unbounded face.  Assign each negative cycle to the smallest positive
        // cycle containing it.
        let mut cells: Vec<(Cycle, Vec<Cycle>)> = Vec::new();
        let mut negatives: Vec<Cycle> = Vec::new();
        for c in cycles {
            if c.area > 0.0 {
                cells.push((c, Vec::new()));
            } else {
                negatives.push(c);
            }
        }
        for neg in negatives {
            // Probe: a point just LEFT of the cycle's first half-edge — i.e.
            // inside the face this CW cycle bounds.  (The midpoint itself lies
            // ON shared boundary; the enclosed disk on the right is a
            // DIFFERENT face.)  The unbounded cycle's left side is outside
            // every positive cell, so it gets dropped naturally.
            let h = neg.hes[0];
            let p = self.verts[he_from(h)];
            let q = self.verts[he_to(h)];
            let d = sub(q, p);
            let len = (d.0 * d.0 + d.1 * d.1).sqrt().max(1e-300);
            let eps = (len * 1e-7).max(1e-12);
            let mid = (
                (p.0 + q.0) * 0.5 - d.1 / len * eps,
                (p.1 + q.1) * 0.5 + d.0 / len * eps,
            );
            let mut best: Option<(usize, f64)> = None;
            for (i, (outer, _)) in cells.iter().enumerate() {
                let poly: Vec<P> = outer.hes.iter().map(|&h| self.verts[he_from(h)]).collect();
                if point_in_polygon(mid, &poly) && outer.area > 0.0 {
                    if best.is_none() || outer.area < best.unwrap().1 {
                        best = Some((i, outer.area));
                    }
                }
            }
            if let Some((i, _)) = best {
                cells[i].1.push(neg);
            }
            // else: the unbounded face — dropped.
        }
        // Materialize.
        cells
            .into_iter()
            .map(|(outer, holes)| {
                let area = outer.area + holes.iter().map(|h| h.area).sum::<f64>();
                Cell {
                    outer: self.loop_steps(&outer.hes),
                    holes: holes.iter().map(|h| self.loop_steps(&h.hes)).collect(),
                    area,
                }
            })
            .collect()
    }

    fn loop_steps(&self, hes: &[usize]) -> Vec<LoopStep> {
        hes.iter()
            .map(|&h| {
                let s = &self.segs[h / 2];
                let fwd = h % 2 == 0;
                let (a, b) = if fwd { (s.a, s.b) } else { (s.b, s.a) };
                LoopStep {
                    tag: s.tag,
                    pts: vec![self.verts[a], self.verts[b]],
                    forward: fwd == s.forward,
                }
            })
            .collect()
    }
}

fn loop_polygon(steps: &[LoopStep]) -> Vec<P> {
    steps.iter().map(|s| s.pts[0]).collect()
}

fn point_in_polygon(p: P, poly: &[P]) -> bool {
    let mut inside = false;
    let n = poly.len();
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = poly[i];
        let (xj, yj) = poly[j];
        if ((yi > p.1) != (yj > p.1))
            && (p.0 < (xj - xi) * (p.1 - yi) / (yj - yi) + xi)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// Minimum distance from `p` to any edge of the closed polygon `poly`.
fn dist_to_polygon(p: P, poly: &[P]) -> f64 {
    let n = poly.len();
    if n < 2 {
        return f64::MAX;
    }
    let mut best = f64::MAX;
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let len2 = dx * dx + dy * dy;
        let t = if len2 > 0.0 {
            (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / len2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let (cx, cy) = (a.0 + dx * t, a.1 + dy * t);
        best = best.min(((p.0 - cx).powi(2) + (p.1 - cy).powi(2)).sqrt());
    }
    best
}

fn interior_point_of(outer: &[P], holes: &[Vec<P>]) -> Option<(f64, f64)> {
    if outer.is_empty() {
        return None;
    }
    // Centroid first; fall back to sampling along a horizontal through the
    // polygon's vertical middle.
    let n = outer.len() as f64;
    let c = outer.iter().fold((0.0, 0.0), |acc, p| (acc.0 + p.0, acc.1 + p.1));
    let c = (c.0 / n, c.1 / n);
    // a scale for the "strictly interior" margin (polygon diagonal)
    let (xmn, xmx) = outer.iter().fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.0), a.1.max(p.0)));
    let (ymn, ymx) = outer.iter().fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.1), a.1.max(p.1)));
    let margin = 1e-6 * ((xmx - xmn).hypot(ymx - ymn)).max(1e-9);
    let ok = |p: P| -> bool {
        point_in_polygon(p, outer)
            && holes.iter().all(|h| !point_in_polygon(p, h))
            // reject points ON the boundary: the vertex-centroid of a
            // non-convex (e.g. L-shaped) cell can land exactly on a concave
            // vertex, which classifies unreliably downstream.
            && dist_to_polygon(p, outer) > margin
            && holes.iter().all(|h| dist_to_polygon(p, h) > margin)
    };
    if ok(c) {
        return Some(c);
    }
    let (ymin, ymax) = outer
        .iter()
        .fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.1), a.1.max(p.1)));
    let (xmin, xmax) = outer
        .iter()
        .fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.0), a.1.max(p.0)));
    for fy in [0.5, 0.37, 0.63, 0.21, 0.79] {
        let y = ymin + (ymax - ymin) * fy;
        for fx in 1..40 {
            let x = xmin + (xmax - xmin) * (fx as f64) / 40.0;
            if ok((x, y)) {
                return Some((x, y));
            }
        }
    }
    None
}

/// Build an arrangement on a **u-periodic** domain (cylinder θ, torus φ) by
/// the canonical seam cut: all chain `u` values are wrapped into `[0, period]`
/// (chains crossing the seam are split at it), and two seam chains at `u=0`
/// and `u=period` spanning `v_range` are added with `seam_tag`.
///
/// The result is the standard STEP representation: the periodic face's wire
/// traverses the seam explicitly; cells left and right of the seam share the
/// same 3D seam edge (identical 3D geometry from `u=0` and `u=period`).
///
/// **Precondition:** adjacent points of every chain differ in `u` by less
/// than `period / 2` (seam-crossing detection unwraps by nearest branch).
/// Projected boundary circles / SSI polylines are densely sampled, so this
/// holds naturally; do not feed 2-point full-period chords.
#[must_use]
pub fn build_periodic_u(
    chains: &[ChainInput],
    period: f64,
    v_range: (f64, f64),
    seam_tag: u32,
    opt: &ArrangementOptions,
) -> Arrangement {
    let mut out: Vec<ChainInput> = Vec::new();
    for ch in chains {
        // 1. Unwrap u continuously along the chain.
        let mut un: Vec<(f64, f64)> = Vec::with_capacity(ch.pts.len());
        for (i, &(u, v)) in ch.pts.iter().enumerate() {
            if i == 0 {
                un.push((u, v));
            } else {
                let prev = un[i - 1].0;
                let mut du = (u - prev) % period;
                if du > period / 2.0 {
                    du -= period;
                } else if du < -period / 2.0 {
                    du += period;
                }
                un.push((prev + du, v));
            }
        }
        // 2. Split wherever a segment crosses a multiple of the period, then
        //    map every piece into [0, period] by its own branch offset.
        let mut piece: Vec<(f64, f64)> = vec![un[0]];
        let mut flush = |piece: &mut Vec<(f64, f64)>| {
            if piece.len() >= 2 {
                // Branch from the minimum u over ALL points: endpoints may sit
                // exactly on the seam while the interior dips into a branch.
                let umin = piece.iter().map(|p| p.0).fold(f64::MAX, f64::min);
                let base = (umin + 1e-9).div_euclid(period);
                let off = -base * period;
                out.push(ChainInput {
                    pts: piece.iter().map(|&(u, v)| (u + off, v)).collect(),
                    tag: ch.tag,
                });
            }
            let last = *piece.last().unwrap();
            piece.clear();
            piece.push(last);
        };
        for w in 1..un.len() {
            let (u0, v0) = un[w - 1];
            let (u1, v1) = un[w];
            // crossing points at integer multiples of the period between u0,u1
            let (lo, hi) = if u0 <= u1 { (u0, u1) } else { (u1, u0) };
            let mut ks: Vec<f64> = Vec::new();
            let mut k = (lo / period).ceil() * period;
            while k < hi {
                if k > lo {
                    ks.push(k);
                }
                k += period;
            }
            if u0 > u1 {
                ks.reverse();
            }
            for kc in ks {
                let t = (kc - u0) / (u1 - u0);
                let vc = v0 + (v1 - v0) * t;
                piece.push((kc, vc));
                flush(&mut piece);
            }
            piece.push((u1, v1));
            // A vertex landing EXACTLY on a branch boundary is also a seam
            // touch — split there so each piece stays within one branch.
            let on_multiple = {
                let r = u1 / period;
                (r - r.round()).abs() * period < 1e-9
            };
            if on_multiple && w + 1 < un.len() {
                flush(&mut piece);
            }
        }
        if piece.len() >= 2 {
            flush(&mut piece);
        }
    }
    // 3. Seam chains at u=0 and u=period.
    out.push(ChainInput {
        pts: vec![(0.0, v_range.0), (0.0, v_range.1)],
        tag: seam_tag,
    });
    out.push(ChainInput {
        pts: vec![(period, v_range.0), (period, v_range.1)],
        tag: seam_tag,
    });
    Arrangement::build(&out, opt)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rect_chain(x0: f64, y0: f64, x1: f64, y1: f64, tag: u32) -> ChainInput {
        ChainInput {
            pts: vec![(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)],
            tag,
        }
    }

    /// A bare rectangle = exactly one CCW cell, no holes.
    #[test]
    fn rectangle_single_cell() {
        let arr = Arrangement::build(
            &[rect_chain(0.0, 0.0, 4.0, 2.0, 1)],
            &ArrangementOptions::default(),
        );
        let cells = arr.cells();
        assert_eq!(cells.len(), 1);
        assert!((cells[0].area - 8.0).abs() < 1e-12, "area {}", cells[0].area);
        assert!(cells[0].holes.is_empty());
        // Interior probe is inside.
        let ip = cells[0].interior_point().unwrap();
        assert!(ip.0 > 0.0 && ip.0 < 4.0 && ip.1 > 0.0 && ip.1 < 2.0);
    }

    /// A closed imprint loop inside the rectangle → 2 cells: the loop interior
    /// and the rectangle-with-hole.
    #[test]
    fn closed_imprint_splits_two_cells() {
        let arr = Arrangement::build(
            &[
                rect_chain(0.0, 0.0, 4.0, 2.0, 1),
                rect_chain(1.0, 0.5, 2.0, 1.5, 2),
            ],
            &ArrangementOptions::default(),
        );
        let mut cells = arr.cells();
        cells.sort_by(|a, b| a.area.total_cmp(&b.area));
        assert_eq!(cells.len(), 2, "cells: {:?}", cells.iter().map(|c| c.area).collect::<Vec<_>>());
        // inner cell: area 1, no holes
        assert!((cells[0].area - 1.0).abs() < 1e-12);
        assert!(cells[0].holes.is_empty());
        // outer cell: 8 − 1 = 7, one hole
        assert!((cells[1].area - 7.0).abs() < 1e-12);
        assert_eq!(cells[1].holes.len(), 1);
        // tags survive
        assert!(cells[0].outer.iter().all(|s| s.tag == 2));
        assert!(cells[1].outer.iter().all(|s| s.tag == 1));
        assert!(cells[1].holes[0].iter().all(|s| s.tag == 2));
        // interior probes don't collide
        let p0 = cells[0].interior_point().unwrap();
        let p1 = cells[1].interior_point().unwrap();
        assert!(p0.0 > 1.0 && p0.0 < 2.0 && p0.1 > 0.5 && p0.1 < 1.5);
        assert!(!(p1.0 > 1.0 && p1.0 < 2.0 && p1.1 > 0.5 && p1.1 < 1.5));
    }

    /// An open imprint chain crossing the rectangle splits it into 2 cells —
    /// the parallel-tube case (two open SSI lines across the lateral face).
    #[test]
    fn open_chain_splits_rectangle() {
        let arr = Arrangement::build(
            &[
                rect_chain(0.0, 0.0, 4.0, 2.0, 1),
                ChainInput { pts: vec![(-0.5, 1.0), (4.5, 1.0)], tag: 7 },
            ],
            &ArrangementOptions::default(),
        );
        let mut cells = arr.cells();
        cells.sort_by(|a, b| a.area.total_cmp(&b.area));
        assert_eq!(cells.len(), 2, "areas: {:?}", cells.iter().map(|c| c.area).collect::<Vec<_>>());
        assert!((cells[0].area - 4.0).abs() < 1e-9);
        assert!((cells[1].area - 4.0).abs() < 1e-9);
        // Both cells' boundaries contain both tags.
        for c in &cells {
            assert!(c.outer.iter().any(|s| s.tag == 1));
            assert!(c.outer.iter().any(|s| s.tag == 7));
        }
    }

    /// T-junction: an open chain ending exactly ON the rectangle edge still
    /// splits cleanly (endpoint-on-interior splitting).
    #[test]
    fn t_junction_splits() {
        let arr = Arrangement::build(
            &[
                rect_chain(0.0, 0.0, 4.0, 2.0, 1),
                ChainInput { pts: vec![(2.0, 0.0), (2.0, 2.0)], tag: 9 },
            ],
            &ArrangementOptions::default(),
        );
        let cells = arr.cells();
        assert_eq!(cells.len(), 2);
        let total: f64 = cells.iter().map(|c| c.area).sum();
        assert!((total - 8.0).abs() < 1e-9);
    }

    /// Crossing imprints: two open chains crossing each other AND the boundary
    /// → 4 cells, conservation of area.
    #[test]
    fn two_crossing_chains_four_cells() {
        let arr = Arrangement::build(
            &[
                rect_chain(0.0, 0.0, 4.0, 2.0, 1),
                ChainInput { pts: vec![(2.0, -0.1), (2.0, 2.1)], tag: 2 },
                ChainInput { pts: vec![(-0.1, 1.0), (4.1, 1.0)], tag: 3 },
            ],
            &ArrangementOptions::default(),
        );
        let cells = arr.cells();
        assert_eq!(cells.len(), 4);
        let total: f64 = cells.iter().map(|c| c.area).sum();
        assert!((total - 8.0).abs() < 1e-9);
        for c in &cells {
            assert!((c.area - 2.0).abs() < 1e-9);
        }
    }

    /// Periodic lateral face of a cylinder, untouched: bottom/top boundary
    /// circles + the two seam chains → exactly one cell of area period×len.
    #[test]
    fn periodic_plain_lateral_one_cell() {
        let p = std::f64::consts::TAU;
        let arr = build_periodic_u(
            &[
                ChainInput { pts: vec![(0.0, 0.0), (p / 3.0, 0.0), (2.0 * p / 3.0, 0.0), (p, 0.0)], tag: 1 },
                ChainInput { pts: vec![(p, 3.0), (2.0 * p / 3.0, 3.0), (p / 3.0, 3.0), (0.0, 3.0)], tag: 2 },
            ],
            p,
            (0.0, 3.0),
            99,
            &ArrangementOptions::default(),
        );
        let cells = arr.cells();
        assert_eq!(cells.len(), 1);
        assert!((cells[0].area - p * 3.0).abs() < 1e-9);
        // The cell's wire contains the seam tag twice (up + down).
        let seam_steps = cells[0].outer.iter().filter(|s| s.tag == 99).count();
        assert!(seam_steps >= 2, "seam traversed twice, got {seam_steps}");
    }

    /// A closed imprint loop CROSSING the seam splits into branch pieces and
    /// the cells still conserve area: 2 half-bites + the big cell.
    #[test]
    fn periodic_hole_across_seam() {
        let p = std::f64::consts::TAU;
        let arr = build_periodic_u(
            &[
                ChainInput { pts: vec![(0.0, 0.0), (p / 3.0, 0.0), (2.0 * p / 3.0, 0.0), (p, 0.0)], tag: 1 },
                ChainInput { pts: vec![(p, 2.0), (2.0 * p / 3.0, 2.0), (p / 3.0, 2.0), (0.0, 2.0)], tag: 2 },
                // diamond centred on the seam (u≈0), crossing it twice
                ChainInput {
                    pts: vec![
                        (-0.3, 1.0),
                        (0.0, 0.7),
                        (0.3, 1.0),
                        (0.0, 1.3),
                        (-0.3, 1.0),
                    ],
                    tag: 5,
                },
            ],
            p,
            (0.0, 2.0),
            99,
            &ArrangementOptions::default(),
        );
        let cells = arr.cells();
        let total: f64 = cells.iter().map(|c| c.area).sum();
        assert!(
            (total - p * 2.0).abs() < 1e-9,
            "area conservation: {total} vs {}",
            p * 2.0
        );
        assert_eq!(cells.len(), 3, "two half-bites + main cell");
        // The two small cells are the diamond halves: combined area = diamond.
        let mut areas: Vec<f64> = cells.iter().map(|c| c.area).collect();
        areas.sort_by(f64::total_cmp);
        assert!((areas[0] + areas[1] - 0.18).abs() < 1e-9, "diamond halves: {areas:?}");
    }

    /// Polyline (multi-point) imprint loop keeps per-step geometry and the
    /// invariant Σ cell areas == domain area.
    #[test]
    fn polyline_loop_area_conservation() {
        // a diamond imprint
        let arr = Arrangement::build(
            &[
                rect_chain(0.0, 0.0, 4.0, 2.0, 1),
                ChainInput {
                    pts: vec![(2.0, 0.4), (2.8, 1.0), (2.0, 1.6), (1.2, 1.0), (2.0, 0.4)],
                    tag: 5,
                },
            ],
            &ArrangementOptions::default(),
        );
        let cells = arr.cells();
        assert_eq!(cells.len(), 2);
        let total: f64 = cells.iter().map(|c| c.area).sum();
        assert!((total - 8.0).abs() < 1e-9, "total {total}");
    }
}
