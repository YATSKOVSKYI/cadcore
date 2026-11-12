//! **Universal filament fuse** — the Boolean union of an arbitrary set of
//! swept filaments (straight legs, torus-fillet elbows, flat end caps) into
//! ONE watertight solid, by *boundary tracing* plus per-face arrangement.
//!
//! Imprinting each crossing as a closed HOLE on the two faces it was traced on
//! is exact only while the crossing sits strictly inside both faces, and real
//! lattices break that all along their outline:
//!
//! * the first row of one layer runs over the **free ends** of the layer
//!   below, so the crossing window straddles a cap rim;
//! * the U-turn connector of one layer lies **parallel** on the outer leg of
//!   the next, so the tubes meet in two ruling lines joined over the elbows
//!   (no closed loop on any single member);
//! * honeycomb bonds **end inside** the row they bridge (a T-junction: the
//!   cap and rim of the bond are swallowed whole).
//!
//! So the fuse computes the real union boundary:
//!
//! ```text
//! 1 trace     for every pair of filaments, march along ∂A ∩ ∂B.  The marcher
//!             knows each filament as a chain of MEMBERS (leg, elbow, cap)
//!             glued along junction circles; when the curve leaves a member it
//!             computes the exact crossing of the junction circle with the
//!             other surface (a shared JOINT) and continues on the neighbour.
//!             The result is a closed loop made of PIECES, each lying on one
//!             face of A and one face of B.
//! 2 chains    every junction circle is sampled ONCE (with all its joints and
//!             the adjacent faces' parametric seam points as vertices) and
//!             handed to both faces that share it; every piece is handed to
//!             both faces it lies on.
//! 3 arrange   each face: DCEL of its boundary circles + pieces; a cell is
//!             kept iff its interior lies outside every OTHER filament's
//!             volume (analytic containment, no ray casting).
//! 4 stitch    kept cells → faces through one [`Assembler`]; shared edges weld
//!             because both sides present byte-identical 3-D points.
//! ```
//!
//! Because the geometry is a true union, the half-space cut
//! ([`crate::boolean::half_space_cut`]) and the plate weld downstream see a
//! clean solid: no interpenetrating tubes, no self-crossing wires.

use std::collections::HashMap;

use cadcore_geom::{Circle3, Plane3};
use cadcore_math::{Point3, Vec3};
use cadcore_topo::{BRep, FaceGeom, FaceNormal, Shell, ShellId, Solid};

use super::assembly::Assembler;
use super::cells::FaceChain;
use super::domain::FaceDomain;
use super::filament::{ElbowSeg, LegSeg};
use super::filament::Filament;
use crate::boolean::Aabb;
use crate::geom::composite::{seeds_between, CompositeTube};
use crate::geom::intersect::intersection_point_near;
use crate::geom::refine::AnalyticSurface;

/// Tunables of the universal fuse.
#[derive(Debug, Clone, Copy)]
pub struct FuseOptions {
    /// Nominal marching step along an intersection curve (mm).  Default 0.02.
    pub step: f64,
    /// Target spacing (mm) of the points kept on every traced piece; the
    /// endpoints are always kept.  The STEP writer interpolates each piece with
    /// a cubic B-spline, whose error grows with the fourth power of the
    /// spacing — 0.03 mm keeps a 0.275 mm-radius filament's trim curves within
    /// a tenth of a micron of both surfaces.  Default 0.03.
    pub sample_spacing: f64,
    /// Uniform samples per junction circle (joints and seam points are added
    /// on top).  Default 48.
    pub circle_samples: usize,
    /// Assembler weld tolerance (mm).  Shared points are byte-identical on
    /// both sides, so this only absorbs float noise — and collapses any
    /// sliver segment shorter than it (OCC rejects edges under 1e-7).
    /// Default 1e-5.
    pub weld: f64,
    /// Hard cap on marching steps per loop.  Default 40 000.
    pub max_steps: usize,
}

impl Default for FuseOptions {
    fn default() -> Self {
        Self {
            step: 0.02,
            sample_spacing: 0.03,
            circle_samples: 48,
            weld: 1e-5,
            max_steps: 40_000,
        }
    }
}

/// What the fuse did — for diagnostics and tests.
#[derive(Debug, Clone, Default)]
pub struct FuseReport {
    /// Closed boundary loops traced.
    pub loops: usize,
    /// Curve pieces (one per face pair along a loop).
    pub pieces: usize,
    /// Joints placed on junction circles.
    pub joints: usize,
    /// Traces started but abandoned (tangency, divergence, budget).
    pub aborted: usize,
    /// Faces emitted.
    pub faces: usize,
    /// Human-readable notes about anything irregular.
    pub notes: Vec<String>,
}

/// Fuse `filaments` into one closed shell of `brep`.  Default options.
pub fn fuse_filaments(brep: &mut BRep, filaments: &[Filament]) -> ShellId {
    fuse_filaments_with(brep, filaments, &FuseOptions::default()).0
}

/// [`fuse_filaments`] with explicit options; also returns the report.
pub fn fuse_filaments_with(
    brep: &mut BRep,
    filaments: &[Filament],
    opts: &FuseOptions,
) -> (ShellId, FuseReport) {
    let mut report = FuseReport::default();
    let bodies: Vec<Body> = filaments.iter().map(Body::new).collect();
    let nb = bodies.len();

    // ── 1. Trace every boundary loop between every pair of filaments ────────
    let mut loops: Vec<TracedLoop> = Vec::new();
    let mut known: PointGrid = PointGrid::new(opts.step * 1.5);
    for i in 0..nb {
        for j in (i + 1)..nb {
            trace_pair(&bodies, i, j, opts, &mut loops, &mut known, &mut report);
        }
        // Self-crossings (a filament weaving back over itself) are not traced
        // here; note them so the caller knows the output interpenetrates there.
        let b = &bodies[i];
        let nl = b.legs.len();
        for la in 0..nl {
            for lb in (la + 2)..nl {
                if !b.member_box(Member::Leg(la)).overlaps(&b.member_box(Member::Leg(lb))) {
                    continue;
                }
                // the tubes interpenetrate only if their AXES come closer than
                // the two radii (boxes of parallel diagonal legs always overlap)
                let (sa, sb) = (&b.legs[la], &b.legs[lb]);
                let seg = |l: &LegSeg| (l.surf.frame.origin, l.surf.frame.origin + l.surf.axis().as_vec() * l.length);
                let (a0, a1) = seg(sa);
                let (b0, b1) = seg(sb);
                if segment_distance(a0, a1, b0, b1) < sa.surf.radius + sb.surf.radius - 1e-6 {
                    report
                        .notes
                        .push(format!("filament {i}: legs {la} and {lb} overlap (self-crossing not fused)"));
                }
            }
        }
    }

    // ── 2. Shared chains: junction circles + pieces ──────────────────────────
    // Park every leg's parametric seam (u = 0) in the widest angular gap left
    // by the pieces on it.  A piece crossing the seam has to be split there on
    // both faces; a piece LYING ALONG the seam (a neighbour's cap plane through
    // this leg's axis puts the boundary on the top ruling line — the everyday
    // woodpile corner) degenerates the periodic DCEL outright.  Rotating the
    // frame about the axis changes nothing geometric.
    let mut bodies = bodies;
    park_seams(&mut bodies, &loops);
    let bodies = bodies;
    let mut junctions: HashMap<JunctionKey, JunctionChain> = HashMap::new();
    for (bi, b) in bodies.iter().enumerate() {
        for l in 0..b.legs.len() {
            for end in 0..2u8 {
                let key = JunctionKey { body: bi, leg: l, end };
                let circle = b.leg_circle(l, end);
                let mut jc = JunctionChain { circle, mandatory: Vec::new() };
                jc.mandatory.extend(b.seam_points_at(l, end));
                junctions.insert(key, jc);
            }
        }
    }
    for lp in &mut loops {
        decimate_loop(lp, opts.sample_spacing);
        for piece in &mut lp.pieces {
            let sa = bodies[piece.a.0].surface(piece.a.1);
            let sb = bodies[piece.b.0].surface(piece.b.1);
            for side in [1.0, -1.0] {
                split_piece_at_seam(&mut piece.pts, &sa, &sb, side, &mut report);
                split_piece_at_seam(&mut piece.pts, &sb, &sa, side, &mut report);
            }
        }
        for j in &lp.joints {
            if let Some(jc) = junctions.get_mut(&j.key) {
                jc.mandatory.push(j.p);
            }
        }
        report.pieces += lp.pieces.len();
        report.joints += lp.joints.len();
    }
    report.loops = loops.len();
    if std::env::var("CADCORE_DUMP_FUSE").map_or(false, |v| v == "pieces") {
        for (li, lp) in loops.iter().enumerate() {
            eprintln!("[fuse][loop {li}] {} pieces, {} joints", lp.pieces.len(), lp.joints.len());
            for pc in &lp.pieces {
                let a = pc.pts[0];
                let b = *pc.pts.last().unwrap();
                eprintln!(
                    "[fuse][loop {li}]   {:?}/{:?} n={} ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                    pc.a, pc.b, pc.pts.len(), a.x, a.y, a.z, b.x, b.y, b.z
                );
            }
            for j in &lp.joints {
                eprintln!("[fuse][loop {li}]   joint {:?} ({:.4},{:.4},{:.4})", j.key, j.p.x, j.p.y, j.p.z);
            }
        }
    }

    let sampled: HashMap<JunctionKey, Vec<Point3>> = junctions
        .iter()
        .map(|(k, jc)| (*k, jc.sample(opts.circle_samples)))
        .collect();

    // ── 3 + 4. Arrange every face and stitch the kept cells ─────────────────
    let mut asm = Assembler::new(brep, opts.weld);
    for (bi, b) in bodies.iter().enumerate() {
        let drop = |p: Point3| {
            bodies
                .iter()
                .enumerate()
                .any(|(k, other)| k != bi && other.tube.contains(p))
        };
        for m in b.members() {
            let domain = b.domain(m);
            let mut chains: Vec<FaceChain> = Vec::new();
            let mut tag = 1u32;
            for key in b.junctions_of(m) {
                let Some(pts) = sampled.get(&JunctionKey { body: bi, leg: key.0, end: key.1 }) else {
                    continue;
                };
                let uv: Vec<(f64, f64)> = pts.iter().map(|&p| domain.uv(p)).collect();
                chains.push(FaceChain { pts: uv, tag });
                tag += 1;
            }
            let mut ptag = 1000u32;
            for lp in &loops {
                for piece in &lp.pieces {
                    if piece.a != (bi, m) && piece.b != (bi, m) {
                        continue;
                    }
                    if piece.pts.len() < 2 {
                        continue;
                    }
                    let uv: Vec<(f64, f64)> = piece.pts.iter().map(|&p| domain.uv(p)).collect();
                    chains.push(FaceChain { pts: uv, tag: ptag });
                    ptag += 1;
                }
            }
            let cells = crate::arrange::cells::arrange_and_classify_halved(&domain, &chains, &drop);
            let geom = match m {
                Member::Leg(l) => FaceGeom::Cylinder(b.legs[l].surf),
                Member::Elbow(e) => FaceGeom::Torus(b.elbows[e].surf),
                Member::Cap(c) => FaceGeom::Plane(b.caps[c].plane),
            };
            if std::env::var("CADCORE_DUMP_FUSE").is_ok() {
                eprintln!(
                    "[fuse] body {bi} {m:?}: chains={} (pieces={}) cells={} kept={}",
                    chains.len(),
                    ptag - 1000,
                    cells.len(),
                    cells.iter().filter(|c| c.keep).count()
                );
            }
            if std::env::var("CADCORE_DUMP_FUSE").map_or(false, |v| v == "cells") {
                for c in &cells {
                    let ip = c.cell.interior_point().map(|(u, v)| domain.lift(u, v));
                    let tags: std::collections::BTreeSet<u32> = c.cell.outer.iter().map(|st| st.tag).collect();
                    eprintln!(
                        "[fuse][cell] body {bi} {m:?}: area={:.4} keep={} holes={} probe={:?} tags={:?}",
                        c.cell.area, c.keep, c.cell.holes.len(),
                        ip.map(|p| ((p.x * 1e4).round() / 1e4, (p.y * 1e4).round() / 1e4, (p.z * 1e4).round() / 1e4)),
                        tags
                    );
                }
            }
            for cell in &cells {
                if cell.keep {
                    asm.emit_cell(&domain, geom.clone(), FaceNormal::Same, cell);
                }
            }
        }
    }

    let faces = asm.faces().to_vec();
    report.faces = faces.len();
    let shell = brep.add_shell(Shell {
        faces: faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = brep.add_solid(Solid {
        shells: vec![shell],
        name: Some("scaffold".into()),
    });
    brep.shells[shell].solid = solid;
    for f in &faces {
        brep.faces[*f].shell = shell;
    }
    (shell, report)
}

// ── Bodies and members ───────────────────────────────────────────────────────

/// One face of a filament solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum Member {
    Leg(usize),
    Elbow(usize),
    /// `Cap(0)` closes the first leg's start, `Cap(1)` the last leg's end.
    Cap(usize),
}

/// A flat end cap.
#[derive(Debug, Clone, Copy)]
struct Cap {
    /// Carrier plane, normal pointing OUT of the filament.
    plane: Plane3,
    centre: Point3,
    radius: f64,
    /// The leg it closes and which end (`0` = start, `1` = end).
    leg: usize,
    end: u8,
}

/// A junction circle is attached to a leg end; elbows and caps refer to it
/// through the leg.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct JunctionKey {
    body: usize,
    leg: usize,
    end: u8,
}

struct Body {
    legs: Vec<LegSeg>,
    elbows: Vec<ElbowSeg>,
    caps: Vec<Cap>,
    tube: CompositeTube,
}

impl Body {
    fn new(f: &Filament) -> Self {
        let n = f.legs.len();
        assert!(n >= 1, "a filament needs at least one leg");
        let first = &f.legs[0];
        let last = &f.legs[n - 1];
        let caps = vec![
            Cap {
                plane: Plane3::from_origin_normal(first.surf.frame.origin, first.surf.axis().flip()),
                centre: first.surf.frame.origin,
                radius: first.surf.radius,
                leg: 0,
                end: 0,
            },
            Cap {
                plane: Plane3::from_origin_normal(
                    last.surf.frame.origin + last.surf.axis().as_vec() * last.length,
                    last.surf.axis(),
                ),
                centre: last.surf.frame.origin + last.surf.axis().as_vec() * last.length,
                radius: last.surf.radius,
                leg: n - 1,
                end: 1,
            },
        ];
        let mut tube = CompositeTube::new();
        for l in &f.legs {
            tube = tube.with_leg(l.surf, l.length);
        }
        for e in &f.elbows {
            tube = tube.with_elbow(e.surf, e.theta_lo, e.theta_hi);
        }
        Self {
            legs: f.legs.clone(),
            elbows: f.elbows.clone(),
            caps,
            tube,
        }
    }

    fn members(&self) -> Vec<Member> {
        let mut out: Vec<Member> = (0..self.legs.len()).map(Member::Leg).collect();
        out.extend((0..self.elbows.len()).map(Member::Elbow));
        out.extend((0..self.caps.len()).map(Member::Cap));
        out
    }

    fn surface(&self, m: Member) -> AnalyticSurface {
        match m {
            Member::Leg(l) => AnalyticSurface::Cylinder(self.legs[l].surf),
            Member::Elbow(e) => AnalyticSurface::Torus(self.elbows[e].surf),
            Member::Cap(c) => AnalyticSurface::Plane(self.caps[c].plane),
        }
    }

    fn domain(&self, m: Member) -> FaceDomain {
        match m {
            Member::Leg(l) => FaceDomain::CylinderBand {
                surf: self.legs[l].surf,
                length: self.legs[l].length,
            },
            Member::Elbow(e) => FaceDomain::TorusPatch {
                surf: self.elbows[e].surf,
                theta_lo: self.elbows[e].theta_lo,
                theta_hi: self.elbows[e].theta_hi,
            },
            Member::Cap(c) => FaceDomain::Plane {
                plane: self.caps[c].plane,
                half_extent: 2.0 * self.caps[c].radius,
            },
        }
    }

    /// Junction circle at a leg end (normal = leg axis).
    fn leg_circle(&self, l: usize, end: u8) -> Circle3 {
        let leg = &self.legs[l];
        let v = if end == 0 { 0.0 } else { leg.length };
        Circle3::new(
            leg.surf.frame.origin + leg.surf.axis().as_vec() * v,
            leg.surf.axis(),
            leg.surf.radius,
        )
    }

    /// The junction circles bounding member `m`, as `(leg, end)` keys.
    fn junctions_of(&self, m: Member) -> Vec<(usize, u8)> {
        match m {
            Member::Leg(l) => vec![(l, 0), (l, 1)],
            Member::Elbow(e) => vec![(e, 1), (e + 1, 0)],
            Member::Cap(c) => vec![(self.caps[c].leg, self.caps[c].end)],
        }
    }

    /// The member on the other side of junction `(leg, end)`.
    fn across(&self, l: usize, end: u8) -> Member {
        if end == 0 {
            if l == 0 {
                Member::Cap(0)
            } else {
                Member::Elbow(l - 1)
            }
        } else if l + 1 == self.legs.len() {
            Member::Cap(1)
        } else {
            Member::Elbow(l)
        }
    }

    /// Parametric seam points of the two faces meeting at a junction: the
    /// leg's `u = 0` (its `frame.x`) and, for an elbow, the torus's `φ = 0`
    /// outer-equator point.  Both must be vertices of the shared circle chain,
    /// or the periodic DCEL of one face splits the circle where the other does
    /// not.
    fn seam_points_at(&self, l: usize, end: u8) -> Vec<Point3> {
        // Each face is emitted as two halves, so BOTH lines that bound a half
        // are seam lines: u = 0 and u = π on the leg, the outer and the inner
        // equator on the elbow.
        let c = self.leg_circle(l, end);
        let x = self.legs[l].surf.frame.x;
        let mut out = vec![c.frame.origin + x * c.radius, c.frame.origin - x * c.radius];
        if let Member::Elbow(e) = self.across(l, end) {
            let t = &self.elbows[e].surf;
            let w = c.frame.origin - t.frame.origin;
            let h = t.frame.z.dot_vec(w);
            let q = w - t.frame.z.as_vec() * h;
            if let Some(ring) = q.try_normalize() {
                out.push(c.frame.origin + ring * t.minor_radius);
                out.push(c.frame.origin - ring * t.minor_radius);
            }
        }
        out
    }

    fn member_box(&self, m: Member) -> Aabb {
        let mut bb = Aabb::empty();
        match m {
            Member::Leg(l) => {
                let leg = &self.legs[l];
                bb.expand(leg.surf.frame.origin);
                bb.expand(leg.surf.frame.origin + leg.surf.axis().as_vec() * leg.length);
                bb.pad(leg.surf.radius)
            }
            Member::Elbow(e) => {
                let el = &self.elbows[e];
                let steps = 8;
                for k in 0..=steps {
                    let th = el.theta_lo + (el.theta_hi - el.theta_lo) * k as f64 / steps as f64;
                    let dir = el.surf.frame.x * th.cos() + el.surf.frame.y * th.sin();
                    bb.expand(el.surf.frame.origin + dir * el.surf.major_radius);
                }
                bb.pad(el.surf.minor_radius)
            }
            Member::Cap(c) => {
                bb.expand(self.caps[c].centre);
                bb.pad(self.caps[c].radius)
            }
        }
    }

    /// Span coordinate of `p` on member `m`: `(value, lo, hi)`.  The member's
    /// real surface is `lo ≤ value ≤ hi`.
    fn span(&self, m: Member, p: Point3) -> (f64, f64, f64) {
        match m {
            Member::Leg(l) => {
                let leg = &self.legs[l];
                (leg.surf.axis().dot_vec(p - leg.surf.frame.origin), 0.0, leg.length)
            }
            Member::Elbow(e) => {
                let el = &self.elbows[e];
                let t = &el.surf;
                let w = p - t.frame.origin;
                let ang = t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w));
                let sweep = el.theta_hi - el.theta_lo;
                let sgn = if sweep < 0.0 { -1.0 } else { 1.0 };
                let mut rel = ((ang - el.theta_lo) * sgn).rem_euclid(std::f64::consts::TAU);
                // values just before the arc wrap to ~2π: bring them back
                if rel > std::f64::consts::PI + 0.5 * sweep.abs() {
                    rel -= std::f64::consts::TAU;
                }
                (rel, 0.0, sweep.abs())
            }
            Member::Cap(c) => {
                let cap = &self.caps[c];
                let w = p - cap.centre;
                let n = cap.plane.normal();
                let radial = (w - n.as_vec() * n.dot_vec(w)).length();
                (radial, -1.0, cap.radius)
            }
        }
    }

    fn inside_span(&self, m: Member, p: Point3, margin: f64) -> bool {
        let (v, lo, hi) = self.span(m, p);
        v > lo + margin && v < hi - margin
    }

    /// If the step `p → q` leaves member `m`, the fraction along the step
    /// where it does and which end (`true` = past `hi`).
    fn exit_of(&self, m: Member, p: Point3, q: Point3) -> Option<(f64, bool)> {
        let (vp, lo, hi) = self.span(m, p);
        let (vq, _, _) = self.span(m, q);
        const EPS: f64 = 1e-9;
        if vq > hi + EPS {
            let f = if (vq - vp).abs() > 1e-15 { (hi - vp) / (vq - vp) } else { 0.0 };
            Some((f.clamp(0.0, 1.0), true))
        } else if vq < lo - EPS && lo > -0.5 {
            let f = if (vq - vp).abs() > 1e-15 { (lo - vp) / (vq - vp) } else { 0.0 };
            Some((f.clamp(0.0, 1.0), false))
        } else {
            None
        }
    }

    /// The member entered when leaving `m` past `hi` (or `lo`), with the
    /// junction circle crossed and its key.
    fn neighbour(&self, m: Member, hi: bool) -> Option<(Member, Circle3, (usize, u8))> {
        match m {
            Member::Leg(l) => {
                let end = if hi { 1 } else { 0 };
                Some((self.across(l, end), self.leg_circle(l, end), (l, end)))
            }
            Member::Elbow(e) => {
                if hi {
                    Some((Member::Leg(e + 1), self.leg_circle(e + 1, 0), (e + 1, 0)))
                } else {
                    Some((Member::Leg(e), self.leg_circle(e, 1), (e, 1)))
                }
            }
            Member::Cap(c) => {
                let cap = &self.caps[c];
                hi.then(|| (Member::Leg(cap.leg), self.leg_circle(cap.leg, cap.end), (cap.leg, cap.end)))
            }
        }
    }

    /// Broad-phase seeds for a member of this body against a member of another.
    fn seeds(&self, m: Member, other: &Body, mo: Member) -> Vec<Point3> {
        match (m, mo) {
            (Member::Cap(_), _) | (_, Member::Cap(_)) => {
                // sample the cap disk(s); the refinement pulls each sample onto
                // the two surfaces.
                let mut out = Vec::new();
                let mut disk = |cap: &Cap| {
                    out.push(cap.centre);
                    let (x, y) = cap.plane.normal().perp_basis();
                    for k in 0..8 {
                        let a = std::f64::consts::TAU * k as f64 / 8.0;
                        for f in [0.45, 0.85] {
                            out.push(
                                cap.centre
                                    + (x.as_vec() * a.cos() + y.as_vec() * a.sin()) * (cap.radius * f),
                            );
                        }
                    }
                };
                if let Member::Cap(c) = m {
                    disk(&self.caps[c]);
                }
                if let Member::Cap(c) = mo {
                    disk(&other.caps[c]);
                }
                out
            }
            _ => seeds_between(&self.member_tube(m), &other.member_tube(mo), 0.4),
        }
    }

    fn member_tube(&self, m: Member) -> CompositeTube {
        match m {
            Member::Leg(l) => CompositeTube::new().with_leg(self.legs[l].surf, self.legs[l].length),
            Member::Elbow(e) => {
                let el = &self.elbows[e];
                CompositeTube::new().with_elbow(el.surf, el.theta_lo, el.theta_hi)
            }
            Member::Cap(_) => CompositeTube::new(),
        }
    }
}

/// Rotate each leg's cylinder frame so its two seam lines — `u = 0` and
/// `u = π`, the rulings every leg face is split along — run where the fewest
/// trim curves do.
///
/// A leg crossed from above and below carries its windows on top and bottom
/// and is free along both flanks, so one line goes down each flank and neither
/// touches a window.  Where no direction clears both lines (windows all round)
/// the lines cross trim curves, which is legitimate — every crossing gets an
/// exact shared vertex — but a line must not run ALONG a curve or graze the
/// end of one, so the direction with the widest clearance still wins.
///
/// The elbow tori are split along their outer / inner equators, and those
/// meet the leg's rim in two fixed points.  A leg seam landing a hair beside
/// one puts two vertices microns apart on the rim (a sliver edge), so the
/// equator direction is preferred whenever it is about as clear as the best.
fn park_seams(bodies: &mut [Body], loops: &[TracedLoop]) {
    use std::f64::consts::{PI, TAU};
    for (bi, body) in bodies.iter_mut().enumerate() {
        for l in 0..body.legs.len() {
            let surf = body.legs[l].surf;
            let angle = |w: Vec3| surf.frame.y.dot_vec(w).atan2(surf.frame.x.dot_vec(w)).rem_euclid(TAU);
            let mut us: Vec<f64> = Vec::new();
            for lp in loops {
                for piece in &lp.pieces {
                    if piece.a != (bi, Member::Leg(l)) && piece.b != (bi, Member::Leg(l)) {
                        continue;
                    }
                    for &p in &piece.pts {
                        us.push(angle(p - surf.frame.origin));
                    }
                }
            }
            us.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            // angular distance from `u` to the nearest trim-curve sample
            let clear = |u: f64| -> f64 {
                if us.is_empty() {
                    return PI;
                }
                let u = u.rem_euclid(TAU);
                let i = us.partition_point(|&x| x < u);
                let hi = if i < us.len() { us[i] } else { us[0] + TAU };
                let lo = if i > 0 { us[i - 1] } else { us[us.len() - 1] - TAU };
                (hi - u).min(u - lo)
            };
            let score = |u: f64| clear(u).min(clear(u + PI));

            // the adjacent elbows' equator directions, seen from this leg
            let mut equators: Vec<f64> = Vec::new();
            let ends = [
                (l.checked_sub(1), surf.frame.origin),
                (Some(l), surf.frame.origin + surf.axis().as_vec() * body.legs[l].length),
            ];
            for (e, rim_c) in ends {
                let Some(e) = e else { continue };
                let Some(el) = body.elbows.get(e) else { continue };
                let w = rim_c - el.surf.frame.origin;
                let w = w - el.surf.frame.z.as_vec() * el.surf.frame.z.dot_vec(w);
                let w = w - surf.axis().as_vec() * surf.axis().dot_vec(w);
                if w.length() > 1e-9 {
                    equators.push(angle(w));
                }
            }

            const N: usize = 720;
            let (mut free_u, mut free_s) = (0.0_f64, -1.0_f64);
            for k in 0..N {
                let u = PI * k as f64 / N as f64; // u and u + π are the same pair
                let s = score(u);
                if s > free_s {
                    free_s = s;
                    free_u = u;
                }
            }
            let eq = equators
                .iter()
                .map(|&u| (score(u), u))
                .max_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            let seam_u = match eq {
                Some((s, u)) if s >= (0.5 * free_s).min(0.05) && s > 1e-3 => u,
                _ => free_u,
            };
            if std::env::var("CADCORE_DUMP_FUSE").is_ok() {
                eprintln!(
                    "[fuse][seam] body {bi} leg {l}: {} piece points, seam_u {:.3} (clearance {:.3} rad; best free {:.3} at {:.3}; equators {:?})",
                    us.len(),
                    seam_u,
                    score(seam_u),
                    free_s,
                    free_u,
                    equators
                );
            }
            if seam_u.abs() < 1e-15 {
                continue;
            }
            let (sn, cs) = seam_u.sin_cos();
            let Some(x) = cadcore_math::UnitVec3::try_from_vec(surf.frame.x * cs + surf.frame.y * sn) else {
                continue;
            };
            let Some(y) = cadcore_math::UnitVec3::try_from_vec(surf.frame.z.as_vec().cross(x.as_vec())) else {
                continue;
            };
            body.legs[l].surf = cadcore_geom::CylSurf {
                frame: cadcore_math::Frame3 { origin: surf.frame.origin, x, y, z: surf.frame.z },
                radius: surf.radius,
            };
        }
    }
}

// ── Traced loops ─────────────────────────────────────────────────────────────

/// One run of the boundary curve lying on face `a` of one body and face `b`
/// of another.  The SAME point list is imprinted on both faces.
#[derive(Debug, Clone)]
struct Piece {
    a: (usize, Member),
    b: (usize, Member),
    pts: Vec<Point3>,
}

/// A point where the boundary curve crosses a junction circle.
#[derive(Debug, Clone, Copy)]
struct Joint {
    key: JunctionKey,
    p: Point3,
}

#[derive(Debug, Clone)]
struct TracedLoop {
    pieces: Vec<Piece>,
    joints: Vec<Joint>,
}

/// Coarse spatial hash of every traced point, for seed de-duplication.
struct PointGrid {
    cell: f64,
    map: HashMap<(i64, i64, i64), Vec<Point3>>,
}

impl PointGrid {
    fn new(cell: f64) -> Self {
        Self { cell: cell.max(1e-6), map: HashMap::new() }
    }
    fn key(&self, p: Point3) -> (i64, i64, i64) {
        (
            (p.x / self.cell).floor() as i64,
            (p.y / self.cell).floor() as i64,
            (p.z / self.cell).floor() as i64,
        )
    }
    fn insert(&mut self, p: Point3) {
        self.map.entry(self.key(p)).or_default().push(p);
    }
    fn near(&self, p: Point3, r: f64) -> bool {
        let k = self.key(p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    if let Some(v) = self.map.get(&(k.0 + dx, k.1 + dy, k.2 + dz)) {
                        if v.iter().any(|&q| (q - p).length() < r) {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }
}

fn trace_pair(
    bodies: &[Body],
    i: usize,
    j: usize,
    opts: &FuseOptions,
    loops: &mut Vec<TracedLoop>,
    known: &mut PointGrid,
    report: &mut FuseReport,
) {
    let (bi, bj) = (&bodies[i], &bodies[j]);
    for mi in bi.members() {
        let box_i = bi.member_box(mi);
        for mj in bj.members() {
            if !box_i.overlaps(&bj.member_box(mj)) {
                continue;
            }
            let si = bi.surface(mi);
            let sj = bj.surface(mj);
            for seed in bi.seeds(mi, bj, mj) {
                let Some(p0) = intersection_point_near(&si, &sj, seed, 1e-12) else {
                    continue;
                };
                // start strictly inside both members: a seed on a junction
                // would make the first step ambiguous.
                if !bi.inside_span(mi, p0, 1e-6) || !bj.inside_span(mj, p0, 1e-6) {
                    continue;
                }
                if known.near(p0, opts.step * 1.5) {
                    continue;
                }
                match trace_boundary(bodies, i, j, mi, mj, p0, opts) {
                    Some(lp) => {
                        for piece in &lp.pieces {
                            for &p in &piece.pts {
                                known.insert(p);
                            }
                        }
                        loops.push(lp);
                    }
                    None => {
                        report.aborted += 1;
                        // remember the failed start so we do not retry it from
                        // a neighbouring seed forever.
                        known.insert(p0);
                    }
                }
            }
        }
    }
}

/// March along `∂A ∩ ∂B` from `start` (strictly inside members `mi`/`mj`)
/// until the curve closes, switching members at junction circles.
fn trace_boundary(
    bodies: &[Body],
    i: usize,
    j: usize,
    mi: Member,
    mj: Member,
    start: Point3,
    opts: &FuseOptions,
) -> Option<TracedLoop> {
    let (bi, bj) = (&bodies[i], &bodies[j]);
    let mut m_i = mi;
    let mut m_j = mj;
    let mut s_i = bi.surface(m_i);
    let mut s_j = bj.surface(m_j);

    let mut pieces: Vec<Piece> = Vec::new();
    let mut joints: Vec<Joint> = Vec::new();
    let mut cur = Piece { a: (i, m_i), b: (j, m_j), pts: vec![start] };
    let mut p = start;
    let mut dir: Option<Vec3> = None;
    let mut travelled = 0.0;
    let mut consecutive_switches = 0usize;

    let dbg = std::env::var("CADCORE_DUMP_TRACE").is_ok();
    let verbose = std::env::var("CADCORE_DUMP_TRACE").map_or(false, |v| v == "verbose");
    let mut min_to_start = f64::MAX;
    let mut switches = 0usize;
    let mut events = 0usize;
    let mut steps_total = 0usize;
    macro_rules! abort {
        ($why:expr) => {{
            if dbg {
                eprintln!("[trace][abort] {} at ({:.4},{:.4},{:.4}) members {:?}/{:?} pieces so far {}", $why, p.x, p.y, p.z, m_i, m_j, pieces.len());
            }
            return None;
        }};
    }
    for _ in 0..opts.max_steps {
        let (Some(na), Some(nb)) = (s_i.normal_at(p), s_j.normal_at(p)) else { abort!("no normal") };
        let cr = na.cross(nb);
        let cl = cr.length();
        if cl < 1e-5 {
            abort!("tangency");
        }
        let mut t = cr * (1.0 / cl);
        if let Some(d) = dir {
            if t.dot(d) < 0.0 {
                t = -t;
            }
        }
        let h = opts.step * cl.clamp(0.15, 1.0);
        // predictor–corrector, halving the step while the corrector struggles
        let mut q = None;
        let mut hh = h;
        for _ in 0..5 {
            let (c, res) = crate::geom::refine::refine_pair(p + t * hh, &s_i, &s_j, 1e-13, 30);
            if res < 1e-9 && (c - p).length() > hh * 0.05 && (c - p).dot(t) > 0.0 {
                q = Some(c);
                break;
            }
            hh *= 0.5;
        }
        let Some(q) = q else { abort!("corrector failed") };

        let ex_i = bi.exit_of(m_i, p, q);
        let ex_j = bj.exit_of(m_j, p, q);
        let first = match (ex_i, ex_j) {
            (None, None) => None,
            (Some((f, hi)), None) => Some((true, f, hi)),
            (None, Some((f, hi))) => Some((false, f, hi)),
            (Some((fi, hi_i)), Some((fj, hi_j))) => {
                if fi <= fj {
                    Some((true, fi, hi_i))
                } else {
                    Some((false, fj, hi_j))
                }
            }
        };

        // closure: the step p → q passes within a step of where we started
        // (tested against the SEGMENT, so an overshoot cannot skip it; and at
        // the nominal step, so a corrector-shortened step cannot shrink the
        // window).
        let seg_dist_to_start = {
            let d = q - p;
            let l2 = d.dot(d);
            let t = if l2 > 1e-30 { ((start - p).dot(d) / l2).clamp(0.0, 1.0) } else { 0.0 };
            (p + d * t - start).length()
        };
        if dbg && travelled > 3.0 * opts.step {
            min_to_start = min_to_start.min(seg_dist_to_start);
        }
        // …and only while marching on the SAME face pair the loop started on:
        // the start lies on that pair, so coming within reach of it on another
        // pair means a junction still lies between here and the start (the
        // start sitting a hair past a rim, say) — closing now would skip that
        // joint and leave the piece ending off its surfaces.
        let on_start_pair = (m_i, m_j) == (mi, mj);
        if travelled > 3.0 * opts.step && seg_dist_to_start < 0.75 * opts.step && first.is_none() && on_start_pair {
            cur.pts.push(start);
            pieces.push(cur);
            return finalize(pieces, joints);
        }
        // A lap is complete but the loop cannot close HERE: the start sits
        // within a step of a junction, so the pass by it happens on another
        // face pair (or on the very step that exits).  Marching on would only
        // circle the same loop until the step budget runs out.  Re-seed in the
        // middle of the longest piece of the lap just walked — as far from
        // every junction as this loop allows — and trace it from there.
        if travelled > 3.0 * opts.step && seg_dist_to_start < 0.75 * opts.step {
            let best = pieces
                .iter()
                .filter(|pc| pc.pts.len() >= 5)
                .max_by(|a, b| a.pts.len().cmp(&b.pts.len()));
            if let Some(pc) = best {
                let mid = pc.pts[pc.pts.len() / 2];
                if (mid - start).length() > 2.0 * opts.step && pc.a.0 == i && pc.b.0 == j {
                    if dbg {
                        eprintln!(
                            "[trace][reseed] lap done without closure at ({:.4},{:.4},{:.4}); restarting from ({:.4},{:.4},{:.4}) on {:?}/{:?}",
                            p.x, p.y, p.z, mid.x, mid.y, mid.z, pc.a.1, pc.b.1
                        );
                    }
                    return trace_boundary(bodies, i, j, pc.a.1, pc.b.1, mid, opts);
                }
            }
        }
        steps_total += 1;
        if verbose && events < 160 && steps_total % 25 == 0 {
            events += 1;
            eprintln!("[trace][step] n={steps_total} ({:.4},{:.4},{:.4}) {:?}/{:?} travelled {:.3}", q.x, q.y, q.z, m_i, m_j, travelled);
        }
        let Some((on_i, frac, hi)) = first else {
            // plain step
            consecutive_switches = 0;
            cur.pts.push(q);
            travelled += (q - p).length();
            dir = Some(q - p);
            p = q;
            continue;
        };

        // ── junction crossing ────────────────────────────────────────────
        consecutive_switches += 1;
        switches += 1;
        if consecutive_switches > 6 {
            abort!("junction ping-pong");
        }
        let approx = p.lerp(q, frac);
        let (body, m_cur, other_surf) = if on_i { (bi, m_i, s_j) } else { (bj, m_j, s_i) };
        let Some((next_m, circle, key)) = body.neighbour(m_cur, hi) else { abort!("no neighbour") };
        let Some(jpt) = joint_on_circle(&circle, &other_surf, approx, hh.max(opts.step)) else { abort!("joint not found") };
        let jkey = JunctionKey { body: if on_i { i } else { j }, leg: key.0, end: key.1 };

        // A step point landing a hair before the joint would leave a sliver
        // edge behind; the joint replaces it.
        if cur.pts.len() >= 2 && (jpt - *cur.pts.last().unwrap()).length() < 0.25 * opts.step {
            cur.pts.pop();
        }
        if (jpt - *cur.pts.last().unwrap()).length() > 1e-12 {
            cur.pts.push(jpt);
        }
        if cur.pts.len() >= 2 {
            pieces.push(cur);
        }
        joints.push(Joint { key: jkey, p: jpt });
        if verbose && events < 160 {
            events += 1;
            eprintln!(
                "[trace][switch] n={steps_total} side {} {:?} → {:?} (hi={hi}) at ({:.4},{:.4},{:.4}) frac {frac:.3}",
                if on_i { "i" } else { "j" }, m_cur, next_m, jpt.x, jpt.y, jpt.z
            );
        }
        if on_i {
            m_i = next_m;
            s_i = bi.surface(m_i);
        } else {
            m_j = next_m;
            s_j = bj.surface(m_j);
        }
        cur = Piece { a: (i, m_i), b: (j, m_j), pts: vec![jpt] };
        travelled += (jpt - p).length();
        p = jpt;

        // new marching direction: enter the member we just switched to
        let (Some(na), Some(nb)) = (s_i.normal_at(p), s_j.normal_at(p)) else { abort!("no normal after switch") };
        let cr = na.cross(nb);
        let cl = cr.length();
        if cl < 1e-5 {
            // Two caps in one plane (filaments starting flush on the same
            // boundary plane) — a coplanar-face Boolean, not a curve to march.
            abort!("tangency after switch (coplanar caps?)");
        }
        let t0 = cr * (1.0 / cl);
        let probe_len = 0.5 * opts.step;
        let entering = |t: Vec3| -> (bool, f64, f64, f64, f64, f64) {
            let (c, res) = crate::geom::refine::refine_pair(p + t * probe_len, &s_i, &s_j, 1e-13, 30);
            let b = if on_i { bi } else { bj };
            let m = if on_i { m_i } else { m_j };
            let (v, lo, hi) = b.span(m, c);
            let (v0, _, _) = b.span(m, p);
            if res > 1e-7 {
                return (false, res, v, lo, hi, v0);
            }
            let ok = v > lo - 1e-9 && v < hi + 1e-9 && {
                // and genuinely moving INTO the span, not along its boundary
                let inward = if lo > -0.5 && (v0 - lo).abs() < (v0 - hi).abs() { v > v0 } else { v < v0 };
                inward || (v0 - v).abs() < 1e-12
            };
            (ok, res, v, lo, hi, v0)
        };
        let (fwd, fres, fv, flo, fhi, fv0) = entering(t0);
        let (bwd, bres, bv, _, _, _) = entering(-t0);
        if verbose && events < 160 {
            events += 1;
            eprintln!(
                "[trace][enter] {:?} span [{flo:.4},{fhi:.4}] v0={fv0:.5}  fwd: ok={fwd} res={fres:.1e} v={fv:.5}  bwd: ok={bwd} res={bres:.1e} v={bv:.5}",
                if on_i { m_i } else { m_j }
            );
        }
        dir = Some(match (fwd, bwd) {
            (true, false) => t0,
            (false, true) => -t0,
            _ => match dir {
                Some(d) if t0.dot(d) < 0.0 => -t0,
                _ => t0,
            },
        });
    }
    if dbg {
        eprintln!(
            "[trace][abort] step budget exhausted near ({:.4},{:.4},{:.4}); start ({:.4},{:.4},{:.4}) members {:?}/{:?}, {} switches, {} pieces, travelled {:.2}, closest approach to start {:.2e}",
            p.x, p.y, p.z, start.x, start.y, start.z, mi, mj, switches, pieces.len(), travelled, min_to_start
        );
    }
    None
}

/// Merge the closing piece into the opening one (same face pair, shared
/// `start` point) and drop degenerate pieces.
fn finalize(mut pieces: Vec<Piece>, joints: Vec<Joint>) -> Option<TracedLoop> {
    if pieces.is_empty() {
        return None;
    }
    if pieces.len() >= 2 {
        let last = pieces.last().unwrap();
        let first = &pieces[0];
        if last.a == first.a && last.b == first.b {
            let last = pieces.pop().unwrap();
            let mut merged = last.pts;
            merged.extend(pieces[0].pts.iter().skip(1).copied());
            pieces[0].pts = merged;
        }
    }
    for piece in &mut pieces {
        piece.pts.dedup_by(|b, a| (*b - *a).length() < 1e-12);
    }
    pieces.retain(|p| p.pts.len() >= 2);
    if pieces.is_empty() {
        return None;
    }
    Some(TracedLoop { pieces, joints })
}

/// The point of `circle` that lies on `surf`, nearest to `approx`.
///
/// `approx` is where the straight step crossed the junction plane — within a
/// sagitta of the true point.  A 1-D golden-section minimisation of the
/// distance to `surf` along the circle, bracketed by `reach` (mm) each way,
/// lands on the exact joint; it is what both the junction chain and the piece
/// will carry.
fn joint_on_circle(circle: &Circle3, surf: &AnalyticSurface, approx: Point3, reach: f64) -> Option<Point3> {
    let theta0 = circle_angle(circle, approx);
    let f = |th: f64| surf.distance(circle.point_at(th));
    let mut delta = (3.0 * reach / circle.radius.max(1e-9)).max(0.05);
    for _ in 0..3 {
        let th = golden_min(theta0 - delta, theta0 + delta, 80, f);
        let q = circle.point_at(th);
        if surf.distance(q) < 1e-7 {
            return Some(q);
        }
        delta *= 4.0;
    }
    None
}

fn circle_angle(c: &Circle3, p: Point3) -> f64 {
    let w = p - c.frame.origin;
    c.frame.y.dot_vec(w).atan2(c.frame.x.dot_vec(w))
}

fn golden_min(mut a: f64, mut b: f64, iters: usize, f: impl Fn(f64) -> f64) -> f64 {
    let gr = (5f64.sqrt() - 1.0) / 2.0;
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..iters {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = f(d);
        }
    }
    0.5 * (a + b)
}

/// Thin a traced loop to roughly `spacing` between kept points, keeping every
/// piece's endpoints (they are joints shared with the junction chains).
fn decimate_loop(lp: &mut TracedLoop, spacing: f64) {
    for piece in &mut lp.pieces {
        let n = piece.pts.len();
        if n <= 3 {
            continue;
        }
        let len: f64 = piece.pts.windows(2).map(|w| (w[1] - w[0]).length()).sum();
        let keep = ((len / spacing.max(1e-6)).ceil() as usize + 1).clamp(3, n);
        if keep >= n {
            continue;
        }
        let mut out: Vec<Point3> = Vec::with_capacity(keep);
        let total_len: f64 = piece.pts.windows(2).map(|w| (w[1] - w[0]).length()).sum();
        let min_gap = 0.3 * total_len / (keep.max(2) - 1) as f64;
        for k in 0..keep {
            let idx = (k as f64 * (n - 1) as f64 / (keep - 1) as f64).round() as usize;
            let p = piece.pts[idx.min(n - 1)];
            let is_last = k + 1 == keep;
            // never let an interior sample sit within a sliver of the endpoint
            if !is_last && out.last().map_or(false, |&q: &Point3| (q - p).length() < min_gap) {
                continue;
            }
            if is_last && out.len() >= 2 && (out[out.len() - 1] - p).length() < min_gap {
                out.pop();
            }
            if out.last().map_or(true, |&q: &Point3| (q - p).length() > 1e-12) {
                out.push(p);
            }
        }
        piece.pts = out;
    }
}

/// Insert the exact point where a piece crosses the parametric seam of its
/// carrier `own` (a cylinder's `u = 0` half-plane, a torus's `φ = 0` outer
/// equator).  The periodic DCEL of the face on `own` has to split the chain
/// there anyway; giving it — and the face on `other`, which shares this very
/// point list — an exact vertex keeps both subdivisions identical.
fn split_piece_at_seam(
    pts: &mut Vec<Point3>,
    own: &AnalyticSurface,
    other: &AnalyticSurface,
    side: f64,
    report: &mut FuseReport,
) {
    // `side` = +1: the seam proper (the u = 0 ruling of a cylinder, the outer
    // equator of a torus); −1: the line half a period away (u = π, the inner
    // equator) — the second boundary of the two half faces every wrapped face
    // is emitted as.
    let (value, guard): (Box<dyn Fn(Point3) -> f64>, Box<dyn Fn(Point3) -> bool>) = match own {
        AnalyticSurface::Cylinder(c) => {
            let c = *c;
            (
                Box::new(move |p: Point3| c.frame.y.dot_vec(p - c.frame.origin)),
                Box::new(move |p: Point3| side * c.frame.x.dot_vec(p - c.frame.origin) > 0.0),
            )
        }
        AnalyticSurface::Torus(t) => {
            let t = *t;
            (
                Box::new(move |p: Point3| t.frame.z.dot_vec(p - t.frame.origin)),
                Box::new(move |p: Point3| {
                    let w = p - t.frame.origin;
                    let h = t.frame.z.dot_vec(w);
                    side * ((w - t.frame.z.as_vec() * h).length() - t.major_radius) > 0.0
                }),
            )
        }
        AnalyticSurface::Plane(_) => return,
    };
    let mut out: Vec<Point3> = Vec::with_capacity(pts.len() + 2);
    let mut k = 0usize;
    while k < pts.len() {
        let a = pts[k];
        out.push(a);
        if k + 1 >= pts.len() {
            break;
        }
        let b = pts[k + 1];
        let (va, vb) = (value(a), value(b));
        if (va > 0.0) == (vb > 0.0) || (va - vb).abs() < 1e-18 {
            k += 1;
            continue;
        }
        let t = va / (va - vb);
        let lerp = a.lerp(b, t);
        if !guard(lerp) {
            k += 1;
            continue;
        }
        let exact = match own {
            AnalyticSurface::Cylinder(c) => {
                // seam ruling line: o + r·x + axis·v
                let base = c.frame.origin + c.frame.x * (side * c.radius);
                let v0 = c.axis().dot_vec(lerp - c.frame.origin);
                let reach = (b - a).length().max(1e-3);
                let f = |v: f64| other.distance(base + c.axis().as_vec() * v);
                let v = golden_min(v0 - reach, v0 + reach, 80, f);
                let q = base + c.axis().as_vec() * v;
                (other.distance(q) < 1e-7).then_some(q)
            }
            AnalyticSurface::Torus(t) => {
                let w = lerp - t.frame.origin;
                let th0 = t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w));
                let reach = ((b - a).length() / t.major_radius.max(1e-9)).max(1e-3);
                let eq = |th: f64| {
                    t.frame.origin + (t.frame.x * th.cos() + t.frame.y * th.sin()) * (t.major_radius + side * t.minor_radius)
                };
                let f = |th: f64| other.distance(eq(th));
                let th = golden_min(th0 - reach, th0 + reach, 80, f);
                let q = eq(th);
                (other.distance(q) < 1e-7).then_some(q)
            }
            AnalyticSurface::Plane(_) => None,
        };
        match exact {
            Some(q) => {
                // An existing INTERIOR vertex a hair off the seam is moved onto
                // it instead of growing a sliver next to it; piece ends are
                // joints shared with junction chains and stay put.
                let near_a = (q - a).length() < 1e-6;
                let near_b = (q - b).length() < 1e-6;
                if near_a {
                    if k > 0 {
                        *out.last_mut().unwrap() = q;
                    }
                } else if near_b {
                    if k + 2 < pts.len() {
                        pts[k + 1] = q;
                    }
                } else {
                    out.push(q);
                }
            }
            None => report
                .notes
                .push(format!("seam crossing near {lerp:?} could not be refined; split left to the DCEL")),
        }
        k += 1;
    }
    *pts = out;
}

/// Minimum distance between segments `a0–a1` and `b0–b1`.
fn segment_distance(a0: Point3, a1: Point3, b0: Point3, b1: Point3) -> f64 {
    let (d1, d2, r) = (a1 - a0, b1 - b0, a0 - b0);
    let (a, e, f) = (d1.dot(d1), d2.dot(d2), d2.dot(r));
    let (mut s, mut t);
    if a <= 1e-18 && e <= 1e-18 {
        return r.length();
    }
    if a <= 1e-18 {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = d1.dot(r);
        if e <= 1e-18 {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = d1.dot(d2);
            let den = a * e - b * b;
            s = if den > 1e-18 { ((b * f - c * e) / den).clamp(0.0, 1.0) } else { 0.0 };
            t = (b * s + f) / e;
            if t < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            }
        }
    }
    ((a0 + d1 * s) - (b0 + d2 * t)).length()
}

// ── Junction chains ──────────────────────────────────────────────────────────

/// A junction circle with the vertices every consumer must carry.
struct JunctionChain {
    circle: Circle3,
    mandatory: Vec<Point3>,
}

impl JunctionChain {
    /// `n` uniform samples merged with the mandatory vertices (a uniform
    /// sample too close to a mandatory one is dropped so no sliver segment
    /// appears), sorted by angle and closed.
    fn sample(&self, n: usize) -> Vec<Point3> {
        use std::f64::consts::TAU;
        let n = n.max(8);
        let step = TAU / n as f64;
        let mut mand: Vec<(f64, Point3)> = self
            .mandatory
            .iter()
            .map(|&p| (circle_angle(&self.circle, p).rem_euclid(TAU), p))
            .collect();
        mand.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        mand.dedup_by(|b, a| (b.1 - a.1).length() < 1e-9);
        let mut all: Vec<(f64, Point3)> = Vec::with_capacity(n + mand.len());
        for k in 0..n {
            let th = step * k as f64;
            let close = mand.iter().any(|(a, _)| {
                let d = (a - th).abs();
                d.min(TAU - d) < 0.3 * step
            });
            if !close {
                all.push((th, self.circle.point_at(th)));
            }
        }
        all.extend(mand);
        all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut pts: Vec<Point3> = all.into_iter().map(|(_, p)| p).collect();
        if let Some(&first) = pts.first() {
            pts.push(first);
        }
        pts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::UnionConfig;
    use crate::validate::validate_shell;
    use crate::Diagnostics;

    fn gate(brep: &BRep, shell: ShellId) -> (usize, usize, usize, f64) {
        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let r = validate_shell(brep, shell, &cfg, &mut diag);
        for w in r.wire_violations.iter().take(6) {
            eprintln!("  wire: {w}");
        }
        for d in r.distance_violations.iter().take(4) {
            eprintln!("  dist: {d}");
        }
        (
            r.wire_violations.len(),
            r.distance_violations.len(),
            r.manifold_violations.len(),
            r.max_edge_deviation,
        )
    }

    fn components(brep: &BRep, shell: ShellId) -> usize {
        let faces = &brep.shells[shell].faces;
        let mut by_edge: HashMap<cadcore_topo::EdgeId, Vec<usize>> = HashMap::new();
        for (i, &fid) in faces.iter().enumerate() {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let st = brep.loops[lid].start;
                let mut c = st;
                loop {
                    by_edge.entry(brep.coedges[c].edge).or_default().push(i);
                    c = brep.coedges[c].next;
                    if c == st {
                        break;
                    }
                }
            }
        }
        let mut parent: Vec<usize> = (0..faces.len()).collect();
        fn find(p: &mut [usize], x: usize) -> usize {
            let mut r = x;
            while p[r] != r {
                r = p[r];
            }
            let mut c = x;
            while p[c] != c {
                let n = p[c];
                p[c] = r;
                c = n;
            }
            r
        }
        for fs in by_edge.values() {
            for w in fs.windows(2) {
                let (a, b) = (find(&mut parent, w[0]), find(&mut parent, w[1]));
                parent[a] = b;
            }
        }
        (0..faces.len()).filter(|&i| find(&mut parent, i) == i).count()
    }

    fn dump_open(brep: &BRep, shell: ShellId, max: usize) {
        let rep = crate::validate::manifold::check(brep, shell);
        for v in rep.violations.iter().take(max) {
            eprintln!("  {v}");
        }
    }

    /// The pair a real 3 mm two-layer scaffold fuses: free ends under rows,
    /// U-turn connectors parallel over outer legs, coaxial corner elbows.
    fn real_pair() -> (Filament, Filament) {
        let r = 0.275;
        let (z0, z1) = (0.35, 0.70);
        let a = Filament::serpentine(
            &[
                Point3::new(0.3, 0.3, z0),
                Point3::new(0.3, 2.7, z0),
                Point3::new(1.5, 2.7, z0),
                Point3::new(1.5, 0.3, z0),
                Point3::new(2.7, 0.3, z0),
                Point3::new(2.7, 2.7, z0),
            ],
            r,
            0.5,
        );
        let b = Filament::serpentine(
            &[
                Point3::new(0.3, 0.3, z1),
                Point3::new(2.7, 0.3, z1),
                Point3::new(2.7, 1.5, z1),
                Point3::new(0.3, 1.5, z1),
                Point3::new(0.3, 2.7, z1),
                Point3::new(2.7, 2.7, z1),
            ],
            r,
            0.5,
        );
        (a, b)
    }

    /// Plain transversal crossings: the simplest case, every window strictly
    /// handled must still come out watertight and fused.
    #[test]
    fn two_crossing_serpentines_watertight() {
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
        let (shell, rep) = fuse_filaments_with(&mut brep, &[a, b], &FuseOptions::default());
        eprintln!("{rep:?}");
        let (w, d, m, dev) = gate(&brep, shell);
        dump_open(&brep, shell, 12);
        assert_eq!(rep.loops, 4, "four transversal crossings");
        assert_eq!(m, 0, "manifold");
        assert_eq!(components(&brep, shell), 1, "one solid");
        assert_eq!(w, 0, "wires");
        assert_eq!(d, 0, "distance (max {:.2} um)", dev * 1e3);
    }

    /// A row running over a FREE END: the window straddles the cap rim, so the
    /// band gets a notch, the cap a bite, and the crosser a window closed by the
    /// cap chord.
    #[test]
    fn row_over_free_end_watertight() {
        let r = 0.275;
        let leg = Filament::serpentine(
            &[Point3::new(0.0, -2.0, 0.0), Point3::new(0.0, 0.0, 0.0)],
            r,
            0.5,
        );
        let row = Filament::serpentine(
            &[Point3::new(-2.0, 0.0, 0.35), Point3::new(2.0, 0.0, 0.35)],
            r,
            0.5,
        );
        let mut brep = BRep::new();
        let (shell, rep) = fuse_filaments_with(&mut brep, &[leg, row], &FuseOptions::default());
        eprintln!("{rep:?}");
        let (w, d, m, dev) = gate(&brep, shell);
        dump_open(&brep, shell, 12);
        assert_eq!(rep.loops, 1);
        assert!(rep.joints >= 2, "the loop must cross the rim twice");
        assert_eq!(m, 0, "manifold");
        assert_eq!(components(&brep, shell), 1);
        assert_eq!(w, 0, "wires");
        assert_eq!(d, 0, "distance (max {:.2} um)", dev * 1e3);
    }

    /// Two parallel legs lying on each other: the surfaces meet in two ruling
    /// lines closed over the caps — no closed loop on any single member.
    #[test]
    fn parallel_legs_watertight() {
        let r = 0.275;
        let a = Filament::serpentine(&[Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 0.0, 0.0)], r, 0.5);
        let b = Filament::serpentine(&[Point3::new(0.5, 0.0, 0.35), Point3::new(2.5, 0.0, 0.35)], r, 0.5);
        let mut brep = BRep::new();
        let (shell, rep) = fuse_filaments_with(&mut brep, &[a, b], &FuseOptions::default());
        eprintln!("{rep:?}");
        let (w, d, m, dev) = gate(&brep, shell);
        dump_open(&brep, shell, 12);
        assert_eq!(m, 0, "manifold");
        assert_eq!(components(&brep, shell), 1);
        assert_eq!(w, 0, "wires");
        assert_eq!(d, 0, "distance (max {:.2} um)", dev * 1e3);
    }

    /// A honeycomb bond ending INSIDE the row it bridges (T-junction): the
    /// bond's cap and the rim end of its band are swallowed whole.
    #[test]
    fn bond_ending_inside_row_watertight() {
        let r = 0.275;
        let row = Filament::serpentine(&[Point3::new(-2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)], r, 0.5);
        let bond = Filament::serpentine(&[Point3::new(0.0, -2.0, 0.3), Point3::new(0.0, 0.05, 0.3)], r, 0.5);
        let mut brep = BRep::new();
        let (shell, rep) = fuse_filaments_with(&mut brep, &[row, bond], &FuseOptions::default());
        eprintln!("{rep:?}");
        let (w, d, m, dev) = gate(&brep, shell);
        dump_open(&brep, shell, 12);
        assert_eq!(m, 0, "manifold");
        assert_eq!(components(&brep, shell), 1);
        assert_eq!(w, 0, "wires");
        assert_eq!(d, 0, "distance (max {:.2} um)", dev * 1e3);
    }

    /// A leg crossing over a U-turn elbow (torus × cylinder).
    #[test]
    fn leg_over_elbow_watertight() {
        let r = 0.275;
        let a = Filament::serpentine(
            &[
                Point3::new(-2.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
            ],
            r,
            0.5,
        );
        // b starts INSIDE a's leg span (x = -1.9, not -2.0): two caps lying in
        // one plane is a separate degenerate case (see `FuseReport::notes`).
        let b = Filament::serpentine(
            &[Point3::new(-1.9, 0.15, 0.35), Point3::new(2.0, 0.15, 0.35)],
            r,
            0.5,
        );
        let mut brep = BRep::new();
        let (shell, rep) = fuse_filaments_with(&mut brep, &[a, b], &FuseOptions::default());
        eprintln!("{rep:?}");
        let (w, d, m, dev) = gate(&brep, shell);
        dump_open(&brep, shell, 12);
        assert_eq!(m, 0, "manifold");
        assert_eq!(components(&brep, shell), 1);
        assert_eq!(w, 0, "wires");
        assert_eq!(d, 0, "distance (max {:.2} um)", dev * 1e3);
    }

    /// Every point of every edge in the result, sampled.
    fn shell_points(brep: &BRep, shell: ShellId) -> Vec<Point3> {
        let mut out = Vec::new();
        for &fid in &brep.shells[shell].faces {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let Some(lp) = brep.loops.get(lid) else { continue };
                let st = lp.start;
                let mut c = st;
                loop {
                    let ce = &brep.coedges[c];
                    out.extend(crate::boolean::aabb::sample_edge(brep, ce.edge));
                    c = ce.next;
                    if c == st {
                        break;
                    }
                }
            }
        }
        out
    }

    /// Both terminals off the fused real pair — the export order the process
    /// needs (sinter, then chop), on a shell that is a TRUE union.
    #[test]
    fn cut_both_terminals_off_the_fused_pair() {
        use crate::boolean::half_space_cut;
        use cadcore_math::UnitVec3;
        let (a, b) = real_pair();
        let mut fused = BRep::new();
        let (shell, _) = fuse_filaments_with(&mut fused, &[a, b], &FuseOptions::default());
        let faces = fused.shells[shell].faces.clone();

        let lo = Plane3::from_origin_normal(Point3::new(0.0, 1.0, 0.0), UnitVec3::Y);
        let first = half_space_cut(&fused, &faces, &lo, 1e-6).expect("first cut");
        let (w, d, m, dev) = gate(&first.brep, first.shell);
        dump_open(&first.brep, first.shell, 12);
        eprintln!("first cut: caps={} dropped={} wire={w} dist={d} manifold={m} maxdev={:.2}um", first.caps, first.dropped, dev * 1e3);
        assert!(first.caps >= 3, "three legs are severed");
        assert_eq!(m, 0, "manifold after the first terminal");
        assert_eq!(w, 0, "wires after the first terminal");
        let deepest = shell_points(&first.brep, first.shell)
            .into_iter()
            .map(|p| lo.signed_distance(p))
            .fold(f64::MAX, f64::min);
        assert!(deepest > -1e-3, "material past the first plane: {deepest:.3e}");

        let faces2 = first.brep.shells[first.shell].faces.clone();
        let hi = Plane3::from_origin_normal(Point3::new(0.0, 2.0, 0.0), -UnitVec3::Y);
        let second = half_space_cut(&first.brep, &faces2, &hi, 1e-6).expect("second cut");
        let (w, d, m, dev) = gate(&second.brep, second.shell);
        dump_open(&second.brep, second.shell, 12);
        eprintln!("second cut: caps={} dropped={} wire={w} dist={d} manifold={m} maxdev={:.2}um", second.caps, second.dropped, dev * 1e3);
        assert_eq!(m, 0, "manifold after the second terminal");
        assert_eq!(w, 0, "wires after the second terminal");
        assert_eq!(components(&second.brep, second.shell), 1, "still one solid");
        let deepest = shell_points(&second.brep, second.shell)
            .into_iter()
            .map(|p| hi.signed_distance(p))
            .fold(f64::MAX, f64::min);
        assert!(deepest > -1e-3, "material past the second plane: {deepest:.3e}");
    }

    /// A plane through the ELBOW bulge (torus × plane sections) — the offset a
    /// user actually types is not tuned to the fillets.
    #[test]
    fn cut_through_the_elbows_of_the_fused_pair() {
        use crate::boolean::half_space_cut;
        use cadcore_math::UnitVec3;
        let (a, b) = real_pair();
        let mut fused = BRep::new();
        let (shell, _) = fuse_filaments_with(&mut fused, &[a, b], &FuseOptions::default());
        let faces = fused.shells[shell].faces.clone();
        // y = 0.85 cuts the elbow tubes (they bulge to y ≈ 0.87 past the 0.8
        // tangent) and the straight Y-legs.
        let lo = Plane3::from_origin_normal(Point3::new(0.0, 0.85, 0.0), UnitVec3::Y);
        let cut = half_space_cut(&fused, &faces, &lo, 1e-6).expect("cut");
        let (w, d, m, dev) = gate(&cut.brep, cut.shell);
        dump_open(&cut.brep, cut.shell, 12);
        eprintln!("elbow cut: caps={} dropped={} wire={w} dist={d} manifold={m} maxdev={:.2}um", cut.caps, cut.dropped, dev * 1e3);
        assert_eq!(m, 0, "manifold");
        assert_eq!(w, 0, "wires");
    }

    /// The real thing: the two-layer woodpile with every outline feature.
    #[test]
    fn real_scaffold_pair_passes_the_gate() {
        let (a, b) = real_pair();
        let mut brep = BRep::new();
        let (shell, rep) = fuse_filaments_with(&mut brep, &[a, b], &FuseOptions::default());
        eprintln!("{rep:?}");
        let (w, d, m, dev) = gate(&brep, shell);
        dump_open(&brep, shell, 20);
        assert_eq!(m, 0, "manifold");
        assert_eq!(components(&brep, shell), 1, "one solid");
        assert_eq!(w, 0, "wires");
        assert_eq!(d, 0, "distance (max {:.2} um)", dev * 1e3);
    }
}


#[cfg(test)]
mod export_tests {
    use super::*;
    use crate::boolean::half_space_cut;
    use cadcore_math::UnitVec3;

    /// A woodpile scaffold as the generator lays it: `layers` alternating
    /// serpentines over an `n_legs`-leg box of `pitch`, legs `len` long.
    pub(crate) fn woodpile(layers: usize, n_legs: usize, pitch: f64, len: f64, r: f64, fillet: f64, dz: f64) -> Vec<Filament> {
        let edge = 0.3;
        let mut out = Vec::new();
        for k in 0..layers {
            let z = dz * (k as f64 + 1.0);
            let mut verts = Vec::new();
            for i in 0..n_legs {
                let a = edge + pitch * i as f64;
                let (lo, hi) = if i % 2 == 0 { (edge, edge + len) } else { (edge + len, edge) };
                if k % 2 == 0 {
                    verts.push(Point3::new(a, lo, z));
                    verts.push(Point3::new(a, hi, z));
                } else {
                    verts.push(Point3::new(lo, a, z));
                    verts.push(Point3::new(hi, a, z));
                }
            }
            out.push(Filament::serpentine(&verts, r, fillet));
        }
        out
    }

    fn gate(brep: &BRep, shell: ShellId) -> (usize, usize, usize) {
        let mut cfg = crate::config::UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = crate::Diagnostics::new();
        let r = crate::validate::validate_shell(brep, shell, &cfg, &mut diag);
        (r.wire_violations.len(), r.distance_violations.len(), r.manifold_violations.len())
    }

    fn write_if_requested(name: &str, step: &str) {
        if let Ok(dir) = std::env::var("CADCORE_STEP_OUT") {
            std::fs::write(format!("{dir}/{name}"), step).expect("write step file");
        }
    }

    /// Fuse the real pair, chop both terminals, write STEP: the text-level
    /// AP214 gate must be clean, and the file goes to `CADCORE_STEP_OUT` for
    /// the external (OCC / SolidWorks) oracle.
    #[test]
    fn fused_and_cut_pair_exports_valid_step() {
        let fils = woodpile(2, 3, 1.2, 2.4, 0.275, 0.5, 0.35);
        let mut fused = BRep::new();
        let (shell, rep) = fuse_filaments_with(&mut fused, &fils, &FuseOptions::default());
        assert_eq!(gate(&fused, shell).2, 0, "fuse manifold");
        eprintln!("fuse: {rep:?}");
        let merged = crate::arrange::merge::merge_shell_edges(&mut fused, shell);
        eprintln!("merged {merged} edges → {} edges", fused.edges.len());
        assert_eq!(gate(&fused, shell).2, 0, "fuse manifold after merge");
        let step = cadcore_step::brep_to_step(&fused).expect("step");
        let v = crate::validate::step_text::manifold_check(&step);
        assert!(v.is_empty(), "fused STEP text gate: {} e.g. {:?}", v.len(), &v[..v.len().min(5)]);
        write_if_requested("u_fused_pair.step", &step);

        let faces = fused.shells[shell].faces.clone();
        let lo = Plane3::from_origin_normal(Point3::new(0.0, 0.85, 0.0), UnitVec3::Y);
        let first = half_space_cut(&fused, &faces, &lo, 1e-6).expect("cut 1");
        let faces2 = first.brep.shells[first.shell].faces.clone();
        let hi = Plane3::from_origin_normal(Point3::new(0.0, 2.15, 0.0), -UnitVec3::Y);
        let mut second = half_space_cut(&first.brep, &faces2, &hi, 1e-6).expect("cut 2");
        let (w, d, m) = gate(&second.brep, second.shell);
        assert_eq!((w, m), (0, 0), "cut gate wire={w} dist={d} manifold={m}");
        let merged = crate::arrange::merge::merge_shell_edges(&mut second.brep, second.shell);
        eprintln!("cut: merged {merged} edges → {} edges", second.brep.edges.len());
        let (w, d, m) = gate(&second.brep, second.shell);
        assert_eq!((w, m), (0, 0), "cut gate after merge wire={w} dist={d} manifold={m}");
        let step = cadcore_step::brep_to_step(&second.brep).expect("step");
        let v = crate::validate::step_text::manifold_check(&step);
        assert!(v.is_empty(), "{} STEP-text violations e.g. {:?}", v.len(), &v[..v.len().min(5)]);
        write_if_requested("u_fused_cut_pair.step", &step);
    }

    /// A bigger woodpile (4 layers × 5 legs, 5 mm) fused and cut on both Y
    /// sides — the shape of a real DIW Forge export.  Also a timing probe.
    #[test]
    fn four_layer_woodpile_fuse_cut_export() {
        let fils = woodpile(4, 5, 1.2, 4.8, 0.275, 0.5, 0.35);
        let t0 = std::time::Instant::now();
        let mut fused = BRep::new();
        let (shell, rep) = fuse_filaments_with(&mut fused, &fils, &FuseOptions::default());
        eprintln!("fuse 4x5: {rep:?} in {} ms", t0.elapsed().as_millis());
        let (w, d, m) = gate(&fused, shell);
        eprintln!("fuse gate wire={w} dist={d} manifold={m}");
        if m != 0 {
            super::real_generated::dump_open_edges(&fused, shell, 40);
        }
        assert_eq!(m, 0, "fuse manifold");
        let mrep = crate::validate::manifold::check(&fused, shell);
        eprintln!("sense warnings: {} e.g. {:?}", mrep.sense_warnings.len(), &mrep.sense_warnings[..mrep.sense_warnings.len().min(3)]);
        assert!(mrep.sense_warnings.is_empty(), "every shared edge must be traversed in opposite senses");
        eprintln!("edges={} vertices={} faces={}", fused.edges.len(), fused.vertices.len(), fused.faces.len());
        let tm = std::time::Instant::now();
        let merged = crate::arrange::merge::merge_shell_edges(&mut fused, shell);
        eprintln!("merged {merged} edges → {} edges in {} ms", fused.edges.len(), tm.elapsed().as_millis());
        let (w, d, m) = gate(&fused, shell);
        assert_eq!((w, d, m), (0, 0, 0), "fuse gate after merge");
        let step = cadcore_step::brep_to_step(&fused).expect("step");
        let v = crate::validate::step_text::manifold_check(&step);
        assert!(v.is_empty(), "fused STEP text gate: {} e.g. {:?}", v.len(), &v[..v.len().min(5)]);
        write_if_requested("u_fused_4x5.step", &step);

        let faces = fused.shells[shell].faces.clone();
        let t1 = std::time::Instant::now();
        let lo = Plane3::from_origin_normal(Point3::new(0.0, 0.9, 0.0), UnitVec3::Y);
        let first = half_space_cut(&fused, &faces, &lo, 1e-6).expect("cut 1");
        let faces2 = first.brep.shells[first.shell].faces.clone();
        let hi = Plane3::from_origin_normal(Point3::new(0.0, 4.5, 0.0), -UnitVec3::Y);
        let mut second = half_space_cut(&first.brep, &faces2, &hi, 1e-6).expect("cut 2");
        eprintln!("two cuts in {} ms", t1.elapsed().as_millis());
        let (w, d, m) = gate(&second.brep, second.shell);
        eprintln!("cut gate wire={w} dist={d} manifold={m}");
        assert_eq!((w, m), (0, 0), "cut gate wire={w} dist={d} manifold={m}");
        let merged = crate::arrange::merge::merge_shell_edges(&mut second.brep, second.shell);
        eprintln!("cut: merged {merged} edges → {} edges", second.brep.edges.len());
        let (w, d, m) = gate(&second.brep, second.shell);
        assert_eq!((w, m), (0, 0), "cut gate after merge wire={w} dist={d} manifold={m}");
        // diagnostics: polyline edges whose points all lie on one circle are
        // arcs that lost their analytic form somewhere
        {
            let b = &second.brep;
            let mut on_circle = 0usize;
            let mut example: Option<String> = None;
            let mut kinds = std::collections::BTreeMap::<&str, usize>::new();
            for (eid, e) in b.edges.iter() {
                let k = match &e.geom { cadcore_topo::EdgeGeom::Line(_) => "line", cadcore_topo::EdgeGeom::Circle(_) => "circle", cadcore_topo::EdgeGeom::Polyline(_) => "poly", _ => "other" };
                *kinds.entry(k).or_default() += 1;
                if let cadcore_topo::EdgeGeom::Polyline(pts) = &e.geom {
                    if pts.len() >= 4 {
                        let (a, m, c) = (pts[0], pts[pts.len() / 2], pts[pts.len() - 1]);
                        let ab = m - a; let ac = c - a; let n = ab.cross(ac);
                        if n.length() > 1e-12 {
                            let n2 = n.dot(n);
                            let to_c = (n.cross(ab) * ac.dot(ac) + ac.cross(n) * ab.dot(ab)) * (1.0 / (2.0 * n2));
                            let centre = a + to_c; let r = to_c.length();
                            let nn = n * (1.0 / n.length());
                            let ok = pts.iter().all(|&p| { let w = p - centre; (nn.dot(w)).abs() < 1e-6 && ((w - nn * nn.dot(w)).length() - r).abs() < 1e-6 });
                            if ok {
                                on_circle += 1;
                                if example.is_none() {
                                    let faces: Vec<String> = b.faces.iter().filter(|(_, f)| { let mut ls = vec![f.outer_loop]; ls.extend(f.inner_loops.iter().copied()); ls.iter().any(|&l| b.loop_coedges(l).map_or(false, |cs| cs.iter().any(|&c| b.coedges[c].edge == eid))) }).map(|(_, f)| format!("{:?}", std::mem::discriminant(&f.geom))).collect();
                                    example = Some(format!("{} pts r={r:.4} centre={centre:?} faces={faces:?}", pts.len()));
                                }
                            }
                        }
                    }
                }
            }
            eprintln!("edge kinds {kinds:?}; polylines lying on a circle: {on_circle} e.g. {example:?}");
        }
        let step = cadcore_step::brep_to_step(&second.brep).expect("step");
        let v = crate::validate::step_text::manifold_check(&step);
        assert!(v.is_empty(), "{} STEP-text violations", v.len());
        write_if_requested("u_fused_cut_4x5.step", &step);
    }
}

#[cfg(test)]
mod real_generated {
    use super::*;
    use crate::arrange::fixtures::gen;
    use cadcore_math::UnitVec3;

    fn filaments() -> Vec<Filament> {
        gen::pair_3mm()
    }

    pub(crate) fn dump_open_edges(brep: &BRep, shell: ShellId, max: usize) {
        use std::collections::HashMap;
        let mut uses: HashMap<cadcore_topo::EdgeId, Vec<cadcore_topo::FaceId>> = HashMap::new();
        for &fid in &brep.shells[shell].faces {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let Some(cs) = brep.loop_coedges(lid) else { continue };
                for c in cs {
                    uses.entry(brep.coedges[c].edge).or_default().push(fid);
                }
            }
        }
        let mut shown = 0;
        for (eid, fs) in uses.iter().filter(|(_, u)| u.len() != 2) {
            if shown >= max { break; }
            shown += 1;
            let e = &brep.edges[*eid];
            let a = brep.vertices[e.v_start].point;
            let b = brep.vertices[e.v_end].point;
            let kinds: Vec<String> = fs.iter().map(|f| match &brep.faces[*f].geom { FaceGeom::Plane(_) => format!("{f:?}/pl"), FaceGeom::Cylinder(_) => format!("{f:?}/cyl"), FaceGeom::Torus(_) => format!("{f:?}/tor"), _ => format!("{f:?}/?") }).collect();
            let kind = match &e.geom { cadcore_topo::EdgeGeom::Line(_) => "line", cadcore_topo::EdgeGeom::Circle(_) => "circle", cadcore_topo::EdgeGeom::Polyline(_) => "poly", _ => "other" };
            eprintln!("  open {kind} ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}) len={:.2e} used {} by {:?}", a.x, a.y, a.z, b.x, b.y, b.z, (b - a).length(), fs.len(), kinds);
        }
    }

    #[test]
    fn real_generated_pair_fuses() {
        let fils = filaments();
        let mut brep = BRep::new();
        let (shell, rep) = fuse_filaments_with(&mut brep, &fils, &FuseOptions::default());
        eprintln!("{rep:?}");
        let mut cfg = crate::config::UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = crate::Diagnostics::new();
        let r = crate::validate::validate_shell(&brep, shell, &cfg, &mut diag);
        eprintln!("gate wire={} dist={} manifold={}", r.wire_violations.len(), r.distance_violations.len(), r.manifold_violations.len());
        for w in r.wire_violations.iter().take(4) { eprintln!("  wire: {w}"); }
        dump_open_edges(&brep, shell, 30);
        assert_eq!(rep.aborted, 0, "every crossing must trace");
        assert_eq!(r.manifold_violations.len(), 0, "manifold");
        assert_eq!(r.wire_violations.len(), 0, "wires");
    }

    /// The pipeline's own trim: offset 0.8 on both Y sides of the real pair.
    #[test]
    fn real_generated_pair_cut_at_the_pipeline_planes() {
        use crate::boolean::half_space_cut;
        use cadcore_geom::Plane3;
        let fils = filaments();
        let mut brep = BRep::new();
        let (shell, _) = fuse_filaments_with(&mut brep, &fils, &FuseOptions::default());
        let faces = brep.shells[shell].faces.clone();
        let gate = |b: &BRep, sh: ShellId| {
            let mut cfg = crate::config::UnionConfig::default();
            cfg.wires_strict = true;
            let mut diag = crate::Diagnostics::new();
            let r = crate::validate::validate_shell(b, sh, &cfg, &mut diag);
            (r.wire_violations.len(), r.distance_violations.len(), r.manifold_violations.len())
        };
        let lo = Plane3::from_origin_normal(Point3::new(0.0, 0.8, 0.0), UnitVec3::Y);
        let c1 = half_space_cut(&brep, &faces, &lo, 1e-6).expect("cut 1");
        let g1 = gate(&c1.brep, c1.shell);
        eprintln!("cut 1: caps={} dropped={} gate (wire,dist,manifold)={g1:?}", c1.caps, c1.dropped);
        dump_open_edges(&c1.brep, c1.shell, 24);
        assert_eq!(g1.2, 0, "manifold after cut 1");
        let f1 = c1.brep.shells[c1.shell].faces.clone();
        let hi = Plane3::from_origin_normal(Point3::new(0.0, 2.2, 0.0), -UnitVec3::Y);
        let mut c2 = half_space_cut(&c1.brep, &f1, &hi, 1e-6).expect("cut 2");
        let g2 = gate(&c2.brep, c2.shell);
        eprintln!("cut 2: caps={} dropped={} gate={g2:?}", c2.caps, c2.dropped);
        dump_open_edges(&c2.brep, c2.shell, 24);
        assert_eq!(g2.2, 0, "manifold after cut 2");
        let merged = crate::arrange::merge::merge_shell_edges(&mut c2.brep, c2.shell);
        let g3 = gate(&c2.brep, c2.shell);
        eprintln!("merged {merged}: gate={g3:?}");
        assert_eq!((g3.0, g3.2), (0, 0), "gate after merge");
    }
}

#[cfg(test)]
mod real_grid6 {
    use super::*;
    use super::real_generated::dump_open_edges;
    use crate::arrange::fixtures::gen;
    use crate::boolean::half_space_cut;

    fn gate(b: &BRep, sh: ShellId) -> (usize, usize, usize) {
        let mut cfg = crate::config::UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = crate::Diagnostics::new();
        let r = crate::validate::validate_shell(b, sh, &cfg, &mut diag);
        for w in r.wire_violations.iter().take(4) { eprintln!("  wire: {w}"); }
        (r.wire_violations.len(), r.distance_violations.len(), r.manifold_violations.len())
    }

    /// Consecutive loop vertices closer than `tol` — a wire that doubles back
    /// on itself (the validator reports it as a self-intersection).
    pub(crate) fn dump_near_duplicates(brep: &BRep, shell: ShellId, tol: f64) {
        for &fid in &brep.shells[shell].faces {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let Some(cs) = brep.loop_coedges(lid) else { continue };
                let mut prev: Option<(cadcore_topo::VertexId, Point3)> = None;
                for c in cs.iter().chain(cs.iter().take(1)) {
                    let ce = &brep.coedges[*c];
                    let e = &brep.edges[ce.edge];
                    let (vs, ve) = if matches!(ce.sense, cadcore_topo::CoEdgeSense::Same) { (e.v_start, e.v_end) } else { (e.v_end, e.v_start) };
                    let ps = brep.vertices[vs].point;
                    if let Some((pv, pp)) = prev {
                        if pv != vs && (pp - ps).length() < tol {
                            eprintln!(
                                "  near-dup {fid:?} {lid:?}: {pv:?} ({:.7},{:.7},{:.7}) → {vs:?} ({:.7},{:.7},{:.7}) d={:.2e}",
                                pp.x, pp.y, pp.z, ps.x, ps.y, ps.z, (pp - ps).length()
                            );
                        }
                    }
                    let pe = brep.vertices[ve].point;
                    if (pe - ps).length() < tol {
                        eprintln!(
                            "  short edge {fid:?} {lid:?} {:?}: {vs:?} ({:.7},{:.7},{:.7}) → {ve:?} ({:.7},{:.7},{:.7}) d={:.2e}",
                            ce.edge, ps.x, ps.y, ps.z, pe.x, pe.y, pe.z, (pe - ps).length()
                        );
                    }
                    prev = Some((ve, pe));
                }
            }
        }
    }

    /// Signed (θ, φ) area of every loop of every torus face: an outer loop
    /// must be counter-clockwise (positive), an inner one clockwise.
    pub(crate) fn dump_torus_winding(brep: &BRep, shell: ShellId) {
        use std::f64::consts::{PI, TAU};
        for &fid in &brep.shells[shell].faces {
            let f = &brep.faces[fid];
            let FaceGeom::Torus(t) = &f.geom else { continue };
            let uv = |p: Point3| -> (f64, f64) {
                let w = p - t.frame.origin;
                let th = t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w));
                let h = t.frame.z.dot_vec(w);
                let q = w - t.frame.z.as_vec() * h;
                let ph = h.atan2(q.length() - t.major_radius);
                (th, ph)
            };
            let mut ls = vec![(f.outer_loop, true)];
            ls.extend(f.inner_loops.iter().map(|&l| (l, false)));
            for (lid, outer) in ls {
                let Some(cs) = brep.loop_coedges(lid) else { continue };
                let mut pts: Vec<(f64, f64)> = Vec::new();
                for c in cs {
                    let ce = &brep.coedges[c];
                    let mut e = crate::boolean::aabb::sample_edge(brep, ce.edge);
                    if !matches!(ce.sense, cadcore_topo::CoEdgeSense::Same) {
                        e.reverse();
                    }
                    for p in e {
                        let (mut a, mut b) = uv(p);
                        if let Some(&(pa, pb)) = pts.last() {
                            while a - pa > PI { a -= TAU; }
                            while a - pa < -PI { a += TAU; }
                            while b - pb > PI { b -= TAU; }
                            while b - pb < -PI { b += TAU; }
                        }
                        pts.push((a, b));
                    }
                }
                let n = pts.len();
                if n < 3 { continue; }
                let mut area = 0.0;
                for i in 0..n {
                    let (a, b) = (pts[i], pts[(i + 1) % n]);
                    area += a.0 * b.1 - b.0 * a.1;
                }
                area *= 0.5;
                let wrong = (outer && area < 0.0) || (!outer && area > 0.0);
                let closes = (pts[0].0 - pts[n - 1].0).abs() < 1e-6 && (pts[0].1 - pts[n - 1].1).abs() < 1e-6;
                if wrong || !closes {
                    let (tl, th) = pts.iter().fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.0), a.1.max(p.0)));
                    let (pl, ph) = pts.iter().fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.1), a.1.max(p.1)));
                    eprintln!(
                        "  torus {fid:?} {lid:?} outer={outer} uv-area={area:.5} closes={closes} θ[{tl:.3}..{th:.3}] φ[{pl:.3}..{ph:.3}] normal={:?} n={n}",
                        f.normal
                    );
                }
            }
        }
    }

    /// Replay any pipeline filament dump: `CADCORE_FILAMENTS=<log>` and
    /// optionally `CADCORE_CUTS=y:1.0:+,y:5.0:-` (axis : offset : kept side).
    #[test]
    fn replay_filament_dump() {
        let Ok(path) = std::env::var("CADCORE_FILAMENTS") else { return };
        let fils = crate::arrange::fixtures::from_dump(&path);
        eprintln!("{} filaments, {} legs", fils.len(), fils.iter().map(|f| f.legs.len()).sum::<usize>());
        let mut brep = BRep::new();
        let t0 = std::time::Instant::now();
        let (shell, rep) = fuse_filaments_with(&mut brep, &fils, &FuseOptions::default());
        eprintln!("fuse loops={} pieces={} joints={} aborted={} faces={} in {} ms", rep.loops, rep.pieces, rep.joints, rep.aborted, rep.faces, t0.elapsed().as_millis());
        for n in rep.notes.iter().take(12) {
            eprintln!("  note: {n}");
        }
        let g0 = gate(&brep, shell);
        eprintln!("fuse gate {g0:?}");
        dump_open_edges(&brep, shell, 30);
        dump_torus_winding(&brep, shell);
        // Independent volume of (union of tubes) ∩ (kept half-spaces), by
        // stratified sampling of the point-membership of the INPUT tubes — a
        // reference for the volume an external CAD reports for the export.
        {
            let tubes: Vec<CompositeTube> = fils
                .iter()
                .map(|f| {
                    let mut t = CompositeTube::new();
                    for l in &f.legs {
                        t = t.with_leg(l.surf, l.length);
                    }
                    for e in &f.elbows {
                        t = t.with_elbow(e.surf, e.theta_lo, e.theta_hi);
                    }
                    t
                })
                .collect();
            let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
            for f in &fils {
                for l in &f.legs {
                    for p in [l.surf.frame.origin, l.surf.frame.origin + l.surf.axis().as_vec() * l.length] {
                        for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                            lo[k] = lo[k].min(c - l.surf.radius);
                            hi[k] = hi[k].max(c + l.surf.radius);
                        }
                    }
                }
            }
            let mut halves: Vec<(usize, f64, f64)> = Vec::new(); // (axis, offset, sign)
            if let Ok(cuts) = std::env::var("CADCORE_CUTS") {
                for c in cuts.split(',') {
                    let parts: Vec<&str> = c.split(':').collect();
                    let ax = match parts[0] { "x" => 0, "y" => 1, _ => 2 };
                    halves.push((ax, parts[1].parse().unwrap(), if parts[2] == "-" { -1.0 } else { 1.0 }));
                }
            }
            let n: [usize; 3] = match std::env::var("CADCORE_REF_FINE") { Ok(_) => [520, 520, 150], Err(_) => [220, 220, 60] };
            let mut inside = 0usize;
            for i in 0..n[0] {
                for j in 0..n[1] {
                    for k in 0..n[2] {
                        let q = [
                            lo[0] + (hi[0] - lo[0]) * (i as f64 + 0.5) / n[0] as f64,
                            lo[1] + (hi[1] - lo[1]) * (j as f64 + 0.5) / n[1] as f64,
                            lo[2] + (hi[2] - lo[2]) * (k as f64 + 0.5) / n[2] as f64,
                        ];
                        if halves.iter().any(|&(ax, off, sg)| sg * (q[ax] - off) < 0.0) {
                            continue;
                        }
                        let p = Point3::new(q[0], q[1], q[2]);
                        if tubes.iter().any(|t| t.contains(p)) {
                            inside += 1;
                        }
                    }
                }
            }
            let cell = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]) / (n[0] * n[1] * n[2]) as f64;
            eprintln!("reference volume (sampled, after cuts) = {:.3} mm³", inside as f64 * cell);
        }
        let mut cur = (brep, shell);
        if let Ok(cuts) = std::env::var("CADCORE_CUTS") {
            for (k, c) in cuts.split(',').enumerate() {
                let parts: Vec<&str> = c.split(':').collect();
                let off: f64 = parts[1].parse().unwrap();
                let (o, n) = match parts[0] {
                    "x" => (Point3::new(off, 0.0, 0.0), cadcore_math::UnitVec3::X),
                    "y" => (Point3::new(0.0, off, 0.0), cadcore_math::UnitVec3::Y),
                    _ => (Point3::new(0.0, 0.0, off), cadcore_math::UnitVec3::Z),
                };
                let n = if parts[2] == "-" { -n } else { n };
                let faces = cur.0.shells[cur.1].faces.clone();
                let t1 = std::time::Instant::now();
                let rc = crate::boolean::half_space_cut_robust(&cur.0, &faces, &Plane3::from_origin_normal(o, n), 1e-6).expect("cut");
                let shift = rc.shift;
                let c = rc.outcome;
                let g = gate(&c.brep, c.shell);
                eprintln!("cut {k} ({c_spec}): shift={shift:+.3} caps={} dropped={} gate={g:?} in {} ms", c.caps, c.dropped, t1.elapsed().as_millis(), c_spec = cuts.split(',').nth(k).unwrap());
                dump_open_edges(&c.brep, c.shell, 30);
                dump_torus_winding(&c.brep, c.shell);
                cur = (c.brep, c.shell);
            }
        }
    }

    #[test]
    fn grid6_fuses_and_cuts() {
        let fils = gen::grid_6mm_4l();
        let mut brep = BRep::new();
        let t0 = std::time::Instant::now();
        let (shell, rep) = fuse_filaments_with(&mut brep, &fils, &FuseOptions::default());
        eprintln!("fuse {rep:?} in {} ms", t0.elapsed().as_millis());
        let g0 = gate(&brep, shell);
        eprintln!("fuse gate {g0:?}");
        dump_near_duplicates(&brep, shell, 1e-4);
        assert_eq!(rep.aborted, 0);
        assert_eq!(g0.2, 0, "fuse manifold");
        let faces = brep.shells[shell].faces.clone();
        let t1 = std::time::Instant::now();
        let lo = Plane3::from_origin_normal(Point3::new(0.0, 1.0, 0.0), cadcore_math::UnitVec3::Y);
        let c1 = half_space_cut(&brep, &faces, &lo, 1e-6).expect("cut 1");
        let g1 = gate(&c1.brep, c1.shell);
        eprintln!("cut 1: caps={} dropped={} gate={g1:?} in {} ms", c1.caps, c1.dropped, t1.elapsed().as_millis());
        dump_open_edges(&c1.brep, c1.shell, 24);
        assert_eq!(g1.2, 0, "manifold after cut 1");
        let f1 = c1.brep.shells[c1.shell].faces.clone();
        let hi = Plane3::from_origin_normal(Point3::new(0.0, 5.0, 0.0), -cadcore_math::UnitVec3::Y);
        let c2 = half_space_cut(&c1.brep, &f1, &hi, 1e-6).expect("cut 2");
        let g2 = gate(&c2.brep, c2.shell);
        eprintln!("cut 2: caps={} dropped={} gate={g2:?}", c2.caps, c2.dropped);
        dump_open_edges(&c2.brep, c2.shell, 24);
        assert_eq!(g2.2, 0, "manifold after cut 2");
    }
}
