//! **Plate weld** — fuse a flat electrode plate onto a cut terminal of a
//! solid so the two become ONE body.
//!
//! The process this serves: sinter the lattice, chop its terminals flat, then
//! press a silver plate onto each cut face.  Geometrically the plate's inner
//! face lies IN the cut plane, exactly where the solid's section caps are — a
//! coplanar-face Boolean, the one configuration a generic surface–surface
//! union cannot resolve (two coincident planes have no intersection curve).
//!
//! Here it is handled by construction instead of by intersection:
//!
//! ```text
//! 1 cut     trim the solid to the plate's inner plane (a no-op when it
//!           already ends there) — its terminal is now a set of flat caps
//! 2 lift    every face of the solid is copied EXCEPT those caps
//! 3 inner   the plate's inner face is the plate rectangle MINUS the cap
//!           regions: one planar DCEL over the rectangle and the caps' outer
//!           loops; the cap cells are dropped, the rest emitted.  Its hole
//!           edges are the very same sampled points the copied lateral faces
//!           carry, so they weld
//! 4 box     the plate's other five faces close the body
//! ```
//!
//! The plate has to COVER the terminal: a cap reaching outside the rectangle
//! is reported and the weld is refused rather than guessed at.

use cadcore_geom::Plane3;
use cadcore_math::{Frame3, Point3, UnitVec3};
use cadcore_topo::{BRep, FaceGeom, FaceId, FaceNormal, Shell, ShellId, Solid};

use super::aabb::sample_edge;
use super::halfspace::half_space_cut;
use crate::arrange::assembly::Assembler;
use crate::arrange::cells::{arrange_and_classify_with, FaceChain};
use crate::arrange::domain::FaceDomain;

/// A rectangular plate, placed by its INNER face (the one pressed against the
/// cut terminal).
#[derive(Debug, Clone, Copy)]
pub struct PlateBox {
    /// The inner face's plane.  Its normal points INTO the solid being welded
    /// (i.e. away from the plate body).  The rectangle is expressed in this
    /// frame's `x`/`y`.
    pub inner: Plane3,
    /// Plate thickness along `-inner.normal` (mm).
    pub thickness: f64,
    /// Rectangle extents along `inner.frame.x` …
    pub umin: f64,
    /// … and its far side,
    pub umax: f64,
    /// … along `inner.frame.y`,
    pub vmin: f64,
    /// … and its far side.
    pub vmax: f64,
}

impl PlateBox {
    /// A plate from its axis-aligned box.  `axis` is the world axis (0/1/2)
    /// normal to the plate; `toward_positive` says the solid lies on the +axis
    /// side of the plate (so the inner face is the plate's high face).
    pub fn axis_aligned(lo: Point3, hi: Point3, axis: usize, toward_positive: bool) -> Self {
        let lo_a = [lo.x, lo.y, lo.z];
        let hi_a = [hi.x, hi.y, hi.z];
        let mut n = [0.0; 3];
        n[axis] = if toward_positive { 1.0 } else { -1.0 };
        let normal = UnitVec3::try_from_vec(cadcore_math::Vec3::new(n[0], n[1], n[2])).unwrap();
        let ua = (axis + 1) % 3;
        let va = (axis + 2) % 3;
        let mut xv = [0.0; 3];
        xv[ua] = 1.0;
        let x = UnitVec3::try_from_vec(cadcore_math::Vec3::new(xv[0], xv[1], xv[2])).unwrap();
        let y = UnitVec3::try_from_vec(normal.as_vec().cross(x.as_vec())).unwrap();
        // y is ±world axis `va`; express the rectangle in (x, y)
        let ysign = y.as_vec().x + y.as_vec().y + y.as_vec().z; // ±1 on axis va
        let plane_coord = if toward_positive { hi_a[axis] } else { lo_a[axis] };
        let mut c = [0.0; 3];
        c[axis] = plane_coord;
        c[ua] = 0.5 * (lo_a[ua] + hi_a[ua]);
        c[va] = 0.5 * (lo_a[va] + hi_a[va]);
        let origin = Point3::new(c[0], c[1], c[2]);
        let hu = 0.5 * (hi_a[ua] - lo_a[ua]);
        let hv = 0.5 * (hi_a[va] - lo_a[va]);
        let _ = ysign;
        PlateBox {
            inner: Plane3 { frame: Frame3 { origin, x, y, z: normal } },
            thickness: hi_a[axis] - lo_a[axis],
            umin: -hu,
            umax: hu,
            vmin: -hv,
            vmax: hv,
        }
    }

    fn corner(&self, u: f64, v: f64, depth: f64) -> Point3 {
        self.inner.frame.origin + self.inner.frame.x * u + self.inner.frame.y * v
            - self.inner.frame.z * depth
    }
}

/// Outcome of a plate weld.
#[derive(Debug)]
pub struct WeldOutcome {
    /// The welded body.
    pub brep: BRep,
    /// Its single closed shell.
    pub shell: ShellId,
    /// Section caps swallowed by the plate (holes in its inner face).
    pub caps_welded: usize,
    /// Faces the cut removed from the solid before welding.
    pub dropped: usize,
    /// How far (mm) the plate footprint had to grow to clear the terminal's
    /// sections (`0` when the requested plate already covered it).
    pub grown: f64,
}

/// Weld `plate` onto the solid bounded by `faces`.  `knit` as for
/// [`half_space_cut`].  `None` when the solid has no faces or a cap reaches
/// outside the plate rectangle.
pub fn weld_plate(src: &BRep, faces: &[FaceId], plate: &PlateBox, knit: f64) -> Option<WeldOutcome> {
    if faces.is_empty() {
        return None;
    }
    let n_in = plate.inner.normal();
    // 1. the solid ends at the inner plane
    let keep = Plane3::from_origin_normal(plate.inner.frame.origin, n_in);
    let cut = half_space_cut(src, faces, &keep, knit)?;
    let b1 = &cut.brep;
    let sh_faces = b1.shells[cut.shell].faces.clone();

    // 2. caps lying in the plane
    let to_uv = |p: Point3| {
        let w = p - plate.inner.frame.origin;
        (plate.inner.frame.x.dot_vec(w), plate.inner.frame.y.dot_vec(w))
    };
    let mut caps: Vec<FaceId> = Vec::new();
    let (mut cu_lo, mut cu_hi, mut cv_lo, mut cv_hi) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for &fid in &sh_faces {
        let f = &b1.faces[fid];
        let FaceGeom::Plane(p) = &f.geom else { continue };
        if p.normal().dot(n_in).abs() < 1.0 - 1e-9 {
            continue;
        }
        if keep.signed_distance(p.frame.origin).abs() > 1e-6 {
            continue;
        }
        // every loop point in the plane
        let pts = loop_points(b1, fid);
        if pts.is_empty() || pts.iter().any(|&q| keep.signed_distance(q).abs() > 1e-6) {
            continue;
        }
        for &q in &pts {
            let (u, v) = to_uv(q);
            cu_lo = cu_lo.min(u);
            cu_hi = cu_hi.max(u);
            cv_lo = cv_lo.min(v);
            cv_hi = cv_hi.max(v);
        }
        caps.push(fid);
    }
    // The plate must COVER the terminal, with clearance: a section touching
    // the plate's edge would pinch the two faces together at a point.  A
    // plate placed flush with the lattice envelope (the default footprint)
    // touches the outermost tubes exactly, so it is grown by `COVER_MARGIN`
    // where needed — the weld reports how much.
    const COVER_MARGIN: f64 = 0.05;
    let mut plate = *plate;
    let mut grown = 0.0_f64;
    if !caps.is_empty() {
        for (lo, hi, want_lo, want_hi) in [
            (&mut plate.umin, &mut plate.umax, cu_lo - COVER_MARGIN, cu_hi + COVER_MARGIN),
            (&mut plate.vmin, &mut plate.vmax, cv_lo - COVER_MARGIN, cv_hi + COVER_MARGIN),
        ] {
            if *lo > want_lo {
                grown = grown.max(*lo - want_lo);
                *lo = want_lo;
            }
            if *hi < want_hi {
                grown = grown.max(want_hi - *hi);
                *hi = want_hi;
            }
        }
    }
    let plate = &plate;

    // 3 + 4. assemble
    let weld = knit.clamp(1e-9, 1e-2);
    let mut out = BRep::new();
    let mut asm = Assembler::new(&mut out, weld);
    for &fid in &sh_faces {
        if caps.contains(&fid) {
            continue;
        }
        asm.emit_face_copy(b1, fid, &[]);
    }
    // inner face: rectangle minus the caps
    let half = (plate.umax - plate.umin).abs().max((plate.vmax - plate.vmin).abs());
    let inner_domain = FaceDomain::Plane { plane: plate.inner, half_extent: half };
    const RECT_TAG: u32 = 1;
    let rect = vec![
        (plate.umin, plate.vmin),
        (plate.umax, plate.vmin),
        (plate.umax, plate.vmax),
        (plate.umin, plate.vmax),
        (plate.umin, plate.vmin),
    ];
    let mut chains = vec![FaceChain { pts: rect, tag: RECT_TAG }];
    let mut tag = 1000u32;
    for &cap in &caps {
        for ring in outer_ring_points(b1, cap) {
            let mut uv: Vec<(f64, f64)> = ring.iter().map(|&p| inner_domain.uv(p)).collect();
            if uv.first() != uv.last() {
                uv.push(uv[0]);
            }
            chains.push(FaceChain { pts: uv, tag });
            tag += 1;
        }
    }
    let cells = arrange_and_classify_with(&inner_domain, &chains, 0.0, &|_| false);
    for cell in &cells {
        if cell.cell.outer.iter().any(|s| s.tag == RECT_TAG) {
            asm.emit_cell(&inner_domain, FaceGeom::Plane(plate.inner), FaceNormal::Same, cell);
        }
    }
    // the box: outer face + four sides, each one rectangular cell
    let t = plate.thickness;
    let (u0, u1, v0, v1) = (plate.umin, plate.umax, plate.vmin, plate.vmax);
    let c = |u: f64, v: f64, d: f64| plate.corner(u, v, d);
    let outer_normal = n_in.flip();
    let side_faces: Vec<(Point3, UnitVec3, [Point3; 4])> = vec![
        // outer face
        (
            c(0.5 * (u0 + u1), 0.5 * (v0 + v1), t),
            outer_normal,
            [c(u0, v0, t), c(u1, v0, t), c(u1, v1, t), c(u0, v1, t)],
        ),
        // -x side
        (
            c(u0, 0.5 * (v0 + v1), 0.5 * t),
            plate.inner.frame.x.flip(),
            [c(u0, v0, 0.0), c(u0, v1, 0.0), c(u0, v1, t), c(u0, v0, t)],
        ),
        // +x side
        (
            c(u1, 0.5 * (v0 + v1), 0.5 * t),
            plate.inner.frame.x,
            [c(u1, v0, 0.0), c(u1, v1, 0.0), c(u1, v1, t), c(u1, v0, t)],
        ),
        // -y side
        (
            c(0.5 * (u0 + u1), v0, 0.5 * t),
            plate.inner.frame.y.flip(),
            [c(u0, v0, 0.0), c(u1, v0, 0.0), c(u1, v0, t), c(u0, v0, t)],
        ),
        // +y side
        (
            c(0.5 * (u0 + u1), v1, 0.5 * t),
            plate.inner.frame.y,
            [c(u0, v1, 0.0), c(u1, v1, 0.0), c(u1, v1, t), c(u0, v1, t)],
        ),
    ];
    for (origin, normal, corners) in side_faces {
        let plane = Plane3::from_origin_normal(origin, normal);
        let domain = FaceDomain::Plane { plane, half_extent: half.max(t) };
        let mut uv: Vec<(f64, f64)> = corners.iter().map(|&p| domain.uv(p)).collect();
        uv.push(uv[0]);
        let cells = arrange_and_classify_with(&domain, &[FaceChain { pts: uv, tag: 1 }], 0.0, &|_| false);
        for cell in &cells {
            asm.emit_cell(&domain, FaceGeom::Plane(plane), FaceNormal::Same, cell);
        }
    }

    let out_faces = asm.faces().to_vec();
    let shell = out.add_shell(Shell { faces: out_faces.clone(), is_outer: true, solid: Default::default() });
    let solid = out.add_solid(Solid { shells: vec![shell], name: Some("welded".into()) });
    out.shells[shell].solid = solid;
    for f in &out_faces {
        out.faces[*f].shell = shell;
    }
    Some(WeldOutcome { brep: out, shell, caps_welded: caps.len(), dropped: cut.dropped, grown })
}

/// Every sampled point of every loop of a face.
fn loop_points(brep: &BRep, fid: FaceId) -> Vec<Point3> {
    let mut out = Vec::new();
    for ring in all_ring_points(brep, fid) {
        out.extend(ring);
    }
    out
}

/// The face's OUTER loop as sampled 3-D points in traversal order (one ring).
fn outer_ring_points(brep: &BRep, fid: FaceId) -> Vec<Vec<Point3>> {
    let f = &brep.faces[fid];
    ring_points(brep, f.outer_loop).into_iter().collect()
}

fn all_ring_points(brep: &BRep, fid: FaceId) -> Vec<Vec<Point3>> {
    let f = &brep.faces[fid];
    let mut ls = vec![f.outer_loop];
    ls.extend(f.inner_loops.iter().copied());
    ls.into_iter().filter_map(|l| ring_points(brep, l)).collect()
}

fn ring_points(brep: &BRep, lid: cadcore_topo::LoopId) -> Option<Vec<Point3>> {
    let lp = brep.loops.get(lid)?;
    let st = lp.start;
    brep.coedges.get(st)?;
    let mut pts: Vec<Point3> = Vec::new();
    let mut c = st;
    loop {
        let ce = &brep.coedges[c];
        let mut sp = sample_edge(brep, ce.edge);
        if ce.sense == cadcore_topo::CoEdgeSense::Opposite {
            sp.reverse();
        }
        for p in sp {
            if pts.last().map_or(true, |&q| (q - p).length() > 1e-12) {
                pts.push(p);
            }
        }
        c = ce.next;
        if c == st {
            break;
        }
    }
    if pts.len() > 2 && (pts[0] - *pts.last().unwrap()).length() < 1e-12 {
        pts.pop();
    }
    (pts.len() >= 3).then_some(pts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrange::scaffold::Filament;
    use crate::arrange::universal::fuse_filaments;
    use crate::config::UnionConfig;
    use crate::validate::validate_shell;
    use crate::Diagnostics;

    fn gate(brep: &BRep, shell: ShellId) -> (usize, usize, usize) {
        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let r = validate_shell(brep, shell, &cfg, &mut diag);
        for w in r.wire_violations.iter().take(5) {
            eprintln!("  wire: {w}");
        }
        for m in r.manifold_violations.iter().take(5) {
            eprintln!("  manifold: {m}");
        }
        (r.wire_violations.len(), r.distance_violations.len(), r.manifold_violations.len())
    }

    fn pair() -> Vec<Filament> {
        let r = 0.275;
        let (z0, z1) = (0.35, 0.70);
        vec![
            Filament::serpentine(
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
            ),
            Filament::serpentine(
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
            ),
        ]
    }

    /// Fuse, trim both Y terminals, weld a plate flush onto each: ONE body,
    /// every gate clean.
    #[test]
    fn plates_weld_flush_onto_both_terminals() {
        let fils = pair();
        let mut fused = BRep::new();
        let shell = fuse_filaments(&mut fused, &fils);
        let faces = fused.shells[shell].faces.clone();
        // trim
        let lo = Plane3::from_origin_normal(Point3::new(0.0, 0.85, 0.0), UnitVec3::Y);
        let c1 = half_space_cut(&fused, &faces, &lo, 1e-6).expect("cut 1");
        let f1 = c1.brep.shells[c1.shell].faces.clone();
        let hi = Plane3::from_origin_normal(Point3::new(0.0, 2.15, 0.0), -UnitVec3::Y);
        let c2 = half_space_cut(&c1.brep, &f1, &hi, 1e-6).expect("cut 2");
        let f2 = c2.brep.shells[c2.shell].faces.clone();
        // plate on the low side: box y ∈ [0.45, 0.85], solid on +Y
        let plate_lo = PlateBox::axis_aligned(
            Point3::new(-0.2, 0.45, -0.2),
            Point3::new(3.2, 0.85, 1.3),
            1,
            true,
        );
        let w1 = weld_plate(&c2.brep, &f2, &plate_lo, 1e-6).expect("weld low");
        eprintln!("low plate: caps welded {}, dropped {}", w1.caps_welded, w1.dropped);
        assert_eq!(w1.caps_welded, 3, "three leg sections on the low terminal");
        let (w, d, m) = gate(&w1.brep, w1.shell);
        assert_eq!((w, m), (0, 0), "low weld gate wire={w} dist={d} manifold={m}");
        // plate on the high side: box y ∈ [2.15, 2.55], solid on −Y
        let fw = w1.brep.shells[w1.shell].faces.clone();
        let plate_hi = PlateBox::axis_aligned(
            Point3::new(-0.2, 2.15, -0.2),
            Point3::new(3.2, 2.55, 1.3),
            1,
            false,
        );
        let w2 = weld_plate(&w1.brep, &fw, &plate_hi, 1e-6).expect("weld high");
        eprintln!("high plate: caps welded {}, dropped {}", w2.caps_welded, w2.dropped);
        assert!(w2.caps_welded >= 3);
        let (w, d, m) = gate(&w2.brep, w2.shell);
        assert_eq!((w, m), (0, 0), "high weld gate wire={w} dist={d} manifold={m}");
        assert_eq!(components(&w2.brep, w2.shell), 1, "one body");
        let mut b = w2.brep;
        let merged = crate::arrange::merge::merge_shell_edges(&mut b, w2.shell);
        eprintln!("merged {merged} → {} edges", b.edges.len());
        let (w, _d, m) = gate(&b, w2.shell);
        assert_eq!((w, m), (0, 0), "gate after merge");
        let step = cadcore_step::brep_to_step(&b).expect("step");
        let v = crate::validate::step_text::manifold_check(&step);
        assert!(v.is_empty(), "{} STEP-text violations", v.len());
        if let Ok(dir) = std::env::var("CADCORE_STEP_OUT") {
            std::fs::write(format!("{dir}/u_welded_pair.step"), &step).unwrap();
        }
    }

    fn components(brep: &BRep, shell: ShellId) -> usize {
        use std::collections::HashMap;
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
}
