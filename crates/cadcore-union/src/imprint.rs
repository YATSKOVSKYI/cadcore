//! Contact imprint: partition a planar face of one solid by the faces of
//! another solid that lie on it.
//!
//! Two solids that only *touch* — an electrode plate pressed flat against the
//! cut terminal of a lattice — share no topology.  A FEM mesher meshes each on
//! its own, the nodes on the two sides of the contact do not coincide, and no
//! current crosses it.  CAD systems answer this with an *imprint*: the outline
//! of every contact patch is drawn onto the other body's face, which is split
//! there.  Both bodies then present the same faces at the contact, and the
//! mesher can make them conformal.  Nothing moves; no volume is added or
//! removed.
//!
//! [`imprint_coplanar`] does exactly that for planar contacts: every *tool*
//! face that lies in the plane of the *target* face (and strictly inside it)
//! becomes
//!
//! * a **patch** of the target solid — a new face on the target's surface,
//!   bounded by copies of the tool face's own edges (same curves, same
//!   parameters, so the two sides are identical, bit for bit), and
//! * an **inner loop** of the target face, sharing those edges with the patch.
//!
//! A hole of a tool face (a region of the plane enclosed by a contact patch but
//! not in contact) becomes a face of its own.  The target solid stays a closed
//! 2-manifold: every new edge is used by exactly two of its faces, in opposite
//! senses.
//!
//! The target face must carry real B-Rep topology; [`build_box_solid`] builds
//! an axis-aligned plate that does (the template boxes of `cadcore-ops` are
//! written from closed-form extents and cannot be split).

use cadcore_geom::{Line3, Plane3};
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{
    BRep, CoEdge, CoEdgeId, CoEdgeSense, Edge, EdgeGeom, EdgeId, Face, FaceExtent, FaceGeom, FaceId, FaceNormal, Loop,
    LoopId, Shell, ShellId, Solid, SolidId, Vertex, VertexId,
};
use std::collections::HashMap;

/// What an imprint did.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ImprintReport {
    /// Contact patches added to the target solid (one per tool face).
    pub patches: usize,
    /// Faces added for holes inside contact patches.
    pub holes: usize,
    /// Area of the contact patches (holes excluded), from the sampled outlines.
    pub area: f64,
}

/// Why an imprint was refused.  Nothing is changed in that case.
#[derive(Debug, Clone, PartialEq)]
pub enum ImprintError {
    /// The target face is not planar.
    TargetNotPlanar,
    /// The target face has no real loops (a closed-form template face).
    TargetHasNoTopology,
    /// A tool face is not planar or does not lie in the target's plane.
    NotCoplanar(FaceId),
    /// A tool face has no real loops.
    ToolHasNoTopology(FaceId),
    /// A tool face reaches the boundary of the target face, or lies outside it:
    /// the patch would touch or cross the target's outline.
    NotInside {
        /// The offending tool face.
        face: FaceId,
        /// Smallest clearance to the target outline (negative = outside).
        clearance: f64,
    },
}

impl std::fmt::Display for ImprintError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ImprintError::TargetNotPlanar => write!(f, "the target face is not planar"),
            ImprintError::TargetHasNoTopology => write!(f, "the target face has no B-Rep loops (template face)"),
            ImprintError::NotCoplanar(id) => write!(f, "tool face {id:?} is not in the target's plane"),
            ImprintError::ToolHasNoTopology(id) => write!(f, "tool face {id:?} has no B-Rep loops"),
            ImprintError::NotInside { face, clearance } => write!(
                f,
                "tool face {face:?} is not strictly inside the target face (clearance {clearance:.2e} mm)"
            ),
        }
    }
}

impl std::error::Error for ImprintError {}

// ── Axis-aligned plate with real topology ─────────────────────────────────────

/// An axis-aligned box `lo .. hi` as a closed solid with explicit B-Rep
/// topology: 8 vertices, 12 line edges, 6 planar faces, each edge used by two
/// faces in opposite senses, every outer loop counter-clockwise about its
/// outward normal.  Unlike a template box its faces can be imprinted.
pub fn build_box_solid(brep: &mut BRep, lo: [f64; 3], hi: [f64; 3], name: Option<String>) -> SolidId {
    let corner = |i: usize| -> Point3 {
        Point3::new(
            if i & 1 == 0 { lo[0] } else { hi[0] },
            if i & 2 == 0 { lo[1] } else { hi[1] },
            if i & 4 == 0 { lo[2] } else { hi[2] },
        )
    };
    let verts: Vec<VertexId> = (0..8).map(|i| brep.add_vertex(Vertex { point: corner(i) })).collect();
    // faces as corner cycles, each counter-clockwise about its outward normal
    let faces: [([usize; 4], UnitVec3); 6] = [
        ([0, 4, 6, 2], -UnitVec3::X),
        ([1, 3, 7, 5], UnitVec3::X),
        ([0, 1, 5, 4], -UnitVec3::Y),
        ([2, 6, 7, 3], UnitVec3::Y),
        ([0, 2, 3, 1], -UnitVec3::Z),
        ([4, 5, 7, 6], UnitVec3::Z),
    ];
    let mut edges: HashMap<(usize, usize), EdgeId> = HashMap::new();
    let mut face_ids = Vec::with_capacity(6);
    for (cyc, normal) in faces {
        let face = brep.add_face(Face {
            geom: FaceGeom::Plane(Plane3::from_origin_normal(corner(cyc[0]), normal)),
            normal: FaceNormal::Same,
            outer_loop: LoopId::default(),
            inner_loops: Vec::new(),
            shell: ShellId::default(),
            extent: FaceExtent::Trimmed,
        });
        let lp = brep.add_loop(Loop { start: CoEdgeId::default(), face });
        let mut ces = Vec::with_capacity(4);
        for k in 0..4 {
            let (a, b) = (cyc[k], cyc[(k + 1) % 4]);
            let key = (a.min(b), a.max(b));
            let (edge, sense) = match edges.get(&key) {
                Some(&e) => (e, if brep.edges[e].v_start == verts[a] { CoEdgeSense::Same } else { CoEdgeSense::Opposite }),
                None => {
                    let (pa, pb) = (corner(a), corner(b));
                    let dir = UnitVec3::try_from_vec(pb - pa).expect("box with zero extent");
                    let e = brep.add_edge(Edge {
                        geom: EdgeGeom::Line(Line3::new(pa, dir)),
                        v_start: verts[a],
                        v_end: verts[b],
                        t_start: 0.0,
                        t_end: (pb - pa).length(),
                        partner: None,
                    });
                    edges.insert(key, e);
                    (e, CoEdgeSense::Same)
                }
            };
            ces.push(brep.add_coedge(CoEdge {
                edge,
                sense,
                next: CoEdgeId::default(),
                prev: CoEdgeId::default(),
                loop_id: lp,
            }));
        }
        brep.patch_coedge_links(&ces);
        brep.loops[lp].start = ces[0];
        brep.faces[face].outer_loop = lp;
        face_ids.push(face);
    }
    let shell = brep.add_shell(Shell { faces: face_ids.clone(), is_outer: true, solid: SolidId::default() });
    let solid = brep.add_solid(Solid { shells: vec![shell], name });
    brep.shells[shell].solid = solid;
    for f in face_ids {
        brep.faces[f].shell = shell;
    }
    solid
}

// ── Imprint ───────────────────────────────────────────────────────────────────

/// Faces of `shell` that lie on `plane` (every sampled boundary point within
/// `tol`) and face the given direction: planar, outward normal within 1e-9 of
/// `facing`.  Used to find what a plate will touch before the plate exists.
pub fn faces_on_plane(brep: &BRep, shell: ShellId, plane: &Plane3, facing: UnitVec3, tol: f64) -> Vec<FaceId> {
    let Some(sh) = brep.shells.get(shell) else { return Vec::new() };
    sh.faces
        .iter()
        .copied()
        .filter(|&f| {
            let Some((_, n)) = outward_plane(brep, f) else { return false };
            n.dot_vec(facing.as_vec()) > 1.0 - 1e-9
                && face_points(brep, f).is_some_and(|pts| !pts.is_empty() && pts.iter().all(|&q| plane.signed_distance(q).abs() <= tol))
        })
        .collect()
}

/// Bounding box of the sampled boundaries of `faces`.
pub fn faces_bounds(brep: &BRep, faces: &[FaceId]) -> Option<([f64; 3], [f64; 3])> {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for &f in faces {
        for p in face_points(brep, f)? {
            for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                lo[k] = lo[k].min(c);
                hi[k] = hi[k].max(c);
            }
        }
    }
    lo[0].is_finite().then_some((lo, hi))
}

/// Faces of `shell` that lie in the plane of `target` (every vertex within
/// `tol`) and face it: planar, outward normal opposite to the target's.
pub fn coplanar_contact_faces(brep: &BRep, shell: ShellId, target: FaceId, tol: f64) -> Vec<FaceId> {
    let Some((plane, n_out)) = outward_plane(brep, target) else { return Vec::new() };
    let Some(sh) = brep.shells.get(shell) else { return Vec::new() };
    sh.faces
        .iter()
        .copied()
        .filter(|&f| {
            let Some((p, n)) = outward_plane(brep, f) else { return false };
            let _ = p;
            n.dot_vec(n_out.as_vec()) < -1.0 + 1e-9
                && face_points(brep, f).is_some_and(|pts| !pts.is_empty() && pts.iter().all(|&q| plane.signed_distance(q).abs() <= tol))
        })
        .collect()
}

/// Imprint the tool faces onto the target face (see the module docs).
///
/// `tol` is the coplanarity tolerance; a tool face must also keep at least
/// `tol` of clearance from the target's outline.  On error nothing changes.
pub fn imprint_coplanar(brep: &mut BRep, target: FaceId, tools: &[FaceId], tol: f64) -> Result<ImprintReport, ImprintError> {
    let (plane, n_out) = outward_plane(brep, target).ok_or(ImprintError::TargetNotPlanar)?;
    if !has_loops(brep, brep.faces[target].outer_loop) {
        return Err(ImprintError::TargetHasNoTopology);
    }
    // the target outline in the plane's (u, v), oriented counter-clockwise about the outward normal
    let frame = plane.frame;
    let uv = |p: Point3| {
        let w = p - frame.origin;
        (frame.x.dot_vec(w), frame.y.dot_vec(w))
    };
    let flip = n_out.dot_vec(frame.z.as_vec()) < 0.0; // outward = −z of the plane frame
    let area_about_out = |poly: &[(f64, f64)]| if flip { -shoelace(poly) } else { shoelace(poly) };
    let outline: Vec<(f64, f64)> = loop_points(brep, brep.faces[target].outer_loop).iter().map(|&p| uv(p)).collect();
    let holes: Vec<Vec<(f64, f64)>> = brep.faces[target]
        .inner_loops
        .iter()
        .map(|&l| loop_points(brep, l).iter().map(|&p| uv(p)).collect())
        .collect();

    // ── validate every tool face before touching anything ─────────────────────
    for &t in tools {
        let Some((_, n)) = outward_plane(brep, t) else { return Err(ImprintError::NotCoplanar(t)) };
        if n.dot_vec(n_out.as_vec()).abs() < 1.0 - 1e-9 {
            return Err(ImprintError::NotCoplanar(t));
        }
        if !has_loops(brep, brep.faces[t].outer_loop) {
            return Err(ImprintError::ToolHasNoTopology(t));
        }
        let pts = face_points(brep, t).unwrap_or_default();
        if pts.iter().any(|&q| plane.signed_distance(q).abs() > tol) {
            return Err(ImprintError::NotCoplanar(t));
        }
        let mut clearance = f64::INFINITY;
        for &q in &pts {
            let p = uv(q);
            let mut c = if point_in_polygon(&outline, p) { dist_to_polygon(&outline, p) } else { -dist_to_polygon(&outline, p) };
            for h in &holes {
                c = c.min(if point_in_polygon(h, p) { -dist_to_polygon(h, p) } else { dist_to_polygon(h, p) });
            }
            clearance = clearance.min(c);
        }
        if !(clearance > tol) {
            return Err(ImprintError::NotInside { face: t, clearance });
        }
    }

    // ── build ────────────────────────────────────────────────────────────────
    let shell = brep.faces[target].shell;
    let (geom, normal) = (brep.faces[target].geom.clone(), brep.faces[target].normal);
    let mut vmap: HashMap<VertexId, VertexId> = HashMap::new();
    let mut emap: HashMap<EdgeId, EdgeId> = HashMap::new();
    let mut report = ImprintReport::default();
    let mut new_faces = Vec::new();
    for &t in tools {
        let tool_loops: Vec<LoopId> =
            std::iter::once(brep.faces[t].outer_loop).chain(brep.faces[t].inner_loops.iter().copied()).collect();
        let patch = brep.add_face(Face {
            geom: geom.clone(),
            normal,
            outer_loop: LoopId::default(),
            inner_loops: Vec::new(),
            shell,
            extent: FaceExtent::Trimmed,
        });
        new_faces.push(patch);
        for (k, &tl) in tool_loops.iter().enumerate() {
            let pts: Vec<(f64, f64)> = loop_points(brep, tl).iter().map(|&p| uv(p)).collect();
            let ccw = area_about_out(&pts) > 0.0;
            let a = shoelace(&pts).abs();
            let steps = loop_steps(brep, tl);
            if k == 0 {
                // the patch's outline (CCW) and the matching hole in the target (CW)
                let patch_loop = copy_loop(brep, &steps, !ccw, patch, &mut vmap, &mut emap);
                brep.faces[patch].outer_loop = patch_loop;
                let hole = copy_loop(brep, &steps, ccw, target, &mut vmap, &mut emap);
                brep.faces[target].inner_loops.push(hole);
                report.area += a;
            } else {
                // a hole of the contact: an inner loop of the patch (CW) and a face of its own (CCW)
                let inner = copy_loop(brep, &steps, ccw, patch, &mut vmap, &mut emap);
                brep.faces[patch].inner_loops.push(inner);
                let face = brep.add_face(Face {
                    geom: geom.clone(),
                    normal,
                    outer_loop: LoopId::default(),
                    inner_loops: Vec::new(),
                    shell,
                    extent: FaceExtent::Trimmed,
                });
                let outer = copy_loop(brep, &steps, !ccw, face, &mut vmap, &mut emap);
                brep.faces[face].outer_loop = outer;
                new_faces.push(face);
                report.holes += 1;
                report.area -= a;
            }
        }
        report.patches += 1;
    }
    if let Some(sh) = brep.shells.get_mut(shell) {
        sh.faces.extend(new_faces);
    }
    Ok(report)
}

// ── helpers ──────────────────────────────────────────────────────────────────

/// The face's carrier plane and its OUTWARD normal (the face normal flag applied).
fn outward_plane(brep: &BRep, f: FaceId) -> Option<(Plane3, UnitVec3)> {
    let face = brep.faces.get(f)?;
    match &face.geom {
        FaceGeom::Plane(p) => Some((*p, match face.normal {
            FaceNormal::Same => p.normal(),
            FaceNormal::Reversed => -p.normal(),
        })),
        _ => None,
    }
}

fn has_loops(brep: &BRep, l: LoopId) -> bool {
    brep.loop_coedges(l).is_some_and(|c| !c.is_empty() && brep.coedges.get(c[0]).is_some())
}

/// The loop's co-edges as (edge, sense) in traversal order.
fn loop_steps(brep: &BRep, l: LoopId) -> Vec<(EdgeId, CoEdgeSense)> {
    brep.loop_coedges(l)
        .unwrap_or_default()
        .into_iter()
        .filter_map(|c| brep.coedges.get(c).map(|ce| (ce.edge, ce.sense)))
        .collect()
}

/// Sampled points of a loop in traversal order (closing point not repeated).
fn loop_points(brep: &BRep, l: LoopId) -> Vec<Point3> {
    let mut out: Vec<Point3> = Vec::new();
    for (e, s) in loop_steps(brep, l) {
        let mut pts = crate::boolean::aabb::sample_edge(brep, e);
        if s == CoEdgeSense::Opposite {
            pts.reverse();
        }
        for p in pts {
            if out.last().map_or(true, |&q| (q - p).length() > 1e-12) {
                out.push(p);
            }
        }
    }
    if out.len() > 1 && (out[0] - *out.last().unwrap()).length() < 1e-12 {
        out.pop();
    }
    out
}

fn face_points(brep: &BRep, f: FaceId) -> Option<Vec<Point3>> {
    let face = brep.faces.get(f)?;
    let mut out = loop_points(brep, face.outer_loop);
    for &l in &face.inner_loops {
        out.extend(loop_points(brep, l));
    }
    Some(out)
}

/// Copy a loop onto `face`, reversing it when asked, with new vertices and
/// edges (shared through the maps, so the two copies of one tool edge — the
/// patch's and the target hole's — are one edge of the target solid).
fn copy_loop(
    brep: &mut BRep,
    steps: &[(EdgeId, CoEdgeSense)],
    reverse: bool,
    face: FaceId,
    vmap: &mut HashMap<VertexId, VertexId>,
    emap: &mut HashMap<EdgeId, EdgeId>,
) -> LoopId {
    let lp = brep.add_loop(Loop { start: CoEdgeId::default(), face });
    let order: Vec<(EdgeId, CoEdgeSense)> = if reverse {
        steps
            .iter()
            .rev()
            .map(|&(e, s)| (e, if s == CoEdgeSense::Same { CoEdgeSense::Opposite } else { CoEdgeSense::Same }))
            .collect()
    } else {
        steps.to_vec()
    };
    let mut ces = Vec::with_capacity(order.len());
    for (e, s) in order {
        let ne = match emap.get(&e) {
            Some(&ne) => ne,
            None => {
                let src = brep.edges[e].clone();
                let mut vert = |v: VertexId, brep: &mut BRep| -> VertexId {
                    *vmap.entry(v).or_insert_with(|| brep.add_vertex(Vertex { point: brep.vertices[v].point }))
                };
                let (vs, ve) = (vert(src.v_start, brep), vert(src.v_end, brep));
                let ne = brep.add_edge(Edge { geom: src.geom, v_start: vs, v_end: ve, t_start: src.t_start, t_end: src.t_end, partner: None });
                emap.insert(e, ne);
                // the copy is written with the source's own curve (see `BRep::edge_twins`)
                let src_root = brep.edge_twins.get(&e).copied().unwrap_or(e);
                brep.edge_twins.insert(ne, src_root);
                ne
            }
        };
        ces.push(brep.add_coedge(CoEdge { edge: ne, sense: s, next: CoEdgeId::default(), prev: CoEdgeId::default(), loop_id: lp }));
    }
    brep.patch_coedge_links(&ces);
    if let Some(&first) = ces.first() {
        brep.loops[lp].start = first;
    }
    lp
}

fn shoelace(p: &[(f64, f64)]) -> f64 {
    let n = p.len();
    (0..n).map(|i| p[i].0 * p[(i + 1) % n].1 - p[(i + 1) % n].0 * p[i].1).sum::<f64>() * 0.5
}

fn point_in_polygon(poly: &[(f64, f64)], p: (f64, f64)) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n.wrapping_sub(1);
    for i in 0..n {
        let (a, b) = (poly[i], poly[j]);
        if (a.1 > p.1) != (b.1 > p.1) && p.0 < (b.0 - a.0) * (p.1 - a.1) / (b.1 - a.1) + a.0 {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn dist_to_polygon(poly: &[(f64, f64)], p: (f64, f64)) -> f64 {
    let n = poly.len();
    (0..n)
        .map(|i| {
            let (a, b) = (poly[i], poly[(i + 1) % n]);
            let (dx, dy) = (b.0 - a.0, b.1 - a.1);
            let l2 = dx * dx + dy * dy;
            let t = if l2 > 0.0 { (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0) } else { 0.0 };
            ((p.0 - a.0 - t * dx).powi(2) + (p.1 - a.1 - t * dy).powi(2)).sqrt()
        })
        .fold(f64::INFINITY, f64::min)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrange::fixtures::gen;
    use crate::{fuse_and_trim, validate_shell, Diagnostics, FuseOptions, UnionConfig};

    fn gates(brep: &BRep, shell: ShellId) -> (usize, usize, usize) {
        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut d = Diagnostics::new();
        let r = validate_shell(brep, shell, &cfg, &mut d);
        for v in r.manifold_violations.iter().chain(&r.wire_violations).chain(&r.distance_violations).take(6) {
            eprintln!("  {v}");
        }
        (r.wire_violations.len(), r.distance_violations.len(), r.manifold_violations.len())
    }

    #[test]
    fn box_solid_is_a_closed_manifold() {
        let mut b = BRep::new();
        let s = build_box_solid(&mut b, [0.0, 0.0, 0.0], [0.3, 4.0, 2.0], Some("plate".into()));
        let sh = b.solids[s].shells[0];
        assert_eq!(b.shells[sh].faces.len(), 6);
        assert_eq!(gates(&b, sh), (0, 0, 0));
    }

    /// A plate pressed onto the cut terminal of a fused pair: every cut face of
    /// the lattice becomes a patch of the plate, with the same edges and the
    /// same area, and both solids stay closed.
    #[test]
    fn plate_on_a_cut_terminal_takes_the_lattice_cut_faces() {
        let fils = gen::pair_3mm();
        let planes = [
            Plane3::from_origin_normal(Point3::new(0.0, 1.0, 0.0), UnitVec3::Y),
            Plane3::from_origin_normal(Point3::new(0.0, 5.0, 0.0), -UnitVec3::Y),
        ];
        let part = fuse_and_trim(&fils, &FuseOptions::default(), &planes, 1.0e-6);
        assert!(part.is_watertight());
        let mut b = part.brep;
        let lattice = part.shell;
        // lattice extents across the plate
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for (_, v) in b.vertices.iter() {
            for k in 0..3 {
                lo[k] = lo[k].min([v.point.x, v.point.y, v.point.z][k]);
                hi[k] = hi[k].max([v.point.x, v.point.y, v.point.z][k]);
            }
        }
        let y = part.trims[0].plane.frame.origin.y;
        let m = 0.05;
        let plate = build_box_solid(&mut b, [lo[0] - m, y - 0.3, lo[2] - m], [hi[0] + m, y, hi[2] + m], Some("electrode".into()));
        let psh = b.solids[plate].shells[0];
        let contact = b.shells[psh].faces.iter().copied().find(|&f| matches!(outward_plane(&b, f), Some((_, n)) if n.dot_vec(UnitVec3::Y.as_vec()) > 0.5)).unwrap();
        let caps = coplanar_contact_faces(&b, lattice, contact, 1e-7);
        assert!(!caps.is_empty(), "the cut faces lie on the plate face");
        let cap_area: f64 = caps.iter().map(|&f| {
            let pts: Vec<(f64, f64)> = loop_points(&b, b.faces[f].outer_loop).iter().map(|p| (p.x, p.z)).collect();
            shoelace(&pts).abs()
        }).sum();
        let rep = imprint_coplanar(&mut b, contact, &caps, 1e-7).unwrap();
        assert_eq!(rep.patches, caps.len());
        assert!((rep.area - cap_area).abs() < 1e-9 * cap_area.max(1.0), "{} vs {}", rep.area, cap_area);
        assert_eq!(b.shells[psh].faces.len(), 6 + caps.len() + rep.holes);
        assert_eq!(gates(&b, psh), (0, 0, 0), "plate after imprint");
        assert_eq!(gates(&b, lattice), (0, 0, 0), "lattice untouched");
        // written to STEP: both solids closed, and the plate's contact edges
        // use the lattice's own curve entities (no second, refitted copy)
        let step = cadcore_step::brep_to_step(&b).unwrap();
        assert!(crate::validate::step_text::manifold_check(&step).is_empty(), "STEP manifold gate");
        let mut in_surface_curve = std::collections::HashSet::new();
        let mut direct = std::collections::HashSet::new();
        for line in step.lines() {
            let refs: Vec<&str> = line.split('#').skip(2).map(|t| t.split(|c: char| !c.is_ascii_digit()).next().unwrap()).collect();
            if line.contains("= SURFACE_CURVE(") {
                in_surface_curve.insert(refs[0].to_string());
            } else if line.contains("= EDGE_CURVE(") {
                direct.insert(refs[2].to_string());
            }
        }
        let shared = direct.intersection(&in_surface_curve).count();
        assert!(shared > 0, "the plate reuses the lattice's curves");
        // every patch has exactly the outline of its lattice face
        for &cap in &caps {
            let mut a: Vec<_> = loop_points(&b, b.faces[cap].outer_loop).iter().map(|p| (p.x.to_bits(), p.z.to_bits())).collect();
            a.sort();
            let found = b.shells[psh].faces.iter().any(|&f| {
                let mut q: Vec<_> = loop_points(&b, b.faces[f].outer_loop).iter().map(|p| (p.x.to_bits(), p.z.to_bits())).collect();
                q.sort();
                q == a
            });
            assert!(found, "patch for {cap:?}");
        }
    }

    /// A tool face that reaches the target's outline is refused and nothing changes.
    #[test]
    fn a_contact_reaching_the_plate_edge_is_refused() {
        let mut b = BRep::new();
        let big = build_box_solid(&mut b, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], None);
        let small = build_box_solid(&mut b, [0.0, 1.0, 0.2], [0.5, 1.5, 0.6], None);
        let bsh = b.solids[big].shells[0];
        let ssh = b.solids[small].shells[0];
        let top = b.shells[bsh].faces.iter().copied().find(|&f| matches!(outward_plane(&b, f), Some((_, n)) if n.dot_vec(UnitVec3::Y.as_vec()) > 0.5)).unwrap();
        let tools = coplanar_contact_faces(&b, ssh, top, 1e-9);
        assert_eq!(tools.len(), 1);
        let before = b.shells[bsh].faces.len();
        assert!(matches!(imprint_coplanar(&mut b, top, &tools, 1e-9), Err(ImprintError::NotInside { .. })));
        assert_eq!(b.shells[bsh].faces.len(), before);
        assert!(b.faces[top].inner_loops.is_empty());
    }
}
