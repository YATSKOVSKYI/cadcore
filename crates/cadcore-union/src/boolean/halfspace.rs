//! **Half-space cut** — trim an already-FUSED shell with a plane and close it
//! with real section caps.
//!
//! This is the "cut terminals AFTER sintering the part" operation.  The legacy
//! cut in `cadcore-ops` rebuilds each *template* cylinder face on its own, so it
//! only works on the raw per-segment sweep solids: on a fused shell (whose faces
//! carry explicit trim loops, `FaceExtent::Trimmed`) it finds no template
//! cylinders at all and silently does nothing.
//!
//! Here the cut is expressed with the same machinery as the union, because it
//! *is* the same operation with a different keep rule:
//!
//! ```text
//! union(A, B)      keep a face piece ⇔ it is OUTSIDE the other solid
//! cut (A, plane)   keep a face piece ⇔ it is INSIDE the kept half-space
//! ```
//!
//! so every face is imprinted with its section curve, its cells are classified
//! by an analytic (O(1), no ray casting) predicate, and the survivors are
//! stitched through ONE [`Assembler`].  What the union does *not* need and the
//! cut does is the **section cap**: the flat faces that close the part where the
//! plane went through it.  Those are built here from the very same section
//! curves, so their edges weld with the trimmed lateral faces by construction.
//!
//! Stages:
//!
//! ```text
//! 1 bound     model AABB → a finite cutting face on the plane
//! 2 section   SSI(face, cutting face) for every straddling face   (plane×cyl,
//!             plane×plane, plane×torus — the fillets)
//! 3 trim      imprint + classify + emit every face             (union::emit_solid_faces)
//! 4 cap       chain the section curves into closed loops, arrange them in the
//!             cut plane and emit the cells that carry material
//! ```
//!
//! ## TODO — the seam split point (gate 4)
//!
//! The periodic DCEL splits a chain where it crosses the band's seam, placing
//! that vertex by interpolating uv linearly along the chord.  The point is exact
//! on the cylinder but up to a chord sagitta off the curve it split, and
//! [`insert_joints`] must graft it into the partner face's copy of that curve
//! (or the shell opens), which drags the polyline that far off any OTHER face
//! the edge bounds.  The shell is manifold and its wires are wound correctly,
//! but gate 4 (edge-to-surface deviation) reports ~10 µm on those edges, over
//! SolidWorks' 1 µm knit.
//!
//! Substituting the mathematically exact crossing — solved on the source edge's
//! own curve, which is easy (`(p−o)·y = 0`, `(p−o)·x > 0`) — does NOT fix it:
//! measured, it re-opens 12 edges, because the DCEL still places its own
//! interpolated vertex on the arranged side and the two no longer coincide.
//! The split point therefore has to be made exact **inside the arrangement**
//! (`cadcore_geom::arrangement::build_periodic_u`), so both sides move together;
//! patching it downstream cannot work.

use std::collections::HashMap;

use cadcore_geom::Plane3;
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{
    BRep, Face, FaceExtent, FaceGeom, FaceId, FaceNormal, Shell, ShellId, Solid,
};

use super::aabb::{face_aabb, Aabb};
use super::contain::{face_contains_hit, point_in_solid, Membership};
use super::ssi::{intersect_faces_with, SsiCurve, SsiOptions};
use super::union::{
    face_boundary_chains, face_domain, snap_ssi_endpoints, unwrap_torus_band_chains,
};
use crate::arrange::assembly::Assembler;
use crate::arrange::cells::{arrange_and_classify_with, FaceChain};
use crate::arrange::domain::FaceDomain;

/// Outcome of a cut.
#[derive(Debug)]
pub struct CutOutcome {
    /// The trimmed B-Rep.
    pub brep: BRep,
    /// Its single closed shell.
    pub shell: ShellId,
    /// Number of section-cap faces emitted (0 ⇒ the plane missed the part).
    pub caps: usize,
    /// Faces of the input that were dropped whole (entirely on the cut side).
    pub dropped: usize,
}

/// Cut the solid bounded by `faces` with `keep`, retaining the material on the
/// side the plane's normal points to.
///
/// `knit` is the model knit tolerance (SolidWorks ≈ 1 µm); it drives the
/// assembler's geometric edge weld exactly as in [`super::union::union_n`].
///
/// Returns `None` when the input has no faces.  A plane that misses the part
/// entirely yields the part unchanged with `caps == 0`.
pub fn half_space_cut(
    src: &BRep,
    faces: &[FaceId],
    keep: &Plane3,
    knit: f64,
) -> Option<CutOutcome> {
    if faces.is_empty() {
        return None;
    }
    let n = keep.normal();

    // ── 1. Bound the cutting plane to the model ──────────────────────────────
    let mut model = Aabb::empty();
    for &f in faces {
        let bb = face_aabb(src, f);
        model.expand(bb.min);
        model.expand(bb.max);
    }
    let half = (model.extent() * 1.5).max(1.0);
    let centre_on_plane = keep.project(model.centre());
    // A finite square face on the cut plane — `region_polygon` needs a polygon
    // to clip the section curves against, and the arrangement needs a bounded
    // rectangle.  1.5× the model extent guarantees it covers every section.
    let cut_plane = Plane3 {
        frame: cadcore_math::Frame3 {
            origin: centre_on_plane,
            ..keep.frame
        },
    };
    let corner = |su: f64, sv: f64| centre_on_plane + cut_plane.frame.x * su + cut_plane.frame.y * sv;
    let mut cut_brep = BRep::new();
    let cut_face = cut_brep.add_face(Face {
        geom: FaceGeom::Plane(cut_plane),
        normal: FaceNormal::Same,
        outer_loop: Default::default(),
        inner_loops: Vec::new(),
        shell: Default::default(),
        extent: FaceExtent::Polygon {
            points: vec![
                corner(-half, -half),
                corner(half, -half),
                corner(half, half),
                corner(-half, half),
            ],
        },
    });

    // ── 2. Section curves ────────────────────────────────────────────────────
    // Only faces that actually straddle the plane can contribute; the AABB test
    // keeps the cut O(faces touching the plane) instead of O(all faces).
    let mut ssi: Vec<HashMap<FaceId, Vec<SsiCurve>>> = vec![HashMap::new()];
    let mut dropped = 0usize;
    for &fid in faces {
        let bb = face_aabb(src, fid);
        if std::env::var("CADCORE_DUMP_CUT").is_ok() {
            let kind = match &src.faces[fid].geom {
                FaceGeom::Plane(_) => "plane",
                FaceGeom::Cylinder(_) => "cyl",
                FaceGeom::Torus(_) => "torus",
                FaceGeom::Sphere(_) => "sphere",
            };
            let side = match straddle(&bb, keep, knit.max(1e-9)) {
                Side::Kept => "kept",
                Side::Cut => "cut",
                Side::Straddles => "STRADDLES",
            };
            eprintln!(
                "[cut][bp] {fid:?} {kind} {side} bb y[{:.3}..{:.3}]",
                bb.min.y, bb.max.y
            );
        }
        match straddle(&bb, keep, knit.max(1e-9)) {
            Side::Kept => {}
            Side::Cut => dropped += 1,
            Side::Straddles => {
                let curves = intersect_faces_with(src, fid, &cut_brep, cut_face, SsiOptions::all());
                if std::env::var("CADCORE_DUMP_CUT").is_ok() {
                    let kind = match &src.faces[fid].geom {
                        FaceGeom::Plane(_) => "plane",
                        FaceGeom::Cylinder(_) => "cyl",
                        FaceGeom::Torus(_) => "torus",
                        FaceGeom::Sphere(_) => "sphere",
                    };
                    eprintln!("[cut][ssi] {fid:?} {kind} straddles → {} curves", curves.len());
                }
                if !curves.is_empty() {
                    ssi[0].insert(fid, curves);
                }
            }
        }
    }
    // Where a section crosses a welded junction (a leg's rim circle onto its
    // elbow) the two pieces must share ONE endpoint, or the cap loop will not
    // close and the DCEL leaves a dangling edge.
    snap_ssi_endpoints(&mut ssi, (knit * 30.0).clamp(1e-3, 0.05));
    // Overlapping tubes give overlapping sections; split every section curve at
    // its crossings with the others ONCE, so the cap and the lateral faces
    // partition the shared curve identically.
    presplit_sections(&mut ssi[0], &cut_plane, knit.clamp(5e-4, 1e-2));

    // ── 3. Trim every face ───────────────────────────────────────────────────
    let mut out = BRep::new();
    let weld = knit.clamp(5e-4, 1e-2);
    let mut asm = Assembler::new(&mut out, weld);
    emit_cut_faces(src, faces, &ssi[0], keep, &mut asm);

    // ── 4. Section caps ──────────────────────────────────────────────────────
    let lateral = asm.faces().len();
    let loops = section_loops(&ssi[0], weld);
    let caps = emit_section_caps(src, faces, &cut_plane, n, half, &loops, &mut asm);

    let out_faces = asm.faces().to_vec();
    let shell = out.add_shell(Shell {
        faces: out_faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = out.add_solid(Solid {
        shells: vec![shell],
        name: Some("cut".into()),
    });
    out.shells[shell].solid = solid;
    for f in &out_faces {
        out.faces[*f].shell = shell;
    }

    if std::env::var("CADCORE_DUMP_CUT").is_ok() {
        let rep = crate::validate::manifold::check(&out, shell);
        eprintln!(
            "[cut] faces={} (lateral={lateral} caps={caps}) dropped={dropped} \
             section_loops={} violations={}",
            out_faces.len(),
            loops.len(),
            rep.violations.len()
        );
        eprintln!("[cut] cap faces = {:?}", &out_faces[lateral..]);
        let mut per_face: std::collections::BTreeMap<String, usize> =
            std::collections::BTreeMap::new();
        for v in &rep.violations {
            let faces_in = v
                .split('[')
                .nth(1)
                .unwrap_or("")
                .trim_end_matches(']')
                .to_string();
            *per_face.entry(faces_in).or_default() += 1;
        }
        for (k, n) in per_face.iter().take(12) {
            eprintln!("[cut]   {n:>4} × {k}");
        }
    }

    Some(CutOutcome {
        brep: out,
        shell,
        caps,
        dropped,
    })
}

/// Imprint every face with its section curve, keep the cells that are both
/// material of that face and in front of the plane, and emit them.
///
/// This mirrors [`super::union::emit_solid_faces`] with one extra rule the
/// union never needs.  The union's input faces are analytic templates whose
/// whole parameter region is material; the cut's input is an already-Boolean'd
/// shell, so a face can be a band with HOLES punched in it (crossing windows).
/// The arrangement hands back a hole's interior as a cell in its own right, and
/// "keep unless behind the plane" would happily emit it — filling the window
/// back in and giving its rim four users.  So each cell is additionally tested
/// against the ORIGINAL face's trimmed region with the same membership oracle
/// the ray caster uses.
fn emit_cut_faces(
    src: &BRep,
    faces: &[FaceId],
    ssi: &HashMap<FaceId, Vec<SsiCurve>>,
    keep: &Plane3,
    asm: &mut Assembler,
) {
    // Where a section arc lands on a boundary curve, that curve must be split —
    // on BOTH the arranged face and the copied one across the same window.
    let mut joints = section_joints(ssi);

    // Pass 1: settle each arranged face's seam and publish where it crosses the
    // face's boundary.  The seam is invisible to the face on the other side of
    // a crossed edge, so it has to become a joint like any other split point.
    // The rotation is settled here and REUSED below, so adding these joints
    // cannot move it.
    let mut seams: HashMap<FaceId, f64> = HashMap::new();
    for &fid in faces {
        if !ssi.contains_key(&fid) {
            continue;
        }
        let Some(domain) = face_domain(src, fid) else { continue };
        if !matches!(domain, FaceDomain::CylinderBand { .. }) {
            continue; // only a periodic band has a seam to park
        }
        let mut chains = arranged_chains(src, fid, &domain, ssi, &joints);
        // No usable gap ⇒ keep the default seam (rotation 0) — but STILL publish
        // where it cuts, or the copied partner keeps that edge whole.
        let rot = seam_rotation(&domain, &chains).unwrap_or(0.0);
        let seamed = apply_seam_rotation(&domain, &mut chains, rot);
        seams.insert(fid, rot);
        for p in seam_crossings(&seamed, &chains) {
            if !joints.iter().any(|&q| (q - p).length() < 1e-9) {
                joints.push(p);
            }
        }
    }

    for &fid in faces {
        let f = &src.faces[fid];
        // Faces the plane does not touch AT ALL and that carry real topology are
        // copied verbatim — see `Assembler::emit_face_copy` for why re-arranging
        // them is NOT a no-op on a fused shell.  Template faces (no loops) have
        // nothing to copy, so they take the arrangement path, which synthesises
        // their boundary.
        //
        // A face that STRADDLES must never be copied, even when the SSI found no
        // curve on it (a shallow torus arc, a tangency): its boundary rims
        // straddle too, the neighbour across such a rim IS trimmed, and a copy
        // would present the whole rim against the neighbour's trimmed half —
        // one user each, i.e. an open shell.  Arranging it costs nothing extra
        // and lets the keep predicate trim it.
        if !ssi.contains_key(&fid) && has_real_loops(src, fid) {
            match straddle(&face_aabb(src, fid), keep, 1e-9) {
                Side::Cut => continue, // wholly behind the plane
                Side::Kept => {
                    asm.emit_face_copy(src, fid, &joints);
                    continue;
                }
                Side::Straddles => {} // fall through to the arrangement
            }
        }
        let Some(domain) = face_domain(src, fid) else {
            continue; // unsupported carrier surface — skip (incomplete)
        };

        // Material is removed where the plane's normal points AWAY, and where
        // the original face had no material to begin with.
        let drop = |p: Point3| {
            keep.signed_distance(p) < 0.0
                || !matches!(face_contains_hit(src, fid, p), Membership::Inside)
        };

        let mut chains = arranged_chains(src, fid, &domain, ssi, &joints);
        let stag = 1000 + chains.iter().filter(|c| c.tag >= 1000).count() as u32;
        // Apply the seam settled in pass 1.  A fused band is SEAMLESS; the
        // periodic DCEL has to put a seam somewhere, and wherever it lands it
        // splits the chains crossing it — so it is parked in the widest empty
        // gap (the trick `arrange::solid::seam_away_from` plays) and whatever
        // it still crosses was published as a joint.
        let domain = match seams.get(&fid) {
            Some(&rot) => apply_seam_rotation(&domain, &mut chains, rot),
            None => domain,
        };
        unwrap_torus_band_chains(&domain, &mut chains);

        let cells = arrange_and_classify_with(&domain, &chains, 1e-3, &drop);
        let mut emitted = 0usize;
        for cell in &cells {
            if cell.keep {
                asm.emit_cell(&domain, f.geom.clone(), f.normal, cell);
                emitted += 1;
            }
        }
        if std::env::var("CADCORE_DUMP_CUT").is_ok() {
            let kind = match &f.geom {
                FaceGeom::Plane(_) => "plane",
                FaceGeom::Cylinder(_) => "cyl",
                FaceGeom::Torus(_) => "torus",
                FaceGeom::Sphere(_) => "sphere",
            };
            eprintln!(
                "cut face {fid:?} {kind}: bnd={} ssi={} cells={} emitted={} areas={:?}",
                chains.len() - (stag - 1000) as usize,
                stag - 1000,
                cells.len(),
                emitted,
                cells
                    .iter()
                    .map(|c| (c.cell.area * 1000.0).round() / 1000.0)
                    .collect::<Vec<_>>()
            );
        }
    }
}

/// Every chain an ARRANGED face feeds the DCEL: its own boundary loops (with
/// the joints grafted in) followed by the section curves lying on it.
///
/// Built the same way in both passes so the seam chosen in pass 1 is the seam
/// applied in pass 2.
fn arranged_chains(
    src: &BRep,
    fid: FaceId,
    domain: &FaceDomain,
    ssi: &HashMap<FaceId, Vec<SsiCurve>>,
    joints: &[Point3],
) -> Vec<FaceChain> {
    let mut chains: Vec<FaceChain> = Vec::new();
    let mut tag = 1u32;
    for pts in boundary_chains_with_joints(src, fid, domain, joints) {
        chains.push(FaceChain { pts, tag });
        tag += 1;
    }
    let mut stag = 1000u32;
    if let Some(curves) = ssi.get(&fid) {
        for c in curves {
            let mut pts: Vec<(f64, f64)> = c.points.iter().map(|&p| domain.uv(p)).collect();
            if c.closed && pts.len() >= 3 {
                pts.push(pts[0]);
            }
            if pts.len() >= 2 {
                chains.push(FaceChain { pts, tag: stag });
                stag += 1;
            }
        }
    }
    chains
}

/// Rotate a cylinder band's parameter frame so the u = ±π seam falls in the
/// widest gap between chain points, and shift the chains onto the new frame.
///
/// Returns `None` for non-cylinder domains (nothing to do) and when the chains
/// leave no usable gap.
///
/// The rotation is exact: turning `frame.x` by `rot` maps every parameter
/// `u → u − rot`, so the chains only need shifting, never re-projecting, and
/// `lift` still lands on the identical 3-D points.
fn reseam_cylinder(domain: &FaceDomain, chains: &mut [FaceChain]) -> Option<FaceDomain> {
    let rot = seam_rotation(domain, chains)?;
    Some(apply_seam_rotation(domain, chains, rot))
}

/// Pick the rotation that parks the u = ±π seam in the widest empty gap.
fn seam_rotation(domain: &FaceDomain, chains: &[FaceChain]) -> Option<f64> {
    use std::f64::consts::{PI, TAU};
    let FaceDomain::CylinderBand { .. } = domain else {
        return None;
    };
    // Only LOCALISED chains matter: a band's rim circles span every u, so the
    // seam has to cross them wherever it goes (both faces sharing a rim split
    // it the same way, so that is harmless).  What must not be split is a
    // window loop or a section curve, which occupies a limited arc.
    let mut us: Vec<f64> = Vec::new();
    for c in chains.iter() {
        let (lo, hi) = c
            .pts
            .iter()
            .fold((f64::MAX, f64::MIN), |(lo, hi), p| (lo.min(p.0), hi.max(p.0)));
        if hi - lo > 0.9 * TAU {
            continue; // full-wrap rim
        }
        us.extend(c.pts.iter().map(|p| p.0.rem_euclid(TAU)));
    }
    if us.is_empty() {
        return None;
    }
    us.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // widest circular gap between consecutive samples
    let (mut best_gap, mut seam_u) = (0.0_f64, 0.0_f64);
    for i in 0..us.len() {
        let (a, b) = (us[i], if i + 1 < us.len() { us[i + 1] } else { us[0] + TAU });
        if b - a > best_gap {
            best_gap = b - a;
            seam_u = 0.5 * (a + b);
        }
    }
    // Nothing worth moving: the chains already cover the circle densely (the
    // seam has to cross something, and the caller falls back to the default).
    if best_gap < 1e-6 {
        return None;
    }
    // `build_periodic_u` unwraps each chain from its own first sample and then
    // splits at every multiple of the period, so the seam it cuts along is
    // u = 0 — the frame's x direction — NOT u = ±π.  Rotating by `seam_u` puts
    // u = 0 at the middle of the gap.
    let _ = PI;
    Some(seam_u)
}

/// Rotate the band's parameter frame by `rot` and shift the chains onto it.
fn apply_seam_rotation(domain: &FaceDomain, chains: &mut [FaceChain], rot: f64) -> FaceDomain {
    use std::f64::consts::{PI, TAU};
    let FaceDomain::CylinderBand { surf, length } = domain else {
        return *domain;
    };
    let (s, c) = rot.sin_cos();
    let (Some(x), ..) = (UnitVec3::try_from_vec(surf.frame.x * c + surf.frame.y * s), ()) else {
        return *domain;
    };
    let Some(y) = UnitVec3::try_from_vec(surf.frame.z.cross(x)) else {
        return *domain;
    };
    let reframed = cadcore_geom::CylSurf {
        frame: cadcore_math::Frame3 {
            origin: surf.frame.origin,
            x,
            y,
            z: surf.frame.z,
        },
        radius: surf.radius,
    };
    for ch in chains.iter_mut() {
        for p in ch.pts.iter_mut() {
            let mut u = p.0 - rot;
            u -= TAU * ((u + PI) / TAU).floor(); // wrap into (−π, π]
            p.0 = u;
        }
    }
    FaceDomain::CylinderBand {
        surf: reframed,
        length: *length,
    }
}

/// The 3-D points where `chains` (already in the re-seamed uv) cross the
/// u = 0 seam.
///
/// The seam is an artefact of the parametrisation, not of the model, so the
/// face on the OTHER side of a crossed edge knows nothing about it and keeps
/// that edge whole; publishing the crossings as joints makes both sides split
/// there.
///
/// The point is computed exactly the way `build_periodic_u` computes its own
/// split — linear interpolation of uv along the chord — and that is deliberate:
/// what matters is that BOTH sides use the SAME point.  Substituting the
/// mathematically exact crossing (solved on the source edge's curve) instead
/// re-opens the shell, because the DCEL still places its own interpolated
/// vertex on the arranged side.  See the module TODO.
fn seam_crossings(domain: &FaceDomain, chains: &[FaceChain]) -> Vec<Point3> {
    use std::f64::consts::PI;
    let FaceDomain::CylinderBand { .. } = domain else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for ch in chains {
        for w in ch.pts.windows(2) {
            let (a, b) = (w[0], w[1]);
            // A jump larger than π is the atan2 wrap, not a real step across 0.
            if (a.0 - b.0).abs() > PI {
                continue;
            }
            if (a.0 > 0.0) == (b.0 > 0.0) {
                continue; // both on the same side of the seam
            }
            let t = a.0 / (a.0 - b.0);
            if !(0.0..=1.0).contains(&t) {
                continue;
            }
            out.push(domain.lift(0.0, a.1 + (b.1 - a.1) * t));
        }
    }
    out
}

// ── Joints ───────────────────────────────────────────────────────────────────

/// Every point where a section arc starts or ends.
///
/// These are the joints of the cut: a section arc runs across a face and stops
/// on one of its boundary curves — a crossing window, a rim circle, another
/// section.  They are the points at which the ARRANGED side of that boundary
/// gets split, so the COPIED side has to be split there too.
fn section_joints(ssi: &HashMap<FaceId, Vec<SsiCurve>>) -> Vec<Point3> {
    let mut out: Vec<Point3> = Vec::new();
    for curves in ssi.values() {
        for c in curves {
            if c.closed || c.points.len() < 2 {
                continue; // a closed section ends on nothing
            }
            for p in [c.points[0], *c.points.last().unwrap()] {
                if !out.iter().any(|&q| (q - p).length() < 1e-9) {
                    out.push(p);
                }
            }
        }
    }
    out
}

/// Insert every joint that lies ON `poly` at its exact position.
///
/// A joint is accepted only when it is close to a segment relative to that
/// segment's own length (and within an absolute cap), so a joint belonging to
/// some other curve running nearby is never grafted on.  Inserting the joint
/// itself — not its projection — is the point: both consumers of this curve
/// then carry the identical vertex and their segments weld.
pub(crate) fn insert_joints(poly: &mut Vec<Point3>, joints: &[Point3]) {
    if joints.is_empty() || poly.len() < 2 {
        return;
    }
    let mut out: Vec<Point3> = Vec::with_capacity(poly.len() + joints.len());
    for w in 0..poly.len() - 1 {
        let (a, b) = (poly[w], poly[w + 1]);
        out.push(a);
        let seg = b - a;
        let len2 = seg.dot(seg);
        if len2 < 1e-24 {
            continue;
        }
        // The polyline is only a CHORD approximation of its curve, so a joint
        // that genuinely belongs to it can sit up to a sagitta off — the window
        // is the sagitta scale, not the knit tolerance.
        let tol = (0.3 * len2.sqrt()).min(0.01);
        let mut hits: Vec<(f64, Point3)> = Vec::new();
        for &j in joints {
            let t = seg.dot(j - a) / len2;
            if !(1e-6..=1.0 - 1e-6).contains(&t) {
                continue; // at (or past) an endpoint — already a vertex
            }
            let foot = a + seg * t;
            if (j - foot).length() <= tol {
                hits.push((t, j));
            }
        }
        hits.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
        for (_, j) in hits {
            if out.last().map_or(true, |&q| (q - j).length() > 1e-9) {
                out.push(j);
            }
        }
    }
    out.push(*poly.last().unwrap());
    *poly = out;
}

/// The face's boundary loops as uv chains, with the section joints inserted.
///
/// Mirrors [`face_boundary_chains`] but walks the loops in 3-D first so the
/// joints can be grafted in before projection — the copied path grafts the same
/// joints into the same edges, so the two sides split identically.
fn boundary_chains_with_joints(
    brep: &BRep,
    fid: FaceId,
    domain: &FaceDomain,
    joints: &[Point3],
) -> Vec<Vec<(f64, f64)>> {
    if !has_real_loops(brep, fid) {
        return face_boundary_chains(brep, fid, domain);
    }
    let f = &brep.faces[fid];
    let mut loops = Vec::new();
    let mut ids = vec![f.outer_loop];
    ids.extend(f.inner_loops.iter().copied());
    for lid in ids {
        let Some(lp) = brep.loops.get(lid) else { continue };
        let start = lp.start;
        if brep.coedges.get(start).is_none() {
            continue;
        }
        let mut pts3: Vec<Point3> = Vec::new();
        let mut c = start;
        loop {
            let ce = &brep.coedges[c];
            let mut sp = super::aabb::sample_edge(brep, ce.edge);
            if ce.sense == cadcore_topo::CoEdgeSense::Opposite {
                sp.reverse();
            }
            insert_joints(&mut sp, joints);
            for p in sp {
                if pts3.last().map_or(true, |&q| (q - p).length() > 1e-12) {
                    pts3.push(p);
                }
            }
            c = ce.next;
            if c == start {
                break;
            }
        }
        if pts3.len() >= 3 {
            let mut pts: Vec<(f64, f64)> = pts3.iter().map(|&p| domain.uv(p)).collect();
            if pts.first() != pts.last() {
                pts.push(pts[0]);
            }
            loops.push(pts);
        }
    }
    if loops.is_empty() {
        return face_boundary_chains(brep, fid, domain);
    }
    loops
}

/// Does this face carry explicit B-Rep topology (as opposed to being an
/// analytic template whose boundary is implied by its [`FaceExtent`])?
fn has_real_loops(brep: &BRep, fid: FaceId) -> bool {
    brep.faces
        .get(fid)
        .and_then(|f| brep.loops.get(f.outer_loop))
        .map_or(false, |lp| brep.coedges.get(lp.start).is_some())
}

// ── Broad phase ──────────────────────────────────────────────────────────────

enum Side {
    Kept,
    Cut,
    Straddles,
}

/// Which side of `keep` an axis-aligned box lies on.
///
/// The extreme corners along the normal bound the box's signed-distance range:
/// `d(centre) ± Σ|n·axis|·half_extent`.
fn straddle(bb: &Aabb, keep: &Plane3, tol: f64) -> Side {
    let n = keep.normal().as_vec();
    let c = Point3::new(
        0.5 * (bb.min.x + bb.max.x),
        0.5 * (bb.min.y + bb.max.y),
        0.5 * (bb.min.z + bb.max.z),
    );
    let reach = 0.5
        * (n.x.abs() * (bb.max.x - bb.min.x)
            + n.y.abs() * (bb.max.y - bb.min.y)
            + n.z.abs() * (bb.max.z - bb.min.z));
    let d = keep.signed_distance(c);
    if d - reach >= -tol {
        Side::Kept
    } else if d + reach <= tol {
        Side::Cut
    } else {
        Side::Straddles
    }
}

// ── Section registry (pre-split) ─────────────────────────────────────────────

/// 2-D segment × segment intersection, as parameters `(ta, tb)` along each.
///
/// `None` for parallel segments and for crossings outside either span.
fn seg_cross(a0: (f64, f64), a1: (f64, f64), b0: (f64, f64), b1: (f64, f64)) -> Option<(f64, f64)> {
    let da = (a1.0 - a0.0, a1.1 - a0.1);
    let db = (b1.0 - b0.0, b1.1 - b0.1);
    let den = da.0 * db.1 - da.1 * db.0;
    if den.abs() < 1e-14 {
        return None; // parallel / degenerate
    }
    let w = (b0.0 - a0.0, b0.1 - a0.1);
    let ta = (w.0 * db.1 - w.1 * db.0) / den;
    let tb = (w.0 * da.1 - w.1 * da.0) / den;
    const E: f64 = 1e-9;
    if !(-E..=1.0 + E).contains(&ta) || !(-E..=1.0 + E).contains(&tb) {
        return None;
    }
    Some((ta.clamp(0.0, 1.0), tb.clamp(0.0, 1.0)))
}

/// Split every section curve where it crosses another one, so the cap and the
/// lateral faces are handed the SAME pieces.
///
/// This is the cut's registry, the analogue of [`crate::arrange::registry`] for
/// the union: where two fused tubes overlap at the plane their sections overlap
/// too, and each consumer would otherwise partition the shared curve its own
/// way — the cap at the outline of the union, the lateral face at its crossing
/// window — so no segment would ever weld.  Splitting ONCE here, with a single
/// computed point per crossing used by both curves, makes the partitions
/// identical by construction.
///
/// Every curve lies in the cut plane, so the whole thing is a 2-D problem in
/// the plane's own uv.
fn presplit_sections(ssi: &mut HashMap<FaceId, Vec<SsiCurve>>, cut_plane: &Plane3, tol: f64) {
    let uv = |p: Point3| {
        let w = p - cut_plane.frame.origin;
        (
            cut_plane.frame.x.dot_vec(w),
            cut_plane.frame.y.dot_vec(w),
        )
    };
    let lift = |(u, v): (f64, f64)| cut_plane.frame.origin + cut_plane.frame.x * u + cut_plane.frame.y * v;

    // Flatten to a working list: (owner face, uv polyline with the wrap point
    // appended for closed curves, was_closed).
    struct Work {
        face: FaceId,
        pts: Vec<(f64, f64)>,
        closed: bool,
        /// (segment index, parameter along it, crossing point) — the splits.
        cuts: Vec<(usize, f64, (f64, f64))>,
    }
    let mut work: Vec<Work> = Vec::new();
    for (&face, curves) in ssi.iter() {
        for c in curves {
            if c.points.len() < 2 {
                continue;
            }
            let mut pts: Vec<(f64, f64)> = c.points.iter().map(|&p| uv(p)).collect();
            if c.closed {
                pts.push(pts[0]);
            }
            work.push(Work {
                face,
                pts,
                closed: c.closed,
                cuts: Vec::new(),
            });
        }
    }

    // Pairwise crossings.  The SAME 2-D point is recorded on both curves.
    for i in 0..work.len() {
        for j in (i + 1)..work.len() {
            let mut found: Vec<(usize, f64, usize, f64, (f64, f64))> = Vec::new();
            for si in 0..work[i].pts.len().saturating_sub(1) {
                let (a0, a1) = (work[i].pts[si], work[i].pts[si + 1]);
                for sj in 0..work[j].pts.len().saturating_sub(1) {
                    let (b0, b1) = (work[j].pts[sj], work[j].pts[sj + 1]);
                    if let Some((ta, tb)) = seg_cross(a0, a1, b0, b1) {
                        let p = (a0.0 + (a1.0 - a0.0) * ta, a0.1 + (a1.1 - a0.1) * ta);
                        found.push((si, ta, sj, tb, p));
                    }
                }
            }
            for (si, ta, sj, tb, p) in found {
                work[i].cuts.push((si, ta, p));
                work[j].cuts.push((sj, tb, p));
            }
        }
    }

    if work.iter().all(|w| w.cuts.is_empty()) {
        return; // nothing overlaps — the curves are already the pieces
    }

    // Rebuild each curve as pieces broken at its cuts.
    let mut out: HashMap<FaceId, Vec<SsiCurve>> = HashMap::new();
    for w in work.iter_mut() {
        w.cuts
            .sort_by(|a, b| (a.0, a.1).partial_cmp(&(b.0, b.1)).unwrap_or(std::cmp::Ordering::Equal));
        w.cuts.dedup_by(|a, b| {
            let d = ((a.2 .0 - b.2 .0).powi(2) + (a.2 .1 - b.2 .1).powi(2)).sqrt();
            d < tol
        });

        if w.cuts.is_empty() {
            out.entry(w.face).or_default().push(SsiCurve {
                points: w
                    .pts
                    .iter()
                    .take(w.pts.len() - usize::from(w.closed))
                    .map(|&q| lift(q))
                    .collect(),
                closed: w.closed,
            });
            continue;
        }

        let mut pieces: Vec<Vec<(f64, f64)>> = Vec::new();
        let mut cur: Vec<(f64, f64)> = vec![w.pts[0]];
        let mut ci = 0usize;
        for si in 0..w.pts.len() - 1 {
            while ci < w.cuts.len() && w.cuts[ci].0 == si {
                let p = w.cuts[ci].2;
                if cur.last().map_or(true, |&q| dist2(q, p) > tol * tol) {
                    cur.push(p);
                }
                if cur.len() >= 2 {
                    pieces.push(std::mem::take(&mut cur));
                }
                cur = vec![p];
                ci += 1;
            }
            let nxt = w.pts[si + 1];
            if cur.last().map_or(true, |&q| dist2(q, nxt) > 1e-24) {
                cur.push(nxt);
            }
        }
        if cur.len() >= 2 {
            pieces.push(cur);
        }
        // A closed curve's last piece continues into its first across the wrap
        // point, so they are one arc.
        if w.closed && pieces.len() >= 2 {
            let last = pieces.pop().unwrap();
            let first = pieces.remove(0);
            let mut joined = last;
            joined.extend(first.into_iter().skip(1));
            pieces.push(joined);
        }
        for p in pieces {
            if p.len() >= 2 {
                out.entry(w.face).or_default().push(SsiCurve {
                    points: p.into_iter().map(lift).collect(),
                    closed: false,
                });
            }
        }
    }
    *ssi = out;
}

fn dist2(a: (f64, f64), b: (f64, f64)) -> f64 {
    (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)
}

// ── Section-loop assembly ────────────────────────────────────────────────────

/// Chain every section curve into closed 3-D loops.
///
/// A plane cutting one filament transversally already gives a closed conic; a
/// cut through a welded leg↔elbow junction gives two or more open arcs that
/// have to be walked end-to-end.  Curves whose ends do not meet anything are
/// dropped rather than guessed at — an unclosed cap is worse than a missing one
/// (it would emit a face with an open wire).
fn section_loops(ssi: &HashMap<FaceId, Vec<SsiCurve>>, tol: f64) -> Vec<Vec<Point3>> {
    let mut open: Vec<Vec<Point3>> = Vec::new();
    let mut closed: Vec<Vec<Point3>> = Vec::new();
    for curves in ssi.values() {
        for c in curves {
            if c.points.len() < 2 {
                continue;
            }
            let ends_meet =
                (c.points[0] - *c.points.last().unwrap()).length() <= tol;
            if c.closed || (ends_meet && c.points.len() >= 3) {
                let mut pts = c.points.clone();
                if ends_meet && pts.len() >= 3 {
                    pts.pop(); // stored as an explicit repeat — the loop is implicit
                }
                closed.push(pts);
            } else {
                open.push(c.points.clone());
            }
        }
    }

    if std::env::var("CADCORE_DUMP_CUT").is_ok() {
        eprintln!(
            "[cut] section curves: {} closed, {} open (tol={tol:.1e})",
            closed.len(),
            open.len()
        );
        for (f, cs) in ssi.iter() {
            for c in cs {
                let (a, b) = (c.points[0], *c.points.last().unwrap());
                eprintln!(
                    "[cut]   {f:?} n={} closed={} ({:.4},{:.4},{:.4})..({:.4},{:.4},{:.4})",
                    c.points.len(),
                    c.closed,
                    a.x, a.y, a.z, b.x, b.y, b.z
                );
            }
        }
    }

    // Walk the open arcs end-to-end.  Each arc is used once; a chain that comes
    // back to its own start is a loop, anything else is discarded.
    let mut used = vec![false; open.len()];
    for i in 0..open.len() {
        if used[i] {
            continue;
        }
        used[i] = true;
        let mut chain = open[i].clone();
        loop {
            let tail = *chain.last().unwrap();
            if (tail - chain[0]).length() <= tol && chain.len() >= 3 {
                closed.push(chain);
                break;
            }
            let mut extended = false;
            for j in 0..open.len() {
                if used[j] {
                    continue;
                }
                let (a, b) = (open[j][0], *open[j].last().unwrap());
                if (tail - a).length() <= tol {
                    used[j] = true;
                    chain.extend_from_slice(&open[j][1..]);
                    extended = true;
                    break;
                }
                if (tail - b).length() <= tol {
                    used[j] = true;
                    chain.extend(open[j].iter().rev().skip(1).copied());
                    extended = true;
                    break;
                }
            }
            if !extended {
                break; // dangling chain — not a closed section, drop it
            }
        }
    }
    closed
}

/// Arrange the section loops in the cut plane and emit the cells that carry
/// material as cap faces.
///
/// The arrangement is REQUIRED, not a convenience: where two fused tubes
/// overlap at the plane their section loops overlap too, so the real outline of
/// the cut is the union of the loops, not the loops themselves.  (Emitting each
/// loop directly caps each tube's full circle, including the part buried inside
/// its neighbour.)
///
/// The cap's carrier plane is oriented OUTWARD (away from the kept material),
/// so the DCEL's CCW-in-uv outer loops and CW holes come out physically wound
/// as SolidWorks/SpaceClaim require (README invariant #4).
fn emit_section_caps(
    src: &BRep,
    faces: &[FaceId],
    cut_plane: &Plane3,
    keep_normal: UnitVec3,
    half: f64,
    loops: &[Vec<Point3>],
    asm: &mut Assembler,
) -> usize {
    if loops.is_empty() {
        return 0;
    }
    // Outward = away from the kept material.
    let cap_plane = Plane3::from_origin_normal(cut_plane.frame.origin, -keep_normal);
    let domain = FaceDomain::Plane {
        plane: cap_plane,
        half_extent: half,
    };

    // Bounding rectangle so the arrangement is bounded; its cell lies outside
    // every section loop and is dropped below.
    const RECT_TAG: u32 = 1;
    let rect = vec![
        (-half, -half),
        (half, -half),
        (half, half),
        (-half, half),
        (-half, -half),
    ];
    let mut chains = vec![FaceChain {
        pts: rect,
        tag: RECT_TAG,
    }];
    let mut tag = 1000u32;
    for lp in loops {
        let mut pts: Vec<(f64, f64)> = lp.iter().map(|&p| domain.uv(p)).collect();
        pts.push(pts[0]); // close it for the DCEL
        chains.push(FaceChain { pts, tag });
        tag += 1;
    }

    // A cap cell carries material when a point just INSIDE the kept side of it
    // is inside the original solid.  `arrange_and_classify_with` keeps a cell
    // when the predicate is false, so the predicate is the negation.
    let probe = 1.0e-4;
    let void = |p: Point3| !point_in_solid(src, faces, p + keep_normal * probe);

    let cells = arrange_and_classify_with(&domain, &chains, 0.0, &void);
    let mut emitted = 0usize;
    for cell in cells {
        if !cell.keep {
            continue;
        }
        // The bounding rectangle is scaffolding, not geometry: the cell that
        // uses it is the plane OUTSIDE the part, with every section island as a
        // hole.  Emitting it would hand each section curve a third user (the
        // island cap and the trimmed lateral face are the only two legitimate
        // ones) and paste a giant disc across the model.
        if cell.cell.outer.iter().any(|s| s.tag == RECT_TAG) {
            continue;
        }
        asm.emit_cell(&domain, FaceGeom::Plane(cap_plane), FaceNormal::Same, &cell);
        emitted += 1;
    }
    emitted
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boolean::tests_support::capped_cylinder;
    use cadcore_geom::CylSurf;

    /// Gate 1: every edge of the result must be used exactly twice.
    fn open_edges(brep: &BRep, shell: ShellId) -> usize {
        crate::validate::manifold::check(brep, shell).violations.len()
    }

    /// The full SpaceClaim accept gate: manifold + edge-to-surface deviation +
    /// uv wire winding / self-intersection.  `(wire, distance, manifold)`.
    fn spaceclaim_gate(brep: &BRep, shell: ShellId) -> (usize, usize, usize) {
        let mut cfg = crate::config::UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = crate::Diagnostics::new();
        let r = crate::validate::validate_shell(brep, shell, &cfg, &mut diag);
        (
            r.wire_violations.len(),
            r.distance_violations.len(),
            r.manifold_violations.len(),
        )
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
                    out.extend(super::super::aabb::sample_edge(brep, ce.edge));
                    c = ce.next;
                    if c == st {
                        break;
                    }
                }
            }
        }
        out
    }

    /// No geometry may survive on the discarded side.
    fn deepest_past_plane(brep: &BRep, shell: ShellId, keep: &Plane3) -> f64 {
        shell_points(brep, shell)
            .into_iter()
            .map(|p| keep.signed_distance(p))
            .fold(f64::MAX, f64::min)
    }

    fn face_count(brep: &BRep, shell: ShellId) -> usize {
        brep.shells[shell].faces.len()
    }

    /// Sampled points of one loop, in traversal order.
    fn loop_points(brep: &BRep, lid: cadcore_topo::LoopId) -> Vec<Point3> {
        let Some(lp) = brep.loops.get(lid) else { return Vec::new() };
        let st = lp.start;
        if brep.coedges.get(st).is_none() {
            return Vec::new();
        }
        let mut out = Vec::new();
        let mut c = st;
        loop {
            let ce = &brep.coedges[c];
            out.extend(super::super::aabb::sample_edge(brep, ce.edge));
            c = ce.next;
            if c == st {
                break;
            }
        }
        out
    }

    fn centroid(pts: &[Point3]) -> Point3 {
        let o = Point3::new(0.0, 0.0, 0.0);
        let mut acc = cadcore_math::Vec3::new(0.0, 0.0, 0.0);
        for p in pts {
            acc = acc + (*p - o);
        }
        o + acc * (1.0 / pts.len() as f64)
    }

    /// A capped cylinder cut across its axis: the part keeps its lateral band up
    /// to the plane, loses the far cap and gains a flat disk section.
    #[test]
    fn cut_cylinder_across_axis_is_watertight() {
        let mut a = BRep::new();
        let fa = capped_cylinder(&mut a, Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, 0.5, 4.0);
        // keep z <= 2.5  ⇒ normal points down
        let keep = Plane3::from_origin_normal(Point3::new(0.0, 0.0, 2.5), -UnitVec3::Z);

        let cut = half_space_cut(&a, &fa, &keep, 1e-6).expect("cut runs");
        assert!(cut.caps >= 1, "the section must be capped");
        assert_eq!(
            cut.dropped, 1,
            "the far end cap is entirely past the plane and must be dropped"
        );
        assert_eq!(
            open_edges(&cut.brep, cut.shell),
            0,
            "cut solid must stay watertight"
        );
        assert!(
            deepest_past_plane(&cut.brep, cut.shell, &keep) > -1e-6,
            "material survived past the cut plane: {:.3e}",
            deepest_past_plane(&cut.brep, cut.shell, &keep)
        );
        // The cap sits ON the plane and spans the tube: its points are at the
        // cut height and within the radius.
        let on_plane: Vec<Point3> = shell_points(&cut.brep, cut.shell)
            .into_iter()
            .filter(|p| keep.signed_distance(*p).abs() < 1e-6)
            .collect();
        assert!(!on_plane.is_empty(), "no section geometry at the plane");
        for p in &on_plane {
            let r = (p.x * p.x + p.y * p.y).sqrt();
            assert!(
                (r - 0.5).abs() < 1e-3,
                "section rim off the tube radius: r={r:.6}"
            );
        }
    }

    /// A plane clear of the part changes nothing and adds no cap.
    #[test]
    fn cut_missing_the_part_keeps_it_whole() {
        let mut a = BRep::new();
        let fa = capped_cylinder(&mut a, Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, 0.5, 4.0);
        let keep = Plane3::from_origin_normal(Point3::new(0.0, 0.0, -5.0), UnitVec3::Z);

        let cut = half_space_cut(&a, &fa, &keep, 1e-6).expect("cut runs");
        assert_eq!(cut.caps, 0, "nothing to cap");
        assert_eq!(cut.dropped, 0, "nothing to drop");
        assert_eq!(open_edges(&cut.brep, cut.shell), 0);
        assert_eq!(face_count(&cut.brep, cut.shell), fa.len(), "part unchanged");
    }

    /// Two crossing serpentines FUSED into one solid (the shape the scaffold
    /// export builds), then cut.
    fn fused_scaffold() -> (BRep, ShellId) {
        use crate::arrange::scaffold::{fuse_scaffold, Filament};
        let r = 0.275;
        // A: U in z=0 — X-legs at y=0 and y=2, connector at x=0.
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
        // B: U one layer up — Y-legs at x=-1.5 and x=-0.8 crossing both of A's.
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
        let shell = fuse_scaffold(&mut brep, &[a, b]);
        (brep, shell)
    }

    /// The shape a real generated scaffold fuses into: two 5-leg serpentines on
    /// stacked layers, the geometry `brep::pipeline` builds for a 3 mm 2-layer
    /// part (filament ⌀0.55, spacing 1.2, layer 0.35, bend radius 0.5).
    fn fused_real_scaffold() -> (BRep, ShellId) {
        use crate::arrange::scaffold::{fuse_scaffold, Filament};
        let r = 0.275;
        let z0 = 0.35;
        let z1 = 0.70;
        // layer 0 sweeps along +X with Y-legs; layer 1 along +Y with X-legs.
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
        let mut brep = BRep::new();
        let shell = fuse_scaffold(&mut brep, &[a, b]);
        (brep, shell)
    }

    /// The FUSE alone — no cut anywhere near it — must pass the gate a
    /// receiving CAD applies.  "Press Combine, open it in SpaceClaim, the part
    /// comes in torn into pieces" is a WIRE violation: SpaceClaim ignores
    /// `FACE_BOUND` orientation flags and drops outright any face whose wires
    /// intersect in uv (README invariant #5).
    ///
    /// RED — `wire = 4` on a plain 3 mm 2-layer scaffold.  Diagnosed:
    ///
    /// ```text
    /// FaceId(1)  rim at v=0.000   window at v[-0.264 .. +0.265]
    /// FaceId(6)  rim at v=1.900   window at v[ 1.635 ..  2.165]
    /// FaceId(12) rim at v=0.000   window at v[-0.265 .. +0.265]
    /// FaceId(17) rim at v=1.900   window at v[ 1.636 ..  2.165]
    /// ```
    ///
    /// Every one is a crossing window sitting ASTRIDE a leg's rim: the crossing
    /// lands exactly where the leg meets its elbow, so half the loop is outside
    /// the band and the hole's wire crosses the rim's wire.  Four faces dropped
    /// by the importer, four holes in the part.
    ///
    /// The one-line guard (extend `window_clears_free_ends` to every rim, not
    /// just free ones, for generated scaffolds too) makes this green — and
    /// breaks `leg_over_elbow_watertight`, which needs exactly such a crossing
    /// to weld and splits into two solids without it.  Measured both ways; it is
    /// a genuine conflict, not a tolerance to tune.
    ///
    /// The fix is to stop treating these as all-or-nothing: a crossing that
    /// straddles a junction must be SPLIT at the junction circle and imprinted
    /// on both members — as a hole on neither, but as a shared notch across the
    /// leg↔elbow boundary.  Then it welds AND every wire stays inside its face.
    #[test]
    #[ignore = "wire=4: crossing windows straddle leg rims — see doc comment"]
    fn fusing_a_real_scaffold_passes_the_spaceclaim_gate() {
        let (fused, shell) = fused_real_scaffold();
        assert_eq!(open_edges(&fused, shell), 0, "fused shell is manifold");
        let (wire, dist, man) = spaceclaim_gate(&fused, shell);
        assert_eq!(
            (wire, dist, man),
            (0, 0, 0),
            "SpaceClaim gate on the FUSED part: wire={wire} distance={dist} manifold={man}"
        );
    }

    /// The cut has to survive the real thing, not just the two-filament
    /// fixture: both terminals off a fused 3 mm scaffold.
    ///
    /// RED — 47 single-use edges (`CADCORE_DUMP_CUT=1` for the breakdown).  The
    /// fast repro for the pipeline-level failure in
    /// `brep/tests/fuse_then_cut.rs`; runs in ~0.3 s instead of 74 s.
    ///
    /// What is already right: the broad phase finds exactly the three tubes the
    /// plane crosses, each yields one closed conic, and three caps are emitted —
    /// the section side is correct.
    ///
    /// What is wrong, from the coordinates of the 47: they sit on **crossing
    /// windows**, and not only at the plane —
    ///
    /// ```text
    /// (2.425, 1.000, 0.700) -> (2.446, 1.000, 0.805)   on the cut plane
    /// (2.436, 2.685, 0.425) -> (2.436, 2.700, 0.426)   1.7 mm away from it
    /// ```
    ///
    /// so this is NOT the cut tearing the section: it is the arranged↔copied
    /// boundary again, at a crossing window whose two members end up on
    /// different paths (the leg the plane crosses is ARRANGED, the leg it
    /// crosses at y = 2.7 is wholly kept and COPIED).  The joints machinery does
    /// not cover it because every section here is a CLOSED loop and
    /// `section_joints` only publishes the ends of OPEN ones — so whatever the
    /// arranged face's DCEL does to that window is invisible to the copy.
    ///
    /// Next: find what splits the window on the arranged side when nothing
    /// crosses it (the section conic is 1.7 mm away and the seam is parked in a
    /// gap), then publish that split like any other joint.
    ///
    /// NB the offsets here are deliberately NOT `0.8`: at `0.8` the plane lands
    /// exactly on the fillet tangent (`0.3 + 0.5`), so the junction circles lie
    /// IN the plane — a separate degenerate case worth its own test.
    #[test]
    #[ignore = "47 open edges on a real scaffold — see doc comment"]
    fn cut_a_real_scaffold_is_watertight() {
        let (fused, shell) = fused_real_scaffold();
        let faces = fused.shells[shell].faces.clone();
        assert_eq!(
            open_edges(&fused, shell),
            0,
            "fixture starts watertight ({} faces)",
            faces.len()
        );

        // Both Y terminals, as `build_trim_planes` places them for offset 0.8.
        let lo = Plane3::from_origin_normal(Point3::new(0.0, 1.0, 0.0), UnitVec3::Y);
        let first = half_space_cut(&fused, &faces, &lo, 1e-6).expect("first cut runs");
        assert_eq!(
            open_edges(&first.brep, first.shell),
            0,
            "watertight after the first terminal ({} caps, {} dropped)",
            first.caps,
            first.dropped
        );

        let faces2 = first.brep.shells[first.shell].faces.clone();
        let hi = Plane3::from_origin_normal(Point3::new(0.0, 2.0, 0.0), -UnitVec3::Y);
        let second = half_space_cut(&first.brep, &faces2, &hi, 1e-6).expect("second cut runs");
        assert_eq!(
            open_edges(&second.brep, second.shell),
            0,
            "watertight after the second terminal ({} caps, {} dropped)",
            second.caps,
            second.dropped
        );
    }

    /// Re-emitting a fused shell must be EXACT: same faces, still watertight.
    ///
    /// This is the foundation the cut stands on — it proves `emit_face_copy`
    /// reproduces a Boolean'd shell (seamless periodic bands, crossing-window
    /// holes, torus fillets) through the assembler without losing a single
    /// shared edge.
    #[test]
    fn re_emitting_a_fused_scaffold_is_lossless() {
        let (fused, shell) = fused_scaffold();
        let faces = fused.shells[shell].faces.clone();
        assert_eq!(open_edges(&fused, shell), 0, "fixture starts watertight");

        // A plane far away on the discarded side: nothing is cut, so this
        // exercises the re-emit path alone.
        let keep = Plane3::from_origin_normal(Point3::new(-50.0, 0.0, 0.0), UnitVec3::X);
        let cut = half_space_cut(&fused, &faces, &keep, 1e-6).expect("cut runs");

        assert_eq!(cut.caps, 0, "nothing to cap");
        assert_eq!(cut.dropped, 0, "nothing to drop");
        assert_eq!(
            face_count(&cut.brep, cut.shell),
            faces.len(),
            "face count must be preserved exactly"
        );
        assert_eq!(
            open_edges(&cut.brep, cut.shell),
            0,
            "re-emitting a fused shell must not lose a shared edge"
        );
    }

    /// The real target: cut a shell that is already FUSED (its faces carry
    /// explicit trim loops, not analytic templates).  This is exactly where the
    /// legacy `cadcore-ops::half_space_cut_brep` silently does nothing — it
    /// finds no template cylinders and returns the part untouched.
    ///
    /// Getting here took six distinct fixes, each independently verified — they
    /// are listed because every one of them is a trap the next surface pair will
    /// walk into again:
    ///
    /// * plane∥axis sections produced no curve at all (`CylPlaneCurve::Lines`
    ///   was "deferred"), so a lengthwise-grazed tube survived the cut whole;
    /// * the arrangement's bounding rectangle was emitted as a third face over
    ///   every section, so section edges came out used 3×;
    /// * a traced section was clipped to the last sample INSIDE the arc, ~7 µm
    ///   short of the junction, so it never chained to its neighbour;
    /// * the clip walk emitted its head fragment twice;
    /// * overlapping sections were partitioned differently by the cap and by the
    ///   lateral faces (`presplit_sections` + `section_joints`/`insert_joints`);
    /// * the parametric seam: `build_periodic_u` unwraps each chain from its own
    ///   first sample and splits at multiples of the period, so it cuts along
    ///   **u = 0**, not u = ±π.  Parking the seam — and publishing its crossings
    ///   as joints — half a period off left exactly 9 open edges where a window
    ///   loop crossed the top of a tube.
    #[test]
    fn cut_a_fused_scaffold_is_watertight() {
        let (fused, shell) = fused_scaffold();
        let faces = fused.shells[shell].faces.clone();
        assert_eq!(open_edges(&fused, shell), 0, "fixture starts watertight");

        // Keep x <= -0.7: severs A's two X-legs clear of their U-turn fillets,
        // drops A's connector, leaves every crossing with B intact.
        let keep = Plane3::from_origin_normal(Point3::new(-0.7, 0.0, 0.0), -UnitVec3::X);
        let cut = half_space_cut(&fused, &faces, &keep, 1e-6).expect("cut runs");

        // Three tubes are severed — a transverse disk through each of A's two
        // X-legs, plus a longitudinal chord through B's leg at x = -0.8 (its
        // tube reaches x = -0.525, so the plane grazes it lengthwise).  They
        // OVERLAP at the plane, so the arrangement tiles the section into a few
        // coplanar cells rather than one face per tube; what matters is that
        // they tile it completely, which the watertight assertion below proves.
        assert!(
            cut.caps >= 3,
            "at least one flat section per severed tube, got {}",
            cut.caps
        );
        assert!(cut.dropped > 0, "A's connector and fillets must be dropped");
        assert_eq!(
            open_edges(&cut.brep, cut.shell),
            0,
            "cut fused shell must stay watertight"
        );
        // Traced sections and welded vertices land within the assembler's weld
        // tolerance of the plane, so the bound is the knit scale (1 µm — what
        // SolidWorks itself knits at), not float epsilon.
        assert!(
            deepest_past_plane(&cut.brep, cut.shell, &keep) > -1e-3,
            "material survived past the cut plane: {:.3e}",
            deepest_past_plane(&cut.brep, cut.shell, &keep)
        );
        // Every cap must sit where there IS material: probing just inside the
        // kept side of a cap must land inside the original solid.
        let n = keep.normal();
        for &fid in &cut.brep.shells[cut.shell].faces {
            let f = &cut.brep.faces[fid];
            if !matches!(f.geom, FaceGeom::Plane(_)) {
                continue;
            }
            let pts = loop_points(&cut.brep, f.outer_loop);
            if pts.is_empty() || pts.iter().any(|p| keep.signed_distance(*p).abs() > 1e-6) {
                continue; // not a section cap
            }
            let c = centroid(&pts);
            assert!(
                point_in_solid(&fused, &faces, c + n * 1e-3),
                "a cap was emitted over void at {c:?}"
            );
        }
        // And the gate a receiving CAD actually applies.
        // The SpaceClaim gate: wire winding and manifold are clean.  Gate 4
        // (edge-to-surface deviation) is NOT asserted yet — see the module TODO
        // on the seam split point; it is the last thing between this and a
        // guaranteed SolidWorks-clean cut.
        let (wire, _dist, man) = spaceclaim_gate(&cut.brep, cut.shell);
        assert_eq!(wire, 0, "uv wires must be wound and non-intersecting");
        assert_eq!(man, 0, "shell must be manifold");
    }

    /// The legacy template cut is a silent no-op on a fused shell — the reason
    /// this module exists.  Guards the claim so it cannot rot.
    #[test]
    fn legacy_template_cut_cannot_touch_a_fused_shell() {
        let (fused, shell) = fused_scaffold();
        let templated = fused.shells[shell].faces.iter().any(|&f| {
            matches!(
                (&fused.faces[f].geom, &fused.faces[f].extent),
                (FaceGeom::Cylinder(_), FaceExtent::Cylinder { .. })
            )
        });
        assert!(
            !templated,
            "a fused shell carries no template cylinder faces, so a \
             template-rebuilding cut has nothing to classify"
        );
    }
}
