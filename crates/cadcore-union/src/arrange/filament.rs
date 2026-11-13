//! The filament model and the closed shell of ONE swept filament.
//!
//! A [`Filament`] is a printed path as the kernel sees it: straight
//! [`LegSeg`] cylinders joined by constant-radius [`ElbowSeg`] torus elbows,
//! with a flat cap on each free end.  [`Filament::serpentine`] builds one from
//! a polyline.  Lists of filaments are what the lattice Booleans consume
//! ([`crate::fuse_filaments_with`], [`crate::fuse_and_trim`],
//! [`crate::LatticeModel`]).
//!
//! [`fuse_poly_filament`] / [`fuse_l_filament`] assemble a single filament as
//! a closed analytic solid on its own (full cylinder and torus faces) — the
//! input the generic solid union ([`crate::boolean::union_n`]) works on.

use cadcore_geom::{Circle3, CylSurf, Plane3, TorusSurf};
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{BRep, Shell, ShellId, Solid};

use super::assembly::Assembler;



/// A leg: finite cylinder + its axial span `[v0, v1]` (the band runs
/// `v0..v1`; junctions/caps sit at the ends).
#[derive(Debug, Clone, Copy)]
pub struct LegSeg {
    /// Carrier cylinder (origin = axial 0).
    pub surf: CylSurf,
    /// Axial length.
    pub length: f64,
}

/// An elbow torus fillet over the ring arc `[theta_lo, theta_hi]`.
#[derive(Debug, Clone, Copy)]
pub struct ElbowSeg {
    /// Carrier torus.
    pub surf: TorusSurf,
    /// Ring-arc start.
    pub theta_lo: f64,
    /// Ring-arc end.
    pub theta_hi: f64,
}

/// One serpentine filament path: `legs[i]` joined to `legs[i+1]` by
/// `elbows[i]`.  The two free ends (first leg incoming, last leg outgoing)
/// are capped.  Each leg is built with `surf.frame.origin` at its incoming
/// end and `surf.axis()` pointing toward the outgoing end.
#[derive(Debug, Clone)]
pub struct Filament {
    /// Straight legs in path order.
    pub legs: Vec<LegSeg>,
    /// Elbows (`legs.len() - 1` of them); `elbows[i]` joins legs i and i+1.
    pub elbows: Vec<ElbowSeg>,
}

impl Filament {
    /// Build a planar serpentine filament from a polyline path: straight legs
    /// joined by constant-radius circular fillets (one per interior vertex).
    /// `r` = filament radius, `fillet` = elbow bend radius.  Each leg runs
    /// `surf.frame.origin` → `+axis`, the convention every fuse in this crate expects.
    pub fn serpentine(verts: &[Point3], r: f64, fillet: f64) -> Filament {
        let fillets = vec![fillet; verts.len().saturating_sub(2)];
        Self::serpentine_with_fillets(verts, r, &fillets)
    }

    /// [`Filament::serpentine`] with a bend radius per interior vertex
    /// (`fillets[k - 1]` for `verts[k]`).  Walls of different layers stacked
    /// on one another must END at the same points along the wall: the tangent
    /// length `R·tan(turn/2)` depends on the turn, and turns that differ by
    /// 0.02° (corner noise of a toolpath) put the two junction planes 80 nm
    /// apart — an edge that short in every crossing.  The caller evens that
    /// out by nudging the radius of one elbow (by ~0.05 %).
    pub fn serpentine_with_fillets(verts: &[Point3], r: f64, fillets: &[f64]) -> Filament {
        let n = verts.len();
        assert!(n >= 2, "a filament needs at least two path vertices");
        assert_eq!(fillets.len(), n - 2, "one bend radius per interior vertex");
        let dir = |a: Point3, b: Point3| UnitVec3::try_from_vec(b - a).unwrap();
        let mut tin = vec![Point3::new(0.0, 0.0, 0.0); n];
        let mut tout = vec![Point3::new(0.0, 0.0, 0.0); n];
        let mut elbows = Vec::new();
        for k in 1..n.saturating_sub(1) {
            let din = dir(verts[k - 1], verts[k]);
            let dout = dir(verts[k], verts[k + 1]);
            let turn = din.dot_vec(dout.as_vec()).clamp(-1.0, 1.0).acos();
            // Tangent length of a fillet of radius `fillet` across a corner
            // deflecting by `turn`: R·tan(turn/2).  (R/tan(turn/2) is only
            // equal at 90°; on the slanted connectors a real generator emits
            // it put the leg start 0.07 mm past the elbow's end, leaving the
            // junction circles of leg and elbow on different circles.)
            let fillet = fillets[k - 1];
            let t = fillet * (turn / 2.0).tan();
            let pin = verts[k] - din.as_vec() * t;
            let pout = verts[k] + dout.as_vec() * t;
            let axis = UnitVec3::try_from_vec(din.as_vec().cross(dout.as_vec())).unwrap();
            let inward = UnitVec3::try_from_vec(axis.as_vec().cross(din.as_vec())).unwrap();
            let centre = pin + inward.as_vec() * fillet;
            let torus = TorusSurf::new(centre, axis, fillet, r);
            let theta = |p: Point3| {
                let w = p - torus.frame.origin;
                let rd = w - torus.frame.z.as_vec() * torus.frame.z.dot_vec(w);
                torus.frame.y.dot_vec(rd).atan2(torus.frame.x.dot_vec(rd))
            };
            let lo = theta(pin);
            let mut d = theta(pout) - lo;
            while d > std::f64::consts::PI { d -= std::f64::consts::TAU; }
            while d < -std::f64::consts::PI { d += std::f64::consts::TAU; }
            elbows.push(ElbowSeg { surf: torus, theta_lo: lo, theta_hi: lo + d });
            tin[k] = pin;
            tout[k] = pout;
        }
        let mut legs = Vec::new();
        for k in 0..n - 1 {
            let start = if k == 0 { verts[0] } else { tout[k] };
            let end = if k + 1 == n - 1 { verts[n - 1] } else { tin[k + 1] };
            legs.push(LegSeg {
                surf: CylSurf::new(start, dir(start, end), r),
                length: (end - start).length(),
            });
        }
        Filament { legs, elbows }
    }

}






/// Build a closed shell for an L-filament: `leg_a` → `elbow` → `leg_b`,
/// capped at the two free ends.  `j0` is the leg_a↔elbow junction axial v on
/// leg_a (its junction end); `j1` likewise on leg_b.  Returns the shell.
#[allow(clippy::too_many_arguments)]
pub fn fuse_l_filament(
    brep: &mut BRep,
    leg_a: LegSeg,
    elbow: ElbowSeg,
    leg_b: LegSeg,
    cap_a_at_start: bool,
    cap_b_at_start: bool,
) -> ShellId {
    // Delegate to the general seamless poly-filament builder.  The L-test
    // builds leg_a origin=free(incoming)→junction and leg_b
    // origin=junction(incoming)→free, i.e. both legs origin=incoming, which
    // is exactly the poly convention; the cap flags are implied by the free
    // ends (first leg incoming, last leg outgoing).
    let _ = (cap_a_at_start, cap_b_at_start);
    fuse_poly_filament(brep, &[leg_a, leg_b], &[elbow])
}
/// Build a closed shell for a polyline filament: `legs[0]` → `elbows[0]` →
/// `legs[1]` → … → `legs[n-1]`, capped at the two free ends.  `elbows[i]`
/// connects `legs[i]` and `legs[i+1]`; its θ-range must be oriented so
/// θ_lo touches leg i and θ_hi touches leg i+1.
///
/// Junction circles are sampled ONCE from each elbow and shared by both the
/// elbow patch and the two mating leg bands; each leg's seam is aligned to
/// its junction sample so all three faces split shared circles identically.
pub fn fuse_poly_filament(
    brep: &mut BRep,
    legs: &[LegSeg],
    elbows: &[ElbowSeg],
) -> ShellId {
    assert_eq!(legs.len(), elbows.len() + 1, "n legs <-> n-1 elbows");
    let n = legs.len();

    // SEAMLESS analytic junction / end circles, built ONCE from leg geometry
    // and shared (shared_circle) between each leg rim and the mating elbow /
    // cap.  j_in[i] = leg i incoming-end circle, j_out[i] = outgoing-end.
    let circ_at = |l: &LegSeg, v: f64| -> Circle3 {
        Circle3::new(l.surf.frame.origin + l.surf.axis().as_vec() * v, l.surf.axis(), l.surf.radius)
    };
    let j_in: Vec<Circle3> = legs.iter().map(|l| circ_at(l, 0.0)).collect();
    let j_out: Vec<Circle3> = legs.iter().map(|l| circ_at(l, l.length)).collect();

    let mut asm = Assembler::new(brep, 1e-7);

    // Per-junction weld groups (see scaffold.rs): leg i owns 2·i / 2·i+1.
    let gin = |i: usize| 2 * i as u32;
    let gout = |i: usize| 2 * i as u32 + 1;
    // Legs: seamless cylinder faces (no windows here), rims = the shared
    // junction / end circles.
    for i in 0..n {
        asm.emit_cylinder_face(legs[i].surf, legs[i].length, j_in[i], j_out[i], &[], gin(i), gout(i));
        if i == 0 {
            let n0 = UnitVec3::try_from_vec(legs[i].surf.axis().as_vec() * -1.0).unwrap();
            let plane = Plane3::from_origin_normal(legs[i].surf.frame.origin, n0);
            asm.emit_disk_cap(plane, j_in[i], gin(i));
        }
        if i + 1 == n {
            let c = legs[i].surf.frame.origin + legs[i].surf.axis().as_vec() * legs[i].length;
            let plane = Plane3::from_origin_normal(c, legs[i].surf.axis());
            asm.emit_disk_cap(plane, j_out[i], gout(i));
        }
    }
    // Elbows: seamless torus-fillet faces; junction_lo = leg i outgoing circle,
    // junction_hi = leg i+1 incoming circle (both shared with the legs).
    for (i, e) in elbows.iter().enumerate() {
        asm.emit_torus_fillet(e.surf, j_out[i], j_in[i + 1], gout(i), gin(i + 1));
    }

    let faces = asm.faces().to_vec();
    let shell = brep.add_shell(Shell {
        faces: faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = brep.add_solid(Solid {
        shells: vec![shell],
        name: Some("filament".into()),
    });
    brep.shells[shell].solid = solid;
    for f in &faces {
        brep.faces[*f].shell = shell;
    }
    shell
}




#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::UnionConfig;
    use crate::validate::validate_shell;
    use crate::Diagnostics;
    use std::collections::HashMap;

    /// θ of a junction point in the torus's OWN frame (never assume
    /// perp_basis axes — derive from the actual frame).
    fn theta_of(t: &TorusSurf, p: Point3) -> f64 {
        let w = p - t.frame.origin;
        let rd = w - t.frame.z.as_vec() * t.frame.z.dot_vec(w);
        t.frame.y.dot_vec(rd).atan2(t.frame.x.dot_vec(rd))
    }

    /// L-filament: leg +X → 90° elbow → leg +Y, capped both ends.
    fn l_filament() -> (BRep, cadcore_topo::ShellId) {
        let r = 0.275;
        // elbow torus: centre (0,0.5,0), axis +Z, major 0.5
        let t = TorusSurf::new(Point3::new(0.0, 0.5, 0.0), UnitVec3::Z, 0.5, r);
        // junctions derived from the actual leg ends, in the torus frame:
        // leg A ends at (0,0,0); leg B starts at (0.5,0.5,0)
        let th_a = theta_of(&t, Point3::new(0.0, 0.0, 0.0));
        let th_b = theta_of(&t, Point3::new(0.5, 0.5, 0.0));
        let elbow = ElbowSeg {
            surf: t,
            theta_lo: th_a,
            theta_hi: th_b,
        };
        let leg_a = LegSeg {
            surf: CylSurf::new(Point3::new(-2.0, 0.0, 0.0), UnitVec3::X, r),
            length: 2.0, // ends at (0,0,0)
        };
        let leg_b = LegSeg {
            surf: CylSurf::new(Point3::new(0.5, 0.5, 0.0), UnitVec3::Y, r),
            length: 2.0,
        };
        let mut brep = BRep::new();
        // leg A junction end is at v=length (x=0); cap at v=0 (start).
        // leg B junction end is at v=0; cap at v=length (end).
        let shell = fuse_l_filament(&mut brep, leg_a, elbow, leg_b, true, false);
        (brep, shell)
    }

    fn edge_uses(brep: &BRep, shell: cadcore_topo::ShellId) -> HashMap<cadcore_topo::EdgeId, usize> {
        let mut uses = HashMap::new();
        for &fid in &brep.shells[shell].faces {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let st = brep.loops[lid].start;
                let mut c = st;
                loop {
                    let ce = &brep.coedges[c];
                    *uses.entry(ce.edge).or_insert(0) += 1;
                    c = ce.next;
                    if c == st {
                        break;
                    }
                }
            }
        }
        uses
    }

    /// Build legs+elbows from a polyline path with 90°-ish corners filleted
    /// at radius `fr`, tube radius `r`.  Each leg origin = its incoming end,
    /// axis toward outgoing; each elbow θ oriented leg i → leg i+1.
    fn filament_from_path(verts: &[Point3], r: f64, fr: f64) -> (Vec<LegSeg>, Vec<ElbowSeg>) {
        let n = verts.len();
        assert!(n >= 2);
        // tangent (trim) length at each interior corner + arc data
        let mut t_in = vec![Point3::new(0.0, 0.0, 0.0); n]; // leg-end before corner
        let mut t_out = vec![Point3::new(0.0, 0.0, 0.0); n]; // leg-start after corner
        let mut elbows: Vec<ElbowSeg> = Vec::new();
        let dir = |a: Point3, b: Point3| UnitVec3::try_from_vec(b - a).unwrap();
        for k in 1..n - 1 {
            let din = dir(verts[k - 1], verts[k]); // into corner
            let dout = dir(verts[k], verts[k + 1]); // out of corner
            let turn = din.dot_vec(dout.as_vec()).clamp(-1.0, 1.0).acos();
            let t = fr / (turn / 2.0).tan(); // trim length (90° → fr)
            let ti = verts[k] - din.as_vec() * t;
            let to = verts[k] + dout.as_vec() * t;
            // arc centre: from ti, perpendicular to din toward dout side
            let axis = UnitVec3::try_from_vec(din.as_vec().cross(dout.as_vec())).unwrap();
            let inn = UnitVec3::try_from_vec(axis.as_vec().cross(din.as_vec())).unwrap();
            let centre = ti + inn.as_vec() * fr;
            let torus = TorusSurf::new(centre, axis, fr, r);
            let th_lo = {
                let w = ti - torus.frame.origin;
                let rd = w - torus.frame.z.as_vec() * torus.frame.z.dot_vec(w);
                torus.frame.y.dot_vec(rd).atan2(torus.frame.x.dot_vec(rd))
            };
            let th_hi = {
                let w = to - torus.frame.origin;
                let rd = w - torus.frame.z.as_vec() * torus.frame.z.dot_vec(w);
                torus.frame.y.dot_vec(rd).atan2(torus.frame.x.dot_vec(rd))
            };
            // shortest signed sweep
            let mut d = th_hi - th_lo;
            while d > std::f64::consts::PI { d -= std::f64::consts::TAU; }
            while d < -std::f64::consts::PI { d += std::f64::consts::TAU; }
            elbows.push(ElbowSeg { surf: torus, theta_lo: th_lo, theta_hi: th_lo + d });
            t_in[k] = ti;
            t_out[k] = to;
        }
        // legs: leg k runs from (corner k start) to (corner k+1 end)
        let mut legs: Vec<LegSeg> = Vec::new();
        for k in 0..n - 1 {
            let start = if k == 0 { verts[0] } else { t_out[k] };
            let end = if k + 1 == n - 1 { verts[n - 1] } else { t_in[k + 1] };
            let axis = dir(start, end);
            let length = (end - start).length();
            legs.push(LegSeg {
                surf: CylSurf::new(start, axis, r),
                length,
            });
        }
        (legs, elbows)
    }

    #[test]
    fn u_filament_is_watertight() {
        // planar U: +X, 90° up, +Y, 90° back, -X
        let verts = [
            Point3::new(-2.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(-2.0, 2.0, 0.0),
        ];
        let (legs, elbows) = filament_from_path(&verts, 0.275, 0.5);
        assert_eq!(legs.len(), 3);
        assert_eq!(elbows.len(), 2);
        let mut brep = BRep::new();
        let shell = fuse_poly_filament(&mut brep, &legs, &elbows);
        // 3 legs + 2 elbows + 2 caps
        assert_eq!(brep.shells[shell].faces.len(), 7, "3 legs + 2 elbows + 2 caps");
        let uses = edge_uses(&brep, shell);
        let open = uses.values().filter(|&&n| n == 1).count();
        assert_eq!(open, 0, "U-filament watertight, got {open} open");
        assert!(uses.values().all(|&n| n == 2));
    }

    #[test]
    fn u_filament_passes_gates() {
        let verts = [
            Point3::new(-2.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(-2.0, 2.0, 0.0),
        ];
        let (legs, elbows) = filament_from_path(&verts, 0.275, 0.5);
        let mut brep = BRep::new();
        let shell = fuse_poly_filament(&mut brep, &legs, &elbows);
        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let report = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(report.manifold_violations.is_empty(), "manifold {:?}",
            &report.manifold_violations[..report.manifold_violations.len().min(4)]);
        assert!(report.distance_violations.is_empty(),
            "distance max {:.2}µm", report.max_edge_deviation * 1e3);
    }

    #[test]
    fn l_filament_is_watertight() {
        let (brep, shell) = l_filament();
        // 2 leg laterals + 1 elbow patch + 2 caps
        assert_eq!(brep.shells[shell].faces.len(), 5, "2 legs + elbow + 2 caps");
        let uses = edge_uses(&brep, shell);
        let open = uses.values().filter(|&&n| n == 1).count();
        assert_eq!(open, 0, "watertight: no open edges, got {open}");
        assert!(uses.values().all(|&n| n == 2), "every edge shared exactly twice");
    }

    #[test]
    fn l_filament_passes_gates() {
        let (brep, shell) = l_filament();
        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let report = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(
            report.manifold_violations.is_empty(),
            "manifold: {:?}",
            &report.manifold_violations[..report.manifold_violations.len().min(4)]
        );
        assert!(
            report.distance_violations.is_empty(),
            "distance (max {:.2} µm)",
            report.max_edge_deviation * 1e3
        );
    }
}
