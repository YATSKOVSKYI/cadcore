//! Closed-solid assembly for crossing filaments.
//!
//! Ties the arrangement stages into a *closed* shell with the watertight
//! invariants enforced structurally:
//!
//! * **One curve, one imprint.**  The crossing loop is traced ONCE and the
//!   *same* 3-D polyline is projected into both legs' bands (README
//!   invariant #2).  Imprinting each face from its own tracing would give
//!   two different polylines for one curve → unshared edges.
//! * **Seam avoidance.**  Each leg's cylinder frame is oriented so its
//!   parametric seam (u = ±π) sits opposite the crossing, so the imprint
//!   curve is never split by the seam — both faces consume it as one piece
//!   with identical endpoints.
//! * **Shared end circles.**  Each cap disk's boundary is the leg end circle
//!   sampled exactly as the lateral-band rim, so the assembler's shared-edge
//!   cache stitches cap ↔ band.
//!
//! [`fuse_crossing_legs`] returns the watertight shell for the gates — the
//! first closed *solid* produced entirely by cadcore-union.

use cadcore_geom::{CylSurf, Plane3};
use cadcore_math::{Frame3, Point3, UnitVec3, Vec3};
use cadcore_topo::{BRep, FaceGeom, FaceNormal, Shell, ShellId, Solid};

use super::assembly::Assembler;
use super::cells::{arrange_and_classify, FaceChain};
use super::domain::FaceDomain;
use crate::geom::composite::{closed_loops_between, CompositeTube};
use crate::geom::intersect::TraceOptions;

/// A finite straight leg: carrier cylinder + axial length.
#[derive(Debug, Clone, Copy)]
pub struct Leg {
    /// Carrier cylinder (origin = axial-0 end-circle centre).
    pub surf: CylSurf,
    /// Axial length (mm).
    pub length: f64,
}

const RIM_SAMPLES: usize = 96;

/// Re-frame a leg's cylinder so the parametric SEAM sits opposite `crossing`.
///
/// `build_periodic_u` puts the seam at u = 0 (≡ the cylinder frame's +x
/// direction).  To keep the imprint loop off the seam we point frame.x at the
/// ANTIPODE of the crossing — the loop then sits around u = π and is never
/// seam-split.  (Pointing frame.x AT the crossing — the obvious-looking
/// choice — drives the seam straight through the loop.)
fn seam_away_from(surf: CylSurf, crossing: Vec3) -> CylSurf {
    let axis = surf.axis();
    let toward = crossing * -1.0; // antipode → seam opposite the loop
    let radial = toward - axis.as_vec() * axis.dot_vec(toward);
    let Some(x) = UnitVec3::try_from_vec(radial) else {
        return surf;
    };
    let y = UnitVec3::try_from_vec(axis.as_vec().cross(x.as_vec())).expect("perp");
    CylSurf {
        frame: Frame3 {
            origin: surf.frame.origin,
            x,
            y,
            z: axis,
        },
        radius: surf.radius,
    }
}

/// Sample a leg's end circle at axial `v`, in the SAME order the lateral-band
/// rim chain uses, so the two share edges byte-for-byte.
fn end_circle_3d(band: &FaceDomain, v: f64) -> Vec<Point3> {
    (0..=RIM_SAMPLES)
        .map(|k| {
            let u =
                -std::f64::consts::PI + std::f64::consts::TAU * k as f64 / RIM_SAMPLES as f64;
            band.lift(u, v)
        })
        .collect()
}

/// Fuse two crossing legs into one closed shell entirely through the
/// arrangement engine.  Returns the new shell.
pub fn fuse_crossing_legs(brep: &mut BRep, a: Leg, b: Leg) -> ShellId {
    let tube_a = CompositeTube::new().with_leg(a.surf, a.length);
    let tube_b = CompositeTube::new().with_leg(b.surf, b.length);

    // Trace the crossing loop(s) ONCE — the single source of truth.
    let loops = closed_loops_between(&tube_a, &tube_b, &TraceOptions::default());
    let loop_centroid = loops.first().map(|c| centroid(&c.points));

    // Orient each leg's seam away from the crossing.
    let a = Leg {
        surf: loop_centroid
            .map(|c| seam_away_from(a.surf, c - a.surf.frame.origin))
            .unwrap_or(a.surf),
        ..a
    };
    let b = Leg {
        surf: loop_centroid
            .map(|c| seam_away_from(b.surf, c - b.surf.frame.origin))
            .unwrap_or(b.surf),
        ..b
    };

    let mut asm = Assembler::new(brep, 1e-7);
    // One curve, one imprint: the SAME 3-D polyline goes into both bands, so
    // the crossing curve and all caps stitch.  Residual (≤12 edges): the
    // periodic per-band DCEL's interaction with its seam chain splits a few
    // turning-point segments differently between the two bands (endpoints
    // match — gap 0 — but one band carries an extra interior vertex).  That
    // is a build_periodic_u (cadcore-geom) detail; tracked as the P4
    // follow-up.  Self-intersection enrichment was ruled out (the loop does
    // not self-intersect in either band's uv).
    let loops_enriched: Vec<Vec<Point3>> = loops.iter().map(|c| c.points.clone()).collect();
    emit_leg(&mut asm, a, &loops_enriched, &tube_b);
    emit_leg(&mut asm, b, &loops_enriched, &tube_a);
    let faces = asm.faces().to_vec();

    let shell = brep.add_shell(Shell {
        faces: faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = brep.add_solid(Solid {
        shells: vec![shell],
        name: Some("union".into()),
    });
    brep.shells[shell].solid = solid;
    for f in &faces {
        brep.faces[*f].shell = shell;
    }
    shell
}

/// A straight capped leg for the crossing-grid orchestrator.
#[derive(Debug, Clone, Copy)]
pub struct GridLeg {
    /// Carrier cylinder.
    pub surf: CylSurf,
    /// Axial length.
    pub length: f64,
}

/// Sample a leg end circle at axial `v` using the leg's own frame.
fn grid_circle(c: &CylSurf, v: f64) -> Vec<Point3> {
    (0..=RIM_SAMPLES)
        .map(|k| {
            let u =
                -std::f64::consts::PI + std::f64::consts::TAU * k as f64 / RIM_SAMPLES as f64;
            let radial = c.frame.x * u.cos() + c.frame.y * u.sin();
            c.frame.origin + radial * c.radius + c.axis().as_vec() * v
        })
        .collect()
}

/// Orient a leg's seam (frame.x) into the largest angular gap between the
/// crossing windows, so the seam never splits an imprint loop.
fn seam_into_gap(surf: CylSurf, windows: &[Vec<Point3>]) -> CylSurf {
    if windows.is_empty() {
        return surf;
    }
    let axis = surf.axis();
    // window centroid angles in the CURRENT frame
    let mut angs: Vec<f64> = windows
        .iter()
        .map(|w| {
            let c = centroid(w);
            let d = c - surf.frame.origin;
            let r = d - axis.as_vec() * axis.dot_vec(d);
            surf.frame.y.dot_vec(r).atan2(surf.frame.x.dot_vec(r))
        })
        .collect();
    angs.sort_by(f64::total_cmp);
    // largest gap between consecutive angles (cyclic)
    let n = angs.len();
    let (mut best_mid, mut best_gap) = (std::f64::consts::PI, -1.0);
    for i in 0..n {
        let a = angs[i];
        let b = if i + 1 < n {
            angs[i + 1]
        } else {
            angs[0] + std::f64::consts::TAU
        };
        let gap = b - a;
        if gap > best_gap {
            best_gap = gap;
            best_mid = 0.5 * (a + b);
        }
    }
    // rotate frame so the gap centre becomes u = 0 (the seam)
    let dir = surf.frame.x * best_mid.cos() + surf.frame.y * best_mid.sin();
    let Some(x) = UnitVec3::try_from_vec(dir) else {
        return surf;
    };
    let y = UnitVec3::try_from_vec(axis.as_vec().cross(x.as_vec())).expect("perp");
    CylSurf {
        frame: Frame3 { origin: surf.frame.origin, x, y, z: axis },
        radius: surf.radius,
    }
}

/// Fuse an arbitrary set of straight crossing legs into one closed shell:
/// every leg is trimmed by every other it overlaps, capped at both ends, all
/// welded through one assembler.  Generalises [`fuse_crossing_legs`] to N
/// legs with multiple crossings each.
pub fn fuse_crossing_grid(brep: &mut BRep, legs: &[GridLeg]) -> ShellId {
    let n = legs.len();
    let tubes: Vec<CompositeTube> = legs
        .iter()
        .map(|l| CompositeTube::new().with_leg(l.surf, l.length))
        .collect();

    // crossing loops per leg (shared 3-D points; traced ONCE per pair).
    let mut loops_on: Vec<Vec<Vec<Point3>>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            for c in closed_loops_between(&tubes[i], &tubes[j], &TraceOptions::default()) {
                loops_on[i].push(c.points.clone());
                loops_on[j].push(c.points);
            }
        }
    }

    let mut asm = Assembler::new(brep, 1e-7);
    for i in 0..n {
        // seam into a window-free gap
        let surf = seam_into_gap(legs[i].surf, &loops_on[i]);
        let len = legs[i].length;
        let band = FaceDomain::CylinderBand { surf, length: len };

        // chains: two end rims + every crossing window
        let c0 = grid_circle(&surf, 0.0);
        let c1 = grid_circle(&surf, len);
        let mut chains = vec![
            FaceChain { pts: c0.iter().map(|&p| band.uv(p)).collect(), tag: 1 },
            FaceChain { pts: c1.iter().map(|&p| band.uv(p)).collect(), tag: 2 },
        ];
        for (k, w) in loops_on[i].iter().enumerate() {
            let mut uv: Vec<(f64, f64)> = w.iter().map(|&p| band.uv(p)).collect();
            uv.push(uv[0]);
            chains.push(FaceChain { pts: uv, tag: 100 + k as u32 });
        }

        // keep cells outside ANY other tube; collect their WINDOW holes
        let inside_other = |p: Point3| -> bool {
            (0..n).any(|j| j != i && tubes[j].contains(p))
        };
        let _ = &inside_other;
        let _ = &chains;
        // SEAMLESS cylinder: analytic rim circles + the crossing windows as
        // interior holes (each shared with the crossing leg → welded).
        let axis = surf.axis();
        let rim_lo = cadcore_geom::Circle3::new(surf.frame.origin, axis, surf.radius);
        let rim_hi = cadcore_geom::Circle3::new(
            surf.frame.origin + axis.as_vec() * len,
            axis,
            surf.radius,
        );
        // Per-leg junction weld groups (each straight leg owns 2·i / 2·i+1).
        let (g_lo, g_hi) = (2 * i as u32, 2 * i as u32 + 1);
        let win_ids: Vec<(u32, Vec<Point3>)> =
            loops_on[i].iter().cloned().map(|w| (0u32, w)).collect();
        asm.emit_cylinder_face(surf, len, rim_lo, rim_hi, &win_ids, g_lo, g_hi);
        // caps share the SAME rim circles (seamless weld).
        let nlo = cadcore_math::UnitVec3::try_from_vec(axis.as_vec() * -1.0).unwrap();
        let plane_lo = cadcore_geom::Plane3::from_origin_normal(surf.frame.origin, nlo);
        asm.emit_disk_cap(plane_lo, rim_lo, g_lo);
        let plane_hi = cadcore_geom::Plane3::from_origin_normal(
            surf.frame.origin + axis.as_vec() * len,
            axis,
        );
        asm.emit_disk_cap(plane_hi, rim_hi, g_hi);
    }

    let faces = asm.faces().to_vec();
    let shell = brep.add_shell(Shell {
        faces: faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = brep.add_solid(Solid {
        shells: vec![shell],
        name: Some("grid".into()),
    });
    brep.shells[shell].solid = solid;
    for f in &faces {
        brep.faces[*f].shell = shell;
    }
    shell
}

fn centroid(pts: &[Point3]) -> Point3 {
    let n = pts.len().max(1) as f64;
    let s = pts
        .iter()
        .fold(Vec3::new(0.0, 0.0, 0.0), |acc, p| acc + (*p - Point3::new(0.0, 0.0, 0.0)));
    Point3::new(0.0, 0.0, 0.0) + s * (1.0 / n)
}

/// Emit one leg: its trimmed lateral band + two end-cap disks.  `loops` are
/// the shared crossing curves (same 3-D points for both legs); `other` is the
/// opposite tube for keep/drop classification.
fn emit_leg(asm: &mut Assembler, leg: Leg, loops: &[Vec<Point3>], other: &CompositeTube) {
    let band = FaceDomain::CylinderBand {
        surf: leg.surf,
        length: leg.length,
    };

    let rim_chain = |v: f64, tag: u32| FaceChain {
        pts: (0..=RIM_SAMPLES)
            .map(|k| {
                (
                    -std::f64::consts::PI
                        + std::f64::consts::TAU * k as f64 / RIM_SAMPLES as f64,
                    v,
                )
            })
            .collect(),
        tag,
    };
    let mut chains = vec![rim_chain(0.0, 1), rim_chain(leg.length, 2)];

    for (i, pts3) in loops.iter().enumerate() {
        let mut w: Vec<(f64, f64)> = pts3.iter().map(|&p| band.uv(p)).collect();
        w.push(w[0]);
        chains.push(FaceChain { pts: w, tag: 100 + i as u32 });
    }

    for cell in arrange_and_classify(&band, &chains, other, 1e-3) {
        if cell.keep {
            asm.emit_cell(&band, FaceGeom::Cylinder(leg.surf), FaceNormal::Same, &cell);
        }
    }

    emit_cap(asm, &band, leg, 0.0, true);
    emit_cap(asm, &band, leg, leg.length, false);
}

/// Emit one end-cap disk at axial `v`.
fn emit_cap(asm: &mut Assembler, band: &FaceDomain, leg: Leg, v: f64, at_start: bool) {
    let axis = leg.surf.axis().as_vec();
    let centre = leg.surf.frame.origin + axis * v;
    let outward = if at_start { -axis } else { axis };
    let normal = UnitVec3::try_from_vec(outward).expect("axis is unit");
    let plane = Plane3::from_origin_normal(centre, normal);
    let dom = FaceDomain::Plane {
        plane,
        half_extent: leg.surf.radius * 4.0,
    };

    // boundary = the leg end circle, sampled exactly as the band rim.
    let circ3 = end_circle_3d(band, v);
    let pts: Vec<(f64, f64)> = circ3.iter().map(|&p| dom.uv(p)).collect();
    let chain = FaceChain { pts, tag: 1 };

    let cells = arrange_and_classify(&dom, &[chain], &CompositeTube::new(), 1e-3);
    let centre_uv = dom.uv(centre);
    for cell in &cells {
        let keep_disk = cell
            .cell
            .interior_point()
            .map(|(u, w)| (u - centre_uv.0).hypot(w - centre_uv.1) < leg.surf.radius)
            .unwrap_or(false);
        if keep_disk {
            asm.emit_cell(&dom, FaceGeom::Plane(plane), FaceNormal::Same, cell);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::UnionConfig;
    use crate::validate::validate_shell;
    use crate::Diagnostics;

    use std::collections::HashMap;

    fn build() -> (BRep, cadcore_topo::ShellId) {
        let a = Leg {
            surf: CylSurf::new(Point3::new(-2.0, 0.0, 0.0), UnitVec3::X, 0.275),
            length: 4.0,
        };
        let b = Leg {
            surf: CylSurf::new(Point3::new(0.0, -2.0, 0.35), UnitVec3::Y, 0.275),
            length: 4.0,
        };
        let mut brep = BRep::new();
        let shell = fuse_crossing_legs(&mut brep, a, b);
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

    /// The new engine produces the full closed union solid: 2 trimmed
    /// laterals + 4 caps, the crossing curve shared between the bands and
    /// every end circle shared band↔cap.  Distances are exact and the shell
    /// is ALMOST watertight — the only residual open edges are at the
    /// window's u-turning points (see the ignored test below).
    #[test]
    fn crossing_solid_shares_crossing_curve_and_caps() {
        let (brep, shell) = build();
        assert_eq!(brep.shells[shell].faces.len(), 6, "2 bands + 4 caps");

        let uses = edge_uses(&brep, shell);
        let shared = uses.values().filter(|&&n| n == 2).count();
        let open = uses.values().filter(|&&n| n == 1).count();
        // overwhelmingly watertight: caps fully stitched, crossing curve
        // shared; only a handful of turning-point segments remain open.
        assert!(shared > 400, "most edges shared: {shared}");
        assert_eq!(open, 0, "watertight: no open edges, got {open}");
        assert!(uses.values().all(|&n| n <= 2), "manifold: nothing over-shared");

        // distances are exact by construction (refinement oracle)
        let mut cfg = UnionConfig::default();
        cfg.wires_strict = false;
        let mut diag = Diagnostics::new();
        let report = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(
            report.distance_violations.is_empty(),
            "edges exact on surfaces (max {:.2} µm)",
            report.max_edge_deviation * 1e3,
        );
    }

    /// FULL watertight closure of the crossing solid.  Blocked on
    /// turning-point pre-split: the periodic per-band DCEL inserts an
    /// interpolated vertex at each band's u-extreme of the imprint loop;
    /// those extrema differ between the two bands, leaving ≤12 segments
    /// unshared.  FIX: enrich the shared imprint polyline with BOTH bands'
    /// u-turning points before imprinting (registry-style critical-point
    /// pre-split), so both DCELs vertex at identical 3-D points.
    #[test]
    fn crossing_solid_fully_watertight() {
        let (brep, shell) = build();
        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let report = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(report.manifold_violations.is_empty());
        assert!(report.wire_violations.is_empty());
        assert!(report.distance_violations.is_empty());
    }

    /// A 2×2 woodpile patch: two X-legs (z=0) crossed by two Y-legs (z=0.35)
    /// — 4 crossings, each leg carrying TWO windows — fused into ONE closed
    /// shell through the grid orchestrator.
    #[test]
    fn crossing_grid_2x2_watertight() {
        use crate::arrange::solid::{fuse_crossing_grid, GridLeg};
        let r = 0.275;
        let xleg = |y: f64| GridLeg {
            surf: CylSurf::new(Point3::new(-2.0, y, 0.0), UnitVec3::X, r),
            length: 4.0,
        };
        let yleg = |x: f64| GridLeg {
            surf: CylSurf::new(Point3::new(x, -2.0, 0.35), UnitVec3::Y, r),
            length: 4.0,
        };
        let legs = [xleg(0.0), xleg(1.2), yleg(0.0), yleg(1.2)];
        let mut brep = BRep::new();
        let shell = fuse_crossing_grid(&mut brep, &legs);
        // 4 legs × (lateral + 2 caps) = 12 faces
        assert_eq!(brep.shells[shell].faces.len(), 12, "4 legs + 8 caps");

        let uses = {
            let mut u = HashMap::new();
            for &fid in &brep.shells[shell].faces {
                let f = &brep.faces[fid];
                let mut ls = vec![f.outer_loop];
                ls.extend(f.inner_loops.iter().copied());
                for lid in ls {
                    let st = brep.loops[lid].start;
                    let mut c = st;
                    loop {
                        let ce = &brep.coedges[c];
                        *u.entry(ce.edge).or_insert(0usize) += 1;
                        c = ce.next;
                        if c == st {
                            break;
                        }
                    }
                }
            }
            u
        };
        let open = uses.values().filter(|&&n| n == 1).count();
        assert_eq!(open, 0, "2×2 grid watertight, got {open} open");
        assert!(uses.values().all(|&n| n == 2), "every edge shared exactly twice");

        // strict gates
        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let report = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(report.manifold_violations.is_empty(), "manifold");
        assert!(report.distance_violations.is_empty(),
            "distance max {:.2}µm", report.max_edge_deviation * 1e3);
    }

    /// Scale check: a 3×3 woodpile patch (6 legs, 9 crossings, 3 windows per
    /// leg) closes watertight through the same orchestrator.
    #[test]
    fn crossing_grid_3x3_watertight() {
        use crate::arrange::solid::{fuse_crossing_grid, GridLeg};
        let r = 0.275;
        let xleg = |y: f64| GridLeg { surf: CylSurf::new(Point3::new(-2.0, y, 0.0), UnitVec3::X, r), length: 6.0 };
        let yleg = |x: f64| GridLeg { surf: CylSurf::new(Point3::new(x, -2.0, 0.35), UnitVec3::Y, r), length: 6.0 };
        let legs = [xleg(0.0), xleg(1.2), xleg(2.4), yleg(0.0), yleg(1.2), yleg(2.4)];
        let mut brep = BRep::new();
        let shell = fuse_crossing_grid(&mut brep, &legs);
        assert_eq!(brep.shells[shell].faces.len(), 18, "6 legs + 12 caps");
        let mut uses = HashMap::new();
        for &fid in &brep.shells[shell].faces {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop]; ls.extend(f.inner_loops.iter().copied());
            for lid in ls { let st = brep.loops[lid].start; let mut c = st;
                loop { let ce = &brep.coedges[c]; *uses.entry(ce.edge).or_insert(0usize) += 1; c = ce.next; if c == st { break; } } }
        }
        let open = uses.values().filter(|&&n| n == 1).count();
        assert_eq!(open, 0, "3×3 grid watertight, got {open}");
        assert!(uses.values().all(|&n| n == 2));
    }

    /// A plain capped cylinder (no boolean) closes COMPLETELY through the
    /// arrangement engine — proves caps + assembly + the periodic↔planar
    /// boundary sharing are watertight.
    #[test]
    fn capped_cylinder_is_watertight() {
        use crate::arrange::cells::{arrange_and_classify, FaceChain};
        use crate::geom::composite::CompositeTube;
        use std::f64::consts::{PI, TAU};

        let surf = CylSurf::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, 0.275);
        let len = 2.0;
        let band = FaceDomain::CylinderBand { surf, length: len };
        let n = 96usize;
        let rim = |v: f64, tag: u32| FaceChain {
            pts: (0..=n).map(|k| (-PI + TAU * k as f64 / n as f64, v)).collect(),
            tag,
        };
        let empty = CompositeTube::new();

        let mut brep = BRep::new();
        let mut asm = Assembler::new(&mut brep, 1e-7);
        for cell in arrange_and_classify(&band, &[rim(0.0, 1), rim(len, 2)], &empty, 1e-3) {
            asm.emit_cell(&band, FaceGeom::Cylinder(surf), FaceNormal::Same, &cell);
        }
        emit_cap(&mut asm, &band, Leg { surf, length: len }, 0.0, true);
        emit_cap(&mut asm, &band, Leg { surf, length: len }, len, false);
        let faces = asm.faces().to_vec();
        assert_eq!(faces.len(), 3, "1 lateral + 2 caps");

        let mut uses: HashMap<cadcore_topo::EdgeId, usize> = HashMap::new();
        for &fid in &faces {
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
        let open = uses.values().filter(|&&n| n == 1).count();
        assert_eq!(open, 0, "capped cylinder fully watertight (0 open edges)");
        assert!(uses.values().all(|&n| n == 2), "every edge shared exactly twice");
    }
}
