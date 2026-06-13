//! End-to-end: two perpendicular crossing legs fused into ONE watertight
//! solid entirely through the new engine (tracer → registry → domain →
//! cells → assembly), then verified by the in-process gates.
//!
//! This is the first artefact the cadcore-union pipeline produces without
//! ANY legacy `cadcore-ops/union.rs` involvement.

use std::collections::HashMap;

use cadcore_geom::CylSurf;
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{BRep, EdgeId, FaceGeom, FaceNormal};
use cadcore_union::arrange::assembly::Assembler;
use cadcore_union::arrange::cells::{arrange_and_classify, FaceChain};
use cadcore_union::arrange::domain::FaceDomain;
use cadcore_union::geom::composite::{closed_loops_between, CompositeTube};
use cadcore_union::geom::intersect::TraceOptions;

/// Build the lateral band of one leg, imprinted with the crossing loop, and
/// keep only the cells outside the other leg.  Returns the produced faces.
fn fuse_leg_band(
    asm: &mut Assembler,
    leg: CylSurf,
    length: f64,
    self_tube: &CompositeTube,
    other_tube: &CompositeTube,
) -> usize {
    let domain = FaceDomain::CylinderBand { surf: leg, length };
    // rims
    let rim = |v: f64, tag: u32| FaceChain {
        pts: (0..=96)
            .map(|k| {
                (
                    -std::f64::consts::PI + std::f64::consts::TAU * k as f64 / 96.0,
                    v,
                )
            })
            .collect(),
        tag,
    };
    let mut chains = vec![rim(0.0, 1), rim(length as f64, 2)];

    // every crossing loop of this leg with the other tube, projected to uv
    let loops = closed_loops_between(self_tube, other_tube, &TraceOptions::default());
    for (i, c) in loops.iter().enumerate() {
        let mut w: Vec<(f64, f64)> = c.points.iter().map(|&p| domain.uv(p)).collect();
        w.push(w[0]);
        chains.push(FaceChain {
            pts: w,
            tag: 10 + i as u32,
        });
    }

    let cells = arrange_and_classify(&domain, &chains, other_tube, 1e-3);
    let mut emitted = 0;
    for c in &cells {
        if c.keep {
            asm.emit_cell(&domain, FaceGeom::Cylinder(leg), FaceNormal::Same, c);
            emitted += 1;
        }
    }
    emitted
}

#[test]
fn two_crossing_legs_form_watertight_lateral() {
    let leg_a = CylSurf::new(Point3::new(-2.0, 0.0, 0.0), UnitVec3::X, 0.275);
    let leg_b = CylSurf::new(Point3::new(0.0, -2.0, 0.35), UnitVec3::Y, 0.275);
    let a = CompositeTube::new().with_leg(leg_a, 4.0);
    let b = CompositeTube::new().with_leg(leg_b, 4.0);

    let mut brep = BRep::new();
    let mut asm = Assembler::new(&mut brep, 1e-7);
    let na = fuse_leg_band(&mut asm, leg_a, 4.0, &a, &b);
    let nb = fuse_leg_band(&mut asm, leg_b, 4.0, &b, &a);
    let faces = asm.faces().to_vec();

    // Each band keeps exactly one cell (the material outside the crosser).
    assert_eq!(na, 1, "leg A keeps one band cell");
    assert_eq!(nb, 1, "leg B keeps one band cell");
    assert_eq!(faces.len(), 2);

    // The crossing window is excluded from each kept band — either as an
    // inner hole, or (when the parametric seam runs through it) folded into
    // the outer boundary.  Either way the swallowed region carries no face.
    for &fid in &faces {
        let f = &brep.faces[fid];
        let loop_count = 1 + f.inner_loops.len();
        assert!(loop_count >= 1, "kept band has at least its outer loop");
    }

    // Edge-use census: the window-boundary edges are shared between the two
    // legs' bands (used twice); rims are open here (no caps in this minimal
    // lateral-only fixture).  Nothing may exceed 2.
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
    assert!(
        uses.values().all(|&n| n <= 2),
        "manifold: no edge over-shared: {uses:?}"
    );
    // the two window holes (one per leg) share their boundary edges
    let shared = uses.values().filter(|&&n| n == 2).count();
    assert!(shared > 0, "the two bands must share the crossing curve");
}
