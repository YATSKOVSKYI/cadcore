//! Per-face cell arrangement and keep/drop classification.
//!
//! A face's parameter rectangle is partitioned by its boundary wires plus
//! the registry pieces that land on it (`cadcore-geom`'s DCEL does the 2-D
//! work, seam handling included).  Each resulting cell is classified by
//! lifting an interior point to 3-D and asking whether it sits inside the
//! OTHER solid — material kept by the union keeps the cell, swallowed
//! regions drop it.

use cadcore_geom::arrangement::{build_periodic_u, Arrangement, ArrangementOptions, Cell, ChainInput};
use cadcore_math::Point3;

use super::domain::FaceDomain;
use crate::geom::composite::CompositeTube;

/// A classified cell.
#[derive(Debug, Clone)]
pub struct ClassifiedCell {
    /// The 2-D cell (loops carry source tags).
    pub cell: Cell,
    /// `true` ⇒ the cell's material survives the union (outside the other
    /// solid).
    pub keep: bool,
}

/// Curve input for one face: uv polyline + stable tag (piece id, boundary
/// id…).  Tags survive into [`Cell`] loop steps, which is how the assembly
/// stage maps cell edges back to shared 3-D pieces.
#[derive(Debug, Clone)]
pub struct FaceChain {
    /// uv polyline.
    pub pts: Vec<(f64, f64)>,
    /// Stable source id.
    pub tag: u32,
}

/// Reserved tag for the seam chains the DCEL inserts on periodic domains.
pub const SEAM_TAG: u32 = u32::MAX;

/// Build the face arrangement and classify every cell against `other`.
///
/// The domain's periodic axis is handled by the DCEL's seam machinery; the
/// open axis is bounded by the caller's boundary chains (the face's own
/// rims must be among `chains`).
pub fn arrange_and_classify(
    domain: &FaceDomain,
    chains: &[FaceChain],
    other: &CompositeTube,
    probe_offset: f64,
) -> Vec<ClassifiedCell> {
    arrange_and_classify_with(domain, chains, probe_offset, &|p| other.contains(p))
}

/// As [`arrange_and_classify`] but with an arbitrary "swallowed" predicate —
/// a cell drops when its lifted interior point satisfies `inside_other`
/// (e.g. inside ANY of several crossing tubes).
pub fn arrange_and_classify_with(
    domain: &FaceDomain,
    chains: &[FaceChain],
    probe_offset: f64,
    inside_other: &dyn Fn(Point3) -> bool,
) -> Vec<ClassifiedCell> {
    // A plane has no periodic axis → use the plain (non-periodic) DCEL, so
    // boundary circles are NOT split by a phantom seam.  A full-donut TorusBand
    // is likewise treated non-periodically: it is cut at BOTH the θ-seam and
    // the φ-seam into a [0,2π]² rectangle whose seam edges weld (a 2-periodic
    // surface can't be handled by the single-seam periodic DCEL).
    let plane = matches!(domain, FaceDomain::Plane { .. } | FaceDomain::TorusBand { .. });
    let (periodic_axis, period) = domain.periodic_axis();
    let open = domain.open_range();

    // The DCEL is periodic in u; when the face is periodic in v (torus
    // patch), swap axes for the arrangement and swap back on lift.
    let swap = periodic_axis == 1 && !plane;
    let to_arr = |&(u, v): &(f64, f64)| if swap { (v, u) } else { (u, v) };
    let from_arr = |(x, y): (f64, f64)| if swap { (y, x) } else { (x, y) };

    let inputs: Vec<ChainInput> = chains
        .iter()
        .map(|c| ChainInput {
            pts: c.pts.iter().map(to_arr).collect(),
            tag: c.tag,
        })
        .collect();

    let arr = if plane {
        Arrangement::build(&inputs, &ArrangementOptions::default())
    } else {
        build_periodic_u(&inputs, period, open, SEAM_TAG, &ArrangementOptions::default())
    };
    // Loop-step points come back in ARRANGEMENT coords; map them back to the
    // domain's (u,v) so consumers (assembly) lift them correctly.  For
    // swapped (v-periodic / torus) domains this un-swaps θ↔φ.
    let unswap_step = |st: &mut cadcore_geom::arrangement::LoopStep| {
        for p in st.pts.iter_mut() {
            *p = from_arr(*p);
        }
    };
    arr.cells()
        .into_iter()
        .filter_map(|mut cell| {
            let (iu, iv) = cell.interior_point()?;
            let (u, v) = from_arr((iu, iv));
            let p: Point3 = domain.lift(u, v);
            // probe slightly inside the material side: ON the surface is
            // ambiguous, so test the surface point itself against the other
            // tube (strict containment ⇒ swallowed).
            let _ = probe_offset;
            let keep = !inside_other(p);
            for st in cell.outer.iter_mut() {
                unswap_step(st);
            }
            for h in cell.holes.iter_mut() {
                for st in h.iter_mut() {
                    unswap_step(st);
                }
            }
            Some(ClassifiedCell { cell, keep })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::composite::closed_loops_between;
    use crate::geom::intersect::TraceOptions;
    use cadcore_geom::CylSurf;
    use cadcore_math::UnitVec3;
    use std::f64::consts::TAU;

    /// The production primitive end-to-end on the new engine: a leg crossed
    /// by a perpendicular leg one layer up.  The lateral face must split
    /// into the kept band (outside the crosser) and the dropped window
    /// (inside it), with the window's boundary being the traced loop.
    #[test]
    fn crossing_window_cell_classification() {
        let leg = CylSurf::new(Point3::new(-2.0, 0.0, 0.0), UnitVec3::X, 0.275);
        let a = CompositeTube::new().with_leg(leg, 4.0);
        let b = CompositeTube::new().with_leg(
            CylSurf::new(Point3::new(0.0, -2.0, 0.35), UnitVec3::Y, 0.275),
            4.0,
        );

        let domain = FaceDomain::CylinderBand {
            surf: leg,
            length: 4.0,
        };

        // boundary rims of the band in uv: v=0 and v=length, full period
        let rim = |v: f64, tag: u32| FaceChain {
            pts: (0..=64)
                .map(|k| (-std::f64::consts::PI + TAU * k as f64 / 64.0, v))
                .collect(),
            tag,
        };

        // the crossing loop, projected into this face's uv (NOTE: the leg's
        // frame origin is at x=-2, the crossing sits at x=0 ⇒ v≈2)
        let loops = closed_loops_between(&a, &b, &TraceOptions::default());
        assert_eq!(loops.len(), 1);
        let mut window: Vec<(f64, f64)> = loops[0].points.iter().map(|&p| domain.uv(p)).collect();
        // close the chain explicitly for the DCEL
        let first = window[0];
        window.push(first);

        let chains = vec![
            rim(0.0, 1),
            rim(4.0, 2),
            FaceChain {
                pts: window,
                tag: 10,
            },
        ];

        let cells = arrange_and_classify(&domain, &chains, &b, 1e-3);
        let kept: Vec<_> = cells.iter().filter(|c| c.keep).collect();
        let dropped: Vec<_> = cells.iter().filter(|c| !c.keep).collect();
        // The window may be split by the parametric seam — counts depend on
        // where frame.x points.  Invariants that always hold:
        assert!(!kept.is_empty(), "kept band exists");
        assert!(!dropped.is_empty(), "dropped window exists");
        // every dropped cell is bounded only by the loop and/or the seam
        for d in &dropped {
            assert!(
                d.cell
                    .outer
                    .iter()
                    .all(|s| s.tag == 10 || s.tag == SEAM_TAG),
                "window bounded by loop/seam only"
            );
        }
        // dropped area ≈ the swallowed window: compare against the loop's
        // own uv area via the shoelace over its (closed) chain
        let window_pts = &chains[2].pts;
        let mut wa = 0.0;
        for i in 0..window_pts.len() - 1 {
            let (x0, y0) = window_pts[i];
            let (x1, y1) = window_pts[i + 1];
            wa += x0 * y1 - x1 * y0;
        }
        let wa = (wa / 2.0).abs();
        let da: f64 = dropped.iter().map(|c| c.cell.area.abs()).sum();
        assert!(
            (da - wa).abs() < 0.05 * wa.max(1e-9),
            "dropped area {da:.4} ≈ window area {wa:.4}"
        );
        // cells tile the whole band
        let total: f64 = cells.iter().map(|c| c.cell.area).sum();
        assert!(
            (total - TAU * 4.0).abs() < 1e-2,
            "cells tile the domain: {total} vs {}",
            TAU * 4.0
        );
    }
}
