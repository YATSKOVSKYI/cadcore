//! Per-face cell arrangement and keep/drop classification.
//!
//! A face's parameter rectangle is partitioned by its boundary wires plus
//! the curve pieces that land on it (`cadcore-geom`'s DCEL does the 2-D
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
    // A face that does NOT wrap around its periodic direction (a torus patch
    // between two crossing windows, a band a previous Boolean split along its
    // seam) is not periodic at all: its chains fit in one continuous window of
    // the parameter, and the periodic DCEL would only add a seam — which, when
    // the face's own boundary runs along that seam (two faces of one torus
    // sharing the equator), coincides with a chain and wrecks the cells.  Such
    // a face is unwrapped onto one branch and arranged as a plain region.
    let unwrapped: Option<Vec<FaceChain>> = if plane {
        None
    } else {
        unwrap_non_wrapping(chains, periodic_axis, period)
    };
    let plane = plane || unwrapped.is_some();
    let chains: &[FaceChain] = unwrapped.as_deref().unwrap_or(chains);

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
            // The DCEL hands back CCW outer loops / CW holes in ITS axes.  On a
            // swapped domain those axes are (v, u), and swapping back mirrors
            // the loops: CCW in (v, u) is CW in (u, v) — the face would be wound
            // INWARD (∂u × ∂v is the outward normal on every carrier here).
            // Reverse the traversal so the loop is CCW in the domain's own uv.
            if swap {
                let rev = |steps: &mut Vec<cadcore_geom::arrangement::LoopStep>| {
                    steps.reverse();
                    for st in steps.iter_mut() {
                        st.pts.reverse();
                        st.forward = !st.forward;
                    }
                };
                rev(&mut cell.outer);
                for h in cell.holes.iter_mut() {
                    rev(h);
                }
            }
            Some(ClassifiedCell { cell, keep })
        })
        .collect()
}

/// Tag of the half-period line [`arrange_and_classify_halved`] adds.
pub const HALF_TAG: u32 = u32::MAX - 7;

/// [`arrange_and_classify_with`] for a face that is to be emitted as TWO
/// non-periodic halves: besides the parametric seam, the face is also cut
/// along the line half a period away (the `u = π` ruling of a cylinder, the
/// inner equator of a torus).
///
/// A face that wraps its surface completely has no boundary to anchor its
/// parameterisation, and every consumer invents one: the STEP reader of an
/// OpenCASCADE-based CAD inserts a seam of its own choosing and splits
/// whatever crosses it — per face, so an edge that a long window curve shares
/// with the neighbouring face is split on one side only and the shell opens.
/// Two half faces, each a plain four-sided patch bounded by real ruling-line
/// edges, leave nothing to invent; it is what the mainstream kernels write
/// for a full cylinder for the same reason.
///
/// The caller must have put a vertex on every chain where it crosses EITHER
/// line (see `universal::split_piece_at_seam`), so the two halves and the
/// faces across the chains subdivide identically.
pub fn arrange_and_classify_halved(
    domain: &FaceDomain,
    chains: &[FaceChain],
    inside_other: &dyn Fn(Point3) -> bool,
) -> Vec<ClassifiedCell> {
    if matches!(domain, FaceDomain::Plane { .. } | FaceDomain::TorusBand { .. }) {
        return arrange_and_classify_with(domain, chains, 0.0, inside_other);
    }
    let (axis, period) = domain.periodic_axis();
    let other = |p: &(f64, f64)| if axis == 0 { p.1 } else { p.0 };
    let (mut lo, mut hi) = domain.open_range();
    for c in chains {
        for p in &c.pts {
            lo = lo.min(other(p));
            hi = hi.max(other(p));
        }
    }
    let half = 0.5 * period;
    let mut all: Vec<FaceChain> = chains.to_vec();
    all.push(FaceChain {
        pts: if axis == 0 { vec![(half, lo), (half, hi)] } else { vec![(lo, half), (hi, half)] },
        tag: HALF_TAG,
    });
    arrange_and_classify_with(domain, &all, 0.0, inside_other)
}

/// Does this set of chains describe a face that does NOT wrap its periodic
/// direction (see [`arrange_and_classify_with`])?
pub fn is_non_wrapping(domain: &FaceDomain, chains: &[FaceChain]) -> bool {
    if matches!(domain, FaceDomain::Plane { .. } | FaceDomain::TorusBand { .. }) {
        return true;
    }
    let (axis, period) = domain.periodic_axis();
    unwrap_non_wrapping(chains, axis, period).is_some()
}

/// If no chain wraps the periodic direction and together they leave a clear
/// gap in it, return the chains shifted onto ONE continuous branch (each chain
/// unwrapped along itself, then moved by whole periods to start just past the
/// gap).  `None` when the face genuinely wraps.
fn unwrap_non_wrapping(chains: &[FaceChain], axis: usize, period: f64) -> Option<Vec<FaceChain>> {
    const MARGIN: f64 = 0.2;
    if chains.is_empty() || !(period > 0.0) {
        return None;
    }
    let get = |p: &(f64, f64)| if axis == 0 { p.0 } else { p.1 };
    // per-chain continuous unwrap
    let mut un: Vec<Vec<f64>> = Vec::with_capacity(chains.len());
    for c in chains {
        let mut vals: Vec<f64> = Vec::with_capacity(c.pts.len());
        for (i, p) in c.pts.iter().enumerate() {
            let x = get(p);
            if i == 0 {
                vals.push(x);
            } else {
                let prev = vals[i - 1];
                let mut d = (x - prev) % period;
                if d > 0.5 * period {
                    d -= period;
                } else if d < -0.5 * period {
                    d += period;
                }
                vals.push(prev + d);
            }
        }
        let (lo, hi) = vals.iter().fold((f64::MAX, f64::MIN), |a, &x| (a.0.min(x), a.1.max(x)));
        if hi - lo > period - MARGIN {
            return None; // this chain wraps (a rim ring)
        }
        un.push(vals);
    }
    // widest circular gap of everything, modulo the period
    let mut all: Vec<f64> = un.iter().flatten().map(|x| x.rem_euclid(period)).collect();
    if all.is_empty() {
        return None;
    }
    all.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let (mut gap, mut base) = (0.0_f64, 0.0_f64);
    for i in 0..all.len() {
        let (a, b) = (all[i], if i + 1 < all.len() { all[i + 1] } else { all[0] + period });
        if b - a > gap {
            gap = b - a;
            base = b.rem_euclid(period); // the occupied stretch starts here
        }
    }
    if gap < MARGIN {
        return None;
    }
    let start = base - 0.5 * gap.min(MARGIN); // tolerance below the first sample
    let out = chains
        .iter()
        .zip(un.iter())
        .map(|(c, vals)| {
            // shift so the chain's first sample falls in [start, start + period)
            let k = ((vals[0] - start) / period).floor();
            let shift = -k * period;
            let pts = c
                .pts
                .iter()
                .zip(vals.iter())
                .map(|(p, &x)| if axis == 0 { (x + shift, p.1) } else { (p.0, x + shift) })
                .collect();
            FaceChain { pts, tag: c.tag }
        })
        .collect();
    Some(out)
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
