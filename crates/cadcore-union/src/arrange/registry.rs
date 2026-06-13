//! The pre-split curve registry — README invariant #2, mechanised.
//!
//! An intersection loop between two filaments crosses face boundaries
//! (member junction circles) of BOTH.  Every face that consumes a piece of
//! the loop must agree *exactly* on where the piece starts and ends; two
//! faces partitioning a shared circle from slightly different points produce
//! unshared arcs and an open shell.
//!
//! The registry therefore splits every loop ONCE, globally:
//!
//! 1. walk the loop, tracking the owner member on each side
//!    ([`CompositeTube::member_of`]);
//! 2. at every owner transition, find the junction circle being crossed and
//!    land the split point exactly on it ∩ the other tube
//!    ([`refine_on_circle`]);
//! 3. emit [`Piece`]s whose endpoints ARE those refined joints — every
//!    consumer (face of either solid) references the same piece by id and
//!    the same joint coordinates by construction.

use cadcore_math::Point3;

use crate::geom::composite::{CompositeTube, Junction, MemberRef};
use crate::geom::intersect::TracedCurve;
use crate::geom::refine::{refine_on_circle, Projectable};
use crate::tolerance::Tolerances;

/// Which solid a member reference belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Side {
    /// First tube passed to [`split_loop`].
    A,
    /// Second tube.
    B,
}

/// A joint between consecutive pieces: lies exactly on a junction circle of
/// one tube AND on the surface of the other.
#[derive(Debug, Clone)]
pub struct Joint {
    /// The refined joint position — THE single source of truth.
    pub point: Point3,
    /// Which tube's junction is crossed.
    pub side: Side,
    /// The junction being crossed.
    pub junction: Junction,
}

/// One piece of an intersection loop: constant owner pair.
#[derive(Debug, Clone)]
pub struct Piece {
    /// Points in loop order; first/last are refined joints (or the loop's
    /// only piece when no junction is crossed).
    pub points: Vec<Point3>,
    /// Owning member on tube A.
    pub owner_a: MemberRef,
    /// Owning member on tube B.
    pub owner_b: MemberRef,
}

/// A split loop: pieces in cyclic order with the joints between them.
/// `joints[k]` sits between `pieces[k]` and `pieces[(k+1) % n]`.
#[derive(Debug, Clone)]
pub struct SplitLoop {
    /// Pieces in cyclic order.
    pub pieces: Vec<Piece>,
    /// Joints between consecutive pieces (same length as `pieces`; empty
    /// when the loop never leaves one owner pair).
    pub joints: Vec<Joint>,
}

/// Split one closed intersection loop at every member-junction crossing of
/// both tubes.  Returns `None` when owner classification fails (a loop point
/// off both tubes — upstream tracing guarantees this never happens for its
/// own output).
pub fn split_loop(
    curve: &TracedCurve,
    a: &CompositeTube,
    b: &CompositeTube,
    tol: &Tolerances,
) -> Option<SplitLoop> {
    let pts = &curve.points;
    let n = pts.len();
    if n < 4 {
        return None;
    }
    let junctions_a = a.junctions(tol.weld);
    let junctions_b = b.junctions(tol.weld);

    // Owner pair per point.
    let owners: Vec<(MemberRef, MemberRef)> = pts
        .iter()
        .map(|&p| Some((a.member_of(p)?, b.member_of(p)?)))
        .collect::<Option<Vec<_>>>()?;

    // Transition indices: owner pair changes between i and i+1 (cyclic).
    let mut cuts: Vec<usize> = (0..n)
        .filter(|&i| owners[i] != owners[(i + 1) % n])
        .collect();
    if cuts.is_empty() {
        // single-owner loop (an ordinary crossing window)
        return Some(SplitLoop {
            pieces: vec![Piece {
                points: pts.clone(),
                owner_a: owners[0].0,
                owner_b: owners[0].1,
            }],
            joints: Vec::new(),
        });
    }
    cuts.sort_unstable();

    // Refine each transition onto its junction circle ∩ other tube.
    let mut joints: Vec<Joint> = Vec::with_capacity(cuts.len());
    for &i in &cuts {
        let j = (i + 1) % n;
        let raw = pts[i] + (pts[j] - pts[i]) * 0.5;
        // which side switched members?
        let (side, m_from, m_to) = if owners[i].0 != owners[j].0 {
            (Side::A, owners[i].0, owners[j].0)
        } else {
            (Side::B, owners[i].1, owners[j].1)
        };
        let (tube, other, juncs): (&CompositeTube, &CompositeTube, &[Junction]) = match side {
            Side::A => (a, b, &junctions_a),
            Side::B => (b, a, &junctions_b),
        };
        let _ = tube;
        let junction = juncs
            .iter()
            .find(|jc| {
                (jc.a == m_from && jc.b == m_to) || (jc.a == m_to && jc.b == m_from)
            })
            .copied()
            // members may hop over a short member (sampling skipped it):
            // fall back to the junction circle nearest to the raw crossing
            .or_else(|| {
                juncs
                    .iter()
                    .min_by(|x, y| {
                        let dx = (x.circle.frame.origin - raw).length();
                        let dy = (y.circle.frame.origin - raw).length();
                        dx.total_cmp(&dy)
                    })
                    .copied()
            })?;
        let point = refine_on_circle(raw, &junction.circle, &[], tol.refine, 60);
        // also land on the other tube (junction circle has the final word)
        let point = {
            let mut q = point;
            for _ in 0..40 {
                let before = q;
                q = other.project_p(q);
                q = refine_on_circle(q, &junction.circle, &[], tol.refine, 4);
                if (q - before).length() < tol.refine {
                    break;
                }
            }
            q
        };
        joints.push(Joint {
            point,
            side,
            junction,
        });
    }

    // Assemble pieces between consecutive cuts.
    let m = cuts.len();
    let mut pieces: Vec<Piece> = Vec::with_capacity(m);
    for k in 0..m {
        let from_cut = cuts[k];
        let to_cut = cuts[(k + 1) % m];
        let start_owner = (from_cut + 1) % n;
        let mut points = vec![joints[k].point];
        let mut i = start_owner;
        loop {
            points.push(pts[i]);
            if i == to_cut {
                break;
            }
            i = (i + 1) % n;
        }
        points.push(joints[(k + 1) % m].point);
        pieces.push(Piece {
            points,
            owner_a: owners[start_owner].0,
            owner_b: owners[start_owner].1,
        });
    }
    // rotate joints so joints[k] follows pieces[k]
    joints.rotate_left(1);
    Some(SplitLoop { pieces, joints })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::composite::closed_loops_between;
    use crate::geom::intersect::TraceOptions;
    use cadcore_geom::{CylSurf, TorusSurf};
    use cadcore_math::UnitVec3;

    fn l_filament(z: f64, r: f64) -> CompositeTube {
        let leg_a = CylSurf::new(Point3::new(-2.0, 0.0, z), UnitVec3::X, r);
        let t = TorusSurf::new(Point3::new(0.0, 0.5, z), UnitVec3::Z, 0.5, r);
        let leg_b = CylSurf::new(Point3::new(0.5, 0.5, z), UnitVec3::Y, r);
        CompositeTube::new()
            .with_leg(leg_a, 2.0)
            .with_leg(leg_b, 2.0)
            .with_elbow_between(t, Point3::new(0.0, 0.0, z), Point3::new(0.5, 0.5, z))
    }

    fn stacked_rotated_l() -> CompositeTube {
        let leg_a = CylSurf::new(Point3::new(0.0, -2.0, 0.35), UnitVec3::Y, 0.275);
        let t = TorusSurf::new(Point3::new(0.5, 0.0, 0.35), UnitVec3::Z, 0.5, 0.275);
        let leg_b = CylSurf::new(Point3::new(0.5, 0.5, 0.35), UnitVec3::X, 0.275);
        CompositeTube::new()
            .with_leg(leg_a, 2.0)
            .with_leg(leg_b, 2.0)
            .with_elbow_between(t, Point3::new(0.0, 0.0, 0.35), Point3::new(0.5, 0.5, 0.35))
    }

    /// Junction discovery on the L: one junction per leg↔elbow contact.
    #[test]
    fn junctions_of_l_filament() {
        let f = l_filament(0.0, 0.275);
        let js = f.junctions(1e-6);
        assert_eq!(js.len(), 2, "{js:?}");
        for j in &js {
            assert!((j.circle.radius - 0.275).abs() < 1e-12);
        }
    }

    /// The corner U-on-U loop splits into pieces whose joints are EXACT:
    /// on the junction circle AND on both tube surfaces.
    #[test]
    fn corner_loop_pieces_have_exact_joints() {
        let a = l_filament(0.0, 0.275);
        let b = stacked_rotated_l();
        let tol = Tolerances::default();
        let loops = closed_loops_between(&a, &b, &TraceOptions::default());
        assert_eq!(loops.len(), 1);
        let split = split_loop(&loops[0], &a, &b, &tol).expect("split");
        assert!(
            split.pieces.len() >= 2,
            "corner loop must cross junctions: {} pieces",
            split.pieces.len()
        );
        assert_eq!(split.pieces.len(), split.joints.len());
        for (k, joint) in split.joints.iter().enumerate() {
            // exactly on the junction circle
            let w = joint.point - joint.junction.circle.frame.origin;
            let planar = joint.junction.circle.frame.z.dot_vec(w).abs();
            let radial = (w.length() - joint.junction.circle.radius).abs();
            assert!(planar < 1e-9 && radial < 1e-9, "joint {k} off circle");
            // on both tube surfaces
            assert!(
                a.distance_p(joint.point) < 1e-7 && b.distance_p(joint.point) < 1e-7,
                "joint {k}: dA={:.2e} dB={:.2e}",
                a.distance_p(joint.point),
                b.distance_p(joint.point)
            );
            // and it IS the shared endpoint of the two adjacent pieces
            let p_this = &split.pieces[k];
            let p_next = &split.pieces[(k + 1) % split.pieces.len()];
            assert_eq!(*p_this.points.last().unwrap(), joint.point);
            assert_eq!(p_next.points[0], joint.point);
        }
        // owner pairs differ across every joint
        for k in 0..split.pieces.len() {
            let p_this = &split.pieces[k];
            let p_next = &split.pieces[(k + 1) % split.pieces.len()];
            assert!(
                p_this.owner_a != p_next.owner_a || p_this.owner_b != p_next.owner_b,
                "joint {k} between identical owner pairs"
            );
        }
    }

    /// A plain leg×leg crossing never leaves one owner pair: one piece, no
    /// joints.
    #[test]
    fn plain_crossing_is_single_piece() {
        let a = CompositeTube::new().with_leg(
            CylSurf::new(Point3::new(-2.0, 0.0, 0.0), UnitVec3::X, 0.275),
            4.0,
        );
        let b = CompositeTube::new().with_leg(
            CylSurf::new(Point3::new(0.0, -2.0, 0.35), UnitVec3::Y, 0.275),
            4.0,
        );
        let tol = Tolerances::default();
        let loops = closed_loops_between(&a, &b, &TraceOptions::default());
        assert_eq!(loops.len(), 1);
        let split = split_loop(&loops[0], &a, &b, &tol).expect("split");
        assert_eq!(split.pieces.len(), 1);
        assert!(split.joints.is_empty());
    }
}
