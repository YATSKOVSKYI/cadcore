//! The **general** Boolean union — `union(A, B)` over two arbitrary B-Reps.
//!
//! Unlike the scaffold-shaped orchestrators in [`crate::arrange`] (which take
//! filament legs + elbows), this layer makes NO assumption about what the
//! input solids are.  It is the Parasolid-grade universal entry point:
//!
//! 1. [`aabb`] — broad phase: which faces of A could meet which faces of B.
//! 2. **SSI** — surface×surface intersection of each candidate pair → curves
//!    (the predictor–corrector tracer in [`crate::geom::intersect`]).
//! 3. **imprint** — split each face by the curves on it ([`crate::arrange`]'s
//!    DCEL arrangement).
//! 4. [`contain`] — classify each face-piece: inside the other solid → drop,
//!    outside → keep (union rule), coincident → keep one copy.
//! 5. **stitch** — feed kept pieces to one assembler; shared SSI edges weld.
//!
//! Built incrementally; this module currently lands the broad phase and the
//! universal point-in-solid classifier (steps 1 and 4), which are the
//! solid-agnostic primitives the rest of the pipeline composes.

pub mod aabb;
pub mod contain;
pub mod ssi;
pub mod union;

#[cfg(test)]
pub(crate) mod tests_support;

pub use aabb::{candidate_pairs, face_aabb, Aabb};
pub use contain::point_in_solid;
pub use ssi::{intersect_faces, SsiCurve};
pub use union::{union, union_many, union_n};

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_math::{Point3, UnitVec3};
    use cadcore_topo::BRep;
    use tests_support::{axis_box, capped_cylinder};

    #[test]
    fn point_in_axis_box() {
        let mut brep = BRep::new();
        let faces = axis_box(&mut brep, Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 2.0));
        // inside
        assert!(point_in_solid(&brep, &faces, Point3::new(1.0, 1.0, 1.0)));
        assert!(point_in_solid(&brep, &faces, Point3::new(0.1, 0.1, 0.1)));
        assert!(point_in_solid(&brep, &faces, Point3::new(1.9, 1.5, 0.3)));
        // outside
        assert!(!point_in_solid(&brep, &faces, Point3::new(3.0, 1.0, 1.0)));
        assert!(!point_in_solid(&brep, &faces, Point3::new(-0.1, 1.0, 1.0)));
        assert!(!point_in_solid(&brep, &faces, Point3::new(1.0, 1.0, 2.5)));
        assert!(!point_in_solid(&brep, &faces, Point3::new(1.0, -0.5, 1.0)));
    }

    #[test]
    fn point_in_offset_box() {
        // a box NOT at the origin and non-cubic, to catch axis/offset bugs
        let mut brep = BRep::new();
        let faces = axis_box(
            &mut brep,
            Point3::new(-3.0, 5.0, -1.0),
            Point3::new(1.0, 6.0, 4.0),
        );
        assert!(point_in_solid(&brep, &faces, Point3::new(-1.0, 5.5, 2.0)));
        assert!(!point_in_solid(&brep, &faces, Point3::new(-1.0, 4.9, 2.0)));
        assert!(!point_in_solid(&brep, &faces, Point3::new(2.0, 5.5, 2.0)));
    }

    #[test]
    fn point_in_capped_cylinder() {
        let mut brep = BRep::new();
        let faces = capped_cylinder(
            &mut brep,
            Point3::new(0.0, 0.0, 0.0),
            UnitVec3::Z,
            1.0,
            4.0,
        );
        // inside
        assert!(point_in_solid(&brep, &faces, Point3::new(0.0, 0.0, 2.0)));
        assert!(point_in_solid(&brep, &faces, Point3::new(0.5, 0.5, 0.1)));
        // outside radially
        assert!(!point_in_solid(&brep, &faces, Point3::new(1.5, 0.0, 2.0)));
        // outside axially
        assert!(!point_in_solid(&brep, &faces, Point3::new(0.0, 0.0, 4.5)));
        assert!(!point_in_solid(&brep, &faces, Point3::new(0.0, 0.0, -0.5)));
    }

    /// Count how many times each edge is used across a shell's faces.
    fn edge_use_counts(brep: &BRep, shell: cadcore_topo::ShellId) -> std::collections::HashMap<cadcore_topo::EdgeId, usize> {
        let mut uses = std::collections::HashMap::new();
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

    /// Two axis boxes overlapping at a corner, unioned → ONE watertight shell
    /// (every edge used exactly twice).
    #[test]
    fn union_two_corner_boxes_watertight() {
        let mut a = BRep::new();
        let fa = axis_box(&mut a, Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 2.0));
        let mut b = BRep::new();
        let fb = axis_box(&mut b, Point3::new(1.0, 1.0, 1.0), Point3::new(3.0, 3.0, 3.0));

        let (out, shell) = union(&a, &fa, &b, &fb, 1e-6);
        let uses = edge_use_counts(&out, shell);
        assert!(!uses.is_empty(), "union produced faces");
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "watertight: every edge used twice ({} open of {})", open, uses.len());
    }

    /// A round peg through a plate: a cylinder passing fully through a box
    /// (caps outside, footprint inside).  The box's top & bottom faces get a
    /// circular hole, the cylinder's lateral splits into kept end-bands —
    /// fused into ONE watertight shell.
    #[test]
    fn union_cylinder_through_box_watertight() {
        let mut a = BRep::new();
        let fa = axis_box(&mut a, Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 2.0));
        let mut b = BRep::new();
        // axis Z, centre (2,2), r=1, from z=-1 to z=3 → pokes through the plate
        let fb = capped_cylinder(&mut b, Point3::new(2.0, 2.0, -1.0), UnitVec3::Z, 1.0, 4.0);

        let (out, shell) = union(&a, &fa, &b, &fb, 1e-6);
        let uses = edge_use_counts(&out, shell);
        assert!(!uses.is_empty(), "union produced faces");
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "watertight: every edge used twice ({} open of {})", open, uses.len());
    }

    /// Two perpendicular cylinders crossing through each other (a +), unequal
    /// radii to dodge the equal-radius tangency.  Each lateral band is carved
    /// by the other; fused into ONE watertight shell.
    #[test]
    fn union_two_crossing_cylinders_watertight() {
        let mut a = BRep::new();
        let fa = capped_cylinder(&mut a, Point3::new(-3.0, 0.0, 0.0), UnitVec3::X, 1.0, 6.0);
        let mut b = BRep::new();
        let fb = capped_cylinder(&mut b, Point3::new(0.0, 0.0, -3.0), UnitVec3::Z, 0.6, 6.0);

        let (out, shell) = union(&a, &fa, &b, &fb, 1e-6);
        let uses = edge_use_counts(&out, shell);
        assert!(!uses.is_empty(), "union produced faces");
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "watertight: every edge used twice ({} open of {})", open, uses.len());
    }

    /// Composability: union of THREE boxes via `union_many` (each intermediate
    /// result, all `Trimmed` faces, feeds the next union) → ONE watertight
    /// shell.  This is the N-solid path the scaffold needs.
    #[test]
    fn union_three_boxes_watertight() {
        let mut a = BRep::new();
        let fa = axis_box(&mut a, Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 2.0));
        let mut b = BRep::new();
        let fb = axis_box(&mut b, Point3::new(1.0, 1.0, 1.0), Point3::new(3.0, 3.0, 3.0));
        let mut c = BRep::new();
        let fc = axis_box(&mut c, Point3::new(2.0, 2.0, 2.0), Point3::new(4.0, 4.0, 4.0));

        let (out, shell) = union_many(&[(&a, &fa), (&b, &fb), (&c, &fc)], 1e-6).expect("3 solids");
        let uses = edge_use_counts(&out, shell);
        assert!(!uses.is_empty(), "union produced faces");
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "watertight: every edge used twice ({} open of {})", open, uses.len());
    }

    /// A 3-cylinder woodpile (axes X/Y/Z through the origin, distinct radii)
    /// via the direct N-way [`union_n`] (behind `union_many`): each original
    /// cylinder is classified against ALL others in ONE pass — no
    /// re-processing of intermediates, so the crowded triple-overlap centre
    /// stays watertight.  This is the scaffold pattern.
    #[test]
    fn union_three_cylinder_woodpile_watertight() {
        let mut a = BRep::new();
        let fa = capped_cylinder(&mut a, Point3::new(-3.0, 0.0, 0.0), UnitVec3::X, 1.0, 6.0);
        let mut b = BRep::new();
        let fb = capped_cylinder(&mut b, Point3::new(0.0, -3.0, 0.0), UnitVec3::Y, 0.8, 6.0);
        let mut c = BRep::new();
        let fc = capped_cylinder(&mut c, Point3::new(0.0, 0.0, -3.0), UnitVec3::Z, 0.6, 6.0);

        let (out, shell) = union_many(&[(&a, &fa), (&b, &fb), (&c, &fc)], 1e-6).expect("3 solids");
        let uses = edge_use_counts(&out, shell);
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "watertight: every edge used twice ({} open of {})", open, uses.len());
    }

    /// Union a scaffold L-filament (two legs + a real **torus** elbow + caps)
    /// with a cylinder crossing one LEG.  The elbow is not cut by the crosser,
    /// so this validates that the engine carries a torus face THROUGH a union
    /// (domain, boundary, emission, welding to its mating legs) AND classifies
    /// a point inside a trimmed cylinder band correctly.
    #[test]
    fn union_filament_with_elbow_carries_torus() {
        use crate::arrange::filament::fuse_poly_filament;
        use crate::arrange::scaffold::Filament;

        // L-filament in z=0: leg along X (x:-2→0), 90° elbow, leg along Y (y:0→2)
        let fil = Filament::serpentine(
            &[
                Point3::new(-2.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
            ],
            0.3,
            0.5,
        );
        let mut a = BRep::new();
        let shell_a = fuse_poly_filament(&mut a, &fil.legs, &fil.elbows);
        let fa: Vec<_> = a.shells[shell_a].faces.clone();
        // sanity: the elbow torus face is present
        assert!(
            fa.iter().any(|&f| matches!(a.faces[f].geom, cadcore_topo::FaceGeom::Torus(_))),
            "filament has a torus elbow"
        );

        // crosser: a cylinder along Z through the X-leg at x=-1 (clear of elbow)
        let mut b = BRep::new();
        let fb = capped_cylinder(&mut b, Point3::new(-1.0, 0.0, -1.0), UnitVec3::Z, 0.2, 2.0);

        let (out, shell) = union(&a, &fa, &b, &fb, 1e-6);
        let uses = edge_use_counts(&out, shell);
        assert!(!uses.is_empty(), "union produced faces");
        // the torus survives into the output
        assert!(
            out.shells[shell].faces.iter().any(|&f| matches!(out.faces[f].geom, cadcore_topo::FaceGeom::Torus(_))),
            "torus elbow carried through the union"
        );
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "watertight: every edge used twice ({} open of {})", open, uses.len());
    }

    /// A thin cylinder passing through a thicker one, axes intersecting at the
    /// thick cylinder's centre (different radii) → ONE watertight shell.
    #[test]
    fn union_thin_cylinder_through_thick_watertight() {
        let mut a = BRep::new();
        let fa = capped_cylinder(&mut a, Point3::new(-2.0, 0.0, 0.0), UnitVec3::X, 0.3, 2.0);
        let mut b = BRep::new();
        let fb = capped_cylinder(&mut b, Point3::new(-1.0, 0.0, -1.0), UnitVec3::Z, 0.2, 2.0);
        let (out, shell) = union(&a, &fa, &b, &fb, 1e-6);
        let uses = edge_use_counts(&out, shell);
        let open = uses.values().filter(|&&n| n != 2).count();
        assert_eq!(open, 0, "thin-through-thick: {} open of {}", open, uses.len());
    }

    /// Union an L-filament with a cylinder crossing its **torus elbow** — the
    /// crosser CUTS the torus (torus×cylinder SSI), not just a leg.
    ///
    /// IGNORED — the torus×cyl SSI + φ-seam pre-split infrastructure is in
    /// place, but cutting the torus is not watertight yet (~498 open of 792):
    /// the coarse marching-squares intersection loop, the `TorusPatch`
    /// arrangement with a window across the φ-seam, and the Trimmed-torus pcurve
    /// emission all need dedicated work.  Torus PASS-THROUGH (elbow not cut) IS
    /// validated — see `union_filament_with_elbow_carries_torus`.
    #[test]
    #[ignore = "torus-cut not watertight yet; needs TorusPatch arrangement + emission work"]
    fn union_cylinder_cuts_torus_elbow() {
        use crate::arrange::filament::fuse_poly_filament;
        use crate::arrange::scaffold::Filament;
        use cadcore_topo::FaceGeom;

        // L-filament in z=0; the elbow centre is near (−0.5, 0.5) with the arc
        // bulging toward (0,0)→up.  Place the crosser through the elbow tube.
        let fil = Filament::serpentine(
            &[
                Point3::new(-2.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 2.0, 0.0),
            ],
            0.3,
            0.5,
        );
        let elbow = fil.elbows[0].clone();
        let mut a = BRep::new();
        let shell_a = fuse_poly_filament(&mut a, &fil.legs, &fil.elbows);
        let fa: Vec<_> = a.shells[shell_a].faces.clone();
        assert!(fa.iter().any(|&f| matches!(a.faces[f].geom, FaceGeom::Torus(_))), "has elbow");

        // aim the crosser through the MID-ARC of the elbow tube (away from the
        // leg junctions at the arc ends), perpendicular to the torus plane.
        let tor = elbow.surf;
        let th = 0.5 * (elbow.theta_lo + elbow.theta_hi);
        let ring = tor.frame.origin
            + (tor.frame.x * th.cos() + tor.frame.y * th.sin()) * tor.major_radius;
        let mut b = BRep::new();
        let base = Point3::new(ring.x, ring.y, -1.0);
        let fb = capped_cylinder(&mut b, base, UnitVec3::Z, 0.2, 2.0);

        let (out, shell) = union(&a, &fa, &b, &fb, 1e-6);
        let uses = edge_use_counts(&out, shell);
        assert!(!uses.is_empty(), "union produced faces");
        let openv: Vec<_> = uses.iter().filter(|(_, &n)| n != 2).map(|(e, &n)| (*e, n)).collect();
        if !openv.is_empty() && std::env::var("CADCORE_DUMP_UNION").is_ok() {
            for (eid, n) in openv.iter().take(12) {
                let e = &out.edges[*eid];
                let g = match &e.geom { cadcore_topo::EdgeGeom::Circle(_) => "circle", cadcore_topo::EdgeGeom::Line(_) => "line", cadcore_topo::EdgeGeom::Ellipse(_) => "ellipse", cadcore_topo::EdgeGeom::Polyline(_) => "poly" };
                let p0 = out.vertices[e.v_start].point; let p1 = out.vertices[e.v_end].point;
                eprintln!("open {}x {}: ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})", n, g, p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
            }
        }
        assert_eq!(openv.len(), 0, "watertight: every edge used twice ({} open of {})", openv.len(), uses.len());
    }

    #[test]
    fn broad_phase_culls_distant_faces() {
        let mut a = BRep::new();
        let fa = axis_box(&mut a, Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
        let mut b = BRep::new();
        // overlapping box
        let fb = axis_box(&mut b, Point3::new(0.5, 0.5, 0.5), Point3::new(1.5, 1.5, 1.5));
        let pairs = candidate_pairs(&a, &fa, &b, &fb, 1e-6);
        assert!(!pairs.is_empty(), "overlapping boxes have candidate pairs");

        // far box: no candidates
        let mut c = BRep::new();
        let fc = axis_box(&mut c, Point3::new(10.0, 10.0, 10.0), Point3::new(11.0, 11.0, 11.0));
        let pairs = candidate_pairs(&a, &fa, &c, &fc, 1e-6);
        assert!(pairs.is_empty(), "distant boxes culled");
    }
}
