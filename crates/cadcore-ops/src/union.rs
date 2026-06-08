//! Analytic-surface B-Rep Boolean **UNION** for DIW woodpile scaffolds.
//!
//! # Goal
//!
//! Fuse the independent, mutually-overlapping analytic sweep solids produced by
//! the filament sweep pipeline into a **single watertight manifold B-Rep solid**,
//! so downstream CAD tools (ANSYS SpaceClaim) never need to run "Combine".
//!
//! This is **not** a generic OCCT-style boolean kernel.  It targets exactly the
//! scaffold geometry:
//!
//! * axis-aligned, **equal-radius cylinders** (the filament straight runs),
//! * **planar end caps**,
//! * (v1) **no** sphere caps and **no** torus corner fillets — solids containing
//!   those are left un-fused and reported.
//!
//! # Method (see the approved plan)
//!
//! 1. broad-phase (uniform grid) → candidate intersecting primitive pairs,
//! 2. surface–surface intersection curves ([`cadcore_geom::intersect`]),
//! 3. **2D parameter-domain** trimming per face (project curves into `(u,v)`,
//!    build a trimming graph, extract loops),
//! 4. classify each region by closed-form primitive containment (drop interior),
//! 5. weld shared edges, stitch one outward-wound shell, validate manifold.
//!
//! Phases 1–5 land incrementally; see the per-phase tasks.

use cadcore_topo::BRep;

/// Fuse every solid in `brep` into a single analytic union solid.
///
/// Returns the number of solids remaining after the union.
///
/// **Phase 0 (current):** identity passthrough — performs only the topological
/// grouping of [`BRep::fuse_solids`] (shells merged into one [`Solid`]); the
/// geometric trim/classify/weld stages are added in later phases.  Output is
/// therefore not yet watertight where filaments overlap.
pub fn union_solids(brep: &mut BRep) -> usize {
    brep.fuse_solids();
    brep.solids.len()
}
