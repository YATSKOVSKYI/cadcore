//! **Lattice Booleans — the one-call API.**
//!
//! Everything a caller needs to turn swept filaments into a finished part:
//!
//! ```text
//! filaments ──fuse──▶ one closed solid ──trim (planes)──▶ flat terminals ──merge──▶ export
//! ```
//!
//! * [`fuse_and_trim`] — the whole chain from filaments;
//! * [`trim_shell`] — the trim + edge-merge half, for a shell that is already
//!   fused (it lives in an arena the caller built);
//! * the pieces themselves are re-exported at the crate root:
//!   [`fuse_filaments_with`](crate::fuse_filaments_with),
//!   [`half_space_cut_robust`](crate::half_space_cut_robust),
//!   [`merge_shell_edges`](crate::merge_shell_edges).
//!
//! The result is a single watertight, manifold shell with planar section caps
//! exactly on the (reported) trim planes.  Mating bodies — electrode plates —
//! are NOT united with it: place them on [`TrimInfo::plane`], the plane each
//! cut actually used.

use cadcore_geom::Plane3;
use cadcore_topo::{BRep, ShellId};

use crate::arrange::merge::merge_shell_edges;
use crate::arrange::scaffold::Filament;
use crate::arrange::universal::{fuse_filaments_with, FuseOptions, FuseReport};
use crate::boolean::half_space_cut_robust;

/// What one trim plane did.
#[derive(Debug, Clone, Copy)]
pub struct TrimInfo {
    /// The plane as requested (material is kept on the side its normal points to).
    pub requested: Plane3,
    /// The plane actually cut along: `requested` moved by [`shift`](Self::shift)
    /// along its normal.
    pub plane: Plane3,
    /// Distance the plane was moved to clear a grazing feature (mm; `0.0` when
    /// the requested position cut cleanly, otherwise ±0.01 … ±0.04).
    pub shift: f64,
    /// Planar section caps created.
    pub caps: usize,
    /// Faces of the input dropped whole.
    pub dropped: usize,
    /// Manifold violations left after this cut (`0` = watertight).
    pub violations: usize,
}

/// A trimmed shell in its own arena.
pub struct TrimOutcome {
    /// The trimmed B-Rep.
    pub brep: BRep,
    /// Its single closed shell.
    pub shell: ShellId,
    /// One entry per plane, in the order given.
    pub trims: Vec<TrimInfo>,
    /// Per-step edges folded back into whole curves by the merge pass.
    pub merged_edges: usize,
}

impl TrimOutcome {
    /// `true` when every cut came out watertight.
    pub fn is_watertight(&self) -> bool {
        self.trims.iter().all(|t| t.violations == 0)
    }
}

/// A fused and trimmed lattice.
pub struct LatticeOutcome {
    /// The finished B-Rep (one solid, one shell).
    pub brep: BRep,
    /// Its shell.
    pub shell: ShellId,
    /// What the fuse did.
    pub fuse: FuseReport,
    /// What each trim plane did (empty when no planes were given).
    pub trims: Vec<TrimInfo>,
    /// Per-step edges folded back into whole curves.
    pub merged_edges: usize,
}

impl LatticeOutcome {
    /// `true` when no trace was abandoned and every cut is watertight.
    pub fn is_watertight(&self) -> bool {
        self.fuse.aborted == 0 && self.trims.iter().all(|t| t.violations == 0)
    }
}

/// Trim the fused `shell` of `src` with every plane in turn (keeping the side
/// each normal points to), then fold the per-step edges into whole curves.
///
/// Each plane goes through [`half_space_cut_robust`]: a plane that grazes a
/// junction rim or an end cap is moved by 10–40 µm and the move is reported in
/// [`TrimInfo::shift`].  `knit` is the model knit tolerance (1e-6 for a shell
/// built by the universal fuse).
///
/// Returns `None` when `planes` is empty or the shell does not exist — the
/// caller then keeps what it has.
pub fn trim_shell(src: &BRep, shell: ShellId, planes: &[Plane3], knit: f64) -> Option<TrimOutcome> {
    if planes.is_empty() {
        return None;
    }
    let mut owned: Option<(BRep, ShellId)> = None;
    let mut trims: Vec<TrimInfo> = Vec::with_capacity(planes.len());
    for requested in planes {
        let (b, sh) = match &owned {
            Some((b, s)) => (b, *s),
            None => (src, shell),
        };
        let faces = b.shells.get(sh)?.faces.clone();
        let robust = half_space_cut_robust(b, &faces, requested, knit)?;
        let n = requested.normal();
        trims.push(TrimInfo {
            requested: *requested,
            plane: Plane3::from_origin_normal(requested.frame.origin + n.as_vec() * robust.shift, n),
            shift: robust.shift,
            caps: robust.outcome.caps,
            dropped: robust.outcome.dropped,
            violations: robust.violations,
        });
        owned = Some((robust.outcome.brep, robust.outcome.shell));
    }
    let (mut brep, shell) = owned?;
    let merged_edges = merge_shell_edges(&mut brep, shell);
    Some(TrimOutcome { brep, shell, trims, merged_edges })
}

/// Fuse `filaments` into one solid and trim it with `planes`.
///
/// The whole "sinter the part, then cut its terminals" operation in one call:
/// the universal fuse ([`fuse_filaments_with`]), a robust planar cut per plane
/// and the edge-merge pass.  With no planes the fused solid is returned with
/// its edges merged.
///
/// ```no_run
/// use cadcore_geom::Plane3;
/// use cadcore_math::{Point3, UnitVec3};
/// use cadcore_union::{fuse_and_trim, Filament, FuseOptions};
///
/// # let filaments: Vec<Filament> = Vec::new();
/// // keep 1 mm ≤ y ≤ 5 mm
/// let planes = [
///     Plane3::from_origin_normal(Point3::new(0.0, 1.0, 0.0), UnitVec3::Y),
///     Plane3::from_origin_normal(Point3::new(0.0, 5.0, 0.0), -UnitVec3::Y),
/// ];
/// let part = fuse_and_trim(&filaments, &FuseOptions::default(), &planes, 1.0e-6);
/// assert!(part.is_watertight());
/// for t in &part.trims {
///     println!("cut at {:?}, moved {} mm", t.plane.frame.origin, t.shift);
/// }
/// ```
pub fn fuse_and_trim(filaments: &[Filament], options: &FuseOptions, planes: &[Plane3], knit: f64) -> LatticeOutcome {
    let mut brep = BRep::new();
    let (shell, fuse) = fuse_filaments_with(&mut brep, filaments, options);
    match trim_shell(&brep, shell, planes, knit) {
        Some(t) => LatticeOutcome { brep: t.brep, shell: t.shell, fuse, trims: t.trims, merged_edges: t.merged_edges },
        None => {
            let merged_edges = merge_shell_edges(&mut brep, shell);
            LatticeOutcome { brep, shell, fuse, trims: Vec::new(), merged_edges }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrange::fixtures::gen;
    use cadcore_math::{Point3, UnitVec3};

    /// The public one-call API on a real generator lattice (6 mm × 4 layers).
    #[test]
    fn fuse_and_trim_grid_is_one_watertight_solid() {
        let fils = gen::grid_6mm_4l();
        let planes = [
            Plane3::from_origin_normal(Point3::new(0.0, 1.0, 0.0), UnitVec3::Y),
            Plane3::from_origin_normal(Point3::new(0.0, 5.0, 0.0), -UnitVec3::Y),
        ];
        let part = fuse_and_trim(&fils, &FuseOptions::default(), &planes, 1.0e-6);
        assert!(part.is_watertight(), "fuse {:?} trims {:?}", part.fuse, part.trims);
        assert_eq!(part.trims.len(), 2);
        assert!(part.trims.iter().all(|t| t.caps > 0));
        assert_eq!(part.brep.solids.len(), 1);
        assert!(part.merged_edges > 0);
        let rep = crate::validate::manifold::check(&part.brep, part.shell);
        assert!(rep.violations.is_empty(), "{:?}", &rep.violations[..rep.violations.len().min(4)]);
        // the caps lie exactly on the planes the cuts report
        for t in &part.trims {
            assert!((t.plane.frame.origin - t.requested.frame.origin).length() <= 0.04 + 1e-12);
        }
    }

    /// No planes → the fused solid, edges merged.
    #[test]
    fn fuse_and_trim_without_planes_returns_the_fused_solid() {
        let fils = gen::pair_3mm();
        let part = fuse_and_trim(&fils, &FuseOptions::default(), &[], 1.0e-6);
        assert!(part.is_watertight());
        assert!(part.trims.is_empty());
        let rep = crate::validate::manifold::check(&part.brep, part.shell);
        assert!(rep.violations.is_empty());
    }
}
