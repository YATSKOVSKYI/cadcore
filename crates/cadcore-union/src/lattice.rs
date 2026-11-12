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
//!   [`crate::fuse_filaments_with`],
//!   [`crate::half_space_cut_robust`],
//!   [`crate::merge_shell_edges`].
//!
//! The result is a single watertight, manifold shell with planar section caps
//! exactly on the (reported) trim planes.  Mating bodies — electrode plates —
//! are NOT united with it: place them on [`TrimInfo::plane`], the plane each
//! cut actually used.

use cadcore_geom::{CylSurf, Plane3, TorusSurf};
use cadcore_math::Point3;
use cadcore_topo::{BRep, ShellId};

use crate::geom::composite::CompositeTube;

use crate::arrange::merge::merge_shell_edges;
use crate::arrange::filament::Filament;
use crate::arrange::universal::{fuse_filaments_with, FuseOptions, FuseReport};

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

/// The same solid as [`fuse_and_trim`] builds, as a **point-membership model**:
/// the union of the filament tubes intersected with the kept half-spaces.
///
/// A volume mesher (voxels, marching tetrahedra) does not need the boundary
/// representation — it needs to ask "is this point material?" — and it must
/// get the answer for the SAME solid the exact B-Rep describes, or the mesh
/// and the exported model drift apart at every trimmed end.  This is that
/// oracle: exact tubes (legs and torus elbows, no chord approximation) and
/// flat trim planes.
///
/// ```no_run
/// use cadcore_union::{fuse_and_trim, Filament, FuseOptions, LatticeModel};
/// # let filaments: Vec<Filament> = Vec::new();
/// # let planes = [];
/// let part = fuse_and_trim(&filaments, &FuseOptions::default(), &planes, 1.0e-6);
/// // membership for exactly the solid that was cut (planes as actually used)
/// let model = LatticeModel::for_trims(&filaments, &part.trims);
/// let (lo, hi) = model.bounds();
/// # let _ = (lo, hi);
/// ```
#[derive(Debug, Clone)]
pub struct LatticeModel {
    tubes: Vec<CompositeTube>,
    keep: Vec<Plane3>,
    lo: [f64; 3],
    hi: [f64; 3],
}

impl LatticeModel {
    /// The union of `filaments` kept on the normal side of every plane in `keep`.
    pub fn new(filaments: &[Filament], keep: &[Plane3]) -> Self {
        Self::inflated(filaments, keep, 0.0)
    }

    /// As [`new`](Self::new), on the planes a trim actually used (a robust cut
    /// may have moved them — see [`TrimInfo::plane`]).
    pub fn for_trims(filaments: &[Filament], trims: &[TrimInfo]) -> Self {
        let keep: Vec<Plane3> = trims.iter().map(|t| t.plane).collect();
        Self::new(filaments, &keep)
    }

    /// As [`new`](Self::new) with every tube radius grown by `extra` (mm).
    ///
    /// A voxel mesher marks a cell conductive by its centre; growing the tubes
    /// by half a cell makes a cell conductive when the cell itself reaches the
    /// filament, which keeps thin inter-layer bonds connected on a coarse
    /// grid.  The trim planes are NOT moved.
    pub fn inflated(filaments: &[Filament], keep: &[Plane3], extra: f64) -> Self {
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        let mut grow = |p: Point3, r: f64| {
            for (k, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                lo[k] = lo[k].min(c - r);
                hi[k] = hi[k].max(c + r);
            }
        };
        let mut tubes = Vec::with_capacity(filaments.len());
        for f in filaments {
            let mut t = CompositeTube::new();
            for l in &f.legs {
                let r = l.surf.radius + extra;
                let surf = CylSurf { frame: l.surf.frame, radius: r };
                grow(surf.frame.origin, r);
                grow(surf.frame.origin + surf.axis().as_vec() * l.length, r);
                t = t.with_leg(surf, l.length);
            }
            for e in &f.elbows {
                let r = e.surf.minor_radius + extra;
                let surf = TorusSurf { frame: e.surf.frame, major_radius: e.surf.major_radius, minor_radius: r };
                for k in 0..=8 {
                    let th = e.theta_lo + (e.theta_hi - e.theta_lo) * k as f64 / 8.0;
                    let spine = surf.frame.origin
                        + (surf.frame.x * th.cos() + surf.frame.y * th.sin()) * surf.major_radius;
                    // the chord sagitta of an eighth of the arc is covered by the pad
                    grow(spine, r + surf.major_radius * 0.02);
                }
                t = t.with_elbow(surf, e.theta_lo, e.theta_hi);
            }
            tubes.push(t);
        }
        Self { tubes, keep: keep.to_vec(), lo, hi }
    }

    /// Is `p` inside the solid (in some tube and on the kept side of every plane)?
    pub fn contains(&self, p: Point3) -> bool {
        self.keep.iter().all(|pl| pl.signed_distance(p) >= 0.0) && self.tubes.iter().any(|t| t.contains(p))
    }

    /// Axis-aligned box enclosing the untrimmed tubes (`lo`, `hi`).
    pub fn bounds(&self) -> (Point3, Point3) {
        (Point3::new(self.lo[0], self.lo[1], self.lo[2]), Point3::new(self.hi[0], self.hi[1], self.hi[2]))
    }

    /// The kept half-spaces.
    pub fn keep_planes(&self) -> &[Plane3] {
        &self.keep
    }

    /// Volume by midpoint sampling on an `n[0] × n[1] × n[2]` grid over
    /// [`bounds`](Self::bounds) — an independent check for the volume a CAD
    /// system reports for the exported B-Rep.
    pub fn sampled_volume(&self, n: [usize; 3]) -> f64 {
        let (lo, hi) = (self.lo, self.hi);
        if n.iter().any(|&k| k == 0) || (0..3).any(|k| hi[k] <= lo[k]) {
            return 0.0;
        }
        let mut inside = 0usize;
        for i in 0..n[0] {
            let x = lo[0] + (hi[0] - lo[0]) * (i as f64 + 0.5) / n[0] as f64;
            for j in 0..n[1] {
                let y = lo[1] + (hi[1] - lo[1]) * (j as f64 + 0.5) / n[1] as f64;
                for k in 0..n[2] {
                    let z = lo[2] + (hi[2] - lo[2]) * (k as f64 + 0.5) / n[2] as f64;
                    if self.contains(Point3::new(x, y, z)) {
                        inside += 1;
                    }
                }
            }
        }
        let cell = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]) / (n[0] * n[1] * n[2]) as f64;
        inside as f64 * cell
    }
}

/// Trim the fused `shell` of `src` with every plane in turn (keeping the side
/// each normal points to), then fold the per-step edges into whole curves.
///
/// Each plane goes through [`crate::half_space_cut_robust`]: a plane that grazes a
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
        // open edges the shell already has (0 for a sound fuse): the cut is
        // judged only on the ones it would add
        let baseline = crate::validate::manifold::check(b, sh).violations.len();
        let robust = crate::boolean::halfspace::half_space_cut_robust_over(b, &faces, requested, knit, baseline)?;
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
    let t_merge = std::time::Instant::now();
    let merged_edges = merge_shell_edges(&mut brep, shell);
    if std::env::var("CADCORE_TIME").is_ok() {
        eprintln!("[cut][time] merge_shell_edges: {} ms", t_merge.elapsed().as_millis());
    }
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

    /// The membership model describes the solid the B-Rep bounds: a point just
    /// inside a section cap is material, one just outside is not, and the
    /// sampled volume is that of the lattice (20.66 mm³ for this fixture).
    #[test]
    fn membership_model_matches_the_trimmed_solid() {
        let fils = gen::grid_6mm_4l();
        let planes = [
            Plane3::from_origin_normal(Point3::new(0.0, 1.0, 0.0), UnitVec3::Y),
            Plane3::from_origin_normal(Point3::new(0.0, 5.0, 0.0), -UnitVec3::Y),
        ];
        let part = fuse_and_trim(&fils, &FuseOptions::default(), &planes, 1.0e-6);
        let model = LatticeModel::for_trims(&fils, &part.trims);
        // the first layer-1 leg runs along y at x = 0.275, z = 0.35
        assert!(model.contains(Point3::new(0.275, 1.001, 0.35)));
        assert!(!model.contains(Point3::new(0.275, 0.999, 0.35)));
        assert!(model.contains(Point3::new(0.275, 4.999, 0.35)));
        assert!(!model.contains(Point3::new(0.275, 5.001, 0.35)));
        assert!(!model.contains(Point3::new(0.9, 3.0, 0.35))); // between two legs
        let v = model.sampled_volume([160, 160, 44]);
        assert!((v - 20.66).abs() < 0.15, "sampled volume {v}");
        // inflating only grows it; the planes stay put
        let fat = LatticeModel::inflated(&fils, model.keep_planes(), 0.05);
        assert!(fat.contains(Point3::new(0.275 + 0.3, 3.0, 0.35)));
        assert!(!fat.contains(Point3::new(0.275, 0.999, 0.35)));
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
