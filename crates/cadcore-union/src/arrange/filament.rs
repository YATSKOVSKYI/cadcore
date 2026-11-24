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






/// Smallest admissible elbow bend radius, as a multiple of the tube radius.
///
/// A centre line bent tighter than the tube's own radius sweeps a *spindle*
/// torus — the inner surface folds through itself and the shell is no longer
/// a solid.  At exactly the tube radius the torus is a *horn* (its hole closes
/// to a point), which is just as unusable.  5 % of margin keeps a real hole.
pub const MIN_BEND_FACTOR: f64 = 1.05;

/// Fraction of a leg that its fillets may consume; the rest stays a cylinder.
const LEG_SHARE: f64 = 0.9;

/// Why a polyline cannot be swept as a filament of straight legs and torus
/// elbows.
#[derive(Debug, Clone, PartialEq)]
pub enum PathError {
    /// Fewer than two distinct vertices.
    Degenerate,
    /// A corner that no elbow of an admissible radius fits, and that cannot be
    /// merged with a neighbour (a reversal, or two tight turns of opposite
    /// sense).  `vertex` indexes the resolved vertex list at the time of failure.
    TightCorner {
        /// Where the corner is.
        at: Point3,
        /// Largest bend radius the corner can host.
        capacity: f64,
        /// Smallest admissible bend radius.
        required: f64,
    },
    /// An elbow of a hand-built filament is bent tighter than the tube allows.
    SpindleElbow {
        /// Index of the elbow.
        elbow: usize,
        /// Its bend radius.
        bend: f64,
        /// The tube radius.
        tube: f64,
    },
    /// A leg of a hand-built filament has no length left between its fillets.
    EmptyLeg {
        /// Index of the leg.
        leg: usize,
    },
}

impl std::fmt::Display for PathError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PathError::Degenerate => write!(f, "path has fewer than two distinct vertices"),
            PathError::TightCorner { at, capacity, required } => write!(
                f,
                "corner at ({:.4}, {:.4}, {:.4}) can host a bend radius of {capacity:.4} at most; the tube needs {required:.4}",
                at.x, at.y, at.z
            ),
            PathError::SpindleElbow { elbow, bend, tube } => {
                write!(f, "elbow {elbow}: bend radius {bend:.4} is not larger than the tube radius {tube:.4} (spindle torus)")
            }
            PathError::EmptyLeg { leg } => write!(f, "leg {leg} has no length left between its elbows"),
        }
    }
}

impl std::error::Error for PathError {}

/// A polyline made sweepable: its vertices and one bend radius per interior
/// vertex, every one of them at least [`MIN_BEND_FACTOR`] × the tube radius.
#[derive(Debug, Clone)]
pub struct ResolvedPath {
    /// Path vertices (sharp corners of the centre line).
    pub verts: Vec<Point3>,
    /// Bend radius at `verts[k + 1]`.
    pub fillets: Vec<f64>,
    /// Pairs of corners replaced by one corner because the leg between them
    /// was too short to host both elbows.
    pub merged: usize,
    /// End stubs removed because they were too short for the elbow next to them.
    pub dropped_ends: usize,
}

/// Turn of the path at an interior vertex: `tan(turn / 2)`, or `None` when the
/// vertex is straight.
fn tan_half_turn(a: Point3, b: Point3, c: Point3) -> Option<f64> {
    let (din, dout) = (b - a, c - b);
    let (li, lo) = (din.length(), dout.length());
    if li < 1.0e-12 || lo < 1.0e-12 {
        return None;
    }
    let turn = (din.dot(dout) / (li * lo)).clamp(-1.0, 1.0).acos();
    (turn > 1.0e-7).then(|| (turn * 0.5).tan())
}

/// Make a polyline sweepable as a filament of radius `r`.
///
/// Every interior vertex gets a torus elbow, all of one bend radius: `fillet`,
/// or the largest radius every corner of the path can host.  What a corner can
/// host follows from its legs, whose tangent offsets `R·tan(turn/2)` may take
/// 90 % of each leg.  A leg between two corners is split evenly
/// between them; only when the tube cannot follow the path on those terms is
/// it shared in proportion to what the two corners need (a hairpin next to a
/// gentle turn takes most of the leg, as it must).
///
/// A corner whose legs cannot host even the smallest admissible radius
/// ([`MIN_BEND_FACTOR`]·`r`) is not a corner the tube can follow, and it must
/// not be allowed to set the radius of the path.  It arises
/// where two turns sit closer together than a bead is wide — a boundary stub
/// of a zig-zag next to a box corner, say.  The two turns are then one bend:
/// the stub's end points are replaced by the intersection of the legs on
/// either side of it, and a single elbow rolls through the compound corner.
/// That elbow is never tighter than the circle inscribed between the stub and
/// its neighbours, so the path stays on the inner side of the original
/// polyline.  A stub at a free end is dropped instead: it lies inside the bead
/// of the leg it leads into.
///
/// Fails with [`PathError::TightCorner`] when a tight corner cannot be merged:
/// a reversal, or two tight turns of opposite sense (an S-jog narrower than
/// the bead).  Such a path is not a single clean sweep.
pub fn resolve_path(verts: &[Point3], r: f64, fillet: f64) -> Result<ResolvedPath, PathError> {
    let r_min = r * MIN_BEND_FACTOR;
    let fillet = fillet.max(r_min);

    // distinct, non-collinear vertices
    let mut v: Vec<Point3> = Vec::with_capacity(verts.len());
    for &p in verts {
        if v.last().map_or(true, |&q: &Point3| (p - q).length() > 1.0e-9) {
            v.push(p);
        }
    }
    let straighten = |v: &mut Vec<Point3>, floor: &mut Vec<f64>| {
        let mut k = 1;
        while k + 1 < v.len() {
            if tan_half_turn(v[k - 1], v[k], v[k + 1]).is_none() && (v[k + 1] - v[k - 1]).length() > (v[k] - v[k - 1]).length() {
                v.remove(k);
                floor.remove(k);
            } else {
                k += 1;
            }
        }
    };
    // lower bound of the bend radius at a merged corner (the inscribed circle)
    let mut floor = vec![0.0f64; v.len()];
    straighten(&mut v, &mut floor);
    if v.len() < 2 {
        return Err(PathError::Degenerate);
    }

    let (mut merged, mut dropped_ends) = (0usize, 0usize);
    // Largest radius leg `j` (v[j] → v[j+1]) lets its corners reach.
    let leg_capacity = |v: &[Point3], j: usize| -> f64 {
        let n = v.len();
        let t0 = if j == 0 { 0.0 } else { tan_half_turn(v[j - 1], v[j], v[j + 1]).unwrap_or(0.0) };
        let t1 = if j + 2 >= n { 0.0 } else { tan_half_turn(v[j], v[j + 1], v[j + 2]).unwrap_or(0.0) };
        if t0 + t1 < 1.0e-12 {
            f64::INFINITY
        } else {
            (v[j + 1] - v[j]).length() * LEG_SHARE / (t0 + t1)
        }
    };

    // ── A. End stubs ────────────────────────────────────────────────────────
    // A lead-in too short for the elbow next to it lies inside the bead of the
    // leg it leads into; it goes first, whatever happens to the rest.
    loop {
        let n = v.len();
        if n >= 3 && leg_capacity(&v, 0) < r_min {
            v.remove(0);
            floor.remove(0);
        } else if n >= 3 && leg_capacity(&v, n - 2) < r_min {
            v.pop();
            floor.pop();
        } else {
            break;
        }
        dropped_ends += 1;
        straighten(&mut v, &mut floor);
    }

    // ── B. The even split ───────────────────────────────────────────────────
    // Ordinarily every corner may use 45 % of each of its legs (90 % of an end
    // leg).  The split is even on purpose: mirrored corners of a pattern then
    // get mirrored elbows, and the walls that different layers stack on one
    // another end at the same points.  If the tube can follow the path on
    // these terms, this is the answer.
    {
        let n = v.len();
        let mut even = f64::INFINITY;
        for k in 1..n.saturating_sub(1) {
            let (din, dout) = (v[k] - v[k - 1], v[k + 1] - v[k]);
            let (li, lo) = (din.length(), dout.length());
            if li < 1.0e-12 || lo < 1.0e-12 {
                continue;
            }
            let dot = (din.dot(dout) / (li * lo)).clamp(-1.0, 1.0);
            if dot > 0.999_999 {
                continue;
            }
            let tan_half = (dot.acos() * 0.5).tan();
            if tan_half <= 1.0e-9 {
                continue;
            }
            let fi = if k == 1 { 0.9 } else { 0.45 };
            let fo = if k + 2 == n { 0.9 } else { 0.45 };
            even = even.min((li * fi).min(lo * fo) / tan_half);
        }
        let bend = fillet.min(even);
        if bend >= r_min {
            let fillets = vec![bend; n.saturating_sub(2)];
            return Ok(ResolvedPath { verts: v, fillets, merged, dropped_ends });
        }
    }

    // ── C. Escalation ───────────────────────────────────────────────────────
    // Some corner cannot host the tube on an even split.  Let each leg be
    // shared in proportion to what its two corners need (a hairpin next to a
    // gentle turn takes most of the leg), and merge the corners of any leg
    // that is too short even so.
    for _ in 0..verts.len().max(4) * 2 {
        let n = v.len();
        if n < 3 {
            break;
        }
        // the leg that starves its corners the most
        let mut worst: Option<(f64, usize)> = None;
        for j in 0..n - 1 {
            let cap = leg_capacity(&v, j);
            if cap < r_min && worst.map_or(true, |(w, _)| cap < w) {
                worst = Some((cap, j));
            }
        }
        let Some((cap, j)) = worst else { break };
        if j == 0 {
            v.remove(0);
            floor.remove(0);
            dropped_ends += 1;
        } else if j + 2 == n {
            v.pop();
            floor.pop();
            dropped_ends += 1;
        } else {
            // interior stub v[j] → v[j+1]: merge its two corners into one
            let fail = PathError::TightCorner { at: v[j], capacity: cap, required: r_min };
            let (a, b, c, d) = (v[j - 1], v[j], v[j + 1], v[j + 2]);
            let (din, stub, dout) = (b - a, c - b, d - c);
            let (li, ls, lo) = (din.length(), stub.length(), dout.length());
            let (ui, us, uo) = (din * (1.0 / li), stub * (1.0 / ls), dout * (1.0 / lo));
            let (n1, n2) = (ui.cross(us), us.cross(uo));
            // same sense of rotation, in one plane
            if n1.dot(n2) <= 0.0 {
                return Err(fail);
            }
            let (a1, a2) = (ui.dot(us).clamp(-1.0, 1.0).acos(), us.dot(uo).clamp(-1.0, 1.0).acos());
            let sum = a1 + a2;
            // the outer legs must still meet ahead of the stub: a compound turn
            // of 180° or more is a hairpin, not a corner
            if sum > std::f64::consts::PI - 0.05 {
                return Err(fail);
            }
            let x = b + ui * (ls * a2.sin() / sum.sin());
            let x_back = c - uo * (ls * a1.sin() / sum.sin());
            if (x - x_back).length() > 1.0e-6 * (1.0 + ls) {
                return Err(fail); // the three segments are not coplanar
            }
            // circle inscribed between the stub and the two outer legs
            let rho = ls / ((a1 * 0.5).tan() + (a2 * 0.5).tan());
            let keep = floor[j].max(floor[j + 1]).max(rho);
            v[j] = x;
            floor[j] = keep;
            v.remove(j + 1);
            floor.remove(j + 1);
            merged += 1;
        }
        straighten(&mut v, &mut floor);
    }
    if v.len() < 2 {
        return Err(PathError::Degenerate);
    }

    // One bend radius for the whole path: the requested one, or the largest
    // every corner can host.  Corners that the design stacks on one another —
    // the walls of a honeycomb bonded to the walls below, U-turn columns —
    // then end at the same points in every layer; with a radius of its own at
    // each corner, two stacked walls would end a few microns apart and leave a
    // sliver at every crossing of their junction circles.  What must NOT set
    // this radius is a stub, and none is left: the loop above has merged or
    // dropped every leg that could not host the smallest admissible elbow.
    let n = v.len();
    let caps: Vec<f64> = (1..n.saturating_sub(1)).map(|k| leg_capacity(&v, k - 1).min(leg_capacity(&v, k))).collect();
    let mut bend = fillet;
    for (i, &cap) in caps.iter().enumerate() {
        if cap < r_min {
            return Err(PathError::TightCorner { at: v[i + 1], capacity: cap, required: r_min });
        }
        bend = bend.min(cap);
    }
    // A merged corner cannot go below the circle inscribed in the compound
    // corner it replaces, or its elbow would swing outside the original path.
    let fillets = caps.iter().enumerate().map(|(i, &cap)| bend.max(floor[i + 1]).min(cap)).collect();
    Ok(ResolvedPath { verts: v, fillets, merged, dropped_ends })
}

impl Filament {
    /// Build a filament from any polyline, resolving its corners with
    /// [`resolve_path`]: every elbow gets `fillet` where it fits and the
    /// largest radius that fits elsewhere, and corners too close together for
    /// the tube become one bend.  Fails when the path cannot be swept cleanly.
    pub fn try_serpentine(verts: &[Point3], r: f64, fillet: f64) -> Result<Filament, PathError> {
        let path = resolve_path(verts, r, fillet)?;
        let f = Self::serpentine_with_fillets(&path.verts, r, &path.fillets);
        f.check()?;
        Ok(f)
    }

    /// Check that the filament is a valid sweep: every elbow bent wider than
    /// the tube (no spindle or horn torus) and every leg of positive length.
    /// The lattice Booleans assume both; a filament that fails here fuses into
    /// a shell that is not a solid.
    pub fn check(&self) -> Result<(), PathError> {
        for (i, e) in self.elbows.iter().enumerate() {
            if e.surf.major_radius <= e.surf.minor_radius * (1.0 + 1.0e-6) {
                return Err(PathError::SpindleElbow { elbow: i, bend: e.surf.major_radius, tube: e.surf.minor_radius });
            }
        }
        for (i, l) in self.legs.iter().enumerate() {
            if !(l.length > 1.0e-9) {
                return Err(PathError::EmptyLeg { leg: i });
            }
        }
        Ok(())
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

    // ── resolve_path ────────────────────────────────────────────────────────

    fn p(x: f64, y: f64) -> Point3 {
        Point3::new(x, y, 0.0)
    }

    /// A path with room everywhere keeps its vertices and gets the requested
    /// radius at every corner.
    #[test]
    fn resolve_path_leaves_a_roomy_path_alone() {
        let verts = [p(0.0, 0.0), p(6.0, 0.0), p(6.0, 1.4), p(0.0, 1.4), p(0.0, 2.8), p(6.0, 2.8)];
        let path = resolve_path(&verts, 0.275, 0.5).unwrap();
        assert_eq!(path.verts.len(), verts.len());
        assert_eq!((path.merged, path.dropped_ends), (0, 0));
        assert!(path.fillets.iter().all(|&f| f == 0.5), "{:?}", path.fillets);
    }

    /// A leg that is a real feature of the design (a connector the tube can
    /// follow) sets ONE radius for the whole path, so stacked corners of
    /// different layers end at the same points.
    #[test]
    fn resolve_path_keeps_one_radius_along_the_path() {
        // legs 6 mm, connectors 1.4 mm, one connector 0.62 mm: its corners fit
        // 0.9·0.62/2 = 0.279, which the tube (r = 0.26) can follow
        let verts = [p(0.0, 0.0), p(6.0, 0.0), p(6.0, 1.4), p(0.0, 1.4), p(0.0, 2.02), p(6.0, 2.02)];
        let r = 0.26;
        let path = resolve_path(&verts, r, 0.29).unwrap();
        assert_eq!(path.merged, 0);
        assert!(path.fillets.iter().all(|&f| (f - 0.279).abs() < 1e-9), "{:?}", path.fillets);
        assert!(path.fillets.iter().all(|&f| f >= r * MIN_BEND_FACTOR));
    }

    /// A stub must NOT set that radius: its corners are merged and every other
    /// corner keeps the requested bend (this was the whole-layer spindle bug —
    /// one 0.37 mm stub clamped all 26 elbows of a layer to 0.167 < r).
    #[test]
    fn resolve_path_does_not_let_a_stub_set_the_radius() {
        let a = 25.0_f64.to_radians();
        let d = |len: f64| (len * a.cos(), len * a.sin());
        let (dx, dy) = d(12.0);
        let verts = [
            p(6.0, 3.0),
            p(6.0, 0.275),
            p(2.9, 0.275),
            p(0.275, 0.275),
            p(0.275, 0.646),
            p(0.275 + dx, 0.646 + dy),
            p(0.275 + dx, 0.646 + dy + 1.5),
        ];
        let path = resolve_path(&verts, 0.275, 0.29).unwrap();
        assert_eq!(path.merged, 1);
        assert!(path.fillets.iter().all(|&f| (f - 0.29).abs() < 1e-12), "{:?}", path.fillets);
    }

    /// The ±25° diamond's box corner: a wall, a 0.37 mm stub up the next wall
    /// and a diagonal.  The stub cannot host the two elbows, so its corners
    /// become one bend — tangent to the wall and to the diagonal, on the inner
    /// side of the original corner.
    #[test]
    fn resolve_path_merges_a_stub_between_two_turns_into_one_bend() {
        let (r, fillet) = (0.275, 0.29);
        let a = 25.0_f64.to_radians();
        let stub_top = p(0.275, 0.646);
        let far = p(0.275 + 14.0 * a.cos(), 0.646 + 14.0 * a.sin());
        let verts = [p(2.9, 0.275), p(0.275, 0.275), stub_top, far];
        let path = resolve_path(&verts, r, fillet).unwrap();
        assert_eq!((path.merged, path.verts.len(), path.fillets.len()), (1, 3, 1));
        // the new corner is where the wall and the diagonal meet
        let x = path.verts[1];
        assert!((x.y - 0.275).abs() < 1e-12, "on the wall line: {x:?}");
        let on_diag = (x - stub_top).cross(far - stub_top).length() / (far - stub_top).length();
        assert!(on_diag < 1e-9, "on the diagonal's line: off by {on_diag:e}");
        // one elbow of the requested radius; its tangent points lie on the
        // ORIGINAL legs, i.e. beyond the old corners as seen from the new one
        let big = path.fillets[0];
        assert!(big >= fillet - 1e-12 && big >= r * MIN_BEND_FACTOR);
        let turn = std::f64::consts::PI - a;
        let t = big * (turn * 0.5).tan();
        assert!(t >= (x - verts[1]).length() && t >= (x - stub_top).length(), "tangent {t} inside the extensions");
        // and the result is a valid sweep
        let f = Filament::try_serpentine(&verts, r, fillet).unwrap();
        assert_eq!((f.legs.len(), f.elbows.len()), (2, 1));
        assert!(f.check().is_ok());
        assert!(f.legs.iter().all(|l| l.length > 0.5));
    }

    /// A lead-in too short for the elbow next to it is dropped, not allowed to
    /// clamp the elbow.
    #[test]
    fn resolve_path_drops_an_end_stub() {
        let verts = [p(0.0, 0.0), p(0.08, 0.0), p(0.08, 5.0), p(1.5, 5.0), p(1.5, 0.0)];
        let path = resolve_path(&verts, 0.275, 0.5).unwrap();
        assert_eq!(path.dropped_ends, 1);
        assert_eq!(path.verts.len(), 4);
        assert!(path.fillets.iter().all(|&f| f >= 0.275 * MIN_BEND_FACTOR));
    }

    /// What no single tube can follow is refused, with the corner named.
    #[test]
    fn resolve_path_refuses_reversals_and_narrow_s_jogs() {
        // S-jog: two opposite turns 0.2 mm apart
        let jog = [p(0.0, 0.0), p(5.0, 0.0), p(5.0, 0.2), p(10.0, 0.2)];
        assert!(matches!(resolve_path(&jog, 0.275, 0.3), Err(PathError::TightCorner { .. })));
        // hairpin narrower than the tube: two same-sense 90° turns 0.3 mm apart
        let pin = [p(0.0, 0.0), p(5.0, 0.0), p(5.0, 0.3), p(0.0, 0.3)];
        assert!(matches!(resolve_path(&pin, 0.275, 0.3), Err(PathError::TightCorner { .. })));
        assert!(Filament::try_serpentine(&pin, 0.275, 0.3).is_err());
    }

    /// `check` catches a hand-built spindle elbow.
    #[test]
    fn check_rejects_a_spindle_elbow() {
        let verts = [p(0.0, 0.0), p(3.0, 0.0), p(3.0, 3.0)];
        let ok = Filament::serpentine(&verts, 0.275, 0.5);
        assert!(ok.check().is_ok());
        let spindle = Filament::serpentine(&verts, 0.275, 0.167);
        assert_eq!(spindle.check(), Err(PathError::SpindleElbow { elbow: 0, bend: 0.167, tube: 0.275 }));
    }

    /// The merged compound corner fuses into a closed solid.
    #[test]
    fn merged_corner_filament_passes_the_gates() {
        let a = 25.0_f64.to_radians();
        let verts = [
            p(8.0, 0.275),
            p(0.275, 0.275),
            p(0.275, 0.646),
            p(0.275 + 9.0 * a.cos(), 0.646 + 9.0 * a.sin()),
        ];
        let f = Filament::try_serpentine(&verts, 0.275, 0.29).unwrap();
        let mut brep = BRep::new();
        let (shell, rep) = crate::fuse_filaments_with(&mut brep, &[f], &Default::default());
        assert_eq!(rep.aborted, 0);
        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let r = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(r.passed(), "wire {:?} manifold {:?}", r.wire_violations, r.manifold_violations);
    }
}
