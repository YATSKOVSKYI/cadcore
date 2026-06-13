//! Composite filament surface — a whole swept filament (leg cylinders +
//! elbow torus fillets) as ONE [`Projectable`].
//!
//! Trim curves at structure corners span junctions: they run off a leg, over
//! an elbow, onto the next leg.  Intersecting such filaments pairwise with
//! single analytic surfaces would force the caller to stitch curve pieces
//! across junctions; projecting onto the *composite* lets the tracer walk the
//! whole curve in one go, switching member surfaces implicitly.
//!
//! Member selection is by smallest projection displacement among the members
//! whose *parametric span* (axial range for legs, arc range for elbows)
//! contains the projection — with a small overhang so curves crossing a
//! junction never fall into a gap between members.

use cadcore_geom::{CylSurf, TorusSurf};
use cadcore_math::{Point3, Vec3};

use super::refine::{AnalyticSurface, Projectable};

/// A leg: finite cylinder span.
#[derive(Debug, Clone, Copy)]
pub struct LegSpan {
    /// Carrier cylinder (origin = axial 0).
    pub surf: CylSurf,
    /// Axial length (mm).
    pub length: f64,
}

/// An elbow: torus fillet over a signed arc range.
#[derive(Debug, Clone, Copy)]
pub struct ElbowSpan {
    /// Carrier torus.
    pub surf: TorusSurf,
    /// Arc start angle (radians, in the torus frame).
    pub theta_lo: f64,
    /// Arc end angle (signed sweep from `theta_lo`).
    pub theta_hi: f64,
}

/// A whole filament as one projectable surface.
#[derive(Debug, Clone, Default)]
pub struct CompositeTube {
    /// Straight legs.
    pub legs: Vec<LegSpan>,
    /// Elbow fillets.
    pub elbows: Vec<ElbowSpan>,
    /// Span overhang tolerance (mm along the axis / arc).  Members of a
    /// well-formed filament MEET exactly at junction circles, so this only
    /// absorbs numerical noise.  Keep it SMALL: a generous overhang extends
    /// legs beyond their true span into phantom surface inside the elbow
    /// region, and trim machinery downstream rightly rejects points there.
    /// Default `1e-3`.
    pub overhang: f64,
}

impl CompositeTube {
    /// New empty composite with the default overhang.
    pub fn new() -> Self {
        Self {
            legs: Vec::new(),
            elbows: Vec::new(),
            overhang: 1e-3,
        }
    }

    /// Add a leg.
    pub fn with_leg(mut self, surf: CylSurf, length: f64) -> Self {
        self.legs.push(LegSpan { surf, length });
        self
    }

    /// Add an elbow.
    pub fn with_elbow(mut self, surf: TorusSurf, theta_lo: f64, theta_hi: f64) -> Self {
        self.elbows.push(ElbowSpan {
            surf,
            theta_lo,
            theta_hi,
        });
        self
    }

    /// Elbow whose arc range is derived from two centreline junction points
    /// (angles computed in the torus's OWN frame; shortest signed sweep).
    pub fn with_elbow_between(self, surf: TorusSurf, start: Point3, end: Point3) -> Self {
        let ang = |p: Point3| {
            let w = p - surf.frame.origin;
            surf.frame.y.dot_vec(w).atan2(surf.frame.x.dot_vec(w))
        };
        let lo = ang(start);
        let mut hi = ang(end);
        let mut d = hi - lo;
        while d > std::f64::consts::PI {
            d -= std::f64::consts::TAU;
        }
        while d < -std::f64::consts::PI {
            d += std::f64::consts::TAU;
        }
        hi = lo + d;
        self.with_elbow(surf, lo, hi)
    }

    /// Nearest in-span member projection: `(projected point, member surface)`.
    fn nearest(&self, p: Point3) -> Option<(Point3, AnalyticSurface)> {
        self.nearest_ref(p).map(|(q, s, _)| (q, s))
    }

    /// [`Self::nearest`] plus the member reference.
    fn nearest_ref(&self, p: Point3) -> Option<(Point3, AnalyticSurface, MemberRef)> {
        let mut best: Option<(f64, Point3, AnalyticSurface, MemberRef)> = None;
        for (li, leg) in self.legs.iter().enumerate() {
            let s = leg.surf;
            let w = p - s.frame.origin;
            let ax = s.axis().dot_vec(w);
            if ax < -self.overhang || ax > leg.length + self.overhang {
                continue;
            }
            let surf = AnalyticSurface::Cylinder(s);
            let q = surf.project(p);
            let d = (q - p).length();
            if best.as_ref().map_or(true, |(bd, _, _, _)| d < *bd) {
                best = Some((d, q, surf, MemberRef::Leg(li)));
            }
        }
        for (ei, elbow) in self.elbows.iter().enumerate() {
            let t = elbow.surf;
            let w = p - t.frame.origin;
            let ang = t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w));
            // angular containment with overhang (overhang in radians via
            // major radius)
            let eps = self.overhang / t.major_radius.max(1e-9);
            let sweep = elbow.theta_hi - elbow.theta_lo;
            let rel = (ang - elbow.theta_lo) * sweep.signum();
            let inside = rel > -eps && rel < sweep.abs() + eps;
            // also allow the wrapped branch
            let rel_w = (ang - elbow.theta_lo + std::f64::consts::TAU * sweep.signum())
                * sweep.signum();
            let inside_w = rel_w > -eps && rel_w < sweep.abs() + eps;
            if !inside && !inside_w {
                continue;
            }
            let surf = AnalyticSurface::Torus(t);
            let q = surf.project(p);
            let d = (q - p).length();
            if best.as_ref().map_or(true, |(bd, _, _, _)| d < *bd) {
                best = Some((d, q, surf, MemberRef::Elbow(ei)));
            }
        }
        if best.is_none() {
            // No member's span contains the point — fall back to the global
            // nearest member.  A projection target must ALWAYS exist:
            // returning the input unchanged would masquerade as
            // "on-surface" (distance 0) and silently corrupt callers.
            for (li, leg) in self.legs.iter().enumerate() {
                let surf = AnalyticSurface::Cylinder(leg.surf);
                let q = surf.project(p);
                let d = (q - p).length();
                if best.as_ref().map_or(true, |(bd, _, _, _)| d < *bd) {
                    best = Some((d, q, surf, MemberRef::Leg(li)));
                }
            }
            for (ei, elbow) in self.elbows.iter().enumerate() {
                let surf = AnalyticSurface::Torus(elbow.surf);
                let q = surf.project(p);
                let d = (q - p).length();
                if best.as_ref().map_or(true, |(bd, _, _, _)| d < *bd) {
                    best = Some((d, q, surf, MemberRef::Elbow(ei)));
                }
            }
        }
        best.map(|(_, q, s, m)| (q, s, m))
    }
}

/// Reference to a member of a [`CompositeTube`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum MemberRef {
    /// Index into `legs`.
    Leg(usize),
    /// Index into `elbows`.
    Elbow(usize),
}

/// A junction between two consecutive members: the shared boundary circle.
#[derive(Debug, Clone, Copy)]
pub struct Junction {
    /// First member.
    pub a: MemberRef,
    /// Second member.
    pub b: MemberRef,
    /// The shared circle (centre on the centreline, normal = tangent there).
    pub circle: cadcore_geom::Circle3,
}

impl CompositeTube {
    /// Member whose surface is nearest to `p` (same selection rule as
    /// projection).
    pub fn member_of(&self, p: Point3) -> Option<MemberRef> {
        self.nearest_ref(p).map(|(_, _, m)| m)
    }

    /// Is `p` strictly inside the tube volume?
    pub fn contains(&self, p: Point3) -> bool {
        for leg in &self.legs {
            let w = p - leg.surf.frame.origin;
            let ax = leg.surf.axis().dot_vec(w);
            if ax < 0.0 || ax > leg.length {
                continue;
            }
            let rad = (w - leg.surf.axis().as_vec() * ax).length();
            if rad < leg.surf.radius {
                return true;
            }
        }
        for e in &self.elbows {
            let t = &e.surf;
            let w = p - t.frame.origin;
            let ang = t.frame.y.dot_vec(w).atan2(t.frame.x.dot_vec(w));
            let eps = self.overhang / t.major_radius.max(1e-9);
            let sweep = e.theta_hi - e.theta_lo;
            let rel = (ang - e.theta_lo) * sweep.signum();
            let rel_w = (ang - e.theta_lo + std::f64::consts::TAU * sweep.signum())
                * sweep.signum();
            let inside_arc = (rel > -eps && rel < sweep.abs() + eps)
                || (rel_w > -eps && rel_w < sweep.abs() + eps);
            if !inside_arc {
                continue;
            }
            let h = t.frame.z.dot_vec(w);
            let a = (w - t.frame.z.as_vec() * h).length();
            if a < 1e-12 {
                continue;
            }
            let d_minor = ((a - t.major_radius).powi(2) + h * h).sqrt();
            if d_minor < t.minor_radius {
                return true;
            }
        }
        false
    }

    /// All member-to-member junction circles, found by matching member
    /// centreline endpoints within `tol`.
    pub fn junctions(&self, tol: f64) -> Vec<Junction> {
        // endpoints: (point, tangent, member)
        let mut ends: Vec<(Point3, Vec3, MemberRef)> = Vec::new();
        for (i, leg) in self.legs.iter().enumerate() {
            let d = leg.surf.axis().as_vec();
            ends.push((leg.surf.frame.origin, d, MemberRef::Leg(i)));
            ends.push((leg.surf.frame.origin + d * leg.length, d, MemberRef::Leg(i)));
        }
        for (i, e) in self.elbows.iter().enumerate() {
            for th in [e.theta_lo, e.theta_hi] {
                let dir = e.surf.frame.x * th.cos() + e.surf.frame.y * th.sin();
                let p = e.surf.frame.origin + dir * e.surf.major_radius;
                // ring tangent at the angle
                let tan = (e.surf.frame.y.as_vec() * th.cos()
                    - e.surf.frame.x.as_vec() * th.sin())
                    * (e.theta_hi - e.theta_lo).signum();
                ends.push((p, tan, MemberRef::Elbow(i)));
            }
        }
        let mut out = Vec::new();
        for i in 0..ends.len() {
            for j in (i + 1)..ends.len() {
                if ends[i].2 == ends[j].2 {
                    continue;
                }
                if (ends[i].0 - ends[j].0).length() > tol {
                    continue;
                }
                let radius = match ends[i].2 {
                    MemberRef::Leg(k) => self.legs[k].surf.radius,
                    MemberRef::Elbow(k) => self.elbows[k].surf.minor_radius,
                };
                let normal = match cadcore_math::UnitVec3::try_from_vec(ends[i].1) {
                    Some(n) => n,
                    None => continue,
                };
                out.push(Junction {
                    a: ends[i].2,
                    b: ends[j].2,
                    circle: cadcore_geom::Circle3::new(ends[i].0, normal, radius),
                });
            }
        }
        out
    }

    /// Centreline samples `(point, tube_radius)` for broad-phase seeding.
    ///
    /// Sampling is adaptive: at most `max_gap·radius` between consecutive
    /// samples, so a shallow crossing (centrelines approaching to barely
    /// under `r_a + r_b`) can never slip between two samples.
    pub fn centreline_samples(&self, max_gap: f64) -> Vec<(Point3, f64)> {
        let mut out = Vec::new();
        for leg in &self.legs {
            let gap = (max_gap * leg.surf.radius).max(1e-6);
            let n = ((leg.length / gap).ceil() as usize).max(2);
            for k in 0..=n {
                let t = leg.length * k as f64 / n as f64;
                out.push((
                    leg.surf.frame.origin + leg.surf.axis().as_vec() * t,
                    leg.surf.radius,
                ));
            }
        }
        for e in &self.elbows {
            let arc_len = (e.theta_hi - e.theta_lo).abs() * e.surf.major_radius;
            let gap = (max_gap * e.surf.minor_radius).max(1e-6);
            let n = ((arc_len / gap).ceil() as usize).max(2);
            for k in 0..=n {
                let th = e.theta_lo + (e.theta_hi - e.theta_lo) * k as f64 / n as f64;
                let dir = e.surf.frame.x * th.cos() + e.surf.frame.y * th.sin();
                out.push((
                    e.surf.frame.origin + dir * e.surf.major_radius,
                    e.surf.minor_radius,
                ));
            }
        }
        out
    }
}

/// Broad-phase seeds for the intersection of two composite tubes: midpoints
/// between close centreline sample pairs (surfaces overlap where centreline
/// distance < r_a + r_b).
pub fn seeds_between(a: &CompositeTube, b: &CompositeTube, max_gap: f64) -> Vec<Point3> {
    let sa = a.centreline_samples(max_gap);
    let sb = b.centreline_samples(max_gap);
    let mut seeds: Vec<Point3> = Vec::new();
    for &(pa, ra) in &sa {
        for &(pb, rb) in &sb {
            let d = (pb - pa).length();
            if d >= ra + rb || d < 1e-9 {
                continue;
            }
            let mid = pa + (pb - pa) * 0.5;
            // cluster: keep seeds at least half a radius apart
            let min_sep = 0.5 * ra.min(rb);
            if seeds.iter().all(|&s| (s - mid).length() > min_sep) {
                seeds.push(mid);
            }
        }
    }
    seeds
}

/// All distinct closed intersection loops between two composite tubes:
/// broad-phase seeds → refine → trace → dedup (a loop found from several
/// seeds is kept once).
pub fn closed_loops_between(
    a: &CompositeTube,
    b: &CompositeTube,
    opts: &crate::geom::intersect::TraceOptions,
) -> Vec<crate::geom::intersect::TracedCurve> {
    use crate::geom::intersect::{intersection_point_near, trace_dyn};
    let dbg = std::env::var("CADCORE_DUMP_TRACER").is_ok();
    let mut loops: Vec<crate::geom::intersect::TracedCurve> = Vec::new();
    let (mut n_seed, mut n_refine_fail, mut n_dup, mut ends) =
        (0usize, 0usize, 0usize, std::collections::BTreeMap::<String, usize>::new());
    for seed in seeds_between(a, b, 0.4) {
        n_seed += 1;
        let Some(p) = intersection_point_near(a, b, seed, opts.tol) else {
            n_refine_fail += 1;
            continue;
        };
        // already on a known loop?
        let dup = loops.iter().any(|l| {
            l.points
                .iter()
                .any(|&q| (q - p).length() < opts.step * 1.5)
        });
        if dup {
            n_dup += 1;
            continue;
        }
        if let Some(c) = trace_dyn(a, b, p, opts, |_| true) {
            *ends.entry(format!("{:?}", c.end)).or_insert(0) += 1;
            if c.closed() {
                loops.push(c);
            } else if dbg && ends.get("Diverged").is_some_and(|&n| n <= 3) {
                let q = *c.points.last().unwrap();
                eprintln!(
                    "[tracer:div] n={} start=({:.3},{:.3},{:.3}) last=({:.3},{:.3},{:.3}) dA={:.2e} dB={:.2e}",
                    c.points.len(), p.x, p.y, p.z, q.x, q.y, q.z,
                    a.distance_p(q), b.distance_p(q)
                );
            }
        }
    }
    if dbg {
        eprintln!(
            "[tracer] seeds={n_seed} refine_fail={n_refine_fail} dup={n_dup} ends={ends:?} loops={}",
            loops.len()
        );
    }
    loops
}

impl Projectable for CompositeTube {
    fn project_p(&self, p: Point3) -> Point3 {
        self.nearest(p).map_or(p, |(q, _)| q)
    }

    fn normal_p(&self, p: Point3) -> Option<Vec3> {
        let (_, s) = self.nearest(p)?;
        s.normal_at(p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::intersect::{trace_closed_loop_dyn, TraceOptions};
    use cadcore_math::UnitVec3;

    /// L-shaped filament in the XY plane: leg along +X, 90° elbow, leg along
    /// +Y — the serpentine corner building block.
    fn l_filament(z: f64, r: f64) -> CompositeTube {
        // leg A: from (-2,0,z) to (0,0,z) along +X
        let leg_a = CylSurf::new(Point3::new(-2.0, 0.0, z), UnitVec3::X, r);
        // elbow: torus centre (0, 0.5, z), axis Z, major 0.5: sweeps from the
        // -Y side (angle -90° = touching leg A at (0,0,z)) to the +X side
        // (angle 0° = touching leg B at (0.5,0.5,z))
        let t = TorusSurf::new(Point3::new(0.0, 0.5, z), UnitVec3::Z, 0.5, r);
        // leg B: from (0.5,0.5,z) to (0.5,2.5,z) along +Y
        let leg_b = CylSurf::new(Point3::new(0.5, 0.5, z), UnitVec3::Y, r);
        CompositeTube::new()
            .with_leg(leg_a, 2.0)
            .with_leg(leg_b, 2.0)
            .with_elbow_between(t, Point3::new(0.0, 0.0, z), Point3::new(0.5, 0.5, z))
    }

    #[test]
    fn projects_onto_nearest_member() {
        let f = l_filament(0.0, 0.275);
        // near leg A
        let q = f.project_p(Point3::new(-1.0, 0.2, 0.1));
        assert!(f.distance_p(q) < 1e-12);
        assert!((q.x - -1.0).abs() < 1e-9, "stays at the leg-A station");
        // near the elbow apex (outer side of the bend)
        let q2 = f.project_p(Point3::new(-0.25, 0.25, 0.0));
        assert!(f.distance_p(q2) < 1e-12);
    }

    #[test]
    fn loops_found_without_manual_seeds() {
        // same U-on-U corner, but the engine finds the loop itself
        let a = l_filament(0.0, 0.275);
        let b = stacked_rotated_l();
        let loops = closed_loops_between(&a, &b, &TraceOptions::default());
        assert_eq!(loops.len(), 1, "exactly one corner loop: {}", loops.len());
        assert!(loops[0].points.len() > 10);
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

    #[test]
    fn corner_bite_loop_spans_junction() {
        // Two L-filaments in layers z=0 and z=0.35: the upper one's corner
        // sits on the lower one's corner — the U-on-U configuration.  The
        // intersection loop crosses leg↔elbow junctions on both.
        let a = l_filament(0.0, 0.275);
        let b = {
            // rotate the L 90°: legs along Y then X, elbow at the same corner
            let leg_a = CylSurf::new(Point3::new(0.0, -2.0, 0.35), UnitVec3::Y, 0.275);
            // turn centre right of the +Y heading at (0,0): (0.5, 0); the
            // 90° clockwise sweep exits at (0.5, 0.5) heading +X
            let t = TorusSurf::new(Point3::new(0.5, 0.0, 0.35), UnitVec3::Z, 0.5, 0.275);
            let leg_b = CylSurf::new(Point3::new(0.5, 0.5, 0.35), UnitVec3::X, 0.275);
            CompositeTube::new()
                .with_leg(leg_a, 2.0)
                .with_leg(leg_b, 2.0)
                .with_elbow_between(
                    t,
                    Point3::new(0.0, 0.0, 0.35),
                    Point3::new(0.5, 0.5, 0.35),
                )
        };
        let c = trace_closed_loop_dyn(
            &a,
            &b,
            Point3::new(0.0, 0.0, 0.18),
            &TraceOptions::default(),
        )
        .expect("corner loop closes");
        assert!(c.points.len() > 10);
        let dev = c
            .points
            .iter()
            .map(|&p| a.distance_p(p).max(b.distance_p(p)))
            .fold(0.0, f64::max);
        assert!(dev < 1e-9, "loop exact on both composites: {dev:.2e}");
    }
}
