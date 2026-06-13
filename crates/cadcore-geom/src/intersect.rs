//! Surface–surface intersection (SSI) curves for the analytic Boolean union.
//!
//! v1 scope (DIW woodpile scaffolds): the only surface pairs that meet are
//!
//! * **cylinder ∩ cylinder** — perpendicular axes, equal radius, offset in the
//!   stacking direction.  This is the woodpile crossing curve.  For perpendicular
//!   axes the intersection has a closed form (no marching/Newton needed): every
//!   returned point lies on *both* surfaces to `f64` precision.
//! * **cylinder ∩ plane** — a filament end-cap plane cutting a crossing
//!   cylinder.  Analytic conic: circle / ellipse / line-pair.
//!
//! Sphere and torus pairs are intentionally **not** handled (out of v1 scope).

use cadcore_math::{Point3, UnitVec3, Vec3};

use crate::{Circle3, CylSurf, Ellipse3, Line3, Plane3, TorusSurf};

const TWO_PI: f64 = std::f64::consts::PI * 2.0;

/// A sampled intersection curve (polyline).  `closed` marks a loop whose last
/// point connects back to the first.
#[derive(Clone, Debug)]
pub struct IntersectionPolyline {
    /// Ordered sample points, all lying on both surfaces.
    pub points: Vec<Point3>,
    /// `true` when the polyline is a closed loop.
    pub closed: bool,
}

/// Analytic conic from a cylinder ∩ plane intersection.
#[derive(Clone, Debug)]
pub enum CylPlaneCurve {
    /// Plane perpendicular to the cylinder axis.
    Circle(Circle3),
    /// Plane oblique to the axis (not parallel).
    Ellipse(Ellipse3),
    /// Plane parallel to the axis: 0, 1 (tangent) or 2 ruling lines.
    Lines(Vec<Line3>),
    /// No intersection.
    None,
}

// ── cylinder ∩ cylinder ───────────────────────────────────────────────────────

/// Intersection of two cylinders with **perpendicular axes** and **equal radius**.
///
/// Returns the closed intersection loop(s):
/// * empty — surfaces do not meet, or inputs fall outside v1 scope
///   (unequal radius / non-perpendicular axes);
/// * one loop — the usual transversal woodpile crossing (partial overlap);
/// * two loops — full mutual penetration (axes nearly intersecting).
///
/// `samples` controls angular resolution around cylinder `a`.
///
/// # Exactness
///
/// A point is generated from an angle `θ` on `a` (so it lies on `a` exactly),
/// then the axial parameter `t` is the exact root of `b`'s quadratic — so the
/// point lies on `b` exactly too.  Because the axes are perpendicular
/// (`da · db = 0`), `b`'s constraint is quadratic in `t`:
/// `t² + 2(w₀·da)·t + (|w₀|² − (w₀·db)² − r²) = 0`.
pub fn cyl_cyl_intersection(a: &CylSurf, b: &CylSurf, samples: usize) -> Vec<IntersectionPolyline> {
    if (a.radius - b.radius).abs() > 1e-9 {
        return vec![]; // unequal radius — out of v1 scope
    }
    let da = a.axis();
    let db = b.axis();
    let dab = da.dot(db);
    // General quadratic below handles ANY axis angle; only reject near-parallel
    // (where the intersection degenerates). NOTE: real scaffold legs are only
    // perpendicular to ~1e-5 (float noise), so a strict perpendicular test would
    // wrongly skip most crossings — leaving untrimmed self-intersecting overlaps.
    let aq = 1.0 - dab * dab;
    // Reject clearly non-perpendicular pairs. Orthogonal/diamond scaffold layers
    // cross at ~90° (|da·db| ≈ noise ~1e-5); same-orientation legs are parallel
    // (|da·db| ≈ 1). The 0.9 cut keeps every real crossing and rejects parallel
    // legs + shallow-angle connector pairs (whose general-quadratic root lands
    // far off-segment and would otherwise be spurious).
    if dab.abs() > 0.9 {
        return vec![];
    }
    let r = a.radius;
    let oa = a.frame.origin;
    let ax = a.frame.x;
    let ay = a.frame.y;
    let ob = b.frame.origin;

    // For a point `P = c(θ) + da·t` on cylinder A, the constraint that it also
    // lies on cylinder B, `|w|² − (w·db)² = r²` with `w = P − ob`, is a quadratic
    // in the axial parameter `t` (general — no perpendicularity assumed):
    //   aq·t² + bq·t + cq = 0,
    //   aq = 1 − (da·db)²,  bq = 2[(w₀·da) − (w₀·db)(da·db)],
    //   cq = |w₀|² − (w₀·db)² − r².
    let coeffs = |theta: f64| -> (f64, f64, Point3) {
        let (s, c) = theta.sin_cos();
        let cpt = oa + ax * (c * r) + ay * (s * r);
        let w0 = cpt - ob;
        let wda = da.dot_vec(w0);
        let wdb = db.dot_vec(w0);
        let bq = 2.0 * (wda - wdb * dab);
        let cq = w0.length_sq() - wdb * wdb - r * r;
        (bq, cq, cpt)
    };
    // Discriminant of the quadratic (≥0 ⇒ the ray meets B).
    let delta = |theta: f64| -> f64 {
        let (bq, cq, _) = coeffs(theta);
        bq * bq - 4.0 * aq * cq
    };
    let branch_point = |theta: f64, plus: bool| -> Point3 {
        let (bq, cq, cpt) = coeffs(theta);
        let d = (bq * bq - 4.0 * aq * cq).max(0.0);
        let sq = d.sqrt();
        let t = if plus {
            (-bq + sq) / (2.0 * aq)
        } else {
            (-bq - sq) / (2.0 * aq)
        };
        cpt + da * t
    };

    let n = samples.max(48);
    let valid: Vec<bool> = (0..n)
        .map(|i| delta(i as f64 / n as f64 * TWO_PI) >= 0.0)
        .collect();
    let n_valid = valid.iter().filter(|&&v| v).count();

    if n_valid == 0 {
        return vec![];
    }

    // Full mutual penetration: Δ ≥ 0 everywhere → two separate loops (± branches).
    if n_valid == n {
        let mut loops = Vec::new();
        for plus in [true, false] {
            let pts: Vec<Point3> = (0..n)
                .map(|i| branch_point(i as f64 / n as f64 * TWO_PI, plus))
                .collect();
            loops.push(IntersectionPolyline {
                points: pts,
                closed: true,
            });
        }
        return loops;
    }

    // Partial overlap: assemble loops from the maximal Δ ≥ 0 arcs.  Each arc
    // [θ_lo, θ_hi] (bracketed by Δ = 0) becomes one loop: +branch θ_lo→θ_hi
    // then −branch θ_hi→θ_lo.
    let theta_of = |i: usize| i as f64 / n as f64 * TWO_PI;
    // Refine a Δ=0 crossing between an invalid sample θ_i and a valid θ_v.
    let refine_zero = |mut lo: f64, mut hi: f64| -> f64 {
        // invariant: delta(lo) < 0, delta(hi) >= 0
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if delta(mid) >= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        hi
    };

    // Find a rotation start at an invalid index so runs don't wrap.
    let start = (0..n).find(|&i| !valid[i]).unwrap();
    let mut loops = Vec::new();
    let mut i = 0usize;
    while i < n {
        let idx = (start + i) % n;
        if !valid[idx] {
            i += 1;
            continue;
        }
        // Start of a valid run at `idx`.
        let run_begin = i;
        let mut j = i;
        while j < n && valid[(start + j) % n] {
            j += 1;
        }
        let run_end = j - 1; // last valid offset
                             // θ window, refined at both ends using the bracketing invalid samples.
        let prev_invalid = theta_of((start + run_begin + n - 1) % n);
        let first_valid = theta_of((start + run_begin) % n);
        let last_valid = theta_of((start + run_end) % n);
        let next_invalid = theta_of((start + run_end + 1) % n);
        // Unwrap so the window is monotonic increasing.
        let mut theta_lo = refine_zero_wrap(prev_invalid, first_valid, &refine_zero);
        let mut theta_hi = refine_zero_wrap(next_invalid, last_valid, &refine_zero);
        // Ensure theta_lo < theta_hi on a continuous branch.
        if theta_hi <= theta_lo {
            theta_hi += TWO_PI;
        }
        let _ = (&mut theta_lo, &mut theta_hi);

        let arc = theta_hi - theta_lo;
        let m = ((arc / TWO_PI) * n as f64).ceil() as usize;
        let m = m.max(6);
        let mut pts = Vec::with_capacity(2 * m);
        // +branch θ_lo → θ_hi (inclusive)
        for k in 0..=m {
            let theta = theta_lo + arc * (k as f64 / m as f64);
            pts.push(branch_point(theta, true));
        }
        // −branch θ_hi → θ_lo (exclusive of both ends — they equal +branch ends)
        for k in (1..m).rev() {
            let theta = theta_lo + arc * (k as f64 / m as f64);
            pts.push(branch_point(theta, false));
        }
        loops.push(IntersectionPolyline {
            points: pts,
            closed: true,
        });

        i = j;
    }
    loops
}

/// Refine a Δ=0 crossing where the two θ samples may straddle the 0/2π seam.
fn refine_zero_wrap(
    theta_invalid: f64,
    theta_valid: f64,
    refine: &impl Fn(f64, f64) -> f64,
) -> f64 {
    // Move onto a common branch: pick the representation of theta_invalid
    // closest to theta_valid.
    let mut ti = theta_invalid;
    while ti - theta_valid > std::f64::consts::PI {
        ti -= TWO_PI;
    }
    while theta_valid - ti > std::f64::consts::PI {
        ti += TWO_PI;
    }
    // `refine` brackets by sign (invalid θ, valid θ) — numeric order is irrelevant.
    refine(ti, theta_valid)
}

// ── torus ∩ cylinder ──────────────────────────────────────────────────────────

/// Intersection of a torus-fillet face (arc `θ ∈ [theta_lo, theta_hi]`, full
/// minor circle `φ ∈ [0, 2π)`) with a cylinder.  Returns the closed loop(s)
/// where the cylinder pokes through the torus tube (the woodpile corner case:
/// a curved U-turn crossing a straight filament of the adjacent layer).
///
/// Traced as the zero contour of `g(θ,φ) = signed_dist(torus(θ,φ), cyl)` via
/// marching squares on a `(θ,φ)` grid, then chained into loops.  Sampled to a
/// polyline (the curve is high-order; no closed form).
pub fn torus_cyl_intersection(
    torus: &TorusSurf,
    cyl: &CylSurf,
    theta_lo: f64,
    theta_hi: f64,
    samples: usize,
) -> Vec<IntersectionPolyline> {
    let m = samples.max(64); // θ steps (along the arc)
    let n = samples.max(48); // φ steps (around the tube, periodic)
    let theta_at = |i: usize| theta_lo + (theta_hi - theta_lo) * (i as f64 / m as f64);
    let phi_at = |j: usize| (j as f64 / n as f64) * TWO_PI;
    let g = |i: usize, j: usize| -> f64 {
        cyl_surface_distance(cyl, torus.point_at(theta_at(i), phi_at(j)))
    };

    // Marching squares → contour segments in (θ,φ).
    let frac = |va: f64, vb: f64| va / (va - vb);
    let mut segs: Vec<([f64; 2], [f64; 2])> = Vec::new();
    for i in 0..m {
        for j in 0..n {
            let (t0, t1) = (theta_at(i), theta_at(i + 1));
            let (p0, p1) = (phi_at(j), phi_at(j + 1));
            let g00 = g(i, j);
            let g10 = g(i + 1, j);
            let g11 = g(i + 1, j + 1);
            let g01 = g(i, j + 1);
            let mut pts: Vec<[f64; 2]> = Vec::with_capacity(2);
            if (g00 < 0.0) != (g10 < 0.0) {
                pts.push([t0 + (t1 - t0) * frac(g00, g10), p0]);
            }
            if (g10 < 0.0) != (g11 < 0.0) {
                pts.push([t1, p0 + (p1 - p0) * frac(g10, g11)]);
            }
            if (g11 < 0.0) != (g01 < 0.0) {
                pts.push([t1 + (t0 - t1) * frac(g11, g01), p1]);
            }
            if (g01 < 0.0) != (g00 < 0.0) {
                pts.push([t0, p1 + (p0 - p1) * frac(g01, g00)]);
            }
            if pts.len() == 2 {
                segs.push((pts[0], pts[1]));
            }
        }
    }
    if segs.is_empty() {
        return vec![];
    }

    // Chain segments into closed polylines by endpoint proximity.
    let tol = ((theta_hi - theta_lo).abs() / m as f64).max(TWO_PI / n as f64) * 0.75;
    let mut used = vec![false; segs.len()];
    let close = |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).abs() < tol && (a[1] - b[1]).abs() < tol;
    let mut loops3d = Vec::new();
    for start in 0..segs.len() {
        if used[start] {
            continue;
        }
        used[start] = true;
        let mut chain = vec![segs[start].0, segs[start].1];
        loop {
            let tail = *chain.last().unwrap();
            let mut found = None;
            for (k, seg) in segs.iter().enumerate() {
                if used[k] {
                    continue;
                }
                if close(tail, seg.0) {
                    found = Some((k, seg.1));
                    break;
                } else if close(tail, seg.1) {
                    found = Some((k, seg.0));
                    break;
                }
            }
            match found {
                Some((k, next)) => {
                    used[k] = true;
                    chain.push(next);
                }
                None => break,
            }
        }
        if chain.len() >= 6 {
            let pts: Vec<Point3> = chain.iter().map(|&[t, p]| torus.point_at(t, p)).collect();
            loops3d.push(IntersectionPolyline {
                points: pts,
                closed: true,
            });
        }
    }
    loops3d
}

// ── cylinder ∩ plane ──────────────────────────────────────────────────────────

/// Intersection of a cylinder and a plane as an analytic conic.
pub fn cyl_plane_intersection(cyl: &CylSurf, plane: &Plane3) -> CylPlaneCurve {
    let d = cyl.axis();
    let n = plane.normal();
    let r = cyl.radius;
    let cos = n.dot(d).abs();

    if cos < 1e-6 {
        // Plane parallel to the axis → ruling lines where |dist(axis, plane)| ≤ r.
        let dist = plane.signed_distance(cyl.frame.origin); // along n
        if dist.abs() > r + 1e-12 {
            return CylPlaneCurve::None;
        }
        // In-plane direction perpendicular to the axis.
        let perp = match UnitVec3::try_from_vec(n.cross(d)) {
            Some(u) => u,
            None => return CylPlaneCurve::None,
        };
        // Foot of the axis projected onto the plane.
        let foot = cyl.frame.origin - n.as_vec() * dist;
        let half = (r * r - dist * dist).max(0.0).sqrt();
        let mut lines = Vec::new();
        if half < 1e-12 {
            lines.push(Line3::new(foot, d)); // tangent
        } else {
            lines.push(Line3::new(foot + perp.as_vec() * half, d));
            lines.push(Line3::new(foot - perp.as_vec() * half, d));
        }
        return CylPlaneCurve::Lines(lines);
    }

    // Centre = axis ∩ plane.
    let s0 = plane.signed_distance(cyl.frame.origin);
    let denom = n.dot_vec(d.as_vec());
    let t = -s0 / denom;
    let center = cyl.frame.origin + d * t;

    if (cos - 1.0).abs() < 1e-9 {
        return CylPlaneCurve::Circle(Circle3::new(center, n, r));
    }

    // Ellipse: semi-minor = r along (n × d) (in plane, ⟂ to axis projection);
    // semi-major = r / cos along the axis projected into the plane.
    let minor_dir = match UnitVec3::try_from_vec(n.cross(d)) {
        Some(u) => u,
        None => return CylPlaneCurve::None,
    };
    let major_dir = match UnitVec3::try_from_vec(n.as_vec().cross(minor_dir.as_vec())) {
        Some(u) => u,
        None => return CylPlaneCurve::None,
    };
    CylPlaneCurve::Ellipse(Ellipse3::new(center, n, major_dir, r / cos, r))
}

/// Sample a [`CylPlaneCurve`] ellipse/circle into a closed polyline.
pub fn sample_cyl_plane(curve: &CylPlaneCurve, samples: usize) -> Option<IntersectionPolyline> {
    let n = samples.max(24);
    match curve {
        CylPlaneCurve::Circle(c) => Some(IntersectionPolyline {
            points: (0..n)
                .map(|i| c.point_at(i as f64 / n as f64 * TWO_PI))
                .collect(),
            closed: true,
        }),
        CylPlaneCurve::Ellipse(e) => Some(IntersectionPolyline {
            points: (0..n)
                .map(|i| e.point_at(i as f64 / n as f64 * TWO_PI))
                .collect(),
            closed: true,
        }),
        _ => None,
    }
}

/// Intersection of two **parallel** overlapping cylinders: the two ruling
/// lines along the axis direction (the "tubes lying on each other" case the
/// transversal `cyl_cyl_intersection` cannot represent — its intersection
/// curve is open, not a closed loop).
///
/// Returns an empty vec when the axes are not parallel, the surfaces do not
/// cross (separated or one inside the other), or they are tangent (a single
/// degenerate line — at tangency the material is continuous and there is no
/// transversal curve to imprint).
pub fn cyl_cyl_parallel_lines(a: &CylSurf, b: &CylSurf) -> Vec<Line3> {
    let da = a.axis();
    let db = b.axis();
    if da.dot(db).abs() < 1.0 - 1e-9 {
        return vec![];
    }
    // Perpendicular offset between the two axes.
    let w = b.frame.origin - a.frame.origin;
    let w_perp: Vec3 = w - da.as_vec() * da.dot_vec(w);
    let d = w_perp.length();
    let (ra, rb) = (a.radius, b.radius);
    let tangent_tol = 1e-9;
    if d < tangent_tol || d >= ra + rb - tangent_tol || d <= (ra - rb).abs() + tangent_tol {
        return vec![]; // coaxial, separated/tangent, or nested/tangent
    }
    // 2D circle-circle intersection in the cross-section plane.
    let Some(e) = UnitVec3::try_from_vec(w_perp) else {
        return vec![];
    };
    let f = da.cross(e); // in-plane, ⟂ e
    let along = (d * d + ra * ra - rb * rb) / (2.0 * d);
    let h = (ra * ra - along * along).max(0.0).sqrt();
    let foot = a.frame.origin + e.as_vec() * along;
    vec![
        Line3::new(foot + f * h, da),
        Line3::new(foot - f * h, da),
    ]
}

/// Distance from a point to a cylinder's surface (signed: + outside, − inside).
#[inline]
pub fn cyl_surface_distance(cyl: &CylSurf, p: Point3) -> f64 {
    let w = p - cyl.frame.origin;
    let axial = cyl.axis().dot_vec(w);
    let radial: Vec3 = w - cyl.axis().as_vec() * axial;
    radial.length() - cyl.radius
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_math::Point3;

    /// Two parallel equal-radius tubes stacked with overlap (HL < D): exactly
    /// two ruling lines, each on BOTH surfaces; tangent/separated → none.
    #[test]
    fn parallel_overlap_two_lines() {
        let r = 0.275;
        let a = CylSurf::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::X, r);
        let b = CylSurf::new(Point3::new(0.0, 0.0, 0.35), UnitVec3::X, r);
        let lines = cyl_cyl_parallel_lines(&a, &b);
        assert_eq!(lines.len(), 2);
        for ln in &lines {
            assert!(ln.direction.dot(UnitVec3::X).abs() > 1.0 - 1e-12);
            for t in [-1.0, 0.0, 2.0] {
                let p = ln.origin + ln.direction * t;
                assert!(cyl_surface_distance(&a, p).abs() < 1e-12, "on A");
                assert!(cyl_surface_distance(&b, p).abs() < 1e-12, "on B");
            }
        }
        // The two lines are symmetric about the plane containing both axes.
        assert!((lines[0].origin.y + lines[1].origin.y).abs() < 1e-12);

        // Tangent (d = 2r) and separated → no transversal lines.
        let t1 = CylSurf::new(Point3::new(0.0, 0.0, 2.0 * r), UnitVec3::X, r);
        assert!(cyl_cyl_parallel_lines(&a, &t1).is_empty(), "tangent");
        let t2 = CylSurf::new(Point3::new(0.0, 0.0, 1.0), UnitVec3::X, r);
        assert!(cyl_cyl_parallel_lines(&a, &t2).is_empty(), "separated");
        // Non-parallel input is rejected (the transversal path handles it).
        let t3 = CylSurf::new(Point3::new(0.0, 0.0, 0.35), UnitVec3::Y, r);
        assert!(cyl_cyl_parallel_lines(&a, &t3).is_empty(), "not parallel");
    }

    fn cyl_x(z: f64, r: f64) -> CylSurf {
        CylSurf::new(Point3::new(0.0, 0.0, z), UnitVec3::X, r)
    }
    fn cyl_y(z: f64, r: f64) -> CylSurf {
        CylSurf::new(Point3::new(0.0, 0.0, z), UnitVec3::Y, r)
    }

    /// Laterally-offset crossing (the real scaffold case: an X-leg at y=2 and a
    /// Y-leg at x=3 cross at (3,2)). The loop must still be found near (3,2).
    #[test]
    fn cyl_cyl_lateral_offset_crossing() {
        let r = 0.275;
        let a = CylSurf::new(Point3::new(0.0, 2.0, 0.0), UnitVec3::X, r); // X-leg at y=2,z=0
        let b = CylSurf::new(Point3::new(3.0, 0.0, 0.35), UnitVec3::Y, r); // Y-leg at x=3,z=0.35
        let loops = cyl_cyl_intersection(&a, &b, 96);
        assert_eq!(loops.len(), 1, "offset crossing must produce one loop");
        for p in &loops[0].points {
            assert!(cyl_surface_distance(&a, *p).abs() < 1e-7, "on A");
            assert!(cyl_surface_distance(&b, *p).abs() < 1e-7, "on B");
            // Loop sits at the crossing (x≈3 on A, y≈2 on B).
            assert!((p.x - 3.0).abs() < 0.4, "x near 3, got {}", p.x);
            assert!((p.y - 2.0).abs() < 0.4, "y near 2, got {}", p.y);
        }
    }

    /// HL=0.35 / D=0.55 woodpile crossing: one closed loop, all points on both.
    #[test]
    fn cyl_cyl_woodpile_crossing_single_loop() {
        let r = 0.275;
        let a = cyl_x(0.0, r);
        let b = cyl_y(0.35, r);
        let loops = cyl_cyl_intersection(&a, &b, 96);
        assert_eq!(loops.len(), 1, "partial overlap → one loop");
        let lp = &loops[0];
        assert!(lp.closed);
        assert!(lp.points.len() >= 12);
        for p in &lp.points {
            assert!(cyl_surface_distance(&a, *p).abs() < 1e-7, "on A");
            assert!(cyl_surface_distance(&b, *p).abs() < 1e-7, "on B");
        }
    }

    /// Axes intersecting (offset 0) with equal radius → two loops (ellipses).
    #[test]
    fn cyl_cyl_full_penetration_two_loops() {
        let r = 0.3;
        let a = cyl_x(0.0, r);
        let b = cyl_y(0.0, r);
        let loops = cyl_cyl_intersection(&a, &b, 96);
        assert_eq!(loops.len(), 2, "full penetration → two loops");
        for lp in &loops {
            for p in &lp.points {
                assert!(cyl_surface_distance(&a, *p).abs() < 1e-7);
                assert!(cyl_surface_distance(&b, *p).abs() < 1e-7);
            }
        }
    }

    /// Cylinders too far apart (no overlap) → empty.
    #[test]
    fn cyl_cyl_no_overlap_empty() {
        let r = 0.2;
        let a = cyl_x(0.0, r);
        let b = cyl_y(1.0, r); // 1.0 > 2r = 0.4
        assert!(cyl_cyl_intersection(&a, &b, 64).is_empty());
    }

    /// Out-of-scope guards: unequal radius / parallel axes → empty.
    #[test]
    fn cyl_cyl_out_of_scope_empty() {
        assert!(cyl_cyl_intersection(&cyl_x(0.0, 0.2), &cyl_y(0.1, 0.3), 64).is_empty());
        assert!(cyl_cyl_intersection(&cyl_x(0.0, 0.2), &cyl_x(0.1, 0.2), 64).is_empty());
    }

    /// Plane perpendicular to axis → circle; oblique → ellipse on the surface.
    #[test]
    fn cyl_plane_circle_and_ellipse() {
        let r = 0.275;
        let cyl = cyl_x(0.0, r);
        // Plane x = 1 (normal X) is perpendicular to the X-axis cylinder.
        let perp = Plane3::from_origin_normal(Point3::new(1.0, 0.0, 0.0), UnitVec3::X);
        match cyl_plane_intersection(&cyl, &perp) {
            CylPlaneCurve::Circle(c) => assert!((c.radius - r).abs() < 1e-12),
            other => panic!("expected circle, got {other:?}"),
        }
        // Oblique plane (normal in XZ at 45°).
        let nrm = UnitVec3::try_from_vec(Vec3::new(1.0, 0.0, 1.0)).unwrap();
        let obl = Plane3::from_origin_normal(Point3::new(0.5, 0.0, 0.0), nrm);
        match cyl_plane_intersection(&cyl, &obl) {
            CylPlaneCurve::Ellipse(e) => {
                let sampled = sample_cyl_plane(&CylPlaneCurve::Ellipse(e), 48).unwrap();
                for p in &sampled.points {
                    assert!(
                        cyl_surface_distance(&cyl, *p).abs() < 1e-7,
                        "ellipse on cyl"
                    );
                }
            }
            other => panic!("expected ellipse, got {other:?}"),
        }
    }

    /// Torus elbow ∩ crossing leg cylinder: a bite loop on the torus, all points
    /// on both surfaces (the perimeter U-turn case).
    #[test]
    fn torus_cyl_bite_loop() {
        let r = 0.275;
        // Torus elbow in the z=0 plane (axis Z), tube radius r, major 0.5.
        let torus = TorusSurf::new(Point3::ORIGIN, UnitVec3::Z, 0.5, r);
        // A straight leg along X just above the elbow, crossing it near θ=0.
        let cyl = CylSurf::new(Point3::new(0.0, 0.0, 0.35), UnitVec3::X, r);
        let loops = torus_cyl_intersection(&torus, &cyl, -0.8, 0.8, 96);
        assert!(!loops.is_empty(), "expected a bite loop");
        let lp = &loops[0];
        assert!(lp.points.len() >= 6);
        for p in &lp.points {
            // On the cylinder (marching-squares sampled — modest tolerance) ...
            assert!(
                cyl_surface_distance(&cyl, *p).abs() < 3e-2,
                "near cyl surface"
            );
            // ... and on the torus tube (distance from spine ≈ minor radius).
            let d = *p - torus.frame.origin;
            let in_plane = (d.x * d.x + d.y * d.y).sqrt();
            let spine_dist = ((in_plane - 0.5).powi(2) + d.z * d.z).sqrt();
            assert!((spine_dist - r).abs() < 3e-2, "on torus tube: {spine_dist}");
        }
    }

    /// Cap plane of an X-filament (normal X) cutting a Y-cylinder → ruling lines.
    #[test]
    fn cyl_plane_parallel_axis_lines() {
        let r = 0.275;
        let ycyl = cyl_y(0.0, r);
        let cap = Plane3::from_origin_normal(Point3::new(0.1, 0.0, 0.0), UnitVec3::X);
        match cyl_plane_intersection(&ycyl, &cap) {
            CylPlaneCurve::Lines(ls) => assert_eq!(ls.len(), 2),
            other => panic!("expected 2 lines, got {other:?}"),
        }
    }
}
