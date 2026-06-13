//! Generic parametric surfaces and the **swept-tube** surface — the foundation
//! of the generalized Boolean union (see `docs/swept-tube-ssi-pipeline.md`).
//!
//! A filament is a circular profile of radius `r` swept along its centerline
//! `C(s)` using a **rotation-minimizing frame** (RMF / parallel transport — NOT
//! Frenet, which is unstable at low curvature / inflections):
//!
//! ```text
//! S(s, θ) = C(s) + r·(cos θ · N(s) + sin θ · B(s))
//! ```
//!
//! Straight runs are a single line segment (≡ a cylinder); U-turns are arc
//! segments (≡ a torus patch) — but BOTH are the same `SweptTubeSurface`, so the
//! generic surface∩surface intersection treats them uniformly.

use cadcore_math::{Point3, UnitVec3, Vec3};

const TAU: f64 = std::f64::consts::PI * 2.0;

/// First partial derivatives of a surface at a parameter point.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceDerivatives {
    /// ∂S/∂u
    pub du: Vec3,
    /// ∂S/∂v
    pub dv: Vec3,
}

/// Diagnostic result for projecting a world point to a swept tube.
#[derive(Clone, Copy, Debug)]
pub struct SweptProjection {
    /// Arc-length parameter of the closest centerline point.
    pub u: f64,
    /// Cross-section angle around the tube.
    pub v: f64,
    /// Closest centerline point.
    pub center: Point3,
    /// Distance from `p` to the closest centerline point.
    pub centerline_distance: f64,
    /// Absolute distance from `p` to the tube surface.
    pub surface_distance: f64,
    /// True when the closest centerline point is on a segment/domain boundary.
    pub clamped: bool,
    /// Projection convergence flag. Analytic line/arc projection is exact.
    pub converged: bool,
    /// Number of refinement iterations used after the initial guess.
    pub iterations: usize,
}

/// A surface evaluable in a `(u, v)` parameter domain.
pub trait ParamSurface {
    /// Surface point at `(u, v)`.
    fn point(&self, u: f64, v: f64) -> Point3;
    /// First partials at `(u, v)`.
    fn derivatives(&self, u: f64, v: f64) -> SurfaceDerivatives;
    /// Unit surface normal (default: `du × dv` normalized).
    fn normal(&self, u: f64, v: f64) -> Vec3 {
        let d = self.derivatives(u, v);
        d.du.cross(d.dv).normalize()
    }
    /// `u` parameter range.
    fn u_domain(&self) -> (f64, f64);
    /// `v` parameter range (default: full circle).
    fn v_domain(&self) -> (f64, f64) {
        (0.0, TAU)
    }
    /// Closest `(u, v)` to a world point (best-effort).
    fn project_point(&self, p: Point3) -> Option<(f64, f64)>;
}

// ── Centerline ─────────────────────────────────────────────────────────────────

/// One analytic segment of a filament centerline, parameterised by **arc
/// length** local to the segment (`0..length`).
#[derive(Clone, Copy, Debug)]
pub enum CenterlineSeg {
    /// Straight run.
    Line {
        /// Segment start.
        p0: Point3,
        /// Unit direction.
        dir: UnitVec3,
        /// Length.
        length: f64,
    },
    /// Circular arc (the U-turn / rounded corner).
    Arc {
        /// Arc centre.
        center: Point3,
        /// Arc plane normal (axis).
        axis: UnitVec3,
        /// 0-angle reference direction in the plane.
        xref: UnitVec3,
        /// Arc radius.
        radius: f64,
        /// Start angle (rad).
        ang0: f64,
        /// End angle (rad); arc runs `ang0 → ang1`.
        ang1: f64,
    },
}

impl CenterlineSeg {
    /// Arc length of the segment.
    pub fn length(&self) -> f64 {
        match *self {
            CenterlineSeg::Line { length, .. } => length,
            CenterlineSeg::Arc {
                radius, ang0, ang1, ..
            } => radius * (ang1 - ang0).abs(),
        }
    }

    /// Point at local arc length `s ∈ [0, length]`.
    pub fn point(&self, s: f64) -> Point3 {
        match *self {
            CenterlineSeg::Line { p0, dir, .. } => p0 + dir * s,
            CenterlineSeg::Arc {
                center,
                axis,
                xref,
                radius,
                ang0,
                ang1,
            } => {
                let yref = axis.as_vec().cross(xref.as_vec());
                let sgn = (ang1 - ang0).signum();
                let a = ang0 + sgn * (s / radius);
                let (sn, cs) = a.sin_cos();
                center + xref.as_vec() * (cs * radius) + yref * (sn * radius)
            }
        }
    }

    /// Unit tangent at local arc length `s`.
    pub fn tangent(&self, s: f64) -> UnitVec3 {
        match *self {
            CenterlineSeg::Line { dir, .. } => dir,
            CenterlineSeg::Arc {
                axis,
                xref,
                radius,
                ang0,
                ang1,
                ..
            } => {
                let yref = axis.as_vec().cross(xref.as_vec());
                let sgn = (ang1 - ang0).signum();
                let a = ang0 + sgn * (s / radius);
                let (sn, cs) = a.sin_cos();
                // d/da of (cos a xref + sin a yref) = (-sin a xref + cos a yref)
                let t = (xref.as_vec() * (-sn) + yref * cs) * sgn;
                UnitVec3::try_from_vec(t).unwrap_or(xref)
            }
        }
    }

    fn project_local(&self, p: Point3) -> (f64, bool) {
        match *self {
            CenterlineSeg::Line { p0, dir, length } => {
                let raw = dir.dot_vec(p - p0);
                let s = raw.clamp(0.0, length);
                (s, raw < 0.0 || raw > length)
            }
            CenterlineSeg::Arc {
                center,
                axis,
                xref,
                radius,
                ang0,
                ang1,
            } => {
                let yref = axis.as_vec().cross(xref.as_vec());
                let q = p - center;
                let raw_ang = yref.dot(q).atan2(xref.dot_vec(q));
                let sweep = (ang1 - ang0).abs();
                if sweep <= 1e-14 || radius <= 0.0 {
                    return (0.0, true);
                }

                let sgn = (ang1 - ang0).signum();
                let progress = if sgn >= 0.0 {
                    (raw_ang - ang0).rem_euclid(TAU)
                } else {
                    (ang0 - raw_ang).rem_euclid(TAU)
                };

                if progress <= sweep {
                    (progress * radius, false)
                } else {
                    let len = self.length();
                    let d0 = (p - self.point(0.0)).length_sq();
                    let d1 = (p - self.point(len)).length_sq();
                    if d0 <= d1 {
                        (0.0, true)
                    } else {
                        (len, true)
                    }
                }
            }
        }
    }
}

// ── RMF station ─────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct Station {
    s: f64,      // global arc length
    c: Point3,   // centerline point
    t: UnitVec3, // tangent
    n: Vec3,     // RMF normal (unit, ⟂ t)
}

/// Double-reflection step (Wang et al. 2008): transport reference vector `r0`
/// from `(x0, t0)` to `(x1, t1)`.  Robust, no per-step rotation matrix.
fn double_reflect(x0: Point3, t0: Vec3, r0: Vec3, x1: Point3, t1: Vec3) -> Vec3 {
    let v1 = x1 - x0;
    let c1 = v1.dot(v1);
    let (r_l, t_l) = if c1 > 1e-18 {
        let rl = r0 - v1 * (2.0 / c1 * v1.dot(r0));
        let tl = t0 - v1 * (2.0 / c1 * v1.dot(t0));
        (rl, tl)
    } else {
        (r0, t0)
    };
    let v2 = t1 - t_l;
    let c2 = v2.dot(v2);
    let r1 = if c2 > 1e-18 {
        r_l - v2 * (2.0 / c2 * v2.dot(r_l))
    } else {
        r_l
    };
    // Re-orthogonalise against t1 and normalise.
    let r1 = r1 - t1 * (t1.dot(r1));
    r1.normalize()
}

// ── SweptTubeSurface ─────────────────────────────────────────────────────────────

/// A circular tube of radius `radius` swept along `centerline` with a
/// rotation-minimizing frame.  `u` = arc length along the centerline,
/// `v` = θ around the cross-section.
#[derive(Clone, Debug)]
pub struct SweptTubeSurface {
    segs: Vec<CenterlineSeg>,
    cum: Vec<f64>, // cumulative arc length, len = segs.len()+1
    radius: f64,
    stations: Vec<Station>, // RMF stations, sorted by s
}

impl SweptTubeSurface {
    /// Build from centerline segments and a constant radius.  Stations for the
    /// RMF are sampled adaptively (denser on arcs).
    pub fn new(segs: Vec<CenterlineSeg>, radius: f64) -> Self {
        assert!(!segs.is_empty(), "empty centerline");
        let mut cum = Vec::with_capacity(segs.len() + 1);
        cum.push(0.0);
        for sg in &segs {
            cum.push(cum.last().unwrap() + sg.length());
        }

        // Sample stations: endpoints of each segment + interior on arcs.
        let mut raw: Vec<(f64, Point3, UnitVec3)> = Vec::new();
        for (k, sg) in segs.iter().enumerate() {
            let base = cum[k];
            let len = sg.length();
            let n = match sg {
                CenterlineSeg::Line { .. } => 1,
                CenterlineSeg::Arc { .. } => 24.max((len / (radius.max(1e-6) * 0.25)) as usize),
            };
            for i in 0..=n {
                let sl = len * (i as f64 / n as f64);
                // Avoid duplicating the shared joint (except the very first).
                if i == 0 && k > 0 {
                    continue;
                }
                raw.push((base + sl, sg.point(sl), sg.tangent(sl)));
            }
        }

        // Seed normal ⟂ first tangent.
        let t0 = raw[0].2;
        let seed = {
            let a = if t0.as_vec().x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            let n = a - t0.as_vec() * t0.dot_vec(a);
            n.normalize()
        };
        let mut stations = Vec::with_capacity(raw.len());
        stations.push(Station {
            s: raw[0].0,
            c: raw[0].1,
            t: t0,
            n: seed,
        });
        for i in 1..raw.len() {
            let prev = stations[i - 1];
            let (s, c, t) = raw[i];
            let n = double_reflect(prev.c, prev.t.as_vec(), prev.n, c, t.as_vec());
            stations.push(Station { s, c, t, n });
        }

        Self {
            segs,
            cum,
            radius,
            stations,
        }
    }

    /// Total centerline length (== `u` upper bound).
    pub fn length(&self) -> f64 {
        *self.cum.last().unwrap()
    }

    /// Tube radius.
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// Locate the centerline segment containing global arc length `u`.
    fn locate(&self, u: f64) -> (usize, f64) {
        let u = u.clamp(0.0, self.length());
        // last segment whose start <= u
        let mut k = 0;
        while k + 1 < self.segs.len() && self.cum[k + 1] <= u {
            k += 1;
        }
        (k, u - self.cum[k])
    }

    /// Centerline point + tangent at global arc length `u`.
    pub fn centerline(&self, u: f64) -> (Point3, UnitVec3) {
        let (k, sl) = self.locate(u);
        (self.segs[k].point(sl), self.segs[k].tangent(sl))
    }

    /// RMF `(T, N, B)` at global arc length `u`.
    pub fn frame(&self, u: f64) -> (UnitVec3, Vec3, Vec3) {
        let (c, t) = self.centerline(u);
        // Transport the nearest station's normal to exact u.
        let j = self.nearest_station(u);
        let st = self.stations[j];
        let n = double_reflect(st.c, st.t.as_vec(), st.n, c, t.as_vec());
        let b = t.as_vec().cross(n);
        (t, n, b)
    }

    /// Closest centerline parameter and point for `p`.
    pub fn closest_centerline_point(&self, p: Point3) -> Option<(f64, Point3, bool)> {
        let mut best: Option<(f64, Point3, bool, f64)> = None;
        for (k, sg) in self.segs.iter().enumerate() {
            let (local, clamped) = sg.project_local(p);
            let u = self.cum[k] + local;
            let c = sg.point(local);
            let d2 = (p - c).length_sq();
            if best.map_or(true, |(_, _, _, best_d2)| d2 < best_d2) {
                best = Some((u, c, clamped, d2));
            }
        }
        best.map(|(u, c, clamped, _)| (u, c, clamped))
    }

    /// Project `p` to the tube and return diagnostics used by SSI refinement.
    pub fn project_point_diagnostics(&self, p: Point3) -> Option<SweptProjection> {
        let (u, center, clamped) = self.closest_centerline_point(p)?;
        let (_t, n, b) = self.frame(u);
        let w = p - center;
        let v = w.dot(b).atan2(w.dot(n)).rem_euclid(TAU);
        let centerline_distance = w.length();
        Some(SweptProjection {
            u,
            v,
            center,
            centerline_distance,
            surface_distance: (centerline_distance - self.radius).abs(),
            clamped,
            converged: true,
            iterations: 0,
        })
    }

    /// Signed distance to the swept tube solid (negative inside, positive outside).
    pub fn signed_distance(&self, p: Point3) -> Option<f64> {
        self.project_point_diagnostics(p)
            .map(|proj| proj.centerline_distance - self.radius)
    }

    fn nearest_station(&self, u: f64) -> usize {
        // binary search on s
        let mut lo = 0usize;
        let mut hi = self.stations.len() - 1;
        while lo < hi {
            let mid = (lo + hi) / 2;
            if self.stations[mid].s < u {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        // pick closer of lo, lo-1
        if lo > 0 && (u - self.stations[lo - 1].s).abs() < (self.stations[lo].s - u).abs() {
            lo - 1
        } else {
            lo
        }
    }
}

impl ParamSurface for SweptTubeSurface {
    fn point(&self, u: f64, v: f64) -> Point3 {
        let (_t, n, b) = self.frame(u);
        let (c, _) = self.centerline(u);
        let (sn, cs) = v.sin_cos();
        c + n * (cs * self.radius) + b * (sn * self.radius)
    }

    fn normal(&self, u: f64, v: f64) -> Vec3 {
        // The tube surface normal is exactly the radial direction from the
        // centerline to the surface point (outward).  This is unambiguous and
        // avoids the `du × dv` sign depending on the parametrisation handedness.
        let (_t, n, b) = self.frame(u);
        let (sn, cs) = v.sin_cos();
        (n * cs + b * sn).normalize()
    }

    fn derivatives(&self, u: f64, v: f64) -> SurfaceDerivatives {
        // dv analytic; du by central finite difference (frame variation included).
        let (_t, n, b) = self.frame(u);
        let (sn, cs) = v.sin_cos();
        let dv = (n * (-sn) + b * cs) * self.radius;
        let h = (self.length() * 1e-4).max(1e-6);
        let up = (u + h).min(self.length());
        let dn = (u - h).max(0.0);
        let du = (self.point(up, v) - self.point(dn, v)) / (up - dn);
        SurfaceDerivatives { du, dv }
    }

    fn u_domain(&self) -> (f64, f64) {
        (0.0, self.length())
    }

    fn project_point(&self, p: Point3) -> Option<(f64, f64)> {
        self.project_point_diagnostics(p)
            .map(|proj| (proj.u, proj.v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line_seg(p0: Point3, p1: Point3) -> CenterlineSeg {
        let d = p1 - p0;
        let length = d.length();
        CenterlineSeg::Line {
            p0,
            dir: UnitVec3::try_from_vec(d).unwrap(),
            length,
        }
    }

    /// Straight swept tube == cylinder: every surface point is exactly `r` from
    /// the axis, and the normal points radially outward.
    #[test]
    fn straight_tube_is_a_cylinder() {
        let r = 0.275;
        let tube = SweptTubeSurface::new(
            vec![line_seg(
                Point3::new(-2.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
            )],
            r,
        );
        for iu in 0..=10 {
            for iv in 0..12 {
                let u = tube.length() * iu as f64 / 10.0;
                let v = TAU * iv as f64 / 12.0;
                let p = tube.point(u, v);
                // distance from X axis
                let radial = (p.y * p.y + p.z * p.z).sqrt();
                assert!((radial - r).abs() < 1e-9, "radial={radial}");
                // normal radially outward (⟂ axis, parallel to (0,y,z))
                let nrm = tube.normal(u, v);
                assert!(nrm.x.abs() < 1e-6, "normal not ⟂ axis: {nrm:?}");
                let radial_dir = Vec3::new(0.0, p.y, p.z).normalize();
                assert!(nrm.dot(radial_dir) > 0.999, "normal not outward");
            }
        }
    }

    /// Arc swept tube ≡ torus patch: every surface point is exactly `r` from the
    /// centerline arc; the frame stays orthonormal (no twist).
    #[test]
    fn arc_tube_points_on_tube() {
        let r = 0.275;
        let rad = 0.5;
        let arc = CenterlineSeg::Arc {
            center: Point3::ORIGIN,
            axis: UnitVec3::Z,
            xref: UnitVec3::X,
            radius: rad,
            ang0: 0.0,
            ang1: std::f64::consts::FRAC_PI_2,
        };
        let tube = SweptTubeSurface::new(vec![arc], r);
        for iu in 0..=16 {
            for iv in 0..16 {
                let u = tube.length() * iu as f64 / 16.0;
                let v = TAU * iv as f64 / 16.0;
                let p = tube.point(u, v);
                let (c, t) = tube.centerline(u);
                // distance to the centerline point == r (cross-section is ⟂ t)
                let w = p - c;
                assert!(
                    (w.length() - r).abs() < 1e-9,
                    "tube radius off: {}",
                    w.length()
                );
                // cross-section is perpendicular to the tangent
                assert!(t.dot_vec(w).abs() < 1e-9, "cross-section not ⟂ tangent");
                // frame orthonormal
                let (tt, nn, bb) = tube.frame(u);
                assert!((nn.length() - 1.0).abs() < 1e-9);
                assert!((bb.length() - 1.0).abs() < 1e-9);
                assert!(tt.dot_vec(nn).abs() < 1e-9);
                assert!(tt.dot_vec(bb).abs() < 1e-9);
                assert!(nn.dot(bb).abs() < 1e-9);
            }
        }
    }

    /// Derivatives are consistent: the normal is perpendicular to both partials.
    #[test]
    fn derivatives_consistent_with_normal() {
        let r = 0.3;
        let tube = SweptTubeSurface::new(
            vec![CenterlineSeg::Arc {
                center: Point3::ORIGIN,
                axis: UnitVec3::Z,
                xref: UnitVec3::X,
                radius: 0.6,
                ang0: 0.0,
                ang1: std::f64::consts::FRAC_PI_2,
            }],
            r,
        );
        for iu in 1..16 {
            for iv in 0..12 {
                let u = tube.length() * iu as f64 / 16.0;
                let v = TAU * iv as f64 / 12.0;
                let d = tube.derivatives(u, v);
                let nrm = d.du.cross(d.dv).normalize();
                assert!(nrm.dot(d.dv).abs() < 1e-6, "normal·dv={}", nrm.dot(d.dv));
                // du has a tangential component but normal must be ⟂ du too
                assert!(nrm.dot(d.du).abs() < 1e-4, "normal·du={}", nrm.dot(d.du));
            }
        }
    }

    /// `project_point` round-trips a surface point back to its `(u, v)`.
    #[test]
    fn project_roundtrip() {
        let r = 0.275;
        let tube = SweptTubeSurface::new(
            vec![line_seg(
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(4.0, 1.0, 0.0),
            )],
            r,
        );
        let (u0, v0) = (2.3, 1.1);
        let p = tube.point(u0, v0);
        let (u1, v1) = tube.project_point(p).unwrap();
        let p2 = tube.point(u1, v1);
        assert!(
            (p - p2).length() < 1e-2,
            "project roundtrip off: {}",
            (p - p2).length()
        );
    }

    #[test]
    fn project_arc_is_exact_not_grid_quantized() {
        let r = 0.275;
        let tube = SweptTubeSurface::new(
            vec![CenterlineSeg::Arc {
                center: Point3::new(3.0, 4.0, 0.0),
                axis: UnitVec3::Z,
                xref: UnitVec3::X,
                radius: 0.5,
                ang0: 2.5,
                ang1: 2.5 - std::f64::consts::FRAC_PI_2,
            }],
            r,
        );
        let u0 = tube.length() * 0.37;
        let v0 = 1.8;
        let p = tube.point(u0, v0);
        let proj = tube.project_point_diagnostics(p).unwrap();
        let p2 = tube.point(proj.u, proj.v);
        assert!(
            (p - p2).length() < 1e-9,
            "arc projection off: {}",
            (p - p2).length()
        );
        assert!(proj.surface_distance < 1e-9);
        assert!(!proj.clamped);
    }
}
