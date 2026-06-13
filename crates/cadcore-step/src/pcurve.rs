//! Parametric **(u,v) curves** ("pcurves") for STEP face boundaries.
//!
//! # Why
//!
//! A robust STEP `ADVANCED_FACE` gives, for every boundary edge, BOTH the 3D
//! edge curve AND its 2D representation in the face surface's parameter space
//! (the *pcurve*).  Importers that only get 3D curves must re-project onto the
//! surface, which is a heuristic that breaks on periodic surfaces (cylinder /
//! torus) at the seam — producing "Unorientable shape" / "Invalid imbrication
//! of wires" / non-manifold sewing in OpenCASCADE/ANSYS.
//!
//! This module provides the *math*: map a 3D point on a supported surface to
//! its `(u,v)` parameters, and build a continuous (seam-unwrapped) pcurve
//! polyline for an edge.  Emission of the STEP entities is in `entities.rs`.
//!
//! Parametrisations match the STEP surface definitions exactly:
//! * `PLANE`:               `u = (P−o)·x`, `v = (P−o)·y`.
//! * `CYLINDRICAL_SURFACE`: `u = θ = atan2((P−o)·y, (P−o)·x)`, `v = (P−o)·axis`.

use std::f64::consts::PI;

use cadcore_geom::{CylSurf, Plane3};
use cadcore_math::Point3;

/// `(u, v)` of a 3D point on a plane (`u` along `frame.x`, `v` along `frame.y`).
pub(crate) fn plane_uv(pl: &Plane3, p: Point3) -> (f64, f64) {
    let l = pl.frame.to_local_point(p);
    (l.x, l.y)
}

/// `(u, v) = (θ, axial)` of a 3D point on a cylinder.  `θ ∈ (−π, π]` (raw,
/// un-unwrapped — use [`unwrap_u`] along an edge to get a continuous pcurve).
pub(crate) fn cyl_uv(c: &CylSurf, p: Point3) -> (f64, f64) {
    let l = c.frame.to_local_point(p);
    (l.y.atan2(l.x), l.z)
}

/// Unwrap a sequence of angular `u` values so consecutive entries never jump by
/// more than π — yields a continuous pcurve across the `u = ±π` seam.
pub(crate) fn unwrap_u(us: &mut [f64]) {
    for i in 1..us.len() {
        while us[i] - us[i - 1] > PI {
            us[i] -= 2.0 * PI;
        }
        while us[i] - us[i - 1] < -PI {
            us[i] += 2.0 * PI;
        }
    }
}

/// Build a continuous 2D pcurve `(u,v)` polyline for an edge sampled as 3D
/// `pts` on a cylinder.  `u` is unwrapped for seam continuity; then the whole
/// curve is shifted by a multiple of 2π so it starts in `[start_u−π, start_u+π]`
/// of the supplied `anchor_u` (lets the caller align both ends of a seam, or
/// keep a hole near the face's chosen branch).
pub(crate) fn cyl_pcurve(c: &CylSurf, pts: &[Point3], anchor_u: Option<f64>) -> Vec<(f64, f64)> {
    let mut uv: Vec<(f64, f64)> = pts.iter().map(|&p| cyl_uv(c, p)).collect();
    let mut us: Vec<f64> = uv.iter().map(|p| p.0).collect();
    unwrap_u(&mut us);
    if let (Some(a), Some(&first)) = (anchor_u, us.first()) {
        let mut shift = 0.0;
        let mut d = first - a;
        while d > PI {
            d -= 2.0 * PI;
            shift -= 2.0 * PI;
        }
        while d < -PI {
            d += 2.0 * PI;
            shift += 2.0 * PI;
        }
        for u in &mut us {
            *u += shift;
        }
    }
    for (slot, u) in uv.iter_mut().zip(us) {
        slot.0 = u;
    }
    uv
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_geom::Plane3;
    use cadcore_math::{Point3, UnitVec3};

    /// Cylinder round-trip: param of a point, fed back to `point_at`, returns it.
    #[test]
    fn cyl_uv_roundtrip() {
        let c = CylSurf::new(Point3::new(1.0, 2.0, 3.0), UnitVec3::Z, 0.275);
        for &th in &[0.0_f64, 0.7, 2.5, -1.1, PI - 0.01] {
            for &z in &[0.0_f64, 1.3, -0.4] {
                let p = c.point_at(th, z);
                let (u, v) = cyl_uv(&c, p);
                let back = c.point_at(u, v);
                assert!(
                    (back - p).length() < 1e-9,
                    "th={th} z={z}: {back:?} vs {p:?}"
                );
            }
        }
    }

    /// Plane round-trip: `frame.origin + u·x + v·y` recovers `(u,v)`.
    #[test]
    fn plane_uv_roundtrip() {
        let pl = Plane3::from_origin_normal(Point3::new(0.0, 0.0, 1.0), UnitVec3::Z);
        let p = pl.frame.origin + pl.frame.x * 0.4 + pl.frame.y * (-0.9);
        let (u, v) = plane_uv(&pl, p);
        assert!(
            (u - 0.4).abs() < 1e-12 && (v + 0.9).abs() < 1e-12,
            "u={u} v={v}"
        );
    }

    /// A near-full-circle pcurve on a cylinder must be monotonic in `u` (no 2π
    /// jump at the ±π seam after unwrapping).
    #[test]
    fn cyl_pcurve_unwraps_seam() {
        let c = CylSurf::new(Point3::ORIGIN, UnitVec3::Z, 1.0);
        // sample θ from 2.0 rad up across the +π seam to -2.0-equivalent.
        let pts: Vec<Point3> = (0..=10)
            .map(|k| {
                let th = 2.0 + (k as f64) * 0.25; // 2.0 .. 4.5 rad, crosses π
                c.point_at(th, 0.0)
            })
            .collect();
        let uv = cyl_pcurve(&c, &pts, None);
        for w in uv.windows(2) {
            assert!(
                (w[1].0 - w[0].0).abs() < PI,
                "seam jump: {} -> {}",
                w[0].0,
                w[1].0
            );
        }
        // and each maps back to the surface
        for (&p, &(u, v)) in pts.iter().zip(&uv) {
            assert!((c.point_at(u, v) - p).length() < 1e-9);
        }
    }
}
