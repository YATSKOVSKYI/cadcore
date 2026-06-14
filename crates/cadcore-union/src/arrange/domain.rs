//! Per-face uv parametrisations.
//!
//! The arrangement works in a face's 2-D parameter rectangle; lifting cell
//! interior points back to 3-D drives keep/drop classification.  Conventions
//! match the STEP writer (cylinder: u = angle about the axis from `frame.x`,
//! v = axial; torus: u = θ along the ring, v = φ around the minor circle).

use cadcore_geom::{CylSurf, Plane3, TorusSurf};
use cadcore_math::Point3;

/// A face's uv domain.
#[derive(Debug, Clone, Copy)]
pub enum FaceDomain {
    /// Cylinder band: u = angle (periodic, 2π), v = axial in `[0, length]`.
    CylinderBand {
        /// Carrier surface.
        surf: CylSurf,
        /// Axial extent.
        length: f64,
    },
    /// Torus patch: u = ring angle θ in `[theta_lo, theta_hi]`,
    /// v = minor-circle angle φ (periodic, 2π).
    TorusPatch {
        /// Carrier surface.
        surf: TorusSurf,
        /// Ring-arc start.
        theta_lo: f64,
        /// Ring-arc end (signed sweep from `theta_lo`).
        theta_hi: f64,
    },
    /// Full-donut torus band: u = ring angle θ (periodic, 2π), v = minor-circle
    /// angle φ in `[0, 2π]` (the "open" axis, bounded by the φ-seam — itself a
    /// shared circle, so the band closes).  The mirror of [`Self::TorusPatch`]:
    /// here θ wraps and φ is cut, which is what stacked donuts need.
    TorusBand {
        /// Carrier surface.
        surf: TorusSurf,
    },
    /// Flat plane: u/v are the in-plane Cartesian coordinates (no periodic
    /// axis).  Used for end-cap disks.
    Plane {
        /// Carrier plane.
        plane: Plane3,
        /// Half-extent of the working rectangle (mm); the cap disk and its
        /// trims live well inside.
        half_extent: f64,
    },
}

impl FaceDomain {
    /// Project a 3-D point near the surface into uv.
    pub fn uv(&self, p: Point3) -> (f64, f64) {
        match self {
            FaceDomain::CylinderBand { surf, .. } => {
                let w = p - surf.frame.origin;
                let u = surf.frame.y.dot_vec(w).atan2(surf.frame.x.dot_vec(w));
                (u, surf.axis().dot_vec(w))
            }
            FaceDomain::TorusPatch {
                surf,
                theta_lo,
                theta_hi,
            } => {
                let w = p - surf.frame.origin;
                let raw = surf.frame.y.dot_vec(w).atan2(surf.frame.x.dot_vec(w));
                // θ is periodic; the patch is a sub-arc that may cross the
                // atan2 branch cut (±π).  Unwrap θ onto the branch nearest the
                // arc's midpoint so the open axis is continuous over the patch.
                let mid = 0.5 * (theta_lo + theta_hi);
                let theta = raw
                    + std::f64::consts::TAU
                        * ((mid - raw) / std::f64::consts::TAU).round();
                let h = surf.frame.z.dot_vec(w);
                let q = w - surf.frame.z.as_vec() * h;
                let a = q.length();
                let phi = h.atan2(a - surf.major_radius);
                (theta, phi)
            }
            FaceDomain::TorusBand { surf } => {
                let w = p - surf.frame.origin;
                let theta = surf.frame.y.dot_vec(w).atan2(surf.frame.x.dot_vec(w));
                let h = surf.frame.z.dot_vec(w);
                let q = w - surf.frame.z.as_vec() * h;
                let a = q.length();
                let phi = h.atan2(a - surf.major_radius).rem_euclid(std::f64::consts::TAU);
                (theta, phi)
            }
            FaceDomain::Plane { plane, .. } => {
                let w = p - plane.frame.origin;
                (plane.frame.x.dot_vec(w), plane.frame.y.dot_vec(w))
            }
        }
    }

    /// Lift uv back onto the surface.
    pub fn lift(&self, u: f64, v: f64) -> Point3 {
        match self {
            FaceDomain::CylinderBand { surf, .. } => {
                let radial = surf.frame.x * u.cos() + surf.frame.y * u.sin();
                surf.frame.origin + radial * surf.radius + surf.axis().as_vec() * v
            }
            FaceDomain::TorusPatch { surf, .. } => {
                let ring_dir = surf.frame.x * u.cos() + surf.frame.y * u.sin();
                let ring = surf.frame.origin + ring_dir * surf.major_radius;
                ring + ring_dir * (surf.minor_radius * v.cos())
                    + surf.frame.z.as_vec() * (surf.minor_radius * v.sin())
            }
            FaceDomain::TorusBand { surf } => {
                // u = θ (periodic), v = φ (open)
                let ring_dir = surf.frame.x * u.cos() + surf.frame.y * u.sin();
                let ring = surf.frame.origin + ring_dir * surf.major_radius;
                ring + ring_dir * (surf.minor_radius * v.cos())
                    + surf.frame.z.as_vec() * (surf.minor_radius * v.sin())
            }
            FaceDomain::Plane { plane, .. } => {
                plane.frame.origin + plane.frame.x * u + plane.frame.y * v
            }
        }
    }

    /// Which uv axis is periodic, with its period: `(axis, period)` where
    /// axis 0 = u, 1 = v.
    pub fn periodic_axis(&self) -> (usize, f64) {
        match self {
            FaceDomain::CylinderBand { .. } => (0, std::f64::consts::TAU),
            FaceDomain::TorusPatch { .. } => (1, std::f64::consts::TAU),
            FaceDomain::TorusBand { .. } => (0, std::f64::consts::TAU), // u = θ wraps
            // a plane has no periodic axis; report a huge non-period so the
            // DCEL seam never triggers within the working rectangle.
            FaceDomain::Plane { half_extent, .. } => (0, 4.0 * half_extent),
        }
    }

    /// The non-periodic axis range `(lo, hi)`.
    pub fn open_range(&self) -> (f64, f64) {
        match self {
            FaceDomain::CylinderBand { length, .. } => (0.0, *length),
            FaceDomain::TorusPatch {
                theta_lo, theta_hi, ..
            } => {
                if theta_lo <= theta_hi {
                    (*theta_lo, *theta_hi)
                } else {
                    (*theta_hi, *theta_lo)
                }
            }
            FaceDomain::TorusBand { .. } => (0.0, std::f64::consts::TAU), // v = φ
            FaceDomain::Plane { half_extent, .. } => (-half_extent, *half_extent),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_math::UnitVec3;

    #[test]
    fn cylinder_uv_roundtrip() {
        let d = FaceDomain::CylinderBand {
            surf: CylSurf::new(Point3::new(1.0, 2.0, 3.0), UnitVec3::X, 0.275),
            length: 4.0,
        };
        for &(u, v) in &[(0.3, 0.5), (-2.7, 3.9), (3.0, 0.0)] {
            let p = d.lift(u, v);
            let (u2, v2) = d.uv(p);
            assert!((u - u2).abs() < 1e-12 && (v - v2).abs() < 1e-12);
        }
    }

    #[test]
    fn torus_uv_roundtrip() {
        let d = FaceDomain::TorusPatch {
            surf: TorusSurf::new(Point3::new(0.0, 0.5, 0.0), UnitVec3::Z, 0.5, 0.275),
            theta_lo: -1.5,
            theta_hi: 0.0,
        };
        for &(u, v) in &[(-0.7, 0.4), (-1.2, -2.0), (-0.1, 3.0)] {
            let p = d.lift(u, v);
            let (u2, v2) = d.uv(p);
            let dv = (v - v2 + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
                - std::f64::consts::PI;
            assert!((u - u2).abs() < 1e-12 && dv.abs() < 1e-12, "u={u} v={v}");
        }
    }
}
