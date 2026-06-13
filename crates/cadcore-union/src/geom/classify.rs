//! Surface-pair taxonomy.
//!
//! The engine must know *upfront* whether a configuration is generic or
//! marginal.  Near-tangent crossings (layer height approaching the filament
//! diameter, HL → 2r) produce intersection curves that flutter, slivers below
//! the receiving CAD's capabilities, and projection sets that converge slowly.
//! Classification turns those from downstream numerical surprises into
//! explicit, testable branches.

use cadcore_geom::CylSurf;

/// Relationship between two cylinder axes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxisRelation {
    /// Axes parallel within the angular margin.
    Parallel,
    /// Axes perpendicular within the angular margin.
    Perpendicular,
    /// Anything in between.
    Oblique,
}

/// Depth class of a cylinder×cylinder crossing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossingDepth {
    /// Surfaces don't reach each other.
    Disjoint,
    /// Axis distance within `tangent_margin` of `r1 + r2`: the intersection
    /// degenerates to (near-)tangency.  Special handling required — generic
    /// trim machinery must refuse this.
    NearTangent,
    /// Generic transversal crossing.
    Transversal,
    /// Axis distance within `tangent_margin` of `|r1 − r2|` (one tube nearly
    /// inside the other) — internal tangency.
    NearInternalTangent,
    /// One tube fully inside the other.
    Contained,
}

/// Classified cylinder×cylinder pair.
#[derive(Debug, Clone, Copy)]
pub struct CylCylClass {
    /// Axis relation.
    pub axes: AxisRelation,
    /// Crossing depth.
    pub depth: CrossingDepth,
    /// Distance between the two axes (closest approach, mm).
    pub axis_distance: f64,
    /// Angle between axes in radians, folded to `[0, π/2]`.
    pub angle: f64,
}

/// Classify a cylinder×cylinder pair.
///
/// `angular_margin` is in radians (how far from 0 / π/2 still counts as
/// parallel / perpendicular); `tangent_margin` in mm.
pub fn classify_cyl_cyl(
    a: &CylSurf,
    b: &CylSurf,
    angular_margin: f64,
    tangent_margin: f64,
) -> CylCylClass {
    let da = a.axis().as_vec();
    let db = b.axis().as_vec();
    let cosang = da.dot(db).abs().clamp(0.0, 1.0);
    let angle = cosang.acos(); // [0, π/2]

    let axes = if angle < angular_margin {
        AxisRelation::Parallel
    } else if (std::f64::consts::FRAC_PI_2 - angle) < angular_margin {
        AxisRelation::Perpendicular
    } else {
        AxisRelation::Oblique
    };

    // Closest distance between the two axis lines.
    let w = b.frame.origin - a.frame.origin;
    let cross = da.cross(db);
    let axis_distance = if cross.length() < 1e-12 {
        // parallel lines
        (w - da * da.dot(w)).length()
    } else {
        (w.dot(cross) / cross.length()).abs()
    };

    let rsum = a.radius + b.radius;
    let rdiff = (a.radius - b.radius).abs();
    let depth = if axis_distance > rsum + tangent_margin {
        CrossingDepth::Disjoint
    } else if axis_distance > rsum - tangent_margin {
        CrossingDepth::NearTangent
    } else if axis_distance > rdiff + tangent_margin {
        CrossingDepth::Transversal
    } else if axis_distance > rdiff - tangent_margin {
        CrossingDepth::NearInternalTangent
    } else {
        CrossingDepth::Contained
    };

    CylCylClass {
        axes,
        depth,
        axis_distance,
        angle,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_math::{Point3, UnitVec3};

    const ANG: f64 = 1e-3;
    const TAN: f64 = 5e-2;

    fn cyl(o: Point3, axis: UnitVec3, r: f64) -> CylSurf {
        CylSurf::new(o, axis, r)
    }

    #[test]
    fn woodpile_perpendicular_transversal() {
        // the production case: HL = 0.35, D = 0.55
        let a = cyl(Point3::new(0.0, 0.0, 0.0), UnitVec3::X, 0.275);
        let b = cyl(Point3::new(0.0, 0.0, 0.35), UnitVec3::Y, 0.275);
        let c = classify_cyl_cyl(&a, &b, ANG, TAN);
        assert_eq!(c.axes, AxisRelation::Perpendicular);
        assert_eq!(c.depth, CrossingDepth::Transversal);
        assert!((c.axis_distance - 0.35).abs() < 1e-12);
    }

    #[test]
    fn hl_at_diameter_is_near_tangent() {
        // HL = 0.55 = D: layers exactly touch — the marginal case the engine
        // must FLAG, not stumble into.
        let a = cyl(Point3::new(0.0, 0.0, 0.0), UnitVec3::X, 0.275);
        let b = cyl(Point3::new(0.0, 0.0, 0.55), UnitVec3::Y, 0.275);
        let c = classify_cyl_cyl(&a, &b, ANG, TAN);
        assert_eq!(c.depth, CrossingDepth::NearTangent);
    }

    #[test]
    fn hl_050_still_transversal() {
        // HL = 0.50 (the user's deepest ladder rung): web = 0.05 — inside the
        // tangent margin?  With margin 0.05 it is right at the boundary; the
        // engine treats it as near-tangent and applies the careful path.
        let a = cyl(Point3::new(0.0, 0.0, 0.0), UnitVec3::X, 0.275);
        let b = cyl(Point3::new(0.0, 0.0, 0.50), UnitVec3::Y, 0.275);
        let c = classify_cyl_cyl(&a, &b, ANG, 2e-2);
        assert_eq!(c.depth, CrossingDepth::Transversal);
        let c2 = classify_cyl_cyl(&a, &b, ANG, 6e-2);
        assert_eq!(c2.depth, CrossingDepth::NearTangent);
    }

    #[test]
    fn parallel_stacked() {
        let a = cyl(Point3::new(0.0, 0.0, 0.0), UnitVec3::X, 0.275);
        let b = cyl(Point3::new(0.0, 0.4, 0.1), UnitVec3::X, 0.275);
        let c = classify_cyl_cyl(&a, &b, ANG, TAN);
        assert_eq!(c.axes, AxisRelation::Parallel);
        assert_eq!(c.depth, CrossingDepth::Transversal);
        assert!((c.axis_distance - (0.4f64.hypot(0.1))).abs() < 1e-12);
    }

    #[test]
    fn oblique_detected() {
        let d = (UnitVec3::X.as_vec() + UnitVec3::Y.as_vec()).normalize();
        let a = cyl(Point3::new(0.0, 0.0, 0.0), UnitVec3::X, 0.275);
        let b = cyl(
            Point3::new(0.0, 0.0, 0.3),
            cadcore_math::UnitVec3::try_from_vec(d).unwrap(),
            0.275,
        );
        let c = classify_cyl_cyl(&a, &b, ANG, TAN);
        assert_eq!(c.axes, AxisRelation::Oblique);
    }

    #[test]
    fn disjoint_layers() {
        let a = cyl(Point3::new(0.0, 0.0, 0.0), UnitVec3::X, 0.275);
        let b = cyl(Point3::new(0.0, 0.0, 1.2), UnitVec3::Y, 0.275);
        let c = classify_cyl_cyl(&a, &b, ANG, TAN);
        assert_eq!(c.depth, CrossingDepth::Disjoint);
    }
}
