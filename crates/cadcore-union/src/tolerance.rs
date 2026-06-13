//! The tolerance ladder.
//!
//! A Boolean engine lives and dies by an explicit, ordered tolerance model.
//! Every stage states which rung it works at; mixing rungs implicitly is how
//! the legacy implementation accumulated µm-level drift that SolidWorks
//! punished by dropping faces.

/// The tolerance model (all lengths in **mm**).
///
/// Two distinct axes:
///
/// * the **accuracy ladder** `refine ≤ model ≤ knit` — how exact geometry is
///   at each stage, ending at the contract with the receiving CAD;
/// * **search radii** (`weld`, `sliver`) — how far the engine *looks* when
///   matching joints or detecting degeneracies.  These are legitimately much
///   larger than `knit`: a joint found within `weld` is subsequently refined
///   to `refine` accuracy, so a wide search radius does not loosen output.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tolerances {
    /// Numerical convergence target for iterative refinement (cyclic
    /// projection, triple points).  Geometry is considered *exact* below
    /// this.  Default `1e-12`.
    pub refine: f64,
    /// Modelling tolerance: two points closer than this are the same point;
    /// used for vertex identity.  Default `1e-9`.
    pub model: f64,
    /// Receiving-CAD knitting tolerance — the contract with the outside
    /// world.  Every emitted edge must lie within this distance of every
    /// face it bounds.  SolidWorks knits at ≈1 µm → default `1e-3`.
    pub knit: f64,
    /// Joint search radius: trim-curve endpoints within this distance are
    /// candidates for the same joint (then refined onto it exactly).  Gaps
    /// larger than this are *errors*, not weld candidates.  Default `2e-2`
    /// (matches the legacy chain tolerance).
    pub weld: f64,
    /// Sliver threshold: faces/regions thinner than this are degenerate and
    /// must be collapsed rather than emitted (near-tangent crossings produce
    /// them; SpaceClaim Combine chokes on webs ≤ 0.15 mm).  Default `5e-2`.
    pub sliver: f64,
}

impl Default for Tolerances {
    fn default() -> Self {
        Self {
            refine: 1e-12,
            model: 1e-9,
            weld: 2e-2,
            knit: 1e-3,
            sliver: 5e-2,
        }
    }
}

impl Tolerances {
    /// Validate the model.  Call once at pipeline start.
    pub fn check(&self) -> Result<(), String> {
        let ladder = [
            ("refine", self.refine),
            ("model", self.model),
            ("knit", self.knit),
        ];
        for w in ladder.windows(2) {
            if w[0].1 > w[1].1 {
                return Err(format!(
                    "accuracy ladder violated: {} ({:.1e}) > {} ({:.1e})",
                    w[0].0, w[0].1, w[1].0, w[1].1
                ));
            }
        }
        if self.weld < self.model {
            return Err(format!(
                "weld search radius ({:.1e}) below model tolerance ({:.1e})",
                self.weld, self.model
            ));
        }
        if self.sliver <= 0.0 {
            return Err("sliver threshold must be positive".into());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_ladder_is_ordered() {
        Tolerances::default().check().unwrap();
    }

    #[test]
    fn inverted_ladder_rejected() {
        let t = Tolerances {
            model: 1.0,
            ..Default::default()
        };
        assert!(t.check().is_err()); // model > knit AND weld < model
    }
}
