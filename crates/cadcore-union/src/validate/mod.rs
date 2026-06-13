//! In-process validation gates.
//!
//! These are Rust ports of the external auditors that were tuned until they
//! exactly reproduced SolidWorks' accept/reject behaviour on real models.
//! They run on the in-memory B-Rep *before* any STEP is written: a gate
//! failure means the output must not ship.
//!
//! | Gate | Module | External ancestor |
//! |------|--------|-------------------|
//! | 1. manifold edge sharing | [`manifold`] | `check_ap214_manifold` / python 2-use check |
//! | 2. uv wire sanity        | [`wires`] | `out/strict_audit.py` |
//! | 3. plane-cap wires       | [`wires`] | `out/audit_planes.py` |
//! | 4. edge-to-surface distance | [`distance`] | `out/dist_to_spline.py` |

pub mod distance;
pub mod manifold;
pub mod step_text;
pub mod wires;

use crate::config::UnionConfig;
use crate::diag::{Diagnostics, Severity};
use cadcore_topo::{BRep, ShellId};

/// Aggregated result of running all live gates on one shell.
#[derive(Debug)]
pub struct ValidationReport {
    /// Gate 1 violations (edge id debug strings + use descriptions).
    pub manifold_violations: Vec<String>,
    /// Gates 2+3 violations (uv wire closure / winding / intersections).
    pub wire_violations: Vec<String>,
    /// Gate 4 violations: `(edge, face, distance_mm)` debug strings.
    pub distance_violations: Vec<String>,
    /// Worst edge-to-surface deviation observed anywhere (mm).
    pub max_edge_deviation: f64,
    /// Whether wire violations count as gate failures (see
    /// `UnionConfig::wires_strict`).
    pub wires_strict: bool,
}

impl ValidationReport {
    /// True when every gate passed.
    pub fn passed(&self) -> bool {
        self.manifold_violations.is_empty()
            && (!self.wires_strict || self.wire_violations.is_empty())
            && self.distance_violations.is_empty()
    }
}

/// Run all live gates on one shell.
pub fn validate_shell(
    brep: &BRep,
    shell: ShellId,
    cfg: &UnionConfig,
    diag: &mut Diagnostics,
) -> ValidationReport {
    let tol = &cfg.tol;
    let mreport = manifold::check(brep, shell);
    let manifold_violations = mreport.violations;
    for w in &mreport.sense_warnings {
        diag.event(Severity::Warning, "gate.manifold.sense", w.clone());
    }
    let wire_violations = wires::check(brep, shell);
    let (distance_violations, max_edge_deviation) = distance::check(brep, shell, tol.knit);

    for v in &manifold_violations {
        diag.event(Severity::Error, "gate.manifold", v.clone());
    }
    let wire_sev = if cfg.wires_strict {
        Severity::Error
    } else {
        Severity::Warning
    };
    for v in &wire_violations {
        diag.event(wire_sev, "gate.wires", v.clone());
    }
    for v in &distance_violations {
        diag.event(Severity::Error, "gate.distance", v.clone());
    }
    diag.count_n("gate.manifold.violations", manifold_violations.len() as u64);
    diag.count_n("gate.wires.violations", wire_violations.len() as u64);
    diag.count_n("gate.distance.violations", distance_violations.len() as u64);

    ValidationReport {
        manifold_violations,
        wire_violations,
        distance_violations,
        max_edge_deviation,
        wires_strict: cfg.wires_strict,
    }
}
