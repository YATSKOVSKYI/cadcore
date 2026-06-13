//! The staged union pipeline.
//!
//! ```text
//! Collect ─ gather faces/surfaces per solid, build spatial index
//! Classify ─ taxonomy of every candidate pair (geom::classify); marginal
//!            configurations flagged BEFORE any numerics run
//! Intersect ─ analytic oracle first, marching SSI fallback        (P2)
//! Refine ─ every curve/joint/cut onto emitted geometry (geom::refine)
//! Arrange ─ per-face uv arrangement, cell classification          (P3)
//! Stitch ─ shared-edge registry, loop assembly, winding rules     (P3)
//! Validate ─ gates 1–4 in-process (validate::*)
//! ```
//!
//! P0 ships the frame plus the `Validate` stage: callers run the legacy union
//! (`cadcore-ops`) and pass its shell through [`UnionPipeline::validate`].
//! Later phases move the construction stages in, one at a time, each behind
//! the same configuration and diagnostics.

use crate::config::UnionConfig;
use crate::diag::{Diagnostics, Severity};
use crate::validate::{validate_shell, ValidationReport};
use cadcore_topo::{BRep, ShellId};

/// Pipeline stages, in execution order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Stage {
    /// Gather faces and build the candidate-pair index.
    Collect,
    /// Pair taxonomy and marginal-case flagging.
    Classify,
    /// Intersection curves (analytic oracle / SSI fallback).
    Intersect,
    /// Refinement onto emitted geometry.
    Refine,
    /// Per-face uv arrangements and cell classification.
    Arrange,
    /// Topology assembly: shared edges, loops, winding.
    Stitch,
    /// In-process gates.
    Validate,
}

/// Result of a pipeline run.
#[derive(Debug)]
pub struct UnionOutcome {
    /// Per-gate findings of the final validation.
    pub report: ValidationReport,
    /// Collected diagnostics.
    pub diag: Diagnostics,
}

impl UnionOutcome {
    /// Output may be shipped (all gates green, no Error events).
    pub fn shippable(&self) -> bool {
        self.report.passed() && self.diag.clean_at(Severity::Error)
    }
}

/// The engine façade.
#[derive(Debug, Default)]
pub struct UnionPipeline {
    /// Engine configuration.
    pub config: UnionConfig,
}

impl UnionPipeline {
    /// New pipeline with the given configuration.
    pub fn new(config: UnionConfig) -> Self {
        Self { config }
    }

    /// P0 entry point: validate an already-built shell (e.g. the output of
    /// the legacy union) against all live gates.
    pub fn validate(&self, brep: &BRep, shell: ShellId) -> UnionOutcome {
        let mut diag = Diagnostics::new();
        if let Err(e) = self.config.tol.check() {
            diag.event(Severity::Error, "config.tolerances", e);
        }
        let report = validate_shell(brep, shell, &self.config, &mut diag);
        diag.count_n(
            "validate.max_edge_dev_nm",
            (report.max_edge_deviation * 1e6) as u64,
        );
        UnionOutcome { report, diag }
    }
}
