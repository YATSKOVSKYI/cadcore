//! Engine configuration.

use crate::tolerance::Tolerances;

/// Which trim families the pipeline is allowed to perform.
///
/// Mirrors the bisection axes that proved invaluable during debugging: any
/// family can be switched off to produce an A/B artefact for an external CAD
/// check without rebuilding the engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PassSet {
    /// Leg×leg crossing windows (transversal woodpile crossings).
    pub crossings: bool,
    /// Elbow-involved trims spanning junctions (mid-U merges).
    pub junction_spans: bool,
    /// Elbow×elbow corner bites (U-on-U).
    pub corners: bool,
    /// Capped butt-end mutual bites.
    pub butt_ends: bool,
}

impl Default for PassSet {
    fn default() -> Self {
        Self {
            crossings: true,
            junction_spans: true,
            corners: true,
            butt_ends: true,
        }
    }
}

/// How much work the engine spends double-checking itself.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CheckLevel {
    /// Gates run only at the end, on the final shell.
    Final,
    /// Gates additionally run after every pipeline stage (slower; pinpoints
    /// the stage that introduced a defect).
    PerStage,
    /// `PerStage` plus per-trim micro-checks (each committed trim immediately
    /// verified for edge sharing and surface distance).  Debugging tier.
    Paranoid,
}

/// Full engine configuration.
#[derive(Debug, Clone)]
pub struct UnionConfig {
    /// The tolerance ladder.
    pub tol: Tolerances,
    /// Enabled trim families.
    pub passes: PassSet,
    /// Self-checking intensity.
    pub checks: CheckLevel,
    /// Capture verbose per-pair diagnostic events (memory cost).
    pub trace_pairs: bool,
    /// Treat uv-wire issues (gates 2+3) as hard failures.
    ///
    /// OFF for legacy-union output: its in-memory loops are placeholders —
    /// the STEP writer derives real wires from `FaceExtent` templates, so
    /// BRep-level wire checks are advisory there.  The new pipeline (P3+)
    /// builds authoritative loops and turns this on.
    pub wires_strict: bool,
}

impl Default for UnionConfig {
    fn default() -> Self {
        Self {
            tol: Tolerances::default(),
            passes: PassSet::default(),
            checks: CheckLevel::Final,
            trace_pairs: false,
            wires_strict: false,
        }
    }
}
