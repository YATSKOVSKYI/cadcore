//! Geometric services of the union engine.
//!
//! * [`classify`] — surface and surface-pair taxonomy with explicit tangency
//!   margins (the engine must *know* when a configuration is marginal, not
//!   discover it via downstream numerical misery).
//! * [`refine`] — the projection oracle: pulls every trim-curve point, joint
//!   and weld onto the analytic surfaces that the STEP writer actually emits.
//! * [`intersect`] — the intersection oracle: predictor–corrector tracer
//!   whose every output point lies on both surfaces to refinement accuracy
//!   (any axis angle, unequal radii, torus pairs; tangency reported, never
//!   guessed).

pub mod classify;
pub mod composite;
pub mod intersect;
pub mod refine;
