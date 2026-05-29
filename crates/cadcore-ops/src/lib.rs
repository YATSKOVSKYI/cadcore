//! # cadcore-ops (local fork)
//!
//! Extended with `boolean::half_space_cut_brep` for solid half-space trimming.

pub mod sweep;
pub mod boolean;

pub use sweep::{
    analytic_path_from_polyline_samples, rounded_path_from_polyline,
    sharp_path_from_polyline_samples, solid_corner_centerline_radius,
    sweep_circle_along_path, sweep_circle_along_polyline,
    sweep_circle_along_rounded_polyline, PathApproxOptions,
    SweepOptions, SweepPathSegment, ClipPlane, clip_polyline,
    clip_polyline_with_radius, sweep_circle_along_path_with_caps,
    sweep_circle_along_polyline_with_caps, build_solid_box,
};

pub use boolean::half_space_cut_brep;
