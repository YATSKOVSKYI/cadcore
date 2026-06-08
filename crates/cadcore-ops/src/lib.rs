//! # cadcore-ops (local fork)
//!
//! Pure-Rust analytic operations for sweep construction and half-space trimming.

pub mod boolean;
pub mod sweep;
pub mod union;

pub use sweep::{
    analytic_path_from_polyline_samples, build_solid_box,
    build_solid_rounded_box, build_solid_rounded_box_full, build_solid_rounded_box_xz,
    clip_polyline, clip_polyline_with_radius,
    rounded_path_from_polyline, sharp_path_from_polyline_samples, solid_corner_centerline_radius,
    sweep_circle_along_path, sweep_circle_along_path_with_caps, sweep_circle_along_polyline,
    sweep_circle_along_polyline_with_caps, sweep_circle_along_rounded_polyline, ClipPlane,
    PathApproxOptions, SweepOptions, SweepPathSegment,
};

pub use boolean::half_space_cut_brep;
pub use union::union_solids;
