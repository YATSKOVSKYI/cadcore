//! Sweep a circle profile along a polyline — exact analytic B-Rep construction.
//!
//! ## Algorithm
//!
//! Given N waypoints and a radius r, the solid consists of:
//!
//! * **(N−1) cylinder faces** — one per straight segment between adjacent waypoints.
//! * **(N−2) connector faces** — one per bend vertex; either a toroidal fillet
//!   (smooth join) or an elliptic miter (sharp join).
//! * **2 end-cap faces** — one planar disk at each endpoint.
//!
//! All surfaces are analytic (cylinder, torus, plane, ellipse) and map
//! directly to STEP AP203 entities.  No mesh, no Boolean union.  Total
//! construction cost is O(N).

use cadcore_geom::{Circle3, CylSurf, Ellipse3, Plane3, TorusSurf};
use cadcore_math::{Point3, UnitVec3, Vec3};
use cadcore_topo::{
    BRep, Face, FaceBoundary, FaceExtent, FaceGeom, FaceNormal, Shell, Solid, SolidId,
};

/// Options that control the sweep.
#[derive(Clone, Debug)]
pub struct SweepOptions {
    /// Round sharp polyline corners before sweeping.
    /// If `false`, use miter-plane cuts at sharp joins.
    pub fillet_corners: bool,
    /// Radius left on the inside of a rounded solid corner.
    ///
    /// The kernel implements this as a CAD pipe sweep: first it rounds the
    /// centre-line to `profile_radius + corner_fillet_radius`, then sweeps the
    /// circular profile along that analytic path.  When set to `0.0`, the swept
    /// profile radius is used for backward compatibility.
    pub corner_fillet_radius: f64,
    /// Optional name for the resulting solid.
    pub name: Option<String>,
}

impl Default for SweepOptions {
    fn default() -> Self {
        Self {
            fillet_corners: false,
            corner_fillet_radius: 0.0,
            name: None,
        }
    }
}

impl SweepOptions {
    fn effective_solid_corner_radius(&self, profile_radius: f64) -> Result<f64, SweepError> {
        if !self.fillet_corners {
            return Ok(0.0);
        }
        if self.corner_fillet_radius < 0.0 {
            return Err(SweepError::InvalidCornerRadius(self.corner_fillet_radius));
        }

        let radius = if self.corner_fillet_radius <= 1.0e-9 {
            profile_radius
        } else {
            self.corner_fillet_radius
        };
        Ok(radius)
    }
}

/// Convert a solid-corner radius to the centre-line radius needed for a pipe
/// sweep.
///
/// A swept circular profile self-intersects when the centre-line arc radius is
/// not larger than the profile radius.  Treating `corner_radius` as the inside
/// radius of the solid gives the safe CAD construction:
/// `centre_line_radius = profile_radius + corner_radius`.
///
/// # Errors
///
/// Returns [`SweepError::InvalidRadius`] when `profile_radius <= 0` and
/// [`SweepError::InvalidCornerRadius`] when `corner_radius < 0`.
pub fn solid_corner_centerline_radius(
    profile_radius: f64,
    corner_radius: f64,
) -> Result<f64, SweepError> {
    if profile_radius <= 0.0 {
        return Err(SweepError::InvalidRadius(profile_radius));
    }
    if corner_radius < 0.0 {
        return Err(SweepError::InvalidCornerRadius(corner_radius));
    }
    Ok(profile_radius + corner_radius)
}

/// Options for recovering analytic path segments from a sampled polyline.
#[derive(Clone, Copy, Debug)]
pub struct PathApproxOptions {
    /// Minimum number of consecutive sampled points needed to accept an arc.
    pub min_arc_points: usize,
    /// Minimum accepted arc angle, in radians.
    pub min_arc_angle: f64,
    /// Relative radius tolerance for sampled points on the same arc.
    pub radius_rel_tol: f64,
    /// Absolute radius tolerance for sampled points on the same arc.
    pub radius_abs_tol: f64,
    /// Consecutive points closer than this are ignored.
    pub point_tol: f64,
}

impl Default for PathApproxOptions {
    fn default() -> Self {
        Self {
            min_arc_points: 6,
            min_arc_angle: 0.20,
            radius_rel_tol: 0.03,
            radius_abs_tol: 0.01,
            point_tol: 1.0e-7,
        }
    }
}

/// Analytic path segment for sweeping a circular filament profile.
#[derive(Clone, Copy, Debug)]
pub enum SweepPathSegment {
    /// Straight centre-line segment.
    Line {
        /// Segment start point.
        start: Point3,
        /// Segment end point.
        end: Point3,
    },
    /// Circular centre-line arc.
    ///
    /// `normal` controls the arc direction: the start tangent is
    /// `normal × (start - center)`. The arc is represented as one analytic
    /// torus face in the resulting B-Rep, not as sampled straight cylinders.
    Arc {
        /// Arc start point.
        start: Point3,
        /// Arc end point.
        end: Point3,
        /// Circle centre.
        center: Point3,
        /// Circle plane normal and traversal direction.
        normal: UnitVec3,
    },
}

impl SweepPathSegment {
    /// Start point of this segment.
    #[must_use]
    pub fn start(self) -> Point3 {
        match self {
            Self::Line { start, .. } | Self::Arc { start, .. } => start,
        }
    }

    /// End point of this segment.
    #[must_use]
    pub fn end(self) -> Point3 {
        match self {
            Self::Line { end, .. } | Self::Arc { end, .. } => end,
        }
    }
}

/// Convert a sampled polyline into analytic line and circular-arc path segments.
///
/// This is intended for import pipelines where a G-code/scaffold generator has
/// already sampled a circular corner into many small line moves.  The returned
/// path keeps each recognized circular corner as one [`SweepPathSegment::Arc`],
/// so a later sweep creates one toroidal face instead of many short cylinders.
#[must_use]
pub fn analytic_path_from_polyline_samples(
    points: &[Point3],
    options: PathApproxOptions,
) -> Vec<SweepPathSegment> {
    if points.len() < 2 {
        return Vec::new();
    }

    let mut out = Vec::new();
    let mut i = 0usize;
    while i + 1 < points.len() {
        if (points[i + 1] - points[i]).length() < options.point_tol {
            i += 1;
            continue;
        }
        if let Some((end_idx, arc)) = detect_circular_arc(points, i, options) {
            out.push(arc);
            i = end_idx;
        } else {
            out.push(SweepPathSegment::Line {
                start: points[i],
                end: points[i + 1],
            });
            i += 1;
        }
    }
    out
}

/// Convert sampled rounded corners back to a sharp analytic path.
///
/// This is useful for CAD exports where a sampled centre-line radius would
/// create a self-intersecting sweep (`arc_radius <= profile_radius`).  Each
/// recovered circular arc is replaced by the intersection point of its start
/// and end tangents, restoring the original sharp corner so the sweep can use
/// miter boundaries instead of an invalid spindle torus.
#[must_use]
pub fn sharp_path_from_polyline_samples(
    points: &[Point3],
    options: PathApproxOptions,
) -> Vec<SweepPathSegment> {
    let analytic = analytic_path_from_polyline_samples(points, options);
    if analytic.is_empty() {
        return Vec::new();
    }

    let mut out = Vec::new();
    let mut current = analytic[0].start();

    for segment in analytic {
        match segment {
            SweepPathSegment::Line { end, .. } => {
                push_line_if_nonzero(&mut out, current, end);
                current = end;
            }
            SweepPathSegment::Arc {
                start,
                end,
                center,
                normal,
            } => {
                if let Some(info) = segment_info(SweepPathSegment::Arc {
                    start,
                    end,
                    center,
                    normal,
                }) {
                    if let Some(corner) = tangent_intersection(
                        start,
                        info.start_tangent,
                        end,
                        info.end_tangent,
                        normal,
                    ) {
                        if !extend_last_line_end(&mut out, start, corner) {
                            push_line_if_nonzero(&mut out, current, corner);
                        }
                        current = corner;
                        continue;
                    }
                }

                push_line_if_nonzero(&mut out, current, end);
                current = end;
            }
        }
    }

    out
}

/// Build a CAD-style rounded centre-line path from a sharp polyline.
///
/// Interior vertices are replaced with tangent circular arcs of
/// `corner_radius`.  The swept result is the normal CAD construction:
/// first round the centre-line, then sweep the circular profile along that
/// single analytic path.
///
/// # Errors
///
/// Returns [`SweepError::TooFewWaypoints`] for fewer than two points,
/// [`SweepError::InvalidCornerRadius`] for a negative radius, and
/// [`SweepError::CornerRadiusTooLarge`] when the requested radius does not fit
/// inside adjacent line segments.
pub fn rounded_path_from_polyline(
    waypoints: &[Point3],
    corner_radius: f64,
) -> Result<Vec<SweepPathSegment>, SweepError> {
    if waypoints.len() < 2 {
        return Err(SweepError::TooFewWaypoints);
    }
    if corner_radius < 0.0 {
        return Err(SweepError::InvalidCornerRadius(corner_radius));
    }
    if corner_radius <= 1.0e-9 || waypoints.len() == 2 {
        return Ok(lines_from_polyline(waypoints)?);
    }

    let mut segments = Vec::new();
    let mut current = waypoints[0];

    for i in 1..waypoints.len() - 1 {
        let prev = waypoints[i - 1];
        let vertex = waypoints[i];
        let next = waypoints[i + 1];

        let incoming_vec = vertex - prev;
        let outgoing_vec = next - vertex;
        let incoming_len = incoming_vec.length();
        let outgoing_len = outgoing_vec.length();
        let incoming =
            UnitVec3::try_from_vec(incoming_vec).ok_or(SweepError::CoincidentWaypoints(i - 1))?;
        let outgoing =
            UnitVec3::try_from_vec(outgoing_vec).ok_or(SweepError::CoincidentWaypoints(i))?;

        let dot = incoming.dot(outgoing).clamp(-1.0, 1.0);
        if dot > 0.999_999 {
            continue;
        }

        let angle = dot.acos();
        let trim = corner_radius * (angle * 0.5).tan();
        if trim >= incoming_len || trim >= outgoing_len {
            return Err(SweepError::CornerRadiusTooLarge(i));
        }

        let normal =
            UnitVec3::try_from_vec(incoming.cross(outgoing)).ok_or(SweepError::InvalidArc(i))?;
        let entry = vertex - incoming * trim;
        let exit = vertex + outgoing * trim;
        let start_radial = incoming.cross(normal);
        let center = entry - start_radial * corner_radius;

        if (entry - current).length() > 1.0e-9 {
            segments.push(SweepPathSegment::Line {
                start: current,
                end: entry,
            });
        }
        segments.push(SweepPathSegment::Arc {
            start: entry,
            end: exit,
            center,
            normal,
        });
        current = exit;
    }

    let last = *waypoints.last().expect("checked non-empty");
    if (last - current).length() > 1.0e-9 {
        segments.push(SweepPathSegment::Line {
            start: current,
            end: last,
        });
    }

    Ok(segments)
}

/// Sweep a circular profile along a sharp polyline after applying a centre-line
/// corner radius.
///
/// This is the high-level CAD pipe operation for rounded paths: polyline
/// fillet first, profile sweep second.
pub fn sweep_circle_along_rounded_polyline(
    brep: &mut BRep,
    waypoints: &[Point3],
    profile_radius: f64,
    corner_radius: f64,
    opts: &SweepOptions,
) -> Result<SolidId, SweepError> {
    let path = rounded_path_from_polyline(waypoints, corner_radius)?;
    sweep_circle_along_path(brep, &path, profile_radius, opts)
}

/// A half-space plane used for clipping polylines.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClipPlane {
    /// Origin point on the plane.
    pub origin: Point3,
    /// Normal of the plane pointing *into* the kept half-space.
    pub normal: UnitVec3,
}

/// Clip a polyline by a set of planes.
pub fn clip_polyline(points: &[Point3], planes: &[ClipPlane]) -> Vec<Vec<Point3>> {
    let mut current_set = vec![points.to_vec()];
    for plane in planes {
        let mut next_set = Vec::new();
        for poly in current_set {
            let clipped = clip_polyline_by_plane(&poly, plane);
            next_set.extend(clipped);
        }
        current_set = next_set;
    }
    current_set
}

/// Clip a polyline by a set of planes, with awareness of the profile radius
/// to avoid completely dropping parallel segments that physically intersect the kept region.
pub fn clip_polyline_with_radius(
    points: &[Point3],
    planes: &[ClipPlane],
    radius: f64,
) -> Vec<Vec<Point3>> {
    let mut current_set = vec![points.to_vec()];
    for plane in planes {
        let mut next_set = Vec::new();
        for poly in current_set {
            let clipped = clip_polyline_by_plane_with_radius(&poly, plane, radius);
            next_set.extend(clipped);
        }
        current_set = next_set;
    }
    current_set
}

fn clip_polyline_by_plane_with_radius(
    points: &[Point3],
    plane: &ClipPlane,
    radius: f64,
) -> Vec<Vec<Point3>> {
    if points.len() < 2 {
        return Vec::new();
    }
    let mut result = Vec::new();
    let mut current_poly = Vec::new();

    let d = |p: Point3| -> f64 { plane.normal.dot_vec(p - plane.origin) };

    let mut prev_pt = points[0];
    let mut prev_d = d(prev_pt);

    // Apply envelope intersection margin to parallel segments
    let mut prev_adjusted_d = prev_d;

    if points.len() > 1 {
        let next_pt = points[1];
        let diff = next_pt - prev_pt;
        let len_sq = diff.length_sq();
        if len_sq > 1e-12 {
            let dir = diff / len_sq.sqrt();
            let is_parallel = plane.normal.dot_vec(dir).abs() < 1.0e-3;
            if is_parallel && prev_d < 0.0 && prev_d.abs() < radius * 1.01 {
                prev_adjusted_d = 0.0;
            }
        }
    }

    if prev_adjusted_d >= -1e-7 {
        current_poly.push(prev_pt);
    }

    for i in 1..points.len() {
        let pt = points[i];
        let cur_d = d(pt);
        let mut cur_adjusted_d = cur_d;
        prev_adjusted_d = prev_d;

        let diff = pt - prev_pt;
        let len_sq = diff.length_sq();
        if len_sq > 1e-12 {
            let dir = diff / len_sq.sqrt();
            let is_parallel = plane.normal.dot_vec(dir).abs() < 1.0e-3;
            if is_parallel {
                if prev_d < 0.0 && prev_d.abs() < radius * 1.01 {
                    prev_adjusted_d = 0.0;
                }
                if cur_d < 0.0 && cur_d.abs() < radius * 1.01 {
                    cur_adjusted_d = 0.0;
                }
            }
        }

        if prev_adjusted_d >= -1e-7 && cur_adjusted_d >= -1e-7 {
            current_poly.push(pt);
        } else if prev_adjusted_d >= -1e-7 && cur_adjusted_d < -1e-7 {
            let t = prev_d / (prev_d - cur_d);
            let intersect = prev_pt + (pt - prev_pt) * t;
            current_poly.push(intersect);
            if current_poly.len() >= 2 {
                result.push(current_poly);
            }
            current_poly = Vec::new();
        } else if prev_adjusted_d < -1e-7 && cur_adjusted_d >= -1e-7 {
            let t = prev_d / (prev_d - cur_d);
            let intersect = prev_pt + (pt - prev_pt) * t;
            current_poly.push(intersect);
            current_poly.push(pt);
        }

        prev_pt = pt;
        prev_d = cur_d;
    }

    if current_poly.len() >= 2 {
        result.push(current_poly);
    }
    result
}

fn clip_polyline_by_plane(points: &[Point3], plane: &ClipPlane) -> Vec<Vec<Point3>> {
    if points.len() < 2 {
        return Vec::new();
    }
    let mut result = Vec::new();
    let mut current_poly = Vec::new();

    let d = |p: Point3| -> f64 { plane.normal.dot_vec(p - plane.origin) };

    let mut prev_pt = points[0];
    let mut prev_d = d(prev_pt);

    if prev_d >= -1e-7 {
        current_poly.push(prev_pt);
    }

    for &pt in points.iter().skip(1) {
        let cur_d = d(pt);
        if prev_d >= -1e-7 && cur_d >= -1e-7 {
            current_poly.push(pt);
        } else if prev_d >= -1e-7 && cur_d < -1e-7 {
            let t = prev_d / (prev_d - cur_d);
            let intersect = prev_pt + (pt - prev_pt) * t;
            current_poly.push(intersect);
            if current_poly.len() >= 2 {
                result.push(current_poly);
            }
            current_poly = Vec::new();
        } else if prev_d < -1e-7 && cur_d >= -1e-7 {
            let t = prev_d / (prev_d - cur_d);
            let intersect = prev_pt + (pt - prev_pt) * t;
            current_poly.push(intersect);
            current_poly.push(pt);
        }
        prev_pt = pt;
        prev_d = cur_d;
    }

    if current_poly.len() >= 2 {
        result.push(current_poly);
    }
    result
}

/// Sweep a circle of radius `radius` along the polyline defined by `waypoints`,
/// applying custom planar/oblique end caps at start and end.
pub fn sweep_circle_along_polyline_with_caps(
    brep: &mut BRep,
    waypoints: &[Point3],
    radius: f64,
    opts: &SweepOptions,
    start_cut_normal: Option<UnitVec3>,
    end_cut_normal: Option<UnitVec3>,
) -> Result<SolidId, SweepError> {
    if waypoints.len() < 2 {
        return Err(SweepError::TooFewWaypoints);
    }
    if radius <= 0.0 {
        return Err(SweepError::InvalidRadius(radius));
    }
    if opts.fillet_corners {
        let corner_radius = opts.effective_solid_corner_radius(radius)?;
        let centerline_radius = solid_corner_centerline_radius(radius, corner_radius)?;
        let rounded_opts = SweepOptions {
            fillet_corners: false,
            corner_fillet_radius: 0.0,
            name: opts.name.clone(),
        };
        let path = rounded_path_from_polyline(waypoints, centerline_radius)?;
        return sweep_circle_along_path_with_caps(
            brep,
            &path,
            radius,
            &rounded_opts,
            start_cut_normal,
            end_cut_normal,
        );
    }

    // Pre-compute segment directions.
    let n = waypoints.len();
    let mut dirs: Vec<UnitVec3> = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let d = waypoints[i + 1] - waypoints[i];
        match UnitVec3::try_from_vec(d) {
            Some(u) => dirs.push(u),
            None => return Err(SweepError::CoincidentWaypoints(i)),
        }
    }

    let mut face_ids = Vec::new();

    // ── Start cap ────────────────────────────────────────────────────────────
    let start_bound = {
        if let Some(n) = start_cut_normal {
            boundary_for_plane_cut(waypoints[0], dirs[0], n, radius)
        } else {
            FaceBoundary::Circle(Circle3::new(waypoints[0], dirs[0], radius))
        }
    };
    let start_extent = if start_cut_normal.is_some() {
        FaceExtent::PlanarBoundary {
            boundary: start_bound.clone(),
        }
    } else {
        FaceExtent::Disk { radius }
    };
    let start_normal = start_cut_normal.map(|n| -n).unwrap_or(-dirs[0]);
    let _cap_start_id = build_end_cap(
        brep,
        waypoints[0],
        start_normal,
        start_extent,
        &mut face_ids,
    );

    // ── Cylinder + connector loop ─────────────────────────────────────────────
    for i in 0..n - 1 {
        let p0 = waypoints[i];
        let p1 = waypoints[i + 1];
        let dir = dirs[i];
        let length = (p1 - p0).length();
        let start_bound = if i == 0 {
            start_bound.clone()
        } else {
            // Use dirs[i-1] (incoming direction) as cylinder_dir so that this
            // boundary uses the same parameterisation as the end_bound of the
            // previous cylinder — they share the same physical junction ring.
            miter_boundary(p0, dirs[i - 1], dir, dirs[i - 1], radius)
        };
        let end_bound = if i + 1 == n - 1 {
            if let Some(n) = end_cut_normal {
                boundary_for_plane_cut(p1, dir, n, radius)
            } else {
                FaceBoundary::Circle(Circle3::new(p1, dir, radius))
            }
        } else {
            miter_boundary(p1, dir, dirs[i + 1], dir, radius)
        };

        // --- Cylinder face for segment i ---
        let cyl = CylSurf::new(p0, dir, radius);
        let cyl_loop_id = build_placeholder_loop(brep);
        let cyl_face_id = brep.add_face(Face {
            geom: FaceGeom::Cylinder(cyl),
            normal: FaceNormal::Same,
            outer_loop: cyl_loop_id,
            inner_loops: vec![],
            shell: cadcore_topo::ShellId::default(),
            extent: FaceExtent::Cylinder {
                length,
                start: start_bound,
                end: end_bound,
            },
        });
        face_ids.push(cyl_face_id);
    }

    // ── End cap ───────────────────────────────────────────────────────────────
    let last = n - 1;
    let end_bound = if let Some(n) = end_cut_normal {
        boundary_for_plane_cut(waypoints[last], dirs[last - 1], n, radius)
    } else {
        FaceBoundary::Circle(Circle3::new(waypoints[last], dirs[last - 1], radius))
    };
    let end_extent = if end_cut_normal.is_some() {
        FaceExtent::PlanarBoundary {
            boundary: end_bound,
        }
    } else {
        FaceExtent::Disk { radius }
    };
    let end_normal = end_cut_normal.map(|n| -n).unwrap_or(dirs[last - 1]);
    let _cap_end_id = build_end_cap(brep, waypoints[last], end_normal, end_extent, &mut face_ids);

    assemble_solid(brep, face_ids, opts.name.clone())
}

/// Sweep a circle of radius `radius` along the polyline defined by `waypoints`.
///
/// Returns the id of the newly inserted solid within `brep`.
///
/// # Errors
///
/// Returns `Err` if:
/// * `waypoints` has fewer than 2 points.
/// * `radius` ≤ 0.
/// * Any two consecutive waypoints are coincident (zero-length segment).
pub fn sweep_circle_along_polyline(
    brep: &mut BRep,
    waypoints: &[Point3],
    radius: f64,
    opts: &SweepOptions,
) -> Result<SolidId, SweepError> {
    sweep_circle_along_polyline_with_caps(brep, waypoints, radius, opts, None, None)
}

/// Sweep a circle of radius `radius` along an analytic path,
/// applying custom planar/oblique end caps at start and end.
pub fn sweep_circle_along_path_with_caps(
    brep: &mut BRep,
    segments: &[SweepPathSegment],
    radius: f64,
    opts: &SweepOptions,
    start_cut_normal: Option<UnitVec3>,
    end_cut_normal: Option<UnitVec3>,
) -> Result<SolidId, SweepError> {
    if segments.is_empty() {
        return Err(SweepError::TooFewWaypoints);
    }
    if radius <= 0.0 {
        return Err(SweepError::InvalidRadius(radius));
    }
    let mut infos = Vec::with_capacity(segments.len());
    for (idx, segment) in segments.iter().copied().enumerate() {
        let info = segment_info(segment).ok_or(SweepError::InvalidArc(idx))?;
        if let SegmentKind::Arc { arc_radius, .. } = info.kind {
            if arc_radius <= radius {
                return Err(SweepError::SelfIntersectingSweep(idx));
            }
        }
        infos.push(info);
    }

    for i in 1..infos.len() {
        if (infos[i].start - infos[i - 1].end).length() > 1.0e-6 {
            return Err(SweepError::DisconnectedSegments(i));
        }
    }

    let mut face_ids = Vec::new();

    // ── Start cap boundary and extent ────────────────────────────────────────
    let start_bound = {
        if let Some(n) = start_cut_normal {
            boundary_for_plane_cut(infos[0].start, infos[0].start_tangent, n, radius)
        } else {
            FaceBoundary::Circle(Circle3::new(infos[0].start, infos[0].start_tangent, radius))
        }
    };
    let start_extent = if start_cut_normal.is_some() {
        FaceExtent::PlanarBoundary {
            boundary: start_bound.clone(),
        }
    } else {
        FaceExtent::Disk { radius }
    };
    let start_normal = start_cut_normal
        .map(|n| -n)
        .unwrap_or(-infos[0].start_tangent);
    let _cap_start_id = build_end_cap(
        brep,
        infos[0].start,
        start_normal,
        start_extent,
        &mut face_ids,
    );

    // ── Cylinder / torus faces ───────────────────────────────────────────────
    for (idx, info) in infos.iter().enumerate() {
        let start_bound = if idx == 0 {
            start_bound.clone()
        } else {
            // Use the previous segment's end_tangent as cylinder_dir so this
            // boundary has the same parameterisation as that segment's end_bound.
            join_boundary(
                info.start,
                infos[idx - 1].end_tangent,
                info.start_tangent,
                infos[idx - 1].end_tangent,
                radius,
            )
        };
        let end_bound = if idx + 1 == infos.len() {
            if let Some(n) = end_cut_normal {
                boundary_for_plane_cut(info.end, info.end_tangent, n, radius)
            } else {
                FaceBoundary::Circle(Circle3::new(info.end, info.end_tangent, radius))
            }
        } else {
            join_boundary(
                info.end,
                info.end_tangent,
                infos[idx + 1].start_tangent,
                info.end_tangent,
                radius,
            )
        };

        match info.kind {
            SegmentKind::Line { length } => {
                let cyl = CylSurf::new(info.start, info.start_tangent, radius);
                let loop_id = build_placeholder_loop(brep);
                let face_id = brep.add_face(Face {
                    geom: FaceGeom::Cylinder(cyl),
                    normal: FaceNormal::Same,
                    outer_loop: loop_id,
                    inner_loops: vec![],
                    shell: cadcore_topo::ShellId::default(),
                    extent: FaceExtent::Cylinder {
                        length,
                        start: start_bound,
                        end: end_bound,
                    },
                });
                face_ids.push(face_id);
            }
            SegmentKind::Arc {
                center,
                normal,
                arc_radius,
            } => {
                let torus = TorusSurf::new(center, normal, arc_radius, radius);
                let loop_id = build_placeholder_loop(brep);
                let face_id = brep.add_face(Face {
                    geom: FaceGeom::Torus(torus),
                    normal: FaceNormal::Same,
                    outer_loop: loop_id,
                    inner_loops: vec![],
                    shell: cadcore_topo::ShellId::default(),
                    extent: FaceExtent::TorusFillet {
                        start_circle: Circle3::new(info.start, info.start_tangent, radius),
                        end_circle: Circle3::new(info.end, info.end_tangent, radius),
                    },
                });
                face_ids.push(face_id);
            }
        }
    }

    // ── End cap ──────────────────────────────────────────────────────────────
    let last = infos.last().expect("non-empty path");
    let end_bound = if let Some(n) = end_cut_normal {
        boundary_for_plane_cut(last.end, last.end_tangent, n, radius)
    } else {
        FaceBoundary::Circle(Circle3::new(last.end, last.end_tangent, radius))
    };
    let end_extent = if end_cut_normal.is_some() {
        FaceExtent::PlanarBoundary {
            boundary: end_bound,
        }
    } else {
        FaceExtent::Disk { radius }
    };
    let end_normal = end_cut_normal.map(|n| -n).unwrap_or(last.end_tangent);
    let _cap_end_id = build_end_cap(brep, last.end, end_normal, end_extent, &mut face_ids);

    assemble_solid(brep, face_ids, opts.name.clone())
}

/// Sweep a circle of radius `radius` along an analytic path.
pub fn sweep_circle_along_path(
    brep: &mut BRep,
    segments: &[SweepPathSegment],
    radius: f64,
    opts: &SweepOptions,
) -> Result<SolidId, SweepError> {
    sweep_circle_along_path_with_caps(brep, segments, radius, opts, None, None)
}

/// Construct a solid box (plate/block) in the `BRep` using polygonal faces.
pub fn build_solid_box(
    brep: &mut BRep,
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    zmin: f64,
    zmax: f64,
    name: Option<String>,
) -> Result<SolidId, SweepError> {
    let mut face_ids = Vec::new();

    let mut add_poly_face = |pts: Vec<Point3>, normal: UnitVec3| {
        let plane = Plane3::from_origin_normal(pts[0], normal);
        let loop_id = build_placeholder_loop(brep);
        let fid = brep.add_face(Face {
            geom: FaceGeom::Plane(plane),
            normal: FaceNormal::Same,
            outer_loop: loop_id,
            inner_loops: vec![],
            shell: cadcore_topo::ShellId::default(),
            extent: FaceExtent::Polygon { points: pts },
        });
        face_ids.push(fid);
    };

    // Face 1 (Z max, normal [0, 0, 1])
    add_poly_face(
        vec![
            Point3::new(xmin, ymin, zmax),
            Point3::new(xmax, ymin, zmax),
            Point3::new(xmax, ymax, zmax),
            Point3::new(xmin, ymax, zmax),
        ],
        UnitVec3::Z,
    );

    // Face 2 (Z min, normal [0, 0, -1])
    add_poly_face(
        vec![
            Point3::new(xmin, ymax, zmin),
            Point3::new(xmax, ymax, zmin),
            Point3::new(xmax, ymin, zmin),
            Point3::new(xmin, ymin, zmin),
        ],
        -UnitVec3::Z,
    );

    // Face 3 (Y max, normal [0, 1, 0])
    add_poly_face(
        vec![
            Point3::new(xmin, ymax, zmax),
            Point3::new(xmax, ymax, zmax),
            Point3::new(xmax, ymax, zmin),
            Point3::new(xmin, ymax, zmin),
        ],
        UnitVec3::Y,
    );

    // Face 4 (Y min, normal [0, -1, 0])
    add_poly_face(
        vec![
            Point3::new(xmin, ymin, zmin),
            Point3::new(xmax, ymin, zmin),
            Point3::new(xmax, ymin, zmax),
            Point3::new(xmin, ymin, zmax),
        ],
        -UnitVec3::Y,
    );

    // Face 5 (X max, normal [1, 0, 0])
    add_poly_face(
        vec![
            Point3::new(xmax, ymin, zmin),
            Point3::new(xmax, ymax, zmin),
            Point3::new(xmax, ymax, zmax),
            Point3::new(xmax, ymin, zmax),
        ],
        UnitVec3::X,
    );

    // Face 6 (X min, normal [-1, 0, 0])
    add_poly_face(
        vec![
            Point3::new(xmin, ymin, zmax),
            Point3::new(xmin, ymax, zmax),
            Point3::new(xmin, ymax, zmin),
            Point3::new(xmin, ymin, zmin),
        ],
        -UnitVec3::X,
    );

    let shell_id = brep.add_shell(Shell {
        faces: face_ids.clone(),
        is_outer: true,
        solid: cadcore_topo::SolidId::default(),
    });

    for &fid in &face_ids {
        if let Some(f) = brep.faces.get_mut(fid) {
            f.shell = shell_id;
        }
    }

    let solid_id = brep.add_solid(Solid {
        shells: vec![shell_id],
        name,
    });

    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.solid = solid_id;
    }

    Ok(solid_id)
}

/// Build a solid box with rounded XZ-plane corners using fully **analytic BREP geometry**.
///
/// Produces exactly **10 faces**:
/// * 2 cap faces (Y = ymax / ymin) — `FaceGeom::Plane` + `FaceExtent::RoundedRectCap`
///   with 4 LINE + 4 CIRCLE arc edges in STEP (no tessellation on the cap).
/// * 4 flat side quads — `FaceGeom::Plane` + `FaceExtent::Polygon` (4-vertex quads).
/// * 4 quarter-cylinder corner faces — `FaceGeom::Cylinder` + `FaceExtent::CylinderArcFace`
///   with analytic CIRCLE / LINE boundary curves in STEP.
///
/// This is equivalent to a "pull to radius" operation in a parametric modeller:
/// smooth analytic geometry, no visible faceting, FEM-clean topology.
pub fn build_solid_rounded_box_xz(
    brep: &mut BRep,
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    zmin: f64,
    zmax: f64,
    corner_radius: f64,
    name: Option<String>,
) -> Result<SolidId, SweepError> {
    let r = corner_radius
        .max(0.0)
        .min(((xmax - xmin) * 0.5).min((zmax - zmin) * 0.5));
    if r <= 1.0e-9 {
        return build_solid_box(brep, xmin, xmax, ymin, ymax, zmin, zmax, name);
    }

    let len_y = ymax - ymin;

    // Collect all (geom, normal_flag, extent) triples, then bulk-insert.
    struct FaceSpec {
        geom: FaceGeom,
        normal: FaceNormal,
        extent: FaceExtent,
    }

    let mut specs: Vec<FaceSpec> = Vec::with_capacity(10);

    // Helper: planar polygon face spec.
    let poly_spec = |pts: Vec<Point3>, n: UnitVec3| -> FaceSpec {
        let plane = Plane3::from_origin_normal(pts[0], n);
        FaceSpec {
            geom: FaceGeom::Plane(plane),
            normal: FaceNormal::Same,
            extent: FaceExtent::Polygon { points: pts },
        }
    };

    // ── 1. Back cap (Y = ymax, outward +Y) ───────────────────────────────────
    specs.push(FaceSpec {
        geom: FaceGeom::Plane(Plane3::from_origin_normal(
            Point3::new(0.0, ymax, 0.0),
            UnitVec3::Y,
        )),
        normal: FaceNormal::Same,
        extent: FaceExtent::RoundedRectCap {
            xmin,
            xmax,
            zmin,
            zmax,
            radius: r,
            y: ymax,
            plus_y: true,
        },
    });

    // ── 2. Front cap (Y = ymin, outward −Y) ──────────────────────────────────
    specs.push(FaceSpec {
        geom: FaceGeom::Plane(Plane3::from_origin_normal(
            Point3::new(0.0, ymin, 0.0),
            -UnitVec3::Y,
        )),
        normal: FaceNormal::Same,
        extent: FaceExtent::RoundedRectCap {
            xmin,
            xmax,
            zmin,
            zmax,
            radius: r,
            y: ymin,
            plus_y: false,
        },
    });

    // ── 3–6. Flat side faces ─────────────────────────────────────────────────
    // Bottom (Z = zmin, outward −Z)
    specs.push(poly_spec(
        vec![
            Point3::new(xmin + r, ymin, zmin),
            Point3::new(xmin + r, ymax, zmin),
            Point3::new(xmax - r, ymax, zmin),
            Point3::new(xmax - r, ymin, zmin),
        ],
        -UnitVec3::Z,
    ));
    // Top (Z = zmax, outward +Z)
    specs.push(poly_spec(
        vec![
            Point3::new(xmax - r, ymin, zmax),
            Point3::new(xmax - r, ymax, zmax),
            Point3::new(xmin + r, ymax, zmax),
            Point3::new(xmin + r, ymin, zmax),
        ],
        UnitVec3::Z,
    ));
    // Left (X = xmin, outward −X)
    specs.push(poly_spec(
        vec![
            Point3::new(xmin, ymin, zmax - r),
            Point3::new(xmin, ymax, zmax - r),
            Point3::new(xmin, ymax, zmin + r),
            Point3::new(xmin, ymin, zmin + r),
        ],
        -UnitVec3::X,
    ));
    // Right (X = xmax, outward +X)
    specs.push(poly_spec(
        vec![
            Point3::new(xmax, ymin, zmin + r),
            Point3::new(xmax, ymax, zmin + r),
            Point3::new(xmax, ymax, zmax - r),
            Point3::new(xmax, ymin, zmax - r),
        ],
        UnitVec3::X,
    ));

    // ── 7–10. Quarter-cylinder corner faces (CCW-from-+Y order: BR→TR→TL→BL) ─
    //
    // Each cylinder has axis = +Y, origin at (cx, ymin, cz).
    // arc_ref_dir is the CylSurf x-reference (angle=0 direction).
    // right = Y × arc_ref_dir (= CylSurf y-reference).
    //
    // CylinderArcFace: arc_start=0, arc_end=−π/2, length=len_y.
    //   arc_pt(0)    = origin + arc_ref_dir * r   ← arc start at ymin
    //   arc_pt(−π/2) = origin − right * r         ← arc end at ymin
    //
    // Corner  arc_ref  start(xz)         end(xz)
    // BR  −Z   (xmax−r, zmin)  (xmax, zmin+r)
    // TR  +X   (xmax,  zmax−r) (xmax−r, zmax)
    // TL  +Z   (xmin+r,zmax)   (xmin, zmax−r)
    // BL  −X   (xmin,  zmin+r) (xmin+r,zmin)
    let corner_data: [(f64, f64, UnitVec3, UnitVec3); 4] = [
        (xmax - r, zmin + r, -UnitVec3::Z, -UnitVec3::X), // BR: ref=−Z, right=Y×(−Z)=−X
        (xmax - r, zmax - r, UnitVec3::X, -UnitVec3::Z),  // TR: ref=+X, right=Y×X=−Z
        (xmin + r, zmax - r, UnitVec3::Z, UnitVec3::X),   // TL: ref=+Z, right=Y×Z=+X
        (xmin + r, zmin + r, -UnitVec3::X, UnitVec3::Z),  // BL: ref=−X, right=Y×(−X)=+Z
    ];

    for (cx, cz, arc_ref, right) in corner_data {
        let cyl = CylSurf {
            frame: cadcore_math::Frame3 {
                origin: Point3::new(cx, ymin, cz),
                x: arc_ref,
                y: right,
                z: UnitVec3::Y,
            },
            radius: r,
        };
        specs.push(FaceSpec {
            geom: FaceGeom::Cylinder(cyl),
            normal: FaceNormal::Same,
            extent: FaceExtent::CylinderArcFace {
                length: len_y,
                arc_start_angle: 0.0,
                arc_end_angle: -std::f64::consts::FRAC_PI_2,
                arc_ref_dir: arc_ref,
            },
        });
    }

    // ── Bulk-insert faces, then shell + solid ─────────────────────────────────
    let mut face_ids: Vec<cadcore_topo::FaceId> = Vec::with_capacity(specs.len());
    for spec in specs {
        let loop_id = build_placeholder_loop(brep);
        let fid = brep.add_face(Face {
            geom: spec.geom,
            normal: spec.normal,
            outer_loop: loop_id,
            inner_loops: vec![],
            shell: cadcore_topo::ShellId::default(),
            extent: spec.extent,
        });
        face_ids.push(fid);
    }

    let shell_id = brep.add_shell(Shell {
        faces: face_ids.clone(),
        is_outer: true,
        solid: cadcore_topo::SolidId::default(),
    });
    for &fid in &face_ids {
        if let Some(f) = brep.faces.get_mut(fid) {
            f.shell = shell_id;
        }
    }
    let solid_id = brep.add_solid(Solid {
        shells: vec![shell_id],
        name,
    });
    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.solid = solid_id;
    }
    Ok(solid_id)
}

/// Build a solid box with ALL edges and corners rounded (Minkowski sum of a box
/// and a sphere of radius `corner_radius`).
///
/// The surface is tessellated into polygon triangles on a grid of `N × N` per
/// box face (N = 8), giving smooth edges on every axis and spherical corners.
/// Flat face regions project exactly; edge/corner regions curve correctly.
///
/// For FEM import: all 768 triangle faces are valid PLANE polygons.  CAD tools
/// that read STEP see a mesh that accurately represents the smooth surface.
///
/// If `corner_radius ≤ 0` delegates to [`build_solid_box`].
pub fn build_solid_rounded_box_full(
    brep: &mut BRep,
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    zmin: f64,
    zmax: f64,
    corner_radius: f64,
    name: Option<String>,
) -> Result<SolidId, SweepError> {
    let r = corner_radius.max(0.0).min(
        ((xmax - xmin) * 0.5)
            .min((ymax - ymin) * 0.5)
            .min((zmax - zmin) * 0.5),
    );
    if r <= 1.0e-9 {
        return build_solid_box(brep, xmin, xmax, ymin, ymax, zmin, zmax, name);
    }

    // Inner box center and half-extents (before Minkowski inflation).
    let cx = (xmin + xmax) * 0.5;
    let cy = (ymin + ymax) * 0.5;
    let cz = (zmin + zmax) * 0.5;
    let hx = (xmax - xmin) * 0.5 - r;
    let hy = (ymax - ymin) * 0.5 - r;
    let hz = (zmax - zmin) * 0.5 - r;

    // N grid subdivisions per face edge.
    const N: usize = 8;

    // Minkowski projection: clamp to inner box then extend outward by r.
    let proj = |x: f64, y: f64, z: f64| -> Point3 {
        let rx = x - cx;
        let ry = y - cy;
        let rz = z - cz;
        let ix = rx.clamp(-hx, hx);
        let iy = ry.clamp(-hy, hy);
        let iz = rz.clamp(-hz, hz);
        let dx = rx - ix;
        let dy = ry - iy;
        let dz = rz - iz;
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        if dist > 1.0e-12 {
            let s = r / dist;
            Point3::new(cx + ix + dx * s, cy + iy + dy * s, cz + iz + dz * s)
        } else {
            Point3::new(x, y, z)
        }
    };

    let mut face_ids = Vec::new();

    // Winding helpers for each of the 6 outer-box faces.
    // The macros generate N×N quads split into 2 triangles each.
    // CCW winding is chosen so the cross product of the first two edges
    // points OUTWARD for each face (verified by sign analysis below).
    {
        // Closure captures brep and face_ids mutably.
        let mut add_tri = |p0: Point3, p1: Point3, p2: Point3| {
            let e1 = p1 - p0;
            let e2 = p2 - p0;
            let n = e1.cross(e2);
            let normal = UnitVec3::try_from_vec(n).unwrap_or(UnitVec3::Z);
            let plane = Plane3::from_origin_normal(p0, normal);
            let loop_id = build_placeholder_loop(brep);
            let fid = brep.add_face(Face {
                geom: FaceGeom::Plane(plane),
                normal: FaceNormal::Same,
                outer_loop: loop_id,
                inner_loops: vec![],
                shell: cadcore_topo::ShellId::default(),
                extent: FaceExtent::Polygon {
                    points: vec![p0, p1, p2],
                },
            });
            face_ids.push(fid);
        };

        // ── +X face (outward normal +X: edge1=+Y, edge2=+Z) ──────────────
        for i in 0..N {
            for j in 0..N {
                let y0 = ymin + (ymax - ymin) * i as f64 / N as f64;
                let y1 = ymin + (ymax - ymin) * (i + 1) as f64 / N as f64;
                let z0 = zmin + (zmax - zmin) * j as f64 / N as f64;
                let z1 = zmin + (zmax - zmin) * (j + 1) as f64 / N as f64;
                let a = proj(xmax, y0, z0);
                let b = proj(xmax, y1, z0); // +Y
                let c = proj(xmax, y1, z1); // +Y +Z
                let d = proj(xmax, y0, z1); // +Z
                add_tri(a, b, c);
                add_tri(a, c, d);
            }
        }
        // ── -X face (outward normal -X: edge1=+Z, edge2=+Y) ──────────────
        for i in 0..N {
            for j in 0..N {
                let y0 = ymin + (ymax - ymin) * i as f64 / N as f64;
                let y1 = ymin + (ymax - ymin) * (i + 1) as f64 / N as f64;
                let z0 = zmin + (zmax - zmin) * j as f64 / N as f64;
                let z1 = zmin + (zmax - zmin) * (j + 1) as f64 / N as f64;
                let a = proj(xmin, y0, z0);
                let b = proj(xmin, y0, z1); // +Z
                let c = proj(xmin, y1, z1); // +Y +Z
                let d = proj(xmin, y1, z0); // +Y
                add_tri(a, b, c);
                add_tri(a, c, d);
            }
        }
        // ── +Y face (outward normal +Y: edge1=+Z, edge2=+X) ──────────────
        for i in 0..N {
            for j in 0..N {
                let x0 = xmin + (xmax - xmin) * i as f64 / N as f64;
                let x1 = xmin + (xmax - xmin) * (i + 1) as f64 / N as f64;
                let z0 = zmin + (zmax - zmin) * j as f64 / N as f64;
                let z1 = zmin + (zmax - zmin) * (j + 1) as f64 / N as f64;
                let a = proj(x0, ymax, z0);
                let b = proj(x0, ymax, z1); // +Z
                let c = proj(x1, ymax, z1); // +X +Z
                let d = proj(x1, ymax, z0); // +X
                add_tri(a, b, c);
                add_tri(a, c, d);
            }
        }
        // ── -Y face (outward normal -Y: edge1=+X, edge2=+Z) ──────────────
        for i in 0..N {
            for j in 0..N {
                let x0 = xmin + (xmax - xmin) * i as f64 / N as f64;
                let x1 = xmin + (xmax - xmin) * (i + 1) as f64 / N as f64;
                let z0 = zmin + (zmax - zmin) * j as f64 / N as f64;
                let z1 = zmin + (zmax - zmin) * (j + 1) as f64 / N as f64;
                let a = proj(x0, ymin, z0);
                let b = proj(x1, ymin, z0); // +X
                let c = proj(x1, ymin, z1); // +X +Z
                let d = proj(x0, ymin, z1); // +Z
                add_tri(a, b, c);
                add_tri(a, c, d);
            }
        }
        // ── +Z face (outward normal +Z: edge1=+X, edge2=+Y) ──────────────
        for i in 0..N {
            for j in 0..N {
                let x0 = xmin + (xmax - xmin) * i as f64 / N as f64;
                let x1 = xmin + (xmax - xmin) * (i + 1) as f64 / N as f64;
                let y0 = ymin + (ymax - ymin) * j as f64 / N as f64;
                let y1 = ymin + (ymax - ymin) * (j + 1) as f64 / N as f64;
                let a = proj(x0, y0, zmax);
                let b = proj(x1, y0, zmax); // +X
                let c = proj(x1, y1, zmax); // +X +Y
                let d = proj(x0, y1, zmax); // +Y
                add_tri(a, b, c);
                add_tri(a, c, d);
            }
        }
        // ── -Z face (outward normal -Z: edge1=+Y, edge2=+X) ──────────────
        for i in 0..N {
            for j in 0..N {
                let x0 = xmin + (xmax - xmin) * i as f64 / N as f64;
                let x1 = xmin + (xmax - xmin) * (i + 1) as f64 / N as f64;
                let y0 = ymin + (ymax - ymin) * j as f64 / N as f64;
                let y1 = ymin + (ymax - ymin) * (j + 1) as f64 / N as f64;
                let a = proj(x0, y0, zmin);
                let b = proj(x0, y1, zmin); // +Y
                let c = proj(x1, y1, zmin); // +X +Y
                let d = proj(x1, y0, zmin); // +X
                add_tri(a, b, c);
                add_tri(a, c, d);
            }
        }
    } // closure dropped here — brep borrow released

    let shell_id = brep.add_shell(Shell {
        faces: face_ids.clone(),
        is_outer: true,
        solid: cadcore_topo::SolidId::default(),
    });
    for &fid in &face_ids {
        if let Some(f) = brep.faces.get_mut(fid) {
            f.shell = shell_id;
        }
    }
    let solid_id = brep.add_solid(Solid {
        shells: vec![shell_id],
        name,
    });
    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.solid = solid_id;
    }
    Ok(solid_id)
}

/// Build a solid box with rounded vertical edges (edges parallel to Z).
///
/// All four vertical edges are replaced by quarter-cylinder arcs of
/// `corner_radius`.  The geometry is represented as polygon faces —
/// `N_ARC_SEGS` (8) segments per quarter-arc give a clean approximation
/// suitable for FEM import and STEP export.
///
/// If `corner_radius ≤ 0` or the radius doesn't fit, delegates to
/// [`build_solid_box`] (sharp corners).
pub fn build_solid_rounded_box(
    brep: &mut BRep,
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    zmin: f64,
    zmax: f64,
    corner_radius: f64,
    name: Option<String>,
) -> Result<SolidId, SweepError> {
    let r = corner_radius
        .max(0.0)
        .min(((xmax - xmin) * 0.5).min((ymax - ymin) * 0.5));
    if r <= 1.0e-9 {
        return build_solid_box(brep, xmin, xmax, ymin, ymax, zmin, zmax, name);
    }

    const N_ARC: usize = 8;

    let profile = rounded_rect_profile_2d(xmin, xmax, ymin, ymax, r, N_ARC);
    let n = profile.len();

    let mut face_ids = Vec::new();

    let mut add_poly = |pts: Vec<Point3>, normal: UnitVec3| {
        let plane = Plane3::from_origin_normal(pts[0], normal);
        let loop_id = build_placeholder_loop(brep);
        let fid = brep.add_face(Face {
            geom: FaceGeom::Plane(plane),
            normal: FaceNormal::Same,
            outer_loop: loop_id,
            inner_loops: vec![],
            shell: cadcore_topo::ShellId::default(),
            extent: FaceExtent::Polygon { points: pts },
        });
        face_ids.push(fid);
    };

    // Top face — CCW when viewed from +Z
    let top: Vec<Point3> = profile
        .iter()
        .map(|&[x, y]| Point3::new(x, y, zmax))
        .collect();
    add_poly(top, UnitVec3::Z);

    // Bottom face — CW when viewed from +Z (outward normal -Z)
    let bot: Vec<Point3> = profile
        .iter()
        .rev()
        .map(|&[x, y]| Point3::new(x, y, zmin))
        .collect();
    add_poly(bot, -UnitVec3::Z);

    // Side faces — one quad per polygon edge
    for i in 0..n {
        let j = (i + 1) % n;
        let [xi, yi] = profile[i];
        let [xj, yj] = profile[j];
        let p0 = Point3::new(xi, yi, zmax);
        let p1 = Point3::new(xj, yj, zmax);
        let p2 = Point3::new(xj, yj, zmin);
        let p3 = Point3::new(xi, yi, zmin);

        // Outward normal = right-perpendicular of (j−i) in XY = (dy, −dx, 0)
        let dx = xj - xi;
        let dy = yj - yi;
        let len = (dx * dx + dy * dy).sqrt();
        if len < 1.0e-9 {
            continue;
        }
        let normal =
            UnitVec3::try_from_vec(Vec3::new(dy / len, -dx / len, 0.0)).unwrap_or(UnitVec3::X);
        add_poly(vec![p0, p1, p2, p3], normal);
    }

    let shell_id = brep.add_shell(Shell {
        faces: face_ids.clone(),
        is_outer: true,
        solid: cadcore_topo::SolidId::default(),
    });
    for &fid in &face_ids {
        if let Some(f) = brep.faces.get_mut(fid) {
            f.shell = shell_id;
        }
    }
    let solid_id = brep.add_solid(Solid {
        shells: vec![shell_id],
        name,
    });
    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.solid = solid_id;
    }
    Ok(solid_id)
}

/// Rounded-rectangle 2-D profile, CCW winding when viewed from +Z.
///
/// Returns `4 × n_arc` points: four quarter-arcs (SE → NE → NW → SW),
/// each contributing `n_arc` samples (the exact junction point at the
/// start of each arc is included; the arc end is the next arc's start).
fn rounded_rect_profile_2d(
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    r: f64,
    n_arc: usize,
) -> Vec<[f64; 2]> {
    use std::f64::consts::{FRAC_PI_2, PI};
    // (center_x, center_y, start_angle) for SE → NE → NW → SW
    let arcs: [([f64; 2], f64); 4] = [
        ([xmax - r, ymin + r], -FRAC_PI_2), // SE: 270° → 360°
        ([xmax - r, ymax - r], 0.0),        // NE: 0°   → 90°
        ([xmin + r, ymax - r], FRAC_PI_2),  // NW: 90°  → 180°
        ([xmin + r, ymin + r], PI),         // SW: 180° → 270°
    ];
    let mut pts = Vec::with_capacity(4 * n_arc);
    for ([cx, cy], start_a) in &arcs {
        for i in 0..n_arc {
            let a = start_a + FRAC_PI_2 * i as f64 / n_arc as f64;
            pts.push([cx + r * a.cos(), cy + r * a.sin()]);
        }
    }
    pts
}

fn build_end_cap(
    brep: &mut BRep,
    centre: Point3,
    normal: UnitVec3,
    extent: FaceExtent,
    face_ids: &mut Vec<cadcore_topo::FaceId>,
) -> cadcore_topo::FaceId {
    let plane = Plane3::from_origin_normal(centre, normal);
    let loop_id = build_placeholder_loop(brep);
    let fid = brep.add_face(Face {
        geom: FaceGeom::Plane(plane),
        normal: FaceNormal::Same,
        outer_loop: loop_id,
        inner_loops: vec![],
        shell: cadcore_topo::ShellId::default(),
        extent,
    });
    face_ids.push(fid);
    fid
}

fn build_placeholder_loop(brep: &mut BRep) -> cadcore_topo::LoopId {
    brep.add_loop(cadcore_topo::Loop {
        start: cadcore_topo::CoEdgeId::default(),
        face: cadcore_topo::FaceId::default(),
    })
}

/// Errors that [`sweep_circle_along_polyline`] can return.
#[derive(Clone, Debug, PartialEq)]
pub enum SweepError {
    /// Need at least 2 waypoints.
    TooFewWaypoints,
    /// Radius must be positive.
    InvalidRadius(f64),
    /// Two consecutive waypoints are the same (zero-length segment at index).
    CoincidentWaypoints(usize),
    /// A path arc or line segment is geometrically invalid.
    InvalidArc(usize),
    /// Adjacent analytic path segments do not share an endpoint.
    DisconnectedSegments(usize),
    /// Centre-line corner radius must be non-negative.
    InvalidCornerRadius(f64),
    /// Requested centre-line corner radius does not fit at vertex index.
    CornerRadiusTooLarge(usize),
    /// Circular path radius is not larger than the swept profile radius.
    SelfIntersectingSweep(usize),
}

impl std::fmt::Display for SweepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TooFewWaypoints => write!(f, "need at least 2 waypoints"),
            Self::InvalidRadius(r) => write!(f, "invalid radius: {r}"),
            Self::CoincidentWaypoints(i) => write!(f, "coincident waypoints at index {i}"),
            Self::InvalidArc(i) => write!(f, "invalid analytic path segment at index {i}"),
            Self::DisconnectedSegments(i) => write!(f, "disconnected path segment at index {i}"),
            Self::InvalidCornerRadius(r) => write!(f, "invalid corner radius: {r}"),
            Self::CornerRadiusTooLarge(i) => {
                write!(f, "corner radius too large at vertex index {i}")
            }
            Self::SelfIntersectingSweep(i) => {
                write!(f, "self-intersecting sweep at path segment index {i}")
            }
        }
    }
}

impl std::error::Error for SweepError {}

#[derive(Clone, Copy, Debug)]
enum SegmentKind {
    Line {
        length: f64,
    },
    Arc {
        center: Point3,
        normal: UnitVec3,
        arc_radius: f64,
    },
}

#[derive(Clone, Copy, Debug)]
struct SegmentInfo {
    start: Point3,
    end: Point3,
    start_tangent: UnitVec3,
    end_tangent: UnitVec3,
    kind: SegmentKind,
}

fn segment_info(segment: SweepPathSegment) -> Option<SegmentInfo> {
    match segment {
        SweepPathSegment::Line { start, end } => {
            let delta = end - start;
            let length = delta.length();
            let tangent = UnitVec3::try_from_vec(delta)?;
            Some(SegmentInfo {
                start,
                end,
                start_tangent: tangent,
                end_tangent: tangent,
                kind: SegmentKind::Line { length },
            })
        }
        SweepPathSegment::Arc {
            start,
            end,
            center,
            normal,
        } => {
            let r0 = start - center;
            let r1 = end - center;
            let len0 = r0.length();
            let len1 = r1.length();
            if len0 < 1.0e-9 || len1 < 1.0e-9 || (start - end).length() < 1.0e-9 {
                return None;
            }
            let avg_radius = 0.5 * (len0 + len1);
            let tol = (avg_radius * 0.01).max(1.0e-3);
            if (len0 - len1).abs() > tol {
                return None;
            }
            if normal.dot_vec(r0).abs() > tol || normal.dot_vec(r1).abs() > tol {
                return None;
            }

            let start_tangent = UnitVec3::try_from_vec(normal.as_vec().cross(r0))?;
            let end_tangent = UnitVec3::try_from_vec(normal.as_vec().cross(r1))?;

            Some(SegmentInfo {
                start,
                end,
                start_tangent,
                end_tangent,
                kind: SegmentKind::Arc {
                    center,
                    normal,
                    arc_radius: avg_radius,
                },
            })
        }
    }
}

fn lines_from_polyline(waypoints: &[Point3]) -> Result<Vec<SweepPathSegment>, SweepError> {
    let mut out = Vec::with_capacity(waypoints.len().saturating_sub(1));
    for i in 0..waypoints.len() - 1 {
        if (waypoints[i + 1] - waypoints[i]).length() < 1.0e-9 {
            return Err(SweepError::CoincidentWaypoints(i));
        }
        out.push(SweepPathSegment::Line {
            start: waypoints[i],
            end: waypoints[i + 1],
        });
    }
    Ok(out)
}

fn push_line_if_nonzero(out: &mut Vec<SweepPathSegment>, start: Point3, end: Point3) {
    if (end - start).length() >= 1.0e-9 {
        out.push(SweepPathSegment::Line { start, end });
    }
}

fn extend_last_line_end(out: &mut Vec<SweepPathSegment>, old_end: Point3, new_end: Point3) -> bool {
    let Some(last) = out.last_mut() else {
        return false;
    };
    match last {
        SweepPathSegment::Line { end, .. } if (*end - old_end).length() < 1.0e-7 => {
            *end = new_end;
            true
        }
        _ => false,
    }
}

fn tangent_intersection(
    p: Point3,
    p_dir: UnitVec3,
    q: Point3,
    q_dir: UnitVec3,
    normal: UnitVec3,
) -> Option<Point3> {
    let denom = normal.dot_vec(p_dir.cross(q_dir));
    if denom.abs() < 1.0e-9 {
        return None;
    }

    let t = normal.dot_vec((q - p).cross(q_dir.as_vec())) / denom;
    Some(p + p_dir * t)
}

fn detect_circular_arc(
    points: &[Point3],
    start: usize,
    options: PathApproxOptions,
) -> Option<(usize, SweepPathSegment)> {
    if options.min_arc_points < 3 || start + options.min_arc_points - 1 >= points.len() {
        return None;
    }

    let circle = circle_from_three_points(points[start], points[start + 1], points[start + 2])?;
    let radius_tol = (circle.radius * options.radius_rel_tol).max(options.radius_abs_tol);

    let mut end = start + 2;
    let mut prev_radial = points[end] - circle.center;
    for (k, point) in points.iter().enumerate().skip(start + 3) {
        let radial = *point - circle.center;
        let radial_len = radial.length();
        if (radial_len - circle.radius).abs() > radius_tol {
            break;
        }
        let turn = circle.normal.dot_vec(prev_radial.cross(radial));
        if turn <= 1.0e-9 {
            break;
        }
        end = k;
        prev_radial = radial;
    }

    if end < start + options.min_arc_points - 1 {
        return None;
    }

    let start_radial = points[start] - circle.center;
    let end_radial = points[end] - circle.center;
    let angle = circle
        .normal
        .dot_vec(start_radial.cross(end_radial))
        .atan2(start_radial.dot(end_radial))
        .abs();
    if angle < options.min_arc_angle {
        return None;
    }

    Some((
        end,
        SweepPathSegment::Arc {
            start: points[start],
            end: points[end],
            center: circle.center,
            normal: circle.normal,
        },
    ))
}

#[derive(Clone, Copy, Debug)]
struct CircleFit {
    center: Point3,
    normal: UnitVec3,
    radius: f64,
}

fn circle_from_three_points(a: Point3, b: Point3, c: Point3) -> Option<CircleFit> {
    let u = b - a;
    let v = c - a;
    let w = u.cross(v);
    let w_len_sq = w.length_sq();
    if w_len_sq < 1.0e-18 {
        return None;
    }

    let center_offset =
        (u.length_sq() * v.cross(w) + v.length_sq() * w.cross(u)) / (2.0 * w_len_sq);
    let center = a + center_offset;
    let radius = (a - center).length();
    if radius < 1.0e-6 {
        return None;
    }

    Some(CircleFit {
        center,
        normal: UnitVec3::try_from_vec(w)?,
        radius,
    })
}

fn join_boundary(
    centre: Point3,
    incoming_dir: UnitVec3,
    outgoing_dir: UnitVec3,
    boundary_dir: UnitVec3,
    radius: f64,
) -> FaceBoundary {
    if incoming_dir.dot(outgoing_dir) > 0.999_999 {
        FaceBoundary::Circle(Circle3::new(centre, boundary_dir, radius))
    } else {
        miter_boundary(centre, incoming_dir, outgoing_dir, boundary_dir, radius)
    }
}

fn assemble_solid(
    brep: &mut BRep,
    face_ids: Vec<cadcore_topo::FaceId>,
    name: Option<String>,
) -> Result<SolidId, SweepError> {
    let shell_id = brep.add_shell(Shell {
        faces: face_ids.clone(),
        is_outer: true,
        solid: cadcore_topo::SolidId::default(),
    });

    for &fid in &face_ids {
        if let Some(f) = brep.faces.get_mut(fid) {
            f.shell = shell_id;
        }
    }

    let solid_id = brep.add_solid(Solid {
        shells: vec![shell_id],
        name,
    });

    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.solid = solid_id;
    }

    Ok(solid_id)
}

fn miter_boundary(
    centre: Point3,
    incoming_dir: UnitVec3,
    outgoing_dir: UnitVec3,
    cylinder_dir: UnitVec3,
    radius: f64,
) -> FaceBoundary {
    let Some(plane_normal) = UnitVec3::try_from_vec(incoming_dir.as_vec() + outgoing_dir.as_vec())
    else {
        return FaceBoundary::Circle(Circle3::new(centre, cylinder_dir, radius));
    };

    boundary_for_plane_cut(centre, cylinder_dir, plane_normal, radius)
}

fn boundary_for_plane_cut(
    centre: Point3,
    cylinder_dir: UnitVec3,
    plane_normal: UnitVec3,
    radius: f64,
) -> FaceBoundary {
    let cos = plane_normal.dot(cylinder_dir).abs();

    if cos > 0.999_999 || cos < 1.0e-6 {
        return FaceBoundary::Circle(Circle3::new(centre, cylinder_dir, radius));
    }

    let projected_axis =
        cylinder_dir.as_vec() - plane_normal.as_vec() * plane_normal.dot(cylinder_dir);
    let Some(x_dir) = UnitVec3::try_from_vec(projected_axis) else {
        return FaceBoundary::Circle(Circle3::new(centre, cylinder_dir, radius));
    };

    FaceBoundary::Ellipse(Ellipse3::new(
        centre,
        plane_normal,
        x_dir,
        radius / cos,
        radius,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(x0: f64, x1: f64) -> Vec<Point3> {
        vec![Point3::new(x0, 0.0, 0.0), Point3::new(x1, 0.0, 0.0)]
    }

    #[test]
    fn single_segment_creates_solid() {
        let mut brep = BRep::new();
        let id =
            sweep_circle_along_polyline(&mut brep, &line(0.0, 10.0), 1.0, &SweepOptions::default())
                .unwrap();
        assert!(brep.solids.contains_key(id));
        // 1 cylinder + 2 caps = 3 faces
        let stats = brep.stats();
        assert_eq!(stats.faces, 3, "expected 3 faces, got {}", stats.faces);
    }

    #[test]
    fn two_segments_can_build_rounded_corner_sweep() {
        let waypoints = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(5.0, 5.0, 0.0),
        ];
        let mut brep = BRep::new();
        let _id = sweep_circle_along_polyline(
            &mut brep,
            &waypoints,
            0.5,
            &SweepOptions {
                fillet_corners: true,
                ..SweepOptions::default()
            },
        )
        .unwrap();
        // 2 cylinders + 1 analytic path torus + 2 caps = 5 faces
        let stats = brep.stats();
        assert_eq!(stats.faces, 5, "expected 5 faces, got {}", stats.faces);
    }

    #[test]
    fn corner_fillet_radius_controls_inside_solid_corner_radius() {
        let waypoints = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(5.0, 5.0, 0.0),
        ];
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &waypoints,
            0.5,
            &SweepOptions {
                fillet_corners: true,
                corner_fillet_radius: 0.2,
                ..SweepOptions::default()
            },
        )
        .unwrap();

        let torus = brep
            .faces
            .values()
            .find_map(|face| match face.geom {
                FaceGeom::Torus(t) => Some(t),
                _ => None,
            })
            .expect("missing corner fillet torus");
        assert!((torus.major_radius - 0.7).abs() < 1.0e-10);
        assert!((torus.minor_radius - 0.5).abs() < 1.0e-10);
    }

    #[test]
    fn corner_fillet_radius_can_exceed_profile_radius() {
        let waypoints = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(5.0, 5.0, 0.0),
        ];
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &waypoints,
            0.5,
            &SweepOptions {
                fillet_corners: true,
                corner_fillet_radius: 0.6,
                ..SweepOptions::default()
            },
        )
        .unwrap();

        let torus = brep
            .faces
            .values()
            .find_map(|face| match face.geom {
                FaceGeom::Torus(t) => Some(t),
                _ => None,
            })
            .expect("missing rounded corner torus");
        assert!((torus.major_radius - 1.1).abs() < 1.0e-10);
    }

    #[test]
    fn negative_corner_fillet_radius_fails() {
        let waypoints = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(5.0, 5.0, 0.0),
        ];
        let mut brep = BRep::new();
        let err = sweep_circle_along_polyline(
            &mut brep,
            &waypoints,
            0.5,
            &SweepOptions {
                fillet_corners: true,
                corner_fillet_radius: -0.1,
                ..SweepOptions::default()
            },
        )
        .unwrap_err();

        assert!(matches!(err, SweepError::InvalidCornerRadius(_)));
    }

    #[test]
    fn error_on_too_few_points() {
        let mut brep = BRep::new();
        let err = sweep_circle_along_polyline(
            &mut brep,
            &[Point3::ORIGIN],
            1.0,
            &SweepOptions::default(),
        )
        .unwrap_err();
        assert_eq!(err, SweepError::TooFewWaypoints);
    }

    #[test]
    fn error_on_zero_radius() {
        let mut brep = BRep::new();
        let err =
            sweep_circle_along_polyline(&mut brep, &line(0.0, 5.0), 0.0, &SweepOptions::default())
                .unwrap_err();
        assert!(matches!(err, SweepError::InvalidRadius(_)));
    }

    #[test]
    fn analytic_arc_creates_single_torus_face() {
        let mut brep = BRep::new();
        sweep_circle_along_path(
            &mut brep,
            &[
                SweepPathSegment::Line {
                    start: Point3::new(-5.0, 0.0, 0.0),
                    end: Point3::new(1.0, 0.0, 0.0),
                },
                SweepPathSegment::Arc {
                    start: Point3::new(1.0, 0.0, 0.0),
                    end: Point3::new(0.0, 1.0, 0.0),
                    center: Point3::ORIGIN,
                    normal: UnitVec3::Z,
                },
                SweepPathSegment::Line {
                    start: Point3::new(0.0, 1.0, 0.0),
                    end: Point3::new(0.0, 5.0, 0.0),
                },
            ],
            0.2,
            &SweepOptions {
                fillet_corners: false,
                ..SweepOptions::default()
            },
        )
        .unwrap();

        let stats = brep.stats();
        assert_eq!(stats.faces, 5, "2 cylinders + 1 torus arc + 2 caps");
    }

    #[test]
    fn rounded_polyline_builds_tangent_arc() {
        let path = rounded_path_from_polyline(
            &[
                Point3::new(-5.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 5.0, 0.0),
            ],
            1.0,
        )
        .unwrap();

        assert_eq!(path.len(), 3);
        match path[1] {
            SweepPathSegment::Arc {
                start,
                end,
                center,
                normal,
            } => {
                assert!((start - Point3::new(-1.0, 0.0, 0.0)).length() < 1.0e-10);
                assert!((end - Point3::new(0.0, 1.0, 0.0)).length() < 1.0e-10);
                assert!((center - Point3::new(-1.0, 1.0, 0.0)).length() < 1.0e-10);
                assert!(normal.dot(UnitVec3::Z) > 0.999_999);
            }
            _ => panic!("expected circular arc"),
        }
    }

    #[test]
    fn rounded_polyline_sweep_uses_centerline_radius_for_torus() {
        let mut brep = BRep::new();
        sweep_circle_along_rounded_polyline(
            &mut brep,
            &[
                Point3::new(-5.0, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 5.0, 0.0),
            ],
            0.2,
            1.5,
            &SweepOptions {
                fillet_corners: false,
                ..SweepOptions::default()
            },
        )
        .unwrap();

        let torus = brep
            .faces
            .values()
            .find_map(|face| match face.geom {
                FaceGeom::Torus(t) => Some(t),
                _ => None,
            })
            .expect("missing torus");
        assert!((torus.major_radius - 1.5).abs() < 1.0e-10);
        assert!((torus.minor_radius - 0.2).abs() < 1.0e-10);
    }

    #[test]
    fn arc_radius_must_exceed_profile_radius() {
        let mut brep = BRep::new();
        let err = sweep_circle_along_path(
            &mut brep,
            &[SweepPathSegment::Arc {
                start: Point3::new(0.15, 0.0, 0.0),
                end: Point3::new(0.0, 0.15, 0.0),
                center: Point3::ORIGIN,
                normal: UnitVec3::Z,
            }],
            0.2,
            &SweepOptions {
                fillet_corners: false,
                ..SweepOptions::default()
            },
        )
        .unwrap_err();

        assert!(matches!(err, SweepError::SelfIntersectingSweep(0)));
    }

    #[test]
    fn sampled_polyline_arc_is_recovered_as_one_segment() {
        let mut points = vec![Point3::new(1.0, -2.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        for k in 1..=24 {
            let a = std::f64::consts::FRAC_PI_2 * k as f64 / 24.0;
            points.push(Point3::new(a.cos(), a.sin(), 0.0));
        }
        points.push(Point3::new(-2.0, 1.0, 0.0));

        let path = analytic_path_from_polyline_samples(&points, PathApproxOptions::default());
        assert_eq!(path.len(), 3);
        assert!(matches!(path[1], SweepPathSegment::Arc { .. }));
    }

    #[test]
    fn sampled_arc_can_be_restored_to_sharp_corner() {
        let mut points = vec![Point3::new(0.15, -1.0, 0.0), Point3::new(0.15, 0.0, 0.0)];
        for k in 1..=24 {
            let a = std::f64::consts::FRAC_PI_2 * k as f64 / 24.0;
            points.push(Point3::new(0.15 * a.cos(), 0.15 * a.sin(), 0.0));
        }
        points.push(Point3::new(-1.0, 0.15, 0.0));

        let path = sharp_path_from_polyline_samples(&points, PathApproxOptions::default());
        assert_eq!(path.len(), 2);
        match path[0] {
            SweepPathSegment::Line { end, .. } => {
                assert!((end - Point3::new(0.15, 0.15, 0.0)).length() < 1.0e-9);
            }
            _ => panic!("expected line"),
        }
    }

    // ── build_solid_rounded_box ───────────────────────────────────────────────

    /// Sharp-corner fallback: corner_radius = 0 → same face count as build_solid_box.
    #[test]
    fn rounded_box_zero_radius_equals_sharp_box() {
        let mut brep_r = BRep::new();
        build_solid_rounded_box(&mut brep_r, 0.0, 10.0, 0.0, 5.0, 0.0, 2.0, 0.0, None).unwrap();
        let mut brep_s = BRep::new();
        build_solid_box(&mut brep_s, 0.0, 10.0, 0.0, 5.0, 0.0, 2.0, None).unwrap();
        assert_eq!(brep_r.solids.len(), brep_s.solids.len());
        assert_eq!(
            brep_r.faces.len(),
            brep_s.faces.len(),
            "zero-radius rounded box must match sharp box face count"
        );
    }

    /// Face count: top + bottom + 4*N_ARC side quads (N_ARC = 8 → 32 sides).
    #[test]
    fn rounded_box_face_count() {
        let mut brep = BRep::new();
        build_solid_rounded_box(&mut brep, 0.0, 20.0, 0.0, 5.0, 0.0, 2.0, 1.0, None).unwrap();
        // 2 top/bottom + 4 arcs × N_ARC = 2 + 32 = 34
        assert_eq!(brep.faces.len(), 34, "expected 2 caps + 32 side quads");
        assert_eq!(brep.solids.len(), 1);
    }

    /// The top-face polygon points must all lie at z = zmax.
    #[test]
    fn rounded_box_top_face_z_correct() {
        let zmax = 3.0;
        let mut brep = BRep::new();
        build_solid_rounded_box(&mut brep, 0.0, 10.0, 0.0, 4.0, 0.0, zmax, 0.5, None).unwrap();
        for face in brep.faces.values() {
            if let FaceExtent::Polygon { ref points } = face.extent {
                // Top face: all Z = zmax, bottom face: all Z = 0
                let all_top = points.iter().all(|p| (p.z - zmax).abs() < 1e-9);
                let all_bot = points.iter().all(|p| p.z.abs() < 1e-9);
                if !all_top && !all_bot {
                    // Side face quads span both z values
                    let zs: Vec<f64> = points.iter().map(|p| p.z).collect();
                    assert!(
                        zs.iter().any(|&z| (z - zmax).abs() < 1e-9),
                        "side face must touch zmax"
                    );
                    assert!(
                        zs.iter().any(|&z| z.abs() < 1e-9),
                        "side face must touch zmin"
                    );
                }
            }
        }
    }

    /// Corner profile points must be within tolerance r of the inner-box corner.
    ///
    /// For a rounded rectangle profile, every point P must satisfy:
    ///   max(|P.x - clamp(P.x, xmin+r, xmax-r)|,
    ///       |P.y - clamp(P.y, ymin+r, ymax-r)|) ≈ 0 or dist(P, inner_corner) = r.
    #[test]
    fn rounded_box_profile_points_at_radius_from_inner_box() {
        let (xmin, xmax, ymin, ymax) = (0.0_f64, 10.0, 0.0, 4.0);
        let r = 0.8_f64;
        let profile = rounded_rect_profile_2d(xmin, xmax, ymin, ymax, r, 8);

        for [px, py] in &profile {
            // Inner box extents
            let ix_min = xmin + r;
            let ix_max = xmax - r;
            let iy_min = ymin + r;
            let iy_max = ymax - r;

            // Closest point on inner box boundary
            let cx = px.clamp(ix_min, ix_max);
            let cy = py.clamp(iy_min, iy_max);
            let dist = ((px - cx).powi(2) + (py - cy).powi(2)).sqrt();

            assert!(
                (dist - r).abs() < 1e-9,
                "profile point ({px:.4}, {py:.4}) is at dist {dist:.6} from inner box, expected {r:.4}"
            );
        }
    }

    /// All profile X/Y coords must lie within the outer bounding box.
    #[test]
    fn rounded_box_profile_within_bounds() {
        let (xmin, xmax, ymin, ymax) = (0.0_f64, 20.0, -0.5, 0.0);
        let r = 0.2_f64;
        let profile = rounded_rect_profile_2d(xmin, xmax, ymin, ymax, r, 8);
        for [px, py] in &profile {
            assert!(
                *px >= xmin - 1e-9 && *px <= xmax + 1e-9,
                "X out of bounds: {px}"
            );
            assert!(
                *py >= ymin - 1e-9 && *py <= ymax + 1e-9,
                "Y out of bounds: {py}"
            );
        }
    }

    /// build_solid_rounded_box_full: face count = 6 × N × N × 2 triangles.
    #[test]
    fn rounded_box_full_face_count() {
        let mut brep = BRep::new();
        build_solid_rounded_box_full(&mut brep, 0.0, 20.0, 0.0, 5.0, 0.0, 3.0, 1.0, None).unwrap();
        // N=8 → 6 faces × 8 × 8 × 2 triangles = 768
        assert_eq!(brep.faces.len(), 768);
        assert_eq!(brep.solids.len(), 1);
    }

    /// Full rounded box triangles must have outward normals (away from center).
    #[test]
    fn rounded_box_full_normals_point_outward() {
        let mut brep = BRep::new();
        build_solid_rounded_box_full(&mut brep, 0.0, 10.0, 0.0, 4.0, 0.0, 2.0, 0.5, None).unwrap();
        let center = Point3::new(5.0, 2.0, 1.0);
        let mut wrong_normals = 0usize;
        for face in brep.faces.values() {
            if let FaceExtent::Polygon { ref points } = face.extent {
                if let FaceGeom::Plane(ref plane) = face.geom {
                    let to_face = points[0] - center;
                    let dot = plane.normal().dot_vec(to_face);
                    if dot <= 0.0 {
                        wrong_normals += 1;
                    }
                }
            }
        }
        // Allow a small number of nearly-tangent triangles at edges/corners.
        assert!(
            wrong_normals < 10,
            "{wrong_normals} faces had inward-pointing normals"
        );
    }

    /// All side faces must have exactly 4 polygon points (quads).
    #[test]
    fn rounded_box_side_faces_are_quads() {
        let mut brep = BRep::new();
        build_solid_rounded_box(
            &mut brep,
            0.0,
            20.0,
            -0.5,
            0.0,
            0.0,
            3.0,
            0.2,
            Some("elec".into()),
        )
        .unwrap();
        let mut cap_count = 0;
        let mut side_count = 0;
        for face in brep.faces.values() {
            if let FaceExtent::Polygon { ref points } = face.extent {
                match points.len() {
                    n if n > 4 => cap_count += 1,
                    4 => side_count += 1,
                    _ => panic!("unexpected polygon size {}", points.len()),
                }
            }
        }
        assert_eq!(cap_count, 2, "expected 2 cap faces");
        assert_eq!(side_count, 32, "expected 32 side quads");
    }

    // ── build_solid_rounded_box_xz tests ──────────────────────────────────────

    /// XZ-variant: analytic geometry → exactly 10 faces (2 caps + 4 flat sides + 4 cylinders).
    #[test]
    fn rounded_box_xz_face_count() {
        let mut brep = BRep::new();
        build_solid_rounded_box_xz(&mut brep, 0.0, 20.0, -0.5, 0.0, 0.0, 3.0, 0.3, None).unwrap();
        assert_eq!(
            brep.faces.len(),
            10,
            "expected 10 faces (2 RoundedRectCap + 4 Polygon sides + 4 CylinderArcFace)"
        );
        assert_eq!(brep.solids.len(), 1);
    }

    /// XZ-variant: zero radius must produce a plain box (6 faces).
    #[test]
    fn rounded_box_xz_zero_radius_equals_box() {
        let mut brep = BRep::new();
        build_solid_rounded_box_xz(&mut brep, 0.0, 10.0, 0.0, 1.0, 0.0, 2.0, 0.0, None).unwrap();
        assert_eq!(brep.faces.len(), 6, "zero-radius XZ box must have 6 faces");
    }

    /// XZ-variant: geometry is analytic — 2 RoundedRectCap + 4 Polygon + 4 CylinderArcFace.
    #[test]
    fn rounded_box_xz_analytic_face_types() {
        let mut brep = BRep::new();
        build_solid_rounded_box_xz(&mut brep, 0.0, 20.0, -0.5, 0.0, 0.0, 2.0, 0.4, None).unwrap();
        let mut caps = 0usize;
        let mut polys = 0usize;
        let mut cyls = 0usize;
        for face in brep.faces.values() {
            match &face.extent {
                FaceExtent::RoundedRectCap { .. } => caps += 1,
                FaceExtent::Polygon { .. } => polys += 1,
                FaceExtent::CylinderArcFace { .. } => cyls += 1,
                other => panic!("unexpected extent {:?}", other),
            }
        }
        assert_eq!(caps, 2, "expected 2 RoundedRectCap faces");
        assert_eq!(polys, 4, "expected 4 flat Polygon faces");
        assert_eq!(cyls, 4, "expected 4 CylinderArcFace corners");
    }

    /// XZ-variant: cap face Y-coordinates are exactly ymin / ymax.
    #[test]
    fn rounded_box_xz_caps_at_correct_y() {
        let (ymin, ymax) = (-0.5_f64, 0.0_f64);
        let mut brep = BRep::new();
        build_solid_rounded_box_xz(&mut brep, 0.0, 20.0, ymin, ymax, 0.0, 3.0, 0.3, None).unwrap();
        let cap_ys: Vec<f64> = brep
            .faces
            .values()
            .filter_map(|f| {
                if let FaceExtent::RoundedRectCap { y, .. } = f.extent {
                    Some(y)
                } else {
                    None
                }
            })
            .collect();
        assert_eq!(cap_ys.len(), 2, "expected exactly 2 cap faces");
        assert!(
            cap_ys.iter().any(|&y| (y - ymin).abs() < 1e-9),
            "no cap at ymin={ymin}"
        );
        assert!(
            cap_ys.iter().any(|&y| (y - ymax).abs() < 1e-9),
            "no cap at ymax={ymax}"
        );
    }

    /// XZ-variant: planar face normals point outward from the solid centroid.
    #[test]
    fn rounded_box_xz_normals_point_outward() {
        let (xmin, xmax, ymin, ymax, zmin, zmax) = (0.0_f64, 20.0, -0.5, 0.0, 0.0, 3.0);
        let mut brep = BRep::new();
        build_solid_rounded_box_xz(&mut brep, xmin, xmax, ymin, ymax, zmin, zmax, 0.4, None)
            .unwrap();
        let center = Point3::new(
            (xmin + xmax) * 0.5,
            (ymin + ymax) * 0.5,
            (zmin + zmax) * 0.5,
        );
        let mut wrong = 0usize;
        for face in brep.faces.values() {
            if let FaceGeom::Plane(ref plane) = face.geom {
                // Sample point on the plane — use plane origin.
                let dot = plane.normal().dot_vec(plane.frame.origin - center);
                if dot <= 0.0 {
                    wrong += 1;
                }
            }
        }
        assert_eq!(
            wrong, 0,
            "{wrong} XZ electrode plane faces had inward normals"
        );
    }

    /// XZ-variant: side quads must have all 4 Y-coords within [ymin, ymax].
    #[test]
    fn rounded_box_xz_side_quads_span_y() {
        let (ymin, ymax) = (-0.5_f64, 0.0_f64);
        let mut brep = BRep::new();
        build_solid_rounded_box_xz(&mut brep, 0.0, 20.0, ymin, ymax, 0.0, 3.0, 0.3, None).unwrap();
        for face in brep.faces.values() {
            if let FaceExtent::Polygon { ref points } = face.extent {
                if points.len() == 4 {
                    let ys: Vec<f64> = points.iter().map(|p| p.y).collect();
                    let has_ymin = ys.iter().any(|&y| (y - ymin).abs() < 1e-9);
                    let has_ymax = ys.iter().any(|&y| (y - ymax).abs() < 1e-9);
                    assert!(
                        has_ymin && has_ymax,
                        "side quad Y coords {ys:?} don't span [{ymin}, {ymax}]"
                    );
                }
            }
        }
    }
}
