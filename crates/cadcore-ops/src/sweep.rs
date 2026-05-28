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
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{
    BRep, Face, FaceBoundary, FaceExtent, FaceGeom, FaceNormal, Shell, Solid, SolidId,
};

/// Options that control the sweep.
#[derive(Clone, Debug)]
pub struct SweepOptions {
    /// Use toroidal fillets at corners (smooth G1 join).
    /// If `false`, use miter-plane cuts (sharp join — faster STEP output).
    pub fillet_corners: bool,
    /// Optional name for the resulting solid.
    pub name: Option<String>,
}

impl Default for SweepOptions {
    fn default() -> Self {
        Self {
            fillet_corners: false,
            name: None,
        }
    }
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
    if waypoints.len() < 2 {
        return Err(SweepError::TooFewWaypoints);
    }
    if radius <= 0.0 {
        return Err(SweepError::InvalidRadius(radius));
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

    // Pre-allocate face ids for shell construction.
    let mut face_ids = Vec::new();

    // ── Start cap ────────────────────────────────────────────────────────────
    let _cap_start_id = build_end_cap(brep, waypoints[0], -dirs[0], radius, &mut face_ids);

    // ── Cylinder + connector loop ─────────────────────────────────────────────
    // For each segment i we build: cylinder face for segment i, then (if not last)
    // a connector at waypoints[i+1].
    for i in 0..n - 1 {
        let p0 = waypoints[i];
        let p1 = waypoints[i + 1];
        let dir = dirs[i];
        let length = (p1 - p0).length();
        let start_bound = if i == 0 {
            FaceBoundary::Circle(Circle3::new(p0, dir, radius))
        } else {
            miter_boundary(p0, dirs[i - 1], dir, dir, radius)
        };
        let end_bound = if i + 1 == n - 1 {
            FaceBoundary::Circle(Circle3::new(p1, dir, radius))
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
            shell: cadcore_topo::ShellId::default(), // patched later
            extent: FaceExtent::Cylinder {
                length,
                start: start_bound,
                end: end_bound,
            },
        });
        face_ids.push(cyl_face_id);

        // --- Connector at interior vertex ---
        if i + 1 < n - 1 {
            let next_dir = dirs[i + 1];
            let vertex = waypoints[i + 1];
            let cos_angle = dir.dot(next_dir);

            if opts.fillet_corners && (1.0 - cos_angle.abs()) > 1e-4 {
                // Build toroidal fillet
                if let Some(torus) = TorusSurf::for_corner(vertex, dir, next_dir, radius) {
                    // Compute the two minor-circle boundaries for the STEP face bounds.
                    //
                    // The torus axis A = dir × next_dir (bend-plane normal).
                    // Both `dir` and `next_dir` lie in the equatorial plane of the torus.
                    //
                    // The spine point at the incoming junction has radial direction:
                    //   r_hat_start = normalize(dir × A)
                    // The spine point at the outgoing junction has radial direction:
                    //   r_hat_end   = normalize(next_dir × A)
                    //
                    // Derivation: axis × r_hat = dir ⟹ r_hat = normalize(dir × axis)
                    // (valid when dir ⊥ axis, which holds since axis = dir × next_dir).
                    let extent = if let Some(taxis) = UnitVec3::try_from_vec(dir.cross(next_dir)) {
                        let rs_vec = dir.cross(taxis);
                        let re_vec = next_dir.cross(taxis);
                        match (
                            UnitVec3::try_from_vec(rs_vec),
                            UnitVec3::try_from_vec(re_vec),
                        ) {
                            (Some(r_s), Some(r_e)) => {
                                let start_pt = torus.frame.origin + r_s * torus.major_radius;
                                let end_pt = torus.frame.origin + r_e * torus.major_radius;
                                FaceExtent::TorusFillet {
                                    start_circle: Circle3::new(start_pt, dir, torus.minor_radius),
                                    end_circle: Circle3::new(end_pt, next_dir, torus.minor_radius),
                                }
                            }
                            _ => FaceExtent::None,
                        }
                    } else {
                        FaceExtent::None
                    };

                    let torus_loop_id = build_placeholder_loop(brep);
                    let torus_id = brep.add_face(Face {
                        geom: FaceGeom::Torus(torus),
                        normal: FaceNormal::Same,
                        outer_loop: torus_loop_id,
                        inner_loops: vec![],
                        shell: cadcore_topo::ShellId::default(),
                        extent,
                    });
                    face_ids.push(torus_id);
                }
            }
            // If fillet_corners=false or nearly straight: miter plane cut.
            // In a miter join the cylinder segments share the elliptic edge.
            // The STEP writer handles the shared edge — no extra face needed.
        }
    }

    // ── End cap ───────────────────────────────────────────────────────────────
    let last = n - 1;
    let _cap_end_id = build_end_cap(brep, waypoints[last], dirs[last - 1], radius, &mut face_ids);

    // ── Assemble shell + solid ────────────────────────────────────────────────
    let shell_id = brep.add_shell(Shell {
        faces: face_ids.clone(),
        is_outer: true,
        solid: cadcore_topo::SolidId::default(), // patched below
    });

    // Patch shell back-ref into faces
    for &fid in &face_ids {
        if let Some(f) = brep.faces.get_mut(fid) {
            f.shell = shell_id;
        }
    }

    let solid_id = brep.add_solid(Solid {
        shells: vec![shell_id],
        name: opts.name.clone(),
    });

    // Patch solid back-ref into shell
    if let Some(sh) = brep.shells.get_mut(shell_id) {
        sh.solid = solid_id;
    }

    Ok(solid_id)
}

/// Sweep a circle of radius `radius` along an analytic path.
///
/// Straight segments become cylindrical faces. Circular arc segments become
/// toroidal faces following the exact centre-line radius. This is the preferred
/// API for rounded scaffold corners because it avoids approximating an arc with
/// many short straight cylinders.
///
/// # Errors
///
/// Returns `Err` if the path is empty, contains invalid geometry, or adjacent
/// segments are not position-continuous.
pub fn sweep_circle_along_path(
    brep: &mut BRep,
    segments: &[SweepPathSegment],
    radius: f64,
    opts: &SweepOptions,
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
    let _cap_start_id = build_end_cap(
        brep,
        infos[0].start,
        -infos[0].start_tangent,
        radius,
        &mut face_ids,
    );

    for (idx, info) in infos.iter().enumerate() {
        let start_bound = if idx == 0 {
            FaceBoundary::Circle(Circle3::new(info.start, info.start_tangent, radius))
        } else {
            join_boundary(
                info.start,
                infos[idx - 1].end_tangent,
                info.start_tangent,
                info.start_tangent,
                radius,
            )
        };
        let end_bound = if idx + 1 == infos.len() {
            FaceBoundary::Circle(Circle3::new(info.end, info.end_tangent, radius))
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

        if idx + 1 < infos.len() {
            let current = info.end_tangent;
            let next = infos[idx + 1].start_tangent;
            if opts.fillet_corners && (1.0 - current.dot(next).abs()) > 1.0e-4 {
                if let Some(connector_id) =
                    build_corner_connector(brep, info.end, current, next, radius)
                {
                    face_ids.push(connector_id);
                }
            }
        }
    }

    let last = infos.last().expect("non-empty path");
    let _cap_end_id = build_end_cap(brep, last.end, last.end_tangent, radius, &mut face_ids);

    assemble_solid(brep, face_ids, opts.name.clone())
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

// ── Internal helpers ─────────────────────────────────────────────────────────

/// Build an end-cap (planar disk) face at `centre` with normal `normal`.
///
/// The outward normal of the start cap points *away from* the solid (−dir),
/// and the end cap outward normal is +dir.
fn build_end_cap(
    brep: &mut BRep,
    centre: Point3,
    normal: UnitVec3,
    radius: f64,
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
        extent: FaceExtent::Disk { radius },
    });
    face_ids.push(fid);
    fid
}

/// Insert a placeholder empty loop (patched with real co-edges later by the
/// STEP writer or a topology validator).
fn build_placeholder_loop(brep: &mut BRep) -> cadcore_topo::LoopId {
    // We insert a Loop with a sentinel co-edge id (the default key).
    // The STEP writer's boundary-stitching pass fills in real co-edges.
    brep.add_loop(cadcore_topo::Loop {
        start: cadcore_topo::CoEdgeId::default(),
        face: cadcore_topo::FaceId::default(),
    })
}

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

fn build_corner_connector(
    brep: &mut BRep,
    vertex: Point3,
    dir: UnitVec3,
    next_dir: UnitVec3,
    radius: f64,
) -> Option<cadcore_topo::FaceId> {
    let torus = TorusSurf::for_corner(vertex, dir, next_dir, radius)?;
    let extent = if let Some(taxis) = UnitVec3::try_from_vec(dir.cross(next_dir)) {
        let rs_vec = dir.cross(taxis);
        let re_vec = next_dir.cross(taxis);
        match (
            UnitVec3::try_from_vec(rs_vec),
            UnitVec3::try_from_vec(re_vec),
        ) {
            (Some(r_s), Some(r_e)) => {
                let start_pt = torus.frame.origin + r_s * torus.major_radius;
                let end_pt = torus.frame.origin + r_e * torus.major_radius;
                FaceExtent::TorusFillet {
                    start_circle: Circle3::new(start_pt, dir, torus.minor_radius),
                    end_circle: Circle3::new(end_pt, next_dir, torus.minor_radius),
                }
            }
            _ => FaceExtent::None,
        }
    } else {
        FaceExtent::None
    };

    let torus_loop_id = build_placeholder_loop(brep);
    Some(brep.add_face(Face {
        geom: FaceGeom::Torus(torus),
        normal: FaceNormal::Same,
        outer_loop: torus_loop_id,
        inner_loops: vec![],
        shell: cadcore_topo::ShellId::default(),
        extent,
    }))
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
    fn two_segments_can_add_legacy_connector() {
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
                name: None,
            },
        )
        .unwrap();
        // 2 cylinders + 1 torus fillet + 2 caps = 5 faces
        let stats = brep.stats();
        assert_eq!(stats.faces, 5, "expected 5 faces, got {}", stats.faces);
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
                name: None,
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
                name: None,
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
                name: None,
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
}
