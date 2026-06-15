//! Low-level STEP entity emission helpers.
#![allow(dead_code)]
//!
//! The STEP writer maintains a monotonically increasing entity counter.
//! Every entity is emitted as `#N = ENTITY_TYPE(args);`.

use std::collections::HashMap;
use std::fmt::Write as FmtWrite;

use cadcore_geom::{Circle3, CylSurf, Ellipse3, Line3, Plane3, SphereSurf, TorusSurf};
use cadcore_math::{Frame3, Point3, UnitVec3, Vec3};

/// Errors that can occur during STEP serialisation.
#[derive(Debug, Clone)]
pub enum StepError {
    /// Fmt error (out of memory / write failure).
    Fmt(std::fmt::Error),
}

impl From<std::fmt::Error> for StepError {
    fn from(e: std::fmt::Error) -> Self {
        Self::Fmt(e)
    }
}

impl std::fmt::Display for StepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fmt(e) => write!(f, "fmt error: {e}"),
        }
    }
}
impl std::error::Error for StepError {}

#[derive(Hash, PartialEq, Eq, Debug, Clone, Copy)]
pub(crate) enum StepCurveKey {
    Line {
        v1: usize,
        v2: usize,
    },
    Circle {
        center: [i64; 3],
        radius_micro: i64,
        normal: [i64; 3],
    },
    Ellipse {
        center: [i64; 3],
        semi_major_micro: i64,
        semi_minor_micro: i64,
        normal: [i64; 3],
    },
    /// Canonical arc edge key: v_start ≤ v_end lexicographically so that
    /// ArcEdge(A→B) and ArcEdge(B→A) map to the same entry.
    ArcEdge {
        v_start: [i64; 3],
        v_end: [i64; 3],
        center: [i64; 3],
        radius_micro: i64,
        axis: [i64; 3],
        x_ref: [i64; 3],
    },
    /// Polyline (degree-1 spline) intersection edge, keyed by its two endpoint
    /// vertex ids (sorted) plus an order-independent centroid discriminator so
    /// two different arcs sharing the same endpoints (the two halves of one
    /// intersection loop) do not collide.  Two faces meeting at a sampled cyl∩cyl
    /// curve share one such EDGE_CURVE.
    Polyline {
        v1: usize,
        v2: usize,
        mid: [i64; 3],
    },
    /// Elliptic arc edge, keyed direction-agnostically by endpoints + ellipse.
    EllipseArc {
        v_start: [i64; 3],
        v_end: [i64; 3],
        center: [i64; 3],
        semi_major_micro: i64,
        semi_minor_micro: i64,
        normal: [i64; 3],
    },
}

pub(crate) fn point_key(p: Point3) -> [i64; 3] {
    [
        (p.x * 100_000.0).round() as i64,
        (p.y * 100_000.0).round() as i64,
        (p.z * 100_000.0).round() as i64,
    ]
}

/// Quantise a unit direction vector to 4 decimal places (±10 000 units).
/// No sign normalisation — the raw direction is significant for arc geometry.
pub(crate) fn dir_key(v: UnitVec3) -> [i64; 3] {
    let vec = v.as_vec();
    [
        (vec.x * 10_000.0).round() as i64,
        (vec.y * 10_000.0).round() as i64,
        (vec.z * 10_000.0).round() as i64,
    ]
}

/// Build a direction-agnostic `ArcEdge` key: whichever of `v_start_pos` /
/// `v_end_pos` is lexicographically smaller becomes the canonical `v_start`.
pub(crate) fn arc_edge_key(
    v_start_pos: [i64; 3],
    v_end_pos: [i64; 3],
    center: [i64; 3],
    radius_micro: i64,
    axis: [i64; 3],
    x_ref: [i64; 3],
) -> StepCurveKey {
    let (v_s, v_e) = if v_start_pos <= v_end_pos {
        (v_start_pos, v_end_pos)
    } else {
        (v_end_pos, v_start_pos)
    };
    
    let mut canonical_axis = axis;
    let mut canonical_xref = x_ref;
    
    let ax = axis[0];
    let ay = axis[1];
    let az = axis[2];
    
    let should_flip = if ax != 0 {
        ax < 0
    } else if ay != 0 {
        ay < 0
    } else {
        az < 0
    };
    
    if should_flip {
        canonical_axis = [-ax, -ay, -az];
        canonical_xref = [-x_ref[0], -x_ref[1], -x_ref[2]];
    }

    StepCurveKey::ArcEdge {
        v_start: v_s,
        v_end: v_e,
        center,
        radius_micro,
        axis: canonical_axis,
        x_ref: canonical_xref,
    }
}

/// Mutable counter + output buffer, shared across all entity emission calls.
pub(crate) struct Ctx {
    pub counter: usize,
    pub out: String,
    pub point_cache: HashMap<[i64; 3], usize>,
    pub vertex_cache: HashMap<[i64; 3], usize>,
    pub edge_cache: HashMap<StepCurveKey, (usize, usize)>, // (ec_id, orig_start_vtx_id)
    /// Lazily-emitted shared 2D parametric representation context (for pcurves).
    pub param_ctx_2d: Option<usize>,
}

impl Ctx {
    pub fn new() -> Self {
        Self {
            counter: 0,
            out: String::with_capacity(8192),
            point_cache: HashMap::new(),
            vertex_cache: HashMap::new(),
            edge_cache: HashMap::new(),
            param_ctx_2d: None,
        }
    }

    /// Allocate the next entity id (1-based) and return it.
    pub fn next_id(&mut self) -> usize {
        self.counter += 1;
        self.counter
    }

    /// Write a STEP line: `#N = TYPE(...);\n`.
    pub fn emit_raw(&mut self, id: usize, line: &str) -> Result<(), StepError> {
        writeln!(self.out, "#{id} = {line};")?;
        Ok(())
    }
}

// ── Point3 → CARTESIAN_POINT ─────────────────────────────────────────────────

pub(crate) fn emit_point(ctx: &mut Ctx, p: Point3, label: &str) -> Result<usize, StepError> {
    let key = point_key(p);
    if let Some(&id) = ctx.point_cache.get(&key) {
        return Ok(id);
    }
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        &format!(
            "CARTESIAN_POINT('{label}',({:.10},{:.10},{:.10}))",
            p.x, p.y, p.z
        ),
    )?;
    ctx.point_cache.insert(key, id);
    Ok(id)
}

pub(crate) fn emit_vertex_point(ctx: &mut Ctx, p: Point3) -> Result<usize, StepError> {
    let key = point_key(p);
    if let Some(&id) = ctx.vertex_cache.get(&key) {
        return Ok(id);
    }
    let cp_id = emit_point(ctx, p, "v")?;
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!("VERTEX_POINT('',#{cp_id})"))?;
    ctx.vertex_cache.insert(key, id);
    Ok(id)
}

// ── Vec3 → DIRECTION ─────────────────────────────────────────────────────────

pub(crate) fn emit_direction(ctx: &mut Ctx, v: Vec3, label: &str) -> Result<usize, StepError> {
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        &format!("DIRECTION('{label}',({:.10},{:.10},{:.10}))", v.x, v.y, v.z),
    )?;
    Ok(id)
}

pub(crate) fn emit_unit_direction(
    ctx: &mut Ctx,
    u: UnitVec3,
    label: &str,
) -> Result<usize, StepError> {
    emit_direction(ctx, u.as_vec(), label)
}

// ── VECTOR (direction + magnitude) ──────────────────────────────────────────

pub(crate) fn emit_vector(
    ctx: &mut Ctx,
    dir_id: usize,
    magnitude: f64,
    label: &str,
) -> Result<usize, StepError> {
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!("VECTOR('{label}',#{dir_id},{magnitude:.10})"))?;
    Ok(id)
}

// ── Frame3 → AXIS2_PLACEMENT_3D ──────────────────────────────────────────────

pub(crate) fn emit_axis2(ctx: &mut Ctx, frame: Frame3, label: &str) -> Result<usize, StepError> {
    let origin_id = emit_point(ctx, frame.origin, label)?;
    let z_id = emit_unit_direction(ctx, frame.z, &format!("{label}_Z"))?;
    let x_id = emit_unit_direction(ctx, frame.x, &format!("{label}_X"))?;
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        &format!("AXIS2_PLACEMENT_3D('{label}',#{origin_id},#{z_id},#{x_id})"),
    )?;
    Ok(id)
}

// ── Plane3 → PLANE ──────────────────────────────────────────────────────────

pub(crate) fn emit_plane(ctx: &mut Ctx, plane: &Plane3) -> Result<usize, StepError> {
    let ax_id = emit_axis2(ctx, plane.frame, "plane")?;
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!("PLANE('',#{ax_id})"))?;
    Ok(id)
}

// ── CylSurf → CYLINDRICAL_SURFACE ────────────────────────────────────────────

pub(crate) fn emit_cylinder(ctx: &mut Ctx, cyl: &CylSurf) -> Result<usize, StepError> {
    let ax_id = emit_axis2(ctx, cyl.frame, "cyl_axis")?;
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        &format!("CYLINDRICAL_SURFACE('',#{ax_id},{:.10})", cyl.radius),
    )?;
    Ok(id)
}

// ── SphereSurf → SPHERICAL_SURFACE ────────────────────────────────────────────

pub(crate) fn emit_sphere(ctx: &mut Ctx, sph: &SphereSurf) -> Result<usize, StepError> {
    let ax_id = emit_axis2(ctx, sph.frame, "sph_axis")?;
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        &format!("SPHERICAL_SURFACE('',#{ax_id},{:.10})", sph.radius),
    )?;
    Ok(id)
}

// ── TorusSurf → TOROIDAL_SURFACE ──────────────────────────────────────────────

pub(crate) fn emit_torus(ctx: &mut Ctx, torus: &TorusSurf) -> Result<usize, StepError> {
    let ax_id = emit_axis2(ctx, torus.frame, "torus_axis")?;
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        &format!(
            "TOROIDAL_SURFACE('',#{ax_id},{:.10},{:.10})",
            torus.major_radius, torus.minor_radius
        ),
    )?;
    Ok(id)
}

// ── RationalBSplineSurface → ( B_SPLINE_SURFACE_WITH_KNOTS … RATIONAL … ) ─────

/// Emit an exact rational-NURBS surface patch as the AP214 *complex* entity
/// `( BOUNDED_SURFACE() B_SPLINE_SURFACE(…) B_SPLINE_SURFACE_WITH_KNOTS(…)
///    GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(…)
///    REPRESENTATION_ITEM('') SURFACE() )`.
///
/// Control points are funnelled through the point cache so they may be shared.
pub(crate) fn emit_rational_bspline_surface(
    ctx: &mut Ctx,
    s: &crate::nurbs::RationalBSplineSurface,
) -> Result<usize, StepError> {
    // 1. Control points (one CARTESIAN_POINT per grid node), [u][v].
    let mut pt_ids: Vec<Vec<usize>> = Vec::with_capacity(s.ctrl.len());
    for row in &s.ctrl {
        let mut ids = Vec::with_capacity(row.len());
        for &p in row {
            ids.push(emit_point(ctx, p, "bsp_cp")?);
        }
        pt_ids.push(ids);
    }

    let fmt_row_pts = |ids: &[usize]| -> String {
        let inner: Vec<String> = ids.iter().map(|id| format!("#{id}")).collect();
        format!("({})", inner.join(","))
    };
    let cp_list: String = pt_ids
        .iter()
        .map(|row| fmt_row_pts(row))
        .collect::<Vec<_>>()
        .join(",");

    let fmt_ints = |v: &[usize]| -> String {
        v.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    };
    let fmt_knots = |v: &[f64]| -> String {
        v.iter()
            .map(|x| format!("{x:.10}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    let fmt_weights = |w: &[Vec<f64>]| -> String {
        w.iter()
            .map(|row| {
                let inner: Vec<String> = row.iter().map(|x| format!("{x:.10}")).collect();
                format!("({})", inner.join(","))
            })
            .collect::<Vec<_>>()
            .join(",")
    };

    let u_closed = if s.u_closed { ".T." } else { ".F." };
    let v_closed = if s.v_closed { ".T." } else { ".F." };

    let id = ctx.next_id();
    let line = format!(
        "( BOUNDED_SURFACE() \
B_SPLINE_SURFACE({u_deg},{v_deg},({cp_list}),.UNSPECIFIED.,{u_closed},{v_closed},.F.) \
B_SPLINE_SURFACE_WITH_KNOTS(({u_mult}),({v_mult}),({u_knots}),({v_knots}),.UNSPECIFIED.) \
GEOMETRIC_REPRESENTATION_ITEM() \
RATIONAL_B_SPLINE_SURFACE(({weights})) \
REPRESENTATION_ITEM('') \
SURFACE() )",
        u_deg = s.u_degree,
        v_deg = s.v_degree,
        u_mult = fmt_ints(&s.u_mult),
        v_mult = fmt_ints(&s.v_mult),
        u_knots = fmt_knots(&s.u_knots),
        v_knots = fmt_knots(&s.v_knots),
        weights = fmt_weights(&s.weights),
    );
    ctx.emit_raw(id, &line)?;
    Ok(id)
}

// ── pcurves: 2D parametric curves + PCURVE + SURFACE_CURVE / SEAM_CURVE ───────

/// Shared 2D parametric representation context (emitted once per file).
fn param_context_2d(ctx: &mut Ctx) -> Result<usize, StepError> {
    if let Some(id) = ctx.param_ctx_2d {
        return Ok(id);
    }
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        "( GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('2D','') )",
    )?;
    ctx.param_ctx_2d = Some(id);
    Ok(id)
}

/// `CARTESIAN_POINT` with two coordinates (a `(u,v)` parameter point).
pub(crate) fn emit_point2d(ctx: &mut Ctx, u: f64, v: f64) -> Result<usize, StepError> {
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!("CARTESIAN_POINT('',({u:.10},{v:.10}))"))?;
    Ok(id)
}

/// A 2D `LINE` through `(u0,v0)` with direction `(du,dv)` (auto-normalised; the
/// VECTOR magnitude carries the length so the parametrisation spans the curve).
pub(crate) fn emit_line2d(
    ctx: &mut Ctx,
    u0: f64,
    v0: f64,
    du: f64,
    dv: f64,
) -> Result<usize, StepError> {
    let len = (du * du + dv * dv).sqrt().max(1e-12);
    let pt = emit_point2d(ctx, u0, v0)?;
    let dir = ctx.next_id();
    ctx.emit_raw(
        dir,
        &format!("DIRECTION('',({:.10},{:.10}))", du / len, dv / len),
    )?;
    let vec = ctx.next_id();
    ctx.emit_raw(vec, &format!("VECTOR('',#{dir},{len:.10})"))?;
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!("LINE('',#{pt},#{vec})"))?;
    Ok(id)
}

/// A 2D degree-1 `B_SPLINE_CURVE_WITH_KNOTS` (a polyline in `(u,v)`).
pub(crate) fn emit_polyline2d(ctx: &mut Ctx, pts: &[(f64, f64)]) -> Result<usize, StepError> {
    let cps: Vec<usize> = pts
        .iter()
        .map(|&(u, v)| emit_point2d(ctx, u, v))
        .collect::<Result<_, _>>()?;
    let refs: String = cps
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    // degree 1: knots 0..(n-1), end multiplicities 2, interior 1.
    let n = cps.len();
    let knots: Vec<String> = (0..n).map(|k| format!("{:.6}", k as f64)).collect();
    let mut mult = vec![1usize; n];
    mult[0] = 2;
    mult[n - 1] = 2;
    // degree-1 knot count must be n+degree+1 = n+2; with first/last mult 2 and
    // interior 1 over (n) distinct knots that is 2+(n-2)+2 = n+2. ✓
    let mults: String = mult
        .iter()
        .map(|m| m.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let knots_s = knots.join(",");
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        &format!(
            "B_SPLINE_CURVE_WITH_KNOTS('',1,({refs}),.UNSPECIFIED.,.F.,.F.,({mults}),({knots_s}),.UNSPECIFIED.)"
        ),
    )?;
    Ok(id)
}

/// A 3D polyline emitted as a **degree-1 `B_SPLINE_CURVE_WITH_KNOTS`**.
///
/// Geometrically identical to a `POLYLINE`, but Parasolid-based importers
/// (SolidWorks) silently DROP faces whose edge curves are `POLYLINE` entities
/// — industrial writers (OCC, Parasolid) always emit degree-1 B-splines for
/// piecewise-linear edge geometry.
pub(crate) fn emit_polyline3d(ctx: &mut Ctx, pt_ids: &[usize]) -> Result<usize, StepError> {
    let refs: String = pt_ids
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let n = pt_ids.len();
    let knots: Vec<String> = (0..n).map(|k| format!("{:.6}", k as f64)).collect();
    let mut mult = vec![1usize; n];
    mult[0] = 2;
    mult[n - 1] = 2;
    let mults: String = mult
        .iter()
        .map(std::string::ToString::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let knots_s = knots.join(",");
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        &format!(
            "B_SPLINE_CURVE_WITH_KNOTS('',1,({refs}),.UNSPECIFIED.,.F.,.F.,({mults}),({knots_s}),.UNSPECIFIED.)"
        ),
    )?;
    Ok(id)
}

/// A `PCURVE` binding a 2D curve to a surface (via a DEFINITIONAL_REPRESENTATION).
pub(crate) fn emit_pcurve(
    ctx: &mut Ctx,
    surface_id: usize,
    curve2d_id: usize,
) -> Result<usize, StepError> {
    let ctx2d = param_context_2d(ctx)?;
    let def = ctx.next_id();
    ctx.emit_raw(
        def,
        &format!("DEFINITIONAL_REPRESENTATION('',(#{curve2d_id}),#{ctx2d})"),
    )?;
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!("PCURVE('',#{surface_id},#{def})"))?;
    Ok(id)
}

/// A `SURFACE_CURVE` (or `SEAM_CURVE` when `seam`) wrapping a 3D curve with its
/// associated pcurve(s).  A seam has exactly two pcurves on the same surface.
pub(crate) fn emit_surface_curve(
    ctx: &mut Ctx,
    c3d_id: usize,
    pcurves: &[usize],
    seam: bool,
) -> Result<usize, StepError> {
    let refs: String = pcurves
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let kind = if seam { "SEAM_CURVE" } else { "SURFACE_CURVE" };
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!("{kind}('',#{c3d_id},({refs}),.CURVE_3D.)"))?;
    Ok(id)
}

// ── Line3 → LINE ─────────────────────────────────────────────────────────────

pub(crate) fn emit_line(ctx: &mut Ctx, line: &Line3) -> Result<usize, StepError> {
    let pt_id = emit_point(ctx, line.origin, "line_pt")?;
    let dir_id = emit_unit_direction(ctx, line.direction, "line_dir")?;
    let vec_id = emit_vector(ctx, dir_id, 1.0, "line_vec")?;
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!("LINE('',#{pt_id},#{vec_id})"))?;
    Ok(id)
}

// ── Circle3 → CIRCLE ─────────────────────────────────────────────────────────

pub(crate) fn emit_circle(ctx: &mut Ctx, circ: &Circle3) -> Result<usize, StepError> {
    let ax_id = emit_axis2(ctx, circ.frame, "circ_ax")?;
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!("CIRCLE('',#{ax_id},{:.10})", circ.radius))?;
    Ok(id)
}

// ── Ellipse3 → ELLIPSE ───────────────────────────────────────────────────────

pub(crate) fn emit_ellipse(ctx: &mut Ctx, ellipse: &Ellipse3) -> Result<usize, StepError> {
    let ax_id = emit_axis2(ctx, ellipse.frame, "ellipse_ax")?;
    let id = ctx.next_id();
    ctx.emit_raw(
        id,
        &format!(
            "ELLIPSE('',#{ax_id},{:.10},{:.10})",
            ellipse.semi_major, ellipse.semi_minor
        ),
    )?;
    Ok(id)
}
