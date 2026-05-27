//! Low-level STEP entity emission helpers.
#![allow(dead_code)]
//!
//! The STEP writer maintains a monotonically increasing entity counter.
//! Every entity is emitted as `#N = ENTITY_TYPE(args);`.

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
    fn from(e: std::fmt::Error) -> Self { Self::Fmt(e) }
}

impl std::fmt::Display for StepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { Self::Fmt(e) => write!(f, "fmt error: {e}") }
    }
}
impl std::error::Error for StepError {}

/// Mutable counter + output buffer, shared across all entity emission calls.
pub(crate) struct Ctx {
    pub counter: usize,
    pub out:     String,
}

impl Ctx {
    pub fn new() -> Self { Self { counter: 0, out: String::with_capacity(8192) } }

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
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!(
        "CARTESIAN_POINT('{label}',({:.10},{:.10},{:.10}))",
        p.x, p.y, p.z
    ))?;
    Ok(id)
}

// ── Vec3 → DIRECTION ─────────────────────────────────────────────────────────

pub(crate) fn emit_direction(ctx: &mut Ctx, v: Vec3, label: &str) -> Result<usize, StepError> {
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!(
        "DIRECTION('{label}',({:.10},{:.10},{:.10}))",
        v.x, v.y, v.z
    ))?;
    Ok(id)
}

pub(crate) fn emit_unit_direction(ctx: &mut Ctx, u: UnitVec3, label: &str) -> Result<usize, StepError> {
    emit_direction(ctx, u.as_vec(), label)
}

// ── VECTOR (direction + magnitude) ──────────────────────────────────────────

pub(crate) fn emit_vector(ctx: &mut Ctx, dir_id: usize, magnitude: f64, label: &str) -> Result<usize, StepError> {
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!(
        "VECTOR('{label}',#{dir_id},{magnitude:.10})"
    ))?;
    Ok(id)
}

// ── Frame3 → AXIS2_PLACEMENT_3D ──────────────────────────────────────────────

pub(crate) fn emit_axis2(ctx: &mut Ctx, frame: Frame3, label: &str) -> Result<usize, StepError> {
    let origin_id = emit_point(ctx, frame.origin, label)?;
    let z_id      = emit_unit_direction(ctx, frame.z, &format!("{label}_Z"))?;
    let x_id      = emit_unit_direction(ctx, frame.x, &format!("{label}_X"))?;
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!(
        "AXIS2_PLACEMENT_3D('{label}',#{origin_id},#{z_id},#{x_id})"
    ))?;
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
    ctx.emit_raw(id, &format!("CYLINDRICAL_SURFACE('',#{ax_id},{:.10})", cyl.radius))?;
    Ok(id)
}

// ── SphereSurf → SPHERICAL_SURFACE ────────────────────────────────────────────

pub(crate) fn emit_sphere(ctx: &mut Ctx, sph: &SphereSurf) -> Result<usize, StepError> {
    let ax_id = emit_axis2(ctx, sph.frame, "sph_axis")?;
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!("SPHERICAL_SURFACE('',#{ax_id},{:.10})", sph.radius))?;
    Ok(id)
}

// ── TorusSurf → TOROIDAL_SURFACE ──────────────────────────────────────────────

pub(crate) fn emit_torus(ctx: &mut Ctx, torus: &TorusSurf) -> Result<usize, StepError> {
    let ax_id = emit_axis2(ctx, torus.frame, "torus_axis")?;
    let id = ctx.next_id();
    ctx.emit_raw(id, &format!(
        "TOROIDAL_SURFACE('',#{ax_id},{:.10},{:.10})",
        torus.major_radius, torus.minor_radius
    ))?;
    Ok(id)
}

// ── Line3 → LINE ─────────────────────────────────────────────────────────────

pub(crate) fn emit_line(ctx: &mut Ctx, line: &Line3) -> Result<usize, StepError> {
    let pt_id  = emit_point(ctx, line.origin, "line_pt")?;
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
    ctx.emit_raw(id, &format!(
        "ELLIPSE('',#{ax_id},{:.10},{:.10})",
        ellipse.semi_major, ellipse.semi_minor
    ))?;
    Ok(id)
}
