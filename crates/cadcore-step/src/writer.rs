//! High-level STEP writer — converts a [`BRep`] to an ISO 10303-21 string.
//!
//! ## Output format
//!
//! ```step
//! ISO-10303-21;
//! HEADER;
//!   FILE_DESCRIPTION(...);
//!   FILE_NAME(...);
//!   FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));
//! ENDSEC;
//! DATA;
//! #1 = APPLICATION_PROTOCOL_DEFINITION(...);
//! ...
//! ENDSEC;
//! END-ISO-10303-21;
//! ```

use std::fmt::Write as FmtWrite;

use cadcore_geom::{Circle3, Ellipse3};
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{BRep, FaceBoundary, FaceExtent, FaceGeom, SolidId};

use crate::entities::{
    emit_cylinder, emit_ellipse, emit_plane, emit_point, emit_sphere, emit_torus,
    emit_unit_direction, Ctx, StepError,
};

/// Builder that serialises a [`BRep`] to a STEP AP203 string.
pub struct StepWriter<'a> {
    brep: &'a BRep,
}

impl<'a> StepWriter<'a> {
    /// Create a new writer for `brep`.
    pub fn new(brep: &'a BRep) -> Self {
        Self { brep }
    }

    /// Serialise the entire B-Rep to a STEP string.
    ///
    /// If `solid_ids` is empty, all solids in the B-Rep are exported.
    /// Otherwise only the listed solids are exported.
    pub fn to_step(&self, solid_ids: &[SolidId]) -> Result<String, StepError> {
        let mut ctx = Ctx::new();

        // ── Preamble ───────────────────────────────────────────────────────
        let mut out = String::with_capacity(64 * 1024);
        out.push_str("ISO-10303-21;\n");
        out.push_str("HEADER;\n");
        out.push_str("  FILE_DESCRIPTION(('cadcore B-Rep export'),'2;1');\n");
        out.push_str("  FILE_NAME('','',(''),(''),'cadcore 0.1','','');\n");
        out.push_str("  FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
        out.push_str("ENDSEC;\n");
        out.push_str("DATA;\n");

        // ── Application context ────────────────────────────────────────────
        let app_ctx_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{app_ctx_id} = APPLICATION_CONTEXT('core data for automotive mechanical design processes');"
        )?;

        let apd_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{apd_id} = APPLICATION_PROTOCOL_DEFINITION('draft international standard','automotive_design',1998,#{app_ctx_id});"
        )?;

        let prod_ctx_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{prod_ctx_id} = PRODUCT_CONTEXT('',#{app_ctx_id},'mechanical');"
        )?;

        // ── Product / shape root ───────────────────────────────────────────
        let prod_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{prod_id} = PRODUCT('cadcore_part','cadcore part','',(#{prod_ctx_id}));"
        )?;

        let prod_def_form_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{prod_def_form_id} = PRODUCT_DEFINITION_FORMATION('','',#{prod_id});"
        )?;

        let prod_def_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{prod_def_id} = PRODUCT_DEFINITION('design','',#{prod_def_form_id},#{prod_ctx_id});"
        )?;

        let shape_def_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{shape_def_id} = PRODUCT_DEFINITION_SHAPE('','',#{prod_def_id});"
        )?;

        // ── Determine which solids to export ──────────────────────────────
        let solids_to_export: Vec<SolidId> = if solid_ids.is_empty() {
            self.brep.solids.keys().collect()
        } else {
            solid_ids.to_vec()
        };

        // ── Emit each solid ────────────────────────────────────────────────
        let mut shape_rep_items: Vec<usize> = Vec::new();

        for solid_id in &solids_to_export {
            if let Some(_solid) = self.brep.solids.get(*solid_id) {
                let adv_faces = self.emit_advanced_faces(&mut ctx, *solid_id)?;
                let af_refs: String = adv_faces
                    .iter()
                    .map(|id| format!("#{id}"))
                    .collect::<Vec<_>>()
                    .join(",");

                let csb_id = ctx.next_id();
                writeln!(ctx.out, "#{csb_id} = CLOSED_SHELL('',({af_refs}));")?;

                let msb_id = ctx.next_id();
                writeln!(
                    ctx.out,
                    "#{msb_id} = MANIFOLD_SOLID_BREP('brep',#{csb_id});"
                )?;
                shape_rep_items.push(msb_id);
            }
        }

        // ── Shape representation ───────────────────────────────────────────
        let geom_ctx_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{geom_ctx_id} = (GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{})) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{},#{},#{})) REPRESENTATION_CONTEXT('','3D'));\n",
            ctx.counter + 1,   // length_unit (will be emitted below — forward ref trick)
            ctx.counter + 2,
            ctx.counter + 3,
            ctx.counter + 4,
        )?;

        // Units
        let unc_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{unc_id} = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-007),#{},('',''));",
            ctx.counter + 1
        )?;
        let lu_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{lu_id} = (LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));"
        )?;
        let au_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{au_id} = (NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));"
        )?;
        let su_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{su_id} = (NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT());"
        )?;

        let items_str: String = shape_rep_items
            .iter()
            .map(|id| format!("#{id}"))
            .collect::<Vec<_>>()
            .join(",");

        let sr_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{sr_id} = SHAPE_REPRESENTATION('',({items_str}),#{geom_ctx_id});"
        )?;

        let srr_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{srr_id} = SHAPE_DEFINITION_REPRESENTATION(#{shape_def_id},#{sr_id});"
        )?;

        // ── Assemble ────────────────────────────────────────────────────────
        out.push_str(&ctx.out);
        out.push_str("ENDSEC;\n");
        out.push_str("END-ISO-10303-21;\n");
        Ok(out)
    }

    // ── Internal helpers ────────────────────────────────────────────────────

    /// Emit the carrier surface entity for a face and return its STEP id.
    fn emit_face_geometry(
        &self,
        ctx: &mut Ctx,
        face: &cadcore_topo::Face,
    ) -> Result<usize, StepError> {
        match &face.geom {
            FaceGeom::Plane(p) => emit_plane(ctx, p),
            FaceGeom::Cylinder(c) => emit_cylinder(ctx, c),
            FaceGeom::Sphere(s) => emit_sphere(ctx, s),
            FaceGeom::Torus(t) => emit_torus(ctx, t),
        }
    }

    /// Emit ADVANCED_FACE entities for all faces of a solid.
    fn emit_advanced_faces(
        &self,
        ctx: &mut Ctx,
        solid_id: SolidId,
    ) -> Result<Vec<usize>, StepError> {
        let solid = match self.brep.solids.get(solid_id) {
            Some(s) => s,
            None => return Ok(vec![]),
        };

        let mut ids = Vec::new();

        for &shell_id in &solid.shells {
            let shell = match self.brep.shells.get(shell_id) {
                Some(s) => s,
                None => continue,
            };
            for &face_id in &shell.faces {
                let face = match self.brep.faces.get(face_id) {
                    Some(f) => f,
                    None => continue,
                };

                let surf_id = self.emit_face_geometry(ctx, face)?;

                // Sense flag: .T. if normal agrees with surface, .F. if reversed
                let sense = match face.normal {
                    cadcore_topo::FaceNormal::Same => ".T.",
                    cadcore_topo::FaceNormal::Reversed => ".F.",
                };

                let bounds = emit_face_bounds(ctx, face)?;

                let af_id = ctx.next_id();
                if bounds.is_empty() {
                    // Fallback: emit with empty bounds (importer may reject, but
                    // prevents a write error for unexpected face types).
                    writeln!(
                        ctx.out,
                        "#{af_id} = ADVANCED_FACE('',(),(#{surf_id}),{sense});"
                    )?;
                } else {
                    let bound_refs: String = bounds
                        .iter()
                        .map(|id| format!("#{id}"))
                        .collect::<Vec<_>>()
                        .join(",");
                    writeln!(
                        ctx.out,
                        "#{af_id} = ADVANCED_FACE('',({bound_refs}),(#{surf_id}),{sense});"
                    )?;
                }
                ids.push(af_id);
            }
        }

        Ok(ids)
    }
}

// ── Face bound emission ────────────────────────────────────────────────────────

/// Emit the `FACE_OUTER_BOUND` (and optional `FACE_BOUND`) entities for `face`.
///
/// Returns the list of bound entity ids in order: outer first, then any inner
/// bounds.  For an unsupported extent variant the list is empty.
fn emit_face_bounds(ctx: &mut Ctx, face: &cadcore_topo::Face) -> Result<Vec<usize>, StepError> {
    match &face.extent {
        // ── Cylinder face: start/end circle or miter ellipse bounds ──
        FaceExtent::Cylinder { start, end, .. } => {
            let cyl = match &face.geom {
                FaceGeom::Cylinder(c) => c,
                _ => return Ok(vec![]),
            };

            let outer = emit_boundary(ctx, start, cyl.frame.x, true, true)?;
            let inner = emit_boundary(ctx, end, cyl.frame.x, false, false)?;
            Ok(vec![outer, inner])
        }

        // ── Planar disk: single outer bound ──────────────────────────────────────
        FaceExtent::Disk { radius } => {
            let plane = match &face.geom {
                FaceGeom::Plane(p) => p,
                _ => return Ok(vec![]),
            };
            let centre = plane.frame.origin;
            let normal = plane.frame.z;
            let x_dir = plane.frame.x;
            let bound = emit_circle_bound(ctx, centre, normal, x_dir, *radius, true, true)?;
            Ok(vec![bound])
        }

        // ── Torus fillet: two minor circles ──────────────────────────────────────
        FaceExtent::TorusFillet {
            start_circle,
            end_circle,
        } => {
            let s_x = start_circle.frame.x;
            let e_x = end_circle.frame.x;
            let outer = emit_circle_bound(
                ctx,
                start_circle.frame.origin,
                start_circle.frame.z,
                s_x,
                start_circle.radius,
                true,
                true,
            )?;
            let inner = emit_circle_bound(
                ctx,
                end_circle.frame.origin,
                end_circle.frame.z,
                e_x,
                end_circle.radius,
                false,
                false,
            )?;
            Ok(vec![outer, inner])
        }

        // ── Planar boundary: single circle or ellipse bound ─────────────────────
        FaceExtent::PlanarBoundary { boundary } => {
            let plane = match &face.geom {
                FaceGeom::Plane(p) => p,
                _ => return Ok(vec![]),
            };
            let outer = emit_boundary(ctx, boundary, plane.frame.x, true, true)?;
            Ok(vec![outer])
        }

        // ── Polygonal flat face: EDGE_LOOP of straight lines ─────────────────────
        FaceExtent::Polygon { points } => {
            if points.len() < 3 {
                return Ok(vec![]);
            }
            let _plane = match &face.geom {
                FaceGeom::Plane(p) => p,
                _ => return Ok(vec![]),
            };

            // Emit vertices
            let mut vtx_ids = Vec::with_capacity(points.len());
            for &pt in points {
                let vp_id = emit_point(ctx, pt, "v")?;
                let vtx_id = ctx.next_id();
                writeln!(ctx.out, "#{vtx_id} = VERTEX_POINT('',#{vp_id});")?;
                vtx_ids.push(vtx_id);
            }

            // Emit edges
            let mut oe_ids = Vec::with_capacity(points.len());
            let n = points.len();
            for i in 0..n {
                let p_start = points[i];
                let p_end = points[(i + 1) % n];
                let dir_vec = p_end - p_start;
                let len = dir_vec.length();
                if len < 1e-7 {
                    continue;
                }
                let dir = match UnitVec3::try_from_vec(dir_vec) {
                    Some(u) => u,
                    None => continue,
                };

                // Line placement
                let lp_id = emit_point(ctx, p_start, "lp")?;
                let ld_id = emit_unit_direction(ctx, dir, "ld")?;
                let line_id = ctx.next_id();
                writeln!(ctx.out, "#{line_id} = LINE('',#{lp_id},#{ld_id});")?;

                let v_start = vtx_ids[i];
                let v_end = vtx_ids[(i + 1) % n];

                let ec_id = ctx.next_id();
                writeln!(
                    ctx.out,
                    "#{ec_id} = EDGE_CURVE('',#{v_start},#{v_end},#{line_id},.T.);"
                )?;

                let oe_id = ctx.next_id();
                writeln!(
                    ctx.out,
                    "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},.T.);"
                )?;
                oe_ids.push(oe_id);
            }

            if oe_ids.is_empty() {
                return Ok(vec![]);
            }

            // EDGE_LOOP
            let el_id = ctx.next_id();
            let oe_refs = oe_ids
                .iter()
                .map(|id| format!("#{id}"))
                .collect::<Vec<_>>()
                .join(",");
            writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',({oe_refs}));")?;

            // FACE_OUTER_BOUND
            let fb_id = ctx.next_id();
            writeln!(ctx.out, "#{fb_id} = FACE_OUTER_BOUND('',#{el_id},.T.);")?;

            Ok(vec![fb_id])
        }

        // ── No extent info ────────────────────────────────────────────────────────
        FaceExtent::None => Ok(vec![]),
    }
}

fn emit_boundary(
    ctx: &mut Ctx,
    boundary: &FaceBoundary,
    fallback_x_dir: UnitVec3,
    outer: bool,
    orient: bool,
) -> Result<usize, StepError> {
    match boundary {
        FaceBoundary::Circle(c) => emit_circle_bound(
            ctx,
            c.frame.origin,
            c.frame.z,
            c.frame.x,
            c.radius,
            outer,
            orient,
        ),
        FaceBoundary::Ellipse(e) => emit_ellipse_bound(ctx, e, fallback_x_dir, outer, orient),
    }
}

/// Emit a complete face-bound for a full 360° circle.
///
/// Generates:
/// ```step
/// #N   = CARTESIAN_POINT(...)    /* circle centre */
/// #N+1 = DIRECTION(...)          /* circle normal */
/// #N+2 = DIRECTION(...)          /* circle x-ref */
/// #N+3 = AXIS2_PLACEMENT_3D(...)
/// #N+4 = CIRCLE(...)
/// #N+5 = CARTESIAN_POINT(...)    /* vertex point at θ=0 */
/// #N+6 = VERTEX_POINT(...)
/// #N+7 = EDGE_CURVE(...)         /* start == end → closed */
/// #N+8 = ORIENTED_EDGE(...)
/// #N+9 = EDGE_LOOP(...)
/// #N+10= FACE_OUTER_BOUND/FACE_BOUND(...)
/// ```
///
/// Returns the id of the `FACE_OUTER_BOUND` or `FACE_BOUND` entity.
fn emit_circle_bound(
    ctx: &mut Ctx,
    centre: Point3,
    normal: UnitVec3,
    x_dir: UnitVec3,
    radius: f64,
    outer: bool,  // true = FACE_OUTER_BOUND, false = FACE_BOUND
    orient: bool, // ORIENTED_EDGE orientation: .T. = forward, .F. = reversed
) -> Result<usize, StepError> {
    // Circle placement
    let cp_id = emit_point(ctx, centre, "c")?;
    let cz_id = emit_unit_direction(ctx, normal, "cn")?;
    let cx_id = emit_unit_direction(ctx, x_dir, "cx")?;
    let cax_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{cax_id} = AXIS2_PLACEMENT_3D('',#{cp_id},#{cz_id},#{cx_id});"
    )?;

    let circ_id = ctx.next_id();
    writeln!(ctx.out, "#{circ_id} = CIRCLE('',#{cax_id},{:.10});", radius)?;

    // Vertex point on the circle at θ=0  (centre + x_dir * radius)
    let vp_world = Point3::new(
        centre.x + x_dir.as_vec().x * radius,
        centre.y + x_dir.as_vec().y * radius,
        centre.z + x_dir.as_vec().z * radius,
    );
    let vpt_id = emit_point(ctx, vp_world, "v")?;
    let vtx_id = ctx.next_id();
    writeln!(ctx.out, "#{vtx_id} = VERTEX_POINT('',#{vpt_id});")?;

    // EDGE_CURVE: start == end (closed full circle)
    let ec_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{ec_id} = EDGE_CURVE('',#{vtx_id},#{vtx_id},#{circ_id},.T.);"
    )?;

    // ORIENTED_EDGE
    let sense = if orient { ".T." } else { ".F." };
    let oe_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});"
    )?;

    // EDGE_LOOP
    let el_id = ctx.next_id();
    writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',(#{oe_id}));")?;

    // FACE_OUTER_BOUND or FACE_BOUND
    let fb_id = ctx.next_id();
    let bound_type = if outer {
        "FACE_OUTER_BOUND"
    } else {
        "FACE_BOUND"
    };
    writeln!(ctx.out, "#{fb_id} = {bound_type}('',#{el_id},.T.);")?;

    Ok(fb_id)
}

fn emit_ellipse_bound(
    ctx: &mut Ctx,
    ellipse: &Ellipse3,
    _fallback_x_dir: UnitVec3,
    outer: bool,
    orient: bool,
) -> Result<usize, StepError> {
    let ellipse_id = emit_ellipse(ctx, ellipse)?;

    let vp_world = Point3::new(
        ellipse.frame.origin.x + ellipse.frame.x.as_vec().x * ellipse.semi_major,
        ellipse.frame.origin.y + ellipse.frame.x.as_vec().y * ellipse.semi_major,
        ellipse.frame.origin.z + ellipse.frame.x.as_vec().z * ellipse.semi_major,
    );
    let vpt_id = emit_point(ctx, vp_world, "v")?;
    let vtx_id = ctx.next_id();
    writeln!(ctx.out, "#{vtx_id} = VERTEX_POINT('',#{vpt_id});")?;

    let ec_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{ec_id} = EDGE_CURVE('',#{vtx_id},#{vtx_id},#{ellipse_id},.T.);"
    )?;

    let sense = if orient { ".T." } else { ".F." };
    let oe_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});"
    )?;

    let el_id = ctx.next_id();
    writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',(#{oe_id}));")?;

    let fb_id = ctx.next_id();
    let bound_type = if outer {
        "FACE_OUTER_BOUND"
    } else {
        "FACE_BOUND"
    };
    writeln!(ctx.out, "#{fb_id} = {bound_type}('',#{el_id},.T.);")?;

    Ok(fb_id)
}

// ── Convenience function ──────────────────────────────────────────────────────

/// Write the entire B-Rep to a STEP string, exporting all solids.
pub fn brep_to_step(brep: &BRep) -> Result<String, StepError> {
    StepWriter::new(brep).to_step(&[])
}

// ── Helper re-exports (used by writer only) ───────────────────────────────────

// Ensure Circle3 is available through the emit helpers.
fn _use_circle3(_c: Circle3) {}

#[cfg(test)]
mod tests {
    use super::*;
    use cadcore_math::Point3;
    use cadcore_ops::{sweep_circle_along_polyline, SweepOptions};

    #[test]
    fn simple_rod_produces_valid_step_header() {
        let mut brep = BRep::new();
        let waypoints = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 10.0)];
        sweep_circle_along_polyline(
            &mut brep,
            &waypoints,
            0.5,
            &SweepOptions {
                fillet_corners: true,
                ..SweepOptions::default()
            },
        )
        .unwrap();

        let step = brep_to_step(&brep).unwrap();
        assert!(step.starts_with("ISO-10303-21;"), "missing ISO header");
        assert!(
            step.contains("CYLINDRICAL_SURFACE"),
            "missing cylinder surface"
        );
        assert!(step.contains("PLANE"), "missing plane (end caps)");
        assert!(
            step.contains("FACE_OUTER_BOUND"),
            "missing face outer bound"
        );
        assert!(step.contains("EDGE_LOOP"), "missing edge loop");
        assert!(step.contains("EDGE_CURVE"), "missing edge curve");
        assert!(step.contains("VERTEX_POINT"), "missing vertex point");
        assert!(step.contains("END-ISO-10303-21;"), "missing footer");
    }

    #[test]
    fn simple_rod_no_empty_face_bounds() {
        let mut brep = BRep::new();
        let waypoints = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 10.0)];
        sweep_circle_along_polyline(
            &mut brep,
            &waypoints,
            0.5,
            &SweepOptions {
                fillet_corners: true,
                ..SweepOptions::default()
            },
        )
        .unwrap();

        let step = brep_to_step(&brep).unwrap();
        // No ADVANCED_FACE with empty bounds ()
        assert!(
            !step.contains("ADVANCED_FACE('',(),"),
            "found ADVANCED_FACE with empty bounds"
        );
    }

    #[test]
    fn bent_rod_includes_torus() {
        let mut brep = BRep::new();
        let waypoints = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(5.0, 5.0, 0.0),
        ];
        sweep_circle_along_polyline(
            &mut brep,
            &waypoints,
            0.5,
            &SweepOptions {
                fillet_corners: true,
                ..SweepOptions::default()
            },
        )
        .unwrap();

        let step = brep_to_step(&brep).unwrap();
        assert!(
            step.contains("TOROIDAL_SURFACE"),
            "missing torus corner fillet"
        );
    }

    #[test]
    fn sharp_bent_rod_uses_miter_ellipses() {
        let mut brep = BRep::new();
        let waypoints = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(5.0, 5.0, 0.0),
        ];
        sweep_circle_along_polyline(
            &mut brep,
            &waypoints,
            0.5,
            &SweepOptions {
                fillet_corners: false,
                ..SweepOptions::default()
            },
        )
        .unwrap();

        let step = brep_to_step(&brep).unwrap();
        assert!(step.contains("ELLIPSE"), "missing miter ellipse");
        assert!(
            !step.contains("TOROIDAL_SURFACE"),
            "unexpected torus for sharp miter"
        );
    }
}
