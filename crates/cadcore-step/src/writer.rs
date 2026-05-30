//! High-level STEP writer — converts a [`BRep`] to an ISO 10303-21 string.
//!
//! Extended with `PartialCylinder` and `PartialDisk` for solid half-space cut output.

use std::fmt::Write as FmtWrite;

use cadcore_geom::Ellipse3;
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{BRep, FaceBoundary, FaceExtent, FaceGeom, SolidId};

use crate::entities::{
    emit_cylinder, emit_ellipse, emit_plane, emit_point, emit_sphere, emit_torus,
    emit_unit_direction, Ctx, StepError, emit_vertex_point, StepCurveKey, point_key,
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
    pub fn to_step(&self, solid_ids: &[SolidId]) -> Result<String, StepError> {
        let mut ctx = Ctx::new();

        let mut out = String::with_capacity(64 * 1024);
        out.push_str("ISO-10303-21;\n");
        out.push_str("HEADER;\n");
        out.push_str("  FILE_DESCRIPTION(('cadcore B-Rep export'),'2;1');\n");
        out.push_str("  FILE_NAME('','',(''),(''),'cadcore 0.1','','');\n");
        out.push_str("  FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
        out.push_str("ENDSEC;\n");
        out.push_str("DATA;\n");

        let app_ctx_id = ctx.next_id();
        writeln!(ctx.out, "#{app_ctx_id} = APPLICATION_CONTEXT('core data for automotive mechanical design processes');")?;
        let apd_id = ctx.next_id();
        writeln!(ctx.out, "#{apd_id} = APPLICATION_PROTOCOL_DEFINITION('draft international standard','automotive_design',1998,#{app_ctx_id});")?;
        let prod_ctx_id = ctx.next_id();
        writeln!(ctx.out, "#{prod_ctx_id} = PRODUCT_CONTEXT('',#{app_ctx_id},'mechanical');")?;
        let prod_id = ctx.next_id();
        writeln!(ctx.out, "#{prod_id} = PRODUCT('cadcore_part','cadcore part','',(#{prod_ctx_id}));")?;
        let prod_def_form_id = ctx.next_id();
        writeln!(ctx.out, "#{prod_def_form_id} = PRODUCT_DEFINITION_FORMATION('','',#{prod_id});")?;
        let prod_def_id = ctx.next_id();
        writeln!(ctx.out, "#{prod_def_id} = PRODUCT_DEFINITION('design','',#{prod_def_form_id},#{prod_ctx_id});")?;
        let shape_def_id = ctx.next_id();
        writeln!(ctx.out, "#{shape_def_id} = PRODUCT_DEFINITION_SHAPE('','',#{prod_def_id});")?;

        let solids_to_export: Vec<SolidId> = if solid_ids.is_empty() {
            self.brep.solids.keys().collect()
        } else {
            solid_ids.to_vec()
        };

        let mut shape_rep_items: Vec<usize> = Vec::new();

        for solid_id in &solids_to_export {
            if let Some(_solid) = self.brep.solids.get(*solid_id) {
                let adv_faces = self.emit_advanced_faces(&mut ctx, *solid_id)?;
                let af_refs: String = adv_faces.iter().map(|id| format!("#{id}")).collect::<Vec<_>>().join(",");

                let csb_id = ctx.next_id();
                writeln!(ctx.out, "#{csb_id} = CLOSED_SHELL('',({af_refs}));")?;
                let msb_id = ctx.next_id();
                writeln!(ctx.out, "#{msb_id} = MANIFOLD_SOLID_BREP('brep',#{csb_id});")?;
                shape_rep_items.push(msb_id);
            }
        }

        let geom_ctx_id = ctx.next_id();
        writeln!(ctx.out,
            "#{geom_ctx_id} = (GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{})) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{},#{},#{})) REPRESENTATION_CONTEXT('','3D'));\n",
            ctx.counter + 1, ctx.counter + 2, ctx.counter + 3, ctx.counter + 4,
        )?;
        let unc_id = ctx.next_id();
        writeln!(ctx.out, "#{unc_id} = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-007),#{},('',''));", ctx.counter + 1)?;
        let lu_id = ctx.next_id();
        writeln!(ctx.out, "#{lu_id} = (LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));")?;
        let au_id = ctx.next_id();
        writeln!(ctx.out, "#{au_id} = (NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));")?;
        let su_id = ctx.next_id();
        writeln!(ctx.out, "#{su_id} = (NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT());")?;

        let items_str: String = shape_rep_items.iter().map(|id| format!("#{id}")).collect::<Vec<_>>().join(",");
        let sr_id = ctx.next_id();
        writeln!(ctx.out, "#{sr_id} = SHAPE_REPRESENTATION('',({items_str}),#{geom_ctx_id});")?;
        let srr_id = ctx.next_id();
        writeln!(ctx.out, "#{srr_id} = SHAPE_DEFINITION_REPRESENTATION(#{shape_def_id},#{sr_id});")?;

        out.push_str(&ctx.out);
        out.push_str("ENDSEC;\n");
        out.push_str("END-ISO-10303-21;\n");
        Ok(out)
    }

    fn emit_face_geometry(&self, ctx: &mut Ctx, face: &cadcore_topo::Face) -> Result<usize, StepError> {
        match &face.geom {
            FaceGeom::Plane(p)    => emit_plane(ctx, p),
            FaceGeom::Cylinder(c) => emit_cylinder(ctx, c),
            FaceGeom::Sphere(s)   => emit_sphere(ctx, s),
            FaceGeom::Torus(t)    => emit_torus(ctx, t),
        }
    }
    fn emit_advanced_faces(&self, ctx: &mut Ctx, solid_id: SolidId) -> Result<Vec<usize>, StepError> {
        let solid = match self.brep.solids.get(solid_id) { Some(s) => s, None => return Ok(vec![]) };
        let mut ids = Vec::new();
        for &shell_id in &solid.shells {
            ctx.point_cache.clear();
            ctx.vertex_cache.clear();
            ctx.edge_cache.clear();

            let shell = match self.brep.shells.get(shell_id) { Some(s) => s, None => continue };
            for &face_id in &shell.faces {
                let face = match self.brep.faces.get(face_id) { Some(f) => f, None => continue };
                let surf_id = self.emit_face_geometry(ctx, face)?;
                let sense = match face.normal {
                    cadcore_topo::FaceNormal::Same     => ".T.",
                    cadcore_topo::FaceNormal::Reversed => ".F.",
                };
                let bounds = emit_face_bounds(ctx, face)?;
                let af_id = ctx.next_id();
                if bounds.is_empty() {
                    writeln!(ctx.out, "#{af_id} = ADVANCED_FACE('',(),(#{surf_id}),{sense});")?;
                } else {
                    let bound_refs: String = bounds.iter().map(|id| format!("#{id}")).collect::<Vec<_>>().join(",");
                    writeln!(ctx.out, "#{af_id} = ADVANCED_FACE('',({bound_refs}),(#{surf_id}),{sense});")?;
                }
                ids.push(af_id);
            }
        }
        Ok(ids)
    }
}

// ── Face bound emission ────────────────────────────────────────────────────────

fn emit_face_bounds(ctx: &mut Ctx, face: &cadcore_topo::Face) -> Result<Vec<usize>, StepError> {
    match &face.extent {
        // ── Full cylinder ─────────────────────────────────────────────────────
        FaceExtent::Cylinder { start, end, .. } => {
            let cyl = match &face.geom { FaceGeom::Cylinder(c) => c, _ => return Ok(vec![]) };
            let outer = emit_boundary(ctx, start, cyl.frame.x, true, true)?;
            let inner = emit_boundary(ctx, end, cyl.frame.x, false, false)?;
            Ok(vec![outer, inner])
        }

        // ── Partial cylinder (chord cut) ──────────────────────────────────────
        FaceExtent::PartialCylinder { length, arc_start_angle, arc_end_angle, arc_ref_dir } => {
            let cyl = match &face.geom { FaceGeom::Cylinder(c) => c, _ => return Ok(vec![]) };
            emit_partial_cylinder_bounds(
                ctx,
                cyl.frame.origin,
                cyl.frame.z,
                cyl.frame.x,
                cyl.radius,
                *length,
                *arc_start_angle,
                *arc_end_angle,
                *arc_ref_dir,
            )
        }

        // ── Planar disk (full circle cap) ─────────────────────────────────────
        FaceExtent::Disk { radius } => {
            let plane = match &face.geom { FaceGeom::Plane(p) => p, _ => return Ok(vec![]) };
            let bound = emit_circle_bound(ctx, plane.frame.origin, plane.frame.z, plane.frame.x, *radius, true, true)?;
            Ok(vec![bound])
        }

        // ── Partial disk (arc + chord, end cap of partial cylinder) ───────────
        FaceExtent::PartialDisk { radius, start_angle, end_angle } => {
            let plane = match &face.geom { FaceGeom::Plane(p) => p, _ => return Ok(vec![]) };
            emit_partial_disk_bound(
                ctx,
                plane.frame.origin,
                plane.frame.z,
                plane.frame.x,
                *radius,
                *start_angle,
                *end_angle,
            )
        }

        // ── Torus fillet arc: boundary consists of two minor circles ───────────
        FaceExtent::TorusFillet { start_circle, end_circle } => {
            let s_x = start_circle.frame.x;
            let e_x = end_circle.frame.x;
            let outer = emit_circle_bound(ctx, start_circle.frame.origin, start_circle.frame.z, s_x, start_circle.radius, true, true)?;
            let inner = emit_circle_bound(ctx, end_circle.frame.origin, end_circle.frame.z, e_x, end_circle.radius, false, false)?;
            Ok(vec![outer, inner])
        }

        // ── Planar boundary (circle or ellipse cap) ───────────────────────────
        FaceExtent::PlanarBoundary { boundary } => {
            let plane = match &face.geom { FaceGeom::Plane(p) => p, _ => return Ok(vec![]) };
            let outer = emit_boundary(ctx, boundary, plane.frame.x, true, true)?;
            Ok(vec![outer])
        }

        // ── Polygonal flat face ───────────────────────────────────────────────
        FaceExtent::Polygon { points } => {
            if points.len() < 3 { return Ok(vec![]); }
            let mut vtx_ids = Vec::with_capacity(points.len());
            for &pt in points {
                let vtx_id = emit_vertex_point(ctx, pt)?;
                vtx_ids.push(vtx_id);
            }
            let mut oe_ids = Vec::with_capacity(points.len());
            let n = points.len();
            for i in 0..n {
                let p_start = points[i];
                let p_end   = points[(i + 1) % n];
                let dir_vec = p_end - p_start;
                if dir_vec.length() < 1e-7 { continue; }
                let dir = match UnitVec3::try_from_vec(dir_vec) { Some(u) => u, None => continue };
                let v_s = vtx_ids[i];
                let v_e = vtx_ids[(i + 1) % n];
                
                let v_min = v_s.min(v_e);
                let v_max = v_s.max(v_e);
                let key = StepCurveKey::Line { v1: v_min, v2: v_max };
                let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
                    pair
                } else {
                    let lp_id   = emit_point(ctx, p_start, "lp")?;
                    let ld_id   = emit_unit_direction(ctx, dir, "ld")?;
                    let line_id = ctx.next_id();
                    writeln!(ctx.out, "#{line_id} = LINE('',#{lp_id},#{ld_id});")?;
                    let ec_id   = ctx.next_id();
                    writeln!(ctx.out, "#{ec_id} = EDGE_CURVE('',#{v_s},#{v_e},#{line_id},.T.);")?;
                    ctx.edge_cache.insert(key, (ec_id, v_s));
                    (ec_id, v_s)
                };
                let orient = orig_start == v_s;
                let sense = if orient { ".T." } else { ".F." };
                let oe_id = ctx.next_id();
                writeln!(ctx.out, "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
                oe_ids.push(oe_id);
            }
            if oe_ids.is_empty() { return Ok(vec![]); }
            let el_id  = ctx.next_id();
            let oe_refs = oe_ids.iter().map(|id| format!("#{id}")).collect::<Vec<_>>().join(",");
            writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',({oe_refs}));")?;
            let fb_id  = ctx.next_id();
            writeln!(ctx.out, "#{fb_id} = FACE_OUTER_BOUND('',#{el_id},.T.);")?;
            Ok(vec![fb_id])
        }

        FaceExtent::None => Ok(vec![]),
    }
}

// ── Partial cylinder bound emission ───────────────────────────────────────────

fn emit_partial_cylinder_bounds(
    ctx:             &mut Ctx,
    origin:          Point3,
    axis:            UnitVec3,
    x_ref:           UnitVec3,
    radius:          f64,
    length:          f64,
    arc_start_angle: f64,
    arc_end_angle:   f64,
    arc_ref_dir:     UnitVec3,
) -> Result<Vec<usize>, StepError> {
    let right = match UnitVec3::try_from_vec(axis.cross(arc_ref_dir)) {
        Some(u) => u,
        None => x_ref,
    };

    let n_segs = 24usize;
    let angle_range = arc_end_angle - arc_start_angle;
    let mut pts_start: Vec<Point3> = Vec::with_capacity(n_segs + 1);
    let mut pts_end:   Vec<Point3> = Vec::with_capacity(n_segs + 1);

    let end = origin + axis.as_vec() * length;

    for i in 0..=n_segs {
        let t = i as f64 / n_segs as f64;
        let angle = arc_start_angle + t * angle_range;
        let local = arc_ref_dir.as_vec() * (radius * angle.cos())
                  + right.as_vec()       * (radius * angle.sin());
        pts_start.push(origin + local);
        pts_end.push(end + local);
    }

    let mut all_pts: Vec<Point3> = Vec::new();
    all_pts.extend_from_slice(&pts_start);
    all_pts.extend(pts_end.iter().rev().cloned());

    if (*all_pts.first().unwrap() - *all_pts.last().unwrap()).length() < 1e-7 {
        all_pts.pop();
    }

    if all_pts.len() < 3 { return Ok(vec![]); }

    let mut vtx_ids = Vec::with_capacity(all_pts.len());
    for &pt in &all_pts {
        let vtx_id = emit_vertex_point(ctx, pt)?;
        vtx_ids.push(vtx_id);
    }

    let mut oe_ids = Vec::new();
    let n = all_pts.len();
    for i in 0..n {
        let p_start = all_pts[i];
        let p_end   = all_pts[(i + 1) % n];
        let dv = p_end - p_start;
        if dv.length() < 1e-7 { continue; }
        let dir = match UnitVec3::try_from_vec(dv) { Some(u) => u, None => continue };
        let v_s = vtx_ids[i];
        let v_e = vtx_ids[(i + 1) % n];
        
        let v_min = v_s.min(v_e);
        let v_max = v_s.max(v_e);
        let key = StepCurveKey::Line { v1: v_min, v2: v_max };
        let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
            pair
        } else {
            let lp_id   = emit_point(ctx, p_start, "lp")?;
            let ld_id   = emit_unit_direction(ctx, dir, "ld")?;
            let line_id = ctx.next_id();
            writeln!(ctx.out, "#{line_id} = LINE('',#{lp_id},#{ld_id});")?;
            let ec_id   = ctx.next_id();
            writeln!(ctx.out, "#{ec_id} = EDGE_CURVE('',#{v_s},#{v_e},#{line_id},.T.);")?;
            ctx.edge_cache.insert(key, (ec_id, v_s));
            (ec_id, v_s)
        };
        let orient = orig_start == v_s;
        let sense = if orient { ".T." } else { ".F." };
        let oe_id = ctx.next_id();
        writeln!(ctx.out, "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
        oe_ids.push(oe_id);
    }
    if oe_ids.is_empty() { return Ok(vec![]); }
    let oe_refs = oe_ids.iter().map(|id| format!("#{id}")).collect::<Vec<_>>().join(",");
    let el_id = ctx.next_id();
    writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',({oe_refs}));")?;
    let fb_id = ctx.next_id();
    writeln!(ctx.out, "#{fb_id} = FACE_OUTER_BOUND('',#{el_id},.T.);")?;
    Ok(vec![fb_id])
}

// ── Partial disk bound emission ───────────────────────────────────────────────

fn emit_partial_disk_bound(
    ctx:         &mut Ctx,
    centre:      Point3,
    normal:      UnitVec3,
    x_ref:       UnitVec3,
    radius:      f64,
    start_angle: f64,
    end_angle:   f64,
) -> Result<Vec<usize>, StepError> {
    let n_segs = 24usize;
    let angle_range = end_angle - start_angle;

    let y_ref = match UnitVec3::try_from_vec(normal.cross(x_ref)) {
        Some(u) => u,
        None => return Ok(vec![]),
    };

    let mut pts: Vec<Point3> = Vec::with_capacity(n_segs + 2);
    for i in 0..=n_segs {
        let t     = i as f64 / n_segs as f64;
        let angle = start_angle + t * angle_range;
        let local = x_ref.as_vec() * (radius * angle.cos())
                  + y_ref.as_vec() * (radius * angle.sin());
        pts.push(centre + local);
    }

    let mut vtx_ids = Vec::with_capacity(pts.len());
    for &pt in &pts {
        let vtx_id = emit_vertex_point(ctx, pt)?;
        vtx_ids.push(vtx_id);
    }

    let mut oe_ids = Vec::new();
    let n = pts.len();
    for i in 0..n {
        let p_start = pts[i];
        let p_end   = pts[(i + 1) % n];
        let dv = p_end - p_start;
        if dv.length() < 1e-7 { continue; }
        let dir = match UnitVec3::try_from_vec(dv) { Some(u) => u, None => continue };
        let v_s = vtx_ids[i];
        let v_e = vtx_ids[(i + 1) % n];
        
        let v_min = v_s.min(v_e);
        let v_max = v_s.max(v_e);
        let key = StepCurveKey::Line { v1: v_min, v2: v_max };
        let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
            pair
        } else {
            let lp_id   = emit_point(ctx, p_start, "lp")?;
            let ld_id   = emit_unit_direction(ctx, dir, "ld")?;
            let line_id = ctx.next_id();
            writeln!(ctx.out, "#{line_id} = LINE('',#{lp_id},#{ld_id});")?;
            let ec_id   = ctx.next_id();
            writeln!(ctx.out, "#{ec_id} = EDGE_CURVE('',#{v_s},#{v_e},#{line_id},.T.);")?;
            ctx.edge_cache.insert(key, (ec_id, v_s));
            (ec_id, v_s)
        };
        let orient = orig_start == v_s;
        let sense = if orient { ".T." } else { ".F." };
        let oe_id = ctx.next_id();
        writeln!(ctx.out, "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
        oe_ids.push(oe_id);
    }
    if oe_ids.is_empty() { return Ok(vec![]); }
    let oe_refs = oe_ids.iter().map(|id| format!("#{id}")).collect::<Vec<_>>().join(",");
    let el_id = ctx.next_id();
    writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',({oe_refs}));")?;
    let fb_id = ctx.next_id();
    writeln!(ctx.out, "#{fb_id} = FACE_OUTER_BOUND('',#{el_id},.T.);")?;
    Ok(vec![fb_id])
}

// ── Standard helpers (copied from original) ───────────────────────────────────

fn emit_boundary(
    ctx: &mut Ctx,
    boundary: &FaceBoundary,
    fallback_x_dir: UnitVec3,
    outer: bool,
    orient: bool,
) -> Result<usize, StepError> {
    match boundary {
        FaceBoundary::Circle(c) => emit_circle_bound(ctx, c.frame.origin, c.frame.z, c.frame.x, c.radius, outer, orient),
        FaceBoundary::Ellipse(e) => emit_ellipse_bound(ctx, e, fallback_x_dir, outer, orient),
    }
}

fn emit_circle_bound(
    ctx:    &mut Ctx,
    centre: Point3,
    normal: UnitVec3,
    _x_dir:  UnitVec3,
    radius: f64,
    outer:  bool,
    orient: bool,
) -> Result<usize, StepError> {
    let key = StepCurveKey::Circle {
        center: point_key(centre),
        radius_micro: (radius * 1_000_000.0).round() as i64,
        normal: normal_key(normal),
    };
    let ec_id = if let Some(&(id, _)) = ctx.edge_cache.get(&key) {
        id
    } else {
        // Deterministically compute x_dir orthogonal to normal
        let n_vec = normal.as_vec();
        let ref_vec = if n_vec.x.abs() > 0.9 {
            cadcore_math::Vec3::new(0.0, 1.0, 0.0)
        } else {
            cadcore_math::Vec3::new(1.0, 0.0, 0.0)
        };
        let ortho = n_vec.cross(ref_vec);
        let x_dir = UnitVec3::try_from_vec(ortho).unwrap_or(normal);

        let cp_id  = emit_point(ctx, centre, "c")?;
        let cz_id  = emit_unit_direction(ctx, normal, "cn")?;
        let cx_id  = emit_unit_direction(ctx, x_dir, "cx")?;
        let cax_id = ctx.next_id();
        writeln!(ctx.out, "#{cax_id} = AXIS2_PLACEMENT_3D('',#{cp_id},#{cz_id},#{cx_id});")?;
        let circ_id = ctx.next_id();
        writeln!(ctx.out, "#{circ_id} = CIRCLE('',#{cax_id},{:.10});", radius)?;
        let vp_world = centre + x_dir.as_vec() * radius;
        let vtx_id = emit_vertex_point(ctx, vp_world)?;
        let ec_id  = ctx.next_id();
        writeln!(ctx.out, "#{ec_id} = EDGE_CURVE('',#{vtx_id},#{vtx_id},#{circ_id},.T.);")?;
        ctx.edge_cache.insert(key, (ec_id, vtx_id));
        ec_id
    };
    let sense  = if orient { ".T." } else { ".F." };
    let oe_id  = ctx.next_id();
    writeln!(ctx.out, "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
    let el_id  = ctx.next_id();
    writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',(#{oe_id}));")?;
    let fb_id  = ctx.next_id();
    let btype  = if outer { "FACE_OUTER_BOUND" } else { "FACE_BOUND" };
    writeln!(ctx.out, "#{fb_id} = {btype}('',#{el_id},.T.);")?;
    Ok(fb_id)
}

fn emit_ellipse_bound(
    ctx: &mut Ctx,
    ellipse: &Ellipse3,
    _fallback_x_dir: UnitVec3,
    outer: bool,
    orient: bool,
) -> Result<usize, StepError> {
    let key = StepCurveKey::Ellipse {
        center: point_key(ellipse.frame.origin),
        semi_major_micro: (ellipse.semi_major * 1_000_000.0).round() as i64,
        semi_minor_micro: (ellipse.semi_minor * 1_000_000.0).round() as i64,
        normal: normal_key(ellipse.frame.z),
    };
    let ec_id = if let Some(&(id, _)) = ctx.edge_cache.get(&key) {
        id
    } else {
        let ellipse_id = emit_ellipse(ctx, ellipse)?;
        let vp_world = ellipse.frame.origin + ellipse.frame.x.as_vec() * ellipse.semi_major;
        let vtx_id = emit_vertex_point(ctx, vp_world)?;
        let ec_id  = ctx.next_id();
        writeln!(ctx.out, "#{ec_id} = EDGE_CURVE('',#{vtx_id},#{vtx_id},#{ellipse_id},.T.);")?;
        ctx.edge_cache.insert(key, (ec_id, vtx_id));
        ec_id
    };
    let sense  = if orient { ".T." } else { ".F." };
    let oe_id  = ctx.next_id();
    writeln!(ctx.out, "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
    let el_id  = ctx.next_id();
    writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',(#{oe_id}));")?;
    let fb_id  = ctx.next_id();
    let btype  = if outer { "FACE_OUTER_BOUND" } else { "FACE_BOUND" };
    writeln!(ctx.out, "#{fb_id} = {btype}('',#{el_id},.T.);")?;
    Ok(fb_id)
}

/// Write the entire B-Rep to a STEP string, exporting all solids.
pub fn brep_to_step(brep: &BRep) -> Result<String, StepError> {
    StepWriter::new(brep).to_step(&[])
}

fn normal_key(v: UnitVec3) -> [i64; 3] {
    let vec = v.as_vec();
    let mut x = (vec.x * 10.0).round() as i64;
    let mut y = (vec.y * 10.0).round() as i64;
    let mut z = (vec.z * 10.0).round() as i64;
    if x != 0 {
        if x < 0 { x = -x; y = -y; z = -z; }
    } else if y != 0 {
        if y < 0 { y = -y; z = -z; }
    } else if z != 0 {
        if z < 0 { z = -z; }
    }
    [x, y, z]
}
