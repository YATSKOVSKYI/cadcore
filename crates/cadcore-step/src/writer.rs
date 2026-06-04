//! High-level STEP writer — converts a [`BRep`] to an ISO 10303-21 string.
//!
//! Extended with `PartialCylinder` and `PartialDisk` for solid half-space cut output.

use std::fmt::Write as FmtWrite;

use cadcore_geom::Ellipse3;
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{BRep, FaceBoundary, FaceExtent, FaceGeom, FaceNormal, SolidId};

use crate::entities::{
    arc_edge_key, dir_key, emit_cylinder, emit_ellipse, emit_plane, emit_point, emit_sphere,
    emit_torus, emit_unit_direction, emit_vertex_point, point_key, Ctx, StepCurveKey, StepError,
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
        writeln!(
            ctx.out,
            "#{prod_ctx_id} = PRODUCT_CONTEXT('',#{app_ctx_id},'mechanical');"
        )?;
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

        let solids_to_export: Vec<SolidId> = if solid_ids.is_empty() {
            self.brep.solids.keys().collect()
        } else {
            solid_ids.to_vec()
        };

        let mut shape_rep_items: Vec<usize> = Vec::new();

        for (i, solid_id) in solids_to_export.iter().enumerate() {
            if let Some(solid) = self.brep.solids.get(*solid_id) {
                let adv_faces = self.emit_advanced_faces(&mut ctx, *solid_id)?;
                let af_refs: String = adv_faces
                    .iter()
                    .map(|id| format!("#{id}"))
                    .collect::<Vec<_>>()
                    .join(",");

                let label = solid
                    .name
                    .as_deref()
                    .filter(|n| !n.is_empty())
                    .map(|n| format!("{n}_{i}"))
                    .unwrap_or_else(|| format!("solid_{i}"));

                let csb_id = ctx.next_id();
                writeln!(ctx.out, "#{csb_id} = CLOSED_SHELL('',({af_refs}));")?;
                let msb_id = ctx.next_id();
                writeln!(
                    ctx.out,
                    "#{msb_id} = MANIFOLD_SOLID_BREP('{label}',#{csb_id});"
                )?;
                shape_rep_items.push(msb_id);
            }
        }

        let geom_ctx_id = ctx.next_id();
        writeln!(ctx.out,
            "#{geom_ctx_id} = (GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{})) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{},#{},#{})) REPRESENTATION_CONTEXT('','3D'));\n",
            ctx.counter + 1, ctx.counter + 2, ctx.counter + 3, ctx.counter + 4,
        )?;
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

        out.push_str(&ctx.out);
        out.push_str("ENDSEC;\n");
        out.push_str("END-ISO-10303-21;\n");
        Ok(out)
    }

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
            // Reset caches at shell boundary only — entities must not bleed across
            // separate solids.  Within the shell the caches persist so adjacent
            // faces share VERTEX_POINT and EDGE_CURVE entities (STEP AP214 rule).
            ctx.point_cache.clear();
            ctx.vertex_cache.clear();
            ctx.edge_cache.clear();

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
                let sense = match face.normal {
                    FaceNormal::Same => ".T.",
                    FaceNormal::Reversed => ".F.",
                };
                let bounds = emit_face_bounds(ctx, face)?;
                let af_id = ctx.next_id();
                if bounds.is_empty() {
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

fn emit_face_bounds(ctx: &mut Ctx, face: &cadcore_topo::Face) -> Result<Vec<usize>, StepError> {
    match &face.extent {
        // ── Full cylinder ─────────────────────────────────────────────────────
        FaceExtent::Cylinder { start, end, .. } => {
            let cyl = match &face.geom {
                FaceGeom::Cylinder(c) => c,
                _ => return Ok(vec![]),
            };
            let outer = emit_boundary(ctx, start, cyl.frame.x, true, true)?;
            let inner = emit_boundary(ctx, end, cyl.frame.x, false, false)?;
            Ok(vec![outer, inner])
        }

        // ── Partial cylinder (chord cut) ──────────────────────────────────────
        FaceExtent::PartialCylinder {
            length,
            arc_start_angle,
            arc_end_angle,
            arc_ref_dir,
        } => {
            let cyl = match &face.geom {
                FaceGeom::Cylinder(c) => c,
                _ => return Ok(vec![]),
            };
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
            let plane = match &face.geom {
                FaceGeom::Plane(p) => p,
                _ => return Ok(vec![]),
            };
            let bound = emit_circle_bound(
                ctx,
                plane.frame.origin,
                plane.frame.z,
                plane.frame.x,
                *radius,
                true,
                true,
            )?;
            Ok(vec![bound])
        }

        // ── Partial disk (arc + chord, end cap of partial cylinder) ───────────
        FaceExtent::PartialDisk {
            radius,
            start_angle,
            end_angle,
        } => {
            let plane = match &face.geom {
                FaceGeom::Plane(p) => p,
                _ => return Ok(vec![]),
            };
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

        // ── Planar boundary (circle or ellipse cap) ───────────────────────────
        FaceExtent::PlanarBoundary { boundary } => {
            let plane = match &face.geom {
                FaceGeom::Plane(p) => p,
                _ => return Ok(vec![]),
            };
            let outer = emit_boundary(ctx, boundary, plane.frame.x, true, true)?;
            Ok(vec![outer])
        }

        // ── Polygonal flat face ───────────────────────────────────────────────
        FaceExtent::Polygon { points } => {
            if points.len() < 3 {
                return Ok(vec![]);
            }
            let mut vtx_ids = Vec::with_capacity(points.len());
            for &pt in points {
                let vtx_id = emit_vertex_point(ctx, pt)?;
                vtx_ids.push(vtx_id);
            }
            let mut oe_ids = Vec::with_capacity(points.len());
            let n = points.len();
            for i in 0..n {
                let p_start = points[i];
                let p_end = points[(i + 1) % n];
                let dir_vec = p_end - p_start;
                if dir_vec.length() < 1e-7 {
                    continue;
                }
                let dir = match UnitVec3::try_from_vec(dir_vec) {
                    Some(u) => u,
                    None => continue,
                };
                let v_s = vtx_ids[i];
                let v_e = vtx_ids[(i + 1) % n];

                let v_min = v_s.min(v_e);
                let v_max = v_s.max(v_e);
                let key = StepCurveKey::Line {
                    v1: v_min,
                    v2: v_max,
                };
                let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
                    pair
                } else {
                    let lp_id = emit_point(ctx, p_start, "lp")?;
                    let ld_id = emit_unit_direction(ctx, dir, "ld")?;
                    let line_id = ctx.next_id();
                    writeln!(ctx.out, "#{line_id} = LINE('',#{lp_id},#{ld_id});")?;
                    let ec_id = ctx.next_id();
                    writeln!(
                        ctx.out,
                        "#{ec_id} = EDGE_CURVE('',#{v_s},#{v_e},#{line_id},.T.);"
                    )?;
                    ctx.edge_cache.insert(key, (ec_id, v_s));
                    (ec_id, v_s)
                };
                let orient = orig_start == v_s;
                let sense = if orient { ".T." } else { ".F." };
                let oe_id = ctx.next_id();
                writeln!(
                    ctx.out,
                    "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});"
                )?;
                oe_ids.push(oe_id);
            }
            if oe_ids.is_empty() {
                return Ok(vec![]);
            }
            let el_id = ctx.next_id();
            let oe_refs = oe_ids
                .iter()
                .map(|id| format!("#{id}"))
                .collect::<Vec<_>>()
                .join(",");
            writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',({oe_refs}));")?;
            let fb_id = ctx.next_id();
            writeln!(ctx.out, "#{fb_id} = FACE_OUTER_BOUND('',#{el_id},.T.);")?;
            Ok(vec![fb_id])
        }

        // ── Rounded-rectangle cap (4 lines + 4 arcs) — analytic bounds ─────
        FaceExtent::RoundedRectCap { xmin, xmax, zmin, zmax, radius, y, plus_y } => {
            emit_rounded_rect_cap_bounds(
                ctx, *xmin, *xmax, *zmin, *zmax, *radius, *y, *plus_y,
            )
        }

        // ── Cylinder arc face (quarter-cylinder fillet) — analytic bounds ────
        FaceExtent::CylinderArcFace {
            length,
            arc_start_angle,
            arc_end_angle,
            arc_ref_dir,
        } => {
            let cyl = match &face.geom {
                FaceGeom::Cylinder(c) => c,
                _ => return Ok(vec![]),
            };
            emit_cylinder_arc_face_bounds(
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

        FaceExtent::None => Ok(vec![]),
    }
}

// ── Partial cylinder bound emission ───────────────────────────────────────────

fn emit_partial_cylinder_bounds(
    ctx: &mut Ctx,
    origin: Point3,
    axis: UnitVec3,
    x_ref: UnitVec3,
    radius: f64,
    length: f64,
    arc_start_angle: f64,
    arc_end_angle: f64,
    arc_ref_dir: UnitVec3,
) -> Result<Vec<usize>, StepError> {
    let right = match UnitVec3::try_from_vec(axis.cross(arc_ref_dir)) {
        Some(u) => u,
        None => x_ref,
    };

    let n_segs = 24usize;
    let angle_range = arc_end_angle - arc_start_angle;
    let mut pts_start: Vec<Point3> = Vec::with_capacity(n_segs + 1);
    let mut pts_end: Vec<Point3> = Vec::with_capacity(n_segs + 1);

    let end = origin + axis.as_vec() * length;

    for i in 0..=n_segs {
        let t = i as f64 / n_segs as f64;
        let angle = arc_start_angle + t * angle_range;
        let local =
            arc_ref_dir.as_vec() * (radius * angle.cos()) + right.as_vec() * (radius * angle.sin());
        pts_start.push(origin + local);
        pts_end.push(end + local);
    }

    let mut all_pts: Vec<Point3> = Vec::new();
    all_pts.extend_from_slice(&pts_start);
    all_pts.extend(pts_end.iter().rev().cloned());

    if (*all_pts.first().unwrap() - *all_pts.last().unwrap()).length() < 1e-7 {
        all_pts.pop();
    }

    if all_pts.len() < 3 {
        return Ok(vec![]);
    }

    let mut vtx_ids = Vec::with_capacity(all_pts.len());
    for &pt in &all_pts {
        let vtx_id = emit_vertex_point(ctx, pt)?;
        vtx_ids.push(vtx_id);
    }

    let mut oe_ids = Vec::new();
    let n = all_pts.len();
    for i in 0..n {
        let p_start = all_pts[i];
        let p_end = all_pts[(i + 1) % n];
        let dv = p_end - p_start;
        if dv.length() < 1e-7 {
            continue;
        }
        let dir = match UnitVec3::try_from_vec(dv) {
            Some(u) => u,
            None => continue,
        };
        let v_s = vtx_ids[i];
        let v_e = vtx_ids[(i + 1) % n];

        let v_min = v_s.min(v_e);
        let v_max = v_s.max(v_e);
        let key = StepCurveKey::Line {
            v1: v_min,
            v2: v_max,
        };
        let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
            pair
        } else {
            let lp_id = emit_point(ctx, p_start, "lp")?;
            let ld_id = emit_unit_direction(ctx, dir, "ld")?;
            let line_id = ctx.next_id();
            writeln!(ctx.out, "#{line_id} = LINE('',#{lp_id},#{ld_id});")?;
            let ec_id = ctx.next_id();
            writeln!(
                ctx.out,
                "#{ec_id} = EDGE_CURVE('',#{v_s},#{v_e},#{line_id},.T.);"
            )?;
            ctx.edge_cache.insert(key, (ec_id, v_s));
            (ec_id, v_s)
        };
        let orient = orig_start == v_s;
        let sense = if orient { ".T." } else { ".F." };
        let oe_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});"
        )?;
        oe_ids.push(oe_id);
    }
    if oe_ids.is_empty() {
        return Ok(vec![]);
    }
    let oe_refs = oe_ids
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let el_id = ctx.next_id();
    writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',({oe_refs}));")?;
    let fb_id = ctx.next_id();
    writeln!(ctx.out, "#{fb_id} = FACE_OUTER_BOUND('',#{el_id},.T.);")?;
    Ok(vec![fb_id])
}

// ── Partial disk bound emission ───────────────────────────────────────────────

fn emit_partial_disk_bound(
    ctx: &mut Ctx,
    centre: Point3,
    normal: UnitVec3,
    x_ref: UnitVec3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
) -> Result<Vec<usize>, StepError> {
    let n_segs = 24usize;
    let angle_range = end_angle - start_angle;

    let y_ref = match UnitVec3::try_from_vec(normal.cross(x_ref)) {
        Some(u) => u,
        None => return Ok(vec![]),
    };

    let mut pts: Vec<Point3> = Vec::with_capacity(n_segs + 2);
    for i in 0..=n_segs {
        let t = i as f64 / n_segs as f64;
        let angle = start_angle + t * angle_range;
        let local =
            x_ref.as_vec() * (radius * angle.cos()) + y_ref.as_vec() * (radius * angle.sin());
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
        let p_end = pts[(i + 1) % n];
        let dv = p_end - p_start;
        if dv.length() < 1e-7 {
            continue;
        }
        let dir = match UnitVec3::try_from_vec(dv) {
            Some(u) => u,
            None => continue,
        };
        let v_s = vtx_ids[i];
        let v_e = vtx_ids[(i + 1) % n];

        let v_min = v_s.min(v_e);
        let v_max = v_s.max(v_e);
        let key = StepCurveKey::Line {
            v1: v_min,
            v2: v_max,
        };
        let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
            pair
        } else {
            let lp_id = emit_point(ctx, p_start, "lp")?;
            let ld_id = emit_unit_direction(ctx, dir, "ld")?;
            let line_id = ctx.next_id();
            writeln!(ctx.out, "#{line_id} = LINE('',#{lp_id},#{ld_id});")?;
            let ec_id = ctx.next_id();
            writeln!(
                ctx.out,
                "#{ec_id} = EDGE_CURVE('',#{v_s},#{v_e},#{line_id},.T.);"
            )?;
            ctx.edge_cache.insert(key, (ec_id, v_s));
            (ec_id, v_s)
        };
        let orient = orig_start == v_s;
        let sense = if orient { ".T." } else { ".F." };
        let oe_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});"
        )?;
        oe_ids.push(oe_id);
    }
    if oe_ids.is_empty() {
        return Ok(vec![]);
    }
    let oe_refs = oe_ids
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
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

fn emit_circle_bound(
    ctx: &mut Ctx,
    centre: Point3,
    normal: UnitVec3,
    _x_dir: UnitVec3,
    radius: f64,
    outer: bool,
    orient: bool,
) -> Result<usize, StepError> {
    let key = StepCurveKey::Circle {
        center: point_key(centre),
        radius_micro: (radius * 1_000_000.0).round() as i64,
        normal: normal_key(normal),
    };
    // Cache second element stores first-user's orient (0=.F., 1=.T.) so that
    // the second user (adjacent face sharing this edge) always gets the
    // opposite sense — required by AP214 manifold B-rep rule.
    let (ec_id, sense) = if let Some(&(id, first_orient)) = ctx.edge_cache.get(&key) {
        let s = if first_orient != 0 { ".F." } else { ".T." };
        (id, s)
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
        let vp_world = centre + x_dir.as_vec() * radius;
        let vtx_id = emit_vertex_point(ctx, vp_world)?;
        let ec_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{ec_id} = EDGE_CURVE('',#{vtx_id},#{vtx_id},#{circ_id},.T.);"
        )?;
        ctx.edge_cache.insert(key, (ec_id, orient as usize));
        let s = if orient { ".T." } else { ".F." };
        (ec_id, s)
    };
    let oe_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});"
    )?;
    let el_id = ctx.next_id();
    writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',(#{oe_id}));")?;
    let fb_id = ctx.next_id();
    let btype = if outer {
        "FACE_OUTER_BOUND"
    } else {
        "FACE_BOUND"
    };
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
    // Same flip logic as emit_circle_bound: second user always uses opposite sense.
    let (ec_id, sense) = if let Some(&(id, first_orient)) = ctx.edge_cache.get(&key) {
        let s = if first_orient != 0 { ".F." } else { ".T." };
        (id, s)
    } else {
        let ellipse_id = emit_ellipse(ctx, ellipse)?;
        let vp_world = ellipse.frame.origin + ellipse.frame.x.as_vec() * ellipse.semi_major;
        let vtx_id = emit_vertex_point(ctx, vp_world)?;
        let ec_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{ec_id} = EDGE_CURVE('',#{vtx_id},#{vtx_id},#{ellipse_id},.T.);"
        )?;
        ctx.edge_cache.insert(key, (ec_id, orient as usize));
        let s = if orient { ".T." } else { ".F." };
        (ec_id, s)
    };
    let oe_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{oe_id} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});"
    )?;
    let el_id = ctx.next_id();
    writeln!(ctx.out, "#{el_id} = EDGE_LOOP('',(#{oe_id}));")?;
    let fb_id = ctx.next_id();
    let btype = if outer {
        "FACE_OUTER_BOUND"
    } else {
        "FACE_BOUND"
    };
    writeln!(ctx.out, "#{fb_id} = {btype}('',#{el_id},.T.);")?;
    Ok(fb_id)
}

/// Write the entire B-Rep to a STEP string, exporting all solids.
pub fn brep_to_step(brep: &BRep) -> Result<String, StepError> {
    StepWriter::new(brep).to_step(&[])
}

/// Write the B-Rep as a **STEP AP214 Assembly**: one root product containing
/// every solid as an independent component, each with its own full analytic
/// `MANIFOLD_SOLID_BREP`.  No Boolean union is required — FEM importers
/// (COMSOL, ANSYS, Abaqus, CalculiX) read the assembly and mesh shared
/// contact interfaces automatically.
///
/// Each component is placed at the identity transform (no translation /
/// rotation), which is correct because cadcore sweeps every solid directly
/// into world coordinates.
pub fn brep_to_step_assembly(brep: &BRep) -> Result<String, StepError> {
    let mut ctx = Ctx::new();
    let mut out = String::with_capacity(256 * 1024);

    // ── STEP header ───────────────────────────────────────────────────────────
    out.push_str("ISO-10303-21;\n");
    out.push_str("HEADER;\n");
    out.push_str("  FILE_DESCRIPTION(('cadcore B-Rep assembly export'),'2;1');\n");
    out.push_str("  FILE_NAME('','',(''),(''),'cadcore 0.1','','');\n");
    out.push_str("  FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
    out.push_str("ENDSEC;\n");
    out.push_str("DATA;\n");

    // ── Shared context entities ───────────────────────────────────────────────
    let app_ctx_id = ctx.next_id();
    writeln!(ctx.out, "#{app_ctx_id} = APPLICATION_CONTEXT('core data for automotive mechanical design processes');")?;
    let apd_id = ctx.next_id();
    writeln!(ctx.out, "#{apd_id} = APPLICATION_PROTOCOL_DEFINITION('draft international standard','automotive_design',1998,#{app_ctx_id});")?;
    let prod_ctx_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{prod_ctx_id} = PRODUCT_CONTEXT('',#{app_ctx_id},'mechanical');"
    )?;

    // ── Geometric context (mm, radians) ──────────────────────────────────────
    // We need the ids before writing the entities because GEOMETRIC_REPRESENTATION_CONTEXT
    // references the uncertainty measure which follows it.
    let geom_ctx_id = ctx.next_id();
    let unc_id = ctx.next_id();
    let lu_id = ctx.next_id();
    let au_id = ctx.next_id();
    let su_id = ctx.next_id();
    writeln!(ctx.out,
        "#{geom_ctx_id} = (GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{unc_id})) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{lu_id},#{au_id},#{su_id})) REPRESENTATION_CONTEXT('','3D'));")?;
    writeln!(
        ctx.out,
        "#{unc_id} = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-007),#{lu_id},('',''));"
    )?;
    writeln!(
        ctx.out,
        "#{lu_id} = (LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));"
    )?;
    writeln!(
        ctx.out,
        "#{au_id} = (NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));"
    )?;
    writeln!(
        ctx.out,
        "#{su_id} = (NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT());"
    )?;

    // Shared identity geometry — used as a template for per-component axis placements.
    // We write TWO shared direction/point entities that are safe to reuse (they
    // carry no instance-specific data), but each component gets its OWN two
    // AXIS2_PLACEMENT_3D instances inside ITEM_DEFINED_TRANSFORMATION so that
    // every NAUO pair references distinct entities, which is required by
    // STEP AP214 §4.4.3 and expected by OpenCASCADE / FreeCAD.
    let origin_pt_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{origin_pt_id} = CARTESIAN_POINT('',(0.0,0.0,0.0));"
    )?;
    let z_dir_id = ctx.next_id();
    writeln!(ctx.out, "#{z_dir_id} = DIRECTION('',(0.0,0.0,1.0));")?;
    let x_dir_id = ctx.next_id();
    writeln!(ctx.out, "#{x_dir_id} = DIRECTION('',(1.0,0.0,0.0));")?;

    // Assembly's own axis placement (used in root SHAPE_REPRESENTATION).
    let asm_ax_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{asm_ax_id} = AXIS2_PLACEMENT_3D('',#{origin_pt_id},#{z_dir_id},#{x_dir_id});"
    )?;

    // ── Root assembly product ─────────────────────────────────────────────────
    let asm_prod_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{asm_prod_id} = PRODUCT('scaffold_assembly','scaffold_assembly','',(#{prod_ctx_id}));"
    )?;
    let asm_pdf_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{asm_pdf_id} = PRODUCT_DEFINITION_FORMATION('','',#{asm_prod_id});"
    )?;
    let asm_pd_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{asm_pd_id} = PRODUCT_DEFINITION('design','',#{asm_pdf_id},#{prod_ctx_id});"
    )?;
    let asm_pds_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{asm_pds_id} = PRODUCT_DEFINITION_SHAPE('','',#{asm_pd_id});"
    )?;
    let asm_sr_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{asm_sr_id} = SHAPE_REPRESENTATION('scaffold_assembly',(#{asm_ax_id}),#{geom_ctx_id});"
    )?;
    let asm_sdr_id = ctx.next_id();
    writeln!(
        ctx.out,
        "#{asm_sdr_id} = SHAPE_DEFINITION_REPRESENTATION(#{asm_pds_id},#{asm_sr_id});"
    )?;

    // ── One component per solid ───────────────────────────────────────────────
    let solid_ids: Vec<SolidId> = brep.solids.keys().collect();

    for (i, &solid_id) in solid_ids.iter().enumerate() {
        let solid = match brep.solids.get(solid_id) {
            Some(s) => s,
            None => continue,
        };

        // Build a unique component label: prefer the solid name, but always
        // append the index so names are guaranteed unique in the STEP file
        // (PRODUCT 'id' must be unique per AP214 §4.4.1, and FreeCAD uses the
        // first field as the tree label).
        let base = solid
            .name
            .as_deref()
            .filter(|n| !n.is_empty())
            .unwrap_or("part");
        let label = format!("{base}_{i}");

        // B-Rep geometry (identical to single-solid export).
        let writer = StepWriter::new(brep);
        let adv_faces = writer.emit_advanced_faces(&mut ctx, solid_id)?;
        if adv_faces.is_empty() {
            continue;
        }
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
            "#{msb_id} = MANIFOLD_SOLID_BREP('{label}',#{csb_id});"
        )?;

        // Component product hierarchy.
        let comp_prod_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{comp_prod_id} = PRODUCT('{label}','{label}','',(#{prod_ctx_id}));"
        )?;
        let comp_pdf_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{comp_pdf_id} = PRODUCT_DEFINITION_FORMATION('','',#{comp_prod_id});"
        )?;
        let comp_pd_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{comp_pd_id} = PRODUCT_DEFINITION('design','',#{comp_pdf_id},#{prod_ctx_id});"
        )?;
        let comp_pds_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{comp_pds_id} = PRODUCT_DEFINITION_SHAPE('','',#{comp_pd_id});"
        )?;
        let comp_sr_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{comp_sr_id} = SHAPE_REPRESENTATION('{label}',(#{msb_id}),#{geom_ctx_id});"
        )?;
        let comp_sdr_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{comp_sdr_id} = SHAPE_DEFINITION_REPRESENTATION(#{comp_pds_id},#{comp_sr_id});"
        )?;

        // NAUO: assembly → component.
        let nauo_id = ctx.next_id();
        writeln!(ctx.out,
            "#{nauo_id} = NEXT_ASSEMBLY_USAGE_OCCURRENCE('{i}','{label}','',#{asm_pd_id},#{comp_pd_id},$);")?;

        // Identity placement: TWO separate (but equal) AXIS2_PLACEMENT_3D
        // instances are required by AP214 — one for the "source" frame, one
        // for the "target" frame.  Reusing the same entity id for both is
        // technically valid per the standard but rejected by OpenCASCADE.
        let ax_src_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{ax_src_id} = AXIS2_PLACEMENT_3D('',#{origin_pt_id},#{z_dir_id},#{x_dir_id});"
        )?;
        let ax_tgt_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{ax_tgt_id} = AXIS2_PLACEMENT_3D('',#{origin_pt_id},#{z_dir_id},#{x_dir_id});"
        )?;
        let xform_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{xform_id} = ITEM_DEFINED_TRANSFORMATION('','',#{ax_src_id},#{ax_tgt_id});"
        )?;
        let rr_id = ctx.next_id();
        writeln!(ctx.out,
            "#{rr_id} = (REPRESENTATION_RELATIONSHIP('','',#{comp_sr_id},#{asm_sr_id}) REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#{xform_id}) SHAPE_REPRESENTATION_RELATIONSHIP());")?;
        let cdsr_id = ctx.next_id();
        writeln!(
            ctx.out,
            "#{cdsr_id} = CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#{rr_id},#{nauo_id});"
        )?;

        let _ = (comp_sdr_id, asm_sdr_id, cdsr_id);
    }

    out.push_str(&ctx.out);
    out.push_str("ENDSEC;\n");
    out.push_str("END-ISO-10303-21;\n");
    Ok(out)
}

// ── Shared-edge helpers for rounded-box prism faces ──────────────────────────
//
// The rounded box (electrode) is a prism: a rounded-rectangle XZ profile
// extruded along Y.  For the shell to be a watertight AP214 manifold, EVERY
// edge must be shared by exactly two faces traversed in OPPOSITE directions.
//
// These helpers funnel every straight and arc edge through `ctx.edge_cache`
// (keyed by geometry, direction-agnostic).  The first face to touch an edge
// emits the EDGE_CURVE and records its start vertex; the second face reuses
// the same EDGE_CURVE and gets the opposite ORIENTED_EDGE sense automatically.
//
// Because every face is wound CCW-as-seen-from-outside, the two faces meeting
// at any edge necessarily traverse it in opposite vertex order — so the cache
// sense logic yields exactly one `.T.` and one `.F.` per edge.  This is the
// invariant the `check_ap214_manifold` test enforces.

/// Append one ORIENTED_EDGE for a straight edge `p_from → p_to`, reusing a
/// shared LINE/EDGE_CURVE from the cache when another face already emitted it.
fn push_shared_line(
    ctx: &mut Ctx,
    p_from: Point3,
    p_to: Point3,
    oe_ids: &mut Vec<usize>,
) -> Result<(), StepError> {
    let dv = p_to - p_from;
    if dv.length() < 1e-9 {
        return Ok(());
    }
    let v_from = emit_vertex_point(ctx, p_from)?;
    let v_to = emit_vertex_point(ctx, p_to)?;
    let key = StepCurveKey::Line {
        v1: v_from.min(v_to),
        v2: v_from.max(v_to),
    };
    let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
        pair
    } else {
        let dir = match UnitVec3::try_from_vec(dv) {
            Some(u) => u,
            None => return Ok(()),
        };
        let lp = emit_point(ctx, p_from, "rbl_p")?;
        let ld = emit_unit_direction(ctx, dir, "rbl_d")?;
        let li = ctx.next_id();
        writeln!(ctx.out, "#{li} = LINE('',#{lp},#{ld});")?;
        let ec = ctx.next_id();
        writeln!(ctx.out, "#{ec} = EDGE_CURVE('',#{v_from},#{v_to},#{li},.T.);")?;
        ctx.edge_cache.insert(key, (ec, v_from));
        (ec, v_from)
    };
    let sense = if orig_start == v_from { ".T." } else { ".F." };
    let oe = ctx.next_id();
    writeln!(ctx.out, "#{oe} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
    oe_ids.push(oe);
    Ok(())
}

/// Append one ORIENTED_EDGE for a 90° arc edge `p_from → p_to` lying on the
/// circle (`centre`, normal `axis`, x-reference `xref`, radius `r`), reusing a
/// shared CIRCLE/EDGE_CURVE from the cache when another face already emitted it.
///
/// The EDGE_CURVE `same_sense` flag is derived from geometry: it is `.T.` when
/// `p_from → p_to` runs counter-clockwise about `axis` (the circle's natural
/// parametrisation direction) and `.F.` otherwise.  This makes the helper
/// agnostic to which corner / which Y-level it is called for.
fn push_shared_arc(
    ctx: &mut Ctx,
    centre: Point3,
    axis: UnitVec3,
    xref: UnitVec3,
    r: f64,
    p_from: Point3,
    p_to: Point3,
    oe_ids: &mut Vec<usize>,
) -> Result<(), StepError> {
    let v_from = emit_vertex_point(ctx, p_from)?;
    let v_to = emit_vertex_point(ctx, p_to)?;

    let r_micro = (r * 1_000_000.0).round() as i64;
    let key = arc_edge_key(
        point_key(p_from),
        point_key(p_to),
        point_key(centre),
        r_micro,
        dir_key(axis),
        dir_key(xref),
    );

    let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
        pair
    } else {
        // Determine the CCW direction about `axis` from (xref, axis × xref).
        let xref_v = xref.as_vec();
        let yref_v = axis.as_vec().cross(xref_v);
        let dot = |a: cadcore_math::Vec3, b: cadcore_math::Vec3| a.x * b.x + a.y * b.y + a.z * b.z;
        let ang = |p: Point3| -> f64 {
            let d = p - centre;
            dot(d, yref_v).atan2(dot(d, xref_v))
        };
        let mut dtheta = ang(p_to) - ang(p_from);
        let two_pi = std::f64::consts::PI * 2.0;
        while dtheta <= -std::f64::consts::PI {
            dtheta += two_pi;
        }
        while dtheta > std::f64::consts::PI {
            dtheta -= two_pi;
        }
        let same_sense = if dtheta >= 0.0 { ".T." } else { ".F." };

        let cp = emit_point(ctx, centre, "rba_c")?;
        let cn = emit_unit_direction(ctx, axis, "rba_n")?;
        let cx = emit_unit_direction(ctx, xref, "rba_x")?;
        let ax = ctx.next_id();
        writeln!(ctx.out, "#{ax} = AXIS2_PLACEMENT_3D('',#{cp},#{cn},#{cx});")?;
        let ci = ctx.next_id();
        writeln!(ctx.out, "#{ci} = CIRCLE('',#{ax},{:.10});", r)?;
        let ec = ctx.next_id();
        writeln!(ctx.out, "#{ec} = EDGE_CURVE('',#{v_from},#{v_to},#{ci},{same_sense});")?;
        ctx.edge_cache.insert(key, (ec, v_from));
        (ec, v_from)
    };
    let sense = if orig_start == v_from { ".T." } else { ".F." };
    let oe = ctx.next_id();
    writeln!(ctx.out, "#{oe} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
    oe_ids.push(oe);
    Ok(())
}

/// Emit an EDGE_LOOP + FACE_OUTER_BOUND from a list of ORIENTED_EDGE ids.
fn finish_outer_bound(ctx: &mut Ctx, oe_ids: &[usize]) -> Result<Vec<usize>, StepError> {
    if oe_ids.is_empty() {
        return Ok(vec![]);
    }
    let oe_refs = oe_ids
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let el = ctx.next_id();
    writeln!(ctx.out, "#{el} = EDGE_LOOP('',({oe_refs}));")?;
    let fb = ctx.next_id();
    writeln!(ctx.out, "#{fb} = FACE_OUTER_BOUND('',#{el},.T.);")?;
    Ok(vec![fb])
}

// ── Cylinder arc face bounds (quarter-cylinder corner of the rounded box) ─────
//
// One quarter-cylinder side face of the prism.  Its boundary is a quad:
//   B_i → T_i  (longitudinal line, ymin → ymax, at the arc start vertex)
//   T_i → T_e  (top arc,  at y = ymax)
//   T_e → B_e  (longitudinal line, ymax → ymin, at the arc end vertex)
//   B_e → B_i  (bottom arc, at y = ymin, traversed backward)
//
// This winding is CCW as seen from outside (radially outward normal), so the
// shared arcs pair oppositely with the caps and the longitudinals pair
// oppositely with the adjacent flat side faces.
fn emit_cylinder_arc_face_bounds(
    ctx: &mut Ctx,
    origin: Point3,
    axis: UnitVec3,
    _x_ref: UnitVec3,
    radius: f64,
    length: f64,
    arc_start_angle: f64,
    arc_end_angle: f64,
    arc_ref_dir: UnitVec3,
) -> Result<Vec<usize>, StepError> {
    let right = {
        let v = axis.as_vec().cross(arc_ref_dir.as_vec());
        match UnitVec3::try_from_vec(v) {
            Some(u) => u,
            None => return Ok(vec![]),
        }
    };

    let arc_pt = |angle: f64, z_off: f64| -> Point3 {
        let local = arc_ref_dir.as_vec() * (radius * angle.cos())
            + right.as_vec() * (radius * angle.sin());
        origin + local + axis.as_vec() * z_off
    };

    let p_s0 = arc_pt(arc_start_angle, 0.0); // B_i  (arc start, ymin)
    let p_e0 = arc_pt(arc_end_angle, 0.0); // B_e  (arc end,   ymin)
    let p_s1 = arc_pt(arc_start_angle, length); // T_i  (arc start, ymax)
    let p_e1 = arc_pt(arc_end_angle, length); // T_e  (arc end,   ymax)

    let centre_lo = origin;
    let centre_hi = origin + axis.as_vec() * length;

    let mut oe_ids: Vec<usize> = Vec::with_capacity(4);
    push_shared_line(ctx, p_s0, p_s1, &mut oe_ids)?; // B_i → T_i  (up)
    push_shared_arc(ctx, centre_hi, axis, arc_ref_dir, radius, p_s1, p_e1, &mut oe_ids)?; // top arc
    push_shared_line(ctx, p_e1, p_e0, &mut oe_ids)?; // T_e → B_e  (down)
    push_shared_arc(ctx, centre_lo, axis, arc_ref_dir, radius, p_e0, p_s0, &mut oe_ids)?; // bottom arc (backward)

    finish_outer_bound(ctx, &oe_ids)
}

// ── Rounded-rectangle cap bounds (4 LINE + 4 CIRCLE arc edges) ───────────────
//
// One end cap of the prism (a flat rounded rectangle at y = `y`).  Traversed
// along the profile so the loop is CCW as seen from outside:
//   * front cap (−Y, `plus_y = false`): profile forward  (segments 0..8)
//   * back  cap (+Y, `plus_y = true` ): profile backward  (segments 7..0)
//
// Every line/arc here is shared with the matching prism side face, so the cap
// is sewn into the shell instead of floating as a loose surface.
fn emit_rounded_rect_cap_bounds(
    ctx: &mut Ctx,
    xmin: f64, xmax: f64,
    zmin: f64, zmax: f64,
    r: f64,
    y: f64,
    plus_y: bool,
) -> Result<Vec<usize>, StepError> {
    use cadcore_math::Vec3;

    // Profile junction points (CCW in XZ as seen from −Y).
    let pts = [
        Point3::new(xmin + r, y, zmin),     // 0
        Point3::new(xmax - r, y, zmin),     // 1
        Point3::new(xmax,     y, zmin + r), // 2
        Point3::new(xmax,     y, zmax - r), // 3
        Point3::new(xmax - r, y, zmax),     // 4
        Point3::new(xmin + r, y, zmax),     // 5
        Point3::new(xmin,     y, zmax - r), // 6
        Point3::new(xmin,     y, zmin + r), // 7
    ];

    // Per-segment geometry: (is_arc, centre, xref).  Centres/xrefs match the
    // quarter-cylinder corner faces so the arc-edge cache keys collide.
    let seg_arc: [bool; 8] = [false, true, false, true, false, true, false, true];
    let centres = [
        Point3::new(xmax - r, y, zmin + r), // seg1 BR
        Point3::new(xmax - r, y, zmax - r), // seg3 TR
        Point3::new(xmin + r, y, zmax - r), // seg5 TL
        Point3::new(xmin + r, y, zmin + r), // seg7 BL
    ];
    let xrefs: [Vec3; 4] = [
        Vec3::new(0.0, 0.0, -1.0), // BR: −Z
        Vec3::new(1.0, 0.0, 0.0),  // TR: +X
        Vec3::new(0.0, 0.0, 1.0),  // TL: +Z
        Vec3::new(-1.0, 0.0, 0.0), // BL: −X
    ];

    let axis_y = UnitVec3::try_from_vec(Vec3::new(0.0, 1.0, 0.0)).unwrap();

    // Map a segment index (the arc index 0..4) from the profile segment id.
    let arc_corner = |seg: usize| -> usize {
        match seg {
            1 => 0, // BR
            3 => 1, // TR
            5 => 2, // TL
            7 => 3, // BL
            _ => 0,
        }
    };

    let mut oe_ids: Vec<usize> = Vec::with_capacity(8);

    // Order + direction of the 8 segments depends on the cap's outward normal.
    let order: Vec<(Point3, Point3, usize)> = if plus_y {
        // +Y cap → profile backward: segment i traversed pts[i+1] → pts[i].
        (0..8)
            .rev()
            .map(|i| (pts[(i + 1) % 8], pts[i], i))
            .collect()
    } else {
        // −Y cap → profile forward: segment i traversed pts[i] → pts[i+1].
        (0..8).map(|i| (pts[i], pts[(i + 1) % 8], i)).collect()
    };

    for (p_from, p_to, seg) in order {
        if seg_arc[seg] {
            let c = arc_corner(seg);
            let xref = match UnitVec3::try_from_vec(xrefs[c]) {
                Some(u) => u,
                None => continue,
            };
            push_shared_arc(ctx, centres[c], axis_y, xref, r, p_from, p_to, &mut oe_ids)?;
        } else {
            push_shared_line(ctx, p_from, p_to, &mut oe_ids)?;
        }
    }

    finish_outer_bound(ctx, &oe_ids)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use cadcore_math::Point3;
    use cadcore_ops::{build_solid_rounded_box_xz, sweep_circle_along_polyline, SweepOptions};
    use cadcore_topo::BRep;
    use super::brep_to_step;

    /// Count occurrences of a STEP keyword in the output.
    fn count(step: &str, keyword: &str) -> usize {
        step.matches(keyword).count()
    }

    /// Collect all EDGE_CURVE ids and the ORIENTED_EDGE ids that reference them.
    /// Returns a map: ec_id → list of (oe_id, sense).
    fn oriented_edge_uses(step: &str) -> std::collections::HashMap<usize, Vec<(usize, &str)>> {
        let mut map: std::collections::HashMap<usize, Vec<(usize, &str)>> =
            std::collections::HashMap::new();
        // Match: #N = ORIENTED_EDGE('',*,*,#M,.T.); or .F.
        let re_pat = "ORIENTED_EDGE";
        for line in step.lines() {
            if !line.contains(re_pat) { continue; }
            // parse #OE = ORIENTED_EDGE('',*,*,#EC,.S.);
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() < 4 { continue; }
            let ec_part = parts[parts.len() - 2].trim();
            let sense_part = parts[parts.len() - 1].trim().trim_end_matches(';');
            let ec_id: usize = ec_part.trim_start_matches('#').parse().ok().unwrap_or(0);
            // OE id from start of line
            let oe_id_str = line.trim().split('=').next().unwrap_or("").trim().trim_start_matches('#');
            let oe_id: usize = oe_id_str.parse().ok().unwrap_or(0);
            if ec_id > 0 && oe_id > 0 {
                map.entry(ec_id).or_default().push((oe_id, sense_part));
            }
        }
        map
    }

    /// For every EDGE_CURVE in the shell, each one should be referenced by
    /// exactly 2 ORIENTED_EDGEs with OPPOSITE senses (.T. and .F.) —
    /// the AP214 manifold closed-shell invariant.
    fn check_ap214_manifold(step: &str) -> Result<(), String> {
        let uses = oriented_edge_uses(step);
        // Collect only EC ids that appear in the DATA section as EDGE_CURVE
        let mut violations = Vec::new();
        for (ec_id, refs) in &uses {
            // Only check if this id is actually an EDGE_CURVE line
            let ec_line = format!("#{ec_id} = EDGE_CURVE");
            if !step.contains(&ec_line) { continue; }
            if refs.len() != 2 {
                violations.push(format!("EC#{ec_id} used {} times (want 2): {:?}", refs.len(), refs));
                continue;
            }
            let s0 = refs[0].1;
            let s1 = refs[1].1;
            let t_count = refs.iter().filter(|(_, s)| s.contains(".T.")).count();
            let f_count = refs.iter().filter(|(_, s)| s.contains(".F.")).count();
            if t_count != 1 || f_count != 1 {
                violations.push(format!(
                    "EC#{ec_id} has wrong senses: {s0} and {s1} (want one .T. one .F.)"
                ));
            }
        }
        if violations.is_empty() {
            Ok(())
        } else {
            Err(violations.join("\n"))
        }
    }

    // ── Straight cylinder + caps ──────────────────────────────────────────────

    /// Single straight cylinder. The cap circles must be shared between the cap
    /// face and the cylinder face, traversed in opposite directions.
    #[test]
    fn straight_cylinder_caps_share_edge_curve() {
        let mut brep = BRep::new();
        sweep_circle_along_polyline(
            &mut brep,
            &[Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 5.0, 0.0)],
            0.3,
            &SweepOptions::default(),
        ).unwrap();
        let step = brep_to_step(&brep).unwrap();
        // 1 cylinder + 2 caps = 3 ADVANCED_FACEs
        assert_eq!(count(&step, "CYLINDRICAL_SURFACE"), 1);
        assert_eq!(count(&step, "ADVANCED_FACE"), 3);
        // Each cap rim circle should appear as a single EDGE_CURVE, referenced
        // by two ORIENTED_EDGEs (cap .T. + cylinder .F., or vice-versa).
        if let Err(e) = check_ap214_manifold(&step) {
            panic!("AP214 manifold violation:\n{e}");
        }
    }

    // ── Rounded box — electrode geometry ─────────────────────────────────────

    /// `build_solid_rounded_box_xz` must export as a closed manifold:
    ///   * 10 ADVANCED_FACEs (2 caps + 4 sides + 4 quarter-cylinders)
    ///   * Every shared arc EDGE_CURVE used by exactly 1 .T. + 1 .F.
    #[test]
    fn rounded_box_xz_is_manifold() {
        let mut brep = BRep::new();
        build_solid_rounded_box_xz(
            &mut brep,
            0.0, 10.0,   // xmin xmax
            0.0,  1.0,   // ymin ymax
            0.0,  8.0,   // zmin zmax
            1.0,         // corner_radius
            Some("electrode".into()),
        ).unwrap();
        let step = brep_to_step(&brep).unwrap();

        assert_eq!(count(&step, "ADVANCED_FACE"), 10,
            "rounded box must produce exactly 10 faces");
        assert_eq!(count(&step, "CYLINDRICAL_SURFACE"), 4,
            "four quarter-cylinder corners");

        if let Err(e) = check_ap214_manifold(&step) {
            panic!("Electrode AP214 manifold violation:\n{e}");
        }
    }

    // ── Half-space cut (trim) — filament end caps ─────────────────────────────

    /// A serpentine filament swept as a SINGLE solid (ellipse miters at the
    /// junctions) and then trimmed by ONE half-space plane must stay a watertight
    /// AP214 manifold.  Regression for the "filament ends not drawn" bug: an
    /// axially-truncated leg used to drop its original miter boundary and emit a
    /// plain circle, while the surviving connector kept its miter ellipse — the
    /// shared junction edge then appeared only once → open shell → CAD tools
    /// rendered the end as broken / cut short on exactly one side.
    #[test]
    fn half_space_cut_serpentine_stays_manifold() {
        use cadcore_math::UnitVec3;
        use cadcore_ops::{
            half_space_cut_brep, sweep_circle_along_path_with_caps, ClipPlane, SweepPathSegment,
        };

        let raw = [
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 20.0, 0.0),
            Point3::new(2.0, 20.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 20.0, 0.0),
        ];
        let segs: Vec<SweepPathSegment> = raw
            .windows(2)
            .map(|w| SweepPathSegment::Line { start: w[0], end: w[1] })
            .collect();

        // A single plane (keeps only one side) leaves connectors alive next to
        // truncated legs — the exact configuration that used to open the shell.
        for &(oy, ny) in &[(1.0_f64, 1.0_f64), (19.0, -1.0)] {
            let mut brep = BRep::new();
            sweep_circle_along_path_with_caps(
                &mut brep,
                &segs,
                0.2,
                &SweepOptions {
                    fillet_corners: false,
                    corner_fillet_radius: 0.0,
                    name: Some("serp".into()),
                },
                None,
                None,
            )
            .unwrap();
            let plane = ClipPlane {
                origin: Point3::new(0.0, oy, 0.0),
                normal: if ny > 0.0 { UnitVec3::Y } else { -UnitVec3::Y },
            };
            half_space_cut_brep(&mut brep, &plane);
            let step = brep_to_step(&brep).unwrap();
            if let Err(e) = check_ap214_manifold(&step) {
                panic!("half-space cut opened the shell (oy={oy}):\n{e}");
            }
        }
    }

    /// The 4 arc EDGE_CURVEs shared between caps and corner cylinders must each
    /// appear exactly once (no duplication).
    #[test]
    fn rounded_box_xz_arc_edges_not_duplicated() {
        let mut brep = BRep::new();
        build_solid_rounded_box_xz(
            &mut brep,
            0.0, 10.0,
            0.0,  2.0,
            0.0,  8.0,
            1.5,
            None,
        ).unwrap();
        let step = brep_to_step(&brep).unwrap();
        let uses = oriented_edge_uses(&step);
        // Every arc EDGE_CURVE (referenced from a CIRCLE) that is shared must
        // have exactly 2 ORIENTED_EDGE references.
        let mut dup_arcs = 0usize;
        for (ec_id, refs) in &uses {
            let ec_line = format!("#{ec_id} = EDGE_CURVE");
            if !step.contains(&ec_line) { continue; }
            if refs.len() > 2 {
                dup_arcs += 1;
            }
        }
        assert_eq!(dup_arcs, 0, "some arc EDGE_CURVEs are referenced more than twice");
    }
}

fn normal_key(v: UnitVec3) -> [i64; 3] {
    let vec = v.as_vec();
    let mut x = (vec.x * 10.0).round() as i64;
    let mut y = (vec.y * 10.0).round() as i64;
    let mut z = (vec.z * 10.0).round() as i64;
    if x != 0 {
        if x < 0 {
            x = -x;
            y = -y;
            z = -z;
        }
    } else if y != 0 {
        if y < 0 {
            y = -y;
            z = -z;
        }
    } else if z != 0 {
        if z < 0 {
            z = -z;
        }
    }
    [x, y, z]
}
