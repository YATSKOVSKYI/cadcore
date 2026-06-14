//! High-level STEP writer — converts a [`BRep`] to an ISO 10303-21 string.
//!
//! Extended with `PartialCylinder` and `PartialDisk` for solid half-space cut output.

use std::collections::HashMap;
use std::fmt::Write as FmtWrite;

use cadcore_geom::Ellipse3;
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{
    BRep, CoEdgeSense, EdgeGeom, EdgeId, FaceBoundary, FaceExtent, FaceGeom, FaceId, FaceNormal,
    LoopId, SolidId,
};

use crate::entities::{
    arc_edge_key, dir_key, emit_cylinder, emit_ellipse, emit_pcurve, emit_plane, emit_point,
    emit_polyline2d, emit_polyline3d, emit_rational_bspline_surface, emit_sphere,
    emit_surface_curve, emit_torus,
    emit_unit_direction, emit_vertex_point, point_key, Ctx, StepCurveKey, StepError,
};
use crate::nurbs::torus_patch_nurbs;
use crate::pcurve::cyl_pcurve;

/// Per-shell context for pcurve emission: each face's emitted surface id, and
/// the faces adjacent to each (topologically-shared) edge.
struct PcurveCtx {
    face_surf: HashMap<FaceId, usize>,
    edge_faces: HashMap<EdgeId, Vec<FaceId>>,
    /// Branch-continuous `(u,v)` chains per `(face, polyline edge)`, computed
    /// by walking every wire with a running cursor: each edge's `u` branch is
    /// chained to the previous edge's end (±2π).  Independently-anchored
    /// per-edge pcurves BREAK a wire that crosses the θ-seam between edges —
    /// importers then see the hole loop open in parameter space and fragment
    /// the face.  Stored in the wire's traversal order.
    wire_uv: HashMap<(FaceId, EdgeId), Vec<(f64, f64)>>,
}

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
            FaceGeom::Torus(t) => {
                // A trimmed-to-an-arc TOROIDAL_SURFACE bounded only by its two
                // minor end-circles is ambiguous (90°/270° band).  Emit the
                // exact rational-NURBS patch spanning only the real arc instead,
                // which removes the ambiguity.
                //
                // EXCEPTION — a torus with WINDOW holes (a crosser cut it): the
                // explicit loop topology already removes the ambiguity, and the
                // NURBS patch's numerically-inverted window pcurves come out
                // "Unorientable" in OCC/SolidWorks.  Emit the exact analytic
                // TOROIDAL_SURFACE instead; its native (θ,φ) pcurves are exact.
                if torus_is_toroidal(face) {
                    emit_torus(ctx, t)
                } else if let Some((lo, hi)) = self.torus_arc_range(face) {
                    let patch = torus_patch_nurbs(t, lo, hi);
                    emit_rational_bspline_surface(ctx, &patch)
                } else {
                    emit_torus(ctx, t)
                }
            }
        }
    }

    /// Recover the major-angle arc range `(theta_lo, theta_hi)` of a torus-elbow
    /// face.  Prefer the two endpoint loops over a global boundary scan: trim
    /// notches and bite holes can add vertices inside the elbow band, while the
    /// endpoint loops are the only topology that defines the NURBS patch span.
    /// Returned range follows the short positive sweep so the patch normal stays
    /// consistent with `FaceNormal::Same`.
    fn torus_arc_range(&self, face: &cadcore_topo::Face) -> Option<(f64, f64)> {
        let t = match &face.geom {
            FaceGeom::Torus(t) => t,
            _ => return None,
        };

        let theta_of = |p: Point3| -> f64 {
            let l = t.frame.to_local_point(p);
            l.y.atan2(l.x)
        };
        let directed_range = |a: f64, b: f64| -> Option<(f64, f64)> {
            let mut d = b - a;
            let pi = std::f64::consts::PI;
            while d > pi {
                d -= 2.0 * pi;
            }
            while d < -pi {
                d += 2.0 * pi;
            }
            if d.abs() < 1e-9 {
                return None;
            }
            Some(if d >= 0.0 { (a, a + d) } else { (a + d, a) })
        };

        // Untrimmed re-homed elbows keep `FaceExtent::TorusFillet`: their geometry
        // lives in the extent (start/end minor circles), and the loop topology is
        // empty — so the boundary scan below would find nothing.  Read the band
        // ends straight from the extent circles.
        if let FaceExtent::TorusFillet {
            start_circle,
            end_circle,
        } = &face.extent
        {
            let a = theta_of(start_circle.frame.origin);
            let b = theta_of(end_circle.frame.origin);
            return directed_range(a, b);
        }

        // Trimmed elbows: explicit topology — recover the band from boundary edges.
        if !matches!(face.extent, FaceExtent::Trimmed) {
            return None;
        }

        let loop_theta = |lp: LoopId| -> Option<f64> {
            let coedges = self.brep.loop_coedges(lp)?;
            let mut exact = Vec::new();
            let mut sampled = Vec::new();
            for ce_id in coedges {
                let ce = self.brep.coedges.get(ce_id)?;
                let edge = self.brep.edges.get(ce.edge)?;
                match &edge.geom {
                    EdgeGeom::Circle(c) if (c.radius - t.minor_radius).abs() < 1e-6 => {
                        exact.push(theta_of(c.frame.origin));
                    }
                    EdgeGeom::Polyline(pts) => {
                        sampled.extend(pts.iter().map(|&p| theta_of(p)));
                    }
                    _ => {
                        if let Some(v) = self.brep.vertices.get(edge.v_start) {
                            sampled.push(theta_of(v.point));
                        }
                        if let Some(v) = self.brep.vertices.get(edge.v_end) {
                            sampled.push(theta_of(v.point));
                        }
                    }
                }
            }
            let src = if exact.len() >= 1 { &exact } else { &sampled };
            if src.is_empty() {
                return None;
            }
            let (sx, sy): (f64, f64) = src
                .iter()
                .fold((0.0, 0.0), |(cx, cy), &th| (cx + th.cos(), cy + th.sin()));
            Some(sy.atan2(sx))
        };

        if let Some(&end_loop) = face.inner_loops.first() {
            if let (Some(a), Some(b)) = (loop_theta(face.outer_loop), loop_theta(end_loop)) {
                if let Some(range) = directed_range(a, b) {
                    return Some(range);
                }
            }
        }

        // Preferred source: the two minor end-circles.  Their centre is the tube
        // centre C(θ) at exactly θ_lo / θ_hi, so it pins the band ends precisely.
        // A notched end boundary is split into several CIRCLE arc-edges that all
        // reuse the SAME minor circle (same `frame.origin`), so this stays exact
        // even after trimming.  Bite holes are polylines, never circles of the
        // minor radius, so they cannot pollute this set.
        let mut circle_thetas: Vec<f64> = Vec::new();
        // Fallback source: every boundary vertex (extremes still bracket the arc).
        let mut vertex_thetas: Vec<f64> = Vec::new();
        for &lp in std::iter::once(&face.outer_loop).chain(face.inner_loops.iter()) {
            let Some(coedges) = self.brep.loop_coedges(lp) else {
                continue;
            };
            for ce_id in coedges {
                let Some(ce) = self.brep.coedges.get(ce_id) else {
                    continue;
                };
                let Some(edge) = self.brep.edges.get(ce.edge) else {
                    continue;
                };
                if let EdgeGeom::Circle(c) = &edge.geom {
                    if (c.radius - t.minor_radius).abs() < 1e-6 {
                        circle_thetas.push(theta_of(c.frame.origin));
                    }
                }
                for vid in [edge.v_start, edge.v_end] {
                    if let Some(v) = self.brep.vertices.get(vid) {
                        vertex_thetas.push(theta_of(v.point));
                    }
                }
            }
        }
        // Use the circle ends if they bracket a non-degenerate arc; else vertices.
        let span_of = |ths: &[f64]| -> Option<(f64, f64)> {
            if ths.len() < 2 {
                return None;
            }
            let pi = std::f64::consts::PI;
            let (sx, sy): (f64, f64) = ths
                .iter()
                .fold((0.0, 0.0), |(cx, cy), &th| (cx + th.cos(), cy + th.sin()));
            let mean = sy.atan2(sx);
            let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
            for &th in ths {
                let mut d = th - mean;
                while d > pi {
                    d -= 2.0 * pi;
                }
                while d < -pi {
                    d += 2.0 * pi;
                }
                lo = lo.min(d);
                hi = hi.max(d);
            }
            Some((mean + lo, mean + hi))
        };
        if let Some((lo, hi)) = span_of(&circle_thetas) {
            if (hi - lo).abs() >= 1e-6 {
                return Some((lo, hi));
            }
        }
        match span_of(&vertex_thetas) {
            Some((lo, hi)) if (hi - lo).abs() >= 1e-9 => Some((lo, hi)),
            _ => None,
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

            // Pass 1: emit every face's carrier surface, recording its id.  This
            // must precede edge emission so a shared edge's SURFACE_CURVE can
            // reference BOTH adjacent faces' surfaces in its pcurves.
            let mut face_surf: HashMap<FaceId, usize> = HashMap::new();
            for &face_id in &shell.faces {
                if let Some(face) = self.brep.faces.get(face_id) {
                    let surf_id = self.emit_face_geometry(ctx, face)?;
                    face_surf.insert(face_id, surf_id);
                }
            }

            // Topological edge→faces adjacency (shared hole edges use one EdgeId).
            let mut edge_faces: HashMap<EdgeId, Vec<FaceId>> = HashMap::new();
            for &face_id in &shell.faces {
                let Some(face) = self.brep.faces.get(face_id) else {
                    continue;
                };
                if !matches!(face.extent, FaceExtent::Trimmed) {
                    continue;
                }
                for &lp in std::iter::once(&face.outer_loop).chain(face.inner_loops.iter()) {
                    if let Some(coedges) = self.brep.loop_coedges(lp) {
                        for ce_id in coedges {
                            if let Some(ce) = self.brep.coedges.get(ce_id) {
                                let v = edge_faces.entry(ce.edge).or_default();
                                if !v.contains(&face_id) {
                                    v.push(face_id);
                                }
                            }
                        }
                    }
                }
            }
            // Pass 1.5: branch-continuous pcurve chains (cursor walk per wire).
            let mut wire_uv: HashMap<(FaceId, EdgeId), Vec<(f64, f64)>> = HashMap::new();
            for &face_id in &shell.faces {
                let Some(face) = self.brep.faces.get(face_id) else {
                    continue;
                };
                if !matches!(face.extent, FaceExtent::Trimmed) {
                    continue;
                }
                let FaceGeom::Cylinder(c) = &face.geom else {
                    continue;
                };
                for &lp in std::iter::once(&face.outer_loop).chain(face.inner_loops.iter()) {
                    let Some(coedges) = self.brep.loop_coedges(lp) else {
                        continue;
                    };
                    let mut cursor: Option<f64> = None;
                    for ce_id in coedges {
                        let Some(ce) = self.brep.coedges.get(ce_id) else {
                            continue;
                        };
                        let Some(edge) = self.brep.edges.get(ce.edge) else {
                            continue;
                        };
                        let (vf, vt) = match ce.sense {
                            CoEdgeSense::Same => (edge.v_start, edge.v_end),
                            CoEdgeSense::Opposite => (edge.v_end, edge.v_start),
                        };
                        let (Some(pf), Some(pt)) =
                            (self.brep.vertices.get(vf), self.brep.vertices.get(vt))
                        else {
                            continue;
                        };
                        let pts: Vec<Point3> = match &edge.geom {
                            EdgeGeom::Polyline(p) => order_polyline(p, pf.point, pt.point),
                            EdgeGeom::Line(_) => vec![pf.point, pt.point],
                            EdgeGeom::Circle(cc) => sample_circle_edge(
                                cc,
                                pf.point,
                                pt.point,
                                edge.v_start == edge.v_end,
                                ce.sense,
                            ),
                            EdgeGeom::Ellipse(e) => sample_ellipse_edge(
                                e,
                                pf.point,
                                pt.point,
                                edge.v_start == edge.v_end,
                                ce.sense,
                            ),
                        };
                        if pts.len() < 2 {
                            continue;
                        }
                        let mut uv: Vec<(f64, f64)> =
                            pts.iter().map(|&p| crate::pcurve::cyl_uv(c, p)).collect();
                        let mut us: Vec<f64> = uv.iter().map(|x| x.0).collect();
                        crate::pcurve::unwrap_u(&mut us);
                        let two_pi = 2.0 * std::f64::consts::PI;
                        let shift = match cursor {
                            Some(cu) => {
                                let mut d = us[0] - cu;
                                let mut s = 0.0;
                                while d > std::f64::consts::PI {
                                    d -= two_pi;
                                    s -= two_pi;
                                }
                                while d < -std::f64::consts::PI {
                                    d += two_pi;
                                    s += two_pi;
                                }
                                s
                            }
                            None => {
                                let mut u0 = us[0];
                                let mut s = 0.0;
                                while u0 < 0.0 {
                                    u0 += two_pi;
                                    s += two_pi;
                                }
                                while u0 >= two_pi {
                                    u0 -= two_pi;
                                    s -= two_pi;
                                }
                                s
                            }
                        };
                        for (k, x) in uv.iter_mut().enumerate() {
                            x.0 = us[k] + shift;
                        }
                        cursor = Some(uv.last().unwrap().0);
                        // Store for polylines AND circle arcs: arcs in mixed
                        // trimmed wires need explicit pcurves too, otherwise
                        // the importer projects them onto its own branch and
                        // the wire opens in parameter space.
                        if matches!(edge.geom, EdgeGeom::Polyline(_) | EdgeGeom::Circle(_)) {
                            wire_uv.insert((face_id, ce.edge), uv);
                        }
                    }
                }
            }
            let pc = PcurveCtx {
                face_surf,
                edge_faces,
                wire_uv,
            };

            // Pass 2: emit each face's bounds + ADVANCED_FACE.
            for &face_id in &shell.faces {
                let face = match self.brep.faces.get(face_id) {
                    Some(f) => f,
                    None => continue,
                };
                let surf_id = match pc.face_surf.get(&face_id) {
                    Some(&s) => s,
                    None => continue,
                };
                let sense = match face.normal {
                    FaceNormal::Same => ".T.",
                    FaceNormal::Reversed => ".F.",
                };
                // Analytic-union faces carry explicit B-Rep topology; walk the
                // real CoEdge/Edge arenas instead of a closed-form template.
                let bounds = if matches!(face.extent, FaceExtent::Trimmed) {
                    self.emit_trimmed_face_bounds(ctx, face, &pc)?
                } else {
                    emit_face_bounds(ctx, face)?
                };
                let af_id = ctx.next_id();
                if bounds.is_empty() {
                    writeln!(
                        ctx.out,
                        "#{af_id} = ADVANCED_FACE('',(),#{surf_id},{sense});"
                    )?;
                } else {
                    let bound_refs: String = bounds
                        .iter()
                        .map(|id| format!("#{id}"))
                        .collect::<Vec<_>>()
                        .join(",");
                    writeln!(
                        ctx.out,
                        "#{af_id} = ADVANCED_FACE('',({bound_refs}),#{surf_id},{sense});"
                    )?;
                }
                ids.push(af_id);
            }
        }
        Ok(ids)
    }

    /// Emit the FACE_OUTER_BOUND (+ FACE_BOUND holes) for a `FaceExtent::Trimmed`
    /// face by walking the **real** CoEdge/Edge topology.
    ///
    /// Every edge is funnelled through the shared-edge cache, so the two faces
    /// meeting at an edge emit one `EDGE_CURVE` referenced by two `ORIENTED_EDGE`s
    /// with opposite senses (the AP214 invariant).
    fn emit_trimmed_face_bounds(
        &self,
        ctx: &mut Ctx,
        face: &cadcore_topo::Face,
        pc: &PcurveCtx,
    ) -> Result<Vec<usize>, StepError> {
        let mut bounds = Vec::new();
        if let Some(b) = self.emit_loop_bound(ctx, face.outer_loop, true, pc)? {
            bounds.push(b);
        }
        for &lp in &face.inner_loops {
            if let Some(b) = self.emit_loop_bound(ctx, lp, false, pc)? {
                bounds.push(b);
            }
        }
        Ok(bounds)
    }

    /// Emit one FACE_(OUTER_)BOUND for a single loop.
    ///
    /// A loop that is a single **closed analytic edge** (full circle / full
    /// ellipse, `v_start == v_end`) is routed through `emit_circle_bound` /
    /// `emit_ellipse_bound` so its EDGE_CURVE shares the geometric cache key with
    /// any adjacent **template** face (e.g. an untouched planar cap) — the two
    /// then reference one EDGE_CURVE with opposite senses.  All other loops are
    /// walked edge-by-edge into ORIENTED_EDGEs.
    fn emit_loop_bound(
        &self,
        ctx: &mut Ctx,
        loop_id: LoopId,
        is_outer: bool,
        pc: &PcurveCtx,
    ) -> Result<Option<usize>, StepError> {
        let coedges = match self.brep.loop_coedges(loop_id) {
            Some(c) => c,
            None => return Ok(None),
        };
        let bound_orientation = self.loop_bound_orientation(loop_id);

        // Full closed circle / ellipse → template-compatible bound.
        if coedges.len() == 1 {
            if let Some(ce) = self.brep.coedges.get(coedges[0]) {
                if let Some(edge) = self.brep.edges.get(ce.edge) {
                    let closed = edge.v_start == edge.v_end;
                    // Match the known-good `FaceExtent::Cylinder` template
                    // convention: outer boundary uses orient=true, inner
                    // boundary uses orient=false (so the tube winds correctly).
                    match &edge.geom {
                        EdgeGeom::Circle(c) if closed => {
                            let mut orient = is_outer;
                            if let Some(face_id) = self.brep.loops.get(loop_id).map(|lp| lp.face) {
                                if let Some(face) = self.brep.faces.get(face_id) {
                                    if let FaceGeom::Torus(t) = &face.geom {
                                        let get_loop_circle_origin =
                                            |brep: &BRep, lp_id: LoopId| -> Option<Point3> {
                                                let co = brep.loop_coedges(lp_id)?;
                                                if co.len() == 1 {
                                                    let ce = brep.coedges.get(co[0])?;
                                                    let edge = brep.edges.get(ce.edge)?;
                                                    if let EdgeGeom::Circle(c) = &edge.geom {
                                                        return Some(c.frame.origin);
                                                    }
                                                }
                                                None
                                            };
                                        let start_origin =
                                            get_loop_circle_origin(self.brep, face.outer_loop);
                                        let end_origin = face
                                            .inner_loops
                                            .first()
                                            .and_then(|&lp| get_loop_circle_origin(self.brep, lp));
                                        if let (Some(so), Some(eo)) = (start_origin, end_origin) {
                                            let local_start = t.frame.to_local_point(so);
                                            let local_end = t.frame.to_local_point(eo);
                                            let theta_start = local_start.y.atan2(local_start.x);
                                            let theta_end = local_end.y.atan2(local_end.x);
                                            let mut d = theta_end - theta_start;
                                            let pi = std::f64::consts::PI;
                                            while d > pi {
                                                d -= 2.0 * pi;
                                            }
                                            while d < -pi {
                                                d += 2.0 * pi;
                                            }
                                            if d < 0.0 {
                                                orient = !is_outer;
                                            }
                                        }
                                    }
                                }
                            }

                            let id = emit_circle_bound_oriented(
                                ctx,
                                c.frame.origin,
                                c.frame.z,
                                c.frame.x,
                                c.radius,
                                is_outer,
                                orient,
                                bound_orientation,
                            )?;
                            return Ok(Some(id));
                        }
                        EdgeGeom::Ellipse(e) if closed => {
                            let id = emit_ellipse_bound_oriented(
                                ctx,
                                e,
                                e.frame.x,
                                is_outer,
                                is_outer,
                                bound_orientation,
                            )?;
                            return Ok(Some(id));
                        }
                        _ => {}
                    }
                }
            }
        }

        // General loop: walk co-edges into shared ORIENTED_EDGEs.
        let mut oes = Vec::new();
        for ce_id in coedges {
            let ce = match self.brep.coedges.get(ce_id) {
                Some(c) => c,
                None => continue,
            };
            let edge = match self.brep.edges.get(ce.edge) {
                Some(e) => e,
                None => continue,
            };
            let (v_from_id, v_to_id) = match ce.sense {
                CoEdgeSense::Same => (edge.v_start, edge.v_end),
                CoEdgeSense::Opposite => (edge.v_end, edge.v_start),
            };
            let p_from = self.brep.vertices.get(v_from_id).map(|v| v.point);
            let p_to = self.brep.vertices.get(v_to_id).map(|v| v.point);
            let (p_from, p_to) = match (p_from, p_to) {
                (Some(a), Some(b)) => (a, b),
                _ => continue,
            };
            match &edge.geom {
                EdgeGeom::Line(_) => push_shared_line(ctx, p_from, p_to, &mut oes)?,
                EdgeGeom::Circle(c) => {
                    self.push_shared_arc_pc(ctx, c, p_from, p_to, ce.edge, pc, &mut oes)?;
                }
                EdgeGeom::Ellipse(e) => push_shared_ellipse(ctx, e, p_from, p_to, &mut oes)?,
                EdgeGeom::Polyline(pts) => {
                    let ordered = order_polyline(pts, p_from, p_to);
                    self.emit_polyline_edge_pc(ctx, ce.edge, &ordered, &mut oes, pc)?;
                }
            }
        }
        if oes.is_empty() {
            return Ok(None);
        }
        let wrapped = finish_bound(ctx, &oes, is_outer, bound_orientation)?;
        Ok(wrapped.into_iter().next())
    }

    /// Compute STEP `FACE_BOUND.orientation` from the loop's signed area in the
    /// face surface parameters.  This deliberately does not touch co-edge
    /// senses: adjacent faces still share each EDGE_CURVE with opposite
    /// ORIENTED_EDGE senses, while STEP gets the per-face loop orientation it
    /// needs for trimming.
    fn loop_bound_orientation(&self, loop_id: LoopId) -> bool {
        let Some(lp) = self.brep.loops.get(loop_id) else {
            return true;
        };
        let Some(face) = self.brep.faces.get(lp.face) else {
            return true;
        };
        let Some(pts) = self.loop_points(loop_id) else {
            return true;
        };
        let Some(mut uv) = self.loop_uv(face, &pts) else {
            return true;
        };
        // Periodicity of the face params: only periodic coords may "wind".
        let (per_u, per_v) = match &face.geom {
            FaceGeom::Cylinder(_) => (true, false),
            FaceGeom::Torus(_) => (true, true),
            _ => (false, false),
        };
        unwrap_loop_uv(&mut uv, per_u, per_v);
        // Net winding (incl. the closing step, mapped to the nearest branch).
        let two_pi = 2.0 * std::f64::consts::PI;
        let winding = |first: f64, last: f64| -> f64 {
            let mut close = first - last;
            while close > std::f64::consts::PI {
                close -= two_pi;
            }
            while close < -std::f64::consts::PI {
                close += two_pi;
            }
            (last - first) + close
        };
        let wu = if per_u { winding(uv[0].0, uv[uv.len() - 1].0) } else { 0.0 };
        let wv = if per_v { winding(uv[0].1, uv[uv.len() - 1].1) } else { 0.0 };
        if wu.abs() > std::f64::consts::PI || wv.abs() > std::f64::consts::PI {
            // BOUNDARY-type loop (winds around the surface period).  Its
            // signed area in (u,v) is degenerate (≈ ±sliver), so derive the
            // orientation from the winding direction + which side of the loop
            // the face material lies on (region-on-the-LEFT ⇔ .T.):
            //   winds +u, material above (+v)  → .T.
            //   winds +v, material at  −u side → .T.
            let (loop_mu, loop_mv) = mean_uv(&uv);
            let (face_mu, face_mv) = self.face_mean_uv(face).unwrap_or((loop_mu, loop_mv));
            if wu.abs() >= wv.abs() {
                return (wu > 0.0) == (face_mv > loop_mv);
            }
            return (wv > 0.0) == (face_mu < loop_mu);
        }
        let area = signed_area_2d(&uv);
        if area.abs() < 1e-9 {
            return true;
        }
        // `.T.` = "the loop, as stored, is correctly wound for its role":
        // OUTER loops must be CCW (area > 0), INNER (hole) loops CW (area < 0).
        // Marking a correctly-CW hole `.F.` makes importers flip it to CCW and
        // read it as an island → "invalid imbrication" / fragmented faces.
        let is_outer = face.outer_loop == loop_id;
        if is_outer {
            area > 0.0
        } else {
            area < 0.0
        }
    }

    /// Mean `(u,v)` over ALL of a face's loops — the "material side" reference
    /// for boundary-loop orientation.  Uses the same per-face projection as
    /// `loop_uv` (cylinder z is absolute; torus θ is band-anchored), so the
    /// means are branch-comparable.
    fn face_mean_uv(&self, face: &cadcore_topo::Face) -> Option<(f64, f64)> {
        let mut acc = (0.0, 0.0);
        let mut n = 0usize;
        for lp in std::iter::once(face.outer_loop).chain(face.inner_loops.iter().copied()) {
            let Some(pts) = self.loop_points(lp) else {
                continue;
            };
            let Some(uv) = self.loop_uv(face, &pts) else {
                continue;
            };
            for (u, v) in uv {
                acc.0 += u;
                acc.1 += v;
                n += 1;
            }
        }
        (n > 0).then(|| (acc.0 / n as f64, acc.1 / n as f64))
    }

    fn loop_points(&self, loop_id: LoopId) -> Option<Vec<Point3>> {
        let coedges = self.brep.loop_coedges(loop_id)?;
        let mut out: Vec<Point3> = Vec::new();
        for ce_id in coedges {
            let ce = self.brep.coedges.get(ce_id)?;
            let edge = self.brep.edges.get(ce.edge)?;
            let (v_from_id, v_to_id) = match ce.sense {
                CoEdgeSense::Same => (edge.v_start, edge.v_end),
                CoEdgeSense::Opposite => (edge.v_end, edge.v_start),
            };
            let p_from = self.brep.vertices.get(v_from_id)?.point;
            let p_to = self.brep.vertices.get(v_to_id)?.point;
            let mut pts = match &edge.geom {
                EdgeGeom::Line(_) => vec![p_from, p_to],
                EdgeGeom::Circle(c) => {
                    sample_circle_edge(c, p_from, p_to, edge.v_start == edge.v_end, ce.sense)
                }
                EdgeGeom::Ellipse(e) => {
                    sample_ellipse_edge(e, p_from, p_to, edge.v_start == edge.v_end, ce.sense)
                }
                EdgeGeom::Polyline(pts) => order_polyline(pts, p_from, p_to),
            };
            if out
                .last()
                .zip(pts.first())
                .map_or(false, |(a, b)| (*a - *b).length() < 1e-8)
            {
                pts.remove(0);
            }
            out.extend(pts);
        }
        (out.len() >= 3).then_some(out)
    }

    fn loop_uv(&self, face: &cadcore_topo::Face, pts: &[Point3]) -> Option<Vec<(f64, f64)>> {
        match &face.geom {
            FaceGeom::Plane(p) => Some(
                pts.iter()
                    .map(|&pt| {
                        let l = p.frame.to_local_point(pt);
                        (l.x, l.y)
                    })
                    .collect(),
            ),
            FaceGeom::Cylinder(c) => Some(
                pts.iter()
                    .map(|&pt| {
                        let l = c.frame.to_local_point(pt);
                        (l.y.atan2(l.x), l.z)
                    })
                    .collect(),
            ),
            FaceGeom::Torus(t) => {
                let (lo, hi) = self.torus_arc_range(face)?;
                Some(
                    pts.iter()
                        .map(|&pt| torus_uv_in_band(t, pt, lo, hi))
                        .collect(),
                )
            }
            FaceGeom::Sphere(_) => None,
        }
    }

    /// Emit one ORIENTED_EDGE for a shared **polyline** edge (a bite-window /
    /// intersection curve), attaching a 2D **pcurve** for every adjacent
    /// cylindrical face so importers do not re-project the curve onto the
    /// periodic surface (the cause of seam-crossing "Unorientable shape" /
    /// "Invalid imbrication of wires" defects).  The 3D curve is wrapped in a
    /// `SURFACE_CURVE` carrying the pcurves; the edge is shared via the cache so
    /// the partner face reuses the same EDGE_CURVE with the opposite sense.
    fn emit_polyline_edge_pc(
        &self,
        ctx: &mut Ctx,
        edge_id: EdgeId,
        pts: &[Point3],
        oe_ids: &mut Vec<usize>,
        pc: &PcurveCtx,
    ) -> Result<(), StepError> {
        if pts.len() < 2 {
            return Ok(());
        }
        let p_from = pts[0];
        let p_to = *pts.last().unwrap();
        let v_from = emit_vertex_point(ctx, p_from)?;
        let v_to = emit_vertex_point(ctx, p_to)?;
        let centroid = {
            let mut acc = Point3::new(0.0, 0.0, 0.0);
            for p in pts {
                acc = Point3::new(acc.x + p.x, acc.y + p.y, acc.z + p.z);
            }
            let n = pts.len() as f64;
            Point3::new(acc.x / n, acc.y / n, acc.z / n)
        };
        let key = StepCurveKey::Polyline {
            v1: v_from.min(v_to),
            v2: v_from.max(v_to),
            mid: point_key(centroid),
        };
        let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
            pair
        } else {
            // 3D piecewise-linear edge curve (degree-1 B-spline; POLYLINE is
            // dropped by SolidWorks' importer).
            let pt_ids: Vec<usize> = pts
                .iter()
                .map(|p| emit_point(ctx, *p, "pl"))
                .collect::<Result<_, _>>()?;
            let pl = emit_polyline3d(ctx, &pt_ids)?;
            // pcurve on each adjacent CYLINDER face (torus pcurves: Phase 5).
            let mut pcurves = Vec::new();
            if let Some(faces) = pc.edge_faces.get(&edge_id) {
                for &f in faces {
                    let (Some(face), Some(&surf_id)) =
                        (self.brep.faces.get(f), pc.face_surf.get(&f))
                    else {
                        continue;
                    };
                    match &face.geom {
                        FaceGeom::Cylinder(c) => {
                            // Prefer the branch-continuous wire chain (cursor
                            // pre-pass) — falls back to per-edge principal
                            // anchoring when absent.  Stored chains are in the
                            // FACE's traversal order; the emitted 3D polyline
                            // is in THIS coedge's order — align directions.
                            let uv: Vec<(f64, f64)> = match pc.wire_uv.get(&(f, edge_id)) {
                                Some(stored) if stored.len() == pts.len() => {
                                    let raw0 = crate::pcurve::cyl_uv(c, pts[0]);
                                    let du = (stored[0].0 - raw0.0).rem_euclid(2.0 * std::f64::consts::PI);
                                    let ang_ok = du < 1e-6 || du > 2.0 * std::f64::consts::PI - 1e-6;
                                    let v_ok = (stored[0].1 - raw0.1).abs() < 1e-6;
                                    if ang_ok && v_ok {
                                        stored.clone()
                                    } else {
                                        stored.iter().rev().copied().collect()
                                    }
                                }
                                _ => cyl_pcurve(c, pts, Some(std::f64::consts::PI)),
                            };
                            let c2d = emit_polyline2d(ctx, &uv)?;
                            pcurves.push(emit_pcurve(ctx, surf_id, c2d)?);
                        }
                        FaceGeom::Torus(t) => {
                            if torus_is_toroidal(face) {
                                // Windowed torus emitted as analytic
                                // TOROIDAL_SURFACE: pcurve is native (θ,φ),
                                // unwrapped for continuity across the seams.
                                if let Some((lo, hi)) = self.torus_arc_range(face) {
                                    let mut uv: Vec<(f64, f64)> =
                                        pts.iter().map(|&p| torus_uv_in_band(t, p, lo, hi)).collect();
                                    unwrap_torus_pcurve(&mut uv);
                                    let c2d = emit_polyline2d(ctx, &uv)?;
                                    pcurves.push(emit_pcurve(ctx, surf_id, c2d)?);
                                }
                            } else if let Some((lo, hi)) = self.torus_arc_range(face) {
                                // Elbow face: reconstruct the same NURBS patch the
                                // surface was emitted as, then invert each 3D
                                // point to its (u,v) for the shared trim edge.
                                let patch = torus_patch_nurbs(t, lo, hi);
                                let uv: Vec<(f64, f64)> =
                                    pts.iter().map(|&p| patch.project_point(p)).collect();
                                let c2d = emit_polyline2d(ctx, &uv)?;
                                pcurves.push(emit_pcurve(ctx, surf_id, c2d)?);
                            }
                        }
                        _ => {}
                    }
                }
            }
            let geom_ref = if pcurves.is_empty() {
                pl
            } else {
                emit_surface_curve(ctx, pl, &pcurves, false)?
            };
            let ec = ctx.next_id();
            writeln!(
                ctx.out,
                "#{ec} = EDGE_CURVE('',#{v_from},#{v_to},#{geom_ref},.T.);"
            )?;
            ctx.edge_cache.insert(key, (ec, v_from));
            (ec, v_from)
        };
        let sense = if orig_start == v_from { ".T." } else { ".F." };
        let oe = ctx.next_id();
        writeln!(ctx.out, "#{oe} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
        oe_ids.push(oe);
        Ok(())
    }
}

impl<'a> StepWriter<'a> {
    /// `push_shared_arc` with explicit PCURVEs for each adjacent trimmed
    /// CYLINDER face (from the branch-continuous `wire_uv` chains).  Arcs in
    /// MIXED wires (notched boundaries, strips) must carry pcurves on the same
    /// branch as the wire's polylines — otherwise the importer projects them
    /// onto its own principal branch and the wire opens in parameter space.
    fn push_shared_arc_pc(
        &self,
        ctx: &mut Ctx,
        c: &cadcore_geom::Circle3,
        p_from: Point3,
        p_to: Point3,
        edge_id: EdgeId,
        pc: &PcurveCtx,
        oe_ids: &mut Vec<usize>,
    ) -> Result<(), StepError> {
        let (centre, axis, xref, r) = (c.frame.origin, c.frame.z, c.frame.x, c.radius);
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
            let xref_v = xref.as_vec();
            let yref_v = axis.as_vec().cross(xref_v);
            let dot =
                |a: cadcore_math::Vec3, b: cadcore_math::Vec3| a.x * b.x + a.y * b.y + a.z * b.z;
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
            // Explicit pcurves on every adjacent trimmed cylinder face.
            let mut pcurves = Vec::new();
            if let Some(faces) = pc.edge_faces.get(&edge_id) {
                for &f in faces {
                    let (Some(face), Some(&surf_id)) =
                        (self.brep.faces.get(f), pc.face_surf.get(&f))
                    else {
                        continue;
                    };
                    let FaceGeom::Cylinder(cyl) = &face.geom else {
                        continue;
                    };
                    let Some(stored) = pc.wire_uv.get(&(f, edge_id)) else {
                        continue;
                    };
                    if stored.len() < 2 {
                        continue;
                    }
                    // Align the stored (face-traversal) chain to THIS edge's
                    // 3D direction (p_from → p_to).
                    let raw0 = crate::pcurve::cyl_uv(cyl, p_from);
                    let near = |s: (f64, f64)| -> bool {
                        let du = (s.0 - raw0.0).rem_euclid(2.0 * std::f64::consts::PI);
                        (du < 1e-5 || du > 2.0 * std::f64::consts::PI - 1e-5)
                            && (s.1 - raw0.1).abs() < 1e-5
                    };
                    let uv: Vec<(f64, f64)> = if near(stored[0]) {
                        stored.clone()
                    } else if near(*stored.last().unwrap()) {
                        stored.iter().rev().copied().collect()
                    } else {
                        continue;
                    };
                    let c2d = emit_polyline2d(ctx, &uv)?;
                    pcurves.push(emit_pcurve(ctx, surf_id, c2d)?);
                }
            }
            let geom_ref = if pcurves.is_empty() {
                ci
            } else {
                emit_surface_curve(ctx, ci, &pcurves, false)?
            };
            let ec = ctx.next_id();
            writeln!(
                ctx.out,
                "#{ec} = EDGE_CURVE('',#{v_from},#{v_to},#{geom_ref},{same_sense});"
            )?;
            ctx.edge_cache.insert(key, (ec, v_from));
            (ec, v_from)
        };
        let sense = if orig_start == v_from { ".T." } else { ".F." };
        let oe = ctx.next_id();
        writeln!(ctx.out, "#{oe} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
        oe_ids.push(oe);
        Ok(())
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
            let torus = match &face.geom {
                FaceGeom::Torus(t) => t,
                _ => return Ok(vec![]),
            };

            let local_start = torus.frame.to_local_point(start_circle.frame.origin);
            let local_end = torus.frame.to_local_point(end_circle.frame.origin);
            let theta_start = local_start.y.atan2(local_start.x);
            let theta_end = local_end.y.atan2(local_end.x);

            let mut d = theta_end - theta_start;
            let pi = std::f64::consts::PI;
            while d > pi {
                d -= 2.0 * pi;
            }
            while d < -pi {
                d += 2.0 * pi;
            }

            let (c_outer, c_inner) = if d >= 0.0 {
                (start_circle, end_circle)
            } else {
                (end_circle, start_circle)
            };

            let outer = emit_circle_bound(
                ctx,
                c_outer.frame.origin,
                c_outer.frame.z,
                c_outer.frame.x,
                c_outer.radius,
                true,
                true,
            )?;
            let inner = emit_circle_bound(
                ctx,
                c_inner.frame.origin,
                c_inner.frame.z,
                c_inner.frame.x,
                c_inner.radius,
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
        FaceExtent::RoundedRectCap {
            xmin,
            xmax,
            zmin,
            zmax,
            radius,
            y,
            plus_y,
        } => emit_rounded_rect_cap_bounds(ctx, *xmin, *xmax, *zmin, *zmax, *radius, *y, *plus_y),

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

        // Handled earlier in `emit_advanced_faces` via the real topology walk
        // (`emit_trimmed_face_bounds`); never reached through this template path.
        FaceExtent::Trimmed => Ok(vec![]),
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
    emit_circle_bound_oriented(ctx, centre, normal, _x_dir, radius, outer, orient, true)
}

fn emit_circle_bound_oriented(
    ctx: &mut Ctx,
    centre: Point3,
    normal: UnitVec3,
    _x_dir: UnitVec3,
    radius: f64,
    outer: bool,
    orient: bool,
    bound_orientation: bool,
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
    let bound_sense = if bound_orientation { ".T." } else { ".F." };
    writeln!(ctx.out, "#{fb_id} = {btype}('',#{el_id},{bound_sense});")?;
    Ok(fb_id)
}

fn emit_ellipse_bound(
    ctx: &mut Ctx,
    ellipse: &Ellipse3,
    _fallback_x_dir: UnitVec3,
    outer: bool,
    orient: bool,
) -> Result<usize, StepError> {
    emit_ellipse_bound_oriented(ctx, ellipse, _fallback_x_dir, outer, orient, true)
}

fn emit_ellipse_bound_oriented(
    ctx: &mut Ctx,
    ellipse: &Ellipse3,
    _fallback_x_dir: UnitVec3,
    outer: bool,
    orient: bool,
    bound_orientation: bool,
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
    let bound_sense = if bound_orientation { ".T." } else { ".F." };
    writeln!(ctx.out, "#{fb_id} = {btype}('',#{el_id},{bound_sense});")?;
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
        writeln!(
            ctx.out,
            "#{ec} = EDGE_CURVE('',#{v_from},#{v_to},#{li},.T.);"
        )?;
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
        writeln!(
            ctx.out,
            "#{ec} = EDGE_CURVE('',#{v_from},#{v_to},#{ci},{same_sense});"
        )?;
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
    finish_bound(ctx, oe_ids, true, true)
}

fn finish_bound(
    ctx: &mut Ctx,
    oe_ids: &[usize],
    outer: bool,
    bound_orientation: bool,
) -> Result<Vec<usize>, StepError> {
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
    let btype = if outer {
        "FACE_OUTER_BOUND"
    } else {
        "FACE_BOUND"
    };
    let sense = if bound_orientation { ".T." } else { ".F." };
    writeln!(ctx.out, "#{fb} = {btype}('',#{el},{sense});")?;
    Ok(vec![fb])
}

/// Emit an EDGE_LOOP + FACE_BOUND (an *inner* hole loop) from ORIENTED_EDGEs.
/// Order a stored polyline so it runs `p_from → p_to` (the co-edge direction).
fn order_polyline(pts: &[Point3], p_from: Point3, _p_to: Point3) -> Vec<Point3> {
    if pts.is_empty() {
        return vec![];
    }
    let d_first = (pts[0] - p_from).length();
    let d_last = (pts[pts.len() - 1] - p_from).length();
    if d_last < d_first {
        pts.iter().rev().copied().collect()
    } else {
        pts.to_vec()
    }
}

fn mean_uv(uv: &[(f64, f64)]) -> (f64, f64) {
    let n = uv.len().max(1) as f64;
    let s = uv
        .iter()
        .fold((0.0, 0.0), |a, &(u, v)| (a.0 + u, a.1 + v));
    (s.0 / n, s.1 / n)
}

fn signed_area_2d(uv: &[(f64, f64)]) -> f64 {
    if uv.len() < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    for i in 0..uv.len() {
        let (x0, y0) = uv[i];
        let (x1, y1) = uv[(i + 1) % uv.len()];
        area += x0 * y1 - x1 * y0;
    }
    area * 0.5
}

fn unwrap_loop_uv(uv: &mut [(f64, f64)], per_u: bool, per_v: bool) {
    // Only PERIODIC axes may be branch-unwrapped.  Unwrapping a non-periodic
    // axis (a cylinder's axial v) corrupts loops whose seam edge spans more
    // than π along it (e.g. a full-length seam on a band > π mm long): the
    // seam step gets wrongly shifted by 2π, mangling the polygon and flipping
    // its signed area → the outer loop is mis-oriented and SolidWorks drops
    // the whole face.
    for i in 1..uv.len() {
        if per_u {
            while uv[i].0 - uv[i - 1].0 > std::f64::consts::PI {
                uv[i].0 -= 2.0 * std::f64::consts::PI;
            }
            while uv[i].0 - uv[i - 1].0 < -std::f64::consts::PI {
                uv[i].0 += 2.0 * std::f64::consts::PI;
            }
        }
        if per_v {
            while uv[i].1 - uv[i - 1].1 > std::f64::consts::PI {
                uv[i].1 -= 2.0 * std::f64::consts::PI;
            }
            while uv[i].1 - uv[i - 1].1 < -std::f64::consts::PI {
                uv[i].1 += 2.0 * std::f64::consts::PI;
            }
        }
    }
}

/// Whether a torus face should be emitted as an analytic `TOROIDAL_SURFACE`
/// (with native (θ,φ) pcurves) rather than the rational-NURBS arc patch: true
/// when it has been cut into windows, where the NURBS patch's numerically
/// inverted pcurves are unreliable and OCC reports "Unorientable".
fn torus_is_toroidal(face: &cadcore_topo::Face) -> bool {
    // Any *cut* torus (Trimmed extent) carries native (θ,φ) loops, so emit it
    // analytically — whether it was windowed (inner loops) or sliced into a
    // φ-band by a coaxial stack (outer loop only).  The plain scaffold fillet
    // keeps FaceExtent::TorusFillet and still rides the rational-NURBS patch.
    matches!(face.extent, FaceExtent::Trimmed)
}

/// Make a torus pcurve continuous: unwrap each point's (θ,φ) onto the branch
/// nearest the previous point, so periodic jumps don't tear the 2-D curve.
fn unwrap_torus_pcurve(uv: &mut [(f64, f64)]) {
    let tau = 2.0 * std::f64::consts::PI;
    for i in 1..uv.len() {
        let du = uv[i].0 - uv[i - 1].0;
        uv[i].0 -= tau * (du / tau).round();
        let dv = uv[i].1 - uv[i - 1].1;
        uv[i].1 -= tau * (dv / tau).round();
    }
}

fn torus_uv_in_band(
    t: &cadcore_geom::TorusSurf,
    p: Point3,
    theta_lo: f64,
    theta_hi: f64,
) -> (f64, f64) {
    let l = t.frame.to_local_point(p);
    let mut theta = l.y.atan2(l.x);
    while theta < theta_lo - std::f64::consts::PI {
        theta += 2.0 * std::f64::consts::PI;
    }
    while theta > theta_hi + std::f64::consts::PI {
        theta -= 2.0 * std::f64::consts::PI;
    }
    if theta < theta_lo {
        let alt = theta + 2.0 * std::f64::consts::PI;
        if (theta_lo..=theta_hi).contains(&alt) {
            theta = alt;
        }
    } else if theta > theta_hi {
        let alt = theta - 2.0 * std::f64::consts::PI;
        if (theta_lo..=theta_hi).contains(&alt) {
            theta = alt;
        }
    }
    let radial = (l.x * l.x + l.y * l.y).sqrt() - t.major_radius;
    let phi = l.z.atan2(radial);
    (theta, phi)
}

fn sample_circle_edge(
    c: &cadcore_geom::Circle3,
    p_from: Point3,
    p_to: Point3,
    closed: bool,
    sense: CoEdgeSense,
) -> Vec<Point3> {
    let mut a0 = circle_param(c, p_from);
    let mut a1 = if closed {
        a0 + 2.0 * std::f64::consts::PI
    } else {
        circle_param(c, p_to)
    };
    if !closed {
        while a1 - a0 > std::f64::consts::PI {
            a1 -= 2.0 * std::f64::consts::PI;
        }
        while a1 - a0 < -std::f64::consts::PI {
            a1 += 2.0 * std::f64::consts::PI;
        }
    }
    if matches!(sense, CoEdgeSense::Opposite) {
        std::mem::swap(&mut a0, &mut a1);
    }
    let n = if closed { 32 } else { 8 };
    (0..=n)
        .map(|i| c.point_at(a0 + (a1 - a0) * i as f64 / n as f64))
        .collect()
}

fn sample_ellipse_edge(
    e: &Ellipse3,
    p_from: Point3,
    p_to: Point3,
    closed: bool,
    sense: CoEdgeSense,
) -> Vec<Point3> {
    let mut a0 = ellipse_param(e, p_from);
    let mut a1 = if closed {
        a0 + 2.0 * std::f64::consts::PI
    } else {
        ellipse_param(e, p_to)
    };
    if !closed {
        while a1 - a0 > std::f64::consts::PI {
            a1 -= 2.0 * std::f64::consts::PI;
        }
        while a1 - a0 < -std::f64::consts::PI {
            a1 += 2.0 * std::f64::consts::PI;
        }
    }
    if matches!(sense, CoEdgeSense::Opposite) {
        std::mem::swap(&mut a0, &mut a1);
    }
    let n = if closed { 32 } else { 8 };
    (0..=n)
        .map(|i| e.point_at(a0 + (a1 - a0) * i as f64 / n as f64))
        .collect()
}

fn circle_param(c: &cadcore_geom::Circle3, p: Point3) -> f64 {
    let d = p - c.frame.origin;
    c.frame.y.dot_vec(d).atan2(c.frame.x.dot_vec(d))
}

fn ellipse_param(e: &Ellipse3, p: Point3) -> f64 {
    let d = p - e.frame.origin;
    let x = e.frame.x.dot_vec(d) / e.semi_major;
    let y = e.frame.y.dot_vec(d) / e.semi_minor;
    y.atan2(x)
}

/// Append one ORIENTED_EDGE for a polyline edge `pts[0] → pts[last]`, reusing a
/// shared POLYLINE/EDGE_CURVE from the cache when the partner face emitted it.
fn push_shared_polyline(
    ctx: &mut Ctx,
    pts: &[Point3],
    oe_ids: &mut Vec<usize>,
) -> Result<(), StepError> {
    if pts.len() < 2 {
        return Ok(());
    }
    let p_from = pts[0];
    let p_to = *pts.last().unwrap();
    let v_from = emit_vertex_point(ctx, p_from)?;
    let v_to = emit_vertex_point(ctx, p_to)?;
    // Order-independent centroid: distinguishes two arcs with the same endpoints.
    let centroid = {
        let mut acc = Point3::new(0.0, 0.0, 0.0);
        for p in pts {
            acc = Point3::new(acc.x + p.x, acc.y + p.y, acc.z + p.z);
        }
        let n = pts.len() as f64;
        Point3::new(acc.x / n, acc.y / n, acc.z / n)
    };
    let key = StepCurveKey::Polyline {
        v1: v_from.min(v_to),
        v2: v_from.max(v_to),
        mid: point_key(centroid),
    };
    let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
        pair
    } else {
        let pt_ids: Vec<usize> = pts
            .iter()
            .map(|p| emit_point(ctx, *p, "pl"))
            .collect::<Result<_, _>>()?;
        let pl = emit_polyline3d(ctx, &pt_ids)?;
        let ec = ctx.next_id();
        writeln!(
            ctx.out,
            "#{ec} = EDGE_CURVE('',#{v_from},#{v_to},#{pl},.T.);"
        )?;
        ctx.edge_cache.insert(key, (ec, v_from));
        (ec, v_from)
    };
    let sense = if orig_start == v_from { ".T." } else { ".F." };
    let oe = ctx.next_id();
    writeln!(ctx.out, "#{oe} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
    oe_ids.push(oe);
    Ok(())
}

/// Append one ORIENTED_EDGE for an elliptic arc `p_from → p_to` on `ellipse`,
/// reusing a shared ELLIPSE/EDGE_CURVE from the cache.  `same_sense` is derived
/// from the ellipse's parametric direction (CCW in the `(x,y)` frame).
fn push_shared_ellipse(
    ctx: &mut Ctx,
    ellipse: &Ellipse3,
    p_from: Point3,
    p_to: Point3,
    oe_ids: &mut Vec<usize>,
) -> Result<(), StepError> {
    let v_from = emit_vertex_point(ctx, p_from)?;
    let v_to = emit_vertex_point(ctx, p_to)?;

    let key = {
        let (vs, ve) = if point_key(p_from) <= point_key(p_to) {
            (point_key(p_from), point_key(p_to))
        } else {
            (point_key(p_to), point_key(p_from))
        };
        StepCurveKey::EllipseArc {
            v_start: vs,
            v_end: ve,
            center: point_key(ellipse.frame.origin),
            semi_major_micro: (ellipse.semi_major * 1_000_000.0).round() as i64,
            semi_minor_micro: (ellipse.semi_minor * 1_000_000.0).round() as i64,
            normal: dir_key(ellipse.frame.z),
        }
    };

    let (ec_id, orig_start) = if let Some(&pair) = ctx.edge_cache.get(&key) {
        pair
    } else {
        // CCW parametric angle in the ellipse frame (normalise axes to a circle).
        let fx = ellipse.frame.x.as_vec();
        let fy = ellipse.frame.y.as_vec();
        let dot = |a: cadcore_math::Vec3, b: cadcore_math::Vec3| a.x * b.x + a.y * b.y + a.z * b.z;
        let ang = |p: Point3| -> f64 {
            let d = p - ellipse.frame.origin;
            (dot(d, fy) / ellipse.semi_minor).atan2(dot(d, fx) / ellipse.semi_major)
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

        let el_id = emit_ellipse(ctx, ellipse)?;
        let ec = ctx.next_id();
        writeln!(
            ctx.out,
            "#{ec} = EDGE_CURVE('',#{v_from},#{v_to},#{el_id},{same_sense});"
        )?;
        ctx.edge_cache.insert(key, (ec, v_from));
        (ec, v_from)
    };
    let sense = if orig_start == v_from { ".T." } else { ".F." };
    let oe = ctx.next_id();
    writeln!(ctx.out, "#{oe} = ORIENTED_EDGE('',*,*,#{ec_id},{sense});")?;
    oe_ids.push(oe);
    Ok(())
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
        let local =
            arc_ref_dir.as_vec() * (radius * angle.cos()) + right.as_vec() * (radius * angle.sin());
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
    push_shared_arc(
        ctx,
        centre_hi,
        axis,
        arc_ref_dir,
        radius,
        p_s1,
        p_e1,
        &mut oe_ids,
    )?; // top arc
    push_shared_line(ctx, p_e1, p_e0, &mut oe_ids)?; // T_e → B_e  (down)
    push_shared_arc(
        ctx,
        centre_lo,
        axis,
        arc_ref_dir,
        radius,
        p_e0,
        p_s0,
        &mut oe_ids,
    )?; // bottom arc (backward)

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
    xmin: f64,
    xmax: f64,
    zmin: f64,
    zmax: f64,
    r: f64,
    y: f64,
    plus_y: bool,
) -> Result<Vec<usize>, StepError> {
    use cadcore_math::Vec3;

    // Profile junction points (CCW in XZ as seen from −Y).
    let pts = [
        Point3::new(xmin + r, y, zmin), // 0
        Point3::new(xmax - r, y, zmin), // 1
        Point3::new(xmax, y, zmin + r), // 2
        Point3::new(xmax, y, zmax - r), // 3
        Point3::new(xmax - r, y, zmax), // 4
        Point3::new(xmin + r, y, zmax), // 5
        Point3::new(xmin, y, zmax - r), // 6
        Point3::new(xmin, y, zmin + r), // 7
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
    use super::{brep_to_step, StepWriter};
    use cadcore_geom::{Circle3, Plane3, TorusSurf};
    use cadcore_math::{Point3, UnitVec3};
    use cadcore_ops::{build_solid_rounded_box_xz, sweep_circle_along_polyline, SweepOptions};
    use cadcore_topo::{
        BRep, CoEdge, CoEdgeId, CoEdgeSense, Edge, EdgeGeom, Face, FaceExtent, FaceGeom,
        FaceNormal, Loop, LoopId, Shell, Solid,
    };

    /// Count occurrences of a STEP keyword in the output.
    fn count(step: &str, keyword: &str) -> usize {
        step.matches(keyword).count()
    }

    fn add_single_circle_loop(brep: &mut BRep, circle: Circle3) -> LoopId {
        let seam = circle.point_at(0.0);
        let v = brep.add_vertex(cadcore_topo::Vertex { point: seam });
        let e = brep.add_edge(Edge {
            geom: EdgeGeom::Circle(circle),
            v_start: v,
            v_end: v,
            t_start: 0.0,
            t_end: std::f64::consts::TAU,
            partner: None,
        });
        let ce = brep.add_coedge(CoEdge {
            edge: e,
            sense: CoEdgeSense::Same,
            next: CoEdgeId::default(),
            prev: CoEdgeId::default(),
            loop_id: LoopId::default(),
        });
        brep.patch_coedge_links(&[ce]);
        brep.add_loop(Loop {
            start: ce,
            face: cadcore_topo::FaceId::default(),
        })
    }

    fn add_polyline_loop(brep: &mut BRep, pts: &[Point3]) -> LoopId {
        let verts: Vec<_> = pts
            .iter()
            .map(|&point| brep.add_vertex(cadcore_topo::Vertex { point }))
            .collect();
        let mut coedges = Vec::with_capacity(pts.len());
        for i in 0..pts.len() {
            let j = (i + 1) % pts.len();
            let edge = brep.add_edge(Edge {
                geom: EdgeGeom::Polyline(vec![pts[i], pts[j]]),
                v_start: verts[i],
                v_end: verts[j],
                t_start: 0.0,
                t_end: 1.0,
                partner: None,
            });
            coedges.push(brep.add_coedge(CoEdge {
                edge,
                sense: CoEdgeSense::Same,
                next: CoEdgeId::default(),
                prev: CoEdgeId::default(),
                loop_id: LoopId::default(),
            }));
        }
        brep.patch_coedge_links(&coedges);
        let loop_id = brep.add_loop(Loop {
            start: coedges[0],
            face: cadcore_topo::FaceId::default(),
        });
        for ce in coedges {
            brep.coedges.get_mut(ce).unwrap().loop_id = loop_id;
        }
        loop_id
    }

    #[test]
    fn trimmed_inner_loop_uses_face_bound_orientation_not_edge_flip() {
        let mut brep = BRep::new();
        let outer = add_polyline_loop(
            &mut brep,
            &[
                Point3::new(-2.0, -2.0, 0.0),
                Point3::new(2.0, -2.0, 0.0),
                Point3::new(2.0, 2.0, 0.0),
                Point3::new(-2.0, 2.0, 0.0),
            ],
        );
        let inner_cw = add_polyline_loop(
            &mut brep,
            &[
                Point3::new(-1.0, -1.0, 0.0),
                Point3::new(-1.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(1.0, -1.0, 0.0),
            ],
        );
        let solid = brep.add_solid(Solid {
            shells: vec![],
            name: Some("bound_orientation".into()),
        });
        let shell = brep.add_shell(Shell {
            faces: vec![],
            is_outer: true,
            solid,
        });
        let face_id = brep.add_face(Face {
            geom: FaceGeom::Plane(Plane3::from_origin_normal(Point3::ORIGIN, UnitVec3::Z)),
            normal: FaceNormal::Same,
            outer_loop: outer,
            inner_loops: vec![inner_cw],
            shell,
            extent: FaceExtent::Trimmed,
        });
        for lp in [outer, inner_cw] {
            brep.loops.get_mut(lp).unwrap().face = face_id;
        }
        brep.shells.get_mut(shell).unwrap().faces = vec![face_id];
        brep.solids.get_mut(solid).unwrap().shells = vec![shell];

        let step = brep_to_step(&brep).unwrap();
        // A CW-stored hole is the CONVENTIONAL winding for an inner bound, so
        // FACE_BOUND.orientation must be .T. ("use as stored").  Marking it
        // .F. makes importers reverse it to CCW and read the hole as an
        // island ("invalid imbrication" / fragmented faces).
        assert!(
            step.lines()
                .any(|line| line.contains("= FACE_BOUND") && line.contains(",.T.)")),
            "CW inner loop is conventionally wound → FACE_BOUND.orientation=.T.:\n{step}"
        );
        assert!(
            !step.lines()
                .any(|line| line.contains("= FACE_BOUND") && line.contains(",.F.)")),
            "no inner bound should need reversal here:\n{step}"
        );
    }

    #[test]
    fn torus_arc_range_uses_endpoint_loops_not_extra_holes() {
        let mut brep = BRep::new();
        let torus = TorusSurf::new(Point3::ORIGIN, UnitVec3::Z, 2.0, 0.25);

        let minor_circle = |theta: f64| {
            let radial = torus.frame.x * theta.cos() + torus.frame.y * theta.sin();
            let center = torus.frame.origin + radial * torus.major_radius;
            let tangent = UnitVec3::try_from_vec(torus.frame.z.as_vec().cross(radial)).unwrap();
            Circle3::new(center, tangent, torus.minor_radius)
        };

        let start = 170.0_f64.to_radians();
        let end = -170.0_f64.to_radians();
        let outer = add_single_circle_loop(&mut brep, minor_circle(start));
        let end_loop = add_single_circle_loop(&mut brep, minor_circle(end));
        let unrelated_hole = add_single_circle_loop(&mut brep, minor_circle(0.0));

        let solid = brep.add_solid(Solid {
            shells: vec![],
            name: Some("torus_range".into()),
        });
        let shell = brep.add_shell(Shell {
            faces: vec![],
            is_outer: true,
            solid,
        });
        let face_id = brep.add_face(Face {
            geom: FaceGeom::Torus(torus),
            normal: FaceNormal::Same,
            outer_loop: outer,
            inner_loops: vec![end_loop, unrelated_hole],
            shell,
            extent: FaceExtent::Trimmed,
        });
        let face = brep.faces.get(face_id).unwrap();

        let (lo, hi) = StepWriter::new(&brep).torus_arc_range(face).unwrap();
        let sweep = hi - lo;
        assert!(
            (sweep - 20.0_f64.to_radians()).abs() < 1e-9,
            "range must follow endpoint loops, got sweep {} deg",
            sweep.to_degrees()
        );
    }

    /// Collect all EDGE_CURVE ids and the ORIENTED_EDGE ids that reference them.
    /// Returns a map: ec_id → list of (oe_id, sense).
    fn oriented_edge_uses(step: &str) -> std::collections::HashMap<usize, Vec<(usize, &str)>> {
        let mut map: std::collections::HashMap<usize, Vec<(usize, &str)>> =
            std::collections::HashMap::new();
        // Match: #N = ORIENTED_EDGE('',*,*,#M,.T.); or .F.
        let re_pat = "ORIENTED_EDGE";
        for line in step.lines() {
            if !line.contains(re_pat) {
                continue;
            }
            // parse #OE = ORIENTED_EDGE('',*,*,#EC,.S.);
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() < 4 {
                continue;
            }
            let ec_part = parts[parts.len() - 2].trim();
            let sense_part = parts[parts.len() - 1].trim().trim_end_matches(';');
            let ec_id: usize = ec_part.trim_start_matches('#').parse().ok().unwrap_or(0);
            // OE id from start of line
            let oe_id_str = line
                .trim()
                .split('=')
                .next()
                .unwrap_or("")
                .trim()
                .trim_start_matches('#');
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
            if !step.contains(&ec_line) {
                continue;
            }
            if refs.len() != 2 {
                violations.push(format!(
                    "EC#{ec_id} used {} times (want 2): {:?}",
                    refs.len(),
                    refs
                ));
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

    // ── Trimmed faces (analytic-union output) ─────────────────────────────────

    /// Build a closed tetrahedron entirely from `FaceExtent::Trimmed` planar
    /// faces wired with real CoEdge/Edge topology.  Exercises the topology-
    /// driven `emit_trimmed_face_bounds` path and its shared-edge cache.
    ///
    /// `polyline_edges`: when true every edge is an `EdgeGeom::Polyline`
    /// (endpoints + a collinear midpoint) so the `POLYLINE` emit path is hit;
    /// otherwise edges are `EdgeGeom::Line`.
    fn build_trimmed_tetra(polyline_edges: bool) -> BRep {
        use cadcore_geom::Line3;
        use cadcore_geom::Plane3;
        use cadcore_math::UnitVec3;
        use cadcore_topo::{
            CoEdge, CoEdgeId, CoEdgeSense, Edge, EdgeGeom, Face, FaceExtent, FaceGeom, FaceId,
            FaceNormal, Loop, LoopId, Shell, Solid, Vertex,
        };
        use std::collections::HashMap;

        let mut brep = BRep::new();
        let pos = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        let vids: Vec<_> = pos
            .iter()
            .map(|&p| brep.add_vertex(Vertex { point: p }))
            .collect();
        let centroid = Point3::new(0.25, 0.25, 0.25);
        let tris = [[1usize, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

        let solid_id = brep.add_solid(Solid {
            shells: vec![],
            name: Some("tetra".into()),
        });
        let shell_id = brep.add_shell(Shell {
            faces: vec![],
            is_outer: true,
            solid: solid_id,
        });

        // Undirected edge -> EdgeId, created on first use (v_start/v_end fixed then).
        let mut edge_map: HashMap<(usize, usize), cadcore_topo::EdgeId> = HashMap::new();
        let mut face_ids = Vec::new();

        for tri in tris {
            // Wind CCW as seen from outside (normal points away from centroid).
            let mut v = tri;
            let (a, b, c) = (pos[v[0]], pos[v[1]], pos[v[2]]);
            let n = (b - a).cross(c - a);
            let fc = Point3::new(
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            );
            let outward = fc - centroid;
            if n.x * outward.x + n.y * outward.y + n.z * outward.z < 0.0 {
                v.swap(1, 2);
            }
            let (a, b, c) = (pos[v[0]], pos[v[1]], pos[v[2]]);
            let normal = UnitVec3::try_from_vec((b - a).cross(c - a)).unwrap();
            let fc = Point3::new(
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            );

            // Three directed edges of the wound triangle.
            let mut coedge_ids = Vec::new();
            for k in 0..3 {
                let vi = v[k];
                let vj = v[(k + 1) % 3];
                let key = (vi.min(vj), vi.max(vj));
                let edge_id = *edge_map.entry(key).or_insert_with(|| {
                    let (s, e) = (key.0, key.1);
                    let (ps, pe) = (pos[s], pos[e]);
                    let geom = if polyline_edges {
                        let mid = Point3::new(
                            (ps.x + pe.x) / 2.0,
                            (ps.y + pe.y) / 2.0,
                            (ps.z + pe.z) / 2.0,
                        );
                        EdgeGeom::Polyline(vec![ps, mid, pe])
                    } else {
                        let dir = UnitVec3::try_from_vec(pe - ps).unwrap();
                        EdgeGeom::Line(Line3::new(ps, dir))
                    };
                    brep.add_edge(Edge {
                        geom,
                        v_start: vids[s],
                        v_end: vids[e],
                        t_start: 0.0,
                        t_end: (pe - ps).length(),
                        partner: None,
                    })
                });
                // Co-edge traverses vi->vj; sense relative to stored edge (key.0->key.1).
                let sense = if vi == key.0 {
                    CoEdgeSense::Same
                } else {
                    CoEdgeSense::Opposite
                };
                let ce = brep.add_coedge(CoEdge {
                    edge: edge_id,
                    sense,
                    next: CoEdgeId::default(),
                    prev: CoEdgeId::default(),
                    loop_id: LoopId::default(),
                });
                coedge_ids.push(ce);
            }
            brep.patch_coedge_links(&coedge_ids);

            let loop_id = brep.add_loop(Loop {
                start: coedge_ids[0],
                face: FaceId::default(),
            });
            let face_id = brep.add_face(Face {
                geom: FaceGeom::Plane(Plane3::from_origin_normal(fc, normal)),
                normal: FaceNormal::Same,
                outer_loop: loop_id,
                inner_loops: vec![],
                shell: shell_id,
                extent: FaceExtent::Trimmed,
            });
            if let Some(lp) = brep.loops.get_mut(loop_id) {
                lp.face = face_id;
            }
            face_ids.push(face_id);
        }

        if let Some(sh) = brep.shells.get_mut(shell_id) {
            sh.faces = face_ids;
        }
        if let Some(sd) = brep.solids.get_mut(solid_id) {
            sd.shells = vec![shell_id];
        }
        brep
    }

    #[test]
    fn trimmed_tetrahedron_line_edges_is_manifold() {
        let brep = build_trimmed_tetra(false);
        let step = brep_to_step(&brep).expect("emit STEP");
        assert_eq!(count(&step, "ADVANCED_FACE"), 4, "4 faces");
        assert_eq!(count(&step, "= EDGE_CURVE"), 6, "6 shared edges");
        check_ap214_manifold(&step).expect("tetra (line edges) must be AP214 manifold");
    }

    #[test]
    fn trimmed_tetrahedron_polyline_edges_is_manifold() {
        let brep = build_trimmed_tetra(true);
        let step = brep_to_step(&brep).expect("emit STEP");
        assert_eq!(count(&step, "ADVANCED_FACE"), 4, "4 faces");
        assert_eq!(count(&step, "= EDGE_CURVE"), 6, "6 shared edges");
        // Piecewise-linear edges are emitted as degree-1 B-splines (POLYLINE
        // entities are dropped by SolidWorks' importer).
        assert!(
            count(&step, "B_SPLINE_CURVE_WITH_KNOTS('',1,") >= 6,
            "deg-1 b-spline edges emitted"
        );
        assert_eq!(count(&step, "POLYLINE"), 0, "no POLYLINE entities");
        check_ap214_manifold(&step).expect("tetra (polyline edges) must be AP214 manifold");
    }

    // ── Straight cylinder + caps ──────────────────────────────────────────────

    #[test]
    fn open_boundary_split_exports_one_shared_intersection_edge() {
        use cadcore_geom::{
            swept_swept_intersection, CenterlineSeg, Line3, Plane3, SweptTubeSurface,
        };
        use cadcore_math::{UnitVec3, Vec3};
        use cadcore_ops::union::{split_face_by_open_edge_at_existing_vertices, OpenBoundaryKeep};
        use cadcore_topo::{
            CoEdge, CoEdgeId, CoEdgeSense, Edge, EdgeGeom, Face, FaceExtent, FaceGeom, FaceId,
            FaceNormal, Loop, LoopId, Shell, Solid, Vertex, VertexId,
        };

        fn add_edge_loop_face(
            brep: &mut BRep,
            shell_id: cadcore_topo::ShellId,
            vids: &[VertexId],
            points: &[Point3],
            normal: UnitVec3,
        ) -> cadcore_topo::FaceId {
            let mut coedges = Vec::new();
            for i in 0..vids.len() {
                let j = (i + 1) % vids.len();
                let dir = UnitVec3::try_from_vec(points[j] - points[i]).unwrap();
                let edge = brep.add_edge(Edge {
                    geom: EdgeGeom::Line(Line3::new(points[i], dir)),
                    v_start: vids[i],
                    v_end: vids[j],
                    t_start: 0.0,
                    t_end: (points[j] - points[i]).length(),
                    partner: None,
                });
                coedges.push(brep.add_coedge(CoEdge {
                    edge,
                    sense: CoEdgeSense::Same,
                    next: CoEdgeId::default(),
                    prev: CoEdgeId::default(),
                    loop_id: LoopId::default(),
                }));
            }
            brep.patch_coedge_links(&coedges);
            let loop_id = brep.add_loop(Loop {
                start: coedges[0],
                face: FaceId::default(),
            });
            for &ce in &coedges {
                brep.coedges.get_mut(ce).unwrap().loop_id = loop_id;
            }
            let face = brep.add_face(Face {
                geom: FaceGeom::Plane(Plane3::from_origin_normal(points[0], normal)),
                normal: FaceNormal::Same,
                outer_loop: loop_id,
                inner_loops: vec![],
                shell: shell_id,
                extent: FaceExtent::Trimmed,
            });
            brep.loops.get_mut(loop_id).unwrap().face = face;
            face
        }

        let mut brep = BRep::new();
        let solid = brep.add_solid(Solid {
            shells: vec![],
            name: Some("open_split".into()),
        });
        let shell = brep.add_shell(Shell {
            faces: vec![],
            is_outer: true,
            solid,
        });

        let radius = 0.275;
        let elbow = SweptTubeSurface::new(
            vec![CenterlineSeg::Arc {
                center: Point3::new(0.9, 19.225, 0.35),
                axis: UnitVec3::try_from_vec(Vec3::new(0.0, 0.0, -1.0)).unwrap(),
                xref: UnitVec3::try_from_vec(Vec3::new(0.0, -1.0, 0.0)).unwrap(),
                radius: 0.5,
                ang0: 1.5708,
                ang1: 3.1416,
            }],
            radius,
        );
        let leg_dir =
            UnitVec3::try_from_vec(Point3::new(0.775, 19.6, 0.7) - Point3::new(19.725, 19.6, 0.7))
                .unwrap();
        let leg = SweptTubeSurface::new(
            vec![CenterlineSeg::Line {
                p0: Point3::new(19.725, 19.6, 0.7),
                dir: leg_dir,
                length: 18.95,
            }],
            radius,
        );
        let curves = swept_swept_intersection(&elbow, &leg, 96, 48);
        let open_count = curves.iter().filter(|c| !c.closed).count();
        assert!(
            open_count > 0,
            "controlled elbow x leg must produce open SSI"
        );
        let open = curves
            .into_iter()
            .find(|c| !c.closed && c.points.len() >= 6)
            .expect("usable open SSI curve");
        let a = *open.points.first().unwrap();
        let b = *open.points.last().unwrap();
        let chord = UnitVec3::try_from_vec(b - a).expect("non-degenerate open SSI");
        let ref_axis = if chord.dot(UnitVec3::Z).abs() < 0.9 {
            UnitVec3::Z
        } else {
            UnitVec3::Y
        };
        let side1 = UnitVec3::try_from_vec(chord.as_vec().cross(ref_axis.as_vec()))
            .unwrap()
            .as_vec()
            * (radius * 3.0);
        let side2 = UnitVec3::try_from_vec(chord.as_vec().cross(side1))
            .unwrap()
            .as_vec()
            * (radius * 3.0);
        let p0_pt = a - side1;
        let p1_pt = b - side1;
        let p2_pt = b + side1;
        let p3_pt = a + side1;
        let r0_pt = a - side2;
        let r1_pt = b - side2;
        let r2_pt = b + side2;
        let r3_pt = a + side2;
        let normal1 = UnitVec3::try_from_vec((p1_pt - p0_pt).cross(p2_pt - p0_pt)).unwrap();
        let normal2 = UnitVec3::try_from_vec((r1_pt - r0_pt).cross(r2_pt - r0_pt)).unwrap();
        let va = brep.add_vertex(Vertex { point: a });
        let vb = brep.add_vertex(Vertex { point: b });
        let p0 = brep.add_vertex(Vertex { point: p0_pt });
        let p1 = brep.add_vertex(Vertex { point: p1_pt });
        let p2 = brep.add_vertex(Vertex { point: p2_pt });
        let p3 = brep.add_vertex(Vertex { point: p3_pt });
        let r0 = brep.add_vertex(Vertex { point: r0_pt });
        let r1 = brep.add_vertex(Vertex { point: r1_pt });
        let r2 = brep.add_vertex(Vertex { point: r2_pt });
        let r3 = brep.add_vertex(Vertex { point: r3_pt });

        let face_xy = add_edge_loop_face(
            &mut brep,
            shell,
            &[p0, p1, p2, p3],
            &[p0_pt, p1_pt, p2_pt, p3_pt],
            normal1,
        );
        let face_xz = add_edge_loop_face(
            &mut brep,
            shell,
            &[r0, r1, r2, r3],
            &[r0_pt, r1_pt, r2_pt, r3_pt],
            normal2,
        );
        brep.shells.get_mut(shell).unwrap().faces = vec![face_xy, face_xz];
        brep.solids.get_mut(solid).unwrap().shells = vec![shell];

        let cutter = brep.add_edge(Edge {
            geom: EdgeGeom::Polyline(open.points),
            v_start: va,
            v_end: vb,
            t_start: 0.0,
            t_end: 1.0,
            partner: None,
        });

        let xy = split_face_by_open_edge_at_existing_vertices(
            &mut brep,
            face_xy,
            cutter,
            OpenBoundaryKeep::StartToEnd,
        )
        .expect("split xy face");
        let xz = split_face_by_open_edge_at_existing_vertices(
            &mut brep,
            face_xz,
            cutter,
            OpenBoundaryKeep::EndToStart,
        )
        .expect("split xz face");

        assert_eq!(xy.boundary_coedges, 3);
        assert_eq!(xy.dropped_boundary_coedges, 3);
        assert_eq!(xz.boundary_coedges, 3);
        assert_eq!(xz.dropped_boundary_coedges, 3);
        assert_eq!(brep.shells[shell].faces.len(), 2, "dropped regions removed");
        assert!(brep.faces.get(face_xy).is_none());
        assert!(brep.faces.get(face_xz).is_none());
        assert_eq!(brep.loop_coedges(xy.loop_id).unwrap().len(), 4);
        assert_eq!(brep.loop_coedges(xz.loop_id).unwrap().len(), 4);

        let xy_open = &brep.coedges[xy.open_coedge];
        let xz_open = &brep.coedges[xz.open_coedge];
        assert_eq!(xy_open.edge, cutter);
        assert_eq!(xz_open.edge, cutter);
        assert_ne!(xy_open.sense, xz_open.sense);

        let step = brep_to_step(&brep).expect("emit STEP");
        assert_eq!(count(&step, "ADVANCED_FACE"), 2);
        assert_eq!(
            count(&step, "B_SPLINE_CURVE_WITH_KNOTS('',1,"),
            1,
            "single shared SSI curve (deg-1 b-spline)"
        );
        assert_eq!(count(&step, "POLYLINE"), 0, "no POLYLINE entities");
        let uses = oriented_edge_uses(&step);
        let shared: Vec<_> = uses
            .iter()
            .filter(|(ec_id, refs)| {
                step.contains(&format!("#{} = EDGE_CURVE", ec_id)) && refs.len() == 2
            })
            .collect();
        assert_eq!(shared.len(), 1, "only the cutter edge is shared");
        let senses: Vec<_> = shared[0].1.iter().map(|(_, sense)| *sense).collect();
        assert!(senses.iter().any(|s| s.starts_with(".T.")));
        assert!(senses.iter().any(|s| s.starts_with(".F.")));
    }

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
        )
        .unwrap();
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
            0.0,
            10.0, // xmin xmax
            0.0,
            1.0, // ymin ymax
            0.0,
            8.0, // zmin zmax
            1.0, // corner_radius
            Some("electrode".into()),
        )
        .unwrap();
        let step = brep_to_step(&brep).unwrap();

        assert_eq!(
            count(&step, "ADVANCED_FACE"),
            10,
            "rounded box must produce exactly 10 faces"
        );
        assert_eq!(
            count(&step, "CYLINDRICAL_SURFACE"),
            4,
            "four quarter-cylinder corners"
        );

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
            .map(|w| SweepPathSegment::Line {
                start: w[0],
                end: w[1],
            })
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
        build_solid_rounded_box_xz(&mut brep, 0.0, 10.0, 0.0, 2.0, 0.0, 8.0, 1.5, None).unwrap();
        let step = brep_to_step(&brep).unwrap();
        let uses = oriented_edge_uses(&step);
        // Every arc EDGE_CURVE (referenced from a CIRCLE) that is shared must
        // have exactly 2 ORIENTED_EDGE references.
        let mut dup_arcs = 0usize;
        for (ec_id, refs) in &uses {
            let ec_line = format!("#{ec_id} = EDGE_CURVE");
            if !step.contains(&ec_line) {
                continue;
            }
            if refs.len() > 2 {
                dup_arcs += 1;
            }
        }
        assert_eq!(
            dup_arcs, 0,
            "some arc EDGE_CURVEs are referenced more than twice"
        );
    }
}

fn normal_key(v: UnitVec3) -> [i64; 3] {
    // Quantise to 1e-4 (≈0.006°).  Coarse enough that the tiny tangent
    // differences (~1e-6) between a corner torus and its adjacent cylinder end
    // map to the SAME key (so their coincident boundary circles share one
    // EDGE_CURVE — watertight rounded corners), yet far finer than any genuine
    // normal difference in a scaffold (caps/circles differ by ≥90°).
    //
    // ⚠ Must NOT use 1e-6 here: at that precision a numerically-noisy component
    // like −1.0e-6 rounds to −1 and triggers the sign-normalisation flip below,
    // giving two identical circles different keys → open shell at every corner.
    let vec = v.as_vec();
    let mut x = (vec.x * 10_000.0).round() as i64;
    let mut y = (vec.y * 10_000.0).round() as i64;
    let mut z = (vec.z * 10_000.0).round() as i64;
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
