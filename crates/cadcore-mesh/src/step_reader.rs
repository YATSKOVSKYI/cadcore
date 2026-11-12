use std::collections::HashMap;
use crate::kernel::{
    math::{Point3, Vec3, UnitVec3, Frame3},
    geom::{Plane3, CylSurf, SphereSurf, TorusSurf, Circle3, Ellipse3, Line3},
    topo::{BRep, Vertex, Edge, EdgeGeom, CoEdge, CoEdgeSense, Loop, Face, FaceGeom, FaceNormal, Shell, Solid, FaceExtent, FaceBoundary},
    topo::{VertexId, EdgeId, CoEdgeId, LoopId, FaceId, ShellId, SolidId},
};

/// Parse raw STEP arguments inside parentheses, handling nested groups.
fn split_args(s: &str) -> Vec<String> {
    let mut args = Vec::new();
    let mut current = String::new();
    let mut depth = 0;
    let mut in_string = false;
    for c in s.chars() {
        match c {
            '\'' => {
                in_string = !in_string;
                current.push(c);
            }
            '(' if !in_string => {
                depth += 1;
                current.push(c);
            }
            ')' if !in_string => {
                depth -= 1;
                current.push(c);
            }
            ',' if !in_string && depth == 0 => {
                args.push(current.trim().to_string());
                current.clear();
            }
            _ => {
                current.push(c);
            }
        }
    }
    if !current.is_empty() {
        args.push(current.trim().to_string());
    }
    args
}

/// Parse reference IDs, e.g. "#123" -> 123.
fn parse_ref(s: &str) -> Result<usize, String> {
    let trimmed = s.trim().trim_start_matches('#');
    trimmed.parse::<usize>().map_err(|e| format!("invalid ref '{}': {}", s, e))
}

/// Parse list of references, e.g. "(#10,#20,#30)" -> [10, 20, 30].
fn parse_ref_list(s: &str) -> Result<Vec<usize>, String> {
    let trimmed = s.trim();
    if !trimmed.starts_with('(') || !trimmed.ends_with(')') {
        return Err(format!("invalid list format: {}", s));
    }
    let inner = &trimmed[1..trimmed.len() - 1];
    if inner.trim().is_empty() {
        return Ok(vec![]);
    }
    inner.split(',').map(|p| parse_ref(p)).collect()
}

/// Parse float coordinates, e.g. "(1.0, 2.0, 3.0)" -> [f64; 3].
fn parse_coords(s: &str) -> Result<[f64; 3], String> {
    let trimmed = s.trim();
    if !trimmed.starts_with('(') || !trimmed.ends_with(')') {
        return Err(format!("invalid coordinate format: {}", s));
    }
    let inner = &trimmed[1..trimmed.len() - 1];
    let parts: Vec<&str> = inner.split(',').map(|p| p.trim()).collect();
    if parts.len() != 3 {
        return Err(format!("expected 3 coordinates, got {}: {}", parts.len(), s));
    }
    let x = parts[0].parse::<f64>().map_err(|e| format!("invalid X: {}", e))?;
    let y = parts[1].parse::<f64>().map_err(|e| format!("invalid Y: {}", e))?;
    let z = parts[2].parse::<f64>().map_err(|e| format!("invalid Z: {}", e))?;
    Ok([x, y, z])
}

#[derive(Debug)]
struct RawEntity {
    name: String,
    args: Vec<String>,
}

/// Parse a STEP AP203/214 string into a cadcore BRep.
pub fn read_step_to_brep(contents: &str) -> Result<BRep, String> {
    let mut raw_entities: HashMap<usize, RawEntity> = HashMap::new();

    // 1. Line-by-line scanning & comment stripping
    let mut in_comment = false;
    let mut buffer = String::new();
    for line in contents.lines() {
        let line_str = line.trim();
        if line_str.is_empty() {
            continue;
        }

        // Extremely simple comment stripper
        let mut clean_line = String::new();
        let chars: Vec<char> = line_str.chars().collect();
        let mut i = 0;
        while i < chars.len() {
            if in_comment {
                if i + 1 < chars.len() && chars[i] == '*' && chars[i + 1] == '/' {
                    in_comment = false;
                    i += 2;
                } else {
                    i += 1;
                }
            } else {
                if i + 1 < chars.len() && chars[i] == '/' && chars[i + 1] == '*' {
                    in_comment = true;
                    i += 2;
                } else {
                    clean_line.push(chars[i]);
                    i += 1;
                }
            }
        }

        buffer.push_str(&clean_line);
        if buffer.ends_with(';') {
            let stmt = buffer.trim_end_matches(';').trim();
            if stmt.starts_with('#') {
                if let Some(eq_idx) = stmt.find('=') {
                    let id_str = stmt[..eq_idx].trim().trim_start_matches('#');
                    if let Ok(id) = id_str.parse::<usize>() {
                        let rhs = stmt[eq_idx + 1..].trim();
                        if let Some(paren_idx) = rhs.find('(') {
                            if rhs.ends_with(')') {
                                let type_name = rhs[..paren_idx].trim().to_string();
                                let args_str = &rhs[paren_idx + 1..rhs.len() - 1];
                                let args = split_args(args_str);
                                raw_entities.insert(id, RawEntity { name: type_name, args });
                            }
                        }
                    }
                }
            }
            buffer.clear();
        }
    }

    let mut brep = BRep::new();

    // Resolvers for geometric primitives
    let mut points: HashMap<usize, Point3> = HashMap::new();
    let mut vectors: HashMap<usize, Vec3> = HashMap::new();
    let mut directions: HashMap<usize, UnitVec3> = HashMap::new();
    let mut frames: HashMap<usize, Frame3> = HashMap::new();

    // Resolve CARTESIAN_POINT
    for (&id, ent) in &raw_entities {
        if ent.name == "CARTESIAN_POINT" {
            if ent.args.len() >= 2 {
                let coords = parse_coords(&ent.args[1])?;
                points.insert(id, Point3::new(coords[0], coords[1], coords[2]));
            }
        }
    }

    // Resolve DIRECTION
    for (&id, ent) in &raw_entities {
        if ent.name == "DIRECTION" {
            if ent.args.len() >= 2 {
                let coords = parse_coords(&ent.args[1])?;
                let v = Vec3::new(coords[0], coords[1], coords[2]);
                vectors.insert(id, v);
                if let Some(u) = UnitVec3::try_from_vec(v) {
                    directions.insert(id, u);
                }
            }
        }
    }

    // Resolve AXIS2_PLACEMENT_3D
    let get_dir = |ref_id: usize, fallback: UnitVec3| -> UnitVec3 {
        directions.get(&ref_id).copied().unwrap_or(fallback)
    };

    for (&id, ent) in &raw_entities {
        if ent.name == "AXIS2_PLACEMENT_3D" {
            if ent.args.len() >= 4 {
                let origin_ref = parse_ref(&ent.args[1])?;
                let origin = points.get(&origin_ref).copied().unwrap_or(Point3::ORIGIN);

                let z_ref = ent.args[2].trim();
                let z = if z_ref == "$" {
                    UnitVec3::Z
                } else {
                    get_dir(parse_ref(z_ref)?, UnitVec3::Z)
                };

                let x_ref = ent.args[3].trim();
                let x = if x_ref == "$" {
                    // Compute orthogonal X
                    let n_vec = z.as_vec();
                    let ref_vec = if n_vec.x.abs() > 0.9 {
                        Vec3::new(0.0, 1.0, 0.0)
                    } else {
                        Vec3::new(1.0, 0.0, 0.0)
                    };
                    UnitVec3::try_from_vec(n_vec.cross(ref_vec)).unwrap_or(UnitVec3::X)
                } else {
                    get_dir(parse_ref(x_ref)?, UnitVec3::X)
                };

                let y = UnitVec3::try_from_vec(z.as_vec().cross(x.as_vec())).unwrap_or(UnitVec3::Y);
                frames.insert(id, Frame3 { origin, x, y, z });
            }
        }
    }

    // Reconstruct B-Rep topology maps
    let mut vertex_map: HashMap<usize, VertexId> = HashMap::new();
    let mut edge_map: HashMap<usize, EdgeId> = HashMap::new();
    let mut coedge_map: HashMap<usize, CoEdgeId> = HashMap::new();
    let mut loop_map: HashMap<usize, LoopId> = HashMap::new();
    let mut face_map: HashMap<usize, FaceId> = HashMap::new();
    let mut shell_map: HashMap<usize, ShellId> = HashMap::new();

    // 1. VERTEX_POINT -> Vertex
    for (&id, ent) in &raw_entities {
        if ent.name == "VERTEX_POINT" {
            if ent.args.len() >= 2 {
                let pt_ref = parse_ref(&ent.args[1])?;
                let pt = points.get(&pt_ref).copied().unwrap_or(Point3::ORIGIN);
                let vid = brep.add_vertex(Vertex { point: pt });
                vertex_map.insert(id, vid);
            }
        }
    }

    // Auxiliary curves for edge geometry
    let mut circles: HashMap<usize, Circle3> = HashMap::new();
    let mut ellipses: HashMap<usize, Ellipse3> = HashMap::new();
    let mut lines: HashMap<usize, Line3> = HashMap::new();

    for (&id, ent) in &raw_entities {
        match ent.name.as_str() {
            "LINE" => {
                if ent.args.len() >= 3 {
                    let pt_ref = parse_ref(&ent.args[1])?;
                    let origin = points.get(&pt_ref).copied().unwrap_or(Point3::ORIGIN);
                    let dir_ref = parse_ref(&ent.args[2])?;
                    let dir = directions.get(&dir_ref).copied().unwrap_or(UnitVec3::X);
                    lines.insert(id, Line3::new(origin, dir));
                }
            }
            "CIRCLE" => {
                if ent.args.len() >= 3 {
                    let frame_ref = parse_ref(&ent.args[1])?;
                    let frame = frames.get(&frame_ref).copied().unwrap_or(Frame3::WORLD);
                    let radius = ent.args[2].parse::<f64>().unwrap_or(1.0);
                    circles.insert(id, Circle3::new(frame.origin, frame.z, radius));
                }
            }
            "ELLIPSE" => {
                if ent.args.len() >= 4 {
                    let frame_ref = parse_ref(&ent.args[1])?;
                    let frame = frames.get(&frame_ref).copied().unwrap_or(Frame3::WORLD);
                    let r_maj = ent.args[2].parse::<f64>().unwrap_or(1.0);
                    let r_min = ent.args[3].parse::<f64>().unwrap_or(1.0);
                    ellipses.insert(id, Ellipse3::new(frame.origin, frame.z, frame.x, r_maj, r_min));
                }
            }
            _ => {}
        }
    }

    // 2. EDGE_CURVE -> Edge
    for (&id, ent) in &raw_entities {
        if ent.name == "EDGE_CURVE" {
            if ent.args.len() >= 5 {
                let v_start_ref = parse_ref(&ent.args[1])?;
                let v_end_ref = parse_ref(&ent.args[2])?;
                let curve_ref = parse_ref(&ent.args[3])?;

                let v_start = *vertex_map.get(&v_start_ref).ok_or_else(|| format!("vertex #{} not found for edge #{}", v_start_ref, id))?;
                let v_end = *vertex_map.get(&v_end_ref).ok_or_else(|| format!("vertex #{} not found for edge #{}", v_end_ref, id))?;

                let mut resolved_curve_ref = curve_ref;
                while let Some(curve_ent) = raw_entities.get(&resolved_curve_ref) {
                    if curve_ent.name == "SURFACE_CURVE" {
                        if curve_ent.args.len() >= 2 {
                            if let Ok(next_ref) = parse_ref(&curve_ent.args[1]) {
                                resolved_curve_ref = next_ref;
                                continue;
                            }
                        }
                    }
                    break;
                }

                let geom = if let Some(line) = lines.get(&resolved_curve_ref) {
                    EdgeGeom::Line(*line)
                } else if let Some(circ) = circles.get(&resolved_curve_ref) {
                    EdgeGeom::Circle(*circ)
                } else if let Some(ell) = ellipses.get(&resolved_curve_ref) {
                    EdgeGeom::Ellipse(*ell)
                } else {
                    // Fallback: approximate as a straight line between start and end vertices
                    let p_start = brep.vertices.get(v_start).unwrap().point;
                    let p_end = brep.vertices.get(v_end).unwrap().point;
                    let d_vec = p_end - p_start;
                    let dir = UnitVec3::try_from_vec(d_vec).unwrap_or(UnitVec3::X);
                    EdgeGeom::Line(Line3::new(p_start, dir))
                };

                let t_start = 0.0;
                let t_end = match &geom {
                    EdgeGeom::Line(l) => {
                        let p_start = brep.vertices.get(v_start).unwrap().point;
                        let p_end = brep.vertices.get(v_end).unwrap().point;
                        l.direction.dot_vec(p_end - p_start)
                    }
                    EdgeGeom::Circle(_) | EdgeGeom::Ellipse(_) => std::f32::consts::TAU as f64,
                    // The reader does not currently parse POLYLINE into a
                    // Polyline edge; parametrise defensively by chord length.
                    EdgeGeom::Polyline(pts) => {
                        pts.windows(2).map(|w| (w[1] - w[0]).length()).sum()
                    }
                };

                let eid = brep.add_edge(Edge {
                    geom,
                    v_start,
                    v_end,
                    t_start,
                    t_end,
                    partner: None,
                });
                edge_map.insert(id, eid);
            }
        }
    }

    // 3. ORIENTED_EDGE -> CoEdge (linking temporarily)
    let mut raw_coedges: HashMap<usize, (usize, bool)> = HashMap::new(); // coedge_id -> (edge_ref, sense)
    for (&id, ent) in &raw_entities {
        if ent.name == "ORIENTED_EDGE" {
            if ent.args.len() >= 5 {
                let edge_ref = parse_ref(&ent.args[3])?;
                let sense_str = ent.args[4].trim();
                let sense = sense_str == ".T.";
                raw_coedges.insert(id, (edge_ref, sense));
            }
        }
    }

    // 4. EDGE_LOOP -> Loop (inserts CoEdges in the Loop)
    for (&id, ent) in &raw_entities {
        if ent.name == "EDGE_LOOP" {
            if ent.args.len() >= 2 {
                let coedge_refs = parse_ref_list(&ent.args[1])?;
                let mut coedge_ids = Vec::new();

                let loop_id = brep.add_loop(Loop {
                    start: CoEdgeId::default(),
                    face: FaceId::default(),
                });

                for &ce_ref in &coedge_refs {
                    let &(edge_ref, sense) = raw_coedges.get(&ce_ref)
                        .ok_or_else(|| format!("oriented edge #{} not found for loop #{}", ce_ref, id))?;
                    let edge = *edge_map.get(&edge_ref)
                        .ok_or_else(|| format!("edge #{} not found for loop #{}", edge_ref, id))?;
                    
                    let ce_sense = if sense { CoEdgeSense::Same } else { CoEdgeSense::Opposite };
                    let ce_id = brep.add_coedge(CoEdge {
                        edge,
                        sense: ce_sense,
                        next: CoEdgeId::default(),
                        prev: CoEdgeId::default(),
                        loop_id,
                    });
                    coedge_map.insert(ce_ref, ce_id);
                    coedge_ids.push(ce_id);

                    // Update edge partner backlink
                    if let Some(ed) = brep.edges.get_mut(edge) {
                        ed.partner = Some(ce_id);
                    }
                }

                if !coedge_ids.is_empty() {
                    brep.patch_coedge_links(&coedge_ids);
                    if let Some(lp) = brep.loops.get_mut(loop_id) {
                        lp.start = coedge_ids[0];
                    }
                }
                loop_map.insert(id, loop_id);
            }
        }
    }

    // FACE_OUTER_BOUND / FACE_BOUND wrapper maps
    let mut bound_loops: HashMap<usize, (usize, bool)> = HashMap::new(); // bound_id -> (loop_ref, is_outer)
    for (&id, ent) in &raw_entities {
        if ent.name == "FACE_OUTER_BOUND" {
            if ent.args.len() >= 3 {
                let loop_ref = parse_ref(&ent.args[1])?;
                bound_loops.insert(id, (loop_ref, true));
            }
        } else if ent.name == "FACE_BOUND" {
            if ent.args.len() >= 3 {
                let loop_ref = parse_ref(&ent.args[1])?;
                bound_loops.insert(id, (loop_ref, false));
            }
        }
    }

    // Resolve carrier surfaces for faces
    let mut carrier_planes: HashMap<usize, Plane3> = HashMap::new();
    let mut carrier_cylinders: HashMap<usize, CylSurf> = HashMap::new();
    let mut carrier_spheres: HashMap<usize, SphereSurf> = HashMap::new();
    let mut carrier_tori: HashMap<usize, TorusSurf> = HashMap::new();

    for (&id, ent) in &raw_entities {
        match ent.name.as_str() {
            "PLANE" => {
                if ent.args.len() >= 2 {
                    let frame_ref = parse_ref(&ent.args[1])?;
                    let frame = frames.get(&frame_ref).copied().unwrap_or(Frame3::WORLD);
                    carrier_planes.insert(id, Plane3 { frame });
                }
            }
            "CYLINDRICAL_SURFACE" => {
                if ent.args.len() >= 3 {
                    let frame_ref = parse_ref(&ent.args[1])?;
                    let frame = frames.get(&frame_ref).copied().unwrap_or(Frame3::WORLD);
                    let radius = ent.args[2].parse::<f64>().unwrap_or(1.0);
                    carrier_cylinders.insert(id, CylSurf::new(frame.origin, frame.z, radius));
                }
            }
            "SPHERICAL_SURFACE" => {
                if ent.args.len() >= 3 {
                    let cp_ref = parse_ref(&ent.args[1])?;
                    let center = points.get(&cp_ref).copied().unwrap_or(Point3::ORIGIN);
                    let radius = ent.args[2].parse::<f64>().unwrap_or(1.0);
                    carrier_spheres.insert(id, SphereSurf::new(center, radius));
                }
            }
            "TOROIDAL_SURFACE" => {
                if ent.args.len() >= 4 {
                    let frame_ref = parse_ref(&ent.args[1])?;
                    let frame = frames.get(&frame_ref).copied().unwrap_or(Frame3::WORLD);
                    let r_maj = ent.args[2].parse::<f64>().unwrap_or(1.0);
                    let r_min = ent.args[3].parse::<f64>().unwrap_or(1.0);
                    carrier_tori.insert(id, TorusSurf::new(frame.origin, frame.z, r_maj, r_min));
                }
            }
            _ => {}
        }
    }

    // Helper: reconstruct FaceExtent for a face based on its geometric boundaries
    fn reconstruct_face_extent(brep: &BRep, geom: &FaceGeom, loops: &[LoopId]) -> FaceExtent {
        if loops.is_empty() {
            return FaceExtent::None;
        }

        match geom {
            FaceGeom::Plane(_p) => {
                // If the loop contains circular/elliptic bounds or polygonal vertices
                let loop_id = loops[0];
                if let Some(coedges) = brep.loop_coedges(loop_id) {
                    if coedges.len() == 1 {
                        // Single circular/elliptic boundary
                        let ce = brep.coedges.get(coedges[0]).unwrap();
                        let edge = brep.edges.get(ce.edge).unwrap();
                        match &edge.geom {
                            EdgeGeom::Circle(c) => FaceExtent::Disk { radius: c.radius },
                            EdgeGeom::Ellipse(el) => FaceExtent::PlanarBoundary {
                                boundary: FaceBoundary::Ellipse(*el),
                            },
                            _ => FaceExtent::None,
                        }
                    } else {
                        // Polygon: gather points of the vertices in order
                        let mut poly_pts = Vec::new();
                        for ce_id in coedges {
                            let ce = brep.coedges.get(ce_id).unwrap();
                            let edge = brep.edges.get(ce.edge).unwrap();
                            let start_vid = if ce.sense == CoEdgeSense::Same { edge.v_start } else { edge.v_end };
                            if let Some(vtx) = brep.vertices.get(start_vid) {
                                poly_pts.push(vtx.point);
                            }
                        }
                        if poly_pts.len() >= 3 {
                            FaceExtent::Polygon { points: poly_pts }
                        } else {
                            FaceExtent::None
                        }
                    }
                } else {
                    FaceExtent::None
                }
            }
            FaceGeom::Cylinder(c) => {
                // Cylinders usually have start and end loops.
                // If it is a full cylinder bounded by two circles
                if loops.len() == 2 {
                    let mut circles_found = Vec::new();
                    for &lp_id in loops {
                        if let Some(coedges) = brep.loop_coedges(lp_id) {
                            if coedges.len() == 1 {
                                let ce = brep.coedges.get(coedges[0]).unwrap();
                                let edge = brep.edges.get(ce.edge).unwrap();
                                if let EdgeGeom::Circle(circ) = &edge.geom {
                                    circles_found.push(FaceBoundary::Circle(*circ));
                                }
                            }
                        }
                    }
                    if circles_found.len() == 2 {
                        // Compute length as distance along cylinder axis between centers
                        let c1 = match &circles_found[0] { FaceBoundary::Circle(ci) => ci.frame.origin, _ => Point3::ORIGIN };
                        let c2 = match &circles_found[1] { FaceBoundary::Circle(ci) => ci.frame.origin, _ => Point3::ORIGIN };
                        let length = c.frame.z.dot_vec(c2 - c1).abs();
                        // Put start/end boundaries correctly based on projection
                        let proj1 = c.frame.z.dot_vec(c1 - c.frame.origin);
                        let (start, end) = if proj1 < 1.0e-5 {
                            (circles_found[0].clone(), circles_found[1].clone())
                        } else {
                            (circles_found[1].clone(), circles_found[0].clone())
                        };
                        return FaceExtent::Cylinder { length, start, end };
                    }
                }
                
                // Otherwise check if it is a partial cylinder (chord cut)
                if loops.len() == 1 {
                    let loop_id = loops[0];
                    if let Some(coedges) = brep.loop_coedges(loop_id) {
                        // Gather vertices to compute length
                        let mut min_z = f64::MAX;
                        let mut max_z = f64::MIN;
                        for ce_id in coedges {
                            let ce = brep.coedges.get(ce_id).unwrap();
                            let edge = brep.edges.get(ce.edge).unwrap();
                            let start_vid = if ce.sense == CoEdgeSense::Same { edge.v_start } else { edge.v_end };
                            if let Some(vtx) = brep.vertices.get(start_vid) {
                                let z_val = c.frame.z.dot_vec(vtx.point - c.frame.origin);
                                min_z = min_z.min(z_val);
                                max_z = max_z.max(z_val);
                            }
                        }
                        if min_z.is_finite() && max_z.is_finite() {
                            let length = max_z - min_z;
                            // Approximate angles as full sector for voxelization safety
                            return FaceExtent::PartialCylinder {
                                length,
                                arc_start_angle: -std::f32::consts::FRAC_PI_2 as f64,
                                arc_end_angle: std::f32::consts::FRAC_PI_2 as f64,
                                arc_ref_dir: c.frame.x,
                            };
                        }
                    }
                }
                FaceExtent::None
            }
            FaceGeom::Sphere(_) | FaceGeom::Torus(_) => FaceExtent::None,
        }
    }

    // 5. ADVANCED_FACE -> Face
    for (&id, ent) in &raw_entities {
        if ent.name == "ADVANCED_FACE" {
            if ent.args.len() >= 4 {
                let bound_refs = parse_ref_list(&ent.args[1])?;
                let surface_ref = parse_ref(ent.args[2].trim().trim_start_matches('(').trim_end_matches(')'))?;
                let sense_str = ent.args[3].trim();
                let sense = if sense_str == ".F." { FaceNormal::Reversed } else { FaceNormal::Same };

                let geom = if let Some(plane) = carrier_planes.get(&surface_ref) {
                    FaceGeom::Plane(*plane)
                } else if let Some(cyl) = carrier_cylinders.get(&surface_ref) {
                    FaceGeom::Cylinder(*cyl)
                } else if let Some(sphere) = carrier_spheres.get(&surface_ref) {
                    FaceGeom::Sphere(*sphere)
                } else if let Some(torus) = carrier_tori.get(&surface_ref) {
                    FaceGeom::Torus(*torus)
                } else {
                    return Err(format!("unknown carrier surface #{} for face #{}", surface_ref, id));
                };

                let mut outer_loop = LoopId::default();
                let mut inner_loops = Vec::new();
                let mut all_loops = Vec::new();

                for &b_ref in &bound_refs {
                    let &(loop_ref, is_outer) = bound_loops.get(&b_ref)
                        .ok_or_else(|| format!("face boundary #{} not found for face #{}", b_ref, id))?;
                    let loop_id = *loop_map.get(&loop_ref)
                        .ok_or_else(|| format!("loop #{} not found for face #{}", loop_ref, id))?;
                    
                    all_loops.push(loop_id);
                    if is_outer {
                        outer_loop = loop_id;
                    } else {
                        inner_loops.push(loop_id);
                    }
                }

                if outer_loop == LoopId::default() && !all_loops.is_empty() {
                    // Fallback to first loop as outer
                    outer_loop = all_loops[0];
                }

                // Reconstruct exact face boundaries and extent
                let extent = reconstruct_face_extent(&brep, &geom, &all_loops);

                let fid = brep.add_face(Face {
                    geom,
                    normal: sense,
                    outer_loop,
                    inner_loops,
                    shell: ShellId::default(),
                    extent,
                });

                // Update loops with their parent face ID
                for &lp_id in &all_loops {
                    if let Some(lp) = brep.loops.get_mut(lp_id) {
                        lp.face = fid;
                    }
                }
                face_map.insert(id, fid);
            }
        }
    }

    // 6. CLOSED_SHELL -> Shell
    for (&id, ent) in &raw_entities {
        if ent.name == "CLOSED_SHELL" {
            if ent.args.len() >= 2 {
                let face_refs = parse_ref_list(&ent.args[1])?;
                let mut faces = Vec::new();
                for &f_ref in &face_refs {
                    let face = *face_map.get(&f_ref)
                        .ok_or_else(|| format!("face #{} not found for shell #{}", f_ref, id))?;
                    faces.push(face);
                }

                let shell_id = brep.add_shell(Shell {
                    faces: faces.clone(),
                    is_outer: true,
                    solid: SolidId::default(),
                });

                // Link faces to shell
                for &fid in &faces {
                    if let Some(f) = brep.faces.get_mut(fid) {
                        f.shell = shell_id;
                    }
                }
                shell_map.insert(id, shell_id);
            }
        }
    }

    // 7. MANIFOLD_SOLID_BREP -> Solid
    for (&id, ent) in &raw_entities {
        if ent.name == "MANIFOLD_SOLID_BREP" {
            if ent.args.len() >= 2 {
                let shell_ref = parse_ref(&ent.args[1])?;
                let shell_id = *shell_map.get(&shell_ref)
                    .ok_or_else(|| format!("shell #{} not found for solid #{}", shell_ref, id))?;

                let name = ent.args[0].trim().trim_matches('\'').to_string();
                let solid_id = brep.add_solid(Solid {
                    shells: vec![shell_id],
                    name: Some(name),
                });

                // Link shell to solid
                if let Some(sh) = brep.shells.get_mut(shell_id) {
                    sh.solid = solid_id;
                }
            }
        }
    }

    // Reconstruct compound / assembly hierarchy (e.g. from NEXT_ASSEMBLY_USAGE_OCCURRENCE)
    // For now we just return the flat BRep populated with solids.
    brep.fuse_solids(); // Merge everything under a single solid for unified meshing convenience.

    Ok(brep)
}
