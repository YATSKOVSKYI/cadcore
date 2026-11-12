use slotmap::Key;
use crate::kernel::{
    math::{Point3, Vec3, UnitVec3, Frame3},
    geom::{Plane3, CylSurf, SphereSurf, TorusSurf},
    topo::{BRep, Face, FaceGeom, FaceNormal, FaceExtent, EdgeGeom, PrismAxis},
};
use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FaceMesh {
    pub face_id: u32,
    pub solid_id: u32,
    pub positions: Vec<f32>,
    pub normals: Vec<f32>,
    pub indices: Vec<u32>,
}

fn to_f32_3(v: Vec3) -> [f32; 3] {
    [v.x as f32, v.y as f32, v.z as f32]
}

fn pt_to_f32_3(p: Point3) -> [f32; 3] {
    [p.x as f32, p.y as f32, p.z as f32]
}

/// Tessellate all faces in a BRep solid and return a flat list of face meshes.
pub fn tessellate_brep(brep: &BRep) -> Vec<FaceMesh> {
    let mut meshes = Vec::new();

    for (fid, face) in &brep.faces {
        let face_id_raw = fid.data().as_ffi() as u32; // Extract numeric ID from slotmap key
        let mut solid_id_raw = 0u32;
        if let Some(shell) = brep.shells.get(face.shell) {
            solid_id_raw = shell.solid.data().as_ffi() as u32;
        }

        let mut positions = Vec::new();
        let mut normals = Vec::new();
        let mut indices = Vec::new();

        // Check face normal sense
        let normal_scale = match face.normal {
            FaceNormal::Same => 1.0f32,
            FaceNormal::Reversed => -1.0f32,
        };

        match &face.extent {
            FaceExtent::Disk { radius } => {
                if let FaceGeom::Plane(plane) = &face.geom {
                    tessellate_disk(plane, *radius, normal_scale, &mut positions, &mut normals, &mut indices);
                }
            }
            FaceExtent::PartialDisk { radius, start_angle, end_angle } => {
                if let FaceGeom::Plane(plane) = &face.geom {
                    tessellate_partial_disk(plane, *radius, *start_angle, *end_angle, normal_scale, &mut positions, &mut normals, &mut indices);
                }
            }
            FaceExtent::Polygon { points } => {
                if let FaceGeom::Plane(plane) = &face.geom {
                    tessellate_polygon(plane, points, normal_scale, &mut positions, &mut normals, &mut indices);
                }
            }
            FaceExtent::Cylinder { length, .. } => {
                if let FaceGeom::Cylinder(cyl) = &face.geom {
                    let z_offset = get_face_min_z(brep, face, cyl);
                    tessellate_cylinder(cyl, *length, z_offset, 0.0, std::f32::consts::TAU as f64, cyl.frame.x, normal_scale, &mut positions, &mut normals, &mut indices);
                }
            }
            FaceExtent::PartialCylinder { length, arc_start_angle, arc_end_angle, arc_ref_dir } => {
                if let FaceGeom::Cylinder(cyl) = &face.geom {
                    let z_offset = get_face_min_z(brep, face, cyl);
                    tessellate_cylinder(cyl, *length, z_offset, *arc_start_angle, *arc_end_angle, *arc_ref_dir, normal_scale, &mut positions, &mut normals, &mut indices);
                }
            }
            FaceExtent::CylinderArcFace { length, arc_start_angle, arc_end_angle, arc_ref_dir } => {
                if let FaceGeom::Cylinder(cyl) = &face.geom {
                    let z_offset = get_face_min_z(brep, face, cyl);
                    tessellate_cylinder(cyl, *length, z_offset, *arc_start_angle, *arc_end_angle, *arc_ref_dir, normal_scale, &mut positions, &mut normals, &mut indices);
                }
            }
            FaceExtent::RoundedRectCap { axis, umin, umax, vmin, vmax, radius, w, plus } => {
                if let FaceGeom::Plane(plane) = &face.geom {
                    tessellate_rounded_rect_cap(plane, *axis, *umin, *umax, *vmin, *vmax, *radius, *w, *plus, normal_scale, &mut positions, &mut normals, &mut indices);
                }
            }
            FaceExtent::TorusFillet { start_circle, end_circle } => {
                if let FaceGeom::Torus(torus) = &face.geom {
                    let c1 = start_circle.frame.origin;
                    let c2 = end_circle.frame.origin;
                    let x1 = torus.frame.x.dot_vec(c1 - torus.frame.origin);
                    let y1 = torus.frame.y.dot_vec(c1 - torus.frame.origin);
                    let theta_start = y1.atan2(x1);

                    let x2 = torus.frame.x.dot_vec(c2 - torus.frame.origin);
                    let y2 = torus.frame.y.dot_vec(c2 - torus.frame.origin);
                    let theta_end = y2.atan2(x2);

                    let mut angle_diff = theta_end - theta_start;
                    while angle_diff > std::f64::consts::PI {
                        angle_diff -= std::f32::consts::TAU as f64;
                    }
                    while angle_diff < -std::f64::consts::PI {
                        angle_diff += std::f32::consts::TAU as f64;
                    }

                    tessellate_torus(torus, theta_start, theta_start + angle_diff, normal_scale, &mut positions, &mut normals, &mut indices);
                }
            }
            _ => {
                // Fallback based on geometry type
                match &face.geom {
                    FaceGeom::Sphere(sphere) => {
                        tessellate_sphere(sphere, normal_scale, &mut positions, &mut normals, &mut indices);
                    }
                    FaceGeom::Plane(plane) => {
                        // Planar fallback: disk of default radius
                        tessellate_disk(plane, 2.0, normal_scale, &mut positions, &mut normals, &mut indices);
                    }
                    FaceGeom::Cylinder(cyl) => {
                        let z_offset = get_face_min_z(brep, face, cyl);
                        tessellate_cylinder(cyl, 10.0, z_offset, 0.0, std::f32::consts::TAU as f64, cyl.frame.x, normal_scale, &mut positions, &mut normals, &mut indices);
                    }
                    FaceGeom::Torus(torus) => {
                        if let Some((start, end)) = get_torus_fillet_angles(brep, face, torus) {
                            tessellate_torus(torus, start, end, normal_scale, &mut positions, &mut normals, &mut indices);
                        } else {
                            println!("Warning: Torus face #{:?} fallback to full torus (origin={:?})", fid, torus.frame.origin);
                            tessellate_torus(torus, 0.0, std::f32::consts::TAU as f64, normal_scale, &mut positions, &mut normals, &mut indices);
                        }
                    }
                }
            }
        }

        if !positions.is_empty() {
            meshes.push(FaceMesh {
                face_id: face_id_raw,
                solid_id: solid_id_raw,
                positions,
                normals,
                indices,
            });
        }
    }

    meshes
}

fn get_face_min_z(brep: &BRep, face: &Face, cyl: &CylSurf) -> f64 {
    let mut min_z = f64::MAX;

    let mut process_loop = |lp_id| {
        if let Some(coedges) = brep.loop_coedges(lp_id) {
            for ce_id in coedges {
                if let Some(ce) = brep.coedges.get(ce_id) {
                    if let Some(edge) = brep.edges.get(ce.edge) {
                        for vid in [edge.v_start, edge.v_end] {
                            if let Some(vtx) = brep.vertices.get(vid) {
                                let z_val = cyl.frame.z.dot_vec(vtx.point - cyl.frame.origin);
                                min_z = min_z.min(z_val);
                            }
                        }
                    }
                }
            }
        }
    };

    process_loop(face.outer_loop);
    for &lp_id in &face.inner_loops {
        process_loop(lp_id);
    }

    if min_z.is_finite() {
        min_z
    } else {
        0.0
    }
}

fn get_torus_fillet_angles(brep: &BRep, face: &Face, torus: &TorusSurf) -> Option<(f64, f64)> {
    let mut centers = Vec::new();

    let mut process_loop = |lp_id| {
        if let Some(coedges) = brep.loop_coedges(lp_id) {
            for ce_id in coedges {
                if let Some(ce) = brep.coedges.get(ce_id) {
                    if let Some(edge) = brep.edges.get(ce.edge) {
                        if let EdgeGeom::Circle(circ) = &edge.geom {
                            let proj = torus.frame.z.dot_vec(circ.frame.origin - torus.frame.origin);
                            let radial_vec = circ.frame.origin - torus.frame.origin - torus.frame.z.as_vec() * proj;
                            let dist_to_axis = radial_vec.length();
                            let z_dist = proj.abs();

                            let radius_ok = (circ.radius - torus.minor_radius).abs() < 1.0e-2;
                            let dist_to_axis_ok = (dist_to_axis - torus.major_radius).abs() < 1.0e-2;
                            let z_dist_ok = z_dist < 1.0e-2;

                            if radius_ok && dist_to_axis_ok && z_dist_ok {
                                centers.push(circ.frame.origin);
                            }
                        }
                    }
                }
            }
        }
    };

    process_loop(face.outer_loop);
    for &lp_id in &face.inner_loops {
        process_loop(lp_id);
    }

    // Deduplicate centers that are very close to each other
    let mut unique_centers: Vec<Point3> = Vec::new();
    for c in centers {
        if !unique_centers.iter().any(|&uc| (uc - c).length() < 1.0e-2) {
            unique_centers.push(c);
        }
    }

    if unique_centers.len() >= 2 {
        let c1 = unique_centers[0];
        let c2 = unique_centers[1];

        let x1 = torus.frame.x.dot_vec(c1 - torus.frame.origin);
        let y1 = torus.frame.y.dot_vec(c1 - torus.frame.origin);
        let theta_start = y1.atan2(x1);

        let x2 = torus.frame.x.dot_vec(c2 - torus.frame.origin);
        let y2 = torus.frame.y.dot_vec(c2 - torus.frame.origin);
        let theta_end = y2.atan2(x2);

        let mut angle_diff = theta_end - theta_start;
        while angle_diff > std::f64::consts::PI {
            angle_diff -= std::f32::consts::TAU as f64;
        }
        while angle_diff < -std::f64::consts::PI {
            angle_diff += std::f32::consts::TAU as f64;
        }

        Some((theta_start, theta_start + angle_diff))
    } else {
        None
    }
}

fn add_vertex(
    pos: [f32; 3],
    norm: [f32; 3],
    positions: &mut Vec<f32>,
    normals: &mut Vec<f32>,
) -> u32 {
    let idx = (positions.len() / 3) as u32;
    positions.extend_from_slice(&pos);
    normals.extend_from_slice(&norm);
    idx
}

fn tessellate_disk(
    plane: &Plane3,
    radius: f64,
    n_scale: f32,
    positions: &mut Vec<f32>,
    normals: &mut Vec<f32>,
    indices: &mut Vec<u32>,
) {
    let n = to_f32_3(plane.frame.z.as_vec() * n_scale as f64);
    let center = pt_to_f32_3(plane.frame.origin);
    let c_idx = add_vertex(center, n, positions, normals);

    let segments = 32usize;
    let mut rim_indices = Vec::with_capacity(segments);

    for i in 0..segments {
        let theta = (std::f32::consts::TAU * i as f32 / segments as f32) as f64;
        let pt = plane.frame.origin + plane.frame.x.as_vec() * (radius * theta.cos()) + plane.frame.y.as_vec() * (radius * theta.sin());
        let r_idx = add_vertex(pt_to_f32_3(pt), n, positions, normals);
        rim_indices.push(r_idx);
    }

    for i in 0..segments {
        let i1 = (i + 1) % segments;
        // CCW winding: center -> rim[i] -> rim[i1]
        indices.push(c_idx);
        if n_scale >= 0.0 {
            indices.push(rim_indices[i]);
            indices.push(rim_indices[i1]);
        } else {
            indices.push(rim_indices[i1]);
            indices.push(rim_indices[i]);
        }
    }
}

fn tessellate_partial_disk(
    plane: &Plane3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
    n_scale: f32,
    positions: &mut Vec<f32>,
    normals: &mut Vec<f32>,
    indices: &mut Vec<u32>,
) {
    let n = to_f32_3(plane.frame.z.as_vec() * n_scale as f64);
    let center = pt_to_f32_3(plane.frame.origin);
    let c_idx = add_vertex(center, n, positions, normals);

    let segments = 16usize;
    let angle_range = end_angle - start_angle;
    let mut rim_indices = Vec::with_capacity(segments + 1);

    for i in 0..=segments {
        let theta = start_angle + (i as f64 / segments as f64) * angle_range;
        let pt = plane.frame.origin + plane.frame.x.as_vec() * (radius * theta.cos()) + plane.frame.y.as_vec() * (radius * theta.sin());
        let r_idx = add_vertex(pt_to_f32_3(pt), n, positions, normals);
        rim_indices.push(r_idx);
    }

    for i in 0..segments {
        indices.push(c_idx);
        if n_scale >= 0.0 {
            indices.push(rim_indices[i]);
            indices.push(rim_indices[i + 1]);
        } else {
            indices.push(rim_indices[i + 1]);
            indices.push(rim_indices[i]);
        }
    }
}

fn tessellate_polygon(
    plane: &Plane3,
    points: &[Point3],
    n_scale: f32,
    positions: &mut Vec<f32>,
    normals: &mut Vec<f32>,
    indices: &mut Vec<u32>,
) {
    if points.len() < 3 {
        return;
    }
    let n = to_f32_3(plane.frame.z.as_vec() * n_scale as f64);
    
    let mut v_indices = Vec::with_capacity(points.len());
    for &pt in points {
        let idx = add_vertex(pt_to_f32_3(pt), n, positions, normals);
        v_indices.push(idx);
    }

    // Convex fan triangulation
    for i in 1..points.len() - 1 {
        indices.push(v_indices[0]);
        if n_scale >= 0.0 {
            indices.push(v_indices[i]);
            indices.push(v_indices[i + 1]);
        } else {
            indices.push(v_indices[i + 1]);
            indices.push(v_indices[i]);
        }
    }
}

/// Cap face of a rounded plate (electrode).  The outline is built directly in
/// world space from the prism's own `(u, v, w)` basis, so it is exact for X-,
/// Y- and Z-normal plates alike, and the fillets are drawn as real arcs instead
/// of being squared off.
#[allow(clippy::too_many_arguments)]
fn tessellate_rounded_rect_cap(
    plane: &Plane3,
    axis: PrismAxis,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    radius: f64,
    w: f64,
    plus: bool,
    n_scale: f32,
    positions: &mut Vec<f32>,
    normals: &mut Vec<f32>,
    indices: &mut Vec<u32>,
) {
    let r = radius.max(0.0).min(((umax - umin) * 0.5).min((vmax - vmin) * 0.5));
    let p = |u: f64, v: f64| axis.point(u, v, w);

    let mut pts: Vec<Point3> = Vec::new();
    if r <= 1.0e-9 {
        pts.extend([p(umin, vmin), p(umax, vmin), p(umax, vmax), p(umin, vmax)]);
    } else {
        // Corner centres in CCW order with the arc's start angle, matching the
        // quarter-cylinder fillets of the solid.
        const SEGS: usize = 6;
        let corners = [
            (umax - r, vmin + r, -std::f64::consts::FRAC_PI_2), // bottom-right
            (umax - r, vmax - r, 0.0),                          // top-right
            (umin + r, vmax - r, std::f64::consts::FRAC_PI_2),  // top-left
            (umin + r, vmin + r, std::f64::consts::PI),         // bottom-left
        ];
        for (cu, cv, a0) in corners {
            for i in 0..=SEGS {
                let a = a0 + (i as f64 / SEGS as f64) * std::f64::consts::FRAC_PI_2;
                pts.push(p(cu + r * a.cos(), cv + r * a.sin()));
            }
        }
    }

    // The profile above is CCW in the uv-plane as seen from −axis; the +axis cap
    // needs the opposite traversal so its triangles keep facing outward.
    if plus {
        pts.reverse();
    }
    tessellate_polygon(plane, &pts, n_scale, positions, normals, indices);
}

fn tessellate_cylinder(
    cyl: &CylSurf,
    length: f64,
    z_offset: f64,
    start_angle: f64,
    end_angle: f64,
    arc_ref_dir: UnitVec3,
    n_scale: f32,
    positions: &mut Vec<f32>,
    normals: &mut Vec<f32>,
    indices: &mut Vec<u32>,
) {
    let right = UnitVec3::try_from_vec(cyl.frame.z.as_vec().cross(arc_ref_dir.as_vec())).unwrap_or(cyl.frame.y);

    let segments = 24usize;
    let angle_range = end_angle - start_angle;
    let is_closed = angle_range.abs() >= std::f64::consts::TAU - 1.0e-5;

    let steps = if is_closed { segments } else { segments + 1 };
    
    // Generate rings at z = z_offset and z = z_offset + length
    let mut ring0 = Vec::with_capacity(steps);
    let mut ring1 = Vec::with_capacity(steps);

    for i in 0..steps {
        let theta = start_angle + (i as f64 / segments as f64) * angle_range;
        let normal_dir = arc_ref_dir.as_vec() * theta.cos() + right.as_vec() * theta.sin();
        let normal_dir_u = UnitVec3::try_from_vec(normal_dir).unwrap_or(cyl.frame.x);

        let n = to_f32_3(normal_dir_u.as_vec() * n_scale as f64);
        
        let p0 = cyl.frame.origin + normal_dir_u.as_vec() * cyl.radius + cyl.frame.z.as_vec() * z_offset;
        let p1 = p0 + cyl.frame.z.as_vec() * length;

        let idx0 = add_vertex(pt_to_f32_3(p0), n, positions, normals);
        let idx1 = add_vertex(pt_to_f32_3(p1), n, positions, normals);

        ring0.push(idx0);
        ring1.push(idx1);
    }

    // Connect rings with triangles (CCW winding order)
    for i in 0..segments {
        let i1 = (i + 1) % steps;
        if !is_closed && i1 == 0 {
            break;
        }

        let a0 = ring0[i];
        let a1 = ring0[i1];
        let b0 = ring1[i];
        let b1 = ring1[i1];

        // Quad: (a0, a1, b1, b0)
        if n_scale >= 0.0 {
            indices.extend([a0, a1, b0,   a1, b1, b0]);
        } else {
            indices.extend([a0, b0, a1,   a1, b0, b1]);
        }
    }
}

fn tessellate_sphere(
    sphere: &SphereSurf,
    n_scale: f32,
    positions: &mut Vec<f32>,
    normals: &mut Vec<f32>,
    indices: &mut Vec<u32>,
) {
    let lat_segments = 12usize;
    let lon_segments = 24usize;

    // Grid of vertices
    let mut grid = Vec::new();

    for j in 0..=lat_segments {
        let phi = -std::f32::consts::FRAC_PI_2 + (j as f32 / lat_segments as f32) * std::f32::consts::PI;
        let (sin_phi, cos_phi) = phi.sin_cos();

        let mut row = Vec::new();
        for i in 0..=lon_segments {
            let theta = (i as f32 / lon_segments as f32) * std::f32::consts::TAU;
            let (sin_theta, cos_theta) = theta.sin_cos();

            let nx = (cos_phi * cos_theta) as f64;
            let ny = (cos_phi * sin_theta) as f64;
            let nz = sin_phi as f64;

            let n_vec = Vec3::new(nx, ny, nz);
            let n = to_f32_3(n_vec * n_scale as f64);
            let pt = sphere.centre + n_vec * sphere.radius;

            let idx = add_vertex(pt_to_f32_3(pt), n, positions, normals);
            row.push(idx);
        }
        grid.push(row);
    }

    // Triangulate grid
    for j in 0..lat_segments {
        for i in 0..lon_segments {
            let a0 = grid[j][i];
            let a1 = grid[j][i + 1];
            let b0 = grid[j + 1][i];
            let b1 = grid[j + 1][i + 1];

            if n_scale >= 0.0 {
                indices.extend([a0, a1, b0,   a1, b1, b0]);
            } else {
                indices.extend([a0, b0, a1,   a1, b0, b1]);
            }
        }
    }
}

fn tessellate_torus(
    torus: &TorusSurf,
    start_theta: f64,
    end_theta: f64,
    n_scale: f32,
    positions: &mut Vec<f32>,
    normals: &mut Vec<f32>,
    indices: &mut Vec<u32>,
) {
    // torus is represented as a grid
    let major_segments = 24usize;
    let minor_segments = 12usize;

    let mut grid = Vec::new();
    let theta_range = end_theta - start_theta;

    for j in 0..=major_segments {
        let theta = start_theta + (j as f32 / major_segments as f32) as f64 * theta_range;
        let (sin_theta, cos_theta) = theta.sin_cos();

        let mut row = Vec::new();
        for i in 0..=minor_segments {
            let phi = (i as f32 / minor_segments as f32) * std::f32::consts::TAU;
            let (sin_phi, cos_phi) = phi.sin_cos();

            // Point on major circle
            let major_pt = torus.frame.x.as_vec() * (torus.major_radius * cos_theta as f64) + torus.frame.y.as_vec() * (torus.major_radius * sin_theta as f64);
            
            // Radial direction vector of the major circle
            let radial_dir = torus.frame.x.as_vec() * cos_theta as f64 + torus.frame.y.as_vec() * sin_theta as f64;
            
            // Normal direction of torus surface
            let norm_dir = radial_dir * cos_phi as f64 + torus.frame.z.as_vec() * sin_phi as f64;
            let norm_dir_u = UnitVec3::try_from_vec(norm_dir).unwrap_or(torus.frame.z);

            let n = to_f32_3(norm_dir_u.as_vec() * n_scale as f64);
            let pt = torus.frame.origin + major_pt + norm_dir_u.as_vec() * torus.minor_radius;

            let idx = add_vertex(pt_to_f32_3(pt), n, positions, normals);
            row.push(idx);
        }
        grid.push(row);
    }

    // Triangulate grid
    for j in 0..major_segments {
        for i in 0..minor_segments {
            let a0 = grid[j][i];
            let a1 = grid[j][i + 1];
            let b0 = grid[j + 1][i];
            let b1 = grid[j + 1][i + 1];

            if n_scale >= 0.0 {
                indices.extend([a0, a1, b0,   a1, b1, b0]);
            } else {
                indices.extend([a0, b0, a1,   a1, b0, b1]);
            }
        }
    }
}
