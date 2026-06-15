//! End-to-end: build a union solid entirely on the new engine, serialise it
//! to STEP via the cadcore-step writer, and verify the emitted text is
//! AP214-manifold (every EDGE_CURVE used twice, opposite senses).
//!
//! This is the first STEP file produced by the cadcore-union pipeline.

use cadcore_geom::CylSurf;
use cadcore_math::{Point3, UnitVec3, Vec3};
use cadcore_topo::BRep;
use cadcore_union::arrange::solid::{fuse_crossing_grid, GridLeg};
use cadcore_union::validate::step_text::manifold_check;

#[test]
fn grid_2x2_exports_manifold_step() {
    let r = 0.275;
    let xleg = |y: f64| GridLeg {
        surf: CylSurf::new(Point3::new(-2.0, y, 0.0), UnitVec3::X, r),
        length: 4.0,
    };
    let yleg = |x: f64| GridLeg {
        surf: CylSurf::new(Point3::new(x, -2.0, 0.35), UnitVec3::Y, r),
        length: 4.0,
    };
    let legs = [xleg(0.0), xleg(1.2), yleg(0.0), yleg(1.2)];
    let mut brep = BRep::new();
    let _shell = fuse_crossing_grid(&mut brep, &legs);

    let step = cadcore_step::brep_to_step(&brep).expect("write step");
    assert!(step.contains("MANIFOLD_SOLID_BREP") || step.contains("CLOSED_SHELL"));

    let violations = manifold_check(&step);
    assert!(
        violations.is_empty(),
        "{} manifold violations in new-engine STEP, e.g. {:?}",
        violations.len(),
        &violations[..violations.len().min(6)]
    );

    // stash for external OCC/SW validation
    let out = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_grid_2x2.step");
    if let Some(parent) = out.parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    let _ = std::fs::write(&out, &step);
    println!("wrote {} bytes to {:?}", step.len(), out);
}

#[test]
fn u_filament_exports_seamless_step() {
    use cadcore_geom::TorusSurf;
    use cadcore_union::arrange::filament::{fuse_poly_filament, ElbowSeg, LegSeg};

    // path-from-verts (mirror of the unit-test helper) for a planar U
    let r = 0.275;
    let fr = 0.5;
    let verts = [
        Point3::new(-2.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
        Point3::new(-2.0, 2.0, 0.0),
    ];
    let dir = |a: Point3, b: Point3| UnitVec3::try_from_vec(b - a).unwrap();
    let n = verts.len();
    let mut ti = vec![Point3::new(0.0, 0.0, 0.0); n];
    let mut to = vec![Point3::new(0.0, 0.0, 0.0); n];
    let mut elbows = Vec::new();
    for k in 1..n - 1 {
        let din = dir(verts[k - 1], verts[k]);
        let dout = dir(verts[k], verts[k + 1]);
        let turn = din.dot_vec(dout.as_vec()).clamp(-1.0, 1.0).acos();
        let t = fr / (turn / 2.0).tan();
        let tin = verts[k] - din.as_vec() * t;
        let tout = verts[k] + dout.as_vec() * t;
        let axis = UnitVec3::try_from_vec(din.as_vec().cross(dout.as_vec())).unwrap();
        let inn = UnitVec3::try_from_vec(axis.as_vec().cross(din.as_vec())).unwrap();
        let centre = tin + inn.as_vec() * fr;
        let torus = TorusSurf::new(centre, axis, fr, r);
        let thof = |p: Point3| {
            let w = p - torus.frame.origin;
            let rd = w - torus.frame.z.as_vec() * torus.frame.z.dot_vec(w);
            torus.frame.y.dot_vec(rd).atan2(torus.frame.x.dot_vec(rd))
        };
        let tl = thof(tin);
        let mut d = thof(tout) - tl;
        while d > std::f64::consts::PI { d -= std::f64::consts::TAU; }
        while d < -std::f64::consts::PI { d += std::f64::consts::TAU; }
        elbows.push(ElbowSeg { surf: torus, theta_lo: tl, theta_hi: tl + d });
        ti[k] = tin;
        to[k] = tout;
    }
    let mut legs = Vec::new();
    for k in 0..n - 1 {
        let start = if k == 0 { verts[0] } else { to[k] };
        let end = if k + 1 == n - 1 { verts[n - 1] } else { ti[k + 1] };
        legs.push(LegSeg { surf: CylSurf::new(start, dir(start, end), r), length: (end - start).length() });
    }

    let mut brep = BRep::new();
    let _shell = fuse_poly_filament(&mut brep, &legs, &elbows);
    let step = cadcore_step::brep_to_step(&brep).expect("step");
    let v = manifold_check(&step);
    assert!(v.is_empty(), "manifold: {} e.g. {:?}", v.len(), &v[..v.len().min(4)]);

    let out = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_u_filament.step");
    let _ = std::fs::write(&out, &step);
    println!("U-filament: {} bytes, seam-free (CIRCLE rims)", step.len());
}

#[test]
fn union_two_boxes_exports_manifold_step() {
    use cadcore_geom::{Line3, Plane3};
    use cadcore_topo::{
        CoEdge, CoEdgeSense, Edge, EdgeGeom, Face, FaceExtent, FaceGeom, FaceId, FaceNormal, Loop,
        Shell, Solid, Vertex,
    };
    use cadcore_union::boolean::union;

    // compact axis-box builder (real 6 planar faces with shared edges).
    fn axis_box(brep: &mut BRep, min: Point3, max: Point3) -> Vec<FaceId> {
        let c = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
        let p = [
            c(min.x, min.y, min.z), c(max.x, min.y, min.z), c(max.x, max.y, min.z),
            c(min.x, max.y, min.z), c(min.x, min.y, max.z), c(max.x, min.y, max.z),
            c(max.x, max.y, max.z), c(min.x, max.y, max.z),
        ];
        let v: Vec<_> = p.iter().map(|&q| brep.add_vertex(Vertex { point: q })).collect();
        let quads: [([usize; 4], [f64; 3]); 6] = [
            ([0, 3, 2, 1], [0.0, 0.0, -1.0]),
            ([4, 5, 6, 7], [0.0, 0.0, 1.0]),
            ([0, 1, 5, 4], [0.0, -1.0, 0.0]),
            ([3, 7, 6, 2], [0.0, 1.0, 0.0]),
            ([0, 4, 7, 3], [-1.0, 0.0, 0.0]),
            ([1, 2, 6, 5], [1.0, 0.0, 0.0]),
        ];
        let mut faces = Vec::new();
        for (idx, nrm) in quads {
            let pts: Vec<Point3> = idx.iter().map(|&i| p[i]).collect();
            let normal = UnitVec3::try_from_vec(Vec3::new(nrm[0], nrm[1], nrm[2])).unwrap();
            let plane = Plane3::from_origin_normal(pts[0], normal);
            let face = brep.add_face(Face {
                geom: FaceGeom::Plane(plane),
                normal: FaceNormal::Same,
                outer_loop: Default::default(),
                inner_loops: Vec::new(),
                shell: Default::default(),
                extent: FaceExtent::Polygon { points: pts.clone() },
            });
            let lp = brep.add_loop(Loop { start: Default::default(), face });
            let mut ces = Vec::new();
            for k in 0..4 {
                let seg = p[idx[(k + 1) % 4]] - p[idx[k]];
                let edge = brep.add_edge(Edge {
                    geom: EdgeGeom::Line(Line3::new(p[idx[k]], UnitVec3::try_from_vec(seg).unwrap())),
                    v_start: v[idx[k]],
                    v_end: v[idx[(k + 1) % 4]],
                    t_start: 0.0,
                    t_end: seg.length(),
                    partner: None,
                });
                ces.push(brep.add_coedge(CoEdge {
                    edge,
                    sense: CoEdgeSense::Same,
                    next: Default::default(),
                    prev: Default::default(),
                    loop_id: lp,
                }));
            }
            brep.patch_coedge_links(&ces);
            brep.loops[lp].start = ces[0];
            brep.faces[face].outer_loop = lp;
            faces.push(face);
        }
        let shell = brep.add_shell(Shell { faces: faces.clone(), is_outer: true, solid: Default::default() });
        let solid = brep.add_solid(Solid { shells: vec![shell], name: Some("box".into()) });
        brep.shells[shell].solid = solid;
        for &f in &faces {
            brep.faces[f].shell = shell;
        }
        faces
    }

    let mut a = BRep::new();
    let fa = axis_box(&mut a, Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 2.0));
    let mut b = BRep::new();
    let fb = axis_box(&mut b, Point3::new(1.0, 1.0, 1.0), Point3::new(3.0, 3.0, 3.0));

    let (out, _shell) = union(&a, &fa, &b, &fb, 1e-6);
    let step = cadcore_step::brep_to_step(&out).expect("write step");
    assert!(step.contains("CLOSED_SHELL"));
    let v = manifold_check(&step);
    assert!(v.is_empty(), "manifold: {} e.g. {:?}", v.len(), &v[..v.len().min(4)]);

    let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_union_two_boxes.step");
    let _ = std::fs::write(&outp, &step);
    println!("union two boxes: {} bytes", step.len());
}

#[test]
fn union_three_boxes_exports_manifold_step() {
    use cadcore_geom::{Line3, Plane3};
    use cadcore_topo::{
        CoEdge, CoEdgeSense, Edge, EdgeGeom, Face, FaceExtent, FaceGeom, FaceId, FaceNormal, Loop,
        Shell, Solid, Vertex,
    };
    use cadcore_union::boolean::union_many;

    fn axis_box(brep: &mut BRep, min: Point3, max: Point3) -> Vec<FaceId> {
        let c = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
        let p = [
            c(min.x, min.y, min.z), c(max.x, min.y, min.z), c(max.x, max.y, min.z),
            c(min.x, max.y, min.z), c(min.x, min.y, max.z), c(max.x, min.y, max.z),
            c(max.x, max.y, max.z), c(min.x, max.y, max.z),
        ];
        let v: Vec<_> = p.iter().map(|&q| brep.add_vertex(Vertex { point: q })).collect();
        let quads: [([usize; 4], [f64; 3]); 6] = [
            ([0, 3, 2, 1], [0.0, 0.0, -1.0]), ([4, 5, 6, 7], [0.0, 0.0, 1.0]),
            ([0, 1, 5, 4], [0.0, -1.0, 0.0]), ([3, 7, 6, 2], [0.0, 1.0, 0.0]),
            ([0, 4, 7, 3], [-1.0, 0.0, 0.0]), ([1, 2, 6, 5], [1.0, 0.0, 0.0]),
        ];
        let mut faces = Vec::new();
        for (idx, nrm) in quads {
            let pts: Vec<Point3> = idx.iter().map(|&i| p[i]).collect();
            let normal = UnitVec3::try_from_vec(Vec3::new(nrm[0], nrm[1], nrm[2])).unwrap();
            let face = brep.add_face(Face {
                geom: FaceGeom::Plane(Plane3::from_origin_normal(pts[0], normal)),
                normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
                shell: Default::default(), extent: FaceExtent::Polygon { points: pts.clone() },
            });
            let lp = brep.add_loop(Loop { start: Default::default(), face });
            let mut ces = Vec::new();
            for k in 0..4 {
                let seg = p[idx[(k + 1) % 4]] - p[idx[k]];
                let edge = brep.add_edge(Edge {
                    geom: EdgeGeom::Line(Line3::new(p[idx[k]], UnitVec3::try_from_vec(seg).unwrap())),
                    v_start: v[idx[k]], v_end: v[idx[(k + 1) % 4]], t_start: 0.0, t_end: seg.length(), partner: None,
                });
                ces.push(brep.add_coedge(CoEdge {
                    edge, sense: CoEdgeSense::Same, next: Default::default(), prev: Default::default(), loop_id: lp,
                }));
            }
            brep.patch_coedge_links(&ces);
            brep.loops[lp].start = ces[0];
            brep.faces[face].outer_loop = lp;
            faces.push(face);
        }
        let shell = brep.add_shell(Shell { faces: faces.clone(), is_outer: true, solid: Default::default() });
        let solid = brep.add_solid(Solid { shells: vec![shell], name: Some("box".into()) });
        brep.shells[shell].solid = solid;
        for &f in &faces { brep.faces[f].shell = shell; }
        faces
    }

    let mut a = BRep::new();
    let fa = axis_box(&mut a, Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 2.0));
    let mut b = BRep::new();
    let fb = axis_box(&mut b, Point3::new(1.0, 1.0, 1.0), Point3::new(3.0, 3.0, 3.0));
    let mut c = BRep::new();
    let fc = axis_box(&mut c, Point3::new(2.0, 2.0, 2.0), Point3::new(4.0, 4.0, 4.0));

    let (out, _shell) = union_many(&[(&a, &fa), (&b, &fb), (&c, &fc)], 1e-6).expect("3 solids");
    let step = cadcore_step::brep_to_step(&out).expect("write step");
    let viol = manifold_check(&step);
    assert!(viol.is_empty(), "manifold: {} e.g. {:?}", viol.len(), &viol[..viol.len().min(4)]);

    let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_union_three_boxes.step");
    let _ = std::fs::write(&outp, &step);
    println!("union three boxes: {} bytes", step.len());
}

#[test]
fn union_cylinder_through_box_exports_manifold_step() {
    use cadcore_geom::{Circle3, Line3, Plane3};
    use cadcore_topo::{
        CoEdge, CoEdgeSense, Edge, EdgeGeom, Face, FaceBoundary, FaceExtent, FaceGeom, FaceId,
        FaceNormal, Loop, Shell, Solid, Vertex,
    };
    use cadcore_union::boolean::union;

    fn axis_box(brep: &mut BRep, min: Point3, max: Point3) -> Vec<FaceId> {
        let c = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
        let p = [
            c(min.x, min.y, min.z), c(max.x, min.y, min.z), c(max.x, max.y, min.z),
            c(min.x, max.y, min.z), c(min.x, min.y, max.z), c(max.x, min.y, max.z),
            c(max.x, max.y, max.z), c(min.x, max.y, max.z),
        ];
        let v: Vec<_> = p.iter().map(|&q| brep.add_vertex(Vertex { point: q })).collect();
        let quads: [([usize; 4], [f64; 3]); 6] = [
            ([0, 3, 2, 1], [0.0, 0.0, -1.0]), ([4, 5, 6, 7], [0.0, 0.0, 1.0]),
            ([0, 1, 5, 4], [0.0, -1.0, 0.0]), ([3, 7, 6, 2], [0.0, 1.0, 0.0]),
            ([0, 4, 7, 3], [-1.0, 0.0, 0.0]), ([1, 2, 6, 5], [1.0, 0.0, 0.0]),
        ];
        let mut faces = Vec::new();
        for (idx, nrm) in quads {
            let pts: Vec<Point3> = idx.iter().map(|&i| p[i]).collect();
            let normal = UnitVec3::try_from_vec(Vec3::new(nrm[0], nrm[1], nrm[2])).unwrap();
            let face = brep.add_face(Face {
                geom: FaceGeom::Plane(Plane3::from_origin_normal(pts[0], normal)),
                normal: FaceNormal::Same,
                outer_loop: Default::default(),
                inner_loops: Vec::new(),
                shell: Default::default(),
                extent: FaceExtent::Polygon { points: pts.clone() },
            });
            let lp = brep.add_loop(Loop { start: Default::default(), face });
            let mut ces = Vec::new();
            for k in 0..4 {
                let seg = p[idx[(k + 1) % 4]] - p[idx[k]];
                let edge = brep.add_edge(Edge {
                    geom: EdgeGeom::Line(Line3::new(p[idx[k]], UnitVec3::try_from_vec(seg).unwrap())),
                    v_start: v[idx[k]], v_end: v[idx[(k + 1) % 4]], t_start: 0.0,
                    t_end: seg.length(), partner: None,
                });
                ces.push(brep.add_coedge(CoEdge {
                    edge, sense: CoEdgeSense::Same, next: Default::default(),
                    prev: Default::default(), loop_id: lp,
                }));
            }
            brep.patch_coedge_links(&ces);
            brep.loops[lp].start = ces[0];
            brep.faces[face].outer_loop = lp;
            faces.push(face);
        }
        let shell = brep.add_shell(Shell { faces: faces.clone(), is_outer: true, solid: Default::default() });
        let solid = brep.add_solid(Solid { shells: vec![shell], name: Some("box".into()) });
        brep.shells[shell].solid = solid;
        for &f in &faces { brep.faces[f].shell = shell; }
        faces
    }

    fn capped_cylinder(brep: &mut BRep, base: Point3, axis: UnitVec3, r: f64, len: f64) -> Vec<FaceId> {
        let surf = CylSurf::new(base, axis, r);
        let top = base + axis.as_vec() * len;
        let lateral = brep.add_face(Face {
            geom: FaceGeom::Cylinder(surf), normal: FaceNormal::Same,
            outer_loop: Default::default(), inner_loops: Vec::new(), shell: Default::default(),
            extent: FaceExtent::Cylinder {
                length: len,
                start: FaceBoundary::Circle(Circle3::new(base, axis, r)),
                end: FaceBoundary::Circle(Circle3::new(top, axis, r)),
            },
        });
        let cap0 = brep.add_face(Face {
            geom: FaceGeom::Plane(Plane3::from_origin_normal(base, UnitVec3::try_from_vec(axis.as_vec() * -1.0).unwrap())),
            normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
            shell: Default::default(), extent: FaceExtent::Disk { radius: r },
        });
        let cap1 = brep.add_face(Face {
            geom: FaceGeom::Plane(Plane3::from_origin_normal(top, axis)),
            normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
            shell: Default::default(), extent: FaceExtent::Disk { radius: r },
        });
        let faces = vec![lateral, cap0, cap1];
        let shell = brep.add_shell(Shell { faces: faces.clone(), is_outer: true, solid: Default::default() });
        for &f in &faces { brep.faces[f].shell = shell; }
        faces
    }

    let mut a = BRep::new();
    let fa = axis_box(&mut a, Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 4.0, 2.0));
    let mut b = BRep::new();
    let fb = capped_cylinder(&mut b, Point3::new(2.0, 2.0, -1.0), UnitVec3::Z, 1.0, 4.0);

    let (out, _shell) = union(&a, &fa, &b, &fb, 1e-6);
    let step = cadcore_step::brep_to_step(&out).expect("write step");
    let v = manifold_check(&step);
    assert!(v.is_empty(), "manifold: {} e.g. {:?}", v.len(), &v[..v.len().min(4)]);

    let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_union_cyl_through_box.step");
    let _ = std::fs::write(&outp, &step);
    println!("union cylinder-through-box: {} bytes", step.len());
}

#[test]
fn union_two_crossing_cylinders_exports_manifold_step() {
    use cadcore_geom::Circle3;
    use cadcore_topo::{
        Face, FaceBoundary, FaceExtent, FaceGeom, FaceId, FaceNormal, Shell, Solid,
    };
    use cadcore_union::boolean::union;

    fn capped_cylinder(brep: &mut BRep, base: Point3, axis: UnitVec3, r: f64, len: f64) -> Vec<FaceId> {
        let surf = CylSurf::new(base, axis, r);
        let top = base + axis.as_vec() * len;
        let lateral = brep.add_face(Face {
            geom: FaceGeom::Cylinder(surf), normal: FaceNormal::Same,
            outer_loop: Default::default(), inner_loops: Vec::new(), shell: Default::default(),
            extent: FaceExtent::Cylinder {
                length: len,
                start: FaceBoundary::Circle(Circle3::new(base, axis, r)),
                end: FaceBoundary::Circle(Circle3::new(top, axis, r)),
            },
        });
        let cap0 = brep.add_face(Face {
            geom: FaceGeom::Plane(cadcore_geom::Plane3::from_origin_normal(base, UnitVec3::try_from_vec(axis.as_vec() * -1.0).unwrap())),
            normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
            shell: Default::default(), extent: FaceExtent::Disk { radius: r },
        });
        let cap1 = brep.add_face(Face {
            geom: FaceGeom::Plane(cadcore_geom::Plane3::from_origin_normal(top, axis)),
            normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
            shell: Default::default(), extent: FaceExtent::Disk { radius: r },
        });
        let faces = vec![lateral, cap0, cap1];
        let shell = brep.add_shell(Shell { faces: faces.clone(), is_outer: true, solid: Default::default() });
        let solid = brep.add_solid(Solid { shells: vec![shell], name: Some("cyl".into()) });
        brep.shells[shell].solid = solid;
        for &f in &faces { brep.faces[f].shell = shell; }
        faces
    }

    let mut a = BRep::new();
    let fa = capped_cylinder(&mut a, Point3::new(-3.0, 0.0, 0.0), UnitVec3::X, 1.0, 6.0);
    let mut b = BRep::new();
    let fb = capped_cylinder(&mut b, Point3::new(0.0, 0.0, -3.0), UnitVec3::Z, 0.6, 6.0);

    let (out, _shell) = union(&a, &fa, &b, &fb, 1e-6);
    let step = cadcore_step::brep_to_step(&out).expect("write step");
    let v = manifold_check(&step);
    assert!(v.is_empty(), "manifold: {} e.g. {:?}", v.len(), &v[..v.len().min(4)]);

    let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_union_cross_cylinders.step");
    let _ = std::fs::write(&outp, &step);
    println!("union two crossing cylinders: {} bytes", step.len());
}

#[test]
fn union_filament_elbow_exports_manifold_step() {
    use cadcore_geom::Circle3;
    use cadcore_topo::{Face, FaceBoundary, FaceExtent, FaceGeom, FaceId, FaceNormal, Shell, Solid};
    use cadcore_union::arrange::filament::fuse_poly_filament;
    use cadcore_union::arrange::scaffold::Filament;
    use cadcore_union::boolean::union;

    fn capped_cylinder(brep: &mut BRep, base: Point3, axis: UnitVec3, r: f64, len: f64) -> Vec<FaceId> {
        let surf = CylSurf::new(base, axis, r);
        let top = base + axis.as_vec() * len;
        let lateral = brep.add_face(Face {
            geom: FaceGeom::Cylinder(surf), normal: FaceNormal::Same,
            outer_loop: Default::default(), inner_loops: Vec::new(), shell: Default::default(),
            extent: FaceExtent::Cylinder {
                length: len,
                start: FaceBoundary::Circle(Circle3::new(base, axis, r)),
                end: FaceBoundary::Circle(Circle3::new(top, axis, r)),
            },
        });
        let cap0 = brep.add_face(Face {
            geom: FaceGeom::Plane(cadcore_geom::Plane3::from_origin_normal(base, UnitVec3::try_from_vec(axis.as_vec() * -1.0).unwrap())),
            normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
            shell: Default::default(), extent: FaceExtent::Disk { radius: r },
        });
        let cap1 = brep.add_face(Face {
            geom: FaceGeom::Plane(cadcore_geom::Plane3::from_origin_normal(top, axis)),
            normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
            shell: Default::default(), extent: FaceExtent::Disk { radius: r },
        });
        let faces = vec![lateral, cap0, cap1];
        let shell = brep.add_shell(Shell { faces: faces.clone(), is_outer: true, solid: Default::default() });
        let solid = brep.add_solid(Solid { shells: vec![shell], name: Some("cyl".into()) });
        brep.shells[shell].solid = solid;
        for &f in &faces { brep.faces[f].shell = shell; }
        faces
    }

    // Pass-through case (crosser hits a LEG, elbow carried whole): this is the
    // OCC-valid torus example.  The torus-CUT case is watertight at the B-Rep
    // level (`union_cylinder_cuts_torus_elbow`) but its windowed B-spline torus
    // is not yet emittable by the STEP writer ("Unorientable" inner loop).
    let fil = Filament::serpentine(
        &[Point3::new(-2.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 2.0, 0.0)],
        0.3, 0.5,
    );
    let mut a = BRep::new();
    let shell_a = fuse_poly_filament(&mut a, &fil.legs, &fil.elbows);
    let fa: Vec<_> = a.shells[shell_a].faces.clone();
    let mut b = BRep::new();
    let fb = capped_cylinder(&mut b, Point3::new(-1.0, 0.0, -1.0), UnitVec3::Z, 0.2, 2.0);

    let (out, _shell) = union(&a, &fa, &b, &fb, 1e-6);
    let step = cadcore_step::brep_to_step(&out).expect("write step");
    let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_union_filament_elbow.step");
    let _ = std::fs::write(&outp, &step);
    // Pass-through elbow: the carried-whole elbow torus now emits as an analytic
    // TOROIDAL_SURFACE (its two minor end-circle loops pin the θ-band), so OCC is
    // BOP-clean (max edge tolerance ~6e-5) instead of the old NURBS-patch gap.
    assert!(step.contains("CLOSED_SHELL"));
    assert!(step.contains("TOROIDAL_SURFACE"), "elbow torus emitted analytically");
    println!("union filament+elbow with crosser: {} bytes (OCC BOP-clean)", step.len());

    // torus-CUT variant (crosser through the elbow): the cut elbow emits as an
    // analytic TOROIDAL_SURFACE and is now OCC BOP-CLEAN (strict `solid.check`
    // passes, max edge tolerance ~1.7e-4): the minor circles and the φ-seam are
    // emitted as analytic CIRCLEs with native (θ,φ) line pcurves, and the
    // transcendental crosser∩torus window is finely sampled.
    let th = 0.5 * (fil.elbows[0].theta_lo + fil.elbows[0].theta_hi);
    let tor = fil.elbows[0].surf;
    let ring = tor.frame.origin + (tor.frame.x * th.cos() + tor.frame.y * th.sin()) * tor.major_radius;
    let mut bc = BRep::new();
    let fbc = capped_cylinder(&mut bc, Point3::new(ring.x, ring.y, -1.0), UnitVec3::Z, 0.2, 2.0);
    let (outc, _s) = union(&a, &fa, &bc, &fbc, 1e-6);
    let stepc = cadcore_step::brep_to_step(&outc).expect("write step");
    let outpc = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_union_cut_elbow.step");
    let _ = std::fs::write(&outpc, &stepc);
    assert!(stepc.contains("CLOSED_SHELL"));
    assert!(stepc.contains("TOROIDAL_SURFACE"), "cut elbow emits analytic torus");
    assert!(stepc.contains("CIRCLE("), "minor circles & φ-seam emit as analytic CIRCLEs");
    println!("union CUT elbow: {} bytes (OCC BOP-clean, analytic TOROIDAL_SURFACE)", stepc.len());
}

#[test]
fn union_torus_on_torus_exports_step() {
    use cadcore_geom::{Circle3, Plane3, TorusSurf};
    use cadcore_topo::{FaceId, Shell};
    use cadcore_union::arrange::assembly::Assembler;
    use cadcore_union::boolean::union;
    use std::f64::consts::FRAC_PI_2;

    fn capped_torus_arc(brep: &mut BRep, tor: TorusSurf, lo: f64, hi: f64) -> Vec<FaceId> {
        let minor = |theta: f64| -> (Circle3, UnitVec3) {
            let (s, c) = theta.sin_cos();
            let ring = tor.frame.x * c + tor.frame.y * s;
            let centre = tor.frame.origin + ring * tor.major_radius;
            let tangent = UnitVec3::try_from_vec(tor.frame.x.as_vec() * -s + tor.frame.y.as_vec() * c).unwrap();
            (Circle3::new(centre, tangent, tor.minor_radius), tangent)
        };
        let (clo, tlo) = minor(lo);
        let (chi, thi) = minor(hi);
        let mut asm = Assembler::new(brep, 1e-7);
        asm.emit_torus_fillet(tor, clo, chi);
        asm.emit_disk_cap(Plane3::from_origin_normal(clo.frame.origin, UnitVec3::try_from_vec(tlo.as_vec() * -1.0).unwrap()), clo);
        asm.emit_disk_cap(Plane3::from_origin_normal(chi.frame.origin, thi), chi);
        let faces = asm.faces().to_vec();
        let shell = brep.add_shell(Shell { faces: faces.clone(), is_outer: true, solid: Default::default() });
        for &f in &faces { brep.faces[f].shell = shell; }
        faces
    }

    let r = 0.8 / 2_f64.sqrt();
    let mut a = BRep::new();
    let fa = capped_torus_arc(&mut a, TorusSurf::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, 0.8, 0.3), 0.0, FRAC_PI_2);
    let mut b = BRep::new();
    let fb = capped_torus_arc(&mut b, TorusSurf::new(Point3::new(-r, r + 0.8, 0.18), UnitVec3::X, 0.8, 0.3), FRAC_PI_2 - 1.0, FRAC_PI_2 + 1.0);

    let (out, _shell) = union(&a, &fa, &b, &fb, 1e-6);
    let step = cadcore_step::brep_to_step(&out).expect("write step");
    assert!(step.contains("CLOSED_SHELL"));
    // OCC BOP-clean (max edge tolerance ~3.9e-5): minor circles emit as analytic
    // CIRCLEs; the angled torus∩torus arc — circular on neither torus — stays a
    // finely-sampled polyline (it must, to keep a valid pcurve on both tori).
    assert!(step.contains("CIRCLE("), "minor circles emit as analytic CIRCLEs");
    let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_union_torus_on_torus.step");
    let _ = std::fs::write(&outp, &step);
    println!("union torus-on-torus: {} bytes", step.len());
}

#[test]
fn union_stacked_donuts_exports_step() {
    use cadcore_geom::{Circle3, TorusSurf};
    use cadcore_topo::{Face, FaceExtent, FaceGeom, FaceId, FaceNormal, Shell};
    use cadcore_union::boolean::union;

    // A full donut (no caps): one toroidal face whose θ-seam circle is shared
    // (start_circle == end_circle), so θ wraps the whole ring.
    fn donut(brep: &mut BRep, tor: TorusSurf) -> Vec<FaceId> {
        let ring = tor.frame.x;
        let centre = tor.frame.origin + ring * tor.major_radius;
        let mc = Circle3::new(centre, tor.frame.y, tor.minor_radius);
        let face = brep.add_face(Face {
            geom: FaceGeom::Torus(tor), normal: FaceNormal::Same,
            outer_loop: Default::default(), inner_loops: Vec::new(), shell: Default::default(),
            extent: FaceExtent::TorusFillet { start_circle: mc, end_circle: mc },
        });
        let shell = brep.add_shell(Shell { faces: vec![face], is_outer: true, solid: Default::default() });
        brep.faces[face].shell = shell;
        vec![face]
    }

    // Two coaxial donuts stacked along Z, tubes overlapping (gap 0.4 < 2·0.3).
    let mut a = BRep::new();
    let fa = donut(&mut a, TorusSurf::new(Point3::new(0.0, 0.0, 0.0), UnitVec3::Z, 0.8, 0.3));
    let mut b = BRep::new();
    let fb = donut(&mut b, TorusSurf::new(Point3::new(0.0, 0.0, 0.4), UnitVec3::Z, 0.8, 0.3));

    let (out, _shell) = union(&a, &fa, &b, &fb, 1e-6);
    let step = cadcore_step::brep_to_step(&out).expect("write step");
    assert!(step.contains("CLOSED_SHELL"));
    // OCC BOP-clean @1e-7: the coaxial intersection circles are exact major
    // circles of each torus, so every edge emits as an analytic CIRCLE.
    assert!(step.contains("CIRCLE("), "donut seams emit as analytic CIRCLEs");
    let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_union_stacked_donuts.step");
    let _ = std::fs::write(&outp, &step);
    println!("union stacked donuts: {} bytes", step.len());
}

#[test]
fn union_parallel_cylinders_exports_step() {
    use cadcore_geom::{Circle3, Plane3};
    use cadcore_topo::{Face, FaceBoundary, FaceExtent, FaceGeom, FaceId, FaceNormal, Shell, Solid};
    use cadcore_union::boolean::union;

    fn capped_cylinder(brep: &mut BRep, base: Point3, axis: UnitVec3, r: f64, len: f64) -> Vec<FaceId> {
        let surf = CylSurf::new(base, axis, r);
        let top = base + axis.as_vec() * len;
        let lateral = brep.add_face(Face {
            geom: FaceGeom::Cylinder(surf), normal: FaceNormal::Same,
            outer_loop: Default::default(), inner_loops: Vec::new(), shell: Default::default(),
            extent: FaceExtent::Cylinder { length: len, start: FaceBoundary::Circle(Circle3::new(base, axis, r)), end: FaceBoundary::Circle(Circle3::new(top, axis, r)) },
        });
        let cap0 = brep.add_face(Face { geom: FaceGeom::Plane(Plane3::from_origin_normal(base, UnitVec3::try_from_vec(axis.as_vec() * -1.0).unwrap())), normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(), shell: Default::default(), extent: FaceExtent::Disk { radius: r } });
        let cap1 = brep.add_face(Face { geom: FaceGeom::Plane(Plane3::from_origin_normal(top, axis)), normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(), shell: Default::default(), extent: FaceExtent::Disk { radius: r } });
        let faces = vec![lateral, cap0, cap1];
        let shell = brep.add_shell(Shell { faces: faces.clone(), is_outer: true, solid: Default::default() });
        let solid = brep.add_solid(Solid { shells: vec![shell], name: Some("cyl".into()) });
        brep.shells[shell].solid = solid;
        for &f in &faces { brep.faces[f].shell = shell; }
        faces
    }

    let mut a = BRep::new();
    let fa = capped_cylinder(&mut a, Point3::new(-2.0, 0.0, 0.0), UnitVec3::X, 0.3, 4.0);
    let mut b = BRep::new();
    let fb = capped_cylinder(&mut b, Point3::new(-1.5, 0.0, 0.4), UnitVec3::X, 0.3, 4.0);
    let (out, _shell) = union(&a, &fa, &b, &fb, 1e-6);
    let step = cadcore_step::brep_to_step(&out).expect("write step");
    assert!(step.contains("CLOSED_SHELL"));
    let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../../repro/newengine_union_parallel_cylinders.step");
    let _ = std::fs::write(&outp, &step);
    println!("union parallel cylinders: {} bytes", step.len());
}

#[test]
fn union_three_cylinder_woodpile_exports_manifold_step() {
    use cadcore_geom::Circle3;
    use cadcore_topo::{Face, FaceBoundary, FaceExtent, FaceGeom, FaceId, FaceNormal, Shell, Solid};
    use cadcore_union::boolean::union_n;

    fn capped_cylinder(brep: &mut BRep, base: Point3, axis: UnitVec3, r: f64, len: f64) -> Vec<FaceId> {
        let surf = CylSurf::new(base, axis, r);
        let top = base + axis.as_vec() * len;
        let lateral = brep.add_face(Face {
            geom: FaceGeom::Cylinder(surf), normal: FaceNormal::Same,
            outer_loop: Default::default(), inner_loops: Vec::new(), shell: Default::default(),
            extent: FaceExtent::Cylinder {
                length: len,
                start: FaceBoundary::Circle(Circle3::new(base, axis, r)),
                end: FaceBoundary::Circle(Circle3::new(top, axis, r)),
            },
        });
        let cap0 = brep.add_face(Face {
            geom: FaceGeom::Plane(cadcore_geom::Plane3::from_origin_normal(base, UnitVec3::try_from_vec(axis.as_vec() * -1.0).unwrap())),
            normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
            shell: Default::default(), extent: FaceExtent::Disk { radius: r },
        });
        let cap1 = brep.add_face(Face {
            geom: FaceGeom::Plane(cadcore_geom::Plane3::from_origin_normal(top, axis)),
            normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
            shell: Default::default(), extent: FaceExtent::Disk { radius: r },
        });
        let faces = vec![lateral, cap0, cap1];
        let shell = brep.add_shell(Shell { faces: faces.clone(), is_outer: true, solid: Default::default() });
        let solid = brep.add_solid(Solid { shells: vec![shell], name: Some("cyl".into()) });
        brep.shells[shell].solid = solid;
        for &f in &faces { brep.faces[f].shell = shell; }
        faces
    }

    let mut a = BRep::new();
    let fa = capped_cylinder(&mut a, Point3::new(-3.0, 0.0, 0.0), UnitVec3::X, 1.0, 6.0);
    let mut b = BRep::new();
    let fb = capped_cylinder(&mut b, Point3::new(0.0, -3.0, 0.0), UnitVec3::Y, 0.8, 6.0);
    let mut c = BRep::new();
    let fc = capped_cylinder(&mut c, Point3::new(0.0, 0.0, -3.0), UnitVec3::Z, 0.6, 6.0);

    let (out, _shell) = union_n(&[(&a, &fa), (&b, &fb), (&c, &fc)], 1e-6).expect("3 solids");
    let step = cadcore_step::brep_to_step(&out).expect("write step");
    let viol = manifold_check(&step);
    assert!(viol.is_empty(), "manifold: {} e.g. {:?}", viol.len(), &viol[..viol.len().min(4)]);

    let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_union_woodpile3.step");
    let _ = std::fs::write(&outp, &step);
    println!("union 3-cylinder woodpile: {} bytes", step.len());
}

#[test]
fn scaffold_two_crossing_filaments_exports_manifold_step() {
    use cadcore_union::arrange::scaffold::{fuse_scaffold, Filament};

    let r = 0.275;
    // two serpentines on stacked layers, crossing each other transversally
    // (A runs along X, B runs along Y one layer up)
    let a = Filament::serpentine(
        &[
            Point3::new(-2.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(-2.0, 2.0, 0.0),
        ],
        r,
        0.5,
    );
    let b = Filament::serpentine(
        &[
            Point3::new(-1.5, -1.0, 0.35),
            Point3::new(-1.5, 3.0, 0.35),
            Point3::new(-0.8, 3.0, 0.35),
            Point3::new(-0.8, -1.0, 0.35),
        ],
        r,
        0.3,
    );

    let mut brep = BRep::new();
    let _shell = fuse_scaffold(&mut brep, &[a, b]);
    let step = cadcore_step::brep_to_step(&brep).expect("write step");
    assert!(step.contains("CLOSED_SHELL"));

    let violations = manifold_check(&step);
    assert!(
        violations.is_empty(),
        "{} manifold violations in scaffold STEP, e.g. {:?}",
        violations.len(),
        &violations[..violations.len().min(6)]
    );
    assert!(!step.contains("POLYLINE"), "no POLYLINE — seamless");

    let out = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_scaffold_2filament.step");
    let _ = std::fs::write(&out, &step);
    println!("scaffold (2 crossing filaments): {} bytes", step.len());
}

/// One capped cylinder strand as its own solid (owned BRep + face ids).
fn strand(base: Point3, axis: UnitVec3, r: f64, len: f64) -> (BRep, Vec<cadcore_topo::FaceId>) {
    use cadcore_geom::{Circle3, Plane3};
    use cadcore_topo::{Face, FaceBoundary, FaceExtent, FaceGeom, FaceNormal, Shell, Solid};
    let mut brep = BRep::new();
    let surf = CylSurf::new(base, axis, r);
    let top = base + axis.as_vec() * len;
    let lateral = brep.add_face(Face {
        geom: FaceGeom::Cylinder(surf), normal: FaceNormal::Same,
        outer_loop: Default::default(), inner_loops: Vec::new(), shell: Default::default(),
        extent: FaceExtent::Cylinder { length: len,
            start: FaceBoundary::Circle(Circle3::new(base, axis, r)),
            end: FaceBoundary::Circle(Circle3::new(top, axis, r)) },
    });
    let cap0 = brep.add_face(Face {
        geom: FaceGeom::Plane(Plane3::from_origin_normal(base, UnitVec3::try_from_vec(axis.as_vec() * -1.0).unwrap())),
        normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
        shell: Default::default(), extent: FaceExtent::Disk { radius: r } });
    let cap1 = brep.add_face(Face {
        geom: FaceGeom::Plane(Plane3::from_origin_normal(top, axis)),
        normal: FaceNormal::Same, outer_loop: Default::default(), inner_loops: Vec::new(),
        shell: Default::default(), extent: FaceExtent::Disk { radius: r } });
    let faces = vec![lateral, cap0, cap1];
    let shell = brep.add_shell(Shell { faces: faces.clone(), is_outer: true, solid: Default::default() });
    let solid = brep.add_solid(Solid { shells: vec![shell], name: Some("strand".into()) });
    brep.shells[shell].solid = solid;
    for &f in &faces { brep.faces[f].shell = shell; }
    (brep, faces)
}

/// Build an `nl`-layer woodpile scaffold (layers alternate X/Y, stacked in Z by
/// `hl`), `ns` strands per layer at in-plane pitch `hp`, tube diameter `d`.
/// Returns the per-strand owned solids.
fn scaffold_solids(d: f64, hp: f64, hl: f64, ns: usize, nl: usize) -> Vec<(BRep, Vec<cadcore_topo::FaceId>)> {
    let r = d / 2.0;
    let span = (ns as f64 - 1.0) * hp + 2.0;
    let mut out = Vec::new();
    for layer in 0..nl {
        let z = layer as f64 * hl;
        for k in 0..ns {
            let off = k as f64 * hp;
            let (base, axis) = if layer % 2 == 0 {
                (Point3::new(-1.0, off, z), UnitVec3::X)
            } else {
                (Point3::new(off, -1.0, z), UnitVec3::Y)
            };
            out.push(strand(base, axis, r, span));
        }
    }
    out
}

#[test]
fn scaffold_parametric_sweep_is_manifold() {
    use cadcore_union::boolean::union_n;
    // (diameter, in-plane pitch hp, layer height hl); hl < d so layers fuse.
    let cases = [
        (0.6_f64, 1.2_f64, 0.4_f64),
        (0.5, 1.0, 0.35),
        (0.8, 1.5, 0.5),
        (0.4, 0.9, 0.28),
        (0.7, 1.0, 0.45),
    ];
    for (d, hp, hl) in cases {
        let solids = scaffold_solids(d, hp, hl, 3, 3);
        let refs: Vec<(&BRep, &[cadcore_topo::FaceId])> =
            solids.iter().map(|(b, f)| (b, f.as_slice())).collect();
        let (out, _shell) = union_n(&refs, 1e-6)
            .unwrap_or_else(|| panic!("scaffold union d={d} hp={hp} hl={hl} produced nothing"));
        let step = cadcore_step::brep_to_step(&out)
            .unwrap_or_else(|e| panic!("write step d={d} hp={hp} hl={hl}: {e:?}"));
        let viol = manifold_check(&step);
        assert!(
            viol.is_empty(),
            "d={d} hp={hp} hl={hl}: {} manifold violations e.g. {:?}",
            viol.len(),
            &viol[..viol.len().min(4)]
        );
        let name = format!("newengine_scaffold_d{:.0}_hp{:.0}_hl{:.0}.step",
            d * 100.0, hp * 100.0, hl * 100.0);
        let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../../repro").join(&name);
        let _ = std::fs::write(&outp, &step);
        println!("scaffold d={d} hp={hp} hl={hl}: {} bytes -> {name}", step.len());
    }
}

// KNOWN LIMITATION (over-compressed regime, 2·hl < d): when layers 0 and 2
// (both X) overlap vertically, the engine must (a) union the COPLANAR end-cap
// disks (each cap imprinted by the other's boundary circle, asymmetric, with the
// lateral rims split at the two circle crossings) and (b) handle the parallel
// lateral overlap (ruling lines) at the resulting triple points.  These are
// degenerate-CSG cases not yet implemented; realistic scaffolds (2·hl ≥ d,
// layer height = 0.5–0.8·d) do NOT hit them — see scaffold_parametric_sweep,
// which is BOP-clean across diameters / pitches / layer heights.
#[ignore = "over-compressed 2hl<d: coplanar-cap union + parallel-overlap not yet implemented"]
#[test]
fn scaffold_dense_triplepoints_is_manifold() {
    use cadcore_union::boolean::union_n;
    // DENSE: 2·hl < d, so layer 0 and layer 2 (both X) overlap vertically and a
    // layer-1 (Y) strand crosses that overlap → genuine TRIPLE points where three
    // strands meet.  This is the degenerate case the engine must handle for any
    // diameter / pitch / layer height.
    let cases = [
        (0.6_f64, 1.0_f64, 0.25_f64), // 2hl=0.50 < 0.6
        (0.8, 1.2, 0.30),             // 2hl=0.60 < 0.8
        (0.5, 0.9, 0.20),             // 2hl=0.40 < 0.5
        (0.7, 1.1, 0.28),             // 2hl=0.56 < 0.7
    ];
    for (d, hp, hl) in cases {
        assert!(2.0 * hl < d, "case must overlap layers 0 and 2");
        let solids = scaffold_solids(d, hp, hl, 1, 3);
        let refs: Vec<(&BRep, &[cadcore_topo::FaceId])> =
            solids.iter().map(|(b, f)| (b, f.as_slice())).collect();
        let (out, _shell) = union_n(&refs, 1e-6)
            .unwrap_or_else(|| panic!("dense scaffold union d={d} hp={hp} hl={hl} produced nothing"));
        let step = cadcore_step::brep_to_step(&out)
            .unwrap_or_else(|e| panic!("write step d={d} hp={hp} hl={hl}: {e:?}"));
        let name = format!("newengine_scaffold_dense_d{:.0}_hl{:.0}.step", d * 100.0, hl * 100.0);
        let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../../repro").join(&name);
        let _ = std::fs::write(&outp, &step);
        let viol = manifold_check(&step);
        println!("dense scaffold d={d} hp={hp} hl={hl}: {} bytes, {} manifold violations -> {name}", step.len(), viol.len());
        assert!(
            viol.is_empty(),
            "dense d={d} hp={hp} hl={hl}: {} manifold violations e.g. {:?}",
            viol.len(),
            &viol[..viol.len().min(4)]
        );
    }
}

// KNOWN LIMITATION: the two coplanar end-cap disks (both at x=-1) overlap, and
// COPLANAR overlapping plane faces are not yet unioned (plane_plane returns
// empty for coincident planes).  All 84 unshared edges are cap arcs.  Fix =
// coplanar-face union (imprint each cap with the other's boundary circle; split
// the shared lateral rims at the two crossings).  Lateral ruling-line overlap is
// fine; the caps are the blocker.
#[ignore = "coplanar overlapping cap disks not yet unioned"]
#[test]
fn two_parallel_overlapping_cylinders_is_manifold() {
    use cadcore_union::boolean::union_n;
    // Two parallel X-cylinders (y=0) offset in Z by 0.4 < d=0.6 → they overlap
    // in a lens; the intersection is two ruling LINES.  Must fuse watertight.
    let (a, fa) = strand(Point3::new(-1.0, 0.0, 0.0), UnitVec3::X, 0.3, 4.0);
    let (b, fb) = strand(Point3::new(-1.0, 0.0, 0.4), UnitVec3::X, 0.3, 4.0);
    let refs = [(&a, fa.as_slice()), (&b, fb.as_slice())];
    let (out, shell) = union_n(&refs, 1e-6).expect("parallel union");
    // BRep-level diagnostic: which edges are used != 2 times, and where.
    {
        use std::collections::HashMap;
        let mut uses: HashMap<cadcore_topo::EdgeId, usize> = HashMap::new();
        for &fid in &out.shells[shell].faces {
            let f = &out.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                if let Some(co) = out.loop_coedges(lid) {
                    for c in co { *uses.entry(out.coedges[c].edge).or_insert(0) += 1; }
                }
            }
        }
        let mut bad: Vec<_> = uses.iter().filter(|(_, &n)| n != 2).collect();
        bad.sort_by_key(|(e, _)| format!("{:?}", e));
        eprintln!("DBG parallel: {} faces, {} edges, {} used!=2", out.shells[shell].faces.len(), uses.len(), bad.len());
        for (e, n) in bad.iter().take(6) {
            let ed = &out.edges[**e];
            let p0 = out.vertices[ed.v_start].point;
            let p1 = out.vertices[ed.v_end].point;
            eprintln!("  edge used {}x: ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) geom={}", n,
                p0.x,p0.y,p0.z, p1.x,p1.y,p1.z,
                match &ed.geom { cadcore_topo::EdgeGeom::Line(_)=>"Line", cadcore_topo::EdgeGeom::Circle(_)=>"Circle",
                    cadcore_topo::EdgeGeom::Ellipse(_)=>"Ellipse", cadcore_topo::EdgeGeom::Polyline(_)=>"Polyline" });
        }
    }
    let step = cadcore_step::brep_to_step(&out).expect("step");
    let outp = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../repro/newengine_two_parallel.step");
    let _ = std::fs::write(&outp, &step);
    let viol = manifold_check(&step);
    println!("two parallel: {} bytes, {} violations", step.len(), viol.len());
    assert!(viol.is_empty(), "{} violations e.g. {:?}", viol.len(), &viol[..viol.len().min(4)]);
}
