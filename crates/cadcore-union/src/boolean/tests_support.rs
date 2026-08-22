//! Shared test fixtures for the `boolean` module: arbitrary primitive solids
//! built as real B-Reps (so they double as the kind of inputs the general
//! union must handle).

#![cfg(test)]

use cadcore_geom::{Circle3, CylSurf, Line3, Plane3};
use cadcore_math::{Point3, UnitVec3, Vec3};
use cadcore_topo::{
    BRep, CoEdge, CoEdgeSense, Edge, EdgeGeom, Face, FaceBoundary, FaceExtent, FaceGeom, FaceId,
    FaceNormal, Loop, Shell, Solid, Vertex,
};

/// Build an axis-aligned box `[min,max]` as a real 6-face B-Rep with shared
/// edges (Polygon extents, outward-wound loops).  Returns the face ids.
pub fn axis_box(brep: &mut BRep, min: Point3, max: Point3) -> Vec<FaceId> {
    let c = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let p = [
        c(min.x, min.y, min.z),
        c(max.x, min.y, min.z),
        c(max.x, max.y, min.z),
        c(min.x, max.y, min.z),
        c(min.x, min.y, max.z),
        c(max.x, min.y, max.z),
        c(max.x, max.y, max.z),
        c(min.x, max.y, max.z),
    ];
    let v: Vec<_> = p.iter().map(|&q| brep.add_vertex(Vertex { point: q })).collect();

    let quads: [([usize; 4], Vec3); 6] = [
        ([0, 3, 2, 1], Vec3::new(0.0, 0.0, -1.0)),
        ([4, 5, 6, 7], Vec3::new(0.0, 0.0, 1.0)),
        ([0, 1, 5, 4], Vec3::new(0.0, -1.0, 0.0)),
        ([3, 7, 6, 2], Vec3::new(0.0, 1.0, 0.0)),
        ([0, 4, 7, 3], Vec3::new(-1.0, 0.0, 0.0)),
        ([1, 2, 6, 5], Vec3::new(1.0, 0.0, 0.0)),
    ];

    let mut faces = Vec::new();
    for (idx, nrm) in quads {
        let pts: Vec<Point3> = idx.iter().map(|&i| p[i]).collect();
        let normal = UnitVec3::try_from_vec(nrm).unwrap();
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
            let a = v[idx[k]];
            let b = v[idx[(k + 1) % 4]];
            let seg = p[idx[(k + 1) % 4]] - p[idx[k]];
            let edge = brep.add_edge(Edge {
                geom: EdgeGeom::Line(Line3::new(p[idx[k]], UnitVec3::try_from_vec(seg).unwrap())),
                v_start: a,
                v_end: b,
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
    let shell = brep.add_shell(Shell {
        faces: faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = brep.add_solid(Solid {
        shells: vec![shell],
        name: Some("box".into()),
    });
    brep.shells[shell].solid = solid;
    for &f in &faces {
        brep.faces[f].shell = shell;
    }
    faces
}

/// An axis-aligned box built the way `cadcore_ops::build_solid_box` builds an
/// electrode plate: six `FaceExtent::Polygon` faces with **placeholder loops**
/// (no vertices, edges or co-edges at all).
///
/// The distinction matters — a solid like this carries its whole boundary in
/// the extent template, so every consumer has to synthesise it.  This is the
/// real electrode input, as opposed to [`axis_box`]'s fully-built topology.
pub fn template_box(brep: &mut BRep, min: Point3, max: Point3) -> Vec<FaceId> {
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let quads: [([Point3; 4], Vec3); 6] = [
        (
            [
                p(min.x, min.y, max.z),
                p(max.x, min.y, max.z),
                p(max.x, max.y, max.z),
                p(min.x, max.y, max.z),
            ],
            Vec3::new(0.0, 0.0, 1.0),
        ),
        (
            [
                p(min.x, max.y, min.z),
                p(max.x, max.y, min.z),
                p(max.x, min.y, min.z),
                p(min.x, min.y, min.z),
            ],
            Vec3::new(0.0, 0.0, -1.0),
        ),
        (
            [
                p(min.x, max.y, max.z),
                p(max.x, max.y, max.z),
                p(max.x, max.y, min.z),
                p(min.x, max.y, min.z),
            ],
            Vec3::new(0.0, 1.0, 0.0),
        ),
        (
            [
                p(min.x, min.y, min.z),
                p(max.x, min.y, min.z),
                p(max.x, min.y, max.z),
                p(min.x, min.y, max.z),
            ],
            Vec3::new(0.0, -1.0, 0.0),
        ),
        (
            [
                p(max.x, min.y, min.z),
                p(max.x, max.y, min.z),
                p(max.x, max.y, max.z),
                p(max.x, min.y, max.z),
            ],
            Vec3::new(1.0, 0.0, 0.0),
        ),
        (
            [
                p(min.x, min.y, max.z),
                p(min.x, max.y, max.z),
                p(min.x, max.y, min.z),
                p(min.x, min.y, min.z),
            ],
            Vec3::new(-1.0, 0.0, 0.0),
        ),
    ];
    let mut faces = Vec::new();
    for (pts, nrm) in quads {
        let normal = UnitVec3::try_from_vec(nrm).unwrap();
        let face = brep.add_face(Face {
            geom: FaceGeom::Plane(Plane3::from_origin_normal(pts[0], normal)),
            normal: FaceNormal::Same,
            outer_loop: Default::default(),
            inner_loops: Vec::new(),
            shell: Default::default(),
            extent: FaceExtent::Polygon {
                points: pts.to_vec(),
            },
        });
        faces.push(face);
    }
    let shell = brep.add_shell(Shell {
        faces: faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = brep.add_solid(Solid {
        shells: vec![shell],
        name: Some("plate".into()),
    });
    brep.shells[shell].solid = solid;
    for &f in &faces {
        brep.faces[f].shell = shell;
    }
    faces
}

/// A capped cylinder solid: lateral cylinder + two disk caps.
pub fn capped_cylinder(
    brep: &mut BRep,
    base: Point3,
    axis: UnitVec3,
    r: f64,
    len: f64,
) -> Vec<FaceId> {
    let surf = CylSurf::new(base, axis, r);
    let top = base + axis.as_vec() * len;
    let lateral = brep.add_face(Face {
        geom: FaceGeom::Cylinder(surf),
        normal: FaceNormal::Same,
        outer_loop: Default::default(),
        inner_loops: Vec::new(),
        shell: Default::default(),
        extent: FaceExtent::Cylinder {
            length: len,
            start: FaceBoundary::Circle(Circle3::new(base, axis, r)),
            end: FaceBoundary::Circle(Circle3::new(top, axis, r)),
        },
    });
    let cap0 = brep.add_face(Face {
        geom: FaceGeom::Plane(Plane3::from_origin_normal(
            base,
            UnitVec3::try_from_vec(axis.as_vec() * -1.0).unwrap(),
        )),
        normal: FaceNormal::Same,
        outer_loop: Default::default(),
        inner_loops: Vec::new(),
        shell: Default::default(),
        extent: FaceExtent::Disk { radius: r },
    });
    let cap1 = brep.add_face(Face {
        geom: FaceGeom::Plane(Plane3::from_origin_normal(top, axis)),
        normal: FaceNormal::Same,
        outer_loop: Default::default(),
        inner_loops: Vec::new(),
        shell: Default::default(),
        extent: FaceExtent::Disk { radius: r },
    });
    let faces = vec![lateral, cap0, cap1];
    let shell = brep.add_shell(Shell {
        faces: faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    for &f in &faces {
        brep.faces[f].shell = shell;
    }
    faces
}
