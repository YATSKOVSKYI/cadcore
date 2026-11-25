//! Geometry catalogue for the release film: every kind of contact the universal fuse handles, as real kernel
//! output, plus a lattice cut by a slanted plane (both halves).  For each case it writes `<name>.step` (the
//! fused / trimmed solid) and `<name>.json` (the raw filament centrelines + radius, for animating the parts
//! before they fuse).
//!
//! cargo run --release -p cadcore-union --example video_catalog -- <out_dir>

use cadcore_geom::Plane3;
use cadcore_math::{Point3, UnitVec3, Vec3};
use cadcore_union::{fuse_and_trim, Filament, FuseOptions};

const R: f64 = 0.5;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn centrelines(fils: &[Filament]) -> String {
    let mut out = Vec::new();
    for f in fils {
        let mut pts = Vec::new();
        for (k, leg) in f.legs.iter().enumerate() {
            let o = leg.surf.frame.origin;
            let a = leg.surf.frame.z.as_vec();
            for q in [o, o + a * leg.length] {
                pts.push(format!("[{:.5},{:.5},{:.5}]", q.x, q.y, q.z));
            }
            if let Some(el) = f.elbows.get(k) {
                let t = &el.surf;
                for s in 1..12 {
                    let th = el.theta_lo + (el.theta_hi - el.theta_lo) * s as f64 / 12.0;
                    let d = t.frame.x.as_vec() * th.cos() + t.frame.y.as_vec() * th.sin();
                    let q = t.frame.origin + d * t.major_radius;
                    pts.push(format!("[{:.5},{:.5},{:.5}]", q.x, q.y, q.z));
                }
            }
        }
        out.push(format!("[{}]", pts.join(",")));
    }
    format!("{{\"r\":{},\"filaments\":[{}]}}", R, out.join(","))
}

fn write(dir: &str, name: &str, fils: &[Filament], planes: &[Plane3]) {
    let t = std::time::Instant::now();
    let o = fuse_and_trim(fils, &FuseOptions::default(), planes, 1.0e-6);
    let step = cadcore_step::brep_to_step(&o.brep).expect("step");
    std::fs::write(format!("{dir}/{name}.step"), step).unwrap();
    std::fs::write(format!("{dir}/{name}.json"), centrelines(fils)).unwrap();
    println!("{name}: {} faces, {:.0} ms", o.brep.faces.len(), t.elapsed().as_secs_f64() * 1e3);
}

fn line(a: Point3, b: Point3) -> Filament {
    Filament::serpentine(&[a, b], R, 0.9)
}

fn main() {
    let dir = std::env::args().nth(1).unwrap_or_else(|| "catalog".into());
    std::fs::create_dir_all(&dir).unwrap();
    let dz = 0.75; // layer pitch < 2R: neighbouring layers interpenetrate (bonded)

    // CROSS: two tubes, one over the other, at right angles
    write(&dir, "cross", &[line(p(-4.0, 0.0, 0.0), p(4.0, 0.0, 0.0)), line(p(0.0, -4.0, dz), p(0.0, 4.0, dz))], &[]);
    // ANGLE: crossing at 35 degrees
    let a = 35f64.to_radians();
    write(&dir, "angle", &[line(p(-4.0, 0.0, 0.0), p(4.0, 0.0, 0.0)),
        line(p(-4.0 * a.cos(), -4.0 * a.sin(), dz), p(4.0 * a.cos(), 4.0 * a.sin(), dz))], &[]);
    // TEE: a tube ends inside the side of another
    write(&dir, "tee", &[line(p(-4.0, 0.0, 0.0), p(4.0, 0.0, 0.0)), line(p(0.0, -4.0, 0.3), p(0.0, 0.05, 0.3))], &[]);
    // ELBOW: a U-turn of one filament laid over a straight one
    write(&dir, "elbow", &[
        Filament::serpentine(&[p(-4.0, -1.2, 0.0), p(1.5, -1.2, 0.0), p(1.5, 1.2, 0.0), p(-4.0, 1.2, 0.0)], R, 1.1),
        line(p(1.2, -4.0, dz), p(1.2, 4.0, dz)),
    ], &[]);
    // SIDE: two parallel tubes fused along their length (bonded walls)
    write(&dir, "side", &[line(p(-4.0, 0.0, 0.0), p(4.0, 0.0, 0.0)), line(p(-3.0, 0.0, dz), p(3.0, 0.0, dz))], &[]);
    // WEAVE: two serpentines crossing each other's legs and elbows
    write(&dir, "weave", &[
        Filament::serpentine(&[p(-4.0, -2.0, 0.0), p(3.0, -2.0, 0.0), p(3.0, 0.0, 0.0), p(-3.0, 0.0, 0.0),
            p(-3.0, 2.0, 0.0), p(4.0, 2.0, 0.0)], R, 0.9),
        Filament::serpentine(&[p(-2.0, -4.0, dz), p(-2.0, 3.0, dz), p(0.0, 3.0, dz), p(0.0, -3.0, dz),
            p(2.0, -3.0, dz), p(2.0, 4.0, dz)], R, 0.9),
    ], &[]);

    // TRIM: a 4-layer woodpile, cut by a slanted plane; both halves (material on the normal side is kept)
    let mut wp = Vec::new();
    for k in 0..4 {
        let z = k as f64 * dz;
        let mut v = Vec::new();
        for i in 0..5 {
            let s = -4.0 + 2.0 * i as f64;
            let (a0, a1) = if i % 2 == 0 { (-5.0, 5.0) } else { (5.0, -5.0) };
            if k % 2 == 0 {
                v.push(p(a0, s, z));
                v.push(p(a1, s, z));
            } else {
                v.push(p(s, a0, z));
                v.push(p(s, a1, z));
            }
        }
        wp.push(Filament::serpentine(&v, R, 0.9));
    }
    let n = UnitVec3::try_from_vec(Vec3::new(1.0, 0.35, 0.25)).unwrap();
    let o = p(0.6, 0.0, 1.0);
    write(&dir, "trim_whole", &wp, &[]);
    write(&dir, "trim_keep", &wp, &[Plane3::from_origin_normal(o, n)]);
    write(&dir, "trim_off", &wp, &[Plane3::from_origin_normal(o, UnitVec3::try_from_vec(-n.as_vec()).unwrap())]);
    std::fs::write(format!("{dir}/trim_plane.json"),
        format!("{{\"o\":[{},{},{}],\"n\":[{},{},{}]}}", o.x, o.y, o.z, n.as_vec().x, n.as_vec().y, n.as_vec().z)).unwrap();
}
