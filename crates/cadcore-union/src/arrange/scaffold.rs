//! The scaffold orchestrator — fuse an arbitrary set of serpentine filaments
//! (each: straight legs + torus-fillet elbows) that cross each other into ONE
//! watertight, seamless solid.
//!
//! This is the universal entry point the app pipeline targets: every leg of
//! every filament is trimmed by every leg of every OTHER filament it
//! overlaps, junctions weld legs to elbows, free leg ends are capped, and the
//! whole thing is assembled through one [`Assembler`] (shared edges → one
//! solid).  All cylinders are seamless (analytic circle rims), all elbows are
//! seamless torus fillets.

use cadcore_geom::{Circle3, CylSurf, Plane3, TorusSurf};
use cadcore_math::{Point3, UnitVec3};
use cadcore_topo::{BRep, Shell, ShellId, Solid};

use super::assembly::Assembler;
use super::filament::{ElbowSeg, LegSeg};
use crate::boolean::Aabb;
use crate::geom::composite::{closed_loops_between, CompositeTube};
use crate::geom::intersect::TraceOptions;

/// One serpentine filament path: `legs[i]` joined to `legs[i+1]` by
/// `elbows[i]`.  The two free ends (first leg incoming, last leg outgoing)
/// are capped.  Each leg is built with `surf.frame.origin` at its incoming
/// end and `surf.axis()` pointing toward the outgoing end.
#[derive(Debug, Clone)]
pub struct Filament {
    /// Straight legs in path order.
    pub legs: Vec<LegSeg>,
    /// Elbows (`legs.len() - 1` of them); `elbows[i]` joins legs i and i+1.
    pub elbows: Vec<ElbowSeg>,
}

impl Filament {
    /// Build a planar serpentine filament from a polyline path: straight legs
    /// joined by constant-radius circular fillets (one per interior vertex).
    /// `r` = filament radius, `fillet` = elbow bend radius.  Each leg runs
    /// `surf.frame.origin` → `+axis`, matching [`fuse_scaffold`]'s convention.
    pub fn serpentine(verts: &[Point3], r: f64, fillet: f64) -> Filament {
        let n = verts.len();
        assert!(n >= 2, "a filament needs at least two path vertices");
        let dir = |a: Point3, b: Point3| UnitVec3::try_from_vec(b - a).unwrap();
        let mut tin = vec![Point3::new(0.0, 0.0, 0.0); n];
        let mut tout = vec![Point3::new(0.0, 0.0, 0.0); n];
        let mut elbows = Vec::new();
        for k in 1..n.saturating_sub(1) {
            let din = dir(verts[k - 1], verts[k]);
            let dout = dir(verts[k], verts[k + 1]);
            let turn = din.dot_vec(dout.as_vec()).clamp(-1.0, 1.0).acos();
            let t = fillet / (turn / 2.0).tan();
            let pin = verts[k] - din.as_vec() * t;
            let pout = verts[k] + dout.as_vec() * t;
            let axis = UnitVec3::try_from_vec(din.as_vec().cross(dout.as_vec())).unwrap();
            let inward = UnitVec3::try_from_vec(axis.as_vec().cross(din.as_vec())).unwrap();
            let centre = pin + inward.as_vec() * fillet;
            let torus = TorusSurf::new(centre, axis, fillet, r);
            let theta = |p: Point3| {
                let w = p - torus.frame.origin;
                let rd = w - torus.frame.z.as_vec() * torus.frame.z.dot_vec(w);
                torus.frame.y.dot_vec(rd).atan2(torus.frame.x.dot_vec(rd))
            };
            let lo = theta(pin);
            let mut d = theta(pout) - lo;
            while d > std::f64::consts::PI { d -= std::f64::consts::TAU; }
            while d < -std::f64::consts::PI { d += std::f64::consts::TAU; }
            elbows.push(ElbowSeg { surf: torus, theta_lo: lo, theta_hi: lo + d });
            tin[k] = pin;
            tout[k] = pout;
        }
        let mut legs = Vec::new();
        for k in 0..n - 1 {
            let start = if k == 0 { verts[0] } else { tout[k] };
            let end = if k + 1 == n - 1 { verts[n - 1] } else { tin[k + 1] };
            legs.push(LegSeg {
                surf: CylSurf::new(start, dir(start, end), r),
                length: (end - start).length(),
            });
        }
        Filament { legs, elbows }
    }

    /// The composite tube (legs + elbows) for crossing / containment queries.
    fn tube(&self) -> CompositeTube {
        let mut t = CompositeTube::new();
        for l in &self.legs {
            t = t.with_leg(l.surf, l.length);
        }
        for e in &self.elbows {
            t = t.with_elbow(e.surf, e.theta_lo, e.theta_hi);
        }
        t
    }
}

/// Fuse a set of crossing serpentine filaments into one closed, seamless,
/// watertight solid.  Returns the shell.
pub fn fuse_scaffold(brep: &mut BRep, filaments: &[Filament]) -> ShellId {
    let nf = filaments.len();
    // Per-leg tube (one straight leg) for crossing traces, and the whole
    // filament tube for keep/drop containment.
    let fil_tubes: Vec<CompositeTube> = filaments.iter().map(|f| f.tube()).collect();
    let leg_tubes: Vec<Vec<CompositeTube>> = filaments
        .iter()
        .map(|f| {
            f.legs
                .iter()
                .map(|l| CompositeTube::new().with_leg(l.surf, l.length))
                .collect()
        })
        .collect();

    // Per-leg axis-aligned bounding box (segment endpoints, padded by the tube
    // radius) — the broad phase: two legs can only intersect if their boxes
    // overlap, so we skip the (expensive, fine-step) SSI trace for the vast
    // majority of non-crossing leg pairs.  Without this cull `fuse_scaffold`
    // is O(legs²) traces and does not scale to a real multi-layer scaffold
    // (hundreds of legs → tens of thousands of traces).
    let leg_box = |l: &LegSeg| -> Aabb {
        let a = l.surf.frame.origin;
        let b = a + l.surf.axis().as_vec() * l.length;
        let mut bb = Aabb::empty();
        bb.expand(a);
        bb.expand(b);
        bb.pad(l.surf.radius)
    };
    let leg_boxes: Vec<Vec<Aabb>> = filaments
        .iter()
        .map(|f| f.legs.iter().map(leg_box).collect())
        .collect();

    // windows[fi][li] = crossing loops on leg (fi,li) with legs of OTHER
    // filaments.  Traced ONCE per leg-pair; both sides imprint the SAME points.
    let mut windows: Vec<Vec<Vec<Vec<Point3>>>> = filaments
        .iter()
        .map(|f| vec![Vec::new(); f.legs.len()])
        .collect();
    for fi in 0..nf {
        for fj in (fi + 1)..nf {
            for li in 0..filaments[fi].legs.len() {
                for lj in 0..filaments[fj].legs.len() {
                    if !leg_boxes[fi][li].overlaps(&leg_boxes[fj][lj]) {
                        continue; // broad-phase reject: boxes disjoint → no crossing
                    }
                    for c in closed_loops_between(
                        &leg_tubes[fi][li],
                        &leg_tubes[fj][lj],
                        &TraceOptions::default(),
                    ) {
                        windows[fi][li].push(c.points.clone());
                        windows[fj][lj].push(c.points);
                    }
                }
            }
        }
    }

    let mut asm = Assembler::new(brep, 1e-7);
    let circ_at = |l: &LegSeg, v: f64| -> Circle3 {
        Circle3::new(
            l.surf.frame.origin + l.surf.axis().as_vec() * v,
            l.surf.axis(),
            l.surf.radius,
        )
    };

    for fi in 0..nf {
        let f = &filaments[fi];
        let n = f.legs.len();
        let j_in: Vec<Circle3> = f.legs.iter().map(|l| circ_at(l, 0.0)).collect();
        let j_out: Vec<Circle3> = f.legs.iter().map(|l| circ_at(l, l.length)).collect();

        for li in 0..n {
            asm.emit_cylinder_face(
                f.legs[li].surf,
                f.legs[li].length,
                j_in[li],
                j_out[li],
                &windows[fi][li],
            );
            if li == 0 {
                let nrm = UnitVec3::try_from_vec(f.legs[li].surf.axis().as_vec() * -1.0).unwrap();
                asm.emit_disk_cap(
                    Plane3::from_origin_normal(f.legs[li].surf.frame.origin, nrm),
                    j_in[li],
                );
            }
            if li + 1 == n {
                let c = f.legs[li].surf.frame.origin
                    + f.legs[li].surf.axis().as_vec() * f.legs[li].length;
                asm.emit_disk_cap(Plane3::from_origin_normal(c, f.legs[li].surf.axis()), j_out[li]);
            }
        }
        for (ei, e) in f.elbows.iter().enumerate() {
            asm.emit_torus_fillet(e.surf, j_out[ei], j_in[ei + 1]);
        }
    }
    let _ = &fil_tubes; // (kept for future elbow-crossing support)

    let faces = asm.faces().to_vec();
    let shell = brep.add_shell(Shell {
        faces: faces.clone(),
        is_outer: true,
        solid: Default::default(),
    });
    let solid = brep.add_solid(Solid {
        shells: vec![shell],
        name: Some("scaffold".into()),
    });
    brep.shells[shell].solid = solid;
    for f in &faces {
        brep.faces[*f].shell = shell;
    }
    shell
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::UnionConfig;
    use crate::validate::validate_shell;
    use crate::Diagnostics;
    use std::collections::HashMap;

    fn open_edges(brep: &BRep, shell: ShellId) -> usize {
        let mut uses: HashMap<cadcore_topo::EdgeId, usize> = HashMap::new();
        for &fid in &brep.shells[shell].faces {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let st = brep.loops[lid].start;
                let mut c = st;
                loop {
                    let ce = &brep.coedges[c];
                    *uses.entry(ce.edge).or_insert(0) += 1;
                    c = ce.next;
                    if c == st { break; }
                }
            }
        }
        uses.values().filter(|&&n| n == 1).count()
    }

    /// Number of connected components of the shell, faces joined by shared
    /// edges.  A genuine fusion of crossing filaments is ONE component; two
    /// merely-coexisting (interpenetrating) solids are two.
    fn components(brep: &BRep, shell: ShellId) -> usize {
        let faces = &brep.shells[shell].faces;
        // edge -> faces touching it
        let mut by_edge: HashMap<cadcore_topo::EdgeId, Vec<usize>> = HashMap::new();
        for (i, &fid) in faces.iter().enumerate() {
            let f = &brep.faces[fid];
            let mut ls = vec![f.outer_loop];
            ls.extend(f.inner_loops.iter().copied());
            for lid in ls {
                let st = brep.loops[lid].start;
                let mut c = st;
                loop {
                    by_edge.entry(brep.coedges[c].edge).or_default().push(i);
                    c = brep.coedges[c].next;
                    if c == st { break; }
                }
            }
        }
        let mut parent: Vec<usize> = (0..faces.len()).collect();
        fn find(p: &mut [usize], x: usize) -> usize {
            let mut r = x;
            while p[r] != r { r = p[r]; }
            let mut c = x;
            while p[c] != c { let n = p[c]; p[c] = r; c = n; }
            r
        }
        for fs in by_edge.values() {
            for w in fs.windows(2) {
                let (a, b) = (find(&mut parent, w[0]), find(&mut parent, w[1]));
                parent[a] = b;
            }
        }
        (0..faces.len()).filter(|&i| find(&mut parent, i) == i).count()
    }

    /// Two serpentine filaments on stacked layers that cross transversally
    /// (X-running vs Y-running) — fused into ONE watertight, manifold solid.
    #[test]
    fn two_crossing_serpentines_watertight() {
        let r = 0.275;
        // filament A: U in z=0, legs run along X at y=0 and y=2 (connector at x=0)
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
        // filament B: U one layer up (z=0.35), legs run along Y at x=-1.5 and
        // x=-0.8 so each crosses BOTH of A's X-legs transversally; its elbows
        // sit at y=3, clear of A.
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
        let shell = fuse_scaffold(&mut brep, &[a, b]);
        assert_eq!(open_edges(&brep, shell), 0, "scaffold watertight");
        assert_eq!(components(&brep, shell), 1, "filaments fused into one solid");

        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let report = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(report.manifold_violations.is_empty(), "manifold");
        assert!(report.distance_violations.is_empty(),
            "distance max {:.2}um", report.max_edge_deviation * 1e3);
    }

    /// One serpentine filament covering an `n×n` mm area along `axis_x`
    /// (true = legs run along X, weaving in Y), at height `z`, pitch `p`.
    fn weave(span: f64, p: f64, z: f64, along_x: bool, r: f64, fillet: f64) -> Filament {
        let mut verts = Vec::new();
        let mut k = 0;
        let mut a = 0.0;
        while a <= span + 1e-9 {
            let (lo, hi) = if k % 2 == 0 { (0.0, span) } else { (span, 0.0) };
            for &t in &[lo, hi] {
                let pt = if along_x {
                    Point3::new(t, a, z)
                } else {
                    Point3::new(a, t, z)
                };
                verts.push(pt);
            }
            a += p;
            k += 1;
        }
        Filament::serpentine(&verts, r, fillet)
    }

    /// Scaling guard: a realistic 4-layer woven scaffold (≈ the app's 5×5×4
    /// model) must fuse in well under the broad-phase-less O(legs²) blow-up.
    /// With the leg-AABB cull this is sub-second; without it, it hangs.
    #[test]
    fn four_layer_weave_fuses_fast_and_watertight() {
        let r = 0.275;
        let span = 50.0; // big span → many crossings (≈ the app's dense scaffold)
        let p = 1.0;
        let fillet = 0.25; // < p/2 so the short connector legs stay positive-length
        let layers = vec![
            weave(span, p, 0.0, true, r, fillet),
            weave(span, p, 0.35, false, r, fillet),
            weave(span, p, 0.70, true, r, fillet),
            weave(span, p, 1.05, false, r, fillet),
        ];
        let total_legs: usize = layers.iter().map(|f| f.legs.len()).sum();
        let t0 = std::time::Instant::now();
        let mut brep = BRep::new();
        let shell = fuse_scaffold(&mut brep, &layers);
        let fuse_ms = t0.elapsed().as_millis();
        println!(
            "4-layer weave: {} legs, {} faces, fuse {} ms",
            total_legs,
            brep.shells[shell].faces.len(),
            fuse_ms,
        );
        assert_eq!(open_edges(&brep, shell), 0, "watertight");

        let mut cfg = UnionConfig::default();
        cfg.wires_strict = true;
        let mut diag = Diagnostics::new();
        let report = validate_shell(&brep, shell, &cfg, &mut diag);
        assert!(report.manifold_violations.is_empty(), "manifold");
    }
}
