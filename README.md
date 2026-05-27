# cadcore

**A pure-Rust CAD geometry kernel.**

Zero C++ dependencies. No OCCT. No global state. Fully parallelisable.

---

## What it does

* **B-Rep topology** — arena-based, typed IDs (`VertexId`, `EdgeId`, `FaceId`, …)
* **Analytic surfaces** — `Plane3`, `CylSurf`, `SphereSurf`, `TorusSurf` — exact, no tessellation
* **O(N) filament sweep** — build a solid B-Rep from a polyline path in linear time (no Boolean union)
* **STEP AP203 export** — pure Rust, direct analytic entity mapping, no temp files

---

## Quick start

```toml
[dependencies]
cadcore = { path = "../cadcore/crates/cadcore" }
```

```rust
use cadcore::{
    math::Point3,
    topo::BRep,
    ops::{sweep_circle_along_polyline, SweepOptions},
    step::brep_to_step,
};

let waypoints = vec![
    Point3::new( 0.0, 0.0, 0.0),
    Point3::new(10.0, 0.0, 0.0),
    Point3::new(10.0, 8.0, 0.0),
];

let mut brep = BRep::new();
sweep_circle_along_polyline(&mut brep, &waypoints, 0.2, &SweepOptions::default())?;
let step_text = brep_to_step(&brep)?;
std::fs::write("output.step", &step_text)?;
```

Run the built-in example:

```
cargo run -p cadcore --example scaffold_rod
```

---

## Crate layout

| Crate            | Contents                                                   |
|------------------|------------------------------------------------------------|
| `cadcore-math`   | `Point3`, `Vec3`, `UnitVec3`, `Mat3`, `Frame3`, `Interval` |
| `cadcore-geom`   | `Line3`, `Circle3`, `Ellipse3`, `Plane3`, `CylSurf`, `TorusSurf`, `SphereSurf` |
| `cadcore-topo`   | Arena B-Rep: `BRep`, `Solid`, `Shell`, `Face`, `Edge`, `Vertex` |
| `cadcore-ops`    | `sweep_circle_along_polyline` — O(N) analytic sweep        |
| `cadcore-step`   | `brep_to_step`, `StepWriter` — STEP AP203 exporter         |
| `cadcore`        | Facade re-exporting all of the above                       |

---

## Design decisions

### Why O(N) sweep instead of Boolean union?

For a filament path of N segments, traditional OCCT does an O(N²) fuse to
join N cylinder solids. `cadcore` instead builds the topology analytically:
each segment is a `CylSurf`, each corner is a `TorusSurf`, end caps are
`Plane3`. The result is the same solid B-Rep at O(N) cost.

### Why exact analytic surfaces?

Every surface type maps directly to a STEP AP203 entity:

| cadcore type  | STEP entity              |
|---------------|--------------------------|
| `Plane3`      | `PLANE`                  |
| `CylSurf`     | `CYLINDRICAL_SURFACE`    |
| `SphereSurf`  | `SPHERICAL_SURFACE`      |
| `TorusSurf`   | `TOROIDAL_SURFACE`       |
| `Circle3`     | `CIRCLE`                 |
| `Ellipse3`    | `ELLIPSE`                |

No mesh approximation, no tessellation error, exact geometry in the exported file.

### Why typed arena IDs?

`slotmap::new_key_type!` gives us distinct Rust types for each entity kind.
You can't accidentally pass a `FaceId` where a `VertexId` is expected.
All IDs are stable across insertions (no pointer invalidation).

---

## Status

`v0.1.0` — foundation complete, STEP output works for filament sweeps.

Roadmap:
- [ ] Full edge-loop stitching (boundary co-edges between adjacent faces)
- [ ] Miter-plane joints (sharp-corner alternative to torus fillets)
- [ ] 3MF export
- [ ] IGES export
- [ ] Parallel sweep for multi-path scaffold generation
