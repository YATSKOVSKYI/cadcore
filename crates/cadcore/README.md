# cadcore

[![crates.io](https://img.shields.io/crates/v/cadcore)](https://crates.io/crates/cadcore)
[![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm%20Noncommercial%201.0.0-blue)](../../LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.78%2B-orange)](https://www.rust-lang.org/)

**A pure-Rust CAD geometry kernel.**

Zero C++ dependencies. No OpenCASCADE. No global state. Fully parallelisable.

This crate is the **facade** — it re-exports everything from all cadcore sub-crates under a single dependency.

---

## Quick start

```toml
[dependencies]
cadcore = "0.1"
```

```rust
use cadcore::{
    math::Point3,
    topo::BRep,
    ops::{sweep_circle_along_polyline, SweepOptions},
    step::brep_to_step,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let waypoints = vec![
        Point3::new( 0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(10.0, 8.0, 0.0),
        Point3::new( 0.0, 8.0, 0.0),
    ];

    let mut brep = BRep::new();
    sweep_circle_along_polyline(&mut brep, &waypoints, 0.2, &SweepOptions::default())?;

    let step = brep_to_step(&brep)?;
    std::fs::write("scaffold_leg.step", &step)?;
    Ok(())
}
```

---

## What cadcore provides

| Feature | Detail |
|---|---|
| **Exact analytic geometry** | Planes, cylinders, spheres, toruses — no mesh approximation |
| **Arena-based B-Rep topology** | `Solid → Shell → Face → Loop → Edge → Vertex` with typed stable IDs |
| **O(N) filament sweep** | Build a solid from a polyline path without Boolean union |
| **Pure-Rust STEP AP203 export** | Direct analytic STEP entities — opens in FreeCAD, CATIA, SolidWorks, Rhino |

---

## Architecture

```
cadcore/
├── crates/
│   ├── cadcore-math    — Point3, Vec3, UnitVec3, Mat3, Frame3, Transform3, Interval
│   ├── cadcore-geom    — Line3, Circle3, Ellipse3 / Plane3, CylSurf, TorusSurf, SphereSurf
│   ├── cadcore-topo    — BRep arena: Solid, Shell, Face, Loop, CoEdge, Edge, Vertex
│   ├── cadcore-ops     — sweep_circle_along_polyline (O(N) analytic construction)
│   ├── cadcore-step    — STEP AP203 writer (pure Rust)
│   └── cadcore         — this facade crate
```

### Dependency graph

```
cadcore-math  (no deps)
     ↑
cadcore-geom
     ↑
cadcore-topo
     ↑
cadcore-ops ──────────────────┐
     ↑                        ↓
cadcore-step          cadcore (facade)
```

---

## Crates

| Crate | crates.io | Description |
|---|---|---|
| [`cadcore-math`](https://crates.io/crates/cadcore-math) | [![](https://img.shields.io/crates/v/cadcore-math)](https://crates.io/crates/cadcore-math) | Math primitives (zero deps) |
| [`cadcore-geom`](https://crates.io/crates/cadcore-geom) | [![](https://img.shields.io/crates/v/cadcore-geom)](https://crates.io/crates/cadcore-geom) | Analytic curves and surfaces |
| [`cadcore-topo`](https://crates.io/crates/cadcore-topo) | [![](https://img.shields.io/crates/v/cadcore-topo)](https://crates.io/crates/cadcore-topo) | B-Rep arena topology |
| [`cadcore-ops`](https://crates.io/crates/cadcore-ops) | [![](https://img.shields.io/crates/v/cadcore-ops)](https://crates.io/crates/cadcore-ops) | Geometric operations |
| [`cadcore-step`](https://crates.io/crates/cadcore-step) | [![](https://img.shields.io/crates/v/cadcore-step)](https://crates.io/crates/cadcore-step) | STEP AP203 writer |

---

## Why not OCCT?

| | OCCT | cadcore |
|---|---|---|
| Language | C++ (FFI required) | Pure Rust |
| Boolean union for N paths | O(N²) fuse | O(N) analytic construction |
| STEP export | Temp file + parse back | Direct string serialisation |
| Binary size | ~200 MB prebuilt | < 1 MB |
| Windows install | Complex setup | `cargo add cadcore` |

---

## Roadmap

- [ ] Full edge-loop stitching (co-edge topology between adjacent faces)
- [ ] Miter-plane joints (sharp corners)
- [ ] STEP import (AP203 / AP214 reader)
- [ ] 3MF export
- [ ] Offset surface operations (wall thickness)
- [ ] Parallel multi-path sweep
- [ ] WASM target (browser-side CAD)

---

## License

[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/)

Free for research, education, personal projects, and non-commercial use.
For commercial licensing, contact **dmytroyatskovskiy@gmail.com**.
