# cadcore

**A pure-Rust CAD geometry kernel.**

[![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm%20Noncommercial%201.0.0-blue)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.78%2B-orange)](https://www.rust-lang.org/)
[![crates.io](https://img.shields.io/crates/v/cadcore)](https://crates.io/crates/cadcore)

Zero C++ dependencies. No OpenCASCADE. No global state. Fully parallelisable.

---

## Overview

`cadcore` is a CAD geometry kernel written entirely in Rust. It provides:

- **Exact analytic geometry** — planes, cylinders, spheres, toruses — no mesh approximation
- **Arena-based B-Rep topology** — `Solid → Shell → Face → Loop → Edge → Vertex` with typed stable IDs
- **O(N) filament sweep** — build a solid from a polyline path without Boolean union operations
- **Pure-Rust STEP AP203 export** — direct analytic surface entities, opens in FreeCAD, CATIA, SolidWorks, Rhino

Designed for computational manufacturing: bioprinting scaffolds, filament-path CAD, lattice structures.

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
    // A simple U-shaped filament path
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

Run the built-in example:

```bash
cargo run -p cadcore --example scaffold_rod
```

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
│   └── cadcore         — facade re-exporting everything
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

## Crate reference

### `cadcore-math`

Zero external dependencies. All geometry starts here.

| Type | Description |
|---|---|
| `Point3` | Location in 3-D space (mm) |
| `Vec3` | Free vector (direction + magnitude) |
| `UnitVec3` | Unit-length vector, enforced at construction |
| `Mat3` | 3×3 column-major matrix (rotation, linear maps) |
| `Frame3` | Right-handed orthonormal frame (origin + 3 axes) |
| `Transform3` | Rigid-body transform: rotation + translation |
| `Interval` | Closed real interval `[lo, hi]` |

```rust
use cadcore::math::{Point3, Vec3, UnitVec3, Frame3};

let origin = Point3::new(1.0, 2.0, 3.0);
let dir    = UnitVec3::try_from_vec(Vec3::new(1.0, 0.0, 0.0)).unwrap();
let frame  = Frame3::from_origin_z(origin, dir);
```

### `cadcore-geom`

Analytic curves and surfaces. All types are `Copy` and parameterised by `f64`.

**Curves**

| Type | STEP entity | Description |
|---|---|---|
| `Line3` | `LINE` | Infinite directed line |
| `Circle3` | `CIRCLE` | Planar circle |
| `Ellipse3` | `ELLIPSE` | Planar ellipse (e.g. miter cross-section) |
| `BezierCubic` | — | Cubic Bézier (approximation use-cases) |

**Surfaces**

| Type | STEP entity | Description |
|---|---|---|
| `Plane3` | `PLANE` | Infinite plane |
| `CylSurf` | `CYLINDRICAL_SURFACE` | Right circular cylinder |
| `SphereSurf` | `SPHERICAL_SURFACE` | Full sphere |
| `TorusSurf` | `TOROIDAL_SURFACE` | Ring torus (filament corner fillet) |
| `ConeSurf` | `CONICAL_SURFACE` | Right circular cone |

```rust
use cadcore::geom::{CylSurf, TorusSurf};
use cadcore::math::{Point3, UnitVec3};

let cyl = CylSurf::new(Point3::ORIGIN, UnitVec3::Z, 0.2);
let p   = cyl.point_at(0.0, 5.0);  // theta=0°, z=5mm
```

### `cadcore-topo`

B-Rep topology stored in typed arenas (`slotmap`). All entity IDs are distinct Rust types — you cannot accidentally pass a `FaceId` where a `VertexId` is expected.

```
Solid
 └─ Shell  (outer + optional inner voids)
     └─ Face  (bounded surface region)
         └─ Loop  (outer boundary + optional holes)
             └─ CoEdge  (directed edge use)
                 └─ Edge  (curve segment + two vertices)
                     └─ Vertex  (3-D point)
```

```rust
use cadcore::topo::BRep;

let mut brep = BRep::new();
// populate via cadcore-ops or manually
let stats = brep.stats();
println!("{} faces, {} solids", stats.faces, stats.solids);
```

### `cadcore-ops`

High-level operations that build B-Rep topology analytically.

#### `sweep_circle_along_polyline`

```rust
pub fn sweep_circle_along_polyline(
    brep:      &mut BRep,
    waypoints: &[Point3],
    radius:    f64,
    opts:      &SweepOptions,
) -> Result<SolidId, SweepError>
```

Builds an exact analytic B-Rep solid by sweeping a circle of given `radius` along a polyline:

- **N−1 cylinder faces** (one per segment)
- **N−2 toroidal fillets** at each interior bend (G1-smooth corners)
- **2 planar end caps**

Total cost: **O(N)**. No Boolean union. No mesh.

```rust
use cadcore::ops::{sweep_circle_along_polyline, SweepOptions};

let opts = SweepOptions {
    fillet_corners: true,               // torus fillets at bends
    name: Some("my_rod".to_string()),
};
let id = sweep_circle_along_polyline(&mut brep, &points, 0.2, &opts)?;
```

### `cadcore-step`

Pure-Rust STEP AP203 writer. Converts a `BRep` to an ISO 10303-21 exchange file.

Every surface maps to a native STEP entity — **no tessellation, no approximation**:

| cadcore | STEP AP203 |
|---|---|
| `Plane3` | `PLANE` |
| `CylSurf` | `CYLINDRICAL_SURFACE` |
| `SphereSurf` | `SPHERICAL_SURFACE` |
| `TorusSurf` | `TOROIDAL_SURFACE` |
| `Circle3` | `CIRCLE` |
| `Ellipse3` | `ELLIPSE` |
| `Point3` | `CARTESIAN_POINT` |
| `Frame3` | `AXIS2_PLACEMENT_3D` |

```rust
use cadcore::step::brep_to_step;

let step_text = brep_to_step(&brep)?;
std::fs::write("output.step", &step_text)?;
```

---

## Why not OCCT?

OpenCASCADE (OCCT) is the de-facto C++ CAD kernel. `cadcore` was built to replace it for filament-path geometry:

| | OCCT | cadcore |
|---|---|---|
| Language | C++ (FFI required) | Pure Rust |
| Thread safety | No global state issues in Rust | ✅ |
| Boolean union for N paths | O(N²) fuse | O(N) analytic construction |
| STEP export | Temp file → parse back | Direct string serialisation |
| Binary size | ~200 MB prebuilt | < 1 MB |
| Windows prebuilt | Complex setup | `cargo add cadcore` |

---

## Roadmap

- [ ] Full edge-loop stitching (co-edge topology between adjacent faces)
- [ ] Miter-plane joints (sharp-corner alternative to torus fillets)
- [ ] STEP import (AP203 / AP214 reader)
- [ ] 3MF export
- [ ] Offset surface operations (wall thickness)
- [ ] Parallel multi-path sweep (scaffold generation)
- [ ] WASM target (browser-side CAD)

---

## License

[PolyForm Noncommercial License 1.0.0](LICENSE)

Free for research, education, personal projects, and non-commercial use.  
For commercial licensing, contact **dmytroyatskovskiy@gmail.com**.

---

## Contributing

Issues and PRs are welcome for bug fixes, documentation, and non-commercial enhancements.  
See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
