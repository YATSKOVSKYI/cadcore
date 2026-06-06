# cadcore

![cadcore banner](https://raw.githubusercontent.com/YATSKOVSKYI/cadcore/master/banner.png)

**The Single-Entry Facade Crate for the Pure-Rust CAD Geometry Kernel.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![crates.io](https://img.shields.io/crates/v/cadcore)](https://crates.io/crates/cadcore)

`cadcore` is a unified facade crate that re-exports all sub-crates of the `cadcore` CAD kernel under a single namespace. 

Instead of adding five different crates to your `Cargo.toml`, you only need to add `cadcore`. It compiles down to a lightweight, thread-safe, and dependency-free library for mathematically exact 3D geometry.

---

## Features

*   **Exact Analytic Representation:** Curves (Line, Circle, Ellipse) and Surfaces (Plane, Cylinder, Sphere, Torus) map directly to STEP entities. Zero polygonal approximations in the core representation.
*   **Safe B-Rep Topology:** Strongly-typed entity IDs using slotmaps. Safe against memory corruptions and logical mismatched ID bugs.
*   **$O(N)$ Filament Sweep:** Build watertight B-Rep tubes along polyline paths in linear time without relying on heavy Boolean union operations.
*   **Direct STEP Export:** Generate high-precision, manifold closed shell STEP AP203/AP214 files directly from code.

---

## Quick Start

Add `cadcore` to your `Cargo.toml`:

```toml
[dependencies]
cadcore = "0.1.22"
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
        Point3::new(15.0, 0.0, 0.0),
        Point3::new(15.0, 5.0, 1.0),
    ];

    let mut brep = BRep::new();
    let opts = SweepOptions::default();
    
    // Sweep a circle profile along the path
    sweep_circle_along_polyline(&mut brep, &waypoints, 0.3, &opts)?;

    // Generate STEP string content
    let step_data = brep_to_step(&brep)?;
    std::fs::write("scaffold.step", &step_data)?;

    Ok(())
}
```

---

## Sub-crate Ecosystem

This facade re-exports the following modular building blocks:

| Module | Crate Name | Purpose |
|---|---|---|
| `cadcore::math` | [`cadcore-math`](https://crates.io/crates/cadcore-math) | Math primitives (points, vectors, transforms). Zero dependencies. |
| `cadcore::geom` | [`cadcore-geom`](https://crates.io/crates/cadcore-geom) | Curves, surfaces, and analytic geometry calculations. |
| `cadcore::topo` | [`cadcore-topo`](https://crates.io/crates/cadcore-topo) | Arena-backed Boundary Representation (B-Rep) topological database. |
| `cadcore::ops` | [`crates/cadcore-ops`](https://crates.io/crates/cadcore-ops) | High-level sweep construction and half-space trimming solvers. |
| `cadcore::step` | [`crates/cadcore-step`](https://crates.io/crates/cadcore-step) | Native STEP AP203 writer. |

---

## Why cadcore is the Best Choice

1.  **Pure Rust:** Say goodbye to complex C++ link errors, dll dependencies, and compilation issues. `cargo add` and build instantly.
2.  **WebAssembly Native:** Runs in browser client applications for real-time CAD model generation and downloads.
3.  **Engineered for Speed:** Building swept paths (tubes, rods, filaments) takes $O(N)$ linear time instead of $O(N^2)$ in traditional solid modelers.
4.  **No Global State:** Safe to run in parallel in multi-threaded environments.

---

## License & Contact

Licensed under the **MIT License** (see [LICENSE](../../LICENSE)). Free for commercial and non-commercial application.

For questions, support, or custom integrations, please contact Dmytro Yatskovskyi at [dmytroyatskovskyi@outlook.com](mailto:dmytroyatskovskyi@outlook.com).
