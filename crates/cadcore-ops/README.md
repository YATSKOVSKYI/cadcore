# cadcore-ops

**High-Level Analytic B-Rep Operations and Sweep Solvers.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![crates.io](https://img.shields.io/crates/v/cadcore-ops)](https://crates.io/crates/cadcore-ops)

`cadcore-ops` contains the core geometric algorithms of the `cadcore` kernel. It features a custom linear-time sweep solver that builds mathematically exact B-Rep solids directly from polyline paths without relying on expensive, failure-prone Boolean union operations.

---

## The O(N) Sweep Advantage

When modeling long tubes, pipes, or 3D-printing filaments, legacy C++ CAD kernels (like OpenCASCADE) must build individual cylinder segments and fuse them together using Boolean operations. This results in $O(N^2)$ complexity, high memory overhead, and frequent failures at self-intersecting sharp bends.

`cadcore-ops` sweeps a circular profile along a path by analytically constructing:
*   **2 planar end caps**
*   **N-1 cylinder faces**
*   **N-2 connector faces** (either G1-smooth torus fillets or exact elliptical miter joins)

This algorithm executes in **linear $O(N)$ time**, uses minimal memory, and guarantees a watertight manifold solid.

---

## Key Operations

### 🔄 `sweep_circle_along_polyline`
Builds an exact analytic B-Rep solid by sweeping a circle of a given radius along a polyline. Ideal for structural framing, pipes, and 3D-printing toolpaths.

### ✂️ `half_space_cut_brep`
Applies sequential flat/oblique cutting planes to B-Rep solids, trimming cylinder segments laterally or axially and capping them with flat disks or partial disks. Excellent for modeling sub-spindle cuts, electrodes, and trimmed structural ends.

---

## Usage

Add `cadcore-ops` to your `Cargo.toml`:

```toml
[dependencies]
cadcore-ops = "0.1.22"
```

```rust
use cadcore_topo::BRep;
use cadcore_math::Point3;
use cadcore_ops::{sweep_circle_along_polyline, SweepOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut brep = BRep::new();
    let waypoints = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(10.0, 10.0, 0.0),
    ];

    let opts = SweepOptions {
        fillet_corners: true,            // G1-smooth torus corner fillets
        corner_fillet_radius: 0.0,       // Defaults to profile radius
        name: Some("filament".to_string()),
    };

    // 3 waypoints -> 2 cylinders + 1 torus fillet + 2 end caps = 5 faces
    let solid_id = sweep_circle_along_polyline(&mut brep, &waypoints, 0.2, &opts)?;
    
    let stats = brep.stats();
    assert_eq!(stats.faces, 5);
    Ok(())
}
```

---

## License

Licensed under the **MIT License** (see [LICENSE](../../LICENSE)). Free for commercial and non-commercial application.
