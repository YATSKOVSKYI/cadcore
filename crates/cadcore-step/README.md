# cadcore-step

**High-Precision, Pure-Rust STEP AP203/AP214 Exporter.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![crates.io](https://img.shields.io/crates/v/cadcore-step)](https://crates.io/crates/cadcore-step)

`cadcore-step` is a native, pure-Rust STEP (ISO 10303-21) writer for the `cadcore` CAD kernel. It translates Boundary Representation (B-Rep) topological shells into standard CAD exchange files.

**Zero FFI. Zero C++ dependencies. Generates exact mathematical shapes (no polygonal approximations).**

---

## Why cadcore-step?

1.  **Exact Analytic Entities:** Exported models map directly to exact mathematical STEP entities (`PLANE`, `CYLINDRICAL_SURFACE`, `TOROIDAL_SURFACE`, etc.). The resulting files contain true analytic curves and surfaces rather than discretized triangular meshes, ensuring infinite zoom resolution and flawless CAD importing.
2.  **Broad CAD Compatibility:** Generates standard AP203 (`CONFIGURATION_CONTROLLED_DESIGN` protocol) or AP214 assembly structures that open out-of-the-box in major CAD engines including **SolidWorks, Autodesk Inventor, Fusion 360, FreeCAD, Rhino, CATIA, and AutoCAD**.
3.  **Sub-Micron Precision:** Coordinates, directions, radii, and transforms are serialized with `:.10` (10 decimal places) precision, preserving watertight manifolds and tiny geometric tolerances down to sub-atomic scales ($10^{-10}$ mm).
4.  **No Temporary Files:** Writes directly to memory strings or output buffers, avoiding intermediate disk reads/writes. Perfect for high-concurrency microservices and WebAssembly runtimes.

---

## Entity Mapping

| cadcore Type | STEP AP203 Entity | Description |
|---|---|---|
| `Plane3` | `PLANE` | Flat boundary plane |
| `CylSurf` | `CYLINDRICAL_SURFACE` | Cylindrical face |
| `SphereSurf` | `SPHERICAL_SURFACE` | Spherical cap |
| `TorusSurf` | `TOROIDAL_SURFACE` | Toroidal blend face |
| `Circle3` | `CIRCLE` | Circular edge boundary |
| `Ellipse3` | `ELLIPSE` | Elliptical edge boundary |
| `Line3` | `LINE` | Linear edge boundary |
| `Point3` | `CARTESIAN_POINT` | Topological vertex point |
| `Vec3` / `UnitVec3` | `DIRECTION` | Normalized direction vector |
| `Frame3` | `AXIS2_PLACEMENT_3D` | Orthonormal coordinate system placement |

---

## Usage

Add `cadcore-step` to your `Cargo.toml`:

```toml
[dependencies]
cadcore-step = "0.1.22"
```

```rust
use cadcore_topo::BRep;
use cadcore_step::brep_to_step;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Initialize and populate BRep (usually via cadcore-ops)
    let mut brep = BRep::new();
    // (populate geometry here...)

    // 2. Serialize directly to STEP string
    let step_text: String = brep_to_step(&brep)?;
    
    // 3. Save to disk
    std::fs::write("model.step", &step_text)?;
    Ok(())
}
```

---

## License & Contact

Licensed under the **MIT License** (see [LICENSE](../../LICENSE)). Free for commercial and non-commercial application.

For questions, support, or custom integrations, please contact Dmytro Yatskovskyi at [dmytroyatskovskyi@outlook.com](mailto:dmytroyatskovskyi@outlook.com).
