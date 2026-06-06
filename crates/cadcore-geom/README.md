# cadcore-geom

**Exact Analytic curves and Surfaces for CAD Kernels.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![crates.io](https://img.shields.io/crates/v/cadcore-geom)](https://crates.io/crates/cadcore-geom)

`cadcore-geom` contains the analytic curves and surfaces representing 3D geometry in the `cadcore` kernel. Every primitive in this crate represents exact mathematical surfaces rather than discretized mesh models.

**Zero Heap Allocations. Zero FFI dependencies. All types are `Copy` and `f64`-based.**

---

## 📈 Curves

| Curve Type | STEP Entity | Description |
|---|---|---|
| `Line3` | `LINE` | Infinite directed 3D line |
| `Circle3` | `CIRCLE` | Planar circle defined by center, normal, and radius |
| `Ellipse3` | `ELLIPSE` | Planar ellipse defined by semi-major/minor axes |
| `BezierCubic` | — | Cubic Bézier curve (useful for path approximations) |

## 🔲 Surfaces

| Surface Type | STEP Entity | Description |
|---|---|---|
| `Plane3` | `PLANE` | Infinite flat plane |
| `CylSurf` | `CYLINDRICAL_SURFACE` | Right circular cylinder with defined axis and radius |
| `SphereSurf` | `SPHERICAL_SURFACE` | Full 3D sphere with defined center and radius |
| `TorusSurf` | `TOROIDAL_SURFACE` | Ring torus (used for G1-smooth corner fillets) |
| `ConeSurf` | `CONICAL_SURFACE` | Right circular cone |

---

## Why cadcore-geom?

1.  **Exact Mathematics:** No chordal deviation or tessellation errors. Surfaces and curves are represented analytically using their true geometric formulas.
2.  **1-to-1 STEP Mapping:** Every type corresponds directly to an ISO 10303-21 (STEP) entity, ensuring flawless translations during exports.
3.  **WASM Friendly:** Designed with zero allocations, making it extremely fast when compiled to WebAssembly for browser applications.

---

## Usage

Add `cadcore-geom` to your `Cargo.toml`:

```toml
[dependencies]
cadcore-geom = "0.1.22"
```

```rust
use cadcore_geom::{CylSurf, Plane3};
use cadcore_math::{Point3, UnitVec3, Frame3};

fn main() {
    // 1. Create a cylindrical surface (axis along Z-axis, radius 0.5 mm)
    let cyl = CylSurf::new(Point3::ORIGIN, UnitVec3::Z, 0.5);
    
    // Evaluate coordinate points on the surface (theta = 0, z = 10 mm)
    let pt = cyl.point_at(0.0, 10.0);
    assert!((pt.x - 0.5).abs() < 1e-10);

    // 2. Project points onto a Plane
    let plane = Plane3::new(Frame3::IDENTITY);
    let projected = plane.project(Point3::new(1.0, 2.0, 3.0));
    assert_eq!(projected, Point3::new(1.0, 2.0, 0.0));
}
```

---

## License & Contact

Licensed under the **MIT License** (see [LICENSE](../../LICENSE)). Free for commercial and non-commercial application.

For questions, support, or custom integrations, please contact Dmytro Yatskovskyi at [dmytroyatskovskyi@outlook.com](mailto:dmytroyatskovskyi@outlook.com).
