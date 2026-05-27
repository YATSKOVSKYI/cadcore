# cadcore-step

[![crates.io](https://img.shields.io/crates/v/cadcore-step)](https://crates.io/crates/cadcore-step)
[![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm%20Noncommercial%201.0.0-blue)](../../LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.78%2B-orange)](https://www.rust-lang.org/)

Pure-Rust STEP AP203 writer for the [cadcore](https://crates.io/crates/cadcore) CAD kernel.

Converts a `BRep` to an ISO 10303-21 exchange file. Every cadcore surface maps directly to a native STEP entity — **no tessellation, no mesh approximation**. The output opens in FreeCAD, SolidWorks, CATIA, Rhino, and any other AP203-compliant CAD tool.

---

## Entity mapping

| cadcore type | STEP AP203 entity |
|---|---|
| `Plane3` | `PLANE` |
| `CylSurf` | `CYLINDRICAL_SURFACE` |
| `SphereSurf` | `SPHERICAL_SURFACE` |
| `TorusSurf` | `TOROIDAL_SURFACE` |
| `ConeSurf` | `CONICAL_SURFACE` |
| `Circle3` | `CIRCLE` |
| `Ellipse3` | `ELLIPSE` |
| `Line3` | `LINE` |
| `Point3` | `CARTESIAN_POINT` |
| `Vec3` | `DIRECTION` |
| `Frame3` | `AXIS2_PLACEMENT_3D` |

---

## Usage

```toml
[dependencies]
cadcore-step = "0.1"
```

```rust
use cadcore_topo::BRep;
use cadcore_step::brep_to_step;

// brep is populated via cadcore-ops or manually
let step_text: String = brep_to_step(&brep)?;
std::fs::write("output.step", &step_text)?;
```

### Example — sweep + export in one shot

```rust
use cadcore_topo::BRep;
use cadcore_math::Point3;
use cadcore_ops::{sweep_circle_along_polyline, SweepOptions};
use cadcore_step::brep_to_step;

let waypoints = vec![
    Point3::new(0.0, 0.0, 0.0),
    Point3::new(10.0, 0.0, 0.0),
    Point3::new(10.0, 8.0, 0.0),
];

let mut brep = BRep::new();
sweep_circle_along_polyline(&mut brep, &waypoints, 0.2, &SweepOptions::default())?;

let step = brep_to_step(&brep)?;
std::fs::write("rod.step", &step)?;
```

---

## Output format

The generated `.step` file conforms to ISO 10303-21 / AP203 (CONFIGURATION_CONTROLLED_DESIGN protocol) and contains:

- `FILE_DESCRIPTION`, `FILE_NAME`, `FILE_SCHEMA` header
- One `PRODUCT` / `SHAPE_DEFINITION_REPRESENTATION` per solid
- Exact surface and curve geometry — no `B_SPLINE_SURFACE_WITH_KNOTS` approximations

---

## Part of cadcore

For the full kernel (sweep + STEP export via the facade), use [`cadcore`](https://crates.io/crates/cadcore):

```toml
[dependencies]
cadcore = "0.1"
```

---

## License

[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/)

Free for research, education, and non-commercial use.
Commercial licensing: **dmytroyatskovskiy@gmail.com**
