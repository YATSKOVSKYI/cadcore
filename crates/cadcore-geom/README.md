# cadcore-geom

[![crates.io](https://img.shields.io/crates/v/cadcore-geom)](https://crates.io/crates/cadcore-geom)
[![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm%20Noncommercial%201.0.0-blue)](../../LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.78%2B-orange)](https://www.rust-lang.org/)

Analytic geometry primitives for the [cadcore](https://crates.io/crates/cadcore) CAD kernel — exact curves and surfaces that map 1-to-1 to STEP AP203 entities.

---

## Curves

| Type | STEP entity | Description |
|---|---|---|
| `Line3` | `LINE` | Infinite directed line |
| `Circle3` | `CIRCLE` | Planar circle with radius |
| `Ellipse3` | `ELLIPSE` | Planar ellipse (semi-major / semi-minor) |
| `BezierCubic` | — | Cubic Bézier (approximation) |

## Surfaces

| Type | STEP entity | Description |
|---|---|---|
| `Plane3` | `PLANE` | Infinite plane |
| `CylSurf` | `CYLINDRICAL_SURFACE` | Right circular cylinder |
| `SphereSurf` | `SPHERICAL_SURFACE` | Full sphere |
| `TorusSurf` | `TOROIDAL_SURFACE` | Ring torus — used as corner fillets in sweeps |
| `ConeSurf` | `CONICAL_SURFACE` | Right circular cone |

All types are `Copy` and parameterised by `f64`. No heap allocation.

---

## Usage

```toml
[dependencies]
cadcore-geom = "0.1"
```

```rust
use cadcore_geom::{CylSurf, TorusSurf, Plane3};
use cadcore_math::{Point3, UnitVec3, Frame3};

// Cylindrical surface: axis along Z, radius 0.2 mm
let cyl = CylSurf::new(Point3::ORIGIN, UnitVec3::Z, 0.2);
let p   = cyl.point_at(0.0, 5.0);   // theta = 0 rad, z = 5 mm

// Toroidal surface: major radius 1.0, tube radius 0.2
let torus = TorusSurf::new(Frame3::IDENTITY, 1.0, 0.2);

// Plane through origin with Z normal
let plane = Plane3::new(Frame3::IDENTITY);
let foot  = plane.project(Point3::new(1.0, 2.0, 3.0));  // (1, 2, 0)
```

---

## Part of cadcore

This crate is a building block. For the full kernel (B-Rep + sweep + STEP export), use the [`cadcore`](https://crates.io/crates/cadcore) facade:

```toml
[dependencies]
cadcore = "0.1"
```

---

## License

[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/)

Free for research, education, and non-commercial use.
Commercial licensing: **dmytroyatskovskiy@gmail.com**
