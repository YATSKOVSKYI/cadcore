# cadcore-math

[![crates.io](https://img.shields.io/crates/v/cadcore-math)](https://crates.io/crates/cadcore-math)
[![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm%20Noncommercial%201.0.0-blue)](../../LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.78%2B-orange)](https://www.rust-lang.org/)

Fundamental math primitives for the [cadcore](https://crates.io/crates/cadcore) CAD kernel.

**Zero external dependencies.** Everything in `cadcore` is built on top of this crate.

---

## Types

| Type | Description |
|---|---|
| `Point3` | Location in 3-D space (mm) |
| `Vec3` | Free vector — direction + magnitude |
| `UnitVec3` | Unit-length vector, enforced at construction |
| `Mat3` | 3×3 column-major matrix (rotation, linear maps) |
| `Frame3` | Right-handed orthonormal frame (origin + 3 axes) |
| `Transform3` | Rigid-body transform: rotation + translation |
| `Interval` | Closed real interval `[lo, hi]` |

Constants: `EPS`, `PI`, `TAU`.

---

## Usage

```toml
[dependencies]
cadcore-math = "0.1"
```

```rust
use cadcore_math::{Point3, Vec3, UnitVec3, Mat3, Frame3};

// Build a right-handed frame from an origin and a Z axis
let origin = Point3::new(1.0, 2.0, 0.0);
let axis   = UnitVec3::try_from_vec(Vec3::new(0.0, 0.0, 1.0)).unwrap();
let frame  = Frame3::from_origin_z(origin, axis);

// Vector arithmetic
let a = Vec3::new(1.0, 0.0, 0.0);
let b = Vec3::new(0.0, 1.0, 0.0);
let c = a.cross(b);       // (0, 0, 1)
let d = a.dot(b);         // 0.0

// 3×3 rotation matrix
let rot = Mat3::rotation_z(std::f64::consts::FRAC_PI_2);  // 90° around Z
let v   = rot * Vec3::new(1.0, 0.0, 0.0);                 // (0, 1, 0)

// Interval utilities
use cadcore_math::Interval;
let i = Interval::new(0.0, 10.0);
assert!(i.contains(5.0));
let (lo, hi) = i.lerp_pair(0.25, 0.75);  // 2.5, 7.5
```

---

## Part of cadcore

This crate is a building block. If you need the full kernel (sweep, STEP export), use the [`cadcore`](https://crates.io/crates/cadcore) facade instead:

```toml
[dependencies]
cadcore = "0.1"
```

---

## License

[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/)

Free for research, education, and non-commercial use.
Commercial licensing: **dmytroyatskovskiy@gmail.com**
