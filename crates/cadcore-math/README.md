# cadcore-math

**Type-Safe, High-Performance 3D Math Primitives for CAD Kernels.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![crates.io](https://img.shields.io/crates/v/cadcore-math)](https://crates.io/crates/cadcore-math)

`cadcore-math` provides the mathematical foundation for the `cadcore` CAD kernel. It is designed to be extremely fast, robust, and zero-allocation, utilizing const-generic arithmetic.

**Zero external dependencies. Pure Rust.**

---

## Key Features

*   **Type Safety at the Compiler Level:** Points (`Point3`) and vectors (`Vec3` and `UnitVec3`) are separate types. The compiler will not let you accidentally add two coordinates or confuse a point location with a direction vector.
*   **Enforced Unit Vectors:** `UnitVec3` guarantees length `1.0` at construction, preventing normalization errors during geometric computations.
*   **Orthonormal Frames:** Built-in `Frame3` makes local-to-world coordinate transformations simple and mathematically exact.
*   **High Performance:** Fully optimized for double-precision float (`f64`) operations.

---

## Math Types

| Type | Description |
|---|---|
| `Point3` | A location in 3D space (measured in millimeters) |
| `Vec3` | A free vector (direction and magnitude) |
| `UnitVec3` | A normalized unit vector (guaranteed direction) |
| `Mat3` | Column-major 3x3 matrix (rotations, linear maps) |
| `Frame3` | Right-handed orthonormal frame (origin + Z, X, Y axes) |
| `Transform3` | Rigid-body transform (rotation + translation) |
| `Interval` | A closed real interval `[lo, hi]` for parameter bounds |

---

## Usage

Add `cadcore-math` to your `Cargo.toml`:

```toml
[dependencies]
cadcore-math = "0.1.22"
```

```rust
use cadcore_math::{Point3, Vec3, UnitVec3, Frame3};

fn main() {
    // Define an origin and a Z-axis direction vector
    let origin = Point3::new(1.0, 2.0, 3.0);
    let z_dir  = UnitVec3::try_from_vec(Vec3::new(0.0, 0.0, 1.0)).unwrap();

    // Construct an orthonormal coordinate frame
    let frame  = Frame3::from_origin_z(origin, z_dir);

    // Transform points between local and world coordinates
    let world_pt = Point3::new(4.0, 5.0, 6.0);
    let local_pt = frame.to_local_point(world_pt);
    let back_pt  = frame.to_world_point(local_pt);

    assert!((world_pt - back_pt).length() < 1e-10);
}
```

---

## License

Licensed under the **MIT License** (see [LICENSE](../../LICENSE)). Free for commercial and non-commercial application.
