# cadcore-ops

[![crates.io](https://img.shields.io/crates/v/cadcore-ops)](https://crates.io/crates/cadcore-ops)
[![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm%20Noncommercial%201.0.0-blue)](../../LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.78%2B-orange)](https://www.rust-lang.org/)

Geometric operations for the [cadcore](https://crates.io/crates/cadcore) CAD kernel — analytic B-Rep construction without Boolean operations.

---

## Operations

### `sweep_circle_along_polyline`

Sweeps a circular cross-section along a polyline path to produce an exact analytic B-Rep solid.

```
2 end-cap planes  +  (N−1) cylinder segments  +  (N−2) torus fillets  =  O(N) total
```

No Boolean union. No mesh. Every surface is stored as an exact analytic entity.

```toml
[dependencies]
cadcore-ops = "0.1"
```

```rust
use cadcore_topo::BRep;
use cadcore_math::Point3;
use cadcore_ops::{sweep_circle_along_polyline, SweepOptions};

let waypoints = vec![
    Point3::new( 0.0,  0.0, 0.0),
    Point3::new(10.0,  0.0, 0.0),
    Point3::new(10.0,  8.0, 0.0),
    Point3::new( 0.0,  8.0, 0.0),
];

let mut brep = BRep::new();
let id = sweep_circle_along_polyline(
    &mut brep,
    &waypoints,
    0.2,   // radius in mm
    &SweepOptions {
        fillet_corners: true,            // G1-smooth torus at each bend
        name: Some("leg".to_string()),
    },
)?;

let stats = brep.stats();
// 4 waypoints → 3 cylinders + 2 torus fillets + 2 end caps = 7 faces
```

---

## SweepOptions

| Field | Type | Default | Description |
|---|---|---|---|
| `fillet_corners` | `bool` | `true` | Insert a `TorusSurf` fillet at each interior bend |
| `name` | `Option<String>` | `None` | Name attached to the `Solid` entity in the B-Rep |

---

## Complexity

| Operation | OCCT Boolean fuse | cadcore sweep |
|---|---|---|
| N path segments | O(N²) | **O(N)** |
| Output | Mesh (tessellation) | Exact analytic B-Rep |
| STEP export | Round-trip via file | Direct string serialisation |

---

## Part of cadcore

For the full kernel (STEP export included), use the [`cadcore`](https://crates.io/crates/cadcore) facade:

```toml
[dependencies]
cadcore = "0.1"
```

---

## License

[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/)

Free for research, education, and non-commercial use.
Commercial licensing: **dmytroyatskovskiy@gmail.com**
