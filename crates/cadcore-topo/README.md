# cadcore-topo

[![crates.io](https://img.shields.io/crates/v/cadcore-topo)](https://crates.io/crates/cadcore-topo)
[![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/license-PolyForm%20Noncommercial%201.0.0-blue)](../../LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.78%2B-orange)](https://www.rust-lang.org/)

Arena-based B-Rep topology for the [cadcore](https://crates.io/crates/cadcore) CAD kernel.

Provides the full **Solid → Shell → Face → Loop → CoEdge → Edge → Vertex** hierarchy with typed stable IDs backed by [`slotmap`](https://crates.io/crates/slotmap).

---

## Topology hierarchy

```
Solid
 └─ Shell  (outer surface + optional inner voids)
     └─ Face  (bounded region of a surface)
         └─ Loop  (outer boundary + optional inner holes)
             └─ CoEdge  (directed use of an edge)
                 └─ Edge  (curve segment between two vertices)
                     └─ Vertex  (3-D point)
```

Every entity is stored in its own typed arena inside `BRep`. IDs are distinct Rust types — the compiler prevents passing a `FaceId` where a `VertexId` is expected.

---

## Usage

```toml
[dependencies]
cadcore-topo = "0.1"
```

```rust
use cadcore_topo::{BRep, Vertex, VertexId};
use cadcore_math::Point3;

let mut brep = BRep::new();

// Add a vertex
let v: VertexId = brep.add_vertex(Vertex {
    point: Point3::new(1.0, 2.0, 3.0),
});

// Query a vertex by ID
let pt = brep[v].point;

// Topology statistics
let stats = brep.stats();
println!("vertices={} edges={} faces={} solids={}",
    stats.vertices, stats.edges, stats.faces, stats.solids);
```

In practice you build `BRep` via high-level ops (see [`cadcore-ops`](https://crates.io/crates/cadcore-ops)) rather than manually.

---

## Part of cadcore

For the full kernel (sweep + STEP export), use the [`cadcore`](https://crates.io/crates/cadcore) facade:

```toml
[dependencies]
cadcore = "0.1"
```

---

## License

[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/)

Free for research, education, and non-commercial use.
Commercial licensing: **dmytroyatskovskiy@gmail.com**
