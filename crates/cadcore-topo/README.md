# cadcore-topo

**Compile-Time Safe, Arena-Based B-Rep Topology for CAD Kernels.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![crates.io](https://img.shields.io/crates/v/cadcore-topo)](https://crates.io/crates/cadcore-topo)

`cadcore-topo` implements the Boundary Representation (B-Rep) topological database for the `cadcore` kernel. 

Instead of traditional C-style pointer structures (which lead to memory leaks or segmentation faults), `cadcore-topo` stores all topological elements in flat, typed arenas backed by [`slotmap`](https://crates.io/crates/slotmap).

---

## The Topological B-Rep Hierarchy

```
Solid
 └─ Shell  (outer surface + optional inner void shells)
     └─ Face  (bounded region of a surface)
         └─ Loop  (outer boundary loop + optional inner hole loops)
             └─ CoEdge  (directed edge occurrence)
                 └─ Edge  (curve segment bounded by two vertices)
                     └─ Vertex  (point in 3D space)
```

---

## Why cadcore-topo?

1.  **Strict Compile-Time Safety:** Every topological ID is a distinct Rust type (`SolidId`, `FaceId`, `EdgeId`, `VertexId`). The compiler makes it impossible to accidentally pass a `FaceId` into a function expecting an `EdgeId`.
2.  **No Dangling Pointers:** Topological references are verified index lookups, completely eliminating pointer corruption issues common in legacy C++ kernels.
3.  **High-Performance Access:** Flat arrays under the hood ensure cache-friendly lookups and high memory efficiency.
4.  **Serialization Ready:** BRep graphs can be easily serialized/deserialized since they use index IDs rather than memory addresses.

---

## Usage

Add `cadcore-topo` to your `Cargo.toml`:

```toml
[dependencies]
cadcore-topo = "0.1.22"
```

```rust
use cadcore_topo::{BRep, Vertex, VertexId};
use cadcore_math::Point3;

fn main() {
    let mut brep = BRep::new();

    // 1. Add vertices to the topological database
    let v1: VertexId = brep.add_vertex(Vertex {
        point: Point3::new(0.0, 0.0, 0.0),
    });
    let v2: VertexId = brep.add_vertex(Vertex {
        point: Point3::new(10.0, 0.0, 0.0),
    });

    // 2. Query vertices and count stats
    let pt1 = brep[v1].point;
    let stats = brep.stats();

    println!("B-Rep contains {} vertices, {} faces.", stats.vertices, stats.faces);
}
```

---

## License

Licensed under the **MIT License** (see [LICENSE](../../LICENSE)). Free for commercial and non-commercial application.
