# cadcore

![cadcore banner](https://raw.githubusercontent.com/YATSKOVSKYI/cadcore/master/banner.png)

**A High-Performance, Pure-Rust CAD Geometry Kernel.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.78%2B-orange)](https://www.rust-lang.org/)
[![crates.io](https://img.shields.io/crates/v/cadcore)](https://crates.io/crates/cadcore)

`cadcore` is a lightweight, mathematically exact 3D CAD geometry kernel written from the ground up in pure Rust. It is designed to replace heavy C++ legacy kernels (like OpenCASCADE) in modern engineering pipelines, computational manufacturing, and 3D printing.

**Zero C++ dependencies. No OpenCASCADE. No global state. Fully thread-safe & parallelizable.**

---

## Why cadcore?

### 🚀 1. Zero C++ Pain (Pure Rust)
Most CAD kernels in Rust are unsafe FFI bindings to 200MB+ C++ libraries. `cadcore` is 100% Rust. It compiles in seconds, has a binary footprint under 1MB, cross-compiles to any target, and runs natively in **WebAssembly (WASM)** for browser-side CAD.

### ⚡ 2. $O(N)$ Sweep Complexity (The OpenCASCADE Killer)
In traditional CAD kernels, sweeping a circular profile along a toolpath with thousands of lines requires expensive Boolean union operations ($O(N^2)$ complexity) that freeze the CPU or fail on self-intersections. `cadcore` uses a custom analytic sweep solver to build exact B-Rep solids directly in **linear $O(N)$ time** without a single Boolean union.

### 📐 3. Mathematically Watertight
`cadcore` represents 3D shapes using exact analytic surfaces (Planes, Cylinders, Spheres, Tori) and curves (Lines, Circles, Ellipses). There are no mesh approximations or chordal deviation errors in the kernel. Corner joins are smoothly blended with G1-smooth torus fillets or exact miter planes.

### 🔒 4. Modern, Safe B-Rep Topology
We replace C-style raw pointers with an arena-based B-Rep topology using typed stable IDs (`slotmap`). The Rust compiler prevents logic bugs at compile-time (e.g., passing a `FaceId` where an `EdgeId` is expected), and memory is managed safely without segmentation faults.

### 📁 5. Direct, High-Precision STEP Export
Generate industry-standard, high-precision STEP AP203/AP214 files directly from memory without temporary files. The exported files are watertight, manifold closed shells that open flawlessly in **SolidWorks, Autodesk Inventor, Fusion 360, FreeCAD, Rhino, and CATIA**.

---

## Use Cases

*   **3D Printing & Slicing:** Building exact mathematical models from G-code or toolpaths for validation and simulation.
*   **Computational Design:** Generating complex lattice structures, gyroids, and structural frames.
*   **Bioprinting Scaffolds:** Designing micro-channel networks and bio-compatible scaffolds with precise filament sweeping.
*   **WASM CAD Tools:** Building lightweight, browser-based CAD designers and configurators without backend rendering.

---

## Quick Start

Add `cadcore` to your `Cargo.toml`:

```toml
[dependencies]
cadcore = "0.1.22"
```

Create a watertight U-shaped filament rod and export it to a STEP file:

```rust
use cadcore::{
    math::Point3,
    topo::BRep,
    ops::{sweep_circle_along_polyline, SweepOptions},
    step::brep_to_step,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Define a U-shaped filament center-line
    let waypoints = vec![
        Point3::new( 0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(10.0, 8.0, 0.0),
        Point3::new( 0.0, 8.0, 0.0),
    ];

    // 2. Initialize the B-Rep database
    let mut brep = BRep::new();

    // 3. Sweep a circle of radius 0.25mm along the path (inserts torus corner fillets)
    let opts = SweepOptions {
        fillet_corners: true,
        corner_fillet_radius: 0.0, // defaults to profile radius
        name: Some("filament_loop".to_string()),
    };
    sweep_circle_along_polyline(&mut brep, &waypoints, 0.25, &opts)?;

    // 4. Export to a watertight STEP file
    let step_content = brep_to_step(&brep)?;
    std::fs::write("filament.step", &step_content)?;

    println!("Watertight STEP model successfully written to filament.step!");
    Ok(())
}
```

---

## Crate Architecture

The kernel is modularly split into six sub-crates:

```
cadcore/
├── crates/
│   ├── cadcore-math    — Points, vectors, matrices, rigid transforms (zero deps)
│   ├── cadcore-geom    — Exact analytic curves and surfaces
│   ├── cadcore-topo    — Arena B-Rep database: Solid → Shell → Face → Loop → Edge
│   ├── cadcore-ops     — High-level sweep and trimming operations
│   ├── cadcore-step    — Watertight STEP AP203/AP214 writer
│   └── cadcore         — Facade crate re-exporting everything under a single namespace
```

### Dependency Flow

```
cadcore-math (no deps)
     ↑
cadcore-geom
     ↑
cadcore-topo
     ↑
cadcore-ops ──────────────────┐
     ↑                        ↓
cadcore-step          cadcore (facade)
```

---

## Comparison: cadcore vs. OpenCASCADE (OCCT)

| Feature | OpenCASCADE (OCCT) | cadcore |
|---|---|---|
| **Language** | C++ (requires FFI bindings) | 100% Pure Rust |
| **Safety** | Raw pointers, potential segfaults | Compile-time ID safety, safe memory |
| **Performance** | $O(N^2)$ Boolean unions for sweeps | **$O(N)$ analytic sweeps** |
| **Binary Footprint** | ~200 MB prebuilt libraries | **< 1 MB** |
| **WASM Support** | Emscripten (large, complex WASM bloat) | **Native WASM target** |
| **Setup & Build** | Painful multi-step C++ toolchain | `cargo add cadcore` in seconds |

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details. It is free for both non-commercial and commercial usage.

---

## Contributions & Support

Contributions, bug reports, and pull requests are welcome! If you have questions, want to discuss integration into your CAD/CAM pipeline, or have commercial/partnership queries, feel free to reach out to Dmytro Yatskovskyi at [dmytroyatskovskyi@outlook.com](mailto:dmytroyatskovskyi@outlook.com) or open an issue/pull request in the [GitHub Repository](https://github.com/YATSKOVSKYI/cadcore).
