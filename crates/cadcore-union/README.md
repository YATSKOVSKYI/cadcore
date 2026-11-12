# cadcore-union

Booleans for the [cadcore](https://github.com/YATSKOVSKYI/cadcore) CAD kernel, in pure Rust — no OpenCASCADE.

## The API

| Call | What it does |
|---|---|
| `fuse_and_trim` | filaments → one fused solid → planar trims → merged edges, in one call |
| `fuse_filaments_with` | universal fuse of swept filaments (legs, torus elbows, flat caps) into one closed shell |
| `trim_shell` | robust planar trims of an already fused shell, then edge merge |
| `half_space_cut`, `half_space_cut_robust` | one planar cut with section caps; the robust form moves a plane that grazes a feature by 10–40 µm and reports the shift |
| `merge_shell_edges` | fold per-step edges back into whole lines, arcs and curves |
| `LatticeModel` | exact point membership of the fused + trimmed solid, for volume meshers |
| `boolean::union_n` | union of analytic solids (boxes, cylinders, swept filaments, coaxial tori) |
| `validate_shell` | the acceptance gates: manifold, wire winding, edge-to-surface distance |

```rust
use cadcore_geom::Plane3;
use cadcore_math::{Point3, UnitVec3};
use cadcore_union::{fuse_and_trim, Filament, FuseOptions};

let filaments: Vec<Filament> = paths.iter().map(|p| Filament::serpentine(p, 0.275, 0.5)).collect();
let keep = [
    Plane3::from_origin_normal(Point3::new(0.0, 1.0, 0.0), UnitVec3::Y),
    Plane3::from_origin_normal(Point3::new(0.0, 5.0, 0.0), -UnitVec3::Y),
];
let part = fuse_and_trim(&filaments, &FuseOptions::default(), &keep, 1.0e-6);
assert!(part.is_watertight());
```

## How it works

* **Fuse** — for every pair of filaments the boundary curve ∂A ∩ ∂B is marched across members (leg / elbow / cap), with exact joints where it crosses a junction circle. Each curve piece is imprinted on both faces it lies on; every face is then arranged in its own parameter plane and the cells outside all other filaments are kept. Shared edges are identical by construction.
* **Cut** — faces are classified against the plane, the straddling ones are imprinted with their section curve, and the section curves are chained into planar caps. Joints are solved once per edge, so every face and the cap present the same point.
* **Faces** — a tube face that wraps its surface is written as two half faces bounded by ruling lines (cylinder) or equators (torus), which is what CAD importers read reliably.

## Verified

Orthogonal-grid, diamond and honeycomb lattices, fused and trimmed at both terminals: one closed manifold solid in Ansys SpaceClaim and in OpenCASCADE-based readers (FreeCAD); volume within 0.01 % of a point-sampling reference; solved in Ansys Maxwell with separately exported electrode plates.

## Limits

* Filament elbows need a bend radius larger than the tube radius (no spindle tori).
* Two filaments ending flush in the same plane (coplanar overlapping caps) are not united; the trace is reported in `FuseReport::aborted`.
* A filament that crosses itself is not fused with itself (reported in `FuseReport::notes`).
* `union_n` does not handle two tori crossing each other transversally.
