# cadcore STEP Writer — Critical Rules

## AP214 Edge-sharing invariant (DO NOT BREAK)

Every `EDGE_CURVE` in a `CLOSED_SHELL` must be used by **exactly two** adjacent
faces, one with `ORIENTED_EDGE .T.` and one with `.F.`.  Violating this makes
SpaceClaim (and strict AP214 validators) report an open shell.

### Global vs per-face cache

The `Ctx` caches (`point_cache`, `vertex_cache`, `edge_cache`) are **cleared
once per shell** (at `emit_advanced_faces` shell boundary), NOT per face.
Per-face clearing would prevent valid edge sharing.

### Circle / Ellipse boundary sharing (cap ↔ cylinder)

`emit_circle_bound` and `emit_ellipse_bound` use `normal_key()` which
**normalises the sign** of the normal vector so that `+Y` and `−Y` map to the
same cache key.  This lets a disk cap face (plane normal = `−dirs[0]`) and the
adjacent cylinder's start boundary (circle normal = `+dirs[0]`) share the same
`EDGE_CURVE`.

**Critical**: on a cache hit the sense must be **flipped** relative to the
first user.  The cache second element stores the first user's `orient` (0=`.F.`,
1=`.T.`).  Second user always gets the opposite:

```rust
// cache hit
let sense = if first_orient != 0 { ".F." } else { ".T." };
// cache miss (first user)
ctx.edge_cache.insert(key, (ec_id, orient as usize));
let sense = if orient { ".T." } else { ".F." };
```

Breaking this rule causes start caps (or end caps) to display as inverted /
inside-out in any CAD tool that interprets face normals.

### Rounded box / electrode = a prism — EVERY edge is shared (lines too)

`build_solid_rounded_box_xz` (the silver electrode) is a rounded-rectangle XZ
profile extruded along Y.  It has 10 faces: 2 `RoundedRectCap` end caps, 4
`FaceExtent::Polygon` flat sides, 4 `CylinderArcFace` quarter-cylinder corners.

For a watertight manifold **every** edge — straight *and* arc — must be shared
by exactly two faces.  The emitters funnel all edges through two helpers in
`writer.rs`:

* `push_shared_line` → `StepCurveKey::Line { v1: min, v2: max }` (keyed by the
  two vertex ids).
* `push_shared_arc`  → `StepCurveKey::ArcEdge` (direction-agnostic canonical
  key: endpoints sorted, plus centre/axis/xref).

Both record the first user's start vertex in the cache; the second user gets the
opposite `ORIENTED_EDGE` sense automatically (`orig_start == current_start`).

**The invariant that makes this work is consistent OUTWARD winding.**  Every
face loop is wound CCW-as-seen-from-outside:

* `−Y` front cap: profile **forward** (segments 0→7).
* `+Y` back cap: profile **backward** (segments 7→0).
* each `CylinderArcFace`: `B_i → T_i → T_e → B_e` (up, top-arc, down,
  bottom-arc-backward) — radially-outward normal.

With consistent outward winding, the two faces meeting at any edge always
traverse it in opposite vertex order, so the cache yields exactly one `.T.` +
one `.F.` per edge.

> ⚠️ HISTORICAL BUG (do not reintroduce): an earlier version deliberately did
> NOT cache the straight longitudinal / cap edges, on the false belief that
> adjacent faces traversed them in the same direction.  That belief was a
> symptom of a *backwards `CylinderArcFace` loop* (it was wound inward).  The
> result was 32 single-use `EDGE_CURVE`s → SpaceClaim saw an **open shell** and
> rendered the electrode as a transparent/broken part with protruding cap edges.
> Fix = correct the winding + share ALL edges, never "skip caching to dodge a
> same-sense clash".

For an arc EDGE_CURVE the `same_sense` flag (5th field) is derived from geometry
in `push_shared_arc`: `.T.` when `p_from → p_to` runs CCW about the circle axis,
else `.F.` (computed via `atan2` on the `(xref, axis × xref)` basis).  Never
hard-code it per corner — that does not generalise across Y-levels / corners.

### Half-space cut (trim) must preserve the UNCUT-end boundary

`half_space_cut_brep` (cadcore-ops `boolean.rs`) rebuilds each cylinder face
independently against the cut plane.  When a cylinder is **axially truncated**
(`CylinderOutcome::AxialCut`), only the end that actually meets the plane may be
replaced with a fresh plain `Circle` (it pairs with the new flat disk cap).  The
**other** (untouched) end MUST keep its *original* `FaceBoundary` — which is
usually a **miter ellipse** shared with a surviving neighbour (e.g. the
serpentine connector that joins this leg).

> ⚠️ HISTORICAL BUG (do not reintroduce): `AxialTruncated` used to emit a plain
> circle at *both* ends.  The kept neighbour still carried its miter ellipse at
> the shared junction, so the two faces no longer shared one EDGE_CURVE → that
> edge was used **once** → open shell.  CAD tools (SpaceClaim) then rendered the
> filament as **cut short / not drawn at exactly one end** — and because layers
> are laid orthogonally, the broken end alternated sides per layer.  Fix =
> `FaceTemplate::AxialTruncated` carries `start`/`end: FaceBoundary` and the
> cut-end-only gets a new circle; the uncut end is `start.clone()`/`end.clone()`.

Covered by `half_space_cut_serpentine_stays_manifold` (writer.rs tests) and the
`two_layer_caps` integration test in `gcode-viewer-rs/crates/brep`.

### Regression tests (writer.rs `mod tests`)

`check_ap214_manifold` parses the STEP output and asserts every `EDGE_CURVE` is
referenced by exactly two `ORIENTED_EDGE`s with opposite senses.  Covered by:

* `straight_cylinder_caps_share_edge_curve` — swept filament + caps.
* `rounded_box_xz_is_manifold` — the electrode (10 faces, 4 cylinders).
* `rounded_box_xz_arc_edges_not_duplicated` — no arc used > 2×.

Run `cargo test -p cadcore-step` after ANY change to the box/cap/cylinder
emitters.
