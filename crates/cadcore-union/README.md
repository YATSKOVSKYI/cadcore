# cadcore-union — the Boolean-union engine

Industrial-grade Boolean union for cadcore, designed as a **separate engine**
(the way Parasolid sits under SolidWorks).  The legacy single-file implementation
(`cadcore-ops/src/union.rs`) stays untouched and in production; this crate grows
alongside it and takes over stage by stage.

## Why this exists

The legacy union was grown case-by-case (woodpile crossings → mid-U elbows →
corner U×U → butt-end bites).  It works — SolidWorks- and Ansys-confirmed — but
it is one 4 000-line file, its diagnostics are env-var printfs, and every new
geometric situation (varying layer heights, shifted filaments, oblique
crossings, tangencies) needs hand-written passes.  This crate is the
generalisation: one staged pipeline that handles *configurations*, not cases.

## Hard-won invariants (violate these and SolidWorks breaks)

These cost weeks to learn.  They are law:

1. **Trim curves must lie on the EMITTED geometry to ≤ 1 µm.**
   SolidWorks knits at ~1 µm.  Curves computed against sampled/marched
   surfaces (SSI) sit 1–5 µm off the analytic surfaces the writer emits, and
   SW silently drops the affected faces ("filaments vanish").  OpenCASCADE
   forgives this — it is NOT an oracle for SW.  Refinement onto emitted
   geometry is a first-class pipeline stage (`geom::refine`).
2. **Every joint is a triple point.**  Where two trim pieces meet on a face
   boundary, the junction point must be refined onto ALL incident surfaces and
   then used as the *single source of truth* everywhere — runs, arcs, strips.
   Two faces partitioning a shared circle from slightly different endpoints
   produce unshared arcs → open shell.
3. **Every edge is used exactly twice, opposite senses** (AP214 manifold
   invariant, see `cadcore/CLAUDE.md`).
4. **Wires must be physically wound** (outer CCW, holes CW in uv as seen along
   the face normal): SW ignores `FACE_BOUND` orientation flags.
5. **Wires must not intersect each other or themselves in uv** — SW drops the
   face; OCC does not care.
6. **Validate in-process, before writing STEP.**  The four gates
   (`validate::*`) are the Rust ports of the external Python auditors that
   exactly reproduced SolidWorks' accept/reject behaviour.

## Architecture

```
lib.rs            public API: plan() / run() over cadcore-topo BReps
config.rs         UnionConfig: tolerance ladder, pass toggles, diag level
tolerance.rs      Tolerances: model / refine / weld / knit / sliver
diag.rs           structured diagnostics (events + counters), no env-var printf
geom/
  classify.rs     surface & surface-pair taxonomy (incl. tangency margins —
                  near-tangent crossings at HL→2r are the hardest inputs)
  refine.rs       projection oracle: AnalyticSurface::project, cyclic point /
                  joint / curve refinement onto emitted geometry   [LIVE]
  intersect.rs    intersection oracle: predictor–corrector tracer, every
                  point exact on both surfaces; tangency reported   [LIVE]
  composite.rs    whole filament (legs+elbows) as ONE projectable;
                  members, junctions, containment, loop discovery   [LIVE]
arrange/
  registry.rs     pre-split: every loop split ONCE at all junction
                  crossings of both solids; joints refined onto
                  circle ∩ other surface = single source of truth   [LIVE]
  domain.rs       per-face uv (cylinder band / torus patch), lift   [LIVE]
  cells.rs        DCEL arrangement + keep/drop classification       [LIVE]
validate/
  manifold.rs     gate 1: edge 2-use (BRep) + step_text (emitted)   [LIVE]
  distance.rs     gate 4: edge-to-surface deviation < knit tol      [LIVE]
  wires.rs        gates 2+3: uv wire closure/winding/intersection   [LIVE]
pipeline.rs       staged orchestration: Collect → Classify → Intersect →
                  Refine → Arrange → Stitch → Validate
```

### Migration plan

| Phase | What moves here | Status |
|-------|-----------------|--------|
| P0 | skeleton, tolerances, diagnostics, refine oracle, gates 1+4 | this commit |
| P1 | gate 2+3 (uv wire audit) in Rust; run all gates after legacy union | next |
| P2 | analytic intersection oracle (cyl×cyl any angle, cyl×torus, torus×torus, tangency-aware) | |
| P3 | per-face uv arrangement (DCEL from `cadcore-geom::arrangement`) + pre-split registry: every intersection curve split at ALL boundary/seam crossings once, shared by both faces | |
| P4 | pipeline replaces legacy `union.rs` behind a config flag, case matrix in CI (HL ladder, offsets, oblique) | |

### Testing doctrine

Small models, full feature sets (the 5.4 mm mini reproduces everything the
15 mm part does at 1/20 the bytes).  Every gate runs on every artefact in CI.
A configuration matrix (HL 0.35–0.50 × spacing × offsets) is the regression
suite — single golden files are not enough.
