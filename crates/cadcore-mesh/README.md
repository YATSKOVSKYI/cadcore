# cadcore-mesh

Tessellation, STEP import and plane sections for the [cadcore](https://github.com/YATSKOVSKYI/cadcore) CAD kernel.

* `mesh_brep` — exact B-Rep to one triangle mesh per face, to a chord tolerance.
* `import_step` — any AP203/AP214 STEP text to triangle meshes (trimmed NURBS, analytic surfaces, periodic faces).
* `read_step_to_brep` — a STEP file written by cadcore back into the kernel B-Rep.
* `section_outline` — the closed outline a plane cuts out of a triangulated solid.

Pure Rust, no OpenCASCADE.
