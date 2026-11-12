//! # cadcore-mesh — from exact B-Rep to triangles (and back from STEP)
//!
//! Everything a viewer, a figure or a mesher needs downstream of the exact
//! model, in the kernel instead of in each application:
//!
//! | Call | What it does |
//! |---|---|
//! | [`mesh_brep`] | exact B-Rep → one triangle mesh per face, to a chord tolerance |
//! | [`import_step`] | any AP203/AP214 STEP text → triangle meshes (trimmed NURBS, analytic surfaces, periodic faces) |
//! | [`read_step_to_brep`] | a STEP file written by cadcore → the kernel B-Rep again |
//! | [`tessellate_brep`] | fast shape-recognising tessellation of template faces (sweeps, caps, plates) |
//! | [`section_outline`] | the closed outline a plane cuts out of a triangulated solid |
//! | [`mesh_volume`] | signed volume of a closed triangle mesh |
//!
//! [`mesh_brep`] is the general path: it serialises the B-Rep with the kernel's
//! own STEP writer and meshes that description, so what is drawn is exactly
//! what a CAD system will read from the exported file — trims, half faces,
//! section caps and all.
//!
//! ```no_run
//! use cadcore_mesh::{mesh_brep, ImportOptions};
//! # let brep = cadcore_topo::BRep::new();
//! let (faces, report) = mesh_brep(&brep, &ImportOptions::default()).unwrap();
//! let triangles: usize = faces.iter().map(|f| f.indices.len() / 3).sum();
//! println!("{} faces, {triangles} triangles", faces.len());
//! # let _ = report;
//! ```

/// The kernel crates under the names the modules of this crate were written
/// against.
pub(crate) mod kernel {
    pub use cadcore_geom as geom;
    pub use cadcore_math as math;
    pub use cadcore_topo as topo;
}

pub mod section;
pub mod step_import;
pub mod step_reader;
pub mod tessellate;

pub use section::{plane_axes, section_outline, SectionOptions, SectionOutline};
pub use step_import::{import_step, mesh_volume, ImportOptions, ImportReport};
pub use step_reader::read_step_to_brep;
pub use tessellate::{tessellate_brep, FaceMesh};

/// Triangulate every face of `brep` to the chord tolerance in `opts`.
///
/// The B-Rep is written with the kernel's STEP writer and that text is meshed
/// by [`import_step`]: the general tessellator understands every face the
/// kernel can emit (trimmed cylinders and tori, planar caps with holes, plates),
/// and the mesh is by construction what an external CAD reads from the export.
///
/// Fails when the B-Rep cannot be serialised or holds no faces.
pub fn mesh_brep(brep: &cadcore_topo::BRep, opts: &ImportOptions) -> Result<(Vec<FaceMesh>, ImportReport), String> {
    let text = cadcore_step::brep_to_step(brep).map_err(|e| format!("STEP serialisation failed: {e}"))?;
    import_step(&text, opts)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A box meshed straight from the kernel B-Rep: closed, exact volume.
    #[test]
    fn mesh_brep_of_a_box_has_its_volume() {
        let mut brep = cadcore_topo::BRep::new();
        cadcore_ops::build_solid_box(&mut brep, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, Some("box".into())).expect("box");
        let (faces, _report) = mesh_brep(&brep, &ImportOptions::default()).unwrap();
        assert_eq!(faces.len(), 6);
        let vol = mesh_volume(&faces).abs();
        assert!((vol - 24.0).abs() < 1e-6, "volume {vol}");
    }
}
