//! Example: build a scaffold filament rod and export it to STEP.
//!
//! Run with:
//! ```
//! cargo run --example scaffold_rod
//! ```

use cadcore::{
    math::Point3,
    ops::{sweep_circle_along_polyline, SweepOptions},
    step::brep_to_step,
    topo::BRep,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Build a simple U-shape path (like one leg of a scaffold)
    let waypoints = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(10.0, 0.0, 0.0),
        Point3::new(10.0, 8.0, 0.0),
        Point3::new(0.0, 8.0, 0.0),
    ];

    let filament_radius = 0.2_f64; // 0.4 mm diameter filament

    let mut brep = BRep::new();
    let opts = SweepOptions {
        fillet_corners: true,
        name: Some("scaffold_leg".to_string()),
    };
    let _solid_id = sweep_circle_along_polyline(&mut brep, &waypoints, filament_radius, &opts)?;

    let stats = brep.stats();
    println!(
        "B-Rep stats: {} faces, {} shells, {} solids",
        stats.faces, stats.shells, stats.solids
    );

    let step = brep_to_step(&brep)?;
    let path = "scaffold_leg.step";
    std::fs::write(path, &step)?;
    println!("Wrote {path} ({} bytes)", step.len());

    // Verify the STEP file contains the expected surface types
    assert!(step.contains("CYLINDRICAL_SURFACE"), "missing cylinder");
    assert!(step.contains("TOROIDAL_SURFACE"), "missing torus fillets");
    assert!(step.contains("PLANE"), "missing end caps");
    println!("All expected surfaces present in STEP output ✓");

    Ok(())
}
