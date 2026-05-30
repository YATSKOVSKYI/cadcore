//! Integration test: a faceted Boolean union must export to a *watertight*
//! STEP `MANIFOLD_SOLID_BREP`.
//!
//! "Watertight" is checked directly on the emitted STEP text: in a closed,
//! oriented manifold shell every `EDGE_CURVE` must be referenced by exactly two
//! `ORIENTED_EDGE`s (the two faces that share it).  A naked edge would be
//! referenced once; a non-manifold edge three or more times.

use std::collections::HashMap;

use cadcore_math::Point3;
use cadcore_ops::{sweep_circle_along_polyline, union_solids, SweepOptions, UnionOptions};
use cadcore_step::brep_to_step;
use cadcore_topo::BRep;

/// Parse the leading `#<n>` token from a string slice.
fn leading_hash(s: &str) -> Option<usize> {
    let s = s.trim_start();
    let s = s.strip_prefix('#')?;
    let digits: String = s.chars().take_while(|c| c.is_ascii_digit()).collect();
    digits.parse().ok()
}

/// Count, for every `EDGE_CURVE` id, how many `ORIENTED_EDGE`s reference it.
fn edge_reference_counts(step: &str) -> (Vec<usize>, HashMap<usize, usize>) {
    let mut edge_ids: Vec<usize> = Vec::new();
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for line in step.lines() {
        let line = line.trim();
        if let Some((lhs, _)) = line.split_once("= EDGE_CURVE(") {
            if let Some(id) = leading_hash(lhs) {
                edge_ids.push(id);
            }
        } else if line.contains("= ORIENTED_EDGE(") {
            // Format: #N = ORIENTED_EDGE('',*,*,#<edge>,.T.);
            if let Some(pos) = line.find("*,*,#") {
                let after = &line[pos + 5..];
                let num: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
                if let Ok(id) = num.parse::<usize>() {
                    *counts.entry(id).or_insert(0) += 1;
                }
            }
        }
    }
    (edge_ids, counts)
}

#[test]
fn perpendicular_cylinder_union_exports_watertight_step() {
    let mut brep = BRep::new();
    let r = 0.5;
    let a = sweep_circle_along_polyline(
        &mut brep,
        &[Point3::new(-2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
        r,
        &SweepOptions::default(),
    )
    .unwrap();
    let b = sweep_circle_along_polyline(
        &mut brep,
        &[Point3::new(0.0, 0.0, -2.0), Point3::new(0.0, 0.0, 2.0)],
        r,
        &SweepOptions::default(),
    )
    .unwrap();

    let _fused = union_solids(&mut brep, a, b, &UnionOptions { facets: 24 }).unwrap();

    let step = brep_to_step(&brep).expect("STEP export");

    // Structural sanity: exactly one closed manifold solid.
    assert_eq!(
        step.matches("MANIFOLD_SOLID_BREP").count(),
        1,
        "expected a single MANIFOLD_SOLID_BREP"
    );
    assert_eq!(
        step.matches("CLOSED_SHELL").count(),
        1,
        "expected a single CLOSED_SHELL"
    );
    assert!(
        step.contains("ADVANCED_FACE"),
        "expected ADVANCED_FACE entities"
    );

    // Watertightness: every EDGE_CURVE shared by exactly two ORIENTED_EDGEs.
    let (edge_ids, counts) = edge_reference_counts(&step);
    assert!(!edge_ids.is_empty(), "no EDGE_CURVE entities emitted");

    let mut naked = 0usize;
    let mut non_manifold = 0usize;
    for id in &edge_ids {
        match counts.get(id).copied().unwrap_or(0) {
            2 => {}
            n if n < 2 => naked += 1,
            _ => non_manifold += 1,
        }
    }
    assert_eq!(
        naked, 0,
        "exported STEP has {naked} naked (free) edges — not watertight"
    );
    assert_eq!(
        non_manifold, 0,
        "exported STEP has {non_manifold} non-manifold edges"
    );
}
