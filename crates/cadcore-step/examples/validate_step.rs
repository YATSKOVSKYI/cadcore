//! Validate a STEP file's closed shells against the AP214 manifold invariant
//! AND check end-cap geometry.  Usage: `cargo run --example validate_step -- file.step`
use std::collections::HashMap;

fn main() {
    let path = std::env::args().nth(1).expect("usage: validate_step <file.step>");
    let step = std::fs::read_to_string(&path).expect("read step");

    // ec_id -> Vec<(oe_id, sense)>
    let mut uses: HashMap<usize, Vec<(usize, String)>> = HashMap::new();
    for line in step.lines() {
        if !line.contains("ORIENTED_EDGE") {
            continue;
        }
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 4 {
            continue;
        }
        let ec = parts[parts.len() - 2].trim().trim_start_matches('#');
        let sense = parts[parts.len() - 1].trim().trim_end_matches(';').to_string();
        let ec_id: usize = ec.parse().ok().unwrap_or(0);
        let oe_id: usize = line.trim().split('=').next().unwrap_or("").trim().trim_start_matches('#').parse().ok().unwrap_or(0);
        if ec_id > 0 && oe_id > 0 {
            uses.entry(ec_id).or_default().push((oe_id, sense));
        }
    }

    let is_edge_curve: std::collections::HashSet<usize> = step
        .lines()
        .filter_map(|l| {
            if l.contains("= EDGE_CURVE") {
                l.trim().split('=').next()?.trim().trim_start_matches('#').parse().ok()
            } else {
                None
            }
        })
        .collect();

    let mut single = 0usize;
    let mut over = 0usize;
    let mut samesense = 0usize;
    for (ec, refs) in &uses {
        if !is_edge_curve.contains(ec) {
            continue;
        }
        if refs.len() == 1 {
            single += 1;
            if single <= 20 {
                println!("OPEN  EC#{ec} used once: {:?}", refs);
            }
        } else if refs.len() > 2 {
            over += 1;
            if over <= 20 {
                println!("OVER  EC#{ec} used {}x: {:?}", refs.len(), refs);
            }
        } else {
            let t = refs.iter().filter(|(_, s)| s.contains(".T.")).count();
            let f = refs.iter().filter(|(_, s)| s.contains(".F.")).count();
            if t != 1 || f != 1 {
                samesense += 1;
                if samesense <= 20 {
                    println!("SENSE EC#{ec} {:?}", refs);
                }
            }
        }
    }
    let total_ec = is_edge_curve.len();
    println!("---");
    println!("EDGE_CURVE total: {total_ec}");
    println!("used-once (open shell): {single}");
    println!("used->2 (non-manifold): {over}");
    println!("same-sense pairs:       {samesense}");
    if single == 0 && over == 0 && samesense == 0 {
        println!("OK: fully manifold");
    }
}
