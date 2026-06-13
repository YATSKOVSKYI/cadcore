//! Gate 1 over emitted STEP text.
//!
//! The legacy union's in-memory loops are placeholders: the writer derives
//! the real wires from `FaceExtent` templates and dedupes edges
//! *geometrically* (two distinct `EdgeId`s with identical curves become one
//! `EDGE_CURVE`).  The authoritative manifold check for that path therefore
//! runs on the text the writer actually produced — same rule as ever: every
//! `EDGE_CURVE` referenced by exactly two `ORIENTED_EDGE`s, opposite senses.

use std::collections::HashMap;

/// Check the AP214 edge-sharing invariant on STEP text.
/// Returns one description per violating `EDGE_CURVE` id.
pub fn manifold_check(step: &str) -> Vec<String> {
    // lines look like:  #123 = ORIENTED_EDGE('',*,*,#456,.T.);
    let mut uses: HashMap<u64, Vec<bool>> = HashMap::new();
    for line in step.lines() {
        let Some(pos) = line.find("ORIENTED_EDGE('',*,*,#") else {
            continue;
        };
        let rest = &line[pos + "ORIENTED_EDGE('',*,*,#".len()..];
        let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
        let Ok(ec) = digits.parse::<u64>() else {
            continue;
        };
        let sense = rest[digits.len()..].starts_with(",.T.");
        uses.entry(ec).or_default().push(sense);
    }
    let mut violations = Vec::new();
    for (ec, senses) in &uses {
        let ok = senses.len() == 2 && senses[0] != senses[1];
        if !ok {
            violations.push(format!(
                "EDGE_CURVE #{ec} used {} times (want 2, opposite senses): {senses:?}",
                senses.len()
            ));
        }
    }
    violations.sort();
    violations
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clean_pair_passes() {
        let s = "#1 = ORIENTED_EDGE('',*,*,#10,.T.);\n#2 = ORIENTED_EDGE('',*,*,#10,.F.);\n";
        assert!(manifold_check(s).is_empty());
    }

    #[test]
    fn single_use_flagged() {
        let s = "#1 = ORIENTED_EDGE('',*,*,#10,.T.);\n";
        let v = manifold_check(s);
        assert_eq!(v.len(), 1);
        assert!(v[0].contains("used 1 times"));
    }

    #[test]
    fn same_sense_flagged() {
        let s = "#1 = ORIENTED_EDGE('',*,*,#10,.T.);\n#2 = ORIENTED_EDGE('',*,*,#10,.T.);\n";
        let v = manifold_check(s);
        assert_eq!(v.len(), 1);
        assert!(v[0].contains("[true, true]"));
    }

    #[test]
    fn triple_use_flagged() {
        let s = "#1 = ORIENTED_EDGE('',*,*,#10,.T.);#9 = X;\n#2 = ORIENTED_EDGE('',*,*,#10,.F.);\n#3 = ORIENTED_EDGE('',*,*,#10,.T.);\n";
        let v = manifold_check(s);
        assert_eq!(v.len(), 1);
        assert!(v[0].contains("used 3 times"));
    }
}
