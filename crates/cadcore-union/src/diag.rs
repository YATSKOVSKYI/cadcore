//! Structured diagnostics.
//!
//! The legacy engine reports through `eprintln!` behind a zoo of environment
//! variables; reconstructing what happened means grepping stderr.  Here every
//! noteworthy fact is a typed [`Event`] with a severity, collected in
//! [`Diagnostics`] and queryable by the caller (the Tauri app can surface
//! warnings in the UI; tests can assert on exact counters).

use std::collections::BTreeMap;

/// Event severity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    /// Counters / traces; no action needed.
    Info,
    /// The engine handled it, but the input geometry is marginal
    /// (near-tangency, sliver collapsed, weld at the large end of the
    /// tolerance).
    Warning,
    /// A trim or stage was skipped; output is valid but less merged than
    /// requested.
    Degraded,
    /// Validation gate failure — the output must not be shipped.
    Error,
}

/// A single diagnostic event.
#[derive(Debug, Clone)]
pub struct Event {
    /// Severity tier.
    pub severity: Severity,
    /// Stable machine-readable code, e.g. `"chain_fail"`, `"gate.manifold"`.
    pub code: &'static str,
    /// Human-readable details.
    pub message: String,
}

/// Diagnostics sink: events plus named counters.
#[derive(Debug, Default)]
pub struct Diagnostics {
    events: Vec<Event>,
    counters: BTreeMap<&'static str, u64>,
}

impl Diagnostics {
    /// New empty sink.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record an event.
    pub fn event(&mut self, severity: Severity, code: &'static str, message: impl Into<String>) {
        self.events.push(Event {
            severity,
            code,
            message: message.into(),
        });
    }

    /// Bump a counter.
    pub fn count(&mut self, code: &'static str) {
        *self.counters.entry(code).or_insert(0) += 1;
    }

    /// Bump a counter by `n`.
    pub fn count_n(&mut self, code: &'static str, n: u64) {
        *self.counters.entry(code).or_insert(0) += n;
    }

    /// Current counter value.
    pub fn counter(&self, code: &str) -> u64 {
        self.counters.get(code).copied().unwrap_or(0)
    }

    /// All events (in order of occurrence).
    pub fn events(&self) -> &[Event] {
        &self.events
    }

    /// Highest severity recorded (None when empty).
    pub fn worst(&self) -> Option<Severity> {
        self.events.iter().map(|e| e.severity).max()
    }

    /// True when no event at or above the given severity exists.
    pub fn clean_at(&self, level: Severity) -> bool {
        self.events.iter().all(|e| e.severity < level)
    }

    /// All counters.
    pub fn counters(&self) -> &BTreeMap<&'static str, u64> {
        &self.counters
    }

    /// One-line summary (counters + worst severity) for logs.
    pub fn summary(&self) -> String {
        let worst = self
            .worst()
            .map_or("clean".to_string(), |s| format!("{s:?}"));
        let counters: Vec<String> = self
            .counters
            .iter()
            .map(|(k, v)| format!("{k}={v}"))
            .collect();
        format!("[union] worst={worst} {}", counters.join(" "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn worst_and_clean() {
        let mut d = Diagnostics::new();
        assert!(d.clean_at(Severity::Warning));
        d.event(Severity::Info, "x", "info");
        d.event(Severity::Degraded, "y", "skipped");
        assert_eq!(d.worst(), Some(Severity::Degraded));
        assert!(d.clean_at(Severity::Error));
        assert!(!d.clean_at(Severity::Degraded));
    }

    #[test]
    fn counters_accumulate() {
        let mut d = Diagnostics::new();
        d.count("pairs");
        d.count("pairs");
        d.count_n("trims", 5);
        assert_eq!(d.counter("pairs"), 2);
        assert_eq!(d.counter("trims"), 5);
        assert_eq!(d.counter("absent"), 0);
    }
}
