//! Opt-in wall-clock phase measurements for performance diagnostics.

use std::time::{Duration, Instant};

/// One named, aggregate wall-clock measurement within an algorithm call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PhaseTiming {
    name: &'static str,
    elapsed: Duration,
}

impl PhaseTiming {
    /// Returns the stable machine-readable phase name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        self.name
    }

    /// Returns aggregate monotonic wall time attributed to this phase.
    #[must_use]
    pub const fn elapsed(self) -> Duration {
        self.elapsed
    }

    /// Returns aggregate phase time in nanoseconds.
    #[must_use]
    pub const fn elapsed_ns(self) -> u128 {
        self.elapsed.as_nanos()
    }
}

/// Named phase timings collected during one opt-in profiled algorithm call.
///
/// Repeated phase names are aggregated in first-observation order. The total
/// covers the complete profiled call, while [`Self::unattributed`] captures
/// timer bookkeeping and work deliberately left outside a named phase.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PhaseTimings {
    phases: Box<[PhaseTiming]>,
    total: Duration,
}

impl PhaseTimings {
    /// Returns named phases in first-observation order.
    #[must_use]
    pub const fn phases(&self) -> &[PhaseTiming] {
        &self.phases
    }

    /// Finds one aggregate phase by its stable name.
    #[must_use]
    pub fn phase(&self, name: &str) -> Option<PhaseTiming> {
        self.phases.iter().copied().find(|phase| phase.name == name)
    }

    /// Returns monotonic wall time for the complete profiled call.
    #[must_use]
    pub const fn total(&self) -> Duration {
        self.total
    }

    /// Returns the sum of all named, non-overlapping phase durations.
    #[must_use]
    pub fn attributed(&self) -> Duration {
        self.phases.iter().fold(Duration::ZERO, |total, phase| {
            total.saturating_add(phase.elapsed)
        })
    }

    /// Returns complete-call time not attributed to a named phase.
    #[must_use]
    pub fn unattributed(&self) -> Duration {
        self.total.saturating_sub(self.attributed())
    }
}

/// An algorithm result paired with opt-in monotonic wall-clock phase timings.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Profiled<T> {
    result: T,
    timings: PhaseTimings,
}

impl<T> Profiled<T> {
    /// Creates a profiled result from an algorithm value and completed timings.
    #[doc(hidden)]
    #[must_use]
    pub const fn new(result: T, timings: PhaseTimings) -> Self {
        Self { result, timings }
    }

    /// Returns the unchanged scientific algorithm result.
    #[must_use]
    pub const fn result(&self) -> &T {
        &self.result
    }

    /// Returns the diagnostic phase measurements.
    #[must_use]
    pub const fn timings(&self) -> &PhaseTimings {
        &self.timings
    }

    /// Separates the unchanged result from its diagnostic timings.
    #[must_use]
    pub fn into_parts(self) -> (T, PhaseTimings) {
        (self.result, self.timings)
    }
}

/// Internal monotonic phase collector shared by algorithm crates.
///
/// This type is public only so workspace crates can implement consistent
/// profiled entry points. End users should consume [`Profiled`] and
/// [`PhaseTimings`] instead.
#[doc(hidden)]
#[derive(Debug)]
pub struct PhaseTimingsBuilder {
    started: Instant,
    phases: Vec<PhaseTiming>,
}

impl PhaseTimingsBuilder {
    /// Starts a new complete-call clock with no named phases.
    #[must_use]
    pub fn new() -> Self {
        Self {
            started: Instant::now(),
            phases: Vec::new(),
        }
    }

    /// Measures one closure and aggregates its elapsed time under `name`.
    pub fn measure<T>(&mut self, name: &'static str, operation: impl FnOnce() -> T) -> T {
        let started = Instant::now();
        let result = operation();
        self.record(name, started.elapsed());
        result
    }

    /// Adds an already measured duration to a named phase.
    pub fn record(&mut self, name: &'static str, elapsed: Duration) {
        if let Some(phase) = self.phases.iter_mut().find(|phase| phase.name == name) {
            phase.elapsed = phase.elapsed.saturating_add(elapsed);
        } else {
            self.phases.push(PhaseTiming { name, elapsed });
        }
    }

    /// Finishes the complete-call clock and freezes the aggregate report.
    #[must_use]
    pub fn finish(self) -> PhaseTimings {
        PhaseTimings {
            phases: self.phases.into_boxed_slice(),
            total: self.started.elapsed(),
        }
    }
}

impl Default for PhaseTimingsBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::PhaseTimingsBuilder;
    use std::time::Duration;

    #[test]
    fn repeated_names_are_aggregated_in_first_observation_order() {
        let mut builder = PhaseTimingsBuilder::new();
        builder.record("search", Duration::from_nanos(7));
        builder.record("traceback", Duration::from_nanos(3));
        builder.record("search", Duration::from_nanos(5));
        let timings = builder.finish();

        assert_eq!(timings.phases().len(), 2);
        assert_eq!(timings.phases()[0].name(), "search");
        assert_eq!(timings.phases()[0].elapsed_ns(), 12);
        assert_eq!(timings.phases()[1].name(), "traceback");
        assert_eq!(timings.attributed(), Duration::from_nanos(15));
        assert!(timings.total() >= timings.attributed());
    }
}
