use core::fmt;

use align_core::{Alphabet, Sequence};

use super::{
    ColorMode, PairwiseMetric, PairwiseRunInfo, ReportOptions, Verbosity, sanitized_label,
};

/// One score-, distance-, or penalty-only native result plus its inputs.
#[derive(Clone, Debug)]
pub struct PairwiseValueRun<'a, A: Alphabet, R> {
    first: &'a Sequence<A>,
    second: &'a Sequence<A>,
    result: R,
    info: PairwiseRunInfo,
}

impl<'a, A: Alphabet, R> PairwiseValueRun<'a, A, R> {
    pub(crate) const fn new(
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
        result: R,
        info: PairwiseRunInfo,
    ) -> Self {
        Self {
            first,
            second,
            result,
            info,
        }
    }

    /// Returns the first input sequence.
    #[must_use]
    pub const fn first(&self) -> &'a Sequence<A> {
        self.first
    }

    /// Returns the second input sequence.
    #[must_use]
    pub const fn second(&self) -> &'a Sequence<A> {
        self.second
    }

    /// Returns the untouched native score-only result.
    #[must_use]
    pub const fn result(&self) -> &R {
        &self.result
    }

    /// Consumes the run and returns the untouched native score-only result.
    #[must_use]
    pub fn into_result(self) -> R {
        self.result
    }

    /// Returns algorithm and objective metadata.
    #[must_use]
    pub const fn info(&self) -> &PairwiseRunInfo {
        &self.info
    }

    /// Creates a human-readable objective summary.
    #[must_use]
    pub fn report(&self) -> PairwiseValueReport<'_, 'a, A, R> {
        PairwiseValueReport {
            run: self,
            options: ReportOptions::new(),
            first_label: String::from("First"),
            second_label: String::from("Second"),
        }
    }
}

/// Borrowed configurable score-, distance-, or penalty-only report.
#[derive(Clone, Debug)]
pub struct PairwiseValueReport<'run, 'input, A: Alphabet, R> {
    run: &'run PairwiseValueRun<'input, A, R>,
    options: ReportOptions,
    first_label: String,
    second_label: String,
}

impl<A: Alphabet, R> PairwiseValueReport<'_, '_, A, R> {
    /// Replaces both input labels.
    #[must_use]
    pub fn labels(mut self, first: impl Into<String>, second: impl Into<String>) -> Self {
        self.first_label = sanitized_label(first);
        self.second_label = sanitized_label(second);
        self
    }

    /// Replaces report verbosity.
    #[must_use]
    pub const fn verbosity(mut self, verbosity: Verbosity) -> Self {
        self.options = self.options.with_verbosity(verbosity);
        self
    }

    /// Retains a color preference for terminal-aware writers.
    ///
    /// Stable [`fmt::Display`] output remains colorless.
    #[must_use]
    pub const fn color(mut self, color: ColorMode) -> Self {
        self.options = self.options.with_color(color);
        self
    }

    /// Returns the complete rendering options.
    #[must_use]
    pub const fn options(&self) -> ReportOptions {
        self.options
    }
}

impl<A: Alphabet, R> fmt::Display for PairwiseValueReport<'_, '_, A, R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let info = self.run.info();
        let mut width = 13_usize;
        if self.options.verbosity() == Verbosity::Detailed {
            for detail in info.details() {
                width = width.max(detail.label().chars().count());
            }
        }
        writeln!(formatter, "{:<width$} : {}", "Algorithm", info.algorithm())?;
        writeln!(formatter, "{:<width$} : {}", "Mode", info.mode())?;
        writeln!(formatter, "{:<width$} : {}", "Scoring", info.scoring())?;
        writeln!(
            formatter,
            "{:<width$} : {}={} {}={}",
            "Sequences",
            self.first_label,
            self.run.first.len(),
            self.second_label,
            self.run.second.len()
        )?;
        write_value(formatter, width, info.metric())?;
        if self.options.verbosity() == Verbosity::Detailed {
            for detail in info.details() {
                write!(
                    formatter,
                    "\n{:<width$} : {}",
                    detail.label(),
                    detail.value()
                )?;
            }
        }
        Ok(())
    }
}

fn write_value(
    formatter: &mut fmt::Formatter<'_>,
    width: usize,
    metric: PairwiseMetric,
) -> fmt::Result {
    match metric {
        PairwiseMetric::Score(value) | PairwiseMetric::Penalty(value) => {
            write!(formatter, "{:<width$} : {value}", metric.label())
        }
        PairwiseMetric::Distance(value) => {
            write!(formatter, "{:<width$} : {value}", metric.label())
        }
    }
}
