use core::fmt;

use align_core::{Alphabet, Sequence};
use align_fast::MyersMatch;

use super::{ColorMode, CoordinateStyle, ReportOptions, Verbosity};

/// One exact Myers distance result plus its input sequences.
#[derive(Clone, Debug)]
pub struct MyersDistanceRun<'a, A: Alphabet> {
    pattern: &'a Sequence<A>,
    text: &'a Sequence<A>,
    distance: usize,
}

impl<'a, A: Alphabet> MyersDistanceRun<'a, A> {
    pub(crate) const fn new(
        pattern: &'a Sequence<A>,
        text: &'a Sequence<A>,
        distance: usize,
    ) -> Self {
        Self {
            pattern,
            text,
            distance,
        }
    }

    /// Returns the pattern.
    #[must_use]
    pub const fn pattern(&self) -> &'a Sequence<A> {
        self.pattern
    }

    /// Returns the searched text.
    #[must_use]
    pub const fn text(&self) -> &'a Sequence<A> {
        self.text
    }

    /// Returns the exact global unit edit distance.
    #[must_use]
    pub const fn distance(&self) -> usize {
        self.distance
    }

    /// Creates a human-readable summary.
    #[must_use]
    pub const fn report(&self) -> MyersDistanceReport<'_, 'a, A> {
        MyersDistanceReport {
            run: self,
            options: ReportOptions::new(),
        }
    }
}

/// Borrowed Myers distance report.
#[derive(Clone, Copy, Debug)]
pub struct MyersDistanceReport<'run, 'sequence, A: Alphabet> {
    run: &'run MyersDistanceRun<'sequence, A>,
    options: ReportOptions,
}

impl<A: Alphabet> MyersDistanceReport<'_, '_, A> {
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

impl<A: Alphabet> fmt::Display for MyersDistanceReport<'_, '_, A> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(formatter, "Algorithm      : Myers bit-vector")?;
        writeln!(formatter, "Mode           : global unit edit distance")?;
        writeln!(formatter, "Pattern length : {}", self.run.pattern.len())?;
        writeln!(formatter, "Text length    : {}", self.run.text.len())?;
        write!(formatter, "Distance       : {}", self.run.distance)
    }
}

/// One exact Myers approximate-search result plus its inputs and threshold.
#[derive(Clone, Debug)]
pub struct MyersSearchRun<'a, A: Alphabet> {
    pattern: &'a Sequence<A>,
    text: &'a Sequence<A>,
    maximum_distance: usize,
    matches: Vec<MyersMatch>,
}

impl<'a, A: Alphabet> MyersSearchRun<'a, A> {
    pub(crate) const fn new(
        pattern: &'a Sequence<A>,
        text: &'a Sequence<A>,
        maximum_distance: usize,
        matches: Vec<MyersMatch>,
    ) -> Self {
        Self {
            pattern,
            text,
            maximum_distance,
            matches,
        }
    }

    /// Returns the pattern.
    #[must_use]
    pub const fn pattern(&self) -> &'a Sequence<A> {
        self.pattern
    }

    /// Returns the searched text.
    #[must_use]
    pub const fn text(&self) -> &'a Sequence<A> {
        self.text
    }

    /// Returns the inclusive maximum accepted unit edit distance.
    #[must_use]
    pub const fn maximum_distance(&self) -> usize {
        self.maximum_distance
    }

    /// Returns matches in stable increasing endpoint order.
    #[must_use]
    pub fn matches(&self) -> &[MyersMatch] {
        &self.matches
    }

    /// Consumes the run and returns the native match vector.
    #[must_use]
    pub fn into_matches(self) -> Vec<MyersMatch> {
        self.matches
    }

    /// Creates a normal human-readable match table.
    #[must_use]
    pub const fn report(&self) -> MyersSearchReport<'_, 'a, A> {
        MyersSearchReport {
            run: self,
            options: ReportOptions::new(),
        }
    }
}

/// Borrowed configurable Myers match-table report.
#[derive(Clone, Copy, Debug)]
pub struct MyersSearchReport<'run, 'sequence, A: Alphabet> {
    run: &'run MyersSearchRun<'sequence, A>,
    options: ReportOptions,
}

impl<A: Alphabet> MyersSearchReport<'_, '_, A> {
    /// Replaces report verbosity.
    #[must_use]
    pub const fn verbosity(mut self, verbosity: Verbosity) -> Self {
        self.options = self.options.with_verbosity(verbosity);
        self
    }

    /// Replaces endpoint coordinate style.
    #[must_use]
    pub const fn coordinates(mut self, coordinates: CoordinateStyle) -> Self {
        self.options = self.options.with_coordinates(coordinates);
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

impl<A: Alphabet> fmt::Display for MyersSearchReport<'_, '_, A> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(formatter, "Algorithm        : Myers bit-vector")?;
        writeln!(formatter, "Mode             : approximate suffix search")?;
        writeln!(formatter, "Pattern length   : {}", self.run.pattern.len())?;
        writeln!(formatter, "Text length      : {}", self.run.text.len())?;
        writeln!(
            formatter,
            "Maximum distance : {}",
            self.run.maximum_distance
        )?;
        write!(formatter, "Matches          : {}", self.run.matches.len())?;
        if self.options.verbosity() == Verbosity::Compact {
            return Ok(());
        }
        formatter.write_str("\n\n#  Text endpoint  Distance")?;
        for (index, matched) in self.run.matches.iter().copied().enumerate() {
            let endpoint = match self.options.coordinates() {
                CoordinateStyle::Human => matched.end().to_string(),
                CoordinateStyle::Rust => format!("0..{}", matched.end()),
            };
            write!(
                formatter,
                "\n{}  {:>13}  {}",
                index + 1,
                endpoint,
                matched.distance()
            )?;
        }
        Ok(())
    }
}
