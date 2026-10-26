use core::fmt;

use align_core::{AlignmentOperation, Alphabet, Score, Sequence};

use super::{
    ColorMode, CoordinateStyle, PairwiseStatistics, ReportError, ReportOptions, Verbosity,
    sanitized_label,
};
use crate::result::PairwiseAlignmentView;

/// Primary scientific objective represented by a pairwise run.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PairwiseMetric {
    /// A maximizing alignment score.
    Score(Score),
    /// A minimizing unit edit distance.
    Distance(usize),
    /// A minimizing nonnegative affine penalty.
    Penalty(Score),
}

impl PairwiseMetric {
    /// Returns the unambiguous field label.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Score(_) => "Score",
            Self::Distance(_) => "Distance",
            Self::Penalty(_) => "Penalty",
        }
    }
}

/// One algorithm-specific line retained for detailed output.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReportField {
    label: &'static str,
    value: String,
}

impl ReportField {
    /// Creates a detailed report field.
    #[must_use]
    pub fn new(label: &'static str, value: impl Into<String>) -> Self {
        Self {
            label,
            value: value.into(),
        }
    }

    /// Returns the stable field label.
    #[must_use]
    pub const fn label(&self) -> &'static str {
        self.label
    }

    /// Returns the formatted field value.
    #[must_use]
    pub fn value(&self) -> &str {
        &self.value
    }
}

/// Algorithm identity and objective metadata for a pairwise run.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PairwiseRunInfo {
    algorithm: &'static str,
    mode: &'static str,
    scoring: String,
    metric: PairwiseMetric,
    details: Vec<ReportField>,
}

impl PairwiseRunInfo {
    pub(crate) fn new(
        algorithm: &'static str,
        mode: &'static str,
        scoring: String,
        metric: PairwiseMetric,
    ) -> Self {
        Self {
            algorithm,
            mode,
            scoring,
            metric,
            details: Vec::new(),
        }
    }

    pub(crate) fn with_details(mut self, details: Vec<ReportField>) -> Self {
        self.details = details;
        self
    }

    /// Returns the public algorithm name.
    #[must_use]
    pub const fn algorithm(&self) -> &'static str {
        self.algorithm
    }

    /// Returns the selected alignment mode.
    #[must_use]
    pub const fn mode(&self) -> &'static str {
        self.mode
    }

    /// Returns the exact scoring or cost convention.
    #[must_use]
    pub fn scoring(&self) -> &str {
        &self.scoring
    }

    /// Returns the unambiguous primary objective.
    #[must_use]
    pub const fn metric(&self) -> PairwiseMetric {
        self.metric
    }

    /// Returns algorithm-specific evidence fields.
    #[must_use]
    pub fn details(&self) -> &[ReportField] {
        &self.details
    }
}

/// One native pairwise result plus the inputs and metadata required to explain it.
#[derive(Clone, Debug)]
pub struct PairwiseRun<'a, A: Alphabet, R> {
    first: &'a Sequence<A>,
    second: &'a Sequence<A>,
    result: R,
    info: PairwiseRunInfo,
}

impl<'a, A: Alphabet, R> PairwiseRun<'a, A, R> {
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

    /// Returns the untouched native algorithm result.
    #[must_use]
    pub const fn result(&self) -> &R {
        &self.result
    }

    /// Consumes the run and returns the untouched native result.
    #[must_use]
    pub fn into_result(self) -> R {
        self.result
    }

    /// Returns algorithm and objective metadata.
    #[must_use]
    pub const fn info(&self) -> &PairwiseRunInfo {
        &self.info
    }
}

impl<'a, A: Alphabet, R: PairwiseAlignmentView> PairwiseRun<'a, A, R> {
    /// Returns exact statistics derived from the canonical traceback.
    #[must_use]
    pub fn statistics(&self) -> PairwiseStatistics {
        PairwiseStatistics::from_alignment(self.result.alignment())
    }

    /// Creates a normal colorless report with 60 columns per block.
    #[must_use]
    pub fn report(&self) -> PairwiseReport<'_, 'a, A, R> {
        PairwiseReport {
            run: self,
            options: ReportOptions::new(),
            first_label: String::from("First"),
            second_label: String::from("Second"),
        }
    }

    /// Creates a report using explicit shared options.
    #[must_use]
    pub fn report_with(&self, options: ReportOptions) -> PairwiseReport<'_, 'a, A, R> {
        PairwiseReport {
            run: self,
            options,
            first_label: String::from("First"),
            second_label: String::from("Second"),
        }
    }
}

impl<A: Alphabet, R: PairwiseAlignmentView> PairwiseAlignmentView for PairwiseRun<'_, A, R> {
    fn alignment(&self) -> &align_core::PairwiseAlignment {
        self.result.alignment()
    }
}

/// Borrowed configurable human-readable pairwise report.
#[derive(Clone, Debug)]
pub struct PairwiseReport<'run, 'sequence, A: Alphabet, R> {
    run: &'run PairwiseRun<'sequence, A, R>,
    options: ReportOptions,
    first_label: String,
    second_label: String,
}

impl<A: Alphabet, R> PairwiseReport<'_, '_, A, R> {
    /// Replaces both human-facing sequence labels.
    #[must_use]
    pub fn labels(mut self, first: impl Into<String>, second: impl Into<String>) -> Self {
        self.first_label = sanitized_label(first);
        self.second_label = sanitized_label(second);
        self
    }

    /// Replaces the positive alignment block width.
    ///
    /// # Errors
    ///
    /// Returns [`ReportError::InvalidBlockWidth`] for zero.
    pub fn block_width(mut self, width: usize) -> Result<Self, ReportError> {
        self.options = self.options.with_block_width(width)?;
        Ok(self)
    }

    /// Replaces the verbosity.
    #[must_use]
    pub const fn verbosity(mut self, verbosity: Verbosity) -> Self {
        self.options = self.options.with_verbosity(verbosity);
        self
    }

    /// Replaces the coordinate representation.
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

impl<A: Alphabet, R: PairwiseAlignmentView> fmt::Display for PairwiseReport<'_, '_, A, R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let alignment = self.run.result.alignment();
        let statistics = PairwiseStatistics::from_alignment(alignment);
        let info = &self.run.info;
        let field_width = report_field_width(info, self.options.verbosity());

        write_field(formatter, field_width, "Algorithm", info.algorithm())?;
        write_field(formatter, field_width, "Mode", info.mode())?;
        write_field(formatter, field_width, "Scoring", info.scoring())?;
        writeln!(
            formatter,
            "{:<field_width$} : {}={} {}={}",
            "Sequences",
            self.first_label,
            self.run.first.len(),
            self.second_label,
            self.run.second.len()
        )?;
        write_ranges(
            formatter,
            field_width,
            alignment,
            &self.first_label,
            &self.second_label,
            self.options.coordinates(),
        )?;
        write_metric(formatter, field_width, info.metric())?;

        if self.options.verbosity() == Verbosity::Compact {
            return Ok(());
        }

        formatter.write_str("\n")?;

        writeln!(
            formatter,
            "{:<field_width$} : {}",
            "Columns",
            statistics.columns()
        )?;
        let identity = statistics.identity().map_or_else(
            || String::from("-"),
            |value| format!("{:.2}%", value * 100.0),
        );
        writeln!(
            formatter,
            "{:<field_width$} : {}/{} ({identity})",
            "Identity",
            statistics.matches(),
            statistics.columns()
        )?;
        writeln!(
            formatter,
            "{:<field_width$} : {}",
            "Matches",
            statistics.matches()
        )?;
        writeln!(
            formatter,
            "{:<field_width$} : {}",
            "Substitutions",
            statistics.substitutions()
        )?;
        writeln!(
            formatter,
            "{:<field_width$} : {} residues, {} runs",
            "Insertions",
            statistics.insertion_residues(),
            statistics.insertion_runs()
        )?;
        writeln!(
            formatter,
            "{:<field_width$} : {} residues, {} runs",
            "Deletions",
            statistics.deletion_residues(),
            statistics.deletion_runs()
        )?;
        write!(formatter, "{:<field_width$} : ", "CIGAR")?;
        write_cigar(formatter, alignment.operations())?;
        formatter.write_str("\n")?;

        if self.options.verbosity() == Verbosity::Detailed {
            for field in info.details() {
                write_field(formatter, field_width, field.label(), field.value())?;
            }
        }

        formatter.write_str("\n")?;
        write_alignment_blocks(
            formatter,
            self.run.first,
            self.run.second,
            alignment,
            &self.first_label,
            &self.second_label,
            self.options,
        )?;
        formatter.write_str("\nLegend: | match, . substitution, blank gap")
    }
}

fn report_field_width(info: &PairwiseRunInfo, verbosity: Verbosity) -> usize {
    let mut width = 13_usize;
    if verbosity == Verbosity::Detailed {
        for field in info.details() {
            width = width.max(field.label().chars().count());
        }
    }
    width
}

fn write_field(
    formatter: &mut fmt::Formatter<'_>,
    width: usize,
    label: &str,
    value: &str,
) -> fmt::Result {
    writeln!(formatter, "{label:<width$} : {value}")
}

fn write_metric(
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

fn write_ranges(
    formatter: &mut fmt::Formatter<'_>,
    width: usize,
    alignment: &align_core::PairwiseAlignment,
    first_label: &str,
    second_label: &str,
    coordinates: CoordinateStyle,
) -> fmt::Result {
    let first = alignment.first_range();
    let second = alignment.second_range();
    write!(formatter, "{:<width$} : {first_label} ", "Range")?;
    write_range(formatter, first, coordinates)?;
    write!(formatter, ", {second_label} ")?;
    write_range(formatter, second, coordinates)?;
    formatter.write_str("\n")
}

fn write_range(
    formatter: &mut fmt::Formatter<'_>,
    range: core::ops::Range<usize>,
    coordinates: CoordinateStyle,
) -> fmt::Result {
    if range.is_empty() {
        return formatter.write_str("-");
    }
    match coordinates {
        CoordinateStyle::Human => write!(formatter, "{}..{}", range.start + 1, range.end),
        CoordinateStyle::Rust => write!(formatter, "{}..{}", range.start, range.end),
    }
}

fn write_cigar(
    formatter: &mut fmt::Formatter<'_>,
    operations: &[AlignmentOperation],
) -> fmt::Result {
    let Some(mut selected) = operations.first().copied() else {
        return formatter.write_str("-");
    };
    let mut count = 0_usize;
    for operation in operations.iter().copied() {
        if operation == selected {
            count += 1;
        } else {
            write!(formatter, "{count}{}", operation_code(selected))?;
            selected = operation;
            count = 1;
        }
    }
    write!(formatter, "{count}{}", operation_code(selected))
}

const fn operation_code(operation: AlignmentOperation) -> char {
    match operation {
        AlignmentOperation::Match => '=',
        AlignmentOperation::Substitution => 'X',
        AlignmentOperation::Deletion => 'D',
        AlignmentOperation::Insertion => 'I',
    }
}

#[allow(clippy::too_many_arguments)]
fn write_alignment_blocks<A: Alphabet>(
    formatter: &mut fmt::Formatter<'_>,
    first: &Sequence<A>,
    second: &Sequence<A>,
    alignment: &align_core::PairwiseAlignment,
    first_label: &str,
    second_label: &str,
    options: ReportOptions,
) -> fmt::Result {
    if alignment.operations().is_empty() {
        return formatter.write_str("(empty alignment)");
    }
    let label_width = first_label
        .chars()
        .count()
        .max(second_label.chars().count());
    let coordinate_width = decimal_width(first.len().max(second.len()).max(1));
    let mut first_position = alignment.first_range().start;
    let mut second_position = alignment.second_range().start;

    for (block_index, operations) in alignment
        .operations()
        .chunks(options.block_width())
        .enumerate()
    {
        if block_index != 0 {
            formatter.write_str("\n\n")?;
        }
        let first_before = first_position;
        let second_before = second_position;
        let mut first_row = String::with_capacity(operations.len());
        let mut second_row = String::with_capacity(operations.len());
        let mut middle = String::with_capacity(operations.len());
        for operation in operations.iter().copied() {
            match operation {
                AlignmentOperation::Match | AlignmentOperation::Substitution => {
                    let first_symbol = first
                        .get(first_position)
                        .expect("validated alignment first range is in bounds");
                    let second_symbol = second
                        .get(second_position)
                        .expect("validated alignment second range is in bounds");
                    first_row.push(char::from(first_symbol.as_ascii()));
                    second_row.push(char::from(second_symbol.as_ascii()));
                    middle.push(if operation == AlignmentOperation::Match {
                        '|'
                    } else {
                        '.'
                    });
                    first_position += 1;
                    second_position += 1;
                }
                AlignmentOperation::Deletion => {
                    let symbol = first
                        .get(first_position)
                        .expect("validated alignment first range is in bounds");
                    first_row.push(char::from(symbol.as_ascii()));
                    second_row.push('-');
                    middle.push(' ');
                    first_position += 1;
                }
                AlignmentOperation::Insertion => {
                    let symbol = second
                        .get(second_position)
                        .expect("validated alignment second range is in bounds");
                    first_row.push('-');
                    second_row.push(char::from(symbol.as_ascii()));
                    middle.push(' ');
                    second_position += 1;
                }
            }
        }
        write_block_row(
            formatter,
            first_label,
            label_width,
            coordinate_width,
            first_before,
            first_position,
            &first_row,
            options.coordinates(),
        )?;
        writeln!(
            formatter,
            "{:<label_width$}  {:>coordinate_width$}  {middle}",
            "", ""
        )?;
        write_block_row(
            formatter,
            second_label,
            label_width,
            coordinate_width,
            second_before,
            second_position,
            &second_row,
            options.coordinates(),
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_block_row(
    formatter: &mut fmt::Formatter<'_>,
    label: &str,
    label_width: usize,
    coordinate_width: usize,
    before: usize,
    after: usize,
    row: &str,
    coordinates: CoordinateStyle,
) -> fmt::Result {
    if before == after {
        return writeln!(
            formatter,
            "{label:<label_width$}  {:>coordinate_width$}  {row}  {:>coordinate_width$}",
            "-", "-"
        );
    }
    let (start, end) = match coordinates {
        CoordinateStyle::Human => (before + 1, after),
        CoordinateStyle::Rust => (before, after),
    };
    writeln!(
        formatter,
        "{label:<label_width$}  {start:>coordinate_width$}  {row}  {end:>coordinate_width$}"
    )
}

fn decimal_width(value: usize) -> usize {
    value.to_string().len()
}
