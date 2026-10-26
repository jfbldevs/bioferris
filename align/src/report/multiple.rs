use core::fmt;

use align_core::Alphabet;
use align_msa::MultipleAlignment;

use super::{
    ColorMode, CoordinateStyle, MultipleStatistics, ReportError, ReportField, ReportOptions,
    Verbosity,
};
use crate::result::MultipleAlignmentView;

/// Algorithm identity and retained evidence labels for an MSA run.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MultipleRunInfo {
    algorithm: &'static str,
    profile: &'static str,
    mode: Option<&'static str>,
    summary: Vec<ReportField>,
    details: Vec<ReportField>,
}

impl MultipleRunInfo {
    pub(crate) fn new(
        algorithm: &'static str,
        profile: &'static str,
        mode: Option<&'static str>,
    ) -> Self {
        Self {
            algorithm,
            profile,
            mode,
            summary: Vec::new(),
            details: Vec::new(),
        }
    }

    pub(crate) fn with_summary(mut self, summary: Vec<ReportField>) -> Self {
        self.summary = summary;
        self
    }

    pub(crate) fn with_details(mut self, details: Vec<ReportField>) -> Self {
        self.details = details;
        self
    }

    /// Returns the algorithm-family name.
    #[must_use]
    pub const fn algorithm(&self) -> &'static str {
        self.algorithm
    }

    /// Returns the exact immutable scientific profile identifier.
    #[must_use]
    pub const fn profile(&self) -> &'static str {
        self.profile
    }

    /// Returns the selected mode when the family defines one.
    #[must_use]
    pub const fn mode(&self) -> Option<&'static str> {
        self.mode
    }

    /// Returns algorithm-specific evidence fields.
    #[must_use]
    pub fn details(&self) -> &[ReportField] {
        &self.details
    }

    /// Returns family-specific fields included in normal output.
    #[must_use]
    pub fn summary(&self) -> &[ReportField] {
        &self.summary
    }
}

/// One native MSA result plus the metadata required for presentation.
#[derive(Clone, Debug)]
pub struct MultipleRun<A: Alphabet, R> {
    result: R,
    info: MultipleRunInfo,
    marker: core::marker::PhantomData<A>,
}

impl<A: Alphabet, R> MultipleRun<A, R> {
    pub(crate) const fn new(result: R, info: MultipleRunInfo) -> Self {
        Self {
            result,
            info,
            marker: core::marker::PhantomData,
        }
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

    /// Returns algorithm and profile metadata.
    #[must_use]
    pub const fn info(&self) -> &MultipleRunInfo {
        &self.info
    }
}

impl<A, R> MultipleRun<A, R>
where
    A: Alphabet,
    R: MultipleAlignmentView<Alphabet = A>,
{
    /// Returns the canonical rectangular MSA.
    #[must_use]
    pub fn alignment(&self) -> &MultipleAlignment<A> {
        self.result.alignment()
    }

    /// Returns exact structural statistics.
    #[must_use]
    pub fn statistics(&self) -> MultipleStatistics {
        MultipleStatistics::from_alignment(self.alignment())
    }

    /// Creates a normal colorless report with 60 columns per block.
    #[must_use]
    pub const fn report(&self) -> MultipleReport<'_, A, R> {
        MultipleReport {
            run: self,
            options: ReportOptions::new(),
        }
    }

    /// Creates a report from shared options.
    #[must_use]
    pub const fn report_with(&self, options: ReportOptions) -> MultipleReport<'_, A, R> {
        MultipleReport { run: self, options }
    }
}

impl<A, R> MultipleAlignmentView for MultipleRun<A, R>
where
    A: Alphabet,
    R: MultipleAlignmentView<Alphabet = A>,
{
    type Alphabet = A;

    fn alignment(&self) -> &MultipleAlignment<A> {
        self.result.alignment()
    }
}

/// Borrowed configurable human-readable MSA report.
#[derive(Clone, Copy, Debug)]
pub struct MultipleReport<'run, A: Alphabet, R> {
    run: &'run MultipleRun<A, R>,
    options: ReportOptions,
}

impl<A: Alphabet, R> MultipleReport<'_, A, R> {
    /// Replaces the positive alignment block width.
    ///
    /// # Errors
    ///
    /// Returns [`ReportError::InvalidBlockWidth`] for zero.
    pub fn block_width(mut self, width: usize) -> Result<Self, ReportError> {
        self.options = self.options.with_block_width(width)?;
        Ok(self)
    }

    /// Replaces report verbosity.
    #[must_use]
    pub const fn verbosity(mut self, verbosity: Verbosity) -> Self {
        self.options = self.options.with_verbosity(verbosity);
        self
    }

    /// Replaces row-coordinate representation.
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

    /// Returns shared report options.
    #[must_use]
    pub const fn options(&self) -> ReportOptions {
        self.options
    }
}

impl<A, R> fmt::Display for MultipleReport<'_, A, R>
where
    A: Alphabet,
    R: MultipleAlignmentView<Alphabet = A>,
{
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let info = self.run.info();
        let alignment = self.run.alignment();
        let statistics = MultipleStatistics::from_alignment(alignment);
        let field_width = multiple_field_width(info, self.options.verbosity());

        writeln!(
            formatter,
            "{:<field_width$} : {}",
            "Algorithm",
            info.algorithm()
        )?;
        writeln!(
            formatter,
            "{:<field_width$} : {}",
            "Profile",
            info.profile()
        )?;
        if let Some(mode) = info.mode() {
            writeln!(formatter, "{:<field_width$} : {mode}", "Mode")?;
        }
        writeln!(
            formatter,
            "{:<field_width$} : {}",
            "Sequences",
            statistics.sequences()
        )?;
        write!(
            formatter,
            "{:<field_width$} : {}",
            "Columns",
            statistics.columns()
        )?;
        if self.options.verbosity() == Verbosity::Compact {
            return Ok(());
        }
        formatter.write_str("\n")?;
        writeln!(
            formatter,
            "{:<field_width$} : {}",
            "Conserved",
            statistics.conserved_columns()
        )?;
        writeln!(
            formatter,
            "{:<field_width$} : {}",
            "Gap characters",
            statistics.gap_characters()
        )?;
        for field in info.summary() {
            writeln!(
                formatter,
                "{:<field_width$} : {}",
                field.label(),
                field.value()
            )?;
        }
        if self.options.verbosity() == Verbosity::Detailed {
            for field in info.details() {
                writeln!(
                    formatter,
                    "{:<field_width$} : {}",
                    field.label(),
                    field.value()
                )?;
            }
        }
        formatter.write_str("\n")?;
        write_msa_blocks(formatter, alignment, self.options)?;
        formatter.write_str("\nLegend: * exact conservation in every sequence")
    }
}

fn multiple_field_width(info: &MultipleRunInfo, verbosity: Verbosity) -> usize {
    let mut width = 14_usize;
    if verbosity == Verbosity::Detailed {
        for field in info.details() {
            width = width.max(field.label().chars().count());
        }
    }
    for field in info.summary() {
        width = width.max(field.label().chars().count());
    }
    width
}

fn write_msa_blocks<A: Alphabet>(
    formatter: &mut fmt::Formatter<'_>,
    alignment: &MultipleAlignment<A>,
    options: ReportOptions,
) -> fmt::Result {
    if alignment.width() == 0 {
        return formatter.write_str("(empty alignment)");
    }
    let label_width = alignment
        .rows()
        .iter()
        .map(|row| row.identifier().chars().count())
        .max()
        .unwrap_or(0);
    let maximum_length = alignment
        .rows()
        .iter()
        .map(|row| row.symbols().iter().flatten().count())
        .max()
        .unwrap_or(0);
    let coordinate_width = maximum_length.max(1).to_string().len();
    let mut positions = vec![0_usize; alignment.sequence_count()];

    for (block_index, start) in (0..alignment.width())
        .step_by(options.block_width())
        .enumerate()
    {
        if block_index != 0 {
            formatter.write_str("\n\n")?;
        }
        let end = start
            .saturating_add(options.block_width())
            .min(alignment.width());
        for (row_index, row) in alignment.rows().iter().enumerate() {
            let before = positions[row_index];
            let mut rendered = String::with_capacity(end - start);
            for cell in &row.symbols()[start..end] {
                if let Some(symbol) = cell {
                    positions[row_index] += 1;
                    rendered.push(char::from(symbol.as_ascii()));
                } else {
                    rendered.push('-');
                }
            }
            write_msa_row(
                formatter,
                row.identifier(),
                label_width,
                coordinate_width,
                before,
                positions[row_index],
                &rendered,
                options.coordinates(),
            )?;
        }
        let mut conservation = String::with_capacity(end - start);
        for column in start..end {
            conservation.push(if exactly_conserved(alignment, column) {
                '*'
            } else {
                ' '
            });
        }
        write!(
            formatter,
            "{:<label_width$}  {:>coordinate_width$}  {conservation}",
            "", ""
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_msa_row(
    formatter: &mut fmt::Formatter<'_>,
    identifier: &str,
    label_width: usize,
    coordinate_width: usize,
    before: usize,
    after: usize,
    rendered: &str,
    coordinates: CoordinateStyle,
) -> fmt::Result {
    if before == after {
        return writeln!(
            formatter,
            "{identifier:<label_width$}  {:>coordinate_width$}  {rendered}  {:>coordinate_width$}",
            "-", "-"
        );
    }
    let (first, last) = match coordinates {
        CoordinateStyle::Human => (before + 1, after),
        CoordinateStyle::Rust => (before, after),
    };
    writeln!(
        formatter,
        "{identifier:<label_width$}  {first:>coordinate_width$}  {rendered}  {last:>coordinate_width$}"
    )
}

fn exactly_conserved<A: Alphabet>(alignment: &MultipleAlignment<A>, column: usize) -> bool {
    let mut selected = None;
    for row in alignment.rows() {
        let Some(symbol) = row.symbols()[column] else {
            return false;
        };
        if let Some(expected) = selected {
            if symbol != expected {
                return false;
            }
        } else {
            selected = Some(symbol);
        }
    }
    selected.is_some()
}
