//! Bounded IUPAC nucleotide searches and explicit count/probability/score matrices.
//!
//! DNA and RNA motifs share A/C/G/T bit sets (U uses the T bit), but retain
//! molecule-specific sequence validation. All reported ranges are zero-based,
//! half-open positions in the original input. Both-strand results are ordered
//! by original start and then forward before reverse, including palindromes.
//!
//! Matrix columns use A/C/G/T order for DNA and A/C/G/U for RNA. Matrix scans
//! require concrete input and accumulate cells in `f64` before returning an
//! `f32` score. Log odds use `ln(probability / background) / ln(2)`, and scan
//! thresholds round once to `f32` before inclusive comparisons. Negative
//! infinity is a meaningful score for zero probabilities.
//! Column and background totals use Neumaier compensated summation, matching
//! the declared `CPython` 3.12 numerical observation profile.
//! Search budgets bound the complete worst-case number of windows and symbol
//! comparisons before scanning, rather than depending on early mismatches.

use std::{error::Error, fmt, marker::PhantomData, ops::Range};

use align_core::Alphabet;

use crate::{Dna, Rna, Sequence, SequenceError, Strand, StrandSelection};

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::Dna {}
    impl Sealed for super::Rna {}
}

/// A supported nucleotide alphabet, sealed to the existing DNA and RNA types.
///
/// Custom alphabets cannot implement this trait because neither their IUPAC
/// masks nor their biological complements are defined by this module.
pub trait MotifAlphabet: Alphabet + sealed::Sealed {}

impl MotifAlphabet for Dna {}
impl MotifAlphabet for Rna {}

/// Interpretation of ambiguity in the searched input, separate from the pattern.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum MotifInputAmbiguity {
    /// Every possible input base belongs to the pattern's set (default).
    #[default]
    Definite,
    /// At least one input base belongs to the pattern's set.
    Possible,
    /// Reject the first ambiguous symbol anywhere in the original input.
    Reject,
}

/// Strand selection, ambiguity semantics and global IUPAC search budgets.
///
/// All overlaps are examined. Zero match budget permits only an empty result;
/// excess matches fail without returning a partial result. Window/comparison
/// budgets conservatively include all selected windows and full pattern widths,
/// even when an early mismatch would reduce the actual work.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MotifSearchOptions {
    /// Selected strands; both by default.
    pub strands: StrandSelection,
    /// Source ambiguity semantics; definite by default.
    pub ambiguity: MotifInputAmbiguity,
    /// Maximum total matches across strands; default 100,000.
    pub max_matches: usize,
    /// Maximum total windows across strands; default 10,000,000.
    pub max_windows: usize,
    /// Maximum total window-symbol comparisons; default 100,000,000.
    pub max_comparisons: usize,
}

impl Default for MotifSearchOptions {
    fn default() -> Self {
        Self {
            strands: StrandSelection::Both,
            ambiguity: MotifInputAmbiguity::Definite,
            max_matches: 100_000,
            max_windows: 10_000_000,
            max_comparisons: 100_000_000,
        }
    }
}

/// One overlapping IUPAC match, retaining its original-input span and strand.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MotifMatch {
    range: Range<usize>,
    strand: Strand,
}

impl MotifMatch {
    /// Returns the zero-based half-open range in the original input.
    #[must_use]
    pub fn range(&self) -> Range<usize> {
        self.range.clone()
    }

    /// Returns the strand matching the pattern's 5′ to 3′ orientation.
    #[must_use]
    pub const fn strand(&self) -> Strand {
        self.strand
    }
}

/// An immutable, validated, nonempty IUPAC nucleotide pattern.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct IupacMotif<A: MotifAlphabet> {
    pattern: String,
    masks: Vec<u8>,
    reverse_masks: Vec<u8>,
    marker: PhantomData<A>,
}

impl<A: MotifAlphabet> IupacMotif<A> {
    /// Validates and compiles a pattern without expanding its possible sequences.
    ///
    /// Lowercase is normalized by the existing sequence validator.
    ///
    /// # Errors
    /// Returns [`MotifError::InvalidPattern`] for an invalid molecule-specific
    /// symbol, or [`MotifError::EmptyPattern`] for zero width.
    pub fn new(pattern: impl AsRef<[u8]>) -> Result<Self, MotifError> {
        let sequence = Sequence::<A>::from_ascii(pattern).map_err(MotifError::InvalidPattern)?;
        if sequence.is_empty() {
            return Err(MotifError::EmptyPattern);
        }
        let masks: Vec<_> = sequence
            .iter()
            .map(|base| base_mask(base.as_ascii()))
            .collect();
        let reverse_masks = masks.iter().rev().copied().map(complement_mask).collect();
        Ok(Self {
            pattern: sequence.to_ascii_string(),
            masks,
            reverse_masks,
            marker: PhantomData,
        })
    }

    /// Returns the normalized uppercase pattern in its original orientation.
    #[must_use]
    pub fn pattern(&self) -> &str {
        &self.pattern
    }

    /// Returns the number of nucleotide positions in the nonempty pattern.
    #[must_use]
    pub fn width(&self) -> usize {
        self.masks.len()
    }

    /// Finds both-strand definite matches using the default search budgets.
    ///
    /// # Errors
    /// Returns [`MotifError::LimitExceeded`] when a search budget is exceeded.
    pub fn find_in(&self, input: &Sequence<A>) -> Result<Vec<MotifMatch>, MotifError> {
        self.find_with(input, MotifSearchOptions::default())
    }

    /// Finds overlapping matches with explicit source ambiguity and strand rules.
    ///
    /// Empty or short inputs have no windows. Reject-mode validation examines
    /// the entire source first, including symbols outside complete windows.
    /// Search leaves the input unchanged and does not allocate its reverse copy.
    ///
    /// # Errors
    /// Returns [`MotifError::AmbiguousInput`] under the reject policy, a budget
    /// error, or an arithmetic error if the predicted work cannot fit `usize`.
    pub fn find_with(
        &self,
        input: &Sequence<A>,
        options: MotifSearchOptions,
    ) -> Result<Vec<MotifMatch>, MotifError> {
        if options.ambiguity == MotifInputAmbiguity::Reject {
            validate_concrete(input)?;
        }
        let windows = search_windows(
            input.len(),
            self.width(),
            options.strands,
            options.max_windows,
            options.max_comparisons,
        )?;
        let mut matches = Vec::new();
        for start in 0..windows {
            for strand in selected_strands(options.strands) {
                let pattern = match strand {
                    Strand::Forward => &self.masks,
                    Strand::Reverse => &self.reverse_masks,
                };
                let found = pattern.iter().enumerate().all(|(offset, pattern_mask)| {
                    let input_mask = base_mask(input.as_slice()[start + offset].as_ascii());
                    match options.ambiguity {
                        MotifInputAmbiguity::Definite | MotifInputAmbiguity::Reject => {
                            input_mask & pattern_mask == input_mask
                        }
                        MotifInputAmbiguity::Possible => input_mask & pattern_mask != 0,
                    }
                });
                if found {
                    check_next_match(matches.len(), options.max_matches)?;
                    matches.push(MotifMatch {
                        range: start..start + self.width(),
                        strand: *strand,
                    });
                }
            }
        }
        Ok(matches)
    }
}

/// A molecule-specific IUPAC DNA motif.
pub type DnaMotif = IupacMotif<Dna>;
/// A molecule-specific IUPAC RNA motif.
pub type RnaMotif = IupacMotif<Rna>;

/// Validates a nonempty IUPAC DNA motif.
///
/// # Errors
/// Returns [`MotifError`] for an empty or invalid pattern.
pub fn dna_motif(pattern: impl AsRef<[u8]>) -> Result<DnaMotif, MotifError> {
    DnaMotif::new(pattern)
}

/// Validates a nonempty IUPAC RNA motif.
///
/// # Errors
/// Returns [`MotifError`] for an empty or invalid pattern.
pub fn rna_motif(pattern: impl AsRef<[u8]>) -> Result<RnaMotif, MotifError> {
    RnaMotif::new(pattern)
}

/// Nonnegative finite count columns in A/C/G/T (or A/C/G/U) order.
///
/// Counts may be fractional, and different columns may have different totals.
/// Zero-total columns can be stored but require positive pseudocount totals
/// before normalization. At least one column is required.
#[derive(Clone, Debug, PartialEq)]
pub struct MotifCounts<A: MotifAlphabet> {
    columns: Vec<[f64; 4]>,
    marker: PhantomData<A>,
}

impl<A: MotifAlphabet> MotifCounts<A> {
    /// Builds count columns from nonempty, equally wide, concrete instances.
    ///
    /// Instances are counted independently and are not aligned or gap-filled.
    ///
    /// # Errors
    /// Returns [`MotifError`] for no/empty instances, unequal widths, ambiguous
    /// bases or an arithmetic overflow. Ambiguity offsets identify the position
    /// within the rejected instance; instances are examined in supplied order.
    pub fn from_instances(instances: &[Sequence<A>]) -> Result<Self, MotifError> {
        let first = instances.first().ok_or(MotifError::InvalidMatrix {
            column: None,
            reason: "at least one motif instance is required",
        })?;
        if first.is_empty() {
            return Err(MotifError::InvalidMatrix {
                column: None,
                reason: "motif instances must be nonempty",
            });
        }
        let mut columns: Vec<[f64; 4]> = vec![[0.0; 4]; first.len()];
        for (instance, sequence) in instances.iter().enumerate() {
            if sequence.len() != first.len() {
                return Err(MotifError::DimensionMismatch {
                    expected: first.len(),
                    actual: sequence.len(),
                    instance: Some(instance),
                });
            }
            validate_concrete(sequence)?;
            for (column, symbol) in sequence.iter().enumerate() {
                if let Some(base) = concrete_index(symbol.as_ascii()) {
                    let incremented = columns[column][base] + 1.0;
                    if !incremented.is_finite()
                        || incremented.to_bits() == columns[column][base].to_bits()
                    {
                        return Err(MotifError::InvalidMatrix {
                            column: Some(column),
                            reason: "instance count cannot be represented exactly as f64",
                        });
                    }
                    columns[column][base] = incremented;
                }
            }
        }
        Self::from_columns(columns)
    }

    /// Validates explicit columns, allowing weighted counts and zero columns.
    ///
    /// # Errors
    /// Returns [`MotifError::InvalidMatrix`] for zero width, negative counts,
    /// nonfinite cells or a column whose finite total overflows.
    pub fn from_columns(columns: Vec<[f64; 4]>) -> Result<Self, MotifError> {
        validate_width(columns.len())?;
        for (column, values) in columns.iter().enumerate() {
            validate_nonnegative(values, Some(column))?;
            finite_sum(values, Some(column))?;
        }
        Ok(Self {
            columns,
            marker: PhantomData,
        })
    }

    /// Returns immutable columns in A/C/G/T or A/C/G/U order.
    #[must_use]
    pub fn columns(&self) -> &[[f64; 4]] {
        &self.columns
    }

    /// Returns the nonzero number of positions.
    #[must_use]
    pub fn width(&self) -> usize {
        self.columns.len()
    }

    /// Adds per-base pseudocounts and normalizes each column independently.
    ///
    /// No pseudocount is inferred from the background or instance count.
    ///
    /// # Errors
    /// Returns [`MotifError::InvalidMatrix`] for invalid pseudocounts, a
    /// zero-total column or nonfinite intermediate arithmetic.
    pub fn normalize(&self, pseudocounts: [f64; 4]) -> Result<MotifWeights<A>, MotifError> {
        validate_nonnegative(&pseudocounts, None)?;
        let mut columns = Vec::with_capacity(self.width());
        for (column, counts) in self.columns.iter().enumerate() {
            let mut weights = [0.0; 4];
            for base in 0..4 {
                weights[base] = counts[base] + pseudocounts[base];
            }
            let total = finite_sum(&weights, Some(column))?;
            if total == 0.0 {
                return Err(MotifError::InvalidMatrix {
                    column: Some(column),
                    reason: "count plus pseudocount column has zero total",
                });
            }
            for probability in &mut weights {
                *probability /= total;
            }
            columns.push(weights);
        }
        Ok(MotifWeights {
            columns,
            marker: PhantomData,
        })
    }

    /// Normalizes with zero pseudocounts, preserving zero-probability cells.
    ///
    /// # Errors
    /// Returns [`MotifError::InvalidMatrix`] for a zero-total column.
    pub fn normalize_default(&self) -> Result<MotifWeights<A>, MotifError> {
        self.normalize([0.0; 4])
    }
}

/// Independently normalized probability columns in nucleotide order.
#[derive(Clone, Debug, PartialEq)]
pub struct MotifWeights<A: MotifAlphabet> {
    columns: Vec<[f64; 4]>,
    marker: PhantomData<A>,
}

impl<A: MotifAlphabet> MotifWeights<A> {
    /// Returns immutable A/C/G/T or A/C/G/U probability columns.
    #[must_use]
    pub fn columns(&self) -> &[[f64; 4]] {
        &self.columns
    }

    /// Returns the nonzero number of positions.
    #[must_use]
    pub fn width(&self) -> usize {
        self.columns.len()
    }

    /// Converts probabilities to base-2 log odds against normalized background.
    ///
    /// Background weights must be finite and strictly positive; their sum need
    /// not already equal one. Zero probabilities produce negative infinity.
    ///
    /// # Errors
    /// Returns [`MotifError::InvalidMatrix`] for an invalid background,
    /// background underflow, ratio overflow or another nonfinite intermediate.
    pub fn log_odds(&self, background: [f64; 4]) -> Result<MotifScores<A>, MotifError> {
        let background = normalize_background(background)?;
        let mut columns = Vec::with_capacity(self.width());
        for (column, weights) in self.columns.iter().enumerate() {
            let mut scores = [0.0; 4];
            for base in 0..4 {
                scores[base] = if weights[base] == 0.0 {
                    f64::NEG_INFINITY
                } else {
                    let ratio = weights[base] / background[base];
                    if !ratio.is_finite() {
                        return Err(MotifError::InvalidMatrix {
                            column: Some(column),
                            reason: "probability to background ratio overflows f64",
                        });
                    }
                    // The mathematical log2 and this quotient can differ in
                    // their final bit. Keep the declared reference evaluation.
                    ratio.ln() / std::f64::consts::LN_2
                };
            }
            columns.push(scores);
        }
        // The background is already normalized; avoid an additional division.
        Ok(MotifScores {
            columns,
            background,
            marker: PhantomData,
        })
    }

    /// Converts to log odds against the uniform nucleotide background.
    ///
    /// # Errors
    /// Returns [`MotifError`] if matrix arithmetic cannot be represented.
    pub fn log_odds_uniform(&self) -> Result<MotifScores<A>, MotifError> {
        self.log_odds([0.25; 4])
    }
}

/// Explicit strand selection, inclusive score threshold and global scan budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MatrixSearchOptions {
    /// Selected strands; both by default.
    pub strands: StrandSelection,
    /// Inclusive threshold, rounded to float32 before comparison; default zero.
    /// Must be finite and within the finite float32 range. Two different f64
    /// thresholds that round to the same f32 therefore select the same hits.
    pub threshold: f64,
    /// Maximum retained hits across strands; default 100,000.
    pub max_matches: usize,
    /// Maximum evaluated windows across strands; default 10,000,000.
    pub max_windows: usize,
    /// Maximum window-symbol lookups across strands; default 100,000,000.
    pub max_comparisons: usize,
}

impl Default for MatrixSearchOptions {
    fn default() -> Self {
        Self {
            strands: StrandSelection::Both,
            threshold: 0.0,
            max_matches: 100_000,
            max_windows: 10_000_000,
            max_comparisons: 100_000_000,
        }
    }
}

/// A scoring-matrix hit in original-input coordinates with its float32 score.
#[derive(Clone, Debug, PartialEq)]
pub struct ScoredMotifMatch {
    range: Range<usize>,
    strand: Strand,
    score: f32,
}

impl ScoredMotifMatch {
    /// Returns the zero-based half-open original-input span.
    #[must_use]
    pub fn range(&self) -> Range<usize> {
        self.range.clone()
    }

    /// Returns the matching strand.
    #[must_use]
    pub const fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns the accumulated score after float32 rounding.
    #[must_use]
    pub const fn score(&self) -> f32 {
        self.score
    }
}

/// Base-2 log-odds score columns and the normalized background that defines them.
///
/// Cells are stored as `f64`; negative infinity is valid. Scans reject ambiguous
/// source symbols globally and return float32 scores after double accumulation.
#[derive(Clone, Debug, PartialEq)]
pub struct MotifScores<A: MotifAlphabet> {
    columns: Vec<[f64; 4]>,
    background: [f64; 4],
    marker: PhantomData<A>,
}

impl<A: MotifAlphabet> MotifScores<A> {
    /// Validates explicit score columns and records normalized background.
    ///
    /// Cells are already scores and are not transformed using the background.
    /// This constructor supports external scoring models and numerical probes.
    ///
    /// # Errors
    /// Returns [`MotifError::InvalidMatrix`] for zero width, NaN, positive
    /// infinity or an invalid background. Negative infinity is accepted.
    pub fn from_columns(columns: Vec<[f64; 4]>, background: [f64; 4]) -> Result<Self, MotifError> {
        validate_width(columns.len())?;
        for (column, values) in columns.iter().enumerate() {
            if values
                .iter()
                .any(|value| value.is_nan() || *value == f64::INFINITY)
            {
                return Err(MotifError::InvalidMatrix {
                    column: Some(column),
                    reason: "scores must be finite or negative infinity",
                });
            }
        }
        Ok(Self {
            columns,
            background: normalize_background(background)?,
            marker: PhantomData,
        })
    }

    /// Returns immutable score columns in A/C/G/T or A/C/G/U order.
    #[must_use]
    pub fn columns(&self) -> &[[f64; 4]] {
        &self.columns
    }

    /// Returns the normalized A/C/G/T or A/C/G/U background probabilities.
    #[must_use]
    pub const fn background(&self) -> [f64; 4] {
        self.background
    }

    /// Returns the nonzero number of positions.
    #[must_use]
    pub fn width(&self) -> usize {
        self.columns.len()
    }

    /// Scores exactly one concrete forward-oriented window.
    ///
    /// Cell addition proceeds left to right in `f64`, then rounds once to `f32`.
    /// A zero-probability cell produces negative infinity.
    ///
    /// # Errors
    /// Returns [`MotifError`] for a width mismatch, ambiguous symbol, arithmetic
    /// overflow or a finite sum outside the finite float32 range.
    pub fn score_window(&self, sequence: &Sequence<A>) -> Result<f32, MotifError> {
        if sequence.len() != self.width() {
            return Err(MotifError::DimensionMismatch {
                expected: self.width(),
                actual: sequence.len(),
                instance: None,
            });
        }
        validate_concrete(sequence)?;
        self.score_at(sequence, 0, Strand::Forward)
    }

    /// Scores all complete overlapping forward windows under a window budget.
    ///
    /// Empty or short input returns an empty vector after concrete validation.
    /// This call bounds output windows; use [`Self::find_with`] when a separate
    /// comparison budget and bounded retained hit count are required.
    ///
    /// # Errors
    /// Returns [`MotifError`] for ambiguous input, excess windows or score
    /// arithmetic errors. Errors discard the locally accumulated score vector.
    pub fn calculate(
        &self,
        sequence: &Sequence<A>,
        max_windows: usize,
    ) -> Result<Vec<f32>, MotifError> {
        validate_concrete(sequence)?;
        let windows = window_count(sequence.len(), self.width());
        check_limit("windows", windows, max_windows)?;
        let mut scores = Vec::new();
        for start in 0..windows {
            scores.push(self.score_at(sequence, start, Strand::Forward)?);
        }
        Ok(scores)
    }

    /// Finds all overlapping windows meeting an inclusive score threshold.
    ///
    /// Both-strand hits use original positions and are ordered by start then
    /// forward before reverse. Reverse scans complement the built scores,
    /// preserving asymmetric-background semantics. Palindromes retain both
    /// strand labels. No successful result is silently truncated by a budget.
    ///
    /// # Errors
    /// Returns [`MotifError`] for invalid options, ambiguous source, a budget
    /// failure or nonrepresentable arithmetic. Validation covers the whole
    /// original input, including symbols outside complete scoring windows.
    pub fn find_with(
        &self,
        sequence: &Sequence<A>,
        options: MatrixSearchOptions,
    ) -> Result<Vec<ScoredMotifMatch>, MotifError> {
        if !options.threshold.is_finite() || options.threshold.abs() > f64::from(f32::MAX) {
            return Err(MotifError::InvalidOptions {
                reason: "threshold must be finite and within the finite float32 range",
            });
        }
        // This one rounding is part of the scan contract, including thresholds
        // just above a score in f64 that round back to that same f32 score.
        #[allow(clippy::cast_possible_truncation)]
        let threshold = options.threshold as f32;
        validate_concrete(sequence)?;
        let windows = search_windows(
            sequence.len(),
            self.width(),
            options.strands,
            options.max_windows,
            options.max_comparisons,
        )?;
        let mut matches = Vec::new();
        for start in 0..windows {
            for strand in selected_strands(options.strands) {
                let score = self.score_at(sequence, start, *strand)?;
                if score >= threshold {
                    check_next_match(matches.len(), options.max_matches)?;
                    matches.push(ScoredMotifMatch {
                        range: start..start + self.width(),
                        strand: *strand,
                        score,
                    });
                }
            }
        }
        Ok(matches)
    }

    fn score_at(
        &self,
        sequence: &Sequence<A>,
        start: usize,
        strand: Strand,
    ) -> Result<f32, MotifError> {
        let mut total = 0.0;
        for offset in 0..self.width() {
            let byte = sequence.as_slice()[start + offset].as_ascii();
            let Some(base) = concrete_index(byte) else {
                return Err(MotifError::AmbiguousInput {
                    offset: start + offset,
                    byte,
                });
            };
            let cell = match strand {
                Strand::Forward => self.columns[offset][base],
                Strand::Reverse => self.columns[self.width() - 1 - offset][3 - base],
            };
            if cell == f64::NEG_INFINITY {
                return Ok(f32::NEG_INFINITY);
            }
            total += cell;
            if !total.is_finite() {
                return Err(MotifError::InvalidMatrix {
                    column: Some(offset),
                    reason: "window score accumulation overflows f64",
                });
            }
        }
        if total.abs() > f64::from(f32::MAX) {
            return Err(MotifError::InvalidMatrix {
                column: None,
                reason: "finite window score is outside the finite float32 range",
            });
        }
        // Rounding to the explicitly selected public scan precision is intended.
        #[allow(clippy::cast_possible_truncation)]
        let score = total as f32;
        Ok(score)
    }
}

/// Invalid motif input, matrix configuration or bounded-search request.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum MotifError {
    /// A pattern contains a symbol rejected by its nucleotide alphabet.
    InvalidPattern(SequenceError),
    /// An IUPAC pattern has no nucleotide positions.
    EmptyPattern,
    /// A workflow requiring concrete bases encountered ambiguity.
    AmbiguousInput {
        /// Zero-based original-input offset, or offset within a count instance.
        offset: usize,
        /// Rejected uppercase ASCII symbol.
        byte: u8,
    },
    /// A count, weight, background or score cannot define the specified model.
    InvalidMatrix {
        /// Zero-based affected column, when the error belongs to one column.
        column: Option<usize>,
        /// Stable explanatory text describing the invalid matrix condition.
        reason: &'static str,
    },
    /// A caller-supplied search option or arithmetic bound is invalid.
    InvalidOptions {
        /// Stable explanatory text describing the invalid option.
        reason: &'static str,
    },
    /// A motif instance or exact scoring window has a different width.
    DimensionMismatch {
        /// Required width.
        expected: usize,
        /// Supplied width.
        actual: usize,
        /// Zero-based count-instance index, absent for a scoring window.
        instance: Option<usize>,
    },
    /// The complete request would exceed a declared resource budget.
    LimitExceeded {
        /// Resource name: `windows`, `comparisons` or `matches`.
        resource: &'static str,
        /// Caller-selected maximum.
        limit: usize,
        /// Required count, or the first count exceeding a match budget.
        observed: usize,
    },
}

impl MotifError {
    /// Returns a stable machine-readable error code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::InvalidPattern(_) => "invalid_motif_pattern",
            Self::EmptyPattern => "empty_motif_pattern",
            Self::AmbiguousInput { .. } => "ambiguous_motif_input",
            Self::InvalidMatrix { .. } => "invalid_motif_matrix",
            Self::InvalidOptions { .. } => "invalid_arguments",
            Self::DimensionMismatch { .. } => "motif_dimension_mismatch",
            Self::LimitExceeded { .. } => "limit_exceeded",
        }
    }
}

impl fmt::Display for MotifError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPattern(source) => write!(formatter, "invalid motif pattern: {source}"),
            Self::EmptyPattern => formatter.write_str("motif pattern must be nonempty"),
            Self::AmbiguousInput { offset, byte } => write!(
                formatter,
                "motif workflow requires concrete bases; found {} at nucleotide {offset}",
                char::from(*byte)
            ),
            Self::InvalidMatrix { column, reason } => match column {
                Some(column) => write!(formatter, "invalid motif matrix column {column}: {reason}"),
                None => write!(formatter, "invalid motif matrix: {reason}"),
            },
            Self::InvalidOptions { reason } => write!(formatter, "invalid motif options: {reason}"),
            Self::DimensionMismatch {
                expected,
                actual,
                instance,
            } => match instance {
                Some(instance) => write!(
                    formatter,
                    "motif instance {instance} has width {actual}, expected {expected}"
                ),
                None => write!(
                    formatter,
                    "motif window has width {actual}, expected {expected}"
                ),
            },
            Self::LimitExceeded {
                resource,
                limit,
                observed,
            } => write!(
                formatter,
                "motif {resource} count {observed} exceeds limit {limit}"
            ),
        }
    }
}

impl Error for MotifError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::InvalidPattern(source) => Some(source),
            _ => None,
        }
    }
}

fn base_mask(byte: u8) -> u8 {
    match byte {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}

fn complement_mask(mask: u8) -> u8 {
    ((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3)
}

fn concrete_index(byte: u8) -> Option<usize> {
    match byte {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

fn validate_concrete<A: MotifAlphabet>(input: &Sequence<A>) -> Result<(), MotifError> {
    for (offset, symbol) in input.iter().enumerate() {
        let byte = symbol.as_ascii();
        if concrete_index(byte).is_none() {
            return Err(MotifError::AmbiguousInput { offset, byte });
        }
    }
    Ok(())
}

fn selected_strands(selection: StrandSelection) -> &'static [Strand] {
    match selection {
        StrandSelection::Both => &[Strand::Forward, Strand::Reverse],
        StrandSelection::Forward => &[Strand::Forward],
        StrandSelection::Reverse => &[Strand::Reverse],
    }
}

fn window_count(length: usize, width: usize) -> usize {
    if length < width {
        0
    } else {
        // Width is nonzero for every constructible motif or matrix.
        length - width + 1
    }
}

fn search_windows(
    length: usize,
    width: usize,
    strands: StrandSelection,
    max_windows: usize,
    max_comparisons: usize,
) -> Result<usize, MotifError> {
    let windows = window_count(length, width);
    let total_windows =
        windows
            .checked_mul(selected_strands(strands).len())
            .ok_or(MotifError::InvalidOptions {
                reason: "total window count overflows usize",
            })?;
    check_limit("windows", total_windows, max_windows)?;
    let comparisons = total_windows
        .checked_mul(width)
        .ok_or(MotifError::InvalidOptions {
            reason: "total comparison count overflows usize",
        })?;
    check_limit("comparisons", comparisons, max_comparisons)?;
    Ok(windows)
}

fn check_limit(resource: &'static str, observed: usize, limit: usize) -> Result<(), MotifError> {
    if observed > limit {
        Err(MotifError::LimitExceeded {
            resource,
            limit,
            observed,
        })
    } else {
        Ok(())
    }
}

fn check_next_match(current: usize, limit: usize) -> Result<(), MotifError> {
    let observed = current.checked_add(1).ok_or(MotifError::InvalidOptions {
        reason: "match count overflows usize",
    })?;
    check_limit("matches", observed, limit)
}

fn validate_width(width: usize) -> Result<(), MotifError> {
    if width == 0 {
        Err(MotifError::InvalidMatrix {
            column: None,
            reason: "matrix must have at least one column",
        })
    } else {
        Ok(())
    }
}

fn validate_nonnegative(values: &[f64; 4], column: Option<usize>) -> Result<(), MotifError> {
    if values
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        Err(MotifError::InvalidMatrix {
            column,
            reason: "counts and pseudocounts must be finite and nonnegative",
        })
    } else {
        Ok(())
    }
}

fn finite_sum(values: &[f64; 4], column: Option<usize>) -> Result<f64, MotifError> {
    let mut total: f64 = 0.0;
    let mut correction = 0.0;
    for value in values {
        let next = total + value;
        if !next.is_finite() {
            return Err(MotifError::InvalidMatrix {
                column,
                reason: "matrix normalization total overflows f64",
            });
        }
        // Neumaier compensation retains low-order terms before the final
        // rounding. Column/background normalization is distinct from the
        // deliberately sequential accumulation used for scoring a window.
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    let corrected = total + correction;
    if corrected.is_finite() {
        Ok(corrected)
    } else {
        Err(MotifError::InvalidMatrix {
            column,
            reason: "matrix normalization total overflows f64",
        })
    }
}

fn normalize_background(mut background: [f64; 4]) -> Result<[f64; 4], MotifError> {
    if background
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(MotifError::InvalidMatrix {
            column: None,
            reason: "background weights must be finite and strictly positive",
        });
    }
    let total = finite_sum(&background, None)?;
    for value in &mut background {
        *value /= total;
        if *value == 0.0 {
            return Err(MotifError::InvalidMatrix {
                column: None,
                reason: "normalized background probability underflows to zero",
            });
        }
    }
    Ok(background)
}

#[cfg(test)]
// The rational/integer examples and selected scan-precision probes require exact
// values; a tolerance would hide a changed arithmetic or rounding contract.
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;
    use serde_json::Value;

    fn dna(sequence: &str) -> Sequence<Dna> {
        Sequence::from_ascii(sequence).expect("test DNA must be valid")
    }

    fn forward(ambiguity: MotifInputAmbiguity) -> MotifSearchOptions {
        MotifSearchOptions {
            strands: StrandSelection::Forward,
            ambiguity,
            ..MotifSearchOptions::default()
        }
    }

    fn summary(matches: &[MotifMatch]) -> Vec<(Range<usize>, Strand)> {
        matches
            .iter()
            .map(|hit| (hit.range(), hit.strand()))
            .collect()
    }

    fn worked_scores() -> MotifScores<Dna> {
        MotifCounts::from_instances(&[dna("AC"), dna("AC"), dna("AC"), dna("CA")])
            .expect("equal concrete instances")
            .normalize([1.0; 4])
            .expect("positive pseudocounts")
            .log_odds_uniform()
            .expect("positive uniform background")
    }

    #[test]
    fn patterns_validate_normalize_and_preserve_molecule() {
        assert_eq!(dna_motif("arn").expect("valid motif").pattern(), "ARN");
        assert_eq!(rna_motif("au").expect("valid RNA").width(), 2);
        assert!(matches!(dna_motif(""), Err(MotifError::EmptyPattern)));
        assert!(matches!(
            dna_motif("AU"),
            Err(MotifError::InvalidPattern(_))
        ));
        assert!(matches!(
            rna_motif("AT"),
            Err(MotifError::InvalidPattern(_))
        ));
        let error = dna_motif("A-").expect_err("gap is invalid");
        assert_eq!(error.code(), "invalid_motif_pattern");
        assert!(error.source().is_some());
    }

    #[test]
    fn worked_ambiguity_cases_and_both_strand_order() {
        let motif = dna_motif("R").expect("valid R");
        let sequence = dna("AGRYN");
        assert_eq!(
            summary(&motif.find_in(&sequence).expect("default search")),
            vec![
                (0..1, Strand::Forward),
                (1..2, Strand::Forward),
                (2..3, Strand::Forward),
                (3..4, Strand::Reverse),
            ]
        );
        let possible = motif
            .find_with(
                &sequence,
                MotifSearchOptions {
                    ambiguity: MotifInputAmbiguity::Possible,
                    ..MotifSearchOptions::default()
                },
            )
            .expect("possible search");
        assert_eq!(
            summary(&possible)[4..],
            [(4..5, Strand::Forward), (4..5, Strand::Reverse)]
        );
    }

    #[test]
    fn all_iupac_pairs_follow_independent_explicit_sets() {
        let symbols: &[(u8, &[u8])] = &[
            (b'A', b"A"),
            (b'C', b"C"),
            (b'G', b"G"),
            (b'T', b"T"),
            (b'R', b"AG"),
            (b'Y', b"CT"),
            (b'S', b"CG"),
            (b'W', b"AT"),
            (b'K', b"GT"),
            (b'M', b"AC"),
            (b'B', b"CGT"),
            (b'D', b"AGT"),
            (b'H', b"ACT"),
            (b'V', b"ACG"),
            (b'N', b"ACGT"),
        ];
        for (pattern, accepted) in symbols {
            let motif = dna_motif([*pattern]).expect("valid pattern");
            for (input, possibilities) in symbols {
                let sequence = Sequence::<Dna>::from_ascii([*input]).expect("valid input");
                let definite = possibilities.iter().all(|base| accepted.contains(base));
                let possible = possibilities.iter().any(|base| accepted.contains(base));
                for (ambiguity, expected) in [
                    (MotifInputAmbiguity::Definite, definite),
                    (MotifInputAmbiguity::Possible, possible),
                ] {
                    assert_eq!(
                        !motif
                            .find_with(&sequence, forward(ambiguity))
                            .expect("bounded pair")
                            .is_empty(),
                        expected,
                        "pattern={} input={} mode={ambiguity:?}",
                        char::from(*pattern),
                        char::from(*input)
                    );
                }
            }
        }
    }

    #[test]
    fn overlap_palindromes_and_original_reverse_positions() {
        let overlap = dna_motif("AA")
            .expect("valid motif")
            .find_with(&dna("AAAA"), forward(MotifInputAmbiguity::Definite))
            .expect("overlap search");
        assert_eq!(
            summary(&overlap),
            vec![
                (0..2, Strand::Forward),
                (1..3, Strand::Forward),
                (2..4, Strand::Forward)
            ]
        );
        let palindrome = dna_motif("AT")
            .expect("valid motif")
            .find_in(&dna("ATAT"))
            .expect("palindrome search");
        assert_eq!(
            summary(&palindrome),
            vec![
                (0..2, Strand::Forward),
                (0..2, Strand::Reverse),
                (2..4, Strand::Forward),
                (2..4, Strand::Reverse)
            ]
        );
        let reverse = dna_motif("AC")
            .expect("valid motif")
            .find_in(&dna("TTGTAA"))
            .expect("reverse search");
        assert_eq!(summary(&reverse), vec![(2..4, Strand::Reverse)]);
    }

    #[test]
    fn every_iupac_reverse_complement_uses_the_declared_symbol_sets() {
        let complements = [
            (b'A', b'T'),
            (b'C', b'G'),
            (b'G', b'C'),
            (b'T', b'A'),
            (b'R', b'Y'),
            (b'Y', b'R'),
            (b'S', b'S'),
            (b'W', b'W'),
            (b'K', b'M'),
            (b'M', b'K'),
            (b'B', b'V'),
            (b'D', b'H'),
            (b'H', b'D'),
            (b'V', b'B'),
            (b'N', b'N'),
        ];
        for (symbol, complement) in complements {
            for (second, second_complement) in complements {
                let motif = dna_motif([symbol, second]).expect("IUPAC pair");
                let source = Sequence::<Dna>::from_ascii([second_complement, complement])
                    .expect("independently supplied complement pair");
                let hits = motif
                    .find_with(
                        &source,
                        MotifSearchOptions {
                            strands: StrandSelection::Reverse,
                            ..MotifSearchOptions::default()
                        },
                    )
                    .expect("reverse definite search");
                assert_eq!(summary(&hits), [(0..2, Strand::Reverse)]);
            }
        }
    }

    #[test]
    fn empty_short_inputs_and_global_rejection() {
        let motif = dna_motif("ACG").expect("valid motif");
        assert!(motif.find_in(&dna("")).expect("empty source").is_empty());
        assert!(motif.find_in(&dna("AC")).expect("short source").is_empty());
        assert!(matches!(
            motif.find_with(&dna("AN"), forward(MotifInputAmbiguity::Reject)),
            Err(MotifError::AmbiguousInput {
                offset: 1,
                byte: b'N'
            })
        ));
        let rna = rna_motif("AU").expect("RNA motif");
        let source = Sequence::<Rna>::from_ascii("AUAU").expect("RNA source");
        assert_eq!(rna.find_in(&source).expect("RNA search").len(), 4);
    }

    #[test]
    fn pattern_budgets_count_both_strands_and_worst_case_comparisons() {
        let motif = dna_motif("AA").expect("motif");
        let input = dna("CCCC");
        assert!(matches!(
            motif.find_with(
                &input,
                MotifSearchOptions {
                    max_windows: 5,
                    ..MotifSearchOptions::default()
                }
            ),
            Err(MotifError::LimitExceeded {
                resource: "windows",
                observed: 6,
                ..
            })
        ));
        assert!(matches!(
            motif.find_with(
                &input,
                MotifSearchOptions {
                    max_comparisons: 11,
                    ..MotifSearchOptions::default()
                }
            ),
            Err(MotifError::LimitExceeded {
                resource: "comparisons",
                observed: 12,
                ..
            })
        ));
        assert!(
            motif
                .find_with(
                    &input,
                    MotifSearchOptions {
                        max_matches: 0,
                        ..MotifSearchOptions::default()
                    }
                )
                .expect("no matches within zero budget")
                .is_empty()
        );
        assert!(matches!(
            motif.find_with(
                &dna("AAAA"),
                MotifSearchOptions {
                    max_matches: 2,
                    ..forward(MotifInputAmbiguity::Definite)
                }
            ),
            Err(MotifError::LimitExceeded {
                resource: "matches",
                observed: 3,
                ..
            })
        ));
        assert!(
            search_windows(usize::MAX, 1, StrandSelection::Both, usize::MAX, usize::MAX).is_err()
        );
        assert!(
            search_windows(
                usize::MAX,
                2,
                StrandSelection::Forward,
                usize::MAX,
                usize::MAX
            )
            .is_err()
        );
    }

    #[test]
    fn counts_weights_and_integer_scores_match_hand_calculation() {
        let counts =
            MotifCounts::<Dna>::from_instances(&[dna("AC"), dna("AC"), dna("AC"), dna("CA")])
                .expect("instances");
        assert_eq!(
            counts.columns(),
            &[[3.0, 1.0, 0.0, 0.0], [1.0, 3.0, 0.0, 0.0]]
        );
        let weights = counts.normalize([1.0; 4]).expect("pseudocounts");
        assert_eq!(
            weights.columns(),
            &[[0.5, 0.25, 0.125, 0.125], [0.25, 0.5, 0.125, 0.125]]
        );
        let scores = weights.log_odds_uniform().expect("log odds");
        assert_eq!(
            scores.columns(),
            &[[1.0, 0.0, -1.0, -1.0], [0.0, 1.0, -1.0, -1.0]]
        );
        assert_eq!(
            scores.calculate(&dna("ACGTAC"), 5).expect("scores"),
            [2.0, -1.0, -2.0, -1.0, 2.0]
        );
        let matches = scores
            .find_with(
                &dna("ACGTAC"),
                MatrixSearchOptions {
                    threshold: 2.0,
                    ..MatrixSearchOptions::default()
                },
            )
            .expect("inclusive scan");
        assert_eq!(
            matches
                .iter()
                .map(|hit| (hit.range(), hit.strand(), hit.score()))
                .collect::<Vec<_>>(),
            vec![
                (0..2, Strand::Forward, 2.0),
                (2..4, Strand::Reverse, 2.0),
                (4..6, Strand::Forward, 2.0)
            ]
        );
    }

    #[test]
    fn unequal_column_totals_weighted_counts_and_rna_work() {
        let counts =
            MotifCounts::<Rna>::from_columns(vec![[0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.0, 10.0]])
                .expect("weighted columns");
        let weights = counts.normalize_default().expect("column normalization");
        assert_eq!(
            weights.columns(),
            &[[0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
        );
        let scores = weights.log_odds_uniform().expect("scores");
        let source = Sequence::<Rna>::from_ascii("AU").expect("RNA");
        assert_eq!(scores.score_window(&source).expect("RNA score"), 3.0);
        assert_eq!(weights.width(), 2);
    }

    #[test]
    fn zero_columns_pseudocounts_and_negative_infinity_are_explicit() {
        let zero = MotifCounts::<Dna>::from_columns(vec![[0.0; 4]]).expect("zero counts retained");
        assert!(matches!(
            zero.normalize_default(),
            Err(MotifError::InvalidMatrix {
                column: Some(0),
                ..
            })
        ));
        assert_eq!(
            zero.normalize([1.0; 4])
                .expect("positive pseudocounts")
                .columns(),
            &[[0.25; 4]]
        );
        let scores = MotifCounts::<Dna>::from_columns(vec![[1.0, 0.0, 0.0, 0.0]])
            .expect("counts")
            .normalize_default()
            .expect("nonzero total")
            .log_odds_uniform()
            .expect("scores");
        assert_eq!(
            scores.columns(),
            &[[2.0, f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY]]
        );
        assert_eq!(
            scores
                .calculate(&dna("AC"), 2)
                .expect("negative infinity score"),
            [2.0, f32::NEG_INFINITY]
        );
        assert_eq!(
            scores
                .find_with(
                    &dna("AC"),
                    MatrixSearchOptions {
                        strands: StrandSelection::Forward,
                        ..MatrixSearchOptions::default()
                    }
                )
                .expect("finite threshold")
                .len(),
            1
        );
    }

    #[test]
    fn asymmetric_background_reverse_scans_complement_existing_scores() {
        let scores = MotifCounts::<Dna>::from_columns(vec![[1.0; 4]])
            .expect("counts")
            .normalize_default()
            .expect("weights")
            .log_odds([4.0, 2.0, 1.0, 1.0])
            .expect("positive background");
        assert_eq!(scores.background(), [0.5, 0.25, 0.125, 0.125]);
        assert_eq!(scores.columns(), &[[-1.0, 0.0, 1.0, 1.0]]);
        let hits = scores
            .find_with(&dna("A"), MatrixSearchOptions::default())
            .expect("reverse hit");
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].strand(), Strand::Reverse);
        assert_eq!(hits[0].score(), 1.0);
    }

    #[test]
    fn malformed_counts_instances_backgrounds_and_scores_fail() {
        assert!(MotifCounts::<Dna>::from_columns(vec![]).is_err());
        assert!(MotifCounts::<Dna>::from_instances(&[]).is_err());
        assert!(MotifCounts::<Dna>::from_instances(&[dna("")]).is_err());
        assert!(matches!(
            MotifCounts::<Dna>::from_instances(&[dna("AC"), dna("A")]),
            Err(MotifError::DimensionMismatch {
                instance: Some(1),
                ..
            })
        ));
        assert!(matches!(
            MotifCounts::<Dna>::from_instances(&[dna("AN")]),
            Err(MotifError::AmbiguousInput { offset: 1, .. })
        ));
        for bad in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(MotifCounts::<Dna>::from_columns(vec![[bad, 0.0, 0.0, 0.0]]).is_err());
        }
        assert!(MotifCounts::<Dna>::from_columns(vec![[f64::MAX; 4]]).is_err());
        let counts = MotifCounts::<Dna>::from_columns(vec![[1.0; 4]]).expect("counts");
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(counts.normalize([bad; 4]).is_err());
        }
        let weights = counts.normalize_default().expect("weights");
        for background in [
            [0.0, 1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0, 1.0],
            [f64::NAN; 4],
            [f64::MAX; 4],
            [f64::from_bits(1), f64::MAX, 1.0, 1.0],
        ] {
            assert!(weights.log_odds(background).is_err());
        }
        for bad in [f64::NAN, f64::INFINITY] {
            assert!(MotifScores::<Dna>::from_columns(vec![[bad; 4]], [1.0; 4]).is_err());
        }
        assert!(MotifScores::<Dna>::from_columns(vec![[f64::NEG_INFINITY; 4]], [1.0; 4]).is_ok());
    }

    #[test]
    fn matrix_validation_is_global_and_budget_failures_do_not_truncate() {
        let scores = worked_scores();
        assert!(
            scores
                .calculate(&dna(""), 0)
                .expect("empty source")
                .is_empty()
        );
        assert!(
            scores
                .calculate(&dna("A"), 0)
                .expect("short source")
                .is_empty()
        );
        assert!(matches!(
            scores.calculate(&dna("N"), 0),
            Err(MotifError::AmbiguousInput { offset: 0, .. })
        ));
        assert!(matches!(
            scores.score_window(&dna("A")),
            Err(MotifError::DimensionMismatch { actual: 1, .. })
        ));
        assert!(matches!(
            scores.calculate(&dna("ACG"), 1),
            Err(MotifError::LimitExceeded {
                resource: "windows",
                observed: 2,
                ..
            })
        ));
        assert!(matches!(
            scores.find_with(
                &dna("ACGTAC"),
                MatrixSearchOptions {
                    max_windows: 9,
                    ..MatrixSearchOptions::default()
                }
            ),
            Err(MotifError::LimitExceeded {
                resource: "windows",
                observed: 10,
                ..
            })
        ));
        assert!(matches!(
            scores.find_with(
                &dna("ACGTAC"),
                MatrixSearchOptions {
                    max_comparisons: 19,
                    ..MatrixSearchOptions::default()
                }
            ),
            Err(MotifError::LimitExceeded {
                resource: "comparisons",
                observed: 20,
                ..
            })
        ));
        assert!(matches!(
            scores.find_with(
                &dna("AC"),
                MatrixSearchOptions {
                    max_matches: 0,
                    ..MatrixSearchOptions::default()
                }
            ),
            Err(MotifError::LimitExceeded {
                resource: "matches",
                observed: 1,
                ..
            })
        ));
        assert!(
            scores
                .find_with(
                    &dna("CC"),
                    MatrixSearchOptions {
                        threshold: 2.0,
                        max_matches: 0,
                        ..MatrixSearchOptions::default()
                    }
                )
                .expect("no hits allowed")
                .is_empty()
        );
        for threshold in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, f64::MAX] {
            assert!(matches!(
                scores.find_with(
                    &dna("AC"),
                    MatrixSearchOptions {
                        threshold,
                        ..MatrixSearchOptions::default()
                    }
                ),
                Err(MotifError::InvalidOptions { .. })
            ));
        }
    }

    #[test]
    fn scan_precision_uses_double_accumulation_then_single_rounding() {
        let scores = MotifScores::<Dna>::from_columns(
            vec![[1.0; 4], [2_f64.powi(-25); 4], [-1.0; 4]],
            [1.0; 4],
        )
        .expect("explicit scores");
        assert_eq!(
            f64::from(scores.score_window(&dna("AAA")).expect("score")),
            2_f64.powi(-25)
        );
        let outside =
            MotifScores::<Dna>::from_columns(vec![[f64::MAX; 4]], [1.0; 4]).expect("finite cells");
        assert!(outside.score_window(&dna("A")).is_err());
        let overflow =
            MotifScores::<Dna>::from_columns(vec![[f64::MAX; 4], [f64::MAX; 4]], [1.0; 4])
                .expect("finite cells");
        assert!(overflow.score_window(&dna("AA")).is_err());
    }

    #[test]
    fn non_power_of_two_log_odds_preserve_the_reference_log_quotient() {
        let scores = MotifCounts::<Dna>::from_columns(vec![[1.0, 2.0, 3.0, 5.0]])
            .expect("unequal concrete counts")
            .normalize_default()
            .expect("positive column total")
            .log_odds([0.3, 0.2, 0.2, 0.3])
            .expect("positive background");
        // Biopython 1.88 / CPython 3.12.3 / NumPy 2.3.3 public-API observation:
        // A score = -0x1.b8f3888b9310ap+0. A direct log2 differs by one bit.
        assert_eq!(scores.columns()[0][0].to_bits(), 0xbffb_8f38_88b9_310a);
    }

    #[test]
    fn thresholds_round_to_scan_precision_before_inclusive_comparison() {
        let scores = MotifScores::<Dna>::from_columns(vec![[0.1; 4]], [1.0; 4])
            .expect("finite explicit scores");
        let score = scores.score_window(&dna("A")).expect("one-window score");
        assert_eq!(score.to_bits(), 0x3dcc_cccd);
        let rounded_back_threshold = f64::from(score) + 1e-12;
        assert!(rounded_back_threshold > f64::from(score));
        let hits = scores
            .find_with(
                &dna("AA"),
                MatrixSearchOptions {
                    strands: StrandSelection::Forward,
                    threshold: rounded_back_threshold,
                    ..MatrixSearchOptions::default()
                },
            )
            .expect("threshold rounds back to the same score");
        assert_eq!(hits.len(), 2);
        let next_score = f32::from_bits(score.to_bits() + 1);
        assert!(
            scores
                .find_with(
                    &dna("AA"),
                    MatrixSearchOptions {
                        strands: StrandSelection::Forward,
                        threshold: f64::from(next_score),
                        ..MatrixSearchOptions::default()
                    },
                )
                .expect("next representable threshold")
                .is_empty()
        );
    }

    fn reference_fixture() -> Value {
        serde_json::from_str(include_str!("../tests/fixtures/motifs-biopython-1.88.json"))
            .expect("captured reference fixture must remain valid JSON")
    }

    fn reference_f64(value: &Value) -> f64 {
        f64::from_bits(
            u64::from_str_radix(value.as_str().expect("binary64 string"), 16)
                .expect("16 hexadecimal binary64 digits"),
        )
    }

    fn reference_column(value: &Value) -> [f64; 4] {
        let values = value.as_array().expect("nucleotide column");
        assert_eq!(values.len(), 4);
        std::array::from_fn(|base| reference_f64(&values[base]))
    }

    fn assert_matrix_bits(actual: &[[f64; 4]], expected: &Value, id: &str, label: &str) -> usize {
        let expected = expected.as_array().expect("matrix columns");
        assert_eq!(actual.len(), expected.len(), "{id} {label} width");
        for (column, (actual, expected)) in actual.iter().zip(expected).enumerate() {
            for base in 0..4 {
                assert_eq!(
                    actual[base].to_bits(),
                    reference_f64(&expected[base]).to_bits(),
                    "{id} {label} column={column} base={base}"
                );
            }
        }
        actual.len() * 4
    }

    fn assert_scan_bits(actual: &[f32], expected: &Value, id: &str, label: &str) -> usize {
        let expected = expected.as_array().expect("scan score array");
        assert_eq!(actual.len(), expected.len(), "{id} {label} windows");
        for (window, (actual, expected)) in actual.iter().zip(expected).enumerate() {
            let bits = u32::from_str_radix(expected.as_str().expect("binary32 string"), 16)
                .expect("8 hexadecimal binary32 digits");
            assert_eq!(actual.to_bits(), bits, "{id} {label} window={window}");
        }
        actual.len()
    }

    fn assert_native_hit_bits(actual: &[ScoredMotifMatch], case: &Value, id: &str) -> usize {
        let expected = case["native_hits"].as_array().expect("native hits");
        assert_eq!(actual.len(), expected.len(), "{id} hit count");
        let input_length = i64::try_from(case["sequence"].as_str().expect("source").len())
            .expect("short captured source");
        for (actual, expected) in actual.iter().zip(expected) {
            let start = usize::try_from(expected["start"].as_u64().expect("hit start"))
                .expect("source position fits usize");
            let end = usize::try_from(expected["end"].as_u64().expect("hit end"))
                .expect("source position fits usize");
            let strand = match expected["strand"].as_str().expect("strand label") {
                "forward" => Strand::Forward,
                "reverse" => Strand::Reverse,
                other => panic!("invalid captured strand {other}"),
            };
            let native = expected["native_position"]
                .as_i64()
                .expect("native position");
            let native_original = if native < 0 {
                input_length + native
            } else {
                native
            };
            assert_eq!(
                usize::try_from(native_original).expect("original position"),
                start
            );
            assert_eq!(actual.range(), start..end, "{id} hit range");
            assert_eq!(actual.strand(), strand, "{id} hit strand");
            let score = u32::from_str_radix(
                expected["score_f32_bits"].as_str().expect("hit binary32"),
                16,
            )
            .expect("hit hexadecimal bits");
            assert_eq!(actual.score().to_bits(), score, "{id} hit score");
        }
        actual.len()
    }

    #[test]
    fn pinned_public_reference_matrix_and_scan_corpus_matches_exact_bits() {
        let fixture = reference_fixture();
        assert_eq!(fixture["reference"]["biopython"], "1.88");
        assert_eq!(fixture["reference"]["numpy"], "2.3.3");
        assert_eq!(fixture["reference"]["python"], "3.12.3");
        let cases = fixture["cases"].as_array().expect("captured matrix cases");
        assert_eq!(cases.len(), 360);
        let (mut errors, mut cells, mut scans, mut hits) = (0, 0, 0, 0);
        for case in cases {
            let id = case["id"].as_str().expect("case identifier");
            let columns = case["counts_f64_bits"]
                .as_array()
                .expect("captured input columns")
                .iter()
                .map(reference_column)
                .collect();
            let counts = MotifCounts::<Dna>::from_columns(columns).expect("finite counts");
            let weights = counts.normalize(reference_column(&case["pseudocounts_f64_bits"]));
            if !case["reference_error"].is_null() {
                assert_eq!(case["reference_error"]["type"], "ZeroDivisionError");
                assert!(
                    matches!(weights, Err(MotifError::InvalidMatrix { .. })),
                    "{id}"
                );
                errors += 1;
                continue;
            }
            let weights = weights.expect("reference normalization succeeds");
            cells += assert_matrix_bits(weights.columns(), &case["weights_f64_bits"], id, "PWM");
            let scores = weights
                .log_odds(reference_column(&case["background_f64_bits"]))
                .expect("positive background");
            cells += assert_matrix_bits(scores.columns(), &case["scores_f64_bits"], id, "PSSM");
            let source = dna(case["sequence"].as_str().expect("source sequence"));
            scans += assert_scan_bits(
                &scores.calculate(&source, 1000).expect("forward scan"),
                &case["forward_scores_f32_bits"],
                id,
                "forward",
            );
            let reverse: Vec<_> = (0..window_count(source.len(), scores.width()))
                .map(|start| {
                    scores
                        .score_at(&source, start, Strand::Reverse)
                        .expect("reverse score")
                })
                .collect();
            scans += assert_scan_bits(&reverse, &case["reverse_scores_f32_bits"], id, "reverse");
            let found = scores
                .find_with(
                    &source,
                    MatrixSearchOptions {
                        threshold: reference_f64(&case["threshold_f64_bits"]),
                        ..MatrixSearchOptions::default()
                    },
                )
                .expect("bounded threshold scan");
            hits += assert_native_hit_bits(&found, case, id);
        }
        assert_eq!((errors, cells, scans, hits), (18, 10_488, 28_158, 9127));
    }

    #[test]
    fn reference_precision_probes_and_documented_one_window_difference() {
        let fixture = reference_fixture();
        for (id, case) in fixture["probes"]
            .as_object()
            .expect("public numerical probes")
        {
            let columns = case["columns"]
                .as_array()
                .expect("explicit probe columns")
                .iter()
                .map(|column| {
                    std::array::from_fn(|base| column[base].as_f64().expect("probe cell"))
                })
                .collect();
            let scores =
                MotifScores::<Dna>::from_columns(columns, [1.0; 4]).expect("finite probe cells");
            let source = dna(case["sequence"].as_str().expect("probe source"));
            assert_scan_bits(
                &scores.calculate(&source, 1000).expect("probe scores"),
                &case["forward_scores_f32_bits"],
                id,
                "forward",
            );
            let found = scores
                .find_with(
                    &source,
                    MatrixSearchOptions {
                        threshold: case["threshold"].as_f64().expect("finite probe threshold"),
                        ..MatrixSearchOptions::default()
                    },
                )
                .expect("valid probe search, including a single window");
            if case["reference_search_error"].is_null() {
                assert_native_hit_bits(&found, case, id);
            } else {
                // NumPy 2.3.3 rejects native search's zero-dimensional result.
                // BioFerris deliberately supports the mathematically valid case.
                assert_eq!(id, "one_window_native_search");
                assert_eq!(case["reference_search_error"]["type"], "ValueError");
                assert_eq!(found.len(), 2);
                assert_eq!(found[0].range(), 0..1);
                assert_eq!(found[1].range(), 0..1);
                assert_eq!(found[0].strand(), Strand::Forward);
                assert_eq!(found[1].strand(), Strand::Reverse);
            }
        }
    }
}
