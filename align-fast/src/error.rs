use core::fmt;

use align_core::{AlignmentOperation, AlignmentStructureError, Score};

/// Failure to configure or execute exact wavefront alignment.
#[derive(Debug)]
#[non_exhaustive]
pub enum WavefrontError {
    /// A gap-affine mismatch penalty must be strictly positive.
    NonPositiveMismatchPenalty {
        /// Supplied mismatch penalty.
        penalty: Score,
    },
    /// A gap-opening penalty must be nonnegative.
    NegativeGapOpeningPenalty {
        /// Supplied gap-opening penalty.
        penalty: Score,
    },
    /// A gap-extension penalty must be strictly positive.
    NonPositiveGapExtensionPenalty {
        /// Supplied gap-extension penalty.
        penalty: Score,
    },
    /// A modeled terminal wavefront penalty must be nonnegative.
    NegativeTerminalPenalty {
        /// Supplied terminal penalty.
        penalty: Score,
    },
    /// Input or derived wavefront dimensions overflowed an index type.
    DimensionOverflow {
        /// First sequence length.
        first_len: usize,
        /// Second sequence length.
        second_len: usize,
    },
    /// A penalty or its conventional negative alignment score was unrepresentable.
    PenaltyOverflow {
        /// Penalty at which arithmetic failed.
        penalty: Score,
    },
    /// A required allocation could not be reserved.
    AllocationFailed {
        /// Name of the allocation.
        buffer: &'static str,
        /// Number of requested additional elements.
        elements: usize,
    },
    /// Stored wavefront history did not admit the selected path.
    TracebackInvariant {
        /// Current nonnegative penalty.
        penalty: Score,
        /// Current matrix diagonal.
        diagonal: isize,
        /// Current horizontal offset when representable.
        offset: isize,
    },
    /// The produced operations failed structural validation.
    InvalidResult(AlignmentStructureError),
    /// Independent replay disagreed with the exact wavefront penalty.
    PenaltyMismatch {
        /// Penalty established by the exact wavefront search.
        expected: Score,
        /// Penalty obtained by independently replaying the operations.
        actual: Score,
    },
}

impl fmt::Display for WavefrontError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonPositiveMismatchPenalty { penalty } => write!(
                formatter,
                "wavefront mismatch penalty {penalty} is not strictly positive"
            ),
            Self::NegativeGapOpeningPenalty { penalty } => write!(
                formatter,
                "wavefront gap-opening penalty {penalty} is negative"
            ),
            Self::NonPositiveGapExtensionPenalty { penalty } => write!(
                formatter,
                "wavefront gap-extension penalty {penalty} is not strictly positive"
            ),
            Self::NegativeTerminalPenalty { penalty } => write!(
                formatter,
                "modeled terminal wavefront penalty {penalty} is negative"
            ),
            Self::DimensionOverflow {
                first_len,
                second_len,
            } => write!(
                formatter,
                "wavefront dimensions overflow for lengths {first_len} and {second_len}"
            ),
            Self::PenaltyOverflow { penalty } => write!(
                formatter,
                "wavefront penalty arithmetic overflowed from {penalty}"
            ),
            Self::AllocationFailed { buffer, elements } => write!(
                formatter,
                "could not reserve {elements} additional elements for {buffer}"
            ),
            Self::TracebackInvariant {
                penalty,
                diagonal,
                offset,
            } => write!(
                formatter,
                "wavefront traceback has no valid predecessor at penalty {penalty}, diagonal {diagonal}, offset {offset}"
            ),
            Self::InvalidResult(error) => write!(formatter, "invalid alignment result: {error}"),
            Self::PenaltyMismatch { expected, actual } => write!(
                formatter,
                "wavefront alignment replay penalty {actual} differs from exact penalty {expected}"
            ),
        }
    }
}

impl std::error::Error for WavefrontError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidResult(error) => Some(error),
            _ => None,
        }
    }
}

/// Failure to configure or execute affine X-drop extension.
#[derive(Debug)]
#[non_exhaustive]
pub enum XDropError {
    /// The X-drop threshold must be nonnegative.
    NegativeThreshold {
        /// Supplied threshold.
        threshold: Score,
    },
    /// A forced seed does not identify one symbol in each input.
    SeedOutOfBounds {
        /// Supplied first-sequence symbol index.
        first_index: usize,
        /// Supplied second-sequence symbol index.
        second_index: usize,
        /// First sequence length.
        first_len: usize,
        /// Second sequence length.
        second_len: usize,
    },
    /// Input or derived buffer dimensions overflowed `usize`.
    DimensionOverflow {
        /// First sequence length.
        first_len: usize,
        /// Second sequence length.
        second_len: usize,
    },
    /// A required allocation could not be reserved.
    AllocationFailed {
        /// Name of the allocation.
        buffer: &'static str,
        /// Number of requested additional elements.
        elements: usize,
    },
    /// A positive candidate made an X-drop score unrepresentable.
    ScoreOverflow {
        /// Consumed prefix length of the first extension view.
        first_prefix: usize,
        /// Consumed prefix length of the second extension view.
        second_prefix: usize,
    },
    /// Stored X-drop history did not admit the selected path.
    TracebackInvariant {
        /// Consumed prefix length of the first extension view.
        first_prefix: usize,
        /// Consumed prefix length of the second extension view.
        second_prefix: usize,
    },
    /// The produced operations failed structural validation.
    InvalidResult(AlignmentStructureError),
}

impl fmt::Display for XDropError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NegativeThreshold { threshold } => {
                write!(formatter, "X-drop threshold {threshold} is negative")
            }
            Self::SeedOutOfBounds {
                first_index,
                second_index,
                first_len,
                second_len,
            } => write!(
                formatter,
                "X-drop seed ({first_index}, {second_index}) is outside sequence lengths {first_len} and {second_len}"
            ),
            Self::DimensionOverflow {
                first_len,
                second_len,
            } => write!(
                formatter,
                "X-drop dimensions overflow for lengths {first_len} and {second_len}"
            ),
            Self::AllocationFailed { buffer, elements } => write!(
                formatter,
                "could not reserve {elements} additional elements for {buffer}"
            ),
            Self::ScoreOverflow {
                first_prefix,
                second_prefix,
            } => write!(
                formatter,
                "X-drop score overflow at extension prefix cell ({first_prefix}, {second_prefix})"
            ),
            Self::TracebackInvariant {
                first_prefix,
                second_prefix,
            } => write!(
                formatter,
                "X-drop history has no valid predecessor at extension prefix cell ({first_prefix}, {second_prefix})"
            ),
            Self::InvalidResult(error) => write!(formatter, "invalid alignment result: {error}"),
        }
    }
}

impl std::error::Error for XDropError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidResult(error) => Some(error),
            _ => None,
        }
    }
}

/// Failure to construct or execute the Myers bit-vector profile exactly.
#[derive(Debug)]
#[non_exhaustive]
pub enum MyersError {
    /// Pattern, text, or derived buffer dimensions overflowed `usize`.
    DimensionOverflow {
        /// Pattern length.
        pattern_len: usize,
        /// Text length.
        text_len: usize,
    },
    /// A required allocation could not be reserved.
    AllocationFailed {
        /// Name of the allocation.
        buffer: &'static str,
        /// Number of requested elements.
        elements: usize,
    },
    /// A distance cannot be represented as a negated [`Score`].
    ScoreOverflow {
        /// Exact nonnegative edit distance.
        distance: usize,
    },
    /// Published bit-vector differences did not admit a valid predecessor.
    TracebackInvariant {
        /// Dynamic-programming row.
        row: usize,
        /// Dynamic-programming column.
        column: usize,
    },
    /// A bottom-cell update contradicted the edit-distance bounds.
    StateInvariant {
        /// One-based text column being advanced.
        column: usize,
    },
    /// The produced operations failed structural validation.
    InvalidResult(AlignmentStructureError),
}

impl fmt::Display for MyersError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DimensionOverflow {
                pattern_len,
                text_len,
            } => write!(
                formatter,
                "Myers dimensions overflow for pattern length {pattern_len} and text length {text_len}"
            ),
            Self::AllocationFailed { buffer, elements } => write!(
                formatter,
                "could not reserve {elements} elements for {buffer}"
            ),
            Self::ScoreOverflow { distance } => write!(
                formatter,
                "edit distance {distance} cannot be represented as a negated alignment score"
            ),
            Self::TracebackInvariant { row, column } => write!(
                formatter,
                "Myers history has no valid predecessor at dynamic-programming cell ({row}, {column})"
            ),
            Self::StateInvariant { column } => write!(
                formatter,
                "Myers bottom score became invalid while advancing text column {column}"
            ),
            Self::InvalidResult(error) => write!(formatter, "invalid alignment result: {error}"),
        }
    }
}

impl std::error::Error for MyersError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidResult(error) => Some(error),
            _ => None,
        }
    }
}

/// Failure while independently replaying a unit-edit alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum UnitEditVerificationError {
    /// The alignment does not consume a complete supplied sequence.
    NonGlobalRange {
        /// Sequence label, either `pattern` or `text`.
        sequence: &'static str,
        /// Observed range start.
        start: usize,
        /// Observed exclusive range end.
        end: usize,
        /// Required sequence length.
        length: usize,
    },
    /// An operation attempted to consume beyond a sequence end.
    SequenceExhausted {
        /// Sequence label, either `pattern` or `text`.
        sequence: &'static str,
        /// Alignment column containing the error.
        column: usize,
    },
    /// `Match` or `Substitution` disagreed with literal symbol equality.
    OperationKindMismatch {
        /// Alignment column containing the error.
        column: usize,
        /// Recorded operation.
        operation: AlignmentOperation,
    },
    /// The independently counted distance cannot be represented by [`Score`].
    ScoreOverflow {
        /// Independently counted edit distance.
        distance: usize,
    },
    /// The reported score is not the negative independently counted distance.
    ScoreMismatch {
        /// Score stored in the result.
        reported: Score,
        /// Required score under the unit-edit profile.
        recomputed: Score,
    },
}

impl fmt::Display for UnitEditVerificationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonGlobalRange {
                sequence,
                start,
                end,
                length,
            } => write!(
                formatter,
                "{sequence} range {start}..{end} is not the required global range 0..{length}"
            ),
            Self::SequenceExhausted { sequence, column } => write!(
                formatter,
                "alignment column {column} reads beyond the {sequence} sequence"
            ),
            Self::OperationKindMismatch { column, operation } => write!(
                formatter,
                "alignment column {column} is recorded as {operation:?}, inconsistent with its symbols"
            ),
            Self::ScoreOverflow { distance } => write!(
                formatter,
                "edit distance {distance} cannot be represented as a negated alignment score"
            ),
            Self::ScoreMismatch {
                reported,
                recomputed,
            } => write!(
                formatter,
                "reported score {reported} differs from unit-edit score {recomputed}"
            ),
        }
    }
}

impl std::error::Error for UnitEditVerificationError {}
