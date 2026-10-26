use core::fmt;

use align_core::{AlignmentOperation, AlignmentStructureError, Score};

/// Failure to compute a pairwise alignment exactly.
#[derive(Debug)]
#[non_exhaustive]
pub enum AlignmentError {
    /// Matrix or traceback dimensions overflowed `usize`.
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
        /// Number of requested elements.
        elements: usize,
    },
    /// An exact score could not be represented by [`Score`].
    ScoreOverflow {
        /// Dynamic-programming row.
        row: usize,
        /// Dynamic-programming column.
        column: usize,
    },
    /// A caller-required affine SIMD backend cannot execute this alignment exactly.
    SimdUnavailable {
        /// Stable diagnostic explaining why exact SIMD execution was rejected.
        reason: &'static str,
    },
    /// The traceback matrix contained an impossible predecessor.
    TracebackInvariant {
        /// Dynamic-programming row.
        row: usize,
        /// Dynamic-programming column.
        column: usize,
    },
    /// The produced operations failed structural validation.
    InvalidResult(AlignmentStructureError),
}

impl fmt::Display for AlignmentError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DimensionOverflow {
                first_len,
                second_len,
            } => write!(
                formatter,
                "alignment dimensions overflow for sequence lengths {first_len} and {second_len}"
            ),
            Self::AllocationFailed { buffer, elements } => {
                write!(
                    formatter,
                    "could not reserve {elements} elements for {buffer}"
                )
            }
            Self::ScoreOverflow { row, column } => write!(
                formatter,
                "alignment score overflow at dynamic-programming cell ({row}, {column})"
            ),
            Self::SimdUnavailable { reason } => {
                write!(formatter, "affine SIMD backend is unavailable: {reason}")
            }
            Self::TracebackInvariant { row, column } => write!(
                formatter,
                "invalid traceback predecessor at dynamic-programming cell ({row}, {column})"
            ),
            Self::InvalidResult(error) => write!(formatter, "invalid alignment result: {error}"),
        }
    }
}

impl std::error::Error for AlignmentError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidResult(error) => Some(error),
            _ => None,
        }
    }
}

/// Failure while independently replaying and scoring an alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum VerificationError {
    /// Needleman–Wunsch must consume both complete input sequences.
    NonGlobalRange {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Observed range start.
        start: usize,
        /// Observed exclusive range end.
        end: usize,
        /// Required sequence length.
        length: usize,
    },
    /// A local alignment range lies outside its source sequence.
    LocalRangeOutOfBounds {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Observed range start.
        start: usize,
        /// Observed exclusive range end.
        end: usize,
        /// Source sequence length.
        length: usize,
    },
    /// A local operation path did not consume its declared source range exactly.
    LocalRangeConsumptionMismatch {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Declared exclusive range end.
        expected_end: usize,
        /// Exclusive position reached by replay.
        actual_end: usize,
    },
    /// An end-gap alignment range lies outside its source sequence.
    EndGapRangeOutOfBounds {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Observed range start.
        start: usize,
        /// Observed exclusive range end.
        end: usize,
        /// Source sequence length.
        length: usize,
    },
    /// A prefix was omitted even though that sequence start is scored.
    OmittedScoredStart {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Observed nonzero range start.
        start: usize,
    },
    /// A suffix was omitted even though that sequence end is scored.
    OmittedScoredEnd {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Observed exclusive range end.
        end: usize,
        /// Required sequence length.
        length: usize,
    },
    /// Both selected ranges begin strictly inside their source sequences.
    EndGapStartNotOnBoundary {
        /// First selected range start.
        first_start: usize,
        /// Second selected range start.
        second_start: usize,
    },
    /// Both selected ranges end strictly before their source sequence ends.
    EndGapEndNotOnBoundary {
        /// First selected range end.
        first_end: usize,
        /// First source sequence length.
        first_length: usize,
        /// Second selected range end.
        second_end: usize,
        /// Second source sequence length.
        second_length: usize,
    },
    /// An end-gap operation path did not consume its declared source range exactly.
    EndGapRangeConsumptionMismatch {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Declared exclusive range end.
        expected_end: usize,
        /// Exclusive position reached by replay.
        actual_end: usize,
    },
    /// A local alignment reported a negative score even though zero is a candidate.
    NegativeLocalScore {
        /// Reported negative score.
        score: Score,
    },
    /// A zero-score local result did not use the canonical empty representation.
    NonCanonicalEmptyLocalAlignment,
    /// A positive local traceback crossed a non-positive intermediate score.
    NonPositiveLocalPrefix {
        /// Alignment column containing the non-positive prefix end.
        column: usize,
        /// Recomputed prefix score.
        score: Score,
    },
    /// An operation attempted to consume beyond a sequence end.
    SequenceExhausted {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Alignment column containing the error.
        column: usize,
    },
    /// `Match` or `Substitution` disagreed with the actual symbol equality.
    OperationKindMismatch {
        /// Alignment column containing the error.
        column: usize,
        /// Recorded operation.
        operation: AlignmentOperation,
    },
    /// Recomputed score exceeded the score representation.
    ScoreOverflow {
        /// Alignment column containing the error.
        column: usize,
    },
    /// Recomputed and reported scores differ.
    ScoreMismatch {
        /// Score stored in the result.
        reported: Score,
        /// Score obtained by independent replay.
        recomputed: Score,
    },
}

impl fmt::Display for VerificationError {
    #[allow(clippy::too_many_lines)]
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
            Self::LocalRangeOutOfBounds {
                sequence,
                start,
                end,
                length,
            } => write!(
                formatter,
                "local {sequence} range {start}..{end} exceeds sequence length {length}"
            ),
            Self::LocalRangeConsumptionMismatch {
                sequence,
                expected_end,
                actual_end,
            } => write!(
                formatter,
                "local {sequence} range ends at {expected_end}, but replay ended at {actual_end}"
            ),
            Self::EndGapRangeOutOfBounds {
                sequence,
                start,
                end,
                length,
            } => write!(
                formatter,
                "end-gap {sequence} range {start}..{end} exceeds sequence length {length}"
            ),
            Self::OmittedScoredStart { sequence, start } => write!(
                formatter,
                "end-gap {sequence} range starts at {start}, but that sequence start is scored"
            ),
            Self::OmittedScoredEnd {
                sequence,
                end,
                length,
            } => write!(
                formatter,
                "end-gap {sequence} range ends at {end}, but that scored end requires {length}"
            ),
            Self::EndGapStartNotOnBoundary {
                first_start,
                second_start,
            } => write!(
                formatter,
                "end-gap ranges start at interior coordinates ({first_start}, {second_start})"
            ),
            Self::EndGapEndNotOnBoundary {
                first_end,
                first_length,
                second_end,
                second_length,
            } => write!(
                formatter,
                "end-gap ranges end at interior coordinates ({first_end}/{first_length}, \
                 {second_end}/{second_length})"
            ),
            Self::EndGapRangeConsumptionMismatch {
                sequence,
                expected_end,
                actual_end,
            } => write!(
                formatter,
                "end-gap {sequence} range ends at {expected_end}, but replay ended at {actual_end}"
            ),
            Self::NegativeLocalScore { score } => write!(
                formatter,
                "local alignment score {score} is negative even though the empty alignment scores zero"
            ),
            Self::NonCanonicalEmptyLocalAlignment => write!(
                formatter,
                "a zero-score local alignment must have no operations and canonical ranges 0..0"
            ),
            Self::NonPositiveLocalPrefix { column, score } => write!(
                formatter,
                "local alignment prefix ending at column {column} has non-positive score {score}"
            ),
            Self::SequenceExhausted { sequence, column } => write!(
                formatter,
                "alignment column {column} reads beyond the {sequence} sequence"
            ),
            Self::OperationKindMismatch { column, operation } => write!(
                formatter,
                "alignment column {column} is recorded as {operation:?}, inconsistent with its symbols"
            ),
            Self::ScoreOverflow { column } => {
                write!(
                    formatter,
                    "score overflow while verifying alignment column {column}"
                )
            }
            Self::ScoreMismatch {
                reported,
                recomputed,
            } => write!(
                formatter,
                "reported score {reported} differs from recomputed score {recomputed}"
            ),
        }
    }
}

impl std::error::Error for VerificationError {}
