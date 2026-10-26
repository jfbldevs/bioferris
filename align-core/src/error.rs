use core::fmt;

/// Error returned when raw sequence data contains an invalid symbol.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SequenceError {
    alphabet: &'static str,
    index: usize,
    byte: u8,
}

impl SequenceError {
    pub(crate) const fn new(alphabet: &'static str, index: usize, byte: u8) -> Self {
        Self {
            alphabet,
            index,
            byte,
        }
    }

    /// Returns the alphabet being parsed.
    #[must_use]
    pub const fn alphabet(&self) -> &'static str {
        self.alphabet
    }

    /// Returns the zero-based byte offset of the invalid input.
    #[must_use]
    pub const fn index(&self) -> usize {
        self.index
    }

    /// Returns the rejected byte.
    #[must_use]
    pub const fn byte(&self) -> u8 {
        self.byte
    }
}

impl fmt::Display for SequenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.byte.is_ascii_graphic() {
            write!(
                formatter,
                "invalid {} symbol '{}' at byte {}",
                self.alphabet,
                char::from(self.byte),
                self.index
            )
        } else {
            write!(
                formatter,
                "invalid {} byte 0x{:02X} at byte {}",
                self.alphabet, self.byte, self.index
            )
        }
    }
}

impl std::error::Error for SequenceError {}

/// Invalid scoring configuration.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ScoringError {
    /// A linear gap reward would make the scoring model biologically misleading.
    PositiveGapScore {
        /// Rejected score.
        score: i64,
    },
    /// An affine gap opening score was positive.
    PositiveAffineGapOpenScore {
        /// Rejected score.
        score: i64,
    },
    /// An affine gap extension score was positive.
    PositiveAffineGapExtensionScore {
        /// Rejected score.
        score: i64,
    },
    /// An affine gap opening was less penalizing than an extension.
    AffineGapOpenLessPenalizingThanExtension {
        /// Rejected opening score.
        open_score: i64,
        /// Rejected extension score.
        extension_score: i64,
    },
    /// A dense substitution table did not contain one score per symbol pair.
    MatrixSizeMismatch {
        /// Alphabet name.
        alphabet: &'static str,
        /// Required number of elements.
        expected: usize,
        /// Supplied number of elements.
        actual: usize,
    },
}

impl fmt::Display for ScoringError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PositiveGapScore { score } => {
                write!(
                    formatter,
                    "linear gap score must be non-positive, got {score}"
                )
            }
            Self::PositiveAffineGapOpenScore { score } => write!(
                formatter,
                "affine gap opening score must be non-positive, got {score}"
            ),
            Self::PositiveAffineGapExtensionScore { score } => write!(
                formatter,
                "affine gap extension score must be non-positive, got {score}"
            ),
            Self::AffineGapOpenLessPenalizingThanExtension {
                open_score,
                extension_score,
            } => write!(
                formatter,
                "affine gap opening score {open_score} must not exceed extension score {extension_score}"
            ),
            Self::MatrixSizeMismatch {
                alphabet,
                expected,
                actual,
            } => write!(
                formatter,
                "{alphabet} substitution matrix requires {expected} scores, got {actual}"
            ),
        }
    }
}

impl std::error::Error for ScoringError {}

/// Structural inconsistency in an alignment result or rendering request.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum AlignmentStructureError {
    /// A declared range has its start after its end.
    InvalidRange {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Inclusive range start.
        start: usize,
        /// Exclusive range end.
        end: usize,
    },
    /// Operations do not consume exactly the declared source ranges.
    ConsumptionMismatch {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Number of symbols declared by the range.
        expected: usize,
        /// Number of symbols consumed by operations.
        actual: usize,
    },
    /// A declared range is outside the sequence supplied for rendering.
    RangeOutOfBounds {
        /// Sequence label, either `first` or `second`.
        sequence: &'static str,
        /// Exclusive end of the requested range.
        end: usize,
        /// Available sequence length.
        length: usize,
    },
    /// Memory for a rendered alignment could not be reserved.
    AllocationFailed {
        /// Requested number of alignment columns.
        columns: usize,
    },
}

impl fmt::Display for AlignmentStructureError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRange {
                sequence,
                start,
                end,
            } => write!(
                formatter,
                "{sequence} alignment range starts at {start}, after its end at {end}"
            ),
            Self::ConsumptionMismatch {
                sequence,
                expected,
                actual,
            } => write!(
                formatter,
                "alignment operations consume {actual} symbols from the {sequence} sequence, expected {expected}"
            ),
            Self::RangeOutOfBounds {
                sequence,
                end,
                length,
            } => write!(
                formatter,
                "{sequence} alignment range ends at {end}, but sequence length is {length}"
            ),
            Self::AllocationFailed { columns } => write!(
                formatter,
                "could not reserve memory for {columns} rendered alignment columns"
            ),
        }
    }
}

impl std::error::Error for AlignmentStructureError {}
