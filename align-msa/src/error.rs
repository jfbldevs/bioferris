use core::fmt;

/// Structural or numeric failure in shared MSA data.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum MsaError {
    /// A record identifier is empty or only whitespace.
    EmptyIdentifier,
    /// A record identifier contains an ASCII control character.
    ControlCharacterInIdentifier {
        /// Invalid identifier supplied by the caller.
        identifier: String,
        /// Zero-based byte position of the control character.
        index: usize,
    },
    /// Two rows use the same identifier.
    DuplicateIdentifier {
        /// Repeated identifier.
        identifier: String,
    },
    /// An alignment or comparison was constructed without rows.
    EmptyCollection,
    /// A row has a different width from the first row.
    RowWidthMismatch {
        /// Zero-based row index.
        row: usize,
        /// Required width.
        expected: usize,
        /// Supplied width.
        actual: usize,
    },
    /// An alignment column consists entirely of gaps.
    AllGapColumn {
        /// Zero-based column index.
        column: usize,
    },
    /// Profile weights do not correspond one-for-one with rows.
    WeightCountMismatch {
        /// Number of alignment rows.
        rows: usize,
        /// Number of supplied weights.
        weights: usize,
    },
    /// A profile contains a zero sequence weight.
    ZeroWeight {
        /// Zero-based weight index.
        index: usize,
    },
    /// A flat distance matrix has the wrong number of entries.
    DistanceMatrixSizeMismatch {
        /// Matrix dimension.
        size: usize,
        /// Required flat length.
        expected: usize,
        /// Supplied flat length.
        actual: usize,
    },
    /// A distance is NaN, infinite, negative, or greater than one.
    InvalidDistance {
        /// Row index.
        row: usize,
        /// Column index.
        column: usize,
        /// Invalid value.
        value: f64,
    },
    /// A distance-matrix diagonal is not zero.
    NonZeroDistanceDiagonal {
        /// Diagonal index.
        index: usize,
        /// Invalid value.
        value: f64,
    },
    /// Symmetric positions in a distance matrix disagree.
    AsymmetricDistance {
        /// First index.
        first: usize,
        /// Second index.
        second: usize,
        /// Forward value.
        forward: f64,
        /// Reverse value.
        reverse: f64,
    },
    /// A guide-tree edge length is not finite.
    NonFiniteBranchLength {
        /// Invalid length.
        length: f64,
    },
    /// A guide-tree leaf index exceeds the declared leaf count.
    LeafIndexOutOfRange {
        /// Invalid leaf index.
        index: usize,
        /// Declared leaf count.
        leaf_count: usize,
    },
    /// A guide-tree leaf occurs more than once.
    DuplicateLeaf {
        /// Repeated leaf index.
        index: usize,
    },
    /// A declared guide-tree leaf does not occur in the tree.
    MissingLeaf {
        /// Missing leaf index.
        index: usize,
    },
    /// Candidate and reference alignments contain different identifiers.
    SequenceSetMismatch,
    /// Candidate and reference rows with the same identifier consume different
    /// ungapped sequences.
    UngappedSequenceMismatch {
        /// Identifier whose biological sequence differs.
        identifier: String,
    },
    /// A reference-column selection mask does not match the reference width.
    ReferenceColumnMaskLengthMismatch {
        /// Reference alignment width.
        expected: usize,
        /// Supplied mask length.
        actual: usize,
    },
}

impl fmt::Display for MsaError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyIdentifier => formatter.write_str("sequence identifier is empty"),
            Self::ControlCharacterInIdentifier { identifier, index } => write!(
                formatter,
                "sequence identifier {identifier:?} contains a control character at byte {index}"
            ),
            Self::DuplicateIdentifier { identifier } => {
                write!(formatter, "duplicate sequence identifier {identifier:?}")
            }
            Self::EmptyCollection => formatter.write_str("at least one sequence is required"),
            Self::RowWidthMismatch {
                row,
                expected,
                actual,
            } => write!(
                formatter,
                "alignment row {row} has width {actual}, expected {expected}"
            ),
            Self::AllGapColumn { column } => {
                write!(formatter, "alignment column {column} contains only gaps")
            }
            Self::WeightCountMismatch { rows, weights } => write!(
                formatter,
                "profile has {rows} rows but {weights} sequence weights"
            ),
            Self::ZeroWeight { index } => {
                write!(formatter, "profile weight {index} is zero")
            }
            Self::DistanceMatrixSizeMismatch {
                size,
                expected,
                actual,
            } => write!(
                formatter,
                "{size} by {size} distance matrix requires {expected} entries, got {actual}"
            ),
            Self::InvalidDistance { row, column, value } => write!(
                formatter,
                "distance ({row},{column}) must be finite and in [0,1], got {value}"
            ),
            Self::NonZeroDistanceDiagonal { index, value } => write!(
                formatter,
                "distance diagonal ({index},{index}) must be zero, got {value}"
            ),
            Self::AsymmetricDistance {
                first,
                second,
                forward,
                reverse,
            } => write!(
                formatter,
                "distance matrix is asymmetric at ({first},{second}): {forward} versus {reverse}"
            ),
            Self::NonFiniteBranchLength { length } => {
                write!(
                    formatter,
                    "guide-tree branch length is not finite: {length}"
                )
            }
            Self::LeafIndexOutOfRange { index, leaf_count } => write!(
                formatter,
                "guide-tree leaf {index} is outside declared range 0..{leaf_count}"
            ),
            Self::DuplicateLeaf { index } => {
                write!(formatter, "guide-tree leaf {index} occurs more than once")
            }
            Self::MissingLeaf { index } => {
                write!(formatter, "guide-tree leaf {index} is missing")
            }
            Self::SequenceSetMismatch => formatter.write_str(
                "candidate and reference alignments contain different sequence identifiers",
            ),
            Self::UngappedSequenceMismatch { identifier } => write!(
                formatter,
                "candidate and reference consume different symbols for {identifier:?}"
            ),
            Self::ReferenceColumnMaskLengthMismatch { expected, actual } => write!(
                formatter,
                "reference column mask has length {actual}, expected {expected}"
            ),
        }
    }
}

impl std::error::Error for MsaError {}
