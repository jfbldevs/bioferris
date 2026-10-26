//! Unified errors for the ergonomic facade.

use core::fmt;

/// A result returned by the ergonomic facade.
pub type Result<T> = core::result::Result<T, Error>;

/// A structured failure from sequence parsing, configuration, or execution.
#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    /// A batch execution requested zero worker threads.
    InvalidThreadCount {
        /// Invalid worker count supplied by the caller.
        value: usize,
    },
    /// A batch worker terminated by panicking.
    BatchWorkerPanicked,
    /// An execution option is not available for the selected scientific profile.
    UnsupportedOption {
        /// Selected algorithm or profile.
        algorithm: &'static str,
        /// Unsupported option name.
        option: &'static str,
    },
    /// A biological sequence contains an invalid symbol.
    Sequence(align_core::SequenceError),
    /// A score or gap configuration is invalid.
    Scoring(align_core::ScoringError),
    /// A canonical alignment result could not be constructed.
    AlignmentStructure(align_core::AlignmentStructureError),
    /// An adaptive affine score could not be negated into a representable penalty.
    AdaptivePenaltyOverflow {
        /// Dense affine score that could not be represented as a positive penalty.
        score: align_core::Score,
    },
    /// An exact dynamic-programming pairwise alignment failed.
    Pairwise(align_pairwise::AlignmentError),
    /// Myers preprocessing or execution failed.
    Myers(align_fast::MyersError),
    /// Wavefront execution failed.
    Wavefront(align_fast::WavefrontError),
    /// X-drop configuration or execution failed.
    XDrop(align_fast::XDropError),
    /// Shared multiple-alignment validation failed.
    MultipleAlignment(align_msa::MsaError),
    /// Clustal Omega execution failed.
    ClustalOmega(align_clustal::ClustalOmegaError),
    /// MUSCLE 5 execution failed.
    Muscle5(align_muscle::Muscle5Error),
    /// MAFFT execution failed.
    Mafft(align_mafft::MafftError),
    /// POASTA graph alignment failed.
    Poasta(align_poa::PoaError),
    /// FASTA input or output failed.
    Fasta(align_formats::FastaError),
    /// Alignment, phylogenetic matrix, graph, structured, or report output failed.
    Output(align_formats::FormatError),
    /// A human-readable report configuration is invalid.
    Report(crate::report::ReportError),
}

impl fmt::Display for Error {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidThreadCount { value } => {
                write!(formatter, "thread count must be positive, got {value}")
            }
            Self::BatchWorkerPanicked => formatter.write_str("an alignment batch worker panicked"),
            Self::UnsupportedOption { algorithm, option } => {
                write!(formatter, "{algorithm} does not support option {option:?}")
            }
            Self::Sequence(error) => write!(formatter, "invalid sequence: {error}"),
            Self::Scoring(error) => write!(formatter, "invalid scoring configuration: {error}"),
            Self::AlignmentStructure(error) => {
                write!(formatter, "invalid alignment structure: {error}")
            }
            Self::AdaptivePenaltyOverflow { score } => {
                write!(
                    formatter,
                    "adaptive affine score {score} cannot be negated into a penalty"
                )
            }
            Self::Pairwise(error) => write!(formatter, "pairwise alignment failed: {error}"),
            Self::Myers(error) => write!(formatter, "Myers alignment failed: {error}"),
            Self::Wavefront(error) => write!(formatter, "wavefront alignment failed: {error}"),
            Self::XDrop(error) => write!(formatter, "X-drop alignment failed: {error}"),
            Self::MultipleAlignment(error) => {
                write!(formatter, "multiple-alignment data is invalid: {error}")
            }
            Self::ClustalOmega(error) => write!(formatter, "Clustal Omega failed: {error}"),
            Self::Muscle5(error) => write!(formatter, "MUSCLE 5 failed: {error}"),
            Self::Mafft(error) => write!(formatter, "MAFFT failed: {error}"),
            Self::Poasta(error) => write!(formatter, "POASTA failed: {error}"),
            Self::Fasta(error) => write!(formatter, "FASTA operation failed: {error}"),
            Self::Output(error) => write!(formatter, "output operation failed: {error}"),
            Self::Report(error) => write!(formatter, "report configuration is invalid: {error}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidThreadCount { .. }
            | Self::BatchWorkerPanicked
            | Self::UnsupportedOption { .. }
            | Self::AdaptivePenaltyOverflow { .. } => None,
            Self::Sequence(error) => Some(error),
            Self::Scoring(error) => Some(error),
            Self::AlignmentStructure(error) => Some(error),
            Self::Pairwise(error) => Some(error),
            Self::Myers(error) => Some(error),
            Self::Wavefront(error) => Some(error),
            Self::XDrop(error) => Some(error),
            Self::MultipleAlignment(error) => Some(error),
            Self::ClustalOmega(error) => Some(error),
            Self::Muscle5(error) => Some(error),
            Self::Mafft(error) => Some(error),
            Self::Poasta(error) => Some(error),
            Self::Fasta(error) => Some(error),
            Self::Output(error) => Some(error),
            Self::Report(error) => Some(error),
        }
    }
}

macro_rules! impl_from {
    ($source:ty, $variant:ident) => {
        impl From<$source> for Error {
            fn from(error: $source) -> Self {
                Self::$variant(error)
            }
        }
    };
}

impl_from!(align_core::SequenceError, Sequence);
impl_from!(align_core::ScoringError, Scoring);
impl_from!(align_core::AlignmentStructureError, AlignmentStructure);
impl_from!(align_pairwise::AlignmentError, Pairwise);
impl_from!(align_fast::MyersError, Myers);
impl_from!(align_fast::WavefrontError, Wavefront);
impl_from!(align_fast::XDropError, XDrop);
impl_from!(align_msa::MsaError, MultipleAlignment);
impl_from!(align_clustal::ClustalOmegaError, ClustalOmega);
impl_from!(align_muscle::Muscle5Error, Muscle5);
impl_from!(align_mafft::MafftError, Mafft);
impl_from!(align_poa::PoaError, Poasta);
impl_from!(align_formats::FastaError, Fasta);
impl_from!(align_formats::FormatError, Output);
impl_from!(crate::report::ReportError, Report);
