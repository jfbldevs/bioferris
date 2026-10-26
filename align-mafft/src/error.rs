use core::fmt;

use align_msa::MsaError;

/// Validation, numeric, or structural failure in a MAFFT-family pipeline.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum MafftError {
    /// Shared MSA input or result validation failed.
    Msa(MsaError),
    /// Checked integer arithmetic overflowed.
    ArithmeticOverflow {
        /// Pipeline operation in which overflow occurred.
        stage: &'static str,
    },
    /// A derived floating-point value was NaN or infinite.
    NonFiniteComputation {
        /// Pipeline operation producing the invalid value.
        stage: &'static str,
    },
    /// A deterministic worker count of zero was requested.
    InvalidThreadCount {
        /// Rejected worker count.
        value: usize,
    },
    /// A supposedly impossible internal state was observed.
    InternalInvariant {
        /// Concise invariant description.
        message: &'static str,
    },
}

impl fmt::Display for MafftError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Msa(error) => write!(formatter, "MAFFT MSA validation failed: {error}"),
            Self::ArithmeticOverflow { stage } => {
                write!(formatter, "integer overflow during MAFFT {stage}")
            }
            Self::NonFiniteComputation { stage } => {
                write!(formatter, "non-finite value during MAFFT {stage}")
            }
            Self::InvalidThreadCount { value } => {
                write!(
                    formatter,
                    "MAFFT thread count must be positive, got {value}"
                )
            }
            Self::InternalInvariant { message } => {
                write!(formatter, "MAFFT internal invariant failed: {message}")
            }
        }
    }
}

impl std::error::Error for MafftError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Msa(error) => Some(error),
            _ => None,
        }
    }
}

impl From<MsaError> for MafftError {
    fn from(error: MsaError) -> Self {
        Self::Msa(error)
    }
}
