use core::fmt;

use align_msa::MsaError;

/// Validation, numeric, or structural failure in the MUSCLE 5 pipeline.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum Muscle5Error {
    /// Shared MSA input or result validation failed.
    Msa(MsaError),
    /// A parameter has a non-finite or out-of-range value.
    InvalidParameter {
        /// Public parameter name.
        parameter: &'static str,
        /// Human-readable constraint.
        requirement: &'static str,
    },
    /// Checked integer arithmetic overflowed.
    ArithmeticOverflow {
        /// Pipeline operation in which overflow occurred.
        stage: &'static str,
    },
    /// A derived floating value was NaN or infinite.
    NonFiniteComputation {
        /// Pipeline operation producing the invalid value.
        stage: &'static str,
    },
    /// A supposedly impossible internal state was observed.
    InternalInvariant {
        /// Concise invariant description.
        message: &'static str,
    },
}

impl fmt::Display for Muscle5Error {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Msa(error) => write!(formatter, "MUSCLE 5 MSA validation failed: {error}"),
            Self::InvalidParameter {
                parameter,
                requirement,
            } => write!(
                formatter,
                "MUSCLE 5 parameter {parameter:?} must {requirement}"
            ),
            Self::ArithmeticOverflow { stage } => {
                write!(formatter, "integer overflow during MUSCLE 5 {stage}")
            }
            Self::NonFiniteComputation { stage } => {
                write!(formatter, "non-finite value during MUSCLE 5 {stage}")
            }
            Self::InternalInvariant { message } => {
                write!(formatter, "MUSCLE 5 internal invariant failed: {message}")
            }
        }
    }
}

impl std::error::Error for Muscle5Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Msa(error) => Some(error),
            _ => None,
        }
    }
}

impl From<MsaError> for Muscle5Error {
    fn from(error: MsaError) -> Self {
        Self::Msa(error)
    }
}
