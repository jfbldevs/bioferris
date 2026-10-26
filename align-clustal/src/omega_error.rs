use core::fmt;

use align_msa::MsaError;

/// Failure while configuring or executing a Clustal Omega reconstruction.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum ClustalOmegaError {
    /// Shared MSA structure validation failed.
    Msa(MsaError),
    /// Parameters for one molecule class were supplied to another alphabet.
    MoleculeProfileMismatch {
        /// Profile required by the alphabet.
        expected: &'static str,
        /// Profile encoded by the supplied parameters.
        actual: &'static str,
    },
    /// A floating-point parameter is outside its documented domain.
    InvalidParameter {
        /// Stable parameter name.
        parameter: &'static str,
        /// Rejected value.
        value: f64,
        /// Human-readable accepted domain.
        requirement: &'static str,
    },
    /// An mBed leaf-cluster limit cannot support binary splitting.
    InvalidClusterSize {
        /// Rejected maximum cluster size.
        value: usize,
    },
    /// A deterministic worker pool cannot contain zero threads.
    InvalidThreadCount {
        /// Rejected worker count.
        value: usize,
    },
    /// Checked dimension, score, distance, or allocation arithmetic overflowed.
    ArithmeticOverflow {
        /// Pipeline stage that detected the overflow.
        stage: &'static str,
    },
    /// A probability calculation produced an invalid finite-domain value.
    InvalidProbability {
        /// Pipeline stage that detected the value.
        stage: &'static str,
        /// Rejected probability.
        value: f64,
    },
    /// A traceback, clustering, HMM, or tree invariant contradicted its producer.
    InternalInvariant {
        /// Stable diagnostic description.
        message: &'static str,
    },
}

impl fmt::Display for ClustalOmegaError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Msa(error) => write!(formatter, "invalid MSA structure: {error}"),
            Self::MoleculeProfileMismatch { expected, actual } => write!(
                formatter,
                "alphabet requires the {expected} Clustal Omega profile, but parameters use {actual}"
            ),
            Self::InvalidParameter {
                parameter,
                value,
                requirement,
            } => write!(
                formatter,
                "Clustal Omega parameter {parameter} must be {requirement}; got {value}"
            ),
            Self::InvalidClusterSize { value } => write!(
                formatter,
                "Clustal Omega mBed maximum cluster size must be at least 2; got {value}"
            ),
            Self::InvalidThreadCount { value } => write!(
                formatter,
                "Clustal Omega worker thread count must be at least 1; got {value}"
            ),
            Self::ArithmeticOverflow { stage } => {
                write!(formatter, "checked arithmetic overflow during {stage}")
            }
            Self::InvalidProbability { stage, value } => write!(
                formatter,
                "invalid probability {value} during Clustal Omega {stage}"
            ),
            Self::InternalInvariant { message } => write!(
                formatter,
                "Clustal Omega internal invariant failed: {message}"
            ),
        }
    }
}

impl std::error::Error for ClustalOmegaError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Msa(error) => Some(error),
            _ => None,
        }
    }
}

impl From<MsaError> for ClustalOmegaError {
    fn from(error: MsaError) -> Self {
        Self::Msa(error)
    }
}
