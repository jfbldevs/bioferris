use core::fmt;
use std::io;

/// Failed alignment or graph serialization.
#[derive(Debug)]
#[non_exhaustive]
pub enum FormatError {
    /// Opening, writing, or flushing the destination failed.
    Io(io::Error),
    /// A pairwise traceback could not be rendered against the supplied inputs.
    Alignment(align_core::AlignmentStructureError),
    /// A caller-supplied output identifier is empty or contains a control byte.
    InvalidIdentifier {
        /// Rejected identifier.
        identifier: String,
    },
    /// An identifier cannot be represented losslessly in a selected format.
    InvalidFormatIdentifier {
        /// Output format requesting the restriction.
        format: &'static str,
        /// Rejected identifier.
        identifier: String,
        /// Stable explanation of the restriction.
        reason: &'static str,
    },
    /// A wrapping or block width was zero.
    InvalidLineWidth {
        /// Output format requesting a positive width.
        format: &'static str,
        /// Rejected width.
        value: usize,
    },
}

impl fmt::Display for FormatError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "I/O error: {error}"),
            Self::Alignment(error) => write!(formatter, "invalid pairwise alignment: {error}"),
            Self::InvalidIdentifier { identifier } => {
                write!(formatter, "invalid output identifier {identifier:?}")
            }
            Self::InvalidFormatIdentifier {
                format,
                identifier,
                reason,
            } => write!(
                formatter,
                "identifier {identifier:?} is invalid for {format}: {reason}"
            ),
            Self::InvalidLineWidth { format, value } => {
                write!(
                    formatter,
                    "{format} line width must be positive, got {value}"
                )
            }
        }
    }
}

impl std::error::Error for FormatError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::Alignment(error) => Some(error),
            Self::InvalidIdentifier { .. }
            | Self::InvalidFormatIdentifier { .. }
            | Self::InvalidLineWidth { .. } => None,
        }
    }
}

impl From<io::Error> for FormatError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<align_core::AlignmentStructureError> for FormatError {
    fn from(error: align_core::AlignmentStructureError) -> Self {
        Self::Alignment(error)
    }
}

pub(crate) fn validate_output_identifier(identifier: &str) -> Result<(), FormatError> {
    if identifier.trim().is_empty() || identifier.bytes().any(|byte| byte.is_ascii_control()) {
        return Err(FormatError::InvalidIdentifier {
            identifier: identifier.to_owned(),
        });
    }
    Ok(())
}
