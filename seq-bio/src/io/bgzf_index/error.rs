use std::{fmt, io};

use super::{BgzfError, BioIoError, IndexError};

/// Preserves the precise BGZF framing or ordinary sequence-index failure.
#[derive(Debug)]
#[non_exhaustive]
pub enum BgzfIndexError {
    /// Filesystem, schema, record, resource or stability failure.
    Index(IndexError),
    /// Invalid framing, integrity, EOF or virtual position.
    Bgzf(BgzfError),
}

impl BgzfIndexError {
    /// Stable machine-readable error category.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            Self::Index(error) => error.code(),
            Self::Bgzf(error) => error.code(),
        }
    }
}

impl fmt::Display for BgzfIndexError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Index(error) => error.fmt(formatter),
            Self::Bgzf(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for BgzfIndexError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Index(error) => Some(error),
            Self::Bgzf(error) => Some(error),
        }
    }
}

impl From<IndexError> for BgzfIndexError {
    fn from(error: IndexError) -> Self {
        Self::Index(error)
    }
}
impl From<BgzfError> for BgzfIndexError {
    fn from(error: BgzfError) -> Self {
        Self::Bgzf(error)
    }
}
impl From<BioIoError> for BgzfIndexError {
    fn from(error: BioIoError) -> Self {
        Self::Index(error.into())
    }
}
impl From<io::Error> for BgzfIndexError {
    fn from(error: io::Error) -> Self {
        Self::Index(error.into())
    }
}
