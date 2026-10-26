//! Checked sequence ranges and literal searches over validated alphabets.

use std::fmt;
use std::iter::FusedIterator;
use std::ops::Range;

use align_core::{Alphabet, Sequence, SequenceError};

/// A checked sequence operation could not be completed.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum SequenceOperationError {
    /// The requested half-open range is reversed or exceeds the sequence.
    InvalidRange {
        /// Inclusive start of the requested range.
        start: usize,
        /// Exclusive end of the requested range.
        end: usize,
        /// Source sequence length.
        length: usize,
    },
    /// A pattern is not valid in the source alphabet.
    InvalidSequence(SequenceError),
    /// The requested combined length exceeds the platform's allocation range.
    LengthOverflow,
}

impl SequenceOperationError {
    /// Stable machine-readable category, independent of display wording.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::InvalidRange { .. } => "invalid_range",
            Self::InvalidSequence(_) => "invalid_sequence",
            Self::LengthOverflow => "length_overflow",
        }
    }
}

impl fmt::Display for SequenceOperationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRange { start, end, length } => write!(
                formatter,
                "range {start}..{end} is invalid for sequence length {length}"
            ),
            Self::InvalidSequence(error) => error.fmt(formatter),
            Self::LengthOverflow => {
                formatter.write_str("combined sequence length exceeds the platform limit")
            }
        }
    }
}

impl std::error::Error for SequenceOperationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidSequence(error) => Some(error),
            _ => None,
        }
    }
}

impl From<SequenceError> for SequenceOperationError {
    fn from(error: SequenceError) -> Self {
        Self::InvalidSequence(error)
    }
}

/// Whether a literal match may share symbols with the previous match.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MatchOverlap {
    /// Resume after the preceding match (the default count policy).
    #[default]
    NonOverlapping,
    /// Include matches starting inside the preceding match.
    Overlapping,
}

/// Checked copying operations and normalized literal pattern search.
///
/// Pattern letters are validated in `A` and normalized to uppercase. Ambiguous
/// letters match themselves, not their biological expansions. Empty patterns
/// match all `len + 1` boundaries. Positions are zero-based symbol offsets.
/// Slices are strict half-open ranges; invalid bounds are errors, not clipped.
///
/// ```
/// use bioferris::prelude::*;
/// let sequence = dna("ATGGCTTAA")?;
/// assert_eq!(sequence.slice(3..6)?.to_string(), "GCT");
/// assert_eq!(sequence.find("gct")?, Some(3));
/// assert_eq!(dna("AAAA")?.count_overlapping("AA")?, 3);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait SequenceExt<A: Alphabet>: Sized {
    /// Copies an in-bounds half-open range, including valid empty ranges.
    ///
    /// # Errors
    /// Returns [`SequenceOperationError::InvalidRange`] if `start > end` or
    /// `end > len`. Result validation errors are propagated defensively.
    fn slice(&self, range: Range<usize>) -> Result<Self, SequenceOperationError>;

    /// Concatenates two sequences of the same alphabet without modifying either.
    ///
    /// # Errors
    /// Returns [`SequenceOperationError::LengthOverflow`] if the result cannot be
    /// represented. Result validation errors are propagated defensively.
    fn concat(&self, other: &Self) -> Result<Self, SequenceOperationError>;

    /// Returns the first literal match, or `None`; an empty pattern matches zero.
    ///
    /// # Errors
    /// Returns [`SequenceOperationError::InvalidSequence`] for an invalid pattern.
    fn find(&self, pattern: impl AsRef<[u8]>) -> Result<Option<usize>, SequenceOperationError> {
        Ok(self.find_all(pattern, MatchOverlap::NonOverlapping)?.next())
    }

    /// Reports whether the normalized literal pattern occurs in this sequence.
    ///
    /// # Errors
    /// Returns [`SequenceOperationError::InvalidSequence`] for an invalid pattern.
    fn contains(&self, pattern: impl AsRef<[u8]>) -> Result<bool, SequenceOperationError> {
        Ok(self.find(pattern)?.is_some())
    }

    /// Counts non-overlapping literal matches, including all empty boundaries.
    ///
    /// # Errors
    /// Returns [`SequenceOperationError::InvalidSequence`] for an invalid pattern.
    fn count(&self, pattern: impl AsRef<[u8]>) -> Result<usize, SequenceOperationError> {
        Ok(self
            .find_all(pattern, MatchOverlap::NonOverlapping)?
            .count())
    }

    /// Counts all literal matches, including overlapping occurrences.
    ///
    /// # Errors
    /// Returns [`SequenceOperationError::InvalidSequence`] for an invalid pattern.
    fn count_overlapping(
        &self,
        pattern: impl AsRef<[u8]>,
    ) -> Result<usize, SequenceOperationError> {
        Ok(self.find_all(pattern, MatchOverlap::Overlapping)?.count())
    }

    /// Lazily iterates match starts in ascending order, with explicit overlap.
    ///
    /// A prefix-table search takes O(n + m) time for a full traversal and O(m)
    /// auxiliary memory. The iterator borrows the source and owns its normalized
    /// pattern; collecting every match is an additional caller-owned allocation.
    ///
    /// # Errors
    /// Returns [`SequenceOperationError::InvalidSequence`] for an invalid pattern.
    fn find_all(
        &self,
        pattern: impl AsRef<[u8]>,
        overlap: MatchOverlap,
    ) -> Result<SequenceMatches<'_, A>, SequenceOperationError>;
}

impl<A: Alphabet> SequenceExt<A> for Sequence<A> {
    fn slice(&self, range: Range<usize>) -> Result<Self, SequenceOperationError> {
        let selected =
            self.as_slice()
                .get(range.clone())
                .ok_or(SequenceOperationError::InvalidRange {
                    start: range.start,
                    end: range.end,
                    length: self.len(),
                })?;
        let bytes: Vec<_> = selected.iter().map(|symbol| symbol.as_ascii()).collect();
        Ok(Self::from_ascii(bytes)?)
    }

    fn concat(&self, other: &Self) -> Result<Self, SequenceOperationError> {
        let length = combined_length(self.len(), other.len())?;
        let mut bytes = Vec::with_capacity(length);
        bytes.extend(
            self.iter()
                .chain(other.iter())
                .map(align_core::Symbol::as_ascii),
        );
        Ok(Self::from_ascii(bytes)?)
    }

    fn find_all(
        &self,
        pattern: impl AsRef<[u8]>,
        overlap: MatchOverlap,
    ) -> Result<SequenceMatches<'_, A>, SequenceOperationError> {
        let pattern = Self::from_ascii(pattern)?;
        let mut prefix = vec![0; pattern.len()];
        let symbols = pattern.as_slice();
        let mut matched = 0;
        for index in 1..symbols.len() {
            while matched > 0 && symbols[index] != symbols[matched] {
                matched = prefix[matched - 1];
            }
            if symbols[index] == symbols[matched] {
                matched += 1;
            }
            prefix[index] = matched;
        }
        Ok(SequenceMatches {
            sequence: self,
            pattern,
            prefix,
            overlap,
            cursor: 0,
            matched: 0,
            exhausted: false,
        })
    }
}

fn combined_length(first: usize, second: usize) -> Result<usize, SequenceOperationError> {
    first
        .checked_add(second)
        .filter(|length| isize::try_from(*length).is_ok())
        .ok_or(SequenceOperationError::LengthOverflow)
}

/// A fused, lazy iterator of zero-based literal match starts.
///
/// Construct with [`SequenceExt::find_all`]. No match-position vector is stored.
#[derive(Debug)]
pub struct SequenceMatches<'a, A: Alphabet> {
    sequence: &'a Sequence<A>,
    pattern: Sequence<A>,
    prefix: Vec<usize>,
    overlap: MatchOverlap,
    cursor: usize,
    matched: usize,
    exhausted: bool,
}

impl<A: Alphabet> Iterator for SequenceMatches<'_, A> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            return None;
        }
        if self.pattern.is_empty() {
            let position = self.cursor;
            if position == self.sequence.len() {
                self.exhausted = true;
            } else {
                self.cursor += 1;
            }
            return Some(position);
        }
        let pattern = self.pattern.as_slice();
        while let Some(symbol) = self.sequence.get(self.cursor) {
            while self.matched > 0 && symbol != pattern[self.matched] {
                self.matched = self.prefix[self.matched - 1];
            }
            if symbol == pattern[self.matched] {
                self.matched += 1;
            }
            self.cursor += 1;
            if self.matched == pattern.len() {
                let position = self.cursor - pattern.len();
                self.matched = match self.overlap {
                    MatchOverlap::NonOverlapping => 0,
                    MatchOverlap::Overlapping => self.prefix[self.matched - 1],
                };
                return Some(position);
            }
        }
        self.exhausted = true;
        None
    }
}

impl<A: Alphabet> FusedIterator for SequenceMatches<'_, A> {}

#[cfg(test)]
mod tests {
    use super::{SequenceOperationError, combined_length};

    #[test]
    fn concatenation_checks_arithmetic_and_allocation_range() {
        assert_eq!(combined_length(0, 0), Ok(0));
        assert_eq!(
            combined_length(usize::MAX, 1),
            Err(SequenceOperationError::LengthOverflow)
        );
        assert_eq!(
            combined_length(isize::MAX as usize, 1),
            Err(SequenceOperationError::LengthOverflow)
        );
    }
}
