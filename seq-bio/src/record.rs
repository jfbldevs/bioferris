use crate::io::{BioIoError, IoErrorKind};
use crate::{Sequence, SequenceExt, SequenceOperationError, prelude::Alphabet};
use std::ops::Range;

/// An immutable biological record with original header/case and optional Phred scores.
///
/// Computational sequences are normalized using the existing alphabet. Source
/// case is retained separately; line wrapping and line endings are not retained.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Record<A: Alphabet> {
    header: String,
    id: String,
    description: Option<String>,
    sequence: Sequence<A>,
    source_sequence: String,
    qualities: Option<Vec<u8>>,
}

/// Where a checked header's identifier and description lie in it.
#[derive(Clone, Debug)]
pub(crate) struct HeaderParts {
    id: Range<usize>,
    description: Option<Range<usize>>,
}

/// Applies the record header checks and locates the identifier and the trimmed
/// description, if any. Readers call it before reading a record's sequence, so
/// a bad header is reported first.
pub(crate) fn split_header(header: &str) -> Result<HeaderParts, BioIoError> {
    // For ASCII text, a control character is exactly an ASCII control byte.
    let control = if header.is_ascii() {
        header
            .bytes()
            .any(|byte| byte.is_ascii_control() && byte != b'\t')
    } else {
        header.chars().any(|c| c.is_control() && c != '\t')
    };
    if control {
        return Err(BioIoError::new(
            IoErrorKind::InvalidHeader,
            "record",
            "header contains control characters",
        ));
    }
    let trimmed = header.trim();
    let start = header.len() - header.trim_start().len();
    let end = trimmed.find(char::is_whitespace).unwrap_or(trimmed.len());
    if end == 0 {
        return Err(BioIoError::new(
            IoErrorKind::InvalidHeader,
            "record",
            "record identifier is empty",
        ));
    }
    let rest = &trimmed[end..];
    let description = rest.trim();
    let description = (!description.is_empty()).then(|| {
        let offset = start + end + (rest.len() - rest.trim_start().len());
        offset..offset + description.len()
    });
    Ok(HeaderParts {
        id: start..start + end,
        description,
    })
}

impl<A: Alphabet> Record<A> {
    /// Validates a record, preserving its header and original sequence case.
    ///
    /// # Errors
    /// Rejects an empty identity, header controls other than tab, invalid sequence
    /// symbols, quality/sequence length mismatch or Phred scores outside 0..=93.
    pub fn from_ascii(
        header: impl Into<String>,
        input: impl AsRef<[u8]>,
        qualities: Option<Vec<u8>>,
    ) -> Result<Self, BioIoError> {
        Self::from_ascii_owned(header.into(), input.as_ref().to_vec(), qualities)
    }

    /// [`Self::from_ascii`] for a sequence buffer the caller owns, such as a
    /// reader's assembled record: the buffer becomes the source text without a
    /// copy, after its growth slack is released.
    pub(crate) fn from_ascii_owned(
        header: String,
        input: Vec<u8>,
        qualities: Option<Vec<u8>>,
    ) -> Result<Self, BioIoError> {
        let parts = split_header(&header)?;
        Self::from_checked_header(header, &parts, input, qualities)
    }

    /// [`Self::from_ascii_owned`] for a header that [`split_header`] already
    /// checked and located.
    pub(crate) fn from_checked_header(
        header: String,
        parts: &HeaderParts,
        mut input: Vec<u8>,
        qualities: Option<Vec<u8>>,
    ) -> Result<Self, BioIoError> {
        let id = header[parts.id.clone()].to_owned();
        let description = parts
            .description
            .clone()
            .map(|range| header[range].to_owned());
        // Shrinking first keeps the peak near twice the record: the buffer and
        // the normalized sequence, both exactly sized.
        input.shrink_to_fit();
        let sequence = Sequence::from_ascii(&input).map_err(|error| {
            BioIoError::new(IoErrorKind::InvalidSequence, "record", error.to_string())
        })?;
        let source_sequence = String::from_utf8(input).map_err(|error| {
            BioIoError::new(IoErrorKind::InvalidSequence, "record", error.to_string())
        })?;
        if let Some(scores) = &qualities {
            if scores.len() != sequence.len() {
                return Err(BioIoError::new(
                    IoErrorKind::LengthMismatch,
                    "record",
                    "quality count differs from sequence length",
                ));
            }
            if scores.iter().any(|score| *score > 93) {
                return Err(BioIoError::new(
                    IoErrorKind::InvalidQuality,
                    "record",
                    "Phred+33 scores must be in 0..=93",
                ));
            }
        }
        Ok(Self {
            header,
            id,
            description,
            sequence,
            source_sequence,
            qualities,
        })
    }

    /// Returns the first whitespace-delimited header token.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }
    /// Returns the complete header without its format prefix or line ending.
    #[must_use]
    pub fn header(&self) -> &str {
        &self.header
    }
    /// Returns trimmed description text after the identifier, if present.
    #[must_use]
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }
    /// Returns the validated uppercase computational sequence.
    #[must_use]
    pub const fn sequence(&self) -> &Sequence<A> {
        &self.sequence
    }
    /// Returns validated source symbols, preserving letter case but not wrapping.
    #[must_use]
    pub fn source_sequence(&self) -> &str {
        &self.source_sequence
    }
    /// Returns decoded Phred scores, if supplied by the source.
    #[must_use]
    pub fn qualities(&self) -> Option<&[u8]> {
        self.qualities.as_deref()
    }

    /// Copies a strict half-open range with matching source case and qualities.
    ///
    /// Header, ID and description are retained verbatim. The original record is
    /// unchanged. This is an unannotated record operation: no feature locations
    /// are inferred or transformed. Empty records remain valid in memory even
    /// when a selected output format cannot serialize them.
    ///
    /// # Errors
    /// Returns [`SequenceOperationError::InvalidRange`] for reversed/out-of-bounds
    /// coordinates. Existing sequence validation errors are propagated defensively.
    pub fn slice(&self, range: Range<usize>) -> Result<Self, SequenceOperationError> {
        let sequence = self.sequence.slice(range.clone())?;
        // Construction guarantees ASCII source and quality/sequence length equality.
        Ok(Self {
            header: self.header.clone(),
            id: self.id.clone(),
            description: self.description.clone(),
            sequence,
            source_sequence: self.source_sequence[range.clone()].to_owned(),
            qualities: self.qualities.as_ref().map(|scores| scores[range].to_vec()),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::{Record, split_header};
    use crate::Dna;

    /// The header rules as written before readers shared `split_header`.
    fn reference_parts(header: &str) -> Result<(String, Option<String>), &'static str> {
        if header.chars().any(|c| c.is_control() && c != '\t') {
            return Err("header contains control characters");
        }
        let trimmed = header.trim();
        let end = trimmed.find(char::is_whitespace).unwrap_or(trimmed.len());
        let id = trimmed[..end].to_owned();
        if id.is_empty() {
            return Err("record identifier is empty");
        }
        let description =
            (!trimmed[end..].trim().is_empty()).then(|| trimmed[end..].trim().to_owned());
        Ok((id, description))
    }

    #[test]
    fn header_parts_follow_the_original_rules() {
        let pieces = [
            "", "a", "id", " ", "\t", "x y", "\u{a0}", "\u{85}", "\u{2003}", "\u{3000}", "\u{7f}",
            "\u{1}", "\u{9f}", "\u{200b}", "µ", "é",
        ];
        for first in pieces {
            for second in pieces {
                for third in pieces {
                    let header = format!("{first}{second}{third}");
                    let parts = split_header(&header)
                        .map(|parts| {
                            (
                                header[parts.id.clone()].to_owned(),
                                parts
                                    .description
                                    .clone()
                                    .map(|range| header[range].to_owned()),
                            )
                        })
                        .map_err(|error| error.to_string());
                    let expected = reference_parts(&header).map_err(|message| {
                        crate::io::BioIoError::new(
                            crate::io::IoErrorKind::InvalidHeader,
                            "record",
                            message,
                        )
                        .to_string()
                    });
                    assert_eq!(parts, expected, "{header:?}");
                }
            }
        }
    }

    /// Header, sequence bytes and optional Phred scores.
    type Case<'a> = (&'a str, &'a [u8], Option<Vec<u8>>);

    #[test]
    fn owned_construction_matches_borrowed_construction() {
        let cases: [Case; 7] = [
            ("seq1 description", b"ACGTacgtN", None),
            ("seq2", b"", None),
            ("seq3", b"ACGT", Some(vec![0, 10, 40, 93])),
            ("seq4", b"ACGX", None),
            ("", b"ACGT", None),
            ("bad\u{7}header", b"ACGT", None),
            ("seq5", b"ACGT", Some(vec![0, 94, 0, 0])),
        ];
        for (header, input, qualities) in cases {
            let mut buffer = Vec::with_capacity(input.len() * 2 + 8);
            buffer.extend_from_slice(input);
            let owned =
                Record::<Dna>::from_ascii_owned(header.to_owned(), buffer, qualities.clone());
            let borrowed = Record::<Dna>::from_ascii(header, input, qualities);
            match (owned, borrowed) {
                (Ok(owned), Ok(borrowed)) => assert_eq!(owned, borrowed),
                (Err(owned), Err(borrowed)) => {
                    assert_eq!(owned.kind(), borrowed.kind());
                    assert_eq!(owned.message(), borrowed.message());
                }
                (owned, borrowed) => panic!("{header:?}: {owned:?} versus {borrowed:?}"),
            }
        }
        let quality_mismatch =
            Record::<Dna>::from_ascii_owned("seq6".to_owned(), b"ACGT".to_vec(), Some(vec![1, 2]));
        assert!(quality_mismatch.is_err());
    }
}
