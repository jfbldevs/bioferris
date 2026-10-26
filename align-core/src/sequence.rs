use core::fmt;
use core::str::FromStr;
use std::sync::OnceLock;

use crate::{Alphabet, Dna, Protein, Rna, SequenceError, Symbol};

/// An owned sequence validated against alphabet `A`.
#[derive(Clone)]
pub struct Sequence<A: Alphabet> {
    symbols: Vec<Symbol<A>>,
    indices: OnceLock<Box<[u8]>>,
}

impl<A: Alphabet> PartialEq for Sequence<A> {
    fn eq(&self, other: &Self) -> bool {
        self.symbols == other.symbols
    }
}

impl<A: Alphabet> Eq for Sequence<A> {}

impl<A: Alphabet> Sequence<A> {
    /// Parses a sequence from ASCII bytes.
    ///
    /// Lowercase symbols are normalized to uppercase. Whitespace, gaps, and FASTA
    /// control text are rejected rather than silently removed.
    ///
    /// # Errors
    ///
    /// Returns [`SequenceError`] at the first byte that is not part of alphabet `A`.
    pub fn from_ascii(input: impl AsRef<[u8]>) -> Result<Self, SequenceError> {
        let input = input.as_ref();
        // One exact-size pass; an invalid byte gets a placeholder that is never
        // returned, and the first one becomes the error.
        let mut invalid = None;
        let symbols = input
            .iter()
            .copied()
            .enumerate()
            .map(|(index, byte)| {
                Symbol::<A>::from_ascii(byte).unwrap_or_else(|| {
                    invalid.get_or_insert((index, byte));
                    Symbol::placeholder()
                })
            })
            .collect::<Vec<_>>();
        if let Some((index, byte)) = invalid {
            return Err(SequenceError::new(A::NAME, index, byte));
        }
        Ok(Self {
            symbols,
            indices: OnceLock::new(),
        })
    }

    /// Returns the number of biological symbols.
    #[must_use]
    pub fn len(&self) -> usize {
        self.symbols.len()
    }

    /// Returns whether the sequence has no symbols.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.symbols.is_empty()
    }

    /// Returns all validated symbols.
    #[must_use]
    pub fn as_slice(&self) -> &[Symbol<A>] {
        &self.symbols
    }

    /// Returns compact stable alphabet indices for byte-oriented kernels.
    ///
    /// The view is initialized lazily and then retained by the immutable sequence.
    #[must_use]
    pub fn as_indices(&self) -> &[u8] {
        self.indices.get_or_init(|| {
            self.symbols
                .iter()
                .map(|symbol| symbol.compact_index())
                .collect::<Vec<_>>()
                .into_boxed_slice()
        })
    }

    /// Returns the symbol at `index`, if present.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<Symbol<A>> {
        self.symbols.get(index).copied()
    }

    /// Iterates over validated symbols by value.
    #[must_use]
    pub fn iter(&self) -> impl ExactSizeIterator<Item = Symbol<A>> + '_ {
        self.symbols.iter().copied()
    }

    /// Returns the normalized sequence as an ASCII string.
    #[must_use]
    pub fn to_ascii_string(&self) -> String {
        self.symbols
            .iter()
            .map(|symbol| char::from(symbol.as_ascii()))
            .collect()
    }
}

impl<A: Alphabet> fmt::Debug for Sequence<A> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_tuple(A::NAME)
            .field(&self.to_ascii_string())
            .finish()
    }
}

impl<A: Alphabet> fmt::Display for Sequence<A> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        for symbol in &self.symbols {
            write!(formatter, "{symbol}")?;
        }
        Ok(())
    }
}

impl<A: Alphabet> FromStr for Sequence<A> {
    type Err = SequenceError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        Self::from_ascii(input)
    }
}

impl<A: Alphabet> TryFrom<&str> for Sequence<A> {
    type Error = SequenceError;

    fn try_from(input: &str) -> Result<Self, Self::Error> {
        Self::from_ascii(input)
    }
}

impl<A: Alphabet> TryFrom<String> for Sequence<A> {
    type Error = SequenceError;

    fn try_from(input: String) -> Result<Self, Self::Error> {
        Self::from_ascii(input)
    }
}

/// A validated DNA sequence.
pub type DnaSequence = Sequence<Dna>;

/// A validated RNA sequence.
pub type RnaSequence = Sequence<Rna>;

/// A validated protein sequence.
pub type ProteinSequence = Sequence<Protein>;

#[cfg(test)]
mod tests {
    use super::{DnaSequence, ProteinSequence, RnaSequence};

    #[test]
    fn empty_sequence_is_valid() {
        let sequence = DnaSequence::from_ascii(b"").expect("empty DNA is valid");
        assert!(sequence.is_empty());
        assert_eq!(sequence.to_string(), "");
    }

    #[test]
    fn sequence_normalizes_lowercase() {
        let sequence = DnaSequence::from_ascii(b"acgtnry").expect("IUB DNA is valid");
        assert_eq!(sequence.to_ascii_string(), "ACGTNRY");
    }

    #[test]
    fn standard_string_conversions_validate_and_normalize() {
        let borrowed = DnaSequence::try_from("acgt").expect("valid borrowed DNA");
        let owned = DnaSequence::try_from(String::from("tgca")).expect("valid owned DNA");
        assert_eq!(borrowed.to_string(), "ACGT");
        assert_eq!(owned.to_string(), "TGCA");
    }

    #[test]
    fn invalid_symbol_reports_exact_byte_offset() {
        let error = RnaSequence::from_ascii(b"AC T").expect_err("whitespace is invalid");
        assert_eq!(error.alphabet(), "RNA");
        assert_eq!(error.index(), 2);
        assert_eq!(error.byte(), b' ');
    }

    #[test]
    fn gap_is_not_a_raw_protein_symbol() {
        let error = ProteinSequence::from_ascii(b"AC-D").expect_err("gap must be rejected");
        assert_eq!(error.index(), 2);
        assert_eq!(error.byte(), b'-');
    }

    #[test]
    fn compact_indices_are_stable_and_do_not_affect_sequence_identity() {
        let sequence = DnaSequence::from_ascii(b"ACGTNRY").expect("IUB DNA is valid");
        let clone_before_initialization = sequence.clone();
        let first = sequence.as_indices();
        assert_eq!(first, &[0, 1, 2, 3, 14, 4, 5]);
        assert!(core::ptr::eq(first, sequence.as_indices()));
        assert_eq!(sequence, clone_before_initialization);

        let clone_after_initialization = sequence.clone();
        assert_eq!(clone_after_initialization.as_indices(), first);
        assert_eq!(sequence, clone_after_initialization);
    }
}
