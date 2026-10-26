use core::fmt;
use core::fmt::Write as _;
use core::marker::PhantomData;

/// A finite biological alphabet whose external representation is uppercase ASCII.
pub trait Alphabet: Copy + Eq + fmt::Debug + Send + Sync + 'static {
    /// Human-readable alphabet name used in diagnostics.
    const NAME: &'static str;

    /// Valid symbols in their stable matrix-index order.
    const SYMBOLS: &'static [u8];

    /// Returns the stable zero-based index for an ASCII symbol.
    ///
    /// Input is normalized to uppercase. Whitespace and gap characters are not
    /// alphabet symbols and therefore return `None`.
    #[must_use]
    fn index_of(ascii: u8) -> Option<usize> {
        // One lookup in a table built at compile time from `SYMBOLS`.
        let index = const { symbol_indices(Self::SYMBOLS) }[usize::from(ascii)];
        (index != NO_SYMBOL).then_some(index)
    }
}

const NO_SYMBOL: usize = usize::MAX;

/// For every byte, the position of its uppercase form in `symbols`, as the
/// first match of a linear search, or [`NO_SYMBOL`].
#[allow(clippy::cast_possible_truncation)] // `byte` stays below 256.
const fn symbol_indices(symbols: &[u8]) -> [usize; 256] {
    let mut table = [NO_SYMBOL; 256];
    let mut byte = 0;
    while byte < table.len() {
        let normalized = (byte as u8).to_ascii_uppercase();
        let mut position = 0;
        while position < symbols.len() {
            if symbols[position] == normalized {
                table[byte] = position;
                break;
            }
            position += 1;
        }
        byte += 1;
    }
    table
}

/// A validated symbol from alphabet `A`.
#[repr(transparent)]
#[derive(Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Symbol<A: Alphabet> {
    index: u8,
    marker: PhantomData<A>,
}

impl<A: Alphabet> Symbol<A> {
    /// Parses and validates one ASCII symbol.
    #[must_use]
    pub fn from_ascii(ascii: u8) -> Option<Self> {
        let index = A::index_of(ascii)?;
        let index = u8::try_from(index).ok()?;
        Some(Self {
            index,
            marker: PhantomData,
        })
    }

    /// Returns the normalized uppercase ASCII representation.
    #[must_use]
    pub fn as_ascii(self) -> u8 {
        A::SYMBOLS[usize::from(self.index)]
    }

    /// Returns the stable index used by substitution matrices.
    #[must_use]
    pub fn alphabet_index(self) -> usize {
        usize::from(self.index)
    }

    pub(crate) const fn compact_index(self) -> u8 {
        self.index
    }

    /// Index 0, which may not be a valid symbol of `A`: only for slots that are
    /// discarded before they can be observed.
    pub(crate) const fn placeholder() -> Self {
        Self {
            index: 0,
            marker: PhantomData,
        }
    }
}

impl<A: Alphabet> fmt::Debug for Symbol<A> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_tuple("Symbol")
            .field(&char::from(self.as_ascii()))
            .finish()
    }
}

impl<A: Alphabet> fmt::Display for Symbol<A> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_char(char::from(self.as_ascii()))
    }
}

/// DNA using the IUB ambiguity alphabet.
///
/// Valid symbols are `A C G T R Y S W K M B D H V N`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Dna;

impl Alphabet for Dna {
    const NAME: &'static str = "DNA";
    const SYMBOLS: &'static [u8] = b"ACGTRYSWKMBDHVN";
}

/// RNA using the IUB ambiguity alphabet.
///
/// Valid symbols are `A C G U R Y S W K M B D H V N`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Rna;

impl Alphabet for Rna {
    const NAME: &'static str = "RNA";
    const SYMBOLS: &'static [u8] = b"ACGURYSWKMBDHVN";
}

/// Protein alphabet used by INSDC sequence records.
///
/// In addition to the 20 canonical amino acids, this alphabet accepts `B` (Asx),
/// `Z` (Glx), `X` (unknown), `J` (Leu/Ile), `U` (selenocysteine), and `O`
/// (pyrrolysine). A termination marker is not an amino-acid symbol and is rejected.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Protein;

impl Alphabet for Protein {
    const NAME: &'static str = "protein";
    const SYMBOLS: &'static [u8] = b"ARNDCQEGHILKMFPSTWYVBZXJUO";
}

#[cfg(test)]
mod tests {
    use super::{Alphabet, Dna, Protein, Rna, Symbol};

    /// The definition the lookup table must reproduce for every byte.
    fn linear_index_of<A: Alphabet>(ascii: u8) -> Option<usize> {
        let normalized = ascii.to_ascii_uppercase();
        A::SYMBOLS
            .iter()
            .position(|candidate| *candidate == normalized)
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct Unusual;

    impl Alphabet for Unusual {
        const NAME: &'static str = "unusual";
        // A repeated symbol keeps its first position; a lowercase one never
        // matches, because input is uppercased first.
        const SYMBOLS: &'static [u8] = b"QAxAQ*";
    }

    fn table_matches_linear_search<A: Alphabet>() {
        for byte in 0..=u8::MAX {
            assert_eq!(
                A::index_of(byte),
                linear_index_of::<A>(byte),
                "{} byte {byte}",
                A::NAME
            );
        }
    }

    #[test]
    fn lookup_table_equals_the_linear_search_for_every_byte() {
        table_matches_linear_search::<Dna>();
        table_matches_linear_search::<Rna>();
        table_matches_linear_search::<Protein>();
        table_matches_linear_search::<Unusual>();
        assert_eq!(Unusual::index_of(b'q'), Some(0));
        assert_eq!(Unusual::index_of(b'a'), Some(1));
        assert_eq!(Unusual::index_of(b'x'), None);
        assert_eq!(Unusual::index_of(b'*'), Some(5));
    }

    #[test]
    fn alphabets_have_no_duplicate_symbols() {
        for symbols in [Dna::SYMBOLS, Rna::SYMBOLS, Protein::SYMBOLS] {
            for (index, symbol) in symbols.iter().enumerate() {
                assert!(!symbols[..index].contains(symbol));
            }
        }
    }

    #[test]
    fn symbols_normalize_ascii_case() {
        let symbol = Symbol::<Dna>::from_ascii(b'g').expect("lowercase DNA symbol is valid");
        assert_eq!(symbol.as_ascii(), b'G');
        assert_eq!(symbol.alphabet_index(), 2);
    }

    #[test]
    fn molecule_specific_bases_are_enforced() {
        assert!(Symbol::<Dna>::from_ascii(b'T').is_some());
        assert!(Symbol::<Dna>::from_ascii(b'U').is_none());
        assert!(Symbol::<Rna>::from_ascii(b'U').is_some());
        assert!(Symbol::<Rna>::from_ascii(b'T').is_none());
    }

    #[test]
    fn extended_protein_symbols_are_supported() {
        for symbol in b"BZXJUO" {
            assert!(Symbol::<Protein>::from_ascii(*symbol).is_some());
        }
        assert!(Symbol::<Protein>::from_ascii(b'*').is_none());
        assert!(Symbol::<Protein>::from_ascii(b'-').is_none());
    }
}
