use core::fmt;
use std::collections::HashSet;

use align_core::{Alphabet, Symbol};

use crate::MsaError;
use crate::record::validate_identifier;

/// One named row of a multiple alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignedRow<A: Alphabet> {
    identifier: String,
    symbols: Box<[Option<Symbol<A>>]>,
}

impl<A: Alphabet> AlignedRow<A> {
    /// Creates an aligned row from typed symbols and explicit gaps.
    ///
    /// # Errors
    ///
    /// Returns [`MsaError`] when the identifier is invalid.
    pub fn new(
        identifier: impl Into<String>,
        symbols: Vec<Option<Symbol<A>>>,
    ) -> Result<Self, MsaError> {
        let identifier = identifier.into();
        validate_identifier(&identifier)?;
        Ok(Self {
            identifier,
            symbols: symbols.into_boxed_slice(),
        })
    }

    /// Returns the row identifier.
    #[must_use]
    pub fn identifier(&self) -> &str {
        &self.identifier
    }

    /// Returns the aligned cells; `None` denotes a gap.
    #[must_use]
    pub const fn symbols(&self) -> &[Option<Symbol<A>>] {
        &self.symbols
    }

    /// Returns the alignment width of this row.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.symbols.len()
    }

    /// Returns whether this row has zero alignment columns.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.symbols.is_empty()
    }

    /// Renders the row with `-` for gaps.
    #[must_use]
    pub fn to_ascii_with_gaps(&self) -> String {
        self.symbols
            .iter()
            .map(|cell| cell.map_or('-', |symbol| char::from(symbol.as_ascii())))
            .collect()
    }

    /// Renders only biological symbols, omitting gaps.
    #[must_use]
    pub fn to_ungapped_ascii(&self) -> String {
        self.symbols
            .iter()
            .filter_map(|cell| *cell)
            .map(|symbol| char::from(symbol.as_ascii()))
            .collect()
    }

    /// Separates the identifier and aligned cells.
    #[must_use]
    pub fn into_parts(self) -> (String, Vec<Option<Symbol<A>>>) {
        (self.identifier, self.symbols.into_vec())
    }
}

/// A validated rectangular multiple-sequence alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MultipleAlignment<A: Alphabet> {
    rows: Box<[AlignedRow<A>]>,
    width: usize,
}

impl<A: Alphabet> MultipleAlignment<A> {
    /// Validates and constructs a multiple alignment.
    ///
    /// # Errors
    ///
    /// Returns [`MsaError`] for an empty collection, duplicate identifiers,
    /// inconsistent row widths, or an all-gap column.
    pub fn try_new(rows: Vec<AlignedRow<A>>) -> Result<Self, MsaError> {
        let Some(first) = rows.first() else {
            return Err(MsaError::EmptyCollection);
        };
        let width = first.len();
        let mut identifiers = HashSet::with_capacity(rows.len());
        for (index, row) in rows.iter().enumerate() {
            if row.len() != width {
                return Err(MsaError::RowWidthMismatch {
                    row: index,
                    expected: width,
                    actual: row.len(),
                });
            }
            if !identifiers.insert(row.identifier()) {
                return Err(MsaError::DuplicateIdentifier {
                    identifier: row.identifier().to_owned(),
                });
            }
        }
        for column in 0..width {
            if rows.iter().all(|row| row.symbols()[column].is_none()) {
                return Err(MsaError::AllGapColumn { column });
            }
        }
        Ok(Self {
            rows: rows.into_boxed_slice(),
            width,
        })
    }

    /// Returns all rows in their stable input order.
    #[must_use]
    pub const fn rows(&self) -> &[AlignedRow<A>] {
        &self.rows
    }

    /// Returns one row by zero-based index.
    #[must_use]
    pub fn row(&self, index: usize) -> Option<&AlignedRow<A>> {
        self.rows.get(index)
    }

    /// Returns one row by identifier.
    #[must_use]
    pub fn row_by_identifier(&self, identifier: &str) -> Option<&AlignedRow<A>> {
        self.rows.iter().find(|row| row.identifier() == identifier)
    }

    /// Returns the number of sequence rows.
    #[must_use]
    pub const fn sequence_count(&self) -> usize {
        self.rows.len()
    }

    /// Returns the shared alignment width.
    #[must_use]
    pub const fn width(&self) -> usize {
        self.width
    }

    /// Iterates over one alignment column, if it exists.
    #[must_use]
    pub fn column(
        &self,
        index: usize,
    ) -> Option<impl ExactSizeIterator<Item = Option<Symbol<A>>> + '_> {
        (index < self.width).then(|| self.rows.iter().map(move |row| row.symbols()[index]))
    }

    /// Consumes the alignment and returns its rows.
    #[must_use]
    pub fn into_rows(self) -> Vec<AlignedRow<A>> {
        self.rows.into_vec()
    }
}

impl<A: Alphabet> fmt::Display for MultipleAlignment<A> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let identifier_width = self
            .rows
            .iter()
            .map(|row| row.identifier().chars().count())
            .max()
            .unwrap_or(0);
        for (index, row) in self.rows.iter().enumerate() {
            if index != 0 {
                formatter.write_str("\n")?;
            }
            write!(
                formatter,
                "{:<identifier_width$}  {}",
                row.identifier(),
                row.to_ascii_with_gaps()
            )?;
        }
        Ok(())
    }
}
