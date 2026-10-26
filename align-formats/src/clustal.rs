use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use align_core::Alphabet;
use align_msa::MultipleAlignment;

use crate::FormatError;

const CLUSTAL_WIDTH: usize = 60;

/// Writes a rectangular MSA in CLUSTAL block format.
///
/// # Errors
///
/// Returns [`FormatError::Io`] when the destination cannot be created or written.
pub fn write_clustal<A: Alphabet>(
    path: impl AsRef<Path>,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FormatError> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_clustal_to(&mut writer, alignment)?;
    writer.flush()?;
    Ok(())
}

/// Writes a rectangular MSA in CLUSTAL block format to a stream.
///
/// A `*` is emitted only for a literal all-row conserved column without gaps.
///
/// # Errors
///
/// Returns [`FormatError::Io`] when the destination cannot be written.
pub fn write_clustal_to<A: Alphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FormatError> {
    writer.write_all(b"CLUSTAL multiple sequence alignment\n\n")?;
    let identifier_width = alignment
        .rows()
        .iter()
        .map(|row| row.identifier().chars().count())
        .max()
        .unwrap_or(0);
    for start in (0..alignment.width()).step_by(CLUSTAL_WIDTH) {
        let end = start.saturating_add(CLUSTAL_WIDTH).min(alignment.width());
        for row in alignment.rows() {
            write!(writer, "{:<identifier_width$}  ", row.identifier())?;
            for cell in &row.symbols()[start..end] {
                let byte = cell.map_or(b'-', align_core::Symbol::as_ascii);
                writer.write_all(&[byte])?;
            }
            writer.write_all(b"\n")?;
        }
        write!(writer, "{:<identifier_width$}  ", "")?;
        for column in start..end {
            writer.write_all(if exactly_conserved(alignment, column) {
                b"*"
            } else {
                b" "
            })?;
        }
        writer.write_all(b"\n\n")?;
    }
    Ok(())
}

fn exactly_conserved<A: Alphabet>(alignment: &MultipleAlignment<A>, column: usize) -> bool {
    let mut selected = None;
    for row in alignment.rows() {
        let Some(symbol) = row.symbols()[column] else {
            return false;
        };
        if selected.is_some_and(|expected| symbol != expected) {
            return false;
        }
        selected = Some(symbol);
    }
    selected.is_some()
}

/// Convenient `.write_clustal(path)` output for an MSA.
pub trait ClustalWrite {
    /// Writes this alignment in CLUSTAL format.
    ///
    /// # Errors
    ///
    /// Returns [`FormatError::Io`] for failed output.
    fn write_clustal(&self, path: impl AsRef<Path>) -> Result<(), FormatError>;
}

impl<A: Alphabet> ClustalWrite for MultipleAlignment<A> {
    fn write_clustal(&self, path: impl AsRef<Path>) -> Result<(), FormatError> {
        write_clustal(path, self)
    }
}
