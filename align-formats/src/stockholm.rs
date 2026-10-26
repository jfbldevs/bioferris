use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use align_core::Alphabet;
use align_msa::MultipleAlignment;

use crate::FormatError;

/// Writes an unwrapped Stockholm 1.0 multiple alignment.
///
/// # Errors
///
/// Returns [`FormatError`] when an identifier contains whitespace or output fails.
pub fn write_stockholm<A: Alphabet>(
    path: impl AsRef<Path>,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FormatError> {
    validate_identifiers(alignment)?;
    let mut writer = BufWriter::new(File::create(path)?);
    write_validated(&mut writer, alignment)?;
    writer.flush()?;
    Ok(())
}

/// Streams an unwrapped Stockholm 1.0 multiple alignment.
///
/// # Errors
///
/// Returns [`FormatError`] when an identifier contains whitespace or output fails.
pub fn write_stockholm_to<A: Alphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FormatError> {
    validate_identifiers(alignment)?;
    write_validated(&mut writer, alignment)
}

fn write_validated<A: Alphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FormatError> {
    writer.write_all(b"# STOCKHOLM 1.0\n")?;
    let width = alignment
        .rows()
        .iter()
        .map(|row| row.identifier().chars().count())
        .max()
        .unwrap_or(0);
    for row in alignment.rows() {
        writeln!(
            writer,
            "{:<width$} {}",
            row.identifier(),
            row.to_ascii_with_gaps()
        )?;
    }
    writer.write_all(b"//\n")?;
    Ok(())
}

fn validate_identifiers<A: Alphabet>(alignment: &MultipleAlignment<A>) -> Result<(), FormatError> {
    for row in alignment.rows() {
        let identifier = row.identifier();
        if identifier.chars().any(char::is_whitespace) {
            return Err(FormatError::InvalidFormatIdentifier {
                format: "Stockholm 1.0",
                identifier: identifier.to_owned(),
                reason: "sequence names must not contain whitespace",
            });
        }
        if identifier.starts_with('#') || identifier == "//" {
            return Err(FormatError::InvalidFormatIdentifier {
                format: "Stockholm 1.0",
                identifier: identifier.to_owned(),
                reason: "sequence names must not collide with markup or the terminator",
            });
        }
    }
    Ok(())
}

/// Convenient Stockholm 1.0 output for a rectangular MSA.
pub trait StockholmWrite {
    /// Writes this alignment as unwrapped Stockholm 1.0.
    ///
    /// # Errors
    ///
    /// Returns [`FormatError`] for invalid identifiers or failed output.
    fn write_stockholm(&self, path: impl AsRef<Path>) -> Result<(), FormatError>;
}

impl<A: Alphabet> StockholmWrite for MultipleAlignment<A> {
    fn write_stockholm(&self, path: impl AsRef<Path>) -> Result<(), FormatError> {
        write_stockholm(path, self)
    }
}
