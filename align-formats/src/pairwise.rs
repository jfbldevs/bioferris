use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use align_core::{Alphabet, PairwiseAlignment, Sequence};

use crate::error::{FormatError, validate_output_identifier};

const FASTA_WIDTH: usize = 80;

/// Writes a pairwise traceback as two aligned FASTA records.
///
/// # Errors
///
/// Returns [`FormatError`] for invalid identifiers, an input/result mismatch, or
/// failed I/O.
#[allow(clippy::too_many_arguments)]
pub fn write_pairwise_fasta<A: Alphabet>(
    path: impl AsRef<Path>,
    first_identifier: &str,
    first: &Sequence<A>,
    second_identifier: &str,
    second: &Sequence<A>,
    alignment: &PairwiseAlignment,
) -> Result<(), FormatError> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_pairwise_fasta_to(
        &mut writer,
        first_identifier,
        first,
        second_identifier,
        second,
        alignment,
    )?;
    writer.flush()?;
    Ok(())
}

/// Writes a pairwise traceback as two aligned FASTA records to a stream.
///
/// # Errors
///
/// Returns [`FormatError`] for invalid identifiers, an input/result mismatch, or
/// failed I/O.
#[allow(clippy::too_many_arguments)]
pub fn write_pairwise_fasta_to<A: Alphabet>(
    mut writer: impl Write,
    first_identifier: &str,
    first: &Sequence<A>,
    second_identifier: &str,
    second: &Sequence<A>,
    alignment: &PairwiseAlignment,
) -> Result<(), FormatError> {
    validate_output_identifier(first_identifier)?;
    validate_output_identifier(second_identifier)?;
    let rendered = alignment.render(first, second)?;
    writeln!(writer, ">{first_identifier}")?;
    write_wrapped(&mut writer, rendered.first().as_bytes())?;
    writeln!(writer, ">{second_identifier}")?;
    write_wrapped(&mut writer, rendered.second().as_bytes())?;
    Ok(())
}

fn write_wrapped(mut writer: impl Write, symbols: &[u8]) -> Result<(), std::io::Error> {
    if symbols.is_empty() {
        writer.write_all(b"\n")?;
        return Ok(());
    }
    for chunk in symbols.chunks(FASTA_WIDTH) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
