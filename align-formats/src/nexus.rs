use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use align_core::{Alphabet, Dna, Protein, Rna};
use align_msa::MultipleAlignment;

use crate::FormatError;

/// Built-in biological alphabets with an exact phylogenetic NEXUS datatype.
pub trait NexusAlphabet: Alphabet {
    /// NEXUS `DATATYPE` token.
    const NEXUS_DATATYPE: &'static str;
}

impl NexusAlphabet for Dna {
    const NEXUS_DATATYPE: &'static str = "DNA";
}

impl NexusAlphabet for Rna {
    const NEXUS_DATATYPE: &'static str = "RNA";
}

impl NexusAlphabet for Protein {
    const NEXUS_DATATYPE: &'static str = "PROTEIN";
}

/// Writes a phylogenetic NEXUS `DATA` block containing one rectangular MSA.
///
/// # Errors
///
/// Returns [`FormatError::Io`] when the destination cannot be created or written.
pub fn write_nexus<A: NexusAlphabet>(
    path: impl AsRef<Path>,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FormatError> {
    validate_identifiers(alignment)?;
    let mut writer = BufWriter::new(File::create(path)?);
    write_validated(&mut writer, alignment)?;
    writer.flush()?;
    Ok(())
}

/// Streams a phylogenetic NEXUS `DATA` block containing one rectangular MSA.
///
/// Taxon labels are always single-quoted and embedded quotes are doubled.
///
/// # Errors
///
/// Returns [`FormatError::Io`] when the destination cannot be written.
pub fn write_nexus_to<A: NexusAlphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FormatError> {
    validate_identifiers(alignment)?;
    write_validated(&mut writer, alignment)
}

fn write_validated<A: NexusAlphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FormatError> {
    writer.write_all(b"#NEXUS\n\nBEGIN DATA;\n")?;
    writeln!(
        writer,
        "  DIMENSIONS NTAX={} NCHAR={};",
        alignment.sequence_count(),
        alignment.width()
    )?;
    writeln!(
        writer,
        "  FORMAT DATATYPE={} MISSING=? GAP=-;",
        A::NEXUS_DATATYPE
    )?;
    writer.write_all(b"  MATRIX\n")?;
    for row in alignment.rows() {
        writer.write_all(b"    ")?;
        write_quoted_label(&mut writer, row.identifier())?;
        writeln!(writer, " {}", row.to_ascii_with_gaps())?;
    }
    writer.write_all(b"  ;\nEND;\n")?;
    Ok(())
}

fn validate_identifiers<A: Alphabet>(alignment: &MultipleAlignment<A>) -> Result<(), FormatError> {
    for row in alignment.rows() {
        if row.identifier().chars().any(char::is_control) {
            return Err(FormatError::InvalidFormatIdentifier {
                format: "NEXUS",
                identifier: row.identifier().to_owned(),
                reason: "taxon labels must not contain control characters",
            });
        }
    }
    Ok(())
}

fn write_quoted_label(mut writer: impl Write, identifier: &str) -> Result<(), FormatError> {
    writer.write_all(b"'")?;
    for part in identifier.split_inclusive('\'') {
        writer.write_all(part.as_bytes())?;
        if part.ends_with('\'') {
            writer.write_all(b"'")?;
        }
    }
    writer.write_all(b"'")?;
    Ok(())
}

/// Convenient phylogenetic NEXUS output for a rectangular MSA.
pub trait NexusWrite {
    /// Writes this alignment as a NEXUS `DATA` block.
    ///
    /// # Errors
    ///
    /// Returns [`FormatError::Io`] for failed output.
    fn write_nexus(&self, path: impl AsRef<Path>) -> Result<(), FormatError>;
}

impl<A: NexusAlphabet> NexusWrite for MultipleAlignment<A> {
    fn write_nexus(&self, path: impl AsRef<Path>) -> Result<(), FormatError> {
        write_nexus(path, self)
    }
}
