use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use align_core::Alphabet;
use align_msa::MultipleAlignment;

use crate::FormatError;

const STRICT_NAME_WIDTH: usize = 10;
const DEFAULT_BLOCK_WIDTH: usize = 60;

/// Placement of aligned sequence data in a PHYLIP document.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum PhylipLayout {
    /// Writes each complete aligned sequence before the next taxon.
    #[default]
    Sequential,
    /// Writes equal-width blocks for every taxon before continuing.
    Interleaved,
}

/// Identifier convention used by a PHYLIP document.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum PhylipNameStyle {
    /// Uses whitespace-delimited names accepted by modern relaxed-PHYLIP readers.
    #[default]
    Relaxed,
    /// Uses the original fixed ten-character taxon field.
    Strict,
}

/// Deterministic PHYLIP serialization options.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PhylipOptions {
    layout: PhylipLayout,
    names: PhylipNameStyle,
    block_width: usize,
}

impl PhylipOptions {
    /// Creates relaxed sequential PHYLIP options.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            layout: PhylipLayout::Sequential,
            names: PhylipNameStyle::Relaxed,
            block_width: DEFAULT_BLOCK_WIDTH,
        }
    }

    /// Replaces the sequence layout.
    #[must_use]
    pub const fn with_layout(mut self, layout: PhylipLayout) -> Self {
        self.layout = layout;
        self
    }

    /// Replaces the taxon-name convention.
    #[must_use]
    pub const fn with_name_style(mut self, names: PhylipNameStyle) -> Self {
        self.names = names;
        self
    }

    /// Replaces the positive interleaved block width.
    ///
    /// # Errors
    ///
    /// Returns [`FormatError::InvalidLineWidth`] for zero.
    pub const fn with_block_width(mut self, block_width: usize) -> Result<Self, FormatError> {
        if block_width == 0 {
            return Err(FormatError::InvalidLineWidth {
                format: "PHYLIP",
                value: block_width,
            });
        }
        self.block_width = block_width;
        Ok(self)
    }

    /// Returns the selected sequence layout.
    #[must_use]
    pub const fn layout(self) -> PhylipLayout {
        self.layout
    }

    /// Returns the selected name convention.
    #[must_use]
    pub const fn name_style(self) -> PhylipNameStyle {
        self.names
    }

    /// Returns the interleaved block width.
    #[must_use]
    pub const fn block_width(self) -> usize {
        self.block_width
    }
}

impl Default for PhylipOptions {
    fn default() -> Self {
        Self::new()
    }
}

/// Writes an MSA as relaxed sequential PHYLIP.
///
/// # Errors
///
/// Returns [`FormatError`] for a nonportable identifier or failed output.
pub fn write_phylip<A: Alphabet>(
    path: impl AsRef<Path>,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FormatError> {
    write_phylip_with_options(path, alignment, PhylipOptions::new())
}

/// Writes an MSA using explicit PHYLIP options.
///
/// # Errors
///
/// Returns [`FormatError`] for a nonportable identifier or failed output.
pub fn write_phylip_with_options<A: Alphabet>(
    path: impl AsRef<Path>,
    alignment: &MultipleAlignment<A>,
    options: PhylipOptions,
) -> Result<(), FormatError> {
    validate_identifiers(alignment, options.name_style())?;
    let mut writer = BufWriter::new(File::create(path)?);
    write_validated(&mut writer, alignment, options)?;
    writer.flush()?;
    Ok(())
}

/// Streams an MSA as relaxed sequential PHYLIP.
///
/// # Errors
///
/// Returns [`FormatError`] for a nonportable identifier or failed output.
pub fn write_phylip_to<A: Alphabet>(
    writer: impl Write,
    alignment: &MultipleAlignment<A>,
) -> Result<(), FormatError> {
    write_phylip_with_options_to(writer, alignment, PhylipOptions::new())
}

/// Streams an MSA using explicit PHYLIP options.
///
/// # Errors
///
/// Returns [`FormatError`] for a nonportable identifier or failed output.
pub fn write_phylip_with_options_to<A: Alphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
    options: PhylipOptions,
) -> Result<(), FormatError> {
    validate_identifiers(alignment, options.name_style())?;
    write_validated(&mut writer, alignment, options)
}

fn write_validated<A: Alphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
    options: PhylipOptions,
) -> Result<(), FormatError> {
    writeln!(
        writer,
        "{} {}",
        alignment.sequence_count(),
        alignment.width()
    )?;
    match options.layout() {
        PhylipLayout::Sequential => write_sequential(&mut writer, alignment, options.name_style()),
        PhylipLayout::Interleaved => write_interleaved(&mut writer, alignment, options),
    }
}

fn write_sequential<A: Alphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
    names: PhylipNameStyle,
) -> Result<(), FormatError> {
    for row in alignment.rows() {
        write_name(&mut writer, row.identifier(), names)?;
        writeln!(writer, "{}", row.to_ascii_with_gaps())?;
    }
    Ok(())
}

fn write_interleaved<A: Alphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
    options: PhylipOptions,
) -> Result<(), FormatError> {
    if alignment.width() == 0 {
        for row in alignment.rows() {
            write_name(&mut writer, row.identifier(), options.name_style())?;
            writeln!(writer)?;
        }
        return Ok(());
    }
    for (block_index, start) in (0..alignment.width())
        .step_by(options.block_width())
        .enumerate()
    {
        if block_index != 0 {
            writeln!(writer)?;
        }
        let end = start
            .saturating_add(options.block_width())
            .min(alignment.width());
        for row in alignment.rows() {
            if block_index == 0 {
                write_name(&mut writer, row.identifier(), options.name_style())?;
            } else {
                write_continuation_indent(&mut writer, alignment, options.name_style())?;
            }
            for cell in &row.symbols()[start..end] {
                writer.write_all(&[cell.map_or(b'-', align_core::Symbol::as_ascii)])?;
            }
            writeln!(writer)?;
        }
    }
    Ok(())
}

fn write_name(
    mut writer: impl Write,
    identifier: &str,
    style: PhylipNameStyle,
) -> Result<(), FormatError> {
    match style {
        PhylipNameStyle::Relaxed => write!(writer, "{identifier} ")?,
        PhylipNameStyle::Strict => write!(writer, "{identifier:<STRICT_NAME_WIDTH$}")?,
    }
    Ok(())
}

fn write_continuation_indent<A: Alphabet>(
    mut writer: impl Write,
    alignment: &MultipleAlignment<A>,
    style: PhylipNameStyle,
) -> Result<(), FormatError> {
    let width = match style {
        PhylipNameStyle::Strict => STRICT_NAME_WIDTH,
        PhylipNameStyle::Relaxed => {
            alignment
                .rows()
                .iter()
                .map(|row| row.identifier().len())
                .max()
                .unwrap_or(0)
                + 1
        }
    };
    write!(writer, "{:width$}", "")?;
    Ok(())
}

fn validate_identifiers<A: Alphabet>(
    alignment: &MultipleAlignment<A>,
    style: PhylipNameStyle,
) -> Result<(), FormatError> {
    for row in alignment.rows() {
        let identifier = row.identifier();
        if !identifier.is_ascii() {
            return identifier_error(identifier, "must contain portable ASCII only");
        }
        if identifier
            .bytes()
            .any(|byte| !byte.is_ascii_graphic() && byte != b' ')
        {
            return identifier_error(identifier, "contains a non-printing character");
        }
        if identifier
            .bytes()
            .any(|byte| matches!(byte, b'(' | b')' | b'[' | b']' | b':' | b';' | b','))
        {
            return identifier_error(identifier, "contains PHYLIP-reserved punctuation");
        }
        match style {
            PhylipNameStyle::Relaxed
                if identifier.bytes().any(|byte| byte.is_ascii_whitespace()) =>
            {
                return identifier_error(identifier, "relaxed names must not contain whitespace");
            }
            PhylipNameStyle::Strict if identifier.len() > STRICT_NAME_WIDTH => {
                return identifier_error(identifier, "strict names must be at most 10 bytes");
            }
            PhylipNameStyle::Strict if identifier.trim() != identifier => {
                return identifier_error(
                    identifier,
                    "strict names must not start or end with whitespace",
                );
            }
            PhylipNameStyle::Relaxed | PhylipNameStyle::Strict => {}
        }
    }
    Ok(())
}

fn identifier_error<T>(identifier: &str, reason: &'static str) -> Result<T, FormatError> {
    Err(FormatError::InvalidFormatIdentifier {
        format: "PHYLIP",
        identifier: identifier.to_owned(),
        reason,
    })
}

/// Convenient PHYLIP output for a rectangular MSA.
pub trait PhylipWrite {
    /// Writes relaxed sequential PHYLIP.
    ///
    /// # Errors
    ///
    /// Returns [`FormatError`] for invalid identifiers or failed output.
    fn write_phylip(&self, path: impl AsRef<Path>) -> Result<(), FormatError>;

    /// Writes PHYLIP with explicit layout and name rules.
    ///
    /// # Errors
    ///
    /// Returns [`FormatError`] for invalid identifiers or failed output.
    fn write_phylip_with(
        &self,
        path: impl AsRef<Path>,
        options: PhylipOptions,
    ) -> Result<(), FormatError>;
}

impl<A: Alphabet> PhylipWrite for MultipleAlignment<A> {
    fn write_phylip(&self, path: impl AsRef<Path>) -> Result<(), FormatError> {
        write_phylip(path, self)
    }

    fn write_phylip_with(
        &self,
        path: impl AsRef<Path>,
        options: PhylipOptions,
    ) -> Result<(), FormatError> {
        write_phylip_with_options(path, self, options)
    }
}
