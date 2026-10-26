use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use align_fast::MyersMatch;

use crate::FormatError;

/// Writes Myers approximate-search matches as stable tab-separated values.
///
/// # Errors
///
/// Returns [`FormatError::Io`] for failed output.
pub fn write_myers_tsv(path: impl AsRef<Path>, matches: &[MyersMatch]) -> Result<(), FormatError> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_myers_tsv_to(&mut writer, matches)?;
    writer.flush()?;
    Ok(())
}

/// Writes Myers approximate-search matches as stable tab-separated values.
///
/// `text_end` is the native zero-based exclusive endpoint.
///
/// # Errors
///
/// Returns [`FormatError::Io`] for failed output.
pub fn write_myers_tsv_to(
    mut writer: impl Write,
    matches: &[MyersMatch],
) -> Result<(), FormatError> {
    writer.write_all(b"rank\ttext_end\tdistance\n")?;
    for (index, matched) in matches.iter().copied().enumerate() {
        writeln!(
            writer,
            "{}\t{}\t{}",
            index + 1,
            matched.end(),
            matched.distance()
        )?;
    }
    Ok(())
}

/// Convenient `.write_tsv(path)` output for Myers match slices.
pub trait MyersTsvWrite {
    /// Writes this match collection as TSV.
    ///
    /// # Errors
    ///
    /// Returns [`FormatError::Io`] for failed output.
    fn write_tsv(&self, path: impl AsRef<Path>) -> Result<(), FormatError>;
}

impl MyersTsvWrite for [MyersMatch] {
    fn write_tsv(&self, path: impl AsRef<Path>) -> Result<(), FormatError> {
        write_myers_tsv(path, self)
    }
}

impl MyersTsvWrite for Vec<MyersMatch> {
    fn write_tsv(&self, path: impl AsRef<Path>) -> Result<(), FormatError> {
        write_myers_tsv(path, self)
    }
}
