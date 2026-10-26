use core::fmt::{self, Write as _};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::FormatError;

/// Streams any deterministic human report to a UTF-8 text file and ensures one
/// final newline.
///
/// # Errors
///
/// Returns [`FormatError::Io`] for failed formatting or output.
pub fn write_text_report(
    path: impl AsRef<Path>,
    report: &impl fmt::Display,
) -> Result<(), FormatError> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_text_report_to(&mut writer, report)?;
    writer.flush()?;
    Ok(())
}

/// Streams any deterministic human report and ensures one final newline.
///
/// # Errors
///
/// Returns [`FormatError::Io`] for failed formatting or output.
pub fn write_text_report_to(
    mut writer: impl Write,
    report: &impl fmt::Display,
) -> Result<(), FormatError> {
    let mut adapter = IoFormatWriter::new(&mut writer);
    if write!(&mut adapter, "{report}").is_err() {
        return Err(adapter
            .error
            .take()
            .unwrap_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "report formatting failed")
            })
            .into());
    }
    drop(adapter);
    writer.write_all(b"\n")?;
    Ok(())
}

struct IoFormatWriter<'a, W> {
    writer: &'a mut W,
    error: Option<io::Error>,
    trailing_newlines: usize,
}

impl<'a, W> IoFormatWriter<'a, W> {
    const fn new(writer: &'a mut W) -> Self {
        Self {
            writer,
            error: None,
            trailing_newlines: 0,
        }
    }

    fn write_all(&mut self, bytes: &[u8]) -> fmt::Result
    where
        W: Write,
    {
        self.writer.write_all(bytes).map_err(|error| {
            self.error = Some(error);
            fmt::Error
        })
    }

    fn flush_pending_newlines(&mut self) -> fmt::Result
    where
        W: Write,
    {
        const NEWLINES: &[u8; 64] = b"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
        while self.trailing_newlines != 0 {
            let count = self.trailing_newlines.min(NEWLINES.len());
            self.write_all(&NEWLINES[..count])?;
            self.trailing_newlines -= count;
        }
        Ok(())
    }
}

impl<W: Write> fmt::Write for IoFormatWriter<'_, W> {
    fn write_str(&mut self, value: &str) -> fmt::Result {
        let bytes = value.as_bytes();
        let Some(last_non_newline) = bytes.iter().rposition(|byte| *byte != b'\n') else {
            self.trailing_newlines = self
                .trailing_newlines
                .checked_add(bytes.len())
                .ok_or(fmt::Error)?;
            return Ok(());
        };
        self.flush_pending_newlines()?;
        self.write_all(&bytes[..=last_non_newline])?;
        self.trailing_newlines = bytes.len() - last_non_newline - 1;
        Ok(())
    }
}
