//! Minimal end-to-end worker for MAFFT 7.526 process benchmarks.

use std::env;
use std::error::Error;
use std::fs;
use std::io::{Error as IoError, ErrorKind};
use std::path::Path;

use align_core::{Protein, ProteinSequence};
use align_mafft_compat::{Mafft7526Alignment, Mafft7526Einsi, Mafft7526Linsi};
use align_msa::SequenceRecord;

fn invalid_data(message: impl Into<String>) -> IoError {
    IoError::new(ErrorKind::InvalidData, message.into())
}

fn read_fasta(path: &Path) -> Result<Vec<SequenceRecord<Protein>>, Box<dyn Error>> {
    let source = fs::read_to_string(path)?;
    let mut records = Vec::new();
    let mut identifier: Option<String> = None;
    let mut sequence = String::new();

    for line in source
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
    {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = identifier.replace(
                header
                    .split_ascii_whitespace()
                    .next()
                    .ok_or_else(|| invalid_data("empty FASTA header"))?
                    .to_owned(),
            ) {
                records.push(SequenceRecord::new(
                    previous,
                    ProteinSequence::from_ascii(sequence.as_bytes())?,
                )?);
                sequence.clear();
            }
        } else if identifier.is_some() {
            sequence.extend(
                line.chars()
                    .filter(|character| !character.is_ascii_whitespace()),
            );
        } else {
            return Err(invalid_data("FASTA sequence occurs before its header").into());
        }
    }

    if let Some(identifier) = identifier {
        records.push(SequenceRecord::new(
            identifier,
            ProteinSequence::from_ascii(sequence.as_bytes())?,
        )?);
    }
    if records.is_empty() {
        return Err(invalid_data("FASTA input contains no records").into());
    }
    Ok(records)
}

fn write_alignment(path: &Path, result: &Mafft7526Alignment<Protein>) -> Result<(), IoError> {
    let mut output = String::new();
    for row in result.alignment().rows() {
        output.push('>');
        output.push_str(row.identifier());
        output.push('\n');
        output.push_str(&row.to_ascii_with_gaps());
        output.push('\n');
    }
    fs::write(path, output)
}

fn main() -> Result<(), Box<dyn Error>> {
    let arguments = env::args().skip(1).collect::<Vec<_>>();
    if arguments.len() != 4 {
        return Err(IoError::new(
            ErrorKind::InvalidInput,
            "usage: process_worker <linsi|einsi> <threads> <input.fasta> <output.fasta>",
        )
        .into());
    }

    let mode = arguments[0].as_str();
    let thread_count = arguments[1].parse::<usize>()?;
    let records = read_fasta(Path::new(&arguments[2]))?;
    let result: Mafft7526Alignment<Protein> = match mode {
        "linsi" => Mafft7526Linsi::<Protein>::new()
            .with_thread_count(thread_count)?
            .align(&records)?,
        "einsi" => Mafft7526Einsi::new()
            .with_thread_count(thread_count)?
            .align(&records)?,
        _ => return Err(invalid_data(format!("unsupported profile {mode:?}")).into()),
    };

    write_alignment(Path::new(&arguments[3]), &result)?;
    eprintln!(
        "profile={} threads={} steps={} engine_score={}",
        result.profile_identifier(),
        result.thread_count(),
        result.progressive_step_count(),
        result.engine_score()
    );
    Ok(())
}
