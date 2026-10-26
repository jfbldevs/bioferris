//! Runs the published protein Clustal Omega profile on one external `BAliBASE` 3 case.
//!
//! The benchmark data is intentionally not vendored. Pass the unaligned `.tfa`
//! input and matching annotated `.xml` reference from the official distribution.

use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::fs;
use std::io::{Error as IoError, ErrorKind};
use std::path::Path;
use std::time::Instant;

use align_clustal::ClustalOmega;
use align_core::{Protein, ProteinSequence, Symbol};
use align_msa::{
    AlignedRow, MultipleAlignment, SequenceRecord, compare_alignments_on_reference_columns,
};
use serde_json::json;

fn invalid_data(message: impl Into<String>) -> IoError {
    IoError::new(ErrorKind::InvalidData, message.into())
}

fn extract_tag<'a>(source: &'a str, tag: &str) -> Result<&'a str, IoError> {
    let opening = format!("<{tag}>");
    let closing = format!("</{tag}>");
    let start = source
        .find(&opening)
        .ok_or_else(|| invalid_data(format!("missing {opening}")))?
        + opening.len();
    let relative_end = source[start..]
        .find(&closing)
        .ok_or_else(|| invalid_data(format!("missing {closing}")))?;
    Ok(&source[start..start + relative_end])
}

fn parse_fasta(path: &Path) -> Result<Vec<SequenceRecord<Protein>>, Box<dyn Error>> {
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

fn parse_reference(path: &Path) -> Result<(MultipleAlignment<Protein>, Vec<bool>), Box<dyn Error>> {
    let source = fs::read_to_string(path)?;
    let mut rows = Vec::new();
    let mut remaining = source.as_str();
    while let Some(start) = remaining.find("<sequence ") {
        remaining = &remaining[start..];
        let end = remaining
            .find("</sequence>")
            .ok_or_else(|| invalid_data("unterminated XML sequence"))?
            + "</sequence>".len();
        let block = &remaining[..end];
        let identifier = extract_tag(block, "seq-name")?.trim();
        let aligned = extract_tag(block, "seq-data")?;
        let symbols = aligned
            .bytes()
            .filter(|byte| !byte.is_ascii_whitespace())
            .map(|byte| {
                if byte == b'-' || byte == b'.' {
                    Ok(None)
                } else {
                    Symbol::<Protein>::from_ascii(byte)
                        .map(Some)
                        .ok_or_else(|| {
                            invalid_data(format!(
                                "invalid protein byte {byte:?} in reference row {identifier}"
                            ))
                        })
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        rows.push(AlignedRow::new(identifier, symbols)?);
        remaining = &remaining[end..];
    }
    let reference = MultipleAlignment::try_new(rows)?;
    let scores = extract_tag(&source, "colsco-data")?
        .split_ascii_whitespace()
        .map(str::parse::<i32>)
        .collect::<Result<Vec<_>, _>>()?;
    if scores.len() != reference.width() {
        return Err(invalid_data(format!(
            "core-column vector has width {}, expected {}",
            scores.len(),
            reference.width()
        ))
        .into());
    }
    Ok((
        reference,
        scores.into_iter().map(|score| score > 0).collect(),
    ))
}

fn validate_sequence_sets(
    records: &[SequenceRecord<Protein>],
    reference: &MultipleAlignment<Protein>,
) -> Result<(), IoError> {
    let reference_sequences = reference
        .rows()
        .iter()
        .map(|row| (row.identifier(), row.to_ungapped_ascii()))
        .collect::<HashMap<_, _>>();
    if reference_sequences.len() != records.len() {
        return Err(invalid_data("input and reference sequence counts differ"));
    }
    for record in records {
        let Some(reference_sequence) = reference_sequences.get(record.identifier()) else {
            return Err(invalid_data(format!(
                "reference is missing input row {}",
                record.identifier()
            )));
        };
        if reference_sequence != &record.sequence().to_ascii_string() {
            return Err(invalid_data(format!(
                "reference sequence differs for {}",
                record.identifier()
            )));
        }
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let arguments = env::args().skip(1).collect::<Vec<_>>();
    if !(2..=3).contains(&arguments.len()) {
        return Err(IoError::new(
            ErrorKind::InvalidInput,
            "usage: clustal_omega_balibase <unaligned.tfa> <reference.xml> [thread-count]",
        )
        .into());
    }
    let thread_count = arguments
        .get(2)
        .map_or(Ok(1), |value| value.parse::<usize>())?;
    let input_path = Path::new(&arguments[0]);
    let reference_path = Path::new(&arguments[1]);
    let records = parse_fasta(input_path)?;
    let (reference, core_columns) = parse_reference(reference_path)?;
    validate_sequence_sets(&records, &reference)?;

    let parameters = ClustalOmega::<Protein>::published()
        .parameters()
        .clone()
        .with_thread_count(thread_count)?;
    let aligner = ClustalOmega::<Protein>::with_parameters(parameters)?;
    let started = Instant::now();
    let result = aligner.align(&records)?;
    let elapsed = started.elapsed();
    let metrics =
        compare_alignments_on_reference_columns(result.alignment(), &reference, &core_columns)?;
    let maximum_input_length = records
        .iter()
        .map(|record| record.sequence().len())
        .max()
        .unwrap_or(0);

    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "profile": result.profile_identifier(),
            "case": input_path.file_stem().and_then(|name| name.to_str()),
            "thread_count": thread_count,
            "sequence_count": records.len(),
            "maximum_input_length": maximum_input_length,
            "reference_width": reference.width(),
            "candidate_width": result.alignment().width(),
            "core_columns": core_columns.iter().filter(|selected| **selected).count(),
            "sp": metrics.sp_score(),
            "tc": metrics.tc_score(),
            "correct_pairs": metrics.correct_pairs(),
            "reference_pairs": metrics.reference_pairs(),
            "correct_columns": metrics.correct_columns(),
            "reference_columns": metrics.reference_columns(),
            "elapsed_ns": elapsed.as_nanos(),
            "estimated_peak_bytes": result.memory_estimate().estimated_peak_bytes()
        }))?
    );
    Ok(())
}
