//! Runs the published protein Clustal Omega profile on one external `HomFam` family.
//!
//! Pass the family's unaligned `*_test-only.vie` file and aligned `*_ref.vie`
//! Homstrad reference. Neither dataset file is vendored with this crate.

use std::collections::HashSet;
use std::env;
use std::error::Error;
use std::fs;
use std::io::{Error as IoError, ErrorKind};
use std::path::Path;
use std::time::Instant;

use align_clustal::ClustalOmega;
use align_core::{Protein, ProteinSequence, Symbol};
use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord, compare_alignments};
use serde_json::json;

fn invalid_data(message: impl Into<String>) -> IoError {
    IoError::new(ErrorKind::InvalidData, message.into())
}

fn fasta_entries(path: &Path) -> Result<Vec<(String, String)>, Box<dyn Error>> {
    let source = fs::read_to_string(path)?;
    let mut entries = Vec::new();
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
                entries.push((previous, core::mem::take(&mut sequence)));
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
        entries.push((identifier, sequence));
    }
    if entries.is_empty() {
        return Err(invalid_data("FASTA input contains no records").into());
    }
    Ok(entries)
}

fn parse_unaligned(path: &Path) -> Result<Vec<SequenceRecord<Protein>>, Box<dyn Error>> {
    fasta_entries(path)?
        .into_iter()
        .map(|(identifier, sequence)| {
            Ok(SequenceRecord::new(
                identifier,
                ProteinSequence::from_ascii(sequence.as_bytes())?,
            )?)
        })
        .collect()
}

fn parse_reference(path: &Path) -> Result<MultipleAlignment<Protein>, Box<dyn Error>> {
    let rows = fasta_entries(path)?
        .into_iter()
        .map(|(identifier, sequence)| {
            let symbols = sequence
                .bytes()
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
            Ok(AlignedRow::new(identifier, symbols)?)
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    Ok(MultipleAlignment::try_new(rows)?)
}

fn embedded_candidate(
    alignment: &MultipleAlignment<Protein>,
    reference: &MultipleAlignment<Protein>,
) -> Result<MultipleAlignment<Protein>, Box<dyn Error>> {
    let identifiers = reference
        .rows()
        .iter()
        .map(align_msa::AlignedRow::identifier)
        .collect::<HashSet<_>>();
    let matched = alignment
        .rows()
        .iter()
        .filter(|row| identifiers.contains(row.identifier()))
        .collect::<Vec<_>>();
    if matched.len() != reference.sequence_count() {
        return Err(invalid_data("candidate is missing an embedded Homstrad sequence").into());
    }
    let retained_columns = (0..alignment.width())
        .filter(|column| matched.iter().any(|row| row.symbols()[*column].is_some()))
        .collect::<Vec<_>>();
    let rows = matched
        .into_iter()
        .map(|row| {
            AlignedRow::new(
                row.identifier(),
                retained_columns
                    .iter()
                    .map(|column| row.symbols()[*column])
                    .collect(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(MultipleAlignment::try_new(rows)?)
}

fn main() -> Result<(), Box<dyn Error>> {
    let arguments = env::args().skip(1).collect::<Vec<_>>();
    if !(2..=3).contains(&arguments.len()) {
        return Err(IoError::new(
            ErrorKind::InvalidInput,
            "usage: clustal_omega_homfam <test-only.vie> <ref.vie> [thread-count]",
        )
        .into());
    }
    let thread_count = arguments
        .get(2)
        .map_or(Ok(1), |value| value.parse::<usize>())?;
    let input_path = Path::new(&arguments[0]);
    let reference_path = Path::new(&arguments[1]);
    let mut records = parse_unaligned(input_path)?;
    let reference = parse_reference(reference_path)?;
    let mut seen = records
        .iter()
        .map(|record| record.identifier().to_owned())
        .collect::<HashSet<_>>();
    for row in reference.rows() {
        if !seen.insert(row.identifier().to_owned()) {
            return Err(invalid_data(format!(
                "duplicate test/reference identifier {}",
                row.identifier()
            ))
            .into());
        }
        records.push(SequenceRecord::new(
            row.identifier(),
            ProteinSequence::from_ascii(row.to_ungapped_ascii().as_bytes())?,
        )?);
    }

    let parameters = ClustalOmega::<Protein>::published()
        .parameters()
        .clone()
        .with_thread_count(thread_count)?;
    let aligner = ClustalOmega::<Protein>::with_parameters(parameters)?;
    let started = Instant::now();
    let result = aligner.align(&records)?;
    let elapsed = started.elapsed();
    let candidate = embedded_candidate(result.alignment(), &reference)?;
    let metrics = compare_alignments(&candidate, &reference)?;

    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "profile": result.profile_identifier(),
            "family": input_path
                .file_name()
                .and_then(|name| name.to_str())
                .map(|name| name.trim_end_matches("_test-only.vie")),
            "thread_count": thread_count,
            "test_sequence_count": records.len() - reference.sequence_count(),
            "reference_sequence_count": reference.sequence_count(),
            "total_sequence_count": records.len(),
            "reference_width": reference.width(),
            "candidate_width": result.alignment().width(),
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
