//! Exact replay of the existing opaque MAFFT 7.526 L-INS-i fixture.

use align_core::{Dna, Protein, Rna, Sequence};
use align_mafft::MafftAlphabet;
use align_mafft_compat::Mafft7526Linsi;
use align_msa::SequenceRecord;
use serde::Deserialize;

const FIXTURE: &str =
    include_str!("../../align-mafft/tests/fixtures/mafft-7.526/l-ins-i-localpair-v1.json");

#[derive(Debug, Deserialize)]
struct Fixture {
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Case {
    id: String,
    alphabet: String,
    sequences: Vec<SequenceValue>,
    reference_alignment: Vec<AlignedValue>,
}

#[derive(Debug, Deserialize)]
struct SequenceValue {
    id: String,
    sequence: String,
}

#[derive(Debug, Deserialize)]
struct AlignedValue {
    id: String,
    aligned: String,
}

fn run_case<A: MafftAlphabet>(case: &Case, threads: usize) -> Vec<(String, String)> {
    let records = case
        .sequences
        .iter()
        .map(|value| {
            SequenceRecord::new(
                value.id.clone(),
                Sequence::<A>::from_ascii(value.sequence.as_bytes()).unwrap(),
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    Mafft7526Linsi::<A>::new()
        .with_thread_count(threads)
        .unwrap()
        .align(&records)
        .unwrap()
        .alignment()
        .rows()
        .iter()
        .map(|row| (row.identifier().to_owned(), row.to_ascii_with_gaps()))
        .collect()
}

#[test]
fn all_opaque_native_linsi_cases_are_exact() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).unwrap();
    assert_eq!(fixture.cases.len(), 7);
    for case in &fixture.cases {
        let actual = match case.alphabet.as_str() {
            "dna" => run_case::<Dna>(case, 1),
            "rna" => run_case::<Rna>(case, 1),
            "protein" => run_case::<Protein>(case, 1),
            other => panic!("unexpected fixture alphabet {other}"),
        };
        let expected = case
            .reference_alignment
            .iter()
            .map(|row| (row.id.clone(), row.aligned.clone()))
            .collect::<Vec<_>>();
        assert_eq!(actual, expected, "{} diverged from MAFFT 7.526", case.id);
    }
}

#[test]
fn divergent_nucleotide_gap_scaling_matches_native_at_one_and_four_workers() {
    let fixture: Fixture = serde_json::from_str(include_str!(
        "fixtures/mafft-7.526/linsi-nucleotide-pair-gap-v1.json"
    ))
    .unwrap();
    assert_eq!(fixture.cases.len(), 2);
    for case in &fixture.cases {
        let expected = case
            .reference_alignment
            .iter()
            .map(|row| (row.id.clone(), row.aligned.clone()))
            .collect::<Vec<_>>();
        for threads in [1, 4] {
            let actual = match case.alphabet.as_str() {
                "dna" => run_case::<Dna>(case, threads),
                "rna" => run_case::<Rna>(case, threads),
                other => panic!("unexpected regression alphabet {other}"),
            };
            assert_eq!(actual, expected, "{} at {threads} workers", case.id);
        }
    }
}
