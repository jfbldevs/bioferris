//! Exhaustive small-alphabet validation of both public MUSCLE 5 modes.

use align_core::{Dna, DnaSequence};
use align_msa::SequenceRecord;
use align_muscle::{Muscle5, Muscle5Parameters};

fn binary_strings(maximum_length: usize) -> Vec<Vec<u8>> {
    let mut strings = vec![Vec::new()];
    for length in 1..=maximum_length {
        for bits in 0..(1_usize << length) {
            strings.push(
                (0..length)
                    .map(|offset| {
                        if bits & (1 << offset) == 0 {
                            b'A'
                        } else {
                            b'C'
                        }
                    })
                    .collect(),
            );
        }
    }
    strings
}

fn records(raw: [&[u8]; 3]) -> [SequenceRecord<Dna>; 3] {
    [
        SequenceRecord::new("a", DnaSequence::from_ascii(raw[0]).unwrap()).unwrap(),
        SequenceRecord::new("b", DnaSequence::from_ascii(raw[1]).unwrap()).unwrap(),
        SequenceRecord::new("c", DnaSequence::from_ascii(raw[2]).unwrap()).unwrap(),
    ]
}

fn assert_replay(result: &align_muscle::Muscle5Alignment<Dna>, raw: [&[u8]; 3]) {
    for (row, source) in result.alignment().rows().iter().zip(raw) {
        assert_eq!(row.to_ungapped_ascii().as_bytes(), source);
        assert_eq!(row.len(), result.alignment().width());
    }
    for column in 0..result.alignment().width() {
        assert!(
            result
                .alignment()
                .column(column)
                .unwrap()
                .any(|cell| cell.is_some())
        );
    }
}

#[test]
fn every_ordered_triple_through_length_two_is_valid_and_deterministic() {
    let strings = binary_strings(2);
    let aligner =
        Muscle5::<Dna>::with_parameters(Muscle5Parameters::published_align().with_iterations(1, 3))
            .unwrap();
    let mut cases = 0_usize;
    for first in &strings {
        for second in &strings {
            for third in &strings {
                let raw = [first.as_slice(), second.as_slice(), third.as_slice()];
                let records = records(raw);
                let first_result = aligner.align(&records).unwrap();
                let second_result = aligner.align(&records).unwrap();
                assert_replay(&first_result, raw);
                assert_eq!(first_result.alignment(), second_result.alignment());
                assert_eq!(
                    first_result.consistent_posterior_checksum(),
                    second_result.consistent_posterior_checksum()
                );
                cases += 1;
            }
        }
    }
    assert_eq!(cases, strings.len().pow(3));
}

#[test]
fn super5_covers_every_ordered_binary_pair_through_length_two() {
    let strings = binary_strings(2);
    let parameters = Muscle5Parameters::published_super5()
        .with_iterations(0, 1)
        .with_super5_limits(0.0, 0.30, 0.10, 2, 16, 2)
        .unwrap();
    let aligner = Muscle5::<Dna>::with_parameters(parameters).unwrap();
    for first in &strings {
        for second in &strings {
            let raw = [first.as_slice(), second.as_slice(), b"A".as_slice()];
            let records = records(raw);
            let result = aligner.align(&records).unwrap();
            assert_replay(&result, raw);
            assert!(result.super5_evidence().is_some());
        }
    }
}
