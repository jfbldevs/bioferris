use align_core::Sequence;
use align_msa::{DistanceMatrix, SequenceRecord};

use crate::matrix::{MoleculeKind, molecule_kind};
use crate::{FftNs1Parameters, MafftAlphabet, MafftError};

#[derive(Clone, Debug)]
pub(crate) struct TupleCollection {
    effective_length: usize,
    codes: Vec<Box<[usize]>>,
    stored_bytes: usize,
}

impl TupleCollection {
    pub(crate) fn from_records<A: MafftAlphabet>(
        records: &[SequenceRecord<A>],
        parameters: FftNs1Parameters,
    ) -> Result<Self, MafftError> {
        let effective_length = records
            .iter()
            .map(|record| record.sequence().len())
            .filter(|length| *length > 0)
            .min()
            .map_or(1, |length| length.min(parameters.tuple_length()));
        let kind = molecule_kind::<A>();
        let radix = match kind {
            MoleculeKind::Protein => 6,
            MoleculeKind::Nucleotide => 4,
        };
        let mut codes = Vec::with_capacity(records.len());
        let mut stored_bytes = 0usize;
        for record in records {
            let mut sequence_codes =
                encode_sequence(record.sequence(), kind, radix, effective_length)?;
            sequence_codes.sort_unstable();
            stored_bytes = stored_bytes
                .checked_add(sequence_codes.len().checked_mul(size_of::<usize>()).ok_or(
                    MafftError::ArithmeticOverflow {
                        stage: "compressed tuple storage",
                    },
                )?)
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "compressed tuple storage",
                })?;
            codes.push(sequence_codes.into_boxed_slice());
        }
        Ok(Self {
            effective_length,
            codes,
            stored_bytes,
        })
    }

    pub(crate) const fn effective_length(&self) -> usize {
        self.effective_length
    }

    pub(crate) const fn stored_bytes(&self) -> usize {
        self.stored_bytes
    }

    pub(crate) fn distance_matrix<A: MafftAlphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<DistanceMatrix, MafftError> {
        let count = records.len();
        let cell_count = count
            .checked_mul(count)
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "6-tuple distance matrix",
            })?;
        let mut distances = vec![0.0; cell_count];
        for first in 0..count {
            for second in (first + 1)..count {
                let first_codes = &self.codes[first];
                let second_codes = &self.codes[second];
                let denominator = first_codes.len().min(second_codes.len());
                let distance = if denominator == 0 {
                    if records[first].sequence() == records[second].sequence() {
                        0.0
                    } else {
                        1.0
                    }
                } else {
                    let shared = multiset_intersection_count(first_codes, second_codes);
                    #[allow(clippy::cast_precision_loss)]
                    let fraction = shared as f64 / denominator as f64;
                    (1.0 - fraction).clamp(0.0, 1.0)
                };
                distances[first * count + second] = distance;
                distances[second * count + first] = distance;
            }
        }
        DistanceMatrix::try_from_flat(count, distances).map_err(Into::into)
    }
}

fn encode_sequence<A: MafftAlphabet>(
    sequence: &Sequence<A>,
    kind: MoleculeKind,
    radix: usize,
    tuple_length: usize,
) -> Result<Vec<usize>, MafftError> {
    if sequence.len() < tuple_length {
        return Ok(Vec::new());
    }
    let capacity = sequence.len() - tuple_length + 1;
    let mut codes = Vec::with_capacity(capacity);
    for window in sequence.as_slice().windows(tuple_length) {
        let mut code = 0usize;
        let mut valid = true;
        for symbol in window {
            let Some(digit) = compressed_symbol(kind, symbol.as_ascii()) else {
                valid = false;
                break;
            };
            code = code
                .checked_mul(radix)
                .and_then(|value| value.checked_add(digit))
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "compressed tuple encoding",
                })?;
        }
        if valid {
            codes.push(code);
        }
    }
    Ok(codes)
}

fn compressed_symbol(kind: MoleculeKind, ascii: u8) -> Option<usize> {
    match kind {
        MoleculeKind::Protein => match ascii {
            b'A' | b'G' | b'P' | b'S' | b'T' => Some(0),
            b'D' | b'E' | b'N' | b'Q' | b'B' | b'Z' => Some(1),
            b'H' | b'K' | b'R' | b'O' => Some(2),
            b'I' | b'L' | b'M' | b'V' | b'J' => Some(3),
            b'F' | b'W' | b'Y' => Some(4),
            b'C' | b'U' => Some(5),
            b'X' => None,
            _ => unreachable!("typed protein symbol must be known"),
        },
        MoleculeKind::Nucleotide => match ascii {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' | b'U' => Some(3),
            _ => None,
        },
    }
}

fn multiset_intersection_count(first: &[usize], second: &[usize]) -> usize {
    let mut first_index = 0;
    let mut second_index = 0;
    let mut shared = 0;
    while first_index < first.len() && second_index < second.len() {
        match first[first_index].cmp(&second[second_index]) {
            core::cmp::Ordering::Less => first_index += 1,
            core::cmp::Ordering::Greater => second_index += 1,
            core::cmp::Ordering::Equal => {
                shared += 1;
                first_index += 1;
                second_index += 1;
            }
        }
    }
    shared
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence};
    use align_msa::SequenceRecord;

    use super::TupleCollection;
    use crate::FftNs1Parameters;

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    #[test]
    fn hand_calculated_multiset_distance_uses_repeated_tuples() {
        let records = vec![dna("a", "AAAAAAA"), dna("b", "AAAAAAC")];
        let tuples =
            TupleCollection::from_records(&records, FftNs1Parameters::paper_2002()).unwrap();
        let distances = tuples.distance_matrix(&records).unwrap();
        // AA...A has two AAAAAA tuples; the second sequence has one.
        assert_eq!(distances.get(0, 1), Some(0.5));
    }

    #[test]
    fn short_collection_uses_one_common_shorter_tuple_length() {
        let records = vec![dna("a", "AC"), dna("b", "AG"), dna("c", "ACGT")];
        let tuples =
            TupleCollection::from_records(&records, FftNs1Parameters::paper_2002()).unwrap();
        assert_eq!(tuples.effective_length(), 2);
        let distances = tuples.distance_matrix(&records).unwrap();
        assert_eq!(distances.get(0, 1), Some(1.0));
        assert_eq!(distances.get(0, 2), Some(0.0));
    }

    #[test]
    fn empty_distance_is_total_and_symmetric() {
        let records = vec![dna("e1", ""), dna("e2", ""), dna("a", "A")];
        let tuples =
            TupleCollection::from_records(&records, FftNs1Parameters::paper_2002()).unwrap();
        let distances = tuples.distance_matrix(&records).unwrap();
        assert_eq!(distances.get(0, 1), Some(0.0));
        assert_eq!(distances.get(0, 2), Some(1.0));
        assert_eq!(distances.get(2, 0), Some(1.0));
    }

    #[test]
    fn protein_dayhoff_classes_make_conservative_tuples_identical() {
        let records = vec![
            SequenceRecord::<Protein>::new("a", ProteinSequence::from_ascii("AGPSTA").unwrap())
                .unwrap(),
            SequenceRecord::<Protein>::new("b", ProteinSequence::from_ascii("TSPGAT").unwrap())
                .unwrap(),
        ];
        let tuples =
            TupleCollection::from_records(&records, FftNs1Parameters::paper_2002()).unwrap();
        assert_eq!(
            tuples.distance_matrix(&records).unwrap().get(0, 1),
            Some(0.0)
        );
    }
}
