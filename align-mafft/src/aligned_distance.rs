use core::mem::size_of;

use align_msa::MultipleAlignment;

use crate::{MafftAlphabet, MafftError};

/// Pairwise Poisson-corrected distances inferred from an FFT-NS-1 alignment.
///
/// Values are finite and nonnegative but are not restricted to one because
/// they represent corrected expected substitutions per comparable site.
#[derive(Clone, Debug, PartialEq)]
pub struct FftNs2DistanceMatrix {
    size: usize,
    values: Box<[f64]>,
    comparable_columns: Box<[usize]>,
    identical_columns: Box<[usize]>,
}

impl FftNs2DistanceMatrix {
    /// Returns the square matrix dimension.
    #[must_use]
    pub const fn size(&self) -> usize {
        self.size
    }

    /// Returns one corrected distance when both indices are in range.
    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> Option<f64> {
        self.offset(row, column).map(|offset| self.values[offset])
    }

    /// Returns all corrected distances in row-major order.
    #[must_use]
    pub const fn as_flat_slice(&self) -> &[f64] {
        &self.values
    }

    /// Returns the pairwise-deletion column count for one matrix cell.
    #[must_use]
    pub fn comparable_columns(&self, row: usize, column: usize) -> Option<usize> {
        self.offset(row, column)
            .map(|offset| self.comparable_columns[offset])
    }

    /// Returns the exact-symbol identity count among comparable columns.
    #[must_use]
    pub fn identical_columns(&self, row: usize, column: usize) -> Option<usize> {
        self.offset(row, column)
            .map(|offset| self.identical_columns[offset])
    }

    pub(crate) fn from_alignment<A: MafftAlphabet>(
        alignment: &MultipleAlignment<A>,
        identity_floor: f64,
    ) -> Result<Self, MafftError> {
        if !identity_floor.is_finite() || !(0.0..=1.0).contains(&identity_floor) {
            return Err(MafftError::InternalInvariant {
                message: "FFT-NS-2 Poisson identity floor is invalid",
            });
        }
        if identity_floor == 0.0 {
            return Err(MafftError::InternalInvariant {
                message: "FFT-NS-2 Poisson identity floor must be positive",
            });
        }
        let size = alignment.sequence_count();
        let cell_count = size
            .checked_mul(size)
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-2 aligned distance dimensions",
            })?;
        let mut values = vec![0.0; cell_count];
        let mut comparable_columns = vec![0usize; cell_count];
        let mut identical_columns = vec![0usize; cell_count];
        let ungapped_lengths = alignment
            .rows()
            .iter()
            .map(|row| row.symbols().iter().filter(|cell| cell.is_some()).count())
            .collect::<Vec<_>>();
        let saturation = -identity_floor.ln();

        for row in 0..size {
            let diagonal = row * size + row;
            comparable_columns[diagonal] = ungapped_lengths[row];
            identical_columns[diagonal] = ungapped_lengths[row];
            for column in (row + 1)..size {
                let mut comparable = 0usize;
                let mut identical = 0usize;
                for (left, right) in alignment.rows()[row]
                    .symbols()
                    .iter()
                    .zip(alignment.rows()[column].symbols())
                {
                    if let (Some(left), Some(right)) = (left, right) {
                        comparable =
                            comparable
                                .checked_add(1)
                                .ok_or(MafftError::ArithmeticOverflow {
                                    stage: "FFT-NS-2 comparable column count",
                                })?;
                        if left == right {
                            identical =
                                identical
                                    .checked_add(1)
                                    .ok_or(MafftError::ArithmeticOverflow {
                                        stage: "FFT-NS-2 identical column count",
                                    })?;
                        }
                    }
                }
                let distance = if comparable == 0 {
                    if ungapped_lengths[row] == 0 && ungapped_lengths[column] == 0 {
                        0.0
                    } else {
                        saturation
                    }
                } else {
                    #[allow(clippy::cast_precision_loss)]
                    let unchanged = identical as f64 / comparable as f64;
                    -unchanged.max(identity_floor).ln()
                };
                if !distance.is_finite() || distance < 0.0 {
                    return Err(MafftError::InternalInvariant {
                        message: "FFT-NS-2 Poisson distance is invalid",
                    });
                }
                let forward = row * size + column;
                let reverse = column * size + row;
                values[forward] = distance;
                values[reverse] = distance;
                comparable_columns[forward] = comparable;
                comparable_columns[reverse] = comparable;
                identical_columns[forward] = identical;
                identical_columns[reverse] = identical;
            }
        }

        Ok(Self {
            size,
            values: values.into_boxed_slice(),
            comparable_columns: comparable_columns.into_boxed_slice(),
            identical_columns: identical_columns.into_boxed_slice(),
        })
    }

    pub(crate) fn resident_bytes(&self) -> Result<usize, MafftError> {
        self.values
            .len()
            .checked_mul(size_of::<f64>())
            .and_then(|bytes| {
                self.comparable_columns
                    .len()
                    .checked_mul(size_of::<usize>())
                    .and_then(|counts| bytes.checked_add(counts))
            })
            .and_then(|bytes| {
                self.identical_columns
                    .len()
                    .checked_mul(size_of::<usize>())
                    .and_then(|counts| bytes.checked_add(counts))
            })
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-2 aligned distance memory accounting",
            })
    }

    fn offset(&self, row: usize, column: usize) -> Option<usize> {
        (row < self.size && column < self.size).then_some(row * self.size + column)
    }
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, Symbol};
    use align_msa::{AlignedRow, MultipleAlignment};

    use super::FftNs2DistanceMatrix;

    fn alignment(rows: &[(&str, &str)]) -> MultipleAlignment<Dna> {
        MultipleAlignment::try_new(
            rows.iter()
                .map(|(identifier, sequence)| {
                    AlignedRow::new(
                        *identifier,
                        sequence
                            .bytes()
                            .map(|byte| {
                                if byte == b'-' {
                                    None
                                } else {
                                    Some(Symbol::<Dna>::from_ascii(byte).unwrap())
                                }
                            })
                            .collect(),
                    )
                    .unwrap()
                })
                .collect(),
        )
        .unwrap()
    }

    #[test]
    fn hand_calculated_poisson_distances_use_pairwise_deletion() {
        let msa = alignment(&[("a", "ACGT"), ("b", "ACGA"), ("c", "A-GA")]);
        let matrix = FftNs2DistanceMatrix::from_alignment(&msa, 1.0e-12).unwrap();
        assert!((matrix.get(0, 1).unwrap() + 0.75_f64.ln()).abs() < 1.0e-12);
        assert!((matrix.get(0, 2).unwrap() + (2.0_f64 / 3.0).ln()).abs() < 1.0e-12);
        assert_eq!(matrix.get(1, 2), Some(0.0));
        assert_eq!(matrix.comparable_columns(0, 2), Some(3));
        assert_eq!(matrix.identical_columns(0, 2), Some(2));
    }

    #[test]
    fn empty_rows_and_nonoverlap_are_totalized_explicitly() {
        let both_empty = alignment(&[("a", ""), ("b", "")]);
        let empty_matrix = FftNs2DistanceMatrix::from_alignment(&both_empty, 1.0e-12).unwrap();
        assert_eq!(empty_matrix.get(0, 1), Some(0.0));

        let nonoverlap = alignment(&[("a", "A-"), ("b", "-C")]);
        let matrix = FftNs2DistanceMatrix::from_alignment(&nonoverlap, 1.0e-12).unwrap();
        assert!((matrix.get(0, 1).unwrap() + 1.0e-12_f64.ln()).abs() < 1.0e-12);
        assert_eq!(matrix.comparable_columns(0, 1), Some(0));
    }

    #[test]
    fn zero_identity_saturates_at_the_configured_finite_distance() {
        let msa = alignment(&[("a", "AAAA"), ("b", "CCCC")]);
        let matrix = FftNs2DistanceMatrix::from_alignment(&msa, 1.0e-12).unwrap();
        assert!((matrix.get(0, 1).unwrap() + 1.0e-12_f64.ln()).abs() < 1.0e-12);
        assert!(matrix.get(0, 1).unwrap().is_finite());
    }
}
