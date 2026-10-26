use std::collections::{HashMap, HashSet};

use align_core::Alphabet;

use crate::{MsaError, MultipleAlignment};

type ResiduePair = (usize, usize, usize, usize);

/// SP and TC agreement of a candidate MSA against an independent reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AlignmentMetrics {
    correct_pairs: usize,
    reference_pairs: usize,
    correct_columns: usize,
    reference_columns: usize,
}

impl AlignmentMetrics {
    /// Returns reference-normalized sum-of-pairs recall in `[0,1]`.
    #[must_use]
    pub fn sp_score(self) -> f64 {
        ratio_or_one(self.correct_pairs, self.reference_pairs)
    }

    /// Returns reference-normalized total-column score in `[0,1]`.
    #[must_use]
    pub fn tc_score(self) -> f64 {
        ratio_or_one(self.correct_columns, self.reference_columns)
    }

    /// Returns correctly recovered reference residue pairs.
    #[must_use]
    pub const fn correct_pairs(self) -> usize {
        self.correct_pairs
    }

    /// Returns all aligned residue pairs in the reference.
    #[must_use]
    pub const fn reference_pairs(self) -> usize {
        self.reference_pairs
    }

    /// Returns exactly recovered reference columns.
    #[must_use]
    pub const fn correct_columns(self) -> usize {
        self.correct_columns
    }

    /// Returns all reference columns.
    #[must_use]
    pub const fn reference_columns(self) -> usize {
        self.reference_columns
    }
}

/// Compares a candidate alignment to a reference by biological coordinates.
///
/// Row order may differ, but identifiers and ungapped biological sequences must be
/// identical. SP is the fraction of reference residue pairs recovered. TC is the
/// fraction of reference columns whose complete coordinate signature is recovered.
/// Empty denominators score one.
///
/// # Errors
///
/// Returns [`MsaError`] if the sequence sets or ungapped sequences differ.
pub fn compare_alignments<A: Alphabet>(
    candidate: &MultipleAlignment<A>,
    reference: &MultipleAlignment<A>,
) -> Result<AlignmentMetrics, MsaError> {
    let selected = vec![true; reference.width()];
    compare_alignments_on_reference_columns(candidate, reference, &selected)
}

/// Compares a candidate MSA on an explicit subset of reference columns.
///
/// This is the core-block form used by structural benchmarks such as `BAliBASE`.
/// Only selected reference columns contribute reference residue pairs and the TC
/// denominator; the complete candidate alignment remains searchable for matching
/// coordinate signatures. Empty selections score one.
///
/// # Errors
///
/// Returns [`MsaError`] if the mask length differs from the reference width, or if
/// the sequence sets or ungapped sequences differ.
pub fn compare_alignments_on_reference_columns<A: Alphabet>(
    candidate: &MultipleAlignment<A>,
    reference: &MultipleAlignment<A>,
    selected_reference_columns: &[bool],
) -> Result<AlignmentMetrics, MsaError> {
    if selected_reference_columns.len() != reference.width() {
        return Err(MsaError::ReferenceColumnMaskLengthMismatch {
            expected: reference.width(),
            actual: selected_reference_columns.len(),
        });
    }
    if candidate.sequence_count() != reference.sequence_count() {
        return Err(MsaError::SequenceSetMismatch);
    }

    let reference_indices: HashMap<&str, usize> = reference
        .rows()
        .iter()
        .enumerate()
        .map(|(index, row)| (row.identifier(), index))
        .collect();
    let mut candidate_by_reference = vec![None; reference.sequence_count()];
    for row in candidate.rows() {
        let Some(&reference_index) = reference_indices.get(row.identifier()) else {
            return Err(MsaError::SequenceSetMismatch);
        };
        if candidate_by_reference[reference_index].is_some() {
            return Err(MsaError::SequenceSetMismatch);
        }
        let reference_row = &reference.rows()[reference_index];
        if row.to_ungapped_ascii() != reference_row.to_ungapped_ascii() {
            return Err(MsaError::UngappedSequenceMismatch {
                identifier: row.identifier().to_owned(),
            });
        }
        candidate_by_reference[reference_index] = Some(row);
    }
    let mut ordered_candidate_rows = Vec::with_capacity(candidate_by_reference.len());
    for row in candidate_by_reference {
        let Some(row) = row else {
            return Err(MsaError::SequenceSetMismatch);
        };
        ordered_candidate_rows.push(row);
    }

    let reference_signatures = column_signatures_reference_order(reference)
        .into_iter()
        .zip(selected_reference_columns)
        .filter_map(|(signature, selected)| selected.then_some(signature))
        .collect::<Vec<_>>();
    let candidate_signatures = column_signatures_rows(&ordered_candidate_rows, candidate.width());

    let reference_pairs = residue_pairs(&reference_signatures);
    let candidate_pairs = residue_pairs(&candidate_signatures);
    let correct_pairs = reference_pairs.intersection(&candidate_pairs).count();
    let reference_columns: HashSet<Vec<Option<usize>>> =
        reference_signatures.iter().cloned().collect();
    let candidate_columns: HashSet<Vec<Option<usize>>> =
        candidate_signatures.iter().cloned().collect();
    let correct_columns = reference_columns.intersection(&candidate_columns).count();

    Ok(AlignmentMetrics {
        correct_pairs,
        reference_pairs: reference_pairs.len(),
        correct_columns,
        reference_columns: reference_signatures.len(),
    })
}

fn column_signatures_reference_order<A: Alphabet>(
    alignment: &MultipleAlignment<A>,
) -> Vec<Vec<Option<usize>>> {
    let rows = alignment.rows().iter().collect::<Vec<_>>();
    column_signatures_rows(&rows, alignment.width())
}

fn column_signatures_rows<A: Alphabet>(
    rows: &[&crate::AlignedRow<A>],
    width: usize,
) -> Vec<Vec<Option<usize>>> {
    let mut coordinates = vec![0_usize; rows.len()];
    let mut signatures = Vec::with_capacity(width);
    for column in 0..width {
        let mut signature = Vec::with_capacity(rows.len());
        for (row_index, row) in rows.iter().enumerate() {
            if row.symbols()[column].is_some() {
                signature.push(Some(coordinates[row_index]));
                coordinates[row_index] += 1;
            } else {
                signature.push(None);
            }
        }
        signatures.push(signature);
    }
    signatures
}

fn residue_pairs(signatures: &[Vec<Option<usize>>]) -> HashSet<ResiduePair> {
    let mut pairs = HashSet::new();
    for signature in signatures {
        for first in 0..signature.len() {
            let Some(first_residue) = signature[first] else {
                continue;
            };
            for (second, second_residue) in signature.iter().enumerate().skip(first + 1) {
                if let Some(second_residue) = second_residue {
                    pairs.insert((first, first_residue, second, *second_residue));
                }
            }
        }
    }
    pairs
}

#[allow(clippy::cast_precision_loss)]
fn ratio_or_one(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        1.0
    } else {
        numerator as f64 / denominator as f64
    }
}
