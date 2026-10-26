use align_core::{AlignmentOperation, Alphabet, PairwiseAlignment, Sequence};

use crate::UnitEditVerificationError;

/// Counts independently replayed operations under unit edit costs.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct UnitEditVerificationReport {
    distance: usize,
    matches: usize,
    substitutions: usize,
    deletions: usize,
    insertions: usize,
}

impl UnitEditVerificationReport {
    /// Returns the independently counted unit edit distance.
    #[must_use]
    pub const fn distance(self) -> usize {
        self.distance
    }

    /// Returns the number of equal diagonal columns.
    #[must_use]
    pub const fn matches(self) -> usize {
        self.matches
    }

    /// Returns the number of unequal diagonal columns.
    #[must_use]
    pub const fn substitutions(self) -> usize {
        self.substitutions
    }

    /// Returns the number of pattern symbols aligned to gaps.
    #[must_use]
    pub const fn deletions(self) -> usize {
        self.deletions
    }

    /// Returns the number of text symbols aligned to gaps.
    #[must_use]
    pub const fn insertions(self) -> usize {
        self.insertions
    }
}

/// Independently replays and validates a global unit-edit alignment.
///
/// # Errors
///
/// Returns [`UnitEditVerificationError`] when ranges are not global, an operation
/// disagrees with the supplied symbols, a sequence is exhausted, the distance is
/// not representable as a negated score, or the reported score differs.
pub fn verify_unit_edit_alignment<A: Alphabet>(
    pattern: &Sequence<A>,
    text: &Sequence<A>,
    alignment: &PairwiseAlignment,
) -> Result<UnitEditVerificationReport, UnitEditVerificationError> {
    let pattern_range = alignment.first_range();
    if pattern_range != (0..pattern.len()) {
        return Err(UnitEditVerificationError::NonGlobalRange {
            sequence: "pattern",
            start: pattern_range.start,
            end: pattern_range.end,
            length: pattern.len(),
        });
    }
    let text_range = alignment.second_range();
    if text_range != (0..text.len()) {
        return Err(UnitEditVerificationError::NonGlobalRange {
            sequence: "text",
            start: text_range.start,
            end: text_range.end,
            length: text.len(),
        });
    }

    let mut pattern_index = 0;
    let mut text_index = 0;
    let mut report = UnitEditVerificationReport {
        distance: 0,
        matches: 0,
        substitutions: 0,
        deletions: 0,
        insertions: 0,
    };

    for (column, operation) in alignment.operations().iter().copied().enumerate() {
        match operation {
            AlignmentOperation::Match | AlignmentOperation::Substitution => {
                let pattern_symbol = pattern.get(pattern_index).ok_or(
                    UnitEditVerificationError::SequenceExhausted {
                        sequence: "pattern",
                        column,
                    },
                )?;
                let text_symbol =
                    text.get(text_index)
                        .ok_or(UnitEditVerificationError::SequenceExhausted {
                            sequence: "text",
                            column,
                        })?;
                let equal = pattern_symbol == text_symbol;
                if equal != matches!(operation, AlignmentOperation::Match) {
                    return Err(UnitEditVerificationError::OperationKindMismatch {
                        column,
                        operation,
                    });
                }
                pattern_index += 1;
                text_index += 1;
                if equal {
                    report.matches += 1;
                } else {
                    report.substitutions += 1;
                    report.distance += 1;
                }
            }
            AlignmentOperation::Deletion => {
                pattern
                    .get(pattern_index)
                    .ok_or(UnitEditVerificationError::SequenceExhausted {
                        sequence: "pattern",
                        column,
                    })?;
                pattern_index += 1;
                report.deletions += 1;
                report.distance += 1;
            }
            AlignmentOperation::Insertion => {
                text.get(text_index)
                    .ok_or(UnitEditVerificationError::SequenceExhausted {
                        sequence: "text",
                        column,
                    })?;
                text_index += 1;
                report.insertions += 1;
                report.distance += 1;
            }
        }
    }

    let recomputed = i64::try_from(report.distance)
        .ok()
        .and_then(i64::checked_neg)
        .ok_or(UnitEditVerificationError::ScoreOverflow {
            distance: report.distance,
        })?;
    if alignment.score() != recomputed {
        return Err(UnitEditVerificationError::ScoreMismatch {
            reported: alignment.score(),
            recomputed,
        });
    }
    Ok(report)
}
