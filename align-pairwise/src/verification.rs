use align_core::{
    AffineGap, AlignmentOperation, Alphabet, LinearGap, PairwiseAlignment, Score, Sequence,
    SubstitutionScorer,
};

use crate::{FreeEndGaps, VerificationError};

/// Counts and independently recomputed score from alignment replay.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct VerificationReport {
    recomputed_score: Score,
    columns: usize,
    matches: usize,
    substitutions: usize,
    deletions: usize,
    insertions: usize,
}

impl VerificationReport {
    /// Returns the independently recomputed score.
    #[must_use]
    pub const fn recomputed_score(self) -> Score {
        self.recomputed_score
    }

    /// Returns the number of alignment columns.
    #[must_use]
    pub const fn columns(self) -> usize {
        self.columns
    }

    /// Returns the number of exact-symbol matches.
    #[must_use]
    pub const fn matches(self) -> usize {
        self.matches
    }

    /// Returns the number of substitutions.
    #[must_use]
    pub const fn substitutions(self) -> usize {
        self.substitutions
    }

    /// Returns the number of deletions (gaps in the second sequence).
    #[must_use]
    pub const fn deletions(self) -> usize {
        self.deletions
    }

    /// Returns the number of insertions (gaps in the first sequence).
    #[must_use]
    pub const fn insertions(self) -> usize {
        self.insertions
    }
}

/// Independently replays an affine alignment with explicit free sequence ends.
///
/// This verifier does not read score or traceback buffers. It validates the
/// boundary-to-boundary range contract implied by the configured free ends and
/// scores each contiguous deletion or insertion with one affine opening followed
/// by zero or more extensions.
///
/// # Errors
///
/// Returns [`VerificationError`] if any range, end permission, operation, consumed
/// symbol, affine run, arithmetic result, or reported score is inconsistent.
#[allow(clippy::too_many_lines)]
pub fn verify_affine_end_gap_alignment<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    alignment: &PairwiseAlignment,
    substitution: &S,
    gap: AffineGap,
    free_ends: FreeEndGaps,
) -> Result<VerificationReport, VerificationError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let first_range = alignment.first_range();
    let second_range = alignment.second_range();

    if first_range.end > first.len() {
        return Err(VerificationError::EndGapRangeOutOfBounds {
            sequence: "first",
            start: first_range.start,
            end: first_range.end,
            length: first.len(),
        });
    }
    if second_range.end > second.len() {
        return Err(VerificationError::EndGapRangeOutOfBounds {
            sequence: "second",
            start: second_range.start,
            end: second_range.end,
            length: second.len(),
        });
    }
    if !free_ends.first_start() && first_range.start != 0 {
        return Err(VerificationError::OmittedScoredStart {
            sequence: "first",
            start: first_range.start,
        });
    }
    if !free_ends.second_start() && second_range.start != 0 {
        return Err(VerificationError::OmittedScoredStart {
            sequence: "second",
            start: second_range.start,
        });
    }
    if !free_ends.first_end() && first_range.end != first.len() {
        return Err(VerificationError::OmittedScoredEnd {
            sequence: "first",
            end: first_range.end,
            length: first.len(),
        });
    }
    if !free_ends.second_end() && second_range.end != second.len() {
        return Err(VerificationError::OmittedScoredEnd {
            sequence: "second",
            end: second_range.end,
            length: second.len(),
        });
    }
    if first_range.start != 0 && second_range.start != 0 {
        return Err(VerificationError::EndGapStartNotOnBoundary {
            first_start: first_range.start,
            second_start: second_range.start,
        });
    }
    if first_range.end != first.len() && second_range.end != second.len() {
        return Err(VerificationError::EndGapEndNotOnBoundary {
            first_end: first_range.end,
            first_length: first.len(),
            second_end: second_range.end,
            second_length: second.len(),
        });
    }

    let mut first_index = first_range.start;
    let mut second_index = second_range.start;
    let mut score: Score = 0;
    let mut matches = 0;
    let mut substitutions = 0;
    let mut deletions = 0;
    let mut insertions = 0;
    let mut previous_operation = None;

    for (column, operation) in alignment.operations().iter().copied().enumerate() {
        let column_score =
            match operation {
                AlignmentOperation::Match | AlignmentOperation::Substitution => {
                    if first_index >= first_range.end {
                        return Err(VerificationError::EndGapRangeConsumptionMismatch {
                            sequence: "first",
                            expected_end: first_range.end,
                            actual_end: first_index.saturating_add(1),
                        });
                    }
                    if second_index >= second_range.end {
                        return Err(VerificationError::EndGapRangeConsumptionMismatch {
                            sequence: "second",
                            expected_end: second_range.end,
                            actual_end: second_index.saturating_add(1),
                        });
                    }
                    let first_symbol = first.get(first_index).ok_or(
                        VerificationError::EndGapRangeOutOfBounds {
                            sequence: "first",
                            start: first_range.start,
                            end: first_range.end,
                            length: first.len(),
                        },
                    )?;
                    let second_symbol = second.get(second_index).ok_or(
                        VerificationError::EndGapRangeOutOfBounds {
                            sequence: "second",
                            start: second_range.start,
                            end: second_range.end,
                            length: second.len(),
                        },
                    )?;
                    let equal = first_symbol == second_symbol;
                    if (operation == AlignmentOperation::Match) != equal {
                        return Err(VerificationError::OperationKindMismatch { column, operation });
                    }
                    first_index += 1;
                    second_index += 1;
                    if equal {
                        matches += 1;
                    } else {
                        substitutions += 1;
                    }
                    substitution.score(first_symbol, second_symbol)
                }
                AlignmentOperation::Deletion => {
                    if first_index >= first_range.end {
                        return Err(VerificationError::EndGapRangeConsumptionMismatch {
                            sequence: "first",
                            expected_end: first_range.end,
                            actual_end: first_index.saturating_add(1),
                        });
                    }
                    first_index += 1;
                    deletions += 1;
                    if previous_operation == Some(AlignmentOperation::Deletion) {
                        gap.extension_score()
                    } else {
                        gap.open_score()
                    }
                }
                AlignmentOperation::Insertion => {
                    if second_index >= second_range.end {
                        return Err(VerificationError::EndGapRangeConsumptionMismatch {
                            sequence: "second",
                            expected_end: second_range.end,
                            actual_end: second_index.saturating_add(1),
                        });
                    }
                    second_index += 1;
                    insertions += 1;
                    if previous_operation == Some(AlignmentOperation::Insertion) {
                        gap.extension_score()
                    } else {
                        gap.open_score()
                    }
                }
            };
        score = score
            .checked_add(column_score)
            .ok_or(VerificationError::ScoreOverflow { column })?;
        previous_operation = Some(operation);
    }

    if first_index != first_range.end {
        return Err(VerificationError::EndGapRangeConsumptionMismatch {
            sequence: "first",
            expected_end: first_range.end,
            actual_end: first_index,
        });
    }
    if second_index != second_range.end {
        return Err(VerificationError::EndGapRangeConsumptionMismatch {
            sequence: "second",
            expected_end: second_range.end,
            actual_end: second_index,
        });
    }
    if score != alignment.score() {
        return Err(VerificationError::ScoreMismatch {
            reported: alignment.score(),
            recomputed: score,
        });
    }

    Ok(VerificationReport {
        recomputed_score: score,
        columns: alignment.columns(),
        matches,
        substitutions,
        deletions,
        insertions,
    })
}

/// Independently replays a purported global alignment and checks its score.
///
/// This verifier does not read the dynamic-programming matrix or traceback data.
///
/// # Errors
///
/// Returns [`VerificationError`] if ranges are not global, operations disagree with
/// input symbols, replay exhausts an input, arithmetic overflows, or scores differ.
#[allow(clippy::too_many_lines)]
pub fn verify_global_alignment<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    alignment: &PairwiseAlignment,
    substitution: &S,
    gap: LinearGap,
) -> Result<VerificationReport, VerificationError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let first_range = alignment.first_range();
    if first_range.start != 0 || first_range.end != first.len() {
        return Err(VerificationError::NonGlobalRange {
            sequence: "first",
            start: first_range.start,
            end: first_range.end,
            length: first.len(),
        });
    }
    let second_range = alignment.second_range();
    if second_range.start != 0 || second_range.end != second.len() {
        return Err(VerificationError::NonGlobalRange {
            sequence: "second",
            start: second_range.start,
            end: second_range.end,
            length: second.len(),
        });
    }

    let mut first_index = 0;
    let mut second_index = 0;
    let mut score: Score = 0;
    let mut matches = 0;
    let mut substitutions = 0;
    let mut deletions = 0;
    let mut insertions = 0;

    for (column, operation) in alignment.operations().iter().copied().enumerate() {
        let column_score = match operation {
            AlignmentOperation::Match | AlignmentOperation::Substitution => {
                let first_symbol =
                    first
                        .get(first_index)
                        .ok_or(VerificationError::SequenceExhausted {
                            sequence: "first",
                            column,
                        })?;
                let second_symbol =
                    second
                        .get(second_index)
                        .ok_or(VerificationError::SequenceExhausted {
                            sequence: "second",
                            column,
                        })?;
                let equal = first_symbol == second_symbol;
                if (operation == AlignmentOperation::Match) != equal {
                    return Err(VerificationError::OperationKindMismatch { column, operation });
                }
                first_index += 1;
                second_index += 1;
                if equal {
                    matches += 1;
                } else {
                    substitutions += 1;
                }
                substitution.score(first_symbol, second_symbol)
            }
            AlignmentOperation::Deletion => {
                if first.get(first_index).is_none() {
                    return Err(VerificationError::SequenceExhausted {
                        sequence: "first",
                        column,
                    });
                }
                first_index += 1;
                deletions += 1;
                gap.score()
            }
            AlignmentOperation::Insertion => {
                if second.get(second_index).is_none() {
                    return Err(VerificationError::SequenceExhausted {
                        sequence: "second",
                        column,
                    });
                }
                second_index += 1;
                insertions += 1;
                gap.score()
            }
        };
        score = score
            .checked_add(column_score)
            .ok_or(VerificationError::ScoreOverflow { column })?;
    }

    if first_index != first.len() {
        return Err(VerificationError::SequenceExhausted {
            sequence: "first",
            column: alignment.columns(),
        });
    }
    if second_index != second.len() {
        return Err(VerificationError::SequenceExhausted {
            sequence: "second",
            column: alignment.columns(),
        });
    }
    if score != alignment.score() {
        return Err(VerificationError::ScoreMismatch {
            reported: alignment.score(),
            recomputed: score,
        });
    }

    Ok(VerificationReport {
        recomputed_score: score,
        columns: alignment.columns(),
        matches,
        substitutions,
        deletions,
        insertions,
    })
}

/// Independently replays a global alignment using affine gap-run scores.
///
/// This verifier does not read Gotoh score or traceback buffers. A deletion or
/// insertion following the same operation extends that gap; every other transition
/// opens a new gap.
///
/// # Errors
///
/// Returns [`VerificationError`] if ranges are not global, operations disagree with
/// input symbols, replay exhausts an input, arithmetic overflows, or scores differ.
#[allow(clippy::too_many_lines)]
pub fn verify_global_affine_alignment<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    alignment: &PairwiseAlignment,
    substitution: &S,
    gap: AffineGap,
) -> Result<VerificationReport, VerificationError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let first_range = alignment.first_range();
    if first_range.start != 0 || first_range.end != first.len() {
        return Err(VerificationError::NonGlobalRange {
            sequence: "first",
            start: first_range.start,
            end: first_range.end,
            length: first.len(),
        });
    }
    let second_range = alignment.second_range();
    if second_range.start != 0 || second_range.end != second.len() {
        return Err(VerificationError::NonGlobalRange {
            sequence: "second",
            start: second_range.start,
            end: second_range.end,
            length: second.len(),
        });
    }

    let mut first_index = 0;
    let mut second_index = 0;
    let mut score: Score = 0;
    let mut matches = 0;
    let mut substitutions = 0;
    let mut deletions = 0;
    let mut insertions = 0;
    let mut previous_operation = None;

    for (column, operation) in alignment.operations().iter().copied().enumerate() {
        let column_score = match operation {
            AlignmentOperation::Match | AlignmentOperation::Substitution => {
                let first_symbol =
                    first
                        .get(first_index)
                        .ok_or(VerificationError::SequenceExhausted {
                            sequence: "first",
                            column,
                        })?;
                let second_symbol =
                    second
                        .get(second_index)
                        .ok_or(VerificationError::SequenceExhausted {
                            sequence: "second",
                            column,
                        })?;
                let equal = first_symbol == second_symbol;
                if (operation == AlignmentOperation::Match) != equal {
                    return Err(VerificationError::OperationKindMismatch { column, operation });
                }
                first_index += 1;
                second_index += 1;
                if equal {
                    matches += 1;
                } else {
                    substitutions += 1;
                }
                substitution.score(first_symbol, second_symbol)
            }
            AlignmentOperation::Deletion => {
                if first.get(first_index).is_none() {
                    return Err(VerificationError::SequenceExhausted {
                        sequence: "first",
                        column,
                    });
                }
                first_index += 1;
                deletions += 1;
                if previous_operation == Some(AlignmentOperation::Deletion) {
                    gap.extension_score()
                } else {
                    gap.open_score()
                }
            }
            AlignmentOperation::Insertion => {
                if second.get(second_index).is_none() {
                    return Err(VerificationError::SequenceExhausted {
                        sequence: "second",
                        column,
                    });
                }
                second_index += 1;
                insertions += 1;
                if previous_operation == Some(AlignmentOperation::Insertion) {
                    gap.extension_score()
                } else {
                    gap.open_score()
                }
            }
        };
        score = score
            .checked_add(column_score)
            .ok_or(VerificationError::ScoreOverflow { column })?;
        previous_operation = Some(operation);
    }

    if first_index != first.len() {
        return Err(VerificationError::SequenceExhausted {
            sequence: "first",
            column: alignment.columns(),
        });
    }
    if second_index != second.len() {
        return Err(VerificationError::SequenceExhausted {
            sequence: "second",
            column: alignment.columns(),
        });
    }
    if score != alignment.score() {
        return Err(VerificationError::ScoreMismatch {
            reported: alignment.score(),
            recomputed: score,
        });
    }

    Ok(VerificationReport {
        recomputed_score: score,
        columns: alignment.columns(),
        matches,
        substitutions,
        deletions,
        insertions,
    })
}

/// Independently replays a purported local alignment and checks its ranges and score.
///
/// This verifier does not read dynamic-programming or traceback buffers. It also
/// checks that every prefix of a positive local traceback remains positive, matching
/// traceback termination at zero.
///
/// # Errors
///
/// Returns [`VerificationError`] if ranges are out of bounds, the canonical empty
/// representation is violated, operations disagree with source symbols, replay does
/// not consume the declared ranges, arithmetic overflows, a positive traceback crosses
/// zero, or reported and recomputed scores differ.
#[allow(clippy::too_many_lines)]
pub fn verify_local_alignment<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    alignment: &PairwiseAlignment,
    substitution: &S,
    gap: LinearGap,
) -> Result<VerificationReport, VerificationError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let first_range = alignment.first_range();
    if first_range.end > first.len() {
        return Err(VerificationError::LocalRangeOutOfBounds {
            sequence: "first",
            start: first_range.start,
            end: first_range.end,
            length: first.len(),
        });
    }
    let second_range = alignment.second_range();
    if second_range.end > second.len() {
        return Err(VerificationError::LocalRangeOutOfBounds {
            sequence: "second",
            start: second_range.start,
            end: second_range.end,
            length: second.len(),
        });
    }
    if alignment.score() < 0 {
        return Err(VerificationError::NegativeLocalScore {
            score: alignment.score(),
        });
    }
    if alignment.score() == 0
        && (!alignment.operations().is_empty() || first_range != (0..0) || second_range != (0..0))
    {
        return Err(VerificationError::NonCanonicalEmptyLocalAlignment);
    }

    let mut first_index = first_range.start;
    let mut second_index = second_range.start;
    let mut score: Score = 0;
    let mut matches = 0;
    let mut substitutions = 0;
    let mut deletions = 0;
    let mut insertions = 0;

    for (column, operation) in alignment.operations().iter().copied().enumerate() {
        let column_score = match operation {
            AlignmentOperation::Match | AlignmentOperation::Substitution => {
                if first_index >= first_range.end {
                    return Err(VerificationError::LocalRangeConsumptionMismatch {
                        sequence: "first",
                        expected_end: first_range.end,
                        actual_end: first_index.saturating_add(1),
                    });
                }
                if second_index >= second_range.end {
                    return Err(VerificationError::LocalRangeConsumptionMismatch {
                        sequence: "second",
                        expected_end: second_range.end,
                        actual_end: second_index.saturating_add(1),
                    });
                }
                let first_symbol =
                    first
                        .get(first_index)
                        .ok_or(VerificationError::LocalRangeOutOfBounds {
                            sequence: "first",
                            start: first_range.start,
                            end: first_range.end,
                            length: first.len(),
                        })?;
                let second_symbol =
                    second
                        .get(second_index)
                        .ok_or(VerificationError::LocalRangeOutOfBounds {
                            sequence: "second",
                            start: second_range.start,
                            end: second_range.end,
                            length: second.len(),
                        })?;
                let equal = first_symbol == second_symbol;
                if (operation == AlignmentOperation::Match) != equal {
                    return Err(VerificationError::OperationKindMismatch { column, operation });
                }
                first_index += 1;
                second_index += 1;
                if equal {
                    matches += 1;
                } else {
                    substitutions += 1;
                }
                substitution.score(first_symbol, second_symbol)
            }
            AlignmentOperation::Deletion => {
                if first_index >= first_range.end {
                    return Err(VerificationError::LocalRangeConsumptionMismatch {
                        sequence: "first",
                        expected_end: first_range.end,
                        actual_end: first_index.saturating_add(1),
                    });
                }
                first_index += 1;
                deletions += 1;
                gap.score()
            }
            AlignmentOperation::Insertion => {
                if second_index >= second_range.end {
                    return Err(VerificationError::LocalRangeConsumptionMismatch {
                        sequence: "second",
                        expected_end: second_range.end,
                        actual_end: second_index.saturating_add(1),
                    });
                }
                second_index += 1;
                insertions += 1;
                gap.score()
            }
        };
        score = score
            .checked_add(column_score)
            .ok_or(VerificationError::ScoreOverflow { column })?;
        if alignment.score() > 0 && score <= 0 {
            return Err(VerificationError::NonPositiveLocalPrefix { column, score });
        }
    }

    if first_index != first_range.end {
        return Err(VerificationError::LocalRangeConsumptionMismatch {
            sequence: "first",
            expected_end: first_range.end,
            actual_end: first_index,
        });
    }
    if second_index != second_range.end {
        return Err(VerificationError::LocalRangeConsumptionMismatch {
            sequence: "second",
            expected_end: second_range.end,
            actual_end: second_index,
        });
    }
    if score != alignment.score() {
        return Err(VerificationError::ScoreMismatch {
            reported: alignment.score(),
            recomputed: score,
        });
    }

    Ok(VerificationReport {
        recomputed_score: score,
        columns: alignment.columns(),
        matches,
        substitutions,
        deletions,
        insertions,
    })
}

/// Independently replays an affine local alignment and checks ranges and score.
///
/// A positive traceback must retain a positive score after every operation,
/// matching Smith-Waterman restart at zero. Contiguous equal-direction gaps are
/// replayed as one opening followed by extensions.
///
/// # Errors
///
/// Returns [`VerificationError`] if ranges are invalid, the zero result is not
/// canonical, operations do not consume the declared ranges, a positive prefix
/// crosses zero, arithmetic overflows, or the score differs.
#[allow(clippy::too_many_lines)]
pub fn verify_local_affine_alignment<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    alignment: &PairwiseAlignment,
    substitution: &S,
    gap: AffineGap,
) -> Result<VerificationReport, VerificationError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let first_range = alignment.first_range();
    if first_range.end > first.len() {
        return Err(VerificationError::LocalRangeOutOfBounds {
            sequence: "first",
            start: first_range.start,
            end: first_range.end,
            length: first.len(),
        });
    }
    let second_range = alignment.second_range();
    if second_range.end > second.len() {
        return Err(VerificationError::LocalRangeOutOfBounds {
            sequence: "second",
            start: second_range.start,
            end: second_range.end,
            length: second.len(),
        });
    }
    if alignment.score() < 0 {
        return Err(VerificationError::NegativeLocalScore {
            score: alignment.score(),
        });
    }
    if alignment.score() == 0
        && (!alignment.operations().is_empty() || first_range != (0..0) || second_range != (0..0))
    {
        return Err(VerificationError::NonCanonicalEmptyLocalAlignment);
    }

    let mut first_index = first_range.start;
    let mut second_index = second_range.start;
    let mut score: Score = 0;
    let mut matches = 0;
    let mut substitutions = 0;
    let mut deletions = 0;
    let mut insertions = 0;
    let mut previous_operation = None;

    for (column, operation) in alignment.operations().iter().copied().enumerate() {
        let column_score = match operation {
            AlignmentOperation::Match | AlignmentOperation::Substitution => {
                if first_index >= first_range.end {
                    return Err(VerificationError::LocalRangeConsumptionMismatch {
                        sequence: "first",
                        expected_end: first_range.end,
                        actual_end: first_index.saturating_add(1),
                    });
                }
                if second_index >= second_range.end {
                    return Err(VerificationError::LocalRangeConsumptionMismatch {
                        sequence: "second",
                        expected_end: second_range.end,
                        actual_end: second_index.saturating_add(1),
                    });
                }
                let first_symbol =
                    first
                        .get(first_index)
                        .ok_or(VerificationError::LocalRangeOutOfBounds {
                            sequence: "first",
                            start: first_range.start,
                            end: first_range.end,
                            length: first.len(),
                        })?;
                let second_symbol =
                    second
                        .get(second_index)
                        .ok_or(VerificationError::LocalRangeOutOfBounds {
                            sequence: "second",
                            start: second_range.start,
                            end: second_range.end,
                            length: second.len(),
                        })?;
                let equal = first_symbol == second_symbol;
                if (operation == AlignmentOperation::Match) != equal {
                    return Err(VerificationError::OperationKindMismatch { column, operation });
                }
                first_index += 1;
                second_index += 1;
                if equal {
                    matches += 1;
                } else {
                    substitutions += 1;
                }
                substitution.score(first_symbol, second_symbol)
            }
            AlignmentOperation::Deletion => {
                if first_index >= first_range.end {
                    return Err(VerificationError::LocalRangeConsumptionMismatch {
                        sequence: "first",
                        expected_end: first_range.end,
                        actual_end: first_index.saturating_add(1),
                    });
                }
                first_index += 1;
                deletions += 1;
                if previous_operation == Some(AlignmentOperation::Deletion) {
                    gap.extension_score()
                } else {
                    gap.open_score()
                }
            }
            AlignmentOperation::Insertion => {
                if second_index >= second_range.end {
                    return Err(VerificationError::LocalRangeConsumptionMismatch {
                        sequence: "second",
                        expected_end: second_range.end,
                        actual_end: second_index.saturating_add(1),
                    });
                }
                second_index += 1;
                insertions += 1;
                if previous_operation == Some(AlignmentOperation::Insertion) {
                    gap.extension_score()
                } else {
                    gap.open_score()
                }
            }
        };
        score = score
            .checked_add(column_score)
            .ok_or(VerificationError::ScoreOverflow { column })?;
        if alignment.score() > 0 && score <= 0 {
            return Err(VerificationError::NonPositiveLocalPrefix { column, score });
        }
        previous_operation = Some(operation);
    }

    if first_index != first_range.end {
        return Err(VerificationError::LocalRangeConsumptionMismatch {
            sequence: "first",
            expected_end: first_range.end,
            actual_end: first_index,
        });
    }
    if second_index != second_range.end {
        return Err(VerificationError::LocalRangeConsumptionMismatch {
            sequence: "second",
            expected_end: second_range.end,
            actual_end: second_index,
        });
    }
    if score != alignment.score() {
        return Err(VerificationError::ScoreMismatch {
            reported: alignment.score(),
            recomputed: score,
        });
    }

    Ok(VerificationReport {
        recomputed_score: score,
        columns: alignment.columns(),
        matches,
        substitutions,
        deletions,
        insertions,
    })
}

#[cfg(test)]
mod tests {
    use align_core::{
        AlignmentOperation, DnaSequence, LinearGap, MatchMismatch, PairwiseAlignment,
    };

    use super::{verify_global_affine_alignment, verify_global_alignment, verify_local_alignment};
    use crate::VerificationError;

    #[test]
    fn verifier_rejects_incorrect_reported_score() {
        let first = DnaSequence::from_ascii(b"A").expect("valid DNA");
        let second = DnaSequence::from_ascii(b"A").expect("valid DNA");
        let alignment = PairwiseAlignment::try_new(99, 0..1, 0..1, vec![AlignmentOperation::Match])
            .expect("structurally valid");
        let result = verify_global_alignment(
            &first,
            &second,
            &alignment,
            &MatchMismatch::new(2, -1),
            LinearGap::new(-2).expect("valid gap"),
        );
        assert_eq!(
            result,
            Err(VerificationError::ScoreMismatch {
                reported: 99,
                recomputed: 2,
            })
        );
    }

    #[test]
    fn verifier_rejects_wrong_operation_kind() {
        let first = DnaSequence::from_ascii(b"A").expect("valid DNA");
        let second = DnaSequence::from_ascii(b"C").expect("valid DNA");
        let alignment = PairwiseAlignment::try_new(-1, 0..1, 0..1, vec![AlignmentOperation::Match])
            .expect("structurally valid");
        let result = verify_global_alignment(
            &first,
            &second,
            &alignment,
            &MatchMismatch::new(2, -1),
            LinearGap::new(-2).expect("valid gap"),
        );
        assert!(matches!(
            result,
            Err(VerificationError::OperationKindMismatch { column: 0, .. })
        ));
    }

    #[test]
    fn local_verifier_rejects_noncanonical_zero_alignment() {
        let first = DnaSequence::from_ascii(b"A").expect("valid DNA");
        let second = DnaSequence::from_ascii(b"C").expect("valid DNA");
        let alignment =
            PairwiseAlignment::try_new(0, 0..1, 0..1, vec![AlignmentOperation::Substitution])
                .expect("structurally valid");
        let result = verify_local_alignment(
            &first,
            &second,
            &alignment,
            &MatchMismatch::new(1, 0),
            LinearGap::new(-1).expect("valid gap"),
        );
        assert_eq!(
            result,
            Err(VerificationError::NonCanonicalEmptyLocalAlignment)
        );
    }

    #[test]
    fn local_verifier_rejects_traceback_crossing_zero() {
        let first = DnaSequence::from_ascii(b"CA").expect("valid DNA");
        let second = DnaSequence::from_ascii(b"TA").expect("valid DNA");
        let alignment = PairwiseAlignment::try_new(
            1,
            0..2,
            0..2,
            vec![AlignmentOperation::Substitution, AlignmentOperation::Match],
        )
        .expect("structurally valid");
        let result = verify_local_alignment(
            &first,
            &second,
            &alignment,
            &MatchMismatch::new(2, -1),
            LinearGap::new(-2).expect("valid gap"),
        );
        assert_eq!(
            result,
            Err(VerificationError::NonPositiveLocalPrefix {
                column: 0,
                score: -1,
            })
        );
    }

    #[test]
    fn local_verifier_rejects_a_range_outside_the_supplied_sequence() {
        let first = DnaSequence::from_ascii(b"A").expect("valid DNA");
        let second = DnaSequence::from_ascii(b"AA").expect("valid DNA");
        let alignment = PairwiseAlignment::try_new(
            4,
            0..2,
            0..2,
            vec![AlignmentOperation::Match, AlignmentOperation::Match],
        )
        .expect("structurally valid without source sequences");
        let result = verify_local_alignment(
            &first,
            &second,
            &alignment,
            &MatchMismatch::new(2, -1),
            LinearGap::new(-2).expect("valid gap"),
        );
        assert_eq!(
            result,
            Err(VerificationError::LocalRangeOutOfBounds {
                sequence: "first",
                start: 0,
                end: 2,
                length: 1,
            })
        );
    }

    #[test]
    fn affine_verifier_opens_and_extends_each_contiguous_gap_run() {
        let first = DnaSequence::from_ascii(b"ACGT").expect("valid DNA");
        let second = DnaSequence::from_ascii(b"AT").expect("valid DNA");
        let alignment = PairwiseAlignment::try_new(
            -1,
            0..4,
            0..2,
            vec![
                AlignmentOperation::Match,
                AlignmentOperation::Deletion,
                AlignmentOperation::Deletion,
                AlignmentOperation::Match,
            ],
        )
        .expect("structurally valid");
        let report = verify_global_affine_alignment(
            &first,
            &second,
            &alignment,
            &MatchMismatch::new(2, -3),
            align_core::AffineGap::new(-4, -1).expect("valid gap"),
        )
        .expect("score replays");
        assert_eq!(report.recomputed_score(), -1);
        assert_eq!(report.deletions(), 2);
    }
}
