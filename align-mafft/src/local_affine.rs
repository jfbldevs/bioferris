use core::mem::size_of;
use core::ops::Range;

use align_core::Sequence;
use align_msa::SequenceRecord;

use crate::matrix::ScientificScorer;
use crate::profile::{WorkingProfile, WorkingRow};
use crate::{MafftAlphabet, MafftError};

const TRACE_RESTART: u8 = 0;
const TRACE_DIAGONAL: u8 = 1;
const TRACE_FIRST_ONLY: u8 = 2;
const TRACE_SECOND_ONLY: u8 = 3;
const TRACE_H_MASK: u8 = 0b11;
const TRACE_FIRST_EXTENDS: u8 = 1 << 2;
const TRACE_SECOND_EXTENDS: u8 = 1 << 3;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum LocalOperation {
    Diagonal,
    FirstOnly,
    SecondOnly,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct LocalAffineMemory {
    pub(crate) score_rows_bytes: usize,
    pub(crate) traceback_bytes: usize,
    pub(crate) operation_bytes: usize,
    pub(crate) profile_bytes: usize,
    pub(crate) estimated_peak_bytes: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct LocalAffinePair<A: MafftAlphabet> {
    pub(crate) profile: WorkingProfile<A>,
    pub(crate) score: f64,
    pub(crate) first_range: Range<usize>,
    pub(crate) second_range: Range<usize>,
    pub(crate) operations: Box<[LocalOperation]>,
    pub(crate) memory: LocalAffineMemory,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum TraceState {
    Best,
    FirstOnly,
    SecondOnly,
}

struct LocalDpResult {
    trace: Vec<u8>,
    columns: usize,
    cell_count: usize,
    best_score: f64,
    best_end: (usize, usize),
}

struct LocalTraceback {
    operations: Vec<LocalOperation>,
    first_range: Range<usize>,
    second_range: Range<usize>,
}

fn checked_product(first: usize, second: usize, stage: &'static str) -> Result<usize, MafftError> {
    first
        .checked_mul(second)
        .ok_or(MafftError::ArithmeticOverflow { stage })
}

fn checked_sum(values: &[usize], stage: &'static str) -> Result<usize, MafftError> {
    values
        .iter()
        .try_fold(0usize, |total, value| total.checked_add(*value))
        .ok_or(MafftError::ArithmeticOverflow { stage })
}

fn scores_agree(first: f64, second: f64) -> bool {
    let scale = first.abs().max(second.abs()).max(1.0);
    (first - second).abs() <= 1.0e-11 * scale
}

fn traceback_index(row: usize, column: usize, columns: usize) -> Result<usize, MafftError> {
    row.checked_mul(columns)
        .and_then(|start| start.checked_add(column))
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "L-INS-i local traceback index",
        })
}

fn replay_local<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
    scorer: &ScientificScorer,
    gap_open: f64,
    gap_extend: f64,
    pair: &LocalAffinePair<A>,
) -> Result<f64, MafftError> {
    let first_sequence = records[first].sequence();
    let second_sequence = records[second].sequence();
    if pair.first_range.end > first_sequence.len()
        || pair.second_range.end > second_sequence.len()
        || pair.first_range.start > pair.first_range.end
        || pair.second_range.start > pair.second_range.end
    {
        return Err(MafftError::InternalInvariant {
            message: "L-INS-i local traceback range is outside its source",
        });
    }
    if pair.score == 0.0 {
        if pair.first_range != (0..0) || pair.second_range != (0..0) || !pair.operations.is_empty()
        {
            return Err(MafftError::InternalInvariant {
                message: "L-INS-i zero-score local result is not canonical",
            });
        }
        return Ok(0.0);
    }

    let mut first_position = pair.first_range.start;
    let mut second_position = pair.second_range.start;
    let mut score = 0.0;
    let mut previous_gap = None;
    for operation in pair.operations.iter().copied() {
        match operation {
            LocalOperation::Diagonal => {
                let first_symbol =
                    first_sequence
                        .get(first_position)
                        .ok_or(MafftError::InternalInvariant {
                            message: "L-INS-i local traceback over-consumed the first sequence",
                        })?;
                let second_symbol =
                    second_sequence
                        .get(second_position)
                        .ok_or(MafftError::InternalInvariant {
                            message: "L-INS-i local traceback over-consumed the second sequence",
                        })?;
                score += scorer.score_indices(
                    first_symbol.alphabet_index(),
                    second_symbol.alphabet_index(),
                );
                first_position += 1;
                second_position += 1;
                previous_gap = None;
            }
            LocalOperation::FirstOnly => {
                if first_position >= first_sequence.len() {
                    return Err(MafftError::InternalInvariant {
                        message: "L-INS-i local gap over-consumed the first sequence",
                    });
                }
                score += if previous_gap == Some(LocalOperation::FirstOnly) {
                    gap_extend
                } else {
                    gap_open
                };
                first_position += 1;
                previous_gap = Some(LocalOperation::FirstOnly);
            }
            LocalOperation::SecondOnly => {
                if second_position >= second_sequence.len() {
                    return Err(MafftError::InternalInvariant {
                        message: "L-INS-i local gap over-consumed the second sequence",
                    });
                }
                score += if previous_gap == Some(LocalOperation::SecondOnly) {
                    gap_extend
                } else {
                    gap_open
                };
                second_position += 1;
                previous_gap = Some(LocalOperation::SecondOnly);
            }
        }
        if !score.is_finite() {
            return Err(MafftError::NonFiniteComputation {
                stage: "L-INS-i local traceback replay",
            });
        }
        if score <= 0.0 {
            return Err(MafftError::InternalInvariant {
                message: "L-INS-i positive local traceback crossed a non-positive prefix",
            });
        }
    }
    if first_position != pair.first_range.end || second_position != pair.second_range.end {
        return Err(MafftError::InternalInvariant {
            message: "L-INS-i local traceback did not consume its declared source ranges",
        });
    }
    if !scores_agree(score, pair.score) {
        return Err(MafftError::InternalInvariant {
            message: "L-INS-i local DP and independent traceback replay disagree",
        });
    }
    Ok(score)
}

fn padded_pair_profile<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
    first_range: &Range<usize>,
    second_range: &Range<usize>,
    operations: &[LocalOperation],
) -> Result<WorkingProfile<A>, MafftError> {
    let first_sequence = records[first].sequence();
    let second_sequence = records[second].sequence();
    let suffix_first =
        first_sequence
            .len()
            .checked_sub(first_range.end)
            .ok_or(MafftError::InternalInvariant {
                message: "L-INS-i local first range exceeds its source",
            })?;
    let suffix_second = second_sequence.len().checked_sub(second_range.end).ok_or(
        MafftError::InternalInvariant {
            message: "L-INS-i local second range exceeds its source",
        },
    )?;
    let width = checked_sum(
        &[
            first_range.start,
            second_range.start,
            operations.len(),
            suffix_first,
            suffix_second,
        ],
        "L-INS-i padded local profile width",
    )?;
    let mut first_cells = Vec::with_capacity(width);
    let mut second_cells = Vec::with_capacity(width);

    for position in 0..first_range.start {
        first_cells.push(first_sequence.get(position));
        second_cells.push(None);
    }
    for position in 0..second_range.start {
        first_cells.push(None);
        second_cells.push(second_sequence.get(position));
    }

    let mut first_position = first_range.start;
    let mut second_position = second_range.start;
    for operation in operations.iter().copied() {
        match operation {
            LocalOperation::Diagonal => {
                first_cells.push(first_sequence.get(first_position));
                second_cells.push(second_sequence.get(second_position));
                first_position += 1;
                second_position += 1;
            }
            LocalOperation::FirstOnly => {
                first_cells.push(first_sequence.get(first_position));
                second_cells.push(None);
                first_position += 1;
            }
            LocalOperation::SecondOnly => {
                first_cells.push(None);
                second_cells.push(second_sequence.get(second_position));
                second_position += 1;
            }
        }
    }
    if first_position != first_range.end || second_position != second_range.end {
        return Err(MafftError::InternalInvariant {
            message: "L-INS-i local operations and ranges disagree while padding flanks",
        });
    }

    for position in first_range.end..first_sequence.len() {
        first_cells.push(first_sequence.get(position));
        second_cells.push(None);
    }
    for position in second_range.end..second_sequence.len() {
        first_cells.push(None);
        second_cells.push(second_sequence.get(position));
    }
    if first_cells.len() != width || second_cells.len() != width {
        return Err(MafftError::InternalInvariant {
            message: "L-INS-i padded local pair profile is not rectangular",
        });
    }
    Ok(WorkingProfile {
        rows: vec![
            WorkingRow {
                input_index: first,
                identifier: records[first].identifier().to_owned(),
                cells: first_cells,
                weight: 1.0,
            },
            WorkingRow {
                input_index: second,
                identifier: records[second].identifier().to_owned(),
                cells: second_cells,
                weight: 1.0,
            },
        ],
        width,
    })
}

fn choose_gap(open: f64, extension: f64) -> (f64, bool) {
    if open >= extension {
        (open, false)
    } else {
        (extension, true)
    }
}

fn choose_cell(diagonal: f64, first_gap: f64, second_gap: f64) -> (f64, u8) {
    if diagonal > 0.0 && diagonal >= first_gap && diagonal >= second_gap {
        (diagonal, TRACE_DIAGONAL)
    } else if first_gap > 0.0 && first_gap >= second_gap {
        (first_gap, TRACE_FIRST_ONLY)
    } else if second_gap > 0.0 {
        (second_gap, TRACE_SECOND_ONLY)
    } else {
        (0.0, TRACE_RESTART)
    }
}

fn fill_local_dp<A: MafftAlphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    scorer: &ScientificScorer,
    gap_open: f64,
    gap_extend: f64,
) -> Result<LocalDpResult, MafftError> {
    let rows = first
        .len()
        .checked_add(1)
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "L-INS-i local DP rows",
        })?;
    let columns = second
        .len()
        .checked_add(1)
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "L-INS-i local DP columns",
        })?;
    let cell_count = checked_product(rows, columns, "L-INS-i local DP cells")?;
    let mut trace = vec![TRACE_RESTART; cell_count];
    // H only needs one row: each cell is overwritten after its previous-row
    // value has been saved as the next diagonal. The vertical affine state also
    // updates in place, while the horizontal state is a single left-to-right
    // scalar. This is the same recurrence and tie order as the former six-row
    // representation with substantially less memory traffic.
    let mut best = vec![0.0; columns];
    let mut first_gap = vec![f64::NEG_INFINITY; columns];
    let mut best_score = 0.0;
    let mut best_end = (0usize, 0usize);
    let first_indices = first.as_indices();
    let second_indices = second.as_indices();

    for (first_position, first_index) in first_indices.iter().copied().enumerate() {
        let row = first_position + 1;
        let trace_row = row * columns;
        let mut diagonal_best = best[0];
        best[0] = 0.0;
        first_gap[0] = f64::NEG_INFINITY;
        let mut second_gap = f64::NEG_INFINITY;
        let substitution_row = scorer.score_row(usize::from(first_index));

        for (second_position, second_index) in second_indices.iter().copied().enumerate() {
            let column = second_position + 1;
            let previous_best = best[column];
            let (next_first_gap, first_extends) =
                choose_gap(previous_best + gap_open, first_gap[column] + gap_extend);
            first_gap[column] = next_first_gap;
            let (next_second_gap, second_extends) =
                choose_gap(best[column - 1] + gap_open, second_gap + gap_extend);
            second_gap = next_second_gap;
            let diagonal = diagonal_best + substitution_row[usize::from(second_index)];
            diagonal_best = previous_best;
            let (cell_score, best_trace) = choose_cell(diagonal, next_first_gap, next_second_gap);
            if !cell_score.is_finite() {
                return Err(MafftError::NonFiniteComputation {
                    stage: "L-INS-i local affine DP",
                });
            }
            best[column] = cell_score;
            let extension_bits = (u8::from(first_extends) * TRACE_FIRST_EXTENDS)
                | (u8::from(second_extends) * TRACE_SECOND_EXTENDS);
            trace[trace_row + column] = best_trace | extension_bits;
            if cell_score > best_score {
                best_score = cell_score;
                best_end = (row, column);
            }
        }
    }
    Ok(LocalDpResult {
        trace,
        columns,
        cell_count,
        best_score,
        best_end,
    })
}

fn trace_local_dp(
    dp: &LocalDpResult,
    operation_capacity: usize,
) -> Result<LocalTraceback, MafftError> {
    let mut operations = Vec::with_capacity(operation_capacity);
    let mut row = dp.best_end.0;
    let mut column = dp.best_end.1;
    let mut state = TraceState::Best;
    while state != TraceState::Best
        || dp.trace[traceback_index(row, column, dp.columns)?] & TRACE_H_MASK != TRACE_RESTART
    {
        if row == 0 && state != TraceState::SecondOnly {
            return Err(MafftError::InternalInvariant {
                message: "L-INS-i local traceback crossed the first boundary",
            });
        }
        if column == 0 && state != TraceState::FirstOnly {
            return Err(MafftError::InternalInvariant {
                message: "L-INS-i local traceback crossed the second boundary",
            });
        }
        let packed = dp.trace[traceback_index(row, column, dp.columns)?];
        match state {
            TraceState::Best => match packed & TRACE_H_MASK {
                TRACE_DIAGONAL => {
                    operations.push(LocalOperation::Diagonal);
                    row -= 1;
                    column -= 1;
                }
                TRACE_FIRST_ONLY => state = TraceState::FirstOnly,
                TRACE_SECOND_ONLY => state = TraceState::SecondOnly,
                _ => {
                    return Err(MafftError::InternalInvariant {
                        message: "L-INS-i local traceback encountered an invalid best state",
                    });
                }
            },
            TraceState::FirstOnly => {
                operations.push(LocalOperation::FirstOnly);
                row -= 1;
                state = if packed & TRACE_FIRST_EXTENDS == 0 {
                    TraceState::Best
                } else {
                    TraceState::FirstOnly
                };
            }
            TraceState::SecondOnly => {
                operations.push(LocalOperation::SecondOnly);
                column -= 1;
                state = if packed & TRACE_SECOND_EXTENDS == 0 {
                    TraceState::Best
                } else {
                    TraceState::SecondOnly
                };
            }
        }
    }
    operations.reverse();
    let (first_range, second_range) = if dp.best_score == 0.0 {
        operations.clear();
        (0..0, 0..0)
    } else {
        (row..dp.best_end.0, column..dp.best_end.1)
    };
    Ok(LocalTraceback {
        operations,
        first_range,
        second_range,
    })
}

fn local_memory<A: MafftAlphabet>(
    dp: &LocalDpResult,
    operations: &Vec<LocalOperation>,
    profile: &WorkingProfile<A>,
) -> Result<LocalAffineMemory, MafftError> {
    let score_rows_bytes = checked_product(
        checked_product(
            dp.columns,
            size_of::<f64>(),
            "L-INS-i local score row bytes",
        )?,
        2,
        "L-INS-i local score rows",
    )?;
    let traceback_bytes = checked_product(
        dp.cell_count,
        size_of::<u8>(),
        "L-INS-i local traceback bytes",
    )?;
    let operation_bytes = checked_product(
        operations.capacity(),
        size_of::<LocalOperation>(),
        "L-INS-i local operation bytes",
    )?;
    let profile_cells = checked_product(profile.width, 2, "L-INS-i local profile cells")?;
    let profile_bytes = checked_product(
        profile_cells,
        size_of::<Option<align_core::Symbol<A>>>(),
        "L-INS-i local profile bytes",
    )?;
    let estimated_peak_bytes = checked_sum(
        &[
            score_rows_bytes,
            traceback_bytes,
            operation_bytes,
            profile_bytes,
        ],
        "L-INS-i local pair peak memory",
    )?;
    Ok(LocalAffineMemory {
        score_rows_bytes,
        traceback_bytes,
        operation_bytes,
        profile_bytes,
        estimated_peak_bytes,
    })
}

pub(crate) fn align_local_affine<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
    scorer: &ScientificScorer,
    gap_open: f64,
    gap_extend: f64,
) -> Result<LocalAffinePair<A>, MafftError> {
    if first >= records.len() || second >= records.len() || first == second {
        return Err(MafftError::InternalInvariant {
            message: "L-INS-i local pair indices are invalid",
        });
    }
    if !gap_open.is_finite() || !gap_extend.is_finite() || gap_open >= 0.0 || gap_extend > 0.0 {
        return Err(MafftError::NonFiniteComputation {
            stage: "L-INS-i local affine parameters",
        });
    }

    let first_sequence = records[first].sequence();
    let second_sequence = records[second].sequence();
    let dp = fill_local_dp(
        first_sequence,
        second_sequence,
        scorer,
        gap_open,
        gap_extend,
    )?;
    let operation_capacity = first_sequence
        .len()
        .checked_add(second_sequence.len())
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "L-INS-i local traceback operation capacity",
        })?;
    let traceback = trace_local_dp(&dp, operation_capacity)?;
    let profile = padded_pair_profile(
        records,
        first,
        second,
        &traceback.first_range,
        &traceback.second_range,
        &traceback.operations,
    )?;
    let memory = local_memory(&dp, &traceback.operations, &profile)?;
    let pair = LocalAffinePair {
        profile,
        score: dp.best_score,
        first_range: traceback.first_range,
        second_range: traceback.second_range,
        operations: traceback.operations.into_boxed_slice(),
        memory,
    };
    replay_local(records, first, second, scorer, gap_open, gap_extend, &pair)?;
    Ok(pair)
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Sequence};
    use align_msa::SequenceRecord;

    use super::{LocalOperation, align_local_affine};
    use crate::matrix::ScientificScorer;

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    struct PathSearch<'a> {
        first: &'a Sequence<Dna>,
        second: &'a Sequence<Dna>,
        first_end: usize,
        second_end: usize,
        scorer: &'a ScientificScorer,
        gap_open: f64,
        gap_extend: f64,
        best: &'a mut f64,
    }

    impl PathSearch<'_> {
        fn enumerate(
            &mut self,
            first_position: usize,
            second_position: usize,
            previous_gap: Option<LocalOperation>,
            score: f64,
        ) {
            if first_position == self.first_end && second_position == self.second_end {
                *self.best = (*self.best).max(score);
                return;
            }
            if first_position < self.first_end && second_position < self.second_end {
                let first_symbol = self.first.get(first_position).unwrap();
                let second_symbol = self.second.get(second_position).unwrap();
                let substitution = self.scorer.score_indices(
                    first_symbol.alphabet_index(),
                    second_symbol.alphabet_index(),
                );
                self.enumerate(
                    first_position + 1,
                    second_position + 1,
                    None,
                    score + substitution,
                );
            }
            if first_position < self.first_end {
                let operation = LocalOperation::FirstOnly;
                let gap = if previous_gap == Some(operation) {
                    self.gap_extend
                } else {
                    self.gap_open
                };
                self.enumerate(
                    first_position + 1,
                    second_position,
                    Some(operation),
                    score + gap,
                );
            }
            if second_position < self.second_end {
                let operation = LocalOperation::SecondOnly;
                let gap = if previous_gap == Some(operation) {
                    self.gap_extend
                } else {
                    self.gap_open
                };
                self.enumerate(
                    first_position,
                    second_position + 1,
                    Some(operation),
                    score + gap,
                );
            }
        }
    }

    fn exhaustive_local_score(
        first: &Sequence<Dna>,
        second: &Sequence<Dna>,
        scorer: &ScientificScorer,
        gap_open: f64,
        gap_extend: f64,
    ) -> f64 {
        let mut best = 0.0_f64;
        for first_start in 0..first.len() {
            for first_end in (first_start + 1)..=first.len() {
                for second_start in 0..second.len() {
                    for second_end in (second_start + 1)..=second.len() {
                        PathSearch {
                            first,
                            second,
                            first_end,
                            second_end,
                            scorer,
                            gap_open,
                            gap_extend,
                            best: &mut best,
                        }
                        .enumerate(first_start, second_start, None, 0.0);
                    }
                }
            }
        }
        best
    }

    fn binary_words(maximum_length: usize) -> Vec<String> {
        let mut words = vec![String::new()];
        for length in 1..=maximum_length {
            for bits in 0..(1usize << length) {
                let word = (0..length)
                    .map(|shift| if bits & (1 << shift) == 0 { 'A' } else { 'C' })
                    .collect();
                words.push(word);
            }
        }
        words
    }

    #[test]
    fn no_positive_hit_is_canonical_and_retains_both_full_sources() {
        let records = [dna("a", ""), dna("b", "AC")];
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let pair = align_local_affine(&records, 0, 1, &scorer, -2.0, -0.1).unwrap();
        assert_eq!(pair.score.to_bits(), 0.0_f64.to_bits());
        assert_eq!(pair.first_range, 0..0);
        assert_eq!(pair.second_range, 0..0);
        assert!(pair.operations.is_empty());
        assert_eq!(pair.profile.width, 2);
        assert!(pair.profile.rows[0].cells.iter().all(Option::is_none));
        assert_eq!(
            pair.profile.rows[1]
                .cells
                .iter()
                .filter_map(|cell| *cell)
                .map(align_core::Symbol::as_ascii)
                .collect::<Vec<_>>(),
            b"AC"
        );
    }

    #[test]
    fn unrelated_flanks_are_outside_the_selected_local_ranges() {
        let records = [dna("a", "TTACGTAA"), dna("b", "GGACGTCC")];
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let pair = align_local_affine(&records, 0, 1, &scorer, -2.0, -0.1).unwrap();
        assert_eq!(pair.first_range, 2..6);
        assert_eq!(pair.second_range, 2..6);
        assert_eq!(pair.operations.as_ref(), [LocalOperation::Diagonal; 4]);
        assert!(pair.score > 0.0);
        for (row, expected) in pair.profile.rows.iter().zip([b"TTACGTAA", b"GGACGTCC"]) {
            assert_eq!(
                row.cells
                    .iter()
                    .filter_map(|cell| *cell)
                    .map(align_core::Symbol::as_ascii)
                    .collect::<Vec<_>>(),
                expected
            );
        }
    }

    #[test]
    fn affine_traceback_can_retain_an_internal_gap() {
        let records = [dna("a", "ACGT"), dna("b", "ACT")];
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let pair = align_local_affine(&records, 0, 1, &scorer, -0.25, -0.05).unwrap();
        assert!(pair.operations.contains(&LocalOperation::FirstOnly));
        assert_eq!(pair.first_range, 0..4);
        assert_eq!(pair.second_range, 0..3);
    }

    #[test]
    fn repeated_runs_are_bitwise_deterministic() {
        let records = [dna("a", "AAAAAA"), dna("b", "AAAA")];
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let first = align_local_affine(&records, 0, 1, &scorer, -2.0, -0.1).unwrap();
        let second = align_local_affine(&records, 0, 1, &scorer, -2.0, -0.1).unwrap();
        assert_eq!(first.score.to_bits(), second.score.to_bits());
        assert_eq!(first.first_range, second.first_range);
        assert_eq!(first.second_range, second.second_range);
        assert_eq!(first.operations, second.operations);
    }

    #[test]
    fn local_affine_scores_match_complete_path_enumeration_through_length_four() {
        let words = binary_words(4);
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let gap_open = -0.7;
        let gap_extend = -0.2;
        for (first_index, first) in words.iter().enumerate() {
            for second in words.iter().skip(first_index) {
                let records = [dna("a", first), dna("b", second)];
                let expected = exhaustive_local_score(
                    records[0].sequence(),
                    records[1].sequence(),
                    &scorer,
                    gap_open,
                    gap_extend,
                );
                let observed =
                    align_local_affine(&records, 0, 1, &scorer, gap_open, gap_extend).unwrap();
                let scale = expected.abs().max(observed.score.abs()).max(1.0);
                assert!(
                    (expected - observed.score).abs() <= 1.0e-11 * scale,
                    "score mismatch for {first:?} versus {second:?}: expected {expected}, observed {}",
                    observed.score
                );
            }
        }
    }
}
