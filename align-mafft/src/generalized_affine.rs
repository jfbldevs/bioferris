use core::mem::size_of;
use core::ops::Range;

use align_core::Sequence;
use align_msa::SequenceRecord;

use crate::matrix::ScientificScorer;
use crate::profile::{WorkingProfile, WorkingRow};
use crate::{MafftAlphabet, MafftError};

const H_RESTART: u16 = 0;
const H_FROM_H: u16 = 1;
const H_FROM_X: u16 = 2;
const H_FROM_Y: u16 = 3;
const H_FROM_Z: u16 = 4;
const H_MASK: u16 = 0b111;
const X_EXTENDS: u16 = 1 << 3;
const Y_EXTENDS: u16 = 1 << 4;
const Z_SHIFT: u32 = 5;
const Z_MASK: u16 = 0b111 << Z_SHIFT;
const Z_OPEN_PAIR: u16 = 1;
const Z_OPEN_FIRST: u16 = 2;
const Z_OPEN_SECOND: u16 = 3;
const Z_EXTEND_PAIR: u16 = 4;
const Z_EXTEND_FIRST: u16 = 5;
const Z_EXTEND_SECOND: u16 = 6;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum GeneralizedOperation {
    Aligned,
    OrdinaryFirstOnly,
    OrdinarySecondOnly,
    SkipPair,
    SkipFirstOnly,
    SkipSecondOnly,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct GeneralizedAffineMemory {
    pub(crate) score_rows_bytes: usize,
    pub(crate) traceback_bytes: usize,
    pub(crate) operation_bytes: usize,
    pub(crate) profile_bytes: usize,
    pub(crate) estimated_peak_bytes: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct GeneralizedAffinePair<A: MafftAlphabet> {
    pub(crate) profile: WorkingProfile<A>,
    pub(crate) score: f64,
    pub(crate) first_range: Range<usize>,
    pub(crate) second_range: Range<usize>,
    pub(crate) operations: Box<[GeneralizedOperation]>,
    pub(crate) skip_opening_count: usize,
    pub(crate) skipped_pair_count: usize,
    pub(crate) skipped_unpaired_residue_count: usize,
    pub(crate) memory: GeneralizedAffineMemory,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum TraceState {
    Aligned,
    OrdinaryFirstOnly,
    OrdinarySecondOnly,
    Skip,
}

struct GeneralizedDpResult {
    trace: Vec<u16>,
    columns: usize,
    cell_count: usize,
    best_score: f64,
    best_end: (usize, usize),
}

struct GeneralizedTraceback {
    operations: Vec<GeneralizedOperation>,
    first_range: Range<usize>,
    second_range: Range<usize>,
}

#[derive(Clone, Copy)]
pub(crate) struct GeneralizedParameters {
    pub(crate) ordinary_open: f64,
    pub(crate) ordinary_extension: f64,
    pub(crate) skip_open: f64,
    pub(crate) paired_skip_extension: f64,
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
            stage: "E-INS-i generalized traceback index",
        })
}

fn choose_gap(open: f64, extension: f64) -> (f64, bool) {
    if open >= extension {
        (open, false)
    } else {
        (extension, true)
    }
}

fn choose_h_predecessor(h: f64, x: f64, y: f64, z: f64) -> (f64, u16) {
    if h > 0.0 && h >= x && h >= y && h >= z {
        (h, H_FROM_H)
    } else if x > 0.0 && x >= y && x >= z {
        (x, H_FROM_X)
    } else if y > 0.0 && y >= z {
        (y, H_FROM_Y)
    } else if z > 0.0 {
        (z, H_FROM_Z)
    } else {
        (0.0, H_RESTART)
    }
}

fn choose_skip(candidates: &[(f64, u16); 6]) -> (f64, u16) {
    let mut best = candidates[0];
    for candidate in candidates.iter().copied().skip(1) {
        if candidate.0 > best.0 {
            best = candidate;
        }
    }
    best
}

#[allow(clippy::too_many_lines)]
fn fill_generalized_dp<A: MafftAlphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    scorer: &ScientificScorer,
    parameters: GeneralizedParameters,
) -> Result<GeneralizedDpResult, MafftError> {
    let rows = first
        .len()
        .checked_add(1)
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "E-INS-i generalized DP rows",
        })?;
    let columns = second
        .len()
        .checked_add(1)
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "E-INS-i generalized DP columns",
        })?;
    let cell_count = checked_product(rows, columns, "E-INS-i generalized DP cells")?;
    let mut trace = vec![0u16; cell_count];
    let mut previous_h = vec![0.0; columns];
    let mut current_h = vec![0.0; columns];
    let mut previous_x = vec![f64::NEG_INFINITY; columns];
    let mut current_x = vec![f64::NEG_INFINITY; columns];
    let mut previous_y = vec![f64::NEG_INFINITY; columns];
    let mut current_y = vec![f64::NEG_INFINITY; columns];
    let mut previous_z = vec![f64::NEG_INFINITY; columns];
    let mut current_z = vec![f64::NEG_INFINITY; columns];
    let mut best_score = 0.0;
    let mut best_end = (0usize, 0usize);

    for row in 1..rows {
        current_h[0] = 0.0;
        current_x[0] = f64::NEG_INFINITY;
        current_y[0] = f64::NEG_INFINITY;
        current_z[0] = f64::NEG_INFINITY;
        for column in 1..columns {
            let (x, x_extends) = choose_gap(
                previous_h[column] + parameters.ordinary_open,
                previous_x[column] + parameters.ordinary_extension,
            );
            current_x[column] = x;
            let (y, y_extends) = choose_gap(
                current_h[column - 1] + parameters.ordinary_open,
                current_y[column - 1] + parameters.ordinary_extension,
            );
            current_y[column] = y;

            let (predecessor, h_trace) = choose_h_predecessor(
                previous_h[column - 1],
                previous_x[column - 1],
                previous_y[column - 1],
                previous_z[column - 1],
            );
            let first_symbol = first.get(row - 1).ok_or(MafftError::InternalInvariant {
                message: "E-INS-i generalized DP lost a first-sequence symbol",
            })?;
            let second_symbol = second
                .get(column - 1)
                .ok_or(MafftError::InternalInvariant {
                    message: "E-INS-i generalized DP lost a second-sequence symbol",
                })?;
            let aligned = predecessor
                + scorer.score_indices(
                    first_symbol.alphabet_index(),
                    second_symbol.alphabet_index(),
                );
            let (h, stored_h_trace) = if aligned > 0.0 {
                (aligned, h_trace)
            } else {
                (0.0, H_RESTART)
            };
            current_h[column] = h;

            let (z, z_trace) = choose_skip(&[
                (
                    previous_h[column - 1]
                        + parameters.skip_open
                        + parameters.paired_skip_extension,
                    Z_OPEN_PAIR,
                ),
                (
                    previous_h[column] + parameters.skip_open + parameters.ordinary_extension,
                    Z_OPEN_FIRST,
                ),
                (
                    current_h[column - 1] + parameters.skip_open + parameters.ordinary_extension,
                    Z_OPEN_SECOND,
                ),
                (
                    previous_z[column - 1] + parameters.paired_skip_extension,
                    Z_EXTEND_PAIR,
                ),
                (
                    previous_z[column] + parameters.ordinary_extension,
                    Z_EXTEND_FIRST,
                ),
                (
                    current_z[column - 1] + parameters.ordinary_extension,
                    Z_EXTEND_SECOND,
                ),
            ]);
            current_z[column] = z;

            if !h.is_finite() {
                return Err(MafftError::NonFiniteComputation {
                    stage: "E-INS-i generalized-affine DP",
                });
            }
            let packed = stored_h_trace
                | (u16::from(x_extends) * X_EXTENDS)
                | (u16::from(y_extends) * Y_EXTENDS)
                | (z_trace << Z_SHIFT);
            trace[traceback_index(row, column, columns)?] = packed;
            if h > best_score {
                best_score = h;
                best_end = (row, column);
            }
        }
        core::mem::swap(&mut previous_h, &mut current_h);
        core::mem::swap(&mut previous_x, &mut current_x);
        core::mem::swap(&mut previous_y, &mut current_y);
        core::mem::swap(&mut previous_z, &mut current_z);
    }

    Ok(GeneralizedDpResult {
        trace,
        columns,
        cell_count,
        best_score,
        best_end,
    })
}

#[allow(clippy::too_many_lines)]
fn trace_generalized_dp(
    dp: &GeneralizedDpResult,
    operation_capacity: usize,
) -> Result<GeneralizedTraceback, MafftError> {
    if dp.best_score == 0.0 {
        return Ok(GeneralizedTraceback {
            operations: Vec::new(),
            first_range: 0..0,
            second_range: 0..0,
        });
    }

    let mut operations = Vec::with_capacity(operation_capacity);
    let mut row = dp.best_end.0;
    let mut column = dp.best_end.1;
    let mut state = TraceState::Aligned;
    loop {
        let packed = dp.trace[traceback_index(row, column, dp.columns)?];
        match state {
            TraceState::Aligned => {
                if row == 0 || column == 0 {
                    return Err(MafftError::InternalInvariant {
                        message: "E-INS-i aligned traceback crossed a sequence boundary",
                    });
                }
                operations.push(GeneralizedOperation::Aligned);
                let predecessor = packed & H_MASK;
                row -= 1;
                column -= 1;
                state = match predecessor {
                    H_RESTART => break,
                    H_FROM_H => TraceState::Aligned,
                    H_FROM_X => TraceState::OrdinaryFirstOnly,
                    H_FROM_Y => TraceState::OrdinarySecondOnly,
                    H_FROM_Z => TraceState::Skip,
                    _ => {
                        return Err(MafftError::InternalInvariant {
                            message: "E-INS-i aligned traceback has an invalid predecessor",
                        });
                    }
                };
            }
            TraceState::OrdinaryFirstOnly => {
                if row == 0 {
                    return Err(MafftError::InternalInvariant {
                        message: "E-INS-i first-only traceback crossed its boundary",
                    });
                }
                operations.push(GeneralizedOperation::OrdinaryFirstOnly);
                row -= 1;
                state = if packed & X_EXTENDS == 0 {
                    TraceState::Aligned
                } else {
                    TraceState::OrdinaryFirstOnly
                };
            }
            TraceState::OrdinarySecondOnly => {
                if column == 0 {
                    return Err(MafftError::InternalInvariant {
                        message: "E-INS-i second-only traceback crossed its boundary",
                    });
                }
                operations.push(GeneralizedOperation::OrdinarySecondOnly);
                column -= 1;
                state = if packed & Y_EXTENDS == 0 {
                    TraceState::Aligned
                } else {
                    TraceState::OrdinarySecondOnly
                };
            }
            TraceState::Skip => match (packed & Z_MASK) >> Z_SHIFT {
                Z_OPEN_PAIR | Z_EXTEND_PAIR => {
                    if row == 0 || column == 0 {
                        return Err(MafftError::InternalInvariant {
                            message: "E-INS-i paired skip crossed a sequence boundary",
                        });
                    }
                    operations.push(GeneralizedOperation::SkipPair);
                    row -= 1;
                    column -= 1;
                    state = if (packed & Z_MASK) >> Z_SHIFT == Z_OPEN_PAIR {
                        TraceState::Aligned
                    } else {
                        TraceState::Skip
                    };
                }
                Z_OPEN_FIRST | Z_EXTEND_FIRST => {
                    if row == 0 {
                        return Err(MafftError::InternalInvariant {
                            message: "E-INS-i first-only skip crossed its boundary",
                        });
                    }
                    operations.push(GeneralizedOperation::SkipFirstOnly);
                    row -= 1;
                    state = if (packed & Z_MASK) >> Z_SHIFT == Z_OPEN_FIRST {
                        TraceState::Aligned
                    } else {
                        TraceState::Skip
                    };
                }
                Z_OPEN_SECOND | Z_EXTEND_SECOND => {
                    if column == 0 {
                        return Err(MafftError::InternalInvariant {
                            message: "E-INS-i second-only skip crossed its boundary",
                        });
                    }
                    operations.push(GeneralizedOperation::SkipSecondOnly);
                    column -= 1;
                    state = if (packed & Z_MASK) >> Z_SHIFT == Z_OPEN_SECOND {
                        TraceState::Aligned
                    } else {
                        TraceState::Skip
                    };
                }
                _ => {
                    return Err(MafftError::InternalInvariant {
                        message: "E-INS-i generalized traceback has an invalid predecessor",
                    });
                }
            },
        }
    }
    operations.reverse();
    Ok(GeneralizedTraceback {
        operations,
        first_range: row..dp.best_end.0,
        second_range: column..dp.best_end.1,
    })
}

fn serialized_operation_width(operations: &[GeneralizedOperation]) -> Result<usize, MafftError> {
    operations
        .iter()
        .try_fold(0usize, |width, operation| {
            let columns = usize::from(*operation == GeneralizedOperation::SkipPair) + 1;
            width.checked_add(columns)
        })
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "E-INS-i generalized serialized operation width",
        })
}

fn padded_pair_profile<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
    first_range: &Range<usize>,
    second_range: &Range<usize>,
    operations: &[GeneralizedOperation],
) -> Result<WorkingProfile<A>, MafftError> {
    let first_sequence = records[first].sequence();
    let second_sequence = records[second].sequence();
    let suffix_first =
        first_sequence
            .len()
            .checked_sub(first_range.end)
            .ok_or(MafftError::InternalInvariant {
                message: "E-INS-i generalized first range exceeds its source",
            })?;
    let suffix_second = second_sequence.len().checked_sub(second_range.end).ok_or(
        MafftError::InternalInvariant {
            message: "E-INS-i generalized second range exceeds its source",
        },
    )?;
    let operation_width = serialized_operation_width(operations)?;
    let width = checked_sum(
        &[
            first_range.start,
            second_range.start,
            operation_width,
            suffix_first,
            suffix_second,
        ],
        "E-INS-i padded pair profile width",
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
            GeneralizedOperation::Aligned => {
                first_cells.push(first_sequence.get(first_position));
                second_cells.push(second_sequence.get(second_position));
                first_position += 1;
                second_position += 1;
            }
            GeneralizedOperation::OrdinaryFirstOnly | GeneralizedOperation::SkipFirstOnly => {
                first_cells.push(first_sequence.get(first_position));
                second_cells.push(None);
                first_position += 1;
            }
            GeneralizedOperation::OrdinarySecondOnly | GeneralizedOperation::SkipSecondOnly => {
                first_cells.push(None);
                second_cells.push(second_sequence.get(second_position));
                second_position += 1;
            }
            GeneralizedOperation::SkipPair => {
                first_cells.push(first_sequence.get(first_position));
                second_cells.push(None);
                first_cells.push(None);
                second_cells.push(second_sequence.get(second_position));
                first_position += 1;
                second_position += 1;
            }
        }
    }
    if first_position != first_range.end || second_position != second_range.end {
        return Err(MafftError::InternalInvariant {
            message: "E-INS-i generalized operations and source ranges disagree",
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
            message: "E-INS-i padded generalized pair profile is not rectangular",
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

#[allow(clippy::too_many_lines)]
fn replay_generalized<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
    scorer: &ScientificScorer,
    parameters: GeneralizedParameters,
    pair: &GeneralizedAffinePair<A>,
) -> Result<f64, MafftError> {
    let first_sequence = records[first].sequence();
    let second_sequence = records[second].sequence();
    if pair.first_range.end > first_sequence.len()
        || pair.second_range.end > second_sequence.len()
        || pair.first_range.start > pair.first_range.end
        || pair.second_range.start > pair.second_range.end
    {
        return Err(MafftError::InternalInvariant {
            message: "E-INS-i generalized traceback range is outside its source",
        });
    }
    if pair.score == 0.0 {
        if pair.first_range != (0..0)
            || pair.second_range != (0..0)
            || !pair.operations.is_empty()
            || pair.skip_opening_count != 0
            || pair.skipped_pair_count != 0
            || pair.skipped_unpaired_residue_count != 0
        {
            return Err(MafftError::InternalInvariant {
                message: "E-INS-i zero-score generalized result is not canonical",
            });
        }
        return Ok(0.0);
    }
    if pair.operations.first() != Some(&GeneralizedOperation::Aligned)
        || pair.operations.last() != Some(&GeneralizedOperation::Aligned)
    {
        return Err(MafftError::InternalInvariant {
            message: "E-INS-i positive generalized path does not start and end aligned",
        });
    }

    let mut first_position = pair.first_range.start;
    let mut second_position = pair.second_range.start;
    let mut score = 0.0;
    let mut previous: Option<GeneralizedOperation> = None;
    let mut skip_opening_count = 0usize;
    let mut skipped_pair_count = 0usize;
    let mut skipped_unpaired_residue_count = 0usize;
    for operation in pair.operations.iter().copied() {
        match operation {
            GeneralizedOperation::Aligned => {
                let first_symbol =
                    first_sequence
                        .get(first_position)
                        .ok_or(MafftError::InternalInvariant {
                            message: "E-INS-i replay over-consumed the first sequence",
                        })?;
                let second_symbol =
                    second_sequence
                        .get(second_position)
                        .ok_or(MafftError::InternalInvariant {
                            message: "E-INS-i replay over-consumed the second sequence",
                        })?;
                score += scorer.score_indices(
                    first_symbol.alphabet_index(),
                    second_symbol.alphabet_index(),
                );
                first_position += 1;
                second_position += 1;
            }
            GeneralizedOperation::OrdinaryFirstOnly => {
                if !matches!(
                    previous,
                    Some(GeneralizedOperation::Aligned | GeneralizedOperation::OrdinaryFirstOnly,)
                ) || first_position >= first_sequence.len()
                {
                    return Err(MafftError::InternalInvariant {
                        message: "E-INS-i replay found an invalid ordinary first-only move",
                    });
                }
                score += if previous == Some(GeneralizedOperation::OrdinaryFirstOnly) {
                    parameters.ordinary_extension
                } else {
                    parameters.ordinary_open
                };
                first_position += 1;
            }
            GeneralizedOperation::OrdinarySecondOnly => {
                if !matches!(
                    previous,
                    Some(GeneralizedOperation::Aligned | GeneralizedOperation::OrdinarySecondOnly,)
                ) || second_position >= second_sequence.len()
                {
                    return Err(MafftError::InternalInvariant {
                        message: "E-INS-i replay found an invalid ordinary second-only move",
                    });
                }
                score += if previous == Some(GeneralizedOperation::OrdinarySecondOnly) {
                    parameters.ordinary_extension
                } else {
                    parameters.ordinary_open
                };
                second_position += 1;
            }
            GeneralizedOperation::SkipPair
            | GeneralizedOperation::SkipFirstOnly
            | GeneralizedOperation::SkipSecondOnly => {
                let previous_is_skip = matches!(
                    previous,
                    Some(
                        GeneralizedOperation::SkipPair
                            | GeneralizedOperation::SkipFirstOnly
                            | GeneralizedOperation::SkipSecondOnly,
                    )
                );
                if previous != Some(GeneralizedOperation::Aligned) && !previous_is_skip {
                    return Err(MafftError::InternalInvariant {
                        message: "E-INS-i replay found an invalid generalized skip transition",
                    });
                }
                if !previous_is_skip {
                    score += parameters.skip_open;
                    skip_opening_count += 1;
                }
                match operation {
                    GeneralizedOperation::SkipPair => {
                        if first_position >= first_sequence.len()
                            || second_position >= second_sequence.len()
                        {
                            return Err(MafftError::InternalInvariant {
                                message: "E-INS-i replay over-consumed a paired skip",
                            });
                        }
                        score += parameters.paired_skip_extension;
                        first_position += 1;
                        second_position += 1;
                        skipped_pair_count += 1;
                    }
                    GeneralizedOperation::SkipFirstOnly => {
                        if first_position >= first_sequence.len() {
                            return Err(MafftError::InternalInvariant {
                                message: "E-INS-i replay over-consumed a first-only skip",
                            });
                        }
                        score += parameters.ordinary_extension;
                        first_position += 1;
                        skipped_unpaired_residue_count += 1;
                    }
                    GeneralizedOperation::SkipSecondOnly => {
                        if second_position >= second_sequence.len() {
                            return Err(MafftError::InternalInvariant {
                                message: "E-INS-i replay over-consumed a second-only skip",
                            });
                        }
                        score += parameters.ordinary_extension;
                        second_position += 1;
                        skipped_unpaired_residue_count += 1;
                    }
                    _ => unreachable!(),
                }
            }
        }
        if !score.is_finite() {
            return Err(MafftError::NonFiniteComputation {
                stage: "E-INS-i generalized traceback replay",
            });
        }
        if score <= 0.0 {
            return Err(MafftError::InternalInvariant {
                message: "E-INS-i positive generalized path crossed a non-positive prefix",
            });
        }
        previous = Some(operation);
    }
    if first_position != pair.first_range.end || second_position != pair.second_range.end {
        return Err(MafftError::InternalInvariant {
            message: "E-INS-i generalized replay did not consume its declared ranges",
        });
    }
    if skip_opening_count != pair.skip_opening_count
        || skipped_pair_count != pair.skipped_pair_count
        || skipped_unpaired_residue_count != pair.skipped_unpaired_residue_count
    {
        return Err(MafftError::InternalInvariant {
            message: "E-INS-i generalized replay statistics disagree",
        });
    }
    if !scores_agree(score, pair.score) {
        return Err(MafftError::InternalInvariant {
            message: "E-INS-i generalized DP and independent replay disagree",
        });
    }
    Ok(score)
}

fn generalized_memory<A: MafftAlphabet>(
    dp: &GeneralizedDpResult,
    operations: &Vec<GeneralizedOperation>,
    profile: &WorkingProfile<A>,
) -> Result<GeneralizedAffineMemory, MafftError> {
    let score_rows_bytes = checked_product(
        checked_product(
            dp.columns,
            size_of::<f64>(),
            "E-INS-i generalized score row bytes",
        )?,
        8,
        "E-INS-i generalized score rows",
    )?;
    let traceback_bytes = checked_product(
        dp.cell_count,
        size_of::<u16>(),
        "E-INS-i generalized traceback bytes",
    )?;
    let operation_bytes = checked_product(
        operations.capacity(),
        size_of::<GeneralizedOperation>(),
        "E-INS-i generalized operation bytes",
    )?;
    let profile_cells = checked_product(profile.width, 2, "E-INS-i generalized profile cells")?;
    let profile_bytes = checked_product(
        profile_cells,
        size_of::<Option<align_core::Symbol<A>>>(),
        "E-INS-i generalized profile bytes",
    )?;
    let estimated_peak_bytes = checked_sum(
        &[
            score_rows_bytes,
            traceback_bytes,
            operation_bytes,
            profile_bytes,
        ],
        "E-INS-i generalized pair peak memory",
    )?;
    Ok(GeneralizedAffineMemory {
        score_rows_bytes,
        traceback_bytes,
        operation_bytes,
        profile_bytes,
        estimated_peak_bytes,
    })
}

pub(crate) fn align_generalized_affine<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    first: usize,
    second: usize,
    scorer: &ScientificScorer,
    parameters: GeneralizedParameters,
) -> Result<GeneralizedAffinePair<A>, MafftError> {
    if first >= records.len() || second >= records.len() || first == second {
        return Err(MafftError::InternalInvariant {
            message: "E-INS-i generalized pair indices are invalid",
        });
    }
    if !parameters.ordinary_open.is_finite()
        || !parameters.ordinary_extension.is_finite()
        || !parameters.skip_open.is_finite()
        || !parameters.paired_skip_extension.is_finite()
        || parameters.ordinary_open >= 0.0
        || parameters.ordinary_extension > 0.0
        || parameters.skip_open >= 0.0
        || parameters.paired_skip_extension > 0.0
    {
        return Err(MafftError::NonFiniteComputation {
            stage: "E-INS-i generalized-affine parameters",
        });
    }
    let first_sequence = records[first].sequence();
    let second_sequence = records[second].sequence();
    let dp = fill_generalized_dp(first_sequence, second_sequence, scorer, parameters)?;
    let operation_capacity = first_sequence
        .len()
        .checked_add(second_sequence.len())
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "E-INS-i generalized traceback operation capacity",
        })?;
    let traceback = trace_generalized_dp(&dp, operation_capacity)?;
    let profile = padded_pair_profile(
        records,
        first,
        second,
        &traceback.first_range,
        &traceback.second_range,
        &traceback.operations,
    )?;
    let memory = generalized_memory(&dp, &traceback.operations, &profile)?;
    let mut skip_opening_count = 0usize;
    let mut skipped_pair_count = 0usize;
    let mut skipped_unpaired_residue_count = 0usize;
    let mut in_skip = false;
    for operation in traceback.operations.iter().copied() {
        match operation {
            GeneralizedOperation::SkipPair => {
                if !in_skip {
                    skip_opening_count += 1;
                }
                skipped_pair_count += 1;
                in_skip = true;
            }
            GeneralizedOperation::SkipFirstOnly | GeneralizedOperation::SkipSecondOnly => {
                if !in_skip {
                    skip_opening_count += 1;
                }
                skipped_unpaired_residue_count += 1;
                in_skip = true;
            }
            _ => in_skip = false,
        }
    }
    let pair = GeneralizedAffinePair {
        profile,
        score: dp.best_score,
        first_range: traceback.first_range,
        second_range: traceback.second_range,
        operations: traceback.operations.into_boxed_slice(),
        skip_opening_count,
        skipped_pair_count,
        skipped_unpaired_residue_count,
        memory,
    };
    replay_generalized(records, first, second, scorer, parameters, &pair)?;
    Ok(pair)
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Sequence};
    use align_msa::SequenceRecord;

    use super::{GeneralizedOperation, GeneralizedParameters, align_generalized_affine};
    use crate::matrix::ScientificScorer;

    fn record(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    const fn parameters(
        ordinary_open: f64,
        ordinary_extension: f64,
        skip_open: f64,
        paired_skip_extension: f64,
    ) -> GeneralizedParameters {
        GeneralizedParameters {
            ordinary_open,
            ordinary_extension,
            skip_open,
            paired_skip_extension,
        }
    }

    fn operation_score(
        previous: GeneralizedOperation,
        next: GeneralizedOperation,
        parameters: GeneralizedParameters,
    ) -> Option<f64> {
        use GeneralizedOperation::{
            Aligned, OrdinaryFirstOnly, OrdinarySecondOnly, SkipFirstOnly, SkipPair, SkipSecondOnly,
        };
        match (previous, next) {
            (Aligned, OrdinaryFirstOnly | OrdinarySecondOnly) => Some(parameters.ordinary_open),
            (OrdinaryFirstOnly, OrdinaryFirstOnly)
            | (OrdinarySecondOnly, OrdinarySecondOnly)
            | (SkipPair | SkipFirstOnly | SkipSecondOnly, SkipFirstOnly | SkipSecondOnly) => {
                Some(parameters.ordinary_extension)
            }
            (Aligned, SkipPair) => Some(parameters.skip_open + parameters.paired_skip_extension),
            (Aligned, SkipFirstOnly | SkipSecondOnly) => {
                Some(parameters.skip_open + parameters.ordinary_extension)
            }
            (SkipPair | SkipFirstOnly | SkipSecondOnly, SkipPair) => {
                Some(parameters.paired_skip_extension)
            }
            (
                Aligned | OrdinaryFirstOnly | OrdinarySecondOnly | SkipPair | SkipFirstOnly
                | SkipSecondOnly,
                Aligned,
            ) => Some(0.0),
            _ => None,
        }
    }

    struct ExhaustiveSearch<'a> {
        first: &'a Sequence<Dna>,
        second: &'a Sequence<Dna>,
        scorer: &'a ScientificScorer,
        parameters: GeneralizedParameters,
        best: f64,
    }

    impl ExhaustiveSearch<'_> {
        fn aligned_score(&self, first: usize, second: usize) -> f64 {
            self.scorer.score_indices(
                self.first.get(first).unwrap().alphabet_index(),
                self.second.get(second).unwrap().alphabet_index(),
            )
        }

        fn visit(
            &mut self,
            first: usize,
            second: usize,
            previous: GeneralizedOperation,
            score: f64,
        ) {
            if previous == GeneralizedOperation::Aligned {
                self.best = self.best.max(score);
            }
            let candidates = [
                (GeneralizedOperation::Aligned, 1usize, 1usize),
                (GeneralizedOperation::OrdinaryFirstOnly, 1, 0),
                (GeneralizedOperation::OrdinarySecondOnly, 0, 1),
                (GeneralizedOperation::SkipPair, 1, 1),
                (GeneralizedOperation::SkipFirstOnly, 1, 0),
                (GeneralizedOperation::SkipSecondOnly, 0, 1),
            ];
            for (next, first_step, second_step) in candidates {
                if first + first_step > self.first.len() || second + second_step > self.second.len()
                {
                    continue;
                }
                let Some(mut contribution) = operation_score(previous, next, self.parameters)
                else {
                    continue;
                };
                if next == GeneralizedOperation::Aligned {
                    contribution += self.aligned_score(first, second);
                }
                self.visit(
                    first + first_step,
                    second + second_step,
                    next,
                    score + contribution,
                );
            }
        }

        fn run(mut self) -> f64 {
            for first in 0..self.first.len() {
                for second in 0..self.second.len() {
                    let score = self.aligned_score(first, second);
                    self.best = self.best.max(score);
                    self.visit(first + 1, second + 1, GeneralizedOperation::Aligned, score);
                }
            }
            self.best.max(0.0)
        }
    }

    fn sequences_through_length_three() -> Vec<String> {
        let mut result = vec![String::new()];
        for length in 1..=3 {
            for bits in 0..(1usize << length) {
                let sequence = (0..length)
                    .map(|position| {
                        if bits & (1 << position) == 0 {
                            'A'
                        } else {
                            'C'
                        }
                    })
                    .collect();
                result.push(sequence);
            }
        }
        result
    }

    #[test]
    fn generalized_scores_match_complete_path_enumeration() {
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let parameters = GeneralizedParameters {
            ordinary_open: -0.45,
            ordinary_extension: -0.1,
            skip_open: -0.7,
            paired_skip_extension: 0.0,
        };
        for first in sequences_through_length_three() {
            for second in sequences_through_length_three() {
                let records = [record("a", &first), record("b", &second)];
                let observed =
                    align_generalized_affine(&records, 0, 1, &scorer, parameters).unwrap();
                let expected = ExhaustiveSearch {
                    first: records[0].sequence(),
                    second: records[1].sequence(),
                    scorer: &scorer,
                    parameters,
                    best: 0.0,
                }
                .run();
                assert!(
                    (observed.score - expected).abs() <= 1.0e-11,
                    "{first:?} versus {second:?}: {} != {}",
                    observed.score,
                    expected
                );
            }
        }
    }

    #[test]
    fn paired_skip_connects_two_aligned_blocks_without_homologizing_the_middle() {
        let records = [record("a", "AACCCCCCAA"), record("b", "AAGGGGGGAA")];
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let pair =
            align_generalized_affine(&records, 0, 1, &scorer, parameters(-2.0, -0.1, -0.2, 0.0))
                .unwrap();
        assert_eq!(pair.first_range, 0..10);
        assert_eq!(pair.second_range, 0..10);
        assert_eq!(pair.skip_opening_count, 1);
        assert_eq!(pair.skipped_pair_count, 6);
        assert_eq!(pair.skipped_unpaired_residue_count, 0);
        assert_eq!(
            pair.operations
                .iter()
                .filter(|operation| **operation == GeneralizedOperation::Aligned)
                .count(),
            4
        );
        assert_eq!(pair.profile.width, 16);
        let match_score = scorer.score_indices(0, 0);
        let expected = 4.0 * match_score - 0.2;
        assert!((pair.score - expected).abs() <= 1.0e-11);
    }

    #[test]
    fn unequal_skip_charges_unpaired_extensions() {
        let records = [record("a", "AACCCAA"), record("b", "AAGGGGGAA")];
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let pair =
            align_generalized_affine(&records, 0, 1, &scorer, parameters(-2.0, -0.1, -0.2, 0.0))
                .unwrap();
        assert_eq!(pair.skip_opening_count, 1);
        assert_eq!(pair.skipped_pair_count, 3);
        assert_eq!(pair.skipped_unpaired_residue_count, 2);
        let match_score = scorer.score_indices(0, 0);
        let expected = 4.0 * match_score - 0.2 - 2.0 * 0.1;
        assert!((pair.score - expected).abs() <= 1.0e-11);
    }

    #[test]
    fn ordinary_affine_gap_remains_distinct_from_a_skip() {
        let records = [record("a", "ACGA"), record("b", "ACA")];
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let pair = align_generalized_affine(
            &records,
            0,
            1,
            &scorer,
            parameters(-0.05, -0.01, -10.0, 0.0),
        )
        .unwrap();
        assert!(
            pair.operations
                .contains(&GeneralizedOperation::OrdinaryFirstOnly)
        );
        assert_eq!(pair.skip_opening_count, 0);
        assert_eq!(
            pair.operations.as_ref(),
            [
                GeneralizedOperation::Aligned,
                GeneralizedOperation::Aligned,
                GeneralizedOperation::OrdinaryFirstOnly,
                GeneralizedOperation::Aligned,
            ]
        );
        let first = records[0].sequence();
        let second = records[1].sequence();
        let expected = scorer.score_indices(
            first.get(0).unwrap().alphabet_index(),
            second.get(0).unwrap().alphabet_index(),
        ) + scorer.score_indices(
            first.get(1).unwrap().alphabet_index(),
            second.get(1).unwrap().alphabet_index(),
        ) + scorer.score_indices(
            first.get(3).unwrap().alphabet_index(),
            second.get(2).unwrap().alphabet_index(),
        ) - 0.05;
        assert!((pair.score - expected).abs() <= 1.0e-11);
    }

    #[test]
    fn ordinary_two_residue_gap_pays_one_open_and_one_extension() {
        let records = [record("a", "ACGGA"), record("b", "ACA")];
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let pair = align_generalized_affine(
            &records,
            0,
            1,
            &scorer,
            parameters(-0.05, -0.01, -10.0, 0.0),
        )
        .unwrap();
        assert_eq!(
            pair.operations.as_ref(),
            [
                GeneralizedOperation::Aligned,
                GeneralizedOperation::Aligned,
                GeneralizedOperation::OrdinaryFirstOnly,
                GeneralizedOperation::OrdinaryFirstOnly,
                GeneralizedOperation::Aligned,
            ]
        );
        let first = records[0].sequence();
        let second = records[1].sequence();
        let expected = scorer.score_indices(
            first.get(0).unwrap().alphabet_index(),
            second.get(0).unwrap().alphabet_index(),
        ) + scorer.score_indices(
            first.get(1).unwrap().alphabet_index(),
            second.get(1).unwrap().alphabet_index(),
        ) + scorer.score_indices(
            first.get(4).unwrap().alphabet_index(),
            second.get(2).unwrap().alphabet_index(),
        ) - 0.05
            - 0.01;
        assert!((pair.score - expected).abs() <= 1.0e-11);
    }

    #[test]
    fn first_row_major_endpoint_wins_an_exact_local_tie() {
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let first = [record("a", "A"), record("b", "AA")];
        let second = [record("a", "AA"), record("b", "A")];
        let first_result =
            align_generalized_affine(&first, 0, 1, &scorer, parameters(-2.0, -0.1, -6.0, 0.0))
                .unwrap();
        let second_result =
            align_generalized_affine(&second, 0, 1, &scorer, parameters(-2.0, -0.1, -6.0, 0.0))
                .unwrap();
        assert_eq!(first_result.first_range, 0..1);
        assert_eq!(first_result.second_range, 0..1);
        assert_eq!(second_result.first_range, 0..1);
        assert_eq!(second_result.second_range, 0..1);
    }

    #[test]
    fn no_positive_hit_is_canonical_and_full_sources_are_retained() {
        let records = [record("a", "AAAA"), record("b", "CCCC")];
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let pair =
            align_generalized_affine(&records, 0, 1, &scorer, parameters(-2.0, -0.1, -6.0, 0.0))
                .unwrap();
        assert_eq!(pair.score.to_bits(), 0.0_f64.to_bits());
        assert_eq!(pair.first_range, 0..0);
        assert_eq!(pair.second_range, 0..0);
        assert!(pair.operations.is_empty());
        assert_eq!(pair.profile.rows[0].cells.iter().flatten().count(), 4);
        assert_eq!(pair.profile.rows[1].cells.iter().flatten().count(), 4);
    }

    #[test]
    fn repeated_runs_are_bitwise_deterministic() {
        let records = [record("a", "AACCCCCCAA"), record("b", "AAGGGGGGAA")];
        let scorer = ScientificScorer::ginsi_2005::<Dna>(-0.1).unwrap();
        let first =
            align_generalized_affine(&records, 0, 1, &scorer, parameters(-2.0, -0.1, -0.2, 0.0))
                .unwrap();
        let second =
            align_generalized_affine(&records, 0, 1, &scorer, parameters(-2.0, -0.1, -0.2, 0.0))
                .unwrap();
        assert_eq!(first.score.to_bits(), second.score.to_bits());
        assert_eq!(first.first_range, second.first_range);
        assert_eq!(first.second_range, second.second_range);
        assert_eq!(first.operations, second.operations);
    }
}
