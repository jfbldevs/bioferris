use core::mem::size_of;

use align_core::Symbol;
use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord};

use crate::Muscle5Error;
use crate::muscle5_emissions::Muscle5Alphabet;
use crate::muscle5_posterior::PosteriorProvider;
use crate::muscle5_tree::Muscle5Node;

pub(crate) const OBJECTIVE_SCALE: f64 = 1_000_000_000.0;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct PositionedSymbol<A: Muscle5Alphabet> {
    pub(crate) symbol: Symbol<A>,
    pub(crate) position: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct WorkingRow<A: Muscle5Alphabet> {
    pub(crate) input_index: usize,
    pub(crate) cells: Vec<Option<PositionedSymbol<A>>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct WorkingProfile<A: Muscle5Alphabet> {
    pub(crate) rows: Vec<WorkingRow<A>>,
}

impl<A: Muscle5Alphabet> WorkingProfile<A> {
    pub(crate) fn single(record: &SequenceRecord<A>, input_index: usize) -> Self {
        Self {
            rows: vec![WorkingRow {
                input_index,
                cells: record
                    .sequence()
                    .iter()
                    .enumerate()
                    .map(|(position, symbol)| Some(PositionedSymbol { symbol, position }))
                    .collect(),
            }],
        }
    }

    pub(crate) fn width(&self) -> usize {
        self.rows.first().map_or(0, |row| row.cells.len())
    }

    pub(crate) fn members(&self) -> Vec<usize> {
        let mut members = self
            .rows
            .iter()
            .map(|row| row.input_index)
            .collect::<Vec<_>>();
        members.sort_unstable();
        members
    }

    pub(crate) fn extract(&self, selected: &[bool]) -> Result<Self, Muscle5Error> {
        let mut rows = self
            .rows
            .iter()
            .filter(|row| selected.get(row.input_index).copied().unwrap_or(false))
            .cloned()
            .collect::<Vec<_>>();
        if rows.is_empty() {
            return Err(Muscle5Error::InternalInvariant {
                message: "refinement extracted an empty profile",
            });
        }
        let width = self.width();
        let retained = (0..width)
            .map(|column| rows.iter().any(|row| row.cells[column].is_some()))
            .collect::<Vec<_>>();
        for row in &mut rows {
            row.cells = row
                .cells
                .iter()
                .zip(&retained)
                .filter_map(|(cell, keep)| keep.then_some(*cell))
                .collect();
        }
        Ok(Self { rows })
    }

    pub(crate) fn from_alignment(
        alignment: &MultipleAlignment<A>,
        global_indices: &[usize],
    ) -> Result<Self, Muscle5Error> {
        if alignment.sequence_count() != global_indices.len() {
            return Err(Muscle5Error::InternalInvariant {
                message: "aligned cluster and global-index map have different sizes",
            });
        }
        let rows = alignment
            .rows()
            .iter()
            .zip(global_indices)
            .map(|(row, input_index)| {
                let mut position = 0_usize;
                let cells = row
                    .symbols()
                    .iter()
                    .map(|cell| {
                        cell.map(|symbol| {
                            let positioned = PositionedSymbol { symbol, position };
                            position += 1;
                            positioned
                        })
                    })
                    .collect::<Vec<_>>();
                WorkingRow {
                    input_index: *input_index,
                    cells,
                }
            })
            .collect();
        Ok(Self { rows })
    }

    pub(crate) fn to_alignment(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<MultipleAlignment<A>, Muscle5Error> {
        if self.rows.len() != records.len() {
            return Err(Muscle5Error::InternalInvariant {
                message: "working profile does not contain every input row",
            });
        }
        let mut ordered = self.rows.iter().collect::<Vec<_>>();
        ordered.sort_by_key(|row| row.input_index);
        let rows = ordered
            .into_iter()
            .enumerate()
            .map(|(expected, row)| {
                if row.input_index != expected {
                    return Err(Muscle5Error::InternalInvariant {
                        message: "working profile input indices are not a permutation",
                    });
                }
                AlignedRow::new(
                    records[expected].identifier(),
                    row.cells
                        .iter()
                        .map(|cell| cell.map(|positioned| positioned.symbol))
                        .collect(),
                )
                .map_err(Muscle5Error::from)
            })
            .collect::<Result<Vec<_>, Muscle5Error>>()?;
        MultipleAlignment::try_new(rows).map_err(Muscle5Error::from)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ProfileOperation {
    Both,
    LeftOnly,
    RightOnly,
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum PairSelection<'a> {
    All,
    Selected(&'a [(usize, usize)]),
}

impl PairSelection<'_> {
    fn includes(self, first: usize, second: usize) -> bool {
        match self {
            Self::All => true,
            Self::Selected(pairs) => {
                let pair = canonical_pair(first, second);
                pairs.binary_search(&pair).is_ok()
            }
        }
    }
}

pub(crate) struct ProfileMerge<A: Muscle5Alphabet> {
    pub(crate) profile: WorkingProfile<A>,
    pub(crate) expected_accuracy: f64,
    pub(crate) operations: Box<[ProfileOperation]>,
    pub(crate) workspace_bytes: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct ProgressiveMergeEvidence {
    pub(crate) left_members: Box<[usize]>,
    pub(crate) right_members: Box<[usize]>,
    pub(crate) result_width: usize,
    pub(crate) expected_accuracy: f64,
    pub(crate) operation_count: usize,
}

pub(crate) struct ProgressiveResult<A: Muscle5Alphabet> {
    pub(crate) profile: WorkingProfile<A>,
    pub(crate) merges: Vec<ProgressiveMergeEvidence>,
    pub(crate) peak_workspace_bytes: usize,
}

/// A refinement candidate and the value of the arrangement it would replace,
/// both taken from one column-score matrix.
pub(crate) struct PartitionRealignment<A: Muscle5Alphabet> {
    pub(crate) candidate: ProfileMerge<A>,
    /// Cross-partition posterior sum of the current arrangement, added along its
    /// path over the candidate's column scores, exactly as the DP adds them.
    pub(crate) current_score: f64,
}

/// Posterior sums for every left/right column pair, plus the size of the
/// residue-to-column maps that the memory estimate counts.
struct ColumnScores {
    values: Vec<f64>,
    map_cells: usize,
}

pub(crate) fn align_profiles<A: Muscle5Alphabet, P: PosteriorProvider + ?Sized>(
    left: WorkingProfile<A>,
    right: WorkingProfile<A>,
    library: &P,
    pair_selection: PairSelection<'_>,
) -> Result<ProfileMerge<A>, Muscle5Error> {
    let scores = column_scores(&left, &right, library, pair_selection)?;
    merge_on_scores(left, right, &scores)
}

/// Realigns the two partitions of `profile` and replays its current
/// arrangement through the same column scores.
///
/// The DP and the replay add the same column scores along a path, in path
/// order. Rounded addition is monotone, so the optimum can never fall below the
/// current arrangement's value. An objective summed in another order, such as
/// row pair by row pair, gives no such guarantee.
pub(crate) fn realign_partition<A: Muscle5Alphabet, P: PosteriorProvider + ?Sized>(
    profile: &WorkingProfile<A>,
    first_partition: &[bool],
    library: &P,
) -> Result<PartitionRealignment<A>, Muscle5Error> {
    let current = partition_operations(profile, first_partition)?;
    let complement = first_partition
        .iter()
        .map(|selected| !selected)
        .collect::<Vec<_>>();
    let left = profile.extract(first_partition)?;
    let right = profile.extract(&complement)?;
    let scores = column_scores(&left, &right, library, PairSelection::All)?;
    let current_score = replay_score(&current, &scores.values, left.width(), right.width())?;
    let candidate = merge_on_scores(left, right, &scores)?;
    Ok(PartitionRealignment {
        candidate,
        current_score,
    })
}

/// The current arrangement of `profile` as a path through the lattice of its
/// two extracted partitions: each column consumes the columns that
/// [`WorkingProfile::extract`] keeps on each side.
fn partition_operations<A: Muscle5Alphabet>(
    profile: &WorkingProfile<A>,
    first_partition: &[bool],
) -> Result<Vec<ProfileOperation>, Muscle5Error> {
    let mut operations = Vec::with_capacity(profile.width());
    for column in 0..profile.width() {
        let (mut left, mut right) = (false, false);
        for row in &profile.rows {
            if row.cells[column].is_none() {
                continue;
            }
            match first_partition.get(row.input_index) {
                Some(true) => left = true,
                Some(false) => right = true,
                None => {
                    return Err(Muscle5Error::InternalInvariant {
                        message: "refinement partition omits a profile row",
                    });
                }
            }
        }
        operations.push(match (left, right) {
            (true, true) => ProfileOperation::Both,
            (true, false) => ProfileOperation::LeftOnly,
            (false, true) => ProfileOperation::RightOnly,
            (false, false) => {
                return Err(Muscle5Error::InternalInvariant {
                    message: "PPP profile contains an all-gap column",
                });
            }
        });
    }
    Ok(operations)
}

fn column_scores<A: Muscle5Alphabet, P: PosteriorProvider + ?Sized>(
    left: &WorkingProfile<A>,
    right: &WorkingProfile<A>,
    library: &P,
    pair_selection: PairSelection<'_>,
) -> Result<ColumnScores, Muscle5Error> {
    validate_profile(left, library)?;
    validate_profile(right, library)?;
    let right_width = right.width();
    let score_count =
        left.width()
            .checked_mul(right_width)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "PPP profile score matrix",
            })?;
    let mut values = vec![0.0_f64; score_count];
    let left_position_columns = position_column_maps(left, library)?;
    let right_position_columns = position_column_maps(right, library)?;
    for (left_row, left_columns) in left.rows.iter().zip(&left_position_columns) {
        for (right_row, right_columns) in right.rows.iter().zip(&right_position_columns) {
            if !pair_selection.includes(left_row.input_index, right_row.input_index) {
                continue;
            }
            let posterior = library
                .posterior(left_row.input_index, right_row.input_index)
                .ok_or(Muscle5Error::InternalInvariant {
                    message: "PPP profile pair is absent from its posterior provider",
                })?;
            for (left_position, left_column) in left_columns.iter().copied().enumerate() {
                for entry in posterior.row(left_position) {
                    let right_column = right_columns[entry.column];
                    values[left_column * right_width + right_column] += entry.probability;
                }
            }
        }
    }
    if values.iter().any(|score| !score.is_finite()) {
        return Err(Muscle5Error::NonFiniteComputation {
            stage: "PPP profile column score",
        });
    }
    let map_cells = left_position_columns
        .iter()
        .chain(&right_position_columns)
        .map(Vec::len)
        .sum::<usize>();
    Ok(ColumnScores { values, map_cells })
}

fn merge_on_scores<A: Muscle5Alphabet>(
    left: WorkingProfile<A>,
    right: WorkingProfile<A>,
    scores: &ColumnScores,
) -> Result<ProfileMerge<A>, Muscle5Error> {
    let left_width = left.width();
    let right_width = right.width();
    let column_scores = &scores.values;
    let lattice_rows = left_width
        .checked_add(1)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "PPP profile DP rows",
        })?;
    let lattice_columns = right_width
        .checked_add(1)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "PPP profile DP columns",
        })?;
    let lattice_cells =
        lattice_rows
            .checked_mul(lattice_columns)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "PPP profile DP lattice",
            })?;
    let mut trace = vec![0_u8; lattice_cells];
    for row in 1..lattice_rows {
        trace[row * lattice_columns] = 1;
    }
    trace[1..lattice_columns].fill(2);
    let mut previous_scores = vec![0.0_f64; lattice_columns];
    let mut current_scores = vec![0.0_f64; lattice_columns];
    for row in 1..lattice_rows {
        current_scores[0] = 0.0;
        for column in 1..lattice_columns {
            let diagonal =
                previous_scores[column - 1] + column_scores[(row - 1) * right_width + column - 1];
            let consume_left = previous_scores[column];
            let consume_right = current_scores[column - 1];
            let index = row * lattice_columns + column;
            if diagonal >= consume_left && diagonal >= consume_right {
                current_scores[column] = diagonal;
                trace[index] = 0;
            } else if consume_left >= consume_right {
                current_scores[column] = consume_left;
                trace[index] = 1;
            } else {
                current_scores[column] = consume_right;
                trace[index] = 2;
            }
        }
        core::mem::swap(&mut previous_scores, &mut current_scores);
    }

    let expected_accuracy = previous_scores[right_width];
    let operations = traceback(&trace, left_width, right_width, lattice_columns)?;
    let replayed = replay_score(&operations, column_scores, left_width, right_width)?;
    if (replayed - expected_accuracy).abs() > 1.0e-9 * expected_accuracy.abs().max(1.0) {
        return Err(Muscle5Error::InternalInvariant {
            message: "PPP profile traceback score disagrees with the DP optimum",
        });
    }
    let profile = merge_rows(left, right, &operations)?;
    let workspace_bytes = column_scores
        .len()
        .checked_mul(size_of::<f64>())
        .and_then(|bytes| bytes.checked_add(2 * lattice_columns * size_of::<f64>()))
        .and_then(|bytes| bytes.checked_add(trace.len() * size_of::<u8>()))
        .and_then(|bytes| bytes.checked_add(scores.map_cells * size_of::<usize>()))
        .and_then(|bytes| bytes.checked_add(operations.capacity() * size_of::<ProfileOperation>()))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "PPP profile DP memory estimate",
        })?;
    Ok(ProfileMerge {
        profile,
        expected_accuracy,
        operations: operations.into_boxed_slice(),
        workspace_bytes,
    })
}

fn position_column_maps<A: Muscle5Alphabet, P: PosteriorProvider + ?Sized>(
    profile: &WorkingProfile<A>,
    library: &P,
) -> Result<Vec<Vec<usize>>, Muscle5Error> {
    profile
        .rows
        .iter()
        .map(|row| {
            let length = library.sequence_length(row.input_index).ok_or(
                Muscle5Error::InternalInvariant {
                    message: "PPP row is absent from its posterior provider",
                },
            )?;
            let mut columns = vec![usize::MAX; length];
            for (column, cell) in row.cells.iter().enumerate() {
                if let Some(cell) = cell {
                    columns[cell.position] = column;
                }
            }
            if columns.contains(&usize::MAX) {
                return Err(Muscle5Error::InternalInvariant {
                    message: "PPP profile omitted an original residue position",
                });
            }
            Ok(columns)
        })
        .collect()
}

pub(crate) fn progressive_alignment<A: Muscle5Alphabet, P: PosteriorProvider + ?Sized>(
    root: &Muscle5Node,
    records: &[SequenceRecord<A>],
    library: &P,
) -> Result<ProgressiveResult<A>, Muscle5Error> {
    let mut merges = Vec::new();
    let mut peak_workspace_bytes = 0_usize;
    let profile = progressive_node(
        root,
        records,
        library,
        &mut merges,
        &mut peak_workspace_bytes,
    )?;
    Ok(ProgressiveResult {
        profile,
        merges,
        peak_workspace_bytes,
    })
}

fn progressive_node<A: Muscle5Alphabet, P: PosteriorProvider + ?Sized>(
    node: &Muscle5Node,
    records: &[SequenceRecord<A>],
    library: &P,
    merges: &mut Vec<ProgressiveMergeEvidence>,
    peak_workspace_bytes: &mut usize,
) -> Result<WorkingProfile<A>, Muscle5Error> {
    match node {
        Muscle5Node::Leaf { sequence_index } => records
            .get(*sequence_index)
            .map(|record| WorkingProfile::single(record, *sequence_index))
            .ok_or(Muscle5Error::InternalInvariant {
                message: "guide-tree leaf is outside the input collection",
            }),
        Muscle5Node::Branch { left, right, .. } => {
            let left_profile =
                progressive_node(left, records, library, merges, peak_workspace_bytes)?;
            let right_profile =
                progressive_node(right, records, library, merges, peak_workspace_bytes)?;
            let left_members = left_profile.members().into_boxed_slice();
            let right_members = right_profile.members().into_boxed_slice();
            let profile_merge =
                align_profiles(left_profile, right_profile, library, PairSelection::All)?;
            *peak_workspace_bytes = (*peak_workspace_bytes).max(profile_merge.workspace_bytes);
            merges.push(ProgressiveMergeEvidence {
                left_members,
                right_members,
                result_width: profile_merge.profile.width(),
                expected_accuracy: profile_merge.expected_accuracy,
                operation_count: profile_merge.operations.len(),
            });
            Ok(profile_merge.profile)
        }
    }
}

pub(crate) fn total_objective<A: Muscle5Alphabet, P: PosteriorProvider + ?Sized>(
    profile: &WorkingProfile<A>,
    library: &P,
) -> Result<f64, Muscle5Error> {
    validate_profile(profile, library)?;
    let mut score = 0.0_f64;
    for first in 0..profile.rows.len() {
        for second in first + 1..profile.rows.len() {
            score += row_pair_objective(&profile.rows[first], &profile.rows[second], library)?;
        }
    }
    finite_objective(score)
}

// Refinement decisions use `realign_partition`; this independent row-pair sum
// remains a cross-check in checked builds and tests.
#[cfg(any(test, debug_assertions))]
pub(crate) fn cross_partition_objective<A: Muscle5Alphabet, P: PosteriorProvider + ?Sized>(
    profile: &WorkingProfile<A>,
    first_partition: &[bool],
    library: &P,
) -> Result<f64, Muscle5Error> {
    validate_profile(profile, library)?;
    let mut score = 0.0_f64;
    for first in 0..profile.rows.len() {
        let first_side = first_partition
            .get(profile.rows[first].input_index)
            .copied()
            .unwrap_or(false);
        for second in first + 1..profile.rows.len() {
            let second_side = first_partition
                .get(profile.rows[second].input_index)
                .copied()
                .unwrap_or(false);
            if first_side != second_side {
                score += row_pair_objective(&profile.rows[first], &profile.rows[second], library)?;
            }
        }
    }
    finite_objective(score)
}

#[allow(clippy::cast_precision_loss)]
pub(crate) fn quantized_objective(score: f64) -> Result<i128, Muscle5Error> {
    let scaled = score * OBJECTIVE_SCALE;
    if !scaled.is_finite() || scaled < i128::MIN as f64 || scaled > i128::MAX as f64 {
        return Err(Muscle5Error::NonFiniteComputation {
            stage: "fixed-precision PPP objective",
        });
    }
    #[allow(clippy::cast_possible_truncation)]
    Ok(scaled.round() as i128)
}

fn row_pair_objective<A: Muscle5Alphabet, P: PosteriorProvider + ?Sized>(
    first: &WorkingRow<A>,
    second: &WorkingRow<A>,
    library: &P,
) -> Result<f64, Muscle5Error> {
    if first.cells.len() != second.cells.len() {
        return Err(Muscle5Error::InternalInvariant {
            message: "PPP objective received non-rectangular rows",
        });
    }
    let posterior = library
        .posterior(first.input_index, second.input_index)
        .ok_or(Muscle5Error::InternalInvariant {
            message: "PPP objective pair is absent from its posterior provider",
        })?;
    let mut score = 0.0_f64;
    for (first_cell, second_cell) in first.cells.iter().zip(&second.cells) {
        if let (Some(first_cell), Some(second_cell)) = (first_cell, second_cell) {
            score += posterior.get(first_cell.position, second_cell.position);
        }
    }
    finite_objective(score)
}

fn finite_objective(score: f64) -> Result<f64, Muscle5Error> {
    if score.is_finite() {
        Ok(score)
    } else {
        Err(Muscle5Error::NonFiniteComputation {
            stage: "PPP alignment objective",
        })
    }
}

fn validate_profile<A: Muscle5Alphabet, P: PosteriorProvider + ?Sized>(
    profile: &WorkingProfile<A>,
    library: &P,
) -> Result<(), Muscle5Error> {
    if profile.rows.is_empty() {
        return Err(Muscle5Error::InternalInvariant {
            message: "PPP profile has no rows",
        });
    }
    let width = profile.width();
    let mut seen = std::collections::HashSet::with_capacity(profile.rows.len());
    for row in &profile.rows {
        if row.cells.len() != width {
            return Err(Muscle5Error::InternalInvariant {
                message: "PPP profile is not rectangular",
            });
        }
        if !seen.insert(row.input_index) {
            return Err(Muscle5Error::InternalInvariant {
                message: "PPP profile repeats an input row",
            });
        }
        let sequence_length =
            library
                .sequence_length(row.input_index)
                .ok_or(Muscle5Error::InternalInvariant {
                    message: "PPP profile row index is outside the posterior provider",
                })?;
        let mut previous = None;
        for cell in row.cells.iter().flatten() {
            if cell.position >= sequence_length
                || previous.is_some_and(|position| position >= cell.position)
            {
                return Err(Muscle5Error::InternalInvariant {
                    message: "PPP profile residue positions are invalid",
                });
            }
            previous = Some(cell.position);
        }
    }
    for column in 0..width {
        if profile.rows.iter().all(|row| row.cells[column].is_none()) {
            return Err(Muscle5Error::InternalInvariant {
                message: "PPP profile contains an all-gap column",
            });
        }
    }
    Ok(())
}

fn traceback(
    trace: &[u8],
    mut left: usize,
    mut right: usize,
    columns: usize,
) -> Result<Vec<ProfileOperation>, Muscle5Error> {
    let mut reverse = Vec::with_capacity(left.saturating_add(right));
    while left != 0 || right != 0 {
        match trace[left * columns + right] {
            0 if left != 0 && right != 0 => {
                reverse.push(ProfileOperation::Both);
                left -= 1;
                right -= 1;
            }
            1 if left != 0 => {
                reverse.push(ProfileOperation::LeftOnly);
                left -= 1;
            }
            2 if right != 0 => {
                reverse.push(ProfileOperation::RightOnly);
                right -= 1;
            }
            _ => {
                return Err(Muscle5Error::InternalInvariant {
                    message: "invalid PPP profile traceback state",
                });
            }
        }
    }
    reverse.reverse();
    Ok(reverse)
}

fn replay_score(
    operations: &[ProfileOperation],
    column_scores: &[f64],
    left_width: usize,
    right_width: usize,
) -> Result<f64, Muscle5Error> {
    let mut left = 0_usize;
    let mut right = 0_usize;
    let mut score = 0.0_f64;
    for operation in operations {
        match operation {
            ProfileOperation::Both => {
                let value = column_scores.get(left * right_width + right).ok_or(
                    Muscle5Error::InternalInvariant {
                        message: "PPP score replay exceeded the column matrix",
                    },
                )?;
                score += value;
                left += 1;
                right += 1;
            }
            ProfileOperation::LeftOnly => left += 1,
            ProfileOperation::RightOnly => right += 1,
        }
    }
    if left != left_width || right != right_width {
        return Err(Muscle5Error::InternalInvariant {
            message: "PPP score replay did not consume both profiles",
        });
    }
    finite_objective(score)
}

fn merge_rows<A: Muscle5Alphabet>(
    left: WorkingProfile<A>,
    right: WorkingProfile<A>,
    operations: &[ProfileOperation],
) -> Result<WorkingProfile<A>, Muscle5Error> {
    let mut rows = Vec::with_capacity(left.rows.len() + right.rows.len());
    for row in left.rows {
        rows.push(WorkingRow {
            input_index: row.input_index,
            cells: replay_row(&row.cells, operations, true)?,
        });
    }
    for row in right.rows {
        rows.push(WorkingRow {
            input_index: row.input_index,
            cells: replay_row(&row.cells, operations, false)?,
        });
    }
    rows.sort_by_key(|row| row.input_index);
    Ok(WorkingProfile { rows })
}

fn replay_row<A: Muscle5Alphabet>(
    source: &[Option<PositionedSymbol<A>>],
    operations: &[ProfileOperation],
    left: bool,
) -> Result<Vec<Option<PositionedSymbol<A>>>, Muscle5Error> {
    let mut result = Vec::with_capacity(operations.len());
    let mut source_column = 0_usize;
    for operation in operations {
        let consumes = matches!(operation, ProfileOperation::Both)
            || (left && matches!(operation, ProfileOperation::LeftOnly))
            || (!left && matches!(operation, ProfileOperation::RightOnly));
        if consumes {
            result.push(source.get(source_column).copied().ok_or(
                Muscle5Error::InternalInvariant {
                    message: "PPP profile merge consumed beyond a source row",
                },
            )?);
            source_column += 1;
        } else {
            result.push(None);
        }
    }
    if source_column != source.len() {
        return Err(Muscle5Error::InternalInvariant {
            message: "PPP profile merge did not consume a source row",
        });
    }
    Ok(result)
}

const fn canonical_pair(first: usize, second: usize) -> (usize, usize) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

#[cfg(test)]
mod tests {
    use align_core::DnaSequence;
    use align_msa::SequenceRecord;

    use super::{
        PairSelection, ProfileOperation, WorkingProfile, WorkingRow, align_profiles,
        cross_partition_objective, progressive_alignment, quantized_objective, realign_partition,
        total_objective,
    };
    use crate::muscle5_posterior::{PosteriorLibrary, SparsePosterior};
    use crate::muscle5_tree::Muscle5Node;

    fn record(identifier: &str, sequence: &str) -> SequenceRecord<align_core::Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    fn hand_library() -> PosteriorLibrary {
        PosteriorLibrary::from_test_upper_triangle(
            &[2, 3, 2],
            vec![
                SparsePosterior::from_dense(2, 3, &[0.9, 0.1, 0.0, 0.0, 0.2, 0.8], 0.0).unwrap(),
                SparsePosterior::from_dense(2, 2, &[0.7, 0.0, 0.0, 0.6], 0.0).unwrap(),
                SparsePosterior::from_dense(3, 2, &[0.5, 0.0, 0.0, 0.1, 0.0, 0.4], 0.0).unwrap(),
            ],
        )
        .unwrap()
    }

    #[test]
    fn singleton_profile_alignment_matches_hand_expected_accuracy() {
        let records = [record("x", "AC"), record("y", "AGC")];
        let library = hand_library();
        let merged = align_profiles(
            WorkingProfile::single(&records[0], 0),
            WorkingProfile::single(&records[1], 1),
            &library,
            PairSelection::All,
        )
        .unwrap();
        assert_eq!(
            merged.operations.as_ref(),
            &[
                ProfileOperation::Both,
                ProfileOperation::RightOnly,
                ProfileOperation::Both
            ]
        );
        assert!((merged.expected_accuracy - 1.7).abs() < 1.0e-12);
        let alignment = merged.profile.to_alignment(&records).unwrap();
        assert_eq!(alignment.row(0).unwrap().to_ascii_with_gaps(), "A-C");
        assert_eq!(alignment.row(1).unwrap().to_ascii_with_gaps(), "AGC");
    }

    #[test]
    fn empty_profile_is_expanded_without_losing_residue_positions() {
        let records = [record("empty", ""), record("full", "AGC")];
        let library = PosteriorLibrary::from_test_upper_triangle(
            &[0, 3],
            vec![SparsePosterior::from_dense(0, 3, &[], 0.0).unwrap()],
        )
        .unwrap();
        let merged = align_profiles(
            WorkingProfile::single(&records[0], 0),
            WorkingProfile::single(&records[1], 1),
            &library,
            PairSelection::All,
        )
        .unwrap();
        assert_eq!(
            merged.operations.as_ref(),
            &[ProfileOperation::RightOnly; 3]
        );
        let alignment = merged.profile.to_alignment(&records).unwrap();
        assert_eq!(alignment.row(0).unwrap().to_ascii_with_gaps(), "---");
        assert_eq!(alignment.row(1).unwrap().to_ascii_with_gaps(), "AGC");
    }

    #[test]
    fn selected_pairs_exclude_unselected_cross_terms() {
        let records = [record("x", "AC"), record("y", "AGC"), record("z", "AT")];
        let library = hand_library();
        let right = align_profiles(
            WorkingProfile::single(&records[1], 1),
            WorkingProfile::single(&records[2], 2),
            &library,
            PairSelection::All,
        )
        .unwrap()
        .profile;
        let selected = [(0, 2)];
        let merged = align_profiles(
            WorkingProfile::single(&records[0], 0),
            right,
            &library,
            PairSelection::Selected(&selected),
        )
        .unwrap();
        assert!((merged.expected_accuracy - 1.3).abs() < 1.0e-12);
    }

    #[test]
    fn objective_is_independently_replayed_and_partitioned() {
        let records = [record("x", "AC"), record("y", "AGC"), record("z", "AT")];
        let library = hand_library();
        let first = align_profiles(
            WorkingProfile::single(&records[0], 0),
            WorkingProfile::single(&records[1], 1),
            &library,
            PairSelection::All,
        )
        .unwrap();
        let complete = align_profiles(
            first.profile,
            WorkingProfile::single(&records[2], 2),
            &library,
            PairSelection::All,
        )
        .unwrap()
        .profile;
        let total = total_objective(&complete, &library).unwrap();
        let cross = cross_partition_objective(&complete, &[true, false, false], &library).unwrap();
        let within = cross_partition_objective(&complete, &[false, true, false], &library).unwrap();
        assert!(total >= cross);
        assert!(within > 0.0);
        assert_eq!(quantized_objective(1.234_567_890_4).unwrap(), 1_234_567_890);
    }

    #[test]
    fn extraction_removes_only_partition_all_gap_columns() {
        let records = [record("x", "AC"), record("y", "AGC")];
        let library = hand_library();
        let profile = align_profiles(
            WorkingProfile::single(&records[0], 0),
            WorkingProfile::single(&records[1], 1),
            &library,
            PairSelection::All,
        )
        .unwrap()
        .profile;
        let extracted = profile.extract(&[true, false, false]).unwrap();
        assert_eq!(extracted.width(), 2);
        assert_eq!(extracted.rows[0].cells[0].unwrap().position, 0);
        assert_eq!(extracted.rows[0].cells[1].unwrap().position, 1);
    }

    #[test]
    fn progressive_traversal_records_postorder_merges() {
        let records = [record("x", "AC"), record("y", "AGC"), record("z", "AT")];
        let library = hand_library();
        let tree = Muscle5Node::Branch {
            left: Box::new(Muscle5Node::Branch {
                left: Box::new(Muscle5Node::Leaf { sequence_index: 0 }),
                right: Box::new(Muscle5Node::Leaf { sequence_index: 1 }),
                left_length: 0.1,
                right_length: 0.1,
            }),
            right: Box::new(Muscle5Node::Leaf { sequence_index: 2 }),
            left_length: 0.2,
            right_length: 0.2,
        };
        let result = progressive_alignment(&tree, &records, &library).unwrap();
        assert_eq!(result.merges.len(), 2);
        assert_eq!(result.merges[0].left_members.as_ref(), &[0]);
        assert_eq!(result.merges[0].right_members.as_ref(), &[1]);
        assert_eq!(result.merges[1].left_members.as_ref(), &[0, 1]);
        assert_eq!(result.merges[1].right_members.as_ref(), &[2]);
        assert!(result.peak_workspace_bytes > 0);
    }

    #[test]
    fn nonrectangular_profile_is_rejected() {
        let records = [record("x", "AC"), record("y", "AGC")];
        let library = hand_library();
        let malformed = WorkingProfile {
            rows: vec![
                WorkingProfile::single(&records[0], 0).rows.remove(0),
                WorkingRow {
                    input_index: 1,
                    cells: vec![],
                },
            ],
        };
        assert!(total_objective(&malformed, &library).is_err());
    }

    #[test]
    fn partition_realignment_cannot_round_below_the_current_arrangement() {
        // These cross pairs make the current alignment's posterior sum round onto
        // adjacent 1e-9 steps: added row pair by row pair it is 2017472177.5000002
        // units, added column by column, as the DP adds it, 2017472177.4999998.
        let records = [record("x", "AC"), record("y", "AC"), record("z", "AC")];
        let library = PosteriorLibrary::from_test_upper_triangle(
            &[2, 2, 2],
            vec![
                SparsePosterior::from_dense(2, 2, &[0.9, 0.0, 0.0, 0.9], 0.0).unwrap(),
                SparsePosterior::from_dense(2, 2, &[0.429_529, 0.0, 0.0, 0.915_817], 0.0).unwrap(),
                SparsePosterior::from_dense(2, 2, &[0.119_858, 0.0, 0.0, 0.552_268_177_5], 0.0)
                    .unwrap(),
            ],
        )
        .unwrap();
        let pair = align_profiles(
            WorkingProfile::single(&records[0], 0),
            WorkingProfile::single(&records[1], 1),
            &library,
            PairSelection::All,
        )
        .unwrap()
        .profile;
        let current = align_profiles(
            pair,
            WorkingProfile::single(&records[2], 2),
            &library,
            PairSelection::All,
        )
        .unwrap()
        .profile;
        let partition = [true, true, false];

        let row_pair = cross_partition_objective(&current, &partition, &library).unwrap();
        let realigned = realign_partition(&current, &partition, &library).unwrap();
        let optimum = realigned.candidate.expected_accuracy;
        // Version 1 compared these two values and reported a decrease.
        assert_eq!(quantized_objective(row_pair).unwrap(), 2_017_472_178);
        assert_eq!(quantized_objective(optimum).unwrap(), 2_017_472_177);
        // Version 2 replays the current path through the DP's column scores: the
        // realignment returns the same alignment with exactly the same score.
        assert_eq!(realigned.candidate.profile, current);
        assert_eq!(realigned.current_score.to_bits(), optimum.to_bits());
    }

    #[test]
    fn partition_realignment_never_places_the_optimum_below_the_current_path() {
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            u32::try_from(state >> 32).unwrap()
        };
        for _ in 0..200 {
            let lengths = [0_usize; 4].map(|_| 1 + usize::try_from(next() % 6).unwrap());
            let sources = "ACGTACGT";
            let records =
                [0, 1, 2, 3].map(|index| record(&format!("r{index}"), &sources[..lengths[index]]));
            let mut matrices = Vec::new();
            for first in 0..4 {
                for second in first + 1..4 {
                    let values = (0..lengths[first] * lengths[second])
                        .map(|_| {
                            let draw = next();
                            if draw % 2 == 0 {
                                0.0
                            } else {
                                f64::from(draw) / f64::from(u32::MAX)
                            }
                        })
                        .collect::<Vec<_>>();
                    matrices.push(
                        SparsePosterior::from_dense(lengths[first], lengths[second], &values, 0.0)
                            .unwrap(),
                    );
                }
            }
            let library = PosteriorLibrary::from_test_upper_triangle(&lengths, matrices).unwrap();
            let merge = |left, right| {
                align_profiles(left, right, &library, PairSelection::All)
                    .unwrap()
                    .profile
            };
            let current = merge(
                merge(
                    WorkingProfile::single(&records[0], 0),
                    WorkingProfile::single(&records[1], 1),
                ),
                merge(
                    WorkingProfile::single(&records[2], 2),
                    WorkingProfile::single(&records[3], 3),
                ),
            );
            for mask in 1_u8..15 {
                let partition = [0, 1, 2, 3].map(|bit| mask & (1 << bit) != 0);
                let realigned = realign_partition(&current, &partition, &library).unwrap();
                let optimum = realigned.candidate.expected_accuracy;
                assert!(realigned.current_score <= optimum);
                assert!(
                    quantized_objective(realigned.current_score).unwrap()
                        <= quantized_objective(optimum).unwrap()
                );
                let row_pair = cross_partition_objective(&current, &partition, &library).unwrap();
                assert!((row_pair - realigned.current_score).abs() <= 1.0e-12 * row_pair.max(1.0));
            }
        }
    }
}
