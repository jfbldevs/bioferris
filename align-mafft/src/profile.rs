use core::mem::size_of;
use core::ops::Range;
use std::time::{Duration, Instant};

use align_core::Symbol;
use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord};

use crate::matrix::ScientificScorer;
use crate::{MafftAlphabet, MafftError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ProfileOperation {
    Diagonal,
    LeftOnly,
    RightOnly,
}

#[derive(Clone, Debug)]
pub(crate) struct WorkingRow<A: MafftAlphabet> {
    pub(crate) input_index: usize,
    pub(crate) identifier: String,
    pub(crate) cells: Vec<Option<Symbol<A>>>,
    pub(crate) weight: f64,
}

#[derive(Clone, Debug)]
pub(crate) struct WorkingProfile<A: MafftAlphabet> {
    pub(crate) rows: Vec<WorkingRow<A>>,
    pub(crate) width: usize,
}

impl<A: MafftAlphabet> WorkingProfile<A> {
    pub(crate) fn leaf(input_index: usize, record: &SequenceRecord<A>, weight: f64) -> Self {
        Self {
            rows: vec![WorkingRow {
                input_index,
                identifier: record.identifier().to_owned(),
                cells: record.sequence().iter().map(Some).collect(),
                weight,
            }],
            width: record.sequence().len(),
        }
    }

    pub(crate) fn into_alignment(mut self) -> Result<MultipleAlignment<A>, MafftError> {
        self.rows.sort_by_key(|row| row.input_index);
        let rows = self
            .rows
            .into_iter()
            .map(|row| AlignedRow::new(row.identifier, row.cells))
            .collect::<Result<Vec<_>, _>>()?;
        MultipleAlignment::try_new(rows).map_err(Into::into)
    }

    pub(crate) fn member_indices(&self) -> Vec<usize> {
        let mut members = self
            .rows
            .iter()
            .map(|row| row.input_index)
            .collect::<Vec<_>>();
        members.sort_unstable();
        members
    }
}

#[derive(Clone, Debug)]
pub(crate) struct PreparedProfile {
    width: usize,
    alphabet_size: usize,
    frequencies: Vec<f64>,
    signals: Vec<Vec<f64>>,
    gap_starts: Vec<f64>,
    gap_ends: Vec<f64>,
}

impl PreparedProfile {
    pub(crate) fn new<A: MafftAlphabet>(
        profile: &WorkingProfile<A>,
        scorer: &ScientificScorer,
    ) -> Result<Self, MafftError> {
        let total_weight: f64 = profile.rows.iter().map(|row| row.weight).sum();
        if !total_weight.is_finite() || total_weight <= 0.0 {
            return Err(MafftError::NonFiniteComputation {
                stage: "profile weight normalization",
            });
        }
        let weights = profile
            .rows
            .iter()
            .map(|row| row.weight / total_weight)
            .collect::<Vec<_>>();
        let alphabet_size = scorer.alphabet_size();
        let frequency_count =
            profile
                .width
                .checked_mul(alphabet_size)
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "profile column frequencies",
                })?;
        let mut frequencies = vec![0.0; frequency_count];
        let mut signals = vec![vec![0.0; profile.width]; scorer.channel_count()];
        let boundary_count =
            profile
                .width
                .checked_add(1)
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "profile gap boundaries",
                })?;
        let mut gap_starts = vec![0.0; boundary_count];
        let mut gap_ends = vec![0.0; boundary_count];

        for (row, weight) in profile.rows.iter().zip(weights) {
            if row.cells.len() != profile.width {
                return Err(MafftError::InternalInvariant {
                    message: "working profile row is not rectangular",
                });
            }
            for (column, cell) in row.cells.iter().copied().enumerate() {
                if let Some(symbol) = cell {
                    frequencies[column * alphabet_size + symbol.alphabet_index()] += weight;
                    let signal = scorer.signal(symbol);
                    for channel in 0..scorer.channel_count() {
                        signals[channel][column] += weight * signal[channel];
                    }
                }
            }
            if profile.width > 0 {
                if row.cells[0].is_none() {
                    gap_starts[0] += weight;
                }
                if row.cells[profile.width - 1].is_none() {
                    gap_ends[profile.width] += weight;
                }
                for boundary in 1..profile.width {
                    if row.cells[boundary - 1].is_some() && row.cells[boundary].is_none() {
                        gap_starts[boundary] += weight;
                    }
                    if row.cells[boundary - 1].is_none() && row.cells[boundary].is_some() {
                        gap_ends[boundary] += weight;
                    }
                }
            }
        }
        Ok(Self {
            width: profile.width,
            alphabet_size,
            frequencies,
            signals,
            gap_starts,
            gap_ends,
        })
    }

    pub(crate) const fn width(&self) -> usize {
        self.width
    }

    pub(crate) fn signals(&self) -> &[Vec<f64>] {
        &self.signals
    }

    fn frequencies(&self, column: usize) -> &[f64] {
        let start = column * self.alphabet_size;
        &self.frequencies[start..start + self.alphabet_size]
    }

    fn gap_cost(&self, start_boundary: usize, end_boundary: usize, gap_open: f64) -> f64 {
        debug_assert!(start_boundary < end_boundary);
        gap_open * (1.0 - (self.gap_starts[start_boundary] + self.gap_ends[end_boundary]) * 0.5)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ProfilePairScores {
    alphabet_size: usize,
    left_width: usize,
    right_width: usize,
    transformed_left: Vec<f64>,
    right_offsets: Vec<usize>,
    right_indices: Vec<usize>,
    right_values: Vec<f64>,
    dense_bonus: Option<Box<[f64]>>,
}

impl ProfilePairScores {
    pub(crate) fn new(
        left: &PreparedProfile,
        right: &PreparedProfile,
        scorer: &ScientificScorer,
    ) -> Result<Self, MafftError> {
        if left.alphabet_size != right.alphabet_size || left.alphabet_size != scorer.alphabet_size()
        {
            return Err(MafftError::InternalInvariant {
                message: "profile scoring alphabet sizes differ",
            });
        }
        let alphabet_size = left.alphabet_size;
        let transformed_count =
            left.width
                .checked_mul(alphabet_size)
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "profile score transform",
                })?;
        let mut transformed_left = vec![0.0; transformed_count];
        let mut observed = Vec::with_capacity(alphabet_size);
        for column in 0..left.width {
            let frequencies = left.frequencies(column);
            observed.clear();
            observed.extend(
                frequencies
                    .iter()
                    .copied()
                    .enumerate()
                    .filter(|(_, frequency)| *frequency != 0.0),
            );
            for target in 0..alphabet_size {
                let mut value = 0.0;
                for (source, frequency) in observed.iter().copied() {
                    value += frequency * scorer.score_indices(source, target);
                }
                transformed_left[column * alphabet_size + target] = value;
            }
        }
        let mut right_offsets = Vec::with_capacity(right.width + 1);
        let mut right_indices = Vec::new();
        let mut right_values = Vec::new();
        right_offsets.push(0);
        for column in 0..right.width {
            for (index, frequency) in right.frequencies(column).iter().copied().enumerate() {
                if frequency != 0.0 {
                    right_indices.push(index);
                    right_values.push(frequency);
                }
            }
            right_offsets.push(right_indices.len());
        }
        Ok(Self {
            alphabet_size,
            left_width: left.width,
            right_width: right.width,
            transformed_left,
            right_offsets,
            right_indices,
            right_values,
            dense_bonus: None,
        })
    }

    pub(crate) fn with_dense_bonus(mut self, bonus: Vec<f64>) -> Result<Self, MafftError> {
        let expected = self.left_width.checked_mul(self.right_width);
        let expected = expected.ok_or(MafftError::ArithmeticOverflow {
            stage: "profile consistency bonus dimensions",
        })?;
        if bonus.len() != expected {
            return Err(MafftError::InternalInvariant {
                message: "profile consistency bonus dimensions differ",
            });
        }
        if bonus.iter().any(|value| !value.is_finite()) {
            return Err(MafftError::NonFiniteComputation {
                stage: "profile consistency bonus",
            });
        }
        self.dense_bonus = Some(bonus.into_boxed_slice());
        Ok(self)
    }

    #[inline]
    fn sparse_homology(&self, left_column: usize, right_column: usize) -> f64 {
        debug_assert!(left_column < self.left_width);
        let left_start = left_column * self.alphabet_size;
        let mut score = 0.0;
        for slot in self.right_offsets[right_column]..self.right_offsets[right_column + 1] {
            score += self.transformed_left[left_start + self.right_indices[slot]]
                * self.right_values[slot];
        }
        score
    }

    #[inline]
    pub(crate) fn homology_without_bonus(&self, left_column: usize, right_column: usize) -> f64 {
        debug_assert!(self.dense_bonus.is_none());
        self.sparse_homology(left_column, right_column) + 0.0
    }

    #[inline]
    pub(crate) fn homology(&self, left_column: usize, right_column: usize) -> f64 {
        self.sparse_homology(left_column, right_column)
            + self.dense_bonus.as_ref().map_or(0.0, |bonus| {
                bonus[left_column * self.right_width + right_column]
            })
    }

    fn fill_homology_row<const DENSE_BONUS: bool>(
        &self,
        left_column: usize,
        right_range: Range<usize>,
        output: &mut [f64],
    ) {
        debug_assert_eq!(right_range.len(), output.len());
        for (score, right_column) in output.iter_mut().zip(right_range) {
            *score = if DENSE_BONUS {
                self.homology(left_column, right_column)
            } else {
                self.homology_without_bonus(left_column, right_column)
            };
        }
    }

    fn record_homology_term_distribution(
        &self,
        right_range: &Range<usize>,
        left_repetitions: usize,
        timings: &mut ProfileDpPhaseTimings,
    ) {
        debug_assert!(right_range.end <= self.right_width);
        for right_column in right_range.clone() {
            let term_count = self.right_offsets[right_column + 1]
                .saturating_sub(self.right_offsets[right_column]);
            timings.homology_term_evaluations = timings
                .homology_term_evaluations
                .saturating_add(term_count.saturating_mul(left_repetitions));
            let bucket = term_count.min(timings.homology_term_cell_counts.len() - 1);
            timings.homology_term_cell_counts[bucket] =
                timings.homology_term_cell_counts[bucket].saturating_add(left_repetitions);
        }
    }

    pub(crate) fn resident_bytes(&self) -> Result<usize, MafftError> {
        let float_bytes = self
            .transformed_left
            .len()
            .checked_add(self.right_values.len())
            .and_then(|elements| {
                elements.checked_add(self.dense_bonus.as_ref().map_or(0, |bonus| bonus.len()))
            })
            .and_then(|elements| elements.checked_mul(size_of::<f64>()))
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "profile score workspace accounting",
            })?;
        let index_bytes = self
            .right_offsets
            .len()
            .checked_add(self.right_indices.len())
            .and_then(|elements| elements.checked_mul(size_of::<usize>()))
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "profile score workspace accounting",
            })?;
        float_bytes
            .checked_add(index_bytes)
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "profile score workspace accounting",
            })
    }
}

#[derive(Clone, Debug)]
pub(crate) struct BlockAlignment {
    pub(crate) score: f64,
    pub(crate) operations: Vec<ProfileOperation>,
    pub(crate) peak_bytes: usize,
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct ProfileDpPhaseTimings {
    pub(crate) block_count: usize,
    pub(crate) cell_count: usize,
    pub(crate) homology_term_evaluations: usize,
    pub(crate) homology_term_cell_counts: [usize; 6],
    pub(crate) workspace_initialization: Duration,
    pub(crate) recurrence_fill: Duration,
    pub(crate) endpoint_traceback: Duration,
    pub(crate) block_accounting: Duration,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(usize)]
enum TraceKind {
    Diagonal = 0,
    LeftRun = 1,
    RightRun = 2,
}

#[derive(Clone, Copy, Debug)]
#[repr(transparent)]
struct Trace(usize);

impl Trace {
    const KIND_MASK: usize = 0b11;

    fn new(kind: TraceKind, origin: usize) -> Self {
        debug_assert!(origin <= usize::MAX >> 2);
        Self((origin << 2) | kind as usize)
    }

    const fn kind(self) -> TraceKind {
        match self.0 & Self::KIND_MASK {
            0 => TraceKind::Diagonal,
            1 => TraceKind::LeftRun,
            2 => TraceKind::RightRun,
            _ => unreachable!(),
        }
    }

    const fn origin(self) -> usize {
        self.0 >> 2
    }
}

#[derive(Clone, Copy, Debug)]
struct PrefixMaximum {
    score: f64,
    origin: usize,
}

impl PrefixMaximum {
    const UNREACHABLE: Self = Self {
        score: f64::NEG_INFINITY,
        origin: 0,
    };
}

#[allow(clippy::too_many_lines)]
pub(crate) fn align_block(
    left: &PreparedProfile,
    right: &PreparedProfile,
    pair_scores: &ProfilePairScores,
    left_range: Range<usize>,
    right_range: Range<usize>,
    gap_open: f64,
) -> Result<BlockAlignment, MafftError> {
    align_block_with_phase_timings::<false, true>(
        left,
        right,
        pair_scores,
        left_range,
        right_range,
        gap_open,
        &mut ProfileDpPhaseTimings::default(),
    )
}

#[allow(clippy::too_many_lines)]
pub(crate) fn align_block_without_bonus(
    left: &PreparedProfile,
    right: &PreparedProfile,
    pair_scores: &ProfilePairScores,
    left_range: Range<usize>,
    right_range: Range<usize>,
    gap_open: f64,
) -> Result<BlockAlignment, MafftError> {
    align_block_with_phase_timings::<false, false>(
        left,
        right,
        pair_scores,
        left_range,
        right_range,
        gap_open,
        &mut ProfileDpPhaseTimings::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn align_block_profiled_without_bonus(
    left: &PreparedProfile,
    right: &PreparedProfile,
    pair_scores: &ProfilePairScores,
    left_range: Range<usize>,
    right_range: Range<usize>,
    gap_open: f64,
    timings: &mut ProfileDpPhaseTimings,
) -> Result<BlockAlignment, MafftError> {
    align_block_with_phase_timings::<true, false>(
        left,
        right,
        pair_scores,
        left_range,
        right_range,
        gap_open,
        timings,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn align_block_with_phase_timings<const PROFILED: bool, const DENSE_BONUS: bool>(
    left: &PreparedProfile,
    right: &PreparedProfile,
    pair_scores: &ProfilePairScores,
    left_range: Range<usize>,
    right_range: Range<usize>,
    gap_open: f64,
    timings: &mut ProfileDpPhaseTimings,
) -> Result<BlockAlignment, MafftError> {
    debug_assert!(left_range.end <= left.width());
    debug_assert!(right_range.end <= right.width());
    let left_length = left_range.len();
    let right_length = right_range.len();
    if PROFILED {
        timings.block_count = timings.block_count.saturating_add(1);
    }
    if left_length == 0 || right_length == 0 {
        let started = PROFILED.then(Instant::now);
        let block = align_empty_block(left, right, left_range, right_range, gap_open);
        finish_profile_dp_phase(started, &mut timings.endpoint_traceback);
        return block;
    }
    let cell_count =
        left_length
            .checked_mul(right_length)
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "profile DP matrix",
            })?;
    if PROFILED {
        timings.cell_count = timings.cell_count.saturating_add(cell_count);
    }
    let workspace_started = PROFILED.then(Instant::now);
    let mut two_rows_back = vec![f64::NEG_INFINITY; right_length];
    let mut previous_row = vec![f64::NEG_INFINITY; right_length];
    let mut current_row = vec![f64::NEG_INFINITY; right_length];
    let mut last_column_scores = Vec::with_capacity(left_length);
    if left_length.max(right_length) > usize::MAX >> 2 {
        return Err(MafftError::ArithmeticOverflow {
            stage: "packed profile traceback origin",
        });
    }
    let mut traces = vec![Trace::new(TraceKind::Diagonal, 0); cell_count];
    let mut vertical_prefix = vec![PrefixMaximum::UNREACHABLE; right_length.saturating_sub(1)];
    finish_profile_dp_phase(workspace_started, &mut timings.workspace_initialization);

    let recurrence_started = PROFILED.then(Instant::now);
    for left_consumed in 1..=left_length {
        if left_consumed >= 3 {
            let origin_left = left_consumed - 2;
            let start_adjustment = gap_open * 0.5 * left.gap_starts[left_range.start + origin_left];
            for origin_right in 1..right_length {
                let origin_score = two_rows_back[origin_right - 1] + start_adjustment;
                update_prefix(
                    &mut vertical_prefix[origin_right - 1],
                    origin_score,
                    origin_left,
                );
            }
        }

        let absolute_left = left_range.start + left_consumed - 1;
        let row_offset = (left_consumed - 1) * right_length;
        pair_scores.fill_homology_row::<DENSE_BONUS>(
            absolute_left,
            right_range.clone(),
            &mut current_row,
        );
        if left_consumed == 1 {
            let homology = current_row[0];
            if !homology.is_finite() {
                return Err(MafftError::NonFiniteComputation {
                    stage: "profile dynamic programming",
                });
            }
            current_row[0] = homology;
            traces[row_offset] = Trace::new(TraceKind::Diagonal, 0);

            let horizontal = PrefixMaximum {
                score: gap_open * 0.5 * right.gap_starts[right_range.start],
                origin: 0,
            };
            for right_consumed in 2..=right_length {
                let absolute_right = right_range.start + right_consumed - 1;
                let homology = current_row[right_consumed - 1];
                let best_score = horizontal.score - gap_open
                    + gap_open * 0.5 * right.gap_ends[absolute_right]
                    + homology;
                if !best_score.is_finite() {
                    return Err(MafftError::NonFiniteComputation {
                        stage: "profile dynamic programming",
                    });
                }
                current_row[right_consumed - 1] = best_score;
                traces[row_offset + right_consumed - 1] =
                    Trace::new(TraceKind::RightRun, horizontal.origin);
            }
        } else {
            let homology = current_row[0];
            let vertical = PrefixMaximum {
                score: gap_open * 0.5 * left.gap_starts[left_range.start],
                origin: 0,
            };
            let mut best_score = f64::NEG_INFINITY;
            let mut best_trace = Trace::new(TraceKind::Diagonal, 0);
            let candidate = vertical.score - gap_open
                + gap_open * 0.5 * left.gap_ends[absolute_left]
                + homology;
            if non_nan_score_is_greater(candidate, best_score) {
                best_score = candidate;
                best_trace = Trace::new(TraceKind::LeftRun, vertical.origin);
            }
            if !best_score.is_finite() {
                return Err(MafftError::NonFiniteComputation {
                    stage: "profile dynamic programming",
                });
            }
            current_row[0] = best_score;
            traces[row_offset] = best_trace;

            let mut horizontal_prefix = PrefixMaximum::UNREACHABLE;
            for right_consumed in 2..=right_length {
                if right_consumed >= 3 {
                    let origin_right = right_consumed - 2;
                    let origin_score = previous_row[origin_right - 1]
                        + gap_open * 0.5 * right.gap_starts[right_range.start + origin_right];
                    update_prefix(&mut horizontal_prefix, origin_score, origin_right);
                }

                let absolute_right = right_range.start + right_consumed - 1;
                let homology = current_row[right_consumed - 1];
                let mut best_score = previous_row[right_consumed - 2] + homology;
                let mut best_trace = Trace::new(TraceKind::Diagonal, left_consumed - 1);

                if left_consumed >= 3 {
                    let vertical = vertical_prefix[right_consumed - 2];
                    let candidate = vertical.score - gap_open
                        + gap_open * 0.5 * left.gap_ends[absolute_left]
                        + homology;
                    if non_nan_score_is_greater(candidate, best_score) {
                        best_score = candidate;
                        best_trace = Trace::new(TraceKind::LeftRun, vertical.origin);
                    }
                }
                if right_consumed >= 3 {
                    let candidate = horizontal_prefix.score - gap_open
                        + gap_open * 0.5 * right.gap_ends[absolute_right]
                        + homology;
                    if non_nan_score_is_greater(candidate, best_score) {
                        best_score = candidate;
                        best_trace = Trace::new(TraceKind::RightRun, horizontal_prefix.origin);
                    }
                }

                if !best_score.is_finite() {
                    return Err(MafftError::NonFiniteComputation {
                        stage: "profile dynamic programming",
                    });
                }
                current_row[right_consumed - 1] = best_score;
                traces[row_offset + right_consumed - 1] = best_trace;
            }
        }
        last_column_scores.push(current_row[right_length - 1]);
        core::mem::swap(&mut two_rows_back, &mut previous_row);
        core::mem::swap(&mut previous_row, &mut current_row);
    }
    finish_profile_dp_phase(recurrence_started, &mut timings.recurrence_fill);

    let traceback_started = PROFILED.then(Instant::now);
    let mut final_score = previous_row[right_length - 1];
    let mut final_left = left_length;
    let mut final_right = right_length;
    let mut suffix = None::<ProfileOperation>;
    for origin_left in 1..left_length {
        let candidate = last_column_scores[origin_left - 1]
            - left.gap_cost(left_range.start + origin_left, left_range.end, gap_open);
        if non_nan_score_is_greater(candidate, final_score) {
            final_score = candidate;
            final_left = origin_left;
            final_right = right_length;
            suffix = Some(ProfileOperation::LeftOnly);
        }
    }
    for origin_right in 1..right_length {
        let candidate = previous_row[origin_right - 1]
            - right.gap_cost(right_range.start + origin_right, right_range.end, gap_open);
        if non_nan_score_is_greater(candidate, final_score) {
            final_score = candidate;
            final_left = left_length;
            final_right = origin_right;
            suffix = Some(ProfileOperation::RightOnly);
        }
    }

    let operation_capacity =
        left_length
            .checked_add(right_length)
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "profile traceback operation capacity",
            })?;
    let mut reverse_operations = Vec::with_capacity(operation_capacity);
    match suffix {
        Some(ProfileOperation::LeftOnly) => {
            reverse_operations.extend(core::iter::repeat_n(
                ProfileOperation::LeftOnly,
                left_length - final_left,
            ));
        }
        Some(ProfileOperation::RightOnly) => {
            reverse_operations.extend(core::iter::repeat_n(
                ProfileOperation::RightOnly,
                right_length - final_right,
            ));
        }
        Some(ProfileOperation::Diagonal) | None => {}
    }

    let mut left_consumed = final_left;
    let mut right_consumed = final_right;
    while left_consumed > 0 && right_consumed > 0 {
        let trace = traces[(left_consumed - 1) * right_length + right_consumed - 1];
        reverse_operations.push(ProfileOperation::Diagonal);
        match trace.kind() {
            TraceKind::Diagonal => {
                left_consumed -= 1;
                right_consumed -= 1;
            }
            TraceKind::LeftRun => {
                reverse_operations.extend(core::iter::repeat_n(
                    ProfileOperation::LeftOnly,
                    left_consumed - trace.origin() - 1,
                ));
                left_consumed = trace.origin();
                right_consumed -= 1;
            }
            TraceKind::RightRun => {
                reverse_operations.extend(core::iter::repeat_n(
                    ProfileOperation::RightOnly,
                    right_consumed - trace.origin() - 1,
                ));
                left_consumed -= 1;
                right_consumed = trace.origin();
            }
        }
    }
    if left_consumed != 0 || right_consumed != 0 {
        return Err(MafftError::InternalInvariant {
            message: "profile traceback did not reach the origin",
        });
    }
    reverse_operations.reverse();
    validate_operation_consumption(&reverse_operations, left_length, right_length)?;
    finish_profile_dp_phase(traceback_started, &mut timings.endpoint_traceback);

    let accounting_started = PROFILED.then(Instant::now);
    let trace_bytes =
        cell_count
            .checked_mul(size_of::<Trace>())
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "profile DP memory accounting",
            })?;
    let score_bytes = right_length
        .checked_mul(3)
        .and_then(|elements| elements.checked_add(left_length))
        .and_then(|elements| elements.checked_mul(size_of::<f64>()))
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "profile DP memory accounting",
        })?;
    let prefix_bytes = right_length
        .saturating_sub(1)
        .checked_mul(size_of::<PrefixMaximum>())
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "profile DP memory accounting",
        })?;
    let peak_bytes = trace_bytes
        .checked_add(score_bytes)
        .and_then(|bytes| bytes.checked_add(prefix_bytes))
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "profile DP memory accounting",
        })?;
    if PROFILED {
        pair_scores.record_homology_term_distribution(&right_range, left_length, timings);
    }
    finish_profile_dp_phase(accounting_started, &mut timings.block_accounting);
    Ok(BlockAlignment {
        score: final_score,
        operations: reverse_operations,
        peak_bytes,
    })
}

#[inline]
fn profile_homology<const DENSE_BONUS: bool>(
    pair_scores: &ProfilePairScores,
    left_column: usize,
    right_column: usize,
) -> f64 {
    if DENSE_BONUS {
        pair_scores.homology(left_column, right_column)
    } else {
        pair_scores.homology_without_bonus(left_column, right_column)
    }
}

#[inline]
fn finish_profile_dp_phase(started: Option<Instant>, elapsed: &mut Duration) {
    if let Some(started) = started {
        *elapsed = elapsed.saturating_add(started.elapsed());
    }
}

fn align_empty_block(
    left: &PreparedProfile,
    right: &PreparedProfile,
    left_range: Range<usize>,
    right_range: Range<usize>,
    gap_open: f64,
) -> Result<BlockAlignment, MafftError> {
    let left_length = left_range.len();
    let right_length = right_range.len();
    let (score, operation, count) = match (left_length, right_length) {
        (0, 0) => (0.0, ProfileOperation::Diagonal, 0),
        (_, 0) => (
            -left.gap_cost(left_range.start, left_range.end, gap_open),
            ProfileOperation::LeftOnly,
            left_length,
        ),
        (0, _) => (
            -right.gap_cost(right_range.start, right_range.end, gap_open),
            ProfileOperation::RightOnly,
            right_length,
        ),
        _ => unreachable!(),
    };
    if !score.is_finite() {
        return Err(MafftError::NonFiniteComputation {
            stage: "empty profile block alignment",
        });
    }
    Ok(BlockAlignment {
        score,
        operations: core::iter::repeat_n(operation, count).collect(),
        peak_bytes: 0,
    })
}

fn update_prefix(slot: &mut PrefixMaximum, score: f64, origin: usize) {
    if non_nan_score_is_greater(score, slot.score)
        || (score.to_bits() == slot.score.to_bits() && origin > slot.origin)
    {
        *slot = PrefixMaximum { score, origin };
    }
}

#[inline]
fn non_nan_score_is_greater(candidate: f64, current: f64) -> bool {
    debug_assert!(!candidate.is_nan());
    debug_assert!(!current.is_nan());
    // For non-NaN values, total order differs from `>` only for signed zero.
    candidate > current || (candidate.to_bits() == 0 && current.to_bits() == (-0.0_f64).to_bits())
}

pub(crate) fn merge_profiles<A: MafftAlphabet>(
    left: WorkingProfile<A>,
    right: WorkingProfile<A>,
    operations: &[ProfileOperation],
) -> Result<WorkingProfile<A>, MafftError> {
    validate_operation_consumption(operations, left.width, right.width)?;
    let mut rows = Vec::with_capacity(left.rows.len() + right.rows.len());
    for row in left.rows {
        let mut source = 0usize;
        let mut cells = Vec::with_capacity(operations.len());
        for operation in operations {
            match operation {
                ProfileOperation::Diagonal | ProfileOperation::LeftOnly => {
                    cells.push(row.cells[source]);
                    source += 1;
                }
                ProfileOperation::RightOnly => cells.push(None),
            }
        }
        rows.push(WorkingRow {
            input_index: row.input_index,
            identifier: row.identifier,
            cells,
            weight: row.weight,
        });
    }
    for row in right.rows {
        let mut source = 0usize;
        let mut cells = Vec::with_capacity(operations.len());
        for operation in operations {
            match operation {
                ProfileOperation::Diagonal | ProfileOperation::RightOnly => {
                    cells.push(row.cells[source]);
                    source += 1;
                }
                ProfileOperation::LeftOnly => cells.push(None),
            }
        }
        rows.push(WorkingRow {
            input_index: row.input_index,
            identifier: row.identifier,
            cells,
            weight: row.weight,
        });
    }
    Ok(WorkingProfile {
        rows,
        width: operations.len(),
    })
}

pub(crate) fn replay_operations_score(
    left: &PreparedProfile,
    right: &PreparedProfile,
    pair_scores: &ProfilePairScores,
    left_range: Range<usize>,
    right_range: Range<usize>,
    operations: &[ProfileOperation],
    gap_open: f64,
) -> Result<f64, MafftError> {
    replay_operations_score_with_bonus::<true>(
        left,
        right,
        pair_scores,
        left_range,
        right_range,
        operations,
        gap_open,
    )
}

pub(crate) fn replay_operations_score_without_bonus(
    left: &PreparedProfile,
    right: &PreparedProfile,
    pair_scores: &ProfilePairScores,
    left_range: Range<usize>,
    right_range: Range<usize>,
    operations: &[ProfileOperation],
    gap_open: f64,
) -> Result<f64, MafftError> {
    replay_operations_score_with_bonus::<false>(
        left,
        right,
        pair_scores,
        left_range,
        right_range,
        operations,
        gap_open,
    )
}

#[allow(clippy::too_many_arguments)]
fn replay_operations_score_with_bonus<const DENSE_BONUS: bool>(
    left: &PreparedProfile,
    right: &PreparedProfile,
    pair_scores: &ProfilePairScores,
    left_range: Range<usize>,
    right_range: Range<usize>,
    operations: &[ProfileOperation],
    gap_open: f64,
) -> Result<f64, MafftError> {
    validate_operation_consumption(operations, left_range.len(), right_range.len())?;
    let mut left_column = left_range.start;
    let mut right_column = right_range.start;
    let mut score = 0.0;
    let mut operation_index = 0usize;
    while operation_index < operations.len() {
        match operations[operation_index] {
            ProfileOperation::Diagonal => {
                score += profile_homology::<DENSE_BONUS>(pair_scores, left_column, right_column);
                left_column += 1;
                right_column += 1;
                operation_index += 1;
            }
            ProfileOperation::LeftOnly => {
                let start = left_column;
                while operation_index < operations.len()
                    && operations[operation_index] == ProfileOperation::LeftOnly
                {
                    left_column += 1;
                    operation_index += 1;
                }
                score -= left.gap_cost(start, left_column, gap_open);
            }
            ProfileOperation::RightOnly => {
                let start = right_column;
                while operation_index < operations.len()
                    && operations[operation_index] == ProfileOperation::RightOnly
                {
                    right_column += 1;
                    operation_index += 1;
                }
                score -= right.gap_cost(start, right_column, gap_open);
            }
        }
    }
    if left_column != left_range.end || right_column != right_range.end {
        return Err(MafftError::InternalInvariant {
            message: "score replay did not consume its profile ranges",
        });
    }
    if !score.is_finite() {
        return Err(MafftError::NonFiniteComputation {
            stage: "profile operation score replay",
        });
    }
    Ok(score)
}

fn validate_operation_consumption(
    operations: &[ProfileOperation],
    expected_left: usize,
    expected_right: usize,
) -> Result<(), MafftError> {
    let left = operations
        .iter()
        .filter(|operation| !matches!(operation, ProfileOperation::RightOnly))
        .count();
    let right = operations
        .iter()
        .filter(|operation| !matches!(operation, ProfileOperation::LeftOnly))
        .count();
    if left != expected_left || right != expected_right {
        return Err(MafftError::InternalInvariant {
            message: "profile operation path consumes the wrong number of columns",
        });
    }
    Ok(())
}

#[cfg(test)]
pub(crate) fn cubic_block_score(
    left: &PreparedProfile,
    right: &PreparedProfile,
    pair_scores: &ProfilePairScores,
    left_range: Range<usize>,
    right_range: Range<usize>,
    gap_open: f64,
) -> f64 {
    let n = left_range.len();
    let m = right_range.len();
    if n == 0 || m == 0 {
        return align_empty_block(left, right, left_range, right_range, gap_open)
            .unwrap()
            .score;
    }
    let mut scores = vec![f64::NEG_INFINITY; n * m];
    for i in 1..=n {
        for j in 1..=m {
            let homology =
                pair_scores.homology(left_range.start + i - 1, right_range.start + j - 1);
            let mut best = if i == 1 && j == 1 {
                0.0
            } else if i > 1 && j > 1 {
                scores[(i - 2) * m + j - 2]
            } else {
                f64::NEG_INFINITY
            };
            if j == 1 && i >= 2 {
                best =
                    best.max(-left.gap_cost(left_range.start, left_range.start + i - 1, gap_open));
            } else if j >= 2 {
                for x in 1..i.saturating_sub(1) {
                    best = best.max(
                        scores[(x - 1) * m + j - 2]
                            - left.gap_cost(
                                left_range.start + x,
                                left_range.start + i - 1,
                                gap_open,
                            ),
                    );
                }
            }
            if i == 1 && j >= 2 {
                best = best.max(-right.gap_cost(
                    right_range.start,
                    right_range.start + j - 1,
                    gap_open,
                ));
            } else if i >= 2 {
                for y in 1..j.saturating_sub(1) {
                    best = best.max(
                        scores[(i - 2) * m + y - 1]
                            - right.gap_cost(
                                right_range.start + y,
                                right_range.start + j - 1,
                                gap_open,
                            ),
                    );
                }
            }
            scores[(i - 1) * m + j - 1] = best + homology;
        }
    }
    let mut best = scores[n * m - 1];
    for x in 1..n {
        best = best.max(
            scores[(x - 1) * m + m - 1]
                - left.gap_cost(left_range.start + x, left_range.end, gap_open),
        );
    }
    for y in 1..m {
        best = best.max(
            scores[(n - 1) * m + y - 1]
                - right.gap_cost(right_range.start + y, right_range.end, gap_open),
        );
    }
    best
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence};
    use align_msa::SequenceRecord;

    use super::{
        PreparedProfile, ProfilePairScores, WorkingProfile, align_block, align_block_without_bonus,
        cubic_block_score,
    };
    use crate::FftNs1Parameters;
    use crate::matrix::ScientificScorer;

    #[test]
    fn non_nan_score_comparison_matches_total_order() {
        let values = [
            f64::NEG_INFINITY,
            -f64::MAX,
            -1.0,
            -f64::MIN_POSITIVE,
            -f64::from_bits(1),
            -0.0,
            0.0,
            f64::from_bits(1),
            f64::MIN_POSITIVE,
            1.0,
            f64::MAX,
            f64::INFINITY,
        ];
        for candidate in values {
            for current in values {
                assert_eq!(
                    super::non_nan_score_is_greater(candidate, current),
                    candidate.total_cmp(&current).is_gt(),
                    "{candidate:?} versus {current:?}",
                );
            }
        }
    }

    fn profile(index: usize, sequence: &str) -> WorkingProfile<Dna> {
        let record = SequenceRecord::new(
            format!("s{index}"),
            DnaSequence::from_ascii(sequence).unwrap(),
        )
        .unwrap();
        WorkingProfile::leaf(index, &record, 1.0)
    }

    #[test]
    fn optimized_recurrence_matches_independent_cubic_oracle() {
        let parameters = FftNs1Parameters::paper_2002();
        let scorer = ScientificScorer::new::<Dna>(parameters).unwrap();
        let values = ["", "A", "C", "AA", "AC", "AGT", "AAAA"];
        for first in values {
            for second in values {
                let first_profile = profile(0, first);
                let second_profile = profile(1, second);
                let prepared_first = PreparedProfile::new(&first_profile, &scorer).unwrap();
                let prepared_second = PreparedProfile::new(&second_profile, &scorer).unwrap();
                let pair_scores =
                    ProfilePairScores::new(&prepared_first, &prepared_second, &scorer).unwrap();
                let observed = align_block(
                    &prepared_first,
                    &prepared_second,
                    &pair_scores,
                    0..first.len(),
                    0..second.len(),
                    parameters.gap_open(),
                )
                .unwrap();
                let expected = cubic_block_score(
                    &prepared_first,
                    &prepared_second,
                    &pair_scores,
                    0..first.len(),
                    0..second.len(),
                    parameters.gap_open(),
                );
                assert!(
                    (observed.score - expected).abs() < 1.0e-10,
                    "{first:?} vs {second:?}: {} != {expected}",
                    observed.score
                );
            }
        }
    }

    #[test]
    fn one_gap_run_costs_one_opening_without_extension() {
        let parameters = FftNs1Parameters::paper_2002();
        let scorer = ScientificScorer::new::<Dna>(parameters).unwrap();
        let first = profile(0, "AAAA");
        let second = profile(1, "A");
        let prepared_first = PreparedProfile::new(&first, &scorer).unwrap();
        let prepared_second = PreparedProfile::new(&second, &scorer).unwrap();
        let pair_scores =
            ProfilePairScores::new(&prepared_first, &prepared_second, &scorer).unwrap();
        let result = align_block(
            &prepared_first,
            &prepared_second,
            &pair_scores,
            0..4,
            0..1,
            parameters.gap_open(),
        )
        .unwrap();
        assert_eq!(result.operations.len(), 4);
        assert!(
            (result.score
                - (scorer.score(
                    DnaSequence::from_ascii("A").unwrap().get(0).unwrap(),
                    DnaSequence::from_ascii("A").unwrap().get(0).unwrap(),
                ) - 2.4))
                .abs()
                < 1.0e-10
        );
    }

    #[test]
    fn bonus_free_dp_specialization_matches_the_general_path_bit_for_bit() {
        let parameters = FftNs1Parameters::paper_2002();
        let scorer = ScientificScorer::new::<Dna>(parameters).unwrap();
        let first = profile(0, "ACGTNRY");
        let second = profile(1, "AGTNRY");
        let prepared_first = PreparedProfile::new(&first, &scorer).unwrap();
        let prepared_second = PreparedProfile::new(&second, &scorer).unwrap();
        let pair_scores =
            ProfilePairScores::new(&prepared_first, &prepared_second, &scorer).unwrap();
        let general = align_block(
            &prepared_first,
            &prepared_second,
            &pair_scores,
            0..prepared_first.width(),
            0..prepared_second.width(),
            parameters.gap_open(),
        )
        .unwrap();
        let specialized = align_block_without_bonus(
            &prepared_first,
            &prepared_second,
            &pair_scores,
            0..prepared_first.width(),
            0..prepared_second.width(),
            parameters.gap_open(),
        )
        .unwrap();
        assert_eq!(specialized.score.to_bits(), general.score.to_bits());
        assert_eq!(specialized.operations, general.operations);
        assert_eq!(specialized.peak_bytes, general.peak_bytes);
    }

    #[test]
    fn sparse_homology_matches_dense_equation_bit_for_bit() {
        let parameters = FftNs1Parameters::paper_2002();
        let scorer = ScientificScorer::new::<Dna>(parameters).unwrap();
        let left = profile(0, "ACGTNRY");
        let right = profile(1, "TGCANRY");
        let prepared_left = PreparedProfile::new(&left, &scorer).unwrap();
        let prepared_right = PreparedProfile::new(&right, &scorer).unwrap();
        let pair_scores = ProfilePairScores::new(&prepared_left, &prepared_right, &scorer).unwrap();

        for left_column in 0..prepared_left.width() {
            for right_column in 0..prepared_right.width() {
                let mut expected = 0.0;
                for target in 0..scorer.alphabet_size() {
                    let mut transformed = 0.0;
                    for (source, frequency) in prepared_left
                        .frequencies(left_column)
                        .iter()
                        .copied()
                        .enumerate()
                    {
                        transformed += frequency * scorer.score_indices(source, target);
                    }
                    expected += transformed * prepared_right.frequencies(right_column)[target];
                }
                assert_eq!(
                    pair_scores.homology(left_column, right_column).to_bits(),
                    expected.to_bits()
                );
            }
        }
    }
}
