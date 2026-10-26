//! Compact affine WFA score and inferred-traceback kernels modeled on WFA2's
//! padded source layout.
#![allow(clippy::inline_always)]

use core::mem::size_of;

use align_core::{AlignmentOperation, Score};
use align_fast_kernels::WavefrontExtension;
use pulp::{Simd, WithSimd};

use crate::{
    WavefrontAffinePenalties, WavefrontError, WavefrontTraceRun, wavefront::count_byte_mismatches,
};

// Keep one billion negative values available for branch-free predecessor
// increments, matching WFA2's score-kernel sentinel convention.
const UNREACHABLE: i32 = i32::MIN / 2;
const MAXIMUM_COMPACT_INCREMENT: Score = 63;
const MAXIMUM_PROFITABLE_COMPACT_INCREMENT: Score = 18;
const MAXIMUM_I32_LANES: usize = 16;
// Exact bounds constrain both sides of each front and need more work to recover
// compact-kernel setup than unconstrained public-API calls.
const MINIMUM_BOUNDED_COMPACT_LEVELS: Score = 29;
const MINIMUM_UNBOUNDED_COMPACT_LEVELS: Score = 20;
// Reuse the measured finite-bound crossover for compact traceback; the exact
// bound keeps these low-score histories narrow.
const MINIMUM_BOUNDED_COMPACT_TRACE_LEVELS: Score = 29;
// Histories below the compact traceback crossover do not retain enough padded
// wavefronts to amortize preserving their buffer-to-level layout across calls.
const MINIMUM_ORDERED_TRACE_WAVEFRONTS: usize = 29;
// Avoid adding a mismatch-count pass to tiny unbounded public-API calls.
const MINIMUM_UNBOUNDED_LENGTH: usize = 128;

#[derive(Debug, Default)]
pub(crate) struct CompactScoreWorkspace {
    slots: Vec<Option<CompactWavefront>>,
    recycled: Vec<CompactWavefront>,
    null_offsets: Vec<i32>,
    trace: CompactTraceWorkspace,
    architecture: Option<pulp::Arch>,
    extension: Option<WavefrontExtension>,
}

impl CompactScoreWorkspace {
    pub(crate) const fn new() -> Self {
        Self {
            slots: Vec::new(),
            recycled: Vec::new(),
            null_offsets: Vec::new(),
            trace: CompactTraceWorkspace::new(),
            architecture: None,
            extension: None,
        }
    }

    pub(crate) fn release(&mut self) {
        self.slots = Vec::new();
        self.recycled = Vec::new();
        self.null_offsets = Vec::new();
        self.trace = CompactTraceWorkspace::new();
    }

    pub(crate) fn retained_bytes(&self) -> usize {
        let wavefronts = self
            .slots
            .iter()
            .filter_map(Option::as_ref)
            .chain(&self.recycled)
            .map(|wavefront| wavefront.offsets.capacity() * size_of::<i32>())
            .sum::<usize>();
        wavefronts
            + self.slots.capacity() * size_of::<Option<CompactWavefront>>()
            + self.recycled.capacity() * size_of::<CompactWavefront>()
            + self.null_offsets.capacity() * size_of::<i32>()
            + self.trace.retained_bytes()
    }
}

#[derive(Clone, Copy)]
struct CompactDimensions {
    first: i32,
    second: i32,
    target_diagonal: i32,
}

impl CompactDimensions {
    fn new(first_len: usize, second_len: usize, padding: usize) -> Option<Self> {
        let first = i32::try_from(first_len).ok()?;
        let second = i32::try_from(second_len).ok()?;
        let padding = i32::try_from(padding).ok()?;
        first.checked_add(padding)?;
        second.checked_add(padding)?;
        Some(Self {
            first,
            second,
            target_diagonal: second.checked_sub(first)?,
        })
    }
}

#[derive(Debug)]
struct CompactWavefront {
    penalty: Score,
    lower: i32,
    width: usize,
    storage_lower: i32,
    storage_width: usize,
    padding: usize,
    offsets: Vec<i32>,
    main_active: Option<(i32, i32)>,
    insertion_active: Option<(i32, i32)>,
    deletion_active: Option<(i32, i32)>,
}

impl CompactWavefront {
    fn new(
        penalty: Score,
        lower: i32,
        width: usize,
        padding: usize,
        first_len: usize,
        second_len: usize,
    ) -> Result<Self, WavefrontError> {
        let mut wavefront = Self {
            penalty,
            lower,
            width: 0,
            storage_lower: lower,
            storage_width: 0,
            padding,
            offsets: Vec::new(),
            main_active: None,
            insertion_active: None,
            deletion_active: None,
        };
        wavefront.reset::<false>(penalty, lower, width, padding, first_len, second_len)?;
        Ok(wavefront)
    }

    #[allow(clippy::too_many_arguments)]
    fn reset<const ALLOW_PADDING_REUSE: bool>(
        &mut self,
        penalty: Score,
        lower: i32,
        width: usize,
        padding: usize,
        first_len: usize,
        second_len: usize,
    ) -> Result<(), WavefrontError> {
        let storage_width = width
            .checked_add(
                padding
                    .checked_mul(2)
                    .ok_or(WavefrontError::DimensionOverflow {
                        first_len,
                        second_len,
                    })?,
            )
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let state_length =
            storage_width
                .checked_mul(3)
                .ok_or(WavefrontError::DimensionOverflow {
                    first_len,
                    second_len,
                })?;
        let padding_i32 =
            i32::try_from(padding).map_err(|_| WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let storage_lower =
            lower
                .checked_sub(padding_i32)
                .ok_or(WavefrontError::DimensionOverflow {
                    first_len,
                    second_len,
                })?;
        let layout_matches = ALLOW_PADDING_REUSE
            && self.width == width
            && self.storage_width == storage_width
            && self.padding == padding
            && self.offsets.len() == state_length;

        if self.offsets.capacity() < state_length {
            self.offsets
                .try_reserve(state_length - self.offsets.len())
                .map_err(|_| WavefrontError::AllocationFailed {
                    buffer: "compact affine WFA offsets",
                    elements: state_length,
                })?;
        }
        self.offsets.resize(state_length, UNREACHABLE);
        self.offsets.truncate(state_length);
        if layout_matches {
            debug_assert!((0..3).all(|state| {
                let base = state * storage_width;
                self.offsets[base..base + padding]
                    .iter()
                    .chain(&self.offsets[base + padding + width..base + storage_width])
                    .all(|offset| *offset == UNREACHABLE)
            }));
        } else {
            for state in 0..3 {
                let base = state * storage_width;
                self.offsets[base..base + padding].fill(UNREACHABLE);
                self.offsets[base + padding + width..base + storage_width].fill(UNREACHABLE);
            }
        }

        self.penalty = penalty;
        self.lower = lower;
        self.width = width;
        self.storage_lower = storage_lower;
        self.storage_width = storage_width;
        self.padding = padding;
        self.main_active = None;
        self.insertion_active = None;
        self.deletion_active = None;
        Ok(())
    }

    fn index(&self, diagonal: i32) -> Option<usize> {
        let relative = diagonal.checked_sub(self.lower)?;
        let index = usize::try_from(relative).ok()?;
        (index < self.width).then_some(self.padding + index)
    }

    fn main(&self, diagonal: i32) -> i32 {
        self.index(diagonal)
            .and_then(|index| self.offsets.get(index))
            .copied()
            .unwrap_or(UNREACHABLE)
    }

    fn insertion(&self, diagonal: i32) -> i32 {
        self.index(diagonal)
            .and_then(|index| self.offsets.get(self.storage_width + index))
            .copied()
            .unwrap_or(UNREACHABLE)
    }

    fn deletion(&self, diagonal: i32) -> i32 {
        self.index(diagonal)
            .and_then(|index| self.offsets.get(2 * self.storage_width + index))
            .copied()
            .unwrap_or(UNREACHABLE)
    }

    fn reaches_target(&self, dimensions: CompactDimensions) -> bool {
        self.main(dimensions.target_diagonal) >= dimensions.second
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
enum CompactMainTrace {
    Start = 0,
    Substitution = 1,
    Deletion = 2,
    Insertion = 3,
}

#[derive(Debug)]
struct CompactTraceWavefront {
    wavefront: CompactWavefront,
}

#[derive(Debug, Default)]
struct CompactTraceWorkspace {
    history: Vec<CompactTraceWavefront>,
    recycled_wavefronts: Vec<CompactWavefront>,
    reuse_padding: bool,
}

impl CompactTraceWorkspace {
    const fn new() -> Self {
        Self {
            history: Vec::new(),
            recycled_wavefronts: Vec::new(),
            reuse_padding: false,
        }
    }

    fn recycle_history(&mut self) {
        if !self.history.is_empty() {
            self.reuse_padding =
                recycle_trace_records(&mut self.history, &mut self.recycled_wavefronts);
        }
    }

    fn retained_bytes(&self) -> usize {
        let history = self
            .history
            .iter()
            .map(|record| record.wavefront.offsets.capacity() * size_of::<i32>())
            .sum::<usize>();
        let recycled_wavefronts = self
            .recycled_wavefronts
            .iter()
            .map(|wavefront| wavefront.offsets.capacity() * size_of::<i32>())
            .sum::<usize>();
        history
            + recycled_wavefronts
            + self.history.capacity() * size_of::<CompactTraceWavefront>()
            + self.recycled_wavefronts.capacity() * size_of::<CompactWavefront>()
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct CompactTracebackSummary {
    pub(crate) penalty: Score,
    pub(crate) wavefronts: usize,
    pub(crate) computed_offsets: usize,
    pub(crate) source_wavefront_lookups: usize,
    pub(crate) source_offset_lookups: usize,
    pub(crate) match_extensions: usize,
    pub(crate) maximum_wavefront_width: usize,
    pub(crate) modeled_peak_bytes: usize,
}

#[derive(Clone, Copy)]
enum State {
    Main,
    Insertion,
    Deletion,
}

enum CompactStep {
    Unsupported,
    Empty,
    Wavefront(CompactWavefront, usize),
}

struct HasVectorLanes;

impl WithSimd for HasVectorLanes {
    type Output = bool;

    #[inline(always)]
    fn with_simd<V: Simd>(self, _simd: V) -> Self::Output {
        !V::IS_SCALAR && V::I32_LANES > 1
    }
}

#[inline(always)]
fn compact_dispatch_bound(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    upper_bound: Score,
) -> Option<Score> {
    let (estimated_bound, minimum_levels) = if upper_bound == Score::MAX {
        if first.len() != second.len() || first.len() < MINIMUM_UNBOUNDED_LENGTH {
            return None;
        }
        let mismatches = Score::try_from(count_byte_mismatches(first, second)).ok()?;
        (
            mismatches
                .checked_mul(penalties.mismatch())
                .unwrap_or(Score::MAX),
            MINIMUM_UNBOUNDED_COMPACT_LEVELS,
        )
    } else {
        (upper_bound, MINIMUM_BOUNDED_COMPACT_LEVELS)
    };
    let minimum_increment = penalties.mismatch().min(penalties.gap_extension());
    let levels =
        estimated_bound / minimum_increment + Score::from(estimated_bound % minimum_increment != 0);
    (levels >= minimum_levels).then_some(estimated_bound)
}

#[inline(always)]
fn compact_traceback_dispatch_bound(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    upper_bound: Score,
) -> Option<Score> {
    if upper_bound == Score::MAX {
        return compact_dispatch_bound(first, second, penalties, upper_bound);
    }
    let minimum_increment = penalties.mismatch().min(penalties.gap_extension());
    let levels =
        upper_bound / minimum_increment + Score::from(upper_bound % minimum_increment != 0);
    (levels >= MINIMUM_BOUNDED_COMPACT_TRACE_LEVELS).then_some(upper_bound)
}

#[inline(always)]
fn compact_configuration(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    upper_bound: Score,
    workspace: &mut CompactScoreWorkspace,
) -> Option<(pulp::Arch, WavefrontExtension, usize)> {
    let padding = compact_padding(penalties)?;
    compact_dispatch_bound(first, second, penalties, upper_bound)?;
    let architecture = *workspace.architecture.get_or_insert_with(pulp::Arch::new);
    let extension = *workspace
        .extension
        .get_or_insert_with(WavefrontExtension::new);
    architecture
        .dispatch(HasVectorLanes)
        .then_some((architecture, extension, padding))
}

#[inline(always)]
fn compact_padding(penalties: WavefrontAffinePenalties) -> Option<usize> {
    let maximum_increment = penalties
        .mismatch()
        .max(
            penalties
                .gap_opening()
                .checked_add(penalties.gap_extension())
                .unwrap_or(Score::MAX),
        )
        .max(penalties.gap_extension());
    (maximum_increment <= MAXIMUM_PROFITABLE_COMPACT_INCREMENT)
        .then(|| usize::try_from(maximum_increment + 1).ok())
        .flatten()
}

struct CompactPenaltyScheduler {
    slots: u64,
    cursor: Score,
    overflowed: bool,
}

impl CompactPenaltyScheduler {
    const fn new() -> Self {
        Self {
            slots: 0,
            cursor: 0,
            overflowed: false,
        }
    }

    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn schedule(&mut self, penalty: Score, increment: Score) {
        let Some(next) = penalty.checked_add(increment) else {
            self.overflowed = true;
            return;
        };
        debug_assert!(next > self.cursor);
        debug_assert!(next - self.cursor <= MAXIMUM_COMPACT_INCREMENT);
        self.slots |= 1_u64 << (next as usize & 63);
    }

    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn pop(&mut self) -> Option<Score> {
        if self.slots == 0 {
            return None;
        }
        let cursor_slot = self.cursor as usize & 63;
        let rotation =
            u32::try_from((cursor_slot + 1) & 63).expect("a compact scheduler rotation fits u32");
        let delta = Score::from(self.slots.rotate_right(rotation).trailing_zeros()) + 1;
        let Some(next) = self.cursor.checked_add(delta) else {
            self.overflowed = true;
            return None;
        };
        debug_assert!(delta <= MAXIMUM_COMPACT_INCREMENT);
        self.cursor = next;
        self.slots &= !(1_u64 << (next as usize & 63));
        Some(next)
    }
}

fn checked_width(
    lower: i32,
    upper: i32,
    first_len: usize,
    second_len: usize,
) -> Result<usize, WavefrontError> {
    if lower > upper {
        return Ok(0);
    }
    upper
        .checked_sub(lower)
        .and_then(|span| span.checked_add(1))
        .and_then(|span| usize::try_from(span).ok())
        .ok_or(WavefrontError::DimensionOverflow {
            first_len,
            second_len,
        })
}

fn state_slice<'a>(
    source: Option<&'a CompactWavefront>,
    state: State,
    first_diagonal: i32,
    length: usize,
    null_offsets: &'a [i32],
) -> Option<&'a [i32]> {
    let Some(source) = source else {
        return null_offsets.get(..length);
    };
    let relative = first_diagonal.checked_sub(source.storage_lower)?;
    let relative = usize::try_from(relative).ok()?;
    let state_offset = match state {
        State::Main => 0,
        State::Insertion => source.storage_width,
        State::Deletion => source.storage_width.checked_mul(2)?,
    };
    let start = state_offset.checked_add(relative)?;
    let end = start.checked_add(length)?;
    source.offsets.get(start..end)
}

fn source_wavefront(
    slots: &[Option<CompactWavefront>],
    current_penalty: Score,
    decrement: Score,
) -> Option<&CompactWavefront> {
    let source_penalty = current_penalty.checked_sub(decrement)?;
    if source_penalty < 0 || slots.is_empty() {
        return None;
    }
    let slot_count = Score::try_from(slots.len()).ok()?;
    let index = usize::try_from(source_penalty % slot_count).ok()?;
    slots
        .get(index)?
        .as_ref()
        .filter(|wavefront| wavefront.penalty == source_penalty)
}

fn trace_source_wavefront<'a>(
    history: &'a [CompactTraceWavefront],
    indices: &[usize; 64],
    slot_count: usize,
    current_penalty: Score,
    decrement: Score,
) -> Option<&'a CompactTraceWavefront> {
    let source_penalty = current_penalty.checked_sub(decrement)?;
    if source_penalty < 0 {
        return None;
    }
    let slot_count = Score::try_from(slot_count).ok()?;
    let slot = usize::try_from(source_penalty % slot_count).ok()?;
    history
        .get(*indices.get(slot)?)
        .filter(|record| record.wavefront.penalty == source_penalty)
}

#[inline(always)]
fn compact_main_trace(substitution: i32, deletion: i32, insertion: i32) -> (i32, CompactMainTrace) {
    if substitution >= deletion && substitution >= insertion {
        (substitution, CompactMainTrace::Substitution)
    } else if deletion >= insertion {
        (deletion, CompactMainTrace::Deletion)
    } else {
        (insertion, CompactMainTrace::Insertion)
    }
}

fn merge_shifted_interval(
    lower: &mut i32,
    upper: &mut i32,
    interval: Option<(i32, i32)>,
    shift: i32,
    dimensions: CompactDimensions,
) {
    let Some((source_lower, source_upper)) = interval else {
        return;
    };
    let shifted_lower = source_lower.saturating_add(shift).max(-dimensions.first);
    let shifted_upper = source_upper.saturating_add(shift).min(dimensions.second);
    if shifted_lower <= shifted_upper {
        *lower = (*lower).min(shifted_lower);
        *upper = (*upper).max(shifted_upper);
    }
}

fn active_interval_from_bounds(lower: i32, upper: i32) -> Option<(i32, i32)> {
    (lower <= upper).then_some((lower, upper))
}

fn intersect_interval(interval: Option<(i32, i32)>, lower: i32, upper: i32) -> Option<(i32, i32)> {
    let (interval_lower, interval_upper) = interval?;
    active_interval_from_bounds(interval_lower.max(lower), interval_upper.min(upper))
}

fn constrain_diagonal_bounds(
    lower: &mut i32,
    upper: &mut i32,
    penalty: Score,
    upper_bound: Score,
    gap_extension: Score,
    dimensions: CompactDimensions,
) {
    if upper_bound == Score::MAX {
        return;
    }
    let remaining = upper_bound - penalty;
    if remaining < 0 {
        *lower = 1;
        *upper = 0;
        return;
    }
    let maximum_shift = i32::try_from(remaining / gap_extension).unwrap_or(i32::MAX);
    *lower = (*lower).max(dimensions.target_diagonal.saturating_sub(maximum_shift));
    *upper = (*upper).min(dimensions.target_diagonal.saturating_add(maximum_shift));
}

#[allow(clippy::inline_always)]
#[inline(always)]
fn extend_equal(first: &[u8], second: &[u8], diagonal: i32, offset: i32) -> i32 {
    let mut horizontal = usize::try_from(offset).expect("a legal compact offset is nonnegative");
    let mut vertical = usize::try_from(offset - diagonal)
        .expect("a legal compact vertical coordinate is nonnegative");
    let remaining = (first.len() - vertical).min(second.len() - horizontal);
    let word_bytes = remaining - remaining % size_of::<u64>();
    let first_words = first[vertical..vertical + word_bytes].chunks_exact(size_of::<u64>());
    let second_words = second[horizontal..horizontal + word_bytes].chunks_exact(size_of::<u64>());
    for (word_index, (first_word, second_word)) in first_words.zip(second_words).enumerate() {
        let first_word = u64::from_le_bytes(
            first_word
                .try_into()
                .expect("an exact compact WFA word chunk has eight bytes"),
        );
        let second_word = u64::from_le_bytes(
            second_word
                .try_into()
                .expect("an exact compact WFA word chunk has eight bytes"),
        );
        let difference = first_word ^ second_word;
        if difference != 0 {
            let equal = usize::try_from(difference.trailing_zeros() / u8::BITS)
                .expect("a compact word prefix fits usize");
            return i32::try_from(horizontal + word_index * size_of::<u64>() + equal)
                .expect("compact WFA dimensions fit i32");
        }
    }
    vertical += word_bytes;
    horizontal += word_bytes;
    while vertical < first.len()
        && horizontal < second.len()
        && first[vertical] == second[horizontal]
    {
        vertical += 1;
        horizontal += 1;
    }
    i32::try_from(horizontal).expect("compact WFA dimensions fit i32")
}

fn active_interval(offsets: &[i32], lower: i32) -> Option<(i32, i32)> {
    let first = offsets.iter().position(|offset| *offset >= 0)?;
    let last = offsets.iter().rposition(|offset| *offset >= 0)?;
    let first = lower.checked_add(i32::try_from(first).ok()?)?;
    let last = lower.checked_add(i32::try_from(last).ok()?)?;
    Some((first, last))
}

#[inline(always)]
fn load_i32s<V: Simd>(simd: V, values: &[i32]) -> V::i32s {
    debug_assert!(values.len() <= V::I32_LANES);
    if values.len() == V::I32_LANES {
        let (vectors, remainder) = V::as_simd_i32s(values);
        debug_assert!(remainder.is_empty());
        vectors[0]
    } else {
        simd.partial_load_i32s(values)
    }
}

#[inline(always)]
fn store_i32s<V: Simd>(simd: V, output: &mut [i32], values: V::i32s) {
    debug_assert!(output.len() <= V::I32_LANES);
    if output.len() == V::I32_LANES {
        let (vectors, remainder) = V::as_mut_simd_i32s(output);
        debug_assert!(remainder.is_empty());
        vectors[0] = values;
    } else {
        simd.partial_store_i32s(output, values);
    }
}

struct RecurrenceOperation<'a> {
    mismatch_main: &'a [i32],
    opening_insertion: &'a [i32],
    extension_insertion: &'a [i32],
    opening_deletion: &'a [i32],
    extension_deletion: &'a [i32],
    main_offsets: &'a mut [i32],
    insertion_offsets: &'a mut [i32],
    deletion_offsets: &'a mut [i32],
    lower: i32,
    dimensions: CompactDimensions,
    main_completion: Option<(i32, i32)>,
}

impl WithSimd for RecurrenceOperation<'_> {
    type Output = ();

    #[inline(always)]
    fn with_simd<V: Simd>(mut self, simd: V) {
        let length = self.main_offsets.len();
        debug_assert_eq!(self.mismatch_main.len(), length);
        debug_assert_eq!(self.opening_insertion.len(), length);
        debug_assert_eq!(self.extension_insertion.len(), length);
        debug_assert_eq!(self.opening_deletion.len(), length);
        debug_assert_eq!(self.extension_deletion.len(), length);
        debug_assert_eq!(self.insertion_offsets.len(), length);
        debug_assert_eq!(self.deletion_offsets.len(), length);
        debug_assert!(V::I32_LANES <= MAXIMUM_I32_LANES);

        let (completion_lower, completion_upper) =
            self.main_completion.unwrap_or((i32::MIN, i32::MAX));
        let vector_length = if V::IS_SCALAR {
            length
        } else {
            length / V::I32_LANES * V::I32_LANES
        };
        if vector_length == 0 {
            for index in 0..length {
                self.compute_scalar(index, completion_lower, completion_upper);
            }
            return;
        }

        let mut lane_values = [0_i32; MAXIMUM_I32_LANES];
        for (index, lane) in lane_values[..V::I32_LANES].iter_mut().enumerate() {
            *lane = i32::try_from(index).expect("a SIMD lane index fits i32");
        }
        let lane_indices = load_i32s(simd, &lane_values[..V::I32_LANES]);
        let zero = simd.splat_i32s(0);
        let one = simd.splat_i32s(1);
        let first_limit = simd.splat_i32s(self.dimensions.first);
        let second_limit = simd.splat_i32s(self.dimensions.second);
        let unreachable = simd.splat_i32s(UNREACHABLE);
        let completion_lower = simd.splat_i32s(completion_lower);
        let completion_upper = simd.splat_i32s(completion_upper);

        let mut start = 0_usize;
        while start < vector_length {
            let end = start + V::I32_LANES;
            let diagonal_base = self
                .lower
                .wrapping_add(i32::try_from(start).expect("compact width fits i32"));
            let diagonals = simd.add_i32s(simd.splat_i32s(diagonal_base), lane_indices);
            let insertion = simd.add_i32s(
                simd.max_i32s(
                    load_i32s(simd, &self.opening_insertion[start..end]),
                    load_i32s(simd, &self.extension_insertion[start..end]),
                ),
                one,
            );
            let deletion = simd.max_i32s(
                load_i32s(simd, &self.opening_deletion[start..end]),
                load_i32s(simd, &self.extension_deletion[start..end]),
            );
            let substitution = simd.add_i32s(load_i32s(simd, &self.mismatch_main[start..end]), one);
            let candidate = simd.max_i32s(substitution, simd.max_i32s(deletion, insertion));
            let vertical = simd.sub_i32s(candidate, diagonals);
            let valid = simd.and_m32s(
                simd.and_m32s(
                    simd.greater_than_or_equal_i32s(candidate, zero),
                    simd.less_than_or_equal_i32s(candidate, second_limit),
                ),
                simd.and_m32s(
                    simd.and_m32s(
                        simd.greater_than_or_equal_i32s(vertical, zero),
                        simd.less_than_or_equal_i32s(vertical, first_limit),
                    ),
                    simd.and_m32s(
                        simd.greater_than_or_equal_i32s(diagonals, completion_lower),
                        simd.less_than_or_equal_i32s(diagonals, completion_upper),
                    ),
                ),
            );
            let main = simd.select_i32s(valid, candidate, unreachable);
            store_i32s(simd, &mut self.insertion_offsets[start..end], insertion);
            store_i32s(simd, &mut self.deletion_offsets[start..end], deletion);
            store_i32s(simd, &mut self.main_offsets[start..end], main);
            start = end;
        }
        let (completion_lower, completion_upper) =
            self.main_completion.unwrap_or((i32::MIN, i32::MAX));
        for index in vector_length..length {
            self.compute_scalar(index, completion_lower, completion_upper);
        }
    }
}

impl RecurrenceOperation<'_> {
    #[inline(always)]
    fn compute_scalar(&mut self, index: usize, completion_lower: i32, completion_upper: i32) {
        let diagonal = self
            .lower
            .wrapping_add(i32::try_from(index).expect("compact width fits i32"));
        let insertion = self.opening_insertion[index].max(self.extension_insertion[index]) + 1;
        let deletion = self.opening_deletion[index].max(self.extension_deletion[index]);
        let substitution = self.mismatch_main[index] + 1;
        let candidate = substitution.max(deletion.max(insertion));
        let vertical = candidate - diagonal;
        let valid = candidate >= 0
            && candidate <= self.dimensions.second
            && vertical >= 0
            && vertical <= self.dimensions.first
            && diagonal >= completion_lower
            && diagonal <= completion_upper;
        self.insertion_offsets[index] = insertion;
        self.deletion_offsets[index] = deletion;
        self.main_offsets[index] = if valid { candidate } else { UNREACHABLE };
    }
}

#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]
fn compute_wavefront<const TRACK_MATCHES: bool, const REUSE_PADDING: bool, V: Simd>(
    first: &[u8],
    second: &[u8],
    penalty: Score,
    penalties: WavefrontAffinePenalties,
    dimensions: CompactDimensions,
    upper_bound: Score,
    prune_main: bool,
    padding: usize,
    mismatch: Option<&CompactWavefront>,
    opening: Option<&CompactWavefront>,
    extension: Option<&CompactWavefront>,
    simd: V,
    extension_kernel: WavefrontExtension,
    recycled: &mut Vec<CompactWavefront>,
    null_offsets: &mut Vec<i32>,
) -> Result<CompactStep, WavefrontError> {
    let mut insertion_lower = i32::MAX;
    let mut insertion_upper = i32::MIN;
    merge_shifted_interval(
        &mut insertion_lower,
        &mut insertion_upper,
        opening.and_then(|wavefront| wavefront.main_active),
        1,
        dimensions,
    );
    merge_shifted_interval(
        &mut insertion_lower,
        &mut insertion_upper,
        extension.and_then(|wavefront| wavefront.insertion_active),
        1,
        dimensions,
    );
    let insertion_active = active_interval_from_bounds(insertion_lower, insertion_upper);

    let mut deletion_lower = i32::MAX;
    let mut deletion_upper = i32::MIN;
    merge_shifted_interval(
        &mut deletion_lower,
        &mut deletion_upper,
        opening.and_then(|wavefront| wavefront.main_active),
        -1,
        dimensions,
    );
    merge_shifted_interval(
        &mut deletion_lower,
        &mut deletion_upper,
        extension.and_then(|wavefront| wavefront.deletion_active),
        -1,
        dimensions,
    );
    let deletion_active = active_interval_from_bounds(deletion_lower, deletion_upper);

    let mut lower = i32::MAX;
    let mut upper = i32::MIN;
    merge_shifted_interval(
        &mut lower,
        &mut upper,
        mismatch.and_then(|wavefront| wavefront.main_active),
        0,
        dimensions,
    );
    merge_shifted_interval(&mut lower, &mut upper, insertion_active, 0, dimensions);
    merge_shifted_interval(&mut lower, &mut upper, deletion_active, 0, dimensions);
    constrain_diagonal_bounds(
        &mut lower,
        &mut upper,
        penalty,
        upper_bound,
        penalties.gap_extension(),
        dimensions,
    );
    if lower > upper {
        return Ok(CompactStep::Empty);
    }
    let insertion_active = intersect_interval(insertion_active, lower, upper);
    let deletion_active = intersect_interval(deletion_active, lower, upper);

    let length = checked_width(lower, upper, first.len(), second.len())?;
    if null_offsets.len() < length {
        null_offsets
            .try_reserve(length - null_offsets.len())
            .map_err(|_| WavefrontError::AllocationFailed {
                buffer: "compact affine WFA null source",
                elements: length,
            })?;
        null_offsets.resize(length, UNREACHABLE);
    }
    let null_offsets = &null_offsets[..length];
    let Some(mismatch_main) = state_slice(mismatch, State::Main, lower, length, null_offsets)
    else {
        return Ok(CompactStep::Unsupported);
    };
    let Some(opening_insertion) =
        state_slice(opening, State::Main, lower - 1, length, null_offsets)
    else {
        return Ok(CompactStep::Unsupported);
    };
    let Some(extension_insertion) =
        state_slice(extension, State::Insertion, lower - 1, length, null_offsets)
    else {
        return Ok(CompactStep::Unsupported);
    };
    let Some(opening_deletion) = state_slice(opening, State::Main, lower + 1, length, null_offsets)
    else {
        return Ok(CompactStep::Unsupported);
    };
    let Some(extension_deletion) =
        state_slice(extension, State::Deletion, lower + 1, length, null_offsets)
    else {
        return Ok(CompactStep::Unsupported);
    };

    let mut current = if let Some(mut wavefront) = recycled.pop() {
        wavefront.reset::<REUSE_PADDING>(
            penalty,
            lower,
            length,
            padding,
            first.len(),
            second.len(),
        )?;
        wavefront
    } else {
        CompactWavefront::new(penalty, lower, length, padding, first.len(), second.len())?
    };
    let main_completion = if !prune_main || upper_bound == Score::MAX {
        None
    } else {
        let remaining = upper_bound - penalty;
        let shift = if remaining >= penalties.gap_opening() {
            i32::try_from((remaining - penalties.gap_opening()) / penalties.gap_extension())
                .unwrap_or(i32::MAX)
        } else {
            0
        };
        Some((
            dimensions.target_diagonal.saturating_sub(shift),
            dimensions.target_diagonal.saturating_add(shift),
        ))
    };

    let storage_width = current.storage_width;
    let logical_start = current.padding;
    let logical_end = logical_start + length;
    let (main_storage, gaps) = current.offsets.split_at_mut(storage_width);
    let (insertion_storage, deletion_storage) = gaps.split_at_mut(storage_width);
    let main_offsets = &mut main_storage[logical_start..logical_end];
    let insertion_offsets = &mut insertion_storage[logical_start..logical_end];
    let deletion_offsets = &mut deletion_storage[logical_start..logical_end];
    RecurrenceOperation {
        mismatch_main,
        opening_insertion,
        extension_insertion,
        opening_deletion,
        extension_deletion,
        main_offsets,
        insertion_offsets,
        deletion_offsets,
        lower,
        dimensions,
        main_completion,
    }
    .with_simd(simd);

    let match_extensions;
    current.main_active = if TRACK_MATCHES {
        let result = extension_kernel.extend_with_matches(first, second, lower, main_offsets);
        match_extensions = result.matches();
        result.range()
    } else {
        match_extensions = 0;
        extension_kernel.extend(first, second, lower, main_offsets)
    }
    .map(|range| (range.lower(), range.upper()));
    debug_assert!(
        match (active_interval(insertion_offsets, lower), insertion_active) {
            (None, _) => true,
            (Some(_), None) => false,
            (Some((actual_lower, actual_upper)), Some((predicted_lower, predicted_upper))) => {
                predicted_lower <= actual_lower && actual_upper <= predicted_upper
            }
        }
    );
    debug_assert!(
        match (active_interval(deletion_offsets, lower), deletion_active) {
            (None, _) => true,
            (Some(_), None) => false,
            (Some((actual_lower, actual_upper)), Some((predicted_lower, predicted_upper))) => {
                predicted_lower <= actual_lower && actual_upper <= predicted_upper
            }
        }
    );
    current.insertion_active = insertion_active;
    current.deletion_active = deletion_active;
    if current.main_active.is_none()
        && current.insertion_active.is_none()
        && current.deletion_active.is_none()
    {
        recycled.push(current);
        Ok(CompactStep::Empty)
    } else {
        Ok(CompactStep::Wavefront(current, match_extensions))
    }
}

fn initial_wavefront<const REUSE_PADDING: bool>(
    first: &[u8],
    second: &[u8],
    padding: usize,
    recycled: &mut Vec<CompactWavefront>,
) -> Result<CompactWavefront, WavefrontError> {
    let mut wavefront = if let Some(mut wavefront) = recycled.pop() {
        wavefront.reset::<REUSE_PADDING>(0, 0, 1, padding, first.len(), second.len())?;
        wavefront
    } else {
        CompactWavefront::new(0, 0, 1, padding, first.len(), second.len())?
    };
    let offset = extend_equal(first, second, 0, 0);
    let index = wavefront.padding;
    let storage_width = wavefront.storage_width;
    wavefront.offsets[index] = offset;
    wavefront.offsets[storage_width + index] = UNREACHABLE;
    wavefront.offsets[2 * storage_width + index] = UNREACHABLE;
    wavefront.main_active = Some((0, 0));
    Ok(wavefront)
}

fn initial_trace_wavefront<const REUSE_PADDING: bool>(
    first: &[u8],
    second: &[u8],
    padding: usize,
    recycled_wavefronts: &mut Vec<CompactWavefront>,
) -> Result<CompactTraceWavefront, WavefrontError> {
    Ok(CompactTraceWavefront {
        wavefront: initial_wavefront::<REUSE_PADDING>(first, second, padding, recycled_wavefronts)?,
    })
}

fn recycle_trace_records(
    history: &mut Vec<CompactTraceWavefront>,
    recycled_wavefronts: &mut Vec<CompactWavefront>,
) -> bool {
    let preserve_layout = history.len() >= MINIMUM_ORDERED_TRACE_WAVEFRONTS;
    if preserve_layout {
        for record in history.drain(..).rev() {
            recycled_wavefronts.push(record.wavefront);
        }
    } else {
        for record in history.drain(..) {
            recycled_wavefronts.push(record.wavefront);
        }
    }
    preserve_layout
}

fn find_compact_trace_index(
    history: &[CompactTraceWavefront],
    mut index: usize,
    penalty: Score,
) -> Option<usize> {
    loop {
        let record = history.get(index)?;
        if record.wavefront.penalty == penalty {
            return Some(index);
        }
        if record.wavefront.penalty < penalty || index == 0 {
            return None;
        }
        index -= 1;
    }
}

fn compact_trace_source_index(
    history: &[CompactTraceWavefront],
    history_index: usize,
    penalty: Score,
    decrement: Score,
) -> Option<usize> {
    let source_penalty = penalty.checked_sub(decrement)?;
    (source_penalty >= 0)
        .then(|| find_compact_trace_index(history, history_index, source_penalty))
        .flatten()
}

fn push_compact_trace_run(
    runs: &mut Vec<WavefrontTraceRun>,
    operation: AlignmentOperation,
    count: usize,
    first_len: usize,
    second_len: usize,
) -> Result<(), WavefrontError> {
    if count == 0 {
        return Ok(());
    }
    if let Some(last) = runs.last_mut().filter(|run| run.operation() == operation) {
        let count = last
            .count()
            .checked_add(count)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        *last = WavefrontTraceRun::from_parts(operation, count);
        return Ok(());
    }
    if runs.len() == runs.capacity() {
        runs.try_reserve(1)
            .map_err(|_| WavefrontError::AllocationFailed {
                buffer: "compact affine WFA traceback runs",
                elements: runs.len().saturating_add(1),
            })?;
    }
    runs.push(WavefrontTraceRun::from_parts(operation, count));
    Ok(())
}

fn compact_match_length(
    dimensions: CompactDimensions,
    diagonal: i32,
    source_offset: i32,
    offset: &mut i32,
    penalty: Score,
) -> Result<usize, WavefrontError> {
    let valid = source_offset >= 0
        && source_offset <= *offset
        && source_offset <= dimensions.second
        && *offset <= dimensions.second
        && source_offset - diagonal >= 0
        && source_offset - diagonal <= dimensions.first
        && *offset - diagonal >= 0
        && *offset - diagonal <= dimensions.first;
    if !valid {
        return Err(WavefrontError::TracebackInvariant {
            penalty,
            diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
            offset: isize::try_from(*offset).expect("i32 fits isize"),
        });
    }
    let matches = usize::try_from(*offset - source_offset).map_err(|_| {
        WavefrontError::DimensionOverflow {
            first_len: usize::try_from(dimensions.first).unwrap_or(usize::MAX),
            second_len: usize::try_from(dimensions.second).unwrap_or(usize::MAX),
        }
    })?;
    *offset = source_offset;
    Ok(matches)
}

#[derive(Clone, Copy)]
enum CompactTracebackState {
    Main,
    Insertion,
    Deletion,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn traceback_compact(
    first: &[u8],
    second: &[u8],
    dimensions: CompactDimensions,
    penalties: WavefrontAffinePenalties,
    terminal_penalty: Score,
    history: &[CompactTraceWavefront],
    runs: &mut Vec<WavefrontTraceRun>,
) -> Result<(), WavefrontError> {
    let mut penalty = terminal_penalty;
    let mut diagonal = dimensions.target_diagonal;
    let mut offset = dimensions.second;
    let mut state = CompactTracebackState::Main;
    let mut history_index = history.len().saturating_sub(1);

    loop {
        history_index = find_compact_trace_index(history, history_index, penalty).ok_or(
            WavefrontError::TracebackInvariant {
                penalty,
                diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                offset: isize::try_from(offset).expect("i32 fits isize"),
            },
        )?;
        let record = &history[history_index];
        match state {
            CompactTracebackState::Main => {
                let main = record.wavefront.main(diagonal);
                if main != offset {
                    return Err(WavefrontError::TracebackInvariant {
                        penalty,
                        diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                        offset: isize::try_from(offset).expect("i32 fits isize"),
                    });
                }
                let mismatch_index = compact_trace_source_index(
                    history,
                    history_index,
                    penalty,
                    penalties.mismatch(),
                );
                let substitution = mismatch_index
                    .map_or(UNREACHABLE, |index| history[index].wavefront.main(diagonal))
                    + 1;
                let deletion = record.wavefront.deletion(diagonal);
                let insertion = record.wavefront.insertion(diagonal);
                let (source_offset, main_trace) = if penalty == 0 {
                    (0, CompactMainTrace::Start)
                } else {
                    compact_main_trace(substitution, deletion, insertion)
                };
                let match_length = compact_match_length(
                    dimensions,
                    diagonal,
                    source_offset,
                    &mut offset,
                    penalty,
                )?;
                push_compact_trace_run(
                    runs,
                    AlignmentOperation::Match,
                    match_length,
                    first.len(),
                    second.len(),
                )?;
                match main_trace {
                    CompactMainTrace::Start => {
                        if penalty != 0 || diagonal != 0 || offset != 0 {
                            return Err(WavefrontError::TracebackInvariant {
                                penalty,
                                diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                                offset: isize::try_from(offset).expect("i32 fits isize"),
                            });
                        }
                        break;
                    }
                    CompactMainTrace::Substitution => {
                        let next_history_index =
                            mismatch_index.ok_or(WavefrontError::TracebackInvariant {
                                penalty,
                                diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                                offset: isize::try_from(offset).expect("i32 fits isize"),
                            })?;
                        push_compact_trace_run(
                            runs,
                            AlignmentOperation::Substitution,
                            1,
                            first.len(),
                            second.len(),
                        )?;
                        history_index = next_history_index;
                        offset -= 1;
                        penalty = penalty.checked_sub(penalties.mismatch()).ok_or(
                            WavefrontError::TracebackInvariant {
                                penalty,
                                diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                                offset: isize::try_from(offset).expect("i32 fits isize"),
                            },
                        )?;
                    }
                    CompactMainTrace::Deletion => state = CompactTracebackState::Deletion,
                    CompactMainTrace::Insertion => state = CompactTracebackState::Insertion,
                }
            }
            CompactTracebackState::Deletion => {
                let gap_open_transition = penalties.gap_opening() + penalties.gap_extension();
                let opening_index = compact_trace_source_index(
                    history,
                    history_index,
                    penalty,
                    gap_open_transition,
                );
                let extension_index = compact_trace_source_index(
                    history,
                    history_index,
                    penalty,
                    penalties.gap_extension(),
                );
                let opening = opening_index.map_or(UNREACHABLE, |index| {
                    history[index].wavefront.main(diagonal + 1)
                });
                let extension = extension_index.map_or(UNREACHABLE, |index| {
                    history[index].wavefront.deletion(diagonal + 1)
                });
                let extends = extension > opening;
                let source_offset = opening.max(extension);
                if record.wavefront.deletion(diagonal) != offset
                    || source_offset != offset
                    || source_offset < 0
                {
                    return Err(WavefrontError::TracebackInvariant {
                        penalty,
                        diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                        offset: isize::try_from(offset).expect("i32 fits isize"),
                    });
                }
                push_compact_trace_run(
                    runs,
                    AlignmentOperation::Deletion,
                    1,
                    first.len(),
                    second.len(),
                )?;
                diagonal += 1;
                let decrement = if extends {
                    penalties.gap_extension()
                } else {
                    gap_open_transition
                };
                penalty =
                    penalty
                        .checked_sub(decrement)
                        .ok_or(WavefrontError::TracebackInvariant {
                            penalty,
                            diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                            offset: isize::try_from(offset).expect("i32 fits isize"),
                        })?;
                history_index = if extends {
                    extension_index
                } else {
                    opening_index
                }
                .ok_or(WavefrontError::TracebackInvariant {
                    penalty,
                    diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                    offset: isize::try_from(offset).expect("i32 fits isize"),
                })?;
                state = if extends {
                    CompactTracebackState::Deletion
                } else {
                    CompactTracebackState::Main
                };
            }
            CompactTracebackState::Insertion => {
                let gap_open_transition = penalties.gap_opening() + penalties.gap_extension();
                let opening_index = compact_trace_source_index(
                    history,
                    history_index,
                    penalty,
                    gap_open_transition,
                );
                let extension_index = compact_trace_source_index(
                    history,
                    history_index,
                    penalty,
                    penalties.gap_extension(),
                );
                let opening = opening_index.map_or(UNREACHABLE, |index| {
                    history[index].wavefront.main(diagonal - 1)
                });
                let extension = extension_index.map_or(UNREACHABLE, |index| {
                    history[index].wavefront.insertion(diagonal - 1)
                });
                let extends = extension > opening;
                let source_offset = opening.max(extension);
                let insertion = source_offset + 1;
                if record.wavefront.insertion(diagonal) != offset
                    || insertion != offset
                    || source_offset < 0
                {
                    return Err(WavefrontError::TracebackInvariant {
                        penalty,
                        diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                        offset: isize::try_from(offset).expect("i32 fits isize"),
                    });
                }
                push_compact_trace_run(
                    runs,
                    AlignmentOperation::Insertion,
                    1,
                    first.len(),
                    second.len(),
                )?;
                offset -= 1;
                diagonal -= 1;
                let decrement = if extends {
                    penalties.gap_extension()
                } else {
                    gap_open_transition
                };
                penalty =
                    penalty
                        .checked_sub(decrement)
                        .ok_or(WavefrontError::TracebackInvariant {
                            penalty,
                            diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                            offset: isize::try_from(offset).expect("i32 fits isize"),
                        })?;
                history_index = if extends {
                    extension_index
                } else {
                    opening_index
                }
                .ok_or(WavefrontError::TracebackInvariant {
                    penalty,
                    diagonal: isize::try_from(diagonal).expect("i32 fits isize"),
                    offset: isize::try_from(offset).expect("i32 fits isize"),
                })?;
                state = if extends {
                    CompactTracebackState::Insertion
                } else {
                    CompactTracebackState::Main
                };
            }
        }
    }
    Ok(())
}

fn recycle_slots(slots: &mut [Option<CompactWavefront>], recycled: &mut Vec<CompactWavefront>) {
    for slot in slots {
        if let Some(wavefront) = slot.take() {
            recycled.push(wavefront);
        }
    }
}

/// Attempts the compact score-only kernel, returning `None` for a safe fallback.
#[inline(always)]
pub(crate) fn compact_affine_penalty(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut CompactScoreWorkspace,
    upper_bound: Score,
    prune_main: bool,
) -> Result<Option<Score>, WavefrontError> {
    let Some((architecture, extension_kernel, padding)) =
        compact_configuration(first, second, penalties, upper_bound, workspace)
    else {
        return Ok(None);
    };
    dispatch_compact_alignment(
        architecture,
        first,
        second,
        penalties,
        workspace,
        upper_bound,
        prune_main,
        extension_kernel,
        padding,
    )
}

/// Attempts the compact full-history traceback kernel, returning `None` for a
/// safe fallback to the legacy exact implementation.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compact_affine_traceback(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut CompactScoreWorkspace,
    upper_bound: Score,
    prune_main: bool,
    runs: &mut Vec<WavefrontTraceRun>,
) -> Result<Option<CompactTracebackSummary>, WavefrontError> {
    if compact_traceback_dispatch_bound(first, second, penalties, upper_bound).is_none() {
        return Ok(None);
    }
    let Some((architecture, extension_kernel, padding)) =
        compact_configuration(first, second, penalties, upper_bound, workspace)
    else {
        return Ok(None);
    };
    dispatch_compact_traceback(
        architecture,
        first,
        second,
        penalties,
        workspace,
        upper_bound,
        prune_main,
        extension_kernel,
        padding,
        runs,
    )
}

#[allow(clippy::too_many_arguments)]
#[inline(never)]
fn dispatch_compact_traceback(
    architecture: pulp::Arch,
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut CompactScoreWorkspace,
    upper_bound: Score,
    prune_main: bool,
    extension_kernel: WavefrontExtension,
    padding: usize,
    runs: &mut Vec<WavefrontTraceRun>,
) -> Result<Option<CompactTracebackSummary>, WavefrontError> {
    architecture.dispatch(CompactTracebackOperation {
        first,
        second,
        penalties,
        workspace,
        upper_bound,
        prune_main,
        extension_kernel,
        padding,
        runs,
    })
}

struct CompactTracebackOperation<'a> {
    first: &'a [u8],
    second: &'a [u8],
    penalties: WavefrontAffinePenalties,
    workspace: &'a mut CompactScoreWorkspace,
    upper_bound: Score,
    prune_main: bool,
    extension_kernel: WavefrontExtension,
    padding: usize,
    runs: &'a mut Vec<WavefrontTraceRun>,
}

impl WithSimd for CompactTracebackOperation<'_> {
    type Output = Result<Option<CompactTracebackSummary>, WavefrontError>;

    #[inline(always)]
    fn with_simd<V: Simd>(self, simd: V) -> Self::Output {
        if self.workspace.trace.reuse_padding {
            compact_affine_traceback_core::<true, _>(
                self.first,
                self.second,
                self.penalties,
                self.workspace,
                self.upper_bound,
                self.prune_main,
                self.extension_kernel,
                self.padding,
                self.runs,
                simd,
            )
        } else {
            compact_affine_traceback_core::<false, _>(
                self.first,
                self.second,
                self.penalties,
                self.workspace,
                self.upper_bound,
                self.prune_main,
                self.extension_kernel,
                self.padding,
                self.runs,
                simd,
            )
        }
    }
}

fn compact_trace_history_bytes(
    history: &[CompactTraceWavefront],
    first_len: usize,
    second_len: usize,
) -> Result<usize, WavefrontError> {
    history.iter().try_fold(0_usize, |total, record| {
        let offset_bytes = record
            .wavefront
            .offsets
            .len()
            .checked_mul(size_of::<i32>())
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        total
            .checked_add(offset_bytes)
            .ok_or(WavefrontError::DimensionOverflow {
                first_len,
                second_len,
            })
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn compact_affine_traceback_core<const REUSE_PADDING: bool, V: Simd>(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut CompactScoreWorkspace,
    upper_bound: Score,
    prune_main: bool,
    extension_kernel: WavefrontExtension,
    padding: usize,
    runs: &mut Vec<WavefrontTraceRun>,
    simd: V,
) -> Result<Option<CompactTracebackSummary>, WavefrontError> {
    let Some(dimensions) = CompactDimensions::new(first.len(), second.len(), padding) else {
        return Ok(None);
    };
    let slot_count = padding;
    let CompactScoreWorkspace {
        trace,
        null_offsets,
        slots: _,
        recycled: _,
        architecture: _,
        extension: _,
    } = workspace;
    trace.recycle_history();
    let CompactTraceWorkspace {
        history,
        recycled_wavefronts,
        reuse_padding,
    } = trace;
    runs.clear();

    if history.capacity() == 0 {
        history
            .try_reserve_exact(1)
            .map_err(|_| WavefrontError::AllocationFailed {
                buffer: "compact affine WFA history",
                elements: 1,
            })?;
    }
    let initial =
        initial_trace_wavefront::<REUSE_PADDING>(first, second, padding, recycled_wavefronts)?;
    let initial_matches = usize::try_from(initial.wavefront.main(0)).map_err(|_| {
        WavefrontError::TracebackInvariant {
            penalty: 0,
            diagonal: 0,
            offset: isize::try_from(initial.wavefront.main(0)).expect("i32 fits isize"),
        }
    })?;
    let initial_reaches_target = initial.wavefront.reaches_target(dimensions);
    history.push(initial);

    let mut indices = [usize::MAX; 64];
    indices[0] = 0;
    let mut wavefronts = 1_usize;
    let mut computed_offsets = 1_usize;
    let mut source_wavefront_lookups = 0_usize;
    let mut source_offset_lookups = 0_usize;
    let mut match_extensions = initial_matches;
    let mut maximum_wavefront_width = 1_usize;

    if initial_reaches_target {
        traceback_compact(first, second, dimensions, penalties, 0, history, runs)?;
        let modeled_peak_bytes = compact_trace_history_bytes(history, first.len(), second.len())?;
        let summary = CompactTracebackSummary {
            penalty: 0,
            wavefronts,
            computed_offsets,
            source_wavefront_lookups,
            source_offset_lookups,
            match_extensions,
            maximum_wavefront_width,
            modeled_peak_bytes,
        };
        *reuse_padding = recycle_trace_records(history, recycled_wavefronts);
        return Ok(Some(summary));
    }

    let gap_open_transition = penalties.gap_opening() + penalties.gap_extension();
    let mut scheduler = CompactPenaltyScheduler::new();
    scheduler.schedule(0, penalties.mismatch());
    scheduler.schedule(0, gap_open_transition);
    while let Some(penalty) = scheduler.pop() {
        if upper_bound != Score::MAX && penalty > upper_bound {
            *reuse_padding = recycle_trace_records(history, recycled_wavefronts);
            runs.clear();
            return Ok(None);
        }
        source_wavefront_lookups =
            source_wavefront_lookups
                .checked_add(3)
                .ok_or(WavefrontError::DimensionOverflow {
                    first_len: first.len(),
                    second_len: second.len(),
                })?;

        let mismatch =
            trace_source_wavefront(history, &indices, slot_count, penalty, penalties.mismatch());
        let opening =
            trace_source_wavefront(history, &indices, slot_count, penalty, gap_open_transition);
        let extension = trace_source_wavefront(
            history,
            &indices,
            slot_count,
            penalty,
            penalties.gap_extension(),
        );
        let step = compute_wavefront::<true, REUSE_PADDING, _>(
            first,
            second,
            penalty,
            penalties,
            dimensions,
            upper_bound,
            prune_main,
            padding,
            mismatch.map(|record| &record.wavefront),
            opening.map(|record| &record.wavefront),
            extension.map(|record| &record.wavefront),
            simd,
            extension_kernel,
            recycled_wavefronts,
            null_offsets,
        )?;
        match step {
            CompactStep::Unsupported => {
                *reuse_padding = recycle_trace_records(history, recycled_wavefronts);
                runs.clear();
                return Ok(None);
            }
            CompactStep::Empty => {}
            CompactStep::Wavefront(current, current_matches) => {
                let width = current.width;
                let reaches_target = current.reaches_target(dimensions);
                let has_main = current.main_active.is_some();
                let has_gap =
                    current.insertion_active.is_some() || current.deletion_active.is_some();

                wavefronts =
                    wavefronts
                        .checked_add(1)
                        .ok_or(WavefrontError::DimensionOverflow {
                            first_len: first.len(),
                            second_len: second.len(),
                        })?;
                computed_offsets = computed_offsets.checked_add(width).ok_or(
                    WavefrontError::DimensionOverflow {
                        first_len: first.len(),
                        second_len: second.len(),
                    },
                )?;
                source_offset_lookups =
                    source_offset_lookups
                        .checked_add(width.checked_mul(5).ok_or(
                            WavefrontError::DimensionOverflow {
                                first_len: first.len(),
                                second_len: second.len(),
                            },
                        )?)
                        .ok_or(WavefrontError::DimensionOverflow {
                            first_len: first.len(),
                            second_len: second.len(),
                        })?;
                match_extensions = match_extensions.checked_add(current_matches).ok_or(
                    WavefrontError::DimensionOverflow {
                        first_len: first.len(),
                        second_len: second.len(),
                    },
                )?;
                maximum_wavefront_width = maximum_wavefront_width.max(width);

                if history.len() == history.capacity() {
                    history
                        .try_reserve(1)
                        .map_err(|_| WavefrontError::AllocationFailed {
                            buffer: "compact affine WFA history",
                            elements: history.len().saturating_add(1),
                        })?;
                }
                history.push(CompactTraceWavefront { wavefront: current });
                let slot_count_score =
                    Score::try_from(slot_count).expect("compact slot count fits score");
                let slot = usize::try_from(penalty % slot_count_score)
                    .expect("a compact nonnegative penalty slot fits usize");
                indices[slot] = history.len() - 1;

                if reaches_target {
                    traceback_compact(
                        first, second, dimensions, penalties, penalty, history, runs,
                    )?;
                    let modeled_peak_bytes =
                        compact_trace_history_bytes(history, first.len(), second.len())?;
                    let summary = CompactTracebackSummary {
                        penalty,
                        wavefronts,
                        computed_offsets,
                        source_wavefront_lookups,
                        source_offset_lookups,
                        match_extensions,
                        maximum_wavefront_width,
                        modeled_peak_bytes,
                    };
                    *reuse_padding = recycle_trace_records(history, recycled_wavefronts);
                    return Ok(Some(summary));
                }
                if has_main {
                    scheduler.schedule(penalty, penalties.mismatch());
                    scheduler.schedule(penalty, gap_open_transition);
                }
                if has_gap {
                    scheduler.schedule(penalty, penalties.gap_extension());
                }
            }
        }
    }

    *reuse_padding = recycle_trace_records(history, recycled_wavefronts);
    runs.clear();
    let _ = scheduler.overflowed;
    Ok(None)
}

// Keep the architecture-specialized compact loop out of the legacy score-path
// caller. Inlining this dispatch can perturb the instruction layout of inputs
// that reject the compact kernel before executing any of its work.
#[allow(clippy::too_many_arguments)]
#[inline(never)]
fn dispatch_compact_alignment(
    architecture: pulp::Arch,
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut CompactScoreWorkspace,
    upper_bound: Score,
    prune_main: bool,
    extension_kernel: WavefrontExtension,
    padding: usize,
) -> Result<Option<Score>, WavefrontError> {
    architecture.dispatch(CompactAlignmentOperation {
        first,
        second,
        penalties,
        workspace,
        upper_bound,
        prune_main,
        extension_kernel,
        padding,
    })
}

struct CompactAlignmentOperation<'a> {
    first: &'a [u8],
    second: &'a [u8],
    penalties: WavefrontAffinePenalties,
    workspace: &'a mut CompactScoreWorkspace,
    upper_bound: Score,
    prune_main: bool,
    extension_kernel: WavefrontExtension,
    padding: usize,
}

impl WithSimd for CompactAlignmentOperation<'_> {
    type Output = Result<Option<Score>, WavefrontError>;

    #[inline(always)]
    fn with_simd<V: Simd>(self, simd: V) -> Self::Output {
        compact_affine_penalty_core(
            self.first,
            self.second,
            self.penalties,
            self.workspace,
            self.upper_bound,
            self.prune_main,
            self.extension_kernel,
            self.padding,
            simd,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn compact_affine_penalty_core<V: Simd>(
    first: &[u8],
    second: &[u8],
    penalties: WavefrontAffinePenalties,
    workspace: &mut CompactScoreWorkspace,
    upper_bound: Score,
    prune_main: bool,
    extension_kernel: WavefrontExtension,
    padding: usize,
    simd: V,
) -> Result<Option<Score>, WavefrontError> {
    let Some(dimensions) = CompactDimensions::new(first.len(), second.len(), padding) else {
        return Ok(None);
    };
    let slot_count = padding;

    let CompactScoreWorkspace {
        slots,
        recycled,
        null_offsets,
        trace: _,
        architecture: _,
        extension: _,
    } = workspace;
    recycle_slots(slots, recycled);
    slots.clear();
    slots.resize_with(slot_count, || None);

    let initial = initial_wavefront::<true>(first, second, padding, recycled)?;
    if initial.reaches_target(dimensions) {
        recycled.push(initial);
        return Ok(Some(0));
    }
    slots[0] = Some(initial);

    let gap_open_transition = penalties.gap_opening() + penalties.gap_extension();
    let mut scheduler = CompactPenaltyScheduler::new();
    scheduler.schedule(0, penalties.mismatch());
    scheduler.schedule(0, gap_open_transition);
    while let Some(penalty) = scheduler.pop() {
        if upper_bound != Score::MAX && penalty > upper_bound {
            recycle_slots(slots, recycled);
            return Ok(None);
        }
        let slot_count_score = Score::try_from(slot_count).expect("compact slot count fits score");
        let slot = usize::try_from(penalty % slot_count_score)
            .expect("a compact nonnegative penalty slot fits usize");
        if let Some(expired) = slots[slot].take() {
            recycled.push(expired);
        }

        let mismatch = source_wavefront(slots, penalty, penalties.mismatch());
        let opening = source_wavefront(slots, penalty, gap_open_transition);
        let extension = source_wavefront(slots, penalty, penalties.gap_extension());
        let step = compute_wavefront::<false, true, _>(
            first,
            second,
            penalty,
            penalties,
            dimensions,
            upper_bound,
            prune_main,
            padding,
            mismatch,
            opening,
            extension,
            simd,
            extension_kernel,
            recycled,
            null_offsets,
        )?;
        match step {
            CompactStep::Unsupported => {
                recycle_slots(slots, recycled);
                return Ok(None);
            }
            CompactStep::Empty => {}
            CompactStep::Wavefront(current, _) => {
                let reaches_target = current.reaches_target(dimensions);
                let has_main = current.main_active.is_some();
                let has_gap =
                    current.insertion_active.is_some() || current.deletion_active.is_some();
                slots[slot] = Some(current);
                if reaches_target {
                    recycle_slots(slots, recycled);
                    return Ok(Some(penalty));
                }
                if has_main {
                    scheduler.schedule(penalty, penalties.mismatch());
                    scheduler.schedule(penalty, gap_open_transition);
                }
                if has_gap {
                    scheduler.schedule(penalty, penalties.gap_extension());
                }
            }
        }
    }
    recycle_slots(slots, recycled);
    let _ = scheduler.overflowed;
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_penalties() -> WavefrontAffinePenalties {
        WavefrontAffinePenalties::new(2, 8, 2).expect("test penalties are valid")
    }

    #[test]
    fn bounded_dispatch_gate_starts_at_twenty_nine_estimated_levels() {
        let first = [b'A'; 32];
        let second = first;
        let penalties = test_penalties();

        assert_eq!(compact_dispatch_bound(&first, &second, penalties, 56), None);
        assert_eq!(
            compact_dispatch_bound(&first, &second, penalties, 57),
            Some(57)
        );
    }

    #[test]
    fn bounded_traceback_dispatch_gate_starts_at_twenty_nine_estimated_levels() {
        let first = [b'A'; 64];
        let second = first;
        let penalties = test_penalties();

        assert_eq!(
            compact_traceback_dispatch_bound(&first, &second, penalties, 56),
            None
        );
        assert_eq!(
            compact_traceback_dispatch_bound(&first, &second, penalties, 57),
            Some(57)
        );
    }

    #[test]
    fn rejected_work_does_not_probe_the_runtime_architecture() {
        let first = [b'A'; 32];
        let mut second = first;
        second[0] = b'C';
        let mut workspace = CompactScoreWorkspace::new();

        assert_eq!(
            compact_affine_penalty(&first, &second, test_penalties(), &mut workspace, 56, true,)
                .expect("rejected compact work is a safe fallback"),
            None
        );
        assert!(workspace.architecture.is_none());
    }

    #[test]
    fn trace_history_recycles_in_next_alignment_order() {
        let mut history = (0..MINIMUM_ORDERED_TRACE_WAVEFRONTS)
            .map(|penalty| CompactTraceWavefront {
                wavefront: CompactWavefront::new(
                    Score::try_from(penalty).expect("test penalty fits score"),
                    0,
                    1,
                    1,
                    1,
                    1,
                )
                .expect("test wavefront is valid"),
            })
            .collect::<Vec<_>>();
        let mut recycled = Vec::new();

        assert!(recycle_trace_records(&mut history, &mut recycled));

        assert!(history.is_empty());
        assert_eq!(recycled.pop().map(|wavefront| wavefront.penalty), Some(0));
        assert_eq!(recycled.pop().map(|wavefront| wavefront.penalty), Some(1));
    }

    #[test]
    fn short_trace_history_keeps_the_original_recycling_order() {
        let mut history = [0, 2, 4]
            .into_iter()
            .map(|penalty| CompactTraceWavefront {
                wavefront: CompactWavefront::new(penalty, 0, 1, 1, 1, 1)
                    .expect("test wavefront is valid"),
            })
            .collect::<Vec<_>>();
        let mut recycled = Vec::new();

        assert!(!recycle_trace_records(&mut history, &mut recycled));

        assert!(history.is_empty());
        assert_eq!(recycled.pop().map(|wavefront| wavefront.penalty), Some(4));
        assert_eq!(recycled.pop().map(|wavefront| wavefront.penalty), Some(2));
        assert_eq!(recycled.pop().map(|wavefront| wavefront.penalty), Some(0));
    }

    #[test]
    fn profitable_increment_boundary_is_eighteen() {
        let supported_transition = WavefrontAffinePenalties::new(4, 16, 2).unwrap();
        let supported_mismatch = WavefrontAffinePenalties::new(18, 6, 2).unwrap();
        let rejected = WavefrontAffinePenalties::new(4, 17, 2).unwrap();
        assert_eq!(compact_padding(supported_transition), Some(19));
        assert_eq!(compact_padding(supported_mismatch), Some(19));
        assert_eq!(compact_padding(rejected), None);

        let first = [b'A'; 256];
        let mut second = first;
        second[0] = b'C';
        let mut workspace = CompactScoreWorkspace::new();
        assert_eq!(
            compact_affine_penalty(&first, &second, rejected, &mut workspace, 100, false,).unwrap(),
            None
        );
        assert!(workspace.architecture.is_none());
    }

    #[test]
    fn unbounded_gate_requires_equal_non_tiny_inputs_and_enough_edits() {
        let first = [b'A'; MINIMUM_UNBOUNDED_LENGTH];
        let mut second = first;
        let penalties = test_penalties();

        for symbol in &mut second[..19] {
            *symbol = b'C';
        }
        assert_eq!(
            compact_dispatch_bound(&first, &second, penalties, Score::MAX),
            None
        );
        second[19] = b'C';
        assert_eq!(
            compact_dispatch_bound(&first, &second, penalties, Score::MAX),
            Some(40)
        );
        assert_eq!(
            compact_dispatch_bound(&first[..127], &second[..127], penalties, Score::MAX),
            None
        );
        assert_eq!(
            compact_dispatch_bound(&first, &second[..127], penalties, Score::MAX),
            None
        );
    }

    #[test]
    fn vector_kernel_preserves_a_known_sparse_substitution_penalty() {
        let architecture = pulp::Arch::new();
        if !architecture.dispatch(HasVectorLanes) {
            return;
        }
        let first = [b'A'; 256];
        let mut second = first;
        for symbol in second.iter_mut().step_by(7).take(32) {
            *symbol = b'C';
        }
        let mut workspace = CompactScoreWorkspace::new();

        assert_eq!(
            compact_affine_penalty(&first, &second, test_penalties(), &mut workspace, 64, true,)
                .expect("compact kernel succeeds"),
            Some(64)
        );
    }
}
