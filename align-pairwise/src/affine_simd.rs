//! Safe runtime-dispatched SIMD backend for the exact affine recurrence.
//!
//! Hot helpers deliberately cross generic Pulp dispatch boundaries. Causal codegen
//! inspection and benchmarks require them to inline so the vector loop has no call
//! instructions.
#![allow(clippy::inline_always)]

use core::mem::size_of;

use align_core::{AffineGap, AlignmentOperation, Alphabet, Score, SubstitutionScorer, Symbol};
use pulp::{Simd, WithSimd};

use crate::AlignmentError;
use crate::gotoh::{
    AffineWorkspace, BestTrace, GapTrace, TraceCell, every_path_score_fits_for,
    reserve_workspace_buffer,
};

/// Versioned exact affine SIMD implementation profile.
pub const AFFINE_SIMD_PREFIX_SCAN_V1: &str = "AFFINE-SIMD-PREFIX-SCAN-v1";

const MAX_I32_LANES: usize = 16;
const AUTO_MINIMUM_CELLS: usize = 4_096;
const AUTO_SCORE_MINIMUM_CELLS: usize = 4_096;
const AVX2_I32_LANES: usize = 8;
const AVX2_I16_LANES: usize = 16;

fn padded_length(length: usize, lanes: usize) -> Option<usize> {
    debug_assert!(lanes.is_power_of_two());
    length
        .checked_add(lanes - 1)
        .map(|rounded| rounded & !(lanes - 1))
}

fn auto_shape_is_eligible(first_len: usize, cells: usize, uses_i16: bool) -> bool {
    cells >= AUTO_MINIMUM_CELLS
        && if uses_i16 {
            first_len >= 7
        } else {
            (6..=AVX2_I32_LANES).contains(&first_len) || first_len >= 11
        }
}

fn auto_score_shape_is_eligible(first_len: usize, cells: usize, uses_i16: bool) -> bool {
    cells >= AUTO_SCORE_MINIMUM_CELLS
        && if uses_i16 {
            first_len >= 10
        } else {
            first_len == AVX2_I32_LANES || first_len >= 11
        }
}

/// Policy controlling the implementation backend used with an affine workspace.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[non_exhaustive]
pub enum AffineKernelPolicy {
    /// Selects SIMD only when the ISA, score width, and calibrated size gate allow it.
    #[default]
    Auto,
    /// Always executes the established scalar implementation.
    Scalar,
    /// Requires the SIMD prefix-scan backend or returns a typed error.
    Simd,
}

/// Implementation backend used by the most recent successful workspace call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum AffineBackend {
    /// Scalar bounded or checked dynamic programming.
    Scalar,
    /// Runtime-dispatched signed-16-bit SIMD prefix scan.
    SimdI16 {
        /// Number of signed 16-bit lanes selected at runtime.
        lanes: usize,
    },
    /// Runtime-dispatched signed-32-bit SIMD prefix scan.
    SimdI32 {
        /// Number of signed 32-bit lanes selected at runtime.
        lanes: usize,
    },
}

/// Exact active-buffer memory model for the AVX2 affine SIMD backends.
///
/// Allocation metadata, allocator rounding, input sequences, and retained excess
/// workspace capacity are excluded. The operation-vector estimate is included
/// because a complete alignment with traceback is the backend's public result.
/// Score and profile storage is padded to the signed-32-bit backing width; the
/// signed-16-bit backend safely reinterprets the same storage as twice as many lanes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AffineSimdMemoryEstimate {
    cells: usize,
    final_score_row_bytes: usize,
    score_column_bytes: usize,
    substitution_profile_bytes: usize,
    second_index_bytes: usize,
    traceback_matrix_bytes: usize,
    maximum_operation_bytes: usize,
    peak_buffer_bytes: usize,
}

/// Exact active-buffer memory model for affine SIMD score-only execution.
///
/// Allocation metadata, allocator rounding, input sequences, traceback storage,
/// operation vectors, and retained excess workspace capacity are excluded. Score
/// and profile storage is padded to the signed-32-bit backing width; the signed-16-
/// bit backend safely reinterprets the same storage as twice as many lanes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(clippy::struct_field_names)]
pub struct AffineSimdScoreMemoryEstimate {
    final_score_row_bytes: usize,
    score_column_bytes: usize,
    substitution_profile_bytes: usize,
    second_index_bytes: usize,
    peak_buffer_bytes: usize,
}

impl AffineSimdScoreMemoryEstimate {
    /// Returns bytes occupied by the signed-64-bit terminal-score row.
    #[must_use]
    pub const fn final_score_row_bytes(self) -> usize {
        self.final_score_row_bytes
    }

    /// Returns bytes occupied by both padded signed-32-bit backing score columns.
    #[must_use]
    pub const fn score_column_bytes(self) -> usize {
        self.score_column_bytes
    }

    /// Returns bytes occupied by the padded alphabet-indexed substitution profile.
    #[must_use]
    pub const fn substitution_profile_bytes(self) -> usize {
        self.substitution_profile_bytes
    }

    /// Returns bytes occupied by cached second-sequence alphabet indices.
    #[must_use]
    pub const fn second_index_bytes(self) -> usize {
        self.second_index_bytes
    }

    /// Returns the sum of all simultaneously active score-only SIMD buffers.
    #[must_use]
    pub const fn peak_buffer_bytes(self) -> usize {
        self.peak_buffer_bytes
    }
}

impl AffineSimdMemoryEstimate {
    /// Returns the number of dynamic-programming cells, including boundaries.
    #[must_use]
    pub const fn cells(self) -> usize {
        self.cells
    }

    /// Returns bytes occupied by the signed-64-bit terminal-score row.
    #[must_use]
    pub const fn final_score_row_bytes(self) -> usize {
        self.final_score_row_bytes
    }

    /// Returns bytes occupied by both padded signed-32-bit backing score columns.
    #[must_use]
    pub const fn score_column_bytes(self) -> usize {
        self.score_column_bytes
    }

    /// Returns bytes occupied by the padded alphabet-indexed substitution profile.
    #[must_use]
    pub const fn substitution_profile_bytes(self) -> usize {
        self.substitution_profile_bytes
    }

    /// Returns bytes occupied by cached second-sequence alphabet indices.
    #[must_use]
    pub const fn second_index_bytes(self) -> usize {
        self.second_index_bytes
    }

    /// Returns bytes occupied by column-major packed traceback cells.
    #[must_use]
    pub const fn traceback_matrix_bytes(self) -> usize {
        self.traceback_matrix_bytes
    }

    /// Returns bytes reserved for the maximum possible operation count.
    #[must_use]
    pub const fn maximum_operation_bytes(self) -> usize {
        self.maximum_operation_bytes
    }

    /// Returns the sum of all simultaneously live algorithm-owned buffers.
    #[must_use]
    pub const fn peak_buffer_bytes(self) -> usize {
        self.peak_buffer_bytes
    }
}

impl AffineWorkspace {
    /// Computes exact active element-buffer requirements for the SIMD backends.
    ///
    /// The alphabet type determines the substitution-profile width. Score and
    /// profile buffers include AVX2 backing-width padding. A reused workspace may
    /// retain larger SIMD or scalar buffers from prior calls; those inactive
    /// retained capacities are intentionally excluded.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if a dimension or byte count
    /// cannot be represented by `usize`.
    pub fn simd_memory_requirement<A: Alphabet>(
        first_len: usize,
        second_len: usize,
    ) -> Result<AffineSimdMemoryEstimate, AlignmentError> {
        let overflow = || AlignmentError::DimensionOverflow {
            first_len,
            second_len,
        };
        let padded_first_len = padded_length(first_len, AVX2_I32_LANES).ok_or_else(overflow)?;
        let rows = padded_first_len.checked_add(1).ok_or_else(overflow)?;
        let columns = second_len.checked_add(1).ok_or_else(overflow)?;
        let trace_rows = first_len.checked_add(1).ok_or_else(overflow)?;
        let cells = trace_rows.checked_mul(columns).ok_or_else(overflow)?;
        let profile_elements = padded_first_len
            .checked_mul(A::SYMBOLS.len())
            .ok_or_else(overflow)?;
        let maximum_operations = first_len.checked_add(second_len).ok_or_else(overflow)?;

        let final_score_row_bytes = columns
            .checked_mul(size_of::<Score>())
            .ok_or_else(overflow)?;
        let score_column_bytes = rows
            .checked_mul(size_of::<i32>())
            .and_then(|bytes| bytes.checked_mul(2))
            .ok_or_else(overflow)?;
        let substitution_profile_bytes = profile_elements
            .checked_mul(size_of::<i32>())
            .ok_or_else(overflow)?;
        let second_index_bytes = second_len
            .checked_mul(size_of::<usize>())
            .ok_or_else(overflow)?;
        let traceback_matrix_bytes = cells.checked_mul(size_of::<u8>()).ok_or_else(overflow)?;
        let maximum_operation_bytes = maximum_operations
            .checked_mul(size_of::<AlignmentOperation>())
            .ok_or_else(overflow)?;
        let peak_buffer_bytes = final_score_row_bytes
            .checked_add(score_column_bytes)
            .and_then(|bytes| bytes.checked_add(substitution_profile_bytes))
            .and_then(|bytes| bytes.checked_add(second_index_bytes))
            .and_then(|bytes| bytes.checked_add(traceback_matrix_bytes))
            .and_then(|bytes| bytes.checked_add(maximum_operation_bytes))
            .ok_or_else(overflow)?;

        Ok(AffineSimdMemoryEstimate {
            cells,
            final_score_row_bytes,
            score_column_bytes,
            substitution_profile_bytes,
            second_index_bytes,
            traceback_matrix_bytes,
            maximum_operation_bytes,
            peak_buffer_bytes,
        })
    }

    /// Computes exact active element-buffer requirements for SIMD score-only calls.
    ///
    /// No traceback or operation-vector bytes are included. A reused workspace may
    /// retain inactive capacity from a prior complete alignment; retained excess is
    /// intentionally excluded.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if a dimension or byte count
    /// cannot be represented by `usize`.
    pub fn simd_score_memory_requirement<A: Alphabet>(
        first_len: usize,
        second_len: usize,
    ) -> Result<AffineSimdScoreMemoryEstimate, AlignmentError> {
        let overflow = || AlignmentError::DimensionOverflow {
            first_len,
            second_len,
        };
        first_len.checked_add(second_len).ok_or_else(overflow)?;
        let padded_first_len = padded_length(first_len, AVX2_I32_LANES).ok_or_else(overflow)?;
        let rows = padded_first_len.checked_add(1).ok_or_else(overflow)?;
        let columns = second_len.checked_add(1).ok_or_else(overflow)?;
        let profile_elements = padded_first_len
            .checked_mul(A::SYMBOLS.len())
            .ok_or_else(overflow)?;
        let final_score_row_bytes = columns
            .checked_mul(size_of::<Score>())
            .ok_or_else(overflow)?;
        let score_column_bytes = rows
            .checked_mul(size_of::<i32>())
            .and_then(|bytes| bytes.checked_mul(2))
            .ok_or_else(overflow)?;
        let substitution_profile_bytes = profile_elements
            .checked_mul(size_of::<i32>())
            .ok_or_else(overflow)?;
        let second_index_bytes = second_len
            .checked_mul(size_of::<usize>())
            .ok_or_else(overflow)?;
        let peak_buffer_bytes = final_score_row_bytes
            .checked_add(score_column_bytes)
            .and_then(|bytes| bytes.checked_add(substitution_profile_bytes))
            .and_then(|bytes| bytes.checked_add(second_index_bytes))
            .ok_or_else(overflow)?;
        Ok(AffineSimdScoreMemoryEstimate {
            final_score_row_bytes,
            score_column_bytes,
            substitution_profile_bytes,
            second_index_bytes,
            peak_buffer_bytes,
        })
    }

    /// Fallibly reserves SIMD scratch capacity for an alphabet and maximum shape.
    ///
    /// This method reserves only the strict SIMD path. An automatic-policy call
    /// that is numerically ineligible for SIMD may still allocate scalar fallback
    /// buffers. Existing allocations and contents are retained.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if dimensions cannot be
    /// represented, or [`AlignmentError::AllocationFailed`] if a reservation fails.
    pub fn reserve_simd<A: Alphabet>(
        &mut self,
        first_len: usize,
        second_len: usize,
    ) -> Result<(), AlignmentError> {
        let estimate = Self::simd_memory_requirement::<A>(first_len, second_len)?;
        let padded_first_len =
            padded_length(first_len, AVX2_I32_LANES).ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let rows = padded_first_len + 1;
        let columns = second_len + 1;
        let profile_elements = padded_first_len * A::SYMBOLS.len();

        reserve_workspace_buffer(
            &mut self.best_scores,
            columns,
            "affine workspace SIMD final score row",
        )?;
        reserve_workspace_buffer(
            &mut self.simd_best_scores,
            rows,
            "affine workspace SIMD best-score column",
        )?;
        reserve_workspace_buffer(
            &mut self.simd_insertion_scores,
            rows,
            "affine workspace SIMD insertion-score column",
        )?;
        reserve_workspace_buffer(
            &mut self.simd_profile,
            profile_elements,
            "affine workspace SIMD substitution profile",
        )?;
        reserve_workspace_buffer(
            &mut self.simd_second_indices,
            second_len,
            "affine workspace SIMD second-sequence indices",
        )?;
        reserve_workspace_buffer(
            &mut self.simd_traceback,
            estimate.cells,
            "affine workspace SIMD traceback matrix",
        )
    }

    /// Fallibly reserves only SIMD score-only scratch for an alphabet and shape.
    ///
    /// This does not reserve traceback or operation storage. Existing allocations
    /// and contents are retained.
    ///
    /// # Errors
    ///
    /// Returns [`AlignmentError::DimensionOverflow`] if dimensions cannot be
    /// represented, or [`AlignmentError::AllocationFailed`] if reservation fails.
    pub fn reserve_simd_score<A: Alphabet>(
        &mut self,
        first_len: usize,
        second_len: usize,
    ) -> Result<(), AlignmentError> {
        Self::simd_score_memory_requirement::<A>(first_len, second_len)?;
        let padded_first_len =
            padded_length(first_len, AVX2_I32_LANES).ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let rows = padded_first_len + 1;
        let columns = second_len + 1;
        let profile_elements = padded_first_len * A::SYMBOLS.len();
        reserve_workspace_buffer(
            &mut self.best_scores,
            columns,
            "affine score-only SIMD final score row",
        )?;
        reserve_workspace_buffer(
            &mut self.simd_best_scores,
            rows,
            "affine score-only SIMD best-score column",
        )?;
        reserve_workspace_buffer(
            &mut self.simd_insertion_scores,
            rows,
            "affine score-only SIMD insertion-score column",
        )?;
        reserve_workspace_buffer(
            &mut self.simd_profile,
            profile_elements,
            "affine score-only SIMD substitution profile",
        )?;
        reserve_workspace_buffer(
            &mut self.simd_second_indices,
            second_len,
            "affine score-only SIMD second-sequence indices",
        )
    }

    fn prepare_simd_score(
        &mut self,
        first_len: usize,
        backing_first_len: usize,
        second_len: usize,
        alphabet_width: usize,
    ) -> Result<(), AlignmentError> {
        debug_assert!(backing_first_len >= first_len);
        let rows = backing_first_len
            .checked_add(1)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let columns = second_len
            .checked_add(1)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len,
                second_len,
            })?;
        let profile_elements = backing_first_len.checked_mul(alphabet_width).ok_or(
            AlignmentError::DimensionOverflow {
                first_len,
                second_len: alphabet_width,
            },
        )?;
        reserve_workspace_buffer(
            &mut self.best_scores,
            columns,
            "affine score-only SIMD final score row",
        )?;
        self.best_scores.resize(columns, 0);
        reserve_workspace_buffer(
            &mut self.simd_best_scores,
            rows,
            "affine score-only SIMD best-score column",
        )?;
        self.simd_best_scores.resize(rows, 0);
        reserve_workspace_buffer(
            &mut self.simd_insertion_scores,
            rows,
            "affine score-only SIMD insertion-score column",
        )?;
        self.simd_insertion_scores.resize(rows, 0);
        reserve_workspace_buffer(
            &mut self.simd_profile,
            profile_elements,
            "affine score-only SIMD substitution profile",
        )?;
        self.simd_profile.resize(profile_elements, 0);
        reserve_workspace_buffer(
            &mut self.simd_second_indices,
            second_len,
            "affine score-only SIMD second-sequence indices",
        )?;
        self.simd_second_indices.resize(second_len, 0);
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SimdBoundaryMode {
    Global,
    EndGap {
        first_start_is_free: bool,
        first_end_is_free: bool,
        second_start_is_free: bool,
        second_end_is_free: bool,
    },
    Local,
}

impl SimdBoundaryMode {
    const fn first_start_is_free(self) -> bool {
        match self {
            Self::Global => false,
            Self::EndGap {
                first_start_is_free,
                ..
            } => first_start_is_free,
            Self::Local => true,
        }
    }

    const fn second_start_is_free(self) -> bool {
        match self {
            Self::Global => false,
            Self::EndGap {
                second_start_is_free,
                ..
            } => second_start_is_free,
            Self::Local => true,
        }
    }

    const fn is_local(self) -> bool {
        matches!(self, Self::Local)
    }

    const fn is_global(self) -> bool {
        match self {
            Self::Global => true,
            Self::EndGap {
                first_start_is_free,
                first_end_is_free,
                second_start_is_free,
                second_end_is_free,
            } => {
                !first_start_is_free
                    && !first_end_is_free
                    && !second_start_is_free
                    && !second_end_is_free
            }
            Self::Local => false,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct SimdEndpoint {
    pub(crate) score: Score,
    pub(crate) row: usize,
    pub(crate) column: usize,
}

impl SimdEndpoint {
    const ORIGIN: Self = Self {
        score: 0,
        row: 0,
        column: 0,
    };

    fn observe_local(&mut self, score: Score, row: usize, column: usize) {
        if score > self.score
            || (score == self.score
                && score > 0
                && (row < self.row || (row == self.row && column < self.column)))
        {
            *self = Self { score, row, column };
        }
    }

    fn observe_end_gap(&mut self, score: Score, row: usize, column: usize) {
        let candidate_consumed = row + column;
        let selected_consumed = self.row + self.column;
        if score > self.score
            || (score == self.score && candidate_consumed > selected_consumed)
            || (score == self.score && candidate_consumed == selected_consumed && row > self.row)
        {
            *self = Self { score, row, column };
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct SimdFillOutcome {
    pub(crate) local_endpoint: SimdEndpoint,
    pub(crate) last_column_endpoint: SimdEndpoint,
    pub(crate) lanes: usize,
}

trait AffineSimdI32: Simd {
    fn shift_lanes_right_with_fill(
        self,
        values: Self::i32s,
        amount: usize,
        fill: i32,
    ) -> Self::i32s;
    fn replace_lane_zero(self, values: Self::i32s, first: i32) -> Self::i32s;
    fn pack_trace_bytes(self, values: Self::i32s) -> [u8; 16];
    fn last_active_i32(self, values: Self::i32s, active: usize) -> i32;
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl AffineSimdI32 for pulp::x86::V3 {
    #[inline(always)]
    fn shift_lanes_right_with_fill(
        self,
        values: Self::i32s,
        amount: usize,
        fill: i32,
    ) -> Self::i32s {
        let indices = match amount {
            1 => pulp::i32x8(7, 0, 1, 2, 3, 4, 5, 6),
            2 => pulp::i32x8(6, 7, 0, 1, 2, 3, 4, 5),
            4 => pulp::i32x8(4, 5, 6, 7, 0, 1, 2, 3),
            _ => unreachable!("the affine prefix scan uses only offsets 1, 2, and 4"),
        };
        let shifted = pulp::cast(
            self.avx2
                ._mm256_permutevar8x32_epi32(pulp::cast(values), pulp::cast(indices)),
        );
        let fill = pulp::cast(self.splat_i32s(fill));
        pulp::cast(match amount {
            1 => self.avx2._mm256_blend_epi32::<0xfe>(fill, shifted),
            2 => self.avx2._mm256_blend_epi32::<0xfc>(fill, shifted),
            4 => self.avx2._mm256_blend_epi32::<0xf0>(fill, shifted),
            _ => unreachable!("the affine prefix scan uses only offsets 1, 2, and 4"),
        })
    }

    #[inline(always)]
    fn replace_lane_zero(self, values: Self::i32s, first: i32) -> Self::i32s {
        pulp::cast(
            self.avx2
                ._mm256_blend_epi32::<0x01>(pulp::cast(values), pulp::cast(self.splat_i32s(first))),
        )
    }

    #[inline(always)]
    fn pack_trace_bytes(self, values: Self::i32s) -> [u8; 16] {
        let byte_indices = pulp::u8x32(
            0, 4, 8, 12, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 0, 4, 8, 12,
            128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128, 128,
        );
        let selected = self
            .avx2
            ._mm256_shuffle_epi8(pulp::cast(values), pulp::cast(byte_indices));
        let lower = self.avx2._mm256_extracti128_si256::<0>(selected);
        let upper = self.avx2._mm256_extracti128_si256::<1>(selected);
        pulp::cast(self.sse2._mm_unpacklo_epi32(lower, upper))
    }

    #[inline(always)]
    fn last_active_i32(self, values: Self::i32s, active: usize) -> i32 {
        debug_assert!((1..=Self::I32_LANES).contains(&active));
        if active == Self::I32_LANES {
            self.avx2._mm256_extract_epi32::<7>(pulp::cast(values))
        } else {
            let lanes: [i32; 8] = pulp::cast(values);
            lanes[active - 1]
        }
    }
}

trait AffineSimdI16: Simd {
    fn shift_lanes_right_with_fill_i16(
        self,
        values: Self::i16s,
        amount: usize,
        fill: i16,
    ) -> Self::i16s;
    fn replace_lane_zero_i16(self, values: Self::i16s, first: i16) -> Self::i16s;
    fn select_i16s(self, mask: Self::m16s, if_true: Self::i16s, if_false: Self::i16s)
    -> Self::i16s;
    fn pack_trace_bytes_i16(self, values: Self::i16s) -> [u8; 16];
    fn last_active_i16(self, values: Self::i16s, active: usize) -> i16;
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl AffineSimdI16 for pulp::x86::V3 {
    #[inline(always)]
    fn shift_lanes_right_with_fill_i16(
        self,
        values: Self::i16s,
        amount: usize,
        fill: i16,
    ) -> Self::i16s {
        let values = pulp::cast(values);
        let fill = pulp::cast(self.splat_i16s(fill));
        let previous_halves = self.avx2._mm256_permute2x128_si256::<0x20>(fill, values);
        pulp::cast(match amount {
            1 => self.avx2._mm256_alignr_epi8::<14>(values, previous_halves),
            2 => self.avx2._mm256_alignr_epi8::<12>(values, previous_halves),
            4 => self.avx2._mm256_alignr_epi8::<8>(values, previous_halves),
            8 => previous_halves,
            _ => unreachable!("the affine prefix scan uses only offsets 1, 2, 4, and 8"),
        })
    }

    #[inline(always)]
    fn replace_lane_zero_i16(self, values: Self::i16s, first: i16) -> Self::i16s {
        pulp::cast(self.avx._mm256_insert_epi16::<0>(pulp::cast(values), first))
    }

    #[inline(always)]
    fn select_i16s(
        self,
        mask: Self::m16s,
        if_true: Self::i16s,
        if_false: Self::i16s,
    ) -> Self::i16s {
        pulp::cast(self.avx2._mm256_blendv_epi8(
            pulp::cast(if_false),
            pulp::cast(if_true),
            pulp::cast(mask),
        ))
    }

    #[inline(always)]
    fn pack_trace_bytes_i16(self, values: Self::i16s) -> [u8; 16] {
        let values = pulp::cast(values);
        let lower = self.avx2._mm256_extracti128_si256::<0>(values);
        let upper = self.avx2._mm256_extracti128_si256::<1>(values);
        pulp::cast(self.sse2._mm_packus_epi16(lower, upper))
    }

    #[inline(always)]
    fn last_active_i16(self, values: Self::i16s, active: usize) -> i16 {
        debug_assert!((1..=Self::I16_LANES).contains(&active));
        if active == Self::I16_LANES {
            let extracted = self.avx2._mm256_extract_epi16::<15>(pulp::cast(values));
            let bits = u16::try_from(extracted).expect("AVX2 extracts an unsigned 16-bit lane");
            i16::from_ne_bytes(bits.to_ne_bytes())
        } else {
            let lanes: [i16; 16] = pulp::cast(values);
            lanes[active - 1]
        }
    }
}

#[derive(Clone, Copy)]
enum SimdScoreWidth {
    I16,
    I32,
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
struct FillContext<'a> {
    first_len: usize,
    gap: AffineGap,
    mode: SimdBoundaryMode,
    score_width: SimdScoreWidth,
    workspace: &'a mut AffineWorkspace,
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl FillContext<'_> {
    #[inline(always)]
    fn run(&mut self, simd: pulp::x86::V3) -> Result<SimdFillOutcome, AlignmentError> {
        match self.score_width {
            SimdScoreWidth::I16 => {
                fill_simd_i16(simd, self.first_len, self.gap, self.mode, self.workspace)
            }
            SimdScoreWidth::I32 => {
                fill_simd_i32(simd, self.first_len, self.gap, self.mode, self.workspace)
            }
        }
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
struct ScoreContext<'a> {
    first_len: usize,
    gap: AffineGap,
    mode: SimdBoundaryMode,
    score_width: SimdScoreWidth,
    workspace: &'a mut AffineWorkspace,
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl ScoreContext<'_> {
    #[inline(always)]
    fn run(&mut self, simd: pulp::x86::V3) -> Result<SimdFillOutcome, AlignmentError> {
        match self.score_width {
            SimdScoreWidth::I16 => {
                fill_score_simd_i16(simd, self.first_len, self.gap, self.mode, self.workspace)
            }
            SimdScoreWidth::I32 => {
                fill_score_simd_i32(simd, self.first_len, self.gap, self.mode, self.workspace)
            }
        }
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
struct V3FillOperation<'a, 'workspace> {
    simd: pulp::x86::V3,
    context: &'a mut FillContext<'workspace>,
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl WithSimd for V3FillOperation<'_, '_> {
    type Output = Result<SimdFillOutcome, AlignmentError>;

    #[inline(always)]
    fn with_simd<V: Simd>(self, _enabled: V) -> Self::Output {
        self.context.run(self.simd)
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
struct V3ScoreOperation<'a, 'workspace> {
    simd: pulp::x86::V3,
    context: &'a mut ScoreContext<'workspace>,
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
impl WithSimd for V3ScoreOperation<'_, '_> {
    type Output = Result<SimdFillOutcome, AlignmentError>;

    #[inline(always)]
    fn with_simd<V: Simd>(self, _enabled: V) -> Self::Output {
        self.context.run(self.simd)
    }
}

pub(crate) fn every_path_score_fits_i32<A, S>(
    substitution: &S,
    gap: AffineGap,
    maximum_operations: usize,
) -> bool
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    every_path_score_fits_for::<A, S>(
        substitution,
        gap,
        maximum_operations,
        i128::from(i32::MIN),
        i128::from(i32::MAX),
    ) && i32_prefix_arithmetic_fits(gap)
}

pub(crate) fn every_path_score_fits_i16<A, S>(
    substitution: &S,
    gap: AffineGap,
    maximum_operations: usize,
) -> bool
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    every_path_score_fits_for::<A, S>(
        substitution,
        gap,
        maximum_operations,
        i128::from(i16::MIN),
        i128::from(i16::MAX),
    ) && i16_prefix_arithmetic_fits(gap)
}

fn i16_prefix_arithmetic_fits(gap: AffineGap) -> bool {
    let Ok(extension) = i16::try_from(gap.extension_score()) else {
        return false;
    };
    let Some(extension_delta) = extension.checked_mul(8) else {
        return false;
    };
    i16::MIN.checked_sub(extension_delta).is_some()
}

fn i32_prefix_arithmetic_fits(gap: AffineGap) -> bool {
    let Ok(extension) = i32::try_from(gap.extension_score()) else {
        return false;
    };
    let Some(extension_delta) = extension.checked_mul(4) else {
        return false;
    };
    i32::MIN.checked_sub(extension_delta).is_some()
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
#[inline(never)]
pub(crate) fn try_fill_simd<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    mode: SimdBoundaryMode,
    maximum_operations: usize,
    workspace: &mut AffineWorkspace,
) -> Result<Option<SimdFillOutcome>, AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let policy = workspace.kernel_policy();
    workspace.clear_last_backend();
    if matches!(policy, AffineKernelPolicy::Scalar) {
        return Ok(None);
    }
    if first.is_empty() || second.is_empty() {
        return if matches!(policy, AffineKernelPolicy::Simd) {
            Err(AlignmentError::SimdUnavailable {
                reason: "empty inputs contain no SIMD recurrence cells",
            })
        } else {
            Ok(None)
        };
    }
    if !every_path_score_fits_i32::<A, S>(substitution, gap, maximum_operations) {
        return if matches!(policy, AffineKernelPolicy::Simd) {
            Err(AlignmentError::SimdUnavailable {
                reason: "not every partial path or prefix-scan value is representable in signed 32-bit lanes",
            })
        } else {
            Ok(None)
        };
    }
    let uses_i16 = every_path_score_fits_i16::<A, S>(substitution, gap, maximum_operations);
    let cells = first.len().saturating_mul(second.len());
    if matches!(policy, AffineKernelPolicy::Auto)
        && !auto_shape_is_eligible(first.len(), cells, uses_i16)
    {
        return Ok(None);
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    let outcome = if let Some(simd) = pulp::x86::V3::try_new() {
        let backing_first_len = padded_length(first.len(), AVX2_I32_LANES).ok_or(
            AlignmentError::DimensionOverflow {
                first_len: first.len(),
                second_len: second.len(),
            },
        )?;
        let rows = first
            .len()
            .checked_add(1)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len: first.len(),
                second_len: second.len(),
            })?;
        let columns = second
            .len()
            .checked_add(1)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len: first.len(),
                second_len: second.len(),
            })?;
        let matrix_cells = rows
            .checked_mul(columns)
            .ok_or(AlignmentError::DimensionOverflow {
                first_len: first.len(),
                second_len: second.len(),
            })?;
        workspace.prepare_simd(
            first.len(),
            backing_first_len,
            second.len(),
            A::SYMBOLS.len(),
            matrix_cells,
        )?;
        let score_width = if uses_i16 {
            let profile_stride = padded_length(first.len(), AVX2_I16_LANES).ok_or(
                AlignmentError::DimensionOverflow {
                    first_len: first.len(),
                    second_len: second.len(),
                },
            )?;
            let profile_elements = profile_stride * A::SYMBOLS.len();
            let narrow_profile: &mut [i16] =
                pulp::bytemuck::cast_slice_mut(&mut workspace.simd_profile);
            build_profile_i16(
                first,
                substitution,
                &mut narrow_profile[..profile_elements],
                profile_stride,
            )?;
            SimdScoreWidth::I16
        } else {
            build_profile_i32(
                first,
                substitution,
                &mut workspace.simd_profile,
                backing_first_len,
            )?;
            SimdScoreWidth::I32
        };
        for (output, symbol) in workspace
            .simd_second_indices
            .iter_mut()
            .zip(second.iter().copied())
        {
            *output = symbol.alphabet_index();
        }
        let mut context = FillContext {
            first_len: first.len(),
            gap,
            mode,
            score_width,
            workspace,
        };
        let outcome = Simd::vectorize(
            simd,
            V3FillOperation {
                simd,
                context: &mut context,
            },
        )?;
        Some((outcome, score_width))
    } else {
        None
    };
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    let outcome = None;

    match outcome {
        None if matches!(policy, AffineKernelPolicy::Simd) => {
            Err(AlignmentError::SimdUnavailable {
                reason: "the current architecture has no enabled SIMD backend",
            })
        }
        None => Ok(None),
        Some((outcome, score_width)) => {
            workspace.set_last_backend(match score_width {
                SimdScoreWidth::I16 => AffineBackend::SimdI16 {
                    lanes: outcome.lanes,
                },
                SimdScoreWidth::I32 => AffineBackend::SimdI32 {
                    lanes: outcome.lanes,
                },
            });
            Ok(Some(outcome))
        }
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
#[inline(always)]
pub(crate) fn try_score_simd<A, S>(
    first: &[Symbol<A>],
    second: &[Symbol<A>],
    substitution: &S,
    gap: AffineGap,
    mode: SimdBoundaryMode,
    maximum_operations: usize,
    workspace: &mut AffineWorkspace,
) -> Result<Option<SimdFillOutcome>, AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let policy = workspace.kernel_policy();
    workspace.clear_last_backend();
    if matches!(policy, AffineKernelPolicy::Scalar) {
        return Ok(None);
    }
    if first.is_empty() || second.is_empty() {
        return if matches!(policy, AffineKernelPolicy::Simd) {
            Err(AlignmentError::SimdUnavailable {
                reason: "empty inputs contain no SIMD recurrence cells",
            })
        } else {
            Ok(None)
        };
    }
    if !every_path_score_fits_i32::<A, S>(substitution, gap, maximum_operations) {
        return if matches!(policy, AffineKernelPolicy::Simd) {
            Err(AlignmentError::SimdUnavailable {
                reason: "not every partial path or prefix-scan value is representable in signed 32-bit lanes",
            })
        } else {
            Ok(None)
        };
    }
    let uses_i16 = every_path_score_fits_i16::<A, S>(substitution, gap, maximum_operations);
    let cells = first.len().saturating_mul(second.len());
    if matches!(policy, AffineKernelPolicy::Auto)
        && !auto_score_shape_is_eligible(first.len(), cells, uses_i16)
    {
        return Ok(None);
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    let outcome = if let Some(simd) = pulp::x86::V3::try_new() {
        let backing_first_len = padded_length(first.len(), AVX2_I32_LANES).ok_or(
            AlignmentError::DimensionOverflow {
                first_len: first.len(),
                second_len: second.len(),
            },
        )?;
        workspace.prepare_simd_score(
            first.len(),
            backing_first_len,
            second.len(),
            A::SYMBOLS.len(),
        )?;
        let (score_width, profile_stride) = if uses_i16 {
            let profile_stride = padded_length(first.len(), AVX2_I16_LANES).ok_or(
                AlignmentError::DimensionOverflow {
                    first_len: first.len(),
                    second_len: second.len(),
                },
            )?;
            let profile_elements = profile_stride * A::SYMBOLS.len();
            let narrow_profile: &mut [i16] =
                pulp::bytemuck::cast_slice_mut(&mut workspace.simd_profile);
            build_profile_i16(
                first,
                substitution,
                &mut narrow_profile[..profile_elements],
                profile_stride,
            )?;
            (SimdScoreWidth::I16, profile_stride)
        } else {
            build_profile_i32(
                first,
                substitution,
                &mut workspace.simd_profile,
                backing_first_len,
            )?;
            (SimdScoreWidth::I32, backing_first_len)
        };
        for (output, symbol) in workspace
            .simd_second_indices
            .iter_mut()
            .zip(second.iter().copied())
        {
            *output = symbol
                .alphabet_index()
                .checked_mul(profile_stride)
                .expect("profile allocation proves every row offset representable");
        }
        let mut context = ScoreContext {
            first_len: first.len(),
            gap,
            mode,
            score_width,
            workspace,
        };
        let outcome = Simd::vectorize(
            simd,
            V3ScoreOperation {
                simd,
                context: &mut context,
            },
        )?;
        Some((outcome, score_width))
    } else {
        None
    };
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    let outcome = None;

    match outcome {
        None if matches!(policy, AffineKernelPolicy::Simd) => {
            Err(AlignmentError::SimdUnavailable {
                reason: "the current architecture has no enabled SIMD backend",
            })
        }
        None => Ok(None),
        Some((outcome, score_width)) => {
            workspace.set_last_backend(match score_width {
                SimdScoreWidth::I16 => AffineBackend::SimdI16 {
                    lanes: outcome.lanes,
                },
                SimdScoreWidth::I32 => AffineBackend::SimdI32 {
                    lanes: outcome.lanes,
                },
            });
            Ok(Some(outcome))
        }
    }
}

fn score_to_i32(score: Score) -> Result<i32, AlignmentError> {
    i32::try_from(score).map_err(|_| AlignmentError::SimdUnavailable {
        reason: "a scorer value contradicted its signed-32-bit preflight bound",
    })
}

fn build_profile_i32<A, S>(
    first: &[Symbol<A>],
    substitution: &S,
    profile: &mut [i32],
    stride: usize,
) -> Result<(), AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let width = A::SYMBOLS.len();
    debug_assert!(stride >= first.len());
    debug_assert_eq!(profile.len(), stride * width);
    profile.fill(0);
    let mut score_row = [0; 256];
    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let cached = width <= score_row.len()
            && substitution.fill_score_row(first_symbol, &mut score_row[..width]);
        for (second_index, second_ascii) in A::SYMBOLS.iter().copied().enumerate() {
            let score = if cached {
                score_row[second_index]
            } else {
                let second_symbol = Symbol::<A>::from_ascii(second_ascii)
                    .expect("alphabet constants contain valid symbols");
                substitution.score(first_symbol, second_symbol)
            };
            profile[second_index * stride + first_index] = score_to_i32(score)?;
        }
    }
    Ok(())
}

fn score_to_i16(score: Score) -> Result<i16, AlignmentError> {
    i16::try_from(score).map_err(|_| AlignmentError::SimdUnavailable {
        reason: "a scorer value contradicted its signed-16-bit preflight bound",
    })
}

fn build_profile_i16<A, S>(
    first: &[Symbol<A>],
    substitution: &S,
    profile: &mut [i16],
    stride: usize,
) -> Result<(), AlignmentError>
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let width = A::SYMBOLS.len();
    debug_assert!(stride >= first.len());
    debug_assert_eq!(profile.len(), stride * width);
    profile.fill(0);
    let mut score_row = [0; 256];
    for (first_index, first_symbol) in first.iter().copied().enumerate() {
        let cached = width <= score_row.len()
            && substitution.fill_score_row(first_symbol, &mut score_row[..width]);
        for (second_index, second_ascii) in A::SYMBOLS.iter().copied().enumerate() {
            let score = if cached {
                score_row[second_index]
            } else {
                let second_symbol = Symbol::<A>::from_ascii(second_ascii)
                    .expect("alphabet constants contain valid symbols");
                substitution.score(first_symbol, second_symbol)
            };
            profile[second_index * stride + first_index] = score_to_i16(score)?;
        }
    }
    Ok(())
}

fn gap_boundary(open: i32, extension: i32, consumed: usize) -> Result<i32, AlignmentError> {
    debug_assert!(consumed > 0);
    let extensions = i128::try_from(consumed - 1).map_err(|_| AlignmentError::SimdUnavailable {
        reason: "a boundary length cannot be represented by the SIMD score preflight",
    })?;
    let score = i128::from(open) + i128::from(extension) * extensions;
    i32::try_from(score).map_err(|_| AlignmentError::SimdUnavailable {
        reason: "an affine boundary score is not representable in signed 32-bit lanes",
    })
}

#[inline(always)]
fn load_partial<V: Simd>(simd: V, values: &[i32]) -> V::i32s {
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
fn store_partial<V: Simd>(simd: V, output: &mut [i32], values: V::i32s) {
    debug_assert!(output.len() <= V::I32_LANES);
    if output.len() == V::I32_LANES {
        let (vectors, remainder) = V::as_mut_simd_i32s(output);
        debug_assert!(remainder.is_empty());
        vectors[0] = values;
    } else {
        simd.partial_store_i32s(output, values);
    }
}

#[inline(always)]
fn store_lanes<V: Simd>(simd: V, values: V::i32s, output: &mut [i32; MAX_I32_LANES]) {
    store_partial(simd, &mut output[..V::I32_LANES], values);
}

#[inline(always)]
fn lane_indices<V: Simd>(simd: V) -> V::i32s {
    let mut lanes = [0_i32; MAX_I32_LANES];
    for (index, lane) in lanes[..V::I32_LANES].iter_mut().enumerate() {
        *lane = i32::try_from(index).expect("SIMD lane count fits i32");
    }
    load_partial(simd, &lanes[..V::I32_LANES])
}

fn shift_in_first<V: AffineSimdI32>(simd: V, values: V::i32s, first: i32) -> V::i32s {
    simd.shift_lanes_right_with_fill(values, 1, first)
}

#[derive(Clone, Copy)]
struct PrefixStepI32 {
    extension_delta: i32,
    inactive_base: i32,
}

impl PrefixStepI32 {
    #[inline(always)]
    fn new(extension: i32, offset: i32) -> Self {
        let extension_delta = extension
            .checked_mul(offset)
            .expect("full-vector score preflight covers prefix extension");
        let inactive_base = i32::MIN
            .checked_sub(extension_delta)
            .expect("non-positive affine extension keeps the scan sentinel representable");
        Self {
            extension_delta,
            inactive_base,
        }
    }
}

#[derive(Clone, Copy)]
struct PrefixScanI32 {
    one: PrefixStepI32,
    two: PrefixStepI32,
    four: PrefixStepI32,
    eight: Option<PrefixStepI32>,
}

impl PrefixScanI32 {
    #[inline(always)]
    fn for_simd<V: Simd>(extension: i32) -> Self {
        Self {
            one: PrefixStepI32::new(extension, 1),
            two: PrefixStepI32::new(extension, 2),
            four: PrefixStepI32::new(extension, 4),
            eight: (V::I32_LANES > 8).then(|| PrefixStepI32::new(extension, 8)),
        }
    }
}

#[inline(always)]
fn prefix_step<const OFFSET: usize, V: AffineSimdI32>(
    simd: V,
    deletion: V::i32s,
    step: PrefixStepI32,
) -> V::i32s {
    let safe_shifted = simd.shift_lanes_right_with_fill(deletion, OFFSET, step.inactive_base);
    let candidate = simd.add_i32s(safe_shifted, simd.splat_i32s(step.extension_delta));
    simd.max_i32s(deletion, candidate)
}

#[inline(always)]
fn prefix_deletion<V: AffineSimdI32>(
    simd: V,
    temporary: V::i32s,
    first_seed: i32,
    open: i32,
    scan: PrefixScanI32,
) -> V::i32s {
    let shifted_temporary = simd.shift_lanes_right_with_fill(temporary, 1, 0);
    let mut deletion = simd.add_i32s(shifted_temporary, simd.splat_i32s(open));
    deletion = simd.replace_lane_zero(deletion, first_seed);
    if V::I32_LANES > 1 {
        deletion = prefix_step::<1, V>(simd, deletion, scan.one);
    }
    if V::I32_LANES > 2 {
        deletion = prefix_step::<2, V>(simd, deletion, scan.two);
    }
    if V::I32_LANES > 4 {
        deletion = prefix_step::<4, V>(simd, deletion, scan.four);
    }
    if V::I32_LANES > 8 {
        deletion = prefix_step::<8, V>(
            simd,
            deletion,
            scan.eight
                .expect("prefix constants cover the selected SIMD lane count"),
        );
    }
    deletion
}

fn gap_boundary_i16(open: i16, extension: i16, consumed: usize) -> Result<i16, AlignmentError> {
    debug_assert!(consumed > 0);
    let extensions = i128::try_from(consumed - 1).map_err(|_| AlignmentError::SimdUnavailable {
        reason: "a boundary length cannot be represented by the SIMD score preflight",
    })?;
    let score = i128::from(open) + i128::from(extension) * extensions;
    i16::try_from(score).map_err(|_| AlignmentError::SimdUnavailable {
        reason: "an affine boundary score is not representable in signed-16-bit lanes",
    })
}

#[inline(always)]
fn load_partial_i16<V: Simd>(simd: V, values: &[i16]) -> V::i16s {
    debug_assert!(values.len() <= V::I16_LANES);
    if values.len() == V::I16_LANES {
        let (vectors, remainder) = V::as_simd_i16s(values);
        debug_assert!(remainder.is_empty());
        vectors[0]
    } else {
        simd.partial_load_i16s(values)
    }
}

#[inline(always)]
fn store_partial_i16<V: Simd>(simd: V, output: &mut [i16], values: V::i16s) {
    debug_assert!(output.len() <= V::I16_LANES);
    if output.len() == V::I16_LANES {
        let (vectors, remainder) = V::as_mut_simd_i16s(output);
        debug_assert!(remainder.is_empty());
        vectors[0] = values;
    } else {
        simd.partial_store_i16s(output, values);
    }
}

#[inline(always)]
fn store_lanes_i16<V: Simd>(simd: V, values: V::i16s, output: &mut [i16; 16]) {
    store_partial_i16(simd, &mut output[..V::I16_LANES], values);
}

#[inline(always)]
fn lane_indices_i16<V: Simd>(simd: V) -> V::i16s {
    let mut lanes = [0_i16; 16];
    for (index, lane) in lanes[..V::I16_LANES].iter_mut().enumerate() {
        *lane = i16::try_from(index).expect("SIMD lane count fits i16");
    }
    load_partial_i16(simd, &lanes[..V::I16_LANES])
}

fn shift_in_first_i16<V: AffineSimdI16>(simd: V, values: V::i16s, first: i16) -> V::i16s {
    simd.shift_lanes_right_with_fill_i16(values, 1, first)
}

#[derive(Clone, Copy)]
struct PrefixStepI16 {
    extension_delta: i16,
    inactive_base: i16,
}

impl PrefixStepI16 {
    #[inline(always)]
    fn new(extension: i16, offset: i16) -> Self {
        let extension_delta = extension
            .checked_mul(offset)
            .expect("signed-16-bit prefix preflight covers the extension delta");
        let inactive_base = i16::MIN
            .checked_sub(extension_delta)
            .expect("non-positive affine extension keeps the scan sentinel representable");
        Self {
            extension_delta,
            inactive_base,
        }
    }
}

#[derive(Clone, Copy)]
struct PrefixScanI16 {
    one: PrefixStepI16,
    two: PrefixStepI16,
    four: PrefixStepI16,
    eight: PrefixStepI16,
}

impl PrefixScanI16 {
    #[inline(always)]
    fn new(extension: i16) -> Self {
        Self {
            one: PrefixStepI16::new(extension, 1),
            two: PrefixStepI16::new(extension, 2),
            four: PrefixStepI16::new(extension, 4),
            eight: PrefixStepI16::new(extension, 8),
        }
    }
}

#[inline(always)]
fn prefix_step_i16<const OFFSET: usize, V: AffineSimdI16>(
    simd: V,
    deletion: V::i16s,
    step: PrefixStepI16,
) -> V::i16s {
    let safe_shifted = simd.shift_lanes_right_with_fill_i16(deletion, OFFSET, step.inactive_base);
    let candidate = simd.add_i16s(safe_shifted, simd.splat_i16s(step.extension_delta));
    simd.max_i16s(deletion, candidate)
}

#[inline(always)]
fn prefix_deletion_i16<V: AffineSimdI16>(
    simd: V,
    temporary: V::i16s,
    first_seed: i16,
    open: i16,
    scan: PrefixScanI16,
) -> V::i16s {
    let shifted_temporary = simd.shift_lanes_right_with_fill_i16(temporary, 1, 0);
    let mut deletion = simd.add_i16s(shifted_temporary, simd.splat_i16s(open));
    deletion = simd.replace_lane_zero_i16(deletion, first_seed);
    if V::I16_LANES > 1 {
        deletion = prefix_step_i16::<1, V>(simd, deletion, scan.one);
    }
    if V::I16_LANES > 2 {
        deletion = prefix_step_i16::<2, V>(simd, deletion, scan.two);
    }
    if V::I16_LANES > 4 {
        deletion = prefix_step_i16::<4, V>(simd, deletion, scan.four);
    }
    if V::I16_LANES > 8 {
        deletion = prefix_step_i16::<8, V>(simd, deletion, scan.eight);
    }
    deletion
}

#[inline(always)]
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_simd_i16<V>(
    simd: V,
    first_len: usize,
    gap: AffineGap,
    mode: SimdBoundaryMode,
    workspace: &mut AffineWorkspace,
) -> Result<SimdFillOutcome, AlignmentError>
where
    V: AffineSimdI16,
{
    let rows = first_len + 1;
    let profile_stride =
        padded_length(first_len, V::I16_LANES).ok_or(AlignmentError::DimensionOverflow {
            first_len,
            second_len: 0,
        })?;
    let AffineWorkspace {
        best_scores: final_scores,
        simd_best_scores: wide_best_scores,
        simd_insertion_scores: wide_insertion_scores,
        simd_profile: wide_profile,
        simd_second_indices: second_indices,
        simd_traceback: traceback,
        ..
    } = workspace;
    let second_len = second_indices.len();
    let narrow_best_scores: &mut [i16] = pulp::bytemuck::cast_slice_mut(wide_best_scores);
    let best_scores = &mut narrow_best_scores[..=profile_stride];
    let narrow_insertion_scores: &mut [i16] = pulp::bytemuck::cast_slice_mut(wide_insertion_scores);
    let insertion_scores = &mut narrow_insertion_scores[..=profile_stride];
    let profile: &[i16] = pulp::bytemuck::cast_slice(wide_profile);

    let open = score_to_i16(gap.open_score())?;
    let extension = score_to_i16(gap.extension_score())?;
    let prefix_scan = PrefixScanI16::new(extension);
    let first_start_is_free = mode.first_start_is_free();
    let second_start_is_free = mode.second_start_is_free();
    let local = mode.is_local();

    best_scores.fill(0);
    insertion_scores.fill(0);
    final_scores.fill(0);

    traceback[0] = TraceCell::default().bits();
    let mut first_boundary = open;
    for row in 1..rows {
        if first_start_is_free {
            best_scores[row] = 0;
            traceback[row] = TraceCell::default().bits();
        } else {
            best_scores[row] = first_boundary;
            traceback[row] = TraceCell::default()
                .with_best(BestTrace::Deletion)
                .with_deletion(if row == 1 {
                    GapTrace::Open
                } else {
                    GapTrace::Extend
                })
                .bits();
            if row < first_len {
                first_boundary += extension;
            }
        }
    }
    final_scores[0] = Score::from(best_scores[first_len]);

    let mut local_endpoint = SimdEndpoint::ORIGIN;
    let mut last_column_endpoint = SimdEndpoint {
        score: if second_start_is_free {
            0
        } else {
            Score::from(gap_boundary_i16(open, extension, second_len)?)
        },
        row: 0,
        column: second_len,
    };
    let indices = lane_indices_i16(simd);
    let zero = simd.splat_i16s(0);
    let nonzero_lane = simd.greater_than_i16s(indices, zero);
    let all_lanes = simd.equal_i16s(indices, indices);
    let no_lanes = simd.equal_i16s(indices, simd.splat_i16s(-1));
    let open_vector = simd.splat_i16s(open);
    let extension_vector = simd.splat_i16s(extension);
    let one_vector = simd.splat_i16s(1);
    let two_vector = simd.splat_i16s(2);
    let three_vector = simd.splat_i16s(3);
    let four_vector = simd.splat_i16s(4);
    let eight_vector = simd.splat_i16s(8);

    let mut top_score = if second_start_is_free { 0 } else { open };
    for (second_index, second_symbol_index) in second_indices.iter().copied().enumerate() {
        let column = second_index + 1;
        let trace_column = column * rows;
        traceback[trace_column] = if second_start_is_free {
            TraceCell::default()
        } else {
            TraceCell::default()
                .with_best(BestTrace::Insertion)
                .with_insertion(if column == 1 {
                    GapTrace::Open
                } else {
                    GapTrace::Extend
                })
        }
        .bits();

        let profile_start = second_symbol_index * profile_stride;
        let profile_end = profile_start + profile_stride;
        let profile = &profile[profile_start..profile_end];
        let mut diagonal_carry = best_scores[0];
        best_scores[0] = top_score;
        insertion_scores[0] = top_score;
        let mut best_carry = top_score;
        let mut deletion_carry = 0_i16;

        let mut first_offset = 0;
        while first_offset < first_len {
            let active = (first_len - first_offset).min(V::I16_LANES);
            let profile_range = first_offset..first_offset + V::I16_LANES;
            let score_range = (first_offset + 1)..=(first_offset + V::I16_LANES);
            let old_best = load_partial_i16(simd, &best_scores[score_range.clone()]);
            let old_insertion = load_partial_i16(simd, &insertion_scores[score_range.clone()]);
            let substitution_scores = load_partial_i16(simd, &profile[profile_range]);
            let diagonal_bases = shift_in_first_i16(simd, old_best, diagonal_carry);
            let diagonal = simd.add_i16s(diagonal_bases, substitution_scores);

            let insertion_open = simd.add_i16s(old_best, open_vector);
            let insertion_extension = simd.add_i16s(old_insertion, extension_vector);
            let insertion_extends = if column > 1 {
                simd.greater_than_or_equal_i16s(insertion_extension, insertion_open)
            } else {
                no_lanes
            };
            let insertion = if column > 1 {
                simd.max_i16s(insertion_open, insertion_extension)
            } else {
                insertion_open
            };

            let mut temporary = simd.max_i16s(diagonal, insertion);
            if local {
                temporary = simd.max_i16s(temporary, zero);
            }

            let first_deletion_open = best_carry + open;
            let first_deletion = if first_offset == 0 {
                first_deletion_open
            } else {
                first_deletion_open.max(deletion_carry + extension)
            };
            let deletion = prefix_deletion_i16(simd, temporary, first_deletion, open, prefix_scan);
            let best = simd.max_i16s(temporary, deletion);

            let previous_deletion = shift_in_first_i16(simd, deletion, deletion_carry);
            let deletion_extension = simd.add_i16s(previous_deletion, extension_vector);
            let extension_allowed = if first_offset == 0 {
                nonzero_lane
            } else {
                all_lanes
            };
            let deletion_extends = simd.and_m16s(
                extension_allowed,
                simd.equal_i16s(deletion, deletion_extension),
            );
            let gap_code =
                simd.select_i16s(simd.equal_i16s(best, deletion), two_vector, three_vector);
            let mut best_code =
                simd.select_i16s(simd.equal_i16s(best, diagonal), one_vector, gap_code);
            if local {
                best_code = simd.select_i16s(simd.greater_than_i16s(best, zero), best_code, zero);
            }
            let deletion_bits = simd.select_i16s(deletion_extends, four_vector, zero);
            let insertion_bits = simd.select_i16s(insertion_extends, eight_vector, zero);
            let trace_bits = simd.or_i16s(simd.or_i16s(best_code, deletion_bits), insertion_bits);

            let next_diagonal_carry = simd.last_active_i16(old_best, active);
            let next_best_carry = simd.last_active_i16(best, active);
            let next_deletion_carry = simd.last_active_i16(deletion, active);
            store_partial_i16(simd, &mut best_scores[score_range.clone()], best);
            store_partial_i16(simd, &mut insertion_scores[score_range], insertion);
            let packed_trace = simd.pack_trace_bytes_i16(trace_bits);
            traceback[(trace_column + first_offset + 1)..=(trace_column + first_offset + active)]
                .copy_from_slice(&packed_trace[..active]);
            if local {
                let mut best_lanes = [0_i16; 16];
                store_lanes_i16(simd, best, &mut best_lanes);
                for (lane, score) in best_lanes[..active].iter().copied().enumerate() {
                    let row = first_offset + lane + 1;
                    local_endpoint.observe_local(Score::from(score), row, column);
                }
            }

            diagonal_carry = next_diagonal_carry;
            best_carry = next_best_carry;
            deletion_carry = next_deletion_carry;
            first_offset += V::I16_LANES;
        }

        final_scores[column] = Score::from(best_scores[first_len]);
        if column == second_len {
            for (row, score) in best_scores.iter().copied().enumerate().take(rows).skip(1) {
                last_column_endpoint.observe_end_gap(Score::from(score), row, column);
            }
        }
        if !second_start_is_free && column < second_len {
            top_score += extension;
        }
    }

    Ok(SimdFillOutcome {
        local_endpoint,
        last_column_endpoint,
        lanes: V::I16_LANES,
    })
}

#[inline(always)]
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn fill_score_simd_i16(
    simd: pulp::x86::V3,
    first_len: usize,
    gap: AffineGap,
    mode: SimdBoundaryMode,
    workspace: &mut AffineWorkspace,
) -> Result<SimdFillOutcome, AlignmentError> {
    if mode.is_global() {
        fill_score_simd_i16_kernel::<true>(simd, first_len, gap, mode, workspace)
    } else {
        fill_score_simd_i16_kernel::<false>(simd, first_len, gap, mode, workspace)
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn fill_score_simd_i16_kernel<const GLOBAL: bool>(
    simd: pulp::x86::V3,
    first_len: usize,
    gap: AffineGap,
    mode: SimdBoundaryMode,
    workspace: &mut AffineWorkspace,
) -> Result<SimdFillOutcome, AlignmentError> {
    let profile_stride =
        padded_length(first_len, AVX2_I16_LANES).ok_or(AlignmentError::DimensionOverflow {
            first_len,
            second_len: 0,
        })?;
    let AffineWorkspace {
        best_scores: final_scores,
        simd_best_scores: wide_best_scores,
        simd_insertion_scores: wide_insertion_scores,
        simd_profile: wide_profile,
        simd_second_indices: second_indices,
        ..
    } = workspace;
    let second_len = second_indices.len();
    let narrow_best_scores: &mut [i16] = pulp::bytemuck::cast_slice_mut(wide_best_scores);
    let best_scores = &mut narrow_best_scores[..profile_stride];
    let narrow_insertion_scores: &mut [i16] = pulp::bytemuck::cast_slice_mut(wide_insertion_scores);
    let insertion_scores = &mut narrow_insertion_scores[..profile_stride];
    let profile: &[i16] = pulp::bytemuck::cast_slice(wide_profile);

    let open = score_to_i16(gap.open_score())?;
    let extension = score_to_i16(gap.extension_score())?;
    let prefix_scan = PrefixScanI16::new(extension);
    let first_start_is_free = !GLOBAL && mode.first_start_is_free();
    let second_start_is_free = !GLOBAL && mode.second_start_is_free();
    let local = !GLOBAL && mode.is_local();

    best_scores.fill(0);
    insertion_scores.fill(0);
    if !GLOBAL {
        final_scores.fill(0);
    }

    let mut first_boundary = open;
    for (offset, score) in best_scores.iter_mut().enumerate().take(first_len) {
        let row = offset + 1;
        if first_start_is_free {
            *score = 0;
        } else {
            *score = first_boundary;
            if row < first_len {
                first_boundary += extension;
            }
        }
    }
    if !GLOBAL {
        final_scores[0] = Score::from(best_scores[first_len - 1]);
    }

    let mut local_endpoint = SimdEndpoint::ORIGIN;
    let mut last_column_endpoint = if GLOBAL {
        SimdEndpoint::ORIGIN
    } else {
        SimdEndpoint {
            score: if second_start_is_free {
                0
            } else {
                Score::from(gap_boundary_i16(open, extension, second_len)?)
            },
            row: 0,
            column: second_len,
        }
    };
    let zero = simd.splat_i16s(0);
    let open_vector = simd.splat_i16s(open);
    let extension_vector = simd.splat_i16s(extension);

    let mut previous_top_score = 0_i16;
    let mut top_score = if second_start_is_free { 0 } else { open };
    for (second_index, profile_start) in second_indices.iter().copied().enumerate() {
        let column = second_index + 1;
        let profile_end = profile_start + profile_stride;
        let profile = &profile[profile_start..profile_end];
        let mut diagonal_carry = previous_top_score;
        let mut best_carry = top_score;
        let mut deletion_carry = 0_i16;

        for (block_index, ((best_chunk, insertion_chunk), profile_chunk)) in best_scores
            .chunks_exact_mut(AVX2_I16_LANES)
            .zip(insertion_scores.chunks_exact_mut(AVX2_I16_LANES))
            .zip(profile.chunks_exact(AVX2_I16_LANES))
            .enumerate()
        {
            let first_offset = block_index * AVX2_I16_LANES;
            let active = if GLOBAL {
                AVX2_I16_LANES
            } else {
                (first_len - first_offset).min(AVX2_I16_LANES)
            };
            let old_best = load_partial_i16(simd, best_chunk);
            let old_insertion = load_partial_i16(simd, insertion_chunk);
            let substitution_scores = load_partial_i16(simd, profile_chunk);
            let diagonal_bases = shift_in_first_i16(simd, old_best, diagonal_carry);
            let diagonal = simd.add_i16s(diagonal_bases, substitution_scores);

            let insertion_open = simd.add_i16s(old_best, open_vector);
            let insertion = if column > 1 {
                let insertion_extension = simd.add_i16s(old_insertion, extension_vector);
                simd.max_i16s(insertion_open, insertion_extension)
            } else {
                insertion_open
            };
            let mut temporary = simd.max_i16s(diagonal, insertion);
            if local {
                temporary = simd.max_i16s(temporary, zero);
            }

            let first_deletion_open = best_carry + open;
            let first_deletion = if first_offset == 0 {
                first_deletion_open
            } else {
                first_deletion_open.max(deletion_carry + extension)
            };
            let deletion = prefix_deletion_i16(simd, temporary, first_deletion, open, prefix_scan);
            let best = simd.max_i16s(temporary, deletion);

            let next_diagonal_carry = simd.last_active_i16(old_best, active);
            let next_best_carry = simd.last_active_i16(best, active);
            let next_deletion_carry = simd.last_active_i16(deletion, active);
            store_partial_i16(simd, best_chunk, best);
            store_partial_i16(simd, insertion_chunk, insertion);
            if local {
                let mut best_lanes = [0_i16; 16];
                store_lanes_i16(simd, best, &mut best_lanes);
                for (lane, score) in best_lanes[..active].iter().copied().enumerate() {
                    local_endpoint.observe_local(
                        Score::from(score),
                        first_offset + lane + 1,
                        column,
                    );
                }
            }

            diagonal_carry = next_diagonal_carry;
            best_carry = next_best_carry;
            deletion_carry = next_deletion_carry;
        }

        if !GLOBAL || column == second_len {
            final_scores[column] = Score::from(best_scores[first_len - 1]);
        }
        if !GLOBAL && column == second_len {
            for (offset, score) in best_scores.iter().copied().enumerate().take(first_len) {
                last_column_endpoint.observe_end_gap(Score::from(score), offset + 1, column);
            }
        }
        previous_top_score = top_score;
        if !second_start_is_free && column < second_len {
            top_score += extension;
        }
    }

    Ok(SimdFillOutcome {
        local_endpoint,
        last_column_endpoint,
        lanes: AVX2_I16_LANES,
    })
}

#[inline(always)]
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_simd_i32<V>(
    simd: V,
    first_len: usize,
    gap: AffineGap,
    mode: SimdBoundaryMode,
    workspace: &mut AffineWorkspace,
) -> Result<SimdFillOutcome, AlignmentError>
where
    V: AffineSimdI32,
{
    let second_len = workspace.simd_second_indices.len();
    let rows = first_len + 1;
    let profile_stride =
        padded_length(first_len, V::I32_LANES).ok_or(AlignmentError::DimensionOverflow {
            first_len,
            second_len,
        })?;

    let open = score_to_i32(gap.open_score())?;
    let extension = score_to_i32(gap.extension_score())?;
    let prefix_scan = PrefixScanI32::for_simd::<V>(extension);
    let first_start_is_free = mode.first_start_is_free();
    let second_start_is_free = mode.second_start_is_free();
    let local = mode.is_local();

    workspace.simd_best_scores.fill(0);
    workspace.simd_insertion_scores.fill(0);
    workspace.best_scores.fill(0);

    workspace.simd_traceback[0] = TraceCell::default().bits();
    let mut first_boundary = open;
    for row in 1..rows {
        if first_start_is_free {
            workspace.simd_best_scores[row] = 0;
            workspace.simd_traceback[row] = TraceCell::default().bits();
        } else {
            workspace.simd_best_scores[row] = first_boundary;
            workspace.simd_traceback[row] = TraceCell::default()
                .with_best(BestTrace::Deletion)
                .with_deletion(if row == 1 {
                    GapTrace::Open
                } else {
                    GapTrace::Extend
                })
                .bits();
            if row < first_len {
                first_boundary += extension;
            }
        }
    }
    workspace.best_scores[0] = Score::from(workspace.simd_best_scores[first_len]);

    let mut local_endpoint = SimdEndpoint::ORIGIN;
    let mut last_column_endpoint = SimdEndpoint {
        score: if second_start_is_free {
            0
        } else {
            Score::from(gap_boundary(open, extension, second_len)?)
        },
        row: 0,
        column: second_len,
    };
    let indices = lane_indices(simd);
    let zero = simd.splat_i32s(0);
    let nonzero_lane = simd.greater_than_i32s(indices, zero);
    let no_lanes = simd.equal_i32s(indices, simd.splat_i32s(-1));
    let open_vector = simd.splat_i32s(open);
    let extension_vector = simd.splat_i32s(extension);
    let one_vector = simd.splat_i32s(1);
    let two_vector = simd.splat_i32s(2);
    let three_vector = simd.splat_i32s(3);
    let four_vector = simd.splat_i32s(4);
    let eight_vector = simd.splat_i32s(8);

    let mut top_score = if second_start_is_free { 0 } else { open };
    for second_index in 0..second_len {
        let column = second_index + 1;
        let trace_column = column * rows;
        workspace.simd_traceback[trace_column] = if second_start_is_free {
            TraceCell::default()
        } else {
            TraceCell::default()
                .with_best(BestTrace::Insertion)
                .with_insertion(if column == 1 {
                    GapTrace::Open
                } else {
                    GapTrace::Extend
                })
        }
        .bits();

        let profile_start = workspace.simd_second_indices[second_index] * profile_stride;
        let profile_end = profile_start + profile_stride;
        let profile = &workspace.simd_profile[profile_start..profile_end];
        let mut diagonal_carry = workspace.simd_best_scores[0];
        workspace.simd_best_scores[0] = top_score;
        workspace.simd_insertion_scores[0] = top_score;
        let mut best_carry = top_score;
        let mut deletion_carry = 0_i32;

        let mut first_offset = 0;
        while first_offset < first_len {
            let active = (first_len - first_offset).min(V::I32_LANES);
            let profile_range = first_offset..first_offset + V::I32_LANES;
            let score_range = (first_offset + 1)..=(first_offset + V::I32_LANES);
            let old_best = load_partial(simd, &workspace.simd_best_scores[score_range.clone()]);
            let old_insertion =
                load_partial(simd, &workspace.simd_insertion_scores[score_range.clone()]);
            let substitution_scores = load_partial(simd, &profile[profile_range]);
            let diagonal_bases = shift_in_first(simd, old_best, diagonal_carry);
            let diagonal = simd.add_i32s(diagonal_bases, substitution_scores);

            let insertion_open = simd.add_i32s(old_best, open_vector);
            let insertion_extension = simd.add_i32s(old_insertion, extension_vector);
            let insertion_extends = if column > 1 {
                simd.greater_than_or_equal_i32s(insertion_extension, insertion_open)
            } else {
                no_lanes
            };
            let insertion = if column > 1 {
                simd.max_i32s(insertion_open, insertion_extension)
            } else {
                insertion_open
            };

            let mut temporary = simd.max_i32s(diagonal, insertion);
            if local {
                temporary = simd.max_i32s(temporary, zero);
            }

            let first_deletion_open = best_carry + open;
            let first_deletion = if first_offset == 0 {
                first_deletion_open
            } else {
                first_deletion_open.max(deletion_carry + extension)
            };
            let deletion = prefix_deletion(simd, temporary, first_deletion, open, prefix_scan);
            let best = simd.max_i32s(temporary, deletion);

            let previous_deletion = shift_in_first(simd, deletion, deletion_carry);
            let deletion_extension = simd.add_i32s(previous_deletion, extension_vector);
            let extension_allowed = if first_offset == 0 {
                nonzero_lane
            } else {
                simd.greater_than_or_equal_i32s(indices, zero)
            };
            let deletion_extends = simd.and_m32s(
                extension_allowed,
                simd.equal_i32s(deletion, deletion_extension),
            );
            let gap_code =
                simd.select_i32s(simd.equal_i32s(best, deletion), two_vector, three_vector);
            let mut best_code =
                simd.select_i32s(simd.equal_i32s(best, diagonal), one_vector, gap_code);
            if local {
                best_code = simd.select_i32s(simd.greater_than_i32s(best, zero), best_code, zero);
            }
            let deletion_bits = simd.select_i32s(deletion_extends, four_vector, zero);
            let insertion_bits = simd.select_i32s(insertion_extends, eight_vector, zero);
            let trace_bits = simd.or_i32s(simd.or_i32s(best_code, deletion_bits), insertion_bits);

            let next_diagonal_carry = simd.last_active_i32(old_best, active);
            let next_best_carry = simd.last_active_i32(best, active);
            let next_deletion_carry = simd.last_active_i32(deletion, active);
            store_partial(
                simd,
                &mut workspace.simd_best_scores[score_range.clone()],
                best,
            );
            store_partial(
                simd,
                &mut workspace.simd_insertion_scores[score_range],
                insertion,
            );
            let packed_trace = simd.pack_trace_bytes(trace_bits);
            workspace.simd_traceback
                [(trace_column + first_offset + 1)..=(trace_column + first_offset + active)]
                .copy_from_slice(&packed_trace[..active]);
            if local {
                let mut best_lanes = [0_i32; MAX_I32_LANES];
                store_lanes(simd, best, &mut best_lanes);
                for (lane, score) in best_lanes[..active].iter().copied().enumerate() {
                    let row = first_offset + lane + 1;
                    local_endpoint.observe_local(Score::from(score), row, column);
                }
            }

            diagonal_carry = next_diagonal_carry;
            best_carry = next_best_carry;
            deletion_carry = next_deletion_carry;
            first_offset += V::I32_LANES;
        }

        workspace.best_scores[column] = Score::from(workspace.simd_best_scores[first_len]);
        if column == second_len {
            for row in 1..rows {
                last_column_endpoint.observe_end_gap(
                    Score::from(workspace.simd_best_scores[row]),
                    row,
                    column,
                );
            }
        }
        if !second_start_is_free && column < second_len {
            top_score += extension;
        }
    }

    Ok(SimdFillOutcome {
        local_endpoint,
        last_column_endpoint,
        lanes: V::I32_LANES,
    })
}

#[inline(always)]
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn fill_score_simd_i32(
    simd: pulp::x86::V3,
    first_len: usize,
    gap: AffineGap,
    mode: SimdBoundaryMode,
    workspace: &mut AffineWorkspace,
) -> Result<SimdFillOutcome, AlignmentError> {
    if mode.is_global() {
        fill_score_simd_i32_kernel::<true>(simd, first_len, gap, mode, workspace)
    } else {
        fill_score_simd_i32_kernel::<false>(simd, first_len, gap, mode, workspace)
    }
}

#[inline(always)]
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
fn fill_score_simd_i32_kernel<const GLOBAL: bool>(
    simd: pulp::x86::V3,
    first_len: usize,
    gap: AffineGap,
    mode: SimdBoundaryMode,
    workspace: &mut AffineWorkspace,
) -> Result<SimdFillOutcome, AlignmentError> {
    let profile_stride =
        padded_length(first_len, AVX2_I32_LANES).ok_or(AlignmentError::DimensionOverflow {
            first_len,
            second_len: 0,
        })?;
    let AffineWorkspace {
        best_scores: final_scores,
        simd_best_scores: wide_best_scores,
        simd_insertion_scores: wide_insertion_scores,
        simd_profile: wide_profile,
        simd_second_indices: second_indices,
        ..
    } = workspace;
    let second_len = second_indices.len();
    let best_scores = &mut wide_best_scores[..profile_stride];
    let insertion_scores = &mut wide_insertion_scores[..profile_stride];
    let profile = &wide_profile[..];

    let open = score_to_i32(gap.open_score())?;
    let extension = score_to_i32(gap.extension_score())?;
    let prefix_scan = PrefixScanI32::for_simd::<pulp::x86::V3>(extension);
    let first_start_is_free = !GLOBAL && mode.first_start_is_free();
    let second_start_is_free = !GLOBAL && mode.second_start_is_free();
    let local = !GLOBAL && mode.is_local();

    best_scores.fill(0);
    insertion_scores.fill(0);
    if !GLOBAL {
        final_scores.fill(0);
    }
    let mut first_boundary = open;
    for (offset, score) in best_scores.iter_mut().enumerate().take(first_len) {
        let row = offset + 1;
        if first_start_is_free {
            *score = 0;
        } else {
            *score = first_boundary;
            if row < first_len {
                first_boundary += extension;
            }
        }
    }
    if !GLOBAL {
        final_scores[0] = Score::from(best_scores[first_len - 1]);
    }

    let mut local_endpoint = SimdEndpoint::ORIGIN;
    let mut last_column_endpoint = if GLOBAL {
        SimdEndpoint::ORIGIN
    } else {
        SimdEndpoint {
            score: if second_start_is_free {
                0
            } else {
                Score::from(gap_boundary(open, extension, second_len)?)
            },
            row: 0,
            column: second_len,
        }
    };
    let zero = simd.splat_i32s(0);
    let open_vector = simd.splat_i32s(open);
    let extension_vector = simd.splat_i32s(extension);

    let mut previous_top_score = 0_i32;
    let mut top_score = if second_start_is_free { 0 } else { open };
    for (second_index, profile_start) in second_indices.iter().copied().enumerate() {
        let column = second_index + 1;
        let profile_end = profile_start + profile_stride;
        let profile = &profile[profile_start..profile_end];
        let mut diagonal_carry = previous_top_score;
        let mut best_carry = top_score;
        let mut deletion_carry = 0_i32;

        for (block_index, ((best_chunk, insertion_chunk), profile_chunk)) in best_scores
            .chunks_exact_mut(AVX2_I32_LANES)
            .zip(insertion_scores.chunks_exact_mut(AVX2_I32_LANES))
            .zip(profile.chunks_exact(AVX2_I32_LANES))
            .enumerate()
        {
            let first_offset = block_index * AVX2_I32_LANES;
            let active = if GLOBAL {
                AVX2_I32_LANES
            } else {
                (first_len - first_offset).min(AVX2_I32_LANES)
            };
            let old_best = load_partial(simd, best_chunk);
            let old_insertion = load_partial(simd, insertion_chunk);
            let substitution_scores = load_partial(simd, profile_chunk);
            let diagonal_bases = shift_in_first(simd, old_best, diagonal_carry);
            let diagonal = simd.add_i32s(diagonal_bases, substitution_scores);

            let insertion_open = simd.add_i32s(old_best, open_vector);
            let insertion = if column > 1 {
                let insertion_extension = simd.add_i32s(old_insertion, extension_vector);
                simd.max_i32s(insertion_open, insertion_extension)
            } else {
                insertion_open
            };
            let mut temporary = simd.max_i32s(diagonal, insertion);
            if local {
                temporary = simd.max_i32s(temporary, zero);
            }

            let first_deletion_open = best_carry + open;
            let first_deletion = if first_offset == 0 {
                first_deletion_open
            } else {
                first_deletion_open.max(deletion_carry + extension)
            };
            let deletion = prefix_deletion(simd, temporary, first_deletion, open, prefix_scan);
            let best = simd.max_i32s(temporary, deletion);

            let next_diagonal_carry = simd.last_active_i32(old_best, active);
            let next_best_carry = simd.last_active_i32(best, active);
            let next_deletion_carry = simd.last_active_i32(deletion, active);
            store_partial(simd, best_chunk, best);
            store_partial(simd, insertion_chunk, insertion);
            if local {
                let mut best_lanes = [0_i32; MAX_I32_LANES];
                store_lanes(simd, best, &mut best_lanes);
                for (lane, score) in best_lanes[..active].iter().copied().enumerate() {
                    local_endpoint.observe_local(
                        Score::from(score),
                        first_offset + lane + 1,
                        column,
                    );
                }
            }

            diagonal_carry = next_diagonal_carry;
            best_carry = next_best_carry;
            deletion_carry = next_deletion_carry;
        }

        if !GLOBAL || column == second_len {
            final_scores[column] = Score::from(best_scores[first_len - 1]);
        }
        if !GLOBAL && column == second_len {
            for (offset, score) in best_scores.iter().copied().enumerate().take(first_len) {
                last_column_endpoint.observe_end_gap(Score::from(score), offset + 1, column);
            }
        }
        previous_top_score = top_score;
        if !second_start_is_free && column < second_len {
            top_score += extension;
        }
    }

    Ok(SimdFillOutcome {
        local_endpoint,
        last_column_endpoint,
        lanes: AVX2_I32_LANES,
    })
}

#[cfg(test)]
mod tests {
    use align_core::{AffineGap, Dna, MatchMismatch, Score};

    use super::every_path_score_fits_i32;

    #[test]
    fn signed_32_bit_preflight_is_stricter_than_score_preflight() {
        let gap = AffineGap::new(-5, -1).expect("valid gap");
        assert!(every_path_score_fits_i32::<Dna, _>(
            &MatchMismatch::new(2, -3),
            gap,
            1_000,
        ));
        assert!(!every_path_score_fits_i32::<Dna, _>(
            &MatchMismatch::new(Score::from(i32::MAX), -1),
            gap,
            2,
        ));
    }
}
