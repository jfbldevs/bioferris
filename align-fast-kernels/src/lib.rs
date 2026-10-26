//! Audited architecture-specific kernels behind safe alignment APIs.
//!
//! This crate is intentionally small: `align-fast` remains `unsafe`-free while
//! the operations that require architecture intrinsics live behind validated,
//! safe slice interfaces here.

#![deny(unsafe_op_in_unsafe_fn)]

use core::mem::size_of;

/// Inclusive diagonal interval containing reachable offsets.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DiagonalRange {
    lower: i32,
    upper: i32,
}

/// Reachable interval and total equal-symbol extensions produced by one call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WavefrontExtensionResult {
    range: Option<DiagonalRange>,
    matches: usize,
}

impl WavefrontExtensionResult {
    /// Returns the inclusive interval containing reachable offsets.
    #[must_use]
    pub const fn range(self) -> Option<DiagonalRange> {
        self.range
    }

    /// Returns the total number of equal symbols consumed across all diagonals.
    #[must_use]
    pub const fn matches(self) -> usize {
        self.matches
    }
}

impl DiagonalRange {
    /// Returns the inclusive lower diagonal.
    #[must_use]
    pub const fn lower(self) -> i32 {
        self.lower
    }

    /// Returns the inclusive upper diagonal.
    #[must_use]
    pub const fn upper(self) -> i32 {
        self.upper
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ExtensionKind {
    Portable,
    #[cfg(target_arch = "x86_64")]
    Avx2,
}

/// Runtime-selected wavefront match-extension kernel.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WavefrontExtension {
    kind: ExtensionKind,
}

impl Default for WavefrontExtension {
    fn default() -> Self {
        Self::new()
    }
}

impl WavefrontExtension {
    /// Detects the best supported extension kernel for the current process.
    #[must_use]
    pub fn new() -> Self {
        #[cfg(target_arch = "x86_64")]
        if std::arch::is_x86_feature_detected!("avx2") {
            return Self {
                kind: ExtensionKind::Avx2,
            };
        }
        Self {
            kind: ExtensionKind::Portable,
        }
    }

    /// Returns whether this instance batches eight diagonals with AVX2.
    #[must_use]
    pub const fn is_vectorized(self) -> bool {
        match self.kind {
            ExtensionKind::Portable => false,
            #[cfg(target_arch = "x86_64")]
            ExtensionKind::Avx2 => true,
        }
    }

    /// Extends every reachable offset over equal bytes in place.
    ///
    /// `lower` is the diagonal corresponding to `offsets[0]`. Negative offsets
    /// are unreachable and remain unchanged. Coordinates outside either input
    /// are also left unchanged, keeping this safe API robust to malformed
    /// callers while the alignment kernel supplies only legal coordinates.
    pub fn extend(
        self,
        first: &[u8],
        second: &[u8],
        lower: i32,
        offsets: &mut [i32],
    ) -> Option<DiagonalRange> {
        self.extend_inner::<false>(first, second, lower, offsets)
            .range
    }

    /// Extends offsets and also counts all consumed equal symbols.
    pub fn extend_with_matches(
        self,
        first: &[u8],
        second: &[u8],
        lower: i32,
        offsets: &mut [i32],
    ) -> WavefrontExtensionResult {
        self.extend_inner::<true>(first, second, lower, offsets)
    }

    fn extend_inner<const TRACK_MATCHES: bool>(
        self,
        first: &[u8],
        second: &[u8],
        lower: i32,
        offsets: &mut [i32],
    ) -> WavefrontExtensionResult {
        let interval_fits = offsets.is_empty()
            || i32::try_from(offsets.len() - 1)
                .ok()
                .and_then(|last| lower.checked_add(last))
                .is_some();
        if !interval_fits
            || i32::try_from(first.len()).is_err()
            || i32::try_from(second.len()).is_err()
        {
            return extend_portable_inner::<TRACK_MATCHES>(first, second, lower, offsets);
        }
        match self.kind {
            ExtensionKind::Portable => {
                extend_portable_inner::<TRACK_MATCHES>(first, second, lower, offsets)
            }
            #[cfg(target_arch = "x86_64")]
            ExtensionKind::Avx2 => {
                // SAFETY: construction selected AVX2 only after runtime feature
                // detection. The callee masks every gather to in-bounds 4-byte
                // windows and handles all remaining lanes with safe slices.
                if first.len().min(second.len()) >= x86_64::MINIMUM_CONTIGUOUS_AVX2_LENGTH {
                    unsafe {
                        x86_64::extend_avx2::<TRACK_MATCHES, true>(first, second, lower, offsets)
                    }
                } else {
                    unsafe {
                        x86_64::extend_avx2::<TRACK_MATCHES, false>(first, second, lower, offsets)
                    }
                }
            }
        }
    }
}

#[inline]
fn extend_one(first: &[u8], second: &[u8], diagonal: i32, offset: i32) -> i32 {
    let Ok(mut horizontal) = usize::try_from(offset) else {
        return offset;
    };
    let Some(vertical_offset) = offset.checked_sub(diagonal) else {
        return offset;
    };
    let Ok(mut vertical) = usize::try_from(vertical_offset) else {
        return offset;
    };
    if vertical > first.len() || horizontal > second.len() {
        return offset;
    }

    let remaining = (first.len() - vertical).min(second.len() - horizontal);
    let word_bytes = remaining - remaining % size_of::<u64>();
    let first_words = first[vertical..vertical + word_bytes].chunks_exact(size_of::<u64>());
    let second_words = second[horizontal..horizontal + word_bytes].chunks_exact(size_of::<u64>());
    for (word_index, (first_word, second_word)) in first_words.zip(second_words).enumerate() {
        let first_word = u64::from_le_bytes(
            first_word
                .try_into()
                .expect("an exact word chunk contains eight bytes"),
        );
        let second_word = u64::from_le_bytes(
            second_word
                .try_into()
                .expect("an exact word chunk contains eight bytes"),
        );
        let difference = first_word ^ second_word;
        if difference != 0 {
            let equal = usize::try_from(difference.trailing_zeros() / u8::BITS)
                .expect("a word prefix fits usize");
            return i32::try_from(horizontal + word_index * size_of::<u64>() + equal)
                .expect("validated wavefront dimensions fit i32");
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
    i32::try_from(horizontal).expect("validated wavefront dimensions fit i32")
}

#[cfg(test)]
fn extend_portable(
    first: &[u8],
    second: &[u8],
    lower: i32,
    offsets: &mut [i32],
) -> Option<DiagonalRange> {
    extend_portable_inner::<false>(first, second, lower, offsets).range()
}

fn extend_portable_inner<const TRACK_MATCHES: bool>(
    first: &[u8],
    second: &[u8],
    lower: i32,
    offsets: &mut [i32],
) -> WavefrontExtensionResult {
    let mut active_lower = i32::MAX;
    let mut active_upper = i32::MIN;
    let mut matches = 0_usize;
    for (index, offset) in offsets.iter_mut().enumerate() {
        if *offset < 0 {
            continue;
        }
        let Some(diagonal) = i32::try_from(index)
            .ok()
            .and_then(|index| lower.checked_add(index))
        else {
            continue;
        };
        let source = *offset;
        *offset = extend_one(first, second, diagonal, source);
        if TRACK_MATCHES {
            matches = matches.saturating_add(
                usize::try_from(*offset - source).expect("an extension length is nonnegative"),
            );
        }
        active_lower = active_lower.min(diagonal);
        active_upper = diagonal;
    }
    WavefrontExtensionResult {
        range: (active_lower <= active_upper).then_some(DiagonalRange {
            lower: active_lower,
            upper: active_upper,
        }),
        matches,
    }
}

#[cfg(target_arch = "x86_64")]
mod x86_64 {
    use std::arch::x86_64::{
        __m256i, _mm256_add_epi32, _mm256_add_epi64, _mm256_and_si256, _mm256_andnot_si256,
        _mm256_castps_si256, _mm256_castsi256_ps, _mm256_cmpeq_epi8, _mm256_cmpeq_epi32,
        _mm256_cmpgt_epi32, _mm256_cvtepi32_ps, _mm256_loadu_si256, _mm256_mask_i32gather_epi32,
        _mm256_min_epi16, _mm256_movemask_epi8, _mm256_movemask_ps, _mm256_sad_epu8,
        _mm256_set1_epi32, _mm256_setr_epi8, _mm256_setr_epi32, _mm256_setzero_si256,
        _mm256_shuffle_epi8, _mm256_srli_epi32, _mm256_storeu_si256, _mm256_sub_epi32,
        _mm256_subs_epu16, _mm256_xor_si256,
    };

    use super::{DiagonalRange, WavefrontExtensionResult, extend_one};

    const LANES: usize = 8;
    pub(super) const MINIMUM_CONTIGUOUS_AVX2_LENGTH: usize = 8_192;

    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn leading_zeros_epi32(value: __m256i) -> __m256i {
        // AVX2 has no integer leading-zero instruction. Keep the most
        // significant nonzero byte and recover its exponent through the exact
        // integer-to-float conversion, matching WFA2's AVX2 technique.
        // SAFETY: this helper is only called from the AVX2 entry point, and all
        // operations below are register-only AVX2 intrinsics.
        unsafe {
            let most_significant = _mm256_andnot_si256(_mm256_srli_epi32::<8>(value), value);
            let exponent =
                _mm256_srli_epi32::<23>(_mm256_castps_si256(_mm256_cvtepi32_ps(most_significant)));
            let leading = _mm256_subs_epu16(_mm256_set1_epi32(158), exponent);
            _mm256_min_epi16(leading, _mm256_set1_epi32(32))
        }
    }

    #[inline]
    #[target_feature(enable = "avx2")]
    #[allow(clippy::cast_ptr_alignment)]
    unsafe fn extend_one_avx2(first: &[u8], second: &[u8], diagonal: i32, offset: i32) -> i32 {
        let Ok(horizontal) = usize::try_from(offset) else {
            return offset;
        };
        let Some(vertical_offset) = offset.checked_sub(diagonal) else {
            return offset;
        };
        let Ok(vertical) = usize::try_from(vertical_offset) else {
            return offset;
        };
        if vertical > first.len() || horizontal > second.len() {
            return offset;
        }

        let remaining = (first.len() - vertical).min(second.len() - horizontal);
        let vector_bytes = remaining / size_of::<__m256i>() * size_of::<__m256i>();
        let mut consumed = 0_usize;
        while consumed < vector_bytes {
            // SAFETY: `vector_bytes` is the complete 32-byte prefix inside both
            // remaining slices. Unaligned loads accept arbitrary byte alignment.
            let first_block = unsafe {
                _mm256_loadu_si256(first.as_ptr().add(vertical + consumed).cast::<__m256i>())
            };
            // SAFETY: the same in-bounds argument applies to `second`.
            let second_block = unsafe {
                _mm256_loadu_si256(second.as_ptr().add(horizontal + consumed).cast::<__m256i>())
            };
            // SAFETY: the function's target-feature contract enables AVX2; these
            // operations use registers only after the checked loads above.
            let equal = unsafe { _mm256_cmpeq_epi8(first_block, second_block) };
            // SAFETY: identical AVX2 target-feature argument; movemask has no
            // memory operand.
            let equal_bits =
                u32::from_ne_bytes(unsafe { _mm256_movemask_epi8(equal).to_ne_bytes() });
            if equal_bits != u32::MAX {
                let matching = usize::try_from((!equal_bits).trailing_zeros())
                    .expect("an AVX2 byte prefix fits usize");
                return i32::try_from(horizontal + consumed + matching)
                    .expect("validated wavefront dimensions fit i32");
            }
            consumed += size_of::<__m256i>();
        }

        extend_one(
            first,
            second,
            diagonal,
            i32::try_from(horizontal + consumed).expect("validated wavefront dimensions fit i32"),
        )
    }

    #[target_feature(enable = "avx2")]
    #[allow(clippy::cast_ptr_alignment)]
    #[allow(clippy::too_many_lines)]
    pub(super) unsafe fn extend_avx2<const TRACK_MATCHES: bool, const CONTIGUOUS_AVX2: bool>(
        first: &[u8],
        second: &[u8],
        lower: i32,
        offsets: &mut [i32],
    ) -> WavefrontExtensionResult {
        // SAFETY: the caller established AVX2 support. Pointer-bearing
        // operations are constrained to full output chunks or masked to
        // four-byte windows proven inside the input slices.
        unsafe {
            let first_len = i32::try_from(first.len()).expect("AVX2 dimensions were validated");
            let second_len = i32::try_from(second.len()).expect("AVX2 dimensions were validated");
            let first_four_limit = _mm256_set1_epi32(first_len.saturating_sub(3));
            let second_four_limit = _mm256_set1_epi32(second_len.saturating_sub(3));
            let negative_one = _mm256_set1_epi32(-1);
            let zero = _mm256_setzero_si256();
            let lane_indices = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
            let reverse_bytes = _mm256_setr_epi8(
                3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12, 3, 2, 1, 0, 7, 6, 5, 4, 11,
                10, 9, 8, 15, 14, 13, 12,
            );
            let mut active_lower = i32::MAX;
            let mut active_upper = i32::MIN;
            let mut matches = 0_usize;
            let mut vector_matches = zero;
            let vector_length = offsets.len() / LANES * LANES;

            for start in (0..vector_length).step_by(LANES) {
                let diagonal_base = lower
                    .checked_add(i32::try_from(start).expect("validated interval fits i32"))
                    .expect("validated interval fits i32");
                let diagonals = _mm256_add_epi32(_mm256_set1_epi32(diagonal_base), lane_indices);
                // SAFETY: `start..start + LANES` is inside `offsets` by construction;
                // the unaligned AVX2 load accepts the slice's natural alignment.
                let current = _mm256_loadu_si256(offsets.as_ptr().add(start).cast::<__m256i>());
                let vertical = _mm256_sub_epi32(current, diagonals);
                let reachable = _mm256_cmpgt_epi32(current, negative_one);
                let vertical_nonnegative = _mm256_cmpgt_epi32(vertical, negative_one);
                let horizontal_has_four = _mm256_cmpgt_epi32(second_four_limit, current);
                let vertical_has_four = _mm256_cmpgt_epi32(first_four_limit, vertical);
                let has_four = _mm256_and_si256(
                    reachable,
                    _mm256_and_si256(
                        vertical_nonnegative,
                        _mm256_and_si256(horizontal_has_four, vertical_has_four),
                    ),
                );
                let reachable_bits =
                    u32::try_from(_mm256_movemask_ps(_mm256_castsi256_ps(reachable)))
                        .expect("an eight-lane mask is nonnegative");
                let has_four_bits =
                    u32::try_from(_mm256_movemask_ps(_mm256_castsi256_ps(has_four)))
                        .expect("an eight-lane mask is nonnegative");

                let mut scalar_bits = reachable_bits & !has_four_bits;
                if has_four_bits != 0 {
                    // SAFETY: `has_four` masks every lane whose four-byte window is
                    // outside either slice. AVX2 masked gathers do not access lanes
                    // whose mask sign bit is clear, and scale 1 indexes bytes.
                    let first_words = _mm256_mask_i32gather_epi32::<1>(
                        zero,
                        first.as_ptr().cast::<i32>(),
                        vertical,
                        has_four,
                    );
                    // SAFETY: identical bound argument for `second` and `current`.
                    let second_words = _mm256_mask_i32gather_epi32::<1>(
                        zero,
                        second.as_ptr().cast::<i32>(),
                        current,
                        has_four,
                    );
                    let difference = _mm256_xor_si256(first_words, second_words);
                    let reversed = _mm256_shuffle_epi8(difference, reverse_bytes);
                    let equal_chars = _mm256_and_si256(
                        has_four,
                        _mm256_srli_epi32::<3>(leading_zeros_epi32(reversed)),
                    );
                    if TRACK_MATCHES {
                        vector_matches =
                            _mm256_add_epi64(vector_matches, _mm256_sad_epu8(equal_chars, zero));
                    }
                    let updated = _mm256_add_epi32(current, equal_chars);
                    // SAFETY: the vector store covers the same in-bounds chunk as
                    // the load above and accepts unaligned output.
                    _mm256_storeu_si256(offsets.as_mut_ptr().add(start).cast::<__m256i>(), updated);
                    let all_four_equal =
                        _mm256_and_si256(has_four, _mm256_cmpeq_epi32(difference, zero));
                    scalar_bits |=
                        u32::try_from(_mm256_movemask_ps(_mm256_castsi256_ps(all_four_equal)))
                            .expect("an eight-lane mask is nonnegative");
                }

                if reachable_bits != 0 {
                    let first_lane = i32::try_from(reachable_bits.trailing_zeros())
                        .expect("an AVX2 lane fits i32");
                    let last_lane = i32::try_from(u32::BITS - 1 - reachable_bits.leading_zeros())
                        .expect("an AVX2 lane fits i32");
                    active_lower = active_lower.min(diagonal_base + first_lane);
                    active_upper = diagonal_base + last_lane;
                }
                while scalar_bits != 0 {
                    let lane = usize::try_from(scalar_bits.trailing_zeros())
                        .expect("an AVX2 lane fits usize");
                    let index = start + lane;
                    let diagonal = diagonal_base + i32::try_from(lane).expect("lane fits i32");
                    let source = offsets[index];
                    offsets[index] = if CONTIGUOUS_AVX2 {
                        extend_one_avx2(first, second, diagonal, source)
                    } else {
                        extend_one(first, second, diagonal, source)
                    };
                    if TRACK_MATCHES {
                        matches = matches.saturating_add(
                            usize::try_from(offsets[index] - source)
                                .expect("an extension length is nonnegative"),
                        );
                    }
                    scalar_bits &= scalar_bits - 1;
                }
            }

            for (index, offset) in offsets[vector_length..].iter_mut().enumerate() {
                if *offset < 0 {
                    continue;
                }
                let absolute_index = vector_length + index;
                let diagonal =
                    lower + i32::try_from(absolute_index).expect("validated interval fits i32");
                let source = *offset;
                *offset = extend_one(first, second, diagonal, source);
                if TRACK_MATCHES {
                    matches = matches.saturating_add(
                        usize::try_from(*offset - source)
                            .expect("an extension length is nonnegative"),
                    );
                }
                active_lower = active_lower.min(diagonal);
                active_upper = diagonal;
            }

            if TRACK_MATCHES {
                let mut vector_match_counts = [0_u64; 4];
                _mm256_storeu_si256(
                    vector_match_counts.as_mut_ptr().cast::<__m256i>(),
                    vector_matches,
                );
                matches = vector_match_counts.iter().fold(matches, |total, count| {
                    total.saturating_add(usize::try_from(*count).unwrap_or(usize::MAX))
                });
            }

            WavefrontExtensionResult {
                range: (active_lower <= active_upper).then_some(DiagonalRange {
                    lower: active_lower,
                    upper: active_upper,
                }),
                matches,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{WavefrontExtension, extend_portable};

    fn compare(first: &[u8], second: &[u8], lower: i32, offsets: &[i32]) {
        let mut expected = offsets.to_vec();
        let mut observed = offsets.to_vec();
        let expected_range = extend_portable(first, second, lower, &mut expected);
        let expected_matches = expected
            .iter()
            .zip(offsets)
            .fold(0_usize, |total, (end, start)| {
                total.saturating_add(usize::try_from(end - start).unwrap_or(0))
            });
        let observed_result =
            WavefrontExtension::new().extend_with_matches(first, second, lower, &mut observed);
        assert_eq!(observed, expected);
        assert_eq!(observed_result.range(), expected_range);
        assert_eq!(observed_result.matches(), expected_matches);
    }

    #[test]
    fn extends_mixed_diagonals_and_boundaries() {
        let first = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let second = b"ACGTACGTTCGTACGTACGTTCGTACGTACGA";
        compare(
            first,
            second,
            -8,
            &[-1, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, -1],
        );
        compare(first, second, 0, &[28, 29, 30, 31, 32, -1, -1, -1]);
        compare(first, second, 0, &[]);
    }

    #[test]
    fn agrees_for_all_small_legal_coordinates() {
        let first = b"ACGTACGTTGCA";
        let second = b"ACCTACGATGCA";
        for lower in -12_i32..=12 {
            let mut offsets = Vec::new();
            for lane in 0_i32..24 {
                let diagonal = lower + lane;
                let minimum = diagonal.max(0);
                let maximum = i32::try_from(second.len())
                    .expect("length fits i32")
                    .min(i32::try_from(first.len()).expect("length fits i32") + diagonal);
                offsets.push(if minimum <= maximum {
                    minimum + (lane % (maximum - minimum + 1))
                } else {
                    -1
                });
            }
            compare(first, second, lower, &offsets);
        }
    }

    #[test]
    fn long_extensions_cross_vector_and_tail_boundaries() {
        let first = vec![b'A'; 8_257];
        let boundaries = [0, 1, 30, 31, 32, 33, 62, 63, 64, 65, 8_191, 8_192, 8_256];

        for start in [0, 5] {
            for mismatch in boundaries {
                let mut second = first.clone();
                second[mismatch] = b'C';
                compare(&first, &second, -7, &[start; 8]);
            }
        }

        compare(&first, &first, -7, &[0; 8]);
        compare(&first, &first, -7, &[5; 8]);
    }

    #[test]
    fn invalid_coordinates_remain_memory_safe_and_unchanged() {
        compare(b"ACGT", b"ACGT", -4, &[99, -1, 0, 1, 2, 3, 4, 99]);
    }
}

#[cfg(test)]
mod explicit_portability_tests {
    use core::fmt::Write as _;

    use super::{DiagonalRange, ExtensionKind, WavefrontExtension, WavefrontExtensionResult};

    // Independent byte-at-a-time oracle: no word loads, vector intrinsics or
    // calls to production extension helpers. Wider coordinates avoid reproducing
    // i32 arithmetic overflow in the expected-value calculation.
    fn scalar_expected(
        first: &[u8],
        second: &[u8],
        lower: i32,
        offsets: &mut [i32],
    ) -> WavefrontExtensionResult {
        let mut active = None;
        let mut matches = 0;
        for (lane, offset) in offsets.iter_mut().enumerate() {
            if *offset < 0 {
                continue;
            }
            let diagonal = i64::from(lower) + i64::try_from(lane).unwrap();
            let Ok(diagonal_i32) = i32::try_from(diagonal) else {
                continue;
            };
            active = Some(active.map_or(
                DiagonalRange {
                    lower: diagonal_i32,
                    upper: diagonal_i32,
                },
                |range: DiagonalRange| DiagonalRange {
                    lower: range.lower,
                    upper: diagonal_i32,
                },
            ));
            let mut horizontal = i64::from(*offset);
            let mut vertical = horizontal - diagonal;
            while let (Ok(a), Ok(b)) = (usize::try_from(vertical), usize::try_from(horizontal)) {
                let (Some(a), Some(b)) = (first.get(a), second.get(b)) else {
                    break;
                };
                if a != b {
                    break;
                }
                vertical += 1;
                horizontal += 1;
                matches += 1;
            }
            *offset = i32::try_from(horizontal).unwrap();
        }
        WavefrontExtensionResult {
            range: active,
            matches,
        }
    }

    fn compare(
        kernel: WavefrontExtension,
        first: &[u8],
        second: &[u8],
        lower: i32,
        offsets: &[i32],
    ) {
        let mut expected = offsets.to_vec();
        let expected_result = scalar_expected(first, second, lower, &mut expected);
        let mut actual = offsets.to_vec();
        let observed = kernel.extend_with_matches(first, second, lower, &mut actual);
        assert_eq!(
            actual, expected,
            "backend={:?}, lower={lower}, offsets={offsets:?}",
            kernel.kind
        );
        assert_eq!(observed, expected_result);
        if kernel.is_vectorized() {
            let mut portable_offsets = offsets.to_vec();
            let portable_result = WavefrontExtension {
                kind: ExtensionKind::Portable,
            }
            .extend_with_matches(first, second, lower, &mut portable_offsets);
            assert_eq!(actual, portable_offsets, "explicit AVX2/Portable agreement");
            assert_eq!(observed, portable_result);
        }
        actual.copy_from_slice(offsets);
        assert_eq!(
            kernel.extend(first, second, lower, &mut actual),
            expected_result.range()
        );
        assert_eq!(actual, expected);
    }

    fn next(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        *state
    }

    fn boundaries(kernel: WavefrontExtension) -> usize {
        let lengths = [0, 1, 2, 3, 4, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65];
        let diagonals = [
            i32::MIN,
            i32::MIN + 1,
            -65,
            -33,
            -8,
            -1,
            0,
            1,
            8,
            33,
            65,
            i32::MAX - 8,
            i32::MAX - 1,
            i32::MAX,
        ];
        let mut cases = 0;
        for length in lengths {
            let first = vec![0x80; length];
            let mut second = first.clone();
            if let Some(last) = second.last_mut() {
                *last ^= 1;
            }
            let length_i32 = i32::try_from(length).unwrap();
            let values = [
                i32::MIN,
                -1,
                0,
                1,
                2,
                length_i32.saturating_sub(1),
                length_i32,
                length_i32 + 1,
                i32::MAX,
            ];
            for count in [0, 1, 7, 8, 9, 15, 16, 17] {
                let offsets: Vec<_> = (0..count).map(|lane| values[lane % values.len()]).collect();
                for lower in diagonals {
                    compare(kernel, &first, &second, lower, &offsets);
                    cases += 1;
                }
            }
        }
        cases
    }

    fn long_unaligned(kernel: WavefrontExtension) -> usize {
        let mut cases = 0;
        for length in [255, 256, 257, 8_191, 8_192, 8_193, 8_257, 16_385] {
            for first_shift in [0, 1, 7] {
                for second_shift in [0, 3] {
                    for mismatch in [
                        None,
                        Some(0),
                        Some(3),
                        Some(7),
                        Some(8),
                        Some(31),
                        Some(32),
                        Some(length / 2),
                        Some(length - 1),
                    ] {
                        let first = vec![0xff; first_shift + length];
                        let mut second = vec![0xff; second_shift + length];
                        if let Some(mismatch) = mismatch {
                            second[second_shift + mismatch] = 0;
                        }
                        compare(
                            kernel,
                            &first[first_shift..],
                            &second[second_shift..],
                            -8,
                            &[0; 17],
                        );
                        cases += 1;
                    }
                }
            }
        }
        cases
    }

    fn seeded_arbitrary_bytes(kernel: WavefrontExtension) -> usize {
        let mut state = 0xdaf1_a55e_5eed_0123;
        for case in 0..2_048 {
            let first_len = usize::try_from(next(&mut state) % 257).unwrap();
            let second_len = usize::try_from(next(&mut state) % 257).unwrap();
            let first: Vec<_> = (0..first_len)
                .map(|_| next(&mut state).to_le_bytes()[4])
                .collect();
            let mut second: Vec<_> = (0..second_len)
                .map(|_| next(&mut state).to_le_bytes()[4])
                .collect();
            if case % 2 == 0 {
                let shared = first_len.min(second_len);
                second[..shared].copy_from_slice(&first[..shared]);
            }
            let count = usize::try_from(next(&mut state) % 41).unwrap();
            let offsets: Vec<_> = (0..count)
                .map(|lane| {
                    if lane % 7 == 0 {
                        -1
                    } else {
                        i32::try_from(next(&mut state) % 260).unwrap()
                    }
                })
                .collect();
            let lower = i32::try_from(next(&mut state) % 81).unwrap() - 40;
            compare(kernel, &first, &second, lower, &offsets);
        }
        2_048
    }

    fn word_mismatches(kernel: WavefrontExtension) -> usize {
        for position in 0..8 {
            for difference in 1..=255_u8 {
                let first = [0xa5; 24];
                let mut second = first;
                second[position] ^= difference;
                compare(kernel, &first, &second, 0, &[0, 1, 2, 3, 4, 5, 6, 7]);
            }
        }
        2_040
    }

    #[test]
    fn explicit_portable_boundaries_and_extreme_offsets_match_scalar_oracle() {
        let portable = WavefrontExtension {
            kind: ExtensionKind::Portable,
        };
        assert!(!portable.is_vectorized());
        assert_eq!(boundaries(portable), 1_904);
        println!("Portable boundary lane: EXECUTED (1904 cases)");
    }

    #[test]
    fn explicit_portable_long_unaligned_extensions_match_scalar_oracle() {
        assert_eq!(
            long_unaligned(WavefrontExtension {
                kind: ExtensionKind::Portable
            }),
            432
        );
    }

    #[test]
    fn explicit_portable_seeded_arbitrary_bytes_match_scalar_oracle() {
        assert_eq!(
            seeded_arbitrary_bytes(WavefrontExtension {
                kind: ExtensionKind::Portable
            }),
            2_048
        );
    }

    #[test]
    fn explicit_portable_every_word_byte_mismatch_matches_scalar_oracle() {
        assert_eq!(
            word_mismatches(WavefrontExtension {
                kind: ExtensionKind::Portable
            }),
            2_040
        );
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn explicit_avx2_lane_matches_scalar_oracle_after_runtime_feature_guard() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            println!("AVX2 lane: NOT EXECUTED (CPU/OS capability unavailable)");
            return;
        }
        // The private kind is constructed only after the same runtime guard as
        // production dispatch. No new unsafe API or unchecked intrinsic is used.
        let avx2 = WavefrontExtension {
            kind: ExtensionKind::Avx2,
        };
        assert!(avx2.is_vectorized());
        let cases = boundaries(avx2)
            + long_unaligned(avx2)
            + seeded_arbitrary_bytes(avx2)
            + word_mismatches(avx2);
        assert_eq!(cases, 6_424);
        println!("AVX2 lane: EXECUTED (6424 cases)");
    }

    #[test]
    fn documented_safe_extension_example_matches_complete_output() {
        let kernel = WavefrontExtension::new();
        let mut offsets = [0];
        let result = kernel.extend_with_matches(b"ACGTACGT", b"ACGTTTGT", 0, &mut offsets);
        let mut observed = format!(
            "Offsets: {offsets:?}\nDiagonal range: {:?}\nMatched bytes: {}\n",
            result.range().map(|range| (range.lower(), range.upper())),
            result.matches()
        );
        let mut invalid = [-1, 99];
        let result = kernel.extend_with_matches(b"ACGTACGT", b"ACGTTTGT", 0, &mut invalid);
        assert_eq!(result.matches(), 0);
        writeln!(observed, "Invalid offsets retained: {invalid:?}").unwrap();
        let guide = include_str!("../../docs/testing/portability.md");
        let expected = guide
            .split_once("<!-- portability-rust-output -->\n```text\n")
            .unwrap()
            .1
            .split_once("\n```")
            .unwrap()
            .0;
        assert_eq!(observed, format!("{expected}\n"));
    }
}
