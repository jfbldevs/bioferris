use core::mem::size_of;

use align_core::{AlignmentOperation, Alphabet, PairwiseAlignment, Sequence, Symbol};

use crate::MyersError;

const WORD_BITS: usize = u64::BITS as usize;
const ADAPTIVE_MIN_BLOCKS: usize = 8;
const ADAPTIVE_MIN_DISTANCE: usize = WORD_BITS / 2;

#[derive(Clone, Copy, Debug)]
struct BlockState {
    pv: u64,
    mv: u64,
}

#[derive(Clone, Copy, Debug, Default)]
struct BlockHistory {
    diagonal: u64,
    deletion: u64,
}

#[derive(Clone, Copy, Debug)]
struct HistoryColumn {
    first_block: usize,
    history_end: usize,
}

/// Reusable working memory for repeated Myers executions.
///
/// The workspace retains frontier and traceback buffers between calls. It never
/// retains references to a pattern or text, so one workspace can be reused with
/// different [`MyersPattern`] values and sequence lengths.
#[derive(Clone, Debug, Default)]
pub struct MyersWorkspace {
    states: Vec<BlockState>,
    block_scores: Vec<usize>,
    history: Vec<BlockHistory>,
    history_columns: Vec<HistoryColumn>,
    history_scratch: Vec<BlockHistory>,
}

impl MyersWorkspace {
    /// Creates an empty workspace that allocates on first use.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            states: Vec::new(),
            block_scores: Vec::new(),
            history: Vec::new(),
            history_columns: Vec::new(),
            history_scratch: Vec::new(),
        }
    }

    /// Releases all retained working buffers.
    pub fn release(&mut self) {
        self.states = Vec::new();
        self.block_scores = Vec::new();
        self.history = Vec::new();
        self.history_columns = Vec::new();
        self.history_scratch = Vec::new();
    }

    /// Returns the bytes retained by element buffers at their current capacities.
    #[must_use]
    pub fn retained_bytes(&self) -> usize {
        self.states.capacity() * size_of::<BlockState>()
            + self.block_scores.capacity() * size_of::<usize>()
            + self.history.capacity() * size_of::<BlockHistory>()
            + self.history_columns.capacity() * size_of::<HistoryColumn>()
            + self.history_scratch.capacity() * size_of::<BlockHistory>()
    }
}

/// One approximate-search result from Myers's free-text-prefix boundary.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MyersMatch {
    end: usize,
    distance: usize,
}

impl MyersMatch {
    /// Returns the exclusive text endpoint, in `1..=text.len()`.
    #[must_use]
    pub const fn end(self) -> usize {
        self.end
    }

    /// Returns the exact unit edit distance at this endpoint.
    #[must_use]
    pub const fn distance(self) -> usize {
        self.distance
    }
}

/// Modeled heap element buffers used by the Myers profile.
///
/// Allocator metadata, `Vec` headers, the returned alignment object, and stack
/// values are excluded. Element buffers are included at their Rust element sizes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MyersMemoryEstimate {
    pattern_symbols: usize,
    peq_words: usize,
    state_words: usize,
    band_score_records: usize,
    history_records: usize,
    history_column_records: usize,
    history_scratch_records: usize,
    operation_capacity: usize,
    search_match_capacity: usize,
}

impl MyersMemoryEstimate {
    /// Number of copied validated pattern symbols.
    #[must_use]
    pub const fn pattern_symbols(self) -> usize {
        self.pattern_symbols
    }

    /// Number of precomputed `Peq` words.
    #[must_use]
    pub const fn peq_words(self) -> usize {
        self.peq_words
    }

    /// Number of mutable `Pv`/`Mv` words needed by distance or search.
    #[must_use]
    pub const fn state_words(self) -> usize {
        self.state_words
    }

    /// Number of bottom-score records used by automatically banded distance.
    #[must_use]
    pub const fn band_score_records(self) -> usize {
        self.band_score_records
    }

    /// Conservative maximum number of two-mask records retained for traceback.
    #[must_use]
    pub const fn history_records(self) -> usize {
        self.history_records
    }

    /// Number of per-column records used by long-pattern banded traceback.
    #[must_use]
    pub const fn history_column_records(self) -> usize {
        self.history_column_records
    }

    /// Number of temporary decision-mask records used by banded traceback.
    #[must_use]
    pub const fn history_scratch_records(self) -> usize {
        self.history_scratch_records
    }

    /// Maximum operation elements reserved by global traceback.
    #[must_use]
    pub const fn operation_capacity(self) -> usize {
        self.operation_capacity
    }

    /// Maximum search-result elements when every text endpoint matches.
    #[must_use]
    pub const fn search_match_capacity(self) -> usize {
        self.search_match_capacity
    }

    /// Persistent preprocessed-pattern element bytes.
    #[must_use]
    pub const fn persistent_bytes(self) -> usize {
        self.pattern_symbols * size_of::<u8>() + self.peq_words * size_of::<u64>()
    }

    /// Peak modeled bytes for distance-only execution, including preprocessing.
    #[must_use]
    pub const fn distance_peak_bytes(self) -> usize {
        self.persistent_bytes()
            + self.state_words * size_of::<u64>()
            + self.band_score_records * size_of::<usize>()
    }

    /// Peak modeled bytes for global alignment, including preprocessing and output.
    #[must_use]
    pub const fn alignment_peak_bytes(self) -> usize {
        self.distance_peak_bytes()
            + self.history_records * size_of::<BlockHistory>()
            + self.history_column_records * size_of::<HistoryColumn>()
            + self.history_scratch_records * size_of::<BlockHistory>()
            + self.operation_capacity * size_of::<AlignmentOperation>()
    }

    /// Peak modeled bytes for search when every endpoint is reported.
    #[must_use]
    pub const fn search_peak_bytes(self) -> usize {
        self.distance_peak_bytes() + self.search_match_capacity * size_of::<MyersMatch>()
    }
}

/// Reusable Myers bit-vector pattern for literal unit edit distance.
///
/// `A` may be DNA, RNA, protein, or another finite [`Alphabet`]. Pattern symbols
/// are compared literally; ambiguity codes are not interpreted as wildcards.
#[derive(Clone, Debug)]
pub struct MyersPattern<A: Alphabet> {
    pattern: Vec<Symbol<A>>,
    peq: Vec<u64>,
    blocks: usize,
    last_mask: u64,
    last_high_bit: u64,
}

impl<A: Alphabet> MyersPattern<A> {
    /// Preprocesses a validated pattern into symbol-equality bit masks.
    ///
    /// # Errors
    ///
    /// Returns [`MyersError`] if derived dimensions overflow or memory cannot be
    /// reserved.
    pub fn new(pattern: &Sequence<A>) -> Result<Self, MyersError> {
        let estimate = Self::memory_requirement(pattern.len(), 0)?;
        let blocks = block_count(pattern.len());

        let mut copied_pattern = Vec::new();
        copied_pattern
            .try_reserve_exact(estimate.pattern_symbols)
            .map_err(|_| MyersError::AllocationFailed {
                buffer: "Myers pattern symbols",
                elements: estimate.pattern_symbols,
            })?;
        copied_pattern.extend_from_slice(pattern.as_slice());

        let mut peq = Vec::new();
        peq.try_reserve_exact(estimate.peq_words)
            .map_err(|_| MyersError::AllocationFailed {
                buffer: "Myers Peq table",
                elements: estimate.peq_words,
            })?;
        peq.resize(estimate.peq_words, 0);

        for (row, symbol) in copied_pattern.iter().copied().enumerate() {
            let block = row / WORD_BITS;
            let bit = row % WORD_BITS;
            let index = symbol.alphabet_index() * blocks + block;
            peq[index] |= 1_u64 << bit;
        }

        let last_width = pattern.len() % WORD_BITS;
        let last_width = if last_width == 0 && !pattern.is_empty() {
            WORD_BITS
        } else {
            last_width
        };
        let (last_mask, last_high_bit) = if last_width == 0 {
            (0, 0)
        } else if last_width == WORD_BITS {
            (u64::MAX, 1_u64 << (WORD_BITS - 1))
        } else {
            ((1_u64 << last_width) - 1, 1_u64 << (last_width - 1))
        };

        Ok(Self {
            pattern: copied_pattern,
            peq,
            blocks,
            last_mask,
            last_high_bit,
        })
    }

    /// Returns the number of symbols in the preprocessed pattern.
    #[must_use]
    pub fn pattern_len(&self) -> usize {
        self.pattern.len()
    }

    /// Returns modeled element-buffer memory for the supplied dimensions.
    ///
    /// # Errors
    ///
    /// Returns [`MyersError::DimensionOverflow`] if a derived element count or byte
    /// total cannot be represented by `usize`.
    pub fn memory_requirement(
        pattern_len: usize,
        text_len: usize,
    ) -> Result<MyersMemoryEstimate, MyersError> {
        let blocks = block_count(pattern_len);
        let peq_words =
            A::SYMBOLS
                .len()
                .checked_mul(blocks)
                .ok_or(MyersError::DimensionOverflow {
                    pattern_len,
                    text_len,
                })?;
        let state_words = blocks.checked_mul(2).ok_or(MyersError::DimensionOverflow {
            pattern_len,
            text_len,
        })?;
        let band_score_records = if blocks >= ADAPTIVE_MIN_BLOCKS {
            blocks
        } else {
            0
        };
        let history_records =
            blocks
                .checked_mul(text_len)
                .ok_or(MyersError::DimensionOverflow {
                    pattern_len,
                    text_len,
                })?;
        let (history_column_records, history_scratch_records) =
            if blocks >= ADAPTIVE_MIN_BLOCKS && text_len != 0 {
                (text_len, blocks)
            } else {
                (0, 0)
            };
        let operation_capacity =
            pattern_len
                .checked_add(text_len)
                .ok_or(MyersError::DimensionOverflow {
                    pattern_len,
                    text_len,
                })?;

        let estimate = MyersMemoryEstimate {
            pattern_symbols: pattern_len,
            peq_words,
            state_words,
            band_score_records,
            history_records,
            history_column_records,
            history_scratch_records,
            operation_capacity,
            search_match_capacity: text_len,
        };

        let dimension_error = || MyersError::DimensionOverflow {
            pattern_len,
            text_len,
        };
        let persistent_bytes = pattern_len
            .checked_mul(size_of::<u8>())
            .and_then(|symbols| {
                peq_words
                    .checked_mul(size_of::<u64>())
                    .and_then(|peq_bytes| symbols.checked_add(peq_bytes))
            })
            .ok_or_else(dimension_error)?;
        let distance_bytes = state_words
            .checked_mul(size_of::<u64>())
            .and_then(|state_bytes| persistent_bytes.checked_add(state_bytes))
            .and_then(|bytes| {
                band_score_records
                    .checked_mul(size_of::<usize>())
                    .and_then(|score_bytes| bytes.checked_add(score_bytes))
            })
            .ok_or_else(dimension_error)?;
        history_records
            .checked_mul(size_of::<BlockHistory>())
            .and_then(|history_bytes| distance_bytes.checked_add(history_bytes))
            .and_then(|bytes| {
                history_column_records
                    .checked_mul(size_of::<HistoryColumn>())
                    .and_then(|column_bytes| bytes.checked_add(column_bytes))
            })
            .and_then(|bytes| {
                history_scratch_records
                    .checked_mul(size_of::<BlockHistory>())
                    .and_then(|scratch_bytes| bytes.checked_add(scratch_bytes))
            })
            .and_then(|bytes| {
                operation_capacity
                    .checked_mul(size_of::<AlignmentOperation>())
                    .and_then(|operation_bytes| bytes.checked_add(operation_bytes))
            })
            .ok_or_else(dimension_error)?;
        text_len
            .checked_mul(size_of::<MyersMatch>())
            .and_then(|match_bytes| distance_bytes.checked_add(match_bytes))
            .ok_or_else(dimension_error)?;

        Ok(estimate)
    }

    /// Computes exact global unit edit distance using bit-vector frontier state.
    ///
    /// Multiword patterns use an automatically widened exact band once the full
    /// frontier would span at least eight machine words.
    ///
    /// # Errors
    ///
    /// Returns [`MyersError`] if working memory cannot be reserved or a state
    /// invariant is violated.
    pub fn global_distance(&self, text: &Sequence<A>) -> Result<usize, MyersError> {
        let mut workspace = MyersWorkspace::new();
        self.global_distance_with_workspace(text, &mut workspace)
    }

    /// Computes exact global unit edit distance with reusable working memory.
    ///
    /// This method is equivalent to [`Self::global_distance`]. Reusing the same
    /// workspace avoids frontier allocation for multiword patterns.
    ///
    /// # Errors
    ///
    /// Returns [`MyersError`] if working memory cannot be reserved or a state
    /// invariant is violated.
    pub fn global_distance_with_workspace(
        &self,
        text: &Sequence<A>,
        workspace: &mut MyersWorkspace,
    ) -> Result<usize, MyersError> {
        if self.blocks >= ADAPTIVE_MIN_BLOCKS {
            return self.global_distance_adaptive_with_workspace(text, workspace);
        }
        self.global_distance_full_with_workspace(text, workspace)
    }

    fn global_distance_full_with_workspace(
        &self,
        text: &Sequence<A>,
        workspace: &mut MyersWorkspace,
    ) -> Result<usize, MyersError> {
        if self.pattern.is_empty() {
            return Ok(text.len());
        }
        if self.blocks == 1 {
            return self.global_distance_one_block(text);
        }

        self.initialize_states(&mut workspace.states)?;
        let mut score = self.pattern.len();
        for (column, symbol) in text.iter().enumerate() {
            let horizontal = self.advance_score_column(&mut workspace.states, symbol, 1);
            score = update_score(score, horizontal)
                .ok_or(MyersError::StateInvariant { column: column + 1 })?;
        }
        Ok(score)
    }

    fn global_distance_adaptive_with_workspace(
        &self,
        text: &Sequence<A>,
        workspace: &mut MyersWorkspace,
    ) -> Result<usize, MyersError> {
        let maximum = self.pattern.len().max(text.len());
        let mut max_distance = ADAPTIVE_MIN_DISTANCE
            .max(maximum / 18)
            .max(self.pattern.len().abs_diff(text.len()))
            .min(maximum);
        loop {
            if let Some(distance) =
                self.global_distance_bounded_with_workspace(text, max_distance, workspace)?
            {
                return Ok(distance);
            }
            if max_distance == maximum {
                return Err(MyersError::StateInvariant { column: text.len() });
            }
            max_distance = max_distance
                .saturating_mul(2)
                .checked_next_power_of_two()
                .unwrap_or(maximum)
                .min(maximum);
        }
    }

    /// Computes global unit edit distance when it does not exceed `max_distance`.
    ///
    /// The bounded frontier avoids visiting blocks that cannot be part of an
    /// alignment within the requested limit. [`None`] means that the exact
    /// distance is greater than `max_distance`.
    ///
    /// # Errors
    ///
    /// Returns [`MyersError`] if working memory cannot be reserved or a state
    /// invariant is violated.
    pub fn global_distance_bounded(
        &self,
        text: &Sequence<A>,
        max_distance: usize,
    ) -> Result<Option<usize>, MyersError> {
        let mut workspace = MyersWorkspace::new();
        self.global_distance_bounded_with_workspace(text, max_distance, &mut workspace)
    }

    /// Computes bounded global unit edit distance with reusable working memory.
    ///
    /// This method is equivalent to [`Self::global_distance_bounded`]. Reusing
    /// the workspace avoids reallocating frontier and block-score buffers.
    ///
    /// # Errors
    ///
    /// Returns [`MyersError`] if working memory cannot be reserved or a state
    /// invariant is violated.
    pub fn global_distance_bounded_with_workspace(
        &self,
        text: &Sequence<A>,
        max_distance: usize,
        workspace: &mut MyersWorkspace,
    ) -> Result<Option<usize>, MyersError> {
        self.global_distance_bounded_kernel::<false>(text, max_distance, workspace)
    }

    #[allow(clippy::too_many_lines)]
    fn global_distance_bounded_kernel<const TRACE: bool>(
        &self,
        text: &Sequence<A>,
        max_distance: usize,
        workspace: &mut MyersWorkspace,
    ) -> Result<Option<usize>, MyersError> {
        let max_distance = max_distance.min(self.pattern.len().max(text.len()));
        let length_difference = self.pattern.len().abs_diff(text.len());
        if length_difference > max_distance {
            return Ok(None);
        }
        if self.pattern.is_empty() {
            return Ok(Some(text.len()));
        }
        if text.is_empty() {
            return Ok(Some(self.pattern.len()));
        }
        if max_distance == 0 {
            return Ok((self.pattern.as_slice() == text.as_slice()).then_some(0));
        }
        if self.blocks == 1 {
            let distance = self.global_distance_one_block(text)?;
            return Ok((distance <= max_distance).then_some(distance));
        }

        self.initialize_states(&mut workspace.states)?;
        initialize_block_scores(&mut workspace.block_scores, self.blocks, self.pattern.len())?;
        if TRACE {
            initialize_banded_history(
                workspace,
                self.blocks,
                self.pattern.len(),
                text.len(),
                max_distance,
            )?;
        }

        let initial_numerator = max_distance
            .checked_add(2)
            .and_then(|value| {
                if self.pattern.len() >= text.len() {
                    value.checked_add(self.pattern.len() - text.len())
                } else {
                    value.checked_sub(text.len() - self.pattern.len())
                }
            })
            .ok_or(MyersError::DimensionOverflow {
                pattern_len: self.pattern.len(),
                text_len: text.len(),
            })?;
        let initial_rows = max_distance.min(initial_numerator / 2).max(1);
        let first_block = 0_usize;
        let last_block = ((initial_rows - 1) / WORD_BITS).min(self.blocks - 1);
        if maximum_band_dimension(self.pattern.len(), text.len()) {
            self.global_distance_bounded_columns::<TRACE, true>(
                text,
                max_distance,
                workspace,
                first_block,
                last_block,
            )
        } else {
            self.global_distance_bounded_columns::<TRACE, false>(
                text,
                max_distance,
                workspace,
                first_block,
                last_block,
            )
        }
    }

    #[allow(clippy::too_many_lines)]
    fn global_distance_bounded_columns<const TRACE: bool, const FAST_ARITHMETIC: bool>(
        &self,
        text: &Sequence<A>,
        max_distance: usize,
        workspace: &mut MyersWorkspace,
        mut first_block: usize,
        mut last_block: usize,
    ) -> Result<Option<usize>, MyersError> {
        let band = GlobalBand {
            pattern_len: self.pattern.len(),
            text_len: text.len(),
            max_distance,
            fast_arithmetic: FAST_ARITHMETIC,
        };
        let all_blocks_full = self.last_mask == u64::MAX;

        for (column_offset, symbol) in text.iter().enumerate() {
            let column = column_offset + 1;
            let peq_start = symbol.alphabet_index() * self.blocks;
            let equalities = &self.peq[peq_start..peq_start + self.blocks];
            let active_end = last_block + 1;
            let active_states = &mut workspace.states[first_block..active_end];
            let active_scores = &mut workspace.block_scores[first_block..active_end];
            let active_equalities = &equalities[first_block..active_end];
            let horizontal = if all_blocks_full {
                advance_banded_blocks::<TRACE, true>(
                    active_states,
                    active_scores,
                    active_equalities,
                    &mut workspace.history_scratch,
                    first_block,
                    self.blocks,
                    self.last_mask,
                    self.last_high_bit,
                )
            } else {
                advance_banded_blocks::<TRACE, false>(
                    active_states,
                    active_scores,
                    active_equalities,
                    &mut workspace.history_scratch,
                    first_block,
                    self.blocks,
                    self.last_mask,
                    self.last_high_bit,
                )
            };

            if last_block + 1 < self.blocks {
                let previous_score = if last_block == 0 {
                    column
                } else {
                    workspace.block_scores[last_block - 1]
                };
                let next_block = last_block + 1;
                if !block_is_below_global_band(
                    workspace.block_scores[last_block],
                    previous_score,
                    block_width(last_block, self.pattern.len()),
                    block_last_row(last_block, self.pattern.len()),
                    column,
                    band,
                ) {
                    let (mask, high_bit) = self.block_mask_and_high_bit(next_block);
                    let horizontal_input = horizontal;
                    let (_, ph, mh) = advance_block(
                        &mut workspace.states[next_block],
                        equalities[next_block],
                        horizontal_input,
                        mask,
                        high_bit,
                    );
                    if TRACE {
                        workspace.history_scratch[next_block] = trace_decisions(
                            workspace.states[next_block],
                            equalities[next_block],
                            horizontal_input,
                            ph,
                            mh,
                            mask,
                        );
                    }
                    workspace.block_scores[next_block] = add_vertical_deltas(
                        workspace.block_scores[last_block],
                        workspace.states[next_block],
                        mask,
                        column,
                    )?;
                    last_block = next_block;
                }
            }

            while first_block <= last_block
                && block_is_above_global_band(
                    workspace.block_scores[first_block],
                    block_width(first_block, self.pattern.len()),
                    block_last_row(first_block, self.pattern.len()),
                    column,
                    band,
                )
            {
                first_block += 1;
            }
            if first_block > last_block {
                return Ok(None);
            }

            while last_block > first_block
                && block_is_below_global_band(
                    workspace.block_scores[last_block],
                    workspace.block_scores[last_block - 1],
                    block_width(last_block, self.pattern.len()),
                    block_last_row(last_block, self.pattern.len()),
                    column,
                    band,
                )
            {
                last_block -= 1;
            }
            if TRACE {
                append_banded_history(workspace, first_block, last_block)?;
            }
        }

        if last_block + 1 == self.blocks && first_block <= last_block {
            let distance = workspace.block_scores[last_block];
            Ok((distance <= max_distance).then_some(distance))
        } else {
            Ok(None)
        }
    }

    /// Reports every approximate-search endpoint within `max_distance`.
    ///
    /// Text prefixes are free, as in the problem stated by Myers. Returned endpoints
    /// are exclusive and ordered increasingly. The artificial column zero is not
    /// reported.
    ///
    /// # Errors
    ///
    /// Returns [`MyersError`] if working or result memory cannot be reserved or a
    /// state invariant is violated.
    pub fn find_all(
        &self,
        text: &Sequence<A>,
        max_distance: usize,
    ) -> Result<Vec<MyersMatch>, MyersError> {
        let mut workspace = MyersWorkspace::new();
        self.find_all_with_workspace(text, max_distance, &mut workspace)
    }

    /// Reports approximate-search endpoints using reusable frontier storage.
    ///
    /// The returned vector remains independently owned; only internal frontier
    /// allocation is retained by `workspace`.
    ///
    /// # Errors
    ///
    /// Returns [`MyersError`] if working or result memory cannot be reserved or a
    /// state invariant is violated.
    pub fn find_all_with_workspace(
        &self,
        text: &Sequence<A>,
        max_distance: usize,
        workspace: &mut MyersWorkspace,
    ) -> Result<Vec<MyersMatch>, MyersError> {
        let mut matches = Vec::new();
        if self.pattern.is_empty() {
            matches
                .try_reserve_exact(text.len())
                .map_err(|_| MyersError::AllocationFailed {
                    buffer: "Myers search matches",
                    elements: text.len(),
                })?;
            matches.extend((1..=text.len()).map(|end| MyersMatch { end, distance: 0 }));
            return Ok(matches);
        }
        if self.blocks == 1 {
            return self.find_all_one_block(text, max_distance);
        }

        self.initialize_states(&mut workspace.states)?;
        let mut score = self.pattern.len();
        for (column, symbol) in text.iter().enumerate() {
            let horizontal = self.advance_score_column(&mut workspace.states, symbol, 0);
            score = update_score(score, horizontal)
                .ok_or(MyersError::StateInvariant { column: column + 1 })?;
            if score <= max_distance {
                matches
                    .try_reserve(1)
                    .map_err(|_| MyersError::AllocationFailed {
                        buffer: "Myers search matches",
                        elements: matches.len().saturating_add(1),
                    })?;
                matches.push(MyersMatch {
                    end: column + 1,
                    distance: score,
                });
            }
        }
        Ok(matches)
    }

    /// Computes one deterministic optimal global unit-edit alignment.
    ///
    /// The score is the negated edit distance. Traceback prefers diagonal, then
    /// deletion, then insertion when several predecessors are optimal.
    ///
    /// # Errors
    ///
    /// Returns [`MyersError`] for dimension or allocation failure, an invalid
    /// bit-vector history, score conversion overflow, or invalid output structure.
    pub fn align_global(&self, text: &Sequence<A>) -> Result<PairwiseAlignment, MyersError> {
        let mut workspace = MyersWorkspace::new();
        self.align_global_with_workspace(text, &mut workspace)
    }

    /// Computes an exact global alignment with reusable frontier and history memory.
    ///
    /// The returned alignment remains independently owned. Reusing `workspace`
    /// avoids reallocating the multiword frontier and traceback history.
    ///
    /// # Errors
    ///
    /// Returns [`MyersError`] for dimension or allocation failure, an invalid
    /// bit-vector history, score conversion overflow, or invalid output structure.
    pub fn align_global_with_workspace(
        &self,
        text: &Sequence<A>,
        workspace: &mut MyersWorkspace,
    ) -> Result<PairwiseAlignment, MyersError> {
        if self.blocks >= ADAPTIVE_MIN_BLOCKS && !text.is_empty() {
            return self.align_global_banded_with_workspace(text, workspace);
        }
        self.align_global_full_with_workspace(text, workspace)
    }

    fn align_global_full_with_workspace(
        &self,
        text: &Sequence<A>,
        workspace: &mut MyersWorkspace,
    ) -> Result<PairwiseAlignment, MyersError> {
        let estimate = Self::memory_requirement(self.pattern.len(), text.len())?;
        self.initialize_states(&mut workspace.states)?;
        workspace.history.clear();
        workspace
            .history
            .try_reserve_exact(estimate.history_records)
            .map_err(|_| MyersError::AllocationFailed {
                buffer: "Myers traceback history",
                elements: estimate.history_records,
            })?;

        let mut distance = self.pattern.len();
        for (column, symbol) in text.iter().enumerate() {
            let horizontal =
                self.advance_history_column(&mut workspace.states, symbol, &mut workspace.history);
            distance = update_score(distance, horizontal)
                .ok_or(MyersError::StateInvariant { column: column + 1 })?;
        }

        self.build_alignment_from_history::<false>(
            text,
            distance,
            estimate.operation_capacity,
            workspace,
        )
    }

    fn align_global_banded_with_workspace(
        &self,
        text: &Sequence<A>,
        workspace: &mut MyersWorkspace,
    ) -> Result<PairwiseAlignment, MyersError> {
        let estimate = Self::memory_requirement(self.pattern.len(), text.len())?;
        if self.pattern.as_slice() == text.as_slice() {
            let mut operations = Vec::new();
            operations
                .try_reserve_exact(self.pattern.len())
                .map_err(|_| MyersError::AllocationFailed {
                    buffer: "Myers alignment operations",
                    elements: self.pattern.len(),
                })?;
            operations.resize(self.pattern.len(), AlignmentOperation::Match);
            return PairwiseAlignment::try_new(0, 0..self.pattern.len(), 0..text.len(), operations)
                .map_err(MyersError::InvalidResult);
        }

        let distance = self.global_distance_adaptive_with_workspace(text, workspace)?;
        let bounded = self.global_distance_bounded_kernel::<true>(text, distance, workspace)?;
        if bounded != Some(distance) {
            return Err(MyersError::StateInvariant { column: text.len() });
        }
        self.build_alignment_from_history::<true>(
            text,
            distance,
            estimate.operation_capacity,
            workspace,
        )
    }

    fn build_alignment_from_history<const BANDED: bool>(
        &self,
        text: &Sequence<A>,
        distance: usize,
        operation_capacity: usize,
        workspace: &MyersWorkspace,
    ) -> Result<PairwiseAlignment, MyersError> {
        let mut operations = Vec::new();
        operations
            .try_reserve_exact(operation_capacity)
            .map_err(|_| MyersError::AllocationFailed {
                buffer: "Myers alignment operations",
                elements: operation_capacity,
            })?;

        let mut row = self.pattern.len();
        let mut column = text.len();
        let mut traceback_distance = 0_usize;
        while row > 0 || column > 0 {
            if row == 0 {
                traceback_distance = traceback_distance
                    .checked_add(1)
                    .ok_or(MyersError::TracebackInvariant { row, column })?;
                operations.push(AlignmentOperation::Insertion);
                column -= 1;
                continue;
            }
            if column == 0 {
                traceback_distance = traceback_distance
                    .checked_add(1)
                    .ok_or(MyersError::TracebackInvariant { row, column })?;
                operations.push(AlignmentOperation::Deletion);
                row -= 1;
                continue;
            }

            let record = if BANDED {
                banded_history_at(&workspace.history, &workspace.history_columns, row, column)
            } else {
                history_at(&workspace.history, self.blocks, row, column)
            }
            .ok_or(MyersError::TracebackInvariant { row, column })?;
            let bit = 1_u64 << ((row - 1) % WORD_BITS);
            let equal = self.pattern[row - 1] == text.as_slice()[column - 1];
            if record.diagonal & bit != 0 {
                operations.push(if equal {
                    AlignmentOperation::Match
                } else {
                    traceback_distance = traceback_distance
                        .checked_add(1)
                        .ok_or(MyersError::TracebackInvariant { row, column })?;
                    AlignmentOperation::Substitution
                });
                row -= 1;
                column -= 1;
            } else if record.deletion & bit != 0 {
                traceback_distance = traceback_distance
                    .checked_add(1)
                    .ok_or(MyersError::TracebackInvariant { row, column })?;
                operations.push(AlignmentOperation::Deletion);
                row -= 1;
            } else {
                traceback_distance = traceback_distance
                    .checked_add(1)
                    .ok_or(MyersError::TracebackInvariant { row, column })?;
                operations.push(AlignmentOperation::Insertion);
                column -= 1;
            }
        }

        if traceback_distance != distance {
            return Err(MyersError::TracebackInvariant { row, column });
        }
        operations.reverse();

        let score = i64::try_from(distance)
            .ok()
            .and_then(i64::checked_neg)
            .ok_or(MyersError::ScoreOverflow { distance })?;
        PairwiseAlignment::try_new(score, 0..self.pattern.len(), 0..text.len(), operations)
            .map_err(MyersError::InvalidResult)
    }

    fn initialize_states(&self, states: &mut Vec<BlockState>) -> Result<(), MyersError> {
        states.clear();
        states
            .try_reserve_exact(self.blocks)
            .map_err(|_| MyersError::AllocationFailed {
                buffer: "Myers Pv/Mv state",
                elements: self.blocks,
            })?;
        for block in 0..self.blocks {
            states.push(BlockState {
                pv: if block + 1 == self.blocks {
                    self.last_mask
                } else {
                    u64::MAX
                },
                mv: 0,
            });
        }
        Ok(())
    }

    fn block_mask_and_high_bit(&self, block: usize) -> (u64, u64) {
        if block + 1 == self.blocks {
            (self.last_mask, self.last_high_bit)
        } else {
            (u64::MAX, 1_u64 << (WORD_BITS - 1))
        }
    }

    fn global_distance_one_block(&self, text: &Sequence<A>) -> Result<usize, MyersError> {
        let mut state = BlockState {
            pv: self.last_mask,
            mv: 0,
        };
        let mut score = self.pattern.len();
        for (column, symbol) in text.iter().enumerate() {
            let eq = self.peq[symbol.alphabet_index()];
            let (horizontal, _, _) =
                advance_block(&mut state, eq, 1, self.last_mask, self.last_high_bit);
            score = update_score(score, horizontal)
                .ok_or(MyersError::StateInvariant { column: column + 1 })?;
        }
        Ok(score)
    }

    fn find_all_one_block(
        &self,
        text: &Sequence<A>,
        max_distance: usize,
    ) -> Result<Vec<MyersMatch>, MyersError> {
        let mut matches = Vec::new();
        let mut state = BlockState {
            pv: self.last_mask,
            mv: 0,
        };
        let mut score = self.pattern.len();
        for (column, symbol) in text.iter().enumerate() {
            let eq = self.peq[symbol.alphabet_index()];
            let (horizontal, _, _) =
                advance_block(&mut state, eq, 0, self.last_mask, self.last_high_bit);
            score = update_score(score, horizontal)
                .ok_or(MyersError::StateInvariant { column: column + 1 })?;
            if score <= max_distance {
                matches
                    .try_reserve(1)
                    .map_err(|_| MyersError::AllocationFailed {
                        buffer: "Myers search matches",
                        elements: matches.len().saturating_add(1),
                    })?;
                matches.push(MyersMatch {
                    end: column + 1,
                    distance: score,
                });
            }
        }
        Ok(matches)
    }

    fn advance_score_column(
        &self,
        states: &mut [BlockState],
        symbol: Symbol<A>,
        top_input: i8,
    ) -> i8 {
        debug_assert_eq!(states.len(), self.blocks);
        let peq_start = symbol.alphabet_index() * self.blocks;
        let equalities = &self.peq[peq_start..peq_start + self.blocks];
        let mut horizontal = top_input;
        for (block, (state, &eq)) in states.iter_mut().zip(equalities).enumerate() {
            let (mask, high_bit) = if block + 1 == self.blocks {
                (self.last_mask, self.last_high_bit)
            } else {
                (u64::MAX, 1_u64 << (WORD_BITS - 1))
            };
            let (output, _, _) = advance_block(state, eq, horizontal, mask, high_bit);
            horizontal = output;
        }
        horizontal
    }

    fn advance_history_column(
        &self,
        states: &mut [BlockState],
        symbol: Symbol<A>,
        history: &mut Vec<BlockHistory>,
    ) -> i8 {
        debug_assert_eq!(states.len(), self.blocks);
        if self.blocks == 0 {
            return 1;
        }
        let peq_start = symbol.alphabet_index() * self.blocks;
        let equalities = &self.peq[peq_start..peq_start + self.blocks];
        let mut horizontal = 1;
        for (block, (state, &eq)) in states.iter_mut().zip(equalities).enumerate() {
            let (mask, high_bit) = if block + 1 == self.blocks {
                (self.last_mask, self.last_high_bit)
            } else {
                (u64::MAX, 1_u64 << (WORD_BITS - 1))
            };
            let horizontal_input = horizontal;
            let (output, ph, mh) = advance_block(state, eq, horizontal_input, mask, high_bit);
            horizontal = output;
            history.push(trace_decisions(*state, eq, horizontal_input, ph, mh, mask));
        }
        horizontal
    }
}

fn initialize_banded_history(
    workspace: &mut MyersWorkspace,
    blocks: usize,
    pattern_len: usize,
    text_len: usize,
    max_distance: usize,
) -> Result<(), MyersError> {
    workspace.history.clear();
    workspace.history_columns.clear();
    workspace.history_scratch.clear();

    workspace
        .history_columns
        .try_reserve_exact(text_len)
        .map_err(|_| MyersError::AllocationFailed {
            buffer: "Myers banded history columns",
            elements: text_len,
        })?;
    workspace
        .history_scratch
        .try_reserve_exact(blocks)
        .map_err(|_| MyersError::AllocationFailed {
            buffer: "Myers banded history scratch",
            elements: blocks,
        })?;
    workspace
        .history_scratch
        .resize(blocks, BlockHistory::default());

    let estimated_rows = max_distance
        .saturating_mul(2)
        .saturating_add(WORD_BITS.saturating_mul(2));
    let estimated_blocks = block_count(estimated_rows).max(1).min(blocks);
    let estimated_records =
        estimated_blocks
            .checked_mul(text_len)
            .ok_or(MyersError::DimensionOverflow {
                pattern_len,
                text_len,
            })?;
    workspace
        .history
        .try_reserve_exact(estimated_records)
        .map_err(|_| MyersError::AllocationFailed {
            buffer: "Myers banded traceback history",
            elements: estimated_records,
        })?;
    Ok(())
}

fn append_banded_history(
    workspace: &mut MyersWorkspace,
    first_block: usize,
    last_block: usize,
) -> Result<(), MyersError> {
    let records = last_block - first_block + 1;
    workspace
        .history
        .try_reserve(records)
        .map_err(|_| MyersError::AllocationFailed {
            buffer: "Myers banded traceback history",
            elements: workspace.history.len().saturating_add(records),
        })?;
    workspace
        .history
        .extend_from_slice(&workspace.history_scratch[first_block..=last_block]);
    workspace.history_columns.push(HistoryColumn {
        first_block,
        history_end: workspace.history.len(),
    });
    Ok(())
}

#[inline]
fn trace_decisions(
    state: BlockState,
    equality: u64,
    horizontal_input: i8,
    horizontal_positive: u64,
    horizontal_negative: u64,
    mask: u64,
) -> BlockHistory {
    let horizontal_above_positive =
        ((horizontal_positive << 1) | u64::from(horizontal_input > 0)) & mask;
    let horizontal_above_negative =
        ((horizontal_negative << 1) | u64::from(horizontal_input < 0)) & mask;
    let vertical_zero = !(state.pv | state.mv);
    let horizontal_above_zero = !(horizontal_above_positive | horizontal_above_negative);
    let diagonal_delta_one =
        (state.pv & horizontal_above_zero) | (vertical_zero & horizontal_above_positive);
    let equality = equality & mask;
    let diagonal = equality | (!equality & diagonal_delta_one & mask);
    let deletion = !diagonal & state.pv & mask;
    BlockHistory { diagonal, deletion }
}

fn initialize_block_scores(
    scores: &mut Vec<usize>,
    blocks: usize,
    pattern_len: usize,
) -> Result<(), MyersError> {
    scores.clear();
    scores
        .try_reserve_exact(blocks)
        .map_err(|_| MyersError::AllocationFailed {
            buffer: "Myers band block scores",
            elements: blocks,
        })?;
    scores.extend((0..blocks).map(|block| block_last_row(block, pattern_len)));
    Ok(())
}

fn block_last_row(block: usize, pattern_len: usize) -> usize {
    (block + 1).saturating_mul(WORD_BITS).min(pattern_len)
}

fn block_width(block: usize, pattern_len: usize) -> usize {
    block_last_row(block, pattern_len) - block * WORD_BITS
}

fn add_vertical_deltas(
    top_score: usize,
    state: BlockState,
    mask: u64,
    column: usize,
) -> Result<usize, MyersError> {
    let positive = (state.pv & mask).count_ones() as usize;
    let negative = (state.mv & mask).count_ones() as usize;
    top_score
        .checked_add(positive)
        .and_then(|score| score.checked_sub(negative))
        .ok_or(MyersError::StateInvariant { column })
}

#[derive(Clone, Copy)]
struct GlobalBand {
    pattern_len: usize,
    text_len: usize,
    max_distance: usize,
    fast_arithmetic: bool,
}

#[inline]
fn maximum_band_dimension(pattern_len: usize, text_len: usize) -> bool {
    pattern_len.max(text_len) <= (usize::MAX - WORD_BITS - 1) / 2
}

#[inline]
fn block_is_above_global_band(
    score: usize,
    width: usize,
    last_row: usize,
    column: usize,
    band: GlobalBand,
) -> bool {
    if band.fast_arithmetic {
        if score >= band.max_distance + width {
            return true;
        }
        let remaining_rows = band.pattern_len - last_row;
        let remaining_text = band.text_len - column;
        return score + remaining_rows > band.max_distance + remaining_text;
    }

    block_is_above_global_band_slow(score, width, last_row, column, band)
}

#[cold]
#[inline(never)]
fn block_is_above_global_band_slow(
    score: usize,
    width: usize,
    last_row: usize,
    column: usize,
    band: GlobalBand,
) -> bool {
    if band
        .max_distance
        .checked_add(width)
        .is_some_and(|threshold| score >= threshold)
    {
        return true;
    }
    if let (Some(left), Some(right)) = (
        last_row
            .checked_add(band.max_distance)
            .and_then(|value| value.checked_add(band.text_len)),
        score
            .checked_add(band.pattern_len)
            .and_then(|value| value.checked_add(column)),
    ) {
        left < right
    } else {
        let left = last_row as u128 + band.max_distance as u128 + band.text_len as u128;
        let right = score as u128 + band.pattern_len as u128 + column as u128;
        left < right
    }
}

#[inline]
fn block_is_below_global_band(
    score: usize,
    previous_score: usize,
    width: usize,
    last_row: usize,
    column: usize,
    band: GlobalBand,
) -> bool {
    if band.fast_arithmetic {
        if score >= band.max_distance + width {
            return true;
        }
        let previous_lower_bound = previous_score.saturating_sub(1);
        let within_block_lower_bound = (score + 1).saturating_sub(width);
        let first_score_lower_bound = previous_lower_bound.max(within_block_lower_bound);
        let remaining_rows = band.pattern_len - last_row;
        let remaining_text = band.text_len - column;
        return first_score_lower_bound + remaining_text + 1
            > band.max_distance + remaining_rows + width;
    }

    block_is_below_global_band_slow(score, previous_score, width, last_row, column, band)
}

#[cold]
#[inline(never)]
fn block_is_below_global_band_slow(
    score: usize,
    previous_score: usize,
    width: usize,
    last_row: usize,
    column: usize,
    band: GlobalBand,
) -> bool {
    if band
        .max_distance
        .checked_add(width)
        .is_some_and(|threshold| score >= threshold)
    {
        return true;
    }
    let previous_lower_bound = previous_score.saturating_sub(1);
    let within_block_lower_bound = score
        .checked_add(1)
        .map_or(usize::MAX, |value| value.saturating_sub(width));
    let first_score_lower_bound = previous_lower_bound.max(within_block_lower_bound);
    if let (Some(left), Some(right)) = (
        last_row
            .checked_add(first_score_lower_bound)
            .and_then(|value| value.checked_add(band.text_len))
            .and_then(|value| value.checked_add(1)),
        band.max_distance
            .checked_add(band.pattern_len)
            .and_then(|value| value.checked_add(column))
            .and_then(|value| value.checked_add(width)),
    ) {
        left > right
    } else {
        let left = last_row as u128 + first_score_lower_bound as u128 + band.text_len as u128 + 1;
        let right =
            band.max_distance as u128 + band.pattern_len as u128 + column as u128 + width as u128;
        left > right
    }
}

fn block_count(length: usize) -> usize {
    length / WORD_BITS + usize::from(length % WORD_BITS != 0)
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn advance_banded_blocks<const TRACE: bool, const FULL: bool>(
    states: &mut [BlockState],
    scores: &mut [usize],
    equalities: &[u64],
    history_scratch: &mut [BlockHistory],
    first_block: usize,
    blocks: usize,
    last_mask: u64,
    last_high_bit: u64,
) -> i8 {
    let mut horizontal = 1_i8;
    for (offset, ((state, score), &equality)) in
        states.iter_mut().zip(scores).zip(equalities).enumerate()
    {
        let block = first_block + offset;
        let (mask, high_bit) = if FULL || block + 1 != blocks {
            (u64::MAX, 1_u64 << (WORD_BITS - 1))
        } else {
            (last_mask, last_high_bit)
        };
        let horizontal_input = horizontal;
        let (output, ph, mh) = advance_block(state, equality, horizontal_input, mask, high_bit);
        if TRACE {
            history_scratch[block] =
                trace_decisions(*state, equality, horizontal_input, ph, mh, mask);
        }
        debug_assert!((-1..=1).contains(&output));
        debug_assert!(output != -1 || *score != 0);
        debug_assert!(output != 1 || *score != usize::MAX);
        *score = score.wrapping_add_signed(isize::from(output));
        horizontal = output;
    }
    horizontal
}

#[inline]
fn advance_block(
    state: &mut BlockState,
    mut eq: u64,
    horizontal_input: i8,
    mask: u64,
    high_bit: u64,
) -> (i8, u64, u64) {
    debug_assert!((-1..=1).contains(&horizontal_input));
    eq &= mask;
    if horizontal_input < 0 {
        eq |= 1;
    }

    let xv = eq | state.mv;
    let xh = (((eq & state.pv).wrapping_add(state.pv)) ^ state.pv) | eq;
    let ph = (state.mv | !(xh | state.pv)) & mask;
    let mh = (state.pv & xh) & mask;
    let horizontal_output = i8::from(ph & high_bit != 0) - i8::from(mh & high_bit != 0);

    let shifted_positive = ((ph << 1) | u64::from(horizontal_input > 0)) & mask;
    let shifted_negative = ((mh << 1) | u64::from(horizontal_input < 0)) & mask;
    state.pv = (shifted_negative | !(xv | shifted_positive)) & mask;
    state.mv = (shifted_positive & xv) & mask;

    (horizontal_output, ph, mh)
}

#[inline]
fn update_score(score: usize, difference: i8) -> Option<usize> {
    match difference {
        -1 => score.checked_sub(1),
        0 => Some(score),
        1 => score.checked_add(1),
        _ => None,
    }
}

fn history_at(
    history: &[BlockHistory],
    blocks: usize,
    row: usize,
    column: usize,
) -> Option<BlockHistory> {
    let block = (row - 1) / WORD_BITS;
    let column_start = (column - 1).checked_mul(blocks)?;
    history.get(column_start.checked_add(block)?).copied()
}

fn banded_history_at(
    history: &[BlockHistory],
    columns: &[HistoryColumn],
    row: usize,
    column: usize,
) -> Option<BlockHistory> {
    let metadata = columns.get(column.checked_sub(1)?)?;
    let block = row.checked_sub(1)? / WORD_BITS;
    let relative_block = block.checked_sub(metadata.first_block)?;
    let history_start = if column == 1 {
        0
    } else {
        columns.get(column - 2)?.history_end
    };
    let index = history_start.checked_add(relative_block)?;
    (index < metadata.history_end)
        .then(|| history.get(index).copied())
        .flatten()
}

#[cfg(test)]
mod tests {
    use align_core::DnaSequence;

    use super::{BlockState, MyersPattern, MyersWorkspace, advance_block};

    fn patterned_dna(length: usize, seed: usize) -> Vec<u8> {
        (0..length)
            .map(|position| b"ACGT"[(position.wrapping_mul(17).wrapping_add(seed)) % 4])
            .collect()
    }

    fn related_text(pattern: &[u8], text_len: usize, seed: usize) -> Vec<u8> {
        let mut text = pattern[..pattern.len().min(text_len)].to_vec();
        text.extend(
            (text.len()..text_len)
                .map(|position| b"TGCA"[(position.wrapping_mul(13).wrapping_add(seed)) % 4]),
        );
        for position in (seed % 19..text.len()).step_by(29 + seed % 7) {
            let symbol = text
                .get_mut(position)
                .expect("generated mutation position is in bounds");
            *symbol = match *symbol {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        text
    }

    fn random_dna(length: usize, state: &mut u64) -> Vec<u8> {
        (0..length)
            .map(|_| {
                *state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                b"ACGT"[(*state >> 62) as usize]
            })
            .collect()
    }

    #[test]
    fn one_row_global_match_advances_to_zero() {
        let mut state = BlockState { pv: 1, mv: 0 };
        let (output, ph, mh) = advance_block(&mut state, 1, 1, 1, 1);
        assert_eq!(output, -1);
        assert_eq!(ph, 0);
        assert_eq!(mh, 1);
    }

    #[test]
    fn peq_preprocessing_crosses_a_word_boundary() {
        let pattern =
            DnaSequence::from_ascii(format!("{}C", "A".repeat(64))).expect("canonical DNA");
        let myers = MyersPattern::new(&pattern).expect("preprocessing succeeds");
        assert_eq!(myers.blocks, 2);
        assert_eq!(myers.peq[0], u64::MAX);
        assert_eq!(myers.peq[myers.blocks + 1], 1);
    }

    #[test]
    fn bounded_global_distance_matches_full_frontier_at_word_boundaries() {
        let mut workspace = MyersWorkspace::new();
        let mut full_workspace = MyersWorkspace::new();
        for pattern_len in [65_usize, 127, 128, 129, 191, 192, 193, 255, 256, 257, 513] {
            for delta in [-17_isize, -1, 0, 1, 17] {
                let text_len = pattern_len.saturating_add_signed(delta);
                for seed in 0..4 {
                    let pattern_bytes = patterned_dna(pattern_len, seed);
                    let text_bytes = related_text(&pattern_bytes, text_len, seed + 1);
                    let pattern =
                        DnaSequence::from_ascii(pattern_bytes).expect("generated DNA is valid");
                    let text = DnaSequence::from_ascii(text_bytes).expect("generated DNA is valid");
                    let myers = MyersPattern::new(&pattern).expect("preprocessing succeeds");
                    let expected = myers
                        .global_distance_full_with_workspace(&text, &mut full_workspace)
                        .expect("full distance succeeds");

                    let mut limits = vec![
                        0,
                        pattern_len.abs_diff(text_len),
                        expected.saturating_sub(1),
                        expected,
                        expected.saturating_add(1),
                        pattern_len.max(text_len),
                        usize::MAX,
                    ];
                    limits.sort_unstable();
                    limits.dedup();
                    for limit in limits {
                        let observed = myers
                            .global_distance_bounded_with_workspace(&text, limit, &mut workspace)
                            .expect("bounded distance succeeds");
                        assert_eq!(
                            observed,
                            (expected <= limit).then_some(expected),
                            "pattern_len={pattern_len}, text_len={text_len}, seed={seed}, limit={limit}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn bounded_global_distance_matches_full_frontier_on_diverse_inputs() {
        let mut state = 0x6a09_e667_f3bc_c909_u64;
        let mut workspace = MyersWorkspace::new();
        let mut full_workspace = MyersWorkspace::new();
        for case in 0..160 {
            state = state.rotate_left(17).wrapping_add(case);
            let pattern_len =
                65 + usize::try_from(state % 512).expect("bounded pattern remainder fits usize");
            state = state.rotate_left(23).wrapping_add(0x9e37_79b9);
            let text_len =
                1 + usize::try_from(state % 640).expect("bounded text remainder fits usize");
            let pattern = DnaSequence::from_ascii(random_dna(pattern_len, &mut state))
                .expect("generated DNA is valid");
            let text = DnaSequence::from_ascii(random_dna(text_len, &mut state))
                .expect("generated DNA is valid");
            let myers = MyersPattern::new(&pattern).expect("preprocessing succeeds");
            let expected = myers
                .global_distance_full_with_workspace(&text, &mut full_workspace)
                .expect("full distance succeeds");

            for limit in [
                pattern_len.abs_diff(text_len),
                expected / 2,
                expected.saturating_sub(1),
                expected,
                expected.saturating_add(17),
            ] {
                let observed = myers
                    .global_distance_bounded_with_workspace(&text, limit, &mut workspace)
                    .expect("bounded distance succeeds");
                assert_eq!(
                    observed,
                    (expected <= limit).then_some(expected),
                    "case={case}, pattern_len={pattern_len}, text_len={text_len}, limit={limit}"
                );
            }
        }
    }

    #[test]
    fn automatic_global_distance_matches_full_frontier_for_long_patterns() {
        let mut state = 0xbb67_ae85_84ca_a73b_u64;
        let mut automatic_workspace = MyersWorkspace::new();
        let mut full_workspace = MyersWorkspace::new();
        for (pattern_len, text_len) in [
            (512, 480),
            (512, 544),
            (1_024, 1_024),
            (2_049, 1_907),
            (4_096, 4_096),
        ] {
            let pattern = DnaSequence::from_ascii(random_dna(pattern_len, &mut state))
                .expect("generated DNA is valid");
            let text = DnaSequence::from_ascii(random_dna(text_len, &mut state))
                .expect("generated DNA is valid");
            let myers = MyersPattern::new(&pattern).expect("preprocessing succeeds");
            let expected = myers
                .global_distance_full_with_workspace(&text, &mut full_workspace)
                .expect("full distance succeeds");
            assert_eq!(
                myers
                    .global_distance_with_workspace(&text, &mut automatic_workspace)
                    .expect("automatic distance succeeds"),
                expected,
                "pattern_len={pattern_len}, text_len={text_len}"
            );
        }
    }
}
