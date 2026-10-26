use core::cell::Cell;
use core::mem::size_of;
use std::collections::BTreeMap;

use align_core::Sequence;
use align_msa::{DistanceMatrix, SequenceRecord};

use crate::Muscle5Error;
use crate::muscle5_emissions::Muscle5Alphabet;
use crate::muscle5_hmm::{
    LONG_X_STATE, LONG_Y_STATE, MATCH_STATE, PreparedHmm, SHORT_X_STATE, SHORT_Y_STATE, STATE_COUNT,
};

const FNV_OFFSET: u64 = 0xCBF2_9CE4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01B3;
const PARALLEL_WORK_THRESHOLD: usize = 65_536;
const MAX_PARALLEL_WORKERS: usize = 24;
const PROBABILITY_FAST_PATH_MIN_LOG: f64 = -600.0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct PosteriorEntry {
    pub(crate) column: usize,
    pub(crate) probability: f64,
}

/// Compressed sparse rows of match posterior probabilities.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct SparsePosterior {
    rows: usize,
    columns: usize,
    offsets: Box<[usize]>,
    entries: Box<[PosteriorEntry]>,
}

impl SparsePosterior {
    pub(crate) fn from_dense(
        rows: usize,
        columns: usize,
        values: &[f64],
        cutoff: f64,
    ) -> Result<Self, Muscle5Error> {
        let expected = rows
            .checked_mul(columns)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "posterior matrix dimensions",
            })?;
        if values.len() != expected {
            return Err(Muscle5Error::InternalInvariant {
                message: "dense posterior shape does not match its dimensions",
            });
        }
        let mut offsets = Vec::with_capacity(rows + 1);
        let mut entries = Vec::new();
        offsets.push(0);
        for row in 0..rows {
            for column in 0..columns {
                let probability = values[row * columns + column];
                if !probability.is_finite() || !(0.0..=1.0 + 1.0e-12).contains(&probability) {
                    return Err(Muscle5Error::NonFiniteComputation {
                        stage: "pair posterior probability",
                    });
                }
                let probability = probability.clamp(0.0, 1.0);
                if probability >= cutoff && probability > 0.0 {
                    entries.push(PosteriorEntry {
                        column,
                        probability,
                    });
                }
            }
            offsets.push(entries.len());
        }
        Ok(Self {
            rows,
            columns,
            offsets: offsets.into_boxed_slice(),
            entries: entries.into_boxed_slice(),
        })
    }

    pub(crate) fn identity(size: usize) -> Self {
        let offsets = (0..=size).collect::<Vec<_>>().into_boxed_slice();
        let entries = (0..size)
            .map(|column| PosteriorEntry {
                column,
                probability: 1.0,
            })
            .collect::<Vec<_>>()
            .into_boxed_slice();
        Self {
            rows: size,
            columns: size,
            offsets,
            entries,
        }
    }

    pub(crate) fn transpose(&self) -> Result<Self, Muscle5Error> {
        let mut counts = vec![0_usize; self.columns];
        for entry in &self.entries {
            counts[entry.column] += 1;
        }
        let mut offsets = Vec::with_capacity(self.columns + 1);
        offsets.push(0);
        for count in counts {
            let next = offsets
                .last()
                .copied()
                .unwrap_or(0_usize)
                .checked_add(count)
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "posterior transpose offsets",
                })?;
            offsets.push(next);
        }
        let mut cursors = offsets[..self.columns].to_vec();
        let mut entries = vec![
            PosteriorEntry {
                column: 0,
                probability: 0.0,
            };
            self.entries.len()
        ];
        for row in 0..self.rows {
            for entry in self.row(row) {
                let destination = cursors[entry.column];
                entries[destination] = PosteriorEntry {
                    column: row,
                    probability: entry.probability,
                };
                cursors[entry.column] += 1;
            }
        }
        Ok(Self {
            rows: self.columns,
            columns: self.rows,
            offsets: offsets.into_boxed_slice(),
            entries: entries.into_boxed_slice(),
        })
    }

    pub(crate) const fn rows(&self) -> usize {
        self.rows
    }

    pub(crate) const fn columns(&self) -> usize {
        self.columns
    }

    pub(crate) fn row(&self, row: usize) -> &[PosteriorEntry] {
        if row >= self.rows {
            return &[];
        }
        &self.entries[self.offsets[row]..self.offsets[row + 1]]
    }

    pub(crate) fn get(&self, row: usize, column: usize) -> f64 {
        self.row(row)
            .binary_search_by_key(&column, |entry| entry.column)
            .map_or(0.0, |index| self.row(row)[index].probability)
    }

    pub(crate) const fn entry_count(&self) -> usize {
        self.entries.len()
    }

    pub(crate) fn storage_bytes(&self) -> usize {
        self.offsets.len() * size_of::<usize>() + self.entries.len() * size_of::<PosteriorEntry>()
    }

    #[allow(clippy::cast_possible_truncation)]
    pub(crate) fn checksum(&self) -> u64 {
        let mut hash = FNV_OFFSET;
        for value in [
            self.rows as u64,
            self.columns as u64,
            self.entries.len() as u64,
        ] {
            hash = hash_u64(hash, value);
        }
        for row in 0..self.rows {
            for entry in self.row(row) {
                hash = hash_u64(hash, row as u64);
                hash = hash_u64(hash, entry.column as u64);
                hash = hash_u64(hash, entry.probability.to_bits());
            }
        }
        hash
    }
}

fn hash_u64(mut hash: u64, value: u64) -> u64 {
    for byte in value.to_le_bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    hash
}

#[derive(Clone, Debug)]
pub(crate) struct PairPosteriorCalculation {
    pub(crate) posterior: SparsePosterior,
    pub(crate) log_partition: f64,
    pub(crate) peak_scratch_bytes: usize,
}

/// Pair posteriors computed earlier in one Super5 run, for exact reuse.
///
/// A pair posterior is a pure function of the two sequences, in this order, and
/// of the HMM and cutoff, so a repeated request can return the stored result
/// bit for bit. The cache answers only for the HMM and cutoff it was made for.
#[derive(Debug)]
pub(crate) struct PairPosteriorCache<A: Muscle5Alphabet> {
    hmm: PreparedHmm<A>,
    cutoff: f64,
    reuse: bool,
    entries: BTreeMap<PairKey, PairPosteriorCalculation>,
    hits: Cell<usize>,
}

/// The symbol indices of the first and the second sequence of a pair.
type PairKey = (Box<[u8]>, Box<[u8]>);

impl<A: Muscle5Alphabet> PairPosteriorCache<A> {
    pub(crate) fn new(hmm: &PreparedHmm<A>, cutoff: f64) -> Self {
        Self {
            hmm: hmm.clone(),
            cutoff,
            reuse: true,
            entries: BTreeMap::new(),
            hits: Cell::new(0),
        }
    }

    /// A cache that never stores anything, so every posterior is recomputed;
    /// tests compare runs with and without reuse.
    #[cfg(test)]
    pub(crate) fn without_reuse(hmm: &PreparedHmm<A>, cutoff: f64) -> Self {
        Self {
            reuse: false,
            ..Self::new(hmm, cutoff)
        }
    }

    /// The number of requests answered from the cache.
    #[cfg(test)]
    pub(crate) fn hits(&self) -> usize {
        self.hits.get()
    }

    fn serves(&self, hmm: &PreparedHmm<A>, cutoff: f64) -> bool {
        self.reuse && self.cutoff.to_bits() == cutoff.to_bits() && self.hmm == *hmm
    }

    fn key(first: &Sequence<A>, second: &Sequence<A>) -> PairKey {
        (first.as_indices().into(), second.as_indices().into())
    }

    fn get(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        hmm: &PreparedHmm<A>,
        cutoff: f64,
    ) -> Option<&PairPosteriorCalculation> {
        if !self.serves(hmm, cutoff) {
            return None;
        }
        let found = self.entries.get(&Self::key(first, second));
        if found.is_some() {
            self.hits.set(self.hits.get() + 1);
        }
        found
    }

    /// Returns the stored posterior of the pair, or computes and stores it.
    pub(crate) fn pair_posterior(
        &mut self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        hmm: &PreparedHmm<A>,
        cutoff: f64,
    ) -> Result<PairPosteriorCalculation, Muscle5Error> {
        if let Some(stored) = self.get(first, second, hmm, cutoff) {
            return Ok(stored.clone());
        }
        let calculation = pair_posterior(first, second, hmm, cutoff)?;
        if self.serves(hmm, cutoff) {
            self.entries
                .insert(Self::key(first, second), calculation.clone());
        }
        Ok(calculation)
    }

    /// Returns the stored posterior of the pair without computing one.
    pub(crate) fn stored(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        hmm: &PreparedHmm<A>,
        cutoff: f64,
    ) -> Option<PairPosteriorCalculation> {
        self.get(first, second, hmm, cutoff).cloned()
    }
}

struct IndexedPairPosterior {
    first: usize,
    second: usize,
    calculation: PairPosteriorCalculation,
}

#[derive(Default)]
struct PairHmmWorkspace {
    forward: Vec<f64>,
    backward: Vec<f64>,
    dense: Vec<f64>,
}

struct ConsistentPair {
    first: usize,
    second: usize,
    forward: SparsePosterior,
    reverse: SparsePosterior,
    accumulator_bytes: usize,
}

pub(crate) fn pair_posterior<A: Muscle5Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    hmm: &PreparedHmm<A>,
    cutoff: f64,
) -> Result<PairPosteriorCalculation, Muscle5Error> {
    pair_posterior_with_workspace(first, second, hmm, cutoff, &mut PairHmmWorkspace::default())
}

fn pair_posterior_with_workspace<A: Muscle5Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    hmm: &PreparedHmm<A>,
    cutoff: f64,
    workspace: &mut PairHmmWorkspace,
) -> Result<PairPosteriorCalculation, Muscle5Error> {
    if concrete_path_log_probability(first, second, hmm) >= PROBABILITY_FAST_PATH_MIN_LOG {
        pair_posterior_probability(first, second, hmm, cutoff, workspace)
    } else {
        pair_posterior_log(first, second, hmm, cutoff, workspace)
    }
}

fn concrete_path_log_probability<A: Muscle5Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    hmm: &PreparedHmm<A>,
) -> f64 {
    let shared = first.len().min(second.len());
    if shared == 0 {
        return if first.is_empty() && second.is_empty() {
            0.0
        } else if second.is_empty() {
            gap_run_log(first.as_slice(), None, SHORT_X_STATE, LONG_X_STATE, hmm)
        } else {
            gap_run_log(second.as_slice(), None, SHORT_Y_STATE, LONG_Y_STATE, hmm)
        };
    }

    let mut path =
        hmm.start_log[MATCH_STATE] + hmm.match_log(first.as_slice()[0], second.as_slice()[0]);
    for position in 1..shared {
        path += hmm.transition_log[MATCH_STATE][MATCH_STATE]
            + hmm.match_log(first.as_slice()[position], second.as_slice()[position]);
    }
    if first.len() > shared {
        path + gap_run_log(
            &first.as_slice()[shared..],
            Some(MATCH_STATE),
            SHORT_X_STATE,
            LONG_X_STATE,
            hmm,
        )
    } else if second.len() > shared {
        path + gap_run_log(
            &second.as_slice()[shared..],
            Some(MATCH_STATE),
            SHORT_Y_STATE,
            LONG_Y_STATE,
            hmm,
        )
    } else {
        path + hmm.end_log[MATCH_STATE]
    }
}

fn gap_run_log<A: Muscle5Alphabet>(
    symbols: &[align_core::Symbol<A>],
    predecessor: Option<usize>,
    short_state: usize,
    long_state: usize,
    hmm: &PreparedHmm<A>,
) -> f64 {
    [short_state, long_state]
        .into_iter()
        .map(|state| {
            let first_transition = predecessor.map_or(hmm.start_log[state], |previous| {
                hmm.transition_log[previous][state]
            });
            let mut path = first_transition + hmm.insert_log(symbols[0]);
            for symbol in &symbols[1..] {
                path += hmm.transition_log[state][state] + hmm.insert_log(*symbol);
            }
            path + hmm.end_log[state]
        })
        .fold(f64::NEG_INFINITY, f64::max)
}

#[allow(clippy::similar_names, clippy::too_many_lines)]
fn pair_posterior_probability<A: Muscle5Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    hmm: &PreparedHmm<A>,
    cutoff: f64,
    workspace: &mut PairHmmWorkspace,
) -> Result<PairPosteriorCalculation, Muscle5Error> {
    let rows = first.len();
    let columns = second.len();
    if rows == 0 && columns == 0 {
        return Ok(PairPosteriorCalculation {
            posterior: SparsePosterior::from_dense(0, 0, &[], cutoff)?,
            log_partition: 0.0,
            peak_scratch_bytes: 0,
        });
    }
    let lattice_rows = rows
        .checked_add(1)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "pair-HMM row count",
        })?;
    let lattice_columns = columns
        .checked_add(1)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "pair-HMM column count",
        })?;
    let lattice_cells = lattice_rows
        .checked_mul(lattice_columns)
        .and_then(|cells| cells.checked_mul(STATE_COUNT))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "pair-HMM lattice size",
        })?;
    workspace.forward.resize(lattice_cells, 0.0);
    workspace.forward[..lattice_cells].fill(0.0);
    workspace.backward.resize(lattice_cells, 0.0);
    workspace.backward[..lattice_cells].fill(0.0);
    let forward = &mut workspace.forward;
    let backward = &mut workspace.backward;
    let first_symbols = first.as_slice();
    let second_symbols = second.as_slice();

    for consumed_first in 0..=rows {
        for consumed_second in 0..=columns {
            if consumed_first == 0 && consumed_second == 0 {
                continue;
            }
            let current = lattice_index(consumed_first, consumed_second, 0, lattice_columns);
            if consumed_first > 0 && consumed_second > 0 {
                let previous =
                    lattice_index(consumed_first - 1, consumed_second - 1, 0, lattice_columns);
                let mut prefix = if consumed_first == 1 && consumed_second == 1 {
                    hmm.start_probability[MATCH_STATE]
                } else {
                    0.0
                };
                for previous_state in 0..STATE_COUNT {
                    prefix += forward[previous + previous_state]
                        * hmm.transition_probability[previous_state][MATCH_STATE];
                }
                forward[current + MATCH_STATE] = prefix
                    * hmm.match_probability(
                        first_symbols[consumed_first - 1],
                        second_symbols[consumed_second - 1],
                    );
            }
            if consumed_first > 0 {
                let previous =
                    lattice_index(consumed_first - 1, consumed_second, 0, lattice_columns);
                let insertion = hmm.insert_probability(first_symbols[consumed_first - 1]);
                for state in [SHORT_X_STATE, LONG_X_STATE] {
                    let mut prefix = if consumed_first == 1 && consumed_second == 0 {
                        hmm.start_probability[state]
                    } else {
                        0.0
                    };
                    prefix += forward[previous + MATCH_STATE]
                        * hmm.transition_probability[MATCH_STATE][state];
                    prefix += forward[previous + state] * hmm.transition_probability[state][state];
                    forward[current + state] = prefix * insertion;
                }
            }
            if consumed_second > 0 {
                let previous =
                    lattice_index(consumed_first, consumed_second - 1, 0, lattice_columns);
                let insertion = hmm.insert_probability(second_symbols[consumed_second - 1]);
                for state in [SHORT_Y_STATE, LONG_Y_STATE] {
                    let mut prefix = if consumed_first == 0 && consumed_second == 1 {
                        hmm.start_probability[state]
                    } else {
                        0.0
                    };
                    prefix += forward[previous + MATCH_STATE]
                        * hmm.transition_probability[MATCH_STATE][state];
                    prefix += forward[previous + state] * hmm.transition_probability[state][state];
                    forward[current + state] = prefix * insertion;
                }
            }
        }
    }

    let terminal = lattice_index(rows, columns, 0, lattice_columns);
    backward[terminal..terminal + STATE_COUNT].copy_from_slice(&hmm.end_probability);
    for consumed_first in (0..=rows).rev() {
        for consumed_second in (0..=columns).rev() {
            if consumed_first == rows && consumed_second == columns {
                continue;
            }
            let current = lattice_index(consumed_first, consumed_second, 0, lattice_columns);
            let match_next = if consumed_first < rows && consumed_second < columns {
                let next = lattice_index(
                    consumed_first + 1,
                    consumed_second + 1,
                    MATCH_STATE,
                    lattice_columns,
                );
                hmm.match_probability(
                    first_symbols[consumed_first],
                    second_symbols[consumed_second],
                ) * backward[next]
            } else {
                0.0
            };
            let short_x_next = insertion_probability_tail(
                consumed_first,
                consumed_second,
                rows,
                lattice_columns,
                SHORT_X_STATE,
                first_symbols,
                backward,
                hmm,
            );
            let short_y_next = insertion_probability_tail(
                consumed_second,
                consumed_first,
                columns,
                lattice_columns,
                SHORT_Y_STATE,
                second_symbols,
                backward,
                hmm,
            );
            let long_x_next = insertion_probability_tail(
                consumed_first,
                consumed_second,
                rows,
                lattice_columns,
                LONG_X_STATE,
                first_symbols,
                backward,
                hmm,
            );
            let long_y_next = insertion_probability_tail(
                consumed_second,
                consumed_first,
                columns,
                lattice_columns,
                LONG_Y_STATE,
                second_symbols,
                backward,
                hmm,
            );

            backward[current + MATCH_STATE] = hmm.transition_probability[MATCH_STATE][MATCH_STATE]
                * match_next
                + hmm.transition_probability[MATCH_STATE][SHORT_X_STATE] * short_x_next
                + hmm.transition_probability[MATCH_STATE][SHORT_Y_STATE] * short_y_next
                + hmm.transition_probability[MATCH_STATE][LONG_X_STATE] * long_x_next
                + hmm.transition_probability[MATCH_STATE][LONG_Y_STATE] * long_y_next;
            for (state, insertion_next) in [
                (SHORT_X_STATE, short_x_next),
                (SHORT_Y_STATE, short_y_next),
                (LONG_X_STATE, long_x_next),
                (LONG_Y_STATE, long_y_next),
            ] {
                backward[current + state] = hmm.transition_probability[state][MATCH_STATE]
                    * match_next
                    + hmm.transition_probability[state][state] * insertion_next;
            }
        }
    }

    let mut partition = 0.0_f64;
    for state in 0..STATE_COUNT {
        partition += forward[terminal + state] * hmm.end_probability[state];
    }
    if !partition.is_finite() || partition <= 0.0 {
        return Err(Muscle5Error::NonFiniteComputation {
            stage: "probability-domain pair-HMM partition function",
        });
    }
    let dense_len = rows
        .checked_mul(columns)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "dense posterior dimensions",
        })?;
    workspace.dense.resize(dense_len, 0.0);
    let dense = &mut workspace.dense;
    for first_position in 1..=rows {
        for second_position in 1..=columns {
            let index = lattice_index(
                first_position,
                second_position,
                MATCH_STATE,
                lattice_columns,
            );
            dense[(first_position - 1) * columns + second_position - 1] =
                (forward[index] / partition * backward[index]).clamp(0.0, 1.0);
        }
    }
    let peak_scratch_bytes = pair_scratch_bytes(lattice_cells, dense_len)?;
    let posterior = SparsePosterior::from_dense(rows, columns, dense, cutoff)?;
    Ok(PairPosteriorCalculation {
        posterior,
        log_partition: partition.ln(),
        peak_scratch_bytes,
    })
}

#[allow(clippy::similar_names, clippy::too_many_lines)]
fn pair_posterior_log<A: Muscle5Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    hmm: &PreparedHmm<A>,
    cutoff: f64,
    workspace: &mut PairHmmWorkspace,
) -> Result<PairPosteriorCalculation, Muscle5Error> {
    let rows = first.len();
    let columns = second.len();
    if rows == 0 && columns == 0 {
        return Ok(PairPosteriorCalculation {
            posterior: SparsePosterior::from_dense(0, 0, &[], cutoff)?,
            log_partition: 0.0,
            peak_scratch_bytes: 0,
        });
    }
    let lattice_rows = rows
        .checked_add(1)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "pair-HMM row count",
        })?;
    let lattice_columns = columns
        .checked_add(1)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "pair-HMM column count",
        })?;
    let lattice_cells = lattice_rows
        .checked_mul(lattice_columns)
        .and_then(|cells| cells.checked_mul(STATE_COUNT))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "pair-HMM lattice size",
        })?;
    workspace.forward.resize(lattice_cells, f64::NEG_INFINITY);
    workspace.forward[..lattice_cells].fill(f64::NEG_INFINITY);
    workspace.backward.resize(lattice_cells, f64::NEG_INFINITY);
    workspace.backward[..lattice_cells].fill(f64::NEG_INFINITY);
    let forward = &mut workspace.forward;
    let backward = &mut workspace.backward;
    let first_symbols = first.as_slice();
    let second_symbols = second.as_slice();

    // The published topology is fixed: M can follow every state, while each
    // short/long insertion can follow only M or itself. Expanding that topology
    // here avoids testing the 12 impossible transitions in every lattice cell.
    for consumed_first in 0..=rows {
        for consumed_second in 0..=columns {
            if consumed_first == 0 && consumed_second == 0 {
                continue;
            }
            let current = lattice_index(consumed_first, consumed_second, 0, lattice_columns);
            if consumed_first > 0 && consumed_second > 0 {
                let previous =
                    lattice_index(consumed_first - 1, consumed_second - 1, 0, lattice_columns);
                let mut prefix = if consumed_first == 1 && consumed_second == 1 {
                    hmm.start_log[MATCH_STATE]
                } else {
                    f64::NEG_INFINITY
                };
                for previous_state in 0..STATE_COUNT {
                    prefix = log_add(
                        prefix,
                        forward[previous + previous_state]
                            + hmm.transition_log[previous_state][MATCH_STATE],
                    );
                }
                forward[current + MATCH_STATE] = prefix
                    + hmm.match_log(
                        first_symbols[consumed_first - 1],
                        second_symbols[consumed_second - 1],
                    );
            }
            if consumed_first > 0 {
                let previous =
                    lattice_index(consumed_first - 1, consumed_second, 0, lattice_columns);
                let insertion = hmm.insert_log(first_symbols[consumed_first - 1]);
                for state in [SHORT_X_STATE, LONG_X_STATE] {
                    let mut prefix = if consumed_first == 1 && consumed_second == 0 {
                        hmm.start_log[state]
                    } else {
                        f64::NEG_INFINITY
                    };
                    prefix = log_add(
                        prefix,
                        forward[previous + MATCH_STATE] + hmm.transition_log[MATCH_STATE][state],
                    );
                    prefix = log_add(
                        prefix,
                        forward[previous + state] + hmm.transition_log[state][state],
                    );
                    forward[current + state] = prefix + insertion;
                }
            }
            if consumed_second > 0 {
                let previous =
                    lattice_index(consumed_first, consumed_second - 1, 0, lattice_columns);
                let insertion = hmm.insert_log(second_symbols[consumed_second - 1]);
                for state in [SHORT_Y_STATE, LONG_Y_STATE] {
                    let mut prefix = if consumed_first == 0 && consumed_second == 1 {
                        hmm.start_log[state]
                    } else {
                        f64::NEG_INFINITY
                    };
                    prefix = log_add(
                        prefix,
                        forward[previous + MATCH_STATE] + hmm.transition_log[MATCH_STATE][state],
                    );
                    prefix = log_add(
                        prefix,
                        forward[previous + state] + hmm.transition_log[state][state],
                    );
                    forward[current + state] = prefix + insertion;
                }
            }
        }
    }

    let terminal = lattice_index(rows, columns, 0, lattice_columns);
    backward[terminal..terminal + STATE_COUNT].copy_from_slice(&hmm.end_log);
    for consumed_first in (0..=rows).rev() {
        for consumed_second in (0..=columns).rev() {
            if consumed_first == rows && consumed_second == columns {
                continue;
            }
            let current = lattice_index(consumed_first, consumed_second, 0, lattice_columns);
            let match_next = if consumed_first < rows && consumed_second < columns {
                let next = lattice_index(
                    consumed_first + 1,
                    consumed_second + 1,
                    MATCH_STATE,
                    lattice_columns,
                );
                (
                    hmm.match_log(
                        first_symbols[consumed_first],
                        second_symbols[consumed_second],
                    ),
                    backward[next],
                )
            } else {
                (0.0, f64::NEG_INFINITY)
            };
            let short_x_next = insertion_tail(
                consumed_first,
                consumed_second,
                rows,
                lattice_columns,
                SHORT_X_STATE,
                first_symbols,
                backward,
                hmm,
            );
            let short_y_next = insertion_tail(
                consumed_second,
                consumed_first,
                columns,
                lattice_columns,
                SHORT_Y_STATE,
                second_symbols,
                backward,
                hmm,
            );
            let long_x_next = insertion_tail(
                consumed_first,
                consumed_second,
                rows,
                lattice_columns,
                LONG_X_STATE,
                first_symbols,
                backward,
                hmm,
            );
            let long_y_next = insertion_tail(
                consumed_second,
                consumed_first,
                columns,
                lattice_columns,
                LONG_Y_STATE,
                second_symbols,
                backward,
                hmm,
            );

            let mut match_suffix =
                transition_tail(hmm.transition_log[MATCH_STATE][MATCH_STATE], match_next);
            match_suffix = log_add(
                match_suffix,
                transition_tail(hmm.transition_log[MATCH_STATE][SHORT_X_STATE], short_x_next),
            );
            match_suffix = log_add(
                match_suffix,
                transition_tail(hmm.transition_log[MATCH_STATE][SHORT_Y_STATE], short_y_next),
            );
            match_suffix = log_add(
                match_suffix,
                transition_tail(hmm.transition_log[MATCH_STATE][LONG_X_STATE], long_x_next),
            );
            match_suffix = log_add(
                match_suffix,
                transition_tail(hmm.transition_log[MATCH_STATE][LONG_Y_STATE], long_y_next),
            );
            backward[current + MATCH_STATE] = match_suffix;

            for (state, insertion_next) in [
                (SHORT_X_STATE, short_x_next),
                (SHORT_Y_STATE, short_y_next),
                (LONG_X_STATE, long_x_next),
                (LONG_Y_STATE, long_y_next),
            ] {
                let match_term =
                    transition_tail(hmm.transition_log[state][MATCH_STATE], match_next);
                let insertion_term =
                    transition_tail(hmm.transition_log[state][state], insertion_next);
                backward[current + state] = log_add(match_term, insertion_term);
            }
        }
    }

    let mut log_partition = f64::NEG_INFINITY;
    for state in 0..STATE_COUNT {
        log_partition = log_add(
            log_partition,
            forward[lattice_index(rows, columns, state, lattice_columns)] + hmm.end_log[state],
        );
    }
    if !log_partition.is_finite() {
        return Err(Muscle5Error::NonFiniteComputation {
            stage: "pair-HMM partition function",
        });
    }

    let dense_len = rows
        .checked_mul(columns)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "dense posterior dimensions",
        })?;
    workspace.dense.resize(dense_len, 0.0);
    let dense = &mut workspace.dense;
    for first_position in 1..=rows {
        for second_position in 1..=columns {
            let index = lattice_index(
                first_position,
                second_position,
                MATCH_STATE,
                lattice_columns,
            );
            let log_probability = forward[index] + backward[index] - log_partition;
            let probability = if log_probability < -745.0 {
                0.0
            } else {
                log_probability.exp().clamp(0.0, 1.0)
            };
            dense[(first_position - 1) * columns + second_position - 1] = probability;
        }
    }
    let peak_scratch_bytes = pair_scratch_bytes(lattice_cells, dense_len)?;
    let posterior = SparsePosterior::from_dense(rows, columns, dense, cutoff)?;
    Ok(PairPosteriorCalculation {
        posterior,
        log_partition,
        peak_scratch_bytes,
    })
}

fn pair_scratch_bytes(lattice_cells: usize, dense_len: usize) -> Result<usize, Muscle5Error> {
    lattice_cells
        .checked_mul(2)
        .and_then(|values| values.checked_add(dense_len))
        .and_then(|values| values.checked_mul(size_of::<f64>()))
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "pair-HMM scratch bytes",
        })
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn insertion_tail<A: Muscle5Alphabet>(
    consumed: usize,
    other_consumed: usize,
    length: usize,
    lattice_columns: usize,
    state: usize,
    symbols: &[align_core::Symbol<A>],
    backward: &[f64],
    hmm: &PreparedHmm<A>,
) -> (f64, f64) {
    if consumed >= length {
        return (0.0, f64::NEG_INFINITY);
    }
    let (next_first, next_second) = if state == SHORT_X_STATE || state == LONG_X_STATE {
        (consumed + 1, other_consumed)
    } else {
        (other_consumed, consumed + 1)
    };
    let next = lattice_index(next_first, next_second, state, lattice_columns);
    (hmm.insert_log(symbols[consumed]), backward[next])
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn insertion_probability_tail<A: Muscle5Alphabet>(
    consumed: usize,
    other_consumed: usize,
    length: usize,
    lattice_columns: usize,
    state: usize,
    symbols: &[align_core::Symbol<A>],
    backward: &[f64],
    hmm: &PreparedHmm<A>,
) -> f64 {
    if consumed >= length {
        return 0.0;
    }
    let (next_first, next_second) = if state == SHORT_X_STATE || state == LONG_X_STATE {
        (consumed + 1, other_consumed)
    } else {
        (other_consumed, consumed + 1)
    };
    let next = lattice_index(next_first, next_second, state, lattice_columns);
    hmm.insert_probability(symbols[consumed]) * backward[next]
}

#[inline]
fn transition_tail(transition: f64, emission_and_tail: (f64, f64)) -> f64 {
    let (emission, tail) = emission_and_tail;
    if tail.is_infinite() && tail.is_sign_negative() {
        f64::NEG_INFINITY
    } else {
        transition + emission + tail
    }
}

#[inline]
fn lattice_index(row: usize, column: usize, state: usize, columns: usize) -> usize {
    (row * columns + column) * STATE_COUNT + state
}

#[inline]
fn log_add(first: f64, second: f64) -> f64 {
    if first.is_infinite() && first.is_sign_negative() {
        return second;
    }
    if second.is_infinite() && second.is_sign_negative() {
        return first;
    }
    let maximum = first.max(second);
    let minimum = first.min(second);
    if maximum.is_infinite() {
        maximum
    } else {
        maximum + (minimum - maximum).exp().ln_1p()
    }
}

#[derive(Clone, Debug)]
pub(crate) struct PosteriorLibrary {
    sequence_count: usize,
    lengths: Box<[usize]>,
    matrices: Vec<Option<SparsePosterior>>,
    pair_hmm_peak_bytes: usize,
    consistency_peak_bytes: usize,
}

pub(crate) trait PosteriorProvider {
    fn posterior(&self, first: usize, second: usize) -> Option<&SparsePosterior>;

    fn sequence_length(&self, sequence_index: usize) -> Option<usize>;
}

impl PosteriorLibrary {
    #[cfg(test)]
    pub(crate) fn from_test_upper_triangle(
        lengths: &[usize],
        upper_triangle: Vec<SparsePosterior>,
    ) -> Result<Self, Muscle5Error> {
        let sequence_count = lengths.len();
        let expected_pairs = sequence_count
            .checked_mul(sequence_count.saturating_sub(1))
            .and_then(|value| value.checked_div(2))
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "test posterior pair count",
            })?;
        if upper_triangle.len() != expected_pairs {
            return Err(Muscle5Error::InternalInvariant {
                message: "test posterior triangle has the wrong size",
            });
        }
        let slot_count =
            sequence_count
                .checked_mul(sequence_count)
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "test posterior slots",
                })?;
        let mut matrices = vec![None; slot_count];
        for (index, length) in lengths.iter().copied().enumerate() {
            matrices[index * sequence_count + index] = Some(SparsePosterior::identity(length));
        }
        let mut supplied = upper_triangle.into_iter();
        for first in 0..sequence_count {
            for second in first + 1..sequence_count {
                let matrix = supplied.next().ok_or(Muscle5Error::InternalInvariant {
                    message: "test posterior triangle ended early",
                })?;
                if matrix.rows() != lengths[first] || matrix.columns() != lengths[second] {
                    return Err(Muscle5Error::InternalInvariant {
                        message: "test posterior matrix has the wrong shape",
                    });
                }
                let reverse = matrix.transpose()?;
                matrices[first * sequence_count + second] = Some(matrix);
                matrices[second * sequence_count + first] = Some(reverse);
            }
        }
        Ok(Self {
            sequence_count,
            lengths: lengths.to_vec().into_boxed_slice(),
            matrices,
            pair_hmm_peak_bytes: 0,
            consistency_peak_bytes: 0,
        })
    }

    /// Computes every pair posterior, taking pairs already in `cache` from it.
    ///
    /// With a cache, only the missing pairs are computed. The pair-HMM scratch
    /// estimate is still that of computing every pair with the same worker
    /// schedule, so the library and its evidence equal those built without one.
    pub(crate) fn build<A: Muscle5Alphabet>(
        records: &[SequenceRecord<A>],
        hmm: &PreparedHmm<A>,
        cutoff: f64,
        cache: Option<&PairPosteriorCache<A>>,
    ) -> Result<Self, Muscle5Error> {
        let sequence_count = records.len();
        let slot_count =
            sequence_count
                .checked_mul(sequence_count)
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "posterior library slots",
                })?;
        let mut matrices = vec![None; slot_count];
        let lengths = records
            .iter()
            .map(|record| record.sequence().len())
            .collect::<Vec<_>>();
        for index in 0..sequence_count {
            matrices[index * sequence_count + index] =
                Some(SparsePosterior::identity(lengths[index]));
        }
        let jobs = unordered_pairs(sequence_count);
        let work_units = pair_work_units(&jobs, &lengths)?;
        let workers = parallel_worker_count(jobs.len(), work_units);
        let (calculations, pair_hmm_peak_bytes) = match cache {
            None => calculate_pair_posteriors(&jobs, records, hmm, cutoff, workers)?,
            Some(cache) => {
                let stored = jobs
                    .iter()
                    .map(|&(first, second)| {
                        cache.stored(
                            records[first].sequence(),
                            records[second].sequence(),
                            hmm,
                            cutoff,
                        )
                    })
                    .collect::<Vec<_>>();
                let missing = jobs
                    .iter()
                    .zip(&stored)
                    .filter(|(_, stored)| stored.is_none())
                    .map(|(job, _)| *job)
                    .collect::<Vec<_>>();
                let (computed, _) = if missing.is_empty() {
                    (Vec::new(), 0)
                } else {
                    calculate_pair_posteriors(&missing, records, hmm, cutoff, workers)?
                };
                let mut computed = computed.into_iter();
                let mut calculations = Vec::with_capacity(jobs.len());
                for (&(first, second), stored) in jobs.iter().zip(stored) {
                    let calculation = match stored {
                        Some(calculation) => calculation,
                        None => {
                            computed
                                .next()
                                .ok_or(Muscle5Error::InternalInvariant {
                                    message: "reused posterior library lost a computed pair",
                                })?
                                .calculation
                        }
                    };
                    calculations.push(IndexedPairPosterior {
                        first,
                        second,
                        calculation,
                    });
                }
                let peak = scheduled_scratch_bytes(&calculations, workers)?;
                (calculations, peak)
            }
        };
        for indexed in calculations {
            let reverse = indexed.calculation.posterior.transpose()?;
            matrices[indexed.first * sequence_count + indexed.second] =
                Some(indexed.calculation.posterior);
            matrices[indexed.second * sequence_count + indexed.first] = Some(reverse);
        }
        Ok(Self {
            sequence_count,
            lengths: lengths.into_boxed_slice(),
            matrices,
            pair_hmm_peak_bytes,
            consistency_peak_bytes: 0,
        })
    }

    pub(crate) fn transform_consistency(
        &mut self,
        iterations: usize,
        cutoff: f64,
    ) -> Result<(), Muscle5Error> {
        for _ in 0..iterations {
            let (next, peak_bytes) = self.consistency_round(cutoff)?;
            *self = next;
            self.consistency_peak_bytes = self.consistency_peak_bytes.max(peak_bytes);
        }
        Ok(())
    }

    fn consistency_round(&self, cutoff: f64) -> Result<(Self, usize), Muscle5Error> {
        let slot_count = self.sequence_count.checked_mul(self.sequence_count).ok_or(
            Muscle5Error::ArithmeticOverflow {
                stage: "consistent posterior slots",
            },
        )?;
        let mut matrices = vec![None; slot_count];
        for index in 0..self.sequence_count {
            matrices[index * self.sequence_count + index] =
                Some(SparsePosterior::identity(self.lengths[index]));
        }
        let jobs = unordered_pairs(self.sequence_count);
        let work_units = pair_work_units(&jobs, &self.lengths)?
            .checked_mul(self.sequence_count)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "consistency work estimate",
            })?;
        let workers = parallel_worker_count(jobs.len(), work_units);
        let (pairs, peak_bytes) = calculate_consistent_pairs(self, &jobs, cutoff, workers)?;
        for pair in pairs {
            matrices[pair.first * self.sequence_count + pair.second] = Some(pair.forward);
            matrices[pair.second * self.sequence_count + pair.first] = Some(pair.reverse);
        }
        Ok((
            Self {
                sequence_count: self.sequence_count,
                lengths: self.lengths.clone(),
                matrices,
                pair_hmm_peak_bytes: self.pair_hmm_peak_bytes,
                consistency_peak_bytes: self.consistency_peak_bytes,
            },
            peak_bytes,
        ))
    }

    pub(crate) fn get(&self, first: usize, second: usize) -> &SparsePosterior {
        self.matrices[first * self.sequence_count + second]
            .as_ref()
            .expect("complete posterior library")
    }

    pub(crate) const fn sequence_count(&self) -> usize {
        self.sequence_count
    }

    #[allow(clippy::cast_possible_truncation)]
    pub(crate) fn checksum(&self) -> u64 {
        let mut hash = hash_u64(FNV_OFFSET, self.sequence_count as u64);
        for first in 0..self.sequence_count {
            for second in first..self.sequence_count {
                hash = hash_u64(hash, first as u64);
                hash = hash_u64(hash, second as u64);
                hash = hash_u64(hash, self.get(first, second).checksum());
            }
        }
        hash
    }

    pub(crate) fn entry_count(&self) -> usize {
        let diagonal = self.lengths.iter().sum::<usize>();
        let off_diagonal = (0..self.sequence_count)
            .flat_map(|first| {
                (first + 1..self.sequence_count)
                    .map(move |second| self.get(first, second).entry_count())
            })
            .sum::<usize>();
        diagonal + 2 * off_diagonal
    }

    pub(crate) fn storage_bytes(&self) -> usize {
        self.matrices
            .iter()
            .filter_map(Option::as_ref)
            .map(SparsePosterior::storage_bytes)
            .sum::<usize>()
            + self.lengths.len() * size_of::<usize>()
            + self.matrices.len() * size_of::<Option<SparsePosterior>>()
    }

    pub(crate) const fn pair_hmm_peak_bytes(&self) -> usize {
        self.pair_hmm_peak_bytes
    }

    pub(crate) const fn consistency_peak_bytes(&self) -> usize {
        self.consistency_peak_bytes
    }
}

fn unordered_pairs(sequence_count: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for first in 0..sequence_count {
        for second in first + 1..sequence_count {
            pairs.push((first, second));
        }
    }
    pairs
}

fn pair_work_units(jobs: &[(usize, usize)], lengths: &[usize]) -> Result<usize, Muscle5Error> {
    jobs.iter().try_fold(0_usize, |total, (first, second)| {
        let pair = lengths[*first].checked_mul(lengths[*second]).ok_or(
            Muscle5Error::ArithmeticOverflow {
                stage: "pair-posterior work estimate",
            },
        )?;
        total
            .checked_add(pair)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "pair-posterior work estimate",
            })
    })
}

fn parallel_worker_count(job_count: usize, work_units: usize) -> usize {
    if job_count < 2 || work_units < PARALLEL_WORK_THRESHOLD {
        return 1;
    }
    std::thread::available_parallelism()
        .map_or(1, usize::from)
        .min(MAX_PARALLEL_WORKERS)
        .min(job_count)
}

/// The pair-HMM scratch estimate of computing these pairs, in this order, with
/// `workers` workers: [`calculate_pair_posteriors`]'s estimate for the same jobs.
fn scheduled_scratch_bytes(
    calculations: &[IndexedPairPosterior],
    workers: usize,
) -> Result<usize, Muscle5Error> {
    let peak = |chunk: &[IndexedPairPosterior]| {
        chunk
            .iter()
            .map(|indexed| indexed.calculation.peak_scratch_bytes)
            .max()
            .unwrap_or(0)
    };
    if workers == 1 {
        return Ok(peak(calculations));
    }
    calculations
        .chunks(calculations.len().div_ceil(workers))
        .try_fold(0_usize, |total, chunk| {
            total
                .checked_add(peak(chunk))
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "parallel pair-HMM scratch estimate",
                })
        })
}

fn calculate_pair_posteriors<A: Muscle5Alphabet>(
    jobs: &[(usize, usize)],
    records: &[SequenceRecord<A>],
    hmm: &PreparedHmm<A>,
    cutoff: f64,
    workers: usize,
) -> Result<(Vec<IndexedPairPosterior>, usize), Muscle5Error> {
    if workers == 1 {
        let mut calculations = Vec::with_capacity(jobs.len());
        let mut peak_scratch_bytes = 0_usize;
        let mut workspace = PairHmmWorkspace::default();
        for &(first, second) in jobs {
            let calculation = pair_posterior_with_workspace(
                records[first].sequence(),
                records[second].sequence(),
                hmm,
                cutoff,
                &mut workspace,
            )?;
            peak_scratch_bytes = peak_scratch_bytes.max(calculation.peak_scratch_bytes);
            calculations.push(IndexedPairPosterior {
                first,
                second,
                calculation,
            });
        }
        return Ok((calculations, peak_scratch_bytes));
    }

    let chunk_size = jobs.len().div_ceil(workers);
    std::thread::scope(|scope| {
        let handles = jobs
            .chunks(chunk_size)
            .map(|chunk| {
                scope.spawn(move || {
                    let mut calculations = Vec::with_capacity(chunk.len());
                    let mut peak_scratch_bytes = 0_usize;
                    let mut workspace = PairHmmWorkspace::default();
                    for &(first, second) in chunk {
                        let calculation = pair_posterior_with_workspace(
                            records[first].sequence(),
                            records[second].sequence(),
                            hmm,
                            cutoff,
                            &mut workspace,
                        )?;
                        peak_scratch_bytes = peak_scratch_bytes.max(calculation.peak_scratch_bytes);
                        calculations.push(IndexedPairPosterior {
                            first,
                            second,
                            calculation,
                        });
                    }
                    Ok::<_, Muscle5Error>((calculations, peak_scratch_bytes))
                })
            })
            .collect::<Vec<_>>();
        let mut calculations = Vec::with_capacity(jobs.len());
        let mut concurrent_scratch_bytes = 0_usize;
        for handle in handles {
            let (mut chunk, chunk_peak) =
                handle
                    .join()
                    .map_err(|_| Muscle5Error::InternalInvariant {
                        message: "parallel pair-HMM worker panicked",
                    })??;
            concurrent_scratch_bytes = concurrent_scratch_bytes.checked_add(chunk_peak).ok_or(
                Muscle5Error::ArithmeticOverflow {
                    stage: "parallel pair-HMM scratch estimate",
                },
            )?;
            calculations.append(&mut chunk);
        }
        Ok((calculations, concurrent_scratch_bytes))
    })
}

fn calculate_consistent_pair(
    library: &PosteriorLibrary,
    first: usize,
    second: usize,
    cutoff: f64,
    accumulator: &mut Vec<f64>,
) -> Result<ConsistentPair, Muscle5Error> {
    let rows = library.lengths[first];
    let columns = library.lengths[second];
    let dense_len = rows
        .checked_mul(columns)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "consistency accumulator dimensions",
        })?;
    accumulator.resize(dense_len, 0.0);
    accumulator[..dense_len].fill(0.0);
    for intermediate in 0..library.sequence_count {
        let left = library.get(first, intermediate);
        let right = library.get(intermediate, second);
        for row in 0..rows {
            for left_entry in left.row(row) {
                for right_entry in right.row(left_entry.column) {
                    accumulator[row * columns + right_entry.column] +=
                        left_entry.probability * right_entry.probability;
                }
            }
        }
    }
    let denominator = f64::from(u32::try_from(library.sequence_count).map_err(|_| {
        Muscle5Error::ArithmeticOverflow {
            stage: "consistency sequence-count conversion",
        }
    })?);
    for probability in accumulator.iter_mut() {
        *probability /= denominator;
    }
    let accumulator_bytes = accumulator.len().checked_mul(size_of::<f64>()).ok_or(
        Muscle5Error::ArithmeticOverflow {
            stage: "consistency accumulator bytes",
        },
    )?;
    let forward = SparsePosterior::from_dense(rows, columns, accumulator, cutoff)?;
    let reverse = forward.transpose()?;
    Ok(ConsistentPair {
        first,
        second,
        forward,
        reverse,
        accumulator_bytes,
    })
}

fn calculate_consistent_pairs(
    library: &PosteriorLibrary,
    jobs: &[(usize, usize)],
    cutoff: f64,
    workers: usize,
) -> Result<(Vec<ConsistentPair>, usize), Muscle5Error> {
    if workers == 1 {
        let mut pairs = Vec::with_capacity(jobs.len());
        let mut peak_bytes = 0_usize;
        let mut accumulator = Vec::new();
        for &(first, second) in jobs {
            let pair = calculate_consistent_pair(library, first, second, cutoff, &mut accumulator)?;
            peak_bytes = peak_bytes.max(pair.accumulator_bytes);
            pairs.push(pair);
        }
        return Ok((pairs, peak_bytes));
    }

    let chunk_size = jobs.len().div_ceil(workers);
    std::thread::scope(|scope| {
        let handles = jobs
            .chunks(chunk_size)
            .map(|chunk| {
                scope.spawn(move || {
                    let mut pairs = Vec::with_capacity(chunk.len());
                    let mut peak_bytes = 0_usize;
                    let mut accumulator = Vec::new();
                    for &(first, second) in chunk {
                        let pair = calculate_consistent_pair(
                            library,
                            first,
                            second,
                            cutoff,
                            &mut accumulator,
                        )?;
                        peak_bytes = peak_bytes.max(pair.accumulator_bytes);
                        pairs.push(pair);
                    }
                    Ok::<_, Muscle5Error>((pairs, peak_bytes))
                })
            })
            .collect::<Vec<_>>();
        let mut pairs = Vec::with_capacity(jobs.len());
        let mut concurrent_accumulator_bytes = 0_usize;
        for handle in handles {
            let (mut chunk, chunk_peak) =
                handle
                    .join()
                    .map_err(|_| Muscle5Error::InternalInvariant {
                        message: "parallel consistency worker panicked",
                    })??;
            concurrent_accumulator_bytes = concurrent_accumulator_bytes
                .checked_add(chunk_peak)
                .ok_or(Muscle5Error::ArithmeticOverflow {
                    stage: "parallel consistency scratch estimate",
                })?;
            pairs.append(&mut chunk);
        }
        Ok((pairs, concurrent_accumulator_bytes))
    })
}

impl PosteriorProvider for PosteriorLibrary {
    fn posterior(&self, first: usize, second: usize) -> Option<&SparsePosterior> {
        if first >= self.sequence_count || second >= self.sequence_count {
            return None;
        }
        self.matrices[first * self.sequence_count + second].as_ref()
    }

    fn sequence_length(&self, sequence_index: usize) -> Option<usize> {
        self.lengths.get(sequence_index).copied()
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum PairOperation {
    Match,
    ConsumeFirst,
    ConsumeSecond,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PairExpectedAlignment {
    operations: Box<[PairOperation]>,
    expected_accuracy: f64,
    expected_error: f64,
}

impl PairExpectedAlignment {
    pub(crate) const fn operations(&self) -> &[PairOperation] {
        &self.operations
    }

    pub(crate) const fn expected_accuracy(&self) -> f64 {
        self.expected_accuracy
    }

    pub(crate) const fn expected_error(&self) -> f64 {
        self.expected_error
    }
}

pub(crate) fn maximum_expected_accuracy(
    posterior: &SparsePosterior,
) -> Result<PairExpectedAlignment, Muscle5Error> {
    let lattice_rows = posterior
        .rows()
        .checked_add(1)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "expected-accuracy row count",
        })?;
    let lattice_columns =
        posterior
            .columns()
            .checked_add(1)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "expected-accuracy column count",
            })?;
    let cells =
        lattice_rows
            .checked_mul(lattice_columns)
            .ok_or(Muscle5Error::ArithmeticOverflow {
                stage: "expected-accuracy lattice",
            })?;
    let mut scores = vec![0.0_f64; cells];
    let mut trace = vec![0_u8; cells];
    for row in 1..lattice_rows {
        trace[row * lattice_columns] = 1;
    }
    trace[1..lattice_columns].fill(2);
    for row in 1..lattice_rows {
        for column in 1..lattice_columns {
            let diagonal = scores[(row - 1) * lattice_columns + column - 1]
                + posterior.get(row - 1, column - 1);
            let consume_first = scores[(row - 1) * lattice_columns + column];
            let consume_second = scores[row * lattice_columns + column - 1];
            let index = row * lattice_columns + column;
            if diagonal >= consume_first && diagonal >= consume_second {
                scores[index] = diagonal;
                trace[index] = 0;
            } else if consume_first >= consume_second {
                scores[index] = consume_first;
                trace[index] = 1;
            } else {
                scores[index] = consume_second;
                trace[index] = 2;
            }
        }
    }

    let mut row = posterior.rows();
    let mut column = posterior.columns();
    let mut reverse = Vec::with_capacity(row + column);
    while row > 0 || column > 0 {
        let direction = trace.get(row * lattice_columns + column).copied().ok_or(
            Muscle5Error::InternalInvariant {
                message: "expected-accuracy traceback left its lattice",
            },
        )?;
        match direction {
            0 => {
                reverse.push(PairOperation::Match);
                row -= 1;
                column -= 1;
            }
            1 => {
                reverse.push(PairOperation::ConsumeFirst);
                row -= 1;
            }
            2 => {
                reverse.push(PairOperation::ConsumeSecond);
                column -= 1;
            }
            _ => {
                return Err(Muscle5Error::InternalInvariant {
                    message: "invalid expected-accuracy traceback state",
                });
            }
        }
    }
    reverse.reverse();
    let expected_accuracy = scores[posterior.rows() * lattice_columns + posterior.columns()];
    let expected_error = if reverse.is_empty() {
        0.0
    } else {
        let column_count = f64::from(u32::try_from(reverse.len()).map_err(|_| {
            Muscle5Error::ArithmeticOverflow {
                stage: "expected-error column-count conversion",
            }
        })?);
        (1.0 - expected_accuracy / column_count).clamp(0.0, 1.0)
    };
    Ok(PairExpectedAlignment {
        operations: reverse.into_boxed_slice(),
        expected_accuracy,
        expected_error,
    })
}

pub(crate) fn expected_error_matrix(
    library: &PosteriorLibrary,
) -> Result<DistanceMatrix, Muscle5Error> {
    let count = library.sequence_count();
    let cells = count
        .checked_mul(count)
        .ok_or(Muscle5Error::ArithmeticOverflow {
            stage: "expected-error matrix dimensions",
        })?;
    let mut values = vec![0.0_f64; cells];
    for first in 0..count {
        for second in first + 1..count {
            let alignment = maximum_expected_accuracy(library.get(first, second))?;
            values[first * count + second] = alignment.expected_error();
            values[second * count + first] = alignment.expected_error();
        }
    }
    DistanceMatrix::try_from_flat(count, values).map_err(Muscle5Error::from)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::float_cmp)]

    use align_core::{Dna, DnaSequence};
    use align_msa::SequenceRecord;

    use crate::muscle5_hmm::HmmParameters;
    use crate::muscle5_rng::SplitMix64;

    use super::{
        PROBABILITY_FAST_PATH_MIN_LOG, PairHmmWorkspace, PairOperation, PairPosteriorCache,
        PosteriorLibrary, SparsePosterior, calculate_consistent_pairs, calculate_pair_posteriors,
        concrete_path_log_probability, maximum_expected_accuracy, pair_posterior,
        pair_posterior_log, pair_posterior_probability, unordered_pairs,
    };

    fn seeded_records(count: usize, length: usize, seed: u64) -> Vec<SequenceRecord<Dna>> {
        let mut rng = SplitMix64::new(seed);
        (0..count)
            .map(|index| {
                let sequence = (0..length)
                    .map(|_| char::from(b"ACGT"[rng.index(4)]))
                    .collect::<String>();
                SequenceRecord::new(
                    format!("s{index}"),
                    DnaSequence::from_ascii(&sequence).unwrap(),
                )
                .unwrap()
            })
            .collect()
    }

    #[test]
    fn a_library_built_from_a_partly_filled_cache_equals_a_fresh_one() {
        // Six sequences of 120 bases: 15 pairs and 216,000 work units, enough
        // for the parallel path when more than one CPU is available.
        let records = seeded_records(6, 120, 11);
        let hmm = HmmParameters::published::<Dna>(0, 0.25)
            .unwrap()
            .prepare::<Dna>();
        let fresh = PosteriorLibrary::build(&records, &hmm, 0.01, None).unwrap();
        let mut cache = PairPosteriorCache::new(&hmm, 0.01);
        for (first, second) in [(0, 1), (2, 4), (3, 5), (1, 0)] {
            cache
                .pair_posterior(
                    records[first].sequence(),
                    records[second].sequence(),
                    &hmm,
                    0.01,
                )
                .unwrap();
        }
        let reused = PosteriorLibrary::build(&records, &hmm, 0.01, Some(&cache)).unwrap();
        // The library asks for pairs in increasing order, so (1, 0) is not used.
        assert_eq!(cache.hits(), 3);
        assert_eq!(format!("{reused:?}"), format!("{fresh:?}"));
        assert_eq!(reused.pair_hmm_peak_bytes(), fresh.pair_hmm_peak_bytes());
    }

    #[test]
    fn a_cache_answers_only_for_its_own_hmm_cutoff_and_pair_order() {
        let records = seeded_records(2, 30, 5);
        let (first, second) = (records[0].sequence(), records[1].sequence());
        let hmm = HmmParameters::published::<Dna>(0, 0.25)
            .unwrap()
            .prepare::<Dna>();
        let other_hmm = HmmParameters::published::<Dna>(1, 0.25)
            .unwrap()
            .prepare::<Dna>();
        assert_ne!(hmm, other_hmm);
        let mut cache = PairPosteriorCache::new(&hmm, 0.01);
        let stored = cache.pair_posterior(first, second, &hmm, 0.01).unwrap();
        assert!(cache.stored(first, second, &hmm, 0.01).is_some());
        assert!(cache.stored(first, second, &other_hmm, 0.01).is_none());
        assert!(cache.stored(first, second, &hmm, 0.02).is_none());
        assert!(cache.stored(second, first, &hmm, 0.01).is_none());
        // A request the cache does not serve is computed, not answered from it.
        let other = cache
            .pair_posterior(first, second, &other_hmm, 0.01)
            .unwrap();
        assert_eq!(
            other.posterior,
            pair_posterior(first, second, &other_hmm, 0.01)
                .unwrap()
                .posterior
        );
        assert_eq!(
            stored.posterior,
            pair_posterior(first, second, &hmm, 0.01).unwrap().posterior
        );
        assert_eq!(cache.hits(), 1);
    }

    #[test]
    fn sparse_transpose_preserves_coordinates_and_checksum_is_stable() {
        let matrix =
            SparsePosterior::from_dense(2, 3, &[0.0, 0.2, 0.3, 0.4, 0.0, 0.5], 0.25).unwrap();
        assert_eq!(matrix.entry_count(), 3);
        assert_eq!(matrix.get(0, 1), 0.0);
        assert_eq!(matrix.get(0, 2), 0.3);
        let transpose = matrix.transpose().unwrap();
        assert_eq!(transpose.rows(), 3);
        assert_eq!(transpose.columns(), 2);
        assert_eq!(transpose.get(2, 0), 0.3);
        assert_eq!(transpose.get(0, 1), 0.4);
        assert_eq!(matrix.checksum(), matrix.clone().checksum());
    }

    #[test]
    fn one_letter_pair_has_unit_match_posterior() {
        let first = DnaSequence::from_ascii("A").unwrap();
        let second = DnaSequence::from_ascii("A").unwrap();
        let hmm = HmmParameters::published::<Dna>(0, 0.25)
            .unwrap()
            .prepare::<Dna>();
        let result = pair_posterior(&first, &second, &hmm, 0.0).unwrap();
        assert!((result.posterior.get(0, 0) - 1.0).abs() < 1.0e-12);
        assert!(result.log_partition.is_finite());
    }

    #[test]
    fn two_by_one_posterior_matches_independent_path_enumeration() {
        let first = DnaSequence::from_ascii("AA").unwrap();
        let second = DnaSequence::from_ascii("A").unwrap();
        let hmm = HmmParameters::published::<Dna>(0, 0.25)
            .unwrap()
            .prepare::<Dna>();
        let result = pair_posterior(&first, &second, &hmm, 0.0).unwrap();

        let short_then_match = 0.02 * 0.25 * 0.60 * 0.22 * 0.62;
        let long_then_match = 0.18 * 0.25 * 0.10 * 0.22 * 0.62;
        let match_then_short = 0.62 * 0.22 * 0.012 * 0.25 * 0.02;
        let match_then_long = 0.62 * 0.22 * 0.008 * 0.25 * 0.18;
        let total = short_then_match + long_then_match + match_then_short + match_then_long;
        let first_match = (match_then_short + match_then_long) / total;
        let second_match = (short_then_match + long_then_match) / total;
        assert!((result.posterior.get(0, 0) - first_match).abs() < 1.0e-12);
        assert!((result.posterior.get(1, 0) - second_match).abs() < 1.0e-12);
        assert!((first_match + second_match - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn empty_nonempty_pair_has_no_match_and_finite_partition() {
        let first = DnaSequence::from_ascii("").unwrap();
        let second = DnaSequence::from_ascii("AAA").unwrap();
        let hmm = HmmParameters::published::<Dna>(0, 0.25)
            .unwrap()
            .prepare::<Dna>();
        let result = pair_posterior(&first, &second, &hmm, 0.0).unwrap();
        assert_eq!(result.posterior.rows(), 0);
        assert_eq!(result.posterior.columns(), 3);
        assert_eq!(result.posterior.entry_count(), 0);
        assert!(result.log_partition.is_finite());
    }

    #[test]
    fn guarded_probability_domain_agrees_with_independent_log_domain() {
        let first = DnaSequence::from_ascii("ACGTTACG").unwrap();
        let second = DnaSequence::from_ascii("AGTCAC").unwrap();
        let hmm = HmmParameters::published::<Dna>(7, 0.25)
            .unwrap()
            .prepare::<Dna>();
        assert!(
            concrete_path_log_probability(&first, &second, &hmm) >= PROBABILITY_FAST_PATH_MIN_LOG
        );
        let probability = pair_posterior_probability(
            &first,
            &second,
            &hmm,
            0.0,
            &mut PairHmmWorkspace::default(),
        )
        .unwrap();
        let logarithmic =
            pair_posterior_log(&first, &second, &hmm, 0.0, &mut PairHmmWorkspace::default())
                .unwrap();
        assert!((probability.log_partition - logarithmic.log_partition).abs() < 1.0e-12);
        for row in 0..first.len() {
            for column in 0..second.len() {
                assert!(
                    (probability.posterior.get(row, column)
                        - logarithmic.posterior.get(row, column))
                    .abs()
                        < 2.0e-12
                );
            }
        }

        let divergent = DnaSequence::from_ascii("A".repeat(140)).unwrap();
        let other = DnaSequence::from_ascii("C".repeat(140)).unwrap();
        assert!(
            concrete_path_log_probability(&divergent, &other, &hmm) < PROBABILITY_FAST_PATH_MIN_LOG
        );
    }

    #[test]
    fn parallel_collectors_preserve_sequential_pair_order_and_bits() {
        let records = [
            SequenceRecord::new("a", DnaSequence::from_ascii("ACGTACGT").unwrap()).unwrap(),
            SequenceRecord::new("b", DnaSequence::from_ascii("ACGTCGT").unwrap()).unwrap(),
            SequenceRecord::new("c", DnaSequence::from_ascii("AGGTACGT").unwrap()).unwrap(),
            SequenceRecord::new("d", DnaSequence::from_ascii("TCGTACGA").unwrap()).unwrap(),
        ];
        let hmm = HmmParameters::published::<Dna>(0, 0.25)
            .unwrap()
            .prepare::<Dna>();
        let jobs = unordered_pairs(records.len());
        let (sequential, sequential_peak) =
            calculate_pair_posteriors(&jobs, &records, &hmm, 0.01, 1).unwrap();
        let (parallel, parallel_peak) =
            calculate_pair_posteriors(&jobs, &records, &hmm, 0.01, 3).unwrap();
        assert_eq!(sequential.len(), parallel.len());
        for (left, right) in sequential.iter().zip(&parallel) {
            assert_eq!((left.first, left.second), (right.first, right.second));
            assert_eq!(left.calculation.posterior, right.calculation.posterior);
            assert_eq!(
                left.calculation.log_partition.to_bits(),
                right.calculation.log_partition.to_bits()
            );
        }
        assert!(parallel_peak >= sequential_peak);

        let library = PosteriorLibrary::build(&records, &hmm, 0.01, None).unwrap();
        let (sequential, sequential_peak) =
            calculate_consistent_pairs(&library, &jobs, 0.01, 1).unwrap();
        let (parallel, parallel_peak) =
            calculate_consistent_pairs(&library, &jobs, 0.01, 3).unwrap();
        for (left, right) in sequential.iter().zip(&parallel) {
            assert_eq!((left.first, left.second), (right.first, right.second));
            assert_eq!(left.forward, right.forward);
            assert_eq!(left.reverse, right.reverse);
        }
        assert!(parallel_peak >= sequential_peak);
    }

    #[test]
    fn expected_accuracy_traceback_uses_documented_ties() {
        let posterior = SparsePosterior::from_dense(2, 2, &[0.9, 0.1, 0.2, 0.8], 0.0).unwrap();
        let alignment = maximum_expected_accuracy(&posterior).unwrap();
        assert_eq!(
            alignment.operations(),
            &[PairOperation::Match, PairOperation::Match]
        );
        assert!((alignment.expected_accuracy() - 1.7).abs() < 1.0e-12);
        assert!((alignment.expected_error() - 0.15).abs() < 1.0e-12);

        let zero = SparsePosterior::from_dense(1, 1, &[0.0], 0.0).unwrap();
        assert_eq!(
            maximum_expected_accuracy(&zero).unwrap().operations(),
            &[PairOperation::Match]
        );
    }

    #[test]
    fn one_consistency_round_matches_direct_scalar_product() {
        let identity = SparsePosterior::identity(1);
        let p01 = SparsePosterior::from_dense(1, 1, &[0.2], 0.0).unwrap();
        let p02 = SparsePosterior::from_dense(1, 1, &[0.4], 0.0).unwrap();
        let p12 = SparsePosterior::from_dense(1, 1, &[0.7], 0.0).unwrap();
        let matrices = vec![
            Some(identity.clone()),
            Some(p01.clone()),
            Some(p02.clone()),
            Some(p01),
            Some(identity.clone()),
            Some(p12.clone()),
            Some(p02),
            Some(p12),
            Some(identity),
        ];
        let library = PosteriorLibrary {
            sequence_count: 3,
            lengths: vec![1, 1, 1].into_boxed_slice(),
            matrices,
            pair_hmm_peak_bytes: 0,
            consistency_peak_bytes: 0,
        };
        let (consistent, _) = library.consistency_round(0.0).unwrap();
        let expected = (0.2 + 0.2 + 0.4 * 0.7) / 3.0;
        assert!((consistent.get(0, 1).get(0, 0) - expected).abs() < 1.0e-12);
    }

    #[test]
    fn library_builds_transposed_pairs_and_transforms_deterministically() {
        let records = [
            SequenceRecord::new("a", DnaSequence::from_ascii("AC").unwrap()).unwrap(),
            SequenceRecord::new("b", DnaSequence::from_ascii("A").unwrap()).unwrap(),
            SequenceRecord::new("c", DnaSequence::from_ascii("AG").unwrap()).unwrap(),
        ];
        let hmm = HmmParameters::published::<Dna>(0, 0.25)
            .unwrap()
            .prepare::<Dna>();
        let mut first = PosteriorLibrary::build(&records, &hmm, 0.0, None).unwrap();
        let mut second = first.clone();
        assert_eq!(first.get(0, 1).transpose().unwrap(), *first.get(1, 0));
        first.transform_consistency(2, 0.0).unwrap();
        second.transform_consistency(2, 0.0).unwrap();
        assert_eq!(first.checksum(), second.checksum());
        assert!(first.storage_bytes() > 0);
    }
}
