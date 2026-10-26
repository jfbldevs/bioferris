use core::cmp::Reverse;
use core::hash::{BuildHasherDefault, Hash, Hasher};
use core::mem::size_of;
use core::ops::Range;
use std::collections::{BinaryHeap, HashMap, HashSet, hash_map::Entry};
use std::time::Instant;

use align_core::{Alphabet, PhaseTimingsBuilder, Profiled, Score, Sequence};
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use pulp::Simd;

use crate::affine::{
    AffinePoaAlignmentMemoryEstimate, GlobalAffinePoa, verify_global_affine_operations_validated,
};
use crate::superbubble::PoastaGraphIndex;
use crate::{
    AffinePoaVerificationReport, PoaAlignedPair, PoaError, PoaGraph, PoaNodeId, PoaOperation,
    PoastaAffineCosts,
};

/// Frozen complete profile reconstructed from van Dijk et al. (2025).
pub const POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_V1: &str =
    "POASTA-2025-global-affine-ASTAR-DFS-superbubble-v1";

const POASTA_ASTAR_VALIDATION_PROFILE: &str = "POASTA-2025-validation-ASTAR-v1";
const POASTA_ASTAR_DFS_VALIDATION_PROFILE: &str = "POASTA-2025-validation-ASTAR-DFS-v1";
const POASTA_ASTAR_PRUNING_VALIDATION_PROFILE: &str = "POASTA-2025-validation-ASTAR-superbubble-v1";
const POASTA_EXACT_DENSE_PROFILE: &str = "POASTA-exact-dense-global-affine-v1";
const POASTA_EXACT_AUTO_PROFILE: &str = "POASTA-exact-auto-global-affine-v1";

type Cost = u64;

/// Exact execution strategy used by the POASTA interface.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PoastaBackend {
    /// Published sparse A* search with the selected accelerations.
    Sparse,
    /// Contiguous exact affine dynamic programming over the complete DAG matrix.
    Dense,
    /// Uses the validated optimized policy while preserving exact results.
    ///
    /// The current policy starts with sparse published search plus safe affine-gap
    /// dominance, then activates a stronger relaxed linear-gap potential only when
    /// the observed search exceeds its lightweight probe budget. The dense engine
    /// remains available explicitly but lost on every retained crossover benchmark.
    Auto,
}

/// Exact feature-isolation controls for POASTA validation and benchmarking.
///
/// Every combination remains an exact sparse A* search. Only the default with
/// both features enabled claims the complete published profile identifier.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PoastaSearchOptions {
    exact_match_dfs: bool,
    superbubble_pruning: bool,
    backend: PoastaBackend,
}

impl PoastaSearchOptions {
    /// Enables every characteristic acceleration in the published algorithm.
    #[must_use]
    pub const fn published() -> Self {
        Self {
            exact_match_dfs: true,
            superbubble_pruning: true,
            backend: PoastaBackend::Sparse,
        }
    }

    /// Enables the validated exact optimized-execution policy.
    #[must_use]
    pub const fn optimized() -> Self {
        Self {
            exact_match_dfs: true,
            superbubble_pruning: false,
            backend: PoastaBackend::Auto,
        }
    }

    /// Selects the exact dense affine recurrence explicitly.
    #[must_use]
    pub const fn dense() -> Self {
        Self {
            exact_match_dfs: true,
            superbubble_pruning: true,
            backend: PoastaBackend::Dense,
        }
    }

    /// Creates an ordinary exact A* baseline for differential validation.
    #[must_use]
    pub const fn astar_baseline() -> Self {
        Self {
            exact_match_dfs: false,
            superbubble_pruning: false,
            backend: PoastaBackend::Sparse,
        }
    }

    /// Returns a copy with exact-match DFS explicitly enabled or disabled.
    #[must_use]
    pub const fn with_exact_match_dfs(mut self, enabled: bool) -> Self {
        self.exact_match_dfs = enabled;
        self
    }

    /// Returns a copy with superbubble pruning explicitly enabled or disabled.
    #[must_use]
    pub const fn with_superbubble_pruning(mut self, enabled: bool) -> Self {
        self.superbubble_pruning = enabled;
        self
    }

    /// Returns a copy using the requested exact execution strategy.
    #[must_use]
    pub const fn with_backend(mut self, backend: PoastaBackend) -> Self {
        self.backend = backend;
        self
    }

    /// Returns whether exact zero-cost matches use the DFS frontier.
    #[must_use]
    pub const fn exact_match_dfs(self) -> bool {
        self.exact_match_dfs
    }

    /// Returns whether superbubble-informed dominance pruning is active.
    #[must_use]
    pub const fn superbubble_pruning(self) -> bool {
        self.superbubble_pruning
    }

    /// Returns the requested execution strategy.
    #[must_use]
    pub const fn backend(self) -> PoastaBackend {
        self.backend
    }
}

impl Default for PoastaSearchOptions {
    fn default() -> Self {
        Self::published()
    }
}

/// Detailed-only work counters for compact exact-affine certificate rows.
///
/// These counters are collected only by the opt-in detailed profiler. Ordinary
/// and coarse-profiled alignment paths specialize the collection branch away.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PoastaSearchWorkCounters {
    source_path_certificates: usize,
    early_diagonal_certificates: usize,
    relaxed_linear_table_attempts: usize,
    relaxed_linear_table_diagonal_certificates: usize,
    relaxed_linear_table_affine_opening_certificates: usize,
    relaxed_linear_table_misses: usize,
    relaxed_linear_table_certificate_cells: usize,
    relaxed_linear_table_miss_cells: usize,
    relaxed_linear_rows: usize,
    relaxed_linear_insertion_cells: usize,
    relaxed_linear_insertion_avx2_cells: usize,
    relaxed_linear_zero_node_successor_rows: usize,
    relaxed_linear_single_node_successor_rows: usize,
    relaxed_linear_multiple_node_successor_rows: usize,
    relaxed_linear_two_node_successor_rows: usize,
    relaxed_linear_three_or_four_node_successor_rows: usize,
    relaxed_linear_five_plus_node_successor_rows: usize,
    relaxed_linear_multiple_node_successor_cells: usize,
    relaxed_linear_multiple_node_successor_cell_visits: usize,
    relaxed_linear_multiple_node_distinct_symbol_cell_visits: usize,
    relaxed_linear_successor_cell_visits: usize,
    relaxed_linear_equivalent_single_transition_rows: usize,
    relaxed_linear_equivalent_single_transition_cells: usize,
    relaxed_linear_equivalent_single_transition_insertion_cells: usize,
    relaxed_linear_equivalent_single_transition_successor_cell_visits: usize,
    relaxed_linear_equivalent_multiple_transition_rows: usize,
    relaxed_linear_equivalent_multiple_transition_cells: usize,
    relaxed_linear_equivalent_multiple_transition_insertion_cells: usize,
    relaxed_linear_equivalent_multiple_transition_successor_cell_visits: usize,
    relaxed_linear_initialized_rows: usize,
    relaxed_linear_reused_rows: usize,
    relaxed_linear_initialized_cells: usize,
    exact_affine_attempts: usize,
    exact_affine_bypass_attempts: usize,
    exact_affine_bypass_certificates: usize,
    exact_affine_bypass_cells: usize,
    exact_affine_after_relaxed_attempts: usize,
    exact_affine_after_relaxed_certificates: usize,
    exact_affine_after_relaxed_cells: usize,
    exact_affine_filled_attempts: usize,
    exact_affine_source_cells: usize,
    exact_affine_single_predecessor_rows: usize,
    exact_affine_single_predecessor_cells: usize,
    exact_affine_single_predecessor_fast_cells: usize,
    exact_affine_single_predecessor_avx2_cells: usize,
    exact_affine_two_predecessor_rows: usize,
    exact_affine_two_predecessor_cells: usize,
    exact_affine_two_predecessor_fast_cells: usize,
    exact_affine_two_predecessor_avx2_cells: usize,
    exact_affine_many_predecessor_rows: usize,
    exact_affine_many_predecessor_cells: usize,
    exact_affine_many_predecessor_cell_visits: usize,
    exact_affine_many_predecessor_avx2_cells: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RelaxedLinearTableOutcome {
    DiagonalCertificate,
    AffineOpeningCertificate,
    Miss,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ExactAffineRoute {
    Bypass,
    AfterRelaxed,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct RelaxedLinearRowWork {
    row_cells: usize,
    insertion_cells: usize,
    insertion_avx2_cells: usize,
    node_successors: usize,
    distinct_node_successor_symbols: usize,
    successor_cell_visits: usize,
    equivalent_multiple_transition: bool,
}

impl PoastaSearchWorkCounters {
    /// Returns exact stored-path certificates found before any heuristic table.
    #[must_use]
    pub const fn source_path_certificates(self) -> usize {
        self.source_path_certificates
    }

    /// Returns exact zero diagonal certificates found before a relaxed table.
    #[must_use]
    pub const fn early_diagonal_certificates(self) -> usize {
        self.early_diagonal_certificates
    }

    /// Returns relaxed-linear tables built for progressive queries.
    #[must_use]
    pub const fn relaxed_linear_table_attempts(self) -> usize {
        self.relaxed_linear_table_attempts
    }

    /// Returns relaxed tables that directly certified a diagonal traceback.
    #[must_use]
    pub const fn relaxed_linear_table_diagonal_certificates(self) -> usize {
        self.relaxed_linear_table_diagonal_certificates
    }

    /// Returns relaxed tables that certified mandatory affine openings.
    #[must_use]
    pub const fn relaxed_linear_table_affine_opening_certificates(self) -> usize {
        self.relaxed_linear_table_affine_opening_certificates
    }

    /// Returns relaxed tables that failed to certify the query.
    #[must_use]
    pub const fn relaxed_linear_table_misses(self) -> usize {
        self.relaxed_linear_table_misses
    }

    /// Returns relaxed cells in tables that supplied a certificate.
    #[must_use]
    pub const fn relaxed_linear_table_certificate_cells(self) -> usize {
        self.relaxed_linear_table_certificate_cells
    }

    /// Returns relaxed cells computed before a failed certificate attempt.
    #[must_use]
    pub const fn relaxed_linear_table_miss_cells(self) -> usize {
        self.relaxed_linear_table_miss_cells
    }

    /// Returns relaxed-linear rows whose insertion recurrence was evaluated.
    #[must_use]
    pub const fn relaxed_linear_rows(self) -> usize {
        self.relaxed_linear_rows
    }

    /// Returns cells evaluated by relaxed-linear insertion recurrences.
    #[must_use]
    pub const fn relaxed_linear_insertion_cells(self) -> usize {
        self.relaxed_linear_insertion_cells
    }

    /// Returns relaxed-linear insertion cells finalized by the AVX2 suffix scan.
    #[must_use]
    pub const fn relaxed_linear_insertion_avx2_cells(self) -> usize {
        self.relaxed_linear_insertion_avx2_cells
    }

    /// Returns relaxed-linear rows whose only outgoing transition is the sink.
    #[must_use]
    pub const fn relaxed_linear_zero_node_successor_rows(self) -> usize {
        self.relaxed_linear_zero_node_successor_rows
    }

    /// Returns relaxed-linear rows with exactly one non-sink successor.
    #[must_use]
    pub const fn relaxed_linear_single_node_successor_rows(self) -> usize {
        self.relaxed_linear_single_node_successor_rows
    }

    /// Returns relaxed-linear rows with more than one non-sink successor.
    #[must_use]
    pub const fn relaxed_linear_multiple_node_successor_rows(self) -> usize {
        self.relaxed_linear_multiple_node_successor_rows
    }

    /// Returns relaxed rows with exactly two non-sink successors.
    #[must_use]
    pub const fn relaxed_linear_two_node_successor_rows(self) -> usize {
        self.relaxed_linear_two_node_successor_rows
    }

    /// Returns relaxed rows with three or four non-sink successors.
    #[must_use]
    pub const fn relaxed_linear_three_or_four_node_successor_rows(self) -> usize {
        self.relaxed_linear_three_or_four_node_successor_rows
    }

    /// Returns relaxed rows with at least five non-sink successors.
    #[must_use]
    pub const fn relaxed_linear_five_plus_node_successor_rows(self) -> usize {
        self.relaxed_linear_five_plus_node_successor_rows
    }

    /// Returns active cells in rows with multiple non-sink successors.
    #[must_use]
    pub const fn relaxed_linear_multiple_node_successor_cells(self) -> usize {
        self.relaxed_linear_multiple_node_successor_cells
    }

    /// Returns successor cells visited by multiple-successor rows.
    #[must_use]
    pub const fn relaxed_linear_multiple_node_successor_cell_visits(self) -> usize {
        self.relaxed_linear_multiple_node_successor_cell_visits
    }

    /// Returns active cells weighted by distinct symbols in multiple-successor rows.
    #[must_use]
    pub const fn relaxed_linear_multiple_node_distinct_symbol_cell_visits(self) -> usize {
        self.relaxed_linear_multiple_node_distinct_symbol_cell_visits
    }

    /// Returns cells visited while merging non-sink successor rows.
    #[must_use]
    pub const fn relaxed_linear_successor_cell_visits(self) -> usize {
        self.relaxed_linear_successor_cell_visits
    }

    /// Returns rows sharing their sole transition with an earlier equivalent row.
    #[must_use]
    pub const fn relaxed_linear_equivalent_single_transition_rows(self) -> usize {
        self.relaxed_linear_equivalent_single_transition_rows
    }

    /// Returns active cells in equivalent single-transition rows.
    #[must_use]
    pub const fn relaxed_linear_equivalent_single_transition_cells(self) -> usize {
        self.relaxed_linear_equivalent_single_transition_cells
    }

    /// Returns insertion cells replaceable in equivalent single-transition rows.
    #[must_use]
    pub const fn relaxed_linear_equivalent_single_transition_insertion_cells(self) -> usize {
        self.relaxed_linear_equivalent_single_transition_insertion_cells
    }

    /// Returns successor-cell visits replaceable in equivalent single-transition rows.
    #[must_use]
    pub const fn relaxed_linear_equivalent_single_transition_successor_cell_visits(self) -> usize {
        self.relaxed_linear_equivalent_single_transition_successor_cell_visits
    }

    /// Returns multiple-transition rows matching an earlier complete transition set and bounds.
    #[must_use]
    pub const fn relaxed_linear_equivalent_multiple_transition_rows(self) -> usize {
        self.relaxed_linear_equivalent_multiple_transition_rows
    }

    /// Returns active cells in exact-equivalent multiple-transition rows.
    #[must_use]
    pub const fn relaxed_linear_equivalent_multiple_transition_cells(self) -> usize {
        self.relaxed_linear_equivalent_multiple_transition_cells
    }

    /// Returns insertion cells in exact-equivalent multiple-transition rows.
    #[must_use]
    pub const fn relaxed_linear_equivalent_multiple_transition_insertion_cells(self) -> usize {
        self.relaxed_linear_equivalent_multiple_transition_insertion_cells
    }

    /// Returns successor visits in exact-equivalent multiple-transition rows.
    #[must_use]
    pub const fn relaxed_linear_equivalent_multiple_transition_successor_cell_visits(
        self,
    ) -> usize {
        self.relaxed_linear_equivalent_multiple_transition_successor_cell_visits
    }

    /// Returns relaxed rows whose length-bound seed was materialized.
    #[must_use]
    pub const fn relaxed_linear_initialized_rows(self) -> usize {
        self.relaxed_linear_initialized_rows
    }

    /// Returns relaxed rows whose unchanged length-bound seed was reused.
    #[must_use]
    pub const fn relaxed_linear_reused_rows(self) -> usize {
        self.relaxed_linear_reused_rows
    }

    /// Returns cells written while materializing relaxed length-bound seeds.
    #[must_use]
    pub const fn relaxed_linear_initialized_cells(self) -> usize {
        self.relaxed_linear_initialized_cells
    }

    /// Returns all exact-affine certificate attempts, including range declines.
    #[must_use]
    pub const fn exact_affine_attempts(self) -> usize {
        self.exact_affine_attempts
    }

    /// Returns exact-affine attempts selected before building a relaxed table.
    #[must_use]
    pub const fn exact_affine_bypass_attempts(self) -> usize {
        self.exact_affine_bypass_attempts
    }

    /// Returns direct exact-affine attempts that supplied a certificate.
    #[must_use]
    pub const fn exact_affine_bypass_certificates(self) -> usize {
        self.exact_affine_bypass_certificates
    }

    /// Returns exact-affine cells evaluated by direct bypass attempts.
    #[must_use]
    pub const fn exact_affine_bypass_cells(self) -> usize {
        self.exact_affine_bypass_cells
    }

    /// Returns exact-affine attempts made after a relaxed-table miss.
    #[must_use]
    pub const fn exact_affine_after_relaxed_attempts(self) -> usize {
        self.exact_affine_after_relaxed_attempts
    }

    /// Returns post-relaxed exact-affine attempts that supplied a certificate.
    #[must_use]
    pub const fn exact_affine_after_relaxed_certificates(self) -> usize {
        self.exact_affine_after_relaxed_certificates
    }

    /// Returns exact-affine cells evaluated after relaxed-table misses.
    #[must_use]
    pub const fn exact_affine_after_relaxed_cells(self) -> usize {
        self.exact_affine_after_relaxed_cells
    }

    /// Returns attempts that allocated and filled compact affine tables.
    #[must_use]
    pub const fn exact_affine_filled_attempts(self) -> usize {
        self.exact_affine_filled_attempts
    }

    /// Returns virtual-source row cells initialized by filled attempts.
    #[must_use]
    pub const fn exact_affine_source_cells(self) -> usize {
        self.exact_affine_source_cells
    }

    /// Returns biological rows using the single-predecessor kernel.
    #[must_use]
    pub const fn exact_affine_single_predecessor_rows(self) -> usize {
        self.exact_affine_single_predecessor_rows
    }

    /// Returns banded cells assigned to the single-predecessor kernel.
    #[must_use]
    pub const fn exact_affine_single_predecessor_cells(self) -> usize {
        self.exact_affine_single_predecessor_cells
    }

    /// Returns single-predecessor cells handled by its contiguous fast region.
    #[must_use]
    pub const fn exact_affine_single_predecessor_fast_cells(self) -> usize {
        self.exact_affine_single_predecessor_fast_cells
    }

    /// Returns single-predecessor cells finalized by the runtime-dispatched AVX2 scan.
    #[must_use]
    pub const fn exact_affine_single_predecessor_avx2_cells(self) -> usize {
        self.exact_affine_single_predecessor_avx2_cells
    }

    /// Returns biological rows using the two-predecessor kernel.
    #[must_use]
    pub const fn exact_affine_two_predecessor_rows(self) -> usize {
        self.exact_affine_two_predecessor_rows
    }

    /// Returns banded cells assigned to the two-predecessor kernel.
    #[must_use]
    pub const fn exact_affine_two_predecessor_cells(self) -> usize {
        self.exact_affine_two_predecessor_cells
    }

    /// Returns two-predecessor cells handled by its contiguous fast region.
    #[must_use]
    pub const fn exact_affine_two_predecessor_fast_cells(self) -> usize {
        self.exact_affine_two_predecessor_fast_cells
    }

    /// Returns two-predecessor cells finalized by the runtime-dispatched AVX2 scan.
    #[must_use]
    pub const fn exact_affine_two_predecessor_avx2_cells(self) -> usize {
        self.exact_affine_two_predecessor_avx2_cells
    }

    /// Returns biological rows using the many-predecessor merge kernel.
    #[must_use]
    pub const fn exact_affine_many_predecessor_rows(self) -> usize {
        self.exact_affine_many_predecessor_rows
    }

    /// Returns banded cells assigned to the many-predecessor kernel.
    #[must_use]
    pub const fn exact_affine_many_predecessor_cells(self) -> usize {
        self.exact_affine_many_predecessor_cells
    }

    /// Returns overlapping predecessor cells visited by many-row merges.
    #[must_use]
    pub const fn exact_affine_many_predecessor_cell_visits(self) -> usize {
        self.exact_affine_many_predecessor_cell_visits
    }

    /// Returns many-predecessor cells finalized by the runtime-dispatched AVX2 scan.
    #[must_use]
    pub const fn exact_affine_many_predecessor_avx2_cells(self) -> usize {
        self.exact_affine_many_predecessor_avx2_cells
    }

    fn filled_exact_affine_cells(&self) -> Result<usize, PoaError> {
        self.exact_affine_source_cells
            .checked_add(self.exact_affine_single_predecessor_cells)
            .and_then(|cells| cells.checked_add(self.exact_affine_two_predecessor_cells))
            .and_then(|cells| cells.checked_add(self.exact_affine_many_predecessor_cells))
            .ok_or(PoaError::ArithmeticOverflow {
                context: "POASTA detailed exact-affine routed-cell counter",
            })
    }

    fn add(target: &mut usize, value: usize, context: &'static str) -> Result<(), PoaError> {
        *target = target
            .checked_add(value)
            .ok_or(PoaError::ArithmeticOverflow { context })?;
        Ok(())
    }

    fn record_source_path_certificate(&mut self) -> Result<(), PoaError> {
        Self::add(
            &mut self.source_path_certificates,
            1,
            "POASTA detailed source-path certificate counter",
        )
    }

    fn record_early_diagonal_certificate(&mut self) -> Result<(), PoaError> {
        Self::add(
            &mut self.early_diagonal_certificates,
            1,
            "POASTA detailed early-diagonal certificate counter",
        )
    }

    fn record_relaxed_linear_table_outcome(
        &mut self,
        cells: usize,
        outcome: RelaxedLinearTableOutcome,
    ) -> Result<(), PoaError> {
        Self::add(
            &mut self.relaxed_linear_table_attempts,
            1,
            "POASTA detailed relaxed-table attempt counter",
        )?;
        match outcome {
            RelaxedLinearTableOutcome::DiagonalCertificate => {
                Self::add(
                    &mut self.relaxed_linear_table_diagonal_certificates,
                    1,
                    "POASTA detailed relaxed-table diagonal-certificate counter",
                )?;
                Self::add(
                    &mut self.relaxed_linear_table_certificate_cells,
                    cells,
                    "POASTA detailed relaxed-table certificate-cell counter",
                )
            }
            RelaxedLinearTableOutcome::AffineOpeningCertificate => {
                Self::add(
                    &mut self.relaxed_linear_table_affine_opening_certificates,
                    1,
                    "POASTA detailed relaxed-table affine-opening counter",
                )?;
                Self::add(
                    &mut self.relaxed_linear_table_certificate_cells,
                    cells,
                    "POASTA detailed relaxed-table certificate-cell counter",
                )
            }
            RelaxedLinearTableOutcome::Miss => {
                Self::add(
                    &mut self.relaxed_linear_table_misses,
                    1,
                    "POASTA detailed relaxed-table miss counter",
                )?;
                Self::add(
                    &mut self.relaxed_linear_table_miss_cells,
                    cells,
                    "POASTA detailed relaxed-table miss-cell counter",
                )
            }
        }
    }

    fn record_exact_affine_route(
        &mut self,
        route: ExactAffineRoute,
        cells: usize,
        certified: bool,
    ) -> Result<(), PoaError> {
        match route {
            ExactAffineRoute::Bypass => {
                Self::add(
                    &mut self.exact_affine_bypass_attempts,
                    1,
                    "POASTA detailed exact-affine bypass-attempt counter",
                )?;
                Self::add(
                    &mut self.exact_affine_bypass_cells,
                    cells,
                    "POASTA detailed exact-affine bypass-cell counter",
                )?;
                if certified {
                    Self::add(
                        &mut self.exact_affine_bypass_certificates,
                        1,
                        "POASTA detailed exact-affine bypass-certificate counter",
                    )?;
                }
            }
            ExactAffineRoute::AfterRelaxed => {
                Self::add(
                    &mut self.exact_affine_after_relaxed_attempts,
                    1,
                    "POASTA detailed post-relaxed exact-affine attempt counter",
                )?;
                Self::add(
                    &mut self.exact_affine_after_relaxed_cells,
                    cells,
                    "POASTA detailed post-relaxed exact-affine cell counter",
                )?;
                if certified {
                    Self::add(
                        &mut self.exact_affine_after_relaxed_certificates,
                        1,
                        "POASTA detailed post-relaxed exact-affine certificate counter",
                    )?;
                }
            }
        }
        Ok(())
    }

    fn record_relaxed_linear_row(&mut self, work: RelaxedLinearRowWork) -> Result<(), PoaError> {
        Self::add(
            &mut self.relaxed_linear_rows,
            1,
            "POASTA detailed relaxed-linear row counter",
        )?;
        Self::add(
            &mut self.relaxed_linear_insertion_cells,
            work.insertion_cells,
            "POASTA detailed relaxed-linear insertion-cell counter",
        )?;
        Self::add(
            &mut self.relaxed_linear_insertion_avx2_cells,
            work.insertion_avx2_cells,
            "POASTA detailed relaxed-linear insertion-AVX2-cell counter",
        )?;
        self.record_relaxed_linear_topology(work.node_successors)?;
        Self::add(
            &mut self.relaxed_linear_successor_cell_visits,
            work.successor_cell_visits,
            "POASTA detailed relaxed-linear successor-cell-visit counter",
        )?;
        if work.node_successors < 2 {
            return Ok(());
        }
        Self::add(
            &mut self.relaxed_linear_multiple_node_successor_cells,
            work.row_cells,
            "POASTA detailed relaxed-linear multiple-successor-cell counter",
        )?;
        Self::add(
            &mut self.relaxed_linear_multiple_node_successor_cell_visits,
            work.successor_cell_visits,
            "POASTA detailed relaxed-linear multiple-successor-visit counter",
        )?;
        let distinct_symbol_cell_visits = work
            .row_cells
            .checked_mul(work.distinct_node_successor_symbols)
            .ok_or(PoaError::ArithmeticOverflow {
                context: "POASTA detailed distinct-symbol successor visits",
            })?;
        Self::add(
            &mut self.relaxed_linear_multiple_node_distinct_symbol_cell_visits,
            distinct_symbol_cell_visits,
            "POASTA detailed distinct-symbol successor-visit counter",
        )?;
        if work.equivalent_multiple_transition {
            Self::add(
                &mut self.relaxed_linear_equivalent_multiple_transition_rows,
                1,
                "POASTA detailed equivalent multiple-transition row counter",
            )?;
            Self::add(
                &mut self.relaxed_linear_equivalent_multiple_transition_cells,
                work.row_cells,
                "POASTA detailed equivalent multiple-transition cell counter",
            )?;
            Self::add(
                &mut self.relaxed_linear_equivalent_multiple_transition_insertion_cells,
                work.insertion_cells,
                "POASTA detailed equivalent multiple-transition insertion-cell counter",
            )?;
            Self::add(
                &mut self.relaxed_linear_equivalent_multiple_transition_successor_cell_visits,
                work.successor_cell_visits,
                "POASTA detailed equivalent multiple-transition successor-visit counter",
            )?;
        }
        Ok(())
    }

    fn record_relaxed_linear_topology(&mut self, node_successors: usize) -> Result<(), PoaError> {
        match node_successors {
            0 => Self::add(
                &mut self.relaxed_linear_zero_node_successor_rows,
                1,
                "POASTA detailed relaxed-linear zero-successor-row counter",
            )?,
            1 => Self::add(
                &mut self.relaxed_linear_single_node_successor_rows,
                1,
                "POASTA detailed relaxed-linear single-successor-row counter",
            )?,
            _ => Self::add(
                &mut self.relaxed_linear_multiple_node_successor_rows,
                1,
                "POASTA detailed relaxed-linear multiple-successor-row counter",
            )?,
        }
        match node_successors {
            2 => Self::add(
                &mut self.relaxed_linear_two_node_successor_rows,
                1,
                "POASTA detailed relaxed-linear two-successor-row counter",
            )?,
            3 | 4 => Self::add(
                &mut self.relaxed_linear_three_or_four_node_successor_rows,
                1,
                "POASTA detailed relaxed-linear three/four-successor-row counter",
            )?,
            5.. => Self::add(
                &mut self.relaxed_linear_five_plus_node_successor_rows,
                1,
                "POASTA detailed relaxed-linear five-plus-successor-row counter",
            )?,
            _ => {}
        }
        Ok(())
    }

    fn record_relaxed_linear_copied_row(
        &mut self,
        row_cells: usize,
        replaceable_insertion_cells: usize,
        node_successors: usize,
        replaceable_successor_cell_visits: usize,
    ) -> Result<(), PoaError> {
        Self::add(
            &mut self.relaxed_linear_rows,
            1,
            "POASTA detailed relaxed-linear row counter",
        )?;
        self.record_relaxed_linear_topology(node_successors)?;
        Self::add(
            &mut self.relaxed_linear_equivalent_single_transition_rows,
            1,
            "POASTA detailed equivalent single-transition row counter",
        )?;
        Self::add(
            &mut self.relaxed_linear_equivalent_single_transition_cells,
            row_cells,
            "POASTA detailed equivalent single-transition cell counter",
        )?;
        Self::add(
            &mut self.relaxed_linear_equivalent_single_transition_insertion_cells,
            replaceable_insertion_cells,
            "POASTA detailed equivalent single-transition insertion-cell counter",
        )?;
        Self::add(
            &mut self.relaxed_linear_equivalent_single_transition_successor_cell_visits,
            replaceable_successor_cell_visits,
            "POASTA detailed equivalent single-transition successor-cell-visit counter",
        )
    }

    fn record_relaxed_linear_setup(
        &mut self,
        initialized_rows: usize,
        reused_rows: usize,
        initialized_cells: usize,
    ) -> Result<(), PoaError> {
        Self::add(
            &mut self.relaxed_linear_initialized_rows,
            initialized_rows,
            "POASTA detailed relaxed-linear initialized-row counter",
        )?;
        Self::add(
            &mut self.relaxed_linear_reused_rows,
            reused_rows,
            "POASTA detailed relaxed-linear reused-row counter",
        )?;
        Self::add(
            &mut self.relaxed_linear_initialized_cells,
            initialized_cells,
            "POASTA detailed relaxed-linear initialized-cell counter",
        )
    }

    fn record_exact_affine_attempt(&mut self) -> Result<(), PoaError> {
        Self::add(
            &mut self.exact_affine_attempts,
            1,
            "POASTA detailed exact-affine attempt counter",
        )
    }

    fn record_exact_affine_source(&mut self, cells: usize) -> Result<(), PoaError> {
        Self::add(
            &mut self.exact_affine_filled_attempts,
            1,
            "POASTA detailed filled exact-affine attempt counter",
        )?;
        Self::add(
            &mut self.exact_affine_source_cells,
            cells,
            "POASTA detailed exact-affine source-cell counter",
        )
    }

    fn record_exact_affine_single(
        &mut self,
        cells: usize,
        fast_cells: usize,
        avx2_cells: usize,
    ) -> Result<(), PoaError> {
        Self::add(
            &mut self.exact_affine_single_predecessor_rows,
            1,
            "POASTA detailed exact-affine single-row counter",
        )?;
        Self::add(
            &mut self.exact_affine_single_predecessor_cells,
            cells,
            "POASTA detailed exact-affine single-cell counter",
        )?;
        Self::add(
            &mut self.exact_affine_single_predecessor_fast_cells,
            fast_cells,
            "POASTA detailed exact-affine single-fast-cell counter",
        )?;
        Self::add(
            &mut self.exact_affine_single_predecessor_avx2_cells,
            avx2_cells,
            "POASTA detailed exact-affine single-AVX2-cell counter",
        )
    }

    fn record_exact_affine_two(
        &mut self,
        cells: usize,
        fast_cells: usize,
        avx2_cells: usize,
    ) -> Result<(), PoaError> {
        Self::add(
            &mut self.exact_affine_two_predecessor_rows,
            1,
            "POASTA detailed exact-affine two-row counter",
        )?;
        Self::add(
            &mut self.exact_affine_two_predecessor_cells,
            cells,
            "POASTA detailed exact-affine two-cell counter",
        )?;
        Self::add(
            &mut self.exact_affine_two_predecessor_fast_cells,
            fast_cells,
            "POASTA detailed exact-affine two-fast-cell counter",
        )?;
        Self::add(
            &mut self.exact_affine_two_predecessor_avx2_cells,
            avx2_cells,
            "POASTA detailed exact-affine two-AVX2-cell counter",
        )
    }

    fn record_exact_affine_many(
        &mut self,
        cells: usize,
        predecessor_cell_visits: usize,
        avx2_cells: usize,
    ) -> Result<(), PoaError> {
        Self::add(
            &mut self.exact_affine_many_predecessor_rows,
            1,
            "POASTA detailed exact-affine many-row counter",
        )?;
        Self::add(
            &mut self.exact_affine_many_predecessor_cells,
            cells,
            "POASTA detailed exact-affine many-cell counter",
        )?;
        Self::add(
            &mut self.exact_affine_many_predecessor_cell_visits,
            predecessor_cell_visits,
            "POASTA detailed exact-affine many-predecessor visit counter",
        )?;
        Self::add(
            &mut self.exact_affine_many_predecessor_avx2_cells,
            avx2_cells,
            "POASTA detailed exact-affine many-AVX2-cell counter",
        )
    }
}

/// Exact counters describing one sparse POASTA search.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PoastaSearchStatistics {
    total_possible_states: usize,
    discovered_states: usize,
    expanded_states: usize,
    astar_expanded_states: usize,
    dfs_expanded_states: usize,
    queue_entries: usize,
    reopened_states: usize,
    stale_entries: usize,
    greedy_exact_matches: usize,
    affine_dominated_states: usize,
    linear_heuristic_cells: usize,
    relaxed_linear_source_bound: Option<Cost>,
    relaxed_linear_certificate: bool,
    exact_diagonal_certificate_cells: usize,
    exact_diagonal_cost: Option<Cost>,
    affine_opening_certificate: bool,
    exact_affine_certificate_cells: usize,
    exact_affine_certificate: bool,
    bubble_tests: usize,
    bubble_pruned_states: usize,
    reached_exit_states: usize,
    superbubbles: usize,
    bubble_membership_entries: usize,
    peak_queue_entries: usize,
    peak_dfs_entries: usize,
}

impl PoastaSearchStatistics {
    /// Returns the finite dense alignment-graph state bound.
    #[must_use]
    pub const fn total_possible_states(self) -> usize {
        self.total_possible_states
    }

    /// Returns the number of distinct sparse states ever discovered.
    #[must_use]
    pub const fn discovered_states(self) -> usize {
        self.discovered_states
    }

    /// Returns states expanded after stale and dominance checks.
    #[must_use]
    pub const fn expanded_states(self) -> usize {
        self.expanded_states
    }

    /// Returns states expanded directly from the A* frontier.
    #[must_use]
    pub const fn astar_expanded_states(self) -> usize {
        self.astar_expanded_states
    }

    /// Returns match states expanded by the depth-first frontier.
    #[must_use]
    pub const fn dfs_expanded_states(self) -> usize {
        self.dfs_expanded_states
    }

    /// Returns all priority-queue entries inserted, including later stale ones.
    #[must_use]
    pub const fn queue_entries(self) -> usize {
        self.queue_entries
    }

    /// Returns lower-cost arrivals that reopened already settled states.
    #[must_use]
    pub const fn reopened_states(self) -> usize {
        self.reopened_states
    }

    /// Returns queue or DFS entries skipped as stale or already settled.
    #[must_use]
    pub const fn stale_entries(self) -> usize {
        self.stale_entries
    }

    /// Returns newly improved exact matches traversed through DFS.
    #[must_use]
    pub const fn greedy_exact_matches(self) -> usize {
        self.greedy_exact_matches
    }

    /// Returns open affine-gap states removed by strict same-cell dominance.
    #[must_use]
    pub const fn affine_dominated_states(self) -> usize {
        self.affine_dominated_states
    }

    /// Returns cells refined by the relaxed linear-gap lower-bound table.
    #[must_use]
    pub const fn linear_heuristic_cells(self) -> usize {
        self.linear_heuristic_cells
    }

    /// Returns the relaxed source-state lower bound when its table was built.
    #[must_use]
    pub const fn relaxed_linear_source_bound(self) -> Option<u64> {
        self.relaxed_linear_source_bound
    }

    /// Returns whether a feasible diagonal path met the relaxed lower bound.
    #[must_use]
    pub const fn used_relaxed_linear_certificate(self) -> bool {
        self.relaxed_linear_certificate
    }

    /// Returns cells evaluated by the exact diagonal-only certificate.
    #[must_use]
    pub const fn exact_diagonal_certificate_cells(self) -> usize {
        self.exact_diagonal_certificate_cells
    }

    /// Returns the best complete diagonal-path cost when it was evaluated.
    #[must_use]
    pub const fn exact_diagonal_cost(self) -> Option<u64> {
        self.exact_diagonal_cost
    }

    /// Returns whether mandatory affine openings certified the diagonal optimum.
    #[must_use]
    pub const fn used_affine_opening_certificate(self) -> bool {
        self.affine_opening_certificate
    }

    /// Returns cells evaluated by the compact exact-affine certificate.
    #[must_use]
    pub const fn exact_affine_certificate_cells(self) -> usize {
        self.exact_affine_certificate_cells
    }

    /// Returns whether compact exact-affine DP supplied the final traceback.
    #[must_use]
    pub const fn used_exact_affine_certificate(self) -> bool {
        self.exact_affine_certificate
    }

    /// Returns node-to-bubble dominance tests performed.
    #[must_use]
    pub const fn bubble_tests(self) -> usize {
        self.bubble_tests
    }

    /// Returns states safely removed by superbubble dominance.
    #[must_use]
    pub const fn bubble_pruned_states(self) -> usize {
        self.bubble_pruned_states
    }

    /// Returns settled match-state observations inserted at bubble exits.
    #[must_use]
    pub const fn reached_exit_states(self) -> usize {
        self.reached_exit_states
    }

    /// Returns superbubbles found by Gartner et al. Algorithm 2.
    #[must_use]
    pub const fn superbubbles(self) -> usize {
        self.superbubbles
    }

    /// Returns output-sensitive node-to-containing-bubble entries.
    #[must_use]
    pub const fn bubble_membership_entries(self) -> usize {
        self.bubble_membership_entries
    }

    /// Returns the maximum number of simultaneous priority-queue entries.
    #[must_use]
    pub const fn peak_queue_entries(self) -> usize {
        self.peak_queue_entries
    }

    /// Returns the maximum simultaneous exact-match DFS depth.
    #[must_use]
    pub const fn peak_dfs_entries(self) -> usize {
        self.peak_dfs_entries
    }
}

/// Checked logical payload model for one sparse POASTA call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(clippy::struct_field_names)]
pub struct PoastaMemoryEstimate {
    graph_index_bytes: usize,
    heuristic_payload_bytes: usize,
    state_payload_bytes: usize,
    frontier_payload_peak_bytes: usize,
    dfs_payload_peak_bytes: usize,
    reached_exit_payload_bytes: usize,
    traceback_output_bytes: usize,
    modeled_peak_bytes: usize,
}

impl PoastaMemoryEstimate {
    /// Returns modeled bytes for augmented adjacency, distances, and bubbles.
    #[must_use]
    pub const fn graph_index_bytes(self) -> usize {
        self.graph_index_bytes
    }

    /// Returns modeled bytes for query-specific admissible-heuristic tables.
    #[must_use]
    pub const fn heuristic_payload_bytes(self) -> usize {
        self.heuristic_payload_bytes
    }

    /// Returns modeled payload bytes at final sparse-map capacity.
    #[must_use]
    pub const fn state_payload_bytes(self) -> usize {
        self.state_payload_bytes
    }

    /// Returns modeled payload bytes at peak A* frontier length.
    #[must_use]
    pub const fn frontier_payload_peak_bytes(self) -> usize {
        self.frontier_payload_peak_bytes
    }

    /// Returns modeled payload bytes at peak exact-match DFS depth.
    #[must_use]
    pub const fn dfs_payload_peak_bytes(self) -> usize {
        self.dfs_payload_peak_bytes
    }

    /// Returns modeled payload bytes for ordered reached-exit observations.
    #[must_use]
    pub const fn reached_exit_payload_bytes(self) -> usize {
        self.reached_exit_payload_bytes
    }

    /// Returns bytes owned by returned operation and aligned-pair slices.
    #[must_use]
    pub const fn traceback_output_bytes(self) -> usize {
        self.traceback_output_bytes
    }

    /// Returns the sum of simultaneous modeled payloads.
    #[must_use]
    pub const fn modeled_peak_bytes(self) -> usize {
        self.modeled_peak_bytes
    }
}

/// One deterministic exact global alignment returned by POASTA.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PoastaAlignment {
    score: Score,
    query_len: usize,
    operations: Box<[PoaOperation]>,
    aligned_pairs: Box<[PoaAlignedPair]>,
    statistics: PoastaSearchStatistics,
    memory: PoastaMemoryEstimate,
    options: PoastaSearchOptions,
    backend: PoastaBackend,
}

impl PoastaAlignment {
    /// Returns the exact maximizing score, equal to the negated POASTA cost.
    #[must_use]
    pub const fn score(&self) -> Score {
        self.score
    }

    /// Returns the complete consumed query range.
    #[must_use]
    pub const fn query_range(&self) -> Range<usize> {
        0..self.query_len
    }

    /// Returns canonical forward traceback operations.
    #[must_use]
    pub const fn operations(&self) -> &[PoaOperation] {
        &self.operations
    }

    /// Returns the diagonal graph-node/query-residue projection.
    #[must_use]
    pub const fn aligned_pairs(&self) -> &[PoaAlignedPair] {
        &self.aligned_pairs
    }

    /// Returns exact search counters for this call.
    #[must_use]
    pub const fn statistics(&self) -> PoastaSearchStatistics {
        self.statistics
    }

    /// Returns the logical owned-payload memory model for this call.
    #[must_use]
    pub const fn memory_estimate(&self) -> PoastaMemoryEstimate {
        self.memory
    }

    /// Returns feature-isolation options used for this result.
    #[must_use]
    pub const fn options(&self) -> PoastaSearchOptions {
        self.options
    }

    /// Returns the exact engine selected for this alignment.
    #[must_use]
    pub const fn backend(&self) -> PoastaBackend {
        self.backend
    }

    pub(crate) fn integration_targets(
        &self,
        query_len: usize,
    ) -> Result<Vec<Option<PoaNodeId>>, PoaError> {
        if query_len != self.query_len {
            return verification("POASTA integration target length differs from the query");
        }
        let mut targets = filled(query_len, None, "POASTA iterative integration targets")?;
        for pair in &self.aligned_pairs {
            if pair.query_index() >= query_len
                || targets[pair.query_index()].replace(pair.node()).is_some()
            {
                return verification(
                    "POASTA aligned-pair mapping is not one-to-one with the query",
                );
            }
        }
        Ok(targets)
    }
}

/// Exact sparse POASTA sequence-to-DAG aligner.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Poasta {
    costs: PoastaAffineCosts,
    options: PoastaSearchOptions,
}

impl Poasta {
    /// Creates the complete published profile with every acceleration enabled.
    #[must_use]
    pub const fn new(costs: PoastaAffineCosts) -> Self {
        Self {
            costs,
            options: PoastaSearchOptions::published(),
        }
    }

    /// Creates the exact profile using the validated optimized-execution policy.
    #[must_use]
    pub const fn optimized(costs: PoastaAffineCosts) -> Self {
        Self {
            costs,
            options: PoastaSearchOptions::optimized(),
        }
    }

    /// Creates an exact feature-isolation configuration.
    #[must_use]
    pub const fn with_options(costs: PoastaAffineCosts, options: PoastaSearchOptions) -> Self {
        Self { costs, options }
    }

    /// Returns the configured published costs.
    #[must_use]
    pub const fn costs(self) -> PoastaAffineCosts {
        self.costs
    }

    /// Returns the configured exact-search features.
    #[must_use]
    pub const fn options(self) -> PoastaSearchOptions {
        self.options
    }

    /// Returns a profile identifier that does not mislabel isolation runs.
    #[must_use]
    pub const fn profile_identifier(self) -> &'static str {
        if matches!(self.options.backend, PoastaBackend::Dense) {
            return POASTA_EXACT_DENSE_PROFILE;
        }
        if matches!(self.options.backend, PoastaBackend::Auto) {
            return POASTA_EXACT_AUTO_PROFILE;
        }
        match (
            self.options.exact_match_dfs,
            self.options.superbubble_pruning,
        ) {
            (true, true) => POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_V1,
            (true, false) => POASTA_ASTAR_DFS_VALIDATION_PROFILE,
            (false, true) => POASTA_ASTAR_PRUNING_VALIDATION_PROFILE,
            (false, false) => POASTA_ASTAR_VALIDATION_PROFILE,
        }
    }

    /// Aligns a complete query globally to one complete DAG path.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError`] for invalid input, checked dimension, allocation, or
    /// cost overflow, an unrepresentable public score, invalid traceback, or
    /// failed independent replay.
    pub fn align<A: Alphabet>(
        &self,
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
    ) -> Result<PoastaAlignment, PoaError> {
        let (topological, _) = graph.validate_with_orders()?;
        self.align_validated_with_report(graph, query, &topological, false, None)
            .map(|(alignment, _)| alignment)
    }

    /// Aligns one query and reports monotonic wall time by POASTA phase.
    ///
    /// This diagnostic entry point returns the same exact alignment as
    /// [`Self::align`]. Repeated phase names are aggregated, and timer reads
    /// occur only on this opt-in path.
    ///
    /// # Errors
    ///
    /// Returns the same [`PoaError`] values as [`Self::align`].
    pub fn align_profiled<A: Alphabet>(
        &self,
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
    ) -> Result<Profiled<PoastaAlignment>, PoaError> {
        let mut timings = PhaseTimingsBuilder::new();
        let (topological, _) =
            timings.measure("input_validation", || graph.validate_with_orders())?;
        let (alignment, _) = self.align_validated_with_report_profiled::<A, false>(
            graph,
            query,
            &topological,
            false,
            None,
            Some(&mut timings),
            None,
            None,
        )?;
        Ok(Profiled::new(alignment, timings.finish()))
    }

    pub(crate) fn align_validated_with_report<A: Alphabet>(
        &self,
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
        topological: &[PoaNodeId],
        eager_relaxed_linear: bool,
        workspace: Option<&mut PoastaWorkspace>,
    ) -> Result<(PoastaAlignment, AffinePoaVerificationReport), PoaError> {
        self.align_validated_with_report_profiled::<A, false>(
            graph,
            query,
            topological,
            eager_relaxed_linear,
            workspace,
            None,
            None,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn align_validated_with_report_profiled<A: Alphabet, const DETAILED: bool>(
        &self,
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
        topological: &[PoaNodeId],
        eager_relaxed_linear: bool,
        workspace: Option<&mut PoastaWorkspace>,
        mut timings: Option<&mut PhaseTimingsBuilder>,
        search_timings: Option<&mut PhaseTimingsBuilder>,
        search_counters: Option<&mut PoastaSearchWorkCounters>,
    ) -> Result<(PoastaAlignment, AffinePoaVerificationReport), PoaError> {
        let total_possible_states = possible_state_count(graph.node_count(), query.len())?;
        let backend = resolve_backend(self.options.backend);
        if matches!(backend, PoastaBackend::Dense) {
            return self.align_dense_profiled(graph, query, total_possible_states, timings);
        }
        let graph_index = measure_phase(&mut timings, "graph_index", || {
            PoastaGraphIndex::build(graph, topological, self.options.superbubble_pruning)
        })?;
        let mut search = measure_phase(&mut timings, "search_setup", || {
            Search::new(
                graph,
                query,
                self.costs,
                self.options,
                &graph_index,
                topological,
                eager_relaxed_linear,
                total_possible_states,
                workspace,
            )
        })?;
        let (goal, optimum_cost) = measure_phase(&mut timings, "search", || {
            search.run::<DETAILED>(search_timings, search_counters)
        })?;
        let (operations, aligned_pairs) =
            measure_phase(&mut timings, "traceback", || search.traceback(goal))?;
        let alignment = measure_phase(
            &mut timings,
            "result_accounting",
            || -> Result<_, PoaError> {
                let score = cost_to_score(optimum_cost)?;
                let certified_state_count = match search.certified_traceback.as_ref() {
                    Some(traceback) => traceback
                        .len()
                        .checked_add(usize::from(search.records.is_empty()))
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "POASTA certified state count",
                        })?,
                    None => 0,
                };
                search.statistics.discovered_states = search
                    .records
                    .len()
                    .checked_add(certified_state_count)
                    .ok_or(PoaError::ArithmeticOverflow {
                        context: "POASTA discovered state count",
                    })?;
                let memory = search.memory_estimate(&operations, &aligned_pairs)?;
                let statistics = search.statistics;
                search.recycle_relaxed_linear_costs();
                Ok(PoastaAlignment {
                    score,
                    query_len: query.len(),
                    operations: operations.into_boxed_slice(),
                    aligned_pairs: aligned_pairs.into_boxed_slice(),
                    statistics,
                    memory,
                    options: self.options,
                    backend,
                })
            },
        )?;
        let verification = measure_phase(&mut timings, "verification", || {
            verify_global_affine_operations_validated(
                alignment.score,
                alignment.query_len,
                &alignment.operations,
                &alignment.aligned_pairs,
                graph,
                query,
                &self.costs,
                self.costs.affine_gap(),
            )
        })?;
        Ok((alignment, verification))
    }

    fn align_dense_profiled<A: Alphabet>(
        &self,
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
        total_possible_states: usize,
        mut timings: Option<&mut PhaseTimingsBuilder>,
    ) -> Result<(PoastaAlignment, AffinePoaVerificationReport), PoaError> {
        let dense = measure_phase(&mut timings, "dense_alignment", || {
            GlobalAffinePoa::from_poasta_costs(self.costs).align(graph, query)
        })?;
        let alignment = measure_phase(
            &mut timings,
            "result_accounting",
            || -> Result<_, PoaError> {
                let (score, query_len, operations, aligned_pairs, dense_memory) =
                    dense.into_parts();
                let memory = dense_memory_estimate(dense_memory, &operations, &aligned_pairs)?;
                let statistics = dense_statistics(total_possible_states);
                Ok(PoastaAlignment {
                    score,
                    query_len,
                    operations,
                    aligned_pairs,
                    statistics,
                    memory,
                    options: self.options,
                    backend: PoastaBackend::Dense,
                })
            },
        )?;
        let verification = measure_phase(&mut timings, "verification", || {
            verify_global_affine_operations_validated(
                alignment.score,
                alignment.query_len,
                &alignment.operations,
                &alignment.aligned_pairs,
                graph,
                query,
                &self.costs,
                self.costs.affine_gap(),
            )
        })?;
        Ok((alignment, verification))
    }
}

fn measure_phase<T>(
    timings: &mut Option<&mut PhaseTimingsBuilder>,
    name: &'static str,
    operation: impl FnOnce() -> T,
) -> T {
    if let Some(timings) = timings.as_deref_mut() {
        timings.measure(name, operation)
    } else {
        operation()
    }
}

fn measure_search_phase<const DETAILED: bool, T>(
    timings: &mut Option<&mut PhaseTimingsBuilder>,
    name: &'static str,
    operation: impl FnOnce() -> T,
) -> T {
    if DETAILED {
        measure_phase(timings, name, operation)
    } else {
        operation()
    }
}

fn search_phase_started<const DETAILED: bool>(
    timings: &Option<&mut PhaseTimingsBuilder>,
) -> Option<Instant> {
    if DETAILED && timings.is_some() {
        Some(Instant::now())
    } else {
        None
    }
}

fn record_search_phase<const DETAILED: bool>(
    timings: Option<&mut PhaseTimingsBuilder>,
    name: &'static str,
    started: Option<Instant>,
) {
    if DETAILED {
        if let (Some(timings), Some(started)) = (timings, started) {
            timings.record(name, started.elapsed());
        }
    }
}

/// Independently replays a POASTA result without reading sparse search records.
///
/// # Errors
///
/// Returns [`PoaError`] when path continuity, complete query consumption, affine
/// gap scoring, aligned-pair projection, or the reported score is invalid.
pub fn verify_poasta_alignment<A: Alphabet>(
    alignment: &PoastaAlignment,
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    costs: PoastaAffineCosts,
) -> Result<AffinePoaVerificationReport, PoaError> {
    graph.validate()?;
    verify_global_affine_operations_validated(
        alignment.score,
        alignment.query_len,
        &alignment.operations,
        &alignment.aligned_pairs,
        graph,
        query,
        &costs,
        costs.affine_gap(),
    )
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[repr(u8)]
enum AffineState {
    Match = 0,
    Deletion = 1,
    Insertion = 2,
}

impl AffineState {
    const fn name(self) -> &'static str {
        match self {
            Self::Match => "M",
            Self::Deletion => "D",
            Self::Insertion => "I",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
struct StateKey {
    cursor: usize,
    query_state: usize,
}

impl StateKey {
    const NONE: Self = Self {
        cursor: usize::MAX,
        query_state: usize::MAX,
    };

    const fn new(cursor: usize, query_prefix: usize, state: AffineState) -> Self {
        Self {
            cursor,
            query_state: query_prefix * 3 + state as usize,
        }
    }

    const fn query_prefix(self) -> usize {
        self.query_state / 3
    }

    const fn state(self) -> AffineState {
        match self.query_state % 3 {
            0 => AffineState::Match,
            1 => AffineState::Deletion,
            _ => AffineState::Insertion,
        }
    }

    const fn is_none(self) -> bool {
        self.cursor == usize::MAX && self.query_state == usize::MAX
    }
}

const UNEXPANDED_COST: Cost = Cost::MAX;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct StateRecord {
    cost: Cost,
    predecessor: StateKey,
    expanded_cost: Cost,
    tie_transition: u8,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct QueueEntry {
    estimated_total: Cost,
    cost: Cost,
    key: StateKey,
}

const MAXIMUM_BUCKET_PRIORITY: usize = 65_536;

#[derive(Debug, Default)]
struct BucketFrontier {
    buckets: Vec<Vec<QueueEntry>>,
    overflow: BinaryHeap<Reverse<QueueEntry>>,
    next_bucket: usize,
    len: usize,
}

impl BucketFrontier {
    #[inline]
    fn try_push(&mut self, entry: QueueEntry) -> Result<(), PoaError> {
        let Ok(priority) = usize::try_from(entry.estimated_total) else {
            return self.try_push_overflow(entry);
        };
        if priority > MAXIMUM_BUCKET_PRIORITY {
            return self.try_push_overflow(entry);
        }
        if priority >= self.buckets.len()
            || self.buckets[priority].len() == self.buckets[priority].capacity()
        {
            return self.try_push_bucket_slow(entry, priority);
        }
        self.buckets[priority].push(entry);
        self.record_bucket_push(priority)
    }

    #[cold]
    #[inline(never)]
    fn try_push_bucket_slow(&mut self, entry: QueueEntry, priority: usize) -> Result<(), PoaError> {
        if priority >= self.buckets.len() {
            let additional = priority + 1 - self.buckets.len();
            self.buckets
                .try_reserve(additional)
                .map_err(|_| PoaError::AllocationFailed {
                    buffer: "POASTA bucket frontier",
                    elements: priority + 1,
                })?;
            self.buckets.resize_with(priority + 1, Vec::new);
        }
        let bucket = &mut self.buckets[priority];
        if bucket.len() == bucket.capacity() {
            bucket
                .try_reserve(1)
                .map_err(|_| PoaError::AllocationFailed {
                    buffer: "POASTA bucket entries",
                    elements: 1,
                })?;
        }
        bucket.push(entry);
        self.record_bucket_push(priority)
    }

    #[inline]
    fn record_bucket_push(&mut self, priority: usize) -> Result<(), PoaError> {
        if self.len == 0 || priority < self.next_bucket {
            self.next_bucket = priority;
        }
        self.len = checked_increment(self.len, "POASTA bucket-frontier length")?;
        Ok(())
    }

    #[cold]
    #[inline(never)]
    fn try_push_overflow(&mut self, entry: QueueEntry) -> Result<(), PoaError> {
        if self.overflow.len() == self.overflow.capacity() {
            self.overflow
                .try_reserve(1)
                .map_err(|_| PoaError::AllocationFailed {
                    buffer: "POASTA overflow frontier",
                    elements: 1,
                })?;
        }
        self.overflow.push(Reverse(entry));
        self.len = checked_increment(self.len, "POASTA bucket-frontier length")?;
        Ok(())
    }

    fn pop(&mut self) -> Option<Reverse<QueueEntry>> {
        while self.next_bucket < self.buckets.len() && self.buckets[self.next_bucket].is_empty() {
            self.next_bucket += 1;
        }
        if self.next_bucket < self.buckets.len() {
            let entry = self.buckets[self.next_bucket].pop()?;
            self.len -= 1;
            return Some(Reverse(entry));
        }
        let entry = self.overflow.pop()?;
        self.len -= 1;
        Some(entry)
    }

    const fn len(&self) -> usize {
        self.len
    }

    fn modeled_payload_bytes(&self) -> Result<usize, PoaError> {
        let bucket_table = checked_bytes(
            self.buckets.capacity(),
            size_of::<Vec<QueueEntry>>(),
            "POASTA bucket-frontier table",
        )?;
        let bucket_entries =
            self.buckets
                .iter()
                .try_fold(0_usize, |total, bucket| -> Result<usize, PoaError> {
                    let bytes = checked_bytes(
                        bucket.capacity(),
                        size_of::<QueueEntry>(),
                        "POASTA bucket-frontier entries",
                    )?;
                    total
                        .checked_add(bytes)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "POASTA bucket-frontier entry payload",
                        })
                })?;
        let overflow = checked_bytes(
            self.overflow.capacity(),
            size_of::<Reverse<QueueEntry>>(),
            "POASTA overflow-frontier entries",
        )?;
        bucket_table
            .checked_add(bucket_entries)
            .and_then(|bytes| bytes.checked_add(overflow))
            .ok_or(PoaError::ArithmeticOverflow {
                context: "POASTA bucket-frontier payload",
            })
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Relaxation {
    Improved,
    TraceOnly,
    Unchanged,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct DfsFrame {
    key: StateKey,
    cost: Cost,
    next_successor: usize,
    entered: bool,
    found_new_exact_match: bool,
    requires_insertion: bool,
}

impl DfsFrame {
    const fn new(key: StateKey, cost: Cost) -> Self {
        Self {
            key,
            cost,
            next_successor: 0,
            entered: false,
            found_new_exact_match: false,
            requires_insertion: false,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ReachedPosition {
    query_prefix: usize,
    cost: Cost,
}

#[derive(Clone, Debug, Default)]
struct ReachedExit {
    positions: Vec<ReachedPosition>,
    prefix_insertion_minimum: Vec<Cost>,
    deletion_range_tree: Vec<Option<Cost>>,
    deletion_tree_base: usize,
    cached_query_len: usize,
    cached_extension: Cost,
    bounds_dirty: bool,
}

impl ReachedExit {
    fn insert(&mut self, query_prefix: usize, cost: Cost) -> Result<bool, PoaError> {
        match self
            .positions
            .binary_search_by_key(&query_prefix, |position| position.query_prefix)
        {
            Ok(index) if cost < self.positions[index].cost => {
                self.positions[index].cost = cost;
                self.bounds_dirty = true;
                Ok(true)
            }
            Ok(_) => Ok(false),
            Err(index) => {
                if self.positions.len() == self.positions.capacity() {
                    self.positions
                        .try_reserve(1)
                        .map_err(|_| PoaError::AllocationFailed {
                            buffer: "POASTA reached bubble exits",
                            elements: 1,
                        })?;
                }
                self.positions
                    .insert(index, ReachedPosition { query_prefix, cost });
                self.bounds_dirty = true;
                Ok(true)
            }
        }
    }

    fn upper_bound(
        &mut self,
        query_prefix: usize,
        maximum_deletion: usize,
        query_len: usize,
        costs: PoastaAffineCosts,
    ) -> Result<Option<Cost>, PoaError> {
        if self.positions.is_empty() || query_prefix > query_len {
            return Ok(None);
        }
        self.ensure_bound_index(query_len, costs)?;
        let mut bound = self
            .positions
            .binary_search_by_key(&query_prefix, |position| position.query_prefix)
            .ok()
            .map(|index| self.positions[index].cost);
        let insertion_end = self
            .positions
            .partition_point(|position| position.query_prefix < query_prefix);
        if insertion_end > 0 {
            let transformed = self.prefix_insertion_minimum[insertion_end - 1];
            let tail = checked_cost_mul(
                query_len - query_prefix,
                self.cached_extension,
                "an implicit insertion suffix",
            )?;
            let insertion = checked_cost_add(
                transformed,
                nonnegative_score(costs.gap_opening()),
                "an implicit insertion bound",
            )?
            .checked_sub(tail)
            .ok_or(PoaError::PoastaCostOverflow {
                context: "an implicit insertion bound subtraction",
            })?;
            bound = optional_min(bound, Some(insertion));
        }

        let deletion_start = self
            .positions
            .partition_point(|position| position.query_prefix <= query_prefix);
        let maximum_position = query_prefix
            .checked_add(maximum_deletion)
            .map_or(query_len, |position| position.min(query_len));
        let deletion_end = self
            .positions
            .partition_point(|position| position.query_prefix <= maximum_position);
        if let Some(transformed) = self.deletion_range_minimum(deletion_start, deletion_end) {
            let prefix = checked_cost_mul(
                query_prefix,
                self.cached_extension,
                "an implicit deletion prefix",
            )?;
            let deletion = checked_cost_add(
                transformed,
                nonnegative_score(costs.gap_opening()),
                "an implicit deletion bound",
            )?
            .checked_sub(prefix)
            .ok_or(PoaError::PoastaCostOverflow {
                context: "an implicit deletion bound subtraction",
            })?;
            bound = optional_min(bound, Some(deletion));
        }
        Ok(bound)
    }

    fn ensure_bound_index(
        &mut self,
        query_len: usize,
        costs: PoastaAffineCosts,
    ) -> Result<(), PoaError> {
        let extension = nonnegative_score(costs.gap_extension());
        if !self.bounds_dirty
            && self.prefix_insertion_minimum.len() == self.positions.len()
            && self.cached_query_len == query_len
            && self.cached_extension == extension
        {
            return Ok(());
        }
        self.prefix_insertion_minimum.clear();
        try_reserve_exact(
            &mut self.prefix_insertion_minimum,
            self.positions.len(),
            "POASTA reached-exit prefix bounds",
        )?;
        let mut prefix_minimum = None;
        for position in &self.positions {
            let remaining = query_len.checked_sub(position.query_prefix).ok_or(
                PoaError::ArithmeticOverflow {
                    context: "a reached-exit query suffix",
                },
            )?;
            let transformed = checked_cost_add(
                position.cost,
                checked_cost_mul(remaining, extension, "a reached-exit insertion transform")?,
                "a reached-exit insertion transform",
            )?;
            let minimum =
                prefix_minimum.map_or(transformed, |current| min_cost(current, transformed));
            prefix_minimum = Some(minimum);
            self.prefix_insertion_minimum.push(minimum);
        }

        self.deletion_tree_base = self.positions.len().checked_next_power_of_two().ok_or(
            PoaError::ArithmeticOverflow {
                context: "the reached-exit range-tree base",
            },
        )?;
        let tree_len =
            self.deletion_tree_base
                .checked_mul(2)
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "the reached-exit range-tree length",
                })?;
        self.deletion_range_tree =
            filled(tree_len, None, "POASTA reached-exit deletion range tree")?;
        for (index, position) in self.positions.iter().enumerate() {
            let transformed = checked_cost_add(
                position.cost,
                checked_cost_mul(
                    position.query_prefix,
                    extension,
                    "a reached-exit deletion transform",
                )?,
                "a reached-exit deletion transform",
            )?;
            self.deletion_range_tree[self.deletion_tree_base + index] = Some(transformed);
        }
        for index in (1..self.deletion_tree_base).rev() {
            self.deletion_range_tree[index] = optional_min(
                self.deletion_range_tree[index * 2],
                self.deletion_range_tree[index * 2 + 1],
            );
        }
        self.cached_query_len = query_len;
        self.cached_extension = extension;
        self.bounds_dirty = false;
        Ok(())
    }

    fn deletion_range_minimum(&self, start: usize, end: usize) -> Option<Cost> {
        if start >= end {
            return None;
        }
        let mut left = self.deletion_tree_base + start;
        let mut right = self.deletion_tree_base + end;
        let mut minimum = None;
        while left < right {
            if left & 1 == 1 {
                minimum = optional_min(minimum, self.deletion_range_tree[left]);
                left += 1;
            }
            if right & 1 == 1 {
                right -= 1;
                minimum = optional_min(minimum, self.deletion_range_tree[right]);
            }
            left /= 2;
            right /= 2;
        }
        minimum
    }
}

#[derive(Default)]
struct StableHasher(u64);

impl Hasher for StableHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        let mut hash = if self.0 == 0 {
            0xcbf2_9ce4_8422_2325
        } else {
            self.0
        };
        for byte in bytes {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
        self.0 = hash;
    }

    fn write_usize(&mut self, value: usize) {
        self.0 = (value as u64).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
}

type StateMap = HashMap<usize, StateRecord, BuildHasherDefault<StableHasher>>;

// Radius 24 retains every exact relaxed-linear certificate in the branch-heavy
// crossover corpus; narrower bands force the slower exact-affine fallback.
const RELAXED_LINEAR_RADIUS: usize = 24;
const RELAXED_LINEAR_DISCOVERY_FACTOR: usize = 2;
// Progressive POA inputs often stay in the same difficulty regime for several
// records. After three failed relaxed certificates, prefer the exact compact DP
// and periodically retry the cheaper certificate so phase changes are detected.
const RELAXED_LINEAR_MISS_THRESHOLD: usize = 3;
const RELAXED_LINEAR_REPROBE_INTERVAL: usize = 8;

#[derive(Default)]
pub(crate) struct PoastaWorkspace {
    costs: Vec<u16>,
    row_bounds: Vec<RelaxedLinearRowBounds>,
    query_len: Option<usize>,
    extension: u16,
    exact_match_costs: Vec<u16>,
    exact_deletion_costs: Vec<u16>,
    exact_row_bands: Vec<ExactAffineRowBand>,
    relaxed_certificate_misses: usize,
    exact_affine_bypasses_since_probe: usize,
}

impl PoastaWorkspace {
    fn modeled_payload_bytes(&self) -> Result<usize, PoaError> {
        let cost_bytes = checked_bytes(
            self.costs.capacity(),
            size_of::<u16>(),
            "POASTA cached relaxed linear heuristic payload",
        )?;
        let row_bytes = checked_bytes(
            self.row_bounds.capacity(),
            size_of::<RelaxedLinearRowBounds>(),
            "POASTA cached relaxed linear row payload",
        )?;
        let exact_cell_capacity = self
            .exact_match_costs
            .capacity()
            .checked_add(self.exact_deletion_costs.capacity())
            .ok_or(PoaError::ArithmeticOverflow {
                context: "POASTA cached exact-affine cell capacity",
            })?;
        let exact_table_bytes = checked_bytes(
            exact_cell_capacity,
            size_of::<u16>(),
            "POASTA cached exact-affine table payload",
        )?;
        let exact_band_bytes = checked_bytes(
            self.exact_row_bands.capacity(),
            size_of::<ExactAffineRowBand>(),
            "POASTA cached exact-affine row-band payload",
        )?;
        checked_add_usize(
            checked_add_usize(
                cost_bytes,
                row_bytes,
                "POASTA cached relaxed linear table payload",
            )?,
            checked_add_usize(
                exact_table_bytes,
                exact_band_bytes,
                "POASTA cached exact-affine total payload",
            )?,
            "POASTA cached heuristic total payload",
        )
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct RelaxedLinearRowBounds {
    minimum: usize,
    maximum: usize,
}

const INVALID_RELAXED_LINEAR_ROW: RelaxedLinearRowBounds = RelaxedLinearRowBounds {
    minimum: usize::MAX,
    maximum: usize::MAX,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ExactAffineRowBand {
    first: usize,
    last: usize,
    offset: usize,
}

const INVALID_EXACT_AFFINE_ROW_BAND: ExactAffineRowBand = ExactAffineRowBand {
    first: usize::MAX,
    last: usize::MAX,
    offset: usize::MAX,
};

// The exact-prefix cache and this inlining decision are a measured code-generation
// pair: neither isolated variant is stable across the retained corpus, while the
// combination improves all fourteen final process cases.
#[inline(always)]
#[allow(
    clippy::inline_always,
    clippy::range_plus_one,
    clippy::too_many_arguments
)]
fn merge_relaxed_successor<const INITIALIZE: bool>(
    row: &mut [u16],
    successor_row: &[u16],
    query_indices: &[u8],
    graph_symbol: u8,
    mismatch: u16,
    extension: u16,
    first_query_prefix: usize,
    last_query_prefix: usize,
) {
    let diagonal_end = (last_query_prefix + 1).min(query_indices.len());
    let cells = &mut row[first_query_prefix..diagonal_end];
    let deletion_bases = &successor_row[first_query_prefix..diagonal_end];
    let diagonal_bases = &successor_row[first_query_prefix + 1..diagonal_end + 1];
    let query_symbols = &query_indices[first_query_prefix..diagonal_end];
    if INITIALIZE {
        for (((cell, &deletion_base), &diagonal_base), &query_symbol) in cells
            .iter_mut()
            .zip(deletion_bases)
            .zip(diagonal_bases)
            .zip(query_symbols)
        {
            let deletion = deletion_base + extension;
            let substitution = if query_symbol == graph_symbol {
                0
            } else {
                mismatch
            };
            let diagonal = diagonal_base + substitution;
            *cell = deletion.min(diagonal);
        }
        if last_query_prefix == query_indices.len() {
            row[last_query_prefix] = successor_row[last_query_prefix] + extension;
        }
    } else {
        for (((cell, &deletion_base), &diagonal_base), &query_symbol) in cells
            .iter_mut()
            .zip(deletion_bases)
            .zip(diagonal_bases)
            .zip(query_symbols)
        {
            let deletion = deletion_base + extension;
            let substitution = if query_symbol == graph_symbol {
                0
            } else {
                mismatch
            };
            let diagonal = diagonal_base + substitution;
            *cell = (*cell).min(deletion.min(diagonal));
        }
        if last_query_prefix == query_indices.len() {
            let terminal = &mut row[last_query_prefix];
            *terminal = (*terminal).min(successor_row[last_query_prefix] + extension);
        }
    }
}

#[inline]
fn relaxed_row_pair(
    values: &mut [u16],
    columns: usize,
    cursor: usize,
    successor: usize,
) -> (&mut [u16], &[u16]) {
    debug_assert_ne!(cursor, successor);
    let row_start = cursor * columns;
    let successor_start = successor * columns;
    if row_start < successor_start {
        let (before_successor, successor_and_after) = values.split_at_mut(successor_start);
        (
            &mut before_successor[row_start..row_start + columns],
            &successor_and_after[..columns],
        )
    } else {
        let (before_row, row_and_after) = values.split_at_mut(row_start);
        (
            &mut row_and_after[..columns],
            &before_row[successor_start..successor_start + columns],
        )
    }
}

struct RelaxedLinearHeuristic {
    columns: usize,
    computed_cells: usize,
    construction_scratch_bytes: usize,
    costs: Vec<u16>,
}

struct DiagonalCertificate {
    path: Vec<StateKey>,
    optimum_cost: Cost,
}

struct DiagonalCertificateAttempt {
    certificate: DiagonalCertificate,
    cells: usize,
    payload_bytes: usize,
}

#[allow(clippy::too_many_lines)]
fn exact_source_path_certificate<A: Alphabet>(
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    topological: &[PoaNodeId],
) -> Result<Option<DiagonalCertificate>, PoaError> {
    let mut source_path_exists = false;
    'paths: for source_path in graph.paths() {
        let candidate = source_path.nodes();
        if candidate.len() != query.len() {
            continue;
        }
        let Some((&first, &last)) = candidate.first().zip(candidate.last()) else {
            continue;
        };
        if !graph
            .node(first)
            .ok_or(PoaError::UnknownNode { node: first })?
            .predecessors()
            .is_empty()
            || !graph
                .node(last)
                .ok_or(PoaError::UnknownNode { node: last })?
                .successors()
                .is_empty()
        {
            continue;
        }
        for (query_prefix, &node_id) in candidate.iter().enumerate() {
            let query_symbol = query.get(query_prefix).ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: query.len(),
            })?;
            if graph
                .node(node_id)
                .ok_or(PoaError::UnknownNode { node: node_id })?
                .symbol()
                != query_symbol
            {
                continue 'paths;
            }
        }
        source_path_exists = true;
        break;
    }
    if !source_path_exists {
        return Ok(None);
    }

    // Once one stored path proves that cost zero is attainable, recover the
    // same zero-cost recombinant path selected by the published traceback:
    // smallest terminal cursor, then smallest predecessor cursor backwards.
    let columns = query
        .len()
        .checked_add(1)
        .ok_or(PoaError::DimensionOverflow {
            graph_nodes: graph.node_count(),
            query_len: query.len(),
        })?;
    let words = columns.checked_add(63).ok_or(PoaError::DimensionOverflow {
        graph_nodes: graph.node_count(),
        query_len: query.len(),
    })? / 64;
    let row_count = graph
        .node_count()
        .checked_add(1)
        .ok_or(PoaError::DimensionOverflow {
            graph_nodes: graph.node_count(),
            query_len: query.len(),
        })?;
    let reachability_elements =
        row_count
            .checked_mul(words)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: query.len(),
            })?;
    let mut reachability = filled(
        reachability_elements,
        0_u64,
        "POASTA exact source-path reachability",
    )?;
    reachability[0] = 1;
    let mask_elements = A::SYMBOLS
        .len()
        .checked_mul(words)
        .ok_or(PoaError::DimensionOverflow {
            graph_nodes: graph.node_count(),
            query_len: query.len(),
        })?;
    let mut symbol_masks = filled(
        mask_elements,
        0_u64,
        "POASTA exact source-path symbol masks",
    )?;
    for (query_index, &symbol) in query.as_indices().iter().enumerate() {
        let query_prefix = query_index + 1;
        symbol_masks[usize::from(symbol) * words + query_prefix / 64] |=
            1_u64 << (query_prefix % 64);
    }

    for &node_id in topological {
        let node = graph
            .node(node_id)
            .ok_or(PoaError::UnknownNode { node: node_id })?;
        let row_start = (node_id.index() + 1) * words;
        if node.predecessors().is_empty() {
            reachability.copy_within(0..words, row_start);
        } else {
            for predecessor in node.predecessors() {
                let predecessor_start = (predecessor.index() + 1) * words;
                for word in 0..words {
                    reachability[row_start + word] |= reachability[predecessor_start + word];
                }
            }
        }
        let mask_start = node.symbol().alphabet_index() * words;
        let mut carry = 0_u64;
        for word in 0..words {
            let unshifted = reachability[row_start + word];
            reachability[row_start + word] =
                ((unshifted << 1) | carry) & symbol_masks[mask_start + word];
            carry = unshifted >> 63;
        }
    }

    let goal_cursor = (1..=graph.node_count())
        .find(|&cursor| {
            let node_id = PoaNodeId::new(cursor - 1);
            graph.node(node_id).is_some_and(|node| {
                node.successors().is_empty()
                    && reachability[cursor * words + query.len() / 64]
                        & (1_u64 << (query.len() % 64))
                        != 0
            })
        })
        .ok_or_else(|| PoaError::InvalidGraph {
            message: "an exact source path is not reachable in its validated graph".to_owned(),
        })?;
    let mut reversed_path = Vec::new();
    try_reserve_exact(
        &mut reversed_path,
        query.len(),
        "POASTA exact source-path certificate",
    )?;
    let mut cursor = goal_cursor;
    let mut query_prefix = query.len();
    while cursor != 0 {
        reversed_path.push(StateKey::new(cursor, query_prefix, AffineState::Match));
        let node_id = PoaNodeId::new(cursor - 1);
        let node = graph
            .node(node_id)
            .ok_or(PoaError::UnknownNode { node: node_id })?;
        let predecessor_prefix =
            query_prefix
                .checked_sub(1)
                .ok_or_else(|| PoaError::InvalidGraph {
                    message: "an exact source path consumes too many graph nodes".to_owned(),
                })?;
        let predecessor = augmented_predecessor_cursors(node.predecessors(), 0)
            .find(|&candidate| {
                reachability[candidate * words + predecessor_prefix / 64]
                    & (1_u64 << (predecessor_prefix % 64))
                    != 0
            })
            .ok_or_else(|| PoaError::InvalidGraph {
                message: "an exact source path has no reachable predecessor".to_owned(),
            })?;
        cursor = predecessor;
        query_prefix = predecessor_prefix;
    }
    if query_prefix != 0 || reversed_path.len() != query.len() {
        return Err(PoaError::InvalidGraph {
            message: "an exact source path does not consume the complete query".to_owned(),
        });
    }
    reversed_path.reverse();
    Ok(Some(DiagonalCertificate {
        path: reversed_path,
        optimum_cost: 0,
    }))
}

fn canonical_zero_diagonal_path<A: Alphabet>(
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    topological: &[PoaNodeId],
    index: &PoastaGraphIndex,
    row_bands: &[ExactAffineRowBand],
    diagonal_costs: &mut [u16],
) -> Result<Vec<StateKey>, PoaError> {
    // Once the optimum is zero, nonzero prefix costs are irrelevant. Reuse
    // the cost table itself as a reverse reachability mask: retain a zero only
    // when the cell can continue through exact matches to a terminal. This
    // reproduces the relaxed certificate's forward successor priority without
    // allocating a second table.
    for cursor in topological
        .iter()
        .rev()
        .map(|node| node.index() + 1)
        .chain(core::iter::once(index.source()))
    {
        let row = row_bands[cursor];
        if row == INVALID_EXACT_AFFINE_ROW_BAND {
            continue;
        }
        for query_prefix in (row.first..=row.last).rev() {
            let cell = row.offset + query_prefix - row.first;
            if diagonal_costs[cell] != 0 {
                continue;
            }
            let viable = if query_prefix == query.len() {
                index.successors(cursor) == [index.sink()]
            } else {
                let query_symbol = query.get(query_prefix).ok_or(PoaError::DimensionOverflow {
                    graph_nodes: graph.node_count(),
                    query_len: query.len(),
                })?;
                index.successors(cursor).iter().copied().any(|successor| {
                    if successor == index.sink() {
                        return false;
                    }
                    let Some(successor_cell) =
                        exact_affine_cell(row_bands[successor], query_prefix + 1)
                    else {
                        return false;
                    };
                    if diagonal_costs[successor_cell] != 0 {
                        return false;
                    }
                    index
                        .node_for_cursor(successor)
                        .and_then(|node| graph.node(node))
                        .is_some_and(|node| node.symbol() == query_symbol)
                })
            };
            if !viable {
                diagonal_costs[cell] = u16::MAX;
            }
        }
    }

    let source_cell =
        exact_affine_cell(row_bands[index.source()], 0).ok_or_else(|| PoaError::InvalidGraph {
            message: "a zero-diagonal layout omits its source cell".to_owned(),
        })?;
    if diagonal_costs[source_cell] != 0 {
        return Err(PoaError::InvalidGraph {
            message: "a zero-cost diagonal optimum has no source-to-sink traceback".to_owned(),
        });
    }
    let mut path = Vec::new();
    try_reserve_exact(
        &mut path,
        query.len(),
        "POASTA canonical zero-diagonal traceback",
    )?;
    let mut cursor = index.source();
    for query_prefix in 0..query.len() {
        let query_symbol = query.get(query_prefix).ok_or(PoaError::DimensionOverflow {
            graph_nodes: graph.node_count(),
            query_len: query.len(),
        })?;
        let successor = index
            .successors(cursor)
            .iter()
            .copied()
            .find(|&candidate| {
                if candidate == index.sink() {
                    return false;
                }
                let Some(cell) = exact_affine_cell(row_bands[candidate], query_prefix + 1) else {
                    return false;
                };
                diagonal_costs[cell] == 0
                    && index
                        .node_for_cursor(candidate)
                        .and_then(|node| graph.node(node))
                        .is_some_and(|node| node.symbol() == query_symbol)
            })
            .ok_or_else(|| PoaError::InvalidGraph {
                message: "a canonical zero-diagonal traceback has no viable successor".to_owned(),
            })?;
        path.push(StateKey::new(
            successor,
            query_prefix + 1,
            AffineState::Match,
        ));
        cursor = successor;
    }
    if index.successors(cursor) != [index.sink()] {
        return Err(PoaError::InvalidGraph {
            message: "a canonical zero-diagonal traceback does not reach the sink".to_owned(),
        });
    }
    Ok(path)
}

#[inline]
fn fill_exact_diagonal_single_predecessor_row(
    row: ExactAffineRowBand,
    predecessor_row: ExactAffineRowBand,
    query_indices: &[u8],
    graph_symbol: u8,
    mismatch: u16,
    diagonal_costs: &mut [u16],
) -> usize {
    if predecessor_row == INVALID_EXACT_AFFINE_ROW_BAND {
        return 0;
    }
    let first = row
        .first
        .max(1)
        .max(predecessor_row.first.saturating_add(1));
    let last = row.last.min(predecessor_row.last.saturating_add(1));
    if first > last {
        return 0;
    }
    let count = last + 1 - first;
    let current_start = first - row.first;
    let predecessor_start = predecessor_row.offset + first - 1 - predecessor_row.first;
    let (prior_costs, current_and_later_costs) = diagonal_costs.split_at_mut(row.offset);
    let current_costs = &mut current_and_later_costs[current_start..current_start + count];
    let predecessor_costs = &prior_costs[predecessor_start..predecessor_start + count];
    let query_symbols = &query_indices[first - 1..last];
    for ((current, &predecessor), &query_symbol) in current_costs
        .iter_mut()
        .zip(predecessor_costs)
        .zip(query_symbols)
    {
        let substitution = if query_symbol == graph_symbol {
            0
        } else {
            mismatch
        };
        *current = predecessor.saturating_add(substitution);
    }
    count
}

#[inline]
fn merge_exact_diagonal_predecessor_row(
    row: ExactAffineRowBand,
    predecessor_row: ExactAffineRowBand,
    diagonal_costs: &mut [u16],
) -> usize {
    if predecessor_row == INVALID_EXACT_AFFINE_ROW_BAND {
        return 0;
    }
    let first = row
        .first
        .max(1)
        .max(predecessor_row.first.saturating_add(1));
    let last = row.last.min(predecessor_row.last.saturating_add(1));
    if first > last {
        return 0;
    }
    let count = last + 1 - first;
    let current_start = first - row.first;
    let predecessor_start = predecessor_row.offset + first - 1 - predecessor_row.first;
    let (prior_costs, current_and_later_costs) = diagonal_costs.split_at_mut(row.offset);
    let current_costs = &mut current_and_later_costs[current_start..current_start + count];
    let predecessor_costs = &prior_costs[predecessor_start..predecessor_start + count];
    for (current, &predecessor) in current_costs.iter_mut().zip(predecessor_costs) {
        *current = (*current).min(predecessor);
    }
    count
}

#[inline]
fn finalize_exact_diagonal_row(
    row: ExactAffineRowBand,
    query_indices: &[u8],
    graph_symbol: u8,
    mismatch: u16,
    diagonal_costs: &mut [u16],
) -> usize {
    let first = row.first.max(1);
    if first > row.last {
        return 0;
    }
    let count = row.last + 1 - first;
    let current_start = row.offset + first - row.first;
    let current_costs = &mut diagonal_costs[current_start..current_start + count];
    let query_symbols = &query_indices[first - 1..row.last];
    for (current, &query_symbol) in current_costs.iter_mut().zip(query_symbols) {
        let substitution = if query_symbol == graph_symbol {
            0
        } else {
            mismatch
        };
        *current = current.saturating_add(substitution);
    }
    count
}

#[allow(clippy::too_many_lines)]
fn exact_diagonal_certificate<A: Alphabet>(
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    topological: &[PoaNodeId],
    index: &PoastaGraphIndex,
    mismatch_cost: Cost,
    cached_prefix_bounds: Option<&[RelaxedLinearRowBounds]>,
    mut workspace: Option<&mut PoastaWorkspace>,
) -> Result<Option<DiagonalCertificateAttempt>, PoaError> {
    if graph.node_count() == 0 {
        return if query.is_empty() {
            Ok(Some(DiagonalCertificateAttempt {
                certificate: DiagonalCertificate {
                    path: Vec::new(),
                    optimum_cost: 0,
                },
                cells: 1,
                payload_bytes: 0,
            }))
        } else {
            Ok(None)
        };
    }
    let maximum_cost = checked_cost_mul(
        query.len(),
        mismatch_cost,
        "the exact diagonal certificate range",
    )?;
    if maximum_cost >= Cost::from(u16::MAX) {
        return Ok(None);
    }
    let mismatch = u16::try_from(mismatch_cost).map_err(|_| PoaError::PoastaCostOverflow {
        context: "the exact diagonal mismatch cost",
    })?;
    let computed_prefix_bounds = if cached_prefix_bounds.is_none() {
        Some(exact_affine_prefix_bounds(
            graph,
            topological,
            index,
            query.len(),
        )?)
    } else {
        None
    };
    let prefix_bounds = cached_prefix_bounds
        .or(computed_prefix_bounds.as_deref())
        .expect("exact prefix bounds are cached or computed");
    let cached_rows = workspace.as_deref_mut().map_or_else(Vec::new, |workspace| {
        core::mem::take(&mut workspace.exact_row_bands)
    });
    let ExactAffineBandLayout {
        rows: row_bands,
        cells: elements,
    } = exact_affine_band_layout(
        cached_rows,
        topological,
        index,
        prefix_bounds,
        query.len(),
        Some(0),
    )?;
    let mut diagonal_costs = match workspace.as_deref_mut() {
        Some(workspace) => prepare_exact_affine_buffer(
            core::mem::take(&mut workspace.exact_match_costs),
            elements,
            "POASTA exact diagonal certificate costs",
        )?,
        None => filled(
            elements,
            u16::MAX,
            "POASTA exact diagonal certificate costs",
        )?,
    };
    diagonal_costs[..elements].fill(u16::MAX);
    diagonal_costs[row_bands[index.source()].offset] = 0;
    for &node_id in topological {
        let cursor = node_id.index() + 1;
        let row = row_bands[cursor];
        if row == INVALID_EXACT_AFFINE_ROW_BAND {
            continue;
        }
        let node = graph
            .node(node_id)
            .ok_or(PoaError::UnknownNode { node: node_id })?;
        let graph_symbol =
            u8::try_from(node.symbol().alphabet_index()).map_err(|_| PoaError::InvalidGraph {
                message: "a graph symbol exceeds the compact alphabet range".to_owned(),
            })?;
        let predecessors = node.predecessors();
        let single_predecessor = match predecessors {
            [] => Some(index.source()),
            [predecessor] => Some(predecessor.index() + 1),
            _ => None,
        };
        if let Some(predecessor) = single_predecessor {
            fill_exact_diagonal_single_predecessor_row(
                row,
                row_bands[predecessor],
                query.as_indices(),
                graph_symbol,
                mismatch,
                &mut diagonal_costs,
            );
        } else {
            for predecessor in predecessors {
                merge_exact_diagonal_predecessor_row(
                    row,
                    row_bands[predecessor.index() + 1],
                    &mut diagonal_costs,
                );
            }
            finalize_exact_diagonal_row(
                row,
                query.as_indices(),
                graph_symbol,
                mismatch,
                &mut diagonal_costs,
            );
        }
    }

    let mut best_goal = None;
    for &node_id in topological {
        let cursor = node_id.index() + 1;
        if index.successors(cursor) != [index.sink()] {
            continue;
        }
        let Some(cell) = exact_affine_cell(row_bands[cursor], query.len()) else {
            continue;
        };
        let cost = diagonal_costs[cell];
        if cost == u16::MAX {
            continue;
        }
        if best_goal.is_none_or(|(best_cost, best_cursor)| {
            cost < best_cost || (cost == best_cost && cursor < best_cursor)
        }) {
            best_goal = Some((cost, cursor));
        }
    }
    let Some((optimum_cost, goal_cursor)) = best_goal else {
        if let Some(workspace) = workspace.as_deref_mut() {
            workspace.exact_match_costs = diagonal_costs;
            workspace.exact_row_bands = row_bands;
        }
        return Ok(None);
    };

    let path = if optimum_cost == 0 && mismatch > 0 {
        canonical_zero_diagonal_path(
            graph,
            query,
            topological,
            index,
            &row_bands,
            &mut diagonal_costs,
        )?
    } else {
        let mut reversed_path = Vec::new();
        try_reserve_exact(
            &mut reversed_path,
            query.len(),
            "POASTA exact diagonal certificate traceback",
        )?;
        let mut cursor = goal_cursor;
        let mut query_prefix = query.len();
        while cursor != index.source() {
            reversed_path.push(StateKey::new(cursor, query_prefix, AffineState::Match));
            let node_id = index
                .node_for_cursor(cursor)
                .ok_or_else(|| PoaError::InvalidGraph {
                    message: "an exact diagonal cursor has no graph node".to_owned(),
                })?;
            let node = graph
                .node(node_id)
                .ok_or(PoaError::UnknownNode { node: node_id })?;
            let current_cost = exact_affine_cost(&diagonal_costs, row_bands[cursor], query_prefix);
            let substitution = if query.as_indices()[query_prefix - 1]
                == u8::try_from(node.symbol().alphabet_index()).map_err(|_| {
                    PoaError::InvalidGraph {
                        message: "a graph symbol exceeds the compact alphabet range".to_owned(),
                    }
                })? {
                0
            } else {
                mismatch
            };
            let predecessor_prefix = query_prefix - 1;
            let predecessor = augmented_predecessor_cursors(node.predecessors(), index.source())
                .find(|&candidate| {
                    exact_affine_cost(&diagonal_costs, row_bands[candidate], predecessor_prefix)
                        .saturating_add(substitution)
                        == current_cost
                })
                .ok_or_else(|| PoaError::InvalidGraph {
                    message: "an exact diagonal traceback has no predecessor".to_owned(),
                })?;
            cursor = predecessor;
            query_prefix = predecessor_prefix;
        }
        if query_prefix != 0 || reversed_path.len() != query.len() {
            return Err(PoaError::InvalidGraph {
                message: "an exact diagonal traceback does not consume the complete query"
                    .to_owned(),
            });
        }
        reversed_path.reverse();
        reversed_path
    };

    let prefix_payload_bytes = checked_bytes(
        prefix_bounds.len(),
        size_of::<RelaxedLinearRowBounds>(),
        "POASTA exact diagonal prefix-bound payload",
    )?;
    let payload_bytes = if workspace.is_some() {
        prefix_payload_bytes
    } else {
        checked_add_usize(
            checked_add_usize(
                checked_bytes(
                    diagonal_costs.capacity(),
                    size_of::<u16>(),
                    "POASTA exact diagonal cost payload",
                )?,
                checked_bytes(
                    row_bands.capacity(),
                    size_of::<ExactAffineRowBand>(),
                    "POASTA exact diagonal row-band payload",
                )?,
                "POASTA exact diagonal table payload",
            )?,
            prefix_payload_bytes,
            "POASTA exact diagonal total payload",
        )?
    };
    if let Some(workspace) = workspace {
        workspace.exact_match_costs = diagonal_costs;
        workspace.exact_row_bands = row_bands;
    }
    Ok(Some(DiagonalCertificateAttempt {
        certificate: DiagonalCertificate {
            path,
            optimum_cost: Cost::from(optimum_cost),
        },
        cells: elements,
        payload_bytes,
    }))
}

impl RelaxedLinearHeuristic {
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn build<A: Alphabet, const DETAILED: bool>(
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
        topological: &[PoaNodeId],
        index: &PoastaGraphIndex,
        costs: PoastaAffineCosts,
        mut workspace: Option<&mut PoastaWorkspace>,
        mut search_timings: Option<&mut PhaseTimingsBuilder>,
        mut work_counters: Option<&mut PoastaSearchWorkCounters>,
    ) -> Result<Option<Self>, PoaError> {
        let setup_started = search_phase_started::<DETAILED>(&search_timings);
        let maximum_steps = index
            .max_remaining_nodes(index.source())
            .and_then(|remaining| remaining.checked_add(query.len()))
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: query.len(),
            })?;
        let maximum_unit =
            nonnegative_score(costs.mismatch()).max(nonnegative_score(costs.gap_extension()));
        let maximum_cost = checked_cost_mul(
            maximum_steps,
            maximum_unit,
            "the relaxed linear heuristic range",
        )?;
        if maximum_cost >= Cost::from(u16::MAX) {
            return Ok(None);
        }
        let mismatch = u16::try_from(nonnegative_score(costs.mismatch())).map_err(|_| {
            PoaError::PoastaCostOverflow {
                context: "the relaxed linear mismatch cost",
            }
        })?;
        let extension = u16::try_from(nonnegative_score(costs.gap_extension())).map_err(|_| {
            PoaError::PoastaCostOverflow {
                context: "the relaxed linear gap cost",
            }
        })?;
        let columns = query
            .len()
            .checked_add(1)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: query.len(),
            })?;
        let cursor_count = index
            .sink()
            .checked_add(1)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: query.len(),
            })?;
        let elements = cursor_count
            .checked_mul(columns)
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: query.len(),
            })?;
        let (mut values, previous_elements) = if let Some(workspace) = workspace.as_deref_mut() {
            let mut values = core::mem::take(&mut workspace.costs);
            let previous_elements = values.len();
            if values.len() < elements {
                let additional = elements - values.len();
                values
                    .try_reserve(additional)
                    .map_err(|_| PoaError::AllocationFailed {
                        buffer: "POASTA relaxed linear heuristic",
                        elements,
                    })?;
                values.resize(elements, 0);
            } else {
                values.truncate(elements);
            }
            (values, previous_elements)
        } else {
            (
                filled(elements, 0_u16, "POASTA relaxed linear heuristic")?,
                0,
            )
        };
        let previous_row_count = if let Some(workspace) = workspace.as_deref_mut() {
            if workspace.query_len != Some(query.len()) || workspace.extension != extension {
                workspace.row_bounds.clear();
            }
            let previous_row_count = workspace.row_bounds.len();
            if previous_row_count < cursor_count {
                try_reserve_exact(
                    &mut workspace.row_bounds,
                    cursor_count - previous_row_count,
                    "POASTA cached relaxed linear rows",
                )?;
                workspace
                    .row_bounds
                    .resize(cursor_count, INVALID_RELAXED_LINEAR_ROW);
            } else {
                workspace.row_bounds.truncate(cursor_count);
            }
            workspace.query_len = Some(query.len());
            workspace.extension = extension;
            previous_row_count
        } else {
            0
        };
        let mut initialized_rows = 0_usize;
        let mut reused_rows = 0_usize;
        for cursor in 0..cursor_count {
            let (minimum, maximum) = if cursor == index.sink() {
                (0, 0)
            } else {
                (
                    index
                        .min_remaining_nodes(cursor)
                        .ok_or_else(|| PoaError::InvalidGraph {
                            message: "a heuristic cursor has no minimum terminal distance"
                                .to_owned(),
                        })?,
                    index
                        .max_remaining_nodes(cursor)
                        .ok_or_else(|| PoaError::InvalidGraph {
                            message: "a heuristic cursor has no maximum terminal distance"
                                .to_owned(),
                        })?,
                )
            };
            let row_start = cursor * columns;
            let bounds = RelaxedLinearRowBounds { minimum, maximum };
            let can_reuse_row = workspace.as_deref_mut().is_some_and(|workspace| {
                let reusable = cursor < previous_row_count
                    && row_start + columns <= previous_elements
                    && workspace.row_bounds[cursor] == bounds;
                workspace.row_bounds[cursor] = bounds;
                reusable
            });
            if can_reuse_row {
                if DETAILED {
                    reused_rows += 1;
                }
                continue;
            }
            if DETAILED {
                initialized_rows += 1;
            }
            if extension == 0 {
                values[row_start..row_start + columns].fill(0);
            } else {
                let left_count = query.len().saturating_sub(maximum);
                let mut length_bound =
                    u16::try_from(left_count).map_err(|_| PoaError::PoastaCostOverflow {
                        context: "the relaxed linear length bound",
                    })? * extension;
                for value in &mut values[row_start..row_start + left_count] {
                    *value = length_bound;
                    length_bound -= extension;
                }
                let right_start = if minimum > 0 {
                    let (right_start, initial_gap) = if minimum > query.len() {
                        (0, minimum - query.len())
                    } else {
                        (query.len() - minimum + 1, 1)
                    };
                    let mut length_bound =
                        u16::try_from(initial_gap).map_err(|_| PoaError::PoastaCostOverflow {
                            context: "the relaxed linear length bound",
                        })? * extension;
                    for value in &mut values[row_start + right_start..row_start + columns] {
                        *value = length_bound;
                        length_bound = length_bound.saturating_add(extension);
                    }
                    right_start
                } else {
                    columns
                };
                values[row_start + left_count..row_start + right_start].fill(0);
            }
        }

        record_search_phase::<DETAILED>(
            search_timings.as_deref_mut(),
            "relaxed_linear_setup",
            setup_started,
        );
        if let Some(counters) = work_counters.as_deref_mut() {
            let initialized_cells =
                initialized_rows
                    .checked_mul(columns)
                    .ok_or(PoaError::ArithmeticOverflow {
                        context: "POASTA detailed relaxed-linear initialized cells",
                    })?;
            counters.record_relaxed_linear_setup(
                initialized_rows,
                reused_rows,
                initialized_cells,
            )?;
        }

        let row_compute_started = search_phase_started::<DETAILED>(&search_timings);
        let mut single_transition_representatives = filled(
            cursor_count,
            usize::MAX,
            "POASTA single-transition row representatives",
        )?;
        let mut multiple_transition_representatives = if DETAILED {
            let mut representatives = HashSet::new();
            representatives
                .try_reserve(cursor_count)
                .map_err(|_| PoaError::AllocationFailed {
                    buffer: "POASTA detailed multiple-transition representatives",
                    elements: cursor_count,
                })?;
            Some(representatives)
        } else {
            None
        };
        let query_indices = query.as_indices();
        let insertion_scan = RelaxedLinearInsertionScan::new(extension);
        let mut computed_cells = 0_usize;
        for cursor in topological
            .iter()
            .rev()
            .map(|node| node.index() + 1)
            .chain(core::iter::once(index.source()))
        {
            let minimum =
                index
                    .min_remaining_nodes(cursor)
                    .ok_or_else(|| PoaError::InvalidGraph {
                        message: "a heuristic cursor has no minimum terminal distance".to_owned(),
                    })?;
            let maximum =
                index
                    .max_remaining_nodes(cursor)
                    .ok_or_else(|| PoaError::InvalidGraph {
                        message: "a heuristic cursor has no maximum terminal distance".to_owned(),
                    })?;
            let bounds = RelaxedLinearRowBounds { minimum, maximum };
            let minimum_query_remaining = minimum.saturating_sub(RELAXED_LINEAR_RADIUS);
            let maximum_query_remaining = maximum
                .saturating_add(RELAXED_LINEAR_RADIUS)
                .min(query.len());
            if minimum_query_remaining > maximum_query_remaining {
                continue;
            }
            let single_transition_representative = match index.successors(cursor) {
                [successor] => {
                    let representative = &mut single_transition_representatives[*successor];
                    if *representative == usize::MAX {
                        *representative = cursor;
                        None
                    } else {
                        debug_assert_eq!(index.min_remaining_nodes(*representative), Some(minimum));
                        debug_assert_eq!(index.max_remaining_nodes(*representative), Some(maximum));
                        Some(*representative)
                    }
                }
                _ => None,
            };
            let first_query_prefix = query.len() - maximum_query_remaining;
            let last_query_prefix = query.len() - minimum_query_remaining;
            let row_cells = last_query_prefix + 1 - first_query_prefix;
            let row_start = cursor * columns;
            let equivalent_multiple_transition = if DETAILED {
                let transitions = index.successors(cursor);
                let has_multiple_node_successors = transitions
                    .iter()
                    .filter(|&&successor| successor != index.sink())
                    .take(2)
                    .count()
                    == 2;
                if has_multiple_node_successors {
                    let representatives = multiple_transition_representatives
                        .as_mut()
                        .expect("detailed relaxed rows have a transition map");
                    !representatives.insert((transitions, bounds))
                } else {
                    false
                }
            } else {
                false
            };
            if let Some(representative) = single_transition_representative {
                let representative_start = representative * columns + first_query_prefix;
                let representative_end = representative_start + row_cells;
                values.copy_within(
                    representative_start..representative_end,
                    row_start + first_query_prefix,
                );
                if DETAILED {
                    if let Some(counters) = work_counters.as_deref_mut() {
                        let [successor] = index.successors(cursor) else {
                            unreachable!("a representative requires one transition");
                        };
                        let node_successors = usize::from(*successor != index.sink());
                        let replaceable_successor_cell_visits = row_cells
                            .checked_mul(node_successors)
                            .ok_or(PoaError::ArithmeticOverflow {
                                context: "POASTA detailed equivalent successor-cell visits",
                            })?;
                        let replaceable_insertion_cells =
                            row_cells.saturating_sub(usize::from(last_query_prefix == query.len()));
                        counters.record_relaxed_linear_copied_row(
                            row_cells,
                            replaceable_insertion_cells,
                            node_successors,
                            replaceable_successor_cell_visits,
                        )?;
                    }
                }
                computed_cells =
                    computed_cells
                        .checked_add(row_cells)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "POASTA relaxed linear heuristic cell count",
                        })?;
                continue;
            }
            let mut initialized = false;
            let mut node_successors = 0_usize;
            let mut node_successor_symbols = [0_u64; 4];
            for &successor in index.successors(cursor) {
                if successor == index.sink() {
                    if !initialized {
                        values[row_start + first_query_prefix..=row_start + last_query_prefix]
                            .fill(u16::MAX);
                        initialized = true;
                    }
                    if (first_query_prefix..=last_query_prefix).contains(&query.len()) {
                        values[row_start + query.len()] = 0;
                    }
                    continue;
                }
                if DETAILED {
                    node_successors =
                        node_successors
                            .checked_add(1)
                            .ok_or(PoaError::ArithmeticOverflow {
                                context: "POASTA detailed relaxed-linear successor count",
                            })?;
                }
                let node =
                    index
                        .node_for_cursor(successor)
                        .ok_or_else(|| PoaError::InvalidGraph {
                            message: "a relaxed-heuristic cursor has no graph node".to_owned(),
                        })?;
                let graph_symbol = u8::try_from(
                    graph
                        .node(node)
                        .ok_or(PoaError::UnknownNode { node })?
                        .symbol()
                        .alphabet_index(),
                )
                .map_err(|_| PoaError::InvalidGraph {
                    message: "a graph symbol exceeds the compact alphabet range".to_owned(),
                })?;
                if DETAILED {
                    let symbol = usize::from(graph_symbol);
                    node_successor_symbols[symbol / 64] |= 1_u64 << (symbol % 64);
                }
                let (row, successor_row) =
                    relaxed_row_pair(&mut values, columns, cursor, successor);
                if initialized {
                    merge_relaxed_successor::<false>(
                        row,
                        successor_row,
                        query_indices,
                        graph_symbol,
                        mismatch,
                        extension,
                        first_query_prefix,
                        last_query_prefix,
                    );
                } else {
                    merge_relaxed_successor::<true>(
                        row,
                        successor_row,
                        query_indices,
                        graph_symbol,
                        mismatch,
                        extension,
                        first_query_prefix,
                        last_query_prefix,
                    );
                    initialized = true;
                }
            }
            if !initialized {
                values[row_start + first_query_prefix..=row_start + last_query_prefix]
                    .fill(u16::MAX);
            }
            let row = &mut values[row_start..row_start + columns];
            // Every successor value already includes its length-distance lower
            // bound. A deletion or insertion adds one extension, while a
            // diagonal consumes one graph and one query symbol together. By
            // reverse-topological and reverse-query induction, the recurrence
            // therefore remains at least the current row's length bound; the
            // explicit per-cell maximum would only repeat that same constraint.
            let refinement_end = if last_query_prefix == query.len() {
                if row[last_query_prefix] == u16::MAX {
                    return Err(PoaError::InvalidGraph {
                        message: "a relaxed-heuristic cell cannot reach the sink".to_owned(),
                    });
                }
                last_query_prefix.checked_sub(1)
            } else {
                Some(last_query_prefix)
            };
            let mut insertion = refinement_end
                .and_then(|end| row.get(end + 1))
                .copied()
                .map_or(u16::MAX, |best| best + extension);
            let mut insertion_cells = 0_usize;
            let mut insertion_avx2_cells = 0_usize;
            if let Some(refinement_end) = refinement_end.filter(|&end| first_query_prefix <= end) {
                insertion_cells = refinement_end + 1 - first_query_prefix;
                let (processed, next_insertion, reachable) = insertion_scan.fill(
                    &mut row[first_query_prefix..=refinement_end],
                    extension,
                    insertion,
                );
                if !reachable {
                    return Err(PoaError::InvalidGraph {
                        message: "a relaxed-heuristic cell cannot reach the sink".to_owned(),
                    });
                }
                insertion_avx2_cells = processed;
                insertion = next_insertion;
                let vector_start = refinement_end + 1 - processed;
                for query_prefix in (first_query_prefix + 1..vector_start).rev() {
                    let best = row[query_prefix].min(insertion);
                    if best == u16::MAX {
                        return Err(PoaError::InvalidGraph {
                            message: "a relaxed-heuristic cell cannot reach the sink".to_owned(),
                        });
                    }
                    row[query_prefix] = best;
                    insertion = best + extension;
                }
                if first_query_prefix < vector_start {
                    let best = row[first_query_prefix].min(insertion);
                    if best == u16::MAX {
                        return Err(PoaError::InvalidGraph {
                            message: "a relaxed-heuristic cell cannot reach the sink".to_owned(),
                        });
                    }
                    row[first_query_prefix] = best;
                }
            }
            if DETAILED {
                if let Some(counters) = work_counters.as_deref_mut() {
                    let successor_cell_visits = row_cells.checked_mul(node_successors).ok_or(
                        PoaError::ArithmeticOverflow {
                            context: "POASTA detailed relaxed-linear successor-cell visits",
                        },
                    )?;
                    let distinct_node_successor_symbols = node_successor_symbols
                        .iter()
                        .map(|word| word.count_ones() as usize)
                        .sum();
                    counters.record_relaxed_linear_row(RelaxedLinearRowWork {
                        row_cells,
                        insertion_cells,
                        insertion_avx2_cells,
                        node_successors,
                        distinct_node_successor_symbols,
                        successor_cell_visits,
                        equivalent_multiple_transition,
                    })?;
                }
            }
            computed_cells = computed_cells
                .checked_add(last_query_prefix + 1 - first_query_prefix)
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "POASTA relaxed linear heuristic cell count",
                })?;
        }
        record_search_phase::<DETAILED>(search_timings, "relaxed_linear_rows", row_compute_started);
        let construction_scratch_bytes = checked_bytes(
            single_transition_representatives.capacity(),
            size_of::<usize>(),
            "POASTA single-transition row representative payload",
        )?;
        Ok(Some(Self {
            columns,
            computed_cells,
            construction_scratch_bytes,
            costs: values,
        }))
    }

    fn lower_bound(&self, key: StateKey) -> Cost {
        Cost::from(self.costs[key.cursor * self.columns + key.query_prefix()])
    }

    fn diagonal_certificate<A: Alphabet>(
        &self,
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
        index: &PoastaGraphIndex,
        mismatch: Cost,
    ) -> Result<Option<DiagonalCertificate>, PoaError> {
        let mut path = Vec::new();
        try_reserve_exact(
            &mut path,
            query.len(),
            "POASTA relaxed-linear certificate path",
        )?;
        let source_bound = self.lower_bound(StateKey::new(index.source(), 0, AffineState::Match));
        let mut cursor = index.source();
        let mut path_cost = 0;
        for query_prefix in 0..query.len() {
            let current_bound =
                self.lower_bound(StateKey::new(cursor, query_prefix, AffineState::Match));
            let query_symbol = query.get(query_prefix).ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len: query.len(),
            })?;
            let remaining = query.len() - query_prefix - 1;
            let mut selected = None;
            for &successor in index.successors(cursor) {
                if successor == index.sink()
                    || !index
                        .min_remaining_nodes(successor)
                        .zip(index.max_remaining_nodes(successor))
                        .is_some_and(|(minimum, maximum)| {
                            minimum <= remaining && remaining <= maximum
                        })
                {
                    continue;
                }
                let node =
                    index
                        .node_for_cursor(successor)
                        .ok_or_else(|| PoaError::InvalidGraph {
                            message: "a certificate cursor has no graph node".to_owned(),
                        })?;
                let transition = if graph
                    .node(node)
                    .ok_or(PoaError::UnknownNode { node })?
                    .symbol()
                    == query_symbol
                {
                    0
                } else {
                    mismatch
                };
                let successor_bound = self.lower_bound(StateKey::new(
                    successor,
                    query_prefix + 1,
                    AffineState::Match,
                ));
                if transition
                    .checked_add(successor_bound)
                    .is_some_and(|candidate| candidate == current_bound)
                {
                    selected = Some((successor, transition));
                    break;
                }
            }
            let Some((successor, transition)) = selected else {
                return Ok(None);
            };
            path_cost =
                checked_cost_add(path_cost, transition, "a relaxed-linear certificate cost")?;
            path.push(StateKey::new(
                successor,
                query_prefix + 1,
                AffineState::Match,
            ));
            cursor = successor;
        }
        if index.successors(cursor) != [index.sink()] || path_cost != source_bound {
            return Ok(None);
        }
        Ok(Some(DiagonalCertificate {
            path,
            optimum_cost: path_cost,
        }))
    }

    const fn cell_count(&self) -> usize {
        self.computed_cells
    }

    fn modeled_payload_bytes(&self) -> Result<usize, PoaError> {
        let table_bytes = checked_bytes(
            self.costs.capacity(),
            size_of::<u16>(),
            "POASTA relaxed linear heuristic payload",
        )?;
        checked_add_usize(
            table_bytes,
            self.construction_scratch_bytes,
            "POASTA relaxed linear construction peak payload",
        )
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct CertifiedStep {
    key: StateKey,
    transition_cost: Cost,
    tie_transition: u8,
}

#[derive(Debug, Eq, PartialEq)]
struct ExactAffineCertificate {
    steps: Vec<CertifiedStep>,
    optimum_cost: Cost,
}

struct ExactAffineCertificateAttempt {
    certificate: Option<ExactAffineCertificate>,
    cells: usize,
    payload_bytes: usize,
}

fn consider_exact_affine_predecessor(
    selected: &mut Option<(StateKey, Cost, u8)>,
    predecessor: StateKey,
    transition_cost: Cost,
    tie_transition: u8,
) {
    if selected
        .as_ref()
        .is_none_or(|(selected_predecessor, _, selected_tie)| {
            (tie_transition, predecessor.cursor) < (*selected_tie, selected_predecessor.cursor)
        })
    {
        *selected = Some((predecessor, transition_cost, tie_transition));
    }
}

fn augmented_predecessor_cursors(
    predecessors: &[PoaNodeId],
    source: usize,
) -> impl Iterator<Item = usize> + '_ {
    predecessors
        .iter()
        .map(|node| node.index() + 1)
        .chain(predecessors.is_empty().then_some(source))
}

fn exact_affine_gap_run_cost(
    length: usize,
    opening: Cost,
    extension: Cost,
) -> Result<Cost, PoaError> {
    if length == 0 {
        return Ok(0);
    }
    checked_cost_add(
        opening,
        checked_cost_mul(
            length - 1,
            extension,
            "an exact-affine upper-bound gap extension",
        )?,
        "an exact-affine upper-bound gap run",
    )
}

fn exact_affine_upper_bound<A: Alphabet>(
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    index: &PoastaGraphIndex,
    mismatch: Cost,
    opening: Cost,
    extension: Cost,
) -> Result<Cost, PoaError> {
    if graph.node_count() == 0 {
        return exact_affine_gap_run_cost(query.len(), opening, extension);
    }
    let minimum_path_nodes =
        index
            .min_remaining_nodes(index.source())
            .ok_or_else(|| PoaError::InvalidGraph {
                message: "a non-empty graph has no source-to-terminal path length".to_owned(),
            })?;
    let mut best = checked_cost_add(
        exact_affine_gap_run_cost(minimum_path_nodes, opening, extension)?,
        exact_affine_gap_run_cost(query.len(), opening, extension)?,
        "an exact-affine guaranteed opposite-gap upper bound",
    )?;
    let query_indices = query.as_indices();
    for path in graph.paths() {
        let Some((&first, &last)) = path.nodes().first().zip(path.nodes().last()) else {
            continue;
        };
        if !graph
            .node(first)
            .ok_or(PoaError::UnknownNode { node: first })?
            .predecessors()
            .is_empty()
            || !graph
                .node(last)
                .ok_or(PoaError::UnknownNode { node: last })?
                .successors()
                .is_empty()
        {
            continue;
        }
        let compared = path.source_symbols().len().min(query_indices.len());
        let mismatches = path.source_symbols()[..compared]
            .iter()
            .zip(&query_indices[..compared])
            .filter(|(graph_symbol, query_symbol)| graph_symbol != query_symbol)
            .count();
        let diagonal =
            checked_cost_mul(mismatches, mismatch, "an exact-affine diagonal upper bound")?;
        let length_difference = path.nodes().len().abs_diff(query.len());
        let diagonal_then_gap = checked_cost_add(
            diagonal,
            exact_affine_gap_run_cost(length_difference, opening, extension)?,
            "an exact-affine diagonal-plus-gap upper bound",
        )?;
        let opposite_gaps = checked_cost_add(
            exact_affine_gap_run_cost(path.nodes().len(), opening, extension)?,
            exact_affine_gap_run_cost(query.len(), opening, extension)?,
            "an exact-affine opposite-gap upper bound",
        )?;
        let candidate = diagonal_then_gap.min(opposite_gaps);
        best = best.min(candidate);
    }
    Ok(best)
}

fn exact_affine_prefix_bounds<A: Alphabet>(
    graph: &PoaGraph<A>,
    topological: &[PoaNodeId],
    index: &PoastaGraphIndex,
    query_len: usize,
) -> Result<Vec<RelaxedLinearRowBounds>, PoaError> {
    let mut bounds = filled(
        index.sink(),
        INVALID_RELAXED_LINEAR_ROW,
        "POASTA exact-affine prefix bounds",
    )?;
    bounds[index.source()] = RelaxedLinearRowBounds {
        minimum: 0,
        maximum: 0,
    };
    for &node_id in topological {
        let node = graph
            .node(node_id)
            .ok_or(PoaError::UnknownNode { node: node_id })?;
        let mut minimum = usize::MAX;
        let mut maximum = 0_usize;
        for predecessor in augmented_predecessor_cursors(node.predecessors(), index.source()) {
            let predecessor_bounds = bounds[predecessor];
            if predecessor_bounds == INVALID_RELAXED_LINEAR_ROW {
                return Err(PoaError::InvalidGraph {
                    message: "an exact-affine predecessor has no prefix bounds".to_owned(),
                });
            }
            minimum = minimum.min(predecessor_bounds.minimum);
            maximum = maximum.max(predecessor_bounds.maximum);
        }
        bounds[node_id.index() + 1] = RelaxedLinearRowBounds {
            minimum: minimum.checked_add(1).ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len,
            })?,
            maximum: maximum.checked_add(1).ok_or(PoaError::DimensionOverflow {
                graph_nodes: graph.node_count(),
                query_len,
            })?,
        };
    }
    Ok(bounds)
}

#[cfg(test)]
const fn distance_to_range(value: usize, bounds: RelaxedLinearRowBounds) -> usize {
    if value < bounds.minimum {
        bounds.minimum - value
    } else if value > bounds.maximum {
        value - bounds.maximum
    } else {
        0
    }
}

fn exact_affine_row_band(
    cursor: usize,
    query_len: usize,
    index: &PoastaGraphIndex,
    prefix_bounds: &[RelaxedLinearRowBounds],
    radius: Option<usize>,
) -> Result<Option<(usize, usize)>, PoaError> {
    let Some(radius) = radius else {
        return Ok(Some((0, query_len)));
    };
    let prefix = prefix_bounds[cursor];
    let suffix = RelaxedLinearRowBounds {
        minimum: index
            .min_remaining_nodes(cursor)
            .ok_or_else(|| PoaError::InvalidGraph {
                message: "an exact-affine cursor has no minimum terminal distance".to_owned(),
            })?,
        maximum: index
            .max_remaining_nodes(cursor)
            .ok_or_else(|| PoaError::InvalidGraph {
                message: "an exact-affine cursor has no maximum terminal distance".to_owned(),
            })?,
    };
    Ok(exact_affine_distance_band(
        prefix, suffix, query_len, radius,
    ))
}

fn exact_affine_distance_band(
    prefix: RelaxedLinearRowBounds,
    suffix: RelaxedLinearRowBounds,
    query_len: usize,
    radius: usize,
) -> Option<(usize, usize)> {
    // Reflect the suffix interval into query-prefix coordinates. Adding its
    // maximum first keeps every coordinate nonnegative even when all suffix
    // lengths exceed the query, while u128 prevents the translated endpoints
    // and radius sums from overflowing usize.
    let offset = suffix.maximum as u128;
    let query_end = query_len as u128 + offset;
    let prefix_lower = prefix.minimum as u128 + offset;
    let prefix_upper = prefix.maximum as u128 + offset;
    let suffix_lower = query_len as u128;
    let suffix_upper = query_len as u128 + (suffix.maximum - suffix.minimum) as u128;
    let radius = radius as u128;

    let minimum_distance = if prefix_upper < suffix_lower {
        suffix_lower - prefix_upper
    } else {
        prefix_lower.saturating_sub(suffix_upper)
    };
    if radius < minimum_distance {
        return None;
    }

    let lower_first = prefix_lower.min(suffix_lower);
    let lower_second = prefix_lower.max(suffix_lower);
    let first = if radius < lower_second - lower_first {
        lower_second - radius
    } else {
        (lower_first + lower_second)
            .saturating_sub(radius)
            .div_ceil(2)
    }
    .max(offset);

    let upper_first = prefix_upper.min(suffix_upper);
    let upper_second = prefix_upper.max(suffix_upper);
    let last = if radius < upper_second - upper_first {
        upper_first + radius
    } else {
        (upper_first + upper_second + radius) / 2
    }
    .min(query_end);
    if first > last {
        return None;
    }
    Some((
        usize::try_from(first - offset).expect("a clamped query prefix fits usize"),
        usize::try_from(last - offset).expect("a clamped query prefix fits usize"),
    ))
}

struct ExactAffineBandLayout {
    rows: Vec<ExactAffineRowBand>,
    cells: usize,
}

fn exact_affine_band_layout(
    mut rows: Vec<ExactAffineRowBand>,
    topological: &[PoaNodeId],
    index: &PoastaGraphIndex,
    prefix_bounds: &[RelaxedLinearRowBounds],
    query_len: usize,
    radius: Option<usize>,
) -> Result<ExactAffineBandLayout, PoaError> {
    let row_count = index.sink();
    if rows.len() < row_count {
        let additional = row_count - rows.len();
        try_reserve_exact(
            &mut rows,
            additional,
            "POASTA compact exact-affine row bands",
        )?;
        rows.resize(row_count, INVALID_EXACT_AFFINE_ROW_BAND);
    } else {
        rows.truncate(row_count);
    }
    rows.fill(INVALID_EXACT_AFFINE_ROW_BAND);
    let source = index.source();
    rows[source] = ExactAffineRowBand {
        first: 0,
        last: query_len,
        offset: 0,
    };
    let mut cells = query_len
        .checked_add(1)
        .ok_or(PoaError::ArithmeticOverflow {
            context: "POASTA compact exact-affine source row",
        })?;
    for &node_id in topological {
        let cursor = node_id.index() + 1;
        let Some((first, last)) =
            exact_affine_row_band(cursor, query_len, index, prefix_bounds, radius)?
        else {
            continue;
        };
        rows[cursor] = ExactAffineRowBand {
            first,
            last,
            offset: cells,
        };
        cells = cells
            .checked_add(last + 1 - first)
            .ok_or(PoaError::ArithmeticOverflow {
                context: "POASTA compact exact-affine band cells",
            })?;
    }
    Ok(ExactAffineBandLayout { rows, cells })
}

#[inline]
const fn exact_affine_cell(row: ExactAffineRowBand, query_prefix: usize) -> Option<usize> {
    if query_prefix < row.first || query_prefix > row.last {
        None
    } else {
        Some(row.offset + query_prefix - row.first)
    }
}

#[inline]
const fn exact_affine_row_cells(row: ExactAffineRowBand) -> usize {
    row.last + 1 - row.first
}

#[inline]
fn exact_affine_single_fast_cells(
    row: ExactAffineRowBand,
    predecessor_row: ExactAffineRowBand,
) -> usize {
    if predecessor_row == INVALID_EXACT_AFFINE_ROW_BAND {
        return 0;
    }
    let first = row.first.max(predecessor_row.first.saturating_add(1));
    let last = row.last.min(predecessor_row.last);
    last.checked_sub(first).map_or(0, |distance| distance + 1)
}

#[inline]
fn exact_affine_two_fast_cells(
    row: ExactAffineRowBand,
    first_predecessor_row: ExactAffineRowBand,
    second_predecessor_row: ExactAffineRowBand,
) -> usize {
    if first_predecessor_row == INVALID_EXACT_AFFINE_ROW_BAND
        || second_predecessor_row == INVALID_EXACT_AFFINE_ROW_BAND
    {
        return 0;
    }
    let first = row
        .first
        .max(first_predecessor_row.first.saturating_add(1))
        .max(second_predecessor_row.first.saturating_add(1));
    let last = row
        .last
        .min(first_predecessor_row.last)
        .min(second_predecessor_row.last);
    last.checked_sub(first).map_or(0, |distance| distance + 1)
}

fn exact_affine_many_predecessor_cell_visits(
    row: ExactAffineRowBand,
    predecessors: &[PoaNodeId],
    row_bands: &[ExactAffineRowBand],
) -> Result<usize, PoaError> {
    predecessors.iter().try_fold(0_usize, |total, predecessor| {
        let predecessor_row = row_bands[predecessor.index() + 1];
        if predecessor_row == INVALID_EXACT_AFFINE_ROW_BAND {
            return Ok(total);
        }
        let first = row.first.max(predecessor_row.first);
        let last = row.last.min(predecessor_row.last);
        let cells = last.checked_sub(first).map_or(0, |distance| distance + 1);
        total
            .checked_add(cells)
            .ok_or(PoaError::ArithmeticOverflow {
                context: "POASTA detailed many-predecessor cell visits",
            })
    })
}

#[inline]
fn exact_affine_cost(costs: &[u16], row: ExactAffineRowBand, query_prefix: usize) -> u16 {
    exact_affine_cell(row, query_prefix).map_or(u16::MAX, |cell| costs[cell])
}

fn exact_affine_insertion_costs_at(
    match_costs: &[u16],
    row: ExactAffineRowBand,
    query_prefix: usize,
    opening: u16,
    extension: u16,
) -> (u16, u16) {
    if row == INVALID_EXACT_AFFINE_ROW_BAND || query_prefix <= row.first || query_prefix > row.last
    {
        return (u16::MAX, u16::MAX);
    }
    let mut insertion = u16::MAX;
    let mut predecessor_insertion = u16::MAX;
    for prefix in row.first + 1..=query_prefix {
        predecessor_insertion = insertion;
        insertion = exact_affine_cost(match_costs, row, prefix - 1)
            .saturating_add(opening)
            .min(insertion.saturating_add(extension));
    }
    (insertion, predecessor_insertion)
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn fill_exact_affine_single_predecessor_cell(
    query_prefix: usize,
    row: ExactAffineRowBand,
    predecessor_row: ExactAffineRowBand,
    query_indices: &[u8],
    graph_symbol: u8,
    mismatch: u16,
    opening: u16,
    extension: u16,
    match_costs: &mut [u16],
    deletion_costs: &mut [u16],
    insertion: &mut u16,
) {
    let cell = row.offset + query_prefix - row.first;
    let predecessor_match = exact_affine_cost(match_costs, predecessor_row, query_prefix);
    let predecessor_deletion = exact_affine_cost(deletion_costs, predecessor_row, query_prefix);
    let predecessor_diagonal = if query_prefix == 0 {
        u16::MAX
    } else {
        exact_affine_cost(match_costs, predecessor_row, query_prefix - 1)
    };
    let deletion = predecessor_match
        .saturating_add(opening)
        .min(predecessor_deletion.saturating_add(extension));
    deletion_costs[cell] = deletion;
    let diagonal = if query_prefix == 0 {
        u16::MAX
    } else {
        let substitution = if query_indices[query_prefix - 1] == graph_symbol {
            0
        } else {
            mismatch
        };
        predecessor_diagonal.saturating_add(substitution)
    };
    *insertion = if query_prefix > row.first {
        match_costs[cell - 1]
            .saturating_add(opening)
            .min(insertion.saturating_add(extension))
    } else {
        u16::MAX
    };
    match_costs[cell] = diagonal.min(deletion).min(*insertion);
}

#[derive(Clone, Copy, Default)]
struct RelaxedLinearInsertionScan {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    v3: Option<pulp::x86::V3>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    extension_steps: [u16; 4],
}

impl RelaxedLinearInsertionScan {
    #[inline]
    fn new(extension: u16) -> Self {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            let Some(two) = extension.checked_mul(2) else {
                return Self::default();
            };
            let Some(four) = extension.checked_mul(4) else {
                return Self::default();
            };
            let Some(eight) = extension.checked_mul(8) else {
                return Self::default();
            };
            Self {
                v3: pulp::x86::V3::try_new(),
                extension_steps: [extension, two, four, eight],
            }
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        {
            let _ = extension;
            Self::default()
        }
    }

    #[inline]
    fn fill(self, current: &mut [u16], extension: u16, insertion: u16) -> (usize, u16, bool) {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if let Some(simd) = self.v3 {
            return fill_relaxed_linear_insertions_v3(
                simd,
                current,
                extension,
                self.extension_steps,
                insertion,
            );
        }
        let _ = (current, extension);
        (0, insertion, true)
    }
}

#[derive(Clone, Copy, Default)]
struct ExactAffineInsertionScan {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    v3: Option<pulp::x86::V3>,
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    extension_steps: [u16; 4],
}

impl ExactAffineInsertionScan {
    #[inline]
    fn new(extension: u16) -> Self {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            let Some(two) = extension.checked_mul(2) else {
                return Self::default();
            };
            let Some(four) = extension.checked_mul(4) else {
                return Self::default();
            };
            let Some(eight) = extension.checked_mul(8) else {
                return Self::default();
            };
            Self {
                v3: pulp::x86::V3::try_new(),
                extension_steps: [extension, two, four, eight],
            }
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        {
            let _ = extension;
            Self::default()
        }
    }

    #[inline]
    fn fill(
        self,
        current_matches: &mut [u16],
        opening: u16,
        extension: u16,
        opening_base: u16,
        insertion: u16,
    ) -> (usize, u16, u16) {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if let Some(simd) = self.v3 {
            return fill_exact_affine_insertions_v3(
                simd,
                current_matches,
                opening,
                extension,
                self.extension_steps,
                opening_base,
                insertion,
            );
        }
        let _ = (current_matches, opening, extension);
        (0, opening_base, insertion)
    }

    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn fill_many_bases(
        self,
        predecessor_matches: &mut [u16],
        predecessor_deletions: &mut [u16],
        query_symbols: &[u8],
        graph_symbol: u8,
        mismatch: u16,
        opening: u16,
        extension: u16,
    ) -> usize {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        if let Some(simd) = self.v3 {
            return fill_exact_affine_many_bases_v3(
                simd,
                predecessor_matches,
                predecessor_deletions,
                query_symbols,
                graph_symbol,
                mismatch,
                opening,
                extension,
            );
        }
        let _ = (
            predecessor_matches,
            predecessor_deletions,
            query_symbols,
            graph_symbol,
            mismatch,
            opening,
            extension,
        );
        0
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
fn reverse_u16x16(simd: pulp::x86::V3, values: pulp::u16x16) -> pulp::u16x16 {
    let byte_indices = pulp::u8x32(
        14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7,
        4, 5, 2, 3, 0, 1,
    );
    let reversed_halves = simd
        .avx2
        ._mm256_shuffle_epi8(pulp::cast(values), pulp::cast(byte_indices));
    pulp::cast(
        simd.avx2
            ._mm256_permute2x128_si256::<0x01>(reversed_halves, reversed_halves),
    )
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
fn fill_relaxed_linear_insertions_v3(
    simd: pulp::x86::V3,
    current: &mut [u16],
    extension: u16,
    extension_steps: [u16; 4],
    mut insertion: u16,
) -> (usize, u16, bool) {
    let extension_vectors = extension_steps.map(|step| simd.splat_u16s(step));
    let maximum = simd.splat_u16s(u16::MAX);
    let (_, vectors) = <pulp::x86::V3 as Simd>::as_mut_rsimd_u16s(current);
    let processed = vectors.len() * <pulp::x86::V3 as Simd>::U16_LANES;
    let mut reachable = true;
    for current in vectors.iter_mut().rev() {
        let mut best = reverse_u16x16(simd, *current);
        let first = exact_affine_first_u16(simd, best).min(insertion);
        let first_bits = i16::from_ne_bytes(first.to_ne_bytes());
        best = pulp::cast(
            simd.avx
                ._mm256_insert_epi16::<0>(pulp::cast(best), first_bits),
        );
        best = exact_affine_prefix_step_u16::<1>(simd, best, extension_vectors[0]);
        best = exact_affine_prefix_step_u16::<2>(simd, best, extension_vectors[1]);
        best = exact_affine_prefix_step_u16::<4>(simd, best, extension_vectors[2]);
        best = exact_affine_prefix_step_u16::<8>(simd, best, extension_vectors[3]);
        let unreachable = simd
            .avx2
            ._mm256_cmpeq_epi16(pulp::cast(best), pulp::cast(maximum));
        reachable &= simd.avx2._mm256_movemask_epi8(unreachable) == 0;
        insertion = exact_affine_last_u16(simd, best).saturating_add(extension);
        *current = reverse_u16x16(simd, best);
    }
    (processed, insertion, reachable)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
fn exact_affine_shift_right_u16(
    simd: pulp::x86::V3,
    values: pulp::u16x16,
    amount: usize,
    fill: u16,
) -> pulp::u16x16 {
    let values = pulp::cast(values);
    let fill = pulp::cast(simd.splat_u16s(fill));
    let previous_halves = simd.avx2._mm256_permute2x128_si256::<0x20>(fill, values);
    pulp::cast(match amount {
        1 => simd.avx2._mm256_alignr_epi8::<14>(values, previous_halves),
        2 => simd.avx2._mm256_alignr_epi8::<12>(values, previous_halves),
        4 => simd.avx2._mm256_alignr_epi8::<8>(values, previous_halves),
        8 => previous_halves,
        _ => unreachable!("the exact-affine prefix scan uses only offsets 1, 2, 4, and 8"),
    })
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
fn exact_affine_saturating_add_u16(
    simd: pulp::x86::V3,
    left: pulp::u16x16,
    right: pulp::u16x16,
) -> pulp::u16x16 {
    pulp::cast(
        simd.avx2
            ._mm256_adds_epu16(pulp::cast(left), pulp::cast(right)),
    )
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
fn exact_affine_prefix_step_u16<const OFFSET: usize>(
    simd: pulp::x86::V3,
    insertion: pulp::u16x16,
    extension_step: pulp::u16x16,
) -> pulp::u16x16 {
    let shifted = exact_affine_shift_right_u16(simd, insertion, OFFSET, u16::MAX);
    let candidate = exact_affine_saturating_add_u16(simd, shifted, extension_step);
    simd.min_u16s(insertion, candidate)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
fn fill_exact_affine_insertions_v3(
    simd: pulp::x86::V3,
    current_matches: &mut [u16],
    opening: u16,
    extension: u16,
    extension_steps: [u16; 4],
    mut opening_base: u16,
    mut insertion: u16,
) -> (usize, u16, u16) {
    let opening_vector = simd.splat_u16s(opening);
    let extension_vectors = extension_steps.map(|step| simd.splat_u16s(step));
    let (vectors, _) = <pulp::x86::V3 as Simd>::as_mut_simd_u16s(current_matches);
    let processed = vectors.len() * <pulp::x86::V3 as Simd>::U16_LANES;
    for current in vectors {
        let base = *current;
        let shifted_base = exact_affine_shift_right_u16(simd, base, 1, opening_base);
        let mut next_insertion =
            exact_affine_saturating_add_u16(simd, shifted_base, opening_vector);
        let first = opening_base
            .saturating_add(opening)
            .min(insertion.saturating_add(extension));
        let first_bits = i16::from_ne_bytes(first.to_ne_bytes());
        next_insertion = pulp::cast(
            simd.avx
                ._mm256_insert_epi16::<0>(pulp::cast(next_insertion), first_bits),
        );
        next_insertion =
            exact_affine_prefix_step_u16::<1>(simd, next_insertion, extension_vectors[0]);
        next_insertion =
            exact_affine_prefix_step_u16::<2>(simd, next_insertion, extension_vectors[1]);
        next_insertion =
            exact_affine_prefix_step_u16::<4>(simd, next_insertion, extension_vectors[2]);
        next_insertion =
            exact_affine_prefix_step_u16::<8>(simd, next_insertion, extension_vectors[3]);
        *current = simd.min_u16s(base, next_insertion);

        opening_base = exact_affine_last_u16(simd, base);
        insertion = exact_affine_last_u16(simd, next_insertion);
    }
    (processed, opening_base, insertion)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
#[allow(clippy::too_many_arguments)]
fn fill_exact_affine_many_bases_v3(
    simd: pulp::x86::V3,
    predecessor_matches: &mut [u16],
    predecessor_deletions: &mut [u16],
    query_symbols: &[u8],
    graph_symbol: u8,
    mismatch: u16,
    opening: u16,
    extension: u16,
) -> usize {
    debug_assert_eq!(predecessor_matches.len(), predecessor_deletions.len());
    debug_assert_eq!(
        predecessor_matches.len().saturating_sub(1),
        query_symbols.len()
    );
    let processed = predecessor_matches.len().saturating_sub(1) / 16 * 16;
    let vector_start = predecessor_matches.len() - processed;
    let opening = simd.splat_u16s(opening);
    let extension = simd.splat_u16s(extension);
    let mismatch = simd.splat_u16s(mismatch);
    let graph_symbol = i8::from_ne_bytes([graph_symbol]);
    for start in (vector_start..predecessor_matches.len()).step_by(16).rev() {
        let end = start + 16;
        let predecessor_match: pulp::u16x16 = pulp::cast(
            <[u16; 16]>::try_from(&predecessor_matches[start..end])
                .expect("an exact-affine AVX2 chunk has sixteen match cells"),
        );
        let predecessor_diagonal: pulp::u16x16 = pulp::cast(
            <[u16; 16]>::try_from(&predecessor_matches[start - 1..end - 1])
                .expect("an exact-affine AVX2 chunk has sixteen diagonal cells"),
        );
        let predecessor_deletion: pulp::u16x16 = pulp::cast(
            <[u16; 16]>::try_from(&predecessor_deletions[start..end])
                .expect("an exact-affine AVX2 chunk has sixteen deletion cells"),
        );
        let query_symbol: pulp::u8x16 = pulp::cast(
            <[u8; 16]>::try_from(&query_symbols[start - 1..end - 1])
                .expect("an exact-affine AVX2 chunk has sixteen query symbols"),
        );
        let equal = simd.sse2._mm_cmpeq_epi8(
            pulp::cast(query_symbol),
            simd.sse2._mm_set1_epi8(graph_symbol),
        );
        let equal = simd.avx2._mm256_cvtepi8_epi16(equal);
        let substitution = pulp::cast(simd.avx2._mm256_andnot_si256(equal, pulp::cast(mismatch)));
        let deletion = simd.min_u16s(
            exact_affine_saturating_add_u16(simd, predecessor_match, opening),
            exact_affine_saturating_add_u16(simd, predecessor_deletion, extension),
        );
        let diagonal = exact_affine_saturating_add_u16(simd, predecessor_diagonal, substitution);
        let next_match = simd.min_u16s(diagonal, deletion);
        predecessor_matches[start..end].copy_from_slice(&pulp::cast::<_, [u16; 16]>(next_match));
        predecessor_deletions[start..end].copy_from_slice(&pulp::cast::<_, [u16; 16]>(deletion));
    }
    processed
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
fn exact_affine_last_u16(simd: pulp::x86::V3, values: pulp::u16x16) -> u16 {
    let extracted = simd.avx2._mm256_extract_epi16::<15>(pulp::cast(values));
    u16::try_from(extracted).expect("AVX2 extracts an unsigned 16-bit lane")
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline]
fn exact_affine_first_u16(simd: pulp::x86::V3, values: pulp::u16x16) -> u16 {
    let extracted = simd.avx2._mm256_extract_epi16::<0>(pulp::cast(values));
    u16::try_from(extracted).expect("AVX2 extracts an unsigned 16-bit lane")
}

#[inline]
#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_lines)]
#[allow(clippy::needless_range_loop)]
fn fill_exact_affine_single_predecessor_row(
    row: ExactAffineRowBand,
    predecessor_row: ExactAffineRowBand,
    query_indices: &[u8],
    graph_symbol: u8,
    mismatch: u16,
    opening: u16,
    extension: u16,
    insertion_scan: ExactAffineInsertionScan,
    match_costs: &mut [u16],
    deletion_costs: &mut [u16],
) -> usize {
    let mut insertion = u16::MAX;
    let fast_first = row.first.max(predecessor_row.first.saturating_add(1));
    let fast_last = row.last.min(predecessor_row.last);
    if predecessor_row == INVALID_EXACT_AFFINE_ROW_BAND || fast_first > fast_last {
        for query_prefix in row.first..=row.last {
            fill_exact_affine_single_predecessor_cell(
                query_prefix,
                row,
                predecessor_row,
                query_indices,
                graph_symbol,
                mismatch,
                opening,
                extension,
                match_costs,
                deletion_costs,
                &mut insertion,
            );
        }
        return 0;
    }

    for query_prefix in row.first..fast_first {
        fill_exact_affine_single_predecessor_cell(
            query_prefix,
            row,
            predecessor_row,
            query_indices,
            graph_symbol,
            mismatch,
            opening,
            extension,
            match_costs,
            deletion_costs,
            &mut insertion,
        );
    }

    let count = fast_last + 1 - fast_first;
    let cell_start = row.offset + fast_first - row.first;
    let predecessor_start = predecessor_row.offset + fast_first - predecessor_row.first;
    let previous_match = if fast_first > row.first {
        match_costs[cell_start - 1]
    } else {
        u16::MAX
    };
    {
        let (prior_matches, current_matches) = match_costs.split_at_mut(row.offset);
        let current_matches =
            &mut current_matches[fast_first - row.first..fast_first - row.first + count];
        let predecessor_matches = &prior_matches[predecessor_start..predecessor_start + count];
        let predecessor_diagonals =
            &prior_matches[predecessor_start - 1..predecessor_start - 1 + count];
        let (prior_deletions, current_deletions) = deletion_costs.split_at_mut(row.offset);
        let current_deletions =
            &mut current_deletions[fast_first - row.first..fast_first - row.first + count];
        let predecessor_deletions = &prior_deletions[predecessor_start..predecessor_start + count];
        let query_symbols = &query_indices[fast_first - 1..fast_last];
        for index in 0..count {
            let deletion = predecessor_matches[index]
                .saturating_add(opening)
                .min(predecessor_deletions[index].saturating_add(extension));
            current_deletions[index] = deletion;
            let substitution = if query_symbols[index] == graph_symbol {
                0
            } else {
                mismatch
            };
            current_matches[index] = predecessor_diagonals[index]
                .saturating_add(substitution)
                .min(deletion);
        }
    }
    let avx2_cells = {
        let current_start = fast_first - row.first;
        let current_end = current_start + count;
        let row_len = row.last + 1 - row.first;
        let current_matches = &mut match_costs[row.offset..row.offset + row_len];
        let mut opening_base = previous_match;
        let mut index = current_start;
        while index < current_end && (opening_base == u16::MAX || insertion == u16::MAX) {
            let base = current_matches[index];
            let next_insertion = opening_base
                .saturating_add(opening)
                .min(insertion.saturating_add(extension));
            current_matches[index] = base.min(next_insertion);
            opening_base = base;
            insertion = next_insertion;
            index += 1;
        }
        let (processed, next_opening_base, next_insertion) = insertion_scan.fill(
            &mut current_matches[index..current_end],
            opening,
            extension,
            opening_base,
            insertion,
        );
        index += processed;
        opening_base = next_opening_base;
        insertion = next_insertion;
        // Since opening >= extension, reopening from an insertion can never
        // beat extending it. The recurrence can therefore open from the prior
        // diagonal/deletion base instead of the finalized match value, removing
        // the match-store dependency from the insertion chain.
        for index in index..current_end {
            let base = current_matches[index];
            let next_insertion = opening_base
                .saturating_add(opening)
                .min(insertion.saturating_add(extension));
            current_matches[index] = base.min(next_insertion);
            opening_base = base;
            insertion = next_insertion;
        }
        processed
    };

    if fast_last < row.last {
        for query_prefix in fast_last + 1..=row.last {
            fill_exact_affine_single_predecessor_cell(
                query_prefix,
                row,
                predecessor_row,
                query_indices,
                graph_symbol,
                mismatch,
                opening,
                extension,
                match_costs,
                deletion_costs,
                &mut insertion,
            );
        }
    }
    avx2_cells
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn fill_exact_affine_two_predecessor_cell(
    query_prefix: usize,
    row: ExactAffineRowBand,
    first_predecessor_row: ExactAffineRowBand,
    second_predecessor_row: ExactAffineRowBand,
    query_indices: &[u8],
    graph_symbol: u8,
    mismatch: u16,
    opening: u16,
    extension: u16,
    match_costs: &mut [u16],
    deletion_costs: &mut [u16],
    insertion: &mut u16,
) {
    let cell = row.offset + query_prefix - row.first;
    let predecessor_match = exact_affine_cost(match_costs, first_predecessor_row, query_prefix)
        .min(exact_affine_cost(
            match_costs,
            second_predecessor_row,
            query_prefix,
        ));
    let predecessor_deletion =
        exact_affine_cost(deletion_costs, first_predecessor_row, query_prefix).min(
            exact_affine_cost(deletion_costs, second_predecessor_row, query_prefix),
        );
    let predecessor_diagonal = if query_prefix == 0 {
        u16::MAX
    } else {
        exact_affine_cost(match_costs, first_predecessor_row, query_prefix - 1).min(
            exact_affine_cost(match_costs, second_predecessor_row, query_prefix - 1),
        )
    };
    let deletion = predecessor_match
        .saturating_add(opening)
        .min(predecessor_deletion.saturating_add(extension));
    deletion_costs[cell] = deletion;
    let diagonal = if query_prefix == 0 {
        u16::MAX
    } else {
        let substitution = if query_indices[query_prefix - 1] == graph_symbol {
            0
        } else {
            mismatch
        };
        predecessor_diagonal.saturating_add(substitution)
    };
    *insertion = if query_prefix > row.first {
        match_costs[cell - 1]
            .saturating_add(opening)
            .min(insertion.saturating_add(extension))
    } else {
        u16::MAX
    };
    match_costs[cell] = diagonal.min(deletion).min(*insertion);
}

#[inline]
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
#[allow(clippy::needless_range_loop)]
fn fill_exact_affine_two_predecessor_row(
    row: ExactAffineRowBand,
    first_predecessor_row: ExactAffineRowBand,
    second_predecessor_row: ExactAffineRowBand,
    query_indices: &[u8],
    graph_symbol: u8,
    mismatch: u16,
    opening: u16,
    extension: u16,
    insertion_scan: ExactAffineInsertionScan,
    match_costs: &mut [u16],
    deletion_costs: &mut [u16],
) -> usize {
    let mut insertion = u16::MAX;
    let fast_first = row
        .first
        .max(first_predecessor_row.first.saturating_add(1))
        .max(second_predecessor_row.first.saturating_add(1));
    let fast_last = row
        .last
        .min(first_predecessor_row.last)
        .min(second_predecessor_row.last);
    if first_predecessor_row == INVALID_EXACT_AFFINE_ROW_BAND
        || second_predecessor_row == INVALID_EXACT_AFFINE_ROW_BAND
        || fast_first > fast_last
    {
        for query_prefix in row.first..=row.last {
            fill_exact_affine_two_predecessor_cell(
                query_prefix,
                row,
                first_predecessor_row,
                second_predecessor_row,
                query_indices,
                graph_symbol,
                mismatch,
                opening,
                extension,
                match_costs,
                deletion_costs,
                &mut insertion,
            );
        }
        return 0;
    }

    for query_prefix in row.first..fast_first {
        fill_exact_affine_two_predecessor_cell(
            query_prefix,
            row,
            first_predecessor_row,
            second_predecessor_row,
            query_indices,
            graph_symbol,
            mismatch,
            opening,
            extension,
            match_costs,
            deletion_costs,
            &mut insertion,
        );
    }

    let count = fast_last + 1 - fast_first;
    let cell_start = row.offset + fast_first - row.first;
    let first_predecessor_start =
        first_predecessor_row.offset + fast_first - first_predecessor_row.first;
    let second_predecessor_start =
        second_predecessor_row.offset + fast_first - second_predecessor_row.first;
    let previous_match = if fast_first > row.first {
        match_costs[cell_start - 1]
    } else {
        u16::MAX
    };
    {
        let (prior_matches, current_matches) = match_costs.split_at_mut(row.offset);
        let current_matches =
            &mut current_matches[fast_first - row.first..fast_first - row.first + count];
        let first_predecessor_matches =
            &prior_matches[first_predecessor_start..first_predecessor_start + count];
        let first_predecessor_diagonals =
            &prior_matches[first_predecessor_start - 1..first_predecessor_start - 1 + count];
        let second_predecessor_matches =
            &prior_matches[second_predecessor_start..second_predecessor_start + count];
        let second_predecessor_diagonals =
            &prior_matches[second_predecessor_start - 1..second_predecessor_start - 1 + count];
        let (prior_deletions, current_deletions) = deletion_costs.split_at_mut(row.offset);
        let current_deletions =
            &mut current_deletions[fast_first - row.first..fast_first - row.first + count];
        let first_predecessor_deletions =
            &prior_deletions[first_predecessor_start..first_predecessor_start + count];
        let second_predecessor_deletions =
            &prior_deletions[second_predecessor_start..second_predecessor_start + count];
        let query_symbols = &query_indices[fast_first - 1..fast_last];
        for index in 0..count {
            let predecessor_match =
                first_predecessor_matches[index].min(second_predecessor_matches[index]);
            let predecessor_deletion =
                first_predecessor_deletions[index].min(second_predecessor_deletions[index]);
            let deletion = predecessor_match
                .saturating_add(opening)
                .min(predecessor_deletion.saturating_add(extension));
            current_deletions[index] = deletion;
            let substitution = if query_symbols[index] == graph_symbol {
                0
            } else {
                mismatch
            };
            let predecessor_diagonal =
                first_predecessor_diagonals[index].min(second_predecessor_diagonals[index]);
            current_matches[index] = predecessor_diagonal
                .saturating_add(substitution)
                .min(deletion);
        }
    }
    let avx2_cells = {
        let current_start = fast_first - row.first;
        let current_end = current_start + count;
        let row_len = row.last + 1 - row.first;
        let current_matches = &mut match_costs[row.offset..row.offset + row_len];
        let mut opening_base = previous_match;
        let mut index = current_start;
        while index < current_end && (opening_base == u16::MAX || insertion == u16::MAX) {
            let base = current_matches[index];
            let next_insertion = opening_base
                .saturating_add(opening)
                .min(insertion.saturating_add(extension));
            current_matches[index] = base.min(next_insertion);
            opening_base = base;
            insertion = next_insertion;
            index += 1;
        }
        let (processed, next_opening_base, next_insertion) = insertion_scan.fill(
            &mut current_matches[index..current_end],
            opening,
            extension,
            opening_base,
            insertion,
        );
        index += processed;
        opening_base = next_opening_base;
        insertion = next_insertion;
        for index in index..current_end {
            let base = current_matches[index];
            let next_insertion = opening_base
                .saturating_add(opening)
                .min(insertion.saturating_add(extension));
            current_matches[index] = base.min(next_insertion);
            opening_base = base;
            insertion = next_insertion;
        }
        processed
    };

    if fast_last < row.last {
        for query_prefix in fast_last + 1..=row.last {
            fill_exact_affine_two_predecessor_cell(
                query_prefix,
                row,
                first_predecessor_row,
                second_predecessor_row,
                query_indices,
                graph_symbol,
                mismatch,
                opening,
                extension,
                match_costs,
                deletion_costs,
                &mut insertion,
            );
        }
    }
    avx2_cells
}

#[inline]
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_exact_affine_many_predecessor_row(
    row: ExactAffineRowBand,
    predecessors: &[PoaNodeId],
    row_bands: &[ExactAffineRowBand],
    query_indices: &[u8],
    graph_symbol: u8,
    mismatch: u16,
    opening: u16,
    extension: u16,
    insertion_scan: ExactAffineInsertionScan,
    match_costs: &mut [u16],
    deletion_costs: &mut [u16],
) -> usize {
    debug_assert!(predecessors.len() > 2);
    let row_len = row.last + 1 - row.first;
    let first_predecessor_diagonal = if row.first == 0 {
        u16::MAX
    } else {
        predecessors
            .iter()
            .map(|predecessor| {
                exact_affine_cost(
                    match_costs,
                    row_bands[predecessor.index() + 1],
                    row.first - 1,
                )
            })
            .min()
            .unwrap_or(u16::MAX)
    };

    {
        let (prior_matches, current_and_later_matches) = match_costs.split_at_mut(row.offset);
        let current_matches = &mut current_and_later_matches[..row_len];
        let (prior_deletions, current_and_later_deletions) =
            deletion_costs.split_at_mut(row.offset);
        let current_deletions = &mut current_and_later_deletions[..row_len];
        current_matches.fill(u16::MAX);
        current_deletions.fill(u16::MAX);

        for predecessor in predecessors {
            let predecessor_row = row_bands[predecessor.index() + 1];
            if predecessor_row == INVALID_EXACT_AFFINE_ROW_BAND {
                continue;
            }
            let first = row.first.max(predecessor_row.first);
            let last = row.last.min(predecessor_row.last);
            if first > last {
                continue;
            }
            let count = last + 1 - first;
            let current_start = first - row.first;
            let predecessor_start = predecessor_row.offset + first - predecessor_row.first;
            let current_matches = &mut current_matches[current_start..current_start + count];
            let current_deletions = &mut current_deletions[current_start..current_start + count];
            let predecessor_matches = &prior_matches[predecessor_start..predecessor_start + count];
            let predecessor_deletions =
                &prior_deletions[predecessor_start..predecessor_start + count];
            for (((current_match, current_deletion), &predecessor_match), &predecessor_deletion) in
                current_matches
                    .iter_mut()
                    .zip(current_deletions)
                    .zip(predecessor_matches)
                    .zip(predecessor_deletions)
            {
                *current_match = (*current_match).min(predecessor_match);
                *current_deletion = (*current_deletion).min(predecessor_deletion);
            }
        }
    }

    let current_matches = &mut match_costs[row.offset..row.offset + row_len];
    let current_deletions = &mut deletion_costs[row.offset..row.offset + row_len];
    // Moving backwards preserves every predecessor-match value until its
    // shifted diagonal consumer has read it, while replacing the row in place
    // with the independent diagonal/deletion base for the insertion scan.
    let base_avx2_cells = insertion_scan.fill_many_bases(
        current_matches,
        current_deletions,
        &query_indices[row.first..row.last],
        graph_symbol,
        mismatch,
        opening,
        extension,
    );
    for index in (1..row_len - base_avx2_cells).rev() {
        let query_prefix = row.first + index;
        let predecessor_match = current_matches[index];
        let predecessor_deletion = current_deletions[index];
        let deletion = predecessor_match
            .saturating_add(opening)
            .min(predecessor_deletion.saturating_add(extension));
        let substitution = if query_indices[query_prefix - 1] == graph_symbol {
            0
        } else {
            mismatch
        };
        current_matches[index] = current_matches[index - 1]
            .saturating_add(substitution)
            .min(deletion);
        current_deletions[index] = deletion;
    }
    let predecessor_match = current_matches[0];
    let predecessor_deletion = current_deletions[0];
    let deletion = predecessor_match
        .saturating_add(opening)
        .min(predecessor_deletion.saturating_add(extension));
    let diagonal = if row.first == 0 {
        u16::MAX
    } else {
        let substitution = if query_indices[row.first - 1] == graph_symbol {
            0
        } else {
            mismatch
        };
        first_predecessor_diagonal.saturating_add(substitution)
    };
    current_matches[0] = diagonal.min(deletion);
    current_deletions[0] = deletion;

    let mut opening_base = current_matches[0];
    let mut insertion = u16::MAX;
    let mut index = 1;
    while index < row_len && (opening_base == u16::MAX || insertion == u16::MAX) {
        let base = current_matches[index];
        insertion = opening_base
            .saturating_add(opening)
            .min(insertion.saturating_add(extension));
        current_matches[index] = base.min(insertion);
        opening_base = base;
        index += 1;
    }
    let (avx2_cells, next_opening_base, next_insertion) = insertion_scan.fill(
        &mut current_matches[index..],
        opening,
        extension,
        opening_base,
        insertion,
    );
    index += avx2_cells;
    opening_base = next_opening_base;
    insertion = next_insertion;
    for base in &mut current_matches[index..] {
        let next_insertion = opening_base
            .saturating_add(opening)
            .min(insertion.saturating_add(extension));
        let next_base = *base;
        *base = next_base.min(next_insertion);
        opening_base = next_base;
        insertion = next_insertion;
    }
    avx2_cells
}

fn prepare_exact_affine_buffer(
    mut buffer: Vec<u16>,
    elements: usize,
    name: &'static str,
) -> Result<Vec<u16>, PoaError> {
    // Preserve the initialized high-water length: diagonal and affine
    // certificates alternate in progressive workloads, and truncating here
    // would force Vec::resize to rewrite already allocated storage on the next
    // larger certificate. Callers initialize and access only `..elements`.
    if buffer.len() < elements {
        buffer
            .try_reserve(elements - buffer.len())
            .map_err(|_| PoaError::AllocationFailed {
                buffer: name,
                elements,
            })?;
        buffer.resize(elements, u16::MAX);
    }
    Ok(buffer)
}

#[allow(clippy::struct_field_names)]
struct ExactAffineTables {
    match_costs: Vec<u16>,
    deletion_costs: Vec<u16>,
}

fn prepare_exact_affine_tables(
    workspace: Option<&mut PoastaWorkspace>,
    elements: usize,
) -> Result<ExactAffineTables, PoaError> {
    let Some(workspace) = workspace else {
        return Ok(ExactAffineTables {
            match_costs: filled(elements, u16::MAX, "POASTA exact compact match costs")?,
            deletion_costs: filled(elements, u16::MAX, "POASTA exact compact deletion costs")?,
        });
    };
    let match_costs = prepare_exact_affine_buffer(
        core::mem::take(&mut workspace.exact_match_costs),
        elements,
        "POASTA exact compact match costs",
    )?;
    let deletion_costs = prepare_exact_affine_buffer(
        core::mem::take(&mut workspace.exact_deletion_costs),
        elements,
        "POASTA exact compact deletion costs",
    )?;
    Ok(ExactAffineTables {
        match_costs,
        deletion_costs,
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
#[cfg(test)]
fn exact_affine_certificate<A: Alphabet>(
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    topological: &[PoaNodeId],
    index: &PoastaGraphIndex,
    costs: PoastaAffineCosts,
    feasible_upper_bound: Option<Cost>,
    cached_prefix_bounds: Option<&[RelaxedLinearRowBounds]>,
    workspace: Option<&mut PoastaWorkspace>,
) -> Result<Option<ExactAffineCertificateAttempt>, PoaError> {
    exact_affine_certificate_profiled::<A, false>(
        graph,
        query,
        topological,
        index,
        costs,
        feasible_upper_bound,
        cached_prefix_bounds,
        workspace,
        None,
        None,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn exact_affine_certificate_profiled<A: Alphabet, const DETAILED: bool>(
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    topological: &[PoaNodeId],
    index: &PoastaGraphIndex,
    costs: PoastaAffineCosts,
    feasible_upper_bound: Option<Cost>,
    cached_prefix_bounds: Option<&[RelaxedLinearRowBounds]>,
    mut workspace: Option<&mut PoastaWorkspace>,
    mut search_timings: Option<&mut PhaseTimingsBuilder>,
    mut search_counters: Option<&mut PoastaSearchWorkCounters>,
) -> Result<Option<ExactAffineCertificateAttempt>, PoaError> {
    let setup_started = search_phase_started::<DETAILED>(&search_timings);
    if DETAILED {
        if let Some(counters) = search_counters.as_deref_mut() {
            counters.record_exact_affine_attempt()?;
        }
    }
    let extension_cost = nonnegative_score(costs.gap_extension());
    let opening_cost = checked_cost_add(
        nonnegative_score(costs.gap_opening()),
        extension_cost,
        "the exact compact affine opening cost",
    )?;
    let mismatch_cost = nonnegative_score(costs.mismatch());
    let maximum_steps = index
        .max_remaining_nodes(index.source())
        .and_then(|remaining| remaining.checked_add(query.len()))
        .ok_or(PoaError::DimensionOverflow {
            graph_nodes: graph.node_count(),
            query_len: query.len(),
        })?;
    let maximum_unit = mismatch_cost.max(opening_cost);
    let maximum_cost = checked_cost_mul(
        maximum_steps,
        maximum_unit,
        "the exact compact affine range",
    )?;
    if maximum_cost >= Cost::from(u16::MAX) {
        record_search_phase::<DETAILED>(
            search_timings.as_deref_mut(),
            "exact_affine_setup",
            setup_started,
        );
        return Ok(None);
    }
    let extension = u16::try_from(extension_cost).map_err(|_| PoaError::PoastaCostOverflow {
        context: "the exact compact affine extension cost",
    })?;
    let opening = u16::try_from(opening_cost).map_err(|_| PoaError::PoastaCostOverflow {
        context: "the exact compact affine opening cost",
    })?;
    let mismatch = u16::try_from(mismatch_cost).map_err(|_| PoaError::PoastaCostOverflow {
        context: "the exact compact affine mismatch cost",
    })?;
    let upper_bound = match feasible_upper_bound {
        Some(upper_bound) => upper_bound,
        None => exact_affine_upper_bound(
            graph,
            query,
            index,
            mismatch_cost,
            opening_cost,
            extension_cost,
        )?,
    };
    record_search_phase::<DETAILED>(
        search_timings.as_deref_mut(),
        "exact_affine_setup",
        setup_started,
    );
    let computed_prefix_bounds = if cached_prefix_bounds.is_none() {
        Some(measure_search_phase::<DETAILED, _>(
            &mut search_timings,
            "exact_prefix_bounds",
            || exact_affine_prefix_bounds(graph, topological, index, query.len()),
        )?)
    } else {
        None
    };
    let prefix_bounds = cached_prefix_bounds
        .or(computed_prefix_bounds.as_deref())
        .expect("exact prefix bounds are cached or computed");
    let setup_started = search_phase_started::<DETAILED>(&search_timings);
    let band_radius = if extension_cost == 0 {
        (upper_bound < opening_cost).then_some(0)
    } else if upper_bound < opening_cost {
        Some(0)
    } else {
        // Any non-zero net length discrepancy requires at least one affine
        // gap run: opening + (symbols - 1) * extension. Charging that mandatory
        // opening tightens the exact band while retaining every path no more
        // expensive than the feasible incumbent.
        let radius = 1_u64
            .checked_add((upper_bound - opening_cost) / extension_cost)
            .ok_or(PoaError::ArithmeticOverflow {
                context: "the exact-affine safe-band radius",
            })?;
        Some(
            usize::try_from(radius).map_err(|_| PoaError::PoastaCostOverflow {
                context: "the exact-affine safe-band radius",
            })?,
        )
    };
    let cached_rows = workspace.as_deref_mut().map_or_else(Vec::new, |workspace| {
        core::mem::take(&mut workspace.exact_row_bands)
    });
    record_search_phase::<DETAILED>(
        search_timings.as_deref_mut(),
        "exact_affine_setup",
        setup_started,
    );
    let ExactAffineBandLayout {
        rows: row_bands,
        cells: elements,
    } = measure_search_phase::<DETAILED, _>(
        &mut search_timings,
        "exact_affine_band_layout",
        || {
            exact_affine_band_layout(
                cached_rows,
                topological,
                index,
                prefix_bounds,
                query.len(),
                band_radius,
            )
        },
    )?;
    let ExactAffineTables {
        mut match_costs,
        mut deletion_costs,
    } = measure_search_phase::<DETAILED, _>(
        &mut search_timings,
        "exact_affine_table_prepare",
        || prepare_exact_affine_tables(workspace.as_deref_mut(), elements),
    )?;
    let query_indices = query.as_indices();
    let insertion_scan = ExactAffineInsertionScan::new(extension);
    let source = index.source();
    let source_row = row_bands[source];
    let fill_started = search_phase_started::<DETAILED>(&search_timings);
    if DETAILED {
        if let Some(counters) = search_counters.as_deref_mut() {
            counters.record_exact_affine_source(exact_affine_row_cells(source_row))?;
        }
    }
    deletion_costs[source_row.offset..=source_row.offset + query.len()].fill(u16::MAX);
    match_costs[source_row.offset] = 0;
    let mut source_insertion = u16::MAX;
    for query_prefix in 1..=query.len() {
        let cell = source_row.offset + query_prefix;
        source_insertion = match_costs[cell - 1]
            .saturating_add(opening)
            .min(source_insertion.saturating_add(extension));
        match_costs[cell] = source_insertion;
    }

    for &node_id in topological {
        let cursor = node_id.index() + 1;
        let row = row_bands[cursor];
        if row == INVALID_EXACT_AFFINE_ROW_BAND {
            continue;
        }
        let node = graph
            .node(node_id)
            .ok_or(PoaError::UnknownNode { node: node_id })?;
        let graph_symbol =
            u8::try_from(node.symbol().alphabet_index()).map_err(|_| PoaError::InvalidGraph {
                message: "a graph symbol exceeds the compact alphabet range".to_owned(),
            })?;
        let predecessors = node.predecessors();
        let single_predecessor = match predecessors {
            [] => Some(index.source()),
            [predecessor] => Some(predecessor.index() + 1),
            _ => None,
        };
        if let Some(predecessor) = single_predecessor {
            let avx2_cells = fill_exact_affine_single_predecessor_row(
                row,
                row_bands[predecessor],
                query_indices,
                graph_symbol,
                mismatch,
                opening,
                extension,
                insertion_scan,
                &mut match_costs,
                &mut deletion_costs,
            );
            if DETAILED {
                if let Some(counters) = search_counters.as_deref_mut() {
                    counters.record_exact_affine_single(
                        exact_affine_row_cells(row),
                        exact_affine_single_fast_cells(row, row_bands[predecessor]),
                        avx2_cells,
                    )?;
                }
            }
        } else if let [first_predecessor, second_predecessor] = predecessors {
            let avx2_cells = fill_exact_affine_two_predecessor_row(
                row,
                row_bands[first_predecessor.index() + 1],
                row_bands[second_predecessor.index() + 1],
                query_indices,
                graph_symbol,
                mismatch,
                opening,
                extension,
                insertion_scan,
                &mut match_costs,
                &mut deletion_costs,
            );
            if DETAILED {
                if let Some(counters) = search_counters.as_deref_mut() {
                    counters.record_exact_affine_two(
                        exact_affine_row_cells(row),
                        exact_affine_two_fast_cells(
                            row,
                            row_bands[first_predecessor.index() + 1],
                            row_bands[second_predecessor.index() + 1],
                        ),
                        avx2_cells,
                    )?;
                }
            }
        } else {
            let avx2_cells = fill_exact_affine_many_predecessor_row(
                row,
                predecessors,
                &row_bands,
                query_indices,
                graph_symbol,
                mismatch,
                opening,
                extension,
                insertion_scan,
                &mut match_costs,
                &mut deletion_costs,
            );
            if DETAILED {
                if let Some(counters) = search_counters.as_deref_mut() {
                    counters.record_exact_affine_many(
                        exact_affine_row_cells(row),
                        exact_affine_many_predecessor_cell_visits(row, predecessors, &row_bands)?,
                        avx2_cells,
                    )?;
                }
            }
        }
    }
    record_search_phase::<DETAILED>(
        search_timings.as_deref_mut(),
        "exact_affine_fill",
        fill_started,
    );

    let goal_started = search_phase_started::<DETAILED>(&search_timings);
    let mut best_goal = None;
    if graph.node_count() == 0 {
        best_goal = Some((match_costs[source_row.offset + query.len()], source));
    } else {
        for &node_id in topological {
            let cursor = node_id.index() + 1;
            if index.successors(cursor) != [index.sink()] {
                continue;
            }
            let Some(cell) = exact_affine_cell(row_bands[cursor], query.len()) else {
                continue;
            };
            let cost = match_costs[cell];
            if best_goal.is_none_or(|(best_cost, best_cursor)| {
                cost < best_cost || (cost == best_cost && cursor < best_cursor)
            }) {
                best_goal = Some((cost, cursor));
            }
        }
    }
    let Some((optimum_cost, goal_cursor)) = best_goal.filter(|(cost, _)| *cost != u16::MAX) else {
        return Err(PoaError::InvalidGraph {
            message: "the exact-affine table has no reachable terminal".to_owned(),
        });
    };
    let table_cell_capacity = match_costs
        .capacity()
        .checked_add(deletion_costs.capacity())
        .ok_or(PoaError::ArithmeticOverflow {
            context: "POASTA exact compact affine cell capacity",
        })?;
    let table_payload_bytes = checked_bytes(
        table_cell_capacity,
        size_of::<u16>(),
        "POASTA exact compact affine payload",
    )?;
    let prefix_payload_bytes = checked_bytes(
        prefix_bounds.len(),
        size_of::<RelaxedLinearRowBounds>(),
        "POASTA exact-affine prefix-bound payload",
    )?;
    let band_payload_bytes = checked_bytes(
        row_bands.capacity(),
        size_of::<ExactAffineRowBand>(),
        "POASTA exact-affine compact row-band payload",
    )?;
    let payload_bytes = if workspace.is_some() {
        prefix_payload_bytes
    } else {
        checked_add_usize(
            checked_add_usize(
                table_payload_bytes,
                band_payload_bytes,
                "POASTA exact-affine compact table payload",
            )?,
            prefix_payload_bytes,
            "POASTA exact-affine table and prefix-bound payload",
        )?
    };
    record_search_phase::<DETAILED>(
        search_timings.as_deref_mut(),
        "exact_affine_goal",
        goal_started,
    );

    let traceback_started = search_phase_started::<DETAILED>(&search_timings);
    let maximum_trace_states = maximum_steps
        .checked_mul(2)
        .and_then(|states| states.checked_add(1))
        .ok_or(PoaError::DimensionOverflow {
            graph_nodes: graph.node_count(),
            query_len: query.len(),
        })?;
    let mut reversed_steps = Vec::new();
    try_reserve_exact(
        &mut reversed_steps,
        maximum_trace_states,
        "POASTA exact compact certificate",
    )?;
    let start = StateKey::new(source, 0, AffineState::Match);
    let mut current = StateKey::new(goal_cursor, query.len(), AffineState::Match);
    while current != start {
        if reversed_steps.len() >= maximum_trace_states {
            return Err(PoaError::AffineTracebackInvariant {
                node: index.node_for_cursor(current.cursor),
                column: current.query_prefix(),
                state: current.state().name(),
            });
        }
        let query_prefix = current.query_prefix();
        let current_row = row_bands[current.cursor];
        let cell = exact_affine_cell(current_row, query_prefix).ok_or_else(|| {
            PoaError::AffineTracebackInvariant {
                node: index.node_for_cursor(current.cursor),
                column: query_prefix,
                state: current.state().name(),
            }
        })?;
        let insertion_costs = (current.state() == AffineState::Insertion).then(|| {
            exact_affine_insertion_costs_at(
                &match_costs,
                current_row,
                query_prefix,
                opening,
                extension,
            )
        });
        let current_cost = match current.state() {
            AffineState::Match => match_costs[cell],
            AffineState::Insertion => {
                insertion_costs
                    .expect("an insertion state has recomputed insertion costs")
                    .0
            }
            AffineState::Deletion => deletion_costs[cell],
        };
        let mut selected: Option<(StateKey, Cost, u8)> = None;

        match current.state() {
            AffineState::Match => {
                if current.cursor != source && query_prefix > 0 {
                    let node_id = index.node_for_cursor(current.cursor).ok_or_else(|| {
                        PoaError::InvalidGraph {
                            message: "an exact-affine cursor has no graph node".to_owned(),
                        }
                    })?;
                    let node = graph
                        .node(node_id)
                        .ok_or(PoaError::UnknownNode { node: node_id })?;
                    let graph_symbol =
                        u8::try_from(node.symbol().alphabet_index()).map_err(|_| {
                            PoaError::InvalidGraph {
                                message: "a graph symbol exceeds the compact alphabet range"
                                    .to_owned(),
                            }
                        })?;
                    let substitution = if query_indices[query_prefix - 1] == graph_symbol {
                        0
                    } else {
                        mismatch
                    };
                    for predecessor in augmented_predecessor_cursors(node.predecessors(), source) {
                        let predecessor_cost = exact_affine_cost(
                            &match_costs,
                            row_bands[predecessor],
                            query_prefix - 1,
                        );
                        if predecessor_cost.saturating_add(substitution) == current_cost {
                            consider_exact_affine_predecessor(
                                &mut selected,
                                StateKey::new(predecessor, query_prefix - 1, AffineState::Match),
                                Cost::from(substitution),
                                0,
                            );
                        }
                    }
                }
                if deletion_costs[cell] == current_cost {
                    consider_exact_affine_predecessor(
                        &mut selected,
                        StateKey::new(current.cursor, query_prefix, AffineState::Deletion),
                        0,
                        1,
                    );
                }
                if selected.is_none() {
                    let (insertion, _) = exact_affine_insertion_costs_at(
                        &match_costs,
                        current_row,
                        query_prefix,
                        opening,
                        extension,
                    );
                    if insertion == current_cost {
                        consider_exact_affine_predecessor(
                            &mut selected,
                            StateKey::new(current.cursor, query_prefix, AffineState::Insertion),
                            0,
                            2,
                        );
                    }
                }
            }
            AffineState::Insertion => {
                let (_, predecessor_insertion) =
                    insertion_costs.expect("an insertion state has recomputed insertion costs");
                if query_prefix > 0
                    && predecessor_insertion.saturating_add(extension) == current_cost
                {
                    consider_exact_affine_predecessor(
                        &mut selected,
                        StateKey::new(current.cursor, query_prefix - 1, AffineState::Insertion),
                        Cost::from(extension),
                        0,
                    );
                }
                if query_prefix > 0
                    && exact_affine_cost(&match_costs, current_row, query_prefix - 1)
                        .saturating_add(opening)
                        == current_cost
                {
                    consider_exact_affine_predecessor(
                        &mut selected,
                        StateKey::new(current.cursor, query_prefix - 1, AffineState::Match),
                        Cost::from(opening),
                        1,
                    );
                }
            }
            AffineState::Deletion => {
                let node_id = index.node_for_cursor(current.cursor).ok_or_else(|| {
                    PoaError::InvalidGraph {
                        message: "an exact-affine cursor has no graph node".to_owned(),
                    }
                })?;
                let node = graph
                    .node(node_id)
                    .ok_or(PoaError::UnknownNode { node: node_id })?;
                for predecessor in augmented_predecessor_cursors(node.predecessors(), source) {
                    if exact_affine_cost(&deletion_costs, row_bands[predecessor], query_prefix)
                        .saturating_add(extension)
                        == current_cost
                    {
                        consider_exact_affine_predecessor(
                            &mut selected,
                            StateKey::new(predecessor, query_prefix, AffineState::Deletion),
                            Cost::from(extension),
                            0,
                        );
                    }
                    if exact_affine_cost(&match_costs, row_bands[predecessor], query_prefix)
                        .saturating_add(opening)
                        == current_cost
                    {
                        consider_exact_affine_predecessor(
                            &mut selected,
                            StateKey::new(predecessor, query_prefix, AffineState::Match),
                            Cost::from(opening),
                            1,
                        );
                    }
                }
            }
        }
        let Some((predecessor, transition_cost, tie_transition)) = selected else {
            return Err(PoaError::AffineTracebackInvariant {
                node: index.node_for_cursor(current.cursor),
                column: current.query_prefix(),
                state: current.state().name(),
            });
        };
        reversed_steps.push(CertifiedStep {
            key: current,
            transition_cost,
            tie_transition,
        });
        current = predecessor;
    }
    reversed_steps.reverse();
    record_search_phase::<DETAILED>(search_timings, "exact_affine_traceback", traceback_started);
    if let Some(workspace) = workspace {
        workspace.exact_match_costs = match_costs;
        workspace.exact_deletion_costs = deletion_costs;
        workspace.exact_row_bands = row_bands;
    }
    Ok(Some(ExactAffineCertificateAttempt {
        certificate: Some(ExactAffineCertificate {
            steps: reversed_steps,
            optimum_cost: Cost::from(optimum_cost),
        }),
        cells: elements,
        payload_bytes,
    }))
}

struct Search<'a, A: Alphabet> {
    graph: &'a PoaGraph<A>,
    query: &'a Sequence<A>,
    costs: PoastaAffineCosts,
    options: PoastaSearchOptions,
    index: &'a PoastaGraphIndex,
    topological: &'a [PoaNodeId],
    columns: usize,
    records: StateMap,
    frontier: BucketFrontier,
    dfs_stack: Vec<DfsFrame>,
    reached_exits: Vec<ReachedExit>,
    relaxed_linear_heuristic: Option<RelaxedLinearHeuristic>,
    certified_traceback: Option<Vec<StateKey>>,
    workspace: Option<&'a mut PoastaWorkspace>,
    relaxed_linear_attempted: bool,
    eager_relaxed_linear: bool,
    certificate_heuristic_payload_bytes: usize,
    exact_diagonal_attempted: bool,
    exact_diagonal_attempt: Option<DiagonalCertificateAttempt>,
    exact_prefix_bounds: Option<Vec<RelaxedLinearRowBounds>>,
    exact_affine_upper_bound: Option<Cost>,
    statistics: PoastaSearchStatistics,
    best_goal: Option<(Cost, StateKey)>,
}

impl<'a, A: Alphabet> Search<'a, A> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        graph: &'a PoaGraph<A>,
        query: &'a Sequence<A>,
        costs: PoastaAffineCosts,
        options: PoastaSearchOptions,
        index: &'a PoastaGraphIndex,
        topological: &'a [PoaNodeId],
        eager_relaxed_linear: bool,
        total_possible_states: usize,
        workspace: Option<&'a mut PoastaWorkspace>,
    ) -> Result<Self, PoaError> {
        let mut records = StateMap::default();
        let eager_certificate =
            eager_relaxed_linear && matches!(options.backend, PoastaBackend::Auto);
        let initial_record_capacity = if eager_certificate {
            0
        } else {
            total_possible_states.min(1024)
        };
        records
            .try_reserve(initial_record_capacity)
            .map_err(|_| PoaError::AllocationFailed {
                buffer: "POASTA sparse state records",
                elements: initial_record_capacity,
            })?;
        let mut reached_exits = Vec::new();
        try_reserve_exact(
            &mut reached_exits,
            index.bubbles().len(),
            "POASTA reached-exit maps",
        )?;
        reached_exits.resize_with(index.bubbles().len(), ReachedExit::default);
        let initial_dfs_capacity = if eager_certificate {
            0
        } else {
            graph.node_count().min(256).saturating_add(1)
        };
        let mut dfs_stack = Vec::new();
        dfs_stack
            .try_reserve(initial_dfs_capacity)
            .map_err(|_| PoaError::AllocationFailed {
                buffer: "POASTA exact-match DFS stack",
                elements: initial_dfs_capacity,
            })?;
        Ok(Self {
            graph,
            query,
            costs,
            options,
            index,
            topological,
            columns: query.len() + 1,
            records,
            frontier: BucketFrontier::default(),
            dfs_stack,
            reached_exits,
            relaxed_linear_heuristic: None,
            certified_traceback: None,
            workspace,
            relaxed_linear_attempted: !matches!(options.backend, PoastaBackend::Auto),
            eager_relaxed_linear,
            certificate_heuristic_payload_bytes: 0,
            exact_diagonal_attempted: false,
            exact_diagonal_attempt: None,
            exact_prefix_bounds: None,
            exact_affine_upper_bound: None,
            statistics: PoastaSearchStatistics {
                total_possible_states,
                discovered_states: 0,
                expanded_states: 0,
                astar_expanded_states: 0,
                dfs_expanded_states: 0,
                queue_entries: 0,
                reopened_states: 0,
                stale_entries: 0,
                greedy_exact_matches: 0,
                affine_dominated_states: 0,
                linear_heuristic_cells: 0,
                relaxed_linear_source_bound: None,
                relaxed_linear_certificate: false,
                exact_diagonal_certificate_cells: 0,
                exact_diagonal_cost: None,
                affine_opening_certificate: false,
                exact_affine_certificate_cells: 0,
                exact_affine_certificate: false,
                bubble_tests: 0,
                bubble_pruned_states: 0,
                reached_exit_states: 0,
                superbubbles: index.bubbles().len(),
                bubble_membership_entries: index.membership_entries(),
                peak_queue_entries: 0,
                peak_dfs_entries: 0,
            },
            best_goal: None,
        })
    }

    fn run<const DETAILED: bool>(
        &mut self,
        mut search_timings: Option<&mut PhaseTimingsBuilder>,
        mut search_counters: Option<&mut PoastaSearchWorkCounters>,
    ) -> Result<(StateKey, Cost), PoaError> {
        let start = StateKey::new(self.index.source(), 0, AffineState::Match);
        if self.eager_relaxed_linear && !self.relaxed_linear_attempted {
            self.activate_relaxed_linear::<DETAILED>(
                search_timings.as_deref_mut(),
                search_counters.as_deref_mut(),
            )?;
        } else {
            self.seed_start_record(start);
            self.push_frontier(start, 0)?;
        }

        let mut frontier_started = search_phase_started::<DETAILED>(&search_timings);
        while let Some(Reverse(entry)) = self.frontier.pop() {
            if self
                .best_goal
                .is_some_and(|(cost, _)| entry.estimated_total > cost)
            {
                break;
            }
            if self.should_activate_relaxed_linear() {
                record_search_phase::<DETAILED>(
                    search_timings.as_deref_mut(),
                    "frontier_search",
                    frontier_started.take(),
                );
                self.activate_relaxed_linear::<DETAILED>(
                    search_timings.as_deref_mut(),
                    search_counters.as_deref_mut(),
                )?;
                frontier_started = search_phase_started::<DETAILED>(&search_timings);
                continue;
            }
            if self.options.exact_match_dfs && entry.key.state() == AffineState::Match {
                self.expand_with_dfs(entry)?;
            } else {
                self.expand_from_astar(entry)?;
            }
        }
        record_search_phase::<DETAILED>(search_timings, "frontier_search", frontier_started.take());

        let Some((cost, goal)) = self.best_goal else {
            return self.traceback_error(StateKey::new(
                self.index.source(),
                self.query.len(),
                AffineState::Match,
            ));
        };
        Ok((goal, cost))
    }

    fn recycle_relaxed_linear_costs(&mut self) {
        let Some(workspace) = self.workspace.as_deref_mut() else {
            return;
        };
        if let Some(heuristic) = self.relaxed_linear_heuristic.take() {
            workspace.costs = heuristic.costs;
        }
    }

    fn seed_start_record(&mut self, start: StateKey) {
        let state_index = self.state_index(start);
        self.records.entry(state_index).or_insert(StateRecord {
            cost: 0,
            predecessor: StateKey::NONE,
            expanded_cost: UNEXPANDED_COST,
            tie_transition: 0,
        });
    }

    fn should_activate_relaxed_linear(&self) -> bool {
        if self.relaxed_linear_attempted {
            return false;
        }
        let path_budget = self.query.len().saturating_add(1);
        let discovery_budget = path_budget.saturating_mul(RELAXED_LINEAR_DISCOVERY_FACTOR);
        self.statistics.expanded_states >= path_budget && self.records.len() >= discovery_budget
    }

    fn prefer_exact_affine_certificate(&self) -> bool {
        self.eager_relaxed_linear
            && self.workspace.as_deref().is_some_and(|workspace| {
                workspace.relaxed_certificate_misses >= RELAXED_LINEAR_MISS_THRESHOLD
                    && workspace.exact_affine_bypasses_since_probe < RELAXED_LINEAR_REPROBE_INTERVAL
            })
    }

    fn should_tighten_exact_affine_bypass(&self) -> bool {
        let generic_rows = self
            .graph
            .nodes()
            .iter()
            .filter(|node| node.predecessors().len() > 2)
            .count();
        generic_rows > self.graph.node_count() / 3
    }

    fn prepare_exact_diagonal_certificate<const DETAILED: bool>(
        &mut self,
        mut search_timings: Option<&mut PhaseTimingsBuilder>,
    ) -> Result<(), PoaError> {
        if self.exact_diagonal_attempted {
            return Ok(());
        }
        self.exact_diagonal_attempted = true;
        if self.exact_prefix_bounds.is_none() {
            self.exact_prefix_bounds = Some(measure_search_phase::<DETAILED, _>(
                &mut search_timings,
                "exact_prefix_bounds",
                || {
                    exact_affine_prefix_bounds(
                        self.graph,
                        self.topological,
                        self.index,
                        self.query.len(),
                    )
                },
            )?);
        }
        self.exact_diagonal_attempt = measure_search_phase::<DETAILED, _>(
            &mut search_timings,
            "exact_diagonal_certificate",
            || {
                exact_diagonal_certificate(
                    self.graph,
                    self.query,
                    self.topological,
                    self.index,
                    nonnegative_score(self.costs.mismatch()),
                    self.exact_prefix_bounds.as_deref(),
                    self.workspace.as_deref_mut(),
                )
            },
        )?;
        if let Some(attempt) = self.exact_diagonal_attempt.as_ref() {
            self.statistics.exact_diagonal_certificate_cells = attempt.cells;
            self.statistics.exact_diagonal_cost = Some(attempt.certificate.optimum_cost);
            self.certificate_heuristic_payload_bytes = self
                .certificate_heuristic_payload_bytes
                .max(attempt.payload_bytes);
        } else {
            let prefix_payload_bytes = checked_bytes(
                self.exact_prefix_bounds.as_ref().map_or(0, Vec::capacity),
                size_of::<RelaxedLinearRowBounds>(),
                "POASTA cached exact prefix-bound payload",
            )?;
            self.certificate_heuristic_payload_bytes = self
                .certificate_heuristic_payload_bytes
                .max(prefix_payload_bytes);
        }
        Ok(())
    }

    fn graph_lengths_equal_query(&self) -> bool {
        self.index
            .min_remaining_nodes(self.index.source())
            .zip(self.index.max_remaining_nodes(self.index.source()))
            .is_some_and(|(minimum, maximum)| {
                minimum == self.query.len() && maximum == self.query.len()
            })
    }

    fn try_early_diagonal_certificate<const DETAILED: bool>(
        &mut self,
        mut search_timings: Option<&mut PhaseTimingsBuilder>,
    ) -> Result<bool, PoaError> {
        self.prepare_exact_diagonal_certificate::<DETAILED>(search_timings.as_deref_mut())?;
        if self
            .exact_diagonal_attempt
            .as_ref()
            .is_none_or(|attempt| attempt.certificate.optimum_cost != 0)
        {
            return Ok(false);
        }
        let attempt = self
            .exact_diagonal_attempt
            .take()
            .expect("a checked early diagonal certificate exists");
        measure_search_phase::<DETAILED, _>(&mut search_timings, "certificate_install", || {
            self.install_diagonal_certificate(attempt.certificate, true)
        })?;
        Ok(true)
    }

    fn record_relaxed_certificate_outcome(&mut self, certified: bool) {
        let Some(workspace) = self.workspace.as_deref_mut() else {
            return;
        };
        workspace.exact_affine_bypasses_since_probe = 0;
        if certified {
            workspace.relaxed_certificate_misses = 0;
        } else {
            workspace.relaxed_certificate_misses =
                workspace.relaxed_certificate_misses.saturating_add(1);
        }
    }

    fn record_exact_affine_bypass(&mut self) {
        if let Some(workspace) = self.workspace.as_deref_mut() {
            workspace.exact_affine_bypasses_since_probe = workspace
                .exact_affine_bypasses_since_probe
                .saturating_add(1);
        }
    }

    fn try_exact_affine_certificate<const DETAILED: bool>(
        &mut self,
        mut search_timings: Option<&mut PhaseTimingsBuilder>,
        search_counters: Option<&mut PoastaSearchWorkCounters>,
    ) -> Result<bool, PoaError> {
        let feasible_upper_bound = self.exact_affine_upper_bound.take();
        let attempt = exact_affine_certificate_profiled::<A, DETAILED>(
            self.graph,
            self.query,
            self.topological,
            self.index,
            self.costs,
            feasible_upper_bound,
            self.exact_prefix_bounds.as_deref(),
            self.workspace.as_deref_mut(),
            search_timings.as_deref_mut(),
            search_counters,
        )?;
        let Some(attempt) = attempt else {
            return Ok(false);
        };
        self.statistics.exact_affine_certificate_cells = attempt.cells;
        self.certificate_heuristic_payload_bytes = self
            .certificate_heuristic_payload_bytes
            .max(attempt.payload_bytes);
        let Some(certificate) = attempt.certificate else {
            return Ok(false);
        };
        measure_search_phase::<DETAILED, _>(&mut search_timings, "certificate_install", || {
            self.install_exact_affine_certificate(&certificate)
        })?;
        Ok(true)
    }

    fn try_affine_opening_certificate<const DETAILED: bool>(
        &mut self,
        mut search_timings: Option<&mut PhaseTimingsBuilder>,
    ) -> Result<bool, PoaError> {
        let start = StateKey::new(self.index.source(), 0, AffineState::Match);
        let Some(relaxed_bound) = self
            .relaxed_linear_heuristic
            .as_ref()
            .map(|heuristic| heuristic.lower_bound(start))
        else {
            return Ok(false);
        };
        self.prepare_exact_diagonal_certificate::<DETAILED>(search_timings.as_deref_mut())?;
        let Some(attempt) = self.exact_diagonal_attempt.take() else {
            return Ok(false);
        };
        let minimum_gap_runs = if self.graph_lengths_equal_query() {
            2
        } else {
            1
        };
        let mandatory_openings = checked_cost_mul(
            minimum_gap_runs,
            nonnegative_score(self.costs.gap_opening()),
            "the mandatory affine-opening certificate margin",
        )?;
        let gapped_lower_bound = checked_cost_add(
            relaxed_bound,
            mandatory_openings,
            "the affine-opening certificate lower bound",
        )?;
        if attempt.certificate.optimum_cost >= gapped_lower_bound {
            self.exact_affine_upper_bound = Some(attempt.certificate.optimum_cost);
            return Ok(false);
        }
        self.statistics.affine_opening_certificate = true;
        measure_search_phase::<DETAILED, _>(&mut search_timings, "certificate_install", || {
            self.install_diagonal_certificate(attempt.certificate, true)
        })?;
        Ok(true)
    }

    #[allow(clippy::too_many_lines)]
    fn activate_relaxed_linear<const DETAILED: bool>(
        &mut self,
        mut search_timings: Option<&mut PhaseTimingsBuilder>,
        mut search_counters: Option<&mut PoastaSearchWorkCounters>,
    ) -> Result<(), PoaError> {
        let policy_started = search_phase_started::<DETAILED>(&search_timings);
        self.relaxed_linear_attempted = true;
        let source_path_eligible = self.eager_relaxed_linear
            && nonnegative_score(self.costs.mismatch()) > 0
            && self.opening_transition()? > 0;
        record_search_phase::<DETAILED>(
            search_timings.as_deref_mut(),
            "policy_selection",
            policy_started,
        );
        if source_path_eligible {
            if let Some(certificate) = measure_search_phase::<DETAILED, _>(
                &mut search_timings,
                "source_path_certificate",
                || exact_source_path_certificate(self.graph, self.query, self.topological),
            )? {
                if DETAILED {
                    if let Some(counters) = search_counters.as_deref_mut() {
                        counters.record_source_path_certificate()?;
                    }
                }
                measure_search_phase::<DETAILED, _>(
                    &mut search_timings,
                    "certificate_install",
                    || self.install_diagonal_certificate(certificate, false),
                )?;
                return Ok(());
            }
        }
        let policy_started = search_phase_started::<DETAILED>(&search_timings);
        let tighten_exact_affine = self.should_tighten_exact_affine_bypass();
        let prefer_exact_affine = self.prefer_exact_affine_certificate();
        record_search_phase::<DETAILED>(
            search_timings.as_deref_mut(),
            "policy_selection",
            policy_started,
        );
        if tighten_exact_affine
            && !prefer_exact_affine
            && nonnegative_score(self.costs.mismatch()) > 0
            && nonnegative_score(self.costs.gap_extension()) > 0
            && self.try_early_diagonal_certificate::<DETAILED>(search_timings.as_deref_mut())?
        {
            if DETAILED {
                if let Some(counters) = search_counters.as_deref_mut() {
                    counters.record_early_diagonal_certificate()?;
                }
            }
            self.record_relaxed_certificate_outcome(true);
            self.frontier = BucketFrontier::default();
            return Ok(());
        }
        if prefer_exact_affine {
            if tighten_exact_affine {
                self.prepare_exact_diagonal_certificate::<DETAILED>(search_timings.as_deref_mut())?;
                self.exact_affine_upper_bound = self
                    .exact_diagonal_attempt
                    .take()
                    .map(|attempt| attempt.certificate.optimum_cost);
            }
            let exact_cells_before = if DETAILED {
                search_counters
                    .as_deref()
                    .map_or(Ok(0), PoastaSearchWorkCounters::filled_exact_affine_cells)?
            } else {
                0
            };
            if self.try_exact_affine_certificate::<DETAILED>(
                search_timings.as_deref_mut(),
                search_counters.as_deref_mut(),
            )? {
                if DETAILED {
                    if let Some(counters) = search_counters.as_deref_mut() {
                        let exact_cells = counters
                            .filled_exact_affine_cells()?
                            .checked_sub(exact_cells_before)
                            .ok_or(PoaError::ArithmeticOverflow {
                                context: "POASTA detailed exact-affine bypass-cell subtraction",
                            })?;
                        counters.record_exact_affine_route(
                            ExactAffineRoute::Bypass,
                            exact_cells,
                            true,
                        )?;
                    }
                }
                self.record_exact_affine_bypass();
                self.frontier = BucketFrontier::default();
                return Ok(());
            }
            if DETAILED {
                if let Some(counters) = search_counters.as_deref_mut() {
                    let exact_cells = counters
                        .filled_exact_affine_cells()?
                        .checked_sub(exact_cells_before)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "POASTA detailed exact-affine bypass-cell subtraction",
                        })?;
                    counters.record_exact_affine_route(
                        ExactAffineRoute::Bypass,
                        exact_cells,
                        false,
                    )?;
                }
            }
        }
        if let Some(heuristic) = RelaxedLinearHeuristic::build::<A, DETAILED>(
            self.graph,
            self.query,
            self.topological,
            self.index,
            self.costs,
            self.workspace.as_deref_mut(),
            search_timings.as_deref_mut(),
            search_counters.as_deref_mut(),
        )? {
            self.statistics.linear_heuristic_cells = heuristic.cell_count();
            self.statistics.relaxed_linear_source_bound = Some(
                heuristic.lower_bound(StateKey::new(self.index.source(), 0, AffineState::Match)),
            );
            self.relaxed_linear_heuristic = Some(heuristic);
        }
        let certificate = measure_search_phase::<DETAILED, _>(
            &mut search_timings,
            "relaxed_linear_certificate",
            || match self.relaxed_linear_heuristic.as_ref() {
                Some(heuristic) => heuristic.diagonal_certificate(
                    self.graph,
                    self.query,
                    self.index,
                    nonnegative_score(self.costs.mismatch()),
                ),
                None => Ok(None),
            },
        )?;
        if let Some(certificate) = certificate {
            if DETAILED {
                if let Some(counters) = search_counters.as_deref_mut() {
                    counters.record_relaxed_linear_table_outcome(
                        self.statistics.linear_heuristic_cells,
                        RelaxedLinearTableOutcome::DiagonalCertificate,
                    )?;
                }
            }
            self.record_relaxed_certificate_outcome(true);
            measure_search_phase::<DETAILED, _>(
                &mut search_timings,
                "certificate_install",
                || self.install_diagonal_certificate(certificate, true),
            )?;
            self.frontier = BucketFrontier::default();
            return Ok(());
        }
        if self.relaxed_linear_heuristic.is_some() {
            let affine_opening =
                self.try_affine_opening_certificate::<DETAILED>(search_timings.as_deref_mut())?;
            if affine_opening {
                if DETAILED {
                    if let Some(counters) = search_counters.as_deref_mut() {
                        counters.record_relaxed_linear_table_outcome(
                            self.statistics.linear_heuristic_cells,
                            RelaxedLinearTableOutcome::AffineOpeningCertificate,
                        )?;
                    }
                }
                self.record_relaxed_certificate_outcome(true);
                self.frontier = BucketFrontier::default();
                return Ok(());
            }
            if DETAILED {
                if let Some(counters) = search_counters.as_deref_mut() {
                    counters.record_relaxed_linear_table_outcome(
                        self.statistics.linear_heuristic_cells,
                        RelaxedLinearTableOutcome::Miss,
                    )?;
                }
            }
            self.record_relaxed_certificate_outcome(false);
            let exact_cells_before = if DETAILED {
                search_counters
                    .as_deref()
                    .map_or(Ok(0), PoastaSearchWorkCounters::filled_exact_affine_cells)?
            } else {
                0
            };
            if self.try_exact_affine_certificate::<DETAILED>(
                search_timings,
                search_counters.as_deref_mut(),
            )? {
                if DETAILED {
                    if let Some(counters) = search_counters.as_deref_mut() {
                        let exact_cells = counters
                            .filled_exact_affine_cells()?
                            .checked_sub(exact_cells_before)
                            .ok_or(PoaError::ArithmeticOverflow {
                                context:
                                    "POASTA detailed post-relaxed exact-affine cell subtraction",
                            })?;
                        counters.record_exact_affine_route(
                            ExactAffineRoute::AfterRelaxed,
                            exact_cells,
                            true,
                        )?;
                    }
                }
                self.frontier = BucketFrontier::default();
                return Ok(());
            }
            if DETAILED {
                if let Some(counters) = search_counters {
                    let exact_cells = counters
                        .filled_exact_affine_cells()?
                        .checked_sub(exact_cells_before)
                        .ok_or(PoaError::ArithmeticOverflow {
                            context: "POASTA detailed post-relaxed exact-affine cell subtraction",
                        })?;
                    counters.record_exact_affine_route(
                        ExactAffineRoute::AfterRelaxed,
                        exact_cells,
                        false,
                    )?;
                }
            }
        }

        let start = StateKey::new(self.index.source(), 0, AffineState::Match);
        self.seed_start_record(start);
        let mut pending = Vec::new();
        try_reserve_exact(
            &mut pending,
            self.records.len(),
            "POASTA relaxed-heuristic reprioritization",
        )?;
        pending.extend(self.records.iter().filter_map(|(&state_index, record)| {
            (record.expanded_cost == UNEXPANDED_COST || record.cost < record.expanded_cost)
                .then_some((self.key_from_state_index(state_index), record.cost))
        }));
        self.frontier = BucketFrontier::default();
        for (key, cost) in pending {
            if !self.is_strictly_dominated_open_gap(key, cost) {
                self.push_frontier(key, cost)?;
            }
        }
        Ok(())
    }

    fn install_diagonal_certificate(
        &mut self,
        certificate: DiagonalCertificate,
        relaxed_linear_certificate: bool,
    ) -> Result<(), PoaError> {
        let DiagonalCertificate { path, optimum_cost } = certificate;
        debug_assert_eq!(path.len(), self.query.len());
        let compact_traceback = self.can_install_compact_certificate();
        let mut predecessor = StateKey::new(self.index.source(), 0, AffineState::Match);
        if compact_traceback {
            let goal = path.last().copied().unwrap_or(predecessor);
            if !self.is_goal(goal) {
                return self.traceback_error(goal);
            }
            self.certified_traceback = Some(path);
            self.statistics.relaxed_linear_certificate = relaxed_linear_certificate;
            self.best_goal = Some((optimum_cost, goal));
            return Ok(());
        }
        let mut cost = 0;
        for (query_prefix, &key) in path.iter().enumerate() {
            if key.state() != AffineState::Match || key.query_prefix() != query_prefix + 1 {
                return self.traceback_error(key);
            }
            let node = self
                .index
                .node_for_cursor(key.cursor)
                .ok_or_else(|| self.traceback_error_value(key))?;
            let graph_symbol = self
                .graph
                .node(node)
                .ok_or(PoaError::UnknownNode { node })?
                .symbol();
            let query_symbol = self
                .query
                .get(query_prefix)
                .ok_or(PoaError::DimensionOverflow {
                    graph_nodes: self.graph.node_count(),
                    query_len: self.query.len(),
                })?;
            if graph_symbol != query_symbol {
                cost = checked_cost_add(
                    cost,
                    nonnegative_score(self.costs.mismatch()),
                    "a certified diagonal mismatch",
                )?;
            }
            self.records.insert(
                self.state_index(key),
                StateRecord {
                    cost,
                    predecessor,
                    expanded_cost: UNEXPANDED_COST,
                    tie_transition: 0,
                },
            );
            predecessor = key;
        }
        let goal = predecessor;
        if cost != optimum_cost || !self.is_goal(goal) {
            return self.traceback_error(goal);
        }
        self.statistics.relaxed_linear_certificate = relaxed_linear_certificate;
        self.best_goal = Some((cost, goal));
        Ok(())
    }

    fn install_exact_affine_certificate(
        &mut self,
        certificate: &ExactAffineCertificate,
    ) -> Result<(), PoaError> {
        let compact_traceback = self.can_install_compact_certificate();
        let mut certified_traceback = Vec::new();
        if compact_traceback {
            try_reserve_exact(
                &mut certified_traceback,
                certificate.steps.len(),
                "POASTA compact certified traceback",
            )?;
        }
        let mut predecessor = StateKey::new(self.index.source(), 0, AffineState::Match);
        let mut cost = 0;
        for &step in &certificate.steps {
            cost = checked_cost_add(
                cost,
                step.transition_cost,
                "a certified exact-affine transition",
            )?;
            if compact_traceback {
                certified_traceback.push(step.key);
            } else {
                self.records.insert(
                    self.state_index(step.key),
                    StateRecord {
                        cost,
                        predecessor,
                        expanded_cost: UNEXPANDED_COST,
                        tie_transition: step.tie_transition,
                    },
                );
            }
            predecessor = step.key;
        }
        if cost != certificate.optimum_cost || !self.is_goal(predecessor) {
            return self.traceback_error(predecessor);
        }
        if compact_traceback {
            self.certified_traceback = Some(certified_traceback);
        }
        self.statistics.exact_affine_certificate = true;
        self.best_goal = Some((cost, predecessor));
        Ok(())
    }

    fn can_install_compact_certificate(&self) -> bool {
        let start = StateKey::new(self.index.source(), 0, AffineState::Match);
        self.records.is_empty()
            || (self.records.len() == 1 && self.records.contains_key(&self.state_index(start)))
    }

    fn expand_from_astar(&mut self, entry: QueueEntry) -> Result<(), PoaError> {
        if !self.claim(entry.key, entry.cost) {
            self.statistics.stale_entries =
                checked_increment(self.statistics.stale_entries, "POASTA stale entry count")?;
            return Ok(());
        }
        if self.is_strictly_dominated_open_gap(entry.key, entry.cost) {
            self.statistics.affine_dominated_states = checked_increment(
                self.statistics.affine_dominated_states,
                "POASTA affine-dominated state count",
            )?;
            return Ok(());
        }
        if self.should_prune(entry.key, entry.cost)? {
            self.statistics.bubble_pruned_states = checked_increment(
                self.statistics.bubble_pruned_states,
                "POASTA pruned state count",
            )?;
            return Ok(());
        }
        self.statistics.expanded_states = checked_increment(
            self.statistics.expanded_states,
            "POASTA expanded state count",
        )?;
        self.statistics.astar_expanded_states = checked_increment(
            self.statistics.astar_expanded_states,
            "POASTA A-star expanded state count",
        )?;
        self.observe_state(entry.key, entry.cost)?;
        self.relax_all_edges(entry.key, entry.cost)
    }

    #[allow(clippy::too_many_lines)]
    fn expand_with_dfs(&mut self, entry: QueueEntry) -> Result<(), PoaError> {
        let mut stack = core::mem::take(&mut self.dfs_stack);
        stack.clear();
        let result = self.expand_with_dfs_stack(entry, &mut stack);
        self.dfs_stack = stack;
        result
    }

    #[allow(clippy::too_many_lines)]
    fn expand_with_dfs_stack(
        &mut self,
        entry: QueueEntry,
        stack: &mut Vec<DfsFrame>,
    ) -> Result<(), PoaError> {
        stack.push(DfsFrame::new(entry.key, entry.cost));
        self.statistics.peak_dfs_entries = self.statistics.peak_dfs_entries.max(stack.len());

        while !stack.is_empty() {
            let frame_index = stack.len() - 1;
            let frame = stack[frame_index];
            if !frame.entered {
                if !self.claim(frame.key, frame.cost) {
                    self.statistics.stale_entries = checked_increment(
                        self.statistics.stale_entries,
                        "POASTA stale DFS entry count",
                    )?;
                    stack.pop();
                    continue;
                }
                if self.should_prune(frame.key, frame.cost)? {
                    self.statistics.bubble_pruned_states = checked_increment(
                        self.statistics.bubble_pruned_states,
                        "POASTA DFS-pruned state count",
                    )?;
                    stack.pop();
                    continue;
                }
                self.statistics.expanded_states = checked_increment(
                    self.statistics.expanded_states,
                    "POASTA expanded state count",
                )?;
                self.statistics.dfs_expanded_states = checked_increment(
                    self.statistics.dfs_expanded_states,
                    "POASTA DFS expanded state count",
                )?;
                self.observe_state(frame.key, frame.cost)?;
                stack[frame_index].entered = true;
            }

            let frame = stack[frame_index];
            let successors = self.index.successors(frame.key.cursor);
            if let Some(&successor) = successors.get(frame.next_successor) {
                stack[frame_index].next_successor =
                    checked_increment(frame.next_successor, "POASTA DFS successor index")?;
                if successor == self.index.sink() {
                    continue;
                }
                let node = self
                    .index
                    .node_for_cursor(successor)
                    .ok_or_else(|| self.traceback_error_value(frame.key))?;
                if frame.key.query_prefix() == self.query.len() {
                    self.relax_deletion_open(frame.key, frame.cost, successor)?;
                    continue;
                }
                let query_symbol = self.query.get(frame.key.query_prefix()).ok_or(
                    PoaError::DimensionOverflow {
                        graph_nodes: self.graph.node_count(),
                        query_len: self.query.len(),
                    },
                )?;
                let target =
                    StateKey::new(successor, frame.key.query_prefix() + 1, AffineState::Match);
                if self
                    .graph
                    .node(node)
                    .ok_or(PoaError::UnknownNode { node })?
                    .symbol()
                    == query_symbol
                {
                    let relaxation = self.relax_record(target, frame.cost, frame.key, 0)?;
                    if relaxation == Relaxation::Improved {
                        self.statistics.greedy_exact_matches = checked_increment(
                            self.statistics.greedy_exact_matches,
                            "POASTA greedy exact-match count",
                        )?;
                        stack[frame_index].found_new_exact_match = true;
                        if stack.len() == stack.capacity() {
                            stack
                                .try_reserve(1)
                                .map_err(|_| PoaError::AllocationFailed {
                                    buffer: "POASTA exact-match DFS stack",
                                    elements: 1,
                                })?;
                        }
                        stack.push(DfsFrame::new(target, frame.cost));
                        self.statistics.peak_dfs_entries =
                            self.statistics.peak_dfs_entries.max(stack.len());
                    } else {
                        stack[frame_index].requires_insertion = true;
                        self.relax_deletion_open(frame.key, frame.cost, successor)?;
                    }
                } else {
                    let cost = checked_cost_add(
                        frame.cost,
                        nonnegative_score(self.costs.mismatch()),
                        "a mismatch transition",
                    )?;
                    self.relax_to_frontier(target, cost, frame.key, 0)?;
                    stack[frame_index].requires_insertion = true;
                    self.relax_deletion_open(frame.key, frame.cost, successor)?;
                }
                continue;
            }

            if frame.key.query_prefix() < self.query.len()
                && (!frame.found_new_exact_match || frame.requires_insertion)
            {
                self.relax_insertion_open(frame.key, frame.cost)?;
            }
            stack.pop();
        }
        Ok(())
    }

    #[allow(clippy::too_many_lines)]
    fn relax_all_edges(&mut self, key: StateKey, cost: Cost) -> Result<(), PoaError> {
        match key.state() {
            AffineState::Match => {
                if key.query_prefix() < self.query.len() {
                    let query_symbol =
                        self.query
                            .get(key.query_prefix())
                            .ok_or(PoaError::DimensionOverflow {
                                graph_nodes: self.graph.node_count(),
                                query_len: self.query.len(),
                            })?;
                    let successors = self.index.successors(key.cursor);
                    for &successor in successors {
                        if successor == self.index.sink() {
                            continue;
                        }
                        let node = self
                            .index
                            .node_for_cursor(successor)
                            .ok_or_else(|| self.traceback_error_value(key))?;
                        let graph_symbol = self
                            .graph
                            .node(node)
                            .ok_or(PoaError::UnknownNode { node })?
                            .symbol();
                        let transition_cost = if graph_symbol == query_symbol {
                            0
                        } else {
                            nonnegative_score(self.costs.mismatch())
                        };
                        self.relax_to_frontier(
                            StateKey::new(successor, key.query_prefix() + 1, AffineState::Match),
                            checked_cost_add(cost, transition_cost, "a diagonal transition")?,
                            key,
                            0,
                        )?;
                    }
                }
                self.relax_gap_open_edges(key, cost)
            }
            AffineState::Insertion => {
                if key.query_prefix() < self.query.len() {
                    self.relax_to_frontier(
                        StateKey::new(key.cursor, key.query_prefix() + 1, AffineState::Insertion),
                        checked_cost_add(
                            cost,
                            nonnegative_score(self.costs.gap_extension()),
                            "an insertion extension",
                        )?,
                        key,
                        0,
                    )?;
                }
                self.relax_to_frontier(
                    StateKey::new(key.cursor, key.query_prefix(), AffineState::Match),
                    cost,
                    key,
                    2,
                )?;
                Ok(())
            }
            AffineState::Deletion => {
                let successors = self.index.successors(key.cursor);
                for &successor in successors {
                    if successor == self.index.sink() {
                        continue;
                    }
                    self.relax_to_frontier(
                        StateKey::new(successor, key.query_prefix(), AffineState::Deletion),
                        checked_cost_add(
                            cost,
                            nonnegative_score(self.costs.gap_extension()),
                            "a deletion extension",
                        )?,
                        key,
                        0,
                    )?;
                }
                self.relax_to_frontier(
                    StateKey::new(key.cursor, key.query_prefix(), AffineState::Match),
                    cost,
                    key,
                    1,
                )?;
                Ok(())
            }
        }
    }

    fn relax_gap_open_edges(&mut self, key: StateKey, cost: Cost) -> Result<(), PoaError> {
        debug_assert_eq!(key.state(), AffineState::Match);
        if key.query_prefix() < self.query.len() {
            self.relax_insertion_open(key, cost)?;
        }
        let successors = self.index.successors(key.cursor);
        for &successor in successors {
            if successor == self.index.sink() {
                continue;
            }
            self.relax_deletion_open(key, cost, successor)?;
        }
        Ok(())
    }

    fn relax_insertion_open(&mut self, key: StateKey, cost: Cost) -> Result<(), PoaError> {
        self.relax_to_frontier(
            StateKey::new(key.cursor, key.query_prefix() + 1, AffineState::Insertion),
            checked_cost_add(cost, self.opening_transition()?, "an insertion opening")?,
            key,
            1,
        )?;
        Ok(())
    }

    fn relax_deletion_open(
        &mut self,
        key: StateKey,
        cost: Cost,
        successor: usize,
    ) -> Result<(), PoaError> {
        self.relax_to_frontier(
            StateKey::new(successor, key.query_prefix(), AffineState::Deletion),
            checked_cost_add(cost, self.opening_transition()?, "a deletion opening")?,
            key,
            1,
        )?;
        Ok(())
    }

    fn opening_transition(&self) -> Result<Cost, PoaError> {
        checked_cost_add(
            nonnegative_score(self.costs.gap_opening()),
            nonnegative_score(self.costs.gap_extension()),
            "the affine opening transition",
        )
    }

    fn relax_to_frontier(
        &mut self,
        key: StateKey,
        cost: Cost,
        predecessor: StateKey,
        tie_transition: u8,
    ) -> Result<Relaxation, PoaError> {
        let relaxation = self.relax_record(key, cost, predecessor, tie_transition)?;
        if relaxation == Relaxation::Improved {
            if self.is_strictly_dominated_open_gap(key, cost) {
                self.statistics.affine_dominated_states = checked_increment(
                    self.statistics.affine_dominated_states,
                    "POASTA affine-dominated state count",
                )?;
            } else {
                self.push_frontier(key, cost)?;
            }
        }
        Ok(relaxation)
    }

    fn relax_record(
        &mut self,
        key: StateKey,
        cost: Cost,
        predecessor: StateKey,
        tie_transition: u8,
    ) -> Result<Relaxation, PoaError> {
        let state_index = self.state_index(key);
        if self.records.len() == self.records.capacity() {
            self.records
                .try_reserve(1)
                .map_err(|_| PoaError::AllocationFailed {
                    buffer: "POASTA sparse state records",
                    elements: 1,
                })?;
        }
        match self.records.entry(state_index) {
            Entry::Occupied(mut entry) => {
                let record = entry.get_mut();
                if cost < record.cost {
                    if record.expanded_cost != UNEXPANDED_COST {
                        self.statistics.reopened_states = checked_increment(
                            self.statistics.reopened_states,
                            "POASTA reopened state count",
                        )?;
                    }
                    record.cost = cost;
                    record.predecessor = predecessor;
                    record.tie_transition = tie_transition;
                    return Ok(Relaxation::Improved);
                }
                if cost == record.cost
                    && (tie_transition, predecessor.cursor)
                        < (record.tie_transition, record.predecessor.cursor)
                {
                    record.predecessor = predecessor;
                    record.tie_transition = tie_transition;
                    return Ok(Relaxation::TraceOnly);
                }
                Ok(Relaxation::Unchanged)
            }
            Entry::Vacant(entry) => {
                entry.insert(StateRecord {
                    cost,
                    predecessor,
                    expanded_cost: UNEXPANDED_COST,
                    tie_transition,
                });
                Ok(Relaxation::Improved)
            }
        }
    }

    fn push_frontier(&mut self, key: StateKey, cost: Cost) -> Result<(), PoaError> {
        let estimated_total =
            checked_cost_add(cost, self.heuristic(key)?, "an A-star priority value")?;
        self.frontier.try_push(QueueEntry {
            estimated_total,
            cost,
            key,
        })?;
        self.statistics.queue_entries =
            checked_increment(self.statistics.queue_entries, "POASTA queue-entry count")?;
        self.statistics.peak_queue_entries =
            self.statistics.peak_queue_entries.max(self.frontier.len());
        Ok(())
    }

    fn heuristic(&self, key: StateKey) -> Result<Cost, PoaError> {
        let minimum = self
            .index
            .min_remaining_nodes(key.cursor)
            .ok_or_else(|| self.traceback_error_value(key))?;
        let maximum = self
            .index
            .max_remaining_nodes(key.cursor)
            .ok_or_else(|| self.traceback_error_value(key))?;
        let query_remaining = self.query.len().checked_sub(key.query_prefix()).ok_or(
            PoaError::DimensionOverflow {
                graph_nodes: self.graph.node_count(),
                query_len: self.query.len(),
            },
        )?;
        let (gap_symbols, insertion_only, deletion_only) = if maximum < query_remaining {
            (query_remaining - maximum, true, false)
        } else if minimum > query_remaining {
            (minimum - query_remaining, false, true)
        } else {
            (0, false, false)
        };
        let length_bound = if gap_symbols == 0 {
            0
        } else {
            let extensions = checked_cost_mul(
                gap_symbols,
                nonnegative_score(self.costs.gap_extension()),
                "the minimum remaining gap extensions",
            )?;
            let opening_is_paid = (key.state() == AffineState::Insertion && insertion_only)
                || (key.state() == AffineState::Deletion && deletion_only);
            if opening_is_paid {
                extensions
            } else {
                checked_cost_add(
                    nonnegative_score(self.costs.gap_opening()),
                    extensions,
                    "the minimum remaining affine gap cost",
                )?
            }
        };
        let linear_bound = self
            .relaxed_linear_heuristic
            .as_ref()
            .map_or(0, |heuristic| heuristic.lower_bound(key));
        Ok(length_bound.max(linear_bound))
    }

    fn claim(&mut self, key: StateKey, cost: Cost) -> bool {
        let state_index = self.state_index(key);
        let Some(record) = self.records.get_mut(&state_index) else {
            return false;
        };
        if record.cost != cost
            || (record.expanded_cost != UNEXPANDED_COST && record.expanded_cost <= cost)
        {
            return false;
        }
        record.expanded_cost = cost;
        true
    }

    fn is_strictly_dominated_open_gap(&self, key: StateKey, cost: Cost) -> bool {
        if !matches!(self.options.backend, PoastaBackend::Auto) || key.state() == AffineState::Match
        {
            return false;
        }
        let match_key = StateKey::new(key.cursor, key.query_prefix(), AffineState::Match);
        let Some(match_record) = self.records.get(&self.state_index(match_key)) else {
            return false;
        };
        match_record
            .cost
            .checked_add(nonnegative_score(self.costs.gap_opening()))
            .is_some_and(|reopening_cost| reopening_cost < cost)
    }

    fn observe_state(&mut self, key: StateKey, cost: Cost) -> Result<(), PoaError> {
        if self.is_goal(key) {
            let replace = self.best_goal.is_none_or(|(best_cost, best_key)| {
                cost < best_cost || (cost == best_cost && key.cursor < best_key.cursor)
            });
            if replace {
                self.best_goal = Some((cost, key));
            }
        }
        if key.state() == AffineState::Match {
            for &bubble in self.index.bubbles_ending_at(key.cursor) {
                if self.reached_exits[bubble].insert(key.query_prefix(), cost)? {
                    self.statistics.reached_exit_states = checked_increment(
                        self.statistics.reached_exit_states,
                        "POASTA reached-exit state count",
                    )?;
                }
            }
        }
        Ok(())
    }

    fn is_goal(&self, key: StateKey) -> bool {
        if key.state() != AffineState::Match || key.query_prefix() != self.query.len() {
            return false;
        }
        if self.graph.node_count() == 0 {
            key.cursor == self.index.source()
        } else {
            key.cursor != self.index.source()
                && self.index.successors(key.cursor) == [self.index.sink()]
        }
    }

    fn should_prune(&mut self, key: StateKey, cost: Cost) -> Result<bool, PoaError> {
        if !self.options.superbubble_pruning {
            return Ok(false);
        }
        for membership in self.index.memberships(key.cursor) {
            self.statistics.bubble_tests =
                checked_increment(self.statistics.bubble_tests, "POASTA bubble-test count")?;
            let bubble = self.index.bubbles()[membership.bubble()];
            let Some(minimum_position) = key
                .query_prefix()
                .checked_add(membership.min_edges_to_exit())
            else {
                continue;
            };
            if minimum_position > self.query.len() {
                continue;
            }
            let maximum_position = key
                .query_prefix()
                .checked_add(membership.max_edges_to_exit())
                .map_or(self.query.len(), |position| position.min(self.query.len()));
            if maximum_position < minimum_position {
                continue;
            }
            let reached = &self.reached_exits[membership.bubble()];
            if reached.positions.is_empty() {
                continue;
            }
            if key.state() != AffineState::Match && open_gap_can_improve(cost, reached, self.costs)
            {
                continue;
            }
            let maximum_deletion = self.index.max_remaining_nodes(bubble.exit()).unwrap_or(0);
            let query_len = self.query.len();
            let reached = &mut self.reached_exits[membership.bubble()];
            if bubble_dominates(
                cost,
                minimum_position,
                maximum_position,
                reached,
                maximum_deletion,
                query_len,
                self.costs,
            )? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn traceback(
        &self,
        goal: StateKey,
    ) -> Result<(Vec<PoaOperation>, Vec<PoaAlignedPair>), PoaError> {
        let maximum_operations = self
            .graph
            .node_count()
            .checked_add(self.query.len())
            .ok_or(PoaError::DimensionOverflow {
                graph_nodes: self.graph.node_count(),
                query_len: self.query.len(),
            })?;
        let mut operations = Vec::new();
        try_reserve_exact(
            &mut operations,
            maximum_operations,
            "POASTA traceback operations",
        )?;
        let start = StateKey::new(self.index.source(), 0, AffineState::Match);
        if let Some(certified_traceback) = &self.certified_traceback {
            let mut predecessor = start;
            for &current in certified_traceback {
                if let Some(operation) = self.traceback_operation(current, predecessor)? {
                    operations.push(operation);
                }
                predecessor = current;
                if operations.len() > maximum_operations {
                    return self.traceback_error(current);
                }
            }
            if predecessor != goal {
                return self.traceback_error(predecessor);
            }
        } else {
            let mut current = goal;
            let mut links = 0_usize;
            while current != start {
                let record = self
                    .records
                    .get(&self.state_index(current))
                    .ok_or_else(|| self.traceback_error_value(current))?;
                let predecessor = record.predecessor;
                if predecessor.is_none() {
                    return self.traceback_error(current);
                }
                if let Some(operation) = self.traceback_operation(current, predecessor)? {
                    operations.push(operation);
                }
                current = predecessor;
                links = checked_increment(links, "POASTA traceback link count")?;
                if links > self.records.len() || operations.len() > maximum_operations {
                    return self.traceback_error(current);
                }
            }
            operations.reverse();
        }

        let aligned_count = operations
            .iter()
            .filter(|operation| matches!(operation, PoaOperation::Aligned { .. }))
            .count();
        let mut aligned_pairs = Vec::new();
        try_reserve_exact(
            &mut aligned_pairs,
            aligned_count,
            "POASTA aligned-pair projection",
        )?;
        aligned_pairs.extend(operations.iter().filter_map(|operation| match *operation {
            PoaOperation::Aligned { node, query_index } => {
                Some(PoaAlignedPair::new(node, query_index))
            }
            PoaOperation::GraphOnly { .. } | PoaOperation::QueryOnly { .. } => None,
        }));
        Ok((operations, aligned_pairs))
    }

    fn traceback_operation(
        &self,
        current: StateKey,
        predecessor: StateKey,
    ) -> Result<Option<PoaOperation>, PoaError> {
        match current.state() {
            AffineState::Match
                if predecessor.state() == AffineState::Match
                    && current.cursor != predecessor.cursor
                    && current.query_prefix() == predecessor.query_prefix().saturating_add(1) =>
            {
                let node = self
                    .index
                    .node_for_cursor(current.cursor)
                    .ok_or_else(|| self.traceback_error_value(current))?;
                Ok(Some(PoaOperation::Aligned {
                    node,
                    query_index: current.query_prefix() - 1,
                }))
            }
            AffineState::Match
                if predecessor.state() != AffineState::Match
                    && current.cursor == predecessor.cursor
                    && current.query_prefix() == predecessor.query_prefix() =>
            {
                Ok(None)
            }
            AffineState::Deletion
                if current.cursor != predecessor.cursor
                    && current.query_prefix() == predecessor.query_prefix() =>
            {
                let node = self
                    .index
                    .node_for_cursor(current.cursor)
                    .ok_or_else(|| self.traceback_error_value(current))?;
                Ok(Some(PoaOperation::GraphOnly { node }))
            }
            AffineState::Insertion
                if current.cursor == predecessor.cursor
                    && current.query_prefix() == predecessor.query_prefix().saturating_add(1) =>
            {
                Ok(Some(PoaOperation::QueryOnly {
                    query_index: current.query_prefix() - 1,
                }))
            }
            AffineState::Match | AffineState::Deletion | AffineState::Insertion => {
                self.traceback_error(current)
            }
        }
    }

    #[allow(clippy::too_many_lines)]
    fn memory_estimate(
        &self,
        operations: &[PoaOperation],
        aligned_pairs: &[PoaAlignedPair],
    ) -> Result<PoastaMemoryEstimate, PoaError> {
        let graph_index_bytes = self.index.modeled_payload_bytes()?;
        let retained_heuristic_payload_bytes = self
            .relaxed_linear_heuristic
            .as_ref()
            .map_or(Ok(0), RelaxedLinearHeuristic::modeled_payload_bytes)?;
        let cached_heuristic_payload_bytes = self
            .workspace
            .as_deref()
            .map_or(Ok(0), PoastaWorkspace::modeled_payload_bytes)?;
        let heuristic_payload_bytes = checked_add_usize(
            checked_add_usize(
                retained_heuristic_payload_bytes,
                cached_heuristic_payload_bytes,
                "POASTA retained and cached heuristic payload",
            )?,
            self.certificate_heuristic_payload_bytes,
            "POASTA heuristic peak payload",
        )?;
        let sparse_state_payload_bytes = checked_bytes(
            self.records.capacity(),
            size_of::<(usize, StateRecord)>(),
            "POASTA sparse state payload",
        )?;
        let certified_traceback_payload_bytes =
            self.certified_traceback
                .as_ref()
                .map_or(Ok(0), |traceback| {
                    checked_bytes(
                        traceback.capacity(),
                        size_of::<StateKey>(),
                        "POASTA certified traceback state payload",
                    )
                })?;
        let state_payload_bytes = checked_add_usize(
            sparse_state_payload_bytes,
            certified_traceback_payload_bytes,
            "POASTA sparse and certified state payload",
        )?;
        let frontier_payload_peak_bytes = self.frontier.modeled_payload_bytes()?;
        let dfs_payload_peak_bytes = checked_bytes(
            self.statistics.peak_dfs_entries,
            size_of::<DfsFrame>(),
            "POASTA DFS payload",
        )?;
        let reached_exit_payload_bytes = self.reached_exits.iter().try_fold(
            0_usize,
            |total, reached| -> Result<usize, PoaError> {
                let position_bytes = checked_bytes(
                    reached.positions.capacity(),
                    size_of::<ReachedPosition>(),
                    "POASTA reached-exit positions",
                )?;
                let insertion_bytes = checked_bytes(
                    reached.prefix_insertion_minimum.capacity(),
                    size_of::<Cost>(),
                    "POASTA reached-exit insertion index",
                )?;
                let deletion_bytes = checked_bytes(
                    reached.deletion_range_tree.capacity(),
                    size_of::<Option<Cost>>(),
                    "POASTA reached-exit deletion index",
                )?;
                total
                    .checked_add(position_bytes)
                    .and_then(|bytes| bytes.checked_add(insertion_bytes))
                    .and_then(|bytes| bytes.checked_add(deletion_bytes))
                    .ok_or(PoaError::ArithmeticOverflow {
                        context: "POASTA reached-exit payload",
                    })
            },
        )?;
        let traceback_output_bytes = checked_add_usize(
            checked_bytes(
                operations.len(),
                size_of::<PoaOperation>(),
                "POASTA operation output",
            )?,
            checked_bytes(
                aligned_pairs.len(),
                size_of::<PoaAlignedPair>(),
                "POASTA aligned-pair output",
            )?,
            "POASTA traceback output",
        )?;
        let modeled_peak_bytes = [
            graph_index_bytes,
            heuristic_payload_bytes,
            state_payload_bytes,
            frontier_payload_peak_bytes,
            dfs_payload_peak_bytes,
            reached_exit_payload_bytes,
            traceback_output_bytes,
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .ok_or(PoaError::ArithmeticOverflow {
            context: "POASTA modeled peak payload",
        })?;
        Ok(PoastaMemoryEstimate {
            graph_index_bytes,
            heuristic_payload_bytes,
            state_payload_bytes,
            frontier_payload_peak_bytes,
            dfs_payload_peak_bytes,
            reached_exit_payload_bytes,
            traceback_output_bytes,
            modeled_peak_bytes,
        })
    }

    fn traceback_error<T>(&self, key: StateKey) -> Result<T, PoaError> {
        Err(self.traceback_error_value(key))
    }

    fn traceback_error_value(&self, key: StateKey) -> PoaError {
        PoaError::AffineTracebackInvariant {
            node: self.index.node_for_cursor(key.cursor),
            column: key.query_prefix(),
            state: key.state().name(),
        }
    }

    fn state_index(&self, key: StateKey) -> usize {
        debug_assert!(key.query_prefix() < self.columns);
        // `possible_state_count` proved the complete product fits before the
        // search began. Sink transitions terminate at their predecessor and
        // therefore never require a stored sink state.
        debug_assert!(key.cursor <= self.graph.node_count());
        (key.cursor * self.columns * 3) + key.query_state
    }

    fn key_from_state_index(&self, state_index: usize) -> StateKey {
        let states_per_cursor = self.columns * 3;
        StateKey {
            cursor: state_index / states_per_cursor,
            query_state: state_index % states_per_cursor,
        }
    }
}

fn open_gap_can_improve(
    candidate_cost: Cost,
    reached: &ReachedExit,
    costs: PoastaAffineCosts,
) -> bool {
    reached.positions.iter().any(|position| {
        position
            .cost
            .checked_add(nonnegative_score(costs.gap_opening()))
            .is_none_or(|threshold| candidate_cost < threshold)
    })
}

fn bubble_dominates(
    candidate_cost: Cost,
    minimum_position: usize,
    maximum_position: usize,
    reached: &mut ReachedExit,
    maximum_deletion: usize,
    query_len: usize,
    costs: PoastaAffineCosts,
) -> Result<bool, PoaError> {
    if !position_is_dominated(
        candidate_cost,
        minimum_position,
        reached,
        maximum_deletion,
        query_len,
        costs,
    )? || !position_is_dominated(
        candidate_cost,
        maximum_position,
        reached,
        maximum_deletion,
        query_len,
        costs,
    )? {
        return Ok(false);
    }
    let start = reached
        .positions
        .partition_point(|position| position.query_prefix < minimum_position);
    let end = reached
        .positions
        .partition_point(|position| position.query_prefix <= maximum_position);
    for position_index in start..end {
        let query_prefix = reached.positions[position_index].query_prefix;
        for candidate in [
            query_prefix.checked_sub(1),
            Some(query_prefix),
            query_prefix.checked_add(1),
        ] {
            let Some(candidate) = candidate else {
                continue;
            };
            if candidate < minimum_position || candidate > maximum_position {
                continue;
            }
            if !position_is_dominated(
                candidate_cost,
                candidate,
                reached,
                maximum_deletion,
                query_len,
                costs,
            )? {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn position_is_dominated(
    candidate_cost: Cost,
    query_prefix: usize,
    reached: &mut ReachedExit,
    maximum_deletion: usize,
    query_len: usize,
    costs: PoastaAffineCosts,
) -> Result<bool, PoaError> {
    let upper_bound = reached.upper_bound(query_prefix, maximum_deletion, query_len, costs)?;
    Ok(upper_bound.is_some_and(|bound| bound <= candidate_cost))
}

#[cfg(test)]
fn implicit_gap_bound(
    base: Cost,
    length: usize,
    costs: PoastaAffineCosts,
) -> Result<Option<Cost>, PoaError> {
    if length == 0 {
        return Ok(Some(base));
    }
    let extensions = checked_cost_mul(
        length,
        nonnegative_score(costs.gap_extension()),
        "an implicit bubble-gap extension",
    )?;
    let gap = checked_cost_add(
        nonnegative_score(costs.gap_opening()),
        extensions,
        "an implicit bubble-gap cost",
    )?;
    checked_cost_add(base, gap, "an implicit bubble-gap upper bound").map(Some)
}

const fn min_cost(left: Cost, right: Cost) -> Cost {
    if left < right { left } else { right }
}

const fn optional_min(left: Option<Cost>, right: Option<Cost>) -> Option<Cost> {
    match (left, right) {
        (Some(left), Some(right)) => Some(min_cost(left, right)),
        (Some(value), None) | (None, Some(value)) => Some(value),
        (None, None) => None,
    }
}

const fn resolve_backend(requested: PoastaBackend) -> PoastaBackend {
    match requested {
        PoastaBackend::Sparse | PoastaBackend::Auto => PoastaBackend::Sparse,
        PoastaBackend::Dense => PoastaBackend::Dense,
    }
}

fn dense_statistics(total_possible_states: usize) -> PoastaSearchStatistics {
    PoastaSearchStatistics {
        total_possible_states,
        discovered_states: total_possible_states,
        expanded_states: total_possible_states,
        astar_expanded_states: 0,
        dfs_expanded_states: 0,
        queue_entries: 0,
        reopened_states: 0,
        stale_entries: 0,
        greedy_exact_matches: 0,
        affine_dominated_states: 0,
        linear_heuristic_cells: 0,
        relaxed_linear_source_bound: None,
        relaxed_linear_certificate: false,
        exact_diagonal_certificate_cells: 0,
        exact_diagonal_cost: None,
        affine_opening_certificate: false,
        exact_affine_certificate_cells: 0,
        exact_affine_certificate: false,
        bubble_tests: 0,
        bubble_pruned_states: 0,
        reached_exit_states: 0,
        superbubbles: 0,
        bubble_membership_entries: 0,
        peak_queue_entries: 0,
        peak_dfs_entries: 0,
    }
}

fn dense_memory_estimate(
    dense: AffinePoaAlignmentMemoryEstimate,
    operations: &[PoaOperation],
    aligned_pairs: &[PoaAlignedPair],
) -> Result<PoastaMemoryEstimate, PoaError> {
    let state_payload_bytes = dense
        .wide_score_matrix_bytes()
        .checked_add(dense.traceback_matrix_bytes())
        .ok_or(PoaError::ArithmeticOverflow {
            context: "POASTA dense state payload",
        })?;
    let traceback_output_bytes = checked_add_usize(
        checked_bytes(
            operations.len(),
            size_of::<PoaOperation>(),
            "POASTA dense operation output",
        )?,
        checked_bytes(
            aligned_pairs.len(),
            size_of::<PoaAlignedPair>(),
            "POASTA dense aligned-pair output",
        )?,
        "POASTA dense traceback output",
    )?;
    Ok(PoastaMemoryEstimate {
        graph_index_bytes: dense.node_index_bytes(),
        heuristic_payload_bytes: 0,
        state_payload_bytes,
        frontier_payload_peak_bytes: 0,
        dfs_payload_peak_bytes: 0,
        reached_exit_payload_bytes: 0,
        traceback_output_bytes,
        modeled_peak_bytes: dense.peak_buffer_bytes(),
    })
}

fn possible_state_count(graph_nodes: usize, query_len: usize) -> Result<usize, PoaError> {
    graph_nodes
        .checked_add(1)
        .and_then(|cursors| {
            query_len
                .checked_add(1)
                .and_then(|columns| cursors.checked_mul(columns))
        })
        .and_then(|cells| cells.checked_mul(3))
        .ok_or(PoaError::DimensionOverflow {
            graph_nodes,
            query_len,
        })
}

fn nonnegative_score(score: Score) -> Cost {
    debug_assert!(score >= 0);
    Cost::from(score.unsigned_abs())
}

fn checked_cost_add(left: Cost, right: Cost, context: &'static str) -> Result<Cost, PoaError> {
    left.checked_add(right)
        .filter(|cost| *cost != UNEXPANDED_COST)
        .ok_or(PoaError::PoastaCostOverflow { context })
}

fn checked_cost_mul(count: usize, cost: Cost, context: &'static str) -> Result<Cost, PoaError> {
    Cost::try_from(count)
        .map_err(|_| PoaError::PoastaCostOverflow { context })?
        .checked_mul(cost)
        .filter(|cost| *cost != UNEXPANDED_COST)
        .ok_or(PoaError::PoastaCostOverflow { context })
}

fn cost_to_score(cost: Cost) -> Result<Score, PoaError> {
    let minimum_magnitude = (Score::MAX as Cost) + 1;
    if cost == minimum_magnitude {
        return Ok(Score::MIN);
    }
    if cost <= Score::MAX as Cost {
        let magnitude = Score::try_from(cost).map_err(|_| PoaError::PoastaScoreOutOfRange {
            cost: u128::from(cost),
        })?;
        return Ok(-magnitude);
    }
    Err(PoaError::PoastaScoreOutOfRange {
        cost: u128::from(cost),
    })
}

fn checked_increment(value: usize, context: &'static str) -> Result<usize, PoaError> {
    value
        .checked_add(1)
        .ok_or(PoaError::ArithmeticOverflow { context })
}

fn checked_bytes(
    elements: usize,
    element_size: usize,
    context: &'static str,
) -> Result<usize, PoaError> {
    elements
        .checked_mul(element_size)
        .ok_or(PoaError::ArithmeticOverflow { context })
}

fn checked_add_usize(left: usize, right: usize, context: &'static str) -> Result<usize, PoaError> {
    left.checked_add(right)
        .ok_or(PoaError::ArithmeticOverflow { context })
}

fn try_reserve_exact<T>(
    values: &mut Vec<T>,
    additional: usize,
    buffer: &'static str,
) -> Result<(), PoaError> {
    values
        .try_reserve_exact(additional)
        .map_err(|_| PoaError::AllocationFailed {
            buffer,
            elements: additional,
        })
}

fn filled<T: Clone>(elements: usize, value: T, buffer: &'static str) -> Result<Vec<T>, PoaError> {
    let mut values = Vec::new();
    try_reserve_exact(&mut values, elements, buffer)?;
    values.resize(elements, value);
    Ok(values)
}

fn verification<T>(message: impl Into<String>) -> Result<T, PoaError> {
    Err(PoaError::VerificationFailed {
        message: message.into(),
    })
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, ProteinSequence, RnaSequence};
    use align_msa::SequenceRecord;

    use super::*;
    use crate::GlobalAffinePoa;

    fn costs() -> PoastaAffineCosts {
        PoastaAffineCosts::new(4, 6, 2).expect("valid costs")
    }

    #[test]
    fn specialized_diagonal_rows_match_the_scalar_recurrence() {
        let query_indices = [0, 1, 2, 3, 0, 1, 2, 3];
        let graph_symbol = 2;
        let mismatch = 4;
        let first_predecessor = ExactAffineRowBand {
            first: 0,
            last: 7,
            offset: 0,
        };
        let second_predecessor = ExactAffineRowBand {
            first: 2,
            last: 7,
            offset: 8,
        };
        let row = ExactAffineRowBand {
            first: 1,
            last: 8,
            offset: 14,
        };
        let predecessor_costs = [0, 4, 8, u16::MAX, 16, 20, 24, 28, 7, 3, 9, 1, u16::MAX, 5];

        let scalar_row = |predecessors: &[ExactAffineRowBand]| {
            let mut costs = predecessor_costs.to_vec();
            costs.resize(row.offset + exact_affine_row_cells(row), u16::MAX);
            for query_prefix in row.first.max(1)..=row.last {
                let predecessor = predecessors
                    .iter()
                    .map(|&predecessor| exact_affine_cost(&costs, predecessor, query_prefix - 1))
                    .min()
                    .unwrap_or(u16::MAX);
                let substitution = if query_indices[query_prefix - 1] == graph_symbol {
                    0
                } else {
                    mismatch
                };
                costs[row.offset + query_prefix - row.first] =
                    predecessor.saturating_add(substitution);
            }
            costs
        };

        let mut single = predecessor_costs.to_vec();
        single.resize(row.offset + exact_affine_row_cells(row), u16::MAX);
        let single_cells = fill_exact_diagonal_single_predecessor_row(
            row,
            first_predecessor,
            &query_indices,
            graph_symbol,
            mismatch,
            &mut single,
        );
        assert_eq!(single_cells, exact_affine_row_cells(row));
        assert_eq!(single, scalar_row(&[first_predecessor]));

        let mut branching = predecessor_costs.to_vec();
        branching.resize(row.offset + exact_affine_row_cells(row), u16::MAX);
        assert_eq!(
            merge_exact_diagonal_predecessor_row(row, first_predecessor, &mut branching),
            8
        );
        assert_eq!(
            merge_exact_diagonal_predecessor_row(row, second_predecessor, &mut branching),
            6
        );
        assert_eq!(
            finalize_exact_diagonal_row(
                row,
                &query_indices,
                graph_symbol,
                mismatch,
                &mut branching,
            ),
            exact_affine_row_cells(row)
        );
        assert_eq!(
            branching,
            scalar_row(&[first_predecessor, second_predecessor])
        );
    }

    #[test]
    fn relaxed_linear_simd_insertion_suffix_matches_scalar_recurrence() {
        for &(len, extension, initial_insertion) in &[
            (16, 2, 17),
            (31, 2, 13),
            (32, 0, u16::MAX),
            (49, 4_000, 63_000),
        ] {
            let bases = (0..len)
                .map(|index| {
                    if index % 13 == 0 {
                        u16::MAX
                    } else {
                        u16::try_from((index * 37 + 11) % 257).expect("test value fits u16")
                    }
                })
                .collect::<Vec<_>>();
            let expected_cells = len / 16 * 16;
            let mut expected = bases.clone();
            let mut expected_insertion = initial_insertion;
            let mut expected_reachable = true;
            for value in expected[len - expected_cells..].iter_mut().rev() {
                let best = (*value).min(expected_insertion);
                expected_reachable &= best != u16::MAX;
                *value = best;
                expected_insertion = best.saturating_add(extension);
            }

            let mut actual = bases;
            let (processed, actual_insertion, actual_reachable) = RelaxedLinearInsertionScan::new(
                extension,
            )
            .fill(&mut actual, extension, initial_insertion);
            if processed == 0 {
                continue;
            }
            assert_eq!(processed, expected_cells);
            assert_eq!(actual, expected);
            assert_eq!(actual_insertion, expected_insertion);
            assert_eq!(actual_reachable, expected_reachable);
        }
    }

    #[test]
    fn relaxed_linear_simd_insertion_suffix_reports_unreachable_cells() {
        let mut values = [u16::MAX; 16];
        let (processed, _, reachable) =
            RelaxedLinearInsertionScan::new(2).fill(&mut values, 2, u16::MAX);
        if processed > 0 {
            assert_eq!(processed, values.len());
            assert!(!reachable);
        }
    }

    #[test]
    fn exact_affine_simd_insertion_prefix_matches_scalar_recurrence() {
        for &(len, opening, extension, initial_base, initial_insertion) in &[
            (16, 8, 2, u16::MAX, u16::MAX),
            (31, 8, 2, 13, 17),
            (32, 6, 0, 7, u16::MAX),
            (48, 5_000, 4_000, 62_000, 63_000),
        ] {
            let bases = (0..len)
                .map(|index| {
                    if index % 13 == 0 {
                        u16::MAX
                    } else {
                        u16::try_from((index * 37 + 11) % 257).expect("test value fits u16")
                    }
                })
                .collect::<Vec<_>>();
            let mut expected = bases.clone();
            let mut expected_base = initial_base;
            let mut expected_insertion = initial_insertion;
            let expected_cells = len / 16 * 16;
            for value in expected.iter_mut().take(expected_cells) {
                let base = *value;
                expected_insertion = expected_base
                    .saturating_add(opening)
                    .min(expected_insertion.saturating_add(extension));
                *value = base.min(expected_insertion);
                expected_base = base;
            }

            let mut actual = bases;
            let (processed, actual_base, actual_insertion) =
                ExactAffineInsertionScan::new(extension).fill(
                    &mut actual,
                    opening,
                    extension,
                    initial_base,
                    initial_insertion,
                );
            if processed == 0 {
                continue;
            }
            assert_eq!(processed, expected_cells);
            assert_eq!(actual, expected);
            assert_eq!(actual_base, expected_base);
            assert_eq!(actual_insertion, expected_insertion);
        }
    }

    #[test]
    fn exact_affine_many_base_avx2_matches_scalar_recurrence() {
        let opening = 6;
        let extension = 2;
        let mismatch = 4;
        let graph_symbol = 2;
        for len in [1_usize, 2, 15, 16, 17, 31, 32, 33, 65] {
            let predecessor_matches = (0..len)
                .map(|index| {
                    if index % 19 == 0 {
                        u16::MAX
                    } else {
                        u16::try_from((index * 37 + 11) % 257).expect("test value fits u16")
                    }
                })
                .collect::<Vec<_>>();
            let predecessor_deletions = (0..len)
                .map(|index| {
                    if index % 23 == 0 {
                        u16::MAX
                    } else {
                        u16::try_from((index * 29 + 17) % 251).expect("test value fits u16")
                    }
                })
                .collect::<Vec<_>>();
            let query_symbols = (0..len.saturating_sub(1))
                .map(|index| u8::try_from(index % 4).expect("test symbol fits u8"))
                .collect::<Vec<_>>();
            let mut expected_matches = predecessor_matches.clone();
            let mut expected_deletions = predecessor_deletions.clone();
            for index in (1..len).rev() {
                let deletion = predecessor_matches[index]
                    .saturating_add(opening)
                    .min(predecessor_deletions[index].saturating_add(extension));
                let substitution = if query_symbols[index - 1] == graph_symbol {
                    0
                } else {
                    mismatch
                };
                expected_matches[index] = predecessor_matches[index - 1]
                    .saturating_add(substitution)
                    .min(deletion);
                expected_deletions[index] = deletion;
            }

            let mut actual_matches = predecessor_matches;
            let mut actual_deletions = predecessor_deletions;
            let processed = ExactAffineInsertionScan::new(extension).fill_many_bases(
                &mut actual_matches,
                &mut actual_deletions,
                &query_symbols,
                graph_symbol,
                mismatch,
                opening,
                extension,
            );
            for index in (1..len - processed).rev() {
                let deletion = actual_matches[index]
                    .saturating_add(opening)
                    .min(actual_deletions[index].saturating_add(extension));
                let substitution = if query_symbols[index - 1] == graph_symbol {
                    0
                } else {
                    mismatch
                };
                actual_matches[index] = actual_matches[index - 1]
                    .saturating_add(substitution)
                    .min(deletion);
                actual_deletions[index] = deletion;
            }
            assert_eq!(actual_matches, expected_matches, "len={len}");
            assert_eq!(actual_deletions, expected_deletions, "len={len}");
            assert_eq!(processed % 16, 0);
            assert!(processed <= len.saturating_sub(1));
        }
    }

    fn dna_record(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(
            identifier,
            DnaSequence::from_ascii(sequence.as_bytes()).expect("valid DNA"),
        )
        .expect("valid record")
    }

    fn dna_graph(sequence: &str) -> PoaGraph<Dna> {
        let source = dna_record("backbone", sequence);
        let mut graph = PoaGraph::new();
        graph
            .append_record(&source, &vec![None; source.sequence().len()])
            .expect("seed graph");
        graph
    }

    fn assert_all_modes_equal_oracle(graph: &PoaGraph<Dna>, query: &DnaSequence) {
        let oracle = GlobalAffinePoa::from_poasta_costs(costs())
            .align(graph, query)
            .expect("oracle alignment");
        for options in [
            PoastaSearchOptions::astar_baseline(),
            PoastaSearchOptions::astar_baseline().with_exact_match_dfs(true),
            PoastaSearchOptions::astar_baseline().with_superbubble_pruning(true),
            PoastaSearchOptions::published(),
            PoastaSearchOptions::dense(),
            PoastaSearchOptions::optimized(),
        ] {
            let alignment = Poasta::with_options(costs(), options)
                .align(graph, query)
                .expect("POASTA alignment");
            assert_eq!(alignment.score(), oracle.score(), "options={options:?}");
            assert_eq!(
                verify_poasta_alignment(&alignment, graph, query, costs())
                    .expect("independent replay")
                    .recomputed_score(),
                oracle.score()
            );
        }
    }

    #[allow(clippy::too_many_lines)]
    fn exact_remaining_cost(
        graph: &PoaGraph<Dna>,
        query: &DnaSequence,
        index: &PoastaGraphIndex,
        key: StateKey,
        memo: &mut HashMap<StateKey, Cost>,
    ) -> Cost {
        if let Some(&cost) = memo.get(&key) {
            return cost;
        }
        let query_prefix = key.query_prefix();
        let terminal = index.successors(key.cursor) == [index.sink()];
        if key.state() == AffineState::Match && query_prefix == query.len() && terminal {
            memo.insert(key, 0);
            return 0;
        }

        let extension = nonnegative_score(costs().gap_extension());
        let opening = nonnegative_score(costs().gap_opening()) + extension;
        let mismatch = nonnegative_score(costs().mismatch());
        let mut minimum = Cost::MAX;
        match key.state() {
            AffineState::Match => {
                if query_prefix < query.len() {
                    let query_symbol = query.get(query_prefix).expect("query symbol");
                    for &successor in index.successors(key.cursor) {
                        if successor == index.sink() {
                            continue;
                        }
                        let node = index.node_for_cursor(successor).expect("graph cursor");
                        let transition =
                            if graph.node(node).expect("graph node").symbol() == query_symbol {
                                0
                            } else {
                                mismatch
                            };
                        let suffix = exact_remaining_cost(
                            graph,
                            query,
                            index,
                            StateKey::new(successor, query_prefix + 1, AffineState::Match),
                            memo,
                        );
                        minimum = minimum.min(transition + suffix);
                    }
                    let suffix = exact_remaining_cost(
                        graph,
                        query,
                        index,
                        StateKey::new(key.cursor, query_prefix + 1, AffineState::Insertion),
                        memo,
                    );
                    minimum = minimum.min(opening + suffix);
                }
                for &successor in index.successors(key.cursor) {
                    if successor == index.sink() {
                        continue;
                    }
                    let suffix = exact_remaining_cost(
                        graph,
                        query,
                        index,
                        StateKey::new(successor, query_prefix, AffineState::Deletion),
                        memo,
                    );
                    minimum = minimum.min(opening + suffix);
                }
            }
            AffineState::Insertion => {
                if query_prefix < query.len() {
                    let suffix = exact_remaining_cost(
                        graph,
                        query,
                        index,
                        StateKey::new(key.cursor, query_prefix + 1, AffineState::Insertion),
                        memo,
                    );
                    minimum = minimum.min(extension + suffix);
                }
                minimum = minimum.min(exact_remaining_cost(
                    graph,
                    query,
                    index,
                    StateKey::new(key.cursor, query_prefix, AffineState::Match),
                    memo,
                ));
            }
            AffineState::Deletion => {
                for &successor in index.successors(key.cursor) {
                    if successor == index.sink() {
                        continue;
                    }
                    let suffix = exact_remaining_cost(
                        graph,
                        query,
                        index,
                        StateKey::new(successor, query_prefix, AffineState::Deletion),
                        memo,
                    );
                    minimum = minimum.min(extension + suffix);
                }
                minimum = minimum.min(exact_remaining_cost(
                    graph,
                    query,
                    index,
                    StateKey::new(key.cursor, query_prefix, AffineState::Match),
                    memo,
                ));
            }
        }
        assert_ne!(minimum, Cost::MAX, "every test state reaches a terminal");
        memo.insert(key, minimum);
        minimum
    }

    #[test]
    fn empty_graph_and_query_boundaries_match_scalar_oracle() {
        let empty_graph = PoaGraph::<Dna>::new();
        for query in ["", "ACG"] {
            assert_all_modes_equal_oracle(
                &empty_graph,
                &DnaSequence::from_ascii(query.as_bytes()).expect("query"),
            );
        }
        let graph = dna_graph("ACG");
        assert_all_modes_equal_oracle(&graph, &DnaSequence::from_ascii(b"").expect("empty query"));
    }

    #[test]
    fn sparse_cost_accumulator_preserves_the_complete_public_score_range() {
        let costs = PoastaAffineCosts::new(1, Score::MAX - 3, 2).expect("valid wide costs");
        let alignment = Poasta::new(costs)
            .align(
                &PoaGraph::<Dna>::new(),
                &DnaSequence::from_ascii(b"AA").expect("query"),
            )
            .expect("the minimum public score is representable");
        assert_eq!(alignment.score(), Score::MIN);
    }

    #[test]
    fn identity_mismatch_and_gap_cases_match_scalar_oracle() {
        let graph = dna_graph("ACGT");
        for query in ["ACGT", "AGGT", "TACGT", "ACGTT", "AGT", "ACT"] {
            assert_all_modes_equal_oracle(
                &graph,
                &DnaSequence::from_ascii(query.as_bytes()).expect("query"),
            );
        }
    }

    #[test]
    fn recombining_bubble_matches_scalar_oracle_in_every_mode() {
        let mut graph = dna_graph("ACGT");
        graph
            .append_record(
                &dna_record("alternate", "ATGT"),
                &[
                    Some(PoaNodeId::new(0)),
                    None,
                    Some(PoaNodeId::new(2)),
                    Some(PoaNodeId::new(3)),
                ],
            )
            .expect("alternate path");
        for query in ["ACGT", "ATGT", "AGGT", "ACT", "ATTGT"] {
            assert_all_modes_equal_oracle(
                &graph,
                &DnaSequence::from_ascii(query.as_bytes()).expect("query"),
            );
        }
    }

    #[test]
    fn detailed_relaxed_counters_find_equivalent_single_transition_rows() {
        let mut graph = dna_graph("ACGT");
        graph
            .append_record(
                &dna_record("alternate", "ATGT"),
                &[
                    Some(PoaNodeId::new(0)),
                    None,
                    Some(PoaNodeId::new(2)),
                    Some(PoaNodeId::new(3)),
                ],
            )
            .expect("recombining alternate path");
        let query = DnaSequence::from_ascii(b"AGGT").expect("query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let index = PoastaGraphIndex::build(&graph, &topological, false).expect("graph index");
        let mut counters = PoastaSearchWorkCounters::default();
        let heuristic = RelaxedLinearHeuristic::build::<_, true>(
            &graph,
            &query,
            &topological,
            &index,
            costs(),
            None,
            None,
            Some(&mut counters),
        )
        .expect("relaxed heuristic range")
        .expect("relaxed heuristic");

        assert!(heuristic.cell_count() > 0);
        assert_eq!(
            counters.relaxed_linear_equivalent_single_transition_rows(),
            1
        );
        assert!(counters.relaxed_linear_equivalent_single_transition_cells() > 0);
        assert!(counters.relaxed_linear_equivalent_single_transition_insertion_cells() > 0);
        assert!(counters.relaxed_linear_equivalent_single_transition_successor_cell_visits() > 0);
    }

    #[test]
    fn detailed_relaxed_counters_find_equivalent_multiple_transition_rows() {
        let mut graph = dna_graph("AA");
        graph
            .append_record(
                &dna_record("first-second-allele", "AC"),
                &[Some(PoaNodeId::new(0)), None],
            )
            .expect("first successor allele");
        graph
            .append_record(
                &dna_record("second-first-allele", "CA"),
                &[None, Some(PoaNodeId::new(1))],
            )
            .expect("second predecessor allele");
        graph
            .append_record(
                &dna_record("both-alleles", "CC"),
                &[Some(PoaNodeId::new(3)), Some(PoaNodeId::new(2))],
            )
            .expect("complete recombination");
        let query = DnaSequence::from_ascii(b"GA").expect("query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let index = PoastaGraphIndex::build(&graph, &topological, false).expect("graph index");
        let mut counters = PoastaSearchWorkCounters::default();
        let heuristic = RelaxedLinearHeuristic::build::<_, true>(
            &graph,
            &query,
            &topological,
            &index,
            costs(),
            None,
            None,
            Some(&mut counters),
        )
        .expect("relaxed heuristic range")
        .expect("relaxed heuristic");

        assert!(heuristic.cell_count() > 0);
        assert!(counters.relaxed_linear_two_node_successor_rows() >= 2);
        assert!(counters.relaxed_linear_equivalent_multiple_transition_rows() >= 1);
        assert!(counters.relaxed_linear_equivalent_multiple_transition_cells() > 0);
        assert!(counters.relaxed_linear_equivalent_multiple_transition_insertion_cells() > 0);
        assert!(counters.relaxed_linear_equivalent_multiple_transition_successor_cell_visits() > 0);
    }

    #[test]
    fn nested_superbubbles_match_the_scalar_oracle_in_every_mode() {
        let mut graph = dna_graph("ACGTA");
        graph
            .append_record(
                &dna_record("inner-alternate", "ACCTA"),
                &[
                    Some(PoaNodeId::new(0)),
                    Some(PoaNodeId::new(1)),
                    None,
                    Some(PoaNodeId::new(3)),
                    Some(PoaNodeId::new(4)),
                ],
            )
            .expect("inner allele");
        graph
            .append_record(
                &dna_record("outer-alternate", "AGA"),
                &[Some(PoaNodeId::new(0)), None, Some(PoaNodeId::new(4))],
            )
            .expect("outer allele");

        for query in ["ACGTA", "ACCTA", "AGA", "AGGTA", "ACTA"] {
            assert_all_modes_equal_oracle(
                &graph,
                &DnaSequence::from_ascii(query.as_bytes()).expect("query"),
            );
        }
        let published = Poasta::new(costs())
            .align(&graph, &DnaSequence::from_ascii(b"ACGTA").expect("query"))
            .expect("published nested-bubble alignment");
        assert!(published.statistics().superbubbles() >= 2);
    }

    #[test]
    fn disconnected_roots_and_sinks_match_the_scalar_oracle() {
        let mut graph = dna_graph("AC");
        let disconnected = dna_record("disconnected", "GT");
        graph
            .append_record(&disconnected, &[None, None])
            .expect("disconnected component");
        for query in ["AC", "GT", "AG", "", "ACGT"] {
            assert_all_modes_equal_oracle(
                &graph,
                &DnaSequence::from_ascii(query.as_bytes()).expect("query"),
            );
        }
    }

    #[test]
    fn typed_rna_and_protein_share_the_same_kernel() {
        let rna_source = SequenceRecord::new(
            "rna",
            RnaSequence::from_ascii(b"ACGU").expect("RNA sequence"),
        )
        .expect("RNA record");
        let mut rna_graph = PoaGraph::new();
        rna_graph
            .append_record(&rna_source, &[None, None, None, None])
            .expect("RNA graph");
        let rna = Poasta::new(costs())
            .align(&rna_graph, rna_source.sequence())
            .expect("RNA alignment");
        assert_eq!(rna.score(), 0);

        let protein_source = SequenceRecord::new(
            "protein",
            ProteinSequence::from_ascii(b"MKWV").expect("protein sequence"),
        )
        .expect("protein record");
        let mut protein_graph = PoaGraph::new();
        protein_graph
            .append_record(&protein_source, &[None, None, None, None])
            .expect("protein graph");
        let protein = Poasta::new(costs())
            .align(&protein_graph, protein_source.sequence())
            .expect("protein alignment");
        assert_eq!(protein.score(), 0);
    }

    #[test]
    fn published_profile_reports_sparse_search_evidence() {
        let graph = dna_graph("ACGTACGT");
        let query = DnaSequence::from_ascii(b"ACGTACGT").expect("query");
        let alignment = Poasta::new(costs())
            .align(&graph, &query)
            .expect("alignment");
        assert_eq!(
            Poasta::new(costs()).profile_identifier(),
            POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_V1
        );
        assert!(alignment.statistics().greedy_exact_matches() >= query.len());
        assert!(alignment.statistics().discovered_states() > 0);
        assert!(alignment.memory_estimate().modeled_peak_bytes() > 0);
    }

    #[test]
    fn exact_match_dfs_strictly_reduces_discovery_and_queueing_on_a_long_identity() {
        let sequence = "ACGT".repeat(64);
        let graph = dna_graph(&sequence);
        let query = DnaSequence::from_ascii(sequence.as_bytes()).expect("query");
        let baseline = Poasta::with_options(costs(), PoastaSearchOptions::astar_baseline())
            .align(&graph, &query)
            .expect("ordinary A-star");
        let dfs = Poasta::with_options(
            costs(),
            PoastaSearchOptions::astar_baseline().with_exact_match_dfs(true),
        )
        .align(&graph, &query)
        .expect("A-star plus exact-match DFS");
        assert_eq!(dfs.score(), baseline.score());
        assert_eq!(
            dfs.statistics().expanded_states(),
            baseline.statistics().expanded_states()
        );
        assert!(
            dfs.statistics().discovered_states() < baseline.statistics().discovered_states(),
            "dfs={:?}, baseline={:?}",
            dfs.statistics(),
            baseline.statistics()
        );
        assert!(
            dfs.statistics().queue_entries() < baseline.statistics().queue_entries(),
            "dfs={:?}, baseline={:?}",
            dfs.statistics(),
            baseline.statistics()
        );
        assert_eq!(dfs.statistics().greedy_exact_matches(), query.len());
    }

    #[test]
    fn optimized_policy_keeps_easy_identity_on_the_lightweight_path() {
        let sequence = "ACGT".repeat(64);
        let graph = dna_graph(&sequence);
        let query = DnaSequence::from_ascii(sequence.as_bytes()).expect("query");
        let published = Poasta::new(costs())
            .align(&graph, &query)
            .expect("published identity");
        let optimized = Poasta::optimized(costs())
            .align(&graph, &query)
            .expect("optimized identity");

        assert_eq!(optimized.score(), published.score());
        assert_eq!(optimized.operations(), published.operations());
        assert_eq!(optimized.statistics().linear_heuristic_cells(), 0);
        assert_eq!(optimized.memory_estimate().heuristic_payload_bytes(), 0);
        assert!(optimized.statistics().queue_entries() <= published.statistics().queue_entries());
    }

    #[test]
    fn optimized_policy_activates_for_a_hard_query_and_reduces_search() {
        let graph = dna_graph(&"A".repeat(64));
        let query = DnaSequence::from_ascii("C".repeat(64).as_bytes()).expect("query");
        let oracle = GlobalAffinePoa::from_poasta_costs(costs())
            .align(&graph, &query)
            .expect("scalar oracle");
        let published = Poasta::new(costs())
            .align(&graph, &query)
            .expect("published hard query");
        let optimized = Poasta::optimized(costs())
            .align(&graph, &query)
            .expect("optimized hard query");

        assert_eq!(optimized.score(), oracle.score());
        assert_eq!(published.score(), oracle.score());
        assert_eq!(optimized.operations(), published.operations());
        assert_eq!(optimized.aligned_pairs(), published.aligned_pairs());
        assert_eq!(
            verify_poasta_alignment(&optimized, &graph, &query, costs())
                .expect("independent replay")
                .recomputed_score(),
            oracle.score()
        );
        assert!(optimized.statistics().linear_heuristic_cells() > 0);
        assert!(!optimized.statistics().used_relaxed_linear_certificate());
        assert!(optimized.statistics().used_exact_affine_certificate());
        assert!(optimized.statistics().exact_affine_certificate_cells() > 0);
        assert!(optimized.memory_estimate().heuristic_payload_bytes() > 0);
        assert!(optimized.statistics().affine_dominated_states() > 0);
        assert!(
            optimized.statistics().expanded_states() < published.statistics().expanded_states(),
            "optimized={:?}, published={:?}",
            optimized.statistics(),
            published.statistics()
        );
    }

    #[test]
    fn progressive_workspace_reuses_the_relaxed_linear_allocation_exactly() {
        let graph = dna_graph(&"A".repeat(64));
        let query = DnaSequence::from_ascii("C".repeat(64).as_bytes()).expect("query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let aligner = Poasta::optimized(costs());
        let mut workspace = PoastaWorkspace::default();

        let (first, _) = aligner
            .align_validated_with_report(&graph, &query, &topological, true, Some(&mut workspace))
            .expect("first reusable alignment");
        let allocation = workspace.costs.as_ptr();
        let capacity = workspace.costs.capacity();
        assert!(capacity > 0);

        let (second, _) = aligner
            .align_validated_with_report(&graph, &query, &topological, true, Some(&mut workspace))
            .expect("second reusable alignment");
        assert_eq!(workspace.costs.as_ptr(), allocation);
        assert_eq!(workspace.costs.capacity(), capacity);
        assert_eq!(workspace.row_bounds.len(), graph.node_count() + 2);
        assert_eq!(second.score(), first.score());
        assert_eq!(second.operations(), first.operations());
        assert_eq!(second.aligned_pairs(), first.aligned_pairs());
        assert_eq!(second.statistics(), first.statistics());
    }

    #[test]
    fn eager_policy_proves_zero_cost_before_building_certificate_tables() {
        let sequence = "ACGT".repeat(64);
        let graph = dna_graph(&sequence);
        let query = DnaSequence::from_ascii(sequence.as_bytes()).expect("query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let published = Poasta::new(costs())
            .align(&graph, &query)
            .expect("published identity");
        let mut workspace = PoastaWorkspace {
            relaxed_certificate_misses: RELAXED_LINEAR_MISS_THRESHOLD,
            ..PoastaWorkspace::default()
        };

        let (optimized, _) = Poasta::optimized(costs())
            .align_validated_with_report(&graph, &query, &topological, true, Some(&mut workspace))
            .expect("eager optimized identity");

        assert_eq!(optimized.score(), published.score());
        assert_eq!(optimized.operations(), published.operations());
        assert_eq!(optimized.aligned_pairs(), published.aligned_pairs());
        assert_eq!(optimized.statistics().linear_heuristic_cells(), 0);
        assert_eq!(optimized.statistics().exact_affine_certificate_cells(), 0);
        assert!(!optimized.statistics().used_relaxed_linear_certificate());
        assert!(!optimized.statistics().used_exact_affine_certificate());
        assert_eq!(optimized.statistics().queue_entries(), 0);
        assert_eq!(optimized.statistics().expanded_states(), 0);
    }

    #[test]
    fn exact_source_path_certificate_rejects_nonterminal_stored_paths() {
        for selected_positions in [0..2, 1..3] {
            let mut graph = dna_graph("TAC");
            let targets = graph.paths()[0].nodes()[selected_positions]
                .iter()
                .copied()
                .map(Some)
                .collect::<Vec<_>>();
            let query_text = targets
                .iter()
                .map(|target| {
                    graph
                        .node(target.expect("selected graph node"))
                        .expect("selected node exists")
                        .symbol()
                        .as_ascii()
                })
                .collect::<Vec<_>>();
            let record = SequenceRecord::new(
                "internal-path",
                DnaSequence::from_ascii(&query_text).expect("valid internal path"),
            )
            .expect("valid internal record");
            graph
                .append_record(&record, &targets)
                .expect("append internal stored path");
            let (topological, _) = graph.validate_with_orders().expect("valid graph");

            assert!(
                exact_source_path_certificate(&graph, record.sequence(), &topological)
                    .expect("internal path is not an invalid graph")
                    .is_none()
            );
            let index = PoastaGraphIndex::build(&graph, &topological, false).expect("graph index");
            let attempt = exact_affine_certificate(
                &graph,
                record.sequence(),
                &topological,
                &index,
                costs(),
                None,
                None,
                None,
            )
            .expect("internal path does not make the exact table invalid")
            .expect("small exact-affine table is representable");
            let certificate = attempt
                .certificate
                .expect("exact-affine table has a reachable certificate");
            let oracle = GlobalAffinePoa::from_poasta_costs(costs())
                .align(&graph, record.sequence())
                .expect("scalar oracle");
            assert_eq!(
                cost_to_score(certificate.optimum_cost).expect("certificate score"),
                oracle.score()
            );
        }
    }

    #[test]
    fn progressive_policy_bypasses_repeatedly_unsuccessful_relaxed_certificates() {
        let graph = dna_graph(&"A".repeat(64));
        let query = DnaSequence::from_ascii("C".repeat(64).as_bytes()).expect("query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let published = Poasta::new(costs())
            .align(&graph, &query)
            .expect("published hard query");
        let mut workspace = PoastaWorkspace {
            relaxed_certificate_misses: RELAXED_LINEAR_MISS_THRESHOLD,
            ..PoastaWorkspace::default()
        };

        let (optimized, _) = Poasta::optimized(costs())
            .align_validated_with_report(&graph, &query, &topological, true, Some(&mut workspace))
            .expect("adaptive exact-affine query");

        assert_eq!(optimized.score(), published.score());
        assert_eq!(optimized.operations(), published.operations());
        assert_eq!(optimized.aligned_pairs(), published.aligned_pairs());
        assert_eq!(optimized.statistics().linear_heuristic_cells(), 0);
        assert!(optimized.statistics().used_exact_affine_certificate());
        assert!(optimized.statistics().exact_affine_certificate_cells() > 0);
        assert_eq!(
            workspace.relaxed_certificate_misses,
            RELAXED_LINEAR_MISS_THRESHOLD
        );
        assert_eq!(workspace.exact_affine_bypasses_since_probe, 1);
    }

    #[test]
    fn progressive_policy_periodically_reprobes_the_relaxed_certificate() {
        let graph = dna_graph(&"A".repeat(64));
        let query = DnaSequence::from_ascii("C".repeat(64).as_bytes()).expect("query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let mut workspace = PoastaWorkspace {
            relaxed_certificate_misses: RELAXED_LINEAR_MISS_THRESHOLD,
            exact_affine_bypasses_since_probe: RELAXED_LINEAR_REPROBE_INTERVAL,
            ..PoastaWorkspace::default()
        };

        let (optimized, _) = Poasta::optimized(costs())
            .align_validated_with_report(&graph, &query, &topological, true, Some(&mut workspace))
            .expect("periodic relaxed-certificate probe");

        assert!(optimized.statistics().linear_heuristic_cells() > 0);
        assert!(!optimized.statistics().used_relaxed_linear_certificate());
        assert!(optimized.statistics().used_exact_affine_certificate());
        assert_eq!(workspace.exact_affine_bypasses_since_probe, 0);
        assert_eq!(
            workspace.relaxed_certificate_misses,
            RELAXED_LINEAR_MISS_THRESHOLD + 1
        );
    }

    #[test]
    fn exact_affine_band_closed_form_matches_exhaustive_length_filtering() {
        let mut graph = dna_graph("ACGTA");
        graph
            .append_record(
                &dna_record("short", "AGA"),
                &[Some(PoaNodeId::new(0)), None, Some(PoaNodeId::new(4))],
            )
            .expect("variable-length path");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let index = PoastaGraphIndex::build(&graph, &topological, false).expect("graph index");

        for query_len in 0..=12 {
            let prefix = exact_affine_prefix_bounds(&graph, &topological, &index, query_len)
                .expect("prefix bounds");
            for cursor in 0..index.sink() {
                let suffix = RelaxedLinearRowBounds {
                    minimum: index.min_remaining_nodes(cursor).expect("minimum suffix"),
                    maximum: index.max_remaining_nodes(cursor).expect("maximum suffix"),
                };
                for radius in 0..=20 {
                    let expected = (0..=query_len)
                        .filter(|&query_prefix| {
                            distance_to_range(query_prefix, prefix[cursor])
                                .saturating_add(distance_to_range(query_len - query_prefix, suffix))
                                <= radius
                        })
                        .fold(None, |bounds, query_prefix| {
                            Some(bounds.map_or((query_prefix, query_prefix), |(first, _)| {
                                (first, query_prefix)
                            }))
                        });
                    assert_eq!(
                        exact_affine_row_band(cursor, query_len, &index, &prefix, Some(radius),)
                            .expect("closed-form exact-affine band"),
                        expected,
                        "query_len={query_len}, cursor={cursor}, radius={radius}"
                    );
                }
            }
        }
    }

    #[test]
    fn exact_affine_distance_band_matches_all_small_interval_combinations() {
        for query_len in 0..=6 {
            for prefix_minimum in 0..=8 {
                for prefix_maximum in prefix_minimum..=8 {
                    let prefix = RelaxedLinearRowBounds {
                        minimum: prefix_minimum,
                        maximum: prefix_maximum,
                    };
                    for suffix_minimum in 0..=8 {
                        for suffix_maximum in suffix_minimum..=8 {
                            let suffix = RelaxedLinearRowBounds {
                                minimum: suffix_minimum,
                                maximum: suffix_maximum,
                            };
                            for radius in 0..=12 {
                                let expected =
                                    (0..=query_len)
                                        .filter(|&query_prefix| {
                                            distance_to_range(query_prefix, prefix).saturating_add(
                                                distance_to_range(query_len - query_prefix, suffix),
                                            ) <= radius
                                        })
                                        .fold(None, |bounds, query_prefix| {
                                            Some(bounds.map_or(
                                                (query_prefix, query_prefix),
                                                |(first, _)| (first, query_prefix),
                                            ))
                                        });
                                assert_eq!(
                                    exact_affine_distance_band(prefix, suffix, query_len, radius,),
                                    expected,
                                    "query_len={query_len}, prefix={prefix:?}, suffix={suffix:?}, radius={radius}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn exact_affine_distance_band_handles_usize_boundary_intervals() {
        let maximum = usize::MAX;
        assert_eq!(
            exact_affine_distance_band(
                RelaxedLinearRowBounds {
                    minimum: maximum - 1,
                    maximum,
                },
                RelaxedLinearRowBounds {
                    minimum: maximum - 1,
                    maximum,
                },
                maximum,
                maximum,
            ),
            Some((0, maximum))
        );
        assert_eq!(
            exact_affine_distance_band(
                RelaxedLinearRowBounds {
                    minimum: maximum,
                    maximum,
                },
                RelaxedLinearRowBounds {
                    minimum: maximum,
                    maximum,
                },
                0,
                maximum - 1,
            ),
            None
        );
    }

    #[test]
    fn exact_affine_buffer_keeps_its_initialized_high_water_length() {
        let original = vec![7_u16; 128];
        let allocation = original.as_ptr();
        let prepared = prepare_exact_affine_buffer(original, 16, "test exact-affine buffer")
            .expect("reused exact-affine buffer");

        assert_eq!(prepared.as_ptr(), allocation);
        assert_eq!(prepared.len(), 128);
        assert!(prepared[..16].iter().all(|&value| value == 7));
    }

    #[test]
    fn diagonal_upper_bound_tightens_exact_affine_without_changing_traceback() {
        let mut graph = dna_graph("AAAAAGCCCCC");
        graph
            .append_record(
                &dna_record("alternate", "TTTTTGGGGGG"),
                &[
                    None,
                    None,
                    None,
                    None,
                    None,
                    Some(PoaNodeId::new(5)),
                    None,
                    None,
                    None,
                    None,
                    None,
                ],
            )
            .expect("recombining alternate path");
        let query = DnaSequence::from_ascii(b"AAAAAGGGGGG").expect("recombinant query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let index = PoastaGraphIndex::build(&graph, &topological, false).expect("graph index");
        let diagonal = exact_diagonal_certificate(
            &graph,
            &query,
            &topological,
            &index,
            nonnegative_score(costs().mismatch()),
            None,
            None,
        )
        .expect("diagonal certificate attempt")
        .expect("recombinant diagonal path");
        assert_eq!(diagonal.certificate.optimum_cost, 0);

        let unbounded = exact_affine_certificate(
            &graph,
            &query,
            &topological,
            &index,
            costs(),
            None,
            None,
            None,
        )
        .expect("unbounded exact-affine attempt")
        .expect("compact exact-affine range");
        let bounded = exact_affine_certificate(
            &graph,
            &query,
            &topological,
            &index,
            costs(),
            Some(diagonal.certificate.optimum_cost),
            None,
            None,
        )
        .expect("bounded exact-affine attempt")
        .expect("compact exact-affine range");

        assert!(bounded.cells < unbounded.cells);
        assert_eq!(bounded.certificate, unbounded.certificate);
    }

    #[test]
    fn exact_prefix_cache_preserves_both_certificate_attempts() {
        let mut graph = dna_graph("AAAAAGCCCCC");
        graph
            .append_record(
                &dna_record("alternate", "TTTTTGGGGGG"),
                &[
                    None,
                    None,
                    None,
                    None,
                    None,
                    Some(PoaNodeId::new(5)),
                    None,
                    None,
                    None,
                    None,
                    None,
                ],
            )
            .expect("recombining alternate path");
        let query = DnaSequence::from_ascii(b"AAAAAGGGGGG").expect("recombinant query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let index = PoastaGraphIndex::build(&graph, &topological, false).expect("graph index");
        let prefix_bounds = exact_affine_prefix_bounds(&graph, &topological, &index, query.len())
            .expect("exact prefix bounds");
        let mismatch = nonnegative_score(costs().mismatch());

        let local_diagonal =
            exact_diagonal_certificate(&graph, &query, &topological, &index, mismatch, None, None)
                .expect("local diagonal attempt")
                .expect("local diagonal certificate");
        let cached_diagonal = exact_diagonal_certificate(
            &graph,
            &query,
            &topological,
            &index,
            mismatch,
            Some(&prefix_bounds),
            None,
        )
        .expect("cached diagonal attempt")
        .expect("cached diagonal certificate");
        assert_eq!(
            cached_diagonal.certificate.path,
            local_diagonal.certificate.path
        );
        assert_eq!(
            cached_diagonal.certificate.optimum_cost,
            local_diagonal.certificate.optimum_cost
        );
        assert_eq!(cached_diagonal.cells, local_diagonal.cells);
        assert_eq!(cached_diagonal.payload_bytes, local_diagonal.payload_bytes);

        let local_affine = exact_affine_certificate(
            &graph,
            &query,
            &topological,
            &index,
            costs(),
            Some(local_diagonal.certificate.optimum_cost),
            None,
            None,
        )
        .expect("local exact-affine attempt")
        .expect("local exact-affine certificate");
        let cached_affine = exact_affine_certificate(
            &graph,
            &query,
            &topological,
            &index,
            costs(),
            Some(cached_diagonal.certificate.optimum_cost),
            Some(&prefix_bounds),
            None,
        )
        .expect("cached exact-affine attempt")
        .expect("cached exact-affine certificate");
        assert_eq!(cached_affine.certificate, local_affine.certificate);
        assert_eq!(cached_affine.cells, local_affine.cells);
        assert_eq!(cached_affine.payload_bytes, local_affine.payload_bytes);
    }

    #[test]
    fn declined_diagonal_still_accounts_for_cached_prefix_bounds() {
        let graph = dna_graph("AA");
        let query = DnaSequence::from_ascii(b"CC").expect("query");
        let wide_costs = PoastaAffineCosts::new(40_000, 6, 2).expect("valid wide costs");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let index = PoastaGraphIndex::build(&graph, &topological, false).expect("graph index");
        let total_possible_states =
            possible_state_count(graph.node_count(), query.len()).expect("state count");
        let mut search = Search::new(
            &graph,
            &query,
            wide_costs,
            PoastaSearchOptions::optimized(),
            &index,
            &topological,
            true,
            total_possible_states,
            None,
        )
        .expect("search");

        search
            .prepare_exact_diagonal_certificate::<false>(None)
            .expect("checked diagonal range");

        assert!(search.exact_diagonal_attempt.is_none());
        let prefix_bounds = search
            .exact_prefix_bounds
            .as_ref()
            .expect("cached exact prefix bounds");
        assert_eq!(
            search.certificate_heuristic_payload_bytes,
            prefix_bounds.capacity() * size_of::<RelaxedLinearRowBounds>()
        );
    }

    #[test]
    fn zero_diagonal_certificate_matches_relaxed_forward_tie_breaking() {
        let mut graph = dna_graph("ATG");
        graph
            .append_record(
                &dna_record("low-root-high-sink", "ACG"),
                &[Some(PoaNodeId::new(0)), None, None],
            )
            .expect("first exact alternate path");
        graph
            .append_record(
                &dna_record("high-root-low-sink", "ACG"),
                &[None, None, Some(PoaNodeId::new(2))],
            )
            .expect("second exact alternate path");
        let query = DnaSequence::from_ascii(b"ACG").expect("query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let index = PoastaGraphIndex::build(&graph, &topological, false).expect("graph index");
        let heuristic = RelaxedLinearHeuristic::build::<_, false>(
            &graph,
            &query,
            &topological,
            &index,
            costs(),
            None,
            None,
            None,
        )
        .expect("relaxed heuristic range")
        .expect("relaxed heuristic");
        let relaxed = heuristic
            .diagonal_certificate(
                &graph,
                &query,
                &index,
                nonnegative_score(costs().mismatch()),
            )
            .expect("relaxed certificate attempt")
            .expect("relaxed zero certificate");
        let exact = exact_diagonal_certificate(
            &graph,
            &query,
            &topological,
            &index,
            nonnegative_score(costs().mismatch()),
            None,
            None,
        )
        .expect("exact diagonal attempt")
        .expect("exact zero certificate")
        .certificate;

        assert_eq!(exact.optimum_cost, 0);
        assert_eq!(exact.path, relaxed.path);
        assert_eq!(
            exact.path.iter().map(|key| key.cursor).collect::<Vec<_>>(),
            vec![1, 4, 5]
        );
    }

    #[test]
    fn diagonal_certificate_handles_zero_mismatch_and_recycles_empty_attempt_buffers() {
        let graph = dna_graph("ACGT");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let index = PoastaGraphIndex::build(&graph, &topological, false).expect("graph index");
        let reversed = DnaSequence::from_ascii(b"TGCA").expect("reversed query");
        let zero_mismatch =
            exact_diagonal_certificate(&graph, &reversed, &topological, &index, 0, None, None)
                .expect("zero-mismatch diagonal attempt")
                .expect("zero-mismatch diagonal certificate");
        assert_eq!(zero_mismatch.certificate.optimum_cost, 0);
        assert_eq!(zero_mismatch.certificate.path.len(), reversed.len());

        let short = DnaSequence::from_ascii(b"ACG").expect("short query");
        let mut workspace = PoastaWorkspace::default();
        assert!(
            exact_diagonal_certificate(
                &graph,
                &short,
                &topological,
                &index,
                nonnegative_score(costs().mismatch()),
                None,
                Some(&mut workspace),
            )
            .expect("infeasible diagonal attempt")
            .is_none()
        );
        assert!(!workspace.exact_match_costs.is_empty());
        assert!(!workspace.exact_row_bands.is_empty());
    }

    #[test]
    fn exact_affine_certificate_preserves_an_ambiguous_published_gapped_traceback() {
        let graph = dna_graph(&"A".repeat(64));
        let query = DnaSequence::from_ascii("C".repeat(70).as_bytes()).expect("query");
        let published = Poasta::new(costs())
            .align(&graph, &query)
            .expect("published gapped query");
        let optimized = Poasta::optimized(costs())
            .align(&graph, &query)
            .expect("optimized gapped query");

        assert_eq!(optimized.score(), published.score());
        assert_eq!(optimized.operations(), published.operations());
        assert_eq!(optimized.aligned_pairs(), published.aligned_pairs());
        assert!(
            optimized
                .operations()
                .iter()
                .any(|operation| matches!(operation, PoaOperation::QueryOnly { .. }))
        );
        assert!(optimized.statistics().used_exact_affine_certificate());
        assert!(optimized.statistics().exact_affine_certificate_cells() > 0);
        assert!(optimized.memory_estimate().heuristic_payload_bytes() > 0);
        assert_eq!(
            verify_poasta_alignment(&optimized, &graph, &query, costs())
                .expect("independent gapped replay")
                .recomputed_score(),
            published.score()
        );
    }

    #[test]
    fn exact_affine_certificate_preserves_a_unique_hard_gapped_traceback() {
        let source = "CGTCCGATGCTACGGCATTCGATCGTACCTGACGTTAGC";
        let graph = dna_graph(source);
        let mut query_ascii = source.as_bytes().to_vec();
        query_ascii.splice(20..20, b"AAAAAA".iter().copied());
        let query = DnaSequence::from_ascii(query_ascii).expect("query");
        let published = Poasta::new(costs())
            .align(&graph, &query)
            .expect("published unique gapped query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let (optimized, _) = Poasta::optimized(costs())
            .align_validated_with_report(&graph, &query, &topological, true, None)
            .expect("eager optimized unique gapped query");

        assert_eq!(optimized.score(), published.score());
        assert_eq!(optimized.operations(), published.operations());
        assert_eq!(optimized.aligned_pairs(), published.aligned_pairs());
        assert!(
            optimized
                .operations()
                .iter()
                .any(|operation| matches!(operation, PoaOperation::QueryOnly { .. }))
        );
        assert!(optimized.statistics().used_exact_affine_certificate());
        assert!(optimized.statistics().exact_affine_certificate_cells() > 0);
        assert_eq!(
            verify_poasta_alignment(&optimized, &graph, &query, costs())
                .expect("independent unique gapped replay")
                .recomputed_score(),
            published.score()
        );
    }

    #[test]
    fn compact_heuristics_decline_unrepresentable_u16_ranges_without_changing_results() {
        let wide_costs = PoastaAffineCosts::new(40_000, 30_001, 30_000).expect("valid wide costs");
        let graph = dna_graph(&"A".repeat(32));
        let query = DnaSequence::from_ascii("C".repeat(32).as_bytes()).expect("query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let index =
            PoastaGraphIndex::build(&graph, &topological, false).expect("compact graph index");

        assert!(
            RelaxedLinearHeuristic::build::<_, false>(
                &graph,
                &query,
                &topological,
                &index,
                wide_costs,
                None,
                None,
                None,
            )
            .expect("checked range gate")
            .is_none()
        );
        assert!(
            exact_affine_certificate(
                &graph,
                &query,
                &topological,
                &index,
                wide_costs,
                None,
                None,
                None,
            )
            .expect("checked exact range gate")
            .is_none()
        );

        let published = Poasta::new(wide_costs)
            .align(&graph, &query)
            .expect("published wide-cost query");
        let optimized = Poasta::optimized(wide_costs)
            .align(&graph, &query)
            .expect("optimized wide-cost query");
        assert_eq!(optimized.score(), published.score());
        assert_eq!(optimized.operations(), published.operations());
        assert_eq!(optimized.aligned_pairs(), published.aligned_pairs());
        assert_eq!(optimized.statistics().linear_heuristic_cells(), 0);
        assert_eq!(optimized.statistics().exact_affine_certificate_cells(), 0);
    }

    #[test]
    fn relaxed_linear_potential_is_admissible_from_every_affine_state() {
        let mut graph = dna_graph("ACGTA");
        graph
            .append_record(
                &dna_record("alternate", "ATCTA"),
                &[
                    Some(PoaNodeId::new(0)),
                    None,
                    None,
                    Some(PoaNodeId::new(3)),
                    Some(PoaNodeId::new(4)),
                ],
            )
            .expect("branched graph");
        let query = DnaSequence::from_ascii(b"AGGTT").expect("query");
        let (topological, _) = graph.validate_with_orders().expect("valid graph");
        let index = PoastaGraphIndex::build(&graph, &topological, true).expect("graph index");
        let heuristic = RelaxedLinearHeuristic::build::<_, false>(
            &graph,
            &query,
            &topological,
            &index,
            costs(),
            None,
            None,
            None,
        )
        .expect("heuristic build")
        .expect("compact cost range");
        assert!(heuristic.cell_count() > 0);

        let mut memo = HashMap::new();
        for cursor in index.source()..index.sink() {
            let remaining_bounds = RelaxedLinearRowBounds {
                minimum: index
                    .min_remaining_nodes(cursor)
                    .expect("minimum terminal distance"),
                maximum: index
                    .max_remaining_nodes(cursor)
                    .expect("maximum terminal distance"),
            };
            for query_prefix in 0..=query.len() {
                let length_distance =
                    distance_to_range(query.len() - query_prefix, remaining_bounds);
                let length_bound = Cost::try_from(length_distance)
                    .expect("test length fits the compact cost")
                    * nonnegative_score(costs().gap_extension());
                for state in [
                    AffineState::Match,
                    AffineState::Deletion,
                    AffineState::Insertion,
                ] {
                    let key = StateKey::new(cursor, query_prefix, state);
                    let lower_bound = heuristic.lower_bound(key);
                    let exact = exact_remaining_cost(&graph, &query, &index, key, &mut memo);
                    assert!(
                        lower_bound >= length_bound,
                        "cursor={cursor}, query_prefix={query_prefix}, state={state:?}, lower_bound={lower_bound}, length_bound={length_bound}"
                    );
                    assert!(
                        lower_bound <= exact,
                        "cursor={cursor}, query_prefix={query_prefix}, state={state:?}, lower_bound={lower_bound}, exact={exact}"
                    );
                }
            }
        }
    }

    #[test]
    fn superbubble_pruning_removes_a_dominated_divergent_allele() {
        let backbone = "AACCCCCCGTT";
        let alternate = "AATTTTTTGTT";
        let mut graph = dna_graph(backbone);
        let targets = (0..backbone.len())
            .map(|index| Some(PoaNodeId::new(index)))
            .collect::<Vec<_>>();
        graph
            .append_record(&dna_record("alternate", alternate), &targets)
            .expect("alternate branch");
        let query = DnaSequence::from_ascii(b"AACCCCCCGTA").expect("query");

        let without_pruning = Poasta::with_options(
            costs(),
            PoastaSearchOptions::astar_baseline().with_exact_match_dfs(true),
        )
        .align(&graph, &query)
        .expect("unpruned alignment");
        let published = Poasta::new(costs())
            .align(&graph, &query)
            .expect("published alignment");

        assert_eq!(published.score(), without_pruning.score());
        assert!(
            published.statistics().bubble_pruned_states() > 0,
            "statistics={:?}",
            published.statistics()
        );
        assert!(
            published.statistics().expanded_states()
                < without_pruning.statistics().expanded_states(),
            "published={:?}, unpruned={:?}",
            published.statistics(),
            without_pruning.statistics()
        );
    }

    #[test]
    fn implicit_gap_bounds_cover_exact_insertion_deletion_and_boundaries() {
        let mut reached = ReachedExit::default();
        reached.insert(2, 1).expect("insert first reached position");
        reached
            .insert(5, 3)
            .expect("insert second reached position");

        assert!(position_is_dominated(1, 2, &mut reached, 8, 8, costs()).unwrap());
        assert!(!position_is_dominated(0, 2, &mut reached, 8, 8, costs()).unwrap());

        // From j=2 to x=4, the implicit insertion costs 6 + 2*2.
        assert!(position_is_dominated(11, 4, &mut reached, 0, 8, costs()).unwrap());
        assert!(!position_is_dominated(10, 4, &mut reached, 0, 8, costs()).unwrap());

        let mut deletion_only = ReachedExit::default();
        deletion_only
            .insert(5, 1)
            .expect("insert deletion evidence");
        assert!(position_is_dominated(13, 2, &mut deletion_only, 3, 8, costs()).unwrap());
        assert!(!position_is_dominated(13, 2, &mut deletion_only, 2, 8, costs()).unwrap());

        let mut empty = ReachedExit::default();
        assert!(!position_is_dominated(Cost::MAX, 0, &mut empty, usize::MAX, 0, costs()).unwrap());
    }

    #[test]
    fn bubble_subset_requires_every_endpoint_and_reached_neighbor_to_be_dominated() {
        let mut reached = ReachedExit::default();
        reached.insert(3, 0).unwrap();
        reached.insert(6, 4).unwrap();
        assert!(bubble_dominates(12, 2, 7, &mut reached, 8, 7, costs()).unwrap());
        assert!(!bubble_dominates(11, 2, 7, &mut reached, 8, 7, costs()).unwrap());
        let mut empty = ReachedExit::default();
        assert!(!bubble_dominates(Cost::MAX, 2, 7, &mut empty, 8, 7, costs()).unwrap());
    }

    #[test]
    fn affine_gap_guard_retains_a_candidate_that_can_extend_an_open_gap() {
        let mut reached = ReachedExit::default();
        reached.insert(4, 10).unwrap();
        assert!(open_gap_can_improve(15, &reached, costs()));
        assert!(!open_gap_can_improve(16, &reached, costs()));
        assert!(!open_gap_can_improve(100, &ReachedExit::default(), costs()));
    }

    #[test]
    fn implicit_gap_overflow_is_a_typed_error() {
        assert_eq!(
            implicit_gap_bound(Cost::MAX, 1, costs()),
            Err(PoaError::PoastaCostOverflow {
                context: "an implicit bubble-gap upper bound"
            })
        );
    }

    #[test]
    fn indexed_implicit_bounds_match_independent_exhaustive_scans() {
        for query_len in 0_usize..=7 {
            for mask in 1_usize..(1_usize << (query_len + 1)) {
                let mut reached = ReachedExit::default();
                for position in 0..=query_len {
                    if mask & (1_usize << position) != 0 {
                        reached
                            .insert(position, ((position * 7 + mask) % 13) as Cost)
                            .unwrap();
                    }
                }
                for query_prefix in 0..=query_len {
                    for maximum_deletion in 0..=query_len {
                        let expected = reached.positions.iter().fold(None, |minimum, position| {
                            let bound = if position.query_prefix == query_prefix {
                                Some(position.cost)
                            } else if position.query_prefix < query_prefix {
                                implicit_gap_bound(
                                    position.cost,
                                    query_prefix - position.query_prefix,
                                    costs(),
                                )
                                .unwrap()
                            } else if position.query_prefix - query_prefix <= maximum_deletion {
                                implicit_gap_bound(
                                    position.cost,
                                    position.query_prefix - query_prefix,
                                    costs(),
                                )
                                .unwrap()
                            } else {
                                None
                            };
                            optional_min(minimum, bound)
                        });
                        let actual = reached
                            .upper_bound(query_prefix, maximum_deletion, query_len, costs())
                            .unwrap();
                        assert_eq!(
                            actual, expected,
                            "query_len={query_len}, mask={mask:#x}, query_prefix={query_prefix}, maximum_deletion={maximum_deletion}"
                        );
                    }
                }
            }
        }
    }
}
