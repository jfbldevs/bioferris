//! Fluent entry points for POASTA graph alignment.

use align_core::{Alphabet, Score, Sequence};
use align_msa::{MultipleAlignment, SequenceRecord};
use align_poa::{
    LEE_2002_POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_ITERATIVE_V1,
    POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_V1,
};
use align_poa::{PoaGraph, Poasta, PoastaAffineCosts, PoastaAlignment, PoastaPoa, PoastaPoaResult};

use crate::report::{MultipleRun, MultipleRunInfo, PoastaQueryRun, ReportField};
use crate::result::MultipleAlignmentView;
use crate::{Error, Result};

const DEFAULT_MISMATCH_COST: Score = 4;
const DEFAULT_GAP_OPEN_COST: Score = 6;
const DEFAULT_GAP_EXTENSION_COST: Score = 2;

/// Discoverable constructor for exact sequence-to-graph and progressive POASTA.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Graph;

impl Graph {
    /// Selects the complete published POASTA acceleration profile.
    #[must_use]
    pub const fn poasta() -> PoastaAligner {
        PoastaAligner::new()
    }
}

/// Reusable POASTA affine-cost configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PoastaAligner {
    mismatch: Score,
    gap_open: Score,
    gap_extension: Score,
}

impl PoastaAligner {
    const fn new() -> Self {
        Self {
            mismatch: DEFAULT_MISMATCH_COST,
            gap_open: DEFAULT_GAP_OPEN_COST,
            gap_extension: DEFAULT_GAP_EXTENSION_COST,
        }
    }

    /// Replaces positive POASTA costs.
    #[must_use]
    pub const fn costs(mut self, mismatch: Score, gap_open: Score, gap_extension: Score) -> Self {
        self.mismatch = mismatch;
        self.gap_open = gap_open;
        self.gap_extension = gap_extension;
        self
    }

    fn costs_checked(self) -> Result<PoastaAffineCosts> {
        PoastaAffineCosts::new(self.mismatch, self.gap_open, self.gap_extension)
            .map_err(Error::from)
    }

    /// Aligns one complete query to a complete path through an existing DAG.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Poasta`] for invalid costs or graph, checked arithmetic,
    /// allocation, search, traceback, or replay failure.
    pub fn align_query<A: Alphabet>(
        &self,
        graph: &PoaGraph<A>,
        query: &Sequence<A>,
    ) -> Result<PoastaAlignment> {
        Poasta::new(self.costs_checked()?)
            .align(graph, query)
            .map_err(Error::from)
    }

    /// Aligns one query to a DAG and retains the graph, query, costs, traceback,
    /// and exact-search evidence for presentation.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align_query`].
    pub fn run_query<'a, A: Alphabet>(
        &self,
        graph: &'a PoaGraph<A>,
        query: &'a Sequence<A>,
    ) -> Result<PoastaQueryRun<'a, A>> {
        let costs = self.costs_checked()?;
        let result = Poasta::new(costs).align(graph, query)?;
        Ok(PoastaQueryRun::new(
            graph,
            query,
            result,
            costs,
            POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_V1,
        ))
    }

    /// Builds a partial-order graph and rectangular MSA in caller record order.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Poasta`] for invalid costs or records, checked arithmetic,
    /// allocation, search, graph incorporation, or source-recovery failure.
    pub fn align<A: Alphabet>(&self, records: &[SequenceRecord<A>]) -> Result<PoastaPoaResult<A>> {
        PoastaPoa::new(self.costs_checked()?)
            .align(records)
            .map_err(Error::from)
    }

    /// Builds a partial-order graph and retains graph-growth and memory evidence
    /// alongside the rectangular alignment report.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn run<A: Alphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<MultipleRun<A, PoastaPoaResult<A>>> {
        let result = self.align(records)?;
        let graph = result.graph();
        let memory = result.memory_estimate();
        let aligned_steps = result
            .steps()
            .iter()
            .filter(|step| !step.alignment().is_seed())
            .count();
        let info = MultipleRunInfo::new(
            "POASTA",
            LEE_2002_POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_ITERATIVE_V1,
            Some("progressive partial-order alignment"),
        )
        .with_summary(vec![
            ReportField::new("Graph nodes", graph.node_count().to_string()),
            ReportField::new("Graph edges", graph.edge_count().to_string()),
            ReportField::new("Graph classes", graph.class_count().to_string()),
            ReportField::new("Graph paths", graph.path_count().to_string()),
        ])
        .with_details(vec![
            ReportField::new("Construction steps", result.steps().len().to_string()),
            ReportField::new("Aligned steps", aligned_steps.to_string()),
            ReportField::new(
                "Peak alignment bytes",
                memory.peak_alignment_bytes().to_string(),
            ),
            ReportField::new(
                "Graph storage bytes",
                memory.graph_storage_bytes().to_string(),
            ),
            ReportField::new(
                "Rectangular bytes",
                memory.rectangular_alignment_bytes().to_string(),
            ),
            ReportField::new(
                "Estimated peak bytes",
                memory.estimated_peak_bytes().to_string(),
            ),
        ]);
        Ok(MultipleRun::new(result, info))
    }
}

impl Default for PoastaAligner {
    fn default() -> Self {
        Self::new()
    }
}

impl<A: Alphabet> MultipleAlignmentView for PoastaPoaResult<A> {
    type Alphabet = A;

    fn alignment(&self) -> &MultipleAlignment<A> {
        self.alignment()
    }
}
