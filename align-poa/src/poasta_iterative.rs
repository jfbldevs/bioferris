use core::mem::size_of;

use align_core::{Alphabet, PhaseTimings, PhaseTimingsBuilder, Profiled};
use align_msa::{MultipleAlignment, SequenceRecord};

use crate::iterative_support::{modeled_alignment_bytes, verify_source_recovery};
use crate::poasta::PoastaWorkspace;
use crate::{
    AffinePoaVerificationReport, PoaError, PoaGraph, PoaNodeId, PoaOperation, PoaPathId, Poasta,
    PoastaAffineCosts, PoastaAlignment, PoastaBackend, PoastaSearchOptions,
    PoastaSearchWorkCounters,
};

/// Immutable identifier of the complete iterative accelerated POASTA profile.
pub const LEE_2002_POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_ITERATIVE_V1: &str =
    "LEE-2002-POASTA-2025-global-affine-ASTAR-DFS-superbubble-iterative-v1";

/// Whether one iterative step seeded an empty graph or ran accelerated POASTA.
#[derive(Clone, Debug, Eq, PartialEq)]
#[allow(clippy::large_enum_variant)]
pub enum PoastaPoaStepAlignment {
    /// The record was appended directly because the graph had no biological node.
    Seed,
    /// The record was aligned, replayed, and then incorporated into the graph.
    Aligned {
        /// Exact accelerated sequence-to-current-graph result.
        alignment: PoastaAlignment,
        /// Independent affine path-and-score replay evidence.
        verification: AffinePoaVerificationReport,
    },
}

impl PoastaPoaStepAlignment {
    /// Returns whether this record directly seeded an empty graph.
    #[must_use]
    pub const fn is_seed(&self) -> bool {
        matches!(self, Self::Seed)
    }

    /// Returns the POASTA alignment when this was not a seed.
    #[must_use]
    pub const fn alignment(&self) -> Option<&PoastaAlignment> {
        match self {
            Self::Seed => None,
            Self::Aligned { alignment, .. } => Some(alignment),
        }
    }

    /// Returns independent replay evidence when this was not a seed.
    #[must_use]
    pub const fn verification(&self) -> Option<AffinePoaVerificationReport> {
        match self {
            Self::Seed => None,
            Self::Aligned { verification, .. } => Some(*verification),
        }
    }
}

/// Auditable evidence for one record incorporated by accelerated POASTA.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PoastaPoaStep {
    record_index: usize,
    path: PoaPathId,
    alignment: PoastaPoaStepAlignment,
    nodes_before: usize,
    nodes_after: usize,
    classes_before: usize,
    classes_after: usize,
    edges_before: usize,
    edges_after: usize,
}

impl PoastaPoaStep {
    /// Returns the stable zero-based input record index.
    #[must_use]
    pub const fn record_index(&self) -> usize {
        self.record_index
    }

    /// Returns the graph path created for this source record.
    #[must_use]
    pub const fn path(&self) -> PoaPathId {
        self.path
    }

    /// Returns seed or accelerated-alignment evidence for this step.
    #[must_use]
    pub const fn alignment(&self) -> &PoastaPoaStepAlignment {
        &self.alignment
    }

    /// Returns the node count immediately before incorporation.
    #[must_use]
    pub const fn nodes_before(&self) -> usize {
        self.nodes_before
    }

    /// Returns the node count immediately after incorporation.
    #[must_use]
    pub const fn nodes_after(&self) -> usize {
        self.nodes_after
    }

    /// Returns the alignment-class count immediately before incorporation.
    #[must_use]
    pub const fn classes_before(&self) -> usize {
        self.classes_before
    }

    /// Returns the alignment-class count immediately after incorporation.
    #[must_use]
    pub const fn classes_after(&self) -> usize {
        self.classes_after
    }

    /// Returns the edge count immediately before incorporation.
    #[must_use]
    pub const fn edges_before(&self) -> usize {
        self.edges_before
    }

    /// Returns the edge count immediately after incorporation.
    #[must_use]
    pub const fn edges_after(&self) -> usize {
        self.edges_after
    }
}

/// Checked memory evidence for one complete iterative POASTA call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(clippy::struct_field_names)]
pub struct PoastaPoaMemoryEstimate {
    peak_alignment_bytes: usize,
    peak_integration_target_bytes: usize,
    graph_storage_bytes: usize,
    rectangular_alignment_bytes: usize,
    step_evidence_bytes: usize,
    estimated_peak_bytes: usize,
}

impl PoastaPoaMemoryEstimate {
    /// Returns the largest sparse-search modeled payload requirement.
    #[must_use]
    pub const fn peak_alignment_bytes(self) -> usize {
        self.peak_alignment_bytes
    }

    /// Returns the largest query-to-node target buffer.
    #[must_use]
    pub const fn peak_integration_target_bytes(self) -> usize {
        self.peak_integration_target_bytes
    }

    /// Returns modeled bytes retained by the final graph.
    #[must_use]
    pub const fn graph_storage_bytes(self) -> usize {
        self.graph_storage_bytes
    }

    /// Returns modeled bytes retained by the rectangular serialization.
    #[must_use]
    pub const fn rectangular_alignment_bytes(self) -> usize {
        self.rectangular_alignment_bytes
    }

    /// Returns modeled bytes retained by all per-record evidence.
    #[must_use]
    pub const fn step_evidence_bytes(self) -> usize {
        self.step_evidence_bytes
    }

    /// Returns the largest modeled simultaneous payload total.
    #[must_use]
    pub const fn estimated_peak_bytes(self) -> usize {
        self.estimated_peak_bytes
    }
}

/// Complete deterministic result of iterative accelerated POASTA.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PoastaPoaResult<A: Alphabet> {
    graph: PoaGraph<A>,
    alignment: MultipleAlignment<A>,
    steps: Box<[PoastaPoaStep]>,
    memory: PoastaPoaMemoryEstimate,
}

/// One iterative POASTA result with coarse and nested search diagnostics.
#[derive(Clone, Debug)]
pub struct PoastaPoaDetailedProfile<A: Alphabet> {
    result: PoastaPoaResult<A>,
    timings: PhaseTimings,
    search_timings: PhaseTimings,
    search_work_counters: PoastaSearchWorkCounters,
}

impl<A: Alphabet> PoastaPoaDetailedProfile<A> {
    /// Returns the unchanged scientific POA result.
    #[must_use]
    pub const fn result(&self) -> &PoastaPoaResult<A> {
        &self.result
    }

    /// Returns the complete-call coarse phase timings.
    #[must_use]
    pub const fn timings(&self) -> &PhaseTimings {
        &self.timings
    }

    /// Returns nested timings whose work occurs inside the coarse search phase.
    #[must_use]
    pub const fn search_timings(&self) -> &PhaseTimings {
        &self.search_timings
    }

    /// Returns deterministic work counters collected inside detailed search.
    #[must_use]
    pub const fn search_work_counters(&self) -> PoastaSearchWorkCounters {
        self.search_work_counters
    }

    /// Separates the scientific result and both timing views.
    #[must_use]
    pub fn into_parts(self) -> (PoastaPoaResult<A>, PhaseTimings, PhaseTimings) {
        (self.result, self.timings, self.search_timings)
    }
}

impl<A: Alphabet> PoastaPoaResult<A> {
    /// Returns the authoritative partial-order graph.
    #[must_use]
    pub const fn graph(&self) -> &PoaGraph<A> {
        &self.graph
    }

    /// Returns the deterministic alignment-class serialization.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        &self.alignment
    }

    /// Returns one construction step per input record.
    #[must_use]
    pub const fn steps(&self) -> &[PoastaPoaStep] {
        &self.steps
    }

    /// Returns checked modeled memory evidence.
    #[must_use]
    pub const fn memory_estimate(&self) -> PoastaPoaMemoryEstimate {
        self.memory
    }

    /// Separates the authoritative graph and rectangular serialization.
    #[must_use]
    pub fn into_graph_and_alignment(self) -> (PoaGraph<A>, MultipleAlignment<A>) {
        (self.graph, self.alignment)
    }
}

/// Iterative POA construction driven by the exact accelerated POASTA kernel.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PoastaPoa {
    kernel: Poasta,
}

impl PoastaPoa {
    /// Creates the complete published accelerated search plus Lee incorporation.
    #[must_use]
    pub const fn new(costs: PoastaAffineCosts) -> Self {
        Self {
            kernel: Poasta::new(costs),
        }
    }

    /// Creates the exact pipeline using the validated optimized-execution policy.
    #[must_use]
    pub const fn optimized(costs: PoastaAffineCosts) -> Self {
        Self {
            kernel: Poasta::optimized(costs),
        }
    }

    /// Creates an exact feature-isolation pipeline for validation and benchmarks.
    #[must_use]
    pub const fn with_options(costs: PoastaAffineCosts, options: PoastaSearchOptions) -> Self {
        Self {
            kernel: Poasta::with_options(costs, options),
        }
    }

    /// Returns the immutable complete iterative profile identifier when applicable.
    #[must_use]
    pub const fn profile_identifier(self) -> &'static str {
        if matches!(self.kernel.options().backend(), PoastaBackend::Sparse)
            && self.kernel.options().exact_match_dfs()
            && self.kernel.options().superbubble_pruning()
        {
            LEE_2002_POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_ITERATIVE_V1
        } else {
            self.kernel.profile_identifier()
        }
    }

    /// Returns the configured sequence-to-graph kernel.
    #[must_use]
    pub const fn kernel(self) -> Poasta {
        self.kernel
    }

    /// Builds a PO-MSA by adding records in caller order with exact POASTA search.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError`] for an empty collection, duplicate identifier,
    /// invalid graph, arithmetic or allocation failure, invalid traceback, or
    /// failed source recovery.
    #[allow(clippy::too_many_lines)]
    pub fn align<A: Alphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<PoastaPoaResult<A>, PoaError> {
        self.align_with_phase_timings::<A, false>(records, None, None, None)
    }

    /// Builds a PO-MSA and reports aggregate monotonic wall time by phase.
    ///
    /// This diagnostic entry point returns the same graph, alignment, and audit
    /// evidence as [`Self::align`]. Sequence-to-graph phases are aggregated over
    /// every non-seed record, and clock reads occur only on this opt-in path.
    ///
    /// # Errors
    ///
    /// Returns the same [`PoaError`] values as [`Self::align`].
    pub fn align_profiled<A: Alphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<Profiled<PoastaPoaResult<A>>, PoaError> {
        let mut timings = PhaseTimingsBuilder::new();
        let result =
            self.align_with_phase_timings::<A, false>(records, Some(&mut timings), None, None)?;
        Ok(Profiled::new(result, timings.finish()))
    }

    /// Builds a PO-MSA and adds nested timings for work inside sparse search.
    ///
    /// The ordinary [`Self::align`] and coarse [`Self::align_profiled`] paths do
    /// not read the inner clocks. Nested durations aggregate every non-seed
    /// record and are diagnostic evidence, not native-comparator timings.
    ///
    /// # Errors
    ///
    /// Returns the same [`PoaError`] values as [`Self::align`].
    pub fn align_detailed_profiled<A: Alphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<PoastaPoaDetailedProfile<A>, PoaError> {
        let mut timings = PhaseTimingsBuilder::new();
        let mut search_timings = PhaseTimingsBuilder::new();
        let mut search_work_counters = PoastaSearchWorkCounters::default();
        let result = self.align_with_phase_timings::<A, true>(
            records,
            Some(&mut timings),
            Some(&mut search_timings),
            Some(&mut search_work_counters),
        )?;
        Ok(PoastaPoaDetailedProfile {
            result,
            timings: timings.finish(),
            search_timings: search_timings.finish(),
            search_work_counters,
        })
    }

    #[allow(clippy::too_many_lines)]
    fn align_with_phase_timings<A: Alphabet, const DETAILED: bool>(
        &self,
        records: &[SequenceRecord<A>],
        mut timings: Option<&mut PhaseTimingsBuilder>,
        mut search_timings: Option<&mut PhaseTimingsBuilder>,
        mut search_work_counters: Option<&mut PoastaSearchWorkCounters>,
    ) -> Result<PoastaPoaResult<A>, PoaError> {
        measure_phase(&mut timings, "input_validation", || {
            validate_records(records)
        })?;
        let mut graph = PoaGraph::new();
        let mut node_order = Vec::new();
        let mut steps = Vec::new();
        steps
            .try_reserve_exact(records.len())
            .map_err(|_| PoaError::AllocationFailed {
                buffer: "iterative POASTA construction steps",
                elements: records.len(),
            })?;
        let mut peak_alignment_bytes = 0_usize;
        let mut peak_integration_target_bytes = 0_usize;
        let mut estimated_peak_bytes = 0_usize;
        let mut retained_graph_bytes = 0_usize;
        let mut retained_step_bytes = 0_usize;
        let mut eager_relaxed_linear = false;
        let mut workspace = PoastaWorkspace::default();

        for (record_index, record) in records.iter().enumerate() {
            let nodes_before = graph.node_count();
            let classes_before = graph.class_count();
            let edges_before = graph.edge_count();

            let (step_alignment, targets, alignment_bytes) = if graph.node_count() == 0 {
                (
                    PoastaPoaStepAlignment::Seed,
                    measure_phase(&mut timings, "integration_targets", || {
                        empty_targets(record.sequence().len())
                    })?,
                    0,
                )
            } else {
                let (alignment, verification) = self
                    .kernel
                    .align_validated_with_report_profiled::<A, DETAILED>(
                        &graph,
                        record.sequence(),
                        &node_order,
                        eager_relaxed_linear,
                        Some(&mut workspace),
                        timings.as_deref_mut(),
                        search_timings.as_deref_mut(),
                        search_work_counters.as_deref_mut(),
                    )?;
                eager_relaxed_linear |= alignment.statistics().linear_heuristic_cells() > 0;
                let alignment_bytes = alignment.memory_estimate().modeled_peak_bytes();
                let targets = measure_phase(&mut timings, "integration_targets", || {
                    alignment.integration_targets(record.sequence().len())
                })?;
                (
                    PoastaPoaStepAlignment::Aligned {
                        alignment,
                        verification,
                    },
                    targets,
                    alignment_bytes,
                )
            };
            peak_alignment_bytes = peak_alignment_bytes.max(alignment_bytes);
            let target_bytes = targets
                .len()
                .checked_mul(size_of::<Option<PoaNodeId>>())
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "iterative POASTA integration target storage",
                })?;
            peak_integration_target_bytes = peak_integration_target_bytes.max(target_bytes);

            let (candidate, path, candidate_node_order) =
                measure_phase(&mut timings, "graph_integration", || {
                    graph.into_appended_record_trusted(record, &targets)
                })?;
            let added_nodes = candidate.node_count().checked_sub(nodes_before).ok_or(
                PoaError::ArithmeticOverflow {
                    context: "iterative POASTA node-count growth",
                },
            )?;
            let added_classes = candidate.class_count().checked_sub(classes_before).ok_or(
                PoaError::ArithmeticOverflow {
                    context: "iterative POASTA class-count growth",
                },
            )?;
            let added_edges = candidate.edge_count().checked_sub(edges_before).ok_or(
                PoaError::ArithmeticOverflow {
                    context: "iterative POASTA edge-count growth",
                },
            )?;
            retained_graph_bytes = retained_graph_bytes
                .checked_add(modeled_graph_append_bytes::<A>(
                    record,
                    added_nodes,
                    added_classes,
                    added_edges,
                )?)
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "incremental iterative POASTA graph storage",
                })?;
            let mutation_peak = retained_graph_bytes
                .checked_add(alignment_bytes)
                .and_then(|bytes| bytes.checked_add(target_bytes))
                .and_then(|bytes| bytes.checked_add(retained_step_bytes))
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "iterative POASTA mutation peak",
                })?;
            estimated_peak_bytes = estimated_peak_bytes.max(mutation_peak);
            graph = candidate;
            node_order = candidate_node_order;
            let step = PoastaPoaStep {
                record_index,
                path,
                alignment: step_alignment,
                nodes_before,
                nodes_after: graph.node_count(),
                classes_before,
                classes_after: graph.class_count(),
                edges_before,
                edges_after: graph.edge_count(),
            };
            retained_step_bytes = retained_step_bytes
                .checked_add(modeled_step_bytes(&step)?)
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "incremental iterative POASTA step storage",
                })?;
            steps.push(step);
        }

        // Every mutation above starts from this locally constructed graph,
        // consumes an independently replayed traceback, and goes through the
        // checked graph integration API. Keep the exhaustive audit in debug
        // builds, where it catches violated construction invariants, without
        // rescanning all paths and occurrences in optimized production runs.
        let alignment = measure_phase(
            &mut timings,
            "alignment_serialization",
            || -> Result<_, PoaError> {
                if cfg!(debug_assertions) {
                    graph.validate()?;
                }
                let class_order = graph.stable_class_topological_order()?;
                let alignment = graph.to_multiple_alignment_with_class_order(&class_order)?;
                verify_source_recovery(&alignment, records)?;
                Ok(alignment)
            },
        )?;
        let graph_storage_bytes = graph.modeled_storage_bytes()?;
        let rectangular_alignment_bytes = modeled_alignment_bytes(&alignment)?;
        let step_evidence_bytes = modeled_steps_bytes(&steps)?;
        debug_assert_eq!(retained_graph_bytes, graph_storage_bytes);
        debug_assert_eq!(retained_step_bytes, step_evidence_bytes);
        let returned_bytes = graph_storage_bytes
            .checked_add(rectangular_alignment_bytes)
            .and_then(|bytes| bytes.checked_add(step_evidence_bytes))
            .ok_or(PoaError::ArithmeticOverflow {
                context: "returned iterative POASTA storage",
            })?;
        estimated_peak_bytes = estimated_peak_bytes.max(returned_bytes);

        Ok(PoastaPoaResult {
            graph,
            alignment,
            steps: steps.into_boxed_slice(),
            memory: PoastaPoaMemoryEstimate {
                peak_alignment_bytes,
                peak_integration_target_bytes,
                graph_storage_bytes,
                rectangular_alignment_bytes,
                step_evidence_bytes,
                estimated_peak_bytes,
            },
        })
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

fn modeled_graph_append_bytes<A: Alphabet>(
    record: &SequenceRecord<A>,
    added_nodes: usize,
    added_classes: usize,
    added_edges: usize,
) -> Result<usize, PoaError> {
    let mut bytes = record.identifier().len();
    for (elements, element_size) in [
        (added_nodes, size_of::<crate::PoaNode<A>>()),
        (added_nodes, size_of::<PoaNodeId>()),
        (added_classes, size_of::<crate::PoaAlignmentClass>()),
        (added_edges, size_of::<crate::PoaEdge>()),
        (added_edges, size_of::<PoaNodeId>() * 2),
        (1, size_of::<crate::PoaPath>()),
        (
            record.sequence().len(),
            size_of::<crate::PoaOccurrence>() + size_of::<PoaNodeId>() + size_of::<u8>(),
        ),
    ] {
        bytes = bytes
            .checked_add(element_size.checked_mul(elements).ok_or(
                PoaError::ArithmeticOverflow {
                    context: "incremental iterative POASTA graph component",
                },
            )?)
            .ok_or(PoaError::ArithmeticOverflow {
                context: "incremental iterative POASTA graph storage",
            })?;
    }
    Ok(bytes)
}

fn validate_records<A: Alphabet>(records: &[SequenceRecord<A>]) -> Result<(), PoaError> {
    if records.is_empty() {
        return Err(PoaError::EmptyCollection);
    }
    let mut identifiers = Vec::new();
    identifiers
        .try_reserve_exact(records.len())
        .map_err(|_| PoaError::AllocationFailed {
            buffer: "iterative POASTA identifier audit",
            elements: records.len(),
        })?;
    identifiers.extend(records.iter().map(SequenceRecord::identifier));
    identifiers.sort_unstable();
    if let Some(duplicate) = identifiers
        .windows(2)
        .find(|pair| pair[0] == pair[1])
        .map(|pair| pair[1])
    {
        return Err(PoaError::DuplicateIdentifier {
            identifier: duplicate.to_owned(),
        });
    }
    Ok(())
}

fn empty_targets(length: usize) -> Result<Vec<Option<PoaNodeId>>, PoaError> {
    let mut targets = Vec::new();
    targets
        .try_reserve_exact(length)
        .map_err(|_| PoaError::AllocationFailed {
            buffer: "iterative POASTA seed targets",
            elements: length,
        })?;
    targets.resize(length, None);
    Ok(targets)
}

fn modeled_steps_bytes(steps: &[PoastaPoaStep]) -> Result<usize, PoaError> {
    let mut bytes = 0_usize;
    for step in steps {
        bytes =
            bytes
                .checked_add(modeled_step_bytes(step)?)
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "iterative POASTA retained evidence",
                })?;
    }
    Ok(bytes)
}

fn modeled_step_bytes(step: &PoastaPoaStep) -> Result<usize, PoaError> {
    let mut bytes = size_of::<PoastaPoaStep>();
    let Some(alignment) = step.alignment.alignment() else {
        return Ok(bytes);
    };
    bytes = bytes
        .checked_add(
            alignment
                .operations()
                .len()
                .checked_mul(size_of::<PoaOperation>())
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "retained POASTA operations",
                })?,
        )
        .and_then(|value| {
            alignment
                .aligned_pairs()
                .len()
                .checked_mul(size_of::<crate::PoaAlignedPair>())
                .and_then(|pairs| value.checked_add(pairs))
        })
        .ok_or(PoaError::ArithmeticOverflow {
            context: "iterative POASTA retained evidence",
        })?;
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, ProteinSequence, RnaSequence};

    use super::*;

    fn costs() -> PoastaAffineCosts {
        PoastaAffineCosts::new(4, 6, 2).expect("valid costs")
    }

    fn dna_record(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(
            identifier,
            DnaSequence::from_ascii(sequence.as_bytes()).expect("DNA sequence"),
        )
        .expect("record")
    }

    #[test]
    fn accelerated_iterative_pipeline_is_valid_deterministic_and_source_preserving() {
        let records = [
            dna_record("first", "ACGT"),
            dna_record("second", "ATGT"),
            dna_record("third", "ACGGT"),
        ];
        let aligner = PoastaPoa::new(costs());
        let first = aligner.align(&records).expect("first run");
        let second = aligner.align(&records).expect("second run");
        assert_eq!(first, second);
        first.graph().validate().expect("valid graph");
        assert_eq!(
            first.alignment(),
            &first
                .graph()
                .to_multiple_alignment()
                .expect("independently audited rectangular alignment")
        );
        assert_eq!(first.graph().path_count(), records.len());
        for (row, source) in first.alignment().rows().iter().zip(&records) {
            assert_eq!(row.identifier(), source.identifier());
            assert_eq!(row.to_ungapped_ascii(), source.sequence().to_string());
        }
        assert!(
            first.steps()[1]
                .alignment()
                .alignment()
                .expect("aligned step")
                .statistics()
                .discovered_states()
                > 0
        );
    }

    #[test]
    fn optimized_affine_certificate_handles_a_band_boundary_traceback() {
        let sequences = [
            "ATAATGATACTTCCGTAATTCGCGGCTCGGCCATAATGAAGGCGGAGAGGTGAGATTCCTATGA",
            "CTAACGATACTTCTGTCCAGCCACAGCTTCTTCACGCAATGGGGGCAGCTCGGAGCGATCCTCGGAGGGGTTCAAGTCCCATGG",
            "ACAATGGTACTTCTGTCCAGCCACAGCCCCGGGTGGAACCGGGATCACGCAATGGGGGCAGCTGCGGGAGATCCCCGCGGGGGTGAGAGGCCCATGG",
            "CAAAAGGTACTTCTGTAACTCGCAGCCCGGTCATAACGAAGACGGGGTTGCTGAAATGCCAATGA",
            "ACAATGATACTTCTCGAAAGAGCTGGGCGGAAGCGCAGGGGTGGGACTCCCATGG",
            "ACAATGATACTTCTGCCGAGTCCCAGACATCGCAAGGACGGCAGCCCGCATGACGGGGGAGGTGAAAAGCCTATGT",
            "ATAATGAGACTTCTGTCAATTCGAAATTGAGCGTTGAAAGCGTCGCACGTCAGATGGCAGTCCGAGCGTATTGCTGGGAAGGTGGAATCCCTATGG",
            "ATAATGATACTCCCTTATCGCCATAGCCCGAGCGTTCAAAGCGCCGCAGGCAATGGGTAGGGCTGCGCGAGAACCGCGCAGGGGTGCAATTCCCGTGA",
            "ACAATGATACTCCTGTTCAGTCACAGTCCGAGCGTGATAGCGTTTTCCAACGCAAAGCCGTCGCAGGCAATGAGAGCAGCCGACAGAGAACCTGCCGGGGGTGAAATTCCCGTGG",
            "ATAATGATACTTCTGCCTAGCCGTAGACCGCGCTTTAAGGGCGTGTATGCAAGCATCGGGGGCAGTCCGGGCAGCGCCGGGAAGGTGAAAATCCTATGG",
            "ATAATGGTACTCCTGTCCAGCCACAGCCCTGCTTAAGATCAGAGTAGTCCATAACATTTGGAAACCATCTGGTTGCAGGATCACGCAATGGGGGCAGCTCGGAGAGATCCTCGGAGGGGTTCGATTCCCGTGG",
            "ACAATGATACTCCTGTACGCAGCTCCGGGCGTTCCCCGGAGGGGTGAGAGTCCCGTGA",
            "ATCATGATACTTCTGCGCAGTCGAGATTCAGTTCCAGAGCGCAGCCCGGGCGTTGCCGGGGTGGTGAAATTCCAATGA",
        ];
        let records = sequences
            .iter()
            .enumerate()
            .map(|(index, sequence)| dna_record(&format!("band-{index:02}"), sequence))
            .collect::<Vec<_>>();

        let result = PoastaPoa::optimized(costs())
            .align(&records)
            .expect("the compact affine band must retain a legal traceback");
        result.graph().validate().expect("valid graph");
        for (row, source) in result.alignment().rows().iter().zip(&records) {
            assert_eq!(row.identifier(), source.identifier());
            assert_eq!(row.to_ungapped_ascii(), source.sequence().to_string());
        }
    }

    #[test]
    #[allow(clippy::too_many_lines)]
    fn optimized_branch_heavy_pipeline_matches_the_published_output_exactly() {
        let mut state = 0x6a09_e667_f3bc_c909_u64;
        let backbone = (0..96)
            .map(|_| {
                state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
                let mut mixed = state;
                mixed = (mixed ^ (mixed >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
                mixed = (mixed ^ (mixed >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
                mixed ^= mixed >> 31;
                b"ACGT"[((mixed >> 62) & 3) as usize]
            })
            .collect::<Vec<_>>();
        let records = (0_usize..16)
            .map(|record_index| {
                let mut sequence = backbone.clone();
                if record_index > 0 {
                    for position in (record_index % 7..sequence.len()).step_by(7) {
                        let alphabet_index = b"ACGT"
                            .iter()
                            .position(|symbol| *symbol == sequence[position])
                            .expect("canonical DNA");
                        sequence[position] = b"ACGT"[(alphabet_index + 1 + record_index % 3) % 4];
                    }
                }
                SequenceRecord::new(
                    format!("branch-{record_index:02}"),
                    DnaSequence::from_ascii(sequence).expect("DNA sequence"),
                )
                .expect("record")
            })
            .collect::<Vec<_>>();

        let published = PoastaPoa::new(costs())
            .align(&records)
            .expect("published branch-heavy run");
        let optimized = PoastaPoa::optimized(costs())
            .align(&records)
            .expect("optimized branch-heavy run");

        assert_eq!(optimized.steps().len(), published.steps().len());
        for (optimized_step, published_step) in optimized.steps().iter().zip(published.steps()) {
            assert_eq!(optimized_step.path(), published_step.path());
            if let (Some(optimized_alignment), Some(published_alignment)) = (
                optimized_step.alignment().alignment(),
                published_step.alignment().alignment(),
            ) {
                assert_eq!(
                    optimized_alignment.operations(),
                    published_alignment.operations(),
                    "operations differ at record {}; optimized statistics: {:?}",
                    optimized_step.record_index(),
                    optimized_alignment.statistics()
                );
            }
            assert_eq!(optimized_step.nodes_before(), published_step.nodes_before());
            assert_eq!(optimized_step.nodes_after(), published_step.nodes_after());
            assert_eq!(
                optimized_step.classes_before(),
                published_step.classes_before()
            );
            assert_eq!(
                optimized_step.classes_after(),
                published_step.classes_after()
            );
            assert_eq!(optimized_step.edges_before(), published_step.edges_before());
            assert_eq!(
                optimized_step.edges_after(),
                published_step.edges_after(),
                "edge counts differ after record {}",
                optimized_step.record_index()
            );
            match (
                optimized_step.alignment().alignment(),
                published_step.alignment().alignment(),
            ) {
                (Some(optimized_alignment), Some(published_alignment)) => {
                    assert_eq!(optimized_alignment.score(), published_alignment.score());
                    assert_eq!(
                        optimized_alignment.operations(),
                        published_alignment.operations(),
                        "operations differ at record {}",
                        optimized_step.record_index()
                    );
                    assert_eq!(
                        optimized_alignment.aligned_pairs(),
                        published_alignment.aligned_pairs()
                    );
                    assert_eq!(
                        optimized_step.alignment().verification(),
                        published_step.alignment().verification()
                    );
                }
                (None, None) => {}
                _ => panic!("published and optimized seed steps differ"),
            }
        }
        assert_eq!(optimized.graph(), published.graph());
        assert_eq!(optimized.alignment(), published.alignment());
        assert!(optimized.steps().iter().any(|step| {
            step.alignment().alignment().is_some_and(|alignment| {
                alignment.statistics().linear_heuristic_cells() > 0
                    && (alignment.statistics().used_relaxed_linear_certificate()
                        || alignment.statistics().used_exact_affine_certificate())
            })
        }));
    }

    #[test]
    fn protein_rna_and_dna_share_the_iterative_pipeline() {
        let dna = [dna_record("dna-a", "ACGT"), dna_record("dna-b", "AGGT")];
        assert_eq!(
            PoastaPoa::new(costs())
                .align(&dna)
                .unwrap()
                .graph()
                .path_count(),
            2
        );

        let rna = [
            SequenceRecord::new("rna-a", RnaSequence::from_ascii(b"ACGU").unwrap()).unwrap(),
            SequenceRecord::new("rna-b", RnaSequence::from_ascii(b"AGGU").unwrap()).unwrap(),
        ];
        assert_eq!(
            PoastaPoa::new(costs())
                .align(&rna)
                .unwrap()
                .graph()
                .path_count(),
            2
        );

        let protein = [
            SequenceRecord::new("protein-a", ProteinSequence::from_ascii(b"MKWV").unwrap())
                .unwrap(),
            SequenceRecord::new("protein-b", ProteinSequence::from_ascii(b"MRWV").unwrap())
                .unwrap(),
        ];
        assert_eq!(
            PoastaPoa::new(costs())
                .align(&protein)
                .unwrap()
                .graph()
                .path_count(),
            2
        );
    }

    #[test]
    fn empty_and_duplicate_collections_are_typed_errors() {
        assert_eq!(
            PoastaPoa::new(costs()).align::<Dna>(&[]),
            Err(PoaError::EmptyCollection)
        );
        let duplicate = [dna_record("same", "A"), dna_record("same", "C")];
        assert_eq!(
            PoastaPoa::new(costs()).align(&duplicate),
            Err(PoaError::DuplicateIdentifier {
                identifier: "same".to_owned()
            })
        );
    }
}
