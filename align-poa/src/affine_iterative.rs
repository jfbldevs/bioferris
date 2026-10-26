use core::mem::size_of;

use align_core::{AffineGap, Alphabet, SubstitutionScorer};
use align_msa::{MultipleAlignment, SequenceRecord};

use crate::affine::verify_global_affine_sequence_graph_alignment_validated;
use crate::lee_2002::{modeled_alignment_bytes, verify_source_recovery};
use crate::{
    AffinePoaVerificationReport, GlobalAffinePoa, GlobalAffineSequenceGraphAlignment,
    PoaAlignedPair, PoaError, PoaGraph, PoaNodeId, PoaOperation, PoaPathId, PoastaAffineCosts,
};

/// Immutable identifier of the iterative global affine composite profile.
pub const LEE_2002_POASTA_2025_GLOBAL_AFFINE_ITERATIVE_V1: &str =
    "LEE-2002-POASTA-2025-global-affine-iterative-v1";

/// Whether one iterative step seeded an empty graph or performed affine alignment.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AffinePoaStepAlignment {
    /// The record was appended directly because the graph had no biological node.
    Seed,
    /// The record was globally aligned and independently replayed before mutation.
    Aligned {
        /// Exact global affine sequence-to-current-graph result.
        alignment: GlobalAffineSequenceGraphAlignment,
        /// Independent replay evidence retained with the step.
        verification: AffinePoaVerificationReport,
    },
}

impl AffinePoaStepAlignment {
    /// Returns whether this record directly seeded an empty biological graph.
    #[must_use]
    pub const fn is_seed(&self) -> bool {
        matches!(self, Self::Seed)
    }

    /// Returns the exact sequence-to-graph alignment when this was not a seed.
    #[must_use]
    pub const fn alignment(&self) -> Option<&GlobalAffineSequenceGraphAlignment> {
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

/// Auditable evidence for one record incorporated by global affine POA.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AffinePoaStep {
    record_index: usize,
    path: PoaPathId,
    alignment: AffinePoaStepAlignment,
    nodes_before: usize,
    nodes_after: usize,
    classes_before: usize,
    classes_after: usize,
    edges_before: usize,
    edges_after: usize,
}

impl AffinePoaStep {
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

    /// Returns seed or exact affine-alignment evidence for this step.
    #[must_use]
    pub const fn alignment(&self) -> &AffinePoaStepAlignment {
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

/// Checked memory evidence for one complete iterative affine POA call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(clippy::struct_field_names)]
pub struct AffinePoaMemoryEstimate {
    peak_alignment_bytes: usize,
    peak_integration_target_bytes: usize,
    graph_storage_bytes: usize,
    rectangular_alignment_bytes: usize,
    step_evidence_bytes: usize,
    estimated_peak_bytes: usize,
}

impl AffinePoaMemoryEstimate {
    /// Returns the largest scalar affine alignment element-buffer requirement.
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

    /// Returns the largest modeled simultaneous element-data total.
    #[must_use]
    pub const fn estimated_peak_bytes(self) -> usize {
        self.estimated_peak_bytes
    }
}

/// Complete deterministic result of iterative global affine POA.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AffinePoaResult<A: Alphabet> {
    graph: PoaGraph<A>,
    alignment: MultipleAlignment<A>,
    steps: Box<[AffinePoaStep]>,
    memory: AffinePoaMemoryEstimate,
}

impl<A: Alphabet> AffinePoaResult<A> {
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
    pub const fn steps(&self) -> &[AffinePoaStep] {
        &self.steps
    }

    /// Returns checked modeled memory evidence.
    #[must_use]
    pub const fn memory_estimate(&self) -> AffinePoaMemoryEstimate {
        self.memory
    }

    /// Separates the authoritative graph and rectangular serialization.
    #[must_use]
    pub fn into_graph_and_alignment(self) -> (PoaGraph<A>, MultipleAlignment<A>) {
        (self.graph, self.alignment)
    }
}

/// Iterative partial-order alignment using the exact global affine scalar oracle.
#[derive(Clone, Debug)]
pub struct AffinePoa<S> {
    substitution: S,
    gap: AffineGap,
}

impl<S> AffinePoa<S> {
    /// Creates the explicit composite profile with generic score-side parameters.
    #[must_use]
    pub const fn new(substitution: S, gap: AffineGap) -> Self {
        Self { substitution, gap }
    }

    /// Returns the immutable scientific profile identifier.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        LEE_2002_POASTA_2025_GLOBAL_AFFINE_ITERATIVE_V1
    }

    /// Returns the configured substitution scorer.
    #[must_use]
    pub const fn substitution(&self) -> &S {
        &self.substitution
    }

    /// Returns the configured length-one-open affine gap score.
    #[must_use]
    pub const fn gap(&self) -> AffineGap {
        self.gap
    }

    /// Builds a PO-MSA by adding records in caller order with exact global affine DP.
    ///
    /// Records encountered while the biological graph is empty seed it directly.
    /// Every other record is globally aligned, independently replayed, and then
    /// incorporated with the Lee 2002 same-symbol fusion rules.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError`] for an empty collection, duplicate identifier, invalid
    /// graph, allocation or arithmetic failure, unrepresentable score, invalid
    /// traceback, or failed independent source recovery.
    #[allow(clippy::too_many_lines)]
    pub fn align<A>(&self, records: &[SequenceRecord<A>]) -> Result<AffinePoaResult<A>, PoaError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        validate_records(records)?;
        let scalar = GlobalAffinePoa::new(&self.substitution, self.gap);
        let mut graph = PoaGraph::new();
        let mut class_order = Vec::new();
        let mut steps = Vec::new();
        steps
            .try_reserve_exact(records.len())
            .map_err(|_| PoaError::AllocationFailed {
                buffer: "affine POA construction steps",
                elements: records.len(),
            })?;
        let mut peak_alignment_bytes = 0_usize;
        let mut peak_integration_target_bytes = 0_usize;
        let mut estimated_peak_bytes = 0_usize;

        for (record_index, record) in records.iter().enumerate() {
            let nodes_before = graph.node_count();
            let classes_before = graph.class_count();
            let edges_before = graph.edge_count();

            let (step_alignment, targets, alignment_bytes) = if graph.node_count() == 0 {
                (
                    AffinePoaStepAlignment::Seed,
                    empty_targets(record.sequence().len())?,
                    0,
                )
            } else {
                let alignment = scalar.align(&graph, record.sequence())?;
                let verification = verify_global_affine_sequence_graph_alignment_validated(
                    &alignment,
                    &graph,
                    record.sequence(),
                    &self.substitution,
                    self.gap,
                )?;
                let alignment_bytes = alignment.memory_estimate().peak_buffer_bytes();
                let targets = alignment.integration_targets(record.sequence().len())?;
                (
                    AffinePoaStepAlignment::Aligned {
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
                    context: "affine POA integration target storage",
                })?;
            peak_integration_target_bytes = peak_integration_target_bytes.max(target_bytes);

            let (candidate, path, _, candidate_class_order) =
                graph.into_appended_record(record, &targets)?;
            let graph_after_bytes = candidate.modeled_storage_bytes()?;
            let retained_step_bytes = modeled_steps_bytes(&steps)?;
            let mutation_peak = graph_after_bytes
                .checked_add(alignment_bytes)
                .and_then(|bytes| bytes.checked_add(target_bytes))
                .and_then(|bytes| bytes.checked_add(retained_step_bytes))
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "iterative affine POA mutation peak",
                })?;
            estimated_peak_bytes = estimated_peak_bytes.max(mutation_peak);
            graph = candidate;
            class_order = candidate_class_order;
            steps.push(AffinePoaStep {
                record_index,
                path,
                alignment: step_alignment,
                nodes_before,
                nodes_after: graph.node_count(),
                classes_before,
                classes_after: graph.class_count(),
                edges_before,
                edges_after: graph.edge_count(),
            });
        }

        let alignment = graph.to_multiple_alignment_with_class_order(&class_order)?;
        verify_source_recovery(&alignment, records)?;
        let graph_storage_bytes = graph.modeled_storage_bytes()?;
        let rectangular_alignment_bytes = modeled_alignment_bytes(&alignment)?;
        let step_evidence_bytes = modeled_steps_bytes(&steps)?;
        let returned_bytes = graph_storage_bytes
            .checked_add(rectangular_alignment_bytes)
            .and_then(|bytes| bytes.checked_add(step_evidence_bytes))
            .ok_or(PoaError::ArithmeticOverflow {
                context: "returned iterative affine POA storage",
            })?;
        estimated_peak_bytes = estimated_peak_bytes.max(returned_bytes);

        Ok(AffinePoaResult {
            graph,
            alignment,
            steps: steps.into_boxed_slice(),
            memory: AffinePoaMemoryEstimate {
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

impl AffinePoa<PoastaAffineCosts> {
    /// Creates the strict published literal-cost iterative profile.
    #[must_use]
    pub const fn from_poasta_costs(costs: PoastaAffineCosts) -> Self {
        Self::new(costs, costs.affine_gap())
    }
}

fn validate_records<A: Alphabet>(records: &[SequenceRecord<A>]) -> Result<(), PoaError> {
    if records.is_empty() {
        return Err(PoaError::EmptyCollection);
    }
    let mut identifiers = Vec::new();
    identifiers
        .try_reserve_exact(records.len())
        .map_err(|_| PoaError::AllocationFailed {
            buffer: "affine POA identifier audit",
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
            buffer: "affine POA seed integration targets",
            elements: length,
        })?;
    targets.resize(length, None);
    Ok(targets)
}

fn modeled_steps_bytes(steps: &[AffinePoaStep]) -> Result<usize, PoaError> {
    let mut bytes = steps.len().checked_mul(size_of::<AffinePoaStep>()).ok_or(
        PoaError::ArithmeticOverflow {
            context: "affine POA step arena",
        },
    )?;
    for step in steps {
        let Some(alignment) = step.alignment.alignment() else {
            continue;
        };
        bytes = bytes
            .checked_add(
                alignment
                    .operations()
                    .len()
                    .checked_mul(size_of::<PoaOperation>())
                    .ok_or(PoaError::ArithmeticOverflow {
                        context: "retained affine POA operations",
                    })?,
            )
            .and_then(|total| {
                total.checked_add(
                    alignment
                        .aligned_pairs()
                        .len()
                        .checked_mul(size_of::<PoaAlignedPair>())?,
                )
            })
            .ok_or(PoaError::ArithmeticOverflow {
                context: "affine POA step evidence",
            })?;
    }
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, MatchMismatch, ProteinSequence, RnaSequence};
    use align_msa::SequenceRecord;

    use super::{
        AffinePoa, AffinePoaStepAlignment, LEE_2002_POASTA_2025_GLOBAL_AFFINE_ITERATIVE_V1,
    };
    use crate::{PoaError, PoastaAffineCosts};

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    fn costs() -> PoastaAffineCosts {
        PoastaAffineCosts::new(4, 6, 2).unwrap()
    }

    #[test]
    fn identical_sequences_seed_then_fuse_with_zero_cost() {
        let records = [dna("first", "ACGT"), dna("second", "ACGT")];
        let result = AffinePoa::from_poasta_costs(costs())
            .align(&records)
            .unwrap();
        assert_eq!(result.graph().node_count(), 4);
        assert_eq!(result.graph().class_count(), 4);
        assert!(
            result
                .graph()
                .edges()
                .iter()
                .all(|edge| edge.support() == 2)
        );
        assert!(result.steps()[0].alignment().is_seed());
        assert_eq!(
            result.steps()[1].alignment().alignment().unwrap().score(),
            0
        );
        assert_eq!(
            result.steps()[1]
                .alignment()
                .verification()
                .unwrap()
                .recomputed_score(),
            0
        );
        assert_eq!(
            result.alignment().row(1).unwrap().to_ascii_with_gaps(),
            "ACGT"
        );
    }

    #[test]
    fn mismatch_forms_an_allele_and_internal_indels_form_branches() {
        let records = [
            dna("reference", "ACGT"),
            dna("allele", "ATGT"),
            dna("deletion", "AGT"),
            dna("insertion", "ACAGT"),
        ];
        let result = AffinePoa::from_poasta_costs(costs())
            .align(&records)
            .unwrap();
        result.graph().validate().unwrap();
        assert_eq!(
            result.steps()[1].alignment().alignment().unwrap().score(),
            -4
        );
        for (row, record) in result.alignment().rows().iter().zip(&records) {
            assert_eq!(row.identifier(), record.identifier());
            assert_eq!(row.to_ungapped_ascii(), record.sequence().to_string());
        }
        assert!(
            result
                .graph()
                .classes()
                .iter()
                .any(|class| class.nodes().len() > 1)
        );
    }

    #[test]
    fn leading_empty_records_remain_seeds_and_later_empty_query_is_aligned() {
        let records = [
            dna("empty-first", ""),
            dna("reference", "ACG"),
            dna("empty-last", ""),
        ];
        let result = AffinePoa::from_poasta_costs(costs())
            .align(&records)
            .unwrap();
        assert!(result.steps()[0].alignment().is_seed());
        assert!(result.steps()[1].alignment().is_seed());
        let last = result.steps()[2].alignment().alignment().unwrap();
        assert_eq!(last.score(), -12);
        assert_eq!(result.alignment().row(0).unwrap().to_ungapped_ascii(), "");
        assert_eq!(result.alignment().row(2).unwrap().to_ungapped_ascii(), "");
    }

    #[test]
    fn protein_dna_and_rna_use_one_iterative_pipeline() {
        let protein = [
            SequenceRecord::new("p1", ProteinSequence::from_ascii("ARND").unwrap()).unwrap(),
            SequenceRecord::new("p2", ProteinSequence::from_ascii("ARXD").unwrap()).unwrap(),
        ];
        let protein_result = AffinePoa::from_poasta_costs(costs())
            .align(&protein)
            .unwrap();
        assert_eq!(
            protein_result.steps()[1]
                .alignment()
                .alignment()
                .unwrap()
                .score(),
            -4
        );

        let rna = [
            SequenceRecord::new("r1", RnaSequence::from_ascii("ACGU").unwrap()).unwrap(),
            SequenceRecord::new("r2", RnaSequence::from_ascii("AGU").unwrap()).unwrap(),
        ];
        let rna_result = AffinePoa::from_poasta_costs(costs()).align(&rna).unwrap();
        assert_eq!(
            rna_result.alignment().row(1).unwrap().to_ungapped_ascii(),
            "AGU"
        );
    }

    #[test]
    fn duplicate_and_empty_collections_are_rejected_before_alignment() {
        let aligner = AffinePoa::from_poasta_costs(costs());
        assert_eq!(aligner.align::<Dna>(&[]), Err(PoaError::EmptyCollection));
        assert_eq!(
            aligner.align(&[dna("same", "A"), dna("same", "C")]),
            Err(PoaError::DuplicateIdentifier {
                identifier: "same".to_owned(),
            })
        );
    }

    #[test]
    fn generic_constructor_and_memory_evidence_are_public_and_deterministic() {
        let records = [dna("a", "ACGT"), dna("b", "ACAGT")];
        let aligner = AffinePoa::new(MatchMismatch::new(0, -4), costs().affine_gap());
        assert_eq!(
            aligner.profile_identifier(),
            LEE_2002_POASTA_2025_GLOBAL_AFFINE_ITERATIVE_V1
        );
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        assert_eq!(first, second);
        assert!(first.memory_estimate().peak_alignment_bytes() > 0);
        assert!(first.memory_estimate().peak_integration_target_bytes() > 0);
        assert!(first.memory_estimate().estimated_peak_bytes() > 0);
        assert!(matches!(
            first.steps()[1].alignment(),
            AffinePoaStepAlignment::Aligned { .. }
        ));
    }
}
