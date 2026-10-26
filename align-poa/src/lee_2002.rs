use core::mem::size_of;

use align_core::{Alphabet, LinearGap, SubstitutionScorer, Symbol};
use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord};

use crate::alignment::verify_sequence_graph_alignment_validated;
use crate::{
    LocalPoa, PoaAlignedPair, PoaError, PoaGraph, PoaOperation, PoaPathId, PoaVerificationReport,
    SequenceGraphAlignment,
};

/// Immutable identifier of the original local linear-gap iterative POA profile.
pub const LEE_2002_POA_LOCAL_LINEAR_V1: &str = "LEE-2002-POA-local-linear-v1";

/// Auditable evidence for one sequence incorporated into the PO-MSA graph.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LeePoaStep {
    record_index: usize,
    path: PoaPathId,
    alignment: SequenceGraphAlignment,
    verification: PoaVerificationReport,
    nodes_before: usize,
    nodes_after: usize,
    classes_before: usize,
    classes_after: usize,
    edges_before: usize,
    edges_after: usize,
}

impl LeePoaStep {
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

    /// Returns the exact sequence-to-current-graph alignment.
    #[must_use]
    pub const fn sequence_graph_alignment(&self) -> &SequenceGraphAlignment {
        &self.alignment
    }

    /// Returns independent replay evidence for the sequence-to-graph alignment.
    #[must_use]
    pub const fn verification(&self) -> PoaVerificationReport {
        self.verification
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

/// Checked memory evidence for one complete iterative POA call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[allow(clippy::struct_field_names)]
pub struct LeePoaMemoryEstimate {
    peak_alignment_bytes: usize,
    graph_storage_bytes: usize,
    rectangular_alignment_bytes: usize,
    step_evidence_bytes: usize,
    estimated_peak_bytes: usize,
}

impl LeePoaMemoryEstimate {
    /// Returns the largest scalar sequence-to-DAG element-buffer requirement.
    #[must_use]
    pub const fn peak_alignment_bytes(self) -> usize {
        self.peak_alignment_bytes
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

    /// Returns modeled bytes retained by per-sequence evidence.
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

/// Complete result of the original 2002 iterative POA profile.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LeePoaResult<A: Alphabet> {
    graph: PoaGraph<A>,
    alignment: MultipleAlignment<A>,
    steps: Box<[LeePoaStep]>,
    memory: LeePoaMemoryEstimate,
}

impl<A: Alphabet> LeePoaResult<A> {
    /// Returns the authoritative partial-order graph.
    #[must_use]
    pub const fn graph(&self) -> &PoaGraph<A> {
        &self.graph
    }

    /// Returns the deterministic rectangular alignment-class serialization.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        &self.alignment
    }

    /// Returns one auditable construction step per input record.
    #[must_use]
    pub const fn steps(&self) -> &[LeePoaStep] {
        &self.steps
    }

    /// Returns checked modeled memory evidence.
    #[must_use]
    pub const fn memory_estimate(&self) -> LeePoaMemoryEstimate {
        self.memory
    }

    /// Separates the authoritative graph and rectangular serialization.
    #[must_use]
    pub fn into_graph_and_alignment(self) -> (PoaGraph<A>, MultipleAlignment<A>) {
        (self.graph, self.alignment)
    }
}

/// Original iterative local linear-gap Partial Order Alignment.
#[derive(Clone, Debug)]
pub struct LeePoa<S> {
    substitution: S,
    gap: LinearGap,
}

impl<S> LeePoa<S> {
    /// Creates the immutable 2002 profile with explicit scoring.
    #[must_use]
    pub const fn new(substitution: S, gap: LinearGap) -> Self {
        Self { substitution, gap }
    }

    /// Returns the immutable scientific profile identifier.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        LEE_2002_POA_LOCAL_LINEAR_V1
    }

    /// Returns the configured substitution scorer.
    #[must_use]
    pub const fn substitution(&self) -> &S {
        &self.substitution
    }

    /// Returns the configured non-positive linear gap score.
    #[must_use]
    pub const fn gap(&self) -> LinearGap {
        self.gap
    }

    /// Builds a PO-MSA by adding records one at a time in caller order.
    ///
    /// The graph is authoritative. A same-symbol allele in an aligned class is
    /// fused according to the 2002 paper; unaligned local prefixes and suffixes
    /// become new graph branches. Input order is intentionally preserved because
    /// the published iterative method is order-sensitive.
    ///
    /// # Errors
    ///
    /// Returns [`PoaError`] for an empty collection, duplicate identifier,
    /// invalid graph or traceback, allocation failure, arithmetic overflow, or
    /// failed source recovery.
    #[allow(clippy::too_many_lines)]
    pub fn align<A>(&self, records: &[SequenceRecord<A>]) -> Result<LeePoaResult<A>, PoaError>
    where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        if records.is_empty() {
            return Err(PoaError::EmptyCollection);
        }
        let mut identifiers = Vec::new();
        identifiers
            .try_reserve_exact(records.len())
            .map_err(|_| PoaError::AllocationFailed {
                buffer: "Lee POA identifier audit",
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

        let local = LocalPoa::new(&self.substitution, self.gap);
        let mut graph = PoaGraph::new();
        let mut node_order = Vec::new();
        let mut class_order = Vec::new();
        let mut steps = Vec::new();
        steps
            .try_reserve_exact(records.len())
            .map_err(|_| PoaError::AllocationFailed {
                buffer: "Lee POA construction steps",
                elements: records.len(),
            })?;
        let mut peak_alignment_bytes = 0_usize;
        let mut estimated_peak_bytes = 0_usize;

        for (record_index, record) in records.iter().enumerate() {
            let nodes_before = graph.node_count();
            let classes_before = graph.class_count();
            let edges_before = graph.edge_count();
            let alignment = if graph.node_count() == 0 || record.sequence().is_empty() {
                SequenceGraphAlignment::canonical_empty(LocalPoa::<&S>::memory_requirement(
                    graph.node_count(),
                    record.sequence().len(),
                )?)
            } else {
                local.align_with_order(&graph, record.sequence(), &node_order)?
            };
            let verification = verify_sequence_graph_alignment_validated(
                &alignment,
                &graph,
                record.sequence(),
                &self.substitution,
                self.gap,
            )?;
            peak_alignment_bytes =
                peak_alignment_bytes.max(alignment.memory_estimate().peak_buffer_bytes());
            let targets = alignment.integration_targets(record.sequence().len())?;
            let target_bytes = targets
                .len()
                .checked_mul(size_of::<Option<crate::PoaNodeId>>())
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "POA integration target storage",
                })?;
            let (candidate, path, candidate_node_order, candidate_class_order) =
                graph.into_appended_record(record, &targets)?;
            let graph_after_bytes = candidate.modeled_storage_bytes()?;
            let retained_step_bytes = modeled_steps_bytes(&steps)?;
            let mutation_peak = graph_after_bytes
                .checked_add(alignment.memory_estimate().peak_buffer_bytes())
                .and_then(|bytes| bytes.checked_add(target_bytes))
                .and_then(|bytes| bytes.checked_add(retained_step_bytes))
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "iterative POA mutation peak",
                })?;
            estimated_peak_bytes = estimated_peak_bytes.max(mutation_peak);
            graph = candidate;
            node_order = candidate_node_order;
            class_order = candidate_class_order;
            steps.push(LeePoaStep {
                record_index,
                path,
                alignment,
                verification,
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
                context: "returned iterative POA storage",
            })?;
        estimated_peak_bytes = estimated_peak_bytes.max(returned_bytes);
        Ok(LeePoaResult {
            graph,
            alignment,
            steps: steps.into_boxed_slice(),
            memory: LeePoaMemoryEstimate {
                peak_alignment_bytes,
                graph_storage_bytes,
                rectangular_alignment_bytes,
                step_evidence_bytes,
                estimated_peak_bytes,
            },
        })
    }
}

pub(crate) fn verify_source_recovery<A: Alphabet>(
    alignment: &MultipleAlignment<A>,
    records: &[SequenceRecord<A>],
) -> Result<(), PoaError> {
    if alignment.sequence_count() != records.len() {
        return verification("rectangular serialization changed the sequence count");
    }
    for (row, record) in alignment.rows().iter().zip(records) {
        if row.identifier() != record.identifier() {
            return verification("rectangular serialization changed input order or identifiers");
        }
        if !row
            .symbols()
            .iter()
            .filter_map(|cell| *cell)
            .eq(record.sequence().iter())
        {
            return verification(format!(
                "rectangular serialization does not recover source {:?}",
                record.identifier()
            ));
        }
    }
    Ok(())
}

pub(crate) fn modeled_alignment_bytes<A: Alphabet>(
    alignment: &MultipleAlignment<A>,
) -> Result<usize, PoaError> {
    let cells = alignment
        .sequence_count()
        .checked_mul(alignment.width())
        .ok_or(PoaError::ArithmeticOverflow {
            context: "rectangular POA cell count",
        })?;
    let mut bytes =
        cells
            .checked_mul(size_of::<Option<Symbol<A>>>())
            .ok_or(PoaError::ArithmeticOverflow {
                context: "rectangular POA cell bytes",
            })?;
    bytes = bytes
        .checked_add(
            alignment
                .sequence_count()
                .checked_mul(size_of::<AlignedRow<A>>())
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "rectangular POA row bytes",
                })?,
        )
        .ok_or(PoaError::ArithmeticOverflow {
            context: "rectangular POA storage",
        })?;
    for row in alignment.rows() {
        bytes = bytes
            .checked_add(row.identifier().len())
            .ok_or(PoaError::ArithmeticOverflow {
                context: "rectangular POA identifiers",
            })?;
    }
    Ok(bytes)
}

fn modeled_steps_bytes(steps: &[LeePoaStep]) -> Result<usize, PoaError> {
    let mut bytes =
        steps
            .len()
            .checked_mul(size_of::<LeePoaStep>())
            .ok_or(PoaError::ArithmeticOverflow {
                context: "Lee POA step arena",
            })?;
    for step in steps {
        bytes = bytes
            .checked_add(
                step.alignment
                    .operations()
                    .len()
                    .checked_mul(size_of::<PoaOperation>())
                    .ok_or(PoaError::ArithmeticOverflow {
                        context: "retained POA operations",
                    })?,
            )
            .and_then(|total| {
                total.checked_add(
                    step.alignment
                        .aligned_pairs()
                        .len()
                        .checked_mul(size_of::<PoaAlignedPair>())?,
                )
            })
            .ok_or(PoaError::ArithmeticOverflow {
                context: "Lee POA step evidence",
            })?;
    }
    Ok(bytes)
}

fn verification<T>(message: impl Into<String>) -> Result<T, PoaError> {
    Err(PoaError::VerificationFailed {
        message: message.into(),
    })
}

#[cfg(test)]
mod tests {
    use align_core::{
        Dna, DnaSequence, LinearGap, MatchMismatch, Protein, ProteinSequence, Rna, RnaSequence,
    };
    use align_msa::SequenceRecord;

    use super::{LEE_2002_POA_LOCAL_LINEAR_V1, LeePoa};
    use crate::PoaError;

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    #[test]
    fn identical_sequences_fuse_into_one_supported_path_shape() {
        let records = [dna("a", "ACGT"), dna("b", "ACGT")];
        let result = LeePoa::new(MatchMismatch::new(2, -1), LinearGap::new(-2).unwrap())
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
        assert_eq!(
            result.alignment().row(0).unwrap().to_ascii_with_gaps(),
            "ACGT"
        );
        assert_eq!(
            result.alignment().row(1).unwrap().to_ascii_with_gaps(),
            "ACGT"
        );
        assert_eq!(result.steps()[1].sequence_graph_alignment().score(), 8);
        assert_eq!(result.steps()[1].verification().recomputed_score(), 8);
    }

    #[test]
    fn mismatching_aligned_residues_form_one_allelic_class() {
        let records = [dna("a", "ACG"), dna("b", "ATG")];
        let result = LeePoa::new(MatchMismatch::new(2, 1), LinearGap::new(-2).unwrap())
            .align(&records)
            .unwrap();
        assert_eq!(result.graph().node_count(), 4);
        assert_eq!(result.graph().class_count(), 3);
        assert_eq!(result.graph().classes()[1].nodes().len(), 2);
        assert_eq!(
            result.alignment().row(0).unwrap().to_ascii_with_gaps(),
            "ACG"
        );
        assert_eq!(
            result.alignment().row(1).unwrap().to_ascii_with_gaps(),
            "ATG"
        );
    }

    #[test]
    fn empty_and_zero_hit_sequences_remain_recoverable() {
        let records = [dna("empty", ""), dna("a", "AA"), dna("c", "CC")];
        let result = LeePoa::new(MatchMismatch::new(2, -3), LinearGap::new(-2).unwrap())
            .align(&records)
            .unwrap();
        assert_eq!(result.alignment().row(0).unwrap().to_ungapped_ascii(), "");
        assert_eq!(result.alignment().row(1).unwrap().to_ungapped_ascii(), "AA");
        assert_eq!(result.alignment().row(2).unwrap().to_ungapped_ascii(), "CC");
        assert_eq!(result.graph().node_count(), 4);
        assert_eq!(result.graph().paths().len(), 3);
        assert!(
            result.steps()[2]
                .sequence_graph_alignment()
                .operations()
                .is_empty()
        );
    }

    #[test]
    fn common_api_supports_protein_dna_and_rna() {
        let scoring = MatchMismatch::new(2, -1);
        let gap = LinearGap::new(-2).unwrap();
        let protein = [
            SequenceRecord::<Protein>::new("p1", ProteinSequence::from_ascii("ARND").unwrap())
                .unwrap(),
            SequenceRecord::new("p2", ProteinSequence::from_ascii("ARNE").unwrap()).unwrap(),
        ];
        let rna = [
            SequenceRecord::<Rna>::new("r1", RnaSequence::from_ascii("ACGU").unwrap()).unwrap(),
            SequenceRecord::new("r2", RnaSequence::from_ascii("AGGU").unwrap()).unwrap(),
        ];
        assert_eq!(
            LeePoa::new(scoring, gap)
                .align(&protein)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            LeePoa::new(scoring, gap)
                .align(&rna)
                .unwrap()
                .alignment()
                .sequence_count(),
            2
        );
        assert_eq!(
            LeePoa::new(scoring, gap).profile_identifier(),
            LEE_2002_POA_LOCAL_LINEAR_V1
        );
    }

    #[test]
    fn duplicate_and_empty_collections_are_rejected() {
        let aligner = LeePoa::new(MatchMismatch::new(2, -1), LinearGap::new(-2).unwrap());
        let empty: Vec<SequenceRecord<Dna>> = Vec::new();
        assert!(matches!(
            aligner.align(&empty),
            Err(PoaError::EmptyCollection)
        ));
        let duplicate = [dna("same", "A"), dna("same", "C")];
        assert!(matches!(
            aligner.align(&duplicate),
            Err(PoaError::DuplicateIdentifier { .. })
        ));
    }
}
