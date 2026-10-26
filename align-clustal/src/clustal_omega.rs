use core::any::TypeId;
use core::marker::PhantomData;
use core::mem::size_of;
use std::collections::HashSet;

use align_core::{Alphabet, PhaseTimingsBuilder, Profiled, Protein, Symbol};
use align_msa::{
    AlignedRow, DistanceMatrix, GuideTree, MsaError, MultipleAlignment, SequenceRecord,
};

use crate::omega_distance::{
    DISTANCE_SCALE, FixedTreeNode, full_distance_matrix, mbed, one_tuple_distance, to_public_tree,
    upgma,
};
use crate::omega_hmm::{WorkingProfile, align_profiles};
use crate::{
    ClustalAlphabet, ClustalOmegaError, ClustalOmegaParameters, OmegaEngine, OmegaGuideTreeMode,
    OmegaMoleculeProfile,
};

/// Effective guide-tree construction used by one Clustal Omega call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OmegaGuideTreeMethod {
    /// An explicitly requested complete one-tuple distance matrix and UPGMA.
    Full,
    /// Default mBed behavior collapsed to a complete matrix for a small input.
    MbedEffectiveFull,
    /// Modified mBed embedding, bounded leaf clusters, and grafted UPGMA trees.
    Mbed,
}

/// Auditable evidence retained from Clustal Omega guide-tree construction.
#[derive(Clone, Debug)]
pub struct OmegaGuideEvidence {
    method: OmegaGuideTreeMethod,
    distances: Option<DistanceMatrix>,
    seed_indices: Box<[usize]>,
    clusters: Box<[Box<[usize]>]>,
}

impl OmegaGuideEvidence {
    /// Returns the guide-tree construction actually used.
    #[must_use]
    pub const fn method(&self) -> OmegaGuideTreeMethod {
        self.method
    }

    /// Returns the complete one-tuple matrix when one was materialized.
    #[must_use]
    pub const fn distances(&self) -> Option<&DistanceMatrix> {
        self.distances.as_ref()
    }

    /// Returns systematic mBed seed indices in deterministic sampling order.
    #[must_use]
    pub const fn seed_indices(&self) -> &[usize] {
        &self.seed_indices
    }

    /// Returns mBed leaf-cluster memberships ordered by minimum input index.
    #[must_use]
    pub fn clusters(&self) -> impl ExactSizeIterator<Item = &[usize]> {
        self.clusters.iter().map(AsRef::as_ref)
    }
}

/// One independently replayable progressive HMM-HMM merge.
#[derive(Clone, Debug, PartialEq)]
pub struct OmegaProgressiveStep {
    left_members: Box<[usize]>,
    right_members: Box<[usize]>,
    engine: OmegaEngine,
    score: f64,
    left_core: core::ops::Range<usize>,
    right_core: core::ops::Range<usize>,
    result_width: usize,
    workspace_bytes: usize,
    posterior_bounds: Option<(f64, f64)>,
}

impl OmegaProgressiveStep {
    /// Returns original input indices below the left child profile.
    #[must_use]
    pub const fn left_members(&self) -> &[usize] {
        &self.left_members
    }

    /// Returns original input indices below the right child profile.
    #[must_use]
    pub const fn right_members(&self) -> &[usize] {
        &self.right_members
    }

    /// Returns the HMM-HMM engine selected for this merge.
    #[must_use]
    pub const fn engine(&self) -> OmegaEngine {
        self.engine
    }

    /// Returns the Viterbi log score or MAC expected-accuracy objective.
    #[must_use]
    pub const fn score(&self) -> f64 {
        self.score
    }

    /// Returns the half-open left child core interval.
    #[must_use]
    pub fn left_core(&self) -> core::ops::Range<usize> {
        self.left_core.clone()
    }

    /// Returns the half-open right child core interval.
    #[must_use]
    pub fn right_core(&self) -> core::ops::Range<usize> {
        self.right_core.clone()
    }

    /// Returns the merged profile width after restoring free terminal columns.
    #[must_use]
    pub const fn result_width(&self) -> usize {
        self.result_width
    }

    /// Returns peak scalar DP workspace for this merge.
    #[must_use]
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Returns minimum and maximum match posterior for MAC, if selected.
    #[must_use]
    pub const fn posterior_bounds(&self) -> Option<(f64, f64)> {
        self.posterior_bounds
    }
}

/// Explicit high-water and returned-storage accounting for one Omega call.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ClustalOmegaMemoryEstimate {
    guide_peak_bytes: usize,
    hmm_peak_bytes: usize,
    retained_evidence_bytes: usize,
    returned_alignment_bytes: usize,
    estimated_peak_bytes: usize,
}

impl ClustalOmegaMemoryEstimate {
    /// Returns peak guide-tree construction workspace.
    #[must_use]
    pub const fn guide_peak_bytes(self) -> usize {
        self.guide_peak_bytes
    }

    /// Returns the conservative peak of concurrent HMM-HMM workspaces.
    #[must_use]
    pub const fn hmm_peak_bytes(self) -> usize {
        self.hmm_peak_bytes
    }

    /// Returns complete matrix or mBed metadata retained in the result.
    #[must_use]
    pub const fn retained_evidence_bytes(self) -> usize {
        self.retained_evidence_bytes
    }

    /// Returns expanded aligned cells and identifiers retained in the result.
    #[must_use]
    pub const fn returned_alignment_bytes(self) -> usize {
        self.returned_alignment_bytes
    }

    /// Returns the conservative algorithm-owned peak estimate.
    #[must_use]
    pub const fn estimated_peak_bytes(self) -> usize {
        self.estimated_peak_bytes
    }
}

/// Complete result of a versioned Clustal Omega reconstruction.
#[derive(Clone, Debug)]
pub struct ClustalOmegaAlignment<A: ClustalAlphabet> {
    alignment: MultipleAlignment<A>,
    guide_tree: GuideTree,
    guide_evidence: OmegaGuideEvidence,
    progressive_steps: Box<[OmegaProgressiveStep]>,
    memory: ClustalOmegaMemoryEstimate,
    profile_identifier: &'static str,
}

impl<A: ClustalAlphabet> ClustalOmegaAlignment<A> {
    /// Returns the validated alignment in original input order.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        &self.alignment
    }

    /// Returns the deterministic rooted UPGMA guide tree.
    #[must_use]
    pub const fn guide_tree(&self) -> &GuideTree {
        &self.guide_tree
    }

    /// Returns guide-distance, mBed seed, and leaf-cluster evidence.
    #[must_use]
    pub const fn guide_evidence(&self) -> &OmegaGuideEvidence {
        &self.guide_evidence
    }

    /// Returns progressive merges in guide-tree postorder.
    #[must_use]
    pub const fn progressive_steps(&self) -> &[OmegaProgressiveStep] {
        &self.progressive_steps
    }

    /// Returns checked memory observations and estimates.
    #[must_use]
    pub const fn memory_estimate(&self) -> ClustalOmegaMemoryEstimate {
        self.memory
    }

    /// Returns the stable scientific reconstruction profile identifier.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        self.profile_identifier
    }
}

/// Deterministic Clustal Omega progressive HMM multiple-sequence aligner.
///
/// The type is deterministic and performs no external I/O. [`Self::published`]
/// selects the documented protein 1.0.2 or nucleotide 1.2.4 profile from `A`.
#[derive(Clone, Debug, PartialEq)]
pub struct ClustalOmega<A: ClustalAlphabet> {
    parameters: ClustalOmegaParameters,
    marker: PhantomData<A>,
}

impl<A: ClustalAlphabet> ClustalOmega<A> {
    /// Creates the versioned published profile compatible with `A`.
    #[must_use]
    pub fn published() -> Self {
        let parameters = if is_protein::<A>() {
            ClustalOmegaParameters::protein_v1_0_2()
        } else {
            ClustalOmegaParameters::nucleotide_v1_2_4()
        };
        Self {
            parameters,
            marker: PhantomData,
        }
    }

    /// Creates an aligner from validated molecule-specific parameters.
    ///
    /// # Errors
    ///
    /// Returns [`ClustalOmegaError`] if the molecule profile conflicts with `A`
    /// or a parameter does not satisfy its documented finite domain.
    pub fn with_parameters(parameters: ClustalOmegaParameters) -> Result<Self, ClustalOmegaError> {
        parameters.validate()?;
        let expected = expected_profile::<A>();
        if parameters.molecule_profile() != expected {
            return Err(ClustalOmegaError::MoleculeProfileMismatch {
                expected: expected.diagnostic_name(),
                actual: parameters.molecule_profile().diagnostic_name(),
            });
        }
        Ok(Self {
            parameters,
            marker: PhantomData,
        })
    }

    /// Returns immutable effective parameters.
    #[must_use]
    pub const fn parameters(&self) -> &ClustalOmegaParameters {
        &self.parameters
    }

    /// Evaluates the Wilbur--Lipman one-tuple distance for two sequences.
    ///
    /// # Errors
    ///
    /// Returns [`ClustalOmegaError`] if checked dynamic-programming arithmetic
    /// overflows.
    #[allow(clippy::cast_precision_loss)]
    pub fn pair_distance(
        &self,
        first: &align_core::Sequence<A>,
        second: &align_core::Sequence<A>,
    ) -> Result<f64, ClustalOmegaError> {
        Ok(
            one_tuple_distance(first, second, self.parameters.one_tuple_gap_penalty())? as f64
                / DISTANCE_SCALE as f64,
        )
    }

    /// Aligns a nonempty collection of typed, uniquely named sequences.
    ///
    /// # Errors
    ///
    /// Returns [`ClustalOmegaError`] for invalid collections, checked arithmetic
    /// overflow, invalid derived probabilities, or a failed replay invariant.
    #[allow(clippy::too_many_lines)]
    pub fn align(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<ClustalOmegaAlignment<A>, ClustalOmegaError> {
        self.align_with_phase_timings(records, None)
    }

    /// Aligns records and reports coarse monotonic wall time by pipeline phase.
    ///
    /// This diagnostic entry point returns the same scientific result as
    /// [`Self::align`]. Timer reads occur only on this opt-in path. Reported
    /// phases exclude FASTA parsing and result serialization performed by a
    /// caller.
    ///
    /// # Errors
    ///
    /// Returns the same [`ClustalOmegaError`] values as [`Self::align`].
    pub fn align_profiled(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<Profiled<ClustalOmegaAlignment<A>>, ClustalOmegaError> {
        let mut timings = PhaseTimingsBuilder::new();
        let result = self.align_with_phase_timings(records, Some(&mut timings))?;
        Ok(Profiled::new(result, timings.finish()))
    }

    #[allow(clippy::too_many_lines)]
    fn align_with_phase_timings(
        &self,
        records: &[SequenceRecord<A>],
        mut timings: Option<&mut PhaseTimingsBuilder>,
    ) -> Result<ClustalOmegaAlignment<A>, ClustalOmegaError> {
        measure_phase(&mut timings, "input_validation", || {
            self.parameters.validate()?;
            validate_records(records)
        })?;
        let sequence_count = records.len();
        let leaf_indices = (0..sequence_count).collect::<Vec<_>>();

        let (fixed_root, guide_tree, guide_evidence, guide_peak_bytes) = measure_phase(
            &mut timings,
            "guide_tree",
            || -> Result<_, ClustalOmegaError> {
                let (fixed_root, guide_evidence, guide_peak_bytes) =
                    match self.parameters.guide_tree_mode() {
                        OmegaGuideTreeMode::Full => {
                            let (matrix, distance_workspace) = full_distance_matrix(
                                records,
                                self.parameters.one_tuple_gap_penalty(),
                                self.parameters.thread_count(),
                            )?;
                            let tree = upgma(&matrix, &leaf_indices)?;
                            let matrix_bytes = matrix.memory_bytes()?;
                            let distances = matrix.to_public()?;
                            let guide_peak = matrix_bytes
                                .checked_add(tree.workspace_bytes.max(distance_workspace))
                                .ok_or(ClustalOmegaError::ArithmeticOverflow {
                                    stage: "full guide-tree peak accounting",
                                })?;
                            (
                                tree.root,
                                OmegaGuideEvidence {
                                    method: OmegaGuideTreeMethod::Full,
                                    distances: Some(distances),
                                    seed_indices: Box::new([]),
                                    clusters: Box::new([]),
                                },
                                guide_peak,
                            )
                        }
                        OmegaGuideTreeMode::Mbed
                            if sequence_count <= self.parameters.maximum_cluster_size() =>
                        {
                            let (matrix, distance_workspace) = full_distance_matrix(
                                records,
                                self.parameters.one_tuple_gap_penalty(),
                                self.parameters.thread_count(),
                            )?;
                            let tree = upgma(&matrix, &leaf_indices)?;
                            let matrix_bytes = matrix.memory_bytes()?;
                            let distances = matrix.to_public()?;
                            let guide_peak = matrix_bytes
                                .checked_add(tree.workspace_bytes.max(distance_workspace))
                                .ok_or(ClustalOmegaError::ArithmeticOverflow {
                                    stage: "effective-full guide-tree peak accounting",
                                })?;
                            (
                                tree.root,
                                OmegaGuideEvidence {
                                    method: OmegaGuideTreeMethod::MbedEffectiveFull,
                                    distances: Some(distances),
                                    seed_indices: Box::new([]),
                                    clusters: Box::new([]),
                                },
                                guide_peak,
                            )
                        }
                        OmegaGuideTreeMode::Mbed => {
                            let result = mbed(
                                records,
                                self.parameters.maximum_cluster_size(),
                                self.parameters.kmeans_seed(),
                                self.parameters.one_tuple_gap_penalty(),
                            )?;
                            let evidence = OmegaGuideEvidence {
                                method: OmegaGuideTreeMethod::Mbed,
                                distances: None,
                                seed_indices: result.seed_indices.into_boxed_slice(),
                                clusters: result
                                    .clusters
                                    .into_iter()
                                    .map(Vec::into_boxed_slice)
                                    .collect::<Vec<_>>()
                                    .into_boxed_slice(),
                            };
                            (result.root, evidence, result.workspace_bytes)
                        }
                    };
                let guide_tree = to_public_tree(fixed_root.clone(), sequence_count)?;
                Ok((fixed_root, guide_tree, guide_evidence, guide_peak_bytes))
            },
        )?;
        let subtree = measure_phase(&mut timings, "progressive_alignment", || {
            align_tree_node(
                &fixed_root,
                records,
                &self.parameters,
                is_protein::<A>(),
                self.parameters.thread_count(),
            )
        })?;
        measure_phase(&mut timings, "finalization", || {
            let hmm_peak_bytes = subtree.peak_workspace;
            let alignment = finalize_alignment(subtree.profile, records)?;
            let returned_alignment_bytes = estimate_alignment_bytes(&alignment)?;
            let retained_evidence_bytes = estimate_evidence_bytes(&guide_evidence)?;
            let estimated_peak_bytes = retained_evidence_bytes
                .checked_add(returned_alignment_bytes)
                .and_then(|resident| resident.checked_add(guide_peak_bytes.max(hmm_peak_bytes)))
                .ok_or(ClustalOmegaError::ArithmeticOverflow {
                    stage: "Clustal Omega peak memory accounting",
                })?;
            let memory = ClustalOmegaMemoryEstimate {
                guide_peak_bytes,
                hmm_peak_bytes,
                retained_evidence_bytes,
                returned_alignment_bytes,
                estimated_peak_bytes,
            };
            Ok(ClustalOmegaAlignment {
                alignment,
                guide_tree,
                guide_evidence,
                progressive_steps: subtree.steps.into_boxed_slice(),
                memory,
                profile_identifier: self.parameters.molecule_profile().identifier(),
            })
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

struct AlignedSubtree<A: ClustalAlphabet> {
    profile: WorkingProfile<A>,
    steps: Vec<OmegaProgressiveStep>,
    peak_workspace: usize,
}

fn align_tree_node<A: ClustalAlphabet>(
    node: &FixedTreeNode,
    records: &[SequenceRecord<A>],
    parameters: &ClustalOmegaParameters,
    is_protein: bool,
    thread_count: usize,
) -> Result<AlignedSubtree<A>, ClustalOmegaError> {
    match node {
        FixedTreeNode::Leaf { sequence_index } => Ok(AlignedSubtree {
            profile: WorkingProfile::single(&records[*sequence_index], *sequence_index),
            steps: Vec::new(),
            peak_workspace: 0,
        }),
        FixedTreeNode::Branch { left, right, .. } => {
            let run_parallel = thread_count > 1
                && matches!(left.as_ref(), FixedTreeNode::Branch { .. })
                && matches!(right.as_ref(), FixedTreeNode::Branch { .. });
            let (mut left, mut right) = if run_parallel {
                let left_threads = thread_count / 2;
                let right_threads = thread_count - left_threads;
                std::thread::scope(|scope| {
                    let left_handle = scope.spawn(|| {
                        align_tree_node(left, records, parameters, is_protein, left_threads)
                    });
                    let right =
                        align_tree_node(right, records, parameters, is_protein, right_threads);
                    let left =
                        left_handle
                            .join()
                            .map_err(|_| ClustalOmegaError::InternalInvariant {
                                message: "a deterministic progressive-alignment worker panicked",
                            })?;
                    Ok::<_, ClustalOmegaError>((left?, right?))
                })?
            } else {
                (
                    align_tree_node(left, records, parameters, is_protein, thread_count)?,
                    align_tree_node(right, records, parameters, is_protein, thread_count)?,
                )
            };
            let child_peak = if run_parallel {
                left.peak_workspace
                    .checked_add(right.peak_workspace)
                    .ok_or(ClustalOmegaError::ArithmeticOverflow {
                        stage: "parallel progressive workspace accounting",
                    })?
            } else {
                left.peak_workspace.max(right.peak_workspace)
            };
            let left_members = left.profile.members();
            let right_members = right.profile.members();
            let result = align_profiles(left.profile, right.profile, parameters, is_protein)?;
            left.steps.append(&mut right.steps);
            left.steps.push(OmegaProgressiveStep {
                left_members: left_members.into_boxed_slice(),
                right_members: right_members.into_boxed_slice(),
                engine: result.observation.engine,
                score: result.observation.score,
                left_core: result.observation.left_core,
                right_core: result.observation.right_core,
                result_width: result.profile.width(),
                workspace_bytes: result.observation.workspace_bytes,
                posterior_bounds: result.observation.posterior_bounds,
            });
            debug_assert_eq!(result.observation.operations.len(), result.profile.width());
            Ok(AlignedSubtree {
                profile: result.profile,
                steps: left.steps,
                peak_workspace: child_peak.max(result.observation.workspace_bytes),
            })
        }
    }
}

fn finalize_alignment<A: ClustalAlphabet>(
    mut profile: WorkingProfile<A>,
    records: &[SequenceRecord<A>],
) -> Result<MultipleAlignment<A>, ClustalOmegaError> {
    profile.rows.sort_by_key(|row| row.sequence_index);
    if profile.rows.len() != records.len() {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "root profile row count differs from input count",
        });
    }
    let mut rows = Vec::with_capacity(records.len());
    for (expected_index, row) in profile.rows.into_iter().enumerate() {
        if row.sequence_index != expected_index {
            return Err(ClustalOmegaError::InternalInvariant {
                message: "root profile does not contain every input index exactly once",
            });
        }
        let ungapped = row
            .cells
            .iter()
            .filter_map(|cell| *cell)
            .collect::<Vec<_>>();
        if ungapped != records[expected_index].sequence().as_slice() {
            return Err(ClustalOmegaError::InternalInvariant {
                message: "root profile does not replay to its original sequence",
            });
        }
        rows.push(AlignedRow::new(
            records[expected_index].identifier(),
            row.cells,
        )?);
    }
    Ok(MultipleAlignment::try_new(rows)?)
}

fn validate_records<A: Alphabet>(records: &[SequenceRecord<A>]) -> Result<(), ClustalOmegaError> {
    if records.is_empty() {
        return Err(ClustalOmegaError::Msa(MsaError::EmptyCollection));
    }
    let mut identifiers = HashSet::with_capacity(records.len());
    for record in records {
        if !identifiers.insert(record.identifier()) {
            return Err(ClustalOmegaError::Msa(MsaError::DuplicateIdentifier {
                identifier: record.identifier().to_owned(),
            }));
        }
    }
    Ok(())
}

fn estimate_alignment_bytes<A: Alphabet>(
    alignment: &MultipleAlignment<A>,
) -> Result<usize, ClustalOmegaError> {
    let cells = alignment
        .sequence_count()
        .checked_mul(alignment.width())
        .and_then(|cells| cells.checked_mul(size_of::<Option<Symbol<A>>>()))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "returned Omega alignment memory accounting",
        })?;
    alignment.rows().iter().try_fold(cells, |bytes, row| {
        bytes
            .checked_add(row.identifier().len())
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "returned Omega alignment memory accounting",
            })
    })
}

fn estimate_evidence_bytes(evidence: &OmegaGuideEvidence) -> Result<usize, ClustalOmegaError> {
    let distance_bytes = evidence.distances.as_ref().map_or(Ok(0), |matrix| {
        matrix
            .as_flat_slice()
            .len()
            .checked_mul(size_of::<f64>())
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "retained Omega distance evidence accounting",
            })
    })?;
    let seed_bytes = evidence
        .seed_indices
        .len()
        .checked_mul(size_of::<usize>())
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "retained mBed seed accounting",
        })?;
    let cluster_members = evidence
        .clusters
        .iter()
        .try_fold(0_usize, |sum, cluster| sum.checked_add(cluster.len()))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "retained mBed cluster accounting",
        })?;
    distance_bytes
        .checked_add(seed_bytes)
        .and_then(|bytes| bytes.checked_add(cluster_members * size_of::<usize>()))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "retained Omega guide evidence accounting",
        })
}

fn expected_profile<A: Alphabet>() -> OmegaMoleculeProfile {
    if is_protein::<A>() {
        OmegaMoleculeProfile::ProteinV1_0_2
    } else {
        OmegaMoleculeProfile::NucleotideV1_2_4
    }
}

fn is_protein<A: Alphabet>() -> bool {
    TypeId::of::<A>() == TypeId::of::<Protein>()
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence};
    use align_msa::SequenceRecord;

    use super::{ClustalOmega, OmegaGuideTreeMethod};
    use crate::{ClustalOmegaParameters, OmegaEngineMode, OmegaGuideTreeMode};

    fn records(values: &[(&str, &str)]) -> Vec<SequenceRecord<Dna>> {
        values
            .iter()
            .map(|(identifier, sequence)| {
                SequenceRecord::new(
                    *identifier,
                    DnaSequence::from_ascii(sequence.as_bytes()).unwrap(),
                )
                .unwrap()
            })
            .collect()
    }

    #[test]
    fn complete_pipeline_preserves_rows_and_exposes_every_merge() {
        let records = records(&[("a", "ACGT"), ("b", "AGT"), ("c", "ACGTT")]);
        let parameters = ClustalOmegaParameters::nucleotide_v1_2_4()
            .with_engine_mode(OmegaEngineMode::Viterbi)
            .with_guide_tree_mode(OmegaGuideTreeMode::Full);
        let result = ClustalOmega::<Dna>::with_parameters(parameters)
            .unwrap()
            .align(&records)
            .unwrap();
        assert_eq!(result.progressive_steps().len(), 2);
        assert_eq!(result.guide_tree().leaf_count(), 3);
        assert_eq!(result.guide_evidence().method(), OmegaGuideTreeMethod::Full);
        assert!(result.guide_evidence().distances().is_some());
        for (row, record) in result.alignment().rows().iter().zip(&records) {
            assert_eq!(row.to_ungapped_ascii(), record.sequence().to_ascii_string());
        }
        assert!(result.memory_estimate().estimated_peak_bytes() > 0);
    }
}
