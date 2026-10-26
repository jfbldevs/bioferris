//! Fluent entry points for multiple-sequence alignment.

use align_clustal::{ClustalAlphabet, ClustalOmega, ClustalOmegaAlignment};
use align_core::{Alphabet, Protein};
use align_mafft::{
    EINSI_PROFILE_IDENTIFIER, EinsiAlignment, FFT_NS1_PROFILE_IDENTIFIER,
    FFT_NS2_PROFILE_IDENTIFIER, FFT_NSI_PROFILE_IDENTIFIER, FftNs1Alignment, FftNs2Alignment,
    FftNsiAlignment, GINSI_PROFILE_IDENTIFIER, GinsiAlignment, LINSI_PROFILE_IDENTIFIER,
    LinsiAlignment, MafftAlphabet, MafftEinsi, MafftFftNs1, MafftFftNs2, MafftFftNsi, MafftGinsi,
    MafftLinsi,
};
use align_mafft_compat::{
    MAFFT_7526_NATIVE_VERSION, MAFFT_7526_PORT_VERSION, Mafft7526Alignment, Mafft7526Einsi,
    Mafft7526Linsi,
};
use align_msa::{MultipleAlignment, SequenceRecord};
use align_muscle::{Muscle5, Muscle5Alignment, Muscle5Alphabet};

use crate::graph::PoastaAligner;
use crate::report::{MultipleRun, MultipleRunInfo, ReportField};
use crate::result::MultipleAlignmentView;
use crate::{Error, Result};

pub use align_muscle::Muscle5Mode;

/// Discoverable constructors for the three MSA families.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Multiple;

impl Multiple {
    /// Selects the molecule-specific published Clustal Omega profile.
    #[must_use]
    pub const fn clustal_omega() -> ClustalOmegaAligner {
        ClustalOmegaAligner { thread_count: None }
    }

    /// Selects the MUSCLE 5 publication-default `Align` pipeline.
    #[must_use]
    pub const fn muscle5() -> Muscle5Aligner {
        Muscle5Aligner {
            mode: Muscle5Mode::Align,
        }
    }

    /// Selects the MUSCLE 5 publication-default `Super5` pipeline.
    #[must_use]
    pub const fn muscle5_super5() -> Muscle5Aligner {
        Muscle5Aligner {
            mode: Muscle5Mode::Super5,
        }
    }

    /// Selects one exact retained MAFFT scientific profile.
    #[must_use]
    pub const fn mafft(strategy: MafftStrategy) -> MafftAligner {
        MafftAligner {
            strategy,
            thread_count: None,
        }
    }

    /// Selects native-output-compatible MAFFT 7.526 L-INS-i.
    ///
    /// This path is intentionally separate from [`Self::mafft`], whose six
    /// profiles remain independent publication reconstructions.
    #[must_use]
    pub const fn mafft_7526_linsi() -> Mafft7526LinsiAligner {
        Mafft7526LinsiAligner { thread_count: None }
    }

    /// Selects native-output-compatible MAFFT 7.526 protein E-INS-i.
    ///
    /// DNA and RNA are excluded until the delegated generalized-affine route
    /// is independently byte-compatible with native MAFFT on adversarial cases.
    #[must_use]
    pub const fn mafft_7526_einsi() -> Mafft7526EinsiAligner {
        Mafft7526EinsiAligner { thread_count: None }
    }
}

/// An algorithm facade that accepts a collection of named sequences.
pub trait RecordsAligner<A: Alphabet> {
    /// Native result retained by this algorithm family.
    type Output;

    /// Aligns `records` without changing the selected scientific profile.
    ///
    /// # Errors
    ///
    /// Returns the structured facade error produced by validation or execution.
    fn align_records(&self, records: &[SequenceRecord<A>]) -> Result<Self::Output>;
}

/// Adds Polars-style `.align_with(...)` chaining to record collections.
pub trait AlignRecordsExt<A: Alphabet> {
    /// Aligns this collection through `aligner`.
    ///
    /// # Errors
    ///
    /// Returns the structured facade error produced by validation or execution.
    fn align_with<M>(&self, aligner: M) -> Result<M::Output>
    where
        M: RecordsAligner<A>;
}

impl<A: Alphabet> AlignRecordsExt<A> for [SequenceRecord<A>] {
    fn align_with<M>(&self, aligner: M) -> Result<M::Output>
    where
        M: RecordsAligner<A>,
    {
        aligner.align_records(self)
    }
}

impl<A: Alphabet> AlignRecordsExt<A> for Vec<SequenceRecord<A>> {
    fn align_with<M>(&self, aligner: M) -> Result<M::Output>
    where
        M: RecordsAligner<A>,
    {
        self.as_slice().align_with(aligner)
    }
}

/// Reusable Clustal Omega configuration.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct ClustalOmegaAligner {
    thread_count: Option<usize>,
}

impl ClustalOmegaAligner {
    /// Replaces the deterministic worker count.
    #[must_use]
    pub const fn threads(mut self, count: usize) -> Self {
        self.thread_count = Some(count);
        self
    }

    /// Aligns a nonempty collection with the published molecule-specific profile.
    ///
    /// # Errors
    ///
    /// Returns [`Error::ClustalOmega`] for invalid records, worker count,
    /// dimensions, arithmetic, allocation, or replay failure.
    pub fn align<A: ClustalAlphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<ClustalOmegaAlignment<A>> {
        let published = ClustalOmega::<A>::published();
        if let Some(thread_count) = self.thread_count {
            let parameters = published
                .parameters()
                .clone()
                .with_thread_count(thread_count)?;
            ClustalOmega::with_parameters(parameters)?
                .align(records)
                .map_err(Error::from)
        } else {
            published.align(records).map_err(Error::from)
        }
    }

    /// Aligns records and retains profile and execution evidence for a readable
    /// report.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn run<A: ClustalAlphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<MultipleRun<A, ClustalOmegaAlignment<A>>> {
        let result = self.align(records)?;
        let memory = result.memory_estimate();
        let info = MultipleRunInfo::new(
            "Clustal Omega",
            result.profile_identifier(),
            Some("published progressive HMM"),
        )
        .with_details(vec![
            ReportField::new(
                "Guide method",
                format!("{:?}", result.guide_evidence().method()),
            ),
            ReportField::new(
                "Progressive merges",
                result.progressive_steps().len().to_string(),
            ),
            ReportField::new(
                "Guide seeds",
                result.guide_evidence().seed_indices().len().to_string(),
            ),
            ReportField::new(
                "Guide clusters",
                result.guide_evidence().clusters().len().to_string(),
            ),
            ReportField::new(
                "Estimated peak bytes",
                memory.estimated_peak_bytes().to_string(),
            ),
        ]);
        Ok(MultipleRun::new(result, info))
    }
}

impl<A: ClustalAlphabet> RecordsAligner<A> for ClustalOmegaAligner {
    type Output = ClustalOmegaAlignment<A>;

    fn align_records(&self, records: &[SequenceRecord<A>]) -> Result<Self::Output> {
        self.align(records)
    }
}

/// Reusable MUSCLE 5 mode selection.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Muscle5Aligner {
    mode: Muscle5Mode,
}

impl Muscle5Aligner {
    /// Selects the publication-default `Align` pipeline.
    #[must_use]
    pub const fn align_mode(mut self) -> Self {
        self.mode = Muscle5Mode::Align;
        self
    }

    /// Selects the publication-default `Super5` pipeline.
    #[must_use]
    pub const fn super5(mut self) -> Self {
        self.mode = Muscle5Mode::Super5;
        self
    }

    /// Returns the configured MUSCLE 5 pipeline.
    #[must_use]
    pub const fn mode(self) -> Muscle5Mode {
        self.mode
    }

    /// Aligns a nonempty collection with MUSCLE 5.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Muscle5`] for invalid input, arithmetic, allocation, or
    /// replay failure.
    pub fn align<A: Muscle5Alphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<Muscle5Alignment<A>> {
        match self.mode {
            Muscle5Mode::Align => Muscle5::<A>::published_align(),
            Muscle5Mode::Super5 => Muscle5::<A>::published_super5(),
        }
        .align(records)
        .map_err(Error::from)
    }

    /// Aligns records and retains MUSCLE mode, refinement, objective, and memory
    /// evidence for presentation.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn run<A: Muscle5Alphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<MultipleRun<A, Muscle5Alignment<A>>> {
        let result = self.align(records)?;
        let mode = muscle_mode_name(result.mode());
        let info = MultipleRunInfo::new(
            "MUSCLE 5",
            "MUSCLE-5.1.0-Nature-2022-default-v2",
            Some(mode),
        )
        .with_details(vec![
            ReportField::new("Perturbation seed", result.perturbation_seed().to_string()),
            ReportField::new(
                "Guide permutation",
                format!("{:?}", result.guide_tree_permutation()),
            ),
            ReportField::new(
                "Progressive merges",
                result.progressive_steps().len().to_string(),
            ),
            ReportField::new(
                "Refinement attempts",
                result.refinement_attempts().len().to_string(),
            ),
            ReportField::new(
                "Final objective",
                format!("{:.9}", result.final_objective_score()),
            ),
            ReportField::new(
                "Objective pairs",
                format!(
                    "{}/{}",
                    result.evaluated_objective_pair_count(),
                    result.possible_objective_pair_count()
                ),
            ),
            ReportField::new(
                "Estimated peak bytes",
                result.memory_estimate().estimated_peak_bytes().to_string(),
            ),
        ]);
        Ok(MultipleRun::new(result, info))
    }
}

impl Default for Muscle5Aligner {
    fn default() -> Self {
        Multiple::muscle5()
    }
}

impl<A: Muscle5Alphabet> RecordsAligner<A> for Muscle5Aligner {
    type Output = Muscle5Alignment<A>;

    fn align_records(&self, records: &[SequenceRecord<A>]) -> Result<Self::Output> {
        self.align(records)
    }
}

/// Reusable selector for native-output-compatible MAFFT 7.526 L-INS-i.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Mafft7526LinsiAligner {
    thread_count: Option<usize>,
}

impl Mafft7526LinsiAligner {
    /// Replaces the deterministic Rayon worker count.
    #[must_use]
    pub const fn threads(mut self, count: usize) -> Self {
        self.thread_count = Some(count);
        self
    }

    /// Aligns protein, DNA, or RNA records with MAFFT 7.526 L-INS-i.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Mafft`] for invalid input, worker-pool failure, or a
    /// source-changing result from the pinned compatibility engine.
    pub fn align<A: MafftAlphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<Mafft7526Alignment<A>> {
        let mut aligner = Mafft7526Linsi::<A>::new();
        if let Some(thread_count) = self.thread_count {
            aligner = aligner.with_thread_count(thread_count)?;
        }
        aligner.align(records).map_err(Error::from)
    }

    /// Aligns records and retains compatibility provenance for reporting.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn run<A: MafftAlphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<MultipleRun<A, Mafft7526Alignment<A>>> {
        let result = self.align(records)?;
        let info = MultipleRunInfo::new(
            "MAFFT 7.526",
            result.profile_identifier(),
            Some("L-INS-i native-output-compatible"),
        )
        .with_details(mafft_7526_details(&result));
        Ok(MultipleRun::new(result, info))
    }
}

impl<A: MafftAlphabet> RecordsAligner<A> for Mafft7526LinsiAligner {
    type Output = Mafft7526Alignment<A>;

    fn align_records(&self, records: &[SequenceRecord<A>]) -> Result<Self::Output> {
        self.align(records)
    }
}

/// Reusable selector for native-output-compatible MAFFT 7.526 protein E-INS-i.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Mafft7526EinsiAligner {
    thread_count: Option<usize>,
}

impl Mafft7526EinsiAligner {
    /// Replaces the deterministic Rayon worker count.
    #[must_use]
    pub const fn threads(mut self, count: usize) -> Self {
        self.thread_count = Some(count);
        self
    }

    /// Aligns protein records with MAFFT 7.526 E-INS-i.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Mafft`] for invalid input, worker-pool failure, or a
    /// source-changing result from the pinned compatibility engine.
    pub fn align(
        &self,
        records: &[SequenceRecord<Protein>],
    ) -> Result<Mafft7526Alignment<Protein>> {
        let mut aligner = Mafft7526Einsi::new();
        if let Some(thread_count) = self.thread_count {
            aligner = aligner.with_thread_count(thread_count)?;
        }
        aligner.align(records).map_err(Error::from)
    }

    /// Aligns protein records and retains compatibility provenance for reporting.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn run(
        &self,
        records: &[SequenceRecord<Protein>],
    ) -> Result<MultipleRun<Protein, Mafft7526Alignment<Protein>>> {
        let result = self.align(records)?;
        let info = MultipleRunInfo::new(
            "MAFFT 7.526",
            result.profile_identifier(),
            Some("E-INS-i protein native-output-compatible"),
        )
        .with_details(mafft_7526_details(&result));
        Ok(MultipleRun::new(result, info))
    }
}

impl RecordsAligner<Protein> for Mafft7526EinsiAligner {
    type Output = Mafft7526Alignment<Protein>;

    fn align_records(&self, records: &[SequenceRecord<Protein>]) -> Result<Self::Output> {
        self.align(records)
    }
}

fn mafft_7526_details<A: MafftAlphabet>(result: &Mafft7526Alignment<A>) -> Vec<ReportField> {
    vec![
        ReportField::new("Native target", MAFFT_7526_NATIVE_VERSION),
        ReportField::new("Rust port", format!("mafft {MAFFT_7526_PORT_VERSION}")),
        ReportField::new(
            "Progressive steps",
            result.progressive_step_count().to_string(),
        ),
        ReportField::new("Engine score", format!("{:.9}", result.engine_score())),
        ReportField::new("Workers", result.thread_count().to_string()),
    ]
}

/// One of the six retained MAFFT scientific profiles.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MafftStrategy {
    /// Single progressive FFT-NS-1 pass.
    FftNs1,
    /// Two progressive FFT-NS-2 passes.
    FftNs2,
    /// Iterative FFT-NS-i refinement.
    FftNsi,
    /// Global-pair G-INS-i refinement.
    Ginsi,
    /// Local-pair L-INS-i refinement.
    Linsi,
    /// Generalized-affine E-INS-i refinement.
    Einsi,
}

impl MafftStrategy {
    /// Returns the immutable scientific-profile identifier.
    #[must_use]
    pub const fn profile_identifier(self) -> &'static str {
        match self {
            Self::FftNs1 => FFT_NS1_PROFILE_IDENTIFIER,
            Self::FftNs2 => FFT_NS2_PROFILE_IDENTIFIER,
            Self::FftNsi => FFT_NSI_PROFILE_IDENTIFIER,
            Self::Ginsi => GINSI_PROFILE_IDENTIFIER,
            Self::Linsi => LINSI_PROFILE_IDENTIFIER,
            Self::Einsi => EINSI_PROFILE_IDENTIFIER,
        }
    }
}

/// Reusable selector for one retained MAFFT profile.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MafftAligner {
    strategy: MafftStrategy,
    thread_count: Option<usize>,
}

impl MafftAligner {
    /// Returns the selected strategy.
    #[must_use]
    pub const fn strategy(self) -> MafftStrategy {
        self.strategy
    }

    /// Returns the immutable scientific-profile identifier.
    #[must_use]
    pub const fn profile_identifier(self) -> &'static str {
        self.strategy.profile_identifier()
    }

    /// Replaces the deterministic refinement worker count for FFT-NS-i or L-INS-i.
    ///
    /// The other retained MAFFT profiles are currently scalar and report a
    /// structured unsupported-option error instead of silently ignoring this.
    #[must_use]
    pub const fn threads(mut self, count: usize) -> Self {
        self.thread_count = Some(count);
        self
    }

    /// Aligns a nonempty collection with the selected MAFFT profile.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Mafft`] for invalid records, dimensions, arithmetic,
    /// allocation, numeric, or replay failure.
    pub fn align<A: MafftAlphabet>(&self, records: &[SequenceRecord<A>]) -> Result<MafftResult<A>> {
        if self.thread_count.is_some()
            && !matches!(self.strategy, MafftStrategy::FftNsi | MafftStrategy::Linsi)
        {
            return Err(Error::UnsupportedOption {
                algorithm: self.strategy.profile_identifier(),
                option: "threads",
            });
        }
        let result = match self.strategy {
            MafftStrategy::FftNs1 => MafftResult::FftNs1(MafftFftNs1::<A>::new().align(records)?),
            MafftStrategy::FftNs2 => MafftResult::FftNs2(MafftFftNs2::<A>::new().align(records)?),
            MafftStrategy::FftNsi => {
                let mut aligner = MafftFftNsi::<A>::new();
                if let Some(thread_count) = self.thread_count {
                    aligner = aligner.with_thread_count(thread_count)?;
                }
                MafftResult::FftNsi(aligner.align(records)?)
            }
            MafftStrategy::Ginsi => MafftResult::Ginsi(MafftGinsi::<A>::new().align(records)?),
            MafftStrategy::Linsi => {
                let mut aligner = MafftLinsi::<A>::new();
                if let Some(thread_count) = self.thread_count {
                    aligner = aligner.with_thread_count(thread_count)?;
                }
                MafftResult::Linsi(aligner.align(records)?)
            }
            MafftStrategy::Einsi => MafftResult::Einsi(MafftEinsi::<A>::new().align(records)?),
        };
        Ok(result)
    }

    /// Aligns records and retains exact MAFFT strategy and refinement evidence for
    /// presentation.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn run<A: MafftAlphabet>(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<MultipleRun<A, MafftResult<A>>> {
        let result = self.align(records)?;
        let details = result.report_details();
        let info = MultipleRunInfo::new(
            "MAFFT",
            result.profile_identifier(),
            Some(mafft_strategy_name(result.strategy())),
        )
        .with_details(details);
        Ok(MultipleRun::new(result, info))
    }
}

impl<A: MafftAlphabet> RecordsAligner<A> for MafftAligner {
    type Output = MafftResult<A>;

    fn align_records(&self, records: &[SequenceRecord<A>]) -> Result<Self::Output> {
        self.align(records)
    }
}

impl<A: Alphabet> RecordsAligner<A> for PoastaAligner {
    type Output = align_poa::PoastaPoaResult<A>;

    fn align_records(&self, records: &[SequenceRecord<A>]) -> Result<Self::Output> {
        self.align(records)
    }
}

/// Result and scientific evidence returned by one selected MAFFT profile.
#[derive(Clone, Debug)]
pub enum MafftResult<A: MafftAlphabet> {
    /// FFT-NS-1 result.
    FftNs1(FftNs1Alignment<A>),
    /// FFT-NS-2 result.
    FftNs2(FftNs2Alignment<A>),
    /// FFT-NS-i result.
    FftNsi(FftNsiAlignment<A>),
    /// G-INS-i result.
    Ginsi(GinsiAlignment<A>),
    /// L-INS-i result.
    Linsi(LinsiAlignment<A>),
    /// E-INS-i result.
    Einsi(EinsiAlignment<A>),
}

impl<A: MafftAlphabet> MafftResult<A> {
    /// Returns the selected strategy.
    #[must_use]
    pub const fn strategy(&self) -> MafftStrategy {
        match self {
            Self::FftNs1(_) => MafftStrategy::FftNs1,
            Self::FftNs2(_) => MafftStrategy::FftNs2,
            Self::FftNsi(_) => MafftStrategy::FftNsi,
            Self::Ginsi(_) => MafftStrategy::Ginsi,
            Self::Linsi(_) => MafftStrategy::Linsi,
            Self::Einsi(_) => MafftStrategy::Einsi,
        }
    }

    /// Returns the immutable scientific-profile identifier.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        self.strategy().profile_identifier()
    }

    /// Returns the canonical rectangular multiple alignment.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        match self {
            Self::FftNs1(result) => result.alignment(),
            Self::FftNs2(result) => result.alignment(),
            Self::FftNsi(result) => result.alignment(),
            Self::Ginsi(result) => result.alignment(),
            Self::Linsi(result) => result.alignment(),
            Self::Einsi(result) => result.alignment(),
        }
    }

    /// Consumes the profile evidence and returns the canonical alignment.
    #[must_use]
    pub fn into_alignment(self) -> MultipleAlignment<A> {
        match self {
            Self::FftNs1(result) => result.into_alignment(),
            Self::FftNs2(result) => result.into_alignment(),
            Self::FftNsi(result) => result.into_alignment(),
            Self::Ginsi(result) => result.into_alignment(),
            Self::Linsi(result) => result.into_alignment(),
            Self::Einsi(result) => result.into_alignment(),
        }
    }

    fn report_details(&self) -> Vec<ReportField> {
        match self {
            Self::FftNs1(result) => vec![
                ReportField::new(
                    "Effective tuple",
                    result.effective_tuple_length().to_string(),
                ),
                ReportField::new("Progressive merges", result.merges().len().to_string()),
                ReportField::new("Root merge score", optional_f64(result.root_merge_score())),
                ReportField::new(
                    "Estimated peak bytes",
                    result.memory_estimate().estimated_peak_bytes().to_string(),
                ),
            ],
            Self::FftNs2(result) => vec![
                ReportField::new(
                    "Second-pass merges",
                    result.second_merges().len().to_string(),
                ),
                ReportField::new(
                    "Second root score",
                    optional_f64(result.second_root_merge_score()),
                ),
                ReportField::new(
                    "Guide tree changed",
                    result.guide_tree_changed().to_string(),
                ),
                ReportField::new(
                    "Estimated peak bytes",
                    result.memory_estimate().estimated_peak_bytes().to_string(),
                ),
            ],
            Self::FftNsi(result) => vec![
                ReportField::new("Refinement cycles", result.cycles().len().to_string()),
                ReportField::new("Accepted cycles", result.accepted_cycle_count().to_string()),
                ReportField::new("Termination", format!("{:?}", result.termination())),
                ReportField::new("Final score", format!("{:.9}", result.final_score())),
                ReportField::new(
                    "Changed initial MSA",
                    result.changed_initial_alignment().to_string(),
                ),
                ReportField::new(
                    "Estimated peak bytes",
                    result.memory_estimate().estimated_peak_bytes().to_string(),
                ),
            ],
            Self::Ginsi(result) => iterative_mafft_details(
                result.pair_evidence().len(),
                result.progressive_merges().len(),
                result.cycles().len(),
                result.accepted_cycle_count(),
                format!("{:?}", result.termination()),
                result.final_score(),
                result.memory_estimate().estimated_peak_bytes(),
            ),
            Self::Linsi(result) => iterative_mafft_details(
                result.pair_evidence().len(),
                result.progressive_merges().len(),
                result.cycles().len(),
                result.accepted_cycle_count(),
                format!("{:?}", result.termination()),
                result.final_score(),
                result.memory_estimate().estimated_peak_bytes(),
            ),
            Self::Einsi(result) => iterative_mafft_details(
                result.pair_evidence().len(),
                result.progressive_merges().len(),
                result.cycles().len(),
                result.accepted_cycle_count(),
                format!("{:?}", result.termination()),
                result.final_score(),
                result.memory_estimate().estimated_peak_bytes(),
            ),
        }
    }
}

fn iterative_mafft_details(
    pair_evidence: usize,
    progressive_merges: usize,
    cycles: usize,
    accepted_cycles: usize,
    termination: String,
    final_score: f64,
    estimated_peak_bytes: usize,
) -> Vec<ReportField> {
    vec![
        ReportField::new("Pair evidence", pair_evidence.to_string()),
        ReportField::new("Progressive merges", progressive_merges.to_string()),
        ReportField::new("Refinement cycles", cycles.to_string()),
        ReportField::new("Accepted cycles", accepted_cycles.to_string()),
        ReportField::new("Termination", termination),
        ReportField::new("Final score", format!("{final_score:.9}")),
        ReportField::new("Estimated peak bytes", estimated_peak_bytes.to_string()),
    ]
}

fn optional_f64(value: Option<f64>) -> String {
    value.map_or_else(|| String::from("-"), |value| format!("{value:.9}"))
}

const fn muscle_mode_name(mode: Muscle5Mode) -> &'static str {
    match mode {
        Muscle5Mode::Align => "Align",
        Muscle5Mode::Super5 => "Super5",
    }
}

const fn mafft_strategy_name(strategy: MafftStrategy) -> &'static str {
    match strategy {
        MafftStrategy::FftNs1 => "FFT-NS-1",
        MafftStrategy::FftNs2 => "FFT-NS-2",
        MafftStrategy::FftNsi => "FFT-NS-i",
        MafftStrategy::Ginsi => "G-INS-i",
        MafftStrategy::Linsi => "L-INS-i",
        MafftStrategy::Einsi => "E-INS-i",
    }
}

impl<A: ClustalAlphabet> MultipleAlignmentView for ClustalOmegaAlignment<A> {
    type Alphabet = A;

    fn alignment(&self) -> &MultipleAlignment<A> {
        self.alignment()
    }
}

impl<A: Muscle5Alphabet> MultipleAlignmentView for Muscle5Alignment<A> {
    type Alphabet = A;

    fn alignment(&self) -> &MultipleAlignment<A> {
        self.alignment()
    }
}

impl<A: MafftAlphabet> MultipleAlignmentView for MafftResult<A> {
    type Alphabet = A;

    fn alignment(&self) -> &MultipleAlignment<A> {
        self.alignment()
    }
}

impl<A: MafftAlphabet> MultipleAlignmentView for Mafft7526Alignment<A> {
    type Alphabet = A;

    fn alignment(&self) -> &MultipleAlignment<A> {
        self.alignment()
    }
}

impl<A: Alphabet> MultipleAlignmentView for &MultipleAlignment<A> {
    type Alphabet = A;

    fn alignment(&self) -> &MultipleAlignment<A> {
        self
    }
}
