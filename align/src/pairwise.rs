//! Fluent entry points for pairwise alignment.

use align_core::{
    AffineGap, Alphabet, LinearGap, MatchMismatch, PairwiseAlignment, Score, Sequence,
};
use align_fast::{
    MyersMatch, MyersPattern, WavefrontAffine, WavefrontAffinePenalties, WavefrontAlignment,
    WavefrontScore, WavefrontWorkspace, XDrop, XDropAlignment, XDropDirection, XDropScore,
    XDropSeed, XDropSeedAlignment, XDropSeedScore, adaptive_affine_decision, adaptive_affine_plan,
};
use align_pairwise::{
    AffineEndGapAligner, AffineLocalAligner, AffineWorkspace, FreeEndGaps, Gotoh, NeedlemanWunsch,
    SmithWaterman,
};

use crate::report::{
    MyersDistanceRun, MyersSearchRun, PairwiseMetric, PairwiseRun, PairwiseRunInfo,
    PairwiseValueRun, ReportField,
};
use crate::{Error, Result};

pub use align_fast::{AdaptiveAffineBackend, AdaptiveAffineDecision};

const DEFAULT_MATCH_SCORE: Score = 2;
const DEFAULT_MISMATCH_SCORE: Score = -1;
const DEFAULT_LINEAR_GAP_SCORE: Score = -2;
const DEFAULT_AFFINE_GAP_OPEN: Score = -5;
const DEFAULT_AFFINE_GAP_EXTENSION: Score = -1;
const DEFAULT_WFA_MISMATCH: Score = 4;
const DEFAULT_WFA_GAP_OPEN: Score = 6;
const DEFAULT_WFA_GAP_EXTENSION: Score = 2;
const DEFAULT_XDROP_THRESHOLD: Score = 20;

/// Discoverable constructors for the seven pairwise algorithm families.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Pairwise;

impl Pairwise {
    /// Selects exact global Needleman--Wunsch with simple match/mismatch scoring.
    #[must_use]
    pub const fn needleman_wunsch() -> LinearAligner {
        LinearAligner::new(LinearAlgorithm::NeedlemanWunsch)
    }

    /// Selects exact local Smith--Waterman with simple match/mismatch scoring.
    #[must_use]
    pub const fn smith_waterman() -> LinearAligner {
        LinearAligner::new(LinearAlgorithm::SmithWaterman)
    }

    /// Selects exact affine Gotoh alignment, globally by default.
    #[must_use]
    pub const fn gotoh() -> GotohAligner {
        GotohAligner::new()
    }

    /// Selects exact Myers unit-edit alignment and search.
    #[must_use]
    pub const fn myers() -> MyersAligner {
        MyersAligner { thread_count: 1 }
    }

    /// Selects exact global affine wavefront alignment.
    #[must_use]
    pub const fn wfa() -> WavefrontAligner {
        WavefrontAligner::new()
    }

    /// Selects an exact global affine identity/WFA/SIMD-Gotoh crossover.
    #[must_use]
    pub const fn adaptive_affine() -> AdaptiveAffineAligner {
        AdaptiveAffineAligner::new()
    }

    /// Selects affine X-drop extension.
    #[must_use]
    pub const fn xdrop() -> XDropAligner {
        XDropAligner::new()
    }
}

/// Linear-gap pairwise algorithm selected by [`LinearAligner`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LinearAlgorithm {
    /// Exact global Needleman--Wunsch.
    NeedlemanWunsch,
    /// Exact local Smith--Waterman.
    SmithWaterman,
}

/// Reusable linear-gap pairwise configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LinearAligner {
    algorithm: LinearAlgorithm,
    match_score: Score,
    mismatch_score: Score,
    gap_score: Score,
    thread_count: usize,
}

impl LinearAligner {
    const fn new(algorithm: LinearAlgorithm) -> Self {
        Self {
            algorithm,
            match_score: DEFAULT_MATCH_SCORE,
            mismatch_score: DEFAULT_MISMATCH_SCORE,
            gap_score: DEFAULT_LINEAR_GAP_SCORE,
            thread_count: 1,
        }
    }

    /// Replaces literal match and mismatch scores.
    #[must_use]
    pub const fn scores(mut self, match_score: Score, mismatch_score: Score) -> Self {
        self.match_score = match_score;
        self.mismatch_score = mismatch_score;
        self
    }

    /// Replaces the per-symbol linear gap score.
    #[must_use]
    pub const fn gap(mut self, score: Score) -> Self {
        self.gap_score = score;
        self
    }

    /// Replaces the worker count used by [`Self::align_many`].
    #[must_use]
    pub const fn threads(mut self, count: usize) -> Self {
        self.thread_count = count;
        self
    }

    /// Returns the selected algorithm.
    #[must_use]
    pub const fn algorithm(self) -> LinearAlgorithm {
        self.algorithm
    }

    /// Returns the configured literal scorer.
    #[must_use]
    pub const fn scoring(self) -> MatchMismatch {
        MatchMismatch::new(self.match_score, self.mismatch_score)
    }

    /// Returns the configured gap after validating it.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Scoring`] when the gap score is positive.
    pub const fn linear_gap(self) -> Result<LinearGap> {
        match LinearGap::new(self.gap_score) {
            Ok(gap) => Ok(gap),
            Err(error) => Err(Error::Scoring(error)),
        }
    }

    /// Computes one deterministic alignment.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] for invalid scoring, dimensions, allocation, arithmetic,
    /// traceback, or result validation.
    pub fn align<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<PairwiseAlignment> {
        let substitution = self.scoring();
        let gap = self.linear_gap()?;
        match self.algorithm {
            LinearAlgorithm::NeedlemanWunsch => {
                NeedlemanWunsch::new(substitution, gap).align(first, second)
            }
            LinearAlgorithm::SmithWaterman => {
                SmithWaterman::new(substitution, gap).align(first, second)
            }
        }
        .map_err(Error::from)
    }

    /// Computes one alignment and retains the inputs and configuration for a
    /// human-readable report.
    ///
    /// # Errors
    ///
    /// Returns the same structured failures as [`Self::align`].
    pub fn run<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
    ) -> Result<PairwiseRun<'a, A, PairwiseAlignment>> {
        let result = self.align(first, second)?;
        let (algorithm, mode) = match self.algorithm {
            LinearAlgorithm::NeedlemanWunsch => ("Needleman-Wunsch", "global linear"),
            LinearAlgorithm::SmithWaterman => ("Smith-Waterman", "local linear"),
        };
        let info = PairwiseRunInfo::new(
            algorithm,
            mode,
            format!(
                "match={} mismatch={} gap={}",
                self.match_score, self.mismatch_score, self.gap_score
            ),
            PairwiseMetric::Score(result.score()),
        );
        Ok(PairwiseRun::new(first, second, result, info))
    }

    /// Computes the selected objective value without a traceback.
    ///
    /// The score-only kernels evaluate the same recurrence as [`Self::align`]
    /// over one rolling row, so the value equals `align(first, second)?.score()`
    /// while memory stays linear in `second.len()` instead of the product of both
    /// lengths.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] for invalid scoring, dimension or score overflow, or a
    /// failed score-row allocation.
    pub fn score<A: Alphabet>(&self, first: &Sequence<A>, second: &Sequence<A>) -> Result<Score> {
        let substitution = self.scoring();
        let gap = self.linear_gap()?;
        match self.algorithm {
            LinearAlgorithm::NeedlemanWunsch => {
                NeedlemanWunsch::new(substitution, gap).score(first, second)
            }
            LinearAlgorithm::SmithWaterman => {
                SmithWaterman::new(substitution, gap).score(first, second)
            }
        }
        .map_err(Error::from)
    }

    /// Computes only the selected objective value and retains reporting context.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::score`].
    pub fn run_score<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
    ) -> Result<PairwiseValueRun<'a, A, Score>> {
        let score = self.score(first, second)?;
        let (algorithm, mode) = match self.algorithm {
            LinearAlgorithm::NeedlemanWunsch => ("Needleman-Wunsch", "global linear score-only"),
            LinearAlgorithm::SmithWaterman => ("Smith-Waterman", "local linear score-only"),
        };
        let info = PairwiseRunInfo::new(
            algorithm,
            mode,
            format!(
                "match={} mismatch={} gap={}",
                self.match_score, self.mismatch_score, self.gap_score
            ),
            PairwiseMetric::Score(score),
        );
        Ok(PairwiseValueRun::new(first, second, score, info))
    }

    /// Aligns an iterator of independent pairs with the same configuration.
    ///
    /// # Errors
    ///
    /// Stops at and returns the first failure.
    pub fn align_many<'a, A, I>(&self, pairs: I) -> Result<Vec<PairwiseAlignment>>
    where
        A: Alphabet,
        I: IntoIterator<Item = (&'a Sequence<A>, &'a Sequence<A>)>,
    {
        align_many_ordered(pairs, self.thread_count, |first, second| {
            self.align(first, second)
        })
    }
}

/// Boundary objective selected for affine Gotoh alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GotohMode {
    /// Consume both complete sequences.
    Global,
    /// Find a maximum-scoring local alignment.
    Local,
    /// Apply explicit free-end permissions.
    EndGaps(FreeEndGaps),
}

/// Reusable affine Gotoh configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GotohAligner {
    mode: GotohMode,
    match_score: Score,
    mismatch_score: Score,
    gap_open: Score,
    gap_extension: Score,
    thread_count: usize,
}

impl GotohAligner {
    const fn new() -> Self {
        Self {
            mode: GotohMode::Global,
            match_score: DEFAULT_MATCH_SCORE,
            mismatch_score: DEFAULT_MISMATCH_SCORE,
            gap_open: DEFAULT_AFFINE_GAP_OPEN,
            gap_extension: DEFAULT_AFFINE_GAP_EXTENSION,
            thread_count: 1,
        }
    }

    /// Replaces literal match and mismatch scores.
    #[must_use]
    pub const fn scores(mut self, match_score: Score, mismatch_score: Score) -> Self {
        self.match_score = match_score;
        self.mismatch_score = mismatch_score;
        self
    }

    /// Replaces affine scores using `open + (length - 1) * extension`.
    #[must_use]
    pub const fn gaps(mut self, open: Score, extension: Score) -> Self {
        self.gap_open = open;
        self.gap_extension = extension;
        self
    }

    /// Replaces the worker count used by [`Self::align_many`].
    #[must_use]
    pub const fn threads(mut self, count: usize) -> Self {
        self.thread_count = count;
        self
    }

    /// Requires both complete sequences.
    #[must_use]
    pub const fn global(mut self) -> Self {
        self.mode = GotohMode::Global;
        self
    }

    /// Finds the best local affine alignment.
    #[must_use]
    pub const fn local(mut self) -> Self {
        self.mode = GotohMode::Local;
        self
    }

    /// Makes all four terminal overhangs free.
    #[must_use]
    pub const fn semiglobal(mut self) -> Self {
        self.mode = GotohMode::EndGaps(FreeEndGaps::semiglobal());
        self
    }

    /// Fits the complete first sequence inside the second.
    #[must_use]
    pub const fn fit_first_in_second(mut self) -> Self {
        self.mode = GotohMode::EndGaps(FreeEndGaps::fit_first_in_second());
        self
    }

    /// Fits the complete second sequence inside the first.
    #[must_use]
    pub const fn fit_second_in_first(mut self) -> Self {
        self.mode = GotohMode::EndGaps(FreeEndGaps::fit_second_in_first());
        self
    }

    /// Aligns a suffix of the first sequence to a prefix of the second.
    #[must_use]
    pub const fn first_suffix_second_prefix(mut self) -> Self {
        self.mode = GotohMode::EndGaps(FreeEndGaps::first_suffix_second_prefix());
        self
    }

    /// Aligns a suffix of the second sequence to a prefix of the first.
    #[must_use]
    pub const fn second_suffix_first_prefix(mut self) -> Self {
        self.mode = GotohMode::EndGaps(FreeEndGaps::second_suffix_first_prefix());
        self
    }

    /// Applies exact four-end permissions.
    #[must_use]
    pub const fn free_ends(mut self, free_ends: FreeEndGaps) -> Self {
        self.mode = GotohMode::EndGaps(free_ends);
        self
    }

    /// Returns the selected boundary objective.
    #[must_use]
    pub const fn mode(self) -> GotohMode {
        self.mode
    }

    /// Returns the configured literal scorer.
    #[must_use]
    pub const fn scoring(self) -> MatchMismatch {
        MatchMismatch::new(self.match_score, self.mismatch_score)
    }

    /// Returns the configured affine gap after validating it.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Scoring`] for positive or inconsistent affine scores.
    pub const fn affine_gap(self) -> Result<AffineGap> {
        match AffineGap::new(self.gap_open, self.gap_extension) {
            Ok(gap) => Ok(gap),
            Err(error) => Err(Error::Scoring(error)),
        }
    }

    /// Computes one deterministic affine alignment.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] for invalid scoring, dimensions, allocation, arithmetic,
    /// traceback, or result validation.
    pub fn align<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<PairwiseAlignment> {
        let substitution = self.scoring();
        let gap = self.affine_gap()?;
        match self.mode {
            GotohMode::Global => Gotoh::new(substitution, gap).align(first, second),
            GotohMode::Local => AffineLocalAligner::new(substitution, gap).align(first, second),
            GotohMode::EndGaps(free_ends) => {
                AffineEndGapAligner::new(substitution, gap, free_ends).align(first, second)
            }
        }
        .map_err(Error::from)
    }

    /// Computes one affine alignment and retains enough context to report its
    /// exact boundary and gap conventions.
    ///
    /// # Errors
    ///
    /// Returns the same structured failures as [`Self::align`].
    pub fn run<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
    ) -> Result<PairwiseRun<'a, A, PairwiseAlignment>> {
        let result = self.align(first, second)?;
        let mut details = vec![ReportField::new(
            "Gap convention",
            "open + (length - 1) * extension",
        )];
        if let GotohMode::EndGaps(free_ends) = self.mode {
            details.push(ReportField::new(
                "Free ends",
                format!(
                    "first_start={} first_end={} second_start={} second_end={}",
                    free_ends.first_start(),
                    free_ends.first_end(),
                    free_ends.second_start(),
                    free_ends.second_end()
                ),
            ));
        }
        let info = PairwiseRunInfo::new(
            "Gotoh",
            gotoh_mode_name(self.mode),
            format!(
                "match={} mismatch={} gap_open={} gap_extend={}",
                self.match_score, self.mismatch_score, self.gap_open, self.gap_extension
            ),
            PairwiseMetric::Score(result.score()),
        )
        .with_details(details);
        Ok(PairwiseRun::new(first, second, result, info))
    }

    /// Computes the objective value. Global mode uses the score-only kernel;
    /// other boundary modes reuse their complete deterministic alignment.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn score<A: Alphabet>(&self, first: &Sequence<A>, second: &Sequence<A>) -> Result<Score> {
        let substitution = self.scoring();
        let gap = self.affine_gap()?;
        match self.mode {
            GotohMode::Global => Gotoh::new(substitution, gap)
                .score(first, second)
                .map_err(Error::from),
            GotohMode::Local | GotohMode::EndGaps(_) => {
                self.align(first, second).map(|alignment| alignment.score())
            }
        }
    }

    /// Computes the affine objective value and retains boundary and scoring
    /// context without requiring callers to format a bare integer.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::score`].
    pub fn run_score<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
    ) -> Result<PairwiseValueRun<'a, A, Score>> {
        let score = self.score(first, second)?;
        let mut details = vec![ReportField::new(
            "Gap convention",
            "open + (length - 1) * extension",
        )];
        if let GotohMode::EndGaps(free_ends) = self.mode {
            details.push(ReportField::new(
                "Free ends",
                format!(
                    "first_start={} first_end={} second_start={} second_end={}",
                    free_ends.first_start(),
                    free_ends.first_end(),
                    free_ends.second_start(),
                    free_ends.second_end()
                ),
            ));
        }
        let info = PairwiseRunInfo::new(
            "Gotoh",
            gotoh_score_mode_name(self.mode),
            format!(
                "match={} mismatch={} gap_open={} gap_extend={}",
                self.match_score, self.mismatch_score, self.gap_open, self.gap_extension
            ),
            PairwiseMetric::Score(score),
        )
        .with_details(details);
        Ok(PairwiseValueRun::new(first, second, score, info))
    }

    /// Aligns an iterator of independent pairs with the same configuration.
    ///
    /// # Errors
    ///
    /// Stops at and returns the first failure.
    pub fn align_many<'a, A, I>(&self, pairs: I) -> Result<Vec<PairwiseAlignment>>
    where
        A: Alphabet,
        I: IntoIterator<Item = (&'a Sequence<A>, &'a Sequence<A>)>,
    {
        align_many_ordered(pairs, self.thread_count, |first, second| {
            self.align(first, second)
        })
    }
}

/// Stateless facade for Myers pattern preprocessing, distance, search, and traceback.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MyersAligner {
    thread_count: usize,
}

impl MyersAligner {
    /// Replaces the worker count used by [`Self::align_many`].
    #[must_use]
    pub const fn threads(mut self, count: usize) -> Self {
        self.thread_count = count;
        self
    }

    /// Preprocesses a pattern for repeated comparisons.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Myers`] for unrepresentable dimensions or allocation failure.
    pub fn prepare<A: Alphabet>(&self, pattern: &Sequence<A>) -> Result<MyersPattern<A>> {
        MyersPattern::new(pattern).map_err(Error::from)
    }

    /// Computes exact global unit edit distance.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Myers`] if preprocessing or execution fails.
    pub fn distance<A: Alphabet>(
        &self,
        pattern: &Sequence<A>,
        text: &Sequence<A>,
    ) -> Result<usize> {
        self.prepare(pattern)?
            .global_distance(text)
            .map_err(Error::from)
    }

    /// Computes exact unit edit distance and retains both inputs for reporting.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::distance`].
    pub fn run_distance<'a, A: Alphabet>(
        &self,
        pattern: &'a Sequence<A>,
        text: &'a Sequence<A>,
    ) -> Result<MyersDistanceRun<'a, A>> {
        let distance = self.distance(pattern, text)?;
        Ok(MyersDistanceRun::new(pattern, text, distance))
    }

    /// Computes one exact global unit-edit traceback.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Myers`] if preprocessing or execution fails.
    pub fn align<A: Alphabet>(
        &self,
        pattern: &Sequence<A>,
        text: &Sequence<A>,
    ) -> Result<PairwiseAlignment> {
        self.prepare(pattern)?
            .align_global(text)
            .map_err(Error::from)
    }

    /// Computes one exact unit-edit traceback and retains reporting context.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn run<'a, A: Alphabet>(
        &self,
        pattern: &'a Sequence<A>,
        text: &'a Sequence<A>,
    ) -> Result<PairwiseRun<'a, A, PairwiseAlignment>> {
        let result = self.align(pattern, text)?;
        let distance = result
            .operations()
            .iter()
            .filter(|operation| !matches!(operation, align_core::AlignmentOperation::Match))
            .count();
        let info = PairwiseRunInfo::new(
            "Myers bit-vector",
            "global unit edit traceback",
            String::from("match=0 mismatch=1 insertion=1 deletion=1"),
            PairwiseMetric::Distance(distance),
        );
        Ok(PairwiseRun::new(pattern, text, result, info))
    }

    /// Reports every suffix endpoint within `maximum_distance` edits.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Myers`] if preprocessing or execution fails.
    pub fn find_all<A: Alphabet>(
        &self,
        pattern: &Sequence<A>,
        text: &Sequence<A>,
        maximum_distance: usize,
    ) -> Result<Vec<MyersMatch>> {
        self.prepare(pattern)?
            .find_all(text, maximum_distance)
            .map_err(Error::from)
    }

    /// Searches every suffix endpoint and retains a presentable result table.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::find_all`].
    pub fn run_find_all<'a, A: Alphabet>(
        &self,
        pattern: &'a Sequence<A>,
        text: &'a Sequence<A>,
        maximum_distance: usize,
    ) -> Result<MyersSearchRun<'a, A>> {
        let matches = self.find_all(pattern, text, maximum_distance)?;
        Ok(MyersSearchRun::new(
            pattern,
            text,
            maximum_distance,
            matches,
        ))
    }

    /// Aligns an iterator of pattern/text pairs with unit edit distance.
    ///
    /// # Errors
    ///
    /// Stops at and returns the first preprocessing or execution failure.
    pub fn align_many<'a, A, I>(&self, pairs: I) -> Result<Vec<PairwiseAlignment>>
    where
        A: Alphabet,
        I: IntoIterator<Item = (&'a Sequence<A>, &'a Sequence<A>)>,
    {
        align_many_ordered(pairs, self.thread_count, |pattern, text| {
            self.align(pattern, text)
        })
    }
}

/// Caller-owned scratch storage for the adaptive exact-affine backend.
#[derive(Debug, Default)]
pub struct AdaptiveAffineWorkspace {
    wavefront: WavefrontWorkspace,
    gotoh: AffineWorkspace,
}

impl AdaptiveAffineWorkspace {
    /// Creates empty WFA and Gotoh workspaces that allocate on first use.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            wavefront: WavefrontWorkspace::new(),
            gotoh: AffineWorkspace::new(),
        }
    }

    /// Releases every allocation retained by either backend.
    pub fn release(&mut self) {
        *self = Self::new();
    }
}

/// Exact adaptive affine penalty together with backend provenance.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AdaptiveAffineScore {
    penalty: Score,
    decision: AdaptiveAffineDecision,
}

impl AdaptiveAffineScore {
    /// Returns the exact nonnegative minimizing penalty.
    #[must_use]
    pub const fn penalty(self) -> Score {
        self.penalty
    }

    /// Returns the conventional maximizing score, `-penalty`.
    #[must_use]
    pub const fn alignment_score(self) -> Score {
        -self.penalty
    }

    /// Returns the exact backend selected for this call.
    #[must_use]
    pub const fn backend(self) -> AdaptiveAffineBackend {
        self.decision.backend()
    }

    /// Returns the complete auditable cost-model decision.
    #[must_use]
    pub const fn decision(self) -> AdaptiveAffineDecision {
        self.decision
    }
}

/// Exact adaptive affine traceback together with backend provenance.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AdaptiveAffineAlignment {
    alignment: PairwiseAlignment,
    penalty: Score,
    decision: AdaptiveAffineDecision,
}

impl AdaptiveAffineAlignment {
    /// Borrows the canonical exact global alignment.
    #[must_use]
    pub const fn alignment(&self) -> &PairwiseAlignment {
        &self.alignment
    }

    /// Consumes this result and returns the canonical alignment.
    #[must_use]
    pub fn into_alignment(self) -> PairwiseAlignment {
        self.alignment
    }

    /// Returns the exact nonnegative minimizing penalty.
    #[must_use]
    pub const fn penalty(&self) -> Score {
        self.penalty
    }

    /// Returns the exact backend selected for this call.
    #[must_use]
    pub const fn backend(&self) -> AdaptiveAffineBackend {
        self.decision.backend()
    }

    /// Returns the complete auditable cost-model decision.
    #[must_use]
    pub const fn decision(&self) -> AdaptiveAffineDecision {
        self.decision
    }
}

/// Reusable exact affine identity/WFA/SIMD-Gotoh crossover configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AdaptiveAffineAligner {
    mismatch: Score,
    gap_open: Score,
    gap_extension: Score,
    thread_count: usize,
}

impl AdaptiveAffineAligner {
    const fn new() -> Self {
        Self {
            mismatch: DEFAULT_WFA_MISMATCH,
            gap_open: DEFAULT_WFA_GAP_OPEN,
            gap_extension: DEFAULT_WFA_GAP_EXTENSION,
            thread_count: 1,
        }
    }

    /// Replaces positive WFA penalties using `open + length * extension`.
    #[must_use]
    pub const fn penalties(
        mut self,
        mismatch: Score,
        gap_open: Score,
        gap_extension: Score,
    ) -> Self {
        self.mismatch = mismatch;
        self.gap_open = gap_open;
        self.gap_extension = gap_extension;
        self
    }

    /// Replaces the worker count used by [`Self::align_many`].
    #[must_use]
    pub const fn threads(mut self, count: usize) -> Self {
        self.thread_count = count;
        self
    }

    fn penalties_checked(self) -> Result<WavefrontAffinePenalties> {
        WavefrontAffinePenalties::new(self.mismatch, self.gap_open, self.gap_extension)
            .map_err(Error::from)
    }

    fn gotoh(penalties: WavefrontAffinePenalties) -> Result<Gotoh<MatchMismatch>> {
        let opening = penalties.gap_penalty(1)?;
        let gap = AffineGap::new(-opening, -penalties.gap_extension())?;
        Ok(Gotoh::new(
            MatchMismatch::new(0, -penalties.mismatch()),
            gap,
        ))
    }

    /// Returns the deterministic backend decision without executing an alignment.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Wavefront`] when the configured penalties are invalid.
    pub fn decision<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<AdaptiveAffineDecision> {
        let penalties = self.penalties_checked()?;
        Ok(adaptive_affine_decision(first, second, penalties))
    }

    /// Computes the exact penalty with fresh scratch storage.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] for invalid penalties, allocation failure, arithmetic
    /// overflow, or an internal backend invariant failure.
    pub fn penalty<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<AdaptiveAffineScore> {
        let mut workspace = AdaptiveAffineWorkspace::new();
        self.penalty_with_workspace(first, second, &mut workspace)
    }

    /// Computes the exact penalty with caller-owned WFA and Gotoh scratch storage.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::penalty`].
    pub fn penalty_with_workspace<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut AdaptiveAffineWorkspace,
    ) -> Result<AdaptiveAffineScore> {
        let penalties = self.penalties_checked()?;
        let plan = adaptive_affine_plan(first, second, penalties);
        let decision = plan.decision();
        let penalty = match decision.backend() {
            AdaptiveAffineBackend::Identity => 0,
            AdaptiveAffineBackend::Wavefront => {
                plan.wavefront_penalty_with_workspace(&mut workspace.wavefront)?
            }
            AdaptiveAffineBackend::Gotoh => {
                let score = Self::gotoh(penalties)?.score_with_workspace(
                    first,
                    second,
                    &mut workspace.gotoh,
                )?;
                score
                    .checked_neg()
                    .ok_or(Error::AdaptivePenaltyOverflow { score })?
            }
        };
        Ok(AdaptiveAffineScore { penalty, decision })
    }

    /// Computes one deterministic exact global affine alignment with fresh scratch.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] for invalid penalties, dimensions, allocation, arithmetic,
    /// traceback, or canonical-result validation failure.
    pub fn align<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<AdaptiveAffineAlignment> {
        let mut workspace = AdaptiveAffineWorkspace::new();
        self.align_with_workspace(first, second, &mut workspace)
    }

    /// Computes one exact alignment with caller-owned WFA and Gotoh scratch storage.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn align_with_workspace<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        workspace: &mut AdaptiveAffineWorkspace,
    ) -> Result<AdaptiveAffineAlignment> {
        let penalties = self.penalties_checked()?;
        let plan = adaptive_affine_plan(first, second, penalties);
        let decision = plan.decision();
        let (alignment, penalty) = match decision.backend() {
            AdaptiveAffineBackend::Identity => {
                let mut operations = Vec::new();
                operations.try_reserve_exact(first.len()).map_err(|_| {
                    align_core::AlignmentStructureError::AllocationFailed {
                        columns: first.len(),
                    }
                })?;
                operations.resize(first.len(), align_core::AlignmentOperation::Match);
                (
                    PairwiseAlignment::try_new(0, 0..first.len(), 0..second.len(), operations)?,
                    0,
                )
            }
            AdaptiveAffineBackend::Wavefront => {
                let result = plan.wavefront_align_with_workspace(&mut workspace.wavefront)?;
                let penalty = result.penalty();
                (result.into_alignment(), penalty)
            }
            AdaptiveAffineBackend::Gotoh => {
                let alignment = Self::gotoh(penalties)?.align_with_workspace(
                    first,
                    second,
                    &mut workspace.gotoh,
                )?;
                let score = alignment.score();
                let penalty = score
                    .checked_neg()
                    .ok_or(Error::AdaptivePenaltyOverflow { score })?;
                (alignment, penalty)
            }
        };
        Ok(AdaptiveAffineAlignment {
            alignment,
            penalty,
            decision,
        })
    }

    /// Computes one alignment and retains backend provenance for reporting.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn run<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
    ) -> Result<PairwiseRun<'a, A, AdaptiveAffineAlignment>> {
        let result = self.align(first, second)?;
        let decision = result.decision();
        let info = PairwiseRunInfo::new(
            "Adaptive exact affine",
            "global identity/WFA/SIMD-Gotoh crossover",
            format!(
                "mismatch={} gap_open={} gap_extend={} (positive WFA penalties)",
                self.mismatch, self.gap_open, self.gap_extension
            ),
            PairwiseMetric::Penalty(result.penalty()),
        )
        .with_details(adaptive_affine_report_fields(decision));
        Ok(PairwiseRun::new(first, second, result, info))
    }

    /// Computes only the exact penalty and retains backend provenance for reporting.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::penalty`].
    pub fn run_penalty<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
    ) -> Result<PairwiseValueRun<'a, A, AdaptiveAffineScore>> {
        let result = self.penalty(first, second)?;
        let info = PairwiseRunInfo::new(
            "Adaptive exact affine",
            "global identity/WFA/SIMD-Gotoh score-only crossover",
            format!(
                "mismatch={} gap_open={} gap_extend={} (positive WFA penalties)",
                self.mismatch, self.gap_open, self.gap_extension
            ),
            PairwiseMetric::Penalty(result.penalty()),
        )
        .with_details(adaptive_affine_report_fields(result.decision()));
        Ok(PairwiseValueRun::new(first, second, result, info))
    }

    /// Aligns an iterator of independent pairs with the same adaptive policy.
    ///
    /// # Errors
    ///
    /// Stops at and returns the first configuration or execution failure.
    pub fn align_many<'a, A, I>(&self, pairs: I) -> Result<Vec<AdaptiveAffineAlignment>>
    where
        A: Alphabet,
        I: IntoIterator<Item = (&'a Sequence<A>, &'a Sequence<A>)>,
    {
        align_many_ordered(pairs, self.thread_count, |first, second| {
            self.align(first, second)
        })
    }
}

fn adaptive_affine_report_fields(decision: AdaptiveAffineDecision) -> Vec<ReportField> {
    vec![
        ReportField::new("Selected backend", decision.backend().as_str()),
        ReportField::new(
            "Estimated WFA levels",
            decision.estimated_wavefront_levels().to_string(),
        ),
        ReportField::new("Dense cells", decision.dense_cells().to_string()),
    ]
}

/// Reusable exact affine wavefront configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct WavefrontAligner {
    mismatch: Score,
    gap_open: Score,
    gap_extension: Score,
    thread_count: usize,
}

impl WavefrontAligner {
    const fn new() -> Self {
        Self {
            mismatch: DEFAULT_WFA_MISMATCH,
            gap_open: DEFAULT_WFA_GAP_OPEN,
            gap_extension: DEFAULT_WFA_GAP_EXTENSION,
            thread_count: 1,
        }
    }

    /// Replaces positive WFA penalties using `open + length * extension`.
    #[must_use]
    pub const fn penalties(
        mut self,
        mismatch: Score,
        gap_open: Score,
        gap_extension: Score,
    ) -> Self {
        self.mismatch = mismatch;
        self.gap_open = gap_open;
        self.gap_extension = gap_extension;
        self
    }

    /// Replaces the worker count used by [`Self::align_many`].
    #[must_use]
    pub const fn threads(mut self, count: usize) -> Self {
        self.thread_count = count;
        self
    }

    fn kernel(self) -> Result<WavefrontAffine> {
        WavefrontAffinePenalties::new(self.mismatch, self.gap_open, self.gap_extension)
            .map(WavefrontAffine::new)
            .map_err(Error::from)
    }

    /// Computes only the exact global affine penalty.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Wavefront`] for invalid penalties or execution failure.
    pub fn penalty<A: Alphabet>(&self, first: &Sequence<A>, second: &Sequence<A>) -> Result<Score> {
        self.kernel()?.penalty(first, second).map_err(Error::from)
    }

    /// Computes a score-only WFA execution summary and presents its primary value
    /// explicitly as a minimizing penalty.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Wavefront`] for invalid penalties or execution failure.
    pub fn run_penalty<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
    ) -> Result<PairwiseValueRun<'a, A, WavefrontScore>> {
        let summary = self.kernel()?.score(first, second)?;
        let info = PairwiseRunInfo::new(
            "WFA",
            "global affine wavefront score-only",
            format!(
                "mismatch={} gap_open={} gap_extend={} (positive penalties)",
                self.mismatch, self.gap_open, self.gap_extension
            ),
            PairwiseMetric::Penalty(summary.penalty()),
        )
        .with_details(wfa_report_fields(summary));
        Ok(PairwiseValueRun::new(first, second, summary, info))
    }

    /// Computes one deterministic exact global affine alignment.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Wavefront`] for invalid penalties or execution failure.
    pub fn align<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<WavefrontAlignment> {
        self.kernel()?.align(first, second).map_err(Error::from)
    }

    /// Computes one WFA traceback and retains its minimizing-penalty semantics and
    /// execution telemetry for reporting.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn run<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
    ) -> Result<PairwiseRun<'a, A, WavefrontAlignment>> {
        let result = self.align(first, second)?;
        let summary = result.summary();
        let info = PairwiseRunInfo::new(
            "WFA",
            "global affine wavefront",
            format!(
                "mismatch={} gap_open={} gap_extend={} (positive penalties)",
                self.mismatch, self.gap_open, self.gap_extension
            ),
            PairwiseMetric::Penalty(result.penalty()),
        )
        .with_details(wfa_report_fields(summary));
        Ok(PairwiseRun::new(first, second, result, info))
    }

    /// Aligns an iterator of independent pairs with the same penalties.
    ///
    /// # Errors
    ///
    /// Stops at and returns the first failure.
    pub fn align_many<'a, A, I>(&self, pairs: I) -> Result<Vec<WavefrontAlignment>>
    where
        A: Alphabet,
        I: IntoIterator<Item = (&'a Sequence<A>, &'a Sequence<A>)>,
    {
        align_many_ordered(pairs, self.thread_count, |first, second| {
            self.align(first, second)
        })
    }
}

/// Reusable affine X-drop configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct XDropAligner {
    match_score: Score,
    mismatch_score: Score,
    gap_open: Score,
    gap_extension: Score,
    threshold: Score,
    thread_count: usize,
}

impl XDropAligner {
    const fn new() -> Self {
        Self {
            match_score: DEFAULT_MATCH_SCORE,
            mismatch_score: DEFAULT_MISMATCH_SCORE,
            gap_open: DEFAULT_AFFINE_GAP_OPEN,
            gap_extension: DEFAULT_AFFINE_GAP_EXTENSION,
            threshold: DEFAULT_XDROP_THRESHOLD,
            thread_count: 1,
        }
    }

    /// Replaces literal match and mismatch scores.
    #[must_use]
    pub const fn scores(mut self, match_score: Score, mismatch_score: Score) -> Self {
        self.match_score = match_score;
        self.mismatch_score = mismatch_score;
        self
    }

    /// Replaces affine scores using `open + (length - 1) * extension`.
    #[must_use]
    pub const fn gaps(mut self, open: Score, extension: Score) -> Self {
        self.gap_open = open;
        self.gap_extension = extension;
        self
    }

    /// Replaces the strict nonnegative X-drop threshold.
    #[must_use]
    pub const fn threshold(mut self, threshold: Score) -> Self {
        self.threshold = threshold;
        self
    }

    /// Replaces the worker count used by [`Self::align_many`].
    #[must_use]
    pub const fn threads(mut self, count: usize) -> Self {
        self.thread_count = count;
        self
    }

    fn kernel(self) -> Result<XDrop<MatchMismatch>> {
        let gap = AffineGap::new(self.gap_open, self.gap_extension)?;
        XDrop::new(
            MatchMismatch::new(self.match_score, self.mismatch_score),
            gap,
            self.threshold,
        )
        .map_err(Error::from)
    }

    /// Scores a forward X-drop extension from the origin.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] for invalid scoring or failed X-drop execution.
    pub fn score<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<XDropScore> {
        self.kernel()?.score(first, second).map_err(Error::from)
    }

    /// Computes a forward score-only X-drop traversal with its endpoints and
    /// stopping evidence.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::score`].
    pub fn run_score<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
    ) -> Result<PairwiseValueRun<'a, A, XDropScore>> {
        let summary = self.score(first, second)?;
        let info = self.xdrop_run_info(
            "forward score-only from origin",
            summary.score(),
            xdrop_score_fields(summary),
        );
        Ok(PairwiseValueRun::new(first, second, summary, info))
    }

    /// Computes a forward X-drop extension from the origin.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] for invalid scoring or failed X-drop execution.
    pub fn align<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
    ) -> Result<XDropAlignment> {
        self.kernel()?.align(first, second).map_err(Error::from)
    }

    /// Computes one forward X-drop traceback and retains traversal evidence for a
    /// detailed report.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align`].
    pub fn run<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
    ) -> Result<PairwiseRun<'a, A, XDropAlignment>> {
        let result = self.align(first, second)?;
        let summary = result.summary();
        let info = self.xdrop_run_info(
            "forward from origin",
            result.score(),
            vec![
                ReportField::new(
                    "First native end",
                    format!("{} (exclusive)", result.first_end()),
                ),
                ReportField::new(
                    "Second native end",
                    format!("{} (exclusive)", result.second_end()),
                ),
                ReportField::new("Evaluated cells", result.evaluated_cells().to_string()),
                ReportField::new("Retained cells", result.retained_cells().to_string()),
                ReportField::new("Antidiagonals", summary.antidiagonals().to_string()),
                ReportField::new(
                    "Maximum frontier",
                    summary.maximum_frontier_width().to_string(),
                ),
                ReportField::new("Termination", xdrop_termination_name(result.termination())),
            ],
        );
        Ok(PairwiseRun::new(first, second, result, info))
    }

    /// Computes an upstream, downstream, or two-way extension around a seed.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] for invalid scoring, an out-of-bounds seed, or failed
    /// X-drop execution.
    pub fn align_from_seed<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        seed: XDropSeed,
        direction: XDropDirection,
    ) -> Result<XDropSeedAlignment> {
        self.kernel()?
            .align_from_seed(first, second, seed, direction)
            .map_err(Error::from)
    }

    /// Computes one forced-seed X-drop traceback and retains the seed, direction,
    /// ranges, and traversal counters for reporting.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::align_from_seed`].
    pub fn run_from_seed<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
        seed: XDropSeed,
        direction: XDropDirection,
    ) -> Result<PairwiseRun<'a, A, XDropSeedAlignment>> {
        let result = self.align_from_seed(first, second, seed, direction)?;
        let summary = result.summary();
        let info = self.xdrop_run_info(
            xdrop_direction_name(direction),
            result.score(),
            vec![
                ReportField::new(
                    "Seed",
                    format!(
                        "first={} second={} (zero-based)",
                        seed.first_index(),
                        seed.second_index()
                    ),
                ),
                ReportField::new(
                    "First native range",
                    format!(
                        "{}..{}",
                        result.first_range().start,
                        result.first_range().end
                    ),
                ),
                ReportField::new(
                    "Second native range",
                    format!(
                        "{}..{}",
                        result.second_range().start,
                        result.second_range().end
                    ),
                ),
                ReportField::new("Evaluated cells", summary.evaluated_cells().to_string()),
                ReportField::new("Retained cells", summary.retained_cells().to_string()),
                ReportField::new("Antidiagonals", summary.antidiagonals().to_string()),
                ReportField::new(
                    "Maximum frontier",
                    summary.maximum_frontier_width().to_string(),
                ),
                ReportField::new(
                    "Upstream termination",
                    summary.upstream().map_or("not requested", |value| {
                        xdrop_termination_name(value.termination())
                    }),
                ),
                ReportField::new(
                    "Downstream termination",
                    summary.downstream().map_or("not requested", |value| {
                        xdrop_termination_name(value.termination())
                    }),
                ),
            ],
        );
        Ok(PairwiseRun::new(first, second, result, info))
    }

    /// Scores an upstream, downstream, or two-way extension around a seed.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] for invalid scoring, an out-of-bounds seed, or failed
    /// X-drop execution.
    pub fn score_from_seed<A: Alphabet>(
        &self,
        first: &Sequence<A>,
        second: &Sequence<A>,
        seed: XDropSeed,
        direction: XDropDirection,
    ) -> Result<XDropSeedScore> {
        self.kernel()?
            .score_from_seed(first, second, seed, direction)
            .map_err(Error::from)
    }

    /// Computes a seeded score-only X-drop traversal and retains directional
    /// evidence for presentation.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`Self::score_from_seed`].
    pub fn run_score_from_seed<'a, A: Alphabet>(
        &self,
        first: &'a Sequence<A>,
        second: &'a Sequence<A>,
        seed: XDropSeed,
        direction: XDropDirection,
    ) -> Result<PairwiseValueRun<'a, A, XDropSeedScore>> {
        let summary = self.score_from_seed(first, second, seed, direction)?;
        let info = self.xdrop_run_info(
            xdrop_direction_score_name(direction),
            summary.score(),
            xdrop_seed_score_fields(summary, seed),
        );
        Ok(PairwiseValueRun::new(first, second, summary, info))
    }

    /// Runs forward X-drop on an iterator of independent pairs.
    ///
    /// # Errors
    ///
    /// Stops at and returns the first failure.
    pub fn align_many<'a, A, I>(&self, pairs: I) -> Result<Vec<XDropAlignment>>
    where
        A: Alphabet,
        I: IntoIterator<Item = (&'a Sequence<A>, &'a Sequence<A>)>,
    {
        align_many_ordered(pairs, self.thread_count, |first, second| {
            self.align(first, second)
        })
    }

    fn xdrop_run_info(
        self,
        mode: &'static str,
        score: Score,
        mut details: Vec<ReportField>,
    ) -> PairwiseRunInfo {
        details.push(ReportField::new("Threshold", self.threshold.to_string()));
        details.push(ReportField::new(
            "Gap convention",
            "open + (length - 1) * extension",
        ));
        PairwiseRunInfo::new(
            "X-drop",
            mode,
            format!(
                "match={} mismatch={} gap_open={} gap_extend={}",
                self.match_score, self.mismatch_score, self.gap_open, self.gap_extension
            ),
            PairwiseMetric::Score(score),
        )
        .with_details(details)
    }
}

fn gotoh_mode_name(mode: GotohMode) -> &'static str {
    match mode {
        GotohMode::Global => "global affine",
        GotohMode::Local => "local affine",
        GotohMode::EndGaps(free_ends) if free_ends == FreeEndGaps::semiglobal() => {
            "semiglobal affine"
        }
        GotohMode::EndGaps(free_ends) if free_ends == FreeEndGaps::fit_first_in_second() => {
            "fit first in second"
        }
        GotohMode::EndGaps(free_ends) if free_ends == FreeEndGaps::fit_second_in_first() => {
            "fit second in first"
        }
        GotohMode::EndGaps(free_ends) if free_ends == FreeEndGaps::first_suffix_second_prefix() => {
            "first suffix to second prefix"
        }
        GotohMode::EndGaps(free_ends) if free_ends == FreeEndGaps::second_suffix_first_prefix() => {
            "second suffix to first prefix"
        }
        GotohMode::EndGaps(_) => "custom free-end affine",
    }
}

fn gotoh_score_mode_name(mode: GotohMode) -> &'static str {
    match mode {
        GotohMode::Global => "global affine score-only",
        GotohMode::Local => "local affine score-only",
        GotohMode::EndGaps(free_ends) if free_ends == FreeEndGaps::semiglobal() => {
            "semiglobal affine score-only"
        }
        GotohMode::EndGaps(free_ends) if free_ends == FreeEndGaps::fit_first_in_second() => {
            "fit first in second score-only"
        }
        GotohMode::EndGaps(free_ends) if free_ends == FreeEndGaps::fit_second_in_first() => {
            "fit second in first score-only"
        }
        GotohMode::EndGaps(free_ends) if free_ends == FreeEndGaps::first_suffix_second_prefix() => {
            "first suffix to second prefix score-only"
        }
        GotohMode::EndGaps(free_ends) if free_ends == FreeEndGaps::second_suffix_first_prefix() => {
            "second suffix to first prefix score-only"
        }
        GotohMode::EndGaps(_) => "custom free-end affine score-only",
    }
}

fn wfa_report_fields(summary: WavefrontScore) -> Vec<ReportField> {
    vec![
        ReportField::new("Gap convention", "open + length * extension"),
        ReportField::new("Alignment score", summary.alignment_score().to_string()),
        ReportField::new("Wavefronts", summary.wavefronts().to_string()),
        ReportField::new("Computed offsets", summary.computed_offsets().to_string()),
        ReportField::new("Match extensions", summary.match_extensions().to_string()),
        ReportField::new(
            "Maximum width",
            summary.maximum_wavefront_width().to_string(),
        ),
        ReportField::new(
            "Modeled peak bytes",
            summary.modeled_peak_bytes().to_string(),
        ),
    ]
}

fn xdrop_score_fields(summary: XDropScore) -> Vec<ReportField> {
    vec![
        ReportField::new(
            "First native end",
            format!("{} (exclusive)", summary.first_end()),
        ),
        ReportField::new(
            "Second native end",
            format!("{} (exclusive)", summary.second_end()),
        ),
        ReportField::new("Evaluated cells", summary.evaluated_cells().to_string()),
        ReportField::new("Retained cells", summary.retained_cells().to_string()),
        ReportField::new("Antidiagonals", summary.antidiagonals().to_string()),
        ReportField::new(
            "Maximum frontier",
            summary.maximum_frontier_width().to_string(),
        ),
        ReportField::new("Termination", xdrop_termination_name(summary.termination())),
    ]
}

fn xdrop_seed_score_fields(summary: XDropSeedScore, seed: XDropSeed) -> Vec<ReportField> {
    let first_range = summary.first_range();
    let second_range = summary.second_range();
    vec![
        ReportField::new(
            "Seed",
            format!(
                "first={} second={} (zero-based)",
                seed.first_index(),
                seed.second_index()
            ),
        ),
        ReportField::new(
            "First native range",
            format!("{}..{}", first_range.start, first_range.end),
        ),
        ReportField::new(
            "Second native range",
            format!("{}..{}", second_range.start, second_range.end),
        ),
        ReportField::new("Evaluated cells", summary.evaluated_cells().to_string()),
        ReportField::new("Retained cells", summary.retained_cells().to_string()),
        ReportField::new("Antidiagonals", summary.antidiagonals().to_string()),
        ReportField::new(
            "Maximum frontier",
            summary.maximum_frontier_width().to_string(),
        ),
        ReportField::new(
            "Upstream termination",
            summary.upstream().map_or("not requested", |value| {
                xdrop_termination_name(value.termination())
            }),
        ),
        ReportField::new(
            "Downstream termination",
            summary.downstream().map_or("not requested", |value| {
                xdrop_termination_name(value.termination())
            }),
        ),
    ]
}

const fn xdrop_direction_name(direction: XDropDirection) -> &'static str {
    match direction {
        XDropDirection::Upstream => "seeded upstream",
        XDropDirection::Downstream => "seeded downstream",
        XDropDirection::Both => "seeded two-way",
    }
}

const fn xdrop_direction_score_name(direction: XDropDirection) -> &'static str {
    match direction {
        XDropDirection::Upstream => "seeded upstream score-only",
        XDropDirection::Downstream => "seeded downstream score-only",
        XDropDirection::Both => "seeded two-way score-only",
    }
}

const fn xdrop_termination_name(termination: align_fast::XDropTermination) -> &'static str {
    match termination {
        align_fast::XDropTermination::MatrixExhausted => "matrix exhausted",
        align_fast::XDropTermination::FrontiersExhausted => "frontiers exhausted",
        _ => "unrecognized termination",
    }
}

impl Default for LinearAligner {
    fn default() -> Self {
        Pairwise::needleman_wunsch()
    }
}

impl Default for GotohAligner {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for MyersAligner {
    fn default() -> Self {
        Pairwise::myers()
    }
}

impl Default for WavefrontAligner {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for XDropAligner {
    fn default() -> Self {
        Self::new()
    }
}

fn align_many_ordered<'a, A, O, I, F>(pairs: I, thread_count: usize, operation: F) -> Result<Vec<O>>
where
    A: Alphabet,
    O: Send,
    I: IntoIterator<Item = (&'a Sequence<A>, &'a Sequence<A>)>,
    F: Fn(&Sequence<A>, &Sequence<A>) -> Result<O> + Sync,
{
    if thread_count == 0 {
        return Err(Error::InvalidThreadCount {
            value: thread_count,
        });
    }
    let pairs = pairs.into_iter().collect::<Vec<_>>();
    if thread_count == 1 || pairs.len() < 2 {
        return pairs
            .into_iter()
            .map(|(first, second)| operation(first, second))
            .collect();
    }

    let worker_count = thread_count.min(pairs.len());
    let chunk_size = pairs.len().div_ceil(worker_count);
    std::thread::scope(|scope| {
        let handles = pairs
            .chunks(chunk_size)
            .map(|chunk| {
                let operation = &operation;
                scope.spawn(move || {
                    chunk
                        .iter()
                        .map(|(first, second)| operation(first, second))
                        .collect::<Result<Vec<_>>>()
                })
            })
            .collect::<Vec<_>>();
        let mut output = Vec::with_capacity(pairs.len());
        for handle in handles {
            output.extend(handle.join().map_err(|_| Error::BatchWorkerPanicked)??);
        }
        Ok(output)
    })
}
