/// Immutable parameters of the `MAFFT-2002-FFT-NS-1-v1` scientific profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FftNs1Parameters {
    tuple_length: usize,
    peak_count: usize,
    window_length: usize,
    segment_threshold_per_column: f64,
    maximum_segment_length: usize,
    gap_open: f64,
    random_match_offset: f64,
}

impl FftNs1Parameters {
    /// Returns the defaults printed by Katoh et al. (2002).
    #[must_use]
    pub const fn paper_2002() -> Self {
        Self {
            tuple_length: 6,
            peak_count: 20,
            window_length: 30,
            segment_threshold_per_column: 0.7,
            maximum_segment_length: 150,
            gap_open: 2.4,
            random_match_offset: 0.06,
        }
    }

    pub(crate) const fn ginsi_2005_pairwise() -> Self {
        Self {
            tuple_length: 6,
            peak_count: 20,
            window_length: 30,
            segment_threshold_per_column: 0.7,
            maximum_segment_length: 150,
            gap_open: GINSI_GAP_OPEN,
            random_match_offset: GINSI_RANDOM_MATCH_OFFSET,
        }
    }

    /// Returns the normal compressed tuple length.
    #[must_use]
    pub const fn tuple_length(self) -> usize {
        self.tuple_length
    }

    /// Returns the number of Fourier-correlation peaks inspected per merge.
    #[must_use]
    pub const fn peak_count(self) -> usize {
        self.peak_count
    }

    /// Returns the local homology sliding-window length.
    #[must_use]
    pub const fn window_length(self) -> usize {
        self.window_length
    }

    /// Returns the strict local score threshold per alignment column.
    #[must_use]
    pub const fn segment_threshold_per_column(self) -> f64 {
        self.segment_threshold_per_column
    }

    /// Returns the maximum retained homologous-segment length.
    #[must_use]
    pub const fn maximum_segment_length(self) -> usize {
        self.maximum_segment_length
    }

    /// Returns the extension-free profile gap-opening cost.
    #[must_use]
    pub const fn gap_open(self) -> f64 {
        self.gap_open
    }

    /// Returns the normalized-matrix random-match offset.
    #[must_use]
    pub const fn random_match_offset(self) -> f64 {
        self.random_match_offset
    }
}

impl Default for FftNs1Parameters {
    fn default() -> Self {
        Self::paper_2002()
    }
}

/// Finite identity floor used when the first MSA has no identical comparable sites.
pub const FFT_NS2_POISSON_IDENTITY_FLOOR: f64 = 1.0e-12;

/// Immutable parameters of the `MAFFT-2002-FFT-NS-2-v1` scientific profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FftNs2Parameters {
    progressive: FftNs1Parameters,
    poisson_identity_floor: f64,
}

impl FftNs2Parameters {
    /// Returns the 2002 two-pass profile and its documented finite-distance resolution.
    #[must_use]
    pub const fn paper_2002() -> Self {
        Self {
            progressive: FftNs1Parameters::paper_2002(),
            poisson_identity_floor: FFT_NS2_POISSON_IDENTITY_FLOOR,
        }
    }

    /// Returns the published constants shared by both progressive passes.
    #[must_use]
    pub const fn progressive(self) -> FftNs1Parameters {
        self.progressive
    }

    /// Returns the explicit floor that keeps zero-identity Poisson distances finite.
    #[must_use]
    pub const fn poisson_identity_floor(self) -> f64 {
        self.poisson_identity_floor
    }
}

impl Default for FftNs2Parameters {
    fn default() -> Self {
        Self::paper_2002()
    }
}

/// Maximum accepted refinement cycles in the published FFT-NS-i profile.
pub const FFT_NSI_MAX_CYCLES: usize = 1_000;

/// Gotoh three-way equalization exponent fixed by the versioned profile.
pub const FFT_NSI_EQUALIZATION_FACTOR: f64 = 1.15;

/// Positive floor for zero-length branches in Gotoh's weighting equations.
pub const FFT_NSI_BRANCH_LENGTH_FLOOR: f64 = 1.0e-12;

/// Immutable parameters of the `MAFFT-2002-FFT-NS-i-v2` scientific profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FftNsiParameters {
    initial: FftNs2Parameters,
    maximum_cycles: usize,
    equalization_factor: f64,
    branch_length_floor: f64,
    thread_count: usize,
}

impl FftNsiParameters {
    /// Returns the 2002 iterative profile and its documented resolutions.
    #[must_use]
    pub const fn paper_2002() -> Self {
        Self {
            initial: FftNs2Parameters::paper_2002(),
            maximum_cycles: FFT_NSI_MAX_CYCLES,
            equalization_factor: FFT_NSI_EQUALIZATION_FACTOR,
            branch_length_floor: FFT_NSI_BRANCH_LENGTH_FLOOR,
            thread_count: 1,
        }
    }

    /// Returns the complete two-pass profile used to construct the initial MSA.
    #[must_use]
    pub const fn initial(self) -> FftNs2Parameters {
        self.initial
    }

    /// Returns the maximum number of strictly improving cycles.
    #[must_use]
    pub const fn maximum_cycles(self) -> usize {
        self.maximum_cycles
    }

    /// Returns Gotoh's three-way equalization exponent.
    #[must_use]
    pub const fn equalization_factor(self) -> f64 {
        self.equalization_factor
    }

    /// Returns the positive floor applied to inferred zero-length branches.
    #[must_use]
    pub const fn branch_length_floor(self) -> f64 {
        self.branch_length_floor
    }

    /// Returns the deterministic candidate-evaluation worker count.
    #[must_use]
    pub const fn thread_count(self) -> usize {
        self.thread_count
    }

    /// Replaces the deterministic candidate-evaluation worker count.
    ///
    /// This execution setting does not change the scientific profile or result
    /// ordering. Candidate results are reduced in canonical partition order.
    ///
    /// # Errors
    ///
    /// Returns [`crate::MafftError::InvalidThreadCount`] when `count` is zero.
    pub const fn with_thread_count(mut self, count: usize) -> Result<Self, crate::MafftError> {
        if count == 0 {
            return Err(crate::MafftError::InvalidThreadCount { value: count });
        }
        self.thread_count = count;
        Ok(self)
    }
}

impl Default for FftNsiParameters {
    fn default() -> Self {
        Self::paper_2002()
    }
}

/// Maximum accepted refinement cycles in the ancestral MAFFT 5.3 G-INS-i profile.
pub const GINSI_MAX_CYCLES: usize = 1_000;

/// Weight of the primary-library importance term in G-INS-i.
pub const GINSI_IMPORTANCE_WEIGHT: f64 = 2.7;

/// Arithmetic-linkage fraction in MAFFT's mixed guide-tree update.
pub const GINSI_MIXED_LINKAGE_FRACTION: f64 = 0.1;

/// Extension-free profile gap-opening cost used by MAFFT 5.3.
pub const GINSI_GAP_OPEN: f64 = 1.53;

/// Random-match offset added after matrix normalization in MAFFT 5.3.
pub const GINSI_RANDOM_MATCH_OFFSET: f64 = 0.123;

/// Immutable parameters of the MAFFT-2005-G-INS-i-v1 scientific profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GinsiParameters {
    maximum_cycles: usize,
    importance_weight: f64,
    mixed_linkage_fraction: f64,
    gap_open: f64,
    random_match_offset: f64,
    equalization_factor: f64,
    branch_length_floor: f64,
}

impl GinsiParameters {
    /// Returns the constants of the ancestral G-INS-i method in MAFFT 5.3.
    #[must_use]
    pub const fn paper_2005() -> Self {
        Self {
            maximum_cycles: GINSI_MAX_CYCLES,
            importance_weight: GINSI_IMPORTANCE_WEIGHT,
            mixed_linkage_fraction: GINSI_MIXED_LINKAGE_FRACTION,
            gap_open: GINSI_GAP_OPEN,
            random_match_offset: GINSI_RANDOM_MATCH_OFFSET,
            equalization_factor: FFT_NSI_EQUALIZATION_FACTOR,
            branch_length_floor: FFT_NSI_BRANCH_LENGTH_FLOOR,
        }
    }

    /// Returns the maximum number of strictly improving refinement cycles.
    #[must_use]
    pub const fn maximum_cycles(self) -> usize {
        self.maximum_cycles
    }

    /// Returns the multiplier of the primary-library importance term.
    #[must_use]
    pub const fn importance_weight(self) -> f64 {
        self.importance_weight
    }

    /// Returns the arithmetic component of MAFFT's mixed linkage.
    #[must_use]
    pub const fn mixed_linkage_fraction(self) -> f64 {
        self.mixed_linkage_fraction
    }

    /// Returns the extension-free profile gap-opening cost.
    #[must_use]
    pub const fn gap_open(self) -> f64 {
        self.gap_open
    }

    /// Returns the normalized-matrix random-match offset.
    #[must_use]
    pub const fn random_match_offset(self) -> f64 {
        self.random_match_offset
    }

    /// Returns Gotoh's three-way equalization exponent.
    #[must_use]
    pub const fn equalization_factor(self) -> f64 {
        self.equalization_factor
    }

    /// Returns the positive floor applied to inferred zero-length branches.
    #[must_use]
    pub const fn branch_length_floor(self) -> f64 {
        self.branch_length_floor
    }
}

impl Default for GinsiParameters {
    fn default() -> Self {
        Self::paper_2005()
    }
}

/// Maximum accepted refinement cycles in the ancestral L-INS-i profile.
pub const LINSI_MAX_CYCLES: usize = GINSI_MAX_CYCLES;

/// Weight of the primary local-library importance term in L-INS-i.
pub const LINSI_IMPORTANCE_WEIGHT: f64 = GINSI_IMPORTANCE_WEIGHT;

/// Arithmetic-linkage fraction in the L-INS-i guide-tree update.
pub const LINSI_MIXED_LINKAGE_FRACTION: f64 = GINSI_MIXED_LINKAGE_FRACTION;

/// Extension-free complete-profile gap-opening cost in L-INS-i.
pub const LINSI_GROUP_GAP_OPEN: f64 = GINSI_GAP_OPEN;

/// Complete-profile normalized-matrix offset in L-INS-i.
pub const LINSI_GROUP_RANDOM_MATCH_OFFSET: f64 = GINSI_RANDOM_MATCH_OFFSET;

/// Signed gap-open contribution in the local affine all-pairs stage.
pub const LINSI_LOCAL_GAP_OPEN: f64 = -2.0;

/// Positive normalized-matrix offset penalty in the local affine all-pairs stage.
///
/// The public MAFFT option is written as `--lep 0.1`, but the value is
/// subtracted from the normalized substitution matrix in the local recurrence.
pub const LINSI_LOCAL_OFFSET_PENALTY: f64 = 0.1;

/// Signed per-residue extension contribution after a local gap opens.
pub const LINSI_LOCAL_GAP_EXTENSION: f64 = -0.1;

/// Immutable parameters of the `MAFFT-2005-L-INS-i-v1` scientific profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinsiParameters {
    common: GinsiParameters,
    local_gap_open: f64,
    local_offset_penalty: f64,
    local_gap_extension: f64,
}

impl LinsiParameters {
    /// Returns the published constants and versioned resolutions of ancestral L-INS-i.
    #[must_use]
    pub const fn paper_2005() -> Self {
        Self {
            common: GinsiParameters::paper_2005(),
            local_gap_open: LINSI_LOCAL_GAP_OPEN,
            local_offset_penalty: LINSI_LOCAL_OFFSET_PENALTY,
            local_gap_extension: LINSI_LOCAL_GAP_EXTENSION,
        }
    }

    pub(crate) const fn common(self) -> GinsiParameters {
        self.common
    }

    /// Returns the maximum number of strictly improving refinement cycles.
    #[must_use]
    pub const fn maximum_cycles(self) -> usize {
        self.common.maximum_cycles()
    }

    /// Returns the multiplier of the primary local-library importance term.
    #[must_use]
    pub const fn importance_weight(self) -> f64 {
        self.common.importance_weight()
    }

    /// Returns the arithmetic component of MAFFT's mixed linkage.
    #[must_use]
    pub const fn mixed_linkage_fraction(self) -> f64 {
        self.common.mixed_linkage_fraction()
    }

    /// Returns the extension-free complete-profile gap-opening cost.
    #[must_use]
    pub const fn group_gap_open(self) -> f64 {
        self.common.gap_open()
    }

    /// Returns the complete-profile normalized-matrix offset.
    #[must_use]
    pub const fn group_random_match_offset(self) -> f64 {
        self.common.random_match_offset()
    }

    /// Returns the signed local affine gap-open contribution.
    #[must_use]
    pub const fn local_gap_open(self) -> f64 {
        self.local_gap_open
    }

    /// Returns the positive offset penalty subtracted from the local matrix.
    #[must_use]
    pub const fn local_offset_penalty(self) -> f64 {
        self.local_offset_penalty
    }

    /// Returns the signed local affine gap-extension contribution.
    #[must_use]
    pub const fn local_gap_extension(self) -> f64 {
        self.local_gap_extension
    }

    /// Returns Gotoh's fixed-tree equalization exponent.
    #[must_use]
    pub const fn equalization_factor(self) -> f64 {
        self.common.equalization_factor()
    }

    /// Returns the positive floor applied to inferred zero tree branches.
    #[must_use]
    pub const fn branch_length_floor(self) -> f64 {
        self.common.branch_length_floor()
    }
}

impl Default for LinsiParameters {
    fn default() -> Self {
        Self::paper_2005()
    }
}

/// Maximum accepted refinement cycles in the ancestral E-INS-i profile.
pub const EINSI_MAX_CYCLES: usize = GINSI_MAX_CYCLES;

/// Weight of the primary generalized-local library in E-INS-i.
pub const EINSI_IMPORTANCE_WEIGHT: f64 = GINSI_IMPORTANCE_WEIGHT;

/// Arithmetic-linkage fraction in the E-INS-i guide-tree update.
pub const EINSI_MIXED_LINKAGE_FRACTION: f64 = GINSI_MIXED_LINKAGE_FRACTION;

/// Extension-free complete-profile gap-opening cost in E-INS-i.
pub const EINSI_GROUP_GAP_OPEN: f64 = GINSI_GAP_OPEN;

/// Complete-profile normalized-matrix offset recommended for ancestral E-INS-i.
pub const EINSI_GROUP_RANDOM_MATCH_OFFSET: f64 = 0.0;

/// Signed ordinary gap-open contribution in the generalized-local pair stage.
pub const EINSI_LOCAL_GAP_OPEN: f64 = LINSI_LOCAL_GAP_OPEN;

/// Positive normalized-matrix penalty subtracted in the generalized-local stage.
pub const EINSI_LOCAL_OFFSET_PENALTY: f64 = LINSI_LOCAL_OFFSET_PENALTY;

/// Signed ordinary and unpaired-skip extension contribution.
pub const EINSI_LOCAL_GAP_EXTENSION: f64 = LINSI_LOCAL_GAP_EXTENSION;

/// Signed contribution for opening one generalized unaligned region.
pub const EINSI_SKIP_GAP_OPEN: f64 = -6.0;

/// Signed contribution for every residue pair consumed without homology.
pub const EINSI_PAIRED_SKIP_EXTENSION: f64 = 0.0;

/// Immutable parameters of the `MAFFT-2008-E-INS-i-v1` scientific profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EinsiParameters {
    common: GinsiParameters,
    local_gap_open: f64,
    local_offset_penalty: f64,
    local_gap_extension: f64,
    skip_gap_open: f64,
    paired_skip_extension: f64,
}

impl EinsiParameters {
    /// Returns the published constants and versioned resolutions of ancestral E-INS-i.
    #[must_use]
    pub const fn paper_2008() -> Self {
        Self {
            common: GinsiParameters {
                maximum_cycles: EINSI_MAX_CYCLES,
                importance_weight: EINSI_IMPORTANCE_WEIGHT,
                mixed_linkage_fraction: EINSI_MIXED_LINKAGE_FRACTION,
                gap_open: EINSI_GROUP_GAP_OPEN,
                random_match_offset: EINSI_GROUP_RANDOM_MATCH_OFFSET,
                equalization_factor: FFT_NSI_EQUALIZATION_FACTOR,
                branch_length_floor: FFT_NSI_BRANCH_LENGTH_FLOOR,
            },
            local_gap_open: EINSI_LOCAL_GAP_OPEN,
            local_offset_penalty: EINSI_LOCAL_OFFSET_PENALTY,
            local_gap_extension: EINSI_LOCAL_GAP_EXTENSION,
            skip_gap_open: EINSI_SKIP_GAP_OPEN,
            paired_skip_extension: EINSI_PAIRED_SKIP_EXTENSION,
        }
    }

    pub(crate) const fn common(self) -> GinsiParameters {
        self.common
    }

    /// Returns the maximum number of strictly improving refinement cycles.
    #[must_use]
    pub const fn maximum_cycles(self) -> usize {
        self.common.maximum_cycles()
    }

    /// Returns the multiplier of the primary generalized-local library.
    #[must_use]
    pub const fn importance_weight(self) -> f64 {
        self.common.importance_weight()
    }

    /// Returns the arithmetic component of MAFFT's mixed linkage.
    #[must_use]
    pub const fn mixed_linkage_fraction(self) -> f64 {
        self.common.mixed_linkage_fraction()
    }

    /// Returns the extension-free complete-profile gap-opening cost.
    #[must_use]
    pub const fn group_gap_open(self) -> f64 {
        self.common.gap_open()
    }

    /// Returns the complete-profile normalized-matrix offset.
    #[must_use]
    pub const fn group_random_match_offset(self) -> f64 {
        self.common.random_match_offset()
    }

    /// Returns the signed ordinary local gap-open contribution.
    #[must_use]
    pub const fn local_gap_open(self) -> f64 {
        self.local_gap_open
    }

    /// Returns the positive offset penalty subtracted from the local matrix.
    #[must_use]
    pub const fn local_offset_penalty(self) -> f64 {
        self.local_offset_penalty
    }

    /// Returns the ordinary and unpaired-skip extension contribution.
    #[must_use]
    pub const fn local_gap_extension(self) -> f64 {
        self.local_gap_extension
    }

    /// Returns the signed contribution for opening a generalized skip.
    #[must_use]
    pub const fn skip_gap_open(self) -> f64 {
        self.skip_gap_open
    }

    /// Returns the signed contribution for one paired skipped residue.
    #[must_use]
    pub const fn paired_skip_extension(self) -> f64 {
        self.paired_skip_extension
    }

    /// Returns Gotoh's fixed-tree equalization exponent.
    #[must_use]
    pub const fn equalization_factor(self) -> f64 {
        self.common.equalization_factor()
    }

    /// Returns the positive floor applied to inferred zero tree branches.
    #[must_use]
    pub const fn branch_length_floor(self) -> f64 {
        self.common.branch_length_floor()
    }
}

impl Default for EinsiParameters {
    fn default() -> Self {
        Self::paper_2008()
    }
}

#[cfg(test)]
mod tests {
    use super::{
        EINSI_GROUP_GAP_OPEN, EINSI_GROUP_RANDOM_MATCH_OFFSET, EINSI_IMPORTANCE_WEIGHT,
        EINSI_LOCAL_GAP_EXTENSION, EINSI_LOCAL_GAP_OPEN, EINSI_LOCAL_OFFSET_PENALTY,
        EINSI_MAX_CYCLES, EINSI_MIXED_LINKAGE_FRACTION, EINSI_PAIRED_SKIP_EXTENSION,
        EINSI_SKIP_GAP_OPEN, EinsiParameters, FFT_NS2_POISSON_IDENTITY_FLOOR,
        FFT_NSI_BRANCH_LENGTH_FLOOR, FFT_NSI_EQUALIZATION_FACTOR, FFT_NSI_MAX_CYCLES,
        FftNs1Parameters, FftNs2Parameters, FftNsiParameters, GINSI_GAP_OPEN,
        GINSI_IMPORTANCE_WEIGHT, GINSI_MAX_CYCLES, GINSI_MIXED_LINKAGE_FRACTION,
        GINSI_RANDOM_MATCH_OFFSET, GinsiParameters, LINSI_GROUP_GAP_OPEN,
        LINSI_GROUP_RANDOM_MATCH_OFFSET, LINSI_IMPORTANCE_WEIGHT, LINSI_LOCAL_GAP_EXTENSION,
        LINSI_LOCAL_GAP_OPEN, LINSI_LOCAL_OFFSET_PENALTY, LINSI_MAX_CYCLES,
        LINSI_MIXED_LINKAGE_FRACTION, LinsiParameters,
    };

    #[test]
    fn paper_profile_preserves_printed_constants() {
        let parameters = FftNs1Parameters::paper_2002();
        assert_eq!(parameters.tuple_length(), 6);
        assert_eq!(parameters.peak_count(), 20);
        assert_eq!(parameters.window_length(), 30);
        assert!((parameters.segment_threshold_per_column() - 0.7).abs() < f64::EPSILON);
        assert_eq!(parameters.maximum_segment_length(), 150);
        assert!((parameters.gap_open() - 2.4).abs() < f64::EPSILON);
        assert!((parameters.random_match_offset() - 0.06).abs() < f64::EPSILON);
    }

    #[test]
    fn two_pass_profile_reuses_published_progressive_constants() {
        let parameters = FftNs2Parameters::paper_2002();
        assert_eq!(parameters.progressive(), FftNs1Parameters::paper_2002());
        assert_eq!(
            parameters.poisson_identity_floor().to_bits(),
            FFT_NS2_POISSON_IDENTITY_FLOOR.to_bits()
        );
    }

    #[test]
    fn iterative_profile_preserves_versioned_constants() {
        let parameters = FftNsiParameters::paper_2002();
        assert_eq!(parameters.initial(), FftNs2Parameters::paper_2002());
        assert_eq!(parameters.maximum_cycles(), FFT_NSI_MAX_CYCLES);
        assert_eq!(
            parameters.equalization_factor().to_bits(),
            FFT_NSI_EQUALIZATION_FACTOR.to_bits()
        );
        assert_eq!(
            parameters.branch_length_floor().to_bits(),
            FFT_NSI_BRANCH_LENGTH_FLOOR.to_bits()
        );
        assert_eq!(parameters.thread_count(), 1);
        assert_eq!(parameters.with_thread_count(8).unwrap().thread_count(), 8);
        assert_eq!(
            parameters.with_thread_count(0),
            Err(crate::MafftError::InvalidThreadCount { value: 0 })
        );
    }

    #[test]
    fn ginsi_profile_preserves_mafft_53_constants() {
        let parameters = GinsiParameters::paper_2005();
        assert_eq!(parameters.maximum_cycles(), GINSI_MAX_CYCLES);
        assert_eq!(
            parameters.importance_weight().to_bits(),
            GINSI_IMPORTANCE_WEIGHT.to_bits()
        );
        assert_eq!(
            parameters.mixed_linkage_fraction().to_bits(),
            GINSI_MIXED_LINKAGE_FRACTION.to_bits()
        );
        assert_eq!(parameters.gap_open().to_bits(), GINSI_GAP_OPEN.to_bits());
        assert_eq!(
            parameters.random_match_offset().to_bits(),
            GINSI_RANDOM_MATCH_OFFSET.to_bits()
        );
        assert_eq!(
            parameters.equalization_factor().to_bits(),
            FFT_NSI_EQUALIZATION_FACTOR.to_bits()
        );
        assert_eq!(
            parameters.branch_length_floor().to_bits(),
            FFT_NSI_BRANCH_LENGTH_FLOOR.to_bits()
        );
        let pairwise = FftNs1Parameters::ginsi_2005_pairwise();
        assert_eq!(pairwise.gap_open().to_bits(), GINSI_GAP_OPEN.to_bits());
        assert_eq!(
            pairwise.random_match_offset().to_bits(),
            GINSI_RANDOM_MATCH_OFFSET.to_bits()
        );
    }

    #[test]
    fn linsi_profile_preserves_published_local_and_shared_constants() {
        let parameters = LinsiParameters::paper_2005();
        assert_eq!(parameters.maximum_cycles(), LINSI_MAX_CYCLES);
        assert_eq!(
            parameters.importance_weight().to_bits(),
            LINSI_IMPORTANCE_WEIGHT.to_bits()
        );
        assert_eq!(
            parameters.mixed_linkage_fraction().to_bits(),
            LINSI_MIXED_LINKAGE_FRACTION.to_bits()
        );
        assert_eq!(
            parameters.group_gap_open().to_bits(),
            LINSI_GROUP_GAP_OPEN.to_bits()
        );
        assert_eq!(
            parameters.group_random_match_offset().to_bits(),
            LINSI_GROUP_RANDOM_MATCH_OFFSET.to_bits()
        );
        assert_eq!(
            parameters.local_gap_open().to_bits(),
            LINSI_LOCAL_GAP_OPEN.to_bits()
        );
        assert_eq!(
            parameters.local_offset_penalty().to_bits(),
            LINSI_LOCAL_OFFSET_PENALTY.to_bits()
        );
        assert_eq!(
            parameters.local_gap_extension().to_bits(),
            LINSI_LOCAL_GAP_EXTENSION.to_bits()
        );
        assert_eq!(
            parameters.equalization_factor().to_bits(),
            FFT_NSI_EQUALIZATION_FACTOR.to_bits()
        );
        assert_eq!(
            parameters.branch_length_floor().to_bits(),
            FFT_NSI_BRANCH_LENGTH_FLOOR.to_bits()
        );
    }

    #[test]
    fn einsi_profile_preserves_ancestral_generalized_constants() {
        let parameters = EinsiParameters::paper_2008();
        assert_eq!(parameters.maximum_cycles(), EINSI_MAX_CYCLES);
        assert_eq!(
            parameters.importance_weight().to_bits(),
            EINSI_IMPORTANCE_WEIGHT.to_bits()
        );
        assert_eq!(
            parameters.mixed_linkage_fraction().to_bits(),
            EINSI_MIXED_LINKAGE_FRACTION.to_bits()
        );
        assert_eq!(
            parameters.group_gap_open().to_bits(),
            EINSI_GROUP_GAP_OPEN.to_bits()
        );
        assert_eq!(
            parameters.group_random_match_offset().to_bits(),
            EINSI_GROUP_RANDOM_MATCH_OFFSET.to_bits()
        );
        assert_eq!(
            parameters.local_gap_open().to_bits(),
            EINSI_LOCAL_GAP_OPEN.to_bits()
        );
        assert_eq!(
            parameters.local_offset_penalty().to_bits(),
            EINSI_LOCAL_OFFSET_PENALTY.to_bits()
        );
        assert_eq!(
            parameters.local_gap_extension().to_bits(),
            EINSI_LOCAL_GAP_EXTENSION.to_bits()
        );
        assert_eq!(
            parameters.skip_gap_open().to_bits(),
            EINSI_SKIP_GAP_OPEN.to_bits()
        );
        assert_eq!(
            parameters.paired_skip_extension().to_bits(),
            EINSI_PAIRED_SKIP_EXTENSION.to_bits()
        );
        assert_eq!(
            parameters.equalization_factor().to_bits(),
            FFT_NSI_EQUALIZATION_FACTOR.to_bits()
        );
        assert_eq!(
            parameters.branch_length_floor().to_bits(),
            FFT_NSI_BRANCH_LENGTH_FLOOR.to_bits()
        );
    }
}
