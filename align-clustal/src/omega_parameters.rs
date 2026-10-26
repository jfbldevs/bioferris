use crate::ClustalOmegaError;

/// Molecule-specific, versioned Clustal Omega reconstruction profile.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OmegaMoleculeProfile {
    /// Protein-only Clustal Omega 1.0.2 default pipeline.
    ProteinV1_0_2,
    /// Explicit nucleotide behavior from the stable 1.2.4 companion profile.
    NucleotideV1_2_4,
}

impl OmegaMoleculeProfile {
    /// Returns the stable scientific profile identifier.
    #[must_use]
    pub const fn identifier(self) -> &'static str {
        match self {
            Self::ProteinV1_0_2 => "CLUSTALO-1.0.2-default-protein-v1",
            Self::NucleotideV1_2_4 => "CLUSTALO-1.2.4-default-nucleotide-v1",
        }
    }

    pub(crate) const fn diagnostic_name(self) -> &'static str {
        match self {
            Self::ProteinV1_0_2 => "protein 1.0.2",
            Self::NucleotideV1_2_4 => "nucleotide 1.2.4",
        }
    }
}

/// Guide-tree construction selected for one alignment.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OmegaGuideTreeMode {
    /// Modified mBed, with an effective full matrix for at most one leaf cluster.
    Mbed,
    /// Materialize every real one-tuple pair distance before UPGMA.
    Full,
}

/// HMM-HMM dynamic-programming engine selection policy.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OmegaEngineMode {
    /// Select MAC when its checked workspace fits, otherwise Viterbi.
    Auto,
    /// Require posterior-probability Maximum Accuracy alignment.
    Mac,
    /// Require the five-state Viterbi alignment.
    Viterbi,
}

/// Validated parameters for the versioned Clustal Omega profiles.
#[derive(Clone, Debug, PartialEq)]
pub struct ClustalOmegaParameters {
    molecule: OmegaMoleculeProfile,
    guide_tree_mode: OmegaGuideTreeMode,
    engine_mode: OmegaEngineMode,
    one_tuple_gap_penalty: u32,
    maximum_cluster_size: usize,
    thread_count: usize,
    kmeans_seed: u64,
    mac_memory_limit_bytes: usize,
    mac_threshold: f64,
    pseudocount_a: f64,
    pseudocount_b: f64,
    pseudocount_c: f64,
    transition_admixture: f64,
    gap_open_admixture: f64,
    gap_extend_admixture: f64,
    column_score_shift: f64,
    delete_open_factor: f64,
    delete_extend_factor: f64,
    insert_open_factor: f64,
    insert_extend_factor: f64,
}

impl ClustalOmegaParameters {
    const DEFAULT_KMEANS_SEED: u64 = 0x6D42_4544_4F4D_4547;
    const TWO_GIBIBYTES: usize = 2 * 1024 * 1024 * 1024;
    const EIGHT_GIBIBYTES: usize = 8 * 1024 * 1024 * 1024;

    /// Returns `CLUSTALO-1.0.2-default-protein-v1` defaults.
    #[must_use]
    pub const fn protein_v1_0_2() -> Self {
        Self::defaults(OmegaMoleculeProfile::ProteinV1_0_2, Self::TWO_GIBIBYTES)
    }

    /// Returns `CLUSTALO-1.2.4-default-nucleotide-v1` defaults.
    #[must_use]
    pub const fn nucleotide_v1_2_4() -> Self {
        Self::defaults(
            OmegaMoleculeProfile::NucleotideV1_2_4,
            Self::EIGHT_GIBIBYTES,
        )
    }

    const fn defaults(molecule: OmegaMoleculeProfile, memory_limit: usize) -> Self {
        Self {
            molecule,
            guide_tree_mode: OmegaGuideTreeMode::Mbed,
            engine_mode: OmegaEngineMode::Auto,
            one_tuple_gap_penalty: 3,
            maximum_cluster_size: 100,
            thread_count: 1,
            kmeans_seed: Self::DEFAULT_KMEANS_SEED,
            mac_memory_limit_bytes: memory_limit,
            mac_threshold: 0.35,
            pseudocount_a: match molecule {
                OmegaMoleculeProfile::ProteinV1_0_2 => 1.0,
                OmegaMoleculeProfile::NucleotideV1_2_4 => 0.55,
            },
            pseudocount_b: 1.5,
            pseudocount_c: 1.0,
            transition_admixture: 1.0,
            gap_open_admixture: 0.15,
            gap_extend_admixture: 1.0,
            column_score_shift: -0.03,
            delete_open_factor: 0.6,
            delete_extend_factor: 0.6,
            insert_open_factor: 0.6,
            insert_extend_factor: 0.6,
        }
    }

    /// Returns the exact molecule/version profile.
    #[must_use]
    pub const fn molecule_profile(&self) -> OmegaMoleculeProfile {
        self.molecule
    }

    /// Returns the guide-tree policy.
    #[must_use]
    pub const fn guide_tree_mode(&self) -> OmegaGuideTreeMode {
        self.guide_tree_mode
    }

    /// Returns the HMM-HMM engine policy.
    #[must_use]
    pub const fn engine_mode(&self) -> OmegaEngineMode {
        self.engine_mode
    }

    /// Returns the fixed penalty for each internal one-tuple gap event.
    #[must_use]
    pub const fn one_tuple_gap_penalty(&self) -> u32 {
        self.one_tuple_gap_penalty
    }

    /// Returns the maximum number of records in an mBed leaf cluster.
    #[must_use]
    pub const fn maximum_cluster_size(&self) -> usize {
        self.maximum_cluster_size
    }

    /// Returns the deterministic progressive-alignment worker count.
    #[must_use]
    pub const fn thread_count(&self) -> usize {
        self.thread_count
    }

    /// Returns the fixed seed used by deterministic k-means++ initialization.
    #[must_use]
    pub const fn kmeans_seed(&self) -> u64 {
        self.kmeans_seed
    }

    /// Returns the checked MAC workspace limit in bytes.
    #[must_use]
    pub const fn mac_memory_limit_bytes(&self) -> usize {
        self.mac_memory_limit_bytes
    }

    /// Returns the MAC posterior greediness threshold.
    #[must_use]
    pub const fn mac_threshold(&self) -> f64 {
        self.mac_threshold
    }

    /// Returns `(a,b,c)` for diversity-dependent emission pseudocounts.
    #[must_use]
    pub const fn emission_pseudocounts(&self) -> (f64, f64, f64) {
        (self.pseudocount_a, self.pseudocount_b, self.pseudocount_c)
    }

    /// Returns `(gapb,gapd,gape)` transition-pseudocount parameters.
    #[must_use]
    pub const fn transition_pseudocounts(&self) -> (f64, f64, f64) {
        (
            self.transition_admixture,
            self.gap_open_admixture,
            self.gap_extend_admixture,
        )
    }

    /// Returns the additive HMM column-score offset in bits.
    #[must_use]
    pub const fn column_score_shift(&self) -> f64 {
        self.column_score_shift
    }

    /// Returns delete-open, delete-extend, insert-open, and insert-extend factors.
    #[must_use]
    pub const fn gap_path_factors(&self) -> (f64, f64, f64, f64) {
        (
            self.delete_open_factor,
            self.delete_extend_factor,
            self.insert_open_factor,
            self.insert_extend_factor,
        )
    }

    /// Replaces the guide-tree construction policy.
    #[must_use]
    pub const fn with_guide_tree_mode(mut self, mode: OmegaGuideTreeMode) -> Self {
        self.guide_tree_mode = mode;
        self
    }

    /// Replaces the HMM-HMM engine policy.
    #[must_use]
    pub const fn with_engine_mode(mut self, mode: OmegaEngineMode) -> Self {
        self.engine_mode = mode;
        self
    }

    /// Replaces the fixed internal-gap penalty in the one-tuple score.
    #[must_use]
    pub const fn with_one_tuple_gap_penalty(mut self, penalty: u32) -> Self {
        self.one_tuple_gap_penalty = penalty;
        self
    }

    /// Replaces the deterministic k-means++ seed.
    #[must_use]
    pub const fn with_kmeans_seed(mut self, seed: u64) -> Self {
        self.kmeans_seed = seed;
        self
    }

    /// Replaces the MAC workspace limit.
    #[must_use]
    pub const fn with_mac_memory_limit_bytes(mut self, bytes: usize) -> Self {
        self.mac_memory_limit_bytes = bytes;
        self
    }

    /// Replaces the mBed leaf-cluster cap.
    ///
    /// # Errors
    ///
    /// Returns [`ClustalOmegaError`] when `size` is below two.
    pub const fn with_maximum_cluster_size(
        mut self,
        size: usize,
    ) -> Result<Self, ClustalOmegaError> {
        if size < 2 {
            return Err(ClustalOmegaError::InvalidClusterSize { value: size });
        }
        self.maximum_cluster_size = size;
        Ok(self)
    }

    /// Replaces the deterministic progressive-alignment worker count.
    ///
    /// # Errors
    ///
    /// Returns [`ClustalOmegaError`] when `count` is zero.
    pub const fn with_thread_count(mut self, count: usize) -> Result<Self, ClustalOmegaError> {
        if count == 0 {
            return Err(ClustalOmegaError::InvalidThreadCount { value: count });
        }
        self.thread_count = count;
        Ok(self)
    }

    /// Replaces the MAC posterior threshold.
    ///
    /// # Errors
    ///
    /// Returns [`ClustalOmegaError`] unless `threshold` is finite and in `[0,1)`.
    pub fn with_mac_threshold(mut self, threshold: f64) -> Result<Self, ClustalOmegaError> {
        validate_finite_range("mac_threshold", threshold, 0.0, 1.0, false)?;
        self.mac_threshold = threshold;
        Ok(self)
    }

    /// Replaces diversity-dependent emission-pseudocount parameters.
    ///
    /// # Errors
    ///
    /// Returns [`ClustalOmegaError`] unless `a` is in `[0,1]`, `b >= 1`, and
    /// `c` is in `[0,3]`.
    pub fn with_emission_pseudocounts(
        mut self,
        a: f64,
        b: f64,
        c: f64,
    ) -> Result<Self, ClustalOmegaError> {
        validate_finite_range("pseudocount_a", a, 0.0, 1.0, true)?;
        validate_finite_minimum("pseudocount_b", b, 1.0)?;
        validate_finite_range("pseudocount_c", c, 0.0, 3.0, true)?;
        self.pseudocount_a = a;
        self.pseudocount_b = b;
        self.pseudocount_c = c;
        Ok(self)
    }

    /// Replaces transition-pseudocount parameters `(gapb,gapd,gape)`.
    ///
    /// # Errors
    ///
    /// Returns [`ClustalOmegaError`] unless all values are finite and nonnegative,
    /// `gapd < 17.48`, and `gape <= 1.5`.
    #[allow(clippy::similar_names)]
    pub fn with_transition_pseudocounts(
        mut self,
        gapb: f64,
        gapd: f64,
        gape: f64,
    ) -> Result<Self, ClustalOmegaError> {
        validate_finite_minimum("transition_admixture", gapb, 0.0)?;
        validate_finite_range("gap_open_admixture", gapd, 0.0, 17.48, false)?;
        validate_finite_range("gap_extend_admixture", gape, 0.0, 1.5, true)?;
        self.transition_admixture = gapb;
        self.gap_open_admixture = gapd;
        self.gap_extend_admixture = gape;
        Ok(self)
    }

    /// Replaces the HMM column-score shift in bits.
    ///
    /// # Errors
    ///
    /// Returns [`ClustalOmegaError`] unless the shift is finite and in `[-1,1]`.
    pub fn with_column_score_shift(mut self, shift: f64) -> Result<Self, ClustalOmegaError> {
        validate_finite_range("column_score_shift", shift, -1.0, 1.0, true)?;
        self.column_score_shift = shift;
        Ok(self)
    }

    /// Replaces the four gap-path log-score factors.
    ///
    /// # Errors
    ///
    /// Returns [`ClustalOmegaError`] unless every value is finite and positive.
    pub fn with_gap_path_factors(
        mut self,
        delete_open: f64,
        delete_extend: f64,
        insert_open: f64,
        insert_extend: f64,
    ) -> Result<Self, ClustalOmegaError> {
        for (name, value) in [
            ("delete_open_factor", delete_open),
            ("delete_extend_factor", delete_extend),
            ("insert_open_factor", insert_open),
            ("insert_extend_factor", insert_extend),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ClustalOmegaError::InvalidParameter {
                    parameter: name,
                    value,
                    requirement: "finite and positive",
                });
            }
        }
        self.delete_open_factor = delete_open;
        self.delete_extend_factor = delete_extend;
        self.insert_open_factor = insert_open;
        self.insert_extend_factor = insert_extend;
        Ok(self)
    }

    pub(crate) fn validate(&self) -> Result<(), ClustalOmegaError> {
        if self.maximum_cluster_size < 2 {
            return Err(ClustalOmegaError::InvalidClusterSize {
                value: self.maximum_cluster_size,
            });
        }
        if self.thread_count == 0 {
            return Err(ClustalOmegaError::InvalidThreadCount {
                value: self.thread_count,
            });
        }
        validate_finite_range("mac_threshold", self.mac_threshold, 0.0, 1.0, false)?;
        validate_finite_range("pseudocount_a", self.pseudocount_a, 0.0, 1.0, true)?;
        validate_finite_minimum("pseudocount_b", self.pseudocount_b, 1.0)?;
        validate_finite_range("pseudocount_c", self.pseudocount_c, 0.0, 3.0, true)?;
        validate_finite_minimum("transition_admixture", self.transition_admixture, 0.0)?;
        validate_finite_range(
            "gap_open_admixture",
            self.gap_open_admixture,
            0.0,
            17.48,
            false,
        )?;
        validate_finite_range(
            "gap_extend_admixture",
            self.gap_extend_admixture,
            0.0,
            1.5,
            true,
        )?;
        validate_finite_range(
            "column_score_shift",
            self.column_score_shift,
            -1.0,
            1.0,
            true,
        )?;
        for (name, value) in [
            ("delete_open_factor", self.delete_open_factor),
            ("delete_extend_factor", self.delete_extend_factor),
            ("insert_open_factor", self.insert_open_factor),
            ("insert_extend_factor", self.insert_extend_factor),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ClustalOmegaError::InvalidParameter {
                    parameter: name,
                    value,
                    requirement: "finite and positive",
                });
            }
        }
        Ok(())
    }
}

fn validate_finite_minimum(
    parameter: &'static str,
    value: f64,
    minimum: f64,
) -> Result<(), ClustalOmegaError> {
    if !value.is_finite() || value < minimum {
        return Err(ClustalOmegaError::InvalidParameter {
            parameter,
            value,
            requirement: "finite and above its documented minimum",
        });
    }
    Ok(())
}

fn validate_finite_range(
    parameter: &'static str,
    value: f64,
    minimum: f64,
    maximum: f64,
    inclusive_maximum: bool,
) -> Result<(), ClustalOmegaError> {
    let maximum_is_valid = if inclusive_maximum {
        value <= maximum
    } else {
        value < maximum
    };
    if !value.is_finite() || value < minimum || !maximum_is_valid {
        return Err(ClustalOmegaError::InvalidParameter {
            parameter,
            value,
            requirement: "finite and inside its documented range",
        });
    }
    Ok(())
}
