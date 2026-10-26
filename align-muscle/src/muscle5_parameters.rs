use crate::Muscle5Error;

/// MUSCLE 5 top-level alignment pipeline.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Muscle5Mode {
    /// Parallel Perturbed `ProbCons` for small and medium sequence sets.
    Align,
    /// Published divide-and-conquer pipeline for large sequence sets.
    Super5,
}

/// One of the four published MUSCLE 5 guide-tree resolutions.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GuideTreePermutation {
    /// Preserve the expected-error UPGMA topology.
    None,
    /// Resolve each three-way split as `((a,b),c)`.
    Abc,
    /// Resolve each three-way split as `((a,c),b)`.
    Acb,
    /// Resolve each three-way split as `((b,c),a)`.
    Bca,
}

impl GuideTreePermutation {
    pub(crate) const ALL: [Self; 4] = [Self::None, Self::Abc, Self::Acb, Self::Bca];

    /// Returns the lower-case label used by MUSCLE ensemble replicate names.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Abc => "abc",
            Self::Acb => "acb",
            Self::Bca => "bca",
        }
    }
}

/// Validated parameters for `MUSCLE-5.1.0-Nature-2022-default-v2`, unchanged from v1.
#[derive(Clone, Debug, PartialEq)]
pub struct Muscle5Parameters {
    mode: Muscle5Mode,
    consistency_iterations: usize,
    refinement_iterations: usize,
    perturbation_seed: u64,
    guide_tree_permutation: GuideTreePermutation,
    perturbation_amplitude: f64,
    posterior_cutoff: f64,
    redundancy_ee_threshold: f64,
    coarse_ee_threshold: f64,
    subdivision_ee_threshold: f64,
    maximum_cluster_size: usize,
    maximum_uclust_rejections: usize,
    ppp_pair_limit: usize,
}

impl Muscle5Parameters {
    /// Creates the publication-default PPP `Align` profile.
    #[must_use]
    pub const fn published_align() -> Self {
        Self::published(Muscle5Mode::Align)
    }

    /// Creates the publication-default `Super5` profile.
    #[must_use]
    pub const fn published_super5() -> Self {
        Self::published(Muscle5Mode::Super5)
    }

    const fn published(mode: Muscle5Mode) -> Self {
        Self {
            mode,
            consistency_iterations: 2,
            refinement_iterations: 100,
            perturbation_seed: 0,
            guide_tree_permutation: GuideTreePermutation::None,
            perturbation_amplitude: 0.25,
            posterior_cutoff: 0.01,
            redundancy_ee_threshold: 0.01,
            coarse_ee_threshold: 0.30,
            subdivision_ee_threshold: 0.10,
            maximum_cluster_size: 500,
            maximum_uclust_rejections: 16,
            ppp_pair_limit: 2_000,
        }
    }

    /// Returns the selected top-level pipeline.
    #[must_use]
    pub const fn mode(&self) -> Muscle5Mode {
        self.mode
    }

    /// Returns the number of probabilistic-consistency transformations.
    #[must_use]
    pub const fn consistency_iterations(&self) -> usize {
        self.consistency_iterations
    }

    /// Returns the number of randomized partition-refinement attempts.
    #[must_use]
    pub const fn refinement_iterations(&self) -> usize {
        self.refinement_iterations
    }

    /// Returns the HMM perturbation seed; zero selects the unperturbed model.
    #[must_use]
    pub const fn perturbation_seed(&self) -> u64 {
        self.perturbation_seed
    }

    /// Returns the selected guide-tree permutation.
    #[must_use]
    pub const fn guide_tree_permutation(&self) -> GuideTreePermutation {
        self.guide_tree_permutation
    }

    /// Returns the published multiplicative HMM perturbation amplitude.
    #[must_use]
    pub const fn perturbation_amplitude(&self) -> f64 {
        self.perturbation_amplitude
    }

    /// Returns the sparse posterior retention threshold.
    #[must_use]
    pub const fn posterior_cutoff(&self) -> f64 {
        self.posterior_cutoff
    }

    /// Returns the Super5 redundancy expected-error threshold.
    #[must_use]
    pub const fn redundancy_ee_threshold(&self) -> f64 {
        self.redundancy_ee_threshold
    }

    /// Returns the Super5 coarse-cluster expected-error threshold.
    #[must_use]
    pub const fn coarse_ee_threshold(&self) -> f64 {
        self.coarse_ee_threshold
    }

    /// Returns the Super5 large-cluster subdivision expected-error threshold.
    #[must_use]
    pub const fn subdivision_ee_threshold(&self) -> f64 {
        self.subdivision_ee_threshold
    }

    /// Returns the maximum number of rows in a PPP cluster.
    #[must_use]
    pub const fn maximum_cluster_size(&self) -> usize {
        self.maximum_cluster_size
    }

    /// Returns the maximum rejected UCLUST-EE candidates per query.
    #[must_use]
    pub const fn maximum_uclust_rejections(&self) -> usize {
        self.maximum_uclust_rejections
    }

    /// Returns the PPP-pa cross-pair limit and PPP-ps sample size.
    #[must_use]
    pub const fn ppp_pair_limit(&self) -> usize {
        self.ppp_pair_limit
    }

    /// Selects PPP `Align` or `Super5` while retaining all other values.
    #[must_use]
    pub const fn with_mode(mut self, mode: Muscle5Mode) -> Self {
        self.mode = mode;
        self
    }

    /// Selects an HMM perturbation seed.
    #[must_use]
    pub const fn with_perturbation_seed(mut self, seed: u64) -> Self {
        self.perturbation_seed = seed;
        self
    }

    /// Selects one published guide-tree resolution.
    #[must_use]
    pub const fn with_guide_tree_permutation(mut self, permutation: GuideTreePermutation) -> Self {
        self.guide_tree_permutation = permutation;
        self
    }

    /// Changes consistency and refinement iteration counts.
    ///
    /// Zero is valid and is useful for controlled ablation tests.
    #[must_use]
    pub const fn with_iterations(
        mut self,
        consistency_iterations: usize,
        refinement_iterations: usize,
    ) -> Self {
        self.consistency_iterations = consistency_iterations;
        self.refinement_iterations = refinement_iterations;
        self
    }

    /// Changes the sparse posterior cutoff.
    ///
    /// # Errors
    ///
    /// Returns [`Muscle5Error`] unless `cutoff` is finite and in `[0,1)`.
    pub fn with_posterior_cutoff(mut self, cutoff: f64) -> Result<Self, Muscle5Error> {
        self.posterior_cutoff = cutoff;
        self.validate()?;
        Ok(self)
    }

    /// Changes all documented Super5 thresholds and hard limits.
    ///
    /// # Errors
    ///
    /// Returns [`Muscle5Error`] for non-finite/out-of-order thresholds or zero
    /// hard limits.
    pub fn with_super5_limits(
        mut self,
        redundancy_ee_threshold: f64,
        coarse_ee_threshold: f64,
        subdivision_ee_threshold: f64,
        maximum_cluster_size: usize,
        maximum_uclust_rejections: usize,
        ppp_pair_limit: usize,
    ) -> Result<Self, Muscle5Error> {
        self.redundancy_ee_threshold = redundancy_ee_threshold;
        self.coarse_ee_threshold = coarse_ee_threshold;
        self.subdivision_ee_threshold = subdivision_ee_threshold;
        self.maximum_cluster_size = maximum_cluster_size;
        self.maximum_uclust_rejections = maximum_uclust_rejections;
        self.ppp_pair_limit = ppp_pair_limit;
        self.validate()?;
        Ok(self)
    }

    pub(crate) fn validate(&self) -> Result<(), Muscle5Error> {
        if !self.perturbation_amplitude.is_finite()
            || !(0.0..1.0).contains(&self.perturbation_amplitude)
        {
            return Err(invalid("perturbation_amplitude", "be finite and in [0,1)"));
        }
        if !self.posterior_cutoff.is_finite() || !(0.0..1.0).contains(&self.posterior_cutoff) {
            return Err(invalid("posterior_cutoff", "be finite and in [0,1)"));
        }
        for (name, value) in [
            ("redundancy_ee_threshold", self.redundancy_ee_threshold),
            ("subdivision_ee_threshold", self.subdivision_ee_threshold),
            ("coarse_ee_threshold", self.coarse_ee_threshold),
        ] {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                return Err(invalid(name, "be finite and in [0,1]"));
            }
        }
        if !(self.redundancy_ee_threshold <= self.subdivision_ee_threshold
            && self.subdivision_ee_threshold <= self.coarse_ee_threshold)
        {
            return Err(invalid(
                "Super5 EE thresholds",
                "satisfy redundancy <= subdivision <= coarse",
            ));
        }
        for (name, value) in [
            ("maximum_cluster_size", self.maximum_cluster_size),
            ("maximum_uclust_rejections", self.maximum_uclust_rejections),
            ("ppp_pair_limit", self.ppp_pair_limit),
        ] {
            if value == 0 {
                return Err(invalid(name, "be at least one"));
            }
        }
        Ok(())
    }
}

const fn invalid(parameter: &'static str, requirement: &'static str) -> Muscle5Error {
    Muscle5Error::InvalidParameter {
        parameter,
        requirement,
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::float_cmp)]

    use super::{GuideTreePermutation, Muscle5Mode, Muscle5Parameters};

    #[test]
    fn published_profiles_expose_paper_defaults() {
        let align = Muscle5Parameters::published_align();
        assert_eq!(align.mode(), Muscle5Mode::Align);
        assert_eq!(align.consistency_iterations(), 2);
        assert_eq!(align.refinement_iterations(), 100);
        assert_eq!(align.perturbation_seed(), 0);
        assert_eq!(align.guide_tree_permutation(), GuideTreePermutation::None);
        assert_eq!(align.perturbation_amplitude(), 0.25);
        assert_eq!(align.posterior_cutoff(), 0.01);
        assert_eq!(align.maximum_cluster_size(), 500);
        assert_eq!(align.maximum_uclust_rejections(), 16);
        assert_eq!(align.ppp_pair_limit(), 2_000);
        assert_eq!(
            Muscle5Parameters::published_super5().mode(),
            Muscle5Mode::Super5
        );
    }

    #[test]
    fn invalid_limits_are_rejected() {
        assert!(
            Muscle5Parameters::published_align()
                .with_posterior_cutoff(1.0)
                .is_err()
        );
        assert!(
            Muscle5Parameters::published_super5()
                .with_super5_limits(0.2, 0.1, 0.15, 500, 16, 2_000)
                .is_err()
        );
    }
}
