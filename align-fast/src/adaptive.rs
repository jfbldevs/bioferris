//! Deterministic cost model for exact affine backend selection.

use align_core::{Alphabet, Score, Sequence};

use crate::wavefront::exact_alignment_upper_bound_unbounded;
use crate::{
    WavefrontAffine, WavefrontAffinePenalties, WavefrontAlignment, WavefrontAlignmentView,
    WavefrontError, WavefrontWorkspace,
};

const CROSSOVER_MIN_LENGTH: usize = 256;
const CROSSOVER_MAX_LENGTH: usize = 8_192;

/// Exact affine backend selected by [`adaptive_affine_decision`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AdaptiveAffineBackend {
    /// The two indexed sequences are identical, so their exact result is trivial.
    Identity,
    /// Use exact affine wavefront alignment.
    Wavefront,
    /// Use exact dense Gotoh alignment, preferably through its SIMD kernel.
    Gotoh,
}

impl AdaptiveAffineBackend {
    /// Returns a stable provenance label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Identity => "identity",
            Self::Wavefront => "wavefront-affine",
            Self::Gotoh => "simd-gotoh",
        }
    }
}

/// Auditable inputs and output of the exact affine crossover cost model.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AdaptiveAffineDecision {
    backend: AdaptiveAffineBackend,
    estimated_wavefront_levels: u128,
    dense_cells: u128,
}

impl AdaptiveAffineDecision {
    /// Returns the selected exact backend.
    #[must_use]
    pub const fn backend(self) -> AdaptiveAffineBackend {
        self.backend
    }

    /// Returns the normalized affine upper-bound level estimate.
    #[must_use]
    pub const fn estimated_wavefront_levels(self) -> u128 {
        self.estimated_wavefront_levels
    }

    /// Returns `first_len * second_len`, the dense dynamic-programming cell count.
    #[must_use]
    pub const fn dense_cells(self) -> u128 {
        self.dense_cells
    }
}

/// One deterministic backend decision tied to the exact analyzed inputs.
///
/// Keeping the input borrow inside the plan allows WFA to reuse the selector's
/// exact upper bound without exposing a caller-supplied pruning bound.
pub struct AdaptiveAffinePlan<'a, A: Alphabet> {
    first: &'a Sequence<A>,
    second: &'a Sequence<A>,
    penalties: WavefrontAffinePenalties,
    decision: AdaptiveAffineDecision,
    wavefront_bound: Option<(Score, bool)>,
}

impl<A: Alphabet> AdaptiveAffinePlan<'_, A> {
    /// Returns the auditable backend decision.
    #[must_use]
    pub const fn decision(&self) -> AdaptiveAffineDecision {
        self.decision
    }

    /// Runs WFA score-only execution against the inputs bound to this plan.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`WavefrontAffine::penalty_with_workspace`].
    pub fn wavefront_penalty_with_workspace(
        &self,
        workspace: &mut WavefrontWorkspace,
    ) -> Result<Score, WavefrontError> {
        let wavefront = WavefrontAffine::new(self.penalties);
        if let Some((upper_bound, prune_main)) = self.wavefront_bound {
            wavefront.penalty_with_workspace_bound(
                self.first,
                self.second,
                workspace,
                upper_bound,
                prune_main,
            )
        } else {
            wavefront.penalty_with_workspace(self.first, self.second, workspace)
        }
    }

    /// Runs owned WFA traceback against the inputs bound to this plan.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`WavefrontAffine::align_with_workspace`].
    pub fn wavefront_align_with_workspace(
        &self,
        workspace: &mut WavefrontWorkspace,
    ) -> Result<WavefrontAlignment, WavefrontError> {
        let wavefront = WavefrontAffine::new(self.penalties);
        if let Some((upper_bound, _)) = self.wavefront_bound {
            wavefront.align_with_workspace_bound(self.first, self.second, workspace, upper_bound)
        } else {
            wavefront.align_with_workspace(self.first, self.second, workspace)
        }
    }

    /// Runs workspace-backed WFA traceback against the inputs bound to this plan.
    ///
    /// # Errors
    ///
    /// Returns the same failures as [`WavefrontAffine::align_in_workspace`].
    pub fn wavefront_align_in_workspace<'w>(
        &self,
        workspace: &'w mut WavefrontWorkspace,
    ) -> Result<WavefrontAlignmentView<'w>, WavefrontError> {
        let wavefront = WavefrontAffine::new(self.penalties);
        if let Some((upper_bound, _)) = self.wavefront_bound {
            wavefront.align_in_workspace_bound(self.first, self.second, workspace, upper_bound)
        } else {
            wavefront.align_in_workspace(self.first, self.second, workspace)
        }
    }
}

/// Analyzes and binds an exact affine backend plan to two input sequences.
///
/// Selection itself is deterministic. Identical inputs use the identity path. Other
/// inputs shorter than 256 symbols, or longer than the bounded 8,192-symbol scan,
/// remain on WFA. Within that interval, Gotoh is selected when the square of the
/// normalized WFA upper-bound estimate exceeds the dense matrix cell count.
#[must_use]
#[inline]
pub fn adaptive_affine_plan<'a, A: Alphabet>(
    first: &'a Sequence<A>,
    second: &'a Sequence<A>,
    penalties: WavefrontAffinePenalties,
) -> AdaptiveAffinePlan<'a, A> {
    let maximum_length = first.len().max(second.len());
    let dense_cells = (first.len() as u128) * (second.len() as u128);
    if first.as_indices() == second.as_indices() {
        return AdaptiveAffinePlan {
            first,
            second,
            penalties,
            decision: AdaptiveAffineDecision {
                backend: AdaptiveAffineBackend::Identity,
                estimated_wavefront_levels: 0,
                dense_cells,
            },
            wavefront_bound: Some((0, true)),
        };
    }
    if !(CROSSOVER_MIN_LENGTH..=CROSSOVER_MAX_LENGTH).contains(&maximum_length) || dense_cells == 0
    {
        return AdaptiveAffinePlan {
            first,
            second,
            penalties,
            decision: AdaptiveAffineDecision {
                backend: AdaptiveAffineBackend::Wavefront,
                estimated_wavefront_levels: 0,
                dense_cells,
            },
            wavefront_bound: None,
        };
    }

    let (upper_bound, prune_main) =
        exact_alignment_upper_bound_unbounded(first.as_indices(), second.as_indices(), penalties);
    let mismatch = u128::from(penalties.mismatch().unsigned_abs());
    let opening = u128::from(penalties.gap_opening().unsigned_abs());
    let extension = u128::from(penalties.gap_extension().unsigned_abs());
    let upper_bound_estimate = u128::try_from(upper_bound).unwrap_or(u128::MAX);
    let minimum_increment = mismatch
        .min(opening.saturating_add(extension))
        .min(extension);
    let estimated_wavefront_levels =
        upper_bound_estimate.saturating_add(minimum_increment - 1) / minimum_increment;
    let estimated_wavefront_cells =
        estimated_wavefront_levels.saturating_mul(estimated_wavefront_levels);
    AdaptiveAffinePlan {
        first,
        second,
        penalties,
        decision: AdaptiveAffineDecision {
            backend: if estimated_wavefront_cells > dense_cells {
                AdaptiveAffineBackend::Gotoh
            } else {
                AdaptiveAffineBackend::Wavefront
            },
            estimated_wavefront_levels,
            dense_cells,
        },
        wavefront_bound: Some((upper_bound, prune_main)),
    }
}

/// Selects an exact affine backend from sequence dimensions, contents, and WFA penalties.
#[must_use]
#[inline]
pub fn adaptive_affine_decision<A: Alphabet>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    penalties: WavefrontAffinePenalties,
) -> AdaptiveAffineDecision {
    adaptive_affine_plan(first, second, penalties).decision()
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, Sequence};

    use super::*;

    #[test]
    fn cost_model_selects_all_three_backends_deterministically() {
        let penalties = WavefrontAffinePenalties::new(4, 6, 2).unwrap();
        let identity = Sequence::<Dna>::from_ascii(vec![b'A'; 256]).unwrap();
        assert_eq!(
            adaptive_affine_decision(&identity, &identity, penalties).backend(),
            AdaptiveAffineBackend::Identity
        );

        let short_first = Sequence::<Dna>::from_ascii(b"ACGTACGT").unwrap();
        let short_second = Sequence::<Dna>::from_ascii(b"ACGTTTGT").unwrap();
        assert_eq!(
            adaptive_affine_decision(&short_first, &short_second, penalties).backend(),
            AdaptiveAffineBackend::Wavefront
        );

        let dense_first = Sequence::<Dna>::from_ascii(vec![b'A'; 256]).unwrap();
        let dense_second = Sequence::<Dna>::from_ascii(vec![b'C'; 256]).unwrap();
        let decision = adaptive_affine_decision(&dense_first, &dense_second, penalties);
        assert_eq!(decision.backend(), AdaptiveAffineBackend::Gotoh);
        assert!(decision.estimated_wavefront_levels() > 0);
        assert_eq!(decision.dense_cells(), 256 * 256);
    }

    #[test]
    fn bound_reusing_plan_matches_the_regular_wavefront_paths() {
        let penalties = WavefrontAffinePenalties::new(4, 6, 2).unwrap();
        let first = Sequence::<Dna>::from_ascii(vec![b'A'; 256]).unwrap();
        let mut second_symbols = vec![b'A'; 256];
        for index in [17, 71, 149, 233] {
            second_symbols[index] = b'C';
        }
        let second = Sequence::<Dna>::from_ascii(second_symbols).unwrap();
        let plan = adaptive_affine_plan(&first, &second, penalties);
        assert_eq!(plan.decision().backend(), AdaptiveAffineBackend::Wavefront);

        let wavefront = WavefrontAffine::new(penalties);
        let expected_penalty = wavefront.penalty(&first, &second).unwrap();
        let expected_alignment = wavefront.align(&first, &second).unwrap();
        let mut workspace = WavefrontWorkspace::new();

        assert_eq!(
            plan.wavefront_penalty_with_workspace(&mut workspace)
                .unwrap(),
            expected_penalty
        );
        let observed = plan.wavefront_align_in_workspace(&mut workspace).unwrap();
        assert_eq!(observed.penalty(), expected_alignment.penalty());
        assert_eq!(
            observed.operations(),
            expected_alignment.alignment().operations()
        );

        let mut third_symbols = vec![b'A'; 256];
        third_symbols[127] = b'G';
        let third = Sequence::<Dna>::from_ascii(third_symbols).unwrap();
        let expected_next = wavefront.penalty(&first, &third).unwrap();
        assert_eq!(
            wavefront
                .penalty_with_workspace(&first, &third, &mut workspace)
                .unwrap(),
            expected_next,
            "a consumed plan bound must not leak into the next workspace call"
        );
    }
}
