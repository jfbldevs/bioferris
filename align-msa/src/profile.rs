use align_core::Alphabet;

use crate::{MsaError, MultipleAlignment};

/// A validated multiple alignment with one positive fixed-point weight per row.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Profile<A: Alphabet> {
    alignment: MultipleAlignment<A>,
    weights: Box<[u32]>,
}

impl<A: Alphabet> Profile<A> {
    /// Creates a weighted alignment profile.
    ///
    /// Weights are algorithm-defined fixed-point values. Their common scale does
    /// not affect weighted means, but every value must be positive.
    ///
    /// # Errors
    ///
    /// Returns [`MsaError`] for a count mismatch or zero weight.
    pub fn try_new(alignment: MultipleAlignment<A>, weights: Vec<u32>) -> Result<Self, MsaError> {
        if alignment.sequence_count() != weights.len() {
            return Err(MsaError::WeightCountMismatch {
                rows: alignment.sequence_count(),
                weights: weights.len(),
            });
        }
        if let Some(index) = weights.iter().position(|weight| *weight == 0) {
            return Err(MsaError::ZeroWeight { index });
        }
        Ok(Self {
            alignment,
            weights: weights.into_boxed_slice(),
        })
    }

    /// Returns the underlying alignment.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        &self.alignment
    }

    /// Returns fixed-point sequence weights in row order.
    #[must_use]
    pub const fn weights(&self) -> &[u32] {
        &self.weights
    }

    /// Separates the alignment and weights.
    #[must_use]
    pub fn into_parts(self) -> (MultipleAlignment<A>, Vec<u32>) {
        (self.alignment, self.weights.into_vec())
    }
}
