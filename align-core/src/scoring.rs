use core::marker::PhantomData;

use crate::{Alphabet, ScoringError, Symbol};

/// Integer score used throughout the alignment workspace.
pub type Score = i64;

/// Scores aligned symbol pairs from one alphabet.
pub trait SubstitutionScorer<A: Alphabet>: Send + Sync {
    /// Returns the score for aligning `first` with `second`.
    fn score(&self, first: Symbol<A>, second: Symbol<A>) -> Score;

    /// Optionally fills one complete score row in [`Alphabet::SYMBOLS`] order.
    ///
    /// Implementations return `true` only after writing exactly `output.len()`
    /// values such that `output[second.alphabet_index()]` equals
    /// `self.score(first, second)` for every symbol in the alphabet. Callers pass
    /// a buffer whose length equals the alphabet size. Returning `false` selects
    /// scalar calls to [`Self::score`] and leaves `output` unspecified.
    ///
    /// This opt-in hook lets dynamically dispatched table scorers remove repeated
    /// dispatch and row-address arithmetic from dynamic-programming inner loops.
    fn fill_score_row(&self, first: Symbol<A>, output: &mut [Score]) -> bool {
        let _ = (first, output);
        false
    }

    /// Returns inclusive lower and upper bounds for every score this scorer emits.
    ///
    /// Returning `None` selects fully checked arithmetic in algorithms. Implementors
    /// may return bounds to enable a preflight-proven fast path, but every subsequent
    /// result from [`Self::score`] must remain inside the declared range.
    #[must_use]
    fn score_bounds(&self) -> Option<(Score, Score)> {
        None
    }
}

impl<A, S> SubstitutionScorer<A> for &S
where
    A: Alphabet,
    S: SubstitutionScorer<A> + ?Sized,
{
    fn score(&self, first: Symbol<A>, second: Symbol<A>) -> Score {
        (**self).score(first, second)
    }

    fn fill_score_row(&self, first: Symbol<A>, output: &mut [Score]) -> bool {
        (**self).fill_score_row(first, output)
    }

    fn score_bounds(&self) -> Option<(Score, Score)> {
        (**self).score_bounds()
    }
}

/// Exact-symbol match and mismatch scoring.
///
/// Ambiguity symbols are compared literally. For example, `N` matches `N`, while
/// `N` mismatches `A`; callers requiring set-overlap semantics should provide a
/// [`DenseSubstitutionMatrix`].
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MatchMismatch {
    match_score: Score,
    mismatch_score: Score,
}

impl MatchMismatch {
    /// Creates literal match/mismatch scoring.
    #[must_use]
    pub const fn new(match_score: Score, mismatch_score: Score) -> Self {
        Self {
            match_score,
            mismatch_score,
        }
    }

    /// Returns the exact-symbol match score.
    #[must_use]
    pub const fn match_score(self) -> Score {
        self.match_score
    }

    /// Returns the non-identical-symbol score.
    #[must_use]
    pub const fn mismatch_score(self) -> Score {
        self.mismatch_score
    }
}

impl<A: Alphabet> SubstitutionScorer<A> for MatchMismatch {
    fn score(&self, first: Symbol<A>, second: Symbol<A>) -> Score {
        if first == second {
            self.match_score
        } else {
            self.mismatch_score
        }
    }

    fn score_bounds(&self) -> Option<(Score, Score)> {
        Some((
            self.match_score.min(self.mismatch_score),
            self.match_score.max(self.mismatch_score),
        ))
    }
}

/// A complete row-major substitution matrix for alphabet `A`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DenseSubstitutionMatrix<A: Alphabet> {
    scores: Box<[Score]>,
    bounds: Option<(Score, Score)>,
    marker: PhantomData<A>,
}

impl<A: Alphabet> DenseSubstitutionMatrix<A> {
    /// Builds a matrix from row-major scores in [`Alphabet::SYMBOLS`] order.
    ///
    /// # Errors
    ///
    /// Returns [`ScoringError::MatrixSizeMismatch`] unless exactly one score is
    /// supplied for every ordered symbol pair.
    pub fn from_flat_scores(scores: Vec<Score>) -> Result<Self, ScoringError> {
        let width = A::SYMBOLS.len();
        let Some(expected) = width.checked_mul(width) else {
            return Err(ScoringError::MatrixSizeMismatch {
                alphabet: A::NAME,
                expected: usize::MAX,
                actual: scores.len(),
            });
        };
        if scores.len() != expected {
            return Err(ScoringError::MatrixSizeMismatch {
                alphabet: A::NAME,
                expected,
                actual: scores.len(),
            });
        }
        let bounds = score_bounds(&scores);
        Ok(Self {
            scores: scores.into_boxed_slice(),
            bounds,
            marker: PhantomData,
        })
    }

    /// Builds a complete matrix by evaluating `score` in stable alphabet order.
    #[must_use]
    pub fn from_score_fn(mut score: impl FnMut(u8, u8) -> Score) -> Self {
        let width = A::SYMBOLS.len();
        let mut scores = Vec::with_capacity(width * width);
        for first in A::SYMBOLS {
            for second in A::SYMBOLS {
                scores.push(score(*first, *second));
            }
        }
        let bounds = score_bounds(&scores);
        Self {
            scores: scores.into_boxed_slice(),
            bounds,
            marker: PhantomData,
        }
    }

    /// Returns the matrix width and height.
    #[must_use]
    pub fn dimension(&self) -> usize {
        A::SYMBOLS.len()
    }

    /// Returns all row-major scores in stable alphabet order.
    #[must_use]
    pub fn as_flat_scores(&self) -> &[Score] {
        &self.scores
    }

    /// Copies one matrix row into `output` in [`Alphabet::SYMBOLS`] order.
    ///
    /// Returns `false` without writing when `output` does not have exactly one
    /// entry per alphabet symbol. This is intended for scorer adapters that opt
    /// into [`SubstitutionScorer::fill_score_row`]; direct matrix scoring remains
    /// on the compiler-optimized scalar lookup path.
    pub fn copy_score_row(&self, first: Symbol<A>, output: &mut [Score]) -> bool {
        let width = A::SYMBOLS.len();
        if output.len() != width {
            return false;
        }
        let row_start = first.alphabet_index() * width;
        output.copy_from_slice(&self.scores[row_start..row_start + width]);
        true
    }
}

impl<A: Alphabet> SubstitutionScorer<A> for DenseSubstitutionMatrix<A> {
    fn score(&self, first: Symbol<A>, second: Symbol<A>) -> Score {
        let width = A::SYMBOLS.len();
        self.scores[first.alphabet_index() * width + second.alphabet_index()]
    }

    fn score_bounds(&self) -> Option<(Score, Score)> {
        self.bounds
    }
}

fn score_bounds(scores: &[Score]) -> Option<(Score, Score)> {
    let minimum = scores.iter().copied().min()?;
    let maximum = scores.iter().copied().max()?;
    Some((minimum, maximum))
}

/// A non-positive score applied once per emitted gap symbol.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LinearGap(Score);

impl LinearGap {
    /// Validates a linear gap score.
    ///
    /// # Errors
    ///
    /// Returns [`ScoringError::PositiveGapScore`] when `score` is positive.
    pub const fn new(score: Score) -> Result<Self, ScoringError> {
        if score > 0 {
            Err(ScoringError::PositiveGapScore { score })
        } else {
            Ok(Self(score))
        }
    }

    /// Returns the score applied to each gap symbol.
    #[must_use]
    pub const fn score(self) -> Score {
        self.0
    }
}

/// A non-positive affine score for each contiguous gap run.
///
/// A gap of positive length `k` scores
/// `open_score + (k - 1) * extension_score`. Requiring
/// `open_score <= extension_score <= 0` preserves the non-negative `u` and `v`
/// restrictions in Gotoh's original distance formulation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct AffineGap {
    open_score: Score,
    extension_score: Score,
}

impl AffineGap {
    /// Validates affine gap scores using the length-one-open convention.
    ///
    /// # Errors
    ///
    /// Returns [`ScoringError`] if either score is positive or if opening a gap is
    /// less penalizing than extending it.
    pub const fn new(open_score: Score, extension_score: Score) -> Result<Self, ScoringError> {
        if open_score > 0 {
            return Err(ScoringError::PositiveAffineGapOpenScore { score: open_score });
        }
        if extension_score > 0 {
            return Err(ScoringError::PositiveAffineGapExtensionScore {
                score: extension_score,
            });
        }
        if open_score > extension_score {
            return Err(ScoringError::AffineGapOpenLessPenalizingThanExtension {
                open_score,
                extension_score,
            });
        }
        Ok(Self {
            open_score,
            extension_score,
        })
    }

    /// Returns the score for the first column in a gap run.
    #[must_use]
    pub const fn open_score(self) -> Score {
        self.open_score
    }

    /// Returns the score for each gap column after the first.
    #[must_use]
    pub const fn extension_score(self) -> Score {
        self.extension_score
    }

    /// Computes the score of one positive-length gap without overflow.
    ///
    /// Returns `None` for length zero or when the exact score does not fit in
    /// [`Score`].
    #[must_use]
    pub fn score_for_length(self, length: usize) -> Option<Score> {
        let extension_count = length.checked_sub(1)?;
        if self.extension_score == 0 {
            return Some(self.open_score);
        }
        let extension_count = Score::try_from(extension_count).ok()?;
        self.extension_score
            .checked_mul(extension_count)?
            .checked_add(self.open_score)
    }
}

#[cfg(test)]
mod tests {
    use crate::{Alphabet, Dna, Sequence, SubstitutionScorer};

    use super::{
        AffineGap, DenseSubstitutionMatrix, LinearGap, MatchMismatch, Score, ScoringError,
    };

    #[test]
    fn positive_gap_score_is_rejected() {
        assert_eq!(
            LinearGap::new(1),
            Err(ScoringError::PositiveGapScore { score: 1 })
        );
        assert_eq!(LinearGap::new(0).expect("zero is valid").score(), 0);
    }

    #[test]
    fn affine_gap_validates_and_scores_complete_runs() {
        let gap = AffineGap::new(-5, -1).expect("valid affine gap");
        assert_eq!(gap.open_score(), -5);
        assert_eq!(gap.extension_score(), -1);
        assert_eq!(gap.score_for_length(0), None);
        assert_eq!(gap.score_for_length(1), Some(-5));
        assert_eq!(gap.score_for_length(4), Some(-8));
        assert_eq!(
            AffineGap::new(-7, 0)
                .expect("zero extension is valid")
                .score_for_length(usize::MAX),
            Some(-7)
        );

        assert_eq!(
            AffineGap::new(1, 0),
            Err(ScoringError::PositiveAffineGapOpenScore { score: 1 })
        );
        assert_eq!(
            AffineGap::new(-1, 1),
            Err(ScoringError::PositiveAffineGapExtensionScore { score: 1 })
        );
        assert_eq!(
            AffineGap::new(-1, -2),
            Err(ScoringError::AffineGapOpenLessPenalizingThanExtension {
                open_score: -1,
                extension_score: -2,
            })
        );
    }

    #[test]
    fn dense_matrix_checks_shape_and_scores_by_index() {
        let error = DenseSubstitutionMatrix::<Dna>::from_flat_scores(vec![0; 3])
            .expect_err("wrong shape must fail");
        assert_eq!(
            error,
            ScoringError::MatrixSizeMismatch {
                alphabet: "DNA",
                expected: 225,
                actual: 3,
            }
        );

        let matrix =
            DenseSubstitutionMatrix::<Dna>::from_score_fn(|a, b| if a == b { 7 } else { -3 });
        let sequence = Sequence::<Dna>::from_ascii(b"AC").expect("valid DNA");
        assert_eq!(
            matrix.score(sequence.get(0).unwrap(), sequence.get(0).unwrap()),
            7
        );
        assert_eq!(
            matrix.score(sequence.get(0).unwrap(), sequence.get(1).unwrap()),
            -3
        );
        assert_eq!(matrix.score_bounds(), Some((-3, 7)));
    }

    #[test]
    fn dense_matrix_row_copy_matches_every_scalar_lookup() {
        let matrix =
            DenseSubstitutionMatrix::<Dna>::from_score_fn(|a, b| i64::from(a) - i64::from(b));
        let sequence = Sequence::<Dna>::from_ascii(Dna::SYMBOLS).expect("all DNA symbols");
        let mut row = vec![Score::MAX; Dna::SYMBOLS.len()];
        for first in sequence.as_slice() {
            assert!(matrix.copy_score_row(*first, &mut row));
            for second in sequence.as_slice() {
                assert_eq!(row[second.alphabet_index()], matrix.score(*first, *second));
            }
        }
        let expected_first_value = row[0];
        assert!(!matrix.copy_score_row(sequence.get(0).unwrap(), &mut row[..1]));
        assert_eq!(row[0], expected_first_value);
        assert!(!matrix.fill_score_row(sequence.get(0).unwrap(), &mut row));
    }

    #[test]
    fn ambiguity_symbols_are_literal_by_default() {
        let scorer = MatchMismatch::new(2, -1);
        let sequence = Sequence::<Dna>::from_ascii(b"NAN").expect("valid DNA");
        assert_eq!(
            scorer.score(sequence.get(0).unwrap(), sequence.get(2).unwrap()),
            2
        );
        assert_eq!(
            scorer.score(sequence.get(0).unwrap(), sequence.get(1).unwrap()),
            -1
        );
        assert_eq!(
            <MatchMismatch as SubstitutionScorer<Dna>>::score_bounds(&scorer),
            Some((-1, 2))
        );
    }
}
