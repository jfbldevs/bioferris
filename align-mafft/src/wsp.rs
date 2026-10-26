use core::mem::size_of;

use align_core::Symbol;
use align_msa::MultipleAlignment;

use crate::matrix::ScientificScorer;
use crate::{FftNsiPairWeights, FftNsiPartition, MafftAlphabet, MafftError};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum GapOrientation {
    FirstResidue,
    SecondResidue,
}

pub(crate) struct WspPairScoreCache {
    size: usize,
    weighted_scores: Box<[f64]>,
    total: f64,
}

impl WspPairScoreCache {
    pub(crate) fn from_alignment<A: MafftAlphabet>(
        alignment: &MultipleAlignment<A>,
        weights: &FftNsiPairWeights,
        scorer: &ScientificScorer,
        gap_open: f64,
    ) -> Result<Self, MafftError> {
        validate_dimensions(alignment, weights)?;
        let matrix_len = alignment
            .sequence_count()
            .checked_mul(alignment.sequence_count())
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-i pair-score cache dimensions",
            })?;
        let mut weighted_scores = vec![0.0; matrix_len];
        let mut total = 0.0;
        for first in 0..alignment.sequence_count() {
            for second in 0..first {
                let contribution = weighted_pair_contribution(
                    alignment, weights, first, second, scorer, gap_open,
                )?;
                weighted_scores[first * alignment.sequence_count() + second] = contribution;
                total += contribution;
            }
        }
        validate_total(total, "FFT-NS-i cached weighted sum-of-pairs objective")?;
        Ok(Self {
            size: alignment.sequence_count(),
            weighted_scores: weighted_scores.into_boxed_slice(),
            total,
        })
    }

    pub(crate) const fn total(&self) -> f64 {
        self.total
    }

    pub(crate) fn score_candidate<A: MafftAlphabet>(
        &self,
        alignment: &MultipleAlignment<A>,
        partition: &FftNsiPartition,
        weights: &FftNsiPairWeights,
        scorer: &ScientificScorer,
        gap_open: f64,
    ) -> Result<f64, MafftError> {
        validate_dimensions(alignment, weights)?;
        if alignment.sequence_count() != self.size {
            return Err(MafftError::InternalInvariant {
                message: "FFT-NS-i candidate and pair-score cache dimensions disagree",
            });
        }
        let sides = partition_sides(partition, self.size)?;
        let mut total = 0.0;
        for first in 0..self.size {
            for second in 0..first {
                let contribution = if sides[first] == sides[second] {
                    self.weighted_scores[first * self.size + second]
                } else {
                    weighted_pair_contribution(alignment, weights, first, second, scorer, gap_open)?
                };
                total += contribution;
            }
        }
        validate_total(total, "FFT-NS-i cached candidate WSP objective")?;
        Ok(total)
    }

    pub(crate) fn resident_bytes(&self) -> Result<usize, MafftError> {
        self.weighted_scores
            .len()
            .checked_mul(size_of::<f64>())
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT-NS-i pair-score cache memory accounting",
            })
    }
}

fn partition_sides(partition: &FftNsiPartition, size: usize) -> Result<Vec<u8>, MafftError> {
    let mut sides = vec![u8::MAX; size];
    for (side, members) in [partition.left_members(), partition.right_members()]
        .into_iter()
        .enumerate()
    {
        for member in members.iter().copied() {
            let label = sides.get_mut(member).ok_or(MafftError::InternalInvariant {
                message: "FFT-NS-i WSP partition member is out of range",
            })?;
            if *label != u8::MAX {
                return Err(MafftError::InternalInvariant {
                    message: "FFT-NS-i WSP partition repeats a member",
                });
            }
            *label = u8::try_from(side).expect("two partition sides fit u8");
        }
    }
    if sides.contains(&u8::MAX) {
        return Err(MafftError::InternalInvariant {
            message: "FFT-NS-i WSP partition does not cover every row",
        });
    }
    Ok(sides)
}

pub(crate) fn weighted_sum_of_pairs<A: MafftAlphabet>(
    alignment: &MultipleAlignment<A>,
    weights: &FftNsiPairWeights,
    scorer: &ScientificScorer,
    gap_open: f64,
) -> Result<f64, MafftError> {
    validate_dimensions(alignment, weights)?;
    let mut total = 0.0;
    for first in 0..alignment.sequence_count() {
        for second in 0..first {
            total +=
                weighted_pair_contribution(alignment, weights, first, second, scorer, gap_open)?;
        }
    }
    validate_total(total, "FFT-NS-i weighted sum-of-pairs objective")?;
    Ok(total)
}

fn validate_dimensions<A: MafftAlphabet>(
    alignment: &MultipleAlignment<A>,
    weights: &FftNsiPairWeights,
) -> Result<(), MafftError> {
    if alignment.sequence_count() != weights.size() {
        return Err(MafftError::InternalInvariant {
            message: "WSP alignment and pair-weight dimensions disagree",
        });
    }
    Ok(())
}

fn weighted_pair_contribution<A: MafftAlphabet>(
    alignment: &MultipleAlignment<A>,
    weights: &FftNsiPairWeights,
    first: usize,
    second: usize,
    scorer: &ScientificScorer,
    gap_open: f64,
) -> Result<f64, MafftError> {
    let pair_weight = weights
        .get(first, second)
        .ok_or(MafftError::InternalInvariant {
            message: "WSP pair weight is missing",
        })?;
    let score = pair_alignment_score(
        alignment.rows()[first].symbols(),
        alignment.rows()[second].symbols(),
        scorer,
        gap_open,
    )?;
    let contribution = pair_weight * score;
    if !contribution.is_finite() {
        return Err(MafftError::NonFiniteComputation {
            stage: "FFT-NS-i weighted pair contribution",
        });
    }
    Ok(contribution)
}

fn validate_total(total: f64, stage: &'static str) -> Result<(), MafftError> {
    if !total.is_finite() {
        return Err(MafftError::NonFiniteComputation { stage });
    }
    Ok(())
}

pub(crate) fn pair_alignment_score<A: MafftAlphabet>(
    first: &[Option<Symbol<A>>],
    second: &[Option<Symbol<A>>],
    scorer: &ScientificScorer,
    gap_open: f64,
) -> Result<f64, MafftError> {
    if first.len() != second.len() {
        return Err(MafftError::InternalInvariant {
            message: "WSP induced pair is not rectangular",
        });
    }
    let mut score = 0.0;
    let mut gap_orientation = None;
    for (first_cell, second_cell) in first.iter().copied().zip(second.iter().copied()) {
        match (first_cell, second_cell) {
            (Some(first_symbol), Some(second_symbol)) => {
                score += scorer.score_indices(
                    first_symbol.alphabet_index(),
                    second_symbol.alphabet_index(),
                );
                gap_orientation = None;
            }
            (Some(_), None) => {
                if gap_orientation != Some(GapOrientation::FirstResidue) {
                    score -= gap_open;
                    gap_orientation = Some(GapOrientation::FirstResidue);
                }
            }
            (None, Some(_)) => {
                if gap_orientation != Some(GapOrientation::SecondResidue) {
                    score -= gap_open;
                    gap_orientation = Some(GapOrientation::SecondResidue);
                }
            }
            (None, None) => {}
        }
    }
    if !score.is_finite() {
        return Err(MafftError::NonFiniteComputation {
            stage: "FFT-NS-i induced pair score",
        });
    }
    Ok(score)
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, Symbol};
    use align_msa::{AlignedRow, GuideTree, GuideTreeEdge, GuideTreeNode, MultipleAlignment};

    use super::{WspPairScoreCache, pair_alignment_score, weighted_sum_of_pairs};
    use crate::matrix::ScientificScorer;
    use crate::refinement::GotohWeighting;
    use crate::{FftNs1Parameters, FftNsiParameters};

    fn dna_cells(text: &str) -> Vec<Option<Symbol<Dna>>> {
        text.bytes()
            .map(|byte| {
                if byte == b'-' {
                    None
                } else {
                    Some(Symbol::from_ascii(byte).unwrap())
                }
            })
            .collect()
    }

    fn alignment(rows: &[(&str, &str)]) -> MultipleAlignment<Dna> {
        MultipleAlignment::try_new(
            rows.iter()
                .map(|(identifier, text)| {
                    AlignedRow::new((*identifier).to_owned(), dna_cells(text)).unwrap()
                })
                .collect(),
        )
        .unwrap()
    }

    fn edge(length: f64, child: GuideTreeNode) -> GuideTreeEdge {
        GuideTreeEdge::new(length, child).unwrap()
    }

    fn four_leaf_tree() -> GuideTree {
        let left = GuideTreeNode::branch(
            edge(0.1, GuideTreeNode::leaf(0)),
            edge(0.2, GuideTreeNode::leaf(1)),
        );
        let right = GuideTreeNode::branch(
            edge(0.3, GuideTreeNode::leaf(2)),
            edge(0.4, GuideTreeNode::leaf(3)),
        );
        GuideTree::try_new(GuideTreeNode::branch(edge(0.5, left), edge(0.6, right)), 4).unwrap()
    }

    #[test]
    fn hand_scored_opposite_gap_orientations_open_twice() {
        let scorer = ScientificScorer::new::<Dna>(FftNs1Parameters::paper_2002()).unwrap();
        let first = dna_cells("A-CG");
        let second = dna_cells("AT-G");
        let expected = scorer.score_indices(0, 0) + scorer.score_indices(2, 2) - 4.8;
        let actual = pair_alignment_score(&first, &second, &scorer, 2.4).unwrap();
        assert!((actual - expected).abs() < 1.0e-12);
    }

    #[test]
    fn gap_gap_deletion_does_not_split_one_induced_gap_run() {
        let scorer = ScientificScorer::new::<Dna>(FftNs1Parameters::paper_2002()).unwrap();
        let first = dna_cells("A-C");
        let second = dna_cells("---");
        let actual = pair_alignment_score(&first, &second, &scorer, 2.4).unwrap();
        assert!((actual + 2.4).abs() < 1.0e-12);
    }

    #[test]
    fn leading_internal_and_terminal_runs_are_each_penalized() {
        let scorer = ScientificScorer::new::<Dna>(FftNs1Parameters::paper_2002()).unwrap();
        let first = dna_cells("-AC-GT");
        let second = dna_cells("T-CAG-");
        let actual = pair_alignment_score(&first, &second, &scorer, 2.4).unwrap();
        let expected = scorer.score_indices(1, 1) + scorer.score_indices(2, 2) - 4.0 * 2.4;
        assert!((actual - expected).abs() < 1.0e-12);
    }

    #[test]
    fn fixed_subgroup_contributions_reproduce_direct_wsp_bit_for_bit() {
        let current = alignment(&[("a", "A-CG"), ("b", "ATCG"), ("c", "ACG-"), ("d", "A-GT")]);
        let candidate = alignment(&[
            ("a", "A--CG"),
            ("b", "AT-CG"),
            ("c", "A-CG-"),
            ("d", "A--GT"),
        ]);
        let parameters = FftNsiParameters::paper_2002();
        let weighting = GotohWeighting::from_guide_tree(&four_leaf_tree(), parameters).unwrap();
        let partition = weighting
            .partitions()
            .iter()
            .find(|partition| partition.public.left_members() == [0, 1])
            .unwrap();
        let scorer = ScientificScorer::new::<Dna>(parameters.initial().progressive()).unwrap();
        let gap_open = parameters.initial().progressive().gap_open();
        let cache = WspPairScoreCache::from_alignment(
            &current,
            weighting.pair_weights(),
            &scorer,
            gap_open,
        )
        .unwrap();
        let current_direct =
            weighted_sum_of_pairs(&current, weighting.pair_weights(), &scorer, gap_open).unwrap();
        assert_eq!(cache.total().to_bits(), current_direct.to_bits());

        let optimized = cache
            .score_candidate(
                &candidate,
                &partition.public,
                weighting.pair_weights(),
                &scorer,
                gap_open,
            )
            .unwrap();
        let direct =
            weighted_sum_of_pairs(&candidate, weighting.pair_weights(), &scorer, gap_open).unwrap();
        assert_eq!(optimized.to_bits(), direct.to_bits());
    }
}
