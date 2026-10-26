#![allow(dead_code)]

use align_core::{
    AffineGap, AlignmentOperation, Alphabet, PairwiseAlignment, Score, Sequence, SubstitutionScorer,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum OracleTermination {
    MatrixExhausted,
    FrontiersExhausted,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct OracleSummary {
    pub(crate) score: Score,
    pub(crate) first_end: usize,
    pub(crate) second_end: usize,
    pub(crate) evaluated_cells: usize,
    pub(crate) retained_cells: usize,
    pub(crate) antidiagonals: usize,
    pub(crate) maximum_frontier_width: usize,
    pub(crate) termination: OracleTermination,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct OracleEndpoint {
    pub(crate) score: Score,
    pub(crate) first_end: usize,
    pub(crate) second_end: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct Cell {
    best: Option<Score>,
    deletion: Option<Score>,
    insertion: Option<Score>,
}

#[derive(Clone, Copy, Debug)]
struct Endpoint {
    score: Score,
    first: usize,
    second: usize,
}

impl Endpoint {
    fn consider(&mut self, score: Score, first: usize, second: usize) {
        let distance = first + second;
        let selected_distance = self.first + self.second;
        if score > self.score
            || (score == self.score
                && (distance > selected_distance
                    || (distance == selected_distance
                        && (second > self.second
                            || (second == self.second && first > self.first)))))
        {
            *self = Self {
                score,
                first,
                second,
            };
        }
    }
}

fn add(base: Option<Score>, increment: Score) -> Option<Score> {
    base.and_then(|score| score.checked_add(increment))
}

fn maximum(first: Option<Score>, second: Option<Score>) -> Option<Score> {
    match (first, second) {
        (Some(first), Some(second)) => Some(first.max(second)),
        (Some(score), None) | (None, Some(score)) => Some(score),
        (None, None) => None,
    }
}

fn matrix_index(first_prefix: usize, second_prefix: usize, columns: usize) -> usize {
    first_prefix * columns + second_prefix
}

/// Dense, full-history oracle written independently from the rolling-frontier implementation.
#[allow(clippy::too_many_lines)]
pub(crate) fn oracle<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    substitution: &S,
    gap: AffineGap,
    threshold: Score,
) -> OracleSummary
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    assert!(threshold >= 0, "the test oracle requires a valid threshold");
    let columns = second.len() + 1;
    let mut matrix = vec![Cell::default(); (first.len() + 1) * columns];
    matrix[0].best = Some(0);

    let maximum_steps = first.len() + second.len();
    let mut retained_by_step = vec![Vec::<(usize, usize)>::new(); maximum_steps + 1];
    let mut active_by_step = vec![None; maximum_steps + 1];
    retained_by_step[0].push((0, 0));
    active_by_step[0] = Some((0, 0));
    let mut endpoint = Endpoint {
        score: 0,
        first: 0,
        second: 0,
    };
    let mut global_best = 0;
    let mut evaluated_cells = 1;
    let mut retained_cells = 1;
    let mut antidiagonals = 1;
    let mut maximum_frontier_width = 1;
    let mut termination = OracleTermination::MatrixExhausted;

    for step in 1..=maximum_steps {
        let mut candidate_first_prefixes = Vec::new();
        if let Some((lower, upper)) = active_by_step[step - 1] {
            for first_prefix in lower..=upper {
                let second_prefix = step - 1 - first_prefix;
                if first_prefix < first.len() {
                    candidate_first_prefixes.push(first_prefix + 1);
                }
                if second_prefix < second.len() {
                    candidate_first_prefixes.push(first_prefix);
                }
            }
        }
        if step >= 2 {
            if let Some((lower, upper)) = active_by_step[step - 2] {
                for first_prefix in lower..=upper {
                    let second_prefix = step - 2 - first_prefix;
                    if first_prefix < first.len() && second_prefix < second.len() {
                        candidate_first_prefixes.push(first_prefix + 1);
                    }
                }
            }
        }
        let Some(lower) = candidate_first_prefixes.iter().copied().min() else {
            termination = OracleTermination::FrontiersExhausted;
            break;
        };
        let upper = candidate_first_prefixes
            .iter()
            .copied()
            .max()
            .expect("a minimum implies a maximum");
        let width = upper - lower + 1;
        evaluated_cells += width;
        antidiagonals += 1;
        maximum_frontier_width = maximum_frontier_width.max(width);
        let prior_global_best = global_best;
        let mut diagonal_best = global_best;

        for first_prefix in lower..=upper {
            let second_prefix = step - first_prefix;
            let index = matrix_index(first_prefix, second_prefix, columns);
            let deletion = if first_prefix > 0 {
                let upper = matrix[matrix_index(first_prefix - 1, second_prefix, columns)];
                maximum(
                    add(upper.best, gap.open_score()),
                    add(upper.deletion, gap.extension_score()),
                )
            } else {
                None
            };
            let insertion = if second_prefix > 0 {
                let left = matrix[matrix_index(first_prefix, second_prefix - 1, columns)];
                maximum(
                    add(left.best, gap.open_score()),
                    add(left.insertion, gap.extension_score()),
                )
            } else {
                None
            };
            let diagonal = if first_prefix > 0 && second_prefix > 0 {
                let predecessor =
                    matrix[matrix_index(first_prefix - 1, second_prefix - 1, columns)];
                add(
                    predecessor.best,
                    substitution.score(
                        first
                            .get(first_prefix - 1)
                            .expect("oracle prefix is in range"),
                        second
                            .get(second_prefix - 1)
                            .expect("oracle prefix is in range"),
                    ),
                )
            } else {
                None
            };
            let best = maximum(maximum(diagonal, deletion), insertion);
            if let Some(score) = best {
                if i128::from(prior_global_best) - i128::from(score) <= i128::from(threshold) {
                    matrix[index] = Cell {
                        best: Some(score),
                        deletion,
                        insertion,
                    };
                    retained_by_step[step].push((first_prefix, second_prefix));
                    retained_cells += 1;
                    diagonal_best = diagonal_best.max(score);
                    endpoint.consider(score, first_prefix, second_prefix);
                }
            }
        }
        global_best = diagonal_best;
        active_by_step[step] = retained_by_step[step]
            .first()
            .zip(retained_by_step[step].last())
            .map(|(&(first, _), &(last, _))| (first, last));
    }

    OracleSummary {
        score: endpoint.score,
        first_end: endpoint.first,
        second_end: endpoint.second,
        evaluated_cells,
        retained_cells,
        antidiagonals,
        maximum_frontier_width,
        termination,
    }
}

/// Full affine prefix matrix with no X-drop pruning or frontier logic.
pub(crate) fn unrestricted_oracle<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    substitution: &S,
    gap: AffineGap,
) -> OracleEndpoint
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let columns = second.len() + 1;
    let mut matrix = vec![Cell::default(); (first.len() + 1) * columns];
    matrix[0].best = Some(0);
    let mut endpoint = Endpoint {
        score: 0,
        first: 0,
        second: 0,
    };
    for first_prefix in 0..=first.len() {
        for second_prefix in 0..=second.len() {
            if first_prefix == 0 && second_prefix == 0 {
                continue;
            }
            let index = matrix_index(first_prefix, second_prefix, columns);
            let deletion = if first_prefix > 0 {
                let upper = matrix[matrix_index(first_prefix - 1, second_prefix, columns)];
                maximum(
                    add(upper.best, gap.open_score()),
                    add(upper.deletion, gap.extension_score()),
                )
            } else {
                None
            };
            let insertion = if second_prefix > 0 {
                let left = matrix[matrix_index(first_prefix, second_prefix - 1, columns)];
                maximum(
                    add(left.best, gap.open_score()),
                    add(left.insertion, gap.extension_score()),
                )
            } else {
                None
            };
            let diagonal = if first_prefix > 0 && second_prefix > 0 {
                let predecessor =
                    matrix[matrix_index(first_prefix - 1, second_prefix - 1, columns)];
                add(
                    predecessor.best,
                    substitution.score(
                        first
                            .get(first_prefix - 1)
                            .expect("oracle prefix is in range"),
                        second
                            .get(second_prefix - 1)
                            .expect("oracle prefix is in range"),
                    ),
                )
            } else {
                None
            };
            let best = maximum(maximum(diagonal, deletion), insertion);
            matrix[index] = Cell {
                best,
                deletion,
                insertion,
            };
            if let Some(score) = best {
                endpoint.consider(score, first_prefix, second_prefix);
            }
        }
    }
    OracleEndpoint {
        score: endpoint.score,
        first_end: endpoint.first,
        second_end: endpoint.second,
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum PriorGap {
    None,
    Deletion,
    Insertion,
}

/// Enumerates every operation path and every possible open endpoint for tiny inputs.
pub(crate) fn exhaustive_path_oracle<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    substitution: &S,
    gap: AffineGap,
) -> OracleEndpoint
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    #[allow(clippy::too_many_arguments)]
    fn visit<A, S>(
        first: &Sequence<A>,
        second: &Sequence<A>,
        substitution: &S,
        gap: AffineGap,
        first_prefix: usize,
        second_prefix: usize,
        score: Score,
        prior_gap: PriorGap,
        endpoint: &mut Endpoint,
    ) where
        A: Alphabet,
        S: SubstitutionScorer<A>,
    {
        endpoint.consider(score, first_prefix, second_prefix);
        if first_prefix < first.len() && second_prefix < second.len() {
            visit(
                first,
                second,
                substitution,
                gap,
                first_prefix + 1,
                second_prefix + 1,
                score
                    + substitution.score(
                        first.get(first_prefix).expect("path prefix is in range"),
                        second.get(second_prefix).expect("path prefix is in range"),
                    ),
                PriorGap::None,
                endpoint,
            );
        }
        if first_prefix < first.len() {
            visit(
                first,
                second,
                substitution,
                gap,
                first_prefix + 1,
                second_prefix,
                score
                    + if matches!(prior_gap, PriorGap::Deletion) {
                        gap.extension_score()
                    } else {
                        gap.open_score()
                    },
                PriorGap::Deletion,
                endpoint,
            );
        }
        if second_prefix < second.len() {
            visit(
                first,
                second,
                substitution,
                gap,
                first_prefix,
                second_prefix + 1,
                score
                    + if matches!(prior_gap, PriorGap::Insertion) {
                        gap.extension_score()
                    } else {
                        gap.open_score()
                    },
                PriorGap::Insertion,
                endpoint,
            );
        }
    }

    let mut endpoint = Endpoint {
        score: 0,
        first: 0,
        second: 0,
    };
    visit(
        first,
        second,
        substitution,
        gap,
        0,
        0,
        0,
        PriorGap::None,
        &mut endpoint,
    );
    OracleEndpoint {
        score: endpoint.score,
        first_end: endpoint.first,
        second_end: endpoint.second,
    }
}

/// Independently replays a returned operation stream and checks its structural labels.
pub(crate) fn replay<A, S>(
    first: &Sequence<A>,
    second: &Sequence<A>,
    alignment: &PairwiseAlignment,
    substitution: &S,
    gap: AffineGap,
) -> Score
where
    A: Alphabet,
    S: SubstitutionScorer<A>,
{
    let first_range = alignment.first_range();
    let second_range = alignment.second_range();
    let mut first_index = first_range.start;
    let mut second_index = second_range.start;
    let mut score = 0_i64;
    let mut preceding_gap = None;

    for &operation in alignment.operations() {
        let increment = match operation {
            AlignmentOperation::Match | AlignmentOperation::Substitution => {
                let first_symbol = first.get(first_index).expect("first replay index is valid");
                let second_symbol = second
                    .get(second_index)
                    .expect("second replay index is valid");
                assert_eq!(
                    matches!(operation, AlignmentOperation::Match),
                    first_symbol == second_symbol,
                    "diagonal operation label must reflect literal equality"
                );
                first_index += 1;
                second_index += 1;
                preceding_gap = None;
                substitution.score(first_symbol, second_symbol)
            }
            AlignmentOperation::Deletion => {
                first_index += 1;
                let increment = if matches!(preceding_gap, Some(AlignmentOperation::Deletion)) {
                    gap.extension_score()
                } else {
                    gap.open_score()
                };
                preceding_gap = Some(operation);
                increment
            }
            AlignmentOperation::Insertion => {
                second_index += 1;
                let increment = if matches!(preceding_gap, Some(AlignmentOperation::Insertion)) {
                    gap.extension_score()
                } else {
                    gap.open_score()
                };
                preceding_gap = Some(operation);
                increment
            }
        };
        score = score
            .checked_add(increment)
            .expect("generated replay score fits");
    }

    assert_eq!(first_index, first_range.end);
    assert_eq!(second_index, second_range.end);
    assert_eq!(score, alignment.score());
    score
}
