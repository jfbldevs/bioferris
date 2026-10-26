use crate::FftNs1Parameters;
use crate::fft::LagCorrelation;
use crate::profile::{PreparedProfile, ProfilePairScores};

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct HomologousSegment {
    pub(crate) left_start: usize,
    pub(crate) right_start: usize,
    pub(crate) length: usize,
    pub(crate) score: f64,
}

impl HomologousSegment {
    fn left_end(self) -> usize {
        self.left_start + self.length
    }

    fn right_end(self) -> usize {
        self.right_start + self.length
    }

    pub(crate) fn center(self) -> (usize, usize) {
        (
            self.left_start + self.length / 2,
            self.right_start + self.length / 2,
        )
    }

    fn key(self) -> (usize, usize, usize) {
        (self.left_start, self.right_start, self.length)
    }
}

pub(crate) fn find_segments(
    left: &PreparedProfile,
    right: &PreparedProfile,
    scores: &ProfilePairScores,
    lags: &[LagCorrelation],
    parameters: FftNs1Parameters,
) -> Vec<HomologousSegment> {
    let window = parameters.window_length();
    #[allow(clippy::cast_precision_loss)]
    let threshold = parameters.segment_threshold_per_column() * window as f64;
    let mut candidates = Vec::new();
    for correlation in lags {
        let (left_start, right_start) = if correlation.lag >= 0 {
            (0usize, correlation.lag.unsigned_abs())
        } else {
            (correlation.lag.unsigned_abs(), 0usize)
        };
        if left_start >= left.width() || right_start >= right.width() {
            continue;
        }
        let overlap = (left.width() - left_start).min(right.width() - right_start);
        if overlap < window {
            continue;
        }
        let diagonal_scores = (0..overlap)
            .map(|offset| scores.homology_without_bonus(left_start + offset, right_start + offset))
            .collect::<Vec<_>>();
        let mut window_sum: f64 = diagonal_scores[..window].iter().sum();
        let mut run_start = None::<usize>;
        let window_count = overlap - window + 1;
        for offset in 0..window_count {
            if offset > 0 {
                window_sum += diagonal_scores[offset + window - 1] - diagonal_scores[offset - 1];
            }
            if window_sum > threshold {
                run_start.get_or_insert(offset);
            } else if let Some(start) = run_start.take() {
                add_merged_run(
                    &mut candidates,
                    left_start,
                    right_start,
                    start,
                    offset - 1,
                    &diagonal_scores,
                    parameters.maximum_segment_length(),
                    window,
                );
            }
        }
        if let Some(start) = run_start {
            add_merged_run(
                &mut candidates,
                left_start,
                right_start,
                start,
                window_count - 1,
                &diagonal_scores,
                parameters.maximum_segment_length(),
                window,
            );
        }
    }

    candidates.sort_by(|first, second| {
        first
            .key()
            .cmp(&second.key())
            .then_with(|| second.score.total_cmp(&first.score))
    });
    candidates.dedup_by(|second, first| {
        if first.key() == second.key() {
            if second.score.total_cmp(&first.score).is_gt() {
                first.score = second.score;
            }
            true
        } else {
            false
        }
    });
    candidates
}

#[allow(clippy::too_many_arguments)]
fn add_merged_run(
    output: &mut Vec<HomologousSegment>,
    base_left: usize,
    base_right: usize,
    first_window: usize,
    last_window: usize,
    diagonal_scores: &[f64],
    maximum_length: usize,
    window: usize,
) {
    let run_length = last_window - first_window + window;
    let mut consumed = 0usize;
    while consumed < run_length {
        let length = (run_length - consumed).min(maximum_length);
        let start = first_window + consumed;
        let score = diagonal_scores[start..start + length].iter().sum();
        output.push(HomologousSegment {
            left_start: base_left + start,
            right_start: base_right + start,
            length,
            score,
        });
        consumed += length;
    }
}

pub(crate) fn select_monotone_chain(candidates: &[HomologousSegment]) -> Vec<HomologousSegment> {
    if candidates.is_empty() {
        return Vec::new();
    }
    let mut ordered = candidates.to_vec();
    ordered.sort_by_key(|segment| segment.key());
    let count = ordered.len();
    let mut totals = vec![0.0; count];
    let mut coverage = vec![0usize; count];
    let mut predecessors = vec![None::<usize>; count];

    for current in 0..count {
        totals[current] = ordered[current].score;
        coverage[current] = ordered[current].length;
        for previous in 0..current {
            if ordered[previous].left_end() > ordered[current].left_start
                || ordered[previous].right_end() > ordered[current].right_start
            {
                continue;
            }
            let candidate_total = totals[previous] + ordered[current].score;
            let candidate_coverage = coverage[previous] + ordered[current].length;
            let replace = candidate_total.total_cmp(&totals[current]).is_gt()
                || (candidate_total.total_cmp(&totals[current]).is_eq()
                    && (candidate_coverage > coverage[current]
                        || (candidate_coverage == coverage[current]
                            && chain_with_tail_is_earlier(
                                &ordered,
                                &predecessors,
                                Some(previous),
                                predecessors[current],
                                current,
                            ))));
            if replace {
                totals[current] = candidate_total;
                coverage[current] = candidate_coverage;
                predecessors[current] = Some(previous);
            }
        }
    }

    let mut best = 0usize;
    for candidate in 1..count {
        let replace = totals[candidate].total_cmp(&totals[best]).is_gt()
            || (totals[candidate].total_cmp(&totals[best]).is_eq()
                && (coverage[candidate] > coverage[best]
                    || (coverage[candidate] == coverage[best]
                        && chain_is_earlier(&ordered, &predecessors, candidate, best))));
        if replace {
            best = candidate;
        }
    }
    let indices = chain_indices(&predecessors, Some(best));
    indices.into_iter().map(|index| ordered[index]).collect()
}

fn chain_with_tail_is_earlier(
    segments: &[HomologousSegment],
    predecessors: &[Option<usize>],
    first_predecessor: Option<usize>,
    second_predecessor: Option<usize>,
    tail: usize,
) -> bool {
    let mut first = chain_indices(predecessors, first_predecessor);
    let mut second = chain_indices(predecessors, second_predecessor);
    first.push(tail);
    second.push(tail);
    compare_chain_keys(segments, &first, &second).is_lt()
}

fn chain_is_earlier(
    segments: &[HomologousSegment],
    predecessors: &[Option<usize>],
    first: usize,
    second: usize,
) -> bool {
    let first = chain_indices(predecessors, Some(first));
    let second = chain_indices(predecessors, Some(second));
    compare_chain_keys(segments, &first, &second).is_lt()
}

fn compare_chain_keys(
    segments: &[HomologousSegment],
    first: &[usize],
    second: &[usize],
) -> core::cmp::Ordering {
    first
        .iter()
        .map(|index| segments[*index].key())
        .cmp(second.iter().map(|index| segments[*index].key()))
}

fn chain_indices(predecessors: &[Option<usize>], mut current: Option<usize>) -> Vec<usize> {
    let mut indices = Vec::new();
    while let Some(index) = current {
        indices.push(index);
        current = predecessors[index];
    }
    indices.reverse();
    indices
}

#[cfg(test)]
mod tests {
    use super::{HomologousSegment, select_monotone_chain};

    fn segment(left: usize, right: usize, length: usize, score: f64) -> HomologousSegment {
        HomologousSegment {
            left_start: left,
            right_start: right,
            length,
            score,
        }
    }

    #[test]
    fn segment_dp_rejects_crossing_high_scores() {
        let candidates = [
            segment(0, 10, 3, 8.0),
            segment(10, 0, 3, 9.0),
            segment(0, 0, 3, 6.0),
            segment(5, 5, 3, 6.0),
        ];
        let selected = select_monotone_chain(&candidates);
        assert_eq!(selected, vec![segment(0, 0, 3, 6.0), segment(5, 5, 3, 6.0)]);
    }

    #[test]
    fn equal_chain_scores_prefer_more_coverage_then_earlier_coordinates() {
        let coverage = [segment(0, 0, 3, 10.0), segment(1, 1, 4, 10.0)];
        assert_eq!(
            select_monotone_chain(&coverage),
            vec![segment(1, 1, 4, 10.0)]
        );

        let coordinates = [
            segment(0, 0, 2, 5.0),
            segment(2, 2, 2, 5.0),
            segment(1, 1, 4, 10.0),
        ];
        assert_eq!(
            select_monotone_chain(&coordinates),
            vec![segment(0, 0, 2, 5.0), segment(2, 2, 2, 5.0)]
        );
    }

    #[test]
    fn center_uses_the_right_middle_for_even_segments() {
        assert_eq!(segment(4, 9, 6, 1.0).center(), (7, 12));
    }
}
