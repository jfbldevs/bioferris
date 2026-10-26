use core::mem::size_of;
use std::sync::OnceLock;

use align_core::{Alphabet, Symbol};
use align_msa::SequenceRecord;

use crate::omega_matrix::{
    SubstitutionModel, add_projected_symbol, background, canonical_size, normalize_positive,
    substitution_prior_into,
};
use crate::{ClustalOmegaError, ClustalOmegaParameters, OmegaEngineMode};

const STATE_COUNT: usize = 5;
const MM: usize = 0;
const MI: usize = 1;
const IM: usize = 2;
const DG: usize = 3;
const GD: usize = 4;
const TRACE_START: u8 = u8::MAX;
const NEGATIVE_INFINITY: f64 = f64::NEG_INFINITY;

/// Engine actually used for one progressive HMM-HMM merge.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OmegaEngine {
    /// Forward/Backward posterior probabilities followed by MAC traceback.
    Mac,
    /// Five-state maximum-probability traceback.
    Viterbi,
}

#[derive(Clone, Debug)]
pub(crate) struct WorkingRow<A: Alphabet> {
    pub(crate) sequence_index: usize,
    pub(crate) cells: Vec<Option<Symbol<A>>>,
}

#[derive(Clone, Debug)]
pub(crate) struct WorkingProfile<A: Alphabet> {
    pub(crate) rows: Vec<WorkingRow<A>>,
}

impl<A: Alphabet> WorkingProfile<A> {
    pub(crate) fn single(record: &SequenceRecord<A>, sequence_index: usize) -> Self {
        Self {
            rows: vec![WorkingRow {
                sequence_index,
                cells: record.sequence().iter().map(Some).collect(),
            }],
        }
    }

    pub(crate) fn width(&self) -> usize {
        self.rows.first().map_or(0, |row| row.cells.len())
    }

    pub(crate) fn minimum_member(&self) -> usize {
        self.rows
            .iter()
            .map(|row| row.sequence_index)
            .min()
            .expect("working profiles always contain a row")
    }

    pub(crate) fn members(&self) -> Vec<usize> {
        let mut members = self
            .rows
            .iter()
            .map(|row| row.sequence_index)
            .collect::<Vec<_>>();
        members.sort_unstable();
        members
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Transition {
    mm: f64,
    mi: f64,
    md: f64,
    im: f64,
    ii: f64,
    dm: f64,
    dd: f64,
}

impl Transition {
    fn into_log2(self) -> LogTransition {
        LogTransition {
            mm: self.mm.log2(),
            mi: self.mi.log2(),
            md: self.md.log2(),
            im: self.im.log2(),
            ii: self.ii.log2(),
            dm: self.dm.log2(),
            dd: self.dd.log2(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct LogTransition {
    mm: f64,
    mi: f64,
    md: f64,
    im: f64,
    ii: f64,
    dm: f64,
    dd: f64,
}

#[derive(Clone, Debug)]
pub(crate) struct ProfileHMM {
    emissions: Vec<f64>,
    alphabet_size: usize,
    probability_transitions: Vec<Transition>,
    transitions: Vec<LogTransition>,
    pub(crate) effective_diversity: Vec<f64>,
    pub(crate) sequence_weights: Vec<f64>,
}

impl ProfileHMM {
    fn len(&self) -> usize {
        self.emissions.len() / self.alphabet_size
    }

    fn emission(&self, position: usize) -> &[f64] {
        &self.emissions[position * self.alphabet_size..(position + 1) * self.alphabet_size]
    }

    fn transition(&self, position: usize) -> LogTransition {
        self.transitions[position.min(self.transitions.len() - 1)]
    }

    fn probability_transition(&self, position: usize) -> Transition {
        self.probability_transitions[position.min(self.probability_transitions.len() - 1)]
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ColumnOperation {
    Both,
    Left,
    Right,
}

#[derive(Clone, Debug)]
pub(crate) struct HmmMergeObservation {
    pub(crate) engine: OmegaEngine,
    pub(crate) score: f64,
    pub(crate) left_core: core::ops::Range<usize>,
    pub(crate) right_core: core::ops::Range<usize>,
    pub(crate) operations: Vec<ColumnOperation>,
    pub(crate) workspace_bytes: usize,
    pub(crate) posterior_bounds: Option<(f64, f64)>,
}

#[derive(Clone, Debug)]
pub(crate) struct HmmMergeResult<A: Alphabet> {
    pub(crate) profile: WorkingProfile<A>,
    pub(crate) observation: HmmMergeObservation,
}

#[derive(Clone, Debug)]
struct CoreAlignment {
    score: f64,
    left_start: usize,
    left_end: usize,
    right_start: usize,
    right_end: usize,
    operations: Vec<ColumnOperation>,
    workspace_bytes: usize,
    posterior_bounds: Option<(f64, f64)>,
}

#[cfg(test)]
fn build_profile_hmm<A: Alphabet>(
    profile: &WorkingProfile<A>,
    parameters: &ClustalOmegaParameters,
    is_protein: bool,
) -> Result<ProfileHMM, ClustalOmegaError> {
    build_profile_hmm_with_model(
        profile,
        parameters,
        is_protein,
        shared_substitution_model(is_protein)?,
    )
}

fn shared_substitution_model(
    is_protein: bool,
) -> Result<&'static SubstitutionModel, ClustalOmegaError> {
    static PROTEIN: OnceLock<Result<SubstitutionModel, ClustalOmegaError>> = OnceLock::new();
    static NUCLEOTIDE: OnceLock<Result<SubstitutionModel, ClustalOmegaError>> = OnceLock::new();
    let model = if is_protein {
        PROTEIN.get_or_init(|| SubstitutionModel::new(true))
    } else {
        NUCLEOTIDE.get_or_init(|| SubstitutionModel::new(false))
    };
    model.as_ref().map_err(Clone::clone)
}

fn build_profile_hmm_with_model<A: Alphabet>(
    profile: &WorkingProfile<A>,
    parameters: &ClustalOmegaParameters,
    is_protein: bool,
    substitution_model: &SubstitutionModel,
) -> Result<ProfileHMM, ClustalOmegaError> {
    debug_assert!(validate_profile(profile).is_ok());
    let weights = henikoff_weights(profile)?;
    let alphabet_size = canonical_size(is_protein);
    let emission_values = profile.width().checked_mul(alphabet_size).ok_or(
        ClustalOmegaError::ArithmeticOverflow {
            stage: "profile HMM emission dimensions",
        },
    )?;
    let mut emissions = Vec::with_capacity(emission_values);
    let mut effective_diversity = Vec::with_capacity(profile.width());
    let mut observed = vec![0.0; alphabet_size];
    let mut substitution = vec![0.0; alphabet_size];
    let mut emission = vec![0.0; alphabet_size];
    let (pseudocount_a, pseudocount_b, pseudocount_c) = parameters.emission_pseudocounts();
    for column in 0..profile.width() {
        observed.fill(0.0);
        for (row, weight) in profile.rows.iter().zip(&weights) {
            if let Some(symbol) = row.cells[column] {
                add_projected_symbol(&mut observed, symbol.as_ascii(), *weight, is_protein)?;
            }
        }
        normalize_positive(&mut observed, "profile observed emissions")?;
        let entropy = observed.iter().fold(0.0, |sum, probability| {
            if *probability <= 0.0 {
                sum
            } else {
                sum - probability * probability.log2()
            }
        });
        let diversity = 2.0_f64.powf(entropy);
        if !diversity.is_finite() || diversity < 1.0 {
            return Err(ClustalOmegaError::InvalidProbability {
                stage: "profile effective diversity",
                value: diversity,
            });
        }
        substitution_prior_into(&observed, substitution_model, &mut substitution)?;
        let scaled_diversity = ((diversity - 1.0) / pseudocount_b).powf(pseudocount_c);
        let tau = pseudocount_a / (1.0 + scaled_diversity);
        if !tau.is_finite() || !(0.0..=1.0).contains(&tau) {
            return Err(ClustalOmegaError::InvalidProbability {
                stage: "emission pseudocount admixture",
                value: tau,
            });
        }
        for ((destination, observed), prior) in
            emission.iter_mut().zip(&observed).zip(&substitution)
        {
            *destination = (1.0 - tau) * observed + tau * prior;
        }
        normalize_positive(&mut emission, "profile pseudocount emissions")?;
        emissions.extend_from_slice(&emission);
        effective_diversity.push(diversity);
    }
    let probability_transitions = profile_transitions(profile, &weights, parameters)?;
    let transitions = probability_transitions
        .iter()
        .copied()
        .map(Transition::into_log2)
        .collect();
    Ok(ProfileHMM {
        emissions,
        alphabet_size,
        probability_transitions,
        transitions,
        effective_diversity,
        sequence_weights: weights,
    })
}

pub(crate) fn align_profiles<A: Alphabet>(
    left: WorkingProfile<A>,
    right: WorkingProfile<A>,
    parameters: &ClustalOmegaParameters,
    is_protein: bool,
) -> Result<HmmMergeResult<A>, ClustalOmegaError> {
    let substitution_model = shared_substitution_model(is_protein)?;
    let left_hmm = build_profile_hmm_with_model(&left, parameters, is_protein, substitution_model)?;
    let right_hmm =
        build_profile_hmm_with_model(&right, parameters, is_protein, substitution_model)?;
    debug_assert_eq!(left_hmm.effective_diversity.len(), left_hmm.len());
    debug_assert_eq!(right_hmm.effective_diversity.len(), right_hmm.len());
    debug_assert_eq!(left_hmm.sequence_weights.len(), left.rows.len());
    debug_assert_eq!(right_hmm.sequence_weights.len(), right.rows.len());
    let mac_workspace = mac_workspace_bytes(left_hmm.len(), right_hmm.len())?;
    let engine = match parameters.engine_mode() {
        OmegaEngineMode::Mac => OmegaEngine::Mac,
        OmegaEngineMode::Viterbi => OmegaEngine::Viterbi,
        OmegaEngineMode::Auto => {
            if mac_workspace <= parameters.mac_memory_limit_bytes() {
                OmegaEngine::Mac
            } else {
                OmegaEngine::Viterbi
            }
        }
    };
    let core = if left_hmm.len() == 0 || right_hmm.len() == 0 {
        CoreAlignment {
            score: 0.0,
            left_start: 0,
            left_end: 0,
            right_start: 0,
            right_end: 0,
            operations: Vec::new(),
            workspace_bytes: 0,
            posterior_bounds: (engine == OmegaEngine::Mac).then_some((0.0, 0.0)),
        }
    } else {
        match engine {
            OmegaEngine::Mac => mac_align_scaled(
                &left_hmm,
                &right_hmm,
                parameters,
                left.minimum_member() <= right.minimum_member(),
            )?,
            OmegaEngine::Viterbi => viterbi_align(&left_hmm, &right_hmm, parameters)?,
        }
    };
    let operations = attach_terminals(
        left.width(),
        right.width(),
        core.left_start,
        core.left_end,
        core.right_start,
        core.right_end,
        &core.operations,
    )?;
    let profile = merge_profiles(left, right, &operations)?;
    let observation = HmmMergeObservation {
        engine,
        score: core.score,
        left_core: core.left_start..core.left_end,
        right_core: core.right_start..core.right_end,
        operations,
        workspace_bytes: core.workspace_bytes,
        posterior_bounds: core.posterior_bounds,
    };
    Ok(HmmMergeResult {
        profile,
        observation,
    })
}

#[allow(clippy::cast_precision_loss)]
fn henikoff_weights<A: Alphabet>(
    profile: &WorkingProfile<A>,
) -> Result<Vec<f64>, ClustalOmegaError> {
    let mut weights = vec![0.0; profile.rows.len()];
    let mut counts = vec![0_usize; A::SYMBOLS.len()];
    for column in 0..profile.width() {
        counts.fill(0);
        for row in &profile.rows {
            if let Some(symbol) = row.cells[column] {
                counts[symbol.alphabet_index()] += 1;
            }
        }
        let distinct = counts.iter().filter(|count| **count != 0).count();
        if distinct == 0 {
            return Err(ClustalOmegaError::InternalInvariant {
                message: "working profile contains an all-gap column",
            });
        }
        for (row_index, row) in profile.rows.iter().enumerate() {
            if let Some(symbol) = row.cells[column] {
                let count = counts[symbol.alphabet_index()];
                weights[row_index] += 1.0 / (distinct * count) as f64;
            }
        }
    }
    let smallest_positive = weights
        .iter()
        .copied()
        .filter(|value| *value > 0.0)
        .min_by(f64::total_cmp);
    match smallest_positive {
        Some(smallest) => {
            for weight in &mut weights {
                if *weight == 0.0 {
                    *weight = smallest;
                }
            }
        }
        None => weights.fill(1.0),
    }
    let total = weights.iter().sum::<f64>();
    if !total.is_finite() || total <= 0.0 {
        return Err(ClustalOmegaError::InvalidProbability {
            stage: "Henikoff sequence-weight normalization",
            value: total,
        });
    }
    let scale = profile.rows.len() as f64 / total;
    for weight in &mut weights {
        *weight *= scale;
    }
    Ok(weights)
}

#[allow(clippy::similar_names)]
fn profile_transitions<A: Alphabet>(
    profile: &WorkingProfile<A>,
    weights: &[f64],
    parameters: &ClustalOmegaParameters,
) -> Result<Vec<Transition>, ClustalOmegaError> {
    let (gapb, gapd, gape) = parameters.transition_pseudocounts();
    let open = gapd * 0.0286;
    let prior_mm = 1.0 - 2.0 * open;
    let prior_extend = gape / (gape - 1.0 + 1.0 / 0.75);
    if prior_mm <= 0.0 || !(0.0..1.0).contains(&prior_extend) && gape != 0.0 {
        return Err(ClustalOmegaError::InvalidProbability {
            stage: "transition pseudocount priors",
            value: prior_mm.min(prior_extend),
        });
    }
    let admixture = gapb.max(1.0e-15);
    let core_bounds = profile
        .rows
        .iter()
        .map(|row| {
            row.cells
                .iter()
                .position(Option::is_some)
                .zip(row.cells.iter().rposition(Option::is_some))
                .map(|(start, end)| (start, end + 1))
        })
        .collect::<Vec<_>>();
    let boundary_transition = prior_transition(prior_mm, open, prior_extend)?;
    let mut transitions = Vec::with_capacity(profile.width() + 1);
    for column in 0..profile.width() {
        let mut observed_mm = 0.0;
        let mut observed_md = 0.0;
        let mut observed_dm = 0.0;
        let mut observed_dd = 0.0;
        for ((row, weight), bounds) in profile.rows.iter().zip(weights).zip(&core_bounds) {
            let Some((start, end)) = *bounds else {
                continue;
            };
            if column < start || column >= end {
                continue;
            }
            let previous_match = column == start || row.cells[column - 1].is_some();
            let current_match = row.cells[column].is_some();
            match (previous_match, current_match) {
                (true, true) => observed_mm += *weight,
                (true, false) => observed_md += *weight,
                (false, true) => observed_dm += *weight,
                (false, false) => observed_dd += *weight,
            }
        }
        let mut m_family = [
            observed_mm + admixture * prior_mm,
            admixture * open,
            observed_md + admixture * open,
        ];
        normalize_positive(&mut m_family, "match-state transition family")?;
        let mut d_family = [
            observed_dm + admixture * (1.0 - prior_extend),
            observed_dd + admixture * prior_extend,
        ];
        normalize_positive(&mut d_family, "delete-state transition family")?;
        let mut i_family = [admixture * (1.0 - prior_extend), admixture * prior_extend];
        normalize_positive(&mut i_family, "insert-state transition family")?;
        transitions.push(Transition {
            mm: m_family[0],
            mi: m_family[1],
            md: m_family[2],
            im: i_family[0],
            ii: i_family[1],
            dm: d_family[0],
            dd: d_family[1],
        });
    }
    transitions.push(boundary_transition);
    Ok(transitions)
}

fn prior_transition(
    prior_mm: f64,
    open: f64,
    prior_extend: f64,
) -> Result<Transition, ClustalOmegaError> {
    let mut m_family = [prior_mm, open, open];
    normalize_positive(&mut m_family, "boundary match-state transition family")?;
    let mut side_family = [1.0 - prior_extend, prior_extend];
    normalize_positive(&mut side_family, "boundary side-state transition family")?;
    Ok(Transition {
        mm: m_family[0],
        mi: m_family[1],
        md: m_family[2],
        im: side_family[0],
        ii: side_family[1],
        dm: side_family[0],
        dd: side_family[1],
    })
}

#[allow(clippy::too_many_lines)]
fn viterbi_align(
    left: &ProfileHMM,
    right: &ProfileHMM,
    parameters: &ClustalOmegaParameters,
) -> Result<CoreAlignment, ClustalOmegaError> {
    let rows = left.len() + 1;
    let columns = right.len() + 1;
    let cells = checked_cells(rows, columns, "Viterbi matrix dimensions")?;
    let state_cells =
        cells
            .checked_mul(STATE_COUNT)
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "Viterbi state matrix dimensions",
            })?;
    let mut scores = vec![NEGATIVE_INFINITY; state_cells];
    let mut trace = vec![TRACE_START; state_cells];
    for row in 0..rows {
        scores[state_index(row, 0, columns, MM)] = 0.0;
    }
    for column in 0..columns {
        scores[state_index(0, column, columns, MM)] = 0.0;
    }
    let emission_scores = column_score_matrix(left, right, parameters)?;
    let factors = parameters.gap_path_factors();
    for row in 1..rows {
        for column in 1..columns {
            let emission = emission_scores[(row - 1) * right.len() + column - 1];
            let previous = (row - 1, column - 1);
            let mm_candidates = [
                (
                    MM,
                    scores[state_index(previous.0, previous.1, columns, MM)]
                        + close_to_match(left, right, row, column, MM),
                ),
                (
                    MI,
                    scores[state_index(previous.0, previous.1, columns, MI)]
                        + close_to_match(left, right, row, column, MI),
                ),
                (
                    IM,
                    scores[state_index(previous.0, previous.1, columns, IM)]
                        + close_to_match(left, right, row, column, IM),
                ),
                (
                    DG,
                    scores[state_index(previous.0, previous.1, columns, DG)]
                        + close_to_match(left, right, row, column, DG),
                ),
                (
                    GD,
                    scores[state_index(previous.0, previous.1, columns, GD)]
                        + close_to_match(left, right, row, column, GD),
                ),
            ];
            let (predecessor, best) = maximum_in_order(&mm_candidates);
            let index = state_index(row, column, columns, MM);
            scores[index] = best + emission;
            trace[index] = u8::try_from(predecessor).expect("five states fit in u8");

            let mi_opening = scores[state_index(row - 1, column, columns, MM)]
                + factors.2 * open_mi(left, right, row, column);
            let mi_extension = scores[state_index(row - 1, column, columns, MI)]
                + factors.3 * extend_mi(left, right, row, column);
            fill_side_state(
                &mut scores,
                &mut trace,
                row,
                column,
                columns,
                MI,
                mi_opening,
                mi_extension,
            );
            let dg_opening =
                scores[state_index(row - 1, column, columns, MM)] + factors.0 * open_dg(left, row);
            let dg_extension = scores[state_index(row - 1, column, columns, DG)]
                + factors.1 * extend_dg(left, row);
            fill_side_state(
                &mut scores,
                &mut trace,
                row,
                column,
                columns,
                DG,
                dg_opening,
                dg_extension,
            );
            let im_opening = scores[state_index(row, column - 1, columns, MM)]
                + factors.2 * open_im(left, right, row, column);
            let im_extension = scores[state_index(row, column - 1, columns, IM)]
                + factors.3 * extend_im(left, right, row, column);
            fill_side_state(
                &mut scores,
                &mut trace,
                row,
                column,
                columns,
                IM,
                im_opening,
                im_extension,
            );
            let gd_opening = scores[state_index(row, column - 1, columns, MM)]
                + factors.0 * open_gd(right, column);
            let gd_extension = scores[state_index(row, column - 1, columns, GD)]
                + factors.1 * extend_gd(right, column);
            fill_side_state(
                &mut scores,
                &mut trace,
                row,
                column,
                columns,
                GD,
                gd_opening,
                gd_extension,
            );
        }
    }
    let (mut row, mut column, score) = best_semiglobal_endpoint(&scores, rows, columns)?;
    let left_end = row;
    let right_end = column;
    let mut state = MM;
    let mut operations = Vec::with_capacity(row + column);
    while row > 0 && column > 0 {
        let index = state_index(row, column, columns, state);
        let predecessor = usize::from(trace[index]);
        match state {
            MM => {
                operations.push(ColumnOperation::Both);
                row -= 1;
                column -= 1;
            }
            MI | DG => {
                operations.push(ColumnOperation::Left);
                row -= 1;
            }
            IM | GD => {
                operations.push(ColumnOperation::Right);
                column -= 1;
            }
            _ => unreachable!("five-state index is exhaustive"),
        }
        if predecessor >= STATE_COUNT {
            return Err(ClustalOmegaError::InternalInvariant {
                message: "Viterbi traceback encountered an unset predecessor",
            });
        }
        state = predecessor;
    }
    operations.reverse();
    replay_core(&operations, row, column, left_end, right_end)?;
    let workspace_bytes = scores
        .len()
        .checked_mul(size_of::<f64>())
        .and_then(|bytes| bytes.checked_add(trace.len() * size_of::<u8>()))
        .and_then(|bytes| {
            emission_scores
                .len()
                .checked_mul(size_of::<f64>())
                .and_then(|emission_bytes| bytes.checked_add(emission_bytes))
        })
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "Viterbi workspace accounting",
        })?;
    Ok(CoreAlignment {
        score,
        left_start: row,
        left_end,
        right_start: column,
        right_end,
        operations,
        workspace_bytes,
        posterior_bounds: None,
    })
}

#[derive(Clone, Copy, Debug)]
struct ProbabilityTransitionWeights {
    mm: f64,
    im: f64,
    dm: f64,
    insert_open_mm: f64,
    insert_open_mi: f64,
    insert_extend_mm: f64,
    insert_extend_ii: f64,
    delete_open_md: f64,
    delete_extend_dd: f64,
}

fn probability_transition_weights(
    profile: &ProfileHMM,
    parameters: &ClustalOmegaParameters,
) -> Result<Vec<ProbabilityTransitionWeights>, ClustalOmegaError> {
    let factors = parameters.gap_path_factors();
    (0..=profile.len())
        .map(|position| {
            let transition = profile.probability_transition(position);
            let weights = ProbabilityTransitionWeights {
                mm: transition.mm,
                im: transition.im,
                dm: transition.dm,
                insert_open_mm: transition.mm.powf(factors.2),
                insert_open_mi: transition.mi.powf(factors.2),
                insert_extend_mm: transition.mm.powf(factors.3),
                insert_extend_ii: transition.ii.powf(factors.3),
                delete_open_md: transition.md.powf(factors.0),
                delete_extend_dd: transition.dd.powf(factors.1),
            };
            let finite = [
                weights.mm,
                weights.im,
                weights.dm,
                weights.insert_open_mm,
                weights.insert_open_mi,
                weights.insert_extend_mm,
                weights.insert_extend_ii,
                weights.delete_open_md,
                weights.delete_extend_dd,
            ]
            .into_iter()
            .all(|value| value.is_finite() && value > 0.0);
            if finite {
                Ok(weights)
            } else {
                Err(ClustalOmegaError::InvalidProbability {
                    stage: "scaled HMM transition weights",
                    value: f64::NAN,
                })
            }
        })
        .collect()
}

#[inline]
fn column_odds(left: &[f64], right: &[f64], inverse_background: &[f64]) -> f64 {
    debug_assert_eq!(left.len(), right.len());
    debug_assert_eq!(left.len(), inverse_background.len());
    match left.len() {
        20 => column_odds_20(left, right, inverse_background),
        4 => column_odds_4(left, right, inverse_background),
        _ => unreachable!("published profiles use four or twenty canonical symbols"),
    }
}

#[inline]
fn column_odds_20(left: &[f64], right: &[f64], inverse_background: &[f64]) -> f64 {
    let mut sums = [0.0_f64; 4];
    macro_rules! accumulate {
        ($offset:expr) => {{
            sums[0] += left[$offset] * right[$offset] * inverse_background[$offset];
            sums[1] += left[$offset + 1] * right[$offset + 1] * inverse_background[$offset + 1];
            sums[2] += left[$offset + 2] * right[$offset + 2] * inverse_background[$offset + 2];
            sums[3] += left[$offset + 3] * right[$offset + 3] * inverse_background[$offset + 3];
        }};
    }
    accumulate!(0);
    accumulate!(4);
    accumulate!(8);
    accumulate!(12);
    accumulate!(16);
    (sums[0] + sums[1]) + (sums[2] + sums[3])
}

#[inline]
fn column_odds_4(left: &[f64], right: &[f64], inverse_background: &[f64]) -> f64 {
    let mut sums = [0.0_f64; 4];
    sums[0] += left[0] * right[0] * inverse_background[0];
    sums[1] += left[1] * right[1] * inverse_background[1];
    sums[2] += left[2] * right[2] * inverse_background[2];
    sums[3] += left[3] * right[3] * inverse_background[3];
    (sums[0] + sums[1]) + (sums[2] + sums[3])
}

fn scaled_partition_log(
    forward_matches: &[f64],
    scales: &[f64],
    left_length: usize,
    right_length: usize,
) -> Result<f64, ClustalOmegaError> {
    let columns = right_length + 1;
    let mut maximum_scale = scales[left_length];
    for scale in scales.iter().take(left_length).skip(1) {
        maximum_scale = maximum_scale.max(*scale);
    }
    let mut sum = 0.0;
    for column in 1..=right_length {
        sum += forward_matches[left_length * columns + column]
            * (scales[left_length] - maximum_scale).exp2();
    }
    for row in 1..left_length {
        sum += forward_matches[row * columns + right_length] * (scales[row] - maximum_scale).exp2();
    }
    let partition = maximum_scale + sum.log2();
    if partition.is_finite() {
        Ok(partition)
    } else {
        Err(ClustalOmegaError::InvalidProbability {
            stage: "scaled Forward semi-global partition",
            value: partition,
        })
    }
}

const SCALE_HIGH_WATER: f64 = 1.0e60;
const SCALE_LOW_WATER: f64 = 1.0e-60;

fn scale_adjustment(maximum: f64) -> Option<f64> {
    if (SCALE_LOW_WATER..=SCALE_HIGH_WATER).contains(&maximum) {
        None
    } else {
        Some(maximum.log2().floor())
    }
}

#[allow(clippy::too_many_lines)]
fn mac_align_scaled(
    left: &ProfileHMM,
    right: &ProfileHMM,
    parameters: &ClustalOmegaParameters,
    prefer_left_gap: bool,
) -> Result<CoreAlignment, ClustalOmegaError> {
    let left_length = left.len();
    let right_length = right.len();
    let rows = left_length + 1;
    let columns = right_length + 1;
    let cells = checked_cells(rows, columns, "scaled Forward matrix dimensions")?;
    let left_weights = probability_transition_weights(left, parameters)?;
    let right_weights = probability_transition_weights(right, parameters)?;
    let emission_cells = checked_cells(
        left_length,
        right_length,
        "HMM column-weight matrix dimensions",
    )?;
    let mut emission_weights = Vec::with_capacity(emission_cells);
    let inverse_background = background(left.alphabet_size == 20)
        .into_iter()
        .map(f64::recip)
        .collect::<Vec<_>>();
    let emission_shift = parameters.column_score_shift().exp2();
    let mut forward_matches = vec![0.0; cells];
    let mut forward_scales = vec![0.0_f64; rows];
    let mut previous_states = vec![[0.0; STATE_COUNT]; columns];
    let mut current_states = vec![[0.0; STATE_COUNT]; columns];
    for cell in &mut previous_states {
        cell[MM] = 1.0;
    }
    forward_matches[..columns].fill(1.0);
    for row in 1..rows {
        current_states.fill([0.0; STATE_COUNT]);
        let previous_scale = forward_scales[row - 1];
        let base_scale = previous_scale.max(0.0);
        let previous_factor = (previous_scale - base_scale).exp2();
        let boundary_factor = (-base_scale).exp2();
        let row_start = row * columns;
        current_states[0][MM] = boundary_factor;
        let mut maximum = boundary_factor;
        let left_emission = left.emission(row - 1);
        for column in 1..columns {
            let close_left = left_weights[row - 1];
            let close_right = right_weights[column - 1];
            let emission = column_odds(
                left_emission,
                right.emission(column - 1),
                &inverse_background,
            ) * emission_shift;
            if !emission.is_finite() || emission <= 0.0 {
                return Err(ClustalOmegaError::InvalidProbability {
                    stage: "scaled HMM-HMM column weight",
                    value: emission,
                });
            }
            emission_weights.push(emission);
            let diagonal_states = previous_states[column - 1];
            let close_sum = close_right.mm
                * (diagonal_states[MM] * close_left.mm
                    + diagonal_states[IM] * close_left.im
                    + diagonal_states[DG] * close_left.dm)
                + close_left.mm
                    * (diagonal_states[MI] * close_right.im + diagonal_states[GD] * close_right.dm);
            current_states[column][MM] = previous_factor * emission * close_sum;
            let scaled_up_match = previous_factor * previous_states[column][MM];
            current_states[column][MI] = scaled_up_match
                * (left_weights[row - 1].insert_open_mm * right_weights[column].insert_open_mi)
                + (previous_factor * previous_states[column][MI])
                    * (left_weights[row - 1].insert_extend_mm
                        * right_weights[column].insert_extend_ii);
            current_states[column][DG] = scaled_up_match * left_weights[row - 1].delete_open_md
                + (previous_factor * previous_states[column][DG])
                    * left_weights[row - 1].delete_extend_dd;
            current_states[column][IM] = current_states[column - 1][MM]
                * (left_weights[row].insert_open_mi * right_weights[column - 1].insert_open_mm)
                + current_states[column - 1][IM]
                    * (left_weights[row].insert_extend_ii
                        * right_weights[column - 1].insert_extend_mm);
            current_states[column][GD] = current_states[column - 1][MM]
                * right_weights[column - 1].delete_open_md
                + current_states[column - 1][GD] * right_weights[column - 1].delete_extend_dd;
            maximum = current_states[column]
                .iter()
                .copied()
                .fold(maximum, f64::max);
        }
        if !maximum.is_finite() || maximum <= 0.0 {
            return Err(ClustalOmegaError::InvalidProbability {
                stage: "scaled Forward row",
                value: maximum,
            });
        }
        let adjustment = scale_adjustment(maximum);
        if let Some(adjustment) = adjustment {
            let multiplier = (-adjustment).exp2();
            for cell in &mut current_states {
                for value in cell {
                    *value *= multiplier;
                }
            }
        }
        for (destination, states) in forward_matches[row_start..row_start + columns]
            .iter_mut()
            .zip(&current_states)
        {
            *destination = states[MM];
        }
        forward_scales[row] = base_scale + adjustment.unwrap_or(0.0);
        core::mem::swap(&mut previous_states, &mut current_states);
    }

    let partition =
        scaled_partition_log(&forward_matches, &forward_scales, left_length, right_length)?;
    let mut minimum_posterior = 1.0_f64;
    let mut maximum_posterior = 0.0_f64;
    previous_states.fill([0.0; STATE_COUNT]);
    current_states.fill([0.0; STATE_COUNT]);
    let mut next = previous_states;
    let mut current = current_states;
    for value in current.iter_mut().skip(1) {
        value[MM] = 1.0;
    }
    let posterior_log_scale = forward_scales[left_length] - partition;
    let posterior_scale = if posterior_log_scale < -1_100.0 {
        0.0
    } else {
        posterior_log_scale.exp2()
    };
    for column in 1..columns {
        let probability =
            (forward_matches[left_length * columns + column] * posterior_scale).clamp(0.0, 1.0);
        if !probability.is_finite() {
            return Err(ClustalOmegaError::InvalidProbability {
                stage: "scaled match-state posterior",
                value: probability,
            });
        }
        forward_matches[left_length * columns + column] = probability;
        minimum_posterior = minimum_posterior.min(probability);
        maximum_posterior = maximum_posterior.max(probability);
    }
    core::mem::swap(&mut next, &mut current);
    let mut next_scale = Some(0.0);

    for row in (0..left_length).rev() {
        current.fill([0.0; STATE_COUNT]);
        let has_endpoint = row > 0;
        let mut base_scale = if has_endpoint { 0.0 } else { NEGATIVE_INFINITY };
        if let Some(scale) = next_scale {
            base_scale = base_scale.max(scale);
        }
        let next_factor = next_scale.map_or(0.0, |scale| (scale - base_scale).exp2());
        let boundary_factor = if has_endpoint {
            (-base_scale).exp2()
        } else {
            0.0
        };
        let mut maximum = 0.0_f64;
        for column in (0..columns).rev() {
            let mut value = [0.0; STATE_COUNT];
            let endpoint = column == right_length && row > 0;
            if endpoint {
                value[MM] += boundary_factor;
            }
            if column < right_length {
                let target = next[column + 1];
                let close_left = left_weights[row];
                let close_right = right_weights[column];
                let emission = emission_weights[row * right_length + column];
                let target_match = target[MM];
                let match_weight = next_factor * target_match * emission;
                let match_left = match_weight * close_left.mm;
                let match_right = match_weight * close_right.mm;
                value[MM] += match_left * close_right.mm;
                value[MI] = match_left * close_right.im
                    + next_factor
                        * next[column][MI]
                        * (left_weights[row].insert_extend_mm
                            * right_weights[column].insert_extend_ii);
                value[DG] = match_right * close_left.dm
                    + next_factor * next[column][DG] * left_weights[row].delete_extend_dd;
                value[IM] = match_right * close_left.im
                    + current[column + 1][IM]
                        * (left_weights[row].insert_extend_ii
                            * right_weights[column].insert_extend_mm);
                value[GD] = match_left * close_right.dm
                    + current[column + 1][GD] * right_weights[column].delete_extend_dd;
            }
            if column > 0 {
                let target = next[column];
                value[MM] += next_factor
                    * (target[MI]
                        * (left_weights[row].insert_open_mm
                            * right_weights[column].insert_open_mi)
                        + target[DG] * left_weights[row].delete_open_md);
            }
            if column < right_length && row > 0 {
                let target = current[column + 1];
                value[MM] += target[IM]
                    * (left_weights[row].insert_open_mi * right_weights[column].insert_open_mm)
                    + target[GD] * right_weights[column].delete_open_md;
            }
            maximum = value.iter().copied().fold(maximum, f64::max);
            current[column] = value;
        }
        if !maximum.is_finite() {
            return Err(ClustalOmegaError::InvalidProbability {
                stage: "scaled Backward row",
                value: maximum,
            });
        }
        let adjustment = if maximum > 0.0 {
            scale_adjustment(maximum)
        } else {
            None
        };
        if let Some(adjustment) = adjustment {
            let multiplier = (-adjustment).exp2();
            for value in &mut current {
                for state in value {
                    *state *= multiplier;
                }
            }
        }
        let current_scale = base_scale + adjustment.unwrap_or(0.0);
        if row > 0 {
            let posterior_log_scale = forward_scales[row] + current_scale - partition;
            let posterior_scale = if posterior_log_scale < -1_100.0 {
                0.0
            } else {
                posterior_log_scale.exp2()
            };
            for column in 1..columns {
                let probability = (forward_matches[row * columns + column]
                    * current[column][MM]
                    * posterior_scale)
                    .clamp(0.0, 1.0);
                if !probability.is_finite() {
                    return Err(ClustalOmegaError::InvalidProbability {
                        stage: "scaled match-state posterior",
                        value: probability,
                    });
                }
                forward_matches[row * columns + column] = probability;
                minimum_posterior = minimum_posterior.min(probability);
                maximum_posterior = maximum_posterior.max(probability);
            }
        }
        core::mem::swap(&mut next, &mut current);
        next_scale = Some(current_scale);
    }

    let threshold = parameters.mac_threshold();
    let mut previous_accuracy = vec![0.0; columns];
    let mut current_accuracy = vec![0.0; columns];
    let mut trace = vec![TRACE_START; cells];
    let mut best = (0.0, 0_usize, 0_usize);
    for row in 1..rows {
        current_accuracy.fill(0.0);
        for column in 1..columns {
            let diagonal =
                previous_accuracy[column - 1] + forward_matches[row * columns + column] - threshold;
            let up = previous_accuracy[column] - threshold / 2.0;
            let left_score = current_accuracy[column - 1] - threshold / 2.0;
            let candidates = if prefer_left_gap {
                [(diagonal, 0_u8), (up, 1), (left_score, 2), (0.0, 3)]
            } else {
                [(diagonal, 0_u8), (left_score, 2), (up, 1), (0.0, 3)]
            };
            let mut selected = candidates[0];
            for candidate in candidates.iter().copied().skip(1) {
                if candidate.0 > selected.0 {
                    selected = candidate;
                }
            }
            current_accuracy[column] = selected.0;
            let index = row * columns + column;
            trace[index] = if selected.0 > 0.0 {
                selected.1
            } else {
                TRACE_START
            };
            let candidate_key = (
                selected.0,
                row + column,
                core::cmp::Reverse(row),
                core::cmp::Reverse(column),
            );
            let best_key = (
                best.0,
                best.1 + best.2,
                core::cmp::Reverse(best.1),
                core::cmp::Reverse(best.2),
            );
            if candidate_key > best_key {
                best = (selected.0, row, column);
            }
        }
        core::mem::swap(&mut previous_accuracy, &mut current_accuracy);
    }
    let mut row = best.1;
    let mut column = best.2;
    let left_end = row;
    let right_end = column;
    let mut operations = Vec::with_capacity(row + column);
    while row > 0 && column > 0 {
        let operation =
            *trace
                .get(row * columns + column)
                .ok_or(ClustalOmegaError::InternalInvariant {
                    message: "scaled MAC traceback index lies outside its matrix",
                })?;
        match operation {
            0 => {
                operations.push(ColumnOperation::Both);
                row -= 1;
                column -= 1;
            }
            1 => {
                operations.push(ColumnOperation::Left);
                row -= 1;
            }
            2 => {
                operations.push(ColumnOperation::Right);
                column -= 1;
            }
            3 | TRACE_START => break,
            _ => {
                return Err(ClustalOmegaError::InternalInvariant {
                    message: "scaled MAC traceback contains an unknown operation",
                });
            }
        }
    }
    operations.reverse();
    replay_core(&operations, row, column, left_end, right_end)?;
    Ok(CoreAlignment {
        score: best.0,
        left_start: row,
        left_end,
        right_start: column,
        right_end,
        operations,
        workspace_bytes: mac_workspace_bytes(left_length, right_length)?,
        posterior_bounds: Some((minimum_posterior, maximum_posterior)),
    })
}

#[allow(clippy::too_many_lines)]
#[cfg(test)]
fn mac_align_logspace(
    left: &ProfileHMM,
    right: &ProfileHMM,
    parameters: &ClustalOmegaParameters,
    prefer_left_gap: bool,
) -> Result<CoreAlignment, ClustalOmegaError> {
    let rows = left.len() + 1;
    let columns = right.len() + 1;
    let cells = checked_cells(rows, columns, "Forward matrix dimensions")?;
    let state_cells =
        cells
            .checked_mul(STATE_COUNT)
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "Forward state matrix dimensions",
            })?;
    let mut forward = vec![NEGATIVE_INFINITY; state_cells];
    for row in 0..rows {
        forward[state_index(row, 0, columns, MM)] = 0.0;
    }
    for column in 0..columns {
        forward[state_index(0, column, columns, MM)] = 0.0;
    }
    let emission_scores = column_score_matrix(left, right, parameters)?;
    let factors = parameters.gap_path_factors();
    for row in 1..rows {
        for column in 1..columns {
            let emission = emission_scores[(row - 1) * right.len() + column - 1];
            let diagonal = (row - 1, column - 1);
            let predecessors = [
                forward[state_index(diagonal.0, diagonal.1, columns, MM)]
                    + close_to_match(left, right, row, column, MM),
                forward[state_index(diagonal.0, diagonal.1, columns, MI)]
                    + close_to_match(left, right, row, column, MI),
                forward[state_index(diagonal.0, diagonal.1, columns, IM)]
                    + close_to_match(left, right, row, column, IM),
                forward[state_index(diagonal.0, diagonal.1, columns, DG)]
                    + close_to_match(left, right, row, column, DG),
                forward[state_index(diagonal.0, diagonal.1, columns, GD)]
                    + close_to_match(left, right, row, column, GD),
            ];
            forward[state_index(row, column, columns, MM)] = logsumexp2(&predecessors) + emission;
            forward[state_index(row, column, columns, MI)] = logadd2(
                forward[state_index(row - 1, column, columns, MM)]
                    + factors.2 * open_mi(left, right, row, column),
                forward[state_index(row - 1, column, columns, MI)]
                    + factors.3 * extend_mi(left, right, row, column),
            );
            forward[state_index(row, column, columns, DG)] = logadd2(
                forward[state_index(row - 1, column, columns, MM)] + factors.0 * open_dg(left, row),
                forward[state_index(row - 1, column, columns, DG)]
                    + factors.1 * extend_dg(left, row),
            );
            forward[state_index(row, column, columns, IM)] = logadd2(
                forward[state_index(row, column - 1, columns, MM)]
                    + factors.2 * open_im(left, right, row, column),
                forward[state_index(row, column - 1, columns, IM)]
                    + factors.3 * extend_im(left, right, row, column),
            );
            forward[state_index(row, column, columns, GD)] = logadd2(
                forward[state_index(row, column - 1, columns, MM)]
                    + factors.0 * open_gd(right, column),
                forward[state_index(row, column - 1, columns, GD)]
                    + factors.1 * extend_gd(right, column),
            );
        }
    }
    let partition = semiglobal_partition(&forward, rows, columns);
    if !partition.is_finite() {
        return Err(ClustalOmegaError::InvalidProbability {
            stage: "Forward semi-global partition",
            value: partition,
        });
    }
    let mut posterior = vec![0.0; cells];
    let row_state_count =
        columns
            .checked_mul(STATE_COUNT)
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "Backward rolling-row dimensions",
            })?;
    let mut next = vec![NEGATIVE_INFINITY; row_state_count];
    let mut current = vec![NEGATIVE_INFINITY; row_state_count];
    let mut minimum_posterior = 1.0_f64;
    let mut maximum_posterior = 0.0_f64;
    for row in (0..rows).rev() {
        current.fill(NEGATIVE_INFINITY);
        for column in (0..columns).rev() {
            let endpoint = (row == rows - 1 && column > 0) || (column == columns - 1 && row > 0);
            let mut mm_outgoing = [NEGATIVE_INFINITY; 6];
            let mut outgoing_count = 0;
            if endpoint {
                mm_outgoing[outgoing_count] = 0.0;
                outgoing_count += 1;
            }
            if row + 1 < rows && column + 1 < columns {
                let emission = emission_scores[row * right.len() + column];
                mm_outgoing[outgoing_count] = next[(column + 1) * STATE_COUNT + MM]
                    + close_to_match(left, right, row + 1, column + 1, MM)
                    + emission;
                outgoing_count += 1;
            }
            if row + 1 < rows && column > 0 {
                mm_outgoing[outgoing_count] = next[column * STATE_COUNT + MI]
                    + factors.2 * open_mi(left, right, row + 1, column);
                outgoing_count += 1;
                mm_outgoing[outgoing_count] =
                    next[column * STATE_COUNT + DG] + factors.0 * open_dg(left, row + 1);
                outgoing_count += 1;
            }
            if column + 1 < columns && row > 0 {
                mm_outgoing[outgoing_count] = current[(column + 1) * STATE_COUNT + IM]
                    + factors.2 * open_im(left, right, row, column + 1);
                outgoing_count += 1;
                mm_outgoing[outgoing_count] = current[(column + 1) * STATE_COUNT + GD]
                    + factors.0 * open_gd(right, column + 1);
                outgoing_count += 1;
            }
            current[column * STATE_COUNT + MM] = logsumexp2(&mm_outgoing[..outgoing_count]);

            if row + 1 < rows && column + 1 < columns {
                let emission = emission_scores[row * right.len() + column];
                current[column * STATE_COUNT + MI] = logadd2(
                    next[(column + 1) * STATE_COUNT + MM]
                        + close_to_match(left, right, row + 1, column + 1, MI)
                        + emission,
                    next[column * STATE_COUNT + MI]
                        + factors.3 * extend_mi(left, right, row + 1, column),
                );
                current[column * STATE_COUNT + DG] = logadd2(
                    next[(column + 1) * STATE_COUNT + MM]
                        + close_to_match(left, right, row + 1, column + 1, DG)
                        + emission,
                    next[column * STATE_COUNT + DG] + factors.1 * extend_dg(left, row + 1),
                );
                current[column * STATE_COUNT + IM] = logadd2(
                    next[(column + 1) * STATE_COUNT + MM]
                        + close_to_match(left, right, row + 1, column + 1, IM)
                        + emission,
                    current[(column + 1) * STATE_COUNT + IM]
                        + factors.3 * extend_im(left, right, row, column + 1),
                );
                current[column * STATE_COUNT + GD] = logadd2(
                    next[(column + 1) * STATE_COUNT + MM]
                        + close_to_match(left, right, row + 1, column + 1, GD)
                        + emission,
                    current[(column + 1) * STATE_COUNT + GD]
                        + factors.1 * extend_gd(right, column + 1),
                );
            }

            if row > 0 && column > 0 {
                let log_probability = forward[state_index(row, column, columns, MM)]
                    + current[column * STATE_COUNT + MM]
                    - partition;
                let probability = if log_probability < -1_100.0 {
                    0.0
                } else {
                    log_probability.exp2().clamp(0.0, 1.0)
                };
                if !probability.is_finite() {
                    return Err(ClustalOmegaError::InvalidProbability {
                        stage: "match-state posterior",
                        value: probability,
                    });
                }
                posterior[row * columns + column] = probability;
                minimum_posterior = minimum_posterior.min(probability);
                maximum_posterior = maximum_posterior.max(probability);
            }
        }
        core::mem::swap(&mut current, &mut next);
    }

    let threshold = parameters.mac_threshold();
    let mut accuracy = vec![0.0; cells];
    let mut trace = vec![TRACE_START; cells];
    let mut best = (0.0, 0_usize, 0_usize);
    for row in 1..rows {
        for column in 1..columns {
            let diagonal = accuracy[(row - 1) * columns + column - 1]
                + posterior[row * columns + column]
                - threshold;
            let up = accuracy[(row - 1) * columns + column] - threshold / 2.0;
            let left_score = accuracy[row * columns + column - 1] - threshold / 2.0;
            let candidates = if prefer_left_gap {
                [(diagonal, 0_u8), (up, 1), (left_score, 2), (0.0, 3)]
            } else {
                [(diagonal, 0_u8), (left_score, 2), (up, 1), (0.0, 3)]
            };
            let mut selected = candidates[0];
            for candidate in candidates.iter().copied().skip(1) {
                if candidate.0 > selected.0 {
                    selected = candidate;
                }
            }
            let index = row * columns + column;
            accuracy[index] = selected.0;
            trace[index] = selected.1;
            let candidate_key = (
                selected.0,
                row + column,
                core::cmp::Reverse(row),
                core::cmp::Reverse(column),
            );
            let best_key = (
                best.0,
                best.1 + best.2,
                core::cmp::Reverse(best.1),
                core::cmp::Reverse(best.2),
            );
            if candidate_key > best_key {
                best = (selected.0, row, column);
            }
        }
    }
    let mut row = best.1;
    let mut column = best.2;
    let left_end = row;
    let right_end = column;
    let mut operations = Vec::with_capacity(row + column);
    while row > 0 && column > 0 {
        let trace_operation =
            *trace
                .get(row * columns + column)
                .ok_or(ClustalOmegaError::InternalInvariant {
                    message: "MAC traceback index lies outside its matrix",
                })?;
        match trace_operation {
            0 => {
                operations.push(ColumnOperation::Both);
                row -= 1;
                column -= 1;
            }
            1 => {
                operations.push(ColumnOperation::Left);
                row -= 1;
            }
            2 => {
                operations.push(ColumnOperation::Right);
                column -= 1;
            }
            3 | TRACE_START => break,
            _ => {
                return Err(ClustalOmegaError::InternalInvariant {
                    message: "MAC traceback contains an unknown operation",
                });
            }
        }
        if accuracy[row * columns + column] <= 0.0 {
            break;
        }
    }
    operations.reverse();
    replay_core(&operations, row, column, left_end, right_end)?;
    let workspace_bytes = mac_workspace_bytes(left.len(), right.len())?;
    Ok(CoreAlignment {
        score: best.0,
        left_start: row,
        left_end,
        right_start: column,
        right_end,
        operations,
        workspace_bytes,
        posterior_bounds: Some((minimum_posterior, maximum_posterior)),
    })
}

#[cfg(test)]
fn column_score(
    left: &ProfileHMM,
    right: &ProfileHMM,
    row: usize,
    column: usize,
    background: &[f64],
    parameters: &ClustalOmegaParameters,
) -> Result<f64, ClustalOmegaError> {
    let odds = left
        .emission(row - 1)
        .iter()
        .zip(right.emission(column - 1))
        .zip(background)
        .try_fold(0.0, |sum, ((left, right), background)| {
            let contribution = left * right / background;
            let next = sum + contribution;
            if next.is_finite() {
                Ok(next)
            } else {
                Err(ClustalOmegaError::InvalidProbability {
                    stage: "HMM-HMM column odds",
                    value: next,
                })
            }
        })?;
    if odds <= 0.0 {
        return Err(ClustalOmegaError::InvalidProbability {
            stage: "HMM-HMM column odds",
            value: odds,
        });
    }
    Ok(odds.log2() + parameters.column_score_shift())
}

fn column_score_matrix(
    left: &ProfileHMM,
    right: &ProfileHMM,
    parameters: &ClustalOmegaParameters,
) -> Result<Vec<f64>, ClustalOmegaError> {
    let cells = checked_cells(
        left.len(),
        right.len(),
        "HMM column-score matrix dimensions",
    )?;
    let mut scores = Vec::with_capacity(cells);
    let inverse_background = background(left.alphabet_size == 20)
        .into_iter()
        .map(f64::recip)
        .collect::<Vec<_>>();
    for row in 0..left.len() {
        for column in 0..right.len() {
            let odds = column_odds(
                left.emission(row),
                right.emission(column),
                &inverse_background,
            );
            if !odds.is_finite() || odds <= 0.0 {
                return Err(ClustalOmegaError::InvalidProbability {
                    stage: "HMM-HMM column odds",
                    value: odds,
                });
            }
            scores.push(odds.log2() + parameters.column_score_shift());
        }
    }
    Ok(scores)
}

fn close_to_match(
    left: &ProfileHMM,
    right: &ProfileHMM,
    row: usize,
    column: usize,
    predecessor: usize,
) -> f64 {
    let left = left.transition(row - 1);
    let right = right.transition(column - 1);
    match predecessor {
        MM => left.mm + right.mm,
        MI => left.mm + right.im,
        IM => left.im + right.mm,
        DG => left.dm + right.mm,
        GD => left.mm + right.dm,
        _ => NEGATIVE_INFINITY,
    }
}

fn open_mi(left: &ProfileHMM, right: &ProfileHMM, row: usize, column: usize) -> f64 {
    left.transition(row - 1).mm + right.transition(column).mi
}

fn extend_mi(left: &ProfileHMM, right: &ProfileHMM, row: usize, column: usize) -> f64 {
    left.transition(row - 1).mm + right.transition(column).ii
}

fn open_im(left: &ProfileHMM, right: &ProfileHMM, row: usize, column: usize) -> f64 {
    left.transition(row).mi + right.transition(column - 1).mm
}

fn extend_im(left: &ProfileHMM, right: &ProfileHMM, row: usize, column: usize) -> f64 {
    left.transition(row).ii + right.transition(column - 1).mm
}

fn open_dg(left: &ProfileHMM, row: usize) -> f64 {
    left.transition(row - 1).md
}

fn extend_dg(left: &ProfileHMM, row: usize) -> f64 {
    left.transition(row - 1).dd
}

fn open_gd(right: &ProfileHMM, column: usize) -> f64 {
    right.transition(column - 1).md
}

fn extend_gd(right: &ProfileHMM, column: usize) -> f64 {
    right.transition(column - 1).dd
}

#[allow(clippy::too_many_arguments)]
fn fill_side_state(
    scores: &mut [f64],
    trace: &mut [u8],
    row: usize,
    column: usize,
    columns: usize,
    state: usize,
    opening: f64,
    extension: f64,
) {
    let index = state_index(row, column, columns, state);
    if extension >= opening {
        scores[index] = extension;
        trace[index] = u8::try_from(state).expect("five states fit in u8");
    } else {
        scores[index] = opening;
        trace[index] = u8::try_from(MM).expect("five states fit in u8");
    }
}

fn maximum_in_order(candidates: &[(usize, f64)]) -> (usize, f64) {
    let mut selected = candidates[0];
    for candidate in candidates.iter().copied().skip(1) {
        if candidate.1 > selected.1 {
            selected = candidate;
        }
    }
    selected
}

fn best_semiglobal_endpoint(
    scores: &[f64],
    rows: usize,
    columns: usize,
) -> Result<(usize, usize, f64), ClustalOmegaError> {
    let mut selected: Option<(f64, usize, core::cmp::Reverse<usize>, usize, usize)> = None;
    let last_row = rows - 1;
    let last_column = columns - 1;
    for column in 1..columns {
        let score = scores[state_index(last_row, column, columns, MM)];
        let candidate = (
            score,
            last_row + column,
            core::cmp::Reverse(last_row),
            last_row,
            column,
        );
        if selected.as_ref().is_none_or(|current| candidate > *current) {
            selected = Some(candidate);
        }
    }
    for row in 1..last_row {
        let score = scores[state_index(row, last_column, columns, MM)];
        let candidate = (
            score,
            row + last_column,
            core::cmp::Reverse(row),
            row,
            last_column,
        );
        if selected.as_ref().is_none_or(|current| candidate > *current) {
            selected = Some(candidate);
        }
    }
    selected
        .map(|(score, _, _, row, column)| (row, column, score))
        .ok_or(ClustalOmegaError::InternalInvariant {
            message: "Viterbi has no legal semi-global endpoint",
        })
}

#[cfg(test)]
fn semiglobal_partition(forward: &[f64], rows: usize, columns: usize) -> f64 {
    let last_row = rows - 1;
    let last_column = columns - 1;
    let mut endpoints = Vec::with_capacity(rows + columns);
    for column in 1..columns {
        endpoints.push(forward[state_index(last_row, column, columns, MM)]);
    }
    for row in 1..last_row {
        endpoints.push(forward[state_index(row, last_column, columns, MM)]);
    }
    logsumexp2(&endpoints)
}

#[cfg(test)]
fn logsumexp2(values: &[f64]) -> f64 {
    let maximum = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .max_by(f64::total_cmp)
        .unwrap_or(NEGATIVE_INFINITY);
    if !maximum.is_finite() {
        return NEGATIVE_INFINITY;
    }
    maximum
        + values
            .iter()
            .filter(|value| value.is_finite())
            .map(|value| (*value - maximum).exp2())
            .sum::<f64>()
            .log2()
}

#[cfg(test)]
fn logadd2(first: f64, second: f64) -> f64 {
    logsumexp2(&[first, second])
}

fn attach_terminals(
    left_width: usize,
    right_width: usize,
    left_start: usize,
    left_end: usize,
    right_start: usize,
    right_end: usize,
    core: &[ColumnOperation],
) -> Result<Vec<ColumnOperation>, ClustalOmegaError> {
    if left_start > left_end
        || left_end > left_width
        || right_start > right_end
        || right_end > right_width
    {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "HMM core interval lies outside a child profile",
        });
    }
    if core.is_empty() {
        return aligned_terminal_block(left_width, right_width, true);
    }
    let mut operations = aligned_terminal_block(left_start, right_start, true)?;
    operations.extend_from_slice(core);
    operations.extend(aligned_terminal_block(
        left_width - left_end,
        right_width - right_end,
        false,
    )?);
    Ok(operations)
}

fn aligned_terminal_block(
    left_count: usize,
    right_count: usize,
    right_align: bool,
) -> Result<Vec<ColumnOperation>, ClustalOmegaError> {
    let width = left_count.max(right_count);
    let mut operations = Vec::with_capacity(width);
    for output in 0..width {
        let left_present = if right_align {
            output >= width - left_count
        } else {
            output < left_count
        };
        let right_present = if right_align {
            output >= width - right_count
        } else {
            output < right_count
        };
        operations.push(match (left_present, right_present) {
            (true, true) => ColumnOperation::Both,
            (true, false) => ColumnOperation::Left,
            (false, true) => ColumnOperation::Right,
            (false, false) => {
                return Err(ClustalOmegaError::InternalInvariant {
                    message: "terminal reattachment produced an empty column",
                });
            }
        });
    }
    Ok(operations)
}

fn merge_profiles<A: Alphabet>(
    left: WorkingProfile<A>,
    right: WorkingProfile<A>,
    operations: &[ColumnOperation],
) -> Result<WorkingProfile<A>, ClustalOmegaError> {
    let left_width = left.width();
    let right_width = right.width();
    let mut left_position = 0;
    let mut right_position = 0;
    let mut left_rows = left
        .rows
        .into_iter()
        .map(|row| {
            (
                row.sequence_index,
                row.cells,
                Vec::with_capacity(operations.len()),
            )
        })
        .collect::<Vec<_>>();
    let mut right_rows = right
        .rows
        .into_iter()
        .map(|row| {
            (
                row.sequence_index,
                row.cells,
                Vec::with_capacity(operations.len()),
            )
        })
        .collect::<Vec<_>>();
    for operation in operations {
        let consume_left = matches!(operation, ColumnOperation::Both | ColumnOperation::Left);
        let consume_right = matches!(operation, ColumnOperation::Both | ColumnOperation::Right);
        if consume_left && left_position >= left_width
            || consume_right && right_position >= right_width
        {
            return Err(ClustalOmegaError::InternalInvariant {
                message: "profile merge operation consumes beyond a child width",
            });
        }
        for (_, source, destination) in &mut left_rows {
            destination.push(consume_left.then(|| source[left_position]).flatten());
        }
        for (_, source, destination) in &mut right_rows {
            destination.push(consume_right.then(|| source[right_position]).flatten());
        }
        left_position += usize::from(consume_left);
        right_position += usize::from(consume_right);
    }
    if left_position != left_width || right_position != right_width {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "profile merge did not consume both child profiles exactly",
        });
    }
    let rows = left_rows
        .into_iter()
        .chain(right_rows)
        .map(|(sequence_index, _, cells)| WorkingRow {
            sequence_index,
            cells,
        })
        .collect::<Vec<_>>();
    let profile = WorkingProfile { rows };
    debug_assert!(validate_profile(&profile).is_ok());
    Ok(profile)
}

fn validate_profile<A: Alphabet>(profile: &WorkingProfile<A>) -> Result<(), ClustalOmegaError> {
    if profile.rows.is_empty() {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "working profile has no rows",
        });
    }
    let width = profile.rows[0].cells.len();
    if profile.rows.iter().any(|row| row.cells.len() != width) {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "working profile rows have different widths",
        });
    }
    for column in 0..width {
        if profile.rows.iter().all(|row| row.cells[column].is_none()) {
            return Err(ClustalOmegaError::InternalInvariant {
                message: "working profile contains an all-gap column",
            });
        }
    }
    Ok(())
}

fn replay_core(
    operations: &[ColumnOperation],
    left_start: usize,
    right_start: usize,
    left_end: usize,
    right_end: usize,
) -> Result<(), ClustalOmegaError> {
    let mut left = left_start;
    let mut right = right_start;
    for operation in operations {
        left += usize::from(matches!(
            operation,
            ColumnOperation::Both | ColumnOperation::Left
        ));
        right += usize::from(matches!(
            operation,
            ColumnOperation::Both | ColumnOperation::Right
        ));
    }
    if left != left_end || right != right_end {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "independent HMM traceback replay disagrees with its endpoint",
        });
    }
    Ok(())
}

fn mac_workspace_bytes(
    left_length: usize,
    right_length: usize,
) -> Result<usize, ClustalOmegaError> {
    let rows = left_length
        .checked_add(1)
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "MAC workspace rows",
        })?;
    let columns = right_length
        .checked_add(1)
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "MAC workspace columns",
        })?;
    let cells = checked_cells(rows, columns, "MAC workspace dimensions")?;
    let forward_match_bytes =
        cells
            .checked_mul(size_of::<f64>())
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "MAC Forward match-state workspace accounting",
            })?;
    let trace_bytes =
        cells
            .checked_mul(size_of::<u8>())
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "MAC traceback workspace accounting",
            })?;
    let emission_bytes = left_length
        .checked_mul(right_length)
        .and_then(|values| values.checked_mul(size_of::<f64>()))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "MAC emission workspace accounting",
        })?;
    let scale_bytes =
        rows.checked_mul(size_of::<f64>())
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "MAC scale workspace accounting",
            })?;
    let state_row_bytes = columns
        .checked_mul(2)
        .and_then(|values| values.checked_mul(size_of::<[f64; STATE_COUNT]>()))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "MAC Forward/Backward row workspace accounting",
        })?;
    let accuracy_bytes =
        columns
            .checked_mul(2 * size_of::<f64>())
            .ok_or(ClustalOmegaError::ArithmeticOverflow {
                stage: "MAC accuracy workspace accounting",
            })?;
    let transition_bytes = left_length
        .checked_add(right_length)
        .and_then(|positions| positions.checked_add(2))
        .and_then(|positions| positions.checked_mul(size_of::<ProbabilityTransitionWeights>()))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "MAC transition workspace accounting",
        })?;
    forward_match_bytes
        .checked_add(trace_bytes)
        .and_then(|bytes| bytes.checked_add(emission_bytes))
        .and_then(|bytes| bytes.checked_add(scale_bytes))
        .and_then(|bytes| bytes.checked_add(state_row_bytes))
        .and_then(|bytes| bytes.checked_add(accuracy_bytes))
        .and_then(|bytes| bytes.checked_add(transition_bytes))
        .ok_or(ClustalOmegaError::ArithmeticOverflow {
            stage: "MAC workspace accounting",
        })
}

fn checked_cells(
    rows: usize,
    columns: usize,
    stage: &'static str,
) -> Result<usize, ClustalOmegaError> {
    rows.checked_mul(columns)
        .ok_or(ClustalOmegaError::ArithmeticOverflow { stage })
}

const fn state_index(row: usize, column: usize, columns: usize, state: usize) -> usize {
    (row * columns + column) * STATE_COUNT + state
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence};
    use align_msa::SequenceRecord;

    use super::{
        ColumnOperation, DG, GD, IM, MI, MM, OmegaEngine, ProfileHMM, WorkingProfile, WorkingRow,
        align_profiles, background, build_profile_hmm, close_to_match, column_score, extend_dg,
        extend_gd, extend_im, extend_mi, mac_align_logspace, mac_align_scaled, open_dg, open_gd,
        open_im, open_mi, profile_transitions, viterbi_align,
    };
    use crate::{ClustalOmegaParameters, OmegaEngineMode};

    fn record(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(
            identifier,
            DnaSequence::from_ascii(sequence.as_bytes()).unwrap(),
        )
        .unwrap()
    }

    fn binary_sequences(maximum_length: usize) -> Vec<String> {
        (1..=maximum_length)
            .flat_map(|length| {
                (0..1_usize << length).map(move |bits| {
                    (0..length)
                        .map(|position| {
                            if bits & (1 << position) == 0 {
                                'A'
                            } else {
                                'C'
                            }
                        })
                        .collect()
                })
            })
            .collect()
    }

    #[derive(Clone, Debug)]
    struct EnumeratedPath {
        score: f64,
        left_start: usize,
        right_start: usize,
        left_end: usize,
        right_end: usize,
        operations: Vec<ColumnOperation>,
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn enumerate_state_paths(
        left: &ProfileHMM,
        right: &ProfileHMM,
        parameters: &ClustalOmegaParameters,
        state: usize,
        row: usize,
        column: usize,
        left_start: usize,
        right_start: usize,
        score: f64,
        operations: &mut Vec<ColumnOperation>,
        paths: &mut Vec<EnumeratedPath>,
    ) {
        let endpoint = state == MM
            && ((row == left.len() && column > 0) || (column == right.len() && row > 0));
        if endpoint {
            paths.push(EnumeratedPath {
                score,
                left_start,
                right_start,
                left_end: row,
                right_end: column,
                operations: operations.clone(),
            });
            return;
        }

        if row < left.len() && column < right.len() {
            let next_row = row + 1;
            let next_column = column + 1;
            let emission = column_score(
                left,
                right,
                next_row,
                next_column,
                &background(left.alphabet_size == 20),
                parameters,
            )
            .unwrap();
            operations.push(ColumnOperation::Both);
            enumerate_state_paths(
                left,
                right,
                parameters,
                MM,
                next_row,
                next_column,
                left_start,
                right_start,
                score + close_to_match(left, right, next_row, next_column, state) + emission,
                operations,
                paths,
            );
            operations.pop();
        }

        let factors = parameters.gap_path_factors();
        if row < left.len() && column > 0 {
            let next_row = row + 1;
            for (next_state, transition) in match state {
                MM => vec![
                    (MI, factors.2 * open_mi(left, right, next_row, column)),
                    (DG, factors.0 * open_dg(left, next_row)),
                ],
                MI => vec![(MI, factors.3 * extend_mi(left, right, next_row, column))],
                DG => vec![(DG, factors.1 * extend_dg(left, next_row))],
                IM | GD => Vec::new(),
                _ => unreachable!("five states are exhaustive"),
            } {
                operations.push(ColumnOperation::Left);
                enumerate_state_paths(
                    left,
                    right,
                    parameters,
                    next_state,
                    next_row,
                    column,
                    left_start,
                    right_start,
                    score + transition,
                    operations,
                    paths,
                );
                operations.pop();
            }
        }
        if column < right.len() && row > 0 {
            let next_column = column + 1;
            for (next_state, transition) in match state {
                MM => vec![
                    (IM, factors.2 * open_im(left, right, row, next_column)),
                    (GD, factors.0 * open_gd(right, next_column)),
                ],
                IM => vec![(IM, factors.3 * extend_im(left, right, row, next_column))],
                GD => vec![(GD, factors.1 * extend_gd(right, next_column))],
                MI | DG => Vec::new(),
                _ => unreachable!("five states are exhaustive"),
            } {
                operations.push(ColumnOperation::Right);
                enumerate_state_paths(
                    left,
                    right,
                    parameters,
                    next_state,
                    row,
                    next_column,
                    left_start,
                    right_start,
                    score + transition,
                    operations,
                    paths,
                );
                operations.pop();
            }
        }
    }

    fn enumerate_all_paths(
        left: &ProfileHMM,
        right: &ProfileHMM,
        parameters: &ClustalOmegaParameters,
    ) -> Vec<EnumeratedPath> {
        let mut paths = Vec::new();
        for row in 0..=left.len() {
            enumerate_state_paths(
                left,
                right,
                parameters,
                MM,
                row,
                0,
                row,
                0,
                0.0,
                &mut Vec::new(),
                &mut paths,
            );
        }
        for column in 1..=right.len() {
            enumerate_state_paths(
                left,
                right,
                parameters,
                MM,
                0,
                column,
                0,
                column,
                0.0,
                &mut Vec::new(),
                &mut paths,
            );
        }
        paths
    }

    #[test]
    fn single_sequence_hmm_probabilities_are_normalized() {
        let record = record("a", "ACGT");
        let profile = WorkingProfile::single(&record, 0);
        let hmm = build_profile_hmm(
            &profile,
            &ClustalOmegaParameters::nucleotide_v1_2_4(),
            false,
        )
        .unwrap();
        assert_eq!(hmm.len(), 4);
        assert!(hmm.emissions.chunks(hmm.alphabet_size).all(|emission| {
            (emission.iter().sum::<f64>() - 1.0).abs() < 1.0e-12
                && emission.iter().all(|value| *value > 0.0)
        }));
        assert_eq!(hmm.sequence_weights, vec![1.0]);
        assert!(hmm.effective_diversity.iter().all(|value| *value >= 1.0));
    }

    #[test]
    fn terminal_gaps_do_not_create_profile_delete_transitions() {
        let full = record("full", "ACG");
        let symbols = full.sequence().as_slice();
        let terminal_profile = WorkingProfile {
            rows: vec![
                WorkingRow {
                    sequence_index: 0,
                    cells: symbols.iter().copied().map(Some).collect(),
                },
                WorkingRow {
                    sequence_index: 1,
                    cells: vec![None, Some(symbols[1]), None],
                },
            ],
        };
        let empty_companion = WorkingProfile {
            rows: vec![
                WorkingRow {
                    sequence_index: 0,
                    cells: symbols.iter().copied().map(Some).collect(),
                },
                WorkingRow {
                    sequence_index: 1,
                    cells: vec![None, None, None],
                },
            ],
        };
        let parameters = ClustalOmegaParameters::nucleotide_v1_2_4();
        let terminal = profile_transitions(&terminal_profile, &[1.0, 1.0], &parameters).unwrap();
        let empty = profile_transitions(&empty_companion, &[1.0, 1.0], &parameters).unwrap();
        assert_eq!(terminal[0], empty[0]);
        assert_ne!(terminal[1], empty[1]);
        assert_eq!(terminal[2], empty[2]);
        assert_eq!(terminal[3], empty[3]);
    }

    #[test]
    fn nucleotide_columns_use_the_published_identity_matrix() {
        let parameters = ClustalOmegaParameters::nucleotide_v1_2_4()
            .with_emission_pseudocounts(1.0, 1.5, 1.0)
            .unwrap();
        let adenine = build_profile_hmm(
            &WorkingProfile::single(&record("a", "A"), 0),
            &parameters,
            false,
        )
        .unwrap();
        let cytosine = build_profile_hmm(
            &WorkingProfile::single(&record("c", "C"), 1),
            &parameters,
            false,
        )
        .unwrap();
        let nucleotide_background = background(false);
        let identity = column_score(
            &adenine,
            &adenine,
            1,
            1,
            &nucleotide_background,
            &parameters,
        )
        .unwrap();
        let mismatch = column_score(
            &adenine,
            &cytosine,
            1,
            1,
            &nucleotide_background,
            &parameters,
        )
        .unwrap();
        let normalizer = 2.0_f64.sqrt() + 3.0;
        let diagonal = 2.0_f64.sqrt() / normalizer;
        let off_diagonal = 1.0 / normalizer;
        let expected_identity =
            (4.0 * (diagonal.powi(2) + 3.0 * off_diagonal.powi(2))).log2() - 0.03;
        let expected_mismatch =
            (4.0 * (2.0 * diagonal * off_diagonal + 2.0 * off_diagonal.powi(2))).log2() - 0.03;
        assert!((identity - expected_identity).abs() < 1.0e-12);
        assert!((mismatch - expected_mismatch).abs() < 1.0e-12);
        assert!(identity > mismatch);
    }

    #[test]
    fn viterbi_score_and_traceback_match_exhaustive_state_path_enumeration() {
        let parameters =
            ClustalOmegaParameters::nucleotide_v1_2_4().with_engine_mode(OmegaEngineMode::Viterbi);
        for (left_text, right_text) in [("ACG", "AG"), ("AC", "CA"), ("AGT", "GT")] {
            let left_profile = WorkingProfile::single(&record("left", left_text), 0);
            let right_profile = WorkingProfile::single(&record("right", right_text), 1);
            let left = build_profile_hmm(&left_profile, &parameters, false).unwrap();
            let right = build_profile_hmm(&right_profile, &parameters, false).unwrap();
            let observed = viterbi_align(&left, &right, &parameters).unwrap();
            let paths = enumerate_all_paths(&left, &right, &parameters);
            let best_score = paths
                .iter()
                .map(|path| path.score)
                .max_by(f64::total_cmp)
                .unwrap();
            assert!((observed.score - best_score).abs() < 1.0e-12);
            assert!(paths.iter().any(|path| {
                (path.score - observed.score).abs() < 1.0e-12
                    && path.left_start == observed.left_start
                    && path.right_start == observed.right_start
                    && path.left_end == observed.left_end
                    && path.right_end == observed.right_end
                    && path.operations == observed.operations
            }));
        }
    }

    #[test]
    fn viterbi_traceback_matches_the_identity_model_terminal_choice() {
        let left_record = record("long", "ACGT");
        let right_record = record("short", "AGT");
        let parameters =
            ClustalOmegaParameters::nucleotide_v1_2_4().with_engine_mode(OmegaEngineMode::Viterbi);
        let result = align_profiles(
            WorkingProfile::single(&left_record, 0),
            WorkingProfile::single(&right_record, 1),
            &parameters,
            false,
        )
        .unwrap();
        assert_eq!(result.observation.engine, OmegaEngine::Viterbi);
        let rows = &result.profile.rows;
        let rendered = rows
            .iter()
            .map(|row| {
                row.cells
                    .iter()
                    .map(|cell| cell.map_or('-', |symbol| char::from(symbol.as_ascii())))
                    .collect::<String>()
            })
            .collect::<Vec<_>>();
        assert_eq!(rendered, vec!["ACGT", "-AGT"]);
    }

    #[test]
    fn scaled_mac_matches_the_independent_log_space_oracle() {
        let parameters =
            ClustalOmegaParameters::nucleotide_v1_2_4().with_engine_mode(OmegaEngineMode::Mac);
        let binary = binary_sequences(4);
        let cases = binary
            .iter()
            .flat_map(|left| {
                binary
                    .iter()
                    .map(move |right| (left.as_str(), right.as_str()))
            })
            .chain([
                ("ACGTAC", "ACGAC"),
                ("ACAC", "CAC"),
                ("AGT", "GT"),
                ("AAAA", "CCCC"),
            ]);
        for (left_text, right_text) in cases {
            let left_profile = WorkingProfile::single(&record("left", left_text), 0);
            let right_profile = WorkingProfile::single(&record("right", right_text), 1);
            let left = build_profile_hmm(&left_profile, &parameters, false).unwrap();
            let right = build_profile_hmm(&right_profile, &parameters, false).unwrap();
            for prefer_left_gap in [false, true] {
                let scaled = mac_align_scaled(&left, &right, &parameters, prefer_left_gap).unwrap();
                let log_space =
                    mac_align_logspace(&left, &right, &parameters, prefer_left_gap).unwrap();
                assert_eq!(scaled.left_start, log_space.left_start);
                assert_eq!(scaled.left_end, log_space.left_end);
                assert_eq!(scaled.right_start, log_space.right_start);
                assert_eq!(scaled.right_end, log_space.right_end);
                assert_eq!(scaled.operations, log_space.operations);
                assert!((scaled.score - log_space.score).abs() < 1.0e-10);
                let scaled_bounds = scaled.posterior_bounds.unwrap();
                let log_bounds = log_space.posterior_bounds.unwrap();
                assert!((scaled_bounds.0 - log_bounds.0).abs() < 1.0e-10);
                assert!((scaled_bounds.1 - log_bounds.1).abs() < 1.0e-10);
            }
        }
    }

    #[test]
    fn mac_posteriors_are_bounded_and_traceback_replays() {
        let first = record("a", "ACGTAC");
        let second = record("b", "ACGAC");
        let parameters =
            ClustalOmegaParameters::nucleotide_v1_2_4().with_engine_mode(OmegaEngineMode::Mac);
        let result = align_profiles(
            WorkingProfile::single(&first, 0),
            WorkingProfile::single(&second, 1),
            &parameters,
            false,
        )
        .unwrap();
        assert_eq!(result.observation.engine, OmegaEngine::Mac);
        let (minimum, maximum) = result.observation.posterior_bounds.unwrap();
        assert!((0.0..=1.0).contains(&minimum));
        assert!((0.0..=1.0).contains(&maximum));
        assert!(minimum <= maximum);
        assert_eq!(
            result.profile.rows[0].cells.len(),
            result.profile.rows[1].cells.len()
        );
    }

    #[test]
    fn empty_profile_rows_are_reattached_without_losing_symbols() {
        let empty = record("empty", "");
        let full = record("full", "ACG");
        let result = align_profiles(
            WorkingProfile::single(&empty, 0),
            WorkingProfile::single(&full, 1),
            &ClustalOmegaParameters::nucleotide_v1_2_4(),
            false,
        )
        .unwrap();
        assert_eq!(result.profile.rows[0].cells, vec![None, None, None]);
        assert_eq!(
            result.profile.rows[1]
                .cells
                .iter()
                .filter_map(|cell| *cell)
                .map(|symbol| char::from(symbol.as_ascii()))
                .collect::<String>(),
            "ACG"
        );
    }
}
