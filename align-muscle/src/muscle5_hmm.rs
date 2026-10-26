// Exact equality is intentional for probabilities constructed from one shared
// value; losing bit symmetry would violate sequence-swap invariance.
#![allow(clippy::float_cmp)]

use core::marker::PhantomData;

use align_core::Symbol;

use crate::Muscle5Error;
use crate::muscle5_emissions::{EmissionModel, Muscle5Alphabet};
use crate::muscle5_rng::{HMM_DOMAIN, SplitMix64, derived_seed};

pub(crate) const STATE_COUNT: usize = 5;
pub(crate) const MATCH_STATE: usize = 0;
pub(crate) const SHORT_X_STATE: usize = 1;
pub(crate) const SHORT_Y_STATE: usize = 2;
pub(crate) const LONG_X_STATE: usize = 3;
pub(crate) const LONG_Y_STATE: usize = 4;

const NORMALIZATION_TOLERANCE: f64 = 1.0e-12;

#[derive(Clone, Debug)]
pub(crate) struct HmmParameters {
    start: [f64; STATE_COUNT],
    end: [f64; STATE_COUNT],
    transition: [[f64; STATE_COUNT]; STATE_COUNT],
    emissions: EmissionModel,
}

impl HmmParameters {
    pub(crate) fn published<A: Muscle5Alphabet>(
        perturbation_seed: u64,
        amplitude: f64,
    ) -> Result<Self, Muscle5Error> {
        let mut start = [0.62, 0.02, 0.02, 0.18, 0.18];
        normalize_distribution(&mut start);
        let mut transition = [[0.0; STATE_COUNT]; STATE_COUNT];
        transition[MATCH_STATE] = [0.96, 0.012, 0.012, 0.008, 0.008];
        transition[SHORT_X_STATE][MATCH_STATE] = 0.60;
        transition[SHORT_X_STATE][SHORT_X_STATE] = 0.40;
        transition[SHORT_Y_STATE][MATCH_STATE] = 0.60;
        transition[SHORT_Y_STATE][SHORT_Y_STATE] = 0.40;
        transition[LONG_X_STATE][MATCH_STATE] = 0.10;
        transition[LONG_X_STATE][LONG_X_STATE] = 0.90;
        transition[LONG_Y_STATE][MATCH_STATE] = 0.10;
        transition[LONG_Y_STATE][LONG_Y_STATE] = 0.90;

        let mut parameters = Self {
            start,
            end: start,
            transition,
            emissions: EmissionModel::for_alphabet::<A>(),
        };
        if perturbation_seed != 0 {
            parameters.perturb(perturbation_seed, amplitude);
        }
        parameters.validate()?;
        Ok(parameters)
    }

    pub(crate) fn prepare<A: Muscle5Alphabet>(&self) -> PreparedHmm<A> {
        let alphabet_size = A::SYMBOLS.len();
        let mut match_probability = Vec::with_capacity(alphabet_size * alphabet_size);
        let mut insert_probability = Vec::with_capacity(alphabet_size);
        for symbol in A::SYMBOLS {
            let symbol = Symbol::<A>::from_ascii(*symbol).expect("alphabet constant is valid");
            insert_probability.push(self.emissions.insert_probability(symbol));
        }
        for left in A::SYMBOLS {
            let left = Symbol::<A>::from_ascii(*left).expect("alphabet constant is valid");
            for right in A::SYMBOLS {
                let right = Symbol::<A>::from_ascii(*right).expect("alphabet constant is valid");
                match_probability.push(self.emissions.match_probability(left, right));
            }
        }
        PreparedHmm {
            start_log: self.start.map(f64::ln),
            end_log: self.end.map(f64::ln),
            transition_log: self.transition.map(|row| row.map(log_probability)),
            start_probability: self.start,
            end_probability: self.end,
            transition_probability: self.transition,
            match_log: match_probability.iter().copied().map(f64::ln).collect(),
            insert_log: insert_probability.iter().copied().map(f64::ln).collect(),
            match_probability: match_probability.into_boxed_slice(),
            insert_probability: insert_probability.into_boxed_slice(),
            alphabet_size,
            marker: PhantomData,
        }
    }

    fn perturb(&mut self, seed: u64, amplitude: f64) {
        let mut rng = SplitMix64::new(derived_seed(seed, HMM_DOMAIN));
        self.start = perturb_symmetric_five(self.start, amplitude, &mut rng);
        self.end = self.start;
        self.transition[MATCH_STATE] =
            perturb_symmetric_five(self.transition[MATCH_STATE], amplitude, &mut rng);

        let short = perturb_pair(
            [
                self.transition[SHORT_X_STATE][MATCH_STATE],
                self.transition[SHORT_X_STATE][SHORT_X_STATE],
            ],
            amplitude,
            &mut rng,
        );
        for state in [SHORT_X_STATE, SHORT_Y_STATE] {
            self.transition[state].fill(0.0);
            self.transition[state][MATCH_STATE] = short[0];
            self.transition[state][state] = short[1];
        }

        let long = perturb_pair(
            [
                self.transition[LONG_X_STATE][MATCH_STATE],
                self.transition[LONG_X_STATE][LONG_X_STATE],
            ],
            amplitude,
            &mut rng,
        );
        for state in [LONG_X_STATE, LONG_Y_STATE] {
            self.transition[state].fill(0.0);
            self.transition[state][MATCH_STATE] = long[0];
            self.transition[state][state] = long[1];
        }

        let canonical_count = self.emissions.canonical_count();
        let mut normalization = 0.0;
        for row in 0..canonical_count {
            for column in row..canonical_count {
                let factor = 1.0 + amplitude * rng.signed_unit_f64();
                let value = self.emissions.joint_mut()[row][column] * factor;
                self.emissions.joint_mut()[row][column] = value;
                self.emissions.joint_mut()[column][row] = value;
                normalization += if row == column { value } else { 2.0 * value };
            }
        }
        for row in 0..canonical_count {
            for value in &mut self.emissions.joint_mut()[row][..canonical_count] {
                *value /= normalization;
            }
        }
        self.emissions.recompute_marginal();
    }

    fn validate(&self) -> Result<(), Muscle5Error> {
        validate_distribution(&self.start, "start transition normalization")?;
        validate_distribution(&self.end, "end transition normalization")?;
        for row in &self.transition {
            validate_distribution(row, "state transition normalization")?;
        }
        if self.start[SHORT_X_STATE] != self.start[SHORT_Y_STATE]
            || self.start[LONG_X_STATE] != self.start[LONG_Y_STATE]
            || self.end != self.start
            || self.transition[SHORT_X_STATE][MATCH_STATE]
                != self.transition[SHORT_Y_STATE][MATCH_STATE]
            || self.transition[SHORT_X_STATE][SHORT_X_STATE]
                != self.transition[SHORT_Y_STATE][SHORT_Y_STATE]
            || self.transition[LONG_X_STATE][MATCH_STATE]
                != self.transition[LONG_Y_STATE][MATCH_STATE]
            || self.transition[LONG_X_STATE][LONG_X_STATE]
                != self.transition[LONG_Y_STATE][LONG_Y_STATE]
        {
            return Err(Muscle5Error::InternalInvariant {
                message: "perturbed HMM lost sequence-swap symmetry",
            });
        }
        if !self.emissions.validate() {
            return Err(Muscle5Error::InternalInvariant {
                message: "HMM emissions are not positive normalized probabilities",
            });
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PreparedHmm<A: Muscle5Alphabet> {
    pub(crate) start_log: [f64; STATE_COUNT],
    pub(crate) end_log: [f64; STATE_COUNT],
    pub(crate) transition_log: [[f64; STATE_COUNT]; STATE_COUNT],
    pub(crate) start_probability: [f64; STATE_COUNT],
    pub(crate) end_probability: [f64; STATE_COUNT],
    pub(crate) transition_probability: [[f64; STATE_COUNT]; STATE_COUNT],
    match_log: Box<[f64]>,
    insert_log: Box<[f64]>,
    match_probability: Box<[f64]>,
    insert_probability: Box<[f64]>,
    alphabet_size: usize,
    marker: PhantomData<A>,
}

impl<A: Muscle5Alphabet> PreparedHmm<A> {
    #[inline]
    pub(crate) fn match_log(&self, left: Symbol<A>, right: Symbol<A>) -> f64 {
        self.match_log[left.alphabet_index() * self.alphabet_size + right.alphabet_index()]
    }

    #[inline]
    pub(crate) fn insert_log(&self, symbol: Symbol<A>) -> f64 {
        self.insert_log[symbol.alphabet_index()]
    }

    #[inline]
    pub(crate) fn match_probability(&self, left: Symbol<A>, right: Symbol<A>) -> f64 {
        self.match_probability[left.alphabet_index() * self.alphabet_size + right.alphabet_index()]
    }

    #[inline]
    pub(crate) fn insert_probability(&self, symbol: Symbol<A>) -> f64 {
        self.insert_probability[symbol.alphabet_index()]
    }

    pub(crate) fn storage_bytes(&self) -> usize {
        (self.match_log.len()
            + self.insert_log.len()
            + self.match_probability.len()
            + self.insert_probability.len()
            + 2 * (STATE_COUNT * STATE_COUNT + 2 * STATE_COUNT))
            * size_of::<f64>()
    }
}

fn perturb_symmetric_five(
    values: [f64; STATE_COUNT],
    amplitude: f64,
    rng: &mut SplitMix64,
) -> [f64; STATE_COUNT] {
    let center = values[MATCH_STATE] * (1.0 + amplitude * rng.signed_unit_f64());
    let short = values[SHORT_X_STATE] * (1.0 + amplitude * rng.signed_unit_f64());
    let long = values[LONG_X_STATE] * (1.0 + amplitude * rng.signed_unit_f64());
    let normalization = center + 2.0 * short + 2.0 * long;
    [
        center / normalization,
        short / normalization,
        short / normalization,
        long / normalization,
        long / normalization,
    ]
}

fn perturb_pair(values: [f64; 2], amplitude: f64, rng: &mut SplitMix64) -> [f64; 2] {
    let first = values[0] * (1.0 + amplitude * rng.signed_unit_f64());
    let second = values[1] * (1.0 + amplitude * rng.signed_unit_f64());
    let normalization = first + second;
    [first / normalization, second / normalization]
}

fn validate_distribution(values: &[f64], stage: &'static str) -> Result<(), Muscle5Error> {
    if values
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
        || (values.iter().sum::<f64>() - 1.0).abs() > NORMALIZATION_TOLERANCE
    {
        return Err(Muscle5Error::NonFiniteComputation { stage });
    }
    Ok(())
}

fn normalize_distribution(values: &mut [f64]) {
    let normalization = values.iter().sum::<f64>();
    for value in values {
        *value /= normalization;
    }
}

fn log_probability(probability: f64) -> f64 {
    if probability == 0.0 {
        f64::NEG_INFINITY
    } else {
        probability.ln()
    }
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence};

    use super::{
        HmmParameters, LONG_X_STATE, LONG_Y_STATE, MATCH_STATE, SHORT_X_STATE, SHORT_Y_STATE,
    };

    #[test]
    fn published_transition_sets_match_the_nature_profile() {
        let hmm = HmmParameters::published::<Protein>(0, 0.25).unwrap();
        assert!((hmm.start.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
        assert!((hmm.start[SHORT_X_STATE] / hmm.start[MATCH_STATE] - 0.02 / 0.62).abs() < 1.0e-15);
        assert!((hmm.start[LONG_X_STATE] / hmm.start[MATCH_STATE] - 0.18 / 0.62).abs() < 1.0e-15);
        assert_eq!(hmm.end, hmm.start);
        assert_eq!(
            hmm.transition[MATCH_STATE],
            [0.96, 0.012, 0.012, 0.008, 0.008]
        );
        assert_eq!(hmm.transition[SHORT_X_STATE][MATCH_STATE], 0.60);
        assert_eq!(hmm.transition[SHORT_X_STATE][SHORT_X_STATE], 0.40);
        assert_eq!(hmm.transition[LONG_X_STATE][MATCH_STATE], 0.10);
        assert_eq!(hmm.transition[LONG_X_STATE][LONG_X_STATE], 0.90);
    }

    #[test]
    fn perturbation_is_deterministic_normalized_and_symmetric() {
        let first = HmmParameters::published::<Dna>(7, 0.25).unwrap();
        let second = HmmParameters::published::<Dna>(7, 0.25).unwrap();
        let other = HmmParameters::published::<Dna>(8, 0.25).unwrap();
        assert_eq!(first.start, second.start);
        assert_eq!(first.transition, second.transition);
        assert_ne!(first.start, other.start);
        assert_eq!(first.start[SHORT_X_STATE], first.start[SHORT_Y_STATE]);
        assert_eq!(first.start[LONG_X_STATE], first.start[LONG_Y_STATE]);
        assert!((first.start.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
        for row in first.transition {
            assert!((row.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
        }
    }

    #[test]
    fn prepared_emissions_preserve_swap_symmetry() {
        let prepared = HmmParameters::published::<Protein>(3, 0.25)
            .unwrap()
            .prepare::<Protein>();
        let symbols = ProteinSequence::from_ascii("AW").unwrap();
        let a = symbols.get(0).unwrap();
        let w = symbols.get(1).unwrap();
        assert_eq!(prepared.match_log(a, w), prepared.match_log(w, a));
        assert!(prepared.insert_log(a).is_finite());

        let dna = HmmParameters::published::<Dna>(0, 0.25)
            .unwrap()
            .prepare::<Dna>();
        let a = DnaSequence::from_ascii("A").unwrap().get(0).unwrap();
        assert!(dna.match_log(a, a) > dna.insert_log(a).mul_add(2.0, 0.0));
        assert!(dna.storage_bytes() > 0);
    }
}
