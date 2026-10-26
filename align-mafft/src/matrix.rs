use core::any::TypeId;

use align_core::{Dna, Protein, Rna, Symbol};

use crate::{FftNs1Parameters, MafftAlphabet, MafftError};

const PROTEIN_CANONICAL: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

// Tables II and III of Jones, Taylor, and Thornton (1992). The printed PAM1
// probabilities are scaled by 100,000 and rounded. Rows are source residues;
// columns are target residues in ARNDCQEGHILKMFPSTWYV order.
const JTT_PAM1_SCALED: [[u32; 20]; 20] = [
    [
        98_759, 27, 24, 42, 12, 23, 66, 129, 5, 19, 28, 22, 11, 6, 99, 264, 267, 1, 4, 193,
    ],
    [
        41, 98_962, 19, 8, 21, 125, 20, 102, 74, 13, 34, 390, 10, 3, 36, 69, 38, 18, 8, 11,
    ],
    [
        43, 23, 98_707, 284, 6, 31, 36, 58, 92, 26, 12, 150, 8, 3, 6, 344, 137, 0, 23, 11,
    ],
    [
        63, 8, 235, 98_932, 2, 21, 478, 95, 24, 6, 6, 17, 4, 1, 6, 40, 25, 1, 15, 21,
    ],
    [
        44, 52, 13, 5, 99_450, 4, 3, 41, 17, 8, 15, 3, 10, 28, 6, 147, 28, 16, 68, 41,
    ],
    [
        43, 154, 33, 27, 2, 98_955, 211, 17, 130, 4, 64, 176, 11, 2, 81, 37, 31, 2, 8, 12,
    ],
    [
        82, 16, 25, 398, 1, 140, 99_042, 83, 6, 6, 9, 102, 4, 2, 10, 21, 19, 2, 2, 31,
    ],
    [
        135, 70, 33, 66, 11, 10, 70, 99_369, 5, 3, 6, 16, 3, 2, 11, 129, 19, 8, 2, 32,
    ],
    [
        17, 164, 171, 53, 15, 223, 15, 15, 98_867, 10, 49, 31, 8, 18, 58, 51, 28, 2, 189, 8,
    ],
    [
        28, 12, 21, 6, 3, 3, 7, 4, 4, 98_722, 212, 12, 113, 31, 5, 28, 149, 2, 10, 630,
    ],
    [
        24, 19, 6, 3, 3, 29, 6, 5, 12, 122, 99_328, 9, 90, 101, 53, 40, 16, 8, 8, 117,
    ],
    [
        28, 334, 108, 14, 1, 122, 107, 20, 12, 11, 13, 99_101, 15, 1, 11, 32, 57, 1, 3, 8,
    ],
    [
        36, 22, 14, 10, 8, 19, 11, 10, 8, 253, 350, 37, 98_845, 18, 8, 19, 123, 3, 6, 201,
    ],
    [
        11, 3, 3, 2, 14, 2, 3, 4, 11, 41, 230, 1, 10, 99_357, 8, 65, 8, 8, 179, 40,
    ],
    [
        150, 36, 5, 7, 3, 66, 12, 16, 26, 5, 97, 13, 4, 6, 99_278, 190, 69, 1, 4, 14,
    ],
    [
        297, 51, 214, 30, 44, 22, 19, 139, 17, 21, 54, 28, 7, 38, 140, 98_548, 278, 4, 20, 27,
    ],
    [
        351, 33, 100, 22, 9, 21, 20, 24, 11, 134, 25, 57, 49, 6, 59, 325, 98_670, 1, 6, 76,
    ],
    [
        7, 65, 1, 3, 23, 7, 7, 41, 3, 7, 49, 5, 5, 22, 4, 21, 5, 99_684, 24, 16,
    ],
    [
        11, 12, 30, 23, 43, 10, 4, 4, 134, 16, 22, 5, 4, 222, 6, 43, 12, 11, 99_377, 11,
    ],
    [
        226, 9, 7, 16, 13, 7, 29, 35, 3, 504, 161, 7, 71, 24, 11, 28, 67, 3, 5, 98_772,
    ],
];

const JTT_FREQUENCIES_PRINTED: [f64; 20] = [
    0.077, 0.051, 0.043, 0.052, 0.020, 0.041, 0.062, 0.074, 0.023, 0.053, 0.091, 0.059, 0.024,
    0.040, 0.051, 0.069, 0.059, 0.014, 0.032, 0.066,
];

// Henikoff and Henikoff (1992) BLOSUM62 in ARNDCQEGHILKMFPSTWYV order.
const BLOSUM62: [[i8; 20]; 20] = [
    [
        4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0,
    ],
    [
        -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3,
    ],
    [
        -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3,
    ],
    [
        -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3,
    ],
    [
        0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
    ],
    [
        -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2,
    ],
    [
        -1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2,
    ],
    [
        0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3,
    ],
    [
        -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3,
    ],
    [
        -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3,
    ],
    [
        -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1,
    ],
    [
        -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2,
    ],
    [
        -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1,
    ],
    [
        -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1,
    ],
    [
        -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2,
    ],
    [
        1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2,
    ],
    [
        0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0,
    ],
    [
        -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3,
    ],
    [
        -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1,
    ],
    [
        0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4,
    ],
];

const BLOSUM62_FREQUENCIES_PRINTED: [f64; 20] = [
    0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074, 0.026, 0.068, 0.099, 0.058, 0.025,
    0.047, 0.039, 0.057, 0.051, 0.013, 0.032, 0.073,
];

// Grantham (1974), Table 1, in JTT residue order.
const GRANTHAM_VOLUME: [f64; 20] = [
    31.0, 124.0, 56.0, 54.0, 55.0, 85.0, 83.0, 3.0, 96.0, 111.0, 111.0, 119.0, 105.0, 132.0, 32.5,
    32.0, 61.0, 170.0, 136.0, 84.0,
];
const GRANTHAM_POLARITY: [f64; 20] = [
    8.1, 10.5, 11.6, 13.0, 5.5, 10.5, 12.3, 9.0, 10.4, 5.2, 4.9, 11.3, 5.7, 5.2, 8.0, 9.2, 8.6,
    5.4, 6.2, 5.9,
];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MoleculeKind {
    Protein,
    Nucleotide,
}

#[derive(Clone, Debug)]
pub(crate) struct ScientificScorer {
    kind: MoleculeKind,
    alphabet_size: usize,
    scores: Box<[f64]>,
    signals: Box<[[f64; 4]]>,
    #[cfg(test)]
    canonical_scores: Box<[f64]>,
    #[cfg(test)]
    canonical_frequencies: Box<[f64]>,
}

impl ScientificScorer {
    pub(crate) fn new<A: MafftAlphabet>(parameters: FftNs1Parameters) -> Result<Self, MafftError> {
        let kind = molecule_kind::<A>();
        let (canonical_scores, canonical_frequencies) = match kind {
            MoleculeKind::Protein => protein_matrix(parameters.random_match_offset())?,
            MoleculeKind::Nucleotide => nucleotide_matrix(parameters.random_match_offset())?,
        };
        Self::from_canonical::<A>(kind, canonical_scores, canonical_frequencies)
    }

    pub(crate) fn ginsi_2005<A: MafftAlphabet>(offset: f64) -> Result<Self, MafftError> {
        let kind = molecule_kind::<A>();
        let (canonical_scores, canonical_frequencies) = match kind {
            MoleculeKind::Protein => blosum62_matrix(offset)?,
            MoleculeKind::Nucleotide => nucleotide_matrix(offset)?,
        };
        Self::from_canonical::<A>(kind, canonical_scores, canonical_frequencies)
    }

    // Test builds retain both owned canonical vectors for equation-level audits.
    #[allow(clippy::needless_pass_by_value)]
    fn from_canonical<A: MafftAlphabet>(
        kind: MoleculeKind,
        canonical_scores: Vec<f64>,
        canonical_frequencies: Vec<f64>,
    ) -> Result<Self, MafftError> {
        let frequency_sum: f64 = canonical_frequencies.iter().sum();
        if !frequency_sum.is_finite() || (frequency_sum - 1.0).abs() > 1.0e-12 {
            return Err(MafftError::NonFiniteComputation {
                stage: "scientific background-frequency normalization",
            });
        }
        let canonical_count = match kind {
            MoleculeKind::Protein => 20,
            MoleculeKind::Nucleotide => 4,
        };
        let alphabet_size = A::SYMBOLS.len();
        let mut scores = vec![0.0; alphabet_size * alphabet_size];
        let mut first_weights = vec![0.0; canonical_count];
        let mut second_weights = vec![0.0; canonical_count];
        for (first, first_ascii) in A::SYMBOLS.iter().copied().enumerate() {
            canonical_weights(kind, first_ascii, &mut first_weights);
            for (second, second_ascii) in A::SYMBOLS.iter().copied().enumerate() {
                canonical_weights(kind, second_ascii, &mut second_weights);
                let mut score = 0.0;
                for a in 0..canonical_count {
                    for b in 0..canonical_count {
                        score += first_weights[a]
                            * second_weights[b]
                            * canonical_scores[a * canonical_count + b];
                    }
                }
                scores[first * alphabet_size + second] = score;
            }
        }

        let protein_volume = standardized(GRANTHAM_VOLUME);
        let protein_polarity = standardized(GRANTHAM_POLARITY);
        let signals = A::SYMBOLS
            .iter()
            .copied()
            .map(|ascii| signal_for_ascii(kind, ascii, &protein_volume, &protein_polarity))
            .collect::<Vec<_>>();
        if scores.iter().any(|score| !score.is_finite())
            || signals.iter().flatten().any(|value| !value.is_finite())
        {
            return Err(MafftError::NonFiniteComputation {
                stage: "scientific scoring-matrix construction",
            });
        }
        Ok(Self {
            kind,
            alphabet_size,
            scores: scores.into_boxed_slice(),
            signals: signals.into_boxed_slice(),
            #[cfg(test)]
            canonical_scores: canonical_scores.into_boxed_slice(),
            #[cfg(test)]
            canonical_frequencies: canonical_frequencies.into_boxed_slice(),
        })
    }

    pub(crate) const fn alphabet_size(&self) -> usize {
        self.alphabet_size
    }

    #[cfg(test)]
    pub(crate) fn score<A: MafftAlphabet>(&self, first: Symbol<A>, second: Symbol<A>) -> f64 {
        self.scores[first.alphabet_index() * self.alphabet_size + second.alphabet_index()]
    }

    pub(crate) fn score_indices(&self, first: usize, second: usize) -> f64 {
        self.scores[first * self.alphabet_size + second]
    }

    pub(crate) fn score_row(&self, first: usize) -> &[f64] {
        let start = first * self.alphabet_size;
        &self.scores[start..start + self.alphabet_size]
    }

    pub(crate) fn signal<A: MafftAlphabet>(&self, symbol: Symbol<A>) -> [f64; 4] {
        self.signals[symbol.alphabet_index()]
    }

    pub(crate) const fn channel_count(&self) -> usize {
        match self.kind {
            MoleculeKind::Protein => 2,
            MoleculeKind::Nucleotide => 4,
        }
    }

    #[cfg(test)]
    pub(crate) fn canonical_moments(&self) -> (f64, f64) {
        let count = self.canonical_frequencies.len();
        let mut identical = 0.0;
        let mut independent = 0.0;
        for first in 0..count {
            identical +=
                self.canonical_frequencies[first] * self.canonical_scores[first * count + first];
            for second in 0..count {
                independent += self.canonical_frequencies[first]
                    * self.canonical_frequencies[second]
                    * self.canonical_scores[first * count + second];
            }
        }
        (identical, independent)
    }
}

pub(crate) fn molecule_kind<A: MafftAlphabet>() -> MoleculeKind {
    if TypeId::of::<A>() == TypeId::of::<Protein>() {
        MoleculeKind::Protein
    } else {
        debug_assert!(
            TypeId::of::<A>() == TypeId::of::<Dna>() || TypeId::of::<A>() == TypeId::of::<Rna>()
        );
        MoleculeKind::Nucleotide
    }
}

pub(crate) fn canonical_weights(kind: MoleculeKind, ascii: u8, output: &mut [f64]) {
    output.fill(0.0);
    match kind {
        MoleculeKind::Protein => protein_weights(ascii, output),
        MoleculeKind::Nucleotide => nucleotide_weights(ascii, output),
    }
}

fn signal_for_ascii(
    kind: MoleculeKind,
    ascii: u8,
    protein_volume: &[f64; 20],
    protein_polarity: &[f64; 20],
) -> [f64; 4] {
    let mut weights = [0.0; 20];
    canonical_weights(kind, ascii, &mut weights);
    match kind {
        MoleculeKind::Protein => {
            let mut signal = [0.0; 4];
            for (index, weight) in weights.iter().copied().enumerate().take(20) {
                signal[0] += weight * protein_volume[index];
                signal[1] += weight * protein_polarity[index];
            }
            signal
        }
        MoleculeKind::Nucleotide => [weights[0], weights[1], weights[2], weights[3]],
    }
}

fn protein_weights(ascii: u8, output: &mut [f64]) {
    if let Some(index) = PROTEIN_CANONICAL.iter().position(|value| *value == ascii) {
        output[index] = 1.0;
        return;
    }
    match ascii {
        b'B' => set_equal(output, &[2, 3]),
        b'Z' => set_equal(output, &[5, 6]),
        b'X' => output[..20].fill(0.05),
        b'J' => set_equal(output, &[9, 10]),
        b'U' => output[4] = 1.0,
        b'O' => output[11] = 1.0,
        _ => unreachable!("typed protein symbol must be known"),
    }
}

fn nucleotide_weights(ascii: u8, output: &mut [f64]) {
    match ascii {
        b'A' => output[0] = 1.0,
        b'C' => output[1] = 1.0,
        b'G' => output[2] = 1.0,
        b'T' | b'U' => output[3] = 1.0,
        b'R' => set_equal(output, &[0, 2]),
        b'Y' => set_equal(output, &[1, 3]),
        b'S' => set_equal(output, &[1, 2]),
        b'W' => set_equal(output, &[0, 3]),
        b'K' => set_equal(output, &[2, 3]),
        b'M' => set_equal(output, &[0, 1]),
        b'B' => set_equal(output, &[1, 2, 3]),
        b'D' => set_equal(output, &[0, 2, 3]),
        b'H' => set_equal(output, &[0, 1, 3]),
        b'V' => set_equal(output, &[0, 1, 2]),
        b'N' => set_equal(output, &[0, 1, 2, 3]),
        _ => unreachable!("typed nucleotide symbol must be known"),
    }
}

fn set_equal(output: &mut [f64], indices: &[usize]) {
    #[allow(clippy::cast_precision_loss)]
    let weight = 1.0 / indices.len() as f64;
    for index in indices {
        output[*index] = weight;
    }
}

fn blosum62_matrix(offset: f64) -> Result<(Vec<f64>, Vec<f64>), MafftError> {
    let frequency_sum: f64 = BLOSUM62_FREQUENCIES_PRINTED.iter().sum();
    let frequencies = BLOSUM62_FREQUENCIES_PRINTED
        .iter()
        .map(|frequency| frequency / frequency_sum)
        .collect::<Vec<_>>();
    let raw = BLOSUM62
        .iter()
        .flat_map(|row| row.iter())
        .map(|score| f64::from(*score))
        .collect::<Vec<_>>();
    normalize_matrix(raw, frequencies, offset)
}

fn protein_matrix(offset: f64) -> Result<(Vec<f64>, Vec<f64>), MafftError> {
    let mut transition = [[0.0; 20]; 20];
    for row in 0..20 {
        let row_sum: u32 = JTT_PAM1_SCALED[row].iter().sum();
        for column in 0..20 {
            transition[row][column] = f64::from(JTT_PAM1_SCALED[row][column]) / f64::from(row_sum);
        }
    }
    let transition_200 = matrix_power_20(transition, 200);
    let frequency_sum: f64 = JTT_FREQUENCIES_PRINTED.iter().sum();
    let frequencies: Vec<f64> = JTT_FREQUENCIES_PRINTED
        .iter()
        .map(|frequency| frequency / frequency_sum)
        .collect();
    let mut raw = vec![0.0; 400];
    for first in 0..20 {
        for second in 0..20 {
            let forward = (transition_200[first][second] / frequencies[second]).ln();
            let reverse = (transition_200[second][first] / frequencies[first]).ln();
            raw[first * 20 + second] = (forward + reverse) * 0.5;
        }
    }
    normalize_matrix(raw, frequencies, offset)
}

fn nucleotide_matrix(offset: f64) -> Result<(Vec<f64>, Vec<f64>), MafftError> {
    // alpha / beta = 2 and (alpha + 2 beta)t = 2 imply beta*t = 0.5.
    let exp_four_beta = (-2.0_f64).exp();
    let exp_alpha_beta = (-3.0_f64).exp();
    let same = 0.25 + 0.25 * exp_four_beta + 0.5 * exp_alpha_beta;
    let transition = 0.25 + 0.25 * exp_four_beta - 0.5 * exp_alpha_beta;
    let transversion = 0.25 - 0.25 * exp_four_beta;
    let mut raw = vec![0.0; 16];
    for first in 0..4 {
        for second in 0..4 {
            let probability = if first == second {
                same
            } else if matches!((first, second), (0, 2) | (2, 0) | (1, 3) | (3, 1)) {
                transition
            } else {
                transversion
            };
            raw[first * 4 + second] = (probability / 0.25).ln();
        }
    }
    normalize_matrix(raw, vec![0.25; 4], offset)
}

fn normalize_matrix(
    raw: Vec<f64>,
    frequencies: Vec<f64>,
    offset: f64,
) -> Result<(Vec<f64>, Vec<f64>), MafftError> {
    let count = frequencies.len();
    let mut average_identical = 0.0;
    let mut average_independent = 0.0;
    for first in 0..count {
        average_identical += frequencies[first] * raw[first * count + first];
        for second in 0..count {
            average_independent +=
                frequencies[first] * frequencies[second] * raw[first * count + second];
        }
    }
    let scale = average_identical - average_independent;
    if !scale.is_finite() || scale <= 0.0 {
        return Err(MafftError::NonFiniteComputation {
            stage: "Katoh similarity-matrix normalization",
        });
    }
    let normalized = raw
        .into_iter()
        .map(|score| (score - average_independent) / scale + offset)
        .collect::<Vec<_>>();
    Ok((normalized, frequencies))
}

fn matrix_power_20(mut base: [[f64; 20]; 20], mut exponent: u32) -> [[f64; 20]; 20] {
    let mut result = [[0.0; 20]; 20];
    for (index, row) in result.iter_mut().enumerate() {
        row[index] = 1.0;
    }
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = multiply_20(&result, &base);
        }
        exponent >>= 1;
        if exponent > 0 {
            base = multiply_20(&base, &base);
        }
    }
    result
}

fn multiply_20(first: &[[f64; 20]; 20], second: &[[f64; 20]; 20]) -> [[f64; 20]; 20] {
    let mut product = [[0.0; 20]; 20];
    for row in 0..20 {
        for (shared, second_row) in second.iter().enumerate() {
            let coefficient = first[row][shared];
            for (column, value) in second_row.iter().copied().enumerate() {
                product[row][column] += coefficient * value;
            }
        }
    }
    product
}

fn standardized(values: [f64; 20]) -> [f64; 20] {
    let mean = values.iter().sum::<f64>() / 20.0;
    let variance = values
        .iter()
        .map(|value| (value - mean) * (value - mean))
        .sum::<f64>()
        / 20.0;
    let standard_deviation = variance.sqrt();
    values.map(|value| (value - mean) / standard_deviation)
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence};

    use super::{BLOSUM62, JTT_PAM1_SCALED, ScientificScorer};
    use crate::{FftNs1Parameters, LINSI_LOCAL_OFFSET_PENALTY};

    #[test]
    fn published_jtt_rows_are_stochastic_up_to_print_rounding() {
        for row in JTT_PAM1_SCALED {
            let sum: u32 = row.iter().sum();
            assert!((99_990..=100_010).contains(&sum), "printed row sum {sum}");
        }
    }

    #[test]
    fn normalized_protein_moments_match_the_paper() {
        let scorer = ScientificScorer::new::<Protein>(FftNs1Parameters::paper_2002()).unwrap();
        let (identical, independent) = scorer.canonical_moments();
        assert!((identical - 1.06).abs() < 1.0e-12, "{identical}");
        assert!((independent - 0.06).abs() < 1.0e-12, "{independent}");
    }

    #[test]
    fn normalized_nucleotide_moments_match_the_paper() {
        let scorer = ScientificScorer::new::<Dna>(FftNs1Parameters::paper_2002()).unwrap();
        let (identical, independent) = scorer.canonical_moments();
        assert!((identical - 1.06).abs() < 1.0e-12, "{identical}");
        assert!((independent - 0.06).abs() < 1.0e-12, "{independent}");
    }

    #[test]
    fn ginsi_blosum62_and_k80_have_published_normalized_moments() {
        for (row_index, row) in BLOSUM62.iter().enumerate() {
            for (column_index, value) in row.iter().enumerate() {
                assert_eq!(*value, BLOSUM62[column_index][row_index]);
            }
        }
        assert_eq!(BLOSUM62[0][0], 4);
        assert_eq!(BLOSUM62[17][17], 11);
        assert_eq!(BLOSUM62[2][17], -4);
        let protein = ScientificScorer::ginsi_2005::<Protein>(0.123).unwrap();
        let dna = ScientificScorer::ginsi_2005::<Dna>(0.123).unwrap();
        let rna = ScientificScorer::ginsi_2005::<Rna>(0.123).unwrap();
        for scorer in [&protein, &dna, &rna] {
            let (identical, independent) = scorer.canonical_moments();
            assert!((identical - 1.123).abs() < 1.0e-12, "{identical}");
            assert!((independent - 0.123).abs() < 1.0e-12, "{independent}");
        }
    }

    #[test]
    fn linsi_local_moments_subtract_the_public_offset_penalty() {
        for scorer in [
            ScientificScorer::ginsi_2005::<Protein>(-LINSI_LOCAL_OFFSET_PENALTY).unwrap(),
            ScientificScorer::ginsi_2005::<Dna>(-LINSI_LOCAL_OFFSET_PENALTY).unwrap(),
            ScientificScorer::ginsi_2005::<Rna>(-LINSI_LOCAL_OFFSET_PENALTY).unwrap(),
        ] {
            let (identical, independent) = scorer.canonical_moments();
            assert!((identical - 0.9).abs() < 1.0e-12, "{identical}");
            assert!((independent + 0.1).abs() < 1.0e-12, "{independent}");
        }
    }

    #[test]
    fn ambiguity_scores_are_canonical_expectations() {
        let scorer = ScientificScorer::new::<Protein>(FftNs1Parameters::paper_2002()).unwrap();
        let sequence = ProteinSequence::from_ascii("BDN").unwrap();
        let b = sequence.get(0).unwrap();
        let d = sequence.get(1).unwrap();
        let n = sequence.get(2).unwrap();
        let expected = (scorer.score(d, d) + scorer.score(n, d)) * 0.5;
        assert!((scorer.score(b, d) - expected).abs() < 1.0e-12);
    }

    #[test]
    fn dna_and_rna_use_the_same_uracil_thymine_state() {
        let dna = ScientificScorer::new::<Dna>(FftNs1Parameters::paper_2002()).unwrap();
        let rna = ScientificScorer::new::<Rna>(FftNs1Parameters::paper_2002()).unwrap();
        let t = DnaSequence::from_ascii("T").unwrap().get(0).unwrap();
        let u = RnaSequence::from_ascii("U").unwrap().get(0).unwrap();
        assert!((dna.signal(t)[3] - rna.signal(u)[3]).abs() < f64::EPSILON);
    }
}
