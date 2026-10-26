//! Published Gonnet observations and explicit ambiguity projections.

use crate::ClustalOmegaError;

pub(crate) const PROTEIN_CANONICAL_SIZE: usize = 20;
pub(crate) const NUCLEOTIDE_CANONICAL_SIZE: usize = 4;

// NCBI BLAST GONNET, canonical order A R N D C Q E G H I L K M F P S T W Y V.
// Retrieved 2026-08-03 from https://ftp.ncbi.nih.gov/blast/matrices/GONNET.
// SHA-256 of the immutable ASCII observation:
// c4d247a9c9aa485b4ed5f1422fb435ae1afe315b233cb8bac8c99724b4125135
const GONNET: [i8; 400] = [
    2, -1, 0, 0, 0, 0, 0, 0, -1, -1, -1, 0, -1, -2, 0, 1, 1, -4, -2, 0, -1, 5, 0, 0, -2, 2, 0, -1,
    1, -2, -2, 3, -2, -3, -1, 0, 0, -2, -2, -2, 0, 0, 4, 2, -2, 1, 1, 0, 1, -3, -3, 1, -2, -3, -1,
    1, 0, -4, -1, -2, 0, 0, 2, 5, -3, 1, 3, 0, 0, -4, -4, 0, -3, -4, -1, 0, 0, -5, -3, -3, 0, -2,
    -2, -3, 12, -2, -3, -2, -1, -1, -2, -3, -1, -1, -3, 0, 0, -1, 0, 0, 0, 2, 1, 1, -2, 3, 2, -1,
    1, -2, -2, 2, -1, -3, 0, 0, 0, -3, -2, -2, 0, 0, 1, 3, -3, 2, 4, -1, 0, -3, -3, 1, -2, -4, 0,
    0, 0, -4, -3, -2, 0, -1, 0, 0, -2, -1, -1, 7, -1, -4, -4, -1, -4, -5, -2, 0, -1, -4, -4, -3,
    -1, 1, 1, 0, -1, 1, 0, -1, 6, -2, -2, 1, -1, 0, -1, 0, 0, -1, 2, -2, -1, -2, -3, -4, -1, -2,
    -3, -4, -2, 4, 3, -2, 2, 1, -3, -2, -1, -2, -1, 3, -1, -2, -3, -4, -2, -2, -3, -4, -2, 3, 4,
    -2, 3, 2, -2, -2, -1, -1, 0, 2, 0, 3, 1, 0, -3, 2, 1, -1, 1, -2, -2, 3, -1, -3, -1, 0, 0, -4,
    -2, -2, -1, -2, -2, -3, -1, -1, -2, -4, -1, 2, 3, -1, 4, 2, -2, -1, -1, -1, 0, 2, -2, -3, -3,
    -4, -1, -3, -4, -5, 0, 1, 2, -3, 2, 7, -4, -3, -2, 4, 5, 0, 0, -1, -1, -1, -3, 0, 0, -2, -1,
    -3, -2, -1, -2, -4, 8, 0, 0, -5, -3, -2, 1, 0, 1, 0, 0, 0, 0, 0, 0, -2, -2, 0, -1, -3, 0, 2, 2,
    -3, -2, -1, 1, 0, 0, 0, 0, 0, 0, -1, 0, -1, -1, 0, -1, -2, 0, 2, 2, -4, -2, 0, -4, -2, -4, -5,
    -1, -3, -4, -4, -1, -2, -1, -4, -1, 4, -5, -3, -4, 14, 4, -3, -2, -2, -1, -3, 0, -2, -3, -4, 2,
    -1, 0, -2, 0, 5, -3, -2, -2, 4, 8, -1, 0, -2, -2, -3, 0, -2, -2, -3, -2, 3, 2, -2, 2, 0, -2,
    -1, 0, -3, -1, 3,
];

// Robinson-Robinson residue observations, used as the explicit fixed database
// background required by Soding's log-odds equation. Integer normalization avoids
// storing rounded probabilities as scientific inputs.
const PROTEIN_BACKGROUND_COUNTS: [u32; 20] = [
    35_155, 23_105, 20_212, 24_161, 8_669, 19_208, 28_354, 33_229, 9_906, 23_161, 40_625, 25_872,
    10_101, 17_367, 23_435, 32_070, 26_311, 5_990, 14_488, 29_012,
];

pub(crate) fn canonical_size(is_protein: bool) -> usize {
    if is_protein {
        PROTEIN_CANONICAL_SIZE
    } else {
        NUCLEOTIDE_CANONICAL_SIZE
    }
}

#[allow(clippy::cast_precision_loss)]
pub(crate) fn background(is_protein: bool) -> Vec<f64> {
    if !is_protein {
        return vec![0.25; NUCLEOTIDE_CANONICAL_SIZE];
    }
    let total = PROTEIN_BACKGROUND_COUNTS
        .iter()
        .map(|value| u64::from(*value))
        .sum::<u64>();
    PROTEIN_BACKGROUND_COUNTS
        .iter()
        .map(|value| f64::from(*value) / total as f64)
        .collect()
}

pub(crate) fn add_projected_symbol(
    counts: &mut [f64],
    ascii: u8,
    weight: f64,
    is_protein: bool,
) -> Result<(), ClustalOmegaError> {
    if is_protein {
        add_projected_protein(counts, ascii, weight)
    } else {
        add_projected_nucleotide(counts, ascii, weight)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct SubstitutionModel {
    size: usize,
    conditional: Vec<f64>,
}

impl SubstitutionModel {
    pub(crate) fn new(is_protein: bool) -> Result<Self, ClustalOmegaError> {
        let size = canonical_size(is_protein);
        let background = background(is_protein);
        let mut conditional = vec![0.0; size * size];
        let mut normalizers = vec![0.0; size];
        for source in 0..size {
            for target in 0..size {
                let score = if is_protein {
                    f64::from(GONNET[target * size + source])
                } else {
                    f64::from(target == source)
                };
                let value = background[target] * 2.0_f64.powf(score / 2.0);
                conditional[source * size + target] = value;
                normalizers[source] += value;
            }
            if !normalizers[source].is_finite() || normalizers[source] <= 0.0 {
                return Err(ClustalOmegaError::InvalidProbability {
                    stage: "substitution pseudocount normalization",
                    value: normalizers[source],
                });
            }
        }
        for source in 0..size {
            for value in &mut conditional[source * size..(source + 1) * size] {
                *value /= normalizers[source];
            }
        }
        Ok(Self { size, conditional })
    }
}

#[cfg(test)]
pub(crate) fn substitution_prior(
    observed: &[f64],
    model: &SubstitutionModel,
) -> Result<Vec<f64>, ClustalOmegaError> {
    let mut prior = vec![0.0; model.size];
    substitution_prior_into(observed, model, &mut prior)?;
    Ok(prior)
}

pub(crate) fn substitution_prior_into(
    observed: &[f64],
    model: &SubstitutionModel,
    prior: &mut [f64],
) -> Result<(), ClustalOmegaError> {
    let size = model.size;
    if observed.len() != size || prior.len() != size {
        return Err(ClustalOmegaError::InternalInvariant {
            message: "emission vector has the wrong canonical alphabet size",
        });
    }
    prior.fill(0.0);
    for (source, observed_frequency) in observed.iter().copied().enumerate() {
        let row = &model.conditional[source * size..(source + 1) * size];
        for (prior, conditional) in prior.iter_mut().zip(row) {
            *prior += observed_frequency * conditional;
        }
    }
    normalize_positive(prior, "substitution pseudocount distribution")
}

pub(crate) fn normalize_positive(
    values: &mut [f64],
    stage: &'static str,
) -> Result<(), ClustalOmegaError> {
    const MINIMUM: f64 = 1.0e-15;
    for value in values.iter_mut() {
        if !value.is_finite() || *value < 0.0 {
            return Err(ClustalOmegaError::InvalidProbability {
                stage,
                value: *value,
            });
        }
        *value = value.max(MINIMUM);
    }
    let total = values.iter().sum::<f64>();
    if !total.is_finite() || total <= 0.0 {
        return Err(ClustalOmegaError::InvalidProbability {
            stage,
            value: total,
        });
    }
    for value in values {
        *value /= total;
    }
    Ok(())
}

fn add_projected_protein(
    counts: &mut [f64],
    ascii: u8,
    weight: f64,
) -> Result<(), ClustalOmegaError> {
    let canonical = match ascii {
        b'A' => Some(0),
        b'R' => Some(1),
        b'N' => Some(2),
        b'D' => Some(3),
        b'C' => Some(4),
        b'Q' => Some(5),
        b'E' => Some(6),
        b'G' => Some(7),
        b'H' => Some(8),
        b'I' => Some(9),
        b'L' => Some(10),
        b'K' => Some(11),
        b'M' => Some(12),
        b'F' => Some(13),
        b'P' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'W' => Some(17),
        b'Y' => Some(18),
        b'V' => Some(19),
        _ => None,
    };
    if let Some(index) = canonical {
        counts[index] += weight;
        return Ok(());
    }
    match ascii {
        b'B' => add_uniform(counts, &[2, 3], weight),
        b'Z' => add_uniform(counts, &[5, 6], weight),
        b'J' => add_uniform(counts, &[9, 10], weight),
        b'U' => add_uniform(counts, &[4], weight),
        b'O' => add_uniform(counts, &[11], weight),
        b'X' => {
            let indices = (0..PROTEIN_CANONICAL_SIZE).collect::<Vec<_>>();
            add_uniform(counts, &indices, weight);
        }
        _ => {
            return Err(ClustalOmegaError::InternalInvariant {
                message: "validated protein symbol lacks an HMM projection",
            });
        }
    }
    Ok(())
}

fn add_projected_nucleotide(
    counts: &mut [f64],
    ascii: u8,
    weight: f64,
) -> Result<(), ClustalOmegaError> {
    let bits = match ascii {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101,
        b'Y' => 0b1010,
        b'S' => 0b0110,
        b'W' => 0b1001,
        b'K' => 0b1100,
        b'M' => 0b0011,
        b'B' => 0b1110,
        b'D' => 0b1101,
        b'H' => 0b1011,
        b'V' => 0b0111,
        b'N' => 0b1111,
        _ => {
            return Err(ClustalOmegaError::InternalInvariant {
                message: "validated nucleotide symbol lacks an HMM projection",
            });
        }
    };
    let indices = (0..NUCLEOTIDE_CANONICAL_SIZE)
        .filter(|index| bits & (1 << index) != 0)
        .collect::<Vec<_>>();
    add_uniform(counts, &indices, weight);
    Ok(())
}

#[allow(clippy::cast_precision_loss)]
fn add_uniform(counts: &mut [f64], indices: &[usize], weight: f64) {
    let contribution = weight / indices.len() as f64;
    for index in indices {
        counts[*index] += contribution;
    }
}

#[cfg(test)]
mod tests {
    use super::{
        NUCLEOTIDE_CANONICAL_SIZE, PROTEIN_CANONICAL_SIZE, SubstitutionModel, add_projected_symbol,
        background, substitution_prior,
    };

    #[test]
    fn ambiguity_projection_conserves_observation_weight() {
        let mut protein = vec![0.0; PROTEIN_CANONICAL_SIZE];
        add_projected_symbol(&mut protein, b'B', 3.0, true).unwrap();
        assert!((protein[2] - 1.5).abs() < f64::EPSILON);
        assert!((protein[3] - 1.5).abs() < f64::EPSILON);
        assert!((protein.iter().sum::<f64>() - 3.0).abs() < f64::EPSILON);

        let mut nucleotide = vec![0.0; NUCLEOTIDE_CANONICAL_SIZE];
        add_projected_symbol(&mut nucleotide, b'N', 4.0, false).unwrap();
        assert_eq!(nucleotide, vec![1.0; 4]);
    }

    #[test]
    fn backgrounds_and_substitution_priors_are_normalized_and_positive() {
        for is_protein in [false, true] {
            let background = background(is_protein);
            assert!((background.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
            let model = SubstitutionModel::new(is_protein).unwrap();
            let prior = substitution_prior(&background, &model).unwrap();
            assert!((prior.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
            assert!(prior.iter().all(|value| value.is_finite() && *value > 0.0));
        }
    }
}
