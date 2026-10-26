use crate::MsaError;

const SYMMETRY_TOLERANCE: f64 = 1.0e-12;

/// A finite symmetric pairwise distance matrix with values in `[0,1]`.
#[derive(Clone, Debug, PartialEq)]
pub struct DistanceMatrix {
    size: usize,
    values: Box<[f64]>,
}

impl DistanceMatrix {
    /// Validates a row-major square distance matrix.
    ///
    /// # Errors
    ///
    /// Returns [`MsaError`] for a shape mismatch, invalid value, nonzero diagonal,
    /// or asymmetric pair.
    pub fn try_from_flat(size: usize, values: Vec<f64>) -> Result<Self, MsaError> {
        let expected = size
            .checked_mul(size)
            .ok_or(MsaError::DistanceMatrixSizeMismatch {
                size,
                expected: usize::MAX,
                actual: values.len(),
            })?;
        if values.len() != expected {
            return Err(MsaError::DistanceMatrixSizeMismatch {
                size,
                expected,
                actual: values.len(),
            });
        }
        for row in 0..size {
            for column in 0..size {
                let value = values[row * size + column];
                if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                    return Err(MsaError::InvalidDistance { row, column, value });
                }
                if row == column && value != 0.0 {
                    return Err(MsaError::NonZeroDistanceDiagonal { index: row, value });
                }
                if column < row {
                    let reverse = values[column * size + row];
                    if (value - reverse).abs() > SYMMETRY_TOLERANCE {
                        return Err(MsaError::AsymmetricDistance {
                            first: row,
                            second: column,
                            forward: value,
                            reverse,
                        });
                    }
                }
            }
        }
        Ok(Self {
            size,
            values: values.into_boxed_slice(),
        })
    }

    /// Returns the matrix dimension.
    #[must_use]
    pub const fn size(&self) -> usize {
        self.size
    }

    /// Returns one distance, if both indices exist.
    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> Option<f64> {
        if row < self.size && column < self.size {
            Some(self.values[row * self.size + column])
        } else {
            None
        }
    }

    /// Returns all row-major values.
    #[must_use]
    pub const fn as_flat_slice(&self) -> &[f64] {
        &self.values
    }
}
