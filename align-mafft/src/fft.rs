use core::mem::size_of;
use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

use rustfft::num_complex::Complex;
use rustfft::{Fft, FftPlanner};

use crate::MafftError;

#[derive(Clone)]
struct TransformPlan {
    forward: Arc<dyn Fft<f64>>,
    inverse: Arc<dyn Fft<f64>>,
    scratch_length: usize,
}

struct TransformPlanCache {
    planner: FftPlanner<f64>,
    plans: HashMap<usize, TransformPlan>,
}

impl TransformPlanCache {
    fn new() -> Self {
        Self {
            planner: FftPlanner::new(),
            plans: HashMap::new(),
        }
    }

    fn get(&mut self, length: usize) -> TransformPlan {
        if let Some(plan) = self.plans.get(&length) {
            return plan.clone();
        }
        let forward = self.planner.plan_fft_forward(length);
        let inverse = self.planner.plan_fft_inverse(length);
        let plan = TransformPlan {
            scratch_length: forward
                .get_inplace_scratch_len()
                .max(inverse.get_inplace_scratch_len()),
            forward,
            inverse,
        };
        self.plans.insert(length, plan.clone());
        plan
    }
}

fn transform_plan(length: usize) -> TransformPlan {
    static CACHE: OnceLock<Mutex<TransformPlanCache>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(TransformPlanCache::new()));
    cache
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .get(length)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct LagCorrelation {
    pub(crate) lag: isize,
    pub(crate) score: f64,
}

#[derive(Clone, Debug)]
pub(crate) struct CorrelationResult {
    pub(crate) values: Vec<LagCorrelation>,
    pub(crate) workspace_bytes: usize,
}

pub(crate) fn correlate_channels(
    first: &[Vec<f64>],
    second: &[Vec<f64>],
) -> Result<CorrelationResult, MafftError> {
    if first.len() != second.len() {
        return Err(MafftError::InternalInvariant {
            message: "FFT profile channel counts differ",
        });
    }
    let first_length = first.first().map_or(0, Vec::len);
    let second_length = second.first().map_or(0, Vec::len);
    if first.iter().any(|channel| channel.len() != first_length)
        || second.iter().any(|channel| channel.len() != second_length)
    {
        return Err(MafftError::InternalInvariant {
            message: "FFT profile channels are ragged",
        });
    }
    if first_length == 0 || second_length == 0 {
        return Ok(CorrelationResult {
            values: Vec::new(),
            workspace_bytes: 0,
        });
    }
    let convolution_length = first_length
        .checked_add(second_length)
        .and_then(|sum| sum.checked_sub(1))
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "FFT convolution length",
        })?;
    let transform_length =
        convolution_length
            .checked_next_power_of_two()
            .ok_or(MafftError::ArithmeticOverflow {
                stage: "FFT transform padding",
            })?;
    let mut accumulated = vec![0.0; convolution_length];
    let mut first_buffer = vec![Complex::default(); transform_length];
    let mut second_buffer = vec![Complex::default(); transform_length];
    let plan = transform_plan(transform_length);
    let mut scratch = vec![Complex::default(); plan.scratch_length];
    #[allow(clippy::cast_precision_loss)]
    let inverse_scale = 1.0 / transform_length as f64;

    for (first_channel, second_channel) in first.iter().zip(second) {
        first_buffer.fill(Complex::default());
        second_buffer.fill(Complex::default());
        for (index, value) in first_channel.iter().rev().copied().enumerate() {
            first_buffer[index].re = value;
        }
        for (index, value) in second_channel.iter().copied().enumerate() {
            second_buffer[index].re = value;
        }
        plan.forward
            .process_with_scratch(&mut first_buffer, &mut scratch);
        plan.forward
            .process_with_scratch(&mut second_buffer, &mut scratch);
        for index in 0..transform_length {
            first_buffer[index] *= second_buffer[index];
        }
        plan.inverse
            .process_with_scratch(&mut first_buffer, &mut scratch);
        for (target, value) in accumulated.iter_mut().zip(&first_buffer) {
            *target += value.re * inverse_scale;
        }
    }

    let center = isize::try_from(first_length - 1).map_err(|_| MafftError::ArithmeticOverflow {
        stage: "FFT lag conversion",
    })?;
    let mut values = Vec::with_capacity(convolution_length);
    for (index, score) in accumulated.into_iter().enumerate() {
        let index = isize::try_from(index).map_err(|_| MafftError::ArithmeticOverflow {
            stage: "FFT lag conversion",
        })?;
        if !score.is_finite() {
            return Err(MafftError::NonFiniteComputation {
                stage: "FFT cross-correlation",
            });
        }
        values.push(LagCorrelation {
            lag: index - center,
            score,
        });
    }
    let workspace_bytes = transform_length
        .checked_mul(2)
        .and_then(|elements| elements.checked_add(plan.scratch_length))
        .and_then(|elements| elements.checked_mul(size_of::<Complex<f64>>()))
        .and_then(|bytes| {
            convolution_length
                .checked_mul(size_of::<f64>())
                .and_then(|accumulated_bytes| bytes.checked_add(accumulated_bytes))
        })
        .ok_or(MafftError::ArithmeticOverflow {
            stage: "FFT workspace accounting",
        })?;
    Ok(CorrelationResult {
        values,
        workspace_bytes,
    })
}

pub(crate) fn best_lags(values: &[LagCorrelation], count: usize) -> Vec<LagCorrelation> {
    let mut selected = values.to_vec();
    selected.sort_by(|first, second| {
        second
            .score
            .total_cmp(&first.score)
            .then_with(|| first.lag.unsigned_abs().cmp(&second.lag.unsigned_abs()))
            .then_with(|| first.lag.cmp(&second.lag))
    });
    selected.truncate(count.min(selected.len()));
    selected
}

#[cfg(test)]
mod tests {
    use super::{best_lags, correlate_channels};

    fn deterministic_signals(channels: usize, length: usize, salt: usize) -> Vec<Vec<f64>> {
        (0..channels)
            .map(|channel| {
                (0..length)
                    .map(|index| {
                        let value = (index * (channel + 3) + salt * 5 + channel) % 11;
                        f64::from(u32::try_from(value).unwrap()) - 5.0
                    })
                    .collect()
            })
            .collect()
    }

    fn direct(first: &[Vec<f64>], second: &[Vec<f64>]) -> Vec<(isize, f64)> {
        let first_length = first[0].len();
        let second_length = second[0].len();
        let minimum_lag = -(isize::try_from(first_length).unwrap() - 1);
        let maximum_lag = isize::try_from(second_length).unwrap() - 1;
        (minimum_lag..=maximum_lag)
            .map(|lag| {
                let mut score = 0.0;
                for (first_channel, second_channel) in first.iter().zip(second) {
                    for (first_index, first_value) in first_channel.iter().copied().enumerate() {
                        let second_index = isize::try_from(first_index).unwrap() + lag;
                        if let Ok(second_index) = usize::try_from(second_index) {
                            if second_index < second_length {
                                score += first_value * second_channel[second_index];
                            }
                        }
                    }
                }
                (lag, score)
            })
            .collect()
    }

    #[test]
    fn fft_matches_direct_cross_correlation() {
        let first = vec![vec![1.0, -2.0, 3.0], vec![0.5, 1.5, -1.0]];
        let second = vec![vec![2.0, 4.0], vec![-1.0, 0.25]];
        let observed = correlate_channels(&first, &second).unwrap().values;
        let expected = direct(&first, &second);
        assert_eq!(observed.len(), expected.len());
        for (observed, (lag, score)) in observed.iter().zip(expected) {
            assert_eq!(observed.lag, lag);
            assert!((observed.score - score).abs() < 1.0e-10);
        }
    }

    #[test]
    fn exhaustive_short_binary_signals_match_direct_correlation() {
        for first_length in 1..=5 {
            for second_length in 1..=5 {
                for first_mask in 0..(1usize << first_length) {
                    for second_mask in 0..(1usize << second_length) {
                        let first = vec![
                            (0..first_length)
                                .map(|index| {
                                    if first_mask & (1 << index) == 0 {
                                        -1.0
                                    } else {
                                        1.0
                                    }
                                })
                                .collect::<Vec<_>>(),
                        ];
                        let second = vec![
                            (0..second_length)
                                .map(|index| {
                                    if second_mask & (1 << index) == 0 {
                                        -1.0
                                    } else {
                                        1.0
                                    }
                                })
                                .collect::<Vec<_>>(),
                        ];
                        let observed = correlate_channels(&first, &second).unwrap().values;
                        let expected = direct(&first, &second);
                        for (observed, (lag, score)) in observed.iter().zip(expected) {
                            assert_eq!(observed.lag, lag);
                            assert!((observed.score - score).abs() < 1.0e-9);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn peak_ties_prefer_small_absolute_then_signed_lag() {
        let values = super::LagCorrelation { lag: 0, score: 0.0 };
        let selected = best_lags(
            &[
                super::LagCorrelation { lag: 2, score: 1.0 },
                super::LagCorrelation {
                    lag: -1,
                    score: 1.0,
                },
                super::LagCorrelation { lag: 1, score: 1.0 },
                values,
            ],
            3,
        );
        assert_eq!(
            selected.iter().map(|value| value.lag).collect::<Vec<_>>(),
            vec![-1, 1, 2]
        );
    }

    #[test]
    fn cached_plans_are_deterministic_across_concurrent_lengths() {
        std::thread::scope(|scope| {
            for (first_length, second_length, channels) in
                [(7, 5, 2), (17, 11, 4), (33, 19, 2), (65, 37, 4)]
            {
                scope.spawn(move || {
                    let first = deterministic_signals(channels, first_length, first_length);
                    let second = deterministic_signals(channels, second_length, second_length + 1);
                    let expected = direct(&first, &second);
                    let baseline = correlate_channels(&first, &second).unwrap();
                    for (observed, (lag, score)) in baseline.values.iter().zip(expected) {
                        assert_eq!(observed.lag, lag);
                        assert!((observed.score - score).abs() < 1.0e-8);
                    }
                    for _ in 0..4 {
                        let repeated = correlate_channels(&first, &second).unwrap();
                        assert_eq!(repeated.values, baseline.values);
                        assert_eq!(repeated.workspace_bytes, baseline.workspace_bytes);
                    }
                });
            }
        });
    }
}
