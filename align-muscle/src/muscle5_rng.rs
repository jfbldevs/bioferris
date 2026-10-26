/// Fixed `SplitMix64` stream used by the versioned MUSCLE 5 profile.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub(crate) const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        value ^ (value >> 31)
    }

    #[allow(clippy::cast_precision_loss)]
    pub(crate) fn unit_f64(&mut self) -> f64 {
        const DENOMINATOR: f64 = 9_007_199_254_740_992.0;
        (self.next_u64() >> 11) as f64 / DENOMINATOR
    }

    pub(crate) fn signed_unit_f64(&mut self) -> f64 {
        self.unit_f64().mul_add(2.0, -1.0)
    }

    pub(crate) fn index(&mut self, upper_exclusive: usize) -> usize {
        debug_assert!(upper_exclusive > 0);
        let upper = u64::try_from(upper_exclusive).unwrap_or(u64::MAX);
        let zone = u64::MAX - u64::MAX % upper;
        loop {
            let value = self.next_u64();
            if value < zone {
                return usize::try_from(value % upper).unwrap_or(0);
            }
        }
    }
}

pub(crate) const HMM_DOMAIN: u64 = 0x484D_4D50_4552_5455;
pub(crate) const REFINE_DOMAIN: u64 = 0x5245_4649_4E45_0001;
pub(crate) const CLUSTER_DOMAIN: u64 = 0x434C_5553_5445_5201;
pub(crate) const SAMPLE_DOMAIN: u64 = 0x5341_4D50_4C45_0001;

#[allow(clippy::cast_possible_truncation)]
pub(crate) const fn derived_seed(seed: u64, domain: u64) -> u64 {
    seed ^ domain.rotate_left((seed & 63) as u32)
}

#[cfg(test)]
mod tests {
    use super::{HMM_DOMAIN, SplitMix64, derived_seed};

    #[test]
    fn stream_has_stable_known_prefix() {
        let mut rng = SplitMix64::new(0);
        assert_eq!(rng.next_u64(), 0xE220_A839_7B1D_CDAF);
        assert_eq!(rng.next_u64(), 0x6E78_9E6A_A1B9_65F4);
        assert_eq!(rng.next_u64(), 0x06C4_5D18_8009_454F);
    }

    #[test]
    fn domains_separate_equal_user_seeds() {
        assert_ne!(derived_seed(7, HMM_DOMAIN), derived_seed(8, HMM_DOMAIN));
        let mut rng = SplitMix64::new(derived_seed(7, HMM_DOMAIN));
        assert!((-1.0..=1.0).contains(&rng.signed_unit_f64()));
        assert!(rng.index(3) < 3);
    }
}
