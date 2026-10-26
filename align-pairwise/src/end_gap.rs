/// Permissions for omitting each unaligned sequence end without a gap score.
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct FreeEndGaps {
    first_start: bool,
    first_end: bool,
    second_start: bool,
    second_end: bool,
}

impl FreeEndGaps {
    /// Creates an exact four-end configuration.
    #[must_use]
    #[allow(clippy::fn_params_excessive_bools)]
    pub const fn new(
        first_start: bool,
        first_end: bool,
        second_start: bool,
        second_end: bool,
    ) -> Self {
        Self {
            first_start,
            first_end,
            second_start,
            second_end,
        }
    }

    /// Requires both complete sequences, reproducing global boundaries.
    #[must_use]
    pub const fn global() -> Self {
        Self::new(false, false, false, false)
    }

    /// Makes all four overhangs free, the symmetric semiglobal preset.
    #[must_use]
    pub const fn semiglobal() -> Self {
        Self::new(true, true, true, true)
    }

    /// Aligns the complete first sequence to a substring of the second.
    #[must_use]
    pub const fn fit_first_in_second() -> Self {
        Self::new(false, false, true, true)
    }

    /// Aligns the complete second sequence to a substring of the first.
    #[must_use]
    pub const fn fit_second_in_first() -> Self {
        Self::new(true, true, false, false)
    }

    /// Aligns a suffix of the first sequence to a prefix of the second.
    #[must_use]
    pub const fn first_suffix_second_prefix() -> Self {
        Self::new(true, false, false, true)
    }

    /// Aligns a suffix of the second sequence to a prefix of the first.
    #[must_use]
    pub const fn second_suffix_first_prefix() -> Self {
        Self::new(false, true, true, false)
    }

    /// Returns whether a first-sequence prefix may be omitted.
    #[must_use]
    pub const fn first_start(self) -> bool {
        self.first_start
    }

    /// Returns whether a first-sequence suffix may be omitted.
    #[must_use]
    pub const fn first_end(self) -> bool {
        self.first_end
    }

    /// Returns whether a second-sequence prefix may be omitted.
    #[must_use]
    pub const fn second_start(self) -> bool {
        self.second_start
    }

    /// Returns whether a second-sequence suffix may be omitted.
    #[must_use]
    pub const fn second_end(self) -> bool {
        self.second_end
    }
}
