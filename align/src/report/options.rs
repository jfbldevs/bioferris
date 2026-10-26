use core::fmt;

/// Default number of alignment columns rendered in one block.
pub const DEFAULT_BLOCK_WIDTH: usize = 60;

/// Amount of scientific and execution detail included in a report.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Verbosity {
    /// Algorithm identity, dimensions, and primary objective only.
    Compact,
    /// Summary statistics plus a block alignment.
    #[default]
    Normal,
    /// Normal output plus retained algorithm-specific evidence.
    Detailed,
}

/// Requested terminal-color behavior.
///
/// Deterministic [`core::fmt::Display`] reports are always colorless. This value
/// is retained for terminal writers and future front ends without changing the
/// stable text contract.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ColorMode {
    /// Never emit ANSI terminal escapes.
    #[default]
    Never,
    /// Emit color only when a terminal-aware writer confirms support.
    Auto,
    /// Request color from a terminal-aware writer.
    Always,
}

/// Coordinate representation used by a report.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum CoordinateStyle {
    /// One-based inclusive coordinates intended for people.
    #[default]
    Human,
    /// Zero-based half-open ranges matching native Rust results.
    Rust,
}

/// Shared immutable rendering options.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReportOptions {
    block_width: usize,
    verbosity: Verbosity,
    color: ColorMode,
    coordinates: CoordinateStyle,
}

impl ReportOptions {
    /// Creates the stable normal report configuration.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            block_width: DEFAULT_BLOCK_WIDTH,
            verbosity: Verbosity::Normal,
            color: ColorMode::Never,
            coordinates: CoordinateStyle::Human,
        }
    }

    /// Replaces the positive block width.
    ///
    /// # Errors
    ///
    /// Returns [`ReportError::InvalidBlockWidth`] when `width` is zero.
    pub const fn with_block_width(mut self, width: usize) -> Result<Self, ReportError> {
        if width == 0 {
            return Err(ReportError::InvalidBlockWidth { value: width });
        }
        self.block_width = width;
        Ok(self)
    }

    /// Replaces the report verbosity.
    #[must_use]
    pub const fn with_verbosity(mut self, verbosity: Verbosity) -> Self {
        self.verbosity = verbosity;
        self
    }

    /// Retains a color preference for terminal-aware writers.
    #[must_use]
    pub const fn with_color(mut self, color: ColorMode) -> Self {
        self.color = color;
        self
    }

    /// Replaces the coordinate representation.
    #[must_use]
    pub const fn with_coordinates(mut self, coordinates: CoordinateStyle) -> Self {
        self.coordinates = coordinates;
        self
    }

    /// Returns the number of alignment columns in one block.
    #[must_use]
    pub const fn block_width(self) -> usize {
        self.block_width
    }

    /// Returns the selected verbosity.
    #[must_use]
    pub const fn verbosity(self) -> Verbosity {
        self.verbosity
    }

    /// Returns the retained terminal-color preference.
    #[must_use]
    pub const fn color(self) -> ColorMode {
        self.color
    }

    /// Returns the selected coordinate representation.
    #[must_use]
    pub const fn coordinates(self) -> CoordinateStyle {
        self.coordinates
    }
}

impl Default for ReportOptions {
    fn default() -> Self {
        Self::new()
    }
}

/// Invalid human-report configuration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum ReportError {
    /// A report block must contain at least one column.
    InvalidBlockWidth {
        /// Rejected width.
        value: usize,
    },
}

impl fmt::Display for ReportError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBlockWidth { value } => {
                write!(
                    formatter,
                    "report block width must be positive, got {value}"
                )
            }
        }
    }
}

impl std::error::Error for ReportError {}
