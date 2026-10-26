//! Specified physical-property estimates for canonical, unmodified peptides.
//!
//! The input is borrowed and must contain at least one canonical amino acid.
//! Average molecular masses describe a linear peptide with free termini; charge
//! and pI use an independent-ionization model with sequence-adjusted terminal
//! parameters. These estimates do not include modifications or structural effects.
//! All results are unrounded binary64 values; floating-point identities observed
//! for a particular reference runtime do not promise identity on every platform.

use std::fmt;

use crate::ProteinSequence;

const AVERAGE_WATER_MASS: f64 = 18.0153;
const PI_LOWER: f64 = 4.05;
const PI_UPPER: f64 = 12.0;
const PI_INTERVAL_WIDTH: f64 = 0.0001;

// Canonical order: A C D E F G H I K L M N P Q R S T V W Y.
// Free amino acid masses from Young et al. (2018), supplementary Table S1.
const AVERAGE_FREE_MASSES: [f64; 20] = [
    89.0932, 121.1582, 133.1027, 147.1293, 165.1891, 75.0666, 155.1546, 131.1729, 146.1876,
    131.1729, 149.2113, 132.1179, 115.1305, 146.1445, 174.2010, 105.0926, 119.1192, 117.1463,
    204.2252, 181.1885,
];

// The original Kyte-Doolittle scale as published by ExPASy ProtScale.
const KYTE_DOOLITTLE: [f64; 20] = [
    1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7,
    4.2, -0.9, -1.3,
];

/// A peptide or numerical option outside the supported property model.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum ProteinPropertyError {
    /// A peptide must contain at least one residue.
    EmptySequence,
    /// A typed protein contains a noncanonical residue at this zero-based offset.
    ///
    /// B/Z/J/X and concrete U/O remain valid in the general protein alphabet,
    /// but this property's canonical-20 parameter tables do not support them.
    UnsupportedResidue {
        /// Zero-based offset in the original validated protein.
        offset: usize,
        /// Unsupported uppercase ASCII residue.
        byte: u8,
    },
    /// The pH must be finite and lie in the inclusive interval 0..14.
    InvalidPh {
        /// Rejected caller-supplied value, including NaN or infinity.
        ph: f64,
    },
}

impl ProteinPropertyError {
    /// Returns the stable machine-readable error identifier.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::EmptySequence => "empty_protein",
            Self::UnsupportedResidue { .. } => "unsupported_protein_residue",
            Self::InvalidPh { .. } => "invalid_ph",
        }
    }
}

impl fmt::Display for ProteinPropertyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySequence => {
                formatter.write_str("protein property analysis requires at least one residue")
            }
            Self::UnsupportedResidue { offset, byte } => write!(
                formatter,
                "protein property analysis requires canonical amino acids; found {} at residue {offset}",
                char::from(*byte)
            ),
            Self::InvalidPh { ph } => write!(
                formatter,
                "protein charge pH must be finite and in 0..=14 (received {ph})"
            ),
        }
    }
}

impl std::error::Error for ProteinPropertyError {}

/// Charge-estimation options for the aggregate property workflow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProteinPropertiesOptions {
    /// Finite pH in inclusive 0..14; the default is 7.0.
    pub ph: f64,
}

impl Default for ProteinPropertiesOptions {
    fn default() -> Self {
        Self { ph: 7.0 }
    }
}

/// Whether a bounded pI calculation approaches an interval boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PiBoundary {
    /// The charge changes sign within the inclusive 4.05..12 interval.
    None,
    /// Charge is already negative at 4.05; the model root is lower.
    Lower,
    /// Charge remains positive at 12; the model root is higher.
    Upper,
}

/// A bounded pI estimate and the status of its numerical bracket.
///
/// Bisection starts with a probe at 7.775 in 4.05..12 and returns the midpoint once interval
/// width is at most 0.0001. A boundary-limited value is not a zero-charge root.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsoelectricPointEstimate {
    value: f64,
    boundary: PiBoundary,
}

impl IsoelectricPointEstimate {
    /// Returns the unrounded bounded midpoint estimate in pH units.
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }

    /// Reports whether the model root lies beyond a bracket boundary.
    #[must_use]
    pub const fn boundary(&self) -> PiBoundary {
        self.boundary
    }

    /// Returns true if the model charge changes sign within inclusive 4.05..12.
    #[must_use]
    pub const fn bracketed(&self) -> bool {
        matches!(self.boundary, PiBoundary::None)
    }
}

/// Immutable results for all five canonical-peptide properties.
///
/// Values are unrounded binary64 estimates. Mass is in Da, charge in elementary
/// charge units, pH/pI are dimensionless, aromaticity is a fraction in 0..1,
/// and GRAVY is the mean Kyte-Doolittle hydropathy value.
#[derive(Clone, Debug, PartialEq)]
pub struct ProteinProperties {
    length: usize,
    molecular_weight: f64,
    ph: f64,
    charge: f64,
    pi_estimate: IsoelectricPointEstimate,
    aromaticity: f64,
    gravy: f64,
}

impl ProteinProperties {
    /// Returns the number of amino acid residues.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.length
    }

    /// Returns false: successful property analysis always has a nonempty input.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// Returns average molecular mass in Da for an unmodified linear peptide.
    #[must_use]
    pub const fn molecular_weight(&self) -> f64 {
        self.molecular_weight
    }

    /// Returns the finite pH used for the charge estimate.
    #[must_use]
    pub const fn ph(&self) -> f64 {
        self.ph
    }

    /// Returns net charge in elementary charge units at [`Self::ph`].
    #[must_use]
    pub const fn charge(&self) -> f64 {
        self.charge
    }

    /// Returns the bounded pI midpoint; see [`Self::pi_estimate`] for its status.
    #[must_use]
    pub const fn isoelectric_point(&self) -> f64 {
        self.pi_estimate.value
    }

    /// Returns the bounded pI value with an explicit bracket/boundary status.
    #[must_use]
    pub const fn pi_estimate(&self) -> IsoelectricPointEstimate {
        self.pi_estimate
    }

    /// Returns the combined fraction of phenylalanine, tryptophan and tyrosine.
    #[must_use]
    pub const fn aromaticity(&self) -> f64 {
        self.aromaticity
    }

    /// Returns the grand average of Kyte-Doolittle hydropathy.
    #[must_use]
    pub const fn gravy(&self) -> f64 {
        self.gravy
    }
}

/// Borrowed, validated input and exact counts for a canonical peptide analysis.
///
/// Construction performs one linear pass and retains a fixed 20-entry count
/// array. Mass and GRAVY reductions follow source-sequence order; other property
/// queries use retained counts and the original terminal residues. The source
/// sequence is neither cloned nor modified.
#[derive(Clone, Debug)]
pub struct ProteinAnalysis<'a> {
    sequence: &'a ProteinSequence,
    counts: [usize; 20],
    n_terminal_pk: f64,
    c_terminal_pk: f64,
}

impl<'a> ProteinAnalysis<'a> {
    /// Validates and borrows a nonempty sequence of canonical amino acids.
    ///
    /// # Errors
    /// Returns [`ProteinPropertyError::EmptySequence`] for empty input or
    /// [`ProteinPropertyError::UnsupportedResidue`] for its first B/Z/J/X/U/O.
    pub fn new(sequence: &'a ProteinSequence) -> Result<Self, ProteinPropertyError> {
        if sequence.is_empty() {
            return Err(ProteinPropertyError::EmptySequence);
        }
        let mut counts = [0; 20];
        let mut first = 0;
        let mut last = 0;
        for (offset, symbol) in sequence.iter().enumerate() {
            let byte = symbol.as_ascii();
            let index = canonical_index(byte)
                .ok_or(ProteinPropertyError::UnsupportedResidue { offset, byte })?;
            counts[index] += 1;
            if offset == 0 {
                first = byte;
            }
            last = byte;
        }
        let n_terminal_pk = match first {
            b'A' => 7.59,
            b'E' => 7.70,
            b'M' => 7.00,
            b'P' => 8.36,
            b'S' => 6.93,
            b'T' => 6.82,
            b'V' => 7.44,
            _ => 7.50,
        };
        let c_terminal_pk = match last {
            b'D' => 4.55,
            b'E' => 4.75,
            _ => 3.55,
        };
        Ok(Self {
            sequence,
            counts,
            n_terminal_pk,
            c_terminal_pk,
        })
    }

    /// Returns the number of validated residues.
    #[must_use]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Returns false because construction rejects empty input.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Calculates average molecular mass in Da with free, unmodified termini.
    ///
    /// Free amino acid masses are summed with a compensated scalar reduction in
    /// source order, then `(length - 1) * 18.0153` is subtracted for peptide bonds.
    ///
    /// # Panics
    /// Only if the private canonical-residue invariant is violated internally;
    /// validated construction and immutable borrowing prevent this through the safe API.
    #[must_use]
    pub fn molecular_weight(&self) -> f64 {
        let mass = compensated_sum(self.sequence.iter().map(|symbol| {
            // Construction checked every byte, and the borrowed source is immutable.
            let index = canonical_index(symbol.as_ascii())
                .expect("ProteinAnalysis retains only validated canonical residues");
            AVERAGE_FREE_MASSES[index]
        }));
        mass - count_as_float(self.len() - 1) * AVERAGE_WATER_MASS
    }

    /// Estimates charge at a finite pH within inclusive 0..14.
    ///
    /// The independent-ionization model includes K/R/H/D/E/C/Y side chains and
    /// one free group at each terminus. Sequence-specific D/E terminal parameters
    /// apply to the C-terminal group, not to those residues' side chains.
    /// Positive contributions are added in N-terminal/K/R/H order and negative
    /// contributions in C-terminal/D/E/C/Y order, then the totals are subtracted.
    ///
    /// # Errors
    /// Returns [`ProteinPropertyError::InvalidPh`] for NaN, infinity or pH outside
    /// the inclusive supported interval. No clamping is performed.
    pub fn charge_at_ph(&self, ph: f64) -> Result<f64, ProteinPropertyError> {
        validate_ph(ph)?;
        Ok(self.charge_unchecked(ph))
    }

    /// Returns the specified bounded pI estimate in pH units.
    ///
    /// Initial bounds are 4.05..12, the first probe is 7.775 and bisection
    /// stops at interval width at most 0.0001. A root outside this interval
    /// produces a boundary-limited value; [`Self::isoelectric_point_estimate`]
    /// makes that limitation explicit.
    #[must_use]
    pub fn isoelectric_point(&self) -> f64 {
        let mut lower = PI_LOWER;
        let mut upper = PI_UPPER;
        let mut midpoint = 7.775;
        while upper - lower > PI_INTERVAL_WIDTH {
            if self.charge_unchecked(midpoint) > 0.0 {
                lower = midpoint;
            } else {
                upper = midpoint;
            }
            midpoint = (lower + upper) / 2.0;
        }
        midpoint
    }

    /// Returns the bounded pI value and whether its interval brackets the root.
    ///
    /// Equality at an endpoint is bracketed. Values outside the interval retain
    /// the same bounded midpoint as [`Self::isoelectric_point`], with a lower or
    /// upper boundary marker rather than an assertion of zero net charge.
    #[must_use]
    pub fn isoelectric_point_estimate(&self) -> IsoelectricPointEstimate {
        let boundary = if self.charge_unchecked(PI_LOWER) < 0.0 {
            PiBoundary::Lower
        } else if self.charge_unchecked(PI_UPPER) > 0.0 {
            PiBoundary::Upper
        } else {
            PiBoundary::None
        };
        IsoelectricPointEstimate {
            value: self.isoelectric_point(),
            boundary,
        }
    }

    /// Returns the sum of F, W and Y residue frequencies, a fraction in 0..1.
    ///
    /// Each count is divided separately by total length before the three terms
    /// are added in F/W/Y order. No histogram or output rounding is applied.
    #[must_use]
    pub fn aromaticity(&self) -> f64 {
        let length = count_as_float(self.len());
        count_as_float(self.counts[4]) / length
            + count_as_float(self.counts[18]) / length
            + count_as_float(self.counts[19]) / length
    }

    /// Returns mean Kyte-Doolittle hydropathy, with a compensated source-order sum.
    ///
    /// # Panics
    /// Only if the private canonical-residue invariant is violated internally;
    /// validated construction and immutable borrowing prevent this through the safe API.
    #[must_use]
    pub fn gravy(&self) -> f64 {
        compensated_sum(self.sequence.iter().map(|symbol| {
            // The invariant is the same as in molecular_weight: immutable checked input.
            let index = canonical_index(symbol.as_ascii())
                .expect("ProteinAnalysis retains only validated canonical residues");
            KYTE_DOOLITTLE[index]
        })) / count_as_float(self.len())
    }

    /// Calculates all five properties, using pH 7.0 for charge.
    #[must_use]
    pub fn properties(&self) -> ProteinProperties {
        self.properties_unchecked(7.0)
    }

    /// Calculates all five properties with explicit charge-estimation options.
    ///
    /// # Errors
    /// Returns [`ProteinPropertyError::InvalidPh`] for a nonfinite pH or a value
    /// outside inclusive 0..14. The retained source is never modified.
    pub fn properties_with(
        &self,
        options: ProteinPropertiesOptions,
    ) -> Result<ProteinProperties, ProteinPropertyError> {
        validate_ph(options.ph)?;
        Ok(self.properties_unchecked(options.ph))
    }

    fn properties_unchecked(&self, ph: f64) -> ProteinProperties {
        ProteinProperties {
            length: self.len(),
            molecular_weight: self.molecular_weight(),
            ph,
            charge: self.charge_unchecked(ph),
            pi_estimate: self.isoelectric_point_estimate(),
            aromaticity: self.aromaticity(),
            gravy: self.gravy(),
        }
    }

    fn charge_unchecked(&self, ph: f64) -> f64 {
        let positive = ordered_sum([
            basic_charge(1, self.n_terminal_pk, ph),
            basic_charge(self.counts[8], 10.0, ph),
            basic_charge(self.counts[14], 12.0, ph),
            basic_charge(self.counts[6], 5.98, ph),
        ]);
        let negative = ordered_sum([
            acidic_charge(1, self.c_terminal_pk, ph),
            acidic_charge(self.counts[2], 4.05, ph),
            acidic_charge(self.counts[3], 4.45, ph),
            acidic_charge(self.counts[1], 9.0, ph),
            acidic_charge(self.counts[19], 10.0, ph),
        ]);
        positive - negative
    }
}

/// Short property-analysis workflows for existing typed proteins.
pub trait ProteinPropertiesExt {
    /// Borrows and validates a canonical-20, nonempty peptide for repeated queries.
    ///
    /// # Errors
    /// Returns [`ProteinPropertyError::EmptySequence`] or the first
    /// [`ProteinPropertyError::UnsupportedResidue`].
    fn analyze_protein(&self) -> Result<ProteinAnalysis<'_>, ProteinPropertyError>;

    /// Calculates all five properties, with charge estimated at pH 7.0.
    ///
    /// # Errors
    /// Returns the same input-validation errors as [`Self::analyze_protein`].
    fn protein_properties(&self) -> Result<ProteinProperties, ProteinPropertyError> {
        Ok(self.analyze_protein()?.properties())
    }

    /// Calculates all five properties with an explicit charge-estimation pH.
    ///
    /// # Errors
    /// Validates input first, then rejects nonfinite pH and values outside 0..14
    /// with [`ProteinPropertyError::InvalidPh`].
    fn protein_properties_with(
        &self,
        options: ProteinPropertiesOptions,
    ) -> Result<ProteinProperties, ProteinPropertyError> {
        self.analyze_protein()?.properties_with(options)
    }
}

impl ProteinPropertiesExt for ProteinSequence {
    fn analyze_protein(&self) -> Result<ProteinAnalysis<'_>, ProteinPropertyError> {
        ProteinAnalysis::new(self)
    }
}

fn canonical_index(byte: u8) -> Option<usize> {
    match byte {
        b'A' => Some(0),
        b'C' => Some(1),
        b'D' => Some(2),
        b'E' => Some(3),
        b'F' => Some(4),
        b'G' => Some(5),
        b'H' => Some(6),
        b'I' => Some(7),
        b'K' => Some(8),
        b'L' => Some(9),
        b'M' => Some(10),
        b'N' => Some(11),
        b'P' => Some(12),
        b'Q' => Some(13),
        b'R' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'V' => Some(17),
        b'W' => Some(18),
        b'Y' => Some(19),
        _ => None,
    }
}

fn validate_ph(ph: f64) -> Result<(), ProteinPropertyError> {
    if !ph.is_finite() || !(0.0..=14.0).contains(&ph) {
        return Err(ProteinPropertyError::InvalidPh { ph });
    }
    Ok(())
}

// Physical-property models explicitly use binary64 counts. Values above 2^53
// undergo IEEE 754 rounding; exact integer sequence length/counts stay retained.
#[allow(clippy::cast_precision_loss)]
fn count_as_float(count: usize) -> f64 {
    count as f64
}

fn basic_charge(count: usize, pk: f64, ph: f64) -> f64 {
    let fraction = 1.0 / (10.0_f64.powf(ph - pk) + 1.0);
    count_as_float(count) * fraction
}

fn acidic_charge(count: usize, pk: f64, ph: f64) -> f64 {
    let fraction = 1.0 / (10.0_f64.powf(pk - ph) + 1.0);
    count_as_float(count) * fraction
}

fn ordered_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().fold(0.0, |sum, value| sum + value)
}

// Neumaier's compensated scalar reduction. Inputs and intermediate totals remain
// finite for all addressable canonical peptides and the supported pH interval.
fn compensated_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum: f64 = 0.0;
    let mut correction = 0.0;
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_same_bits(actual: f64, expected: f64) {
        assert_eq!(actual.to_bits(), expected.to_bits());
    }

    fn peptide(text: &str) -> ProteinSequence {
        ProteinSequence::from_ascii(text).unwrap()
    }

    #[test]
    fn rejects_empty_and_first_noncanonical_residue() {
        assert!(matches!(
            ProteinAnalysis::new(&peptide("")),
            Err(ProteinPropertyError::EmptySequence)
        ));
        for byte in b"BZJXUO" {
            let sequence = peptide(&format!("MA{}C", char::from(*byte)));
            assert!(matches!(
                ProteinAnalysis::new(&sequence),
                Err(ProteinPropertyError::UnsupportedResidue { offset: 2, byte: actual })
                    if actual == *byte
            ));
        }
    }

    #[test]
    fn rejects_invalid_ph_and_preserves_valid_endpoints() {
        let sequence = peptide("INGAR");
        let analysis = sequence.analyze_protein().unwrap();
        for ph in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.001, 14.001] {
            assert!(matches!(
                analysis.charge_at_ph(ph),
                Err(ProteinPropertyError::InvalidPh { .. })
            ));
            assert!(
                analysis
                    .properties_with(ProteinPropertiesOptions { ph })
                    .is_err()
            );
        }
        for ph in [-0.0, 0.0, 14.0] {
            assert!(analysis.charge_at_ph(ph).unwrap().is_finite());
        }
    }

    #[test]
    fn independent_mass_and_composition_examples() {
        let sequence = peptide("AGC");
        let analysis = sequence.analyze_protein().unwrap();
        assert!((analysis.molecular_weight() - 249.2874).abs() < 1e-10);
        assert_same_bits(analysis.aromaticity(), 0.0);
        assert_same_bits(analysis.gravy(), 1.3);
        let aromatic = peptide("FWY");
        let properties = aromatic.protein_properties().unwrap();
        assert_eq!(properties.len(), 3);
        assert!(!properties.is_empty());
        assert_same_bits(properties.aromaticity(), 1.0);
        assert!((properties.gravy() - 0.2).abs() < 1e-15);
    }

    #[test]
    fn published_rounded_charge_and_pi_examples() {
        let ingar = peptide("INGAR");
        let analysis = ingar.analyze_protein().unwrap();
        assert_eq!(format!("{:.2}", analysis.isoelectric_point()), "9.75");
        assert_eq!(
            format!("{:.2}", analysis.charge_at_ph(7.0).unwrap()),
            "0.76"
        );
        let peter = peptide("PETER");
        let analysis = peter.analyze_protein().unwrap();
        assert_eq!(format!("{:.2}", analysis.isoelectric_point()), "4.53");
        assert_eq!(
            format!("{:.2}", analysis.charge_at_ph(4.53).unwrap()),
            "0.00"
        );
    }

    #[test]
    fn reports_both_unbracketed_boundaries_without_changing_solver() {
        for (text, boundary) in [("DDDDDD", PiBoundary::Lower), ("RRRRRR", PiBoundary::Upper)] {
            let sequence = peptide(text);
            let analysis = sequence.analyze_protein().unwrap();
            let estimate = analysis.isoelectric_point_estimate();
            assert_eq!(estimate.boundary(), boundary);
            assert!(!estimate.bracketed());
            assert_eq!(
                estimate.value().to_bits(),
                analysis.isoelectric_point().to_bits()
            );
            assert!(analysis.charge_at_ph(estimate.value()).unwrap().abs() > 0.1);
        }
        let sequence = peptide("INGAR");
        assert!(
            sequence
                .analyze_protein()
                .unwrap()
                .isoelectric_point_estimate()
                .bracketed()
        );
    }

    #[test]
    fn charge_decreases_with_ph_for_a_mixed_peptide() {
        let sequence = peptide("ACDEFGHIKLMNPQRSTVWY");
        let analysis = sequence.analyze_protein().unwrap();
        let mut previous = analysis.charge_at_ph(0.0).unwrap();
        for step in 1..=140 {
            let current = analysis.charge_at_ph(f64::from(step) / 10.0).unwrap();
            assert!(current <= previous);
            previous = current;
        }
    }

    #[test]
    fn aggregate_and_individual_queries_agree() {
        let sequence = peptide("MKWVTFISLLFLFSSAYS");
        let analysis = sequence.analyze_protein().unwrap();
        assert!(!analysis.is_empty());
        let properties = sequence
            .protein_properties_with(ProteinPropertiesOptions { ph: 4.53 })
            .unwrap();
        assert_same_bits(properties.ph(), 4.53);
        assert_same_bits(properties.molecular_weight(), analysis.molecular_weight());
        assert_same_bits(properties.charge(), analysis.charge_at_ph(4.53).unwrap());
        assert_eq!(
            properties.pi_estimate(),
            analysis.isoelectric_point_estimate()
        );
        assert_same_bits(properties.aromaticity(), analysis.aromaticity());
        assert_same_bits(properties.gravy(), analysis.gravy());
        assert_eq!(
            analysis.properties(),
            sequence.protein_properties().unwrap()
        );
        assert_same_bits(analysis.properties().ph(), 7.0);
    }

    #[test]
    fn error_codes_are_stable() {
        assert_eq!(ProteinPropertyError::EmptySequence.code(), "empty_protein");
        assert_eq!(
            ProteinPropertyError::UnsupportedResidue {
                offset: 0,
                byte: b'X'
            }
            .code(),
            "unsupported_protein_residue"
        );
        assert_eq!(
            ProteinPropertyError::InvalidPh { ph: 15.0 }.code(),
            "invalid_ph"
        );
    }
}
