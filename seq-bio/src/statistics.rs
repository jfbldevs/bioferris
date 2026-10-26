//! Count-first sequence summaries and explicitly bounded forward k-mers.

use std::collections::BTreeMap;
use std::fmt;

use align_core::{Alphabet, Dna, Protein, Rna, Sequence};
use serde::{Deserialize, Serialize};

use crate::Record;
use crate::io::BioIoError;

/// How nucleotide ambiguity contributes to GC content.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum GcPolicy {
    /// Count G/C/S over A/C/G/T/U/S/W, excluding other ambiguous symbols.
    #[default]
    Remove,
    /// Count G/C/S over every symbol, without resolving ambiguous symbols.
    Ignore,
    /// Weight every IUPAC symbol by its equally likely concrete expansions.
    Weighted,
}

/// How an ambiguous symbol affects a forward k-mer window.
#[derive(Clone, Copy, Debug, Default, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum KmerAmbiguity {
    /// Skip windows containing nucleotide ambiguity or protein B/Z/J/X.
    #[default]
    Skip,
    /// Keep ambiguity letters as literal symbols, without expanding them.
    Literal,
}

/// Explicit k-mer parameters and aggregate allocation/work budgets.
///
/// Budgets apply across all records, including skipped windows. The key-byte
/// budget covers stored key text, not map nodes or the reusable O(k) scratch key.
/// Distinct-key count separately bounds map-node overhead. Zero budgets are valid.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KmerOptions {
    /// Forward overlapping window length; operations require 1..=64.
    pub k: usize,
    /// Skip ambiguous windows by default, or retain their literal symbols.
    pub ambiguity: KmerAmbiguity,
    /// Maximum considered windows, including skipped windows.
    ///
    /// Each record's full window count is checked before scanning its symbols;
    /// an error reports the attempted aggregate including that entire record.
    pub max_windows: u64,
    /// Maximum distinct stored keys.
    pub max_distinct: usize,
    /// Maximum aggregate bytes of stored key text.
    pub max_key_bytes: usize,
}

impl KmerOptions {
    /// Selects `k`, skip ambiguity, one million windows, 100,000 keys and 8 MiB.
    ///
    /// Construction does not validate `k`; the operation rejects an invalid
    /// length before consuming input, including an empty input iterator.
    #[must_use]
    pub const fn new(k: usize) -> Self {
        Self {
            k,
            ambiguity: KmerAmbiguity::Skip,
            max_windows: 1_000_000,
            max_distinct: 100_000,
            max_key_bytes: 8 * 1024 * 1024,
        }
    }

    fn validate(self) -> Result<(), StatisticsError> {
        if !(1..=64).contains(&self.k) {
            return Err(StatisticsError::InvalidOptions(format!(
                "k-mer length must be in 1..=64 (received {})",
                self.k
            )));
        }
        Ok(())
    }
}

/// Summary policies; defaults calculate remove-policy GC without k-mers.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StatisticsOptions {
    /// GC ambiguity policy; it has no effect for proteins, whose GC is absent.
    pub gc_policy: GcPolicy,
    /// Optional bounded k-mer counting, disabled by default.
    pub kmers: Option<KmerOptions>,
}

impl StatisticsOptions {
    /// Validates configuration without opening or consuming sequence input.
    ///
    /// K-mer lengths must be in 1..=64. All budgets, including zero, are valid
    /// configurations; their enforcement depends on the actual input. Summary
    /// operations call this automatically before consuming their iterators.
    ///
    /// # Errors
    /// Returns [`StatisticsError::InvalidOptions`] for an invalid k-mer length.
    pub fn validate(self) -> Result<(), StatisticsError> {
        self.kmers.map_or(Ok(()), KmerOptions::validate)
    }
}

/// Aggregate nucleotide GC content, derived after counting all input symbols.
///
/// This is a readout, not validated input to a calculation. Counts converted to
/// floating point can lose integer precision above 2^53.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct GcContent {
    /// Policy used for the numerator and denominator.
    pub policy: GcPolicy,
    /// GC-equivalent count, including fractional contributions when weighted.
    pub gc_equivalent: f64,
    /// Number of symbols retained in the denominator by the selected policy.
    pub denominator: u64,
    /// GC-equivalent count divided by denominator, or zero for denominator zero.
    pub fraction: f64,
}

/// Deterministically ordered counts of overlapping, forward, per-record k-mers.
///
/// Keys never cross record boundaries or undergo reverse-complement
/// canonicalization. This is a readout, not validated calculation input.
#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct KmerCounts {
    /// Length of each counted key.
    pub k: usize,
    /// Ambiguity handling applied to every considered window.
    pub ambiguity: KmerAmbiguity,
    /// All considered windows, including skipped windows.
    pub windows: u64,
    /// Windows added to the key counts.
    pub counted: u64,
    /// Windows skipped because of ambiguity.
    pub skipped: u64,
    /// Uppercase literal keys in lexicographic order, with positive counts.
    pub counts: BTreeMap<String, u64>,
}

/// Aggregate sequence statistics without retained input sequences or records.
///
/// Public fields are readouts, not accepted back as computational input. Empty
/// input has no length extrema; one empty sequence has extrema of zero. Protein
/// GC is absent. A summary is returned only after every input succeeds.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct SequenceStatistics {
    /// Number of sequences, including empty sequences.
    pub records: u64,
    /// Total symbol count across all sequences.
    pub length: u64,
    /// Shortest input sequence length, or `None` for no sequences.
    pub minimum_length: Option<usize>,
    /// Longest input sequence length, or `None` for no sequences.
    pub maximum_length: Option<usize>,
    /// Positive counts of observed uppercase symbols, sorted by symbol.
    pub composition: BTreeMap<char, u64>,
    /// Non-ACGT/U nucleotides, including S/W, or protein B/Z/J/X symbols.
    pub ambiguous: u64,
    /// Aggregate nucleotide GC, or `None` for a protein alphabet.
    pub gc: Option<GcContent>,
    /// Optional bounded k-mer counts.
    pub kmers: Option<KmerCounts>,
}

/// A summary failed; no successful partial result is returned.
#[derive(Debug)]
#[non_exhaustive]
pub enum StatisticsError {
    /// A parameter is outside the supported statistics contract.
    InvalidOptions(String),
    /// An explicit aggregate work or stored-key budget would be exceeded.
    LimitExceeded {
        /// Stable resource: `kmer_windows`, `kmer_distinct` or `kmer_key_bytes`.
        resource: &'static str,
        /// Configured inclusive limit.
        limit: u64,
        /// Attempted aggregate value before the rejected work or insertion.
        observed: u64,
    },
    /// A count, length, or stored-key byte total cannot be represented.
    CountOverflow,
    /// A streaming record reader failed; its context and cause are retained.
    Input(BioIoError),
}

impl StatisticsError {
    /// Returns a stable machine code, preserving a wrapped reader's code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::InvalidOptions(_) => "invalid_statistics_options",
            Self::LimitExceeded { .. } => "limit_exceeded",
            Self::CountOverflow => "count_overflow",
            Self::Input(error) => error.code(),
        }
    }
}

impl fmt::Display for StatisticsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOptions(message) => formatter.write_str(message),
            Self::LimitExceeded {
                resource,
                limit,
                observed,
            } => write!(
                formatter,
                "{resource} limit {limit} exceeded: observed {observed}"
            ),
            Self::CountOverflow => formatter.write_str("sequence statistics count overflow"),
            Self::Input(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for StatisticsError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Input(error) => Some(error),
            _ => None,
        }
    }
}

impl From<BioIoError> for StatisticsError {
    fn from(error: BioIoError) -> Self {
        Self::Input(error)
    }
}

mod sealed {
    pub trait Sealed {
        const NUCLEOTIDE: bool;
        fn is_ambiguous(symbol: u8) -> bool;
    }
}

/// A supported statistics alphabet: DNA, RNA or stop-free extended protein.
///
/// This trait is sealed so custom alphabets cannot silently inherit unsuitable
/// ambiguity or GC semantics. It does not change the existing alignment types.
pub trait StatisticsAlphabet: Alphabet + sealed::Sealed {}

impl StatisticsAlphabet for Dna {}
impl StatisticsAlphabet for Rna {}
impl StatisticsAlphabet for Protein {}

impl sealed::Sealed for Dna {
    const NUCLEOTIDE: bool = true;

    fn is_ambiguous(symbol: u8) -> bool {
        !matches!(symbol, b'A' | b'C' | b'G' | b'T')
    }
}

impl sealed::Sealed for Rna {
    const NUCLEOTIDE: bool = true;

    fn is_ambiguous(symbol: u8) -> bool {
        !matches!(symbol, b'A' | b'C' | b'G' | b'U')
    }
}

impl sealed::Sealed for Protein {
    const NUCLEOTIDE: bool = false;

    fn is_ambiguous(symbol: u8) -> bool {
        matches!(symbol, b'B' | b'Z' | b'J' | b'X')
    }
}

/// Count-first composition, GC and optional bounded k-mer sequence summaries.
///
/// ```
/// use bioferris::prelude::*;
/// let summary = dna("ACGTSWN")?.statistics()?;
/// assert_eq!(summary.length, 7);
/// assert_eq!(summary.ambiguous, 3);
/// assert_eq!(summary.gc.as_ref().map(|gc| gc.fraction), Some(0.5));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait StatisticsExt {
    /// Summarizes one sequence with remove-policy GC and no k-mer calculation.
    ///
    /// # Errors
    /// Returns [`StatisticsError::CountOverflow`] if a count cannot be represented.
    fn statistics(&self) -> Result<SequenceStatistics, StatisticsError> {
        self.statistics_with(StatisticsOptions::default())
    }

    /// Summarizes one sequence with explicit GC policy and optional bounded k-mers.
    ///
    /// # Errors
    /// Rejects invalid k-mer options, exceeded budgets and count overflow. An
    /// allocator failure is not a recoverable statistics error.
    fn statistics_with(
        &self,
        options: StatisticsOptions,
    ) -> Result<SequenceStatistics, StatisticsError>;
}

impl<A: StatisticsAlphabet> StatisticsExt for Sequence<A> {
    fn statistics_with(
        &self,
        options: StatisticsOptions,
    ) -> Result<SequenceStatistics, StatisticsError> {
        summarize_sequences(std::iter::once(self), options)
    }
}

/// Summarizes borrowed sequences without collecting or copying them.
///
/// K-mer windows never cross sequence boundaries. Fractions are calculated from
/// aggregate symbol counts, not averaged sequence fractions. Invalid options
/// are rejected before consuming the iterator. No partial result is returned.
///
/// # Errors
/// Returns [`StatisticsError`] for invalid options, exceeded aggregate budgets
/// or count overflow. Stops at the first error without requesting further input.
pub fn summarize_sequences<'a, A: StatisticsAlphabet>(
    sequences: impl IntoIterator<Item = &'a Sequence<A>>,
    options: StatisticsOptions,
) -> Result<SequenceStatistics, StatisticsError> {
    let mut accumulator = Accumulator::new(options)?;
    for sequence in sequences {
        accumulator.add(sequence)?;
    }
    accumulator.finish::<A>()
}

/// Summarizes a streaming record reader without retaining records or qualities.
///
/// Uses the validated sequence, not source case, headers or quality scores.
/// Every reader result is consumed on success, allowing late parse/decompression
/// failures to surface. Invalid options are rejected before consuming input.
///
/// # Errors
/// Returns [`StatisticsError::Input`] with the original reader error and context,
/// or a statistics option/budget/count error. The first failure stops consumption
/// and drops the partial accumulator; no successful partial summary is returned.
pub fn summarize_records<A: StatisticsAlphabet>(
    records: impl IntoIterator<Item = Result<Record<A>, BioIoError>>,
    options: StatisticsOptions,
) -> Result<SequenceStatistics, StatisticsError> {
    let mut accumulator = Accumulator::new(options)?;
    for record in records {
        accumulator.add(record?.sequence())?;
    }
    accumulator.finish::<A>()
}

struct Accumulator {
    summary: SequenceStatistics,
    gc_policy: GcPolicy,
    kmers: Option<KmerAccumulator>,
}

impl Accumulator {
    fn new(options: StatisticsOptions) -> Result<Self, StatisticsError> {
        options.validate()?;
        let kmers = options.kmers.map(KmerAccumulator::new).transpose()?;
        Ok(Self {
            summary: SequenceStatistics {
                records: 0,
                length: 0,
                minimum_length: None,
                maximum_length: None,
                composition: BTreeMap::new(),
                ambiguous: 0,
                gc: None,
                kmers: None,
            },
            gc_policy: options.gc_policy,
            kmers,
        })
    }

    fn add<A: StatisticsAlphabet>(
        &mut self,
        sequence: &Sequence<A>,
    ) -> Result<(), StatisticsError> {
        let length = sequence.len();
        self.summary.records = checked_count(self.summary.records, 1)?;
        self.summary.length = checked_count(self.summary.length, count_from_usize(length)?)?;
        if let Some(kmers) = &mut self.kmers {
            kmers.add(sequence)?;
        }
        self.summary.minimum_length = Some(
            self.summary
                .minimum_length
                .map_or(length, |old| old.min(length)),
        );
        self.summary.maximum_length = Some(
            self.summary
                .maximum_length
                .map_or(length, |old| old.max(length)),
        );
        for symbol in sequence.iter() {
            let byte = symbol.as_ascii();
            let count = self
                .summary
                .composition
                .entry(char::from(byte))
                .or_default();
            *count = checked_count(*count, 1)?;
            if A::is_ambiguous(byte) {
                self.summary.ambiguous = checked_count(self.summary.ambiguous, 1)?;
            }
        }
        Ok(())
    }

    fn finish<A: StatisticsAlphabet>(mut self) -> Result<SequenceStatistics, StatisticsError> {
        if A::NUCLEOTIDE {
            self.summary.gc = Some(gc_content(&self.summary, self.gc_policy)?);
        }
        self.summary.kmers = self.kmers.map(|kmers| kmers.result);
        Ok(self.summary)
    }
}

struct KmerAccumulator {
    result: KmerCounts,
    options: KmerOptions,
    key_bytes: usize,
    scratch: String,
}

impl KmerAccumulator {
    fn new(options: KmerOptions) -> Result<Self, StatisticsError> {
        options.validate()?;
        Ok(Self {
            result: KmerCounts {
                k: options.k,
                ambiguity: options.ambiguity,
                windows: 0,
                counted: 0,
                skipped: 0,
                counts: BTreeMap::new(),
            },
            options,
            key_bytes: 0,
            scratch: String::new(),
        })
    }

    fn add<A: StatisticsAlphabet>(
        &mut self,
        sequence: &Sequence<A>,
    ) -> Result<(), StatisticsError> {
        let windows = match sequence.len().checked_sub(self.options.k) {
            Some(difference) => checked_count(count_from_usize(difference)?, 1)?,
            None => 0,
        };
        let total_windows = checked_count(self.result.windows, windows)?;
        check_limit("kmer_windows", self.options.max_windows, total_windows)?;
        self.result.windows = total_windows;
        for window in sequence.as_slice().windows(self.options.k) {
            if self.options.ambiguity == KmerAmbiguity::Skip
                && window
                    .iter()
                    .any(|symbol| A::is_ambiguous(symbol.as_ascii()))
            {
                self.result.skipped = checked_count(self.result.skipped, 1)?;
                continue;
            }
            self.scratch.clear();
            for symbol in window {
                self.scratch.push(char::from(symbol.as_ascii()));
            }
            self.count_key()?;
        }
        Ok(())
    }

    fn count_key(&mut self) -> Result<(), StatisticsError> {
        if let Some(count) = self.result.counts.get_mut(self.scratch.as_str()) {
            *count = checked_count(*count, 1)?;
        } else {
            let distinct = self
                .result
                .counts
                .len()
                .checked_add(1)
                .ok_or(StatisticsError::CountOverflow)?;
            let key_bytes = self
                .key_bytes
                .checked_add(self.options.k)
                .ok_or(StatisticsError::CountOverflow)?;
            check_limit(
                "kmer_distinct",
                count_from_usize(self.options.max_distinct)?,
                count_from_usize(distinct)?,
            )?;
            check_limit(
                "kmer_key_bytes",
                count_from_usize(self.options.max_key_bytes)?,
                count_from_usize(key_bytes)?,
            )?;
            self.result.counts.insert(self.scratch.clone(), 1);
            self.key_bytes = key_bytes;
        }
        self.result.counted = checked_count(self.result.counted, 1)?;
        Ok(())
    }
}

fn checked_count(current: u64, added: u64) -> Result<u64, StatisticsError> {
    current
        .checked_add(added)
        .ok_or(StatisticsError::CountOverflow)
}

fn count_from_usize(value: usize) -> Result<u64, StatisticsError> {
    u64::try_from(value).map_err(|_| StatisticsError::CountOverflow)
}

fn check_limit(resource: &'static str, limit: u64, observed: u64) -> Result<(), StatisticsError> {
    if observed > limit {
        Err(StatisticsError::LimitExceeded {
            resource,
            limit,
            observed,
        })
    } else {
        Ok(())
    }
}

fn gc_content(
    summary: &SequenceStatistics,
    policy: GcPolicy,
) -> Result<GcContent, StatisticsError> {
    let mut denominator = 0;
    let mut gc_equivalent = 0.0;
    for (&symbol, &count) in &summary.composition {
        if policy == GcPolicy::Remove && !matches!(symbol, 'A' | 'C' | 'G' | 'T' | 'U' | 'S' | 'W')
        {
            continue;
        }
        denominator = checked_count(denominator, count)?;
        let weight = if policy == GcPolicy::Weighted {
            match symbol {
                'G' | 'C' | 'S' => 1.0,
                'B' | 'V' => 2.0 / 3.0,
                'D' | 'H' => 1.0 / 3.0,
                'R' | 'Y' | 'K' | 'M' | 'N' => 0.5,
                _ => 0.0,
            }
        } else {
            f64::from(u8::from(matches!(symbol, 'G' | 'C' | 'S')))
        };
        gc_equivalent += count_as_f64(count) * weight;
    }
    let fraction = if denominator == 0 {
        0.0
    } else {
        gc_equivalent / count_as_f64(denominator)
    };
    Ok(GcContent {
        policy,
        gc_equivalent,
        denominator,
        fraction,
    })
}

// Counts remain exact u64 values. Only documented floating GC readouts can lose
// integer precision above 2^53; this conversion does not affect count arithmetic.
#[allow(clippy::cast_precision_loss)]
fn count_as_f64(value: u64) -> f64 {
    value as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn configuration_validation_needs_no_sequence_input() {
        assert!(StatisticsOptions::default().validate().is_ok());
        for k in [1, 64] {
            let options = StatisticsOptions {
                gc_policy: GcPolicy::Weighted,
                kmers: Some(KmerOptions {
                    max_windows: 0,
                    max_distinct: 0,
                    max_key_bytes: 0,
                    ..KmerOptions::new(k)
                }),
            };
            assert!(options.validate().is_ok());
        }
        for k in [0, 65, usize::MAX] {
            let options = StatisticsOptions {
                kmers: Some(KmerOptions::new(k)),
                ..StatisticsOptions::default()
            };
            assert!(matches!(
                options.validate(),
                Err(StatisticsError::InvalidOptions(_))
            ));
        }
    }

    #[test]
    fn count_overflow_is_checked_before_wrapping() {
        assert!(matches!(
            checked_count(u64::MAX, 1),
            Err(StatisticsError::CountOverflow)
        ));
        let sequence = Sequence::<Dna>::from_ascii("A").unwrap();
        let mut records = Accumulator::new(StatisticsOptions::default()).unwrap();
        records.summary.records = u64::MAX;
        assert!(matches!(
            records.add(&sequence),
            Err(StatisticsError::CountOverflow)
        ));
        let mut length = Accumulator::new(StatisticsOptions::default()).unwrap();
        length.summary.length = u64::MAX;
        assert!(matches!(
            length.add(&sequence),
            Err(StatisticsError::CountOverflow)
        ));
        let mut composition = Accumulator::new(StatisticsOptions::default()).unwrap();
        composition.summary.composition.insert('A', u64::MAX);
        assert!(matches!(
            composition.add(&sequence),
            Err(StatisticsError::CountOverflow)
        ));
        let mut ambiguous = Accumulator::new(StatisticsOptions::default()).unwrap();
        ambiguous.summary.ambiguous = u64::MAX;
        assert!(matches!(
            ambiguous.add(&Sequence::<Dna>::from_ascii("N").unwrap()),
            Err(StatisticsError::CountOverflow)
        ));
    }

    #[test]
    fn kmer_count_overflow_is_checked_independently() {
        let sequence = Sequence::<Dna>::from_ascii("A").unwrap();
        let options = KmerOptions {
            max_windows: u64::MAX,
            ..KmerOptions::new(1)
        };
        let mut windows = KmerAccumulator::new(options).unwrap();
        windows.result.windows = u64::MAX;
        assert!(matches!(
            windows.add(&sequence),
            Err(StatisticsError::CountOverflow)
        ));
        let mut counts = KmerAccumulator::new(options).unwrap();
        counts.result.counts.insert("A".to_owned(), u64::MAX);
        assert!(matches!(
            counts.add(&sequence),
            Err(StatisticsError::CountOverflow)
        ));
        let mut counted = KmerAccumulator::new(options).unwrap();
        counted.result.counted = u64::MAX;
        assert!(matches!(
            counted.add(&sequence),
            Err(StatisticsError::CountOverflow)
        ));
        let mut skipped = KmerAccumulator::new(options).unwrap();
        skipped.result.skipped = u64::MAX;
        assert!(matches!(
            skipped.add(&Sequence::<Dna>::from_ascii("N").unwrap()),
            Err(StatisticsError::CountOverflow)
        ));
    }

    #[test]
    fn stored_key_budget_is_checked_before_insertion() {
        let options = KmerOptions {
            max_key_bytes: 1,
            ..KmerOptions::new(2)
        };
        let mut kmers = KmerAccumulator::new(options).unwrap();
        let sequence = Sequence::<Dna>::from_ascii("AC").unwrap();
        assert!(matches!(
            kmers.add(&sequence),
            Err(StatisticsError::LimitExceeded {
                resource: "kmer_key_bytes",
                limit: 1,
                observed: 2
            })
        ));
        assert!(kmers.result.counts.is_empty());
        assert_eq!(kmers.key_bytes, 0);
    }
}
