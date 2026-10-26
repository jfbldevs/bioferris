use core::{fmt, ops::Range};

use align_core::{DnaSequence, ProteinSequence, RnaSequence, SequenceError};

use crate::{GeneticCode, ReadingFrame, sequence::complement_byte, translation::codon_outcome};

const FRAMES: [ReadingFrame; 3] = [ReadingFrame::Zero, ReadingFrame::One, ReadingFrame::Two];

/// Translation direction relative to the original nucleotide input.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Strand {
    /// Read the original input from left to right.
    Forward,
    /// Read the reverse complement from left to right.
    Reverse,
}

/// Strand selection for ORF discovery; each selected strand uses all three frames.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum StrandSelection {
    /// Forward frames first, then reverse-complement frames (default).
    #[default]
    Both,
    /// Only the original input's three frames.
    Forward,
    /// Only the reverse complement's three frames.
    Reverse,
}

/// How an ORF candidate begins within a stop-delimited interval.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum OrfStartPolicy {
    /// Earliest ATG/AUG; recode its first residue to M (default).
    #[default]
    Atg,
    /// Earliest initiator from the selected NCBI code; recode its first residue to M.
    GeneticCode,
    /// First non-stop codon after a stop or frame boundary; do not recode it.
    AnySense,
}

/// Whether an ORF candidate must have a complete terminal stop codon.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum OrfStopPolicy {
    /// Emit only candidates terminated by a stop (default).
    #[default]
    Require,
    /// Also emit a final unterminated candidate, excluding incomplete codons.
    AllowPartial,
}

/// Explicit biological policies and result-count budget for ORF discovery.
///
/// Only concrete DNA/RNA is accepted, even in skipped or incomplete regions.
/// Each stop-delimited interval retains only its earliest accepted start;
/// nested starts are not enumerated. Candidates are not gene predictions.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OrfOptions {
    /// NCBI code 1, 2 or 11; default is code 1.
    pub genetic_code: GeneticCode,
    /// Selected strands; default is both.
    pub strands: StrandSelection,
    /// Initiation rule; default is ATG/AUG only.
    pub start_policy: OrfStartPolicy,
    /// Termination rule; default requires a stop.
    pub stop_policy: OrfStopPolicy,
    /// Inclusive minimum residue count, excluding stop; must be positive (default 1).
    pub min_amino_acids: usize,
    /// Maximum retained ORFs (default 100,000); excess results fail, never truncate.
    ///
    /// Zero is valid and permits only an empty result. This bounds the result
    /// count, not total process memory; input and peptide storage remain linear.
    pub max_orfs: usize,
}

impl Default for OrfOptions {
    fn default() -> Self {
        Self {
            genetic_code: GeneticCode::Standard,
            strands: StrandSelection::Both,
            start_policy: OrfStartPolicy::Atg,
            stop_policy: OrfStopPolicy::Require,
            min_amino_acids: 1,
            max_orfs: 100_000,
        }
    }
}

/// Complete-codon translation of one oriented frame, including `*` stop markers.
///
/// All ranges refer to the original input. The frame offset is measured in the
/// oriented sequence, so it is not the genomic coordinate modulo three.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FrameTranslation {
    strand: Strand,
    reading_frame: ReadingFrame,
    genetic_code: GeneticCode,
    range: Range<usize>,
    peptide: String,
    skipped_leading_bases: usize,
    ignored_trailing_bases: usize,
}

impl FrameTranslation {
    /// Returns the translation direction.
    #[must_use]
    pub const fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns offset 0/1/2 in the original or reverse-complement sequence.
    #[must_use]
    pub const fn reading_frame(&self) -> ReadingFrame {
        self.reading_frame
    }

    /// Returns the NCBI code used for every complete codon.
    #[must_use]
    pub const fn genetic_code(&self) -> GeneticCode {
        self.genetic_code
    }

    /// Returns the complete-codon span as zero-based half-open original coordinates.
    ///
    /// Empty forward frames are anchored at the actual skipped-base count;
    /// reverse frames map that empty boundary back into the original input.
    #[must_use]
    pub fn range(&self) -> Range<usize> {
        self.range.clone()
    }

    /// Returns translated residues, retaining definite stops and IUPAC B/Z/J/X.
    #[must_use]
    pub fn peptide(&self) -> &str {
        &self.peptide
    }

    /// Returns actual leading bases skipped, capped at the input length.
    #[must_use]
    pub const fn skipped_leading_bases(&self) -> usize {
        self.skipped_leading_bases
    }

    /// Returns excluded incomplete trailing bases (zero, one or two).
    #[must_use]
    pub const fn ignored_trailing_bases(&self) -> usize {
        self.ignored_trailing_bases
    }
}

/// One nonempty ORF candidate, with its terminal stop stored separately.
///
/// The coding range excludes the stop. Peptide order follows translation
/// direction even though all ranges use ascending original-input coordinates.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OpenReadingFrame {
    strand: Strand,
    reading_frame: ReadingFrame,
    genetic_code: GeneticCode,
    range: Range<usize>,
    stop_range: Option<Range<usize>>,
    peptide: String,
}

impl OpenReadingFrame {
    /// Returns the translation direction.
    #[must_use]
    pub const fn strand(&self) -> Strand {
        self.strand
    }

    /// Returns offset 0/1/2 in the selected oriented sequence.
    #[must_use]
    pub const fn reading_frame(&self) -> ReadingFrame {
        self.reading_frame
    }

    /// Returns the NCBI code used to identify and translate this candidate.
    #[must_use]
    pub const fn genetic_code(&self) -> GeneticCode {
        self.genetic_code
    }

    /// Returns zero-based half-open original-input coordinates, excluding stop.
    #[must_use]
    pub fn range(&self) -> Range<usize> {
        self.range.clone()
    }

    /// Returns the terminal stop's original-input range, or `None` at an open end.
    #[must_use]
    pub fn stop_range(&self) -> Option<Range<usize>> {
        self.stop_range.clone()
    }

    /// Returns a stop-free peptide; start-based policies recode its first residue to M.
    #[must_use]
    pub fn peptide(&self) -> &str {
        &self.peptide
    }

    /// Constructs the existing validated protein type for downstream alignment.
    ///
    /// # Errors
    /// Propagates [`SequenceError`] from protein validation. Scanner-produced
    /// peptides contain only concrete amino acids and therefore validate.
    pub fn to_protein(&self) -> Result<ProteinSequence, SequenceError> {
        ProteinSequence::from_ascii(self.peptide.as_bytes())
    }
}

/// A rejected ORF policy/input or an exceeded result budget; no partial results escape.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum OrfError {
    /// A zero minimum length is not meaningful for the nonempty-candidate contract.
    InvalidMinimumLength,
    /// The first IUPAC ambiguous base in the original input, before strand selection.
    AmbiguousBase {
        /// Zero-based original-input nucleotide offset.
        offset: usize,
        /// Uppercase ASCII ambiguity symbol.
        byte: u8,
    },
    /// Another qualifying candidate would exceed the requested result budget.
    ResultLimitExceeded {
        /// Caller-specified maximum number of ORFs.
        limit: usize,
    },
    /// The selected code has no codon that always ends translation, so ORF
    /// boundaries are undefined: NCBI tables 27, 28 and 31.
    ConditionalStops {
        /// The selected genetic code.
        genetic_code: GeneticCode,
    },
}

impl fmt::Display for OrfError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidMinimumLength => {
                formatter.write_str("minimum ORF length must be at least one amino acid")
            }
            Self::AmbiguousBase { offset, byte } => write!(
                formatter,
                "ORF scanning requires concrete bases; found {} at nucleotide {offset}",
                char::from(*byte)
            ),
            Self::ResultLimitExceeded { limit } => {
                write!(formatter, "ORF result count exceeds limit {limit}")
            }
            Self::ConditionalStops { genetic_code } => write!(
                formatter,
                "genetic code {} has no codon that always ends translation, so ORF boundaries are undefined",
                genetic_code.id()
            ),
        }
    }
}

impl std::error::Error for OrfError {}

/// Six-frame translation and explicitly configured ORF discovery on typed DNA/RNA.
///
/// Existing single-frame translation defaults remain unchanged. Six-frame
/// translation excludes and reports incomplete codons; ORF scanning accepts only
/// unambiguous bases. These synchronous workflows retain their returned products.
///
/// ```
/// use bioferris::prelude::*;
/// let sequence = dna("ATGAAATAA")?;
/// let frames = sequence.translate_six_frames();
/// assert_eq!(frames[0].peptide(), "MK*");
/// let orfs = sequence.find_orfs()?;
/// assert_eq!(orfs.len(), 1);
/// assert_eq!(orfs[0].range(), 0..6);
/// assert_eq!(orfs[0].stop_range(), Some(6..9));
/// assert_eq!(orfs[0].peptide(), "MK");
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait ReadingFramesExt {
    /// Translates all six frames using NCBI table 1, retaining stops.
    ///
    /// Order is forward 0/1/2 then reverse 0/1/2. Empty/short inputs still produce
    /// six entries; actual skipped-leading and ignored-trailing counts are explicit.
    fn translate_six_frames(&self) -> [FrameTranslation; 6] {
        self.translate_six_frames_with(GeneticCode::Standard)
    }

    /// Translates all six frames with the selected code and no initiation recoding.
    fn translate_six_frames_with(&self, genetic_code: GeneticCode) -> [FrameTranslation; 6];

    /// Finds earliest-ATG/AUG, stop-terminated candidates in both strands, using code 1.
    ///
    /// # Errors
    /// Returns [`OrfError`] for ambiguous input or more than 100,000 results.
    fn find_orfs(&self) -> Result<Vec<OpenReadingFrame>, OrfError> {
        self.find_orfs_with(OrfOptions::default())
    }

    /// Finds candidates under explicit initiation, termination, strand and count policies.
    ///
    /// Order is forward frames 0/1/2 then reverse frames 0/1/2, each in translation
    /// direction. Only the first accepted start per stop-delimited interval is
    /// retained. The main range excludes stop; no circular wrapping or annotation
    /// interpretation is performed. Minimum length excludes stop and is inclusive.
    ///
    /// # Errors
    /// Returns [`OrfError::InvalidMinimumLength`] before scanning for a zero minimum,
    /// [`OrfError::ConditionalStops`] before scanning with NCBI table 27, 28 or 31,
    /// [`OrfError::AmbiguousBase`] at the first ambiguous original-input symbol,
    /// or [`OrfError::ResultLimitExceeded`] without returning a truncated result.
    fn find_orfs_with(&self, options: OrfOptions) -> Result<Vec<OpenReadingFrame>, OrfError>;
}

fn source_range(range: Range<usize>, length: usize, strand: Strand) -> Range<usize> {
    match strand {
        Strand::Forward => range,
        Strand::Reverse => length - range.end..length - range.start,
    }
}

fn translate_codons(bytes: &[u8], code: GeneticCode, recode_start: bool) -> String {
    let mut peptide = String::with_capacity(bytes.len() / 3);
    for (index, codon) in bytes.chunks_exact(3).enumerate() {
        let residue = if index == 0 && recode_start {
            b'M'
        } else {
            codon_outcome([codon[0], codon[1], codon[2]], code).residue
        };
        peptide.push(char::from(residue));
    }
    peptide
}

fn reverse_bytes(bytes: &[u8], adenine_complement: u8) -> Vec<u8> {
    bytes
        .iter()
        .rev()
        .map(|&byte| complement_byte(byte, adenine_complement))
        .collect()
}

fn six_frames(bytes: &[u8], adenine_complement: u8, code: GeneticCode) -> [FrameTranslation; 6] {
    let reverse = reverse_bytes(bytes, adenine_complement);
    std::array::from_fn(|index| {
        let (oriented, strand) = if index < 3 {
            (bytes, Strand::Forward)
        } else {
            (reverse.as_slice(), Strand::Reverse)
        };
        let frame = FRAMES[index % 3];
        let start = frame.offset().min(oriented.len());
        let trailing = (oriented.len() - start) % 3;
        let end = oriented.len() - trailing;
        FrameTranslation {
            strand,
            reading_frame: frame,
            genetic_code: code,
            range: source_range(start..end, oriented.len(), strand),
            peptide: translate_codons(&oriented[start..end], code, false),
            skipped_leading_bases: start,
            ignored_trailing_bases: trailing,
        }
    })
}

struct OrfScanner<'a> {
    bytes: &'a [u8],
    strand: Strand,
    options: OrfOptions,
    results: &'a mut Vec<OpenReadingFrame>,
}

impl OrfScanner<'_> {
    fn scan(&mut self) -> Result<(), OrfError> {
        for frame in FRAMES {
            self.scan_frame(frame)?;
        }
        Ok(())
    }

    fn scan_frame(&mut self, frame: ReadingFrame) -> Result<(), OrfError> {
        let bytes = self.bytes;
        let offset = frame.offset().min(bytes.len());
        let mut candidate = None;
        for (index, triplet) in bytes[offset..].chunks_exact(3).enumerate() {
            let position = offset + index * 3;
            let codon = [triplet[0], triplet[1], triplet[2]];
            let outcome = codon_outcome(codon, self.options.genetic_code);
            if outcome.residue == b'*' {
                if let Some(start) = candidate.take() {
                    self.emit(frame, start, position, true)?;
                }
            } else if candidate.is_none()
                && match self.options.start_policy {
                    OrfStartPolicy::Atg => codon == *b"ATG" || codon == *b"AUG",
                    OrfStartPolicy::GeneticCode => outcome.all_start,
                    OrfStartPolicy::AnySense => true,
                }
            {
                candidate = Some(position);
            }
        }
        if self.options.stop_policy == OrfStopPolicy::AllowPartial {
            if let Some(start) = candidate {
                let end = bytes.len() - (bytes.len() - offset) % 3;
                self.emit(frame, start, end, false)?;
            }
        }
        Ok(())
    }

    fn emit(
        &mut self,
        frame: ReadingFrame,
        start: usize,
        end: usize,
        has_stop: bool,
    ) -> Result<(), OrfError> {
        if (end - start) / 3 < self.options.min_amino_acids {
            return Ok(());
        }
        if self.results.len() >= self.options.max_orfs {
            return Err(OrfError::ResultLimitExceeded {
                limit: self.options.max_orfs,
            });
        }
        self.results.push(OpenReadingFrame {
            strand: self.strand,
            reading_frame: frame,
            genetic_code: self.options.genetic_code,
            range: source_range(start..end, self.bytes.len(), self.strand),
            stop_range: has_stop.then(|| source_range(end..end + 3, self.bytes.len(), self.strand)),
            peptide: translate_codons(
                &self.bytes[start..end],
                self.options.genetic_code,
                self.options.start_policy != OrfStartPolicy::AnySense,
            ),
        });
        Ok(())
    }
}

fn scan_orfs(
    bytes: &[u8],
    adenine_complement: u8,
    options: OrfOptions,
) -> Result<Vec<OpenReadingFrame>, OrfError> {
    if options.min_amino_acids == 0 {
        return Err(OrfError::InvalidMinimumLength);
    }
    if options.genetic_code.has_conditional_stops() {
        return Err(OrfError::ConditionalStops {
            genetic_code: options.genetic_code,
        });
    }
    if let Some((offset, &byte)) = bytes
        .iter()
        .enumerate()
        .find(|(_, byte)| !matches!(byte, b'A' | b'C' | b'G' | b'T' | b'U'))
    {
        return Err(OrfError::AmbiguousBase { offset, byte });
    }
    let mut results = Vec::new();
    if options.strands != StrandSelection::Reverse {
        OrfScanner {
            bytes,
            strand: Strand::Forward,
            options,
            results: &mut results,
        }
        .scan()?;
    }
    if options.strands != StrandSelection::Forward {
        let reverse = reverse_bytes(bytes, adenine_complement);
        OrfScanner {
            bytes: &reverse,
            strand: Strand::Reverse,
            options,
            results: &mut results,
        }
        .scan()?;
    }
    Ok(results)
}

macro_rules! reading_frame_operations {
    ($sequence:ty, $adenine_complement:expr) => {
        impl ReadingFramesExt for $sequence {
            fn translate_six_frames_with(
                &self,
                genetic_code: GeneticCode,
            ) -> [FrameTranslation; 6] {
                let bytes: Vec<_> = self.iter().map(|symbol| symbol.as_ascii()).collect();
                six_frames(&bytes, $adenine_complement, genetic_code)
            }

            fn find_orfs_with(
                &self,
                options: OrfOptions,
            ) -> Result<Vec<OpenReadingFrame>, OrfError> {
                let bytes: Vec<_> = self.iter().map(|symbol| symbol.as_ascii()).collect();
                scan_orfs(&bytes, $adenine_complement, options)
            }
        }
    };
}

reading_frame_operations!(DnaSequence, b'T');
reading_frame_operations!(RnaSequence, b'U');
