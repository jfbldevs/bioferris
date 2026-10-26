use core::fmt;

use align_core::{ProteinSequence, SequenceError};

// NCBI translation tables, transcribed from gc.prt version 4.6 (SHA-256
// 2aecbdd09ecc5e0d233c105090a593737298ef93d88f6bc0c1bef26df6051a04), downloaded on
// 2026-10-22 from https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt.
// Biological table data, not reference implementation source. Each codon
// position is ordered T, C, A, G. `amino_acids` is gc.prt's `ncbieaa`; `starts`
// is its `sncbieaa`, where `M` marks an initiator and `*` a codon that can end
// translation. In tables 27, 28 and 31 some `*` codons also code amino acids.
const STOP_BIT: u32 = 1 << 26;

struct CodeTable {
    id: u16,
    name: &'static str,
    amino_acids: &'static [u8; 64],
    starts: &'static [u8; 64],
}

const TABLES: [CodeTable; 27] = [
    CodeTable {
        id: 1,
        name: "Standard",
        amino_acids: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**--*----M---------------M----------------------------",
    },
    CodeTable {
        id: 2,
        name: "Vertebrate Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        starts: b"----------**--------------------MMMM----------**---M------------",
    },
    CodeTable {
        id: 3,
        name: "Yeast Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**----------------------MM---------------M------------",
    },
    CodeTable {
        id: 4,
        name: "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--MM------**-------M------------MMMM---------------M------------",
    },
    CodeTable {
        id: 5,
        name: "Invertebrate Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        starts: b"---M------**--------------------MMMM---------------M------------",
    },
    CodeTable {
        id: 6,
        name: "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
        amino_acids: b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    CodeTable {
        id: 9,
        name: "Echinoderm Mitochondrial; Flatworm Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M---------------M------------",
    },
    CodeTable {
        id: 10,
        name: "Euplotid Nuclear",
        amino_acids: b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M----------------------------",
    },
    CodeTable {
        id: 11,
        name: "Bacterial, Archaeal and Plant Plastid",
        amino_acids: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**--*----M------------MMMM---------------M------------",
    },
    CodeTable {
        id: 12,
        name: "Alternative Yeast Nuclear",
        amino_acids: b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*----M---------------M----------------------------",
    },
    CodeTable {
        id: 13,
        name: "Ascidian Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        starts: b"---M------**----------------------MM---------------M------------",
    },
    CodeTable {
        id: 14,
        name: "Alternative Flatworm Mitochondrial",
        amino_acids: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"-----------*-----------------------M----------------------------",
    },
    CodeTable {
        id: 15,
        name: "Blepharisma Macronuclear",
        amino_acids: b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------*---*--------------------M----------------------------",
    },
    CodeTable {
        id: 16,
        name: "Chlorophycean Mitochondrial",
        amino_acids: b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------*---*--------------------M----------------------------",
    },
    CodeTable {
        id: 21,
        name: "Trematode Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M---------------M------------",
    },
    CodeTable {
        id: 22,
        name: "Scenedesmus obliquus Mitochondrial",
        amino_acids: b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"------*---*---*--------------------M----------------------------",
    },
    CodeTable {
        id: 23,
        name: "Thraustochytrium Mitochondrial",
        amino_acids: b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--*-------**--*-----------------M--M---------------M------------",
    },
    CodeTable {
        id: 24,
        name: "Rhabdopleuridae Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M------**-------M---------------M---------------M------------",
    },
    CodeTable {
        id: 25,
        name: "Candidate Division SR1 and Gracilibacteria",
        amino_acids: b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**-----------------------M---------------M------------",
    },
    CodeTable {
        id: 26,
        name: "Pachysolen tannophilus Nuclear",
        amino_acids: b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*----M---------------M----------------------------",
    },
    CodeTable {
        id: 27,
        name: "Karyorelict Nuclear",
        amino_acids: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    CodeTable {
        id: 28,
        name: "Condylostoma Nuclear",
        amino_acids: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*--------------------M----------------------------",
    },
    CodeTable {
        id: 29,
        name: "Mesodinium Nuclear",
        amino_acids: b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    CodeTable {
        id: 30,
        name: "Peritrich Nuclear",
        amino_acids: b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    CodeTable {
        id: 31,
        name: "Blastocrithidia Nuclear",
        amino_acids: b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M----------------------------",
    },
    CodeTable {
        id: 32,
        name: "Balanophoraceae Plastid",
        amino_acids: b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------*---*----M------------MMMM---------------M------------",
    },
    CodeTable {
        id: 33,
        name: "Cephalodiscidae Mitochondrial",
        amino_acids: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M-------*-------M---------------M---------------M------------",
    },
];

/// Every current NCBI genetic code, with its amino-acid assignments, initiators
/// and stops.
///
/// Selecting a table does not infer an organism or apply context-dependent
/// recoding, RNA editing, biological frameshifts or incomplete-stop completion.
/// In tables 27, 28 and 31, the codons that can end translation also code amino
/// acids: they translate as amino acids, and only a complete CDS may end at one,
/// as its terminal codon. ORF discovery refuses those three tables.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[non_exhaustive]
pub enum GeneticCode {
    /// NCBI table 1, Standard: initiators TTG, CTG, ATG; stops TAA, TAG, TGA.
    ///
    /// This is the default.
    #[default]
    Standard,
    /// NCBI table 2, Vertebrate Mitochondrial: initiators ATT, ATC, ATA, ATG, GTG;
    /// stops TAA, TAG, AGA, AGG.
    VertebrateMitochondrial,
    /// NCBI table 3, Yeast Mitochondrial: initiators ATA, ATG, GTG; stops TAA, TAG.
    YeastMitochondrial,
    /// NCBI table 4, Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate
    /// Mitochondrial; Mycoplasma; Spiroplasma: initiators TTA, TTG, CTG, ATT, ATC, ATA,
    /// ATG, GTG; stops TAA, TAG.
    MoldProtozoanCoelenterateMitochondrial,
    /// NCBI table 5, Invertebrate Mitochondrial: initiators TTG, ATT, ATC, ATA, ATG,
    /// GTG; stops TAA, TAG.
    InvertebrateMitochondrial,
    /// NCBI table 6, Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear:
    /// initiators ATG; stops TGA.
    CiliateDasycladaceanHexamitaNuclear,
    /// NCBI table 9, Echinoderm Mitochondrial; Flatworm Mitochondrial: initiators ATG,
    /// GTG; stops TAA, TAG.
    EchinodermFlatwormMitochondrial,
    /// NCBI table 10, Euplotid Nuclear: initiators ATG; stops TAA, TAG.
    EuplotidNuclear,
    /// NCBI table 11, Bacterial, Archaeal and Plant Plastid: initiators TTG, CTG, ATT,
    /// ATC, ATA, ATG, GTG; stops TAA, TAG, TGA.
    BacterialArchaealPlastid,
    /// NCBI table 12, Alternative Yeast Nuclear: initiators CTG, ATG; stops TAA, TAG,
    /// TGA.
    AlternativeYeastNuclear,
    /// NCBI table 13, Ascidian Mitochondrial: initiators TTG, ATA, ATG, GTG; stops TAA,
    /// TAG.
    AscidianMitochondrial,
    /// NCBI table 14, Alternative Flatworm Mitochondrial: initiators ATG; stops TAG.
    AlternativeFlatwormMitochondrial,
    /// NCBI table 15, Blepharisma Macronuclear: initiators ATG; stops TAA, TGA.
    BlepharismaMacronuclear,
    /// NCBI table 16, Chlorophycean Mitochondrial: initiators ATG; stops TAA, TGA.
    ChlorophyceanMitochondrial,
    /// NCBI table 21, Trematode Mitochondrial: initiators ATG, GTG; stops TAA, TAG.
    TrematodeMitochondrial,
    /// NCBI table 22, Scenedesmus obliquus Mitochondrial: initiators ATG; stops TCA,
    /// TAA, TGA.
    ScenedesmusObliquusMitochondrial,
    /// NCBI table 23, Thraustochytrium Mitochondrial: initiators ATT, ATG, GTG; stops
    /// TTA, TAA, TAG, TGA.
    ThraustochytriumMitochondrial,
    /// NCBI table 24, Rhabdopleuridae Mitochondrial: initiators TTG, CTG, ATG, GTG;
    /// stops TAA, TAG.
    RhabdopleuridaeMitochondrial,
    /// NCBI table 25, Candidate Division SR1 and Gracilibacteria: initiators TTG, ATG,
    /// GTG; stops TAA, TAG.
    CandidateDivisionSr1Gracilibacteria,
    /// NCBI table 26, Pachysolen tannophilus Nuclear: initiators CTG, ATG; stops TAA,
    /// TAG, TGA.
    PachysolenTannophilusNuclear,
    /// NCBI table 27, Karyorelict Nuclear: initiators ATG; TGA code amino acids and can
    /// also end translation.
    KaryorelictNuclear,
    /// NCBI table 28, Condylostoma Nuclear: initiators ATG; TAA, TAG, TGA code amino
    /// acids and can also end translation.
    CondylostomaNuclear,
    /// NCBI table 29, Mesodinium Nuclear: initiators ATG; stops TGA.
    MesodiniumNuclear,
    /// NCBI table 30, Peritrich Nuclear: initiators ATG; stops TGA.
    PeritrichNuclear,
    /// NCBI table 31, Blastocrithidia Nuclear: initiators ATG; TAA, TAG code amino
    /// acids and can also end translation.
    BlastocrithidiaNuclear,
    /// NCBI table 32, Balanophoraceae Plastid: initiators TTG, CTG, ATT, ATC, ATA, ATG,
    /// GTG; stops TAA, TGA.
    BalanophoraceaePlastid,
    /// NCBI table 33, Cephalodiscidae Mitochondrial: initiators TTG, CTG, ATG, GTG;
    /// stops TAG.
    CephalodiscidaeMitochondrial,
}

impl GeneticCode {
    /// Every supported table, in ascending NCBI identifier order.
    pub const ALL: [Self; 27] = [
        Self::Standard,
        Self::VertebrateMitochondrial,
        Self::YeastMitochondrial,
        Self::MoldProtozoanCoelenterateMitochondrial,
        Self::InvertebrateMitochondrial,
        Self::CiliateDasycladaceanHexamitaNuclear,
        Self::EchinodermFlatwormMitochondrial,
        Self::EuplotidNuclear,
        Self::BacterialArchaealPlastid,
        Self::AlternativeYeastNuclear,
        Self::AscidianMitochondrial,
        Self::AlternativeFlatwormMitochondrial,
        Self::BlepharismaMacronuclear,
        Self::ChlorophyceanMitochondrial,
        Self::TrematodeMitochondrial,
        Self::ScenedesmusObliquusMitochondrial,
        Self::ThraustochytriumMitochondrial,
        Self::RhabdopleuridaeMitochondrial,
        Self::CandidateDivisionSr1Gracilibacteria,
        Self::PachysolenTannophilusNuclear,
        Self::KaryorelictNuclear,
        Self::CondylostomaNuclear,
        Self::MesodiniumNuclear,
        Self::PeritrichNuclear,
        Self::BlastocrithidiaNuclear,
        Self::BalanophoraceaePlastid,
        Self::CephalodiscidaeMitochondrial,
    ];

    const fn table(self) -> &'static CodeTable {
        // Variants are declared in the order of `TABLES`.
        &TABLES[self as usize]
    }

    /// Returns the NCBI translation-table identifier.
    #[must_use]
    pub const fn id(self) -> u16 {
        self.table().id
    }

    /// Returns the NCBI name of the selected table.
    #[must_use]
    pub const fn name(self) -> &'static str {
        self.table().name
    }

    /// Returns whether some codons both code an amino acid and can end
    /// translation, which holds for NCBI tables 27, 28 and 31 only.
    #[must_use]
    pub fn has_conditional_stops(self) -> bool {
        let table = self.table();
        table
            .starts
            .iter()
            .zip(table.amino_acids)
            .any(|(start, residue)| *start == b'*' && *residue != b'*')
    }

    const fn assignments(self) -> &'static [u8; 64] {
        self.table().amino_acids
    }

    const fn is_start(self, index: usize) -> bool {
        self.table().starts[index] == b'M'
    }

    const fn can_end_translation(self, index: usize) -> bool {
        self.table().starts[index] == b'*'
    }
}

impl TryFrom<u16> for GeneticCode {
    type Error = TranslationError;

    /// Selects an NCBI table by identifier without inferring a fallback.
    ///
    /// # Errors
    /// Returns [`TranslationError::UnsupportedGeneticCode`] for identifiers that
    /// are not current NCBI tables, such as 0, 7, 8, 17 to 20 and 34 or more.
    fn try_from(id: u16) -> Result<Self, Self::Error> {
        Self::ALL
            .into_iter()
            .find(|code| code.id() == id)
            .ok_or(TranslationError::UnsupportedGeneticCode { id })
    }
}

/// A zero-based forward offset into the original, fully validated input.
///
/// This does not select a strand. Reverse-complement translation is an explicit
/// separate operation. Offsets beyond the input length are rejected at execution.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum ReadingFrame {
    /// Start with the first nucleotide, skipping no leading bases.
    #[default]
    Zero,
    /// Skip exactly one leading base.
    One,
    /// Skip exactly two leading bases.
    Two,
}

impl ReadingFrame {
    /// Returns the leading-base count: 0, 1 or 2.
    #[must_use]
    pub const fn offset(self) -> usize {
        match self {
            Self::Zero => 0,
            Self::One => 1,
            Self::Two => 2,
        }
    }
}

impl TryFrom<usize> for ReadingFrame {
    type Error = TranslationError;

    /// Selects a forward offset without clipping or modulo conversion.
    ///
    /// # Errors
    /// Returns [`TranslationError::InvalidReadingFrame`] for values above two.
    fn try_from(frame: usize) -> Result<Self, Self::Error> {
        match frame {
            0 => Ok(Self::Zero),
            1 => Ok(Self::One),
            2 => Ok(Self::Two),
            _ => Err(TranslationError::InvalidReadingFrame { frame }),
        }
    }
}

/// Stop handling for a selected genetic-code translation.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum TranslationMode {
    /// Translate all complete codons, retaining definite stops as `*`.
    #[default]
    Full,
    /// Stop before the first definite stop; do not include its `*`.
    ToStop,
    /// Require a complete CDS; recode the start to M and omit the terminal stop.
    Cds,
}

/// Handling of one or two trailing bases outside a complete codon.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum PartialCodonPolicy {
    /// Reject rather than silently discard biological input (default).
    #[default]
    Reject,
    /// Ignore trailing bases and report their count in the result.
    ///
    /// This never overrides complete-CDS length validation.
    Ignore,
}

/// Explicit policies for table 1 translation; no other genetic code is selected.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct TranslationOptions {
    /// Full translation, stop at first definite stop, or complete CDS.
    pub mode: TranslationMode,
    /// Reject incomplete codons by default, or explicitly ignore trailing bases.
    pub partial_codons: PartialCodonPolicy,
}

/// Explicit genetic code, forward frame and translation policies.
///
/// Defaults are table 1, frame zero, full translation and rejected trailing
/// partial codons. Existing [`TranslationOptions`] and standard-code methods
/// retain their original meaning; conversion from those options selects these
/// same code/frame defaults.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct TranslationConfig {
    /// Explicit NCBI genetic code; no organism or annotation is inferred.
    pub genetic_code: GeneticCode,
    /// Strict zero-based forward offset into the validated input.
    pub reading_frame: ReadingFrame,
    /// Full translation, first definite stop or complete CDS validation.
    pub mode: TranslationMode,
    /// Handling of incomplete trailing codons in the selected suffix.
    pub partial_codons: PartialCodonPolicy,
}

impl From<TranslationOptions> for TranslationConfig {
    fn from(options: TranslationOptions) -> Self {
        Self {
            genetic_code: GeneticCode::Standard,
            reading_frame: ReadingFrame::Zero,
            mode: options.mode,
            partial_codons: options.partial_codons,
        }
    }
}

/// A translation result, possibly containing `*` stop markers.
///
/// Unlike the historical `ProteinSequence`, this is a translation product, not
/// necessarily an ungapped stop-free protein. Ambiguous residues use B, Z, J or X.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Translation {
    peptide: String,
    ignored_trailing_bases: usize,
    genetic_code: GeneticCode,
    reading_frame: ReadingFrame,
}

impl Translation {
    /// Returns the translated ASCII residues and any retained `*` markers.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.peptide
    }

    /// Returns the number of translated residues, counting retained stops.
    #[must_use]
    pub fn len(&self) -> usize {
        self.peptide.len()
    }

    /// Returns whether no residues or stops were emitted.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.peptide.is_empty()
    }

    /// Returns the selected suffix's incomplete trailing-base count (0, 1 or 2).
    ///
    /// This does not count complete codons beyond an early `ToStop` termination.
    #[must_use]
    pub const fn ignored_trailing_bases(&self) -> usize {
        self.ignored_trailing_bases
    }

    /// Returns the genetic-code assignments and initiation set used.
    #[must_use]
    pub const fn genetic_code(&self) -> GeneticCode {
        self.genetic_code
    }

    /// Returns the selected forward reading frame.
    #[must_use]
    pub const fn reading_frame(&self) -> ReadingFrame {
        self.reading_frame
    }

    /// Returns the leading bases skipped by the explicit frame selection.
    ///
    /// This is independent of ignored incomplete trailing bases and always
    /// equals the selected frame's offset for a successful translation.
    #[must_use]
    pub const fn skipped_leading_bases(&self) -> usize {
        self.reading_frame.offset()
    }

    /// Validates the product as a historical stop-free protein for alignment.
    ///
    /// # Errors
    ///
    /// Returns [`SequenceError`] at the first `*`; no stops are silently stripped.
    pub fn to_protein(&self) -> Result<ProteinSequence, SequenceError> {
        ProteinSequence::from_ascii(self.peptide.as_bytes())
    }
}

impl fmt::Display for Translation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.peptide)
    }
}

/// A structurally invalid translation request on an already validated sequence.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum TranslationError {
    /// One or two trailing bases do not form a complete codon.
    PartialCodon {
        /// Nucleotide count in the selected suffix, after the frame offset.
        length: usize,
        /// Number of trailing bases outside complete codons.
        remainder: usize,
    },
    /// A CDS needs at least one start codon and a distinct terminal stop codon.
    CdsTooShort {
        /// Nucleotide count in the selected suffix, after the frame offset.
        length: usize,
    },
    /// The initial codon is not universally a table 1 start.
    InvalidStartCodon {
        /// Original DNA or RNA codon.
        codon: [u8; 3],
    },
    /// The final codon is not universally a stop in the selected code.
    MissingTerminalStop {
        /// Original DNA or RNA codon.
        codon: [u8; 3],
    },
    /// A complete CDS contains a definite stop before its terminal codon.
    InternalStop {
        /// Zero-based stop-codon offset in the original input, before slicing.
        offset: usize,
    },
    /// The requested identifier is not a current NCBI translation table.
    UnsupportedGeneticCode {
        /// Requested NCBI translation-table identifier.
        id: u16,
    },
    /// A numeric selector is not one of the forward frames 0, 1 and 2.
    InvalidReadingFrame {
        /// Requested zero-based forward offset.
        frame: usize,
    },
    /// A valid frame selector lies beyond the actual sequence length.
    FrameOutOfBounds {
        /// Requested zero-based forward offset.
        frame: usize,
        /// Original, fully validated nucleotide count.
        length: usize,
    },
    /// The first codon is not universally an initiator in a nonstandard code.
    ///
    /// Standard-code failures retain the historical [`Self::InvalidStartCodon`]
    /// variant and display wording.
    InvalidStartCodonForCode {
        /// Original DNA or RNA codon in the selected suffix.
        codon: [u8; 3],
        /// Nonstandard genetic code whose initiation set was used.
        genetic_code: GeneticCode,
    },
}

impl fmt::Display for TranslationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::PartialCodon { length, remainder } => write!(
                formatter,
                "sequence length {length} leaves {remainder} incomplete codon bases"
            ),
            Self::CdsTooShort { length } => {
                write!(
                    formatter,
                    "CDS length {length} needs distinct start and stop codons"
                )
            }
            Self::InvalidStartCodon { codon } => {
                write!(
                    formatter,
                    "invalid table 1 start codon {}",
                    String::from_utf8_lossy(codon)
                )
            }
            Self::MissingTerminalStop { codon } => {
                write!(
                    formatter,
                    "expected terminal stop, got {}",
                    String::from_utf8_lossy(codon)
                )
            }
            Self::InternalStop { offset } => {
                write!(
                    formatter,
                    "CDS contains an internal stop at nucleotide {offset}"
                )
            }
            Self::UnsupportedGeneticCode { id } => write!(
                formatter,
                "unsupported genetic code {id}; supported NCBI tables are 1-6, 9-16 and 21-33"
            ),
            Self::InvalidReadingFrame { frame } => {
                write!(
                    formatter,
                    "invalid reading frame {frame}; expected 0, 1 or 2"
                )
            }
            Self::FrameOutOfBounds { frame, length } => write!(
                formatter,
                "reading frame {frame} exceeds sequence length {length}"
            ),
            Self::InvalidStartCodonForCode {
                codon,
                genetic_code,
            } => write!(
                formatter,
                "invalid table {} start codon {}",
                genetic_code.id(),
                String::from_utf8_lossy(codon)
            ),
        }
    }
}

impl std::error::Error for TranslationError {}

pub(crate) fn translate_standard(
    bytes: &[u8],
    options: TranslationOptions,
) -> Result<Translation, TranslationError> {
    translate_configured(bytes, options.into())
}

pub(crate) fn translate_configured(
    input: &[u8],
    config: TranslationConfig,
) -> Result<Translation, TranslationError> {
    let frame = config.reading_frame.offset();
    let bytes = input
        .get(frame..)
        .ok_or(TranslationError::FrameOutOfBounds {
            frame,
            length: input.len(),
        })?;
    let remainder = bytes.len() % 3;
    let cds = config.mode == TranslationMode::Cds;
    if remainder != 0 && (cds || config.partial_codons == PartialCodonPolicy::Reject) {
        return Err(TranslationError::PartialCodon {
            length: bytes.len(),
            remainder,
        });
    }
    if cds {
        if bytes.len() < 6 {
            return Err(TranslationError::CdsTooShort {
                length: bytes.len(),
            });
        }
        let first = [bytes[0], bytes[1], bytes[2]];
        if !codon_outcome(first, config.genetic_code).all_start {
            return Err(if config.genetic_code == GeneticCode::Standard {
                TranslationError::InvalidStartCodon { codon: first }
            } else {
                TranslationError::InvalidStartCodonForCode {
                    codon: first,
                    genetic_code: config.genetic_code,
                }
            });
        }
        let end = bytes.len();
        let last = [bytes[end - 3], bytes[end - 2], bytes[end - 1]];
        if !codon_outcome(last, config.genetic_code).all_end_translation {
            return Err(TranslationError::MissingTerminalStop { codon: last });
        }
    }

    let mut peptide = String::with_capacity(bytes.len() / 3);
    let codon_count = bytes.len() / 3;
    for (index, codon) in bytes.chunks_exact(3).enumerate() {
        if cds && index + 1 == codon_count {
            break;
        }
        let residue = if cds && index == 0 {
            b'M'
        } else {
            codon_outcome([codon[0], codon[1], codon[2]], config.genetic_code).residue
        };
        if residue == b'*' {
            if cds {
                return Err(TranslationError::InternalStop {
                    offset: frame + index * 3,
                });
            }
            if config.mode == TranslationMode::ToStop {
                break;
            }
        }
        peptide.push(char::from(residue));
    }
    Ok(Translation {
        peptide,
        ignored_trailing_bases: remainder,
        genetic_code: config.genetic_code,
        reading_frame: config.reading_frame,
    })
}

pub(crate) struct CodonOutcome {
    pub(crate) residue: u8,
    pub(crate) all_start: bool,
    /// Every concrete codon can end translation: a stop, or in tables 27, 28
    /// and 31 also a codon that otherwise codes an amino acid.
    pub(crate) all_end_translation: bool,
}

// Bits follow NCBI's T,C,A,G order. This function is called only with validated
// IUPAC symbols; U and T denote the same concrete base for the genetic code.
fn base_mask(base: u8) -> u8 {
    match base {
        b'T' | b'U' => 1,
        b'C' => 2,
        b'A' => 4,
        b'G' => 8,
        b'R' => 12,
        b'Y' => 3,
        b'S' => 10,
        b'W' => 5,
        b'K' => 9,
        b'M' => 6,
        b'B' => 11,
        b'D' => 13,
        b'H' => 7,
        b'V' => 14,
        b'N' => 15,
        _ => 0,
    }
}

const fn amino_bit(residue: u8) -> u32 {
    1 << (residue - b'A')
}

pub(crate) fn codon_outcome(codon: [u8; 3], genetic_code: GeneticCode) -> CodonOutcome {
    let masks = codon.map(base_mask);
    let mut outcomes = 0_u32;
    let mut all_start = true;
    let mut all_end_translation = true;
    for first in 0..4 {
        for second in 0..4 {
            for third in 0..4 {
                if masks[0] & (1 << first) == 0
                    || masks[1] & (1 << second) == 0
                    || masks[2] & (1 << third) == 0
                {
                    continue;
                }
                let index = first * 16 + second * 4 + third;
                let residue = genetic_code.assignments()[index];
                outcomes |= if residue == b'*' {
                    STOP_BIT
                } else {
                    amino_bit(residue)
                };
                all_start &= genetic_code.is_start(index);
                all_end_translation &= genetic_code.can_end_translation(index);
            }
        }
    }
    CodonOutcome {
        residue: resolve_outcomes(outcomes),
        all_start,
        all_end_translation,
    }
}

fn resolve_outcomes(outcomes: u32) -> u8 {
    if outcomes == STOP_BIT {
        return b'*';
    }
    for residue in b'A'..=b'Z' {
        if outcomes == amino_bit(residue) {
            return residue;
        }
    }
    if outcomes == amino_bit(b'D') | amino_bit(b'N') {
        b'B'
    } else if outcomes == amino_bit(b'E') | amino_bit(b'Q') {
        b'Z'
    } else if outcomes == amino_bit(b'I') | amino_bit(b'L') {
        b'J'
    } else {
        b'X'
    }
}
