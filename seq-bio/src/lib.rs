//! `BioFerris`: Rust-native bioinformatics, from sequences to workflows.
//!
//! Typed IUPAC DNA/RNA translation supports every current NCBI code with explicit
//! forward frames, six-frame workflows and concrete-base ORF discovery, not the
//! complete Biopython API. Input types normalize case and reject gaps. Single-frame
//! translation rejects partial codons by default; six-frame translation reports
//! excluded leading/trailing bases explicitly.
//!
//! ```
//! use bioferris::prelude::*;
//!
//! let dna = DnaSequence::try_from("ATGGCTTAA")?;
//! assert_eq!(dna.reverse_complement()?.to_string(), "TTAAGCCAT");
//! let rna = dna.transcribe()?;
//! let translation = rna.translate_standard(TranslationOptions {
//!     mode: TranslationMode::Cds,
//!     ..TranslationOptions::default()
//! })?;
//! let protein = translation.to_protein()?;
//! assert_eq!(protein.to_string(), "MA");
//! let run = Pairwise::gotoh().run(&protein, &protein)?;
//! assert_eq!(run.result().first_range(), 0..2);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
#![forbid(unsafe_code)]
#![deny(nonstandard_style)]

mod alignment_analysis;
mod motifs;
mod operations;
mod protein_properties;
mod publication;
mod reading_frames;
mod record;
mod sequence;
mod statistics;
mod translation;

pub use alignment_analysis::{
    AlignmentAmbiguityPolicy, AlignmentAnalysis, AlignmentAnalysisAlphabet, AlignmentAnalysisError,
    AlignmentAnalysisLimits, AlignmentAnalysisRow, AlignmentAnnotationPolicy, AlignmentCell,
    AlignmentCellCounts, AlignmentColumnStatistics, AlignmentConsensus, AlignmentGapRuns,
    AlignmentPairCounts, AlignmentRowStatistics, AlignmentSlice, AlignmentStatistics,
    AlignmentStatisticsOptions, ConsensusColumn, ConsensusGapPolicy, ConsensusMissingPolicy,
    ConsensusOptions, ConsensusStatus, ConsensusTiePolicy, IdentityDenominator, RowCoordinates,
};
pub use io::annotated::{
    CdsReferencePolicy, CdsReferenceStatus, CdsTranslation, CdsTranslationError,
    CdsTranslationOptions,
};
pub use motifs::{
    DnaMotif, IupacMotif, MatrixSearchOptions, MotifAlphabet, MotifCounts, MotifError,
    MotifInputAmbiguity, MotifMatch, MotifScores, MotifSearchOptions, MotifWeights, RnaMotif,
    ScoredMotifMatch, dna_motif, rna_motif,
};
pub use operations::{MatchOverlap, SequenceExt, SequenceMatches, SequenceOperationError};
pub use protein_properties::{
    IsoelectricPointEstimate, PiBoundary, ProteinAnalysis, ProteinProperties, ProteinPropertiesExt,
    ProteinPropertiesOptions, ProteinPropertyError,
};
pub use reading_frames::{
    FrameTranslation, OpenReadingFrame, OrfError, OrfOptions, OrfStartPolicy, OrfStopPolicy,
    ReadingFramesExt, Strand, StrandSelection,
};
pub use record::Record;

/// The unchanged alignment facade, including pairwise, multiple, and graph APIs.
pub use align;
pub use align_core::{
    Dna, DnaSequence, Protein, ProteinSequence, Rna, RnaSequence, Sequence, SequenceError,
};
pub use sequence::{
    ConfiguredTranslationExt, DnaSequenceExt, NucleotideSequenceExt, RnaSequenceExt,
};
pub use statistics::{
    GcContent, GcPolicy, KmerAmbiguity, KmerCounts, KmerOptions, SequenceStatistics,
    StatisticsAlphabet, StatisticsError, StatisticsExt, StatisticsOptions, summarize_records,
    summarize_sequences,
};
pub use translation::{
    GeneticCode, PartialCodonPolicy, ReadingFrame, Translation, TranslationConfig,
    TranslationError, TranslationMode, TranslationOptions,
};

/// Parses typed IUPAC DNA, normalizing lowercase to uppercase.
///
/// # Errors
///
/// Returns [`SequenceError`] with the first invalid byte and its zero-based offset.
pub fn dna(input: impl AsRef<[u8]>) -> Result<DnaSequence, SequenceError> {
    DnaSequence::from_ascii(input)
}

/// Parses typed IUPAC RNA, normalizing lowercase to uppercase.
///
/// # Errors
///
/// Returns [`SequenceError`] for invalid symbols, including T instead of U.
pub fn rna(input: impl AsRef<[u8]>) -> Result<RnaSequence, SequenceError> {
    RnaSequence::from_ascii(input)
}

/// Parses a typed stop-free protein for the existing alignment algorithms.
///
/// # Errors
///
/// Returns [`SequenceError`] for invalid symbols, including `*` and gaps.
pub fn protein(input: impl AsRef<[u8]>) -> Result<ProteinSequence, SequenceError> {
    ProteinSequence::from_ascii(input)
}

/// File readers, validated records and the unchanged legacy format interfaces.
pub mod io;

/// Typed sequences, new biological operations, and the existing alignment API.
pub mod prelude {
    pub use crate::{
        AlignmentAmbiguityPolicy, AlignmentAnalysis, AlignmentAnalysisAlphabet,
        AlignmentAnalysisError, AlignmentAnalysisLimits, AlignmentAnalysisRow,
        AlignmentAnnotationPolicy, AlignmentCell, AlignmentCellCounts, AlignmentColumnStatistics,
        AlignmentConsensus, AlignmentGapRuns, AlignmentPairCounts, AlignmentRowStatistics,
        AlignmentSlice, AlignmentStatistics, AlignmentStatisticsOptions, ConsensusColumn,
        ConsensusGapPolicy, ConsensusMissingPolicy, ConsensusOptions, ConsensusStatus,
        ConsensusTiePolicy, DnaMotif, IdentityDenominator, IupacMotif, MatrixSearchOptions,
        MotifAlphabet, MotifCounts, MotifError, MotifInputAmbiguity, MotifMatch, MotifScores,
        MotifSearchOptions, MotifWeights, RnaMotif, RowCoordinates, ScoredMotifMatch, dna_motif,
        rna_motif,
    };
    pub use align::prelude::{
        AdaptiveAffineBackend, AdaptiveAffineDecision, AdaptiveAffineWorkspace, AffineGap,
        AlignRecordsExt, AlignedRow, AlignmentOperation, Alphabet, ClustalWrite, ColorMode,
        CoordinateStyle, DenseSubstitutionMatrix, Dna, DnaSequence, FastaWrite, FreeEndGaps,
        GotohMode, Graph, LinearAlgorithm, LinearGap, Mafft7526Alignment, Mafft7526EinsiAligner,
        Mafft7526LinsiAligner, Mafft7526Profile, MafftStrategy, MatchMismatch, Multiple,
        MultipleAlignment, MultipleAlignmentExt, MultipleAlignmentView, MultipleRun, Muscle5Mode,
        MyersDistanceRun, MyersMatch, MyersPattern, MyersSearchRun, MyersTsvWrite, NexusAlphabet,
        NexusWrite, Pairwise, PairwiseAlignment, PairwiseAlignmentExt, PairwiseAlignmentView,
        PairwiseRun, PairwiseValueRun, PhylipLayout, PhylipNameStyle, PhylipOptions, PhylipWrite,
        PoaGraph, PoaGraphWrite, PoastaAligner, PoastaQueryRun, Protein, ProteinSequence,
        RecordsAligner, ReportOptions, Rna, RnaSequence, Score, Sequence, SequenceRecord,
        StockholmWrite, SubstitutionScorer, Symbol, Verbosity, XDropDirection, XDropSeed,
        parse_fasta, read_fasta, write_clustal, write_fasta, write_myers_tsv, write_nexus,
        write_pairwise_fasta, write_phylip, write_phylip_with_options, write_poa_dot,
        write_poa_gfa, write_stockholm, write_text_report,
    };
    // Do not shadow Rust's generic Result with an alignment-only error alias
    // in workflows that also perform translation and file I/O.
    pub use align::{Error as AlignmentError, Result as AlignmentResult};

    pub use crate::{
        CdsReferencePolicy, CdsReferenceStatus, CdsTranslation, CdsTranslationError,
        CdsTranslationOptions, ConfiguredTranslationExt, DnaSequenceExt, FrameTranslation,
        GcContent, GcPolicy, GeneticCode, IsoelectricPointEstimate, KmerAmbiguity, KmerCounts,
        KmerOptions, MatchOverlap, NucleotideSequenceExt, OpenReadingFrame, OrfError, OrfOptions,
        OrfStartPolicy, OrfStopPolicy, PartialCodonPolicy, PiBoundary, ProteinAnalysis,
        ProteinProperties, ProteinPropertiesExt, ProteinPropertiesOptions, ProteinPropertyError,
        ReadingFrame, ReadingFramesExt, Record, RnaSequenceExt, SequenceExt, SequenceMatches,
        SequenceOperationError, SequenceStatistics, StatisticsAlphabet, StatisticsError,
        StatisticsExt, StatisticsOptions, Strand, StrandSelection, Translation, TranslationConfig,
        TranslationError, TranslationMode, TranslationOptions, dna, protein, rna,
        summarize_records, summarize_sequences,
    };
}
