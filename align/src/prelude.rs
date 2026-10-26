//! Common imports for the ergonomic facade.

pub use align_core::{
    AffineGap, AlignmentOperation, Alphabet, DenseSubstitutionMatrix, Dna, DnaSequence, LinearGap,
    MatchMismatch, PairwiseAlignment, Protein, ProteinSequence, Rna, RnaSequence, Score, Sequence,
    SubstitutionScorer, Symbol,
};
pub use align_fast::{MyersMatch, MyersPattern, XDropDirection, XDropSeed};
pub use align_formats::{
    ClustalWrite, FastaWrite, MyersTsvWrite, NexusAlphabet, NexusWrite, PhylipLayout,
    PhylipNameStyle, PhylipOptions, PhylipWrite, PoaGraphWrite, StockholmWrite, parse_fasta,
    read_fasta, write_clustal, write_fasta, write_myers_tsv, write_nexus, write_pairwise_fasta,
    write_phylip, write_phylip_with_options, write_poa_dot, write_poa_gfa, write_stockholm,
    write_text_report,
};
pub use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord};
pub use align_pairwise::FreeEndGaps;
pub use align_poa::PoaGraph;

pub use crate::graph::{Graph, PoastaAligner};
pub use crate::multiple::{
    AlignRecordsExt, Mafft7526EinsiAligner, Mafft7526LinsiAligner, MafftStrategy, Multiple,
    Muscle5Mode, RecordsAligner,
};
pub use crate::pairwise::{
    AdaptiveAffineBackend, AdaptiveAffineDecision, AdaptiveAffineWorkspace, GotohMode,
    LinearAlgorithm, Pairwise,
};
pub use crate::report::{
    ColorMode, CoordinateStyle, MultipleRun, MyersDistanceRun, MyersSearchRun, PairwiseRun,
    PairwiseValueRun, PoastaQueryRun, ReportOptions, Verbosity,
};
pub use crate::result::{
    MultipleAlignmentExt, MultipleAlignmentView, PairwiseAlignmentExt, PairwiseAlignmentView,
};
pub use crate::{Error, Mafft7526Alignment, Mafft7526Profile, Result};
