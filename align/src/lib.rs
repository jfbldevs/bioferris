//! A small, Rust-native facade over the validated alignment kernels.
//!
//! Start with [`prelude`] and select an algorithm through [`Pairwise`],
//! [`Multiple`], or [`Graph`]. The facade preserves the scientific contracts of
//! the underlying crates while hiding workspaces and implementation backends.
//! Completed runs can be written directly as human text, versioned JSON, standard
//! alignment matrices, tabular search results, or graph files.
//!
//! ```
//! use align::prelude::*;
//!
//! let target = DnaSequence::try_from("GATTACA")?;
//! let query = DnaSequence::try_from("GCATGCT")?;
//! let run = Pairwise::gotoh()
//!     .scores(2, -1)
//!     .gaps(-5, -1)
//!     .run(&target, &query)?;
//!
//! assert_eq!(run.result().first_range(), 0..target.len());
//! let report = format!("{}", run.report().labels("Target", "Query"));
//! assert!(report.contains("Algorithm     : Gotoh"));
//! assert!(report.contains("Identity"));
//! # Ok::<(), align::Error>(())
//! ```
#![forbid(unsafe_code)]

pub mod error;
pub mod graph;
pub mod multiple;
pub mod output;
pub mod pairwise;
pub mod prelude;
pub mod report;
pub mod result;

/// Direct access to the detailed kernel APIs and diagnostics.
pub mod advanced {
    pub use align_clustal as clustal;
    pub use align_core as core;
    pub use align_fast as fast;
    pub use align_formats as formats;
    pub use align_mafft as mafft;
    pub use align_mafft_compat as mafft_compat;
    pub use align_msa as msa;
    pub use align_muscle as muscle;
    pub use align_pairwise as pairwise;
    pub use align_poa as poa;
}

pub use align_mafft_compat::{Mafft7526Alignment, Mafft7526Profile};
pub use error::{Error, Result};
pub use graph::{Graph, PoastaAligner};
pub use multiple::{
    AlignRecordsExt, ClustalOmegaAligner, Mafft7526EinsiAligner, Mafft7526LinsiAligner,
    MafftAligner, MafftResult, MafftStrategy, Multiple, Muscle5Aligner, Muscle5Mode,
    RecordsAligner,
};
pub use pairwise::{
    AdaptiveAffineAligner, AdaptiveAffineAlignment, AdaptiveAffineBackend, AdaptiveAffineDecision,
    AdaptiveAffineScore, AdaptiveAffineWorkspace, GotohAligner, GotohMode, LinearAlgorithm,
    LinearAligner, MyersAligner, Pairwise, WavefrontAligner, XDropAligner,
};
pub use report::{
    ColorMode, CoordinateStyle, MultipleReport, MultipleRun, MultipleRunInfo, MultipleStatistics,
    MyersDistanceReport, MyersDistanceRun, MyersSearchReport, MyersSearchRun, PairwiseMetric,
    PairwiseReport, PairwiseRun, PairwiseRunInfo, PairwiseStatistics, PairwiseValueReport,
    PairwiseValueRun, PoastaQueryReport, PoastaQueryRun, ReportError, ReportField, ReportOptions,
    Verbosity,
};
pub use result::{
    MultipleAlignmentExt, MultipleAlignmentView, PairwiseAlignmentExt, PairwiseAlignmentView,
};
