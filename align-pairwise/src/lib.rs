//! Exact pairwise alignment algorithms reconstructed from scientific descriptions.
//!
//! The implemented algorithms are global Needleman--Wunsch, local Smith--Waterman,
//! and global, local, and free-end Gotoh recurrences. Exact affine objectives
//! provide scalar and runtime-dispatched AVX2 complete-result and score-only
//! execution paths.
#![forbid(unsafe_code)]

mod affine_modes;
mod affine_score;
mod affine_simd;
mod end_gap;
mod error;
mod gotoh;
mod needleman_wunsch;
mod smith_waterman;
mod verification;

pub use affine_modes::{
    AFFINE_END_GAP_V1, AffineEndGapAligner, AffineLocalAligner,
    SMITH_WATERMAN_GOTOH_LOCAL_AFFINE_V1,
};
pub use affine_score::{AFFINE_SCORE_ONLY_V1, AffineScoreMemoryEstimate};
pub use affine_simd::{
    AFFINE_SIMD_PREFIX_SCAN_V1, AffineBackend, AffineKernelPolicy, AffineSimdMemoryEstimate,
    AffineSimdScoreMemoryEstimate,
};
pub use end_gap::FreeEndGaps;
pub use error::{AlignmentError, VerificationError};
pub use gotoh::{AffineWorkspace, Gotoh, GotohMemoryEstimate};
pub use needleman_wunsch::{MemoryEstimate, NeedlemanWunsch};
pub use smith_waterman::{SmithWaterman, SmithWatermanMemoryEstimate};
pub use verification::{
    VerificationReport, verify_affine_end_gap_alignment, verify_global_affine_alignment,
    verify_global_alignment, verify_local_affine_alignment, verify_local_alignment,
};
