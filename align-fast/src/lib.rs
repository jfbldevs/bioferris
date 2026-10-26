//! Accelerated sequence-alignment algorithms reconstructed from publications.
//!
//! Implemented profiles include exact global affine WFA, affine X-drop seed
//! extension, and Myers's exact bit-vector unit edit distance.
#![forbid(unsafe_code)]

mod adaptive;
mod compact_wavefront;
mod error;
mod myers;
mod verification;
mod wavefront;
mod xdrop;

pub use adaptive::{
    AdaptiveAffineBackend, AdaptiveAffineDecision, AdaptiveAffinePlan, adaptive_affine_decision,
    adaptive_affine_plan,
};
pub use error::{MyersError, UnitEditVerificationError, WavefrontError, XDropError};
pub use myers::{MyersMatch, MyersMemoryEstimate, MyersPattern, MyersWorkspace};
pub use verification::{UnitEditVerificationReport, verify_unit_edit_alignment};
pub use wavefront::{
    WavefrontAffine, WavefrontAffinePenalties, WavefrontAlignment, WavefrontAlignmentView,
    WavefrontMemoryEstimate, WavefrontScore, WavefrontTraceRun, WavefrontTracebackView,
    WavefrontWorkspace,
};
pub use xdrop::{
    XDrop, XDropAlignment, XDropDirection, XDropMemoryEstimate, XDropScore, XDropSeed,
    XDropSeedAlignment, XDropSeedScore, XDropTermination,
};
