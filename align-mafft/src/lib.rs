//! Independent MAFFT-family multiple-sequence alignment algorithms.
//!
//! The public profiles reconstruct the 2002 FFT-NS-1, FFT-NS-2, and FFT-NS-i
//! pipelines plus the 2005 G-INS-i/L-INS-i and 2008 E-INS-i pipelines from
//! their scientific descriptions. They perform no external I/O and have no
//! runtime dependency on another aligner.
#![forbid(unsafe_code)]
// Public memory diagnostics retain explicit unit suffixes as part of the stable API.
#![allow(clippy::struct_field_names)]

mod aligned_distance;
mod alphabet;
mod einsi;
mod error;
mod fft;
mod fft_ns1;
mod fft_ns2;
mod fft_nsi;
mod generalized_affine;
mod ginsi;
mod linsi;
mod local_affine;
mod matrix;
mod parameters;
mod profile;
mod refinement;
mod segment;
mod tuple;
mod upgma;
mod wsp;

pub use aligned_distance::FftNs2DistanceMatrix;
pub use alphabet::MafftAlphabet;
pub use einsi::{
    EINSI_PROFILE_IDENTIFIER, EinsiAlignment, EinsiCycle, EinsiLibraryStatistics,
    EinsiMemoryEstimate, EinsiMerge, EinsiPairEvidence, EinsiPairWeights, EinsiPartition,
    EinsiTermination, MafftEinsi,
};
pub use error::MafftError;
pub use fft_ns1::{
    FFT_NS1_PROFILE_IDENTIFIER, FftNs1Alignment, FftNs1Anchor, FftNs1MemoryEstimate, FftNs1Merge,
    MafftFftNs1,
};
pub use fft_ns2::{
    FFT_NS2_PROFILE_IDENTIFIER, FftNs2Alignment, FftNs2Anchor, FftNs2MemoryEstimate, FftNs2Merge,
    MafftFftNs2,
};
pub use fft_nsi::{
    FFT_NSI_PROFILE_IDENTIFIER, FftNsiAlignment, FftNsiCandidatePhaseTimings, FftNsiCycle,
    FftNsiDetailedProfile, FftNsiMemoryEstimate, FftNsiTermination, MafftFftNsi,
};
pub use ginsi::{
    GINSI_PROFILE_IDENTIFIER, GinsiAlignment, GinsiCycle, GinsiLibraryStatistics,
    GinsiMemoryEstimate, GinsiMerge, GinsiPairEvidence, GinsiPairWeights, GinsiPartition,
    GinsiTermination, MafftGinsi,
};
pub use linsi::{
    LINSI_PROFILE_IDENTIFIER, LinsiAlignment, LinsiCycle, LinsiLibraryStatistics,
    LinsiMemoryEstimate, LinsiMerge, LinsiPairEvidence, LinsiPairWeights, LinsiPartition,
    LinsiTermination, MafftLinsi,
};
pub use parameters::{
    EINSI_GROUP_GAP_OPEN, EINSI_GROUP_RANDOM_MATCH_OFFSET, EINSI_IMPORTANCE_WEIGHT,
    EINSI_LOCAL_GAP_EXTENSION, EINSI_LOCAL_GAP_OPEN, EINSI_LOCAL_OFFSET_PENALTY, EINSI_MAX_CYCLES,
    EINSI_MIXED_LINKAGE_FRACTION, EINSI_PAIRED_SKIP_EXTENSION, EINSI_SKIP_GAP_OPEN,
    EinsiParameters, FFT_NS2_POISSON_IDENTITY_FLOOR, FFT_NSI_BRANCH_LENGTH_FLOOR,
    FFT_NSI_EQUALIZATION_FACTOR, FFT_NSI_MAX_CYCLES, FftNs1Parameters, FftNs2Parameters,
    FftNsiParameters, GINSI_GAP_OPEN, GINSI_IMPORTANCE_WEIGHT, GINSI_MAX_CYCLES,
    GINSI_MIXED_LINKAGE_FRACTION, GINSI_RANDOM_MATCH_OFFSET, GinsiParameters, LINSI_GROUP_GAP_OPEN,
    LINSI_GROUP_RANDOM_MATCH_OFFSET, LINSI_IMPORTANCE_WEIGHT, LINSI_LOCAL_GAP_EXTENSION,
    LINSI_LOCAL_GAP_OPEN, LINSI_LOCAL_OFFSET_PENALTY, LINSI_MAX_CYCLES,
    LINSI_MIXED_LINKAGE_FRACTION, LinsiParameters,
};
pub use refinement::{FftNsiPairWeights, FftNsiPartition};
