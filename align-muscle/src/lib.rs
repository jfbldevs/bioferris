//! MUSCLE 5 multiple-sequence alignment reconstructed from its scientific description.
//!
//! The public profile reconstructs the MUSCLE 5 `Align`, `Super5`, and ensemble
//! pipelines. It performs no external I/O and has no runtime dependency on
//! another aligner.
#![forbid(unsafe_code)]
// Public memory diagnostics retain explicit unit suffixes as part of the stable API.
#![allow(clippy::struct_field_names)]

mod muscle5;
mod muscle5_emissions;
mod muscle5_error;
mod muscle5_hmm;
mod muscle5_parameters;
mod muscle5_posterior;
mod muscle5_profile;
mod muscle5_rng;
mod muscle5_super5;
mod muscle5_tree;
pub use muscle5::{
    Muscle5, Muscle5Alignment, Muscle5Ensemble, Muscle5EnsembleMember, Muscle5MemoryEstimate,
    Muscle5ProgressiveStep, Muscle5RefinementAttempt,
};
pub use muscle5_emissions::Muscle5Alphabet;
pub use muscle5_error::Muscle5Error;
pub use muscle5_parameters::{GuideTreePermutation, Muscle5Mode, Muscle5Parameters};
pub use muscle5_super5::{
    Muscle5ProfileAlignmentMethod, Muscle5RedundancyAssignment, Muscle5Super5Cluster,
    Muscle5Super5Evidence, Muscle5Super5Merge,
};
