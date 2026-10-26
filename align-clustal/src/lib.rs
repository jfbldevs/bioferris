//! Clustal Omega multiple-sequence alignment reconstructed from its publications.
//!
//! The versioned profiles accept typed protein, DNA, and RNA records. Independent
//! guide distances and progressive subtrees can be evaluated concurrently while
//! preserving deterministic output.
#![forbid(unsafe_code)]
// Public memory diagnostics retain explicit unit suffixes as part of the stable API.
#![allow(clippy::struct_field_names)]

mod alphabet;
mod clustal_omega;
mod omega_distance;
mod omega_error;
mod omega_hmm;
mod omega_matrix;
mod omega_parameters;

pub use alphabet::ClustalAlphabet;
pub use clustal_omega::{
    ClustalOmega, ClustalOmegaAlignment, ClustalOmegaMemoryEstimate, OmegaGuideEvidence,
    OmegaGuideTreeMethod, OmegaProgressiveStep,
};
pub use omega_error::ClustalOmegaError;
pub use omega_hmm::OmegaEngine;
pub use omega_parameters::{
    ClustalOmegaParameters, OmegaEngineMode, OmegaGuideTreeMode, OmegaMoleculeProfile,
};
