//! POASTA sequence-to-graph and iterative partial-order alignment.
//!
//! [`POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_V1`] implements the
//! paper's exact sparse A* search, exact-match DFS, and superbubble pruning;
//! [`LEE_2002_POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_ITERATIVE_V1`]
//! integrates it with deterministic graph incorporation. The graph is authoritative; a
//! rectangular multiple alignment is a deterministic serialization of its
//! alignment classes.
#![forbid(unsafe_code)]

// The dense affine recurrence is both an independent oracle and an explicitly
// selectable exact backend. Automatic execution does not silently fall back to it.
#[allow(dead_code)]
mod affine;
#[cfg(test)]
#[allow(dead_code, clippy::unused_self)]
mod affine_iterative;
#[cfg(test)]
#[allow(dead_code)]
mod alignment;
mod error;
mod graph;
mod iterative_support;
#[cfg(test)]
#[allow(dead_code, clippy::unused_self)]
mod lee_2002;
mod operation;
mod poasta;
mod poasta_iterative;
mod superbubble;

#[cfg(test)]
mod properties;

pub use affine::{AffinePoaVerificationReport, PoastaAffineCosts};
#[cfg(test)]
pub(crate) use affine::{
    GlobalAffinePoa, GlobalAffineSequenceGraphAlignment,
    verify_global_affine_sequence_graph_alignment,
};
#[cfg(test)]
#[allow(unused_imports)]
pub(crate) use affine_iterative::{
    AffinePoa, AffinePoaMemoryEstimate, AffinePoaResult, AffinePoaStep, AffinePoaStepAlignment,
    LEE_2002_POASTA_2025_GLOBAL_AFFINE_ITERATIVE_V1,
};
#[cfg(test)]
pub(crate) use alignment::{
    LocalPoa, PoaVerificationReport, SequenceGraphAlignment, verify_sequence_graph_alignment,
};
pub use error::PoaError;
pub use graph::{
    PoaAlignmentClass, PoaClassId, PoaEdge, PoaGraph, PoaNode, PoaNodeId, PoaOccurrence, PoaPath,
    PoaPathId,
};
#[cfg(test)]
pub(crate) use lee_2002::LeePoa;
pub use operation::{PoaAlignedPair, PoaOperation};
pub use poasta::{
    POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_V1, Poasta, PoastaAlignment, PoastaBackend,
    PoastaMemoryEstimate, PoastaSearchOptions, PoastaSearchStatistics, PoastaSearchWorkCounters,
    verify_poasta_alignment,
};
pub use poasta_iterative::{
    LEE_2002_POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_ITERATIVE_V1, PoastaPoa,
    PoastaPoaDetailedProfile, PoastaPoaMemoryEstimate, PoastaPoaResult, PoastaPoaStep,
    PoastaPoaStepAlignment,
};
