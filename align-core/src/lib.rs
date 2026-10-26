//! Shared, dependency-free primitives for sequence alignment.
//!
//! The crate validates biological alphabets at construction time and keeps gaps out
//! of unaligned sequences. Alignment algorithms operate on typed [`Sequence`] values
//! and return an operation-based [`PairwiseAlignment`].
#![forbid(unsafe_code)]

mod alignment;
mod alphabet;
mod error;
mod phase;
mod scoring;
mod sequence;

pub use alignment::{AlignmentOperation, PairwiseAlignment, RenderedAlignment};
pub use alphabet::{Alphabet, Dna, Protein, Rna, Symbol};
pub use error::{AlignmentStructureError, ScoringError, SequenceError};
pub use phase::{PhaseTiming, PhaseTimings, PhaseTimingsBuilder, Profiled};
pub use scoring::{
    AffineGap, DenseSubstitutionMatrix, LinearGap, MatchMismatch, Score, SubstitutionScorer,
};
pub use sequence::{DnaSequence, ProteinSequence, RnaSequence, Sequence};
