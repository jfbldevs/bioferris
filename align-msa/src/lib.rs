//! Shared, algorithm-independent multiple-sequence-alignment primitives.
//!
//! The crate owns validated named inputs, rectangular gapped alignments, weighted
//! profiles, guide trees, distance matrices, and reference-comparison metrics.
//! Algorithm pipelines remain in their family crates.
#![forbid(unsafe_code)]

mod alignment;
mod distance;
mod error;
mod metrics;
mod profile;
mod record;
mod tree;

pub use alignment::{AlignedRow, MultipleAlignment};
pub use distance::DistanceMatrix;
pub use error::MsaError;
pub use metrics::{AlignmentMetrics, compare_alignments, compare_alignments_on_reference_columns};
pub use profile::Profile;
pub use record::SequenceRecord;
pub use tree::{GuideTree, GuideTreeEdge, GuideTreeNode};
