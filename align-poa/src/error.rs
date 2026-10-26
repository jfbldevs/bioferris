use core::fmt;

use align_core::Score;
use align_msa::MsaError;

use crate::{PoaNodeId, PoaPathId};

/// Failure reported by partial-order alignment or graph validation.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum PoaError {
    /// An iterative alignment was requested without any input records.
    EmptyCollection,
    /// Two input records use the same identifier.
    DuplicateIdentifier {
        /// Repeated identifier.
        identifier: String,
    },
    /// A matrix dimension, element count, or byte count cannot fit in `usize`.
    DimensionOverflow {
        /// Number of graph nodes.
        graph_nodes: usize,
        /// Query sequence length.
        query_len: usize,
    },
    /// A dynamic allocation could not be reserved.
    AllocationFailed {
        /// Logical buffer being allocated.
        buffer: &'static str,
        /// Requested element count.
        elements: usize,
    },
    /// A graph counter or modeled byte total overflowed `usize`.
    ArithmeticOverflow {
        /// Operation whose exact result could not be represented.
        context: &'static str,
    },
    /// A POASTA affine cost parameter is negative.
    NegativeAffineCost {
        /// Name of the invalid parameter.
        parameter: &'static str,
        /// Supplied negative cost.
        cost: Score,
    },
    /// The affine gap-opening cost is not strictly greater than extension.
    AffineGapOpeningNotGreaterThanExtension {
        /// Supplied non-negative gap-opening cost.
        opening: Score,
        /// Supplied non-negative gap-extension cost.
        extension: Score,
    },
    /// The affine opening transition `opening + extension` is unrepresentable.
    AffineCostOverflow {
        /// Supplied gap-opening cost.
        opening: Score,
        /// Supplied gap-extension cost.
        extension: Score,
    },
    /// An exact POASTA search-cost operation overflowed its checked accumulator.
    PoastaCostOverflow {
        /// Operation whose exact result could not be represented.
        context: &'static str,
    },
    /// The negated exact POASTA optimum cannot fit in the public score type.
    PoastaScoreOutOfRange {
        /// Exact non-negative optimum accumulated by sparse search.
        cost: u128,
    },
    /// The exact positive score at one DP cell cannot fit in the score type.
    ScoreOverflow {
        /// Stable graph node whose row was being evaluated.
        node: PoaNodeId,
        /// One-based query column.
        column: usize,
    },
    /// An exact affine optimum cannot be represented by the public score type.
    AffineScoreOutOfRange {
        /// Exact mathematical optimum accumulated in the scalar oracle.
        score: i128,
    },
    /// Traceback metadata does not describe a legal predecessor.
    TracebackInvariant {
        /// Stable graph node at the invalid traceback cell.
        node: PoaNodeId,
        /// Query column at the invalid traceback cell.
        column: usize,
    },
    /// Affine traceback metadata does not describe a legal predecessor.
    AffineTracebackInvariant {
        /// Stable graph node, or `None` for the virtual source.
        node: Option<PoaNodeId>,
        /// Consumed query prefix length.
        column: usize,
        /// Affine state being reconstructed.
        state: &'static str,
    },
    /// An alignment operation refers to a missing graph node.
    UnknownNode {
        /// Invalid stable node identifier.
        node: PoaNodeId,
    },
    /// A path identifier is outside the graph's stable path arena.
    UnknownPath {
        /// Invalid stable path identifier.
        path: PoaPathId,
    },
    /// A graph invariant was violated.
    InvalidGraph {
        /// Precise invariant diagnostic.
        message: String,
    },
    /// Independent score or traceback replay failed.
    VerificationFailed {
        /// Precise verification diagnostic.
        message: String,
    },
    /// A shared multiple-alignment primitive rejected the serialization.
    InvalidMultipleAlignment(MsaError),
}

impl fmt::Display for PoaError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyCollection => formatter.write_str("at least one sequence is required"),
            Self::DuplicateIdentifier { identifier } => {
                write!(formatter, "duplicate sequence identifier {identifier:?}")
            }
            Self::DimensionOverflow {
                graph_nodes,
                query_len,
            } => write!(
                formatter,
                "POA dimensions overflow for {graph_nodes} graph nodes and query length {query_len}"
            ),
            Self::AllocationFailed { buffer, elements } => {
                write!(
                    formatter,
                    "could not allocate {elements} elements for {buffer}"
                )
            }
            Self::ArithmeticOverflow { context } => {
                write!(
                    formatter,
                    "POA arithmetic overflow while computing {context}"
                )
            }
            Self::NegativeAffineCost { parameter, cost } => {
                write!(formatter, "affine {parameter} cost {cost} is negative")
            }
            Self::AffineGapOpeningNotGreaterThanExtension { opening, extension } => write!(
                formatter,
                "affine gap-opening cost {opening} is not greater than extension cost {extension}"
            ),
            Self::AffineCostOverflow { opening, extension } => write!(
                formatter,
                "affine opening transition {opening} + {extension} is unrepresentable"
            ),
            Self::PoastaCostOverflow { context } => {
                write!(formatter, "POASTA cost overflow while computing {context}")
            }
            Self::PoastaScoreOutOfRange { cost } => write!(
                formatter,
                "exact POASTA optimum cost {cost} cannot be negated into the public score type"
            ),
            Self::ScoreOverflow { node, column } => write!(
                formatter,
                "POA score overflow at graph node {} and query column {column}",
                node.index()
            ),
            Self::AffineScoreOutOfRange { score } => write!(
                formatter,
                "exact affine POA optimum {score} does not fit the public score type"
            ),
            Self::TracebackInvariant { node, column } => write!(
                formatter,
                "invalid POA traceback at graph node {} and query column {column}",
                node.index()
            ),
            Self::AffineTracebackInvariant {
                node,
                column,
                state,
            } => match node {
                Some(node) => write!(
                    formatter,
                    "invalid affine POA traceback in state {state} at graph node {} and query column {column}",
                    node.index()
                ),
                None => write!(
                    formatter,
                    "invalid affine POA traceback in state {state} at the virtual source and query column {column}"
                ),
            },
            Self::UnknownNode { node } => {
                write!(formatter, "unknown POA graph node {}", node.index())
            }
            Self::UnknownPath { path } => {
                write!(formatter, "unknown POA source path {}", path.index())
            }
            Self::InvalidGraph { message } => write!(formatter, "invalid POA graph: {message}"),
            Self::VerificationFailed { message } => {
                write!(formatter, "POA verification failed: {message}")
            }
            Self::InvalidMultipleAlignment(error) => {
                write!(formatter, "invalid rectangular POA serialization: {error}")
            }
        }
    }
}

impl std::error::Error for PoaError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidMultipleAlignment(error) => Some(error),
            _ => None,
        }
    }
}

impl From<MsaError> for PoaError {
    fn from(value: MsaError) -> Self {
        Self::InvalidMultipleAlignment(value)
    }
}
