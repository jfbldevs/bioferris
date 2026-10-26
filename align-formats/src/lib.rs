//! Input and output formats for validated biological sequences, alignments, and
//! partial-order graphs.
//!
//! FASTA is deliberately separate from the scientific kernels: parsing cannot
//! change an algorithm's alphabet, scoring objective, or traceback semantics.
//! CLUSTAL, PHYLIP, NEXUS, Stockholm, TSV, GFA, DOT, and human-text writers
//! serialize native typed results; they never recover data by parsing presentation
//! output.
#![forbid(unsafe_code)]

mod clustal;
mod error;
mod fasta;
mod myers;
mod nexus;
mod pairwise;
mod phylip;
mod poa;
mod stockholm;
mod text;

pub use clustal::{ClustalWrite, write_clustal, write_clustal_to};
pub use error::FormatError;

pub use fasta::{
    FastaError, FastaWrite, parse_fasta, read_fasta, write_alignment_fasta,
    write_alignment_fasta_to, write_fasta, write_fasta_to,
};
pub use myers::{MyersTsvWrite, write_myers_tsv, write_myers_tsv_to};
pub use nexus::{NexusAlphabet, NexusWrite, write_nexus, write_nexus_to};
pub use pairwise::{write_pairwise_fasta, write_pairwise_fasta_to};
pub use phylip::{
    PhylipLayout, PhylipNameStyle, PhylipOptions, PhylipWrite, write_phylip, write_phylip_to,
    write_phylip_with_options, write_phylip_with_options_to,
};
pub use poa::{PoaGraphWrite, write_poa_dot, write_poa_dot_to, write_poa_gfa, write_poa_gfa_to};
pub use stockholm::{StockholmWrite, write_stockholm, write_stockholm_to};
pub use text::{write_text_report, write_text_report_to};
