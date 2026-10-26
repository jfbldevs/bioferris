//! Builds a DNA PO-MSA with the exact accelerated POASTA profile.

use std::error::Error;

use align_core::{Dna, DnaSequence};
use align_msa::SequenceRecord;
use align_poa::{PoastaAffineCosts, PoastaPoa};

fn main() -> Result<(), Box<dyn Error>> {
    let records = [
        SequenceRecord::<Dna>::new("reference", DnaSequence::from_ascii("ACGTACGT")?)?,
        SequenceRecord::new("allele", DnaSequence::from_ascii("ACGTTCGT")?)?,
        SequenceRecord::new("deletion", DnaSequence::from_ascii("ACGTCGT")?)?,
        SequenceRecord::new("insertion", DnaSequence::from_ascii("ACGTAACGT")?)?,
    ];
    let aligner = PoastaPoa::new(PoastaAffineCosts::new(4, 6, 2)?);
    let result = aligner.align(&records)?;

    println!("profile: {}", aligner.profile_identifier());
    println!(
        "graph: {} nodes, {} classes, {} edges, {} paths",
        result.graph().node_count(),
        result.graph().class_count(),
        result.graph().edge_count(),
        result.graph().path_count()
    );
    for step in result.steps() {
        if let Some(alignment) = step.alignment().alignment() {
            let statistics = alignment.statistics();
            println!(
                "record {}: score={}, expanded={}, pruned={}, exact_dfs={}",
                step.record_index(),
                alignment.score(),
                statistics.expanded_states(),
                statistics.bubble_pruned_states(),
                statistics.greedy_exact_matches()
            );
        }
    }
    for row in result.alignment().rows() {
        println!("{:<12} {}", row.identifier(), row.to_ascii_with_gaps());
    }
    println!(
        "modeled peak: {} bytes",
        result.memory_estimate().estimated_peak_bytes()
    );
    Ok(())
}
