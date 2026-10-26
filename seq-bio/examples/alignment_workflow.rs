//! Myers distance, a MAFFT alignment and a POASTA graph query in one Rust program.

use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let target = dna("ACGT")?;
    let query = dna("AGT")?;
    let distance = Pairwise::myers().distance(&target, &query)?;
    println!("Edit distance: {distance}");

    let records = vec![
        SequenceRecord::new("reference", target)?,
        SequenceRecord::new("deletion", query.clone())?,
    ];
    let multiple = Multiple::mafft(MafftStrategy::FftNs1).align(&records)?;
    println!("MSA records: {}", multiple.alignment().sequence_count());

    let graph = Graph::poasta().align(&records)?;
    let aligned_query = Graph::poasta().run_query(graph.graph(), &query)?;
    println!("Graph paths: {}", graph.graph().path_count());
    println!("{}", aligned_query.report());
    Ok(())
}
