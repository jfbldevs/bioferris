//! Builds a POASTA graph and prints one query-to-graph traceback.

use align::prelude::*;

fn main() -> align::Result<()> {
    let records = [
        SequenceRecord::new("reference", DnaSequence::try_from("ACGTACGT")?)?,
        SequenceRecord::new("allele", DnaSequence::try_from("ACGTTCGT")?)?,
        SequenceRecord::new("deletion", DnaSequence::try_from("ACGTCGT")?)?,
    ];
    let graph_run = Graph::poasta().run(&records)?;
    let query = DnaSequence::try_from("ACGTTACGT")?;
    let query_run = Graph::poasta().run_query(graph_run.result().graph(), &query)?;
    println!("{}", query_run.report().verbosity(Verbosity::Detailed));
    Ok(())
}
