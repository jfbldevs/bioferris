//! Builds a partial-order DNA alignment through the ergonomic Rust facade.

use align::prelude::*;

fn main() -> align::Result<()> {
    let records = [
        SequenceRecord::new("reference", DnaSequence::try_from("ACGTACGT")?)?,
        SequenceRecord::new("allele", DnaSequence::try_from("ACGTTCGT")?)?,
        SequenceRecord::new("deletion", DnaSequence::try_from("ACGTCGT")?)?,
    ];

    let run = Graph::poasta().run(&records)?;
    println!("{}", run.report());
    Ok(())
}
