//! Runs a retained MAFFT profile through the ergonomic Rust facade.

use align::prelude::*;

fn main() -> align::Result<()> {
    let records = [
        SequenceRecord::new("human", ProteinSequence::try_from("MKTAYIAKQRQ")?)?,
        SequenceRecord::new("mouse", ProteinSequence::try_from("MKTAYIAKQKQ")?)?,
        SequenceRecord::new("fish", ProteinSequence::try_from("MKTAYIAQRQ")?)?,
    ];

    let run = Multiple::mafft(MafftStrategy::Linsi).run(&records)?;
    println!("{}", run.report());
    Ok(())
}
