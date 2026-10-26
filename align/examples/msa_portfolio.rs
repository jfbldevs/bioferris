//! Prints compact summaries for Clustal Omega, MUSCLE 5, and six MAFFT profiles.

use align::prelude::*;

fn main() -> align::Result<()> {
    let records = [
        SequenceRecord::new("human", ProteinSequence::try_from("MKTAYIAKQRQ")?)?,
        SequenceRecord::new("mouse", ProteinSequence::try_from("MKTAYIAKQKQ")?)?,
        SequenceRecord::new("fish", ProteinSequence::try_from("MKTAYIAQRQ")?)?,
    ];

    println!(
        "{}\n",
        Multiple::clustal_omega()
            .run(&records)?
            .report()
            .verbosity(Verbosity::Compact)
    );
    println!(
        "{}\n",
        Multiple::muscle5()
            .run(&records)?
            .report()
            .verbosity(Verbosity::Compact)
    );
    for strategy in [
        MafftStrategy::FftNs1,
        MafftStrategy::FftNs2,
        MafftStrategy::FftNsi,
        MafftStrategy::Ginsi,
        MafftStrategy::Linsi,
        MafftStrategy::Einsi,
    ] {
        println!(
            "{}\n",
            Multiple::mafft(strategy)
                .run(&records)?
                .report()
                .verbosity(Verbosity::Compact)
        );
    }
    Ok(())
}
