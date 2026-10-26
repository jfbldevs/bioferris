//! Prints a two-way seeded X-drop extension with detailed traversal evidence.

use align::prelude::*;

fn main() -> align::Result<()> {
    let target = DnaSequence::try_from("TTGATTACAGG")?;
    let query = DnaSequence::try_from("CCGACTACATT")?;
    let run = Pairwise::xdrop()
        .scores(2, -3)
        .gaps(-6, -1)
        .threshold(8)
        .run_from_seed(&target, &query, XDropSeed::new(5, 5), XDropDirection::Both)?;
    println!(
        "{}",
        run.report()
            .labels("Target", "Query")
            .verbosity(Verbosity::Detailed)
    );
    Ok(())
}
