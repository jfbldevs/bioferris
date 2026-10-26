//! Composition, explicit GC policies, bounded k-mers and streaming summaries.

use std::io::Cursor;

use bioferris::io::FastaReader;
use bioferris::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sequence = dna("acgnsw")?;
    let basic = sequence.statistics()?;
    println!("DNA: {sequence}");
    println!("Length: {}; ambiguous: {}", basic.length, basic.ambiguous);
    println!("Composition: {:?}", basic.composition);
    for (name, policy) in [
        ("remove", GcPolicy::Remove),
        ("ignore", GcPolicy::Ignore),
        ("weighted", GcPolicy::Weighted),
    ] {
        let summary = sequence.statistics_with(StatisticsOptions {
            gc_policy: policy,
            ..StatisticsOptions::default()
        })?;
        if let Some(gc) = summary.gc {
            println!(
                "GC ({name}): {:.6} / {} = {:.6}",
                gc.gc_equivalent, gc.denominator, gc.fraction
            );
        }
    }
    for (name, ambiguity) in [
        ("skip", KmerAmbiguity::Skip),
        ("literal", KmerAmbiguity::Literal),
    ] {
        let summary = sequence.statistics_with(StatisticsOptions {
            kmers: Some(KmerOptions {
                ambiguity,
                ..KmerOptions::new(2)
            }),
            ..StatisticsOptions::default()
        })?;
        if let Some(kmers) = summary.kmers {
            println!("2-mers ({name}): {:?}", kmers.counts);
            println!(
                "Windows: {}; counted: {}; skipped: {}",
                kmers.windows, kmers.counted, kmers.skipped
            );
        }
    }

    let fasta = ">short\nG\n>long\nAAAA\n";
    let records = FastaReader::<_, Dna>::new(Cursor::new(fasta));
    let summary = summarize_records(records, StatisticsOptions::default())?;
    println!(
        "Stream: {} records; {} symbols",
        summary.records, summary.length
    );
    if let Some(gc) = summary.gc {
        println!("Stream GC: {:.6}", gc.fraction);
    }
    let peptide = protein("MABZXUO")?.statistics()?;
    println!(
        "Protein: {} symbols; {} ambiguous; GC absent: {}",
        peptide.length,
        peptide.ambiguous,
        peptide.gc.is_none()
    );
    let empty = summarize_sequences(
        std::iter::empty::<&DnaSequence>(),
        StatisticsOptions::default(),
    )?;
    println!(
        "No records: {}; minimum: {:?}",
        empty.records, empty.minimum_length
    );

    let limited = StatisticsOptions {
        kmers: Some(KmerOptions {
            max_windows: 1,
            ..KmerOptions::new(2)
        }),
        ..StatisticsOptions::default()
    };
    if let Err(error) = sequence.statistics_with(limited) {
        println!("Window limit: {}", error.code());
    }
    Ok(())
}
