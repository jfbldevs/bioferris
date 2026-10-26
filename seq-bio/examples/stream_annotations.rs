//! Stream `GenBank` records one at a time, skipping records that fail to parse.

use bioferris::Dna;
use bioferris::io::annotated::{RecordErrorPolicy, genbank};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // The documented three-record file; pass a path to read another one.
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "seq-bio/tests/fixtures/annotated-stream.gb".to_owned());
    let records = genbank::<Dna>(&path)?.with_record_errors(RecordErrorPolicy::Continue);
    for result in records {
        match result {
            Ok(entry) => println!(
                "{}: {} bp, {} feature(s)",
                entry.record().id(),
                entry.record().sequence().len(),
                entry.features().len()
            ),
            Err(error) => println!(
                "Skipped record {} at line {}: {}",
                error.record().unwrap_or_default(),
                error.line().unwrap_or_default(),
                error.code()
            ),
        }
    }
    Ok(())
}
