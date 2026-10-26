//! Stream BLAST XML reports one complete query at a time, as for files too large
//! to collect into one document.

use bioferris::io::blast::{BlastXmlFormat, BlastXmlOptions, BlastXmlReader};
use bioferris::io::open_input;

const REPORTS: [(BlastXmlFormat, &str); 4] = [
    (
        BlastXmlFormat::Legacy,
        "seq-bio/tests/fixtures/blast-example.xml",
    ),
    (
        BlastXmlFormat::Legacy,
        "seq-bio/tests/fixtures/blast-ncbi-2.9.0-nucleotide.xml",
    ),
    (
        BlastXmlFormat::Xml2,
        "seq-bio/tests/fixtures/blast-ncbi-2.9.0-nucleotide-xml2.xml",
    ),
    (
        BlastXmlFormat::Legacy,
        "seq-bio/tests/fixtures/blast-ncbi-2.9.0-translated.xml",
    ),
];

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for (format, path) in REPORTS {
        let name = path.rsplit('/').next().unwrap_or(path);
        let options = BlastXmlOptions {
            format,
            ..BlastXmlOptions::default()
        };
        // Each item is one complete query; the reader keeps no earlier ones.
        for query in BlastXmlReader::new(open_input(path)?, options)? {
            let query = query?;
            let hsps: usize = query.hits().iter().map(|hit| hit.hsps().len()).sum();
            println!(
                "{name}: {} with {} hit(s), {hsps} HSP(s); first hit {}",
                query.id(),
                query.hits().len(),
                query.hits().first().map_or("none", |hit| hit.id())
            );
        }
    }

    // A query larger than its source budget fails instead of being buffered.
    let mut options = BlastXmlOptions::default();
    options.limits.max_query_bytes = 1024;
    let mut reader = BlastXmlReader::new(open_input(REPORTS[1].1)?, options)?;
    if let Some(Err(error)) = reader.next() {
        println!("Query over a 1 KiB budget: {}: {error}", error.code());
    }
    Ok(())
}
