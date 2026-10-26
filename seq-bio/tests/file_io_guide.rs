//! Literal inputs and Rust results documented in the file-I/O guide.

fn documented_block(marker: &str) -> &'static str {
    let document = include_str!("../../docs/api/file-io-v1.md");
    let tail = document
        .split_once(&format!("<!-- {marker} -->"))
        .unwrap()
        .1;
    let fence = tail.find("```").unwrap();
    let content = tail[fence + 3..].split_once('\n').unwrap().1;
    content.split_once("\n```").unwrap().0
}

#[test]
fn documented_fastq_input_matches_its_fixture() {
    assert_eq!(
        documented_block("io-fastq-input"),
        include_str!("fixtures/io-example.fastq").trim_end()
    );
}

#[test]
fn documented_rust_record_results_match_the_public_api() {
    use bioferris::io::{FastaReader, FastqReader, one_record};
    use bioferris::prelude::*;
    let record = one_record(FastqReader::<_, Dna>::new(std::io::Cursor::new(
        include_str!("fixtures/io-example.fastq"),
    )))
    .unwrap();
    let observed = format!(
        "ID: {}\nDescription: {}\nSource: {}\nDNA: {}\nPhred: {:?}",
        record.id(),
        record.description().unwrap_or(""),
        record.source_sequence(),
        record.sequence(),
        record.qualities().unwrap_or(&[])
    );
    assert_eq!(documented_block("io-rust-record-output"), observed);
    let input = ">gene description\nATGGCTTAA\n>second\nATGTAA\n";
    let lines: Vec<_> = FastaReader::<_, Dna>::new(std::io::Cursor::new(input))
        .map(|record| {
            let record = record.unwrap();
            format!(
                "{}: {}",
                record.id(),
                record.sequence().translate_cds().unwrap()
            )
        })
        .collect();
    assert_eq!(documented_block("io-rust-fasta-output"), lines.join("\n"));
}
