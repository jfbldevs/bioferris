//! Read reported BLAST coordinates without guessing missing biological values.

use std::io::Cursor;

use bioferris::io::blast::{BlastTabularOptions, parse_blast_tabular};

const INPUT: &str = "q\ts\t100\t4\t0\t0\t1\t4\t8\t5\t1e-20\t42\n";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let document = parse_blast_tabular(Cursor::new(INPUT), BlastTabularOptions::default())?;
    let query = &document.queries()[0];
    let hit = &query.hits()[0];
    let hsp = &hit.hsps()[0];
    let query_coordinates = hsp.query_coordinates().ok_or("missing query coordinates")?;
    let hit_coordinates = hsp.hit_coordinates().ok_or("missing subject coordinates")?;

    println!("Query: {} ({} hit)", query.id(), query.hits().len());
    println!("Query range: {:?}", query_coordinates.range());
    println!(
        "Subject: {} {:?} (reported {}..{})",
        hit.id(),
        hit_coordinates.range(),
        hit_coordinates.reported_start(),
        hit_coordinates.reported_end()
    );
    println!("Reported strand: {:?}", hsp.hit_strand());
    println!(
        "E-value source: {}",
        hsp.field("evalue").ok_or("missing evalue")?
    );

    let minimal = parse_blast_tabular(
        Cursor::new("q\ts\n"),
        BlastTabularOptions {
            fields: Some(vec!["qseqid".into(), "sseqid".into()]),
            ..BlastTabularOptions::default()
        },
    )?;
    println!(
        "Omitted e-value: {:?}",
        minimal.queries()[0].hits()[0].hsps()[0].evalue()
    );

    let mut limited = BlastTabularOptions::default();
    limited.limits.max_queries = 0;
    if let Err(error) = parse_blast_tabular(Cursor::new(INPUT), limited) {
        println!("Rejected zero query budget: {}", error.code());
    }
    Ok(())
}
