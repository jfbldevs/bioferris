//! Lossless compressed records, explicit finalization and verified BGZF retrieval.

use std::io::{Cursor, Read, Write};
use std::path::PathBuf;

use bioferris::io::bgzf::{
    BgzfOutputError, BgzfReader, BgzfVirtualOffset, BgzfWriteOptions, BgzfWriter, write_bgzf_new,
};
use bioferris::io::bgzf_index::{BgzfIndexOptions, BgzfSequenceIndex};
use bioferris::io::gzip::{GzipWriteOptions, write_gzip_new};
use bioferris::io::index::IndexedFormat;
use bioferris::io::{FastqReader, one_record, open_input, write_record_fastq};
use bioferris::prelude::*;

fn unique_directory() -> std::io::Result<PathBuf> {
    let timestamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(std::io::Error::other)?
        .as_nanos();
    for attempt in 0..128 {
        let directory = std::env::temp_dir().join(format!(
            "bioferris-compressed-io-{}-{timestamp}-{attempt}",
            std::process::id()
        ));
        match std::fs::create_dir(&directory) {
            Ok(()) => return Ok(directory),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => return Err(error),
        }
    }
    Err(std::io::Error::other(
        "no unique example directory available",
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let directory = unique_directory()?;
    let gzip_path = directory.join("read.fastq.gz");
    let bgzf_path = directory.join("read.fastq.bgz");
    let index_path = directory.join("read.index.json");
    let record = Record::<Dna>::from_ascii("read demo", "aCgT", Some(vec![0, 20, 40, 93]))?;
    println!("Input case: {}", record.source_sequence());

    write_gzip_new(&gzip_path, GzipWriteOptions::default(), |output| {
        write_record_fastq(output, &record)
    })?;
    let decoded = one_record(FastqReader::<_, Dna>::new(open_input(&gzip_path)?))?;
    println!(
        "Gzip record: {} {} {:?}",
        decoded.id(),
        decoded.source_sequence(),
        decoded.qualities().ok_or("missing qualities")?
    );

    write_bgzf_new(&bgzf_path, BgzfWriteOptions::default(), |output| {
        write_record_fastq(output, &record)
    })?;
    let index = BgzfSequenceIndex::<Dna>::build(
        &bgzf_path,
        IndexedFormat::FastqSanger,
        BgzfIndexOptions::default(),
    )?;
    println!("BGZF IDs: {}", index.ids().collect::<Vec<_>>().join(","));
    index.save_new(&index_path)?;
    let mut reopened =
        BgzfSequenceIndex::<Dna>::open(&bgzf_path, &index_path, BgzfIndexOptions::default())?;
    let found = reopened.get("read")?.ok_or("missing expected record")?;
    println!(
        "BGZF record: {} {} {:?}",
        found.id(),
        found.source_sequence(),
        found.qualities().ok_or("missing qualities")?
    );
    println!("Normalized: {}", found.sequence());
    println!("Missing ID: {}", reopened.get("absent")?.is_none());

    let refused = write_bgzf_new(&bgzf_path, BgzfWriteOptions::default(), |output| {
        write_record_fastq(output, &record)
    });
    println!(
        "Overwrite refused: {}",
        matches!(refused, Err(BgzfOutputError::DestinationExists(_)))
    );
    if let Err(error) = BgzfSequenceIndex::<Dna>::build(
        &gzip_path,
        IndexedFormat::FastqSanger,
        BgzfIndexOptions::default(),
    ) {
        println!("Ordinary gzip as BGZF index source: {}", error.code());
    }

    let mut writer = BgzfWriter::with_options(
        Vec::new(),
        BgzfWriteOptions {
            block_size: 4,
            ..Default::default()
        },
    )?;
    writer.write_all(b"ACGTACGT")?;
    let mut reader = BgzfReader::new(Cursor::new(writer.finish()?));
    reader.seek_virtual(BgzfVirtualOffset::new(0, 2)?)?;
    let mut suffix = String::new();
    reader.read_to_string(&mut suffix)?;
    println!("BGZF suffix: {suffix}");
    println!("BGZF EOF: {:?}", reader.finish()?);
    if let Err(error) = BgzfVirtualOffset::new(0, 65_536) {
        println!("Unrepresentable offset: {}", error.code());
    }

    // Remove only the three files this run created in its own new directory.
    // Never recursively remove a workspace or clean unrelated/user input files.
    drop(reopened);
    drop(index);
    std::fs::remove_file(gzip_path)?;
    std::fs::remove_file(bgzf_path)?;
    std::fs::remove_file(index_path)?;
    std::fs::remove_dir(directory)?;
    Ok(())
}
