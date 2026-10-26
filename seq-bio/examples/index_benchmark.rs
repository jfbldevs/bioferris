//! Small reproducible phase measurements, not a general indexing speed claim.

use std::fs::{self, OpenOptions};
use std::hint::black_box;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use bioferris::io::index::{IndexOptions, IndexedFormat, SequenceIndex};
use bioferris::prelude::Dna;

const RECORDS: usize = 10_000;
const SYMBOLS: usize = 256;
const LOOKUPS: usize = 1_000;

struct Directory(PathBuf);

impl Directory {
    fn create() -> io::Result<Self> {
        for attempt in 0..128 {
            let path = std::env::temp_dir().join(format!(
                "bioferris-index-benchmark-{}-{attempt}",
                std::process::id()
            ));
            match fs::create_dir(&path) {
                Ok(()) => return Ok(Self(path)),
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => {}
                Err(error) => return Err(error),
            }
        }
        Err(io::Error::other("cannot reserve a benchmark directory"))
    }
}

impl Drop for Directory {
    fn drop(&mut self) {
        // This directory was exclusively created above and contains only this
        // benchmark's generated source and sidecars, never user input files.
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn measure(format: IndexedFormat) -> Result<(), Box<dyn std::error::Error>> {
    let directory = Directory::create()?;
    let source = directory.0.join("records.txt");
    let mut writer = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&source)?,
    );
    let sequence = "ACGT".repeat(SYMBOLS / 4);
    let qualities = "I".repeat(SYMBOLS);
    for record in 0..RECORDS {
        let prefix = if format == IndexedFormat::Fasta {
            '>'
        } else {
            '@'
        };
        writeln!(writer, "{prefix}record{record:05}\n{sequence}")?;
        if format == IndexedFormat::FastqSanger {
            writeln!(writer, "+\n{qualities}")?;
        }
    }
    writer.flush()?;
    drop(writer);
    for iteration in 1..=3 {
        let options = IndexOptions::default();
        let started = Instant::now();
        let index = SequenceIndex::<Dna>::build(&source, format, options)?;
        let build = started.elapsed();
        let sidecar = directory.0.join(format!("index-{iteration}.json"));
        index.save_new(&sidecar)?;
        let source_bytes = index.source_len();
        drop(index);
        let started = Instant::now();
        let mut index = SequenceIndex::<Dna>::open(&source, &sidecar, options)?;
        let open = started.elapsed();
        let started = Instant::now();
        let mut total_symbols = 0usize;
        for lookup in 0..LOOKUPS {
            let id = format!("record{:05}", (lookup * 7919) % RECORDS);
            let record = index
                .get(&id)?
                .ok_or_else(|| io::Error::other("generated record was not indexed"))?;
            total_symbols += black_box(record.sequence().len());
        }
        let retrieval = started.elapsed();
        if total_symbols != LOOKUPS * SYMBOLS || index.len() != RECORDS {
            return Err(io::Error::other("benchmark correctness check failed").into());
        }
        println!(
            "{} {iteration} {RECORDS} {SYMBOLS} {source_bytes} {} {:.3} {:.3} {LOOKUPS} {:.3}",
            format.name(),
            fs::metadata(&sidecar)?.len(),
            build.as_secs_f64() * 1000.0,
            open.as_secs_f64() * 1000.0,
            retrieval.as_secs_f64() * 1_000_000.0 / f64::from(u32::try_from(LOOKUPS)?),
        );
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!(
        "format iteration records symbols_per_record source_bytes index_bytes build_ms verified_open_ms lookups mean_lookup_us"
    );
    measure(IndexedFormat::Fasta)?;
    measure(IndexedFormat::FastqSanger)
}
