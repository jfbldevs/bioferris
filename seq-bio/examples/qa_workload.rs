//! Bounded, synthetic QA workloads; measurements are host-specific, not limits.

use std::error::Error;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use bioferris::Dna;
use bioferris::io::ReadLimits;
use bioferris::io::bgzf::{BgzfWriteOptions, write_bgzf_new};
use bioferris::io::bgzf_index::{BgzfIndexOptions, BgzfSequenceIndex};
use bioferris::io::blast::{BlastReadLimits, BlastTabularOptions, BlastTabularReader};
use bioferris::io::index::{IndexOptions, IndexedFormat, SequenceIndex};
use serde_json::{Value, json};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

struct Options {
    directory: PathBuf,
    records: usize,
    symbols: usize,
}

impl Options {
    fn parse(arguments: impl IntoIterator<Item = String>) -> Result<Self> {
        let mut arguments = arguments.into_iter();
        let mut directory = None;
        let mut records = None;
        let mut symbols = None;
        while let Some(option) = arguments.next() {
            let value = arguments.next().ok_or("each option needs a value")?;
            match option.as_str() {
                "--directory" if directory.is_none() => directory = Some(PathBuf::from(value)),
                "--records" if records.is_none() => records = Some(value.parse::<usize>()?),
                "--symbols" if symbols.is_none() => symbols = Some(value.parse::<usize>()?),
                _ => return Err(format!("unknown or duplicate option: {option}").into()),
            }
        }
        let options = Self {
            directory: directory
                .ok_or("--directory is required (an existing trusted directory)")?,
            records: records.unwrap_or(2_048),
            symbols: symbols.unwrap_or(256),
        };
        if !(1..=100_000).contains(&options.records)
            || !(1..=4_096).contains(&options.symbols)
            || options
                .records
                .checked_mul(options.symbols)
                .is_none_or(|size| size > 16 * 1024 * 1024)
        {
            return Err(
                "workload bounds: 1..=100000 records, 1..=4096 symbols, product <=16777216".into(),
            );
        }
        if !fs::metadata(&options.directory)?.is_dir() {
            return Err("--directory must be an existing trusted directory".into());
        }
        Ok(options)
    }
}

/// Linux cumulative process high-water RSS, in the kernel's KiB units.
/// Unavailable on other platforms or when procfs cannot be read.
#[must_use]
pub fn peak_rss_kib() -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let status = fs::read_to_string("/proc/self/status").ok()?;
        status.lines().find_map(|line| {
            line.strip_prefix("VmHWM:")?
                .split_whitespace()
                .next()?
                .parse()
                .ok()
        })
    }
    #[cfg(not(target_os = "linux"))]
    {
        None
    }
}

fn measured<T>(
    metrics: &mut Vec<Value>,
    name: &str,
    operation: impl FnOnce() -> Result<T>,
) -> Result<T> {
    let started = Instant::now();
    let result = operation()?;
    metrics.push(json!({
        "phase": name,
        "elapsed_us": started.elapsed().as_micros(),
        "process_peak_rss_kib": peak_rss_kib(),
    }));
    Ok(result)
}

fn create(path: &Path) -> io::Result<BufWriter<File>> {
    File::options()
        .write(true)
        .create_new(true)
        .open(path)
        .map(BufWriter::new)
}

fn index_options(records: usize) -> IndexOptions {
    IndexOptions {
        read_limits: ReadLimits {
            max_line_bytes: 8_192,
            max_record_symbols: 4_096,
            max_total_bytes: 32 * 1024 * 1024,
        },
        max_records: records,
        max_id_bytes: 2 * 1024 * 1024,
        max_index_bytes: 64 * 1024 * 1024,
    }
}

fn plain_workload(options: &Options, sequence: &str, metrics: &mut Vec<Value>) -> Result<PathBuf> {
    let source = options.directory.join("sequences.fa");
    let sidecar = options.directory.join("sequences.index.json");
    measured(metrics, "source_write", || {
        let mut output = create(&source)?;
        for number in 0..options.records {
            writeln!(output, ">record_{number:06} synthetic\n{sequence}")?;
        }
        output.flush()?;
        Ok(())
    })?;
    let index_options = index_options(options.records);
    let index = measured(metrics, "plain_build", || {
        Ok(SequenceIndex::<Dna>::build(
            &source,
            IndexedFormat::Fasta,
            index_options,
        )?)
    })?;
    if index.len() != options.records {
        return Err("plain index record count mismatch".into());
    }
    measured(metrics, "plain_save", || Ok(index.save_new(&sidecar)?))?;
    drop(index);
    let mut index = measured(metrics, "plain_open", || {
        Ok(SequenceIndex::<Dna>::open(
            &source,
            &sidecar,
            index_options,
        )?)
    })?;
    measured(metrics, "plain_lookup", || {
        for number in 0..options.records {
            let id = format!("record_{number:06}");
            let record = index.get(&id)?.ok_or("plain record missing")?;
            if record.source_sequence() != sequence || record.header() != format!("{id} synthetic")
            {
                return Err("plain record content mismatch".into());
            }
        }
        if index.get("absent")?.is_some() {
            return Err("plain missing-ID lookup returned a record".into());
        }
        Ok(())
    })?;
    drop(index);
    Ok(source)
}

fn bgzf_workload(
    options: &Options,
    sequence: &str,
    source: &Path,
    metrics: &mut Vec<Value>,
) -> Result<PathBuf> {
    let compressed = options.directory.join("sequences.fa.bgz");
    let compressed_sidecar = options.directory.join("sequences.bgzf.index.json");
    measured(metrics, "bgzf_write", || {
        write_bgzf_new(&compressed, BgzfWriteOptions::default(), |output| {
            io::copy(&mut File::open(source)?, output).map(|_| ())
        })?;
        Ok(())
    })?;
    let bgzf_options = BgzfIndexOptions {
        index: index_options(options.records),
        max_compressed_bytes: 32 * 1024 * 1024,
        max_blocks: 4_096,
    };
    let index = measured(metrics, "bgzf_build", || {
        Ok(BgzfSequenceIndex::<Dna>::build(
            &compressed,
            IndexedFormat::Fasta,
            bgzf_options,
        )?)
    })?;
    if index.len() != options.records {
        return Err("BGZF index record count mismatch".into());
    }
    measured(metrics, "bgzf_save", || {
        Ok(index.save_new(&compressed_sidecar)?)
    })?;
    drop(index);
    let mut index = measured(metrics, "bgzf_open", || {
        Ok(BgzfSequenceIndex::<Dna>::open(
            &compressed,
            &compressed_sidecar,
            bgzf_options,
        )?)
    })?;
    measured(metrics, "bgzf_lookup", || {
        for number in 0..options.records {
            let id = format!("record_{number:06}");
            let record = index.get(&id)?.ok_or("BGZF record missing")?;
            if record.source_sequence() != sequence || record.header() != format!("{id} synthetic")
            {
                return Err("BGZF record content mismatch".into());
            }
        }
        if index.get("absent")?.is_some() {
            return Err("BGZF missing-ID lookup returned a record".into());
        }
        Ok(())
    })?;
    drop(index);
    Ok(compressed)
}

fn blast_workload(options: &Options, metrics: &mut Vec<Value>) -> Result<()> {
    let blast = options.directory.join("blast.tsv");
    measured(metrics, "blast_write", || {
        let mut output = create(&blast)?;
        for number in 0..options.records {
            writeln!(output, "q{number:06}\ts{number:06}")?;
        }
        output.flush()?;
        Ok(())
    })?;
    measured(metrics, "blast_read", || {
        let reader = BlastTabularReader::new(
            BufReader::new(File::open(&blast)?),
            BlastTabularOptions {
                fields: Some(vec!["qseqid".into(), "sseqid".into()]),
                limits: BlastReadLimits {
                    max_queries: options.records,
                    max_hits_per_query: 1,
                    max_hsps_per_query: 1,
                    max_query_bytes: 128,
                    ..BlastReadLimits::default()
                },
                ..BlastTabularOptions::default()
            },
        )?;
        let mut count = 0;
        for (number, query) in reader.enumerate() {
            let query = query?;
            if query.id() != format!("q{number:06}")
                || query.hits().len() != 1
                || query.hits()[0].id() != format!("s{number:06}")
                || query.hits()[0].hsps().len() != 1
            {
                return Err("BLAST query/hit/HSP content mismatch".into());
            }
            count += 1;
        }
        if count != options.records {
            return Err("BLAST query count mismatch".into());
        }
        Ok(())
    })?;
    Ok(())
}

/// Run real library operations on bounded synthetic data in caller-owned space.
/// Files are create-only and deliberately retained for inspection, even on error.
///
/// # Errors
/// Rejects invalid arguments/bounds, unavailable or occupied paths, library I/O
/// and validation failures, and any unexpected record/query content or count.
pub fn run(arguments: impl IntoIterator<Item = String>) -> Result<Value> {
    let options = Options::parse(arguments)?;
    let sequence = "aCgT"
        .chars()
        .cycle()
        .take(options.symbols)
        .collect::<String>();
    let mut metrics = Vec::new();
    let source = plain_workload(&options, &sequence, &mut metrics)?;
    let compressed = bgzf_workload(&options, &sequence, &source, &mut metrics)?;
    blast_workload(&options, &mut metrics)?;
    Ok(json!({
        "kind": "bioferris-qa-workload",
        "schema_version": 1,
        "records": options.records,
        "symbols_per_record": options.symbols,
        "source_bytes": fs::metadata(&source)?.len(),
        "bgzf_bytes": fs::metadata(&compressed)?.len(),
        "verified_plain_records": options.records,
        "verified_bgzf_records": options.records,
        "verified_blast_queries": options.records,
        "process_peak_rss_kib": peak_rss_kib(),
        "phases": metrics,
    }))
}

// The integration harness imports this example to run the identical workload in
// deadline-controlled subprocesses; its entry point is used only by Cargo.
#[allow(dead_code)]
fn main() -> std::process::ExitCode {
    match run(std::env::args().skip(1)) {
        Ok(report) => {
            println!(
                "{}",
                serde_json::to_string_pretty(&report).expect("serializable QA report")
            );
            std::process::ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("qa_workload: {error}");
            std::process::ExitCode::from(2)
        }
    }
}
