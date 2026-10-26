//! Synthetic boundedness checks, isolated process deadlines and host RSS reports.

use std::cell::Cell;
use std::fs::{self, File};
use std::io::{self, BufReader, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, ExitStatus, Stdio};
use std::rc::Rc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

use bioferris::Dna;
use bioferris::io::bgzf::{BgzfWriteOptions, write_bgzf_new};
use bioferris::io::bgzf_index::{BgzfIndexOptions, BgzfSequenceIndex};
use bioferris::io::blast::{BlastReadLimits, BlastTabularOptions, BlastTabularReader};
use bioferris::io::index::{IndexError, IndexOptions, IndexedFormat, SequenceIndex};
use bioferris::io::{FastaReader, IoErrorKind, ReadLimits};
use serde_json::{Value, json};

#[path = "../examples/qa_workload.rs"]
mod qa_workload;

struct Directory(PathBuf);

impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "bioferris-resource-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
}

impl Drop for Directory {
    fn drop(&mut self) {
        // Only a directory this test created successfully and owns exclusively.
        fs::remove_dir_all(&self.0).unwrap();
    }
}

struct OwnedChild(Child);

impl Drop for OwnedChild {
    fn drop(&mut self) {
        if self.0.try_wait().ok().flatten().is_none() {
            let _ = self.0.kill();
            let _ = self.0.wait();
        }
    }
}

struct Outcome {
    status: ExitStatus,
    timed_out: bool,
    elapsed: Duration,
    sampled_peak_rss_kib: Option<u64>,
    stdout: String,
    stderr: String,
}

fn sample_rss(pid: u32) -> Option<u64> {
    #[cfg(target_os = "linux")]
    {
        let status = fs::read_to_string(format!("/proc/{pid}/status")).ok()?;
        status
            .lines()
            .filter_map(|line| {
                if line.starts_with("VmHWM:") || line.starts_with("VmRSS:") {
                    line.split_whitespace().nth(1)?.parse::<u64>().ok()
                } else {
                    None
                }
            })
            .max()
    }
    #[cfg(not(target_os = "linux"))]
    {
        let _ = pid;
        None
    }
}

fn spawn_worker(directory: &Path, mode: &str, records: usize, symbols: usize) -> OwnedChild {
    let stdout = File::options()
        .write(true)
        .create_new(true)
        .open(directory.join("stdout.log"))
        .unwrap();
    let stderr = File::options()
        .write(true)
        .create_new(true)
        .open(directory.join("stderr.log"))
        .unwrap();
    let child = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "resource_worker",
            "--nocapture",
            "--test-threads=1",
        ])
        .env("BIOFERRIS_QA_WORKER_MODE", mode)
        .env("BIOFERRIS_QA_WORKER_DIRECTORY", directory)
        .env("BIOFERRIS_QA_WORKER_RECORDS", records.to_string())
        .env("BIOFERRIS_QA_WORKER_SYMBOLS", symbols.to_string())
        .stdin(Stdio::null())
        .stdout(stdout)
        .stderr(stderr)
        .spawn()
        .unwrap();
    OwnedChild(child)
}

fn wait_with_deadline(mut child: OwnedChild, directory: &Path, timeout: Duration) -> Outcome {
    let started = Instant::now();
    let mut peak = None;
    let mut timed_out = false;
    let status = loop {
        if let Some(rss) = sample_rss(child.0.id()) {
            peak = Some(peak.unwrap_or(0).max(rss));
        }
        if let Some(status) = child.0.try_wait().unwrap() {
            break status;
        }
        if started.elapsed() >= timeout {
            timed_out = true;
            child.0.kill().unwrap();
            break child.0.wait().unwrap();
        }
        std::thread::sleep(Duration::from_millis(2));
    };
    Outcome {
        status,
        timed_out,
        elapsed: started.elapsed(),
        sampled_peak_rss_kib: peak,
        stdout: fs::read_to_string(directory.join("stdout.log")).unwrap(),
        stderr: fs::read_to_string(directory.join("stderr.log")).unwrap(),
    }
}

fn report(outcome: &Outcome) -> Value {
    assert!(
        !outcome.timed_out,
        "worker exceeded deadline: {}",
        outcome.stderr
    );
    assert!(
        outcome.status.success(),
        "worker failed: {}\n{}",
        outcome.stdout,
        outcome.stderr
    );
    let line = outcome
        .stdout
        .lines()
        .find_map(|line| {
            line.split_once("BIOFERRIS_QA_REPORT=")
                .map(|(_, report)| report)
        })
        .expect("child report");
    serde_json::from_str(line).unwrap()
}

#[test]
fn bounded_real_workloads_finish_and_report_separate_phases() {
    for (records, symbols) in [(32, 64), (2_048, 256)] {
        let directory = Directory::new();
        let child = spawn_worker(&directory.0, "workload", records, symbols);
        let outcome = wait_with_deadline(child, &directory.0, Duration::from_secs(30));
        let report = report(&outcome);
        assert_eq!(report["verified_plain_records"], records);
        assert_eq!(report["verified_bgzf_records"], records);
        assert_eq!(report["verified_blast_queries"], records);
        let phases = report["phases"].as_array().unwrap();
        assert_eq!(phases.len(), 12);
        for name in [
            "plain_build",
            "plain_open",
            "plain_lookup",
            "bgzf_build",
            "bgzf_open",
            "bgzf_lookup",
            "blast_read",
        ] {
            assert!(phases.iter().any(|phase| phase["phase"] == name));
        }
        // Measurements are recorded, not compared with a brittle host-specific
        // memory ceiling. API byte/count budgets are tested independently below.
        eprintln!(
            "resource_workload: records={records} symbols={symbols} elapsed_ms={} sampled_peak_rss_kib={:?} worker_peak_rss_kib={}",
            outcome.elapsed.as_millis(),
            outcome.sampled_peak_rss_kib,
            report["process_peak_rss_kib"]
        );
        if cfg!(target_os = "linux") && Path::new("/proc/self/status").exists() {
            assert!(
                report["process_peak_rss_kib"]
                    .as_u64()
                    .is_some_and(|rss| rss > 0)
            );
        }
    }
}

#[test]
fn deadline_terminates_and_reaps_a_ready_owned_worker() {
    let directory = Directory::new();
    let mut child = spawn_worker(&directory.0, "wait", 1, 1);
    let startup = Instant::now();
    loop {
        if fs::read_to_string(directory.0.join("stdout.log"))
            .unwrap()
            .contains("BIOFERRIS_QA_READY")
        {
            break;
        }
        assert!(
            child.0.try_wait().unwrap().is_none(),
            "worker exited before readiness"
        );
        assert!(
            startup.elapsed() < Duration::from_secs(10),
            "worker startup deadline"
        );
        std::thread::sleep(Duration::from_millis(2));
    }
    let outcome = wait_with_deadline(child, &directory.0, Duration::from_millis(50));
    assert!(outcome.timed_out);
    assert!(!outcome.status.success());
    assert!(outcome.elapsed < Duration::from_secs(5));
    assert!(outcome.stdout.contains("BIOFERRIS_QA_READY"));
}

#[test]
fn adversarial_inputs_reject_within_budgets_in_an_isolated_process() {
    let directory = Directory::new();
    let child = spawn_worker(&directory.0, "adversarial", 1, 1);
    let outcome = wait_with_deadline(child, &directory.0, Duration::from_secs(10));
    let report = report(&outcome);
    assert_eq!(report["plain_index_error"], "limit_exceeded");
    assert_eq!(report["bgzf_index_error"], "limit_exceeded");
    assert_eq!(report["blast_error"], "limit_exceeded");
    assert_eq!(report["fasta_error"], "limit_exceeded");
    assert!(report["fasta_transport_bytes"].as_u64().unwrap() <= 128);
    assert_eq!(report["blast_successful_queries"], 0);
    assert_eq!(report["sidecars_created"], 0);
    eprintln!("adversarial_limits: {report}");
}

#[test]
fn workload_argument_bounds_are_checked_before_creating_files() {
    let directory = Directory::new();
    for (records, symbols) in [
        ("0", "4"),
        ("1", "4097"),
        ("100001", "1"),
        ("100000", "4096"),
    ] {
        let error = qa_workload::run([
            "--directory".into(),
            directory.0.to_string_lossy().into_owned(),
            "--records".into(),
            records.into(),
            "--symbols".into(),
            symbols.into(),
        ])
        .unwrap_err();
        assert_eq!(
            error.to_string(),
            "workload bounds: 1..=100000 records, 1..=4096 symbols, product <=16777216"
        );
    }
    assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 0);
}

#[test]
fn workload_files_are_create_only_and_not_implicitly_cleaned() {
    let directory = Directory::new();
    let path = directory.0.join("sequences.fa");
    fs::write(&path, b"sentinel source").unwrap();
    let error = qa_workload::run([
        "--directory".into(),
        directory.0.to_string_lossy().into_owned(),
        "--records".into(),
        "2".into(),
        "--symbols".into(),
        "4".into(),
    ])
    .unwrap_err();
    assert_eq!(
        error.downcast_ref::<io::Error>().unwrap().kind(),
        io::ErrorKind::AlreadyExists
    );
    assert_eq!(fs::read(path).unwrap(), b"sentinel source");
    assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 1);
}

fn fence<'a>(guide: &'a str, marker: &str) -> &'a str {
    let after = guide.split_once(marker).expect("guide marker").1;
    let fenced = after.split_once("```").expect("opening fence").1;
    let body = fenced.split_once('\n').expect("fence language").1;
    body.split_once("```").expect("closing fence").0
}

#[test]
fn documented_literal_inputs_and_report_counts_match_the_real_workload() {
    let guide = include_str!("../../docs/testing/faults-and-resources.md");
    let directory = Directory::new();
    let actual = qa_workload::run([
        "--directory".into(),
        directory.0.to_string_lossy().into_owned(),
        "--records".into(),
        "2".into(),
        "--symbols".into(),
        "4".into(),
    ])
    .unwrap();
    for (name, marker) in [
        ("sequences.fa", "<!-- qa-workload-fasta -->"),
        ("blast.tsv", "<!-- qa-workload-blast -->"),
    ] {
        assert_eq!(
            fs::read_to_string(directory.0.join(name)).unwrap(),
            fence(guide, marker)
        );
    }
    let expected: Value =
        serde_json::from_str(fence(guide, "<!-- qa-workload-report -->")).unwrap();
    for (key, value) in expected.as_object().unwrap() {
        if key != "phases" && key != "process_peak_rss_kib" {
            assert_eq!(&actual[key], value, "documented report field {key}");
        }
    }
    let expected_phases = expected["phases"].as_array().unwrap();
    let actual_phases = actual["phases"].as_array().unwrap();
    assert_eq!(actual_phases.len(), expected_phases.len());
    for (actual, expected) in actual_phases.iter().zip(expected_phases) {
        assert_eq!(actual["phase"], expected["phase"]);
        assert!(actual["elapsed_us"].as_u64().is_some());
        assert!(
            actual["process_peak_rss_kib"].is_null()
                || actual["process_peak_rss_kib"].as_u64().is_some()
        );
    }
    let error = qa_workload::run([
        "--directory".into(),
        directory.0.to_string_lossy().into_owned(),
        "--records".into(),
        "0".into(),
        "--symbols".into(),
        "4".into(),
    ])
    .unwrap_err();
    assert_eq!(
        format!("qa_workload: {error}\n"),
        fence(guide, "<!-- qa-workload-error -->")
    );
}

fn adversarial_blast() -> String {
    // Thousands of HSP/hit rows share one query: reject its fifth retained hit
    // rather than materializing this complete query.
    let mut rows = Vec::new();
    for number in 0..4_096 {
        writeln!(rows, "q\ts{number}").unwrap();
    }
    let mut blast = BlastTabularReader::new(
        Cursor::new(rows),
        BlastTabularOptions {
            fields: Some(vec!["qseqid".into(), "sseqid".into()]),
            limits: BlastReadLimits {
                max_hits_per_query: 4,
                max_hsps_per_query: 4,
                ..BlastReadLimits::default()
            },
            ..BlastTabularOptions::default()
        },
    )
    .unwrap();
    let error = blast.next().unwrap().unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::LimitExceeded);
    assert!(blast.next().is_none(), "error must fuse the reader");
    error.code().to_owned()
}

struct Counted<R> {
    inner: R,
    count: Rc<Cell<usize>>,
}

impl<R: Read> Read for Counted<R> {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        let read = self.inner.read(bytes)?;
        self.count.set(self.count.get() + read);
        Ok(read)
    }
}

fn adversarial_fasta() -> (String, usize) {
    // A virtual 1MiB unterminated physical line is generated on demand, not
    // allocated. Fixed-size transport read-ahead is counted separately.
    let count = Rc::new(Cell::new(0));
    let transport = Counted {
        inner: Cursor::new(b">id\n").chain(io::repeat(b'A').take(1024 * 1024)),
        count: Rc::clone(&count),
    };
    let mut fasta = FastaReader::<_, Dna>::with_limits(
        BufReader::with_capacity(64, transport),
        ReadLimits {
            max_line_bytes: 32,
            max_record_symbols: 64,
            max_total_bytes: 128,
        },
    );
    let error = fasta.next().unwrap().unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::LimitExceeded);
    assert!(fasta.next().is_none());
    assert!(
        count.get() <= 128,
        "bounded parser must not drain its virtual 1MiB input"
    );
    (error.code().to_owned(), count.get())
}

fn adversarial(directory: &Path) -> Value {
    const FASTA: &[u8] = b">first\naCgT\n>second\nNN\n";
    let source = directory.join("bounded.fa");
    File::options()
        .write(true)
        .create_new(true)
        .open(&source)
        .unwrap()
        .write_all(FASTA)
        .unwrap();
    let plain_error = SequenceIndex::<Dna>::build(
        &source,
        IndexedFormat::Fasta,
        IndexOptions {
            max_records: 1,
            ..IndexOptions::default()
        },
    )
    .unwrap_err();
    assert!(matches!(
        plain_error,
        IndexError::LimitExceeded {
            resource: "index_records",
            limit: 1,
            ..
        }
    ));
    let compressed = directory.join("bounded.fa.bgz");
    write_bgzf_new(&compressed, BgzfWriteOptions::default(), |output| {
        output.write_all(FASTA)
    })
    .unwrap();
    let bgzf_error = BgzfSequenceIndex::<Dna>::build(
        &compressed,
        IndexedFormat::Fasta,
        BgzfIndexOptions {
            max_compressed_bytes: 1,
            ..BgzfIndexOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(bgzf_error.code(), "limit_exceeded");

    let blast_error = adversarial_blast();
    let (fasta_error, transport_bytes) = adversarial_fasta();
    let sidecars = fs::read_dir(directory)
        .unwrap()
        .filter_map(Result::ok)
        .filter(|entry| {
            entry
                .path()
                .extension()
                .is_some_and(|extension| extension == "json")
        })
        .count();
    json!({
        "plain_index_error": plain_error.code(),
        "bgzf_index_error": bgzf_error.code(),
        "blast_error": blast_error,
        "fasta_error": fasta_error,
        "fasta_transport_bytes": transport_bytes,
        "blast_successful_queries": 0,
        "sidecars_created": sidecars,
        "process_peak_rss_kib": qa_workload::peak_rss_kib(),
    })
}

#[test]
fn resource_worker() {
    let Some(mode) = std::env::var_os("BIOFERRIS_QA_WORKER_MODE") else {
        return;
    };
    let directory = PathBuf::from(std::env::var_os("BIOFERRIS_QA_WORKER_DIRECTORY").unwrap());
    if mode == "wait" {
        println!("BIOFERRIS_QA_READY");
        io::stdout().flush().unwrap();
        std::thread::sleep(Duration::from_secs(30));
        panic!("parent should terminate the ready child before this point");
    }
    let report = if mode == "adversarial" {
        adversarial(&directory)
    } else {
        assert_eq!(mode, "workload");
        qa_workload::run([
            "--directory".into(),
            directory.to_string_lossy().into_owned(),
            "--records".into(),
            std::env::var("BIOFERRIS_QA_WORKER_RECORDS").unwrap(),
            "--symbols".into(),
            std::env::var("BIOFERRIS_QA_WORKER_SYMBOLS").unwrap(),
        ])
        .unwrap()
    };
    println!("BIOFERRIS_QA_REPORT={report}");
}
