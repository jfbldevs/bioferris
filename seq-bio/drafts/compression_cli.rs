//! Actual-process compression adapter checks; promote only after public wiring.

use std::fs;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

use serde_json::Value;

struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "bioferris-compression-cli-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().unwrap().to_owned()
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}
fn binary() -> PathBuf {
    if let Some(path) = std::env::var_os("BIOFERRIS_COMPRESSION_BINARY") {
        return PathBuf::from(path);
    }
    match option_env!("CARGO_BIN_EXE_bioferris") {
        Some(path) => PathBuf::from(path),
        None => panic!("set BIOFERRIS_COMPRESSION_BINARY for the isolated adapter harness"),
    }
}
fn run(args: &[&str]) -> Output {
    Command::new(binary()).args(args).output().unwrap()
}
fn success(output: &Output) -> Value {
    assert!(
        output.status.success(),
        "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    let value: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["schema_version"], 1);
    assert_eq!(value["ok"], true);
    value["result"].clone()
}
fn failure(output: &Output, code: &str) -> Value {
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stderr.is_empty());
    let value: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(value["schema_version"], 1);
    assert_eq!(value["ok"], false);
    assert_eq!(value["error"]["code"], code);
    value["error"].clone()
}
fn compress(directory: &Directory, compression: &str, input: &[u8]) -> String {
    let source = directory.path("input.raw");
    let output = directory.path(compression);
    fs::write(&source, input).unwrap();
    success(&run(&[
        "compress-file",
        "--compression",
        compression,
        "--input",
        &source,
        "--output",
        &output,
        "--json",
    ]));
    output
}

#[test]
fn gzip_and_bgzf_preserve_all_bytes_and_report_exact_counts() {
    let input: Vec<u8> = (0..160_000)
        .map(|index| u8::try_from(index % 256).unwrap())
        .collect();
    for compression in ["gzip", "bgzf"] {
        let directory = Directory::new();
        let source = directory.path("source");
        let target = directory.path("target");
        fs::write(&source, &input).unwrap();
        let result = success(&run(&[
            "compress-file",
            "--compression",
            compression,
            "--input",
            &source,
            "--output",
            &target,
            "--json",
        ]));
        assert_eq!(result["input_bytes"], input.len());
        assert_eq!(result["output_bytes"], fs::metadata(&target).unwrap().len());
        let encoded = fs::read(&target).unwrap();
        let mut decoded = Vec::new();
        flate2::read::MultiGzDecoder::new(encoded.as_slice())
            .read_to_end(&mut decoded)
            .unwrap();
        assert_eq!(decoded, input);
        assert_eq!(encoded[3], if compression == "gzip" { 0 } else { 4 });
        let second = directory.path("second");
        success(&run(&[
            "compress-file",
            "--compression",
            compression,
            "--input",
            &source,
            "--output",
            &second,
            "--json",
        ]));
        assert_eq!(encoded, fs::read(second).unwrap());
    }
}

#[test]
fn stdin_empty_and_exact_zero_limits_and_human_report_work() {
    for compression in ["gzip", "bgzf"] {
        let directory = Directory::new();
        let target = directory.path("target");
        let output = Command::new(binary())
            .args([
                "compress-file",
                "--compression",
                compression,
                "--input",
                "-",
                "--output",
                &target,
                "--max-input-bytes",
                "0",
                "--json",
            ])
            .stdin(Stdio::null())
            .output()
            .unwrap();
        assert_eq!(success(&output)["input_bytes"], 0);
        let source = directory.path("source");
        fs::write(&source, b"A").unwrap();
        let human = directory.path("human");
        let output = run(&[
            "compress-file",
            "--compression",
            compression,
            "--input",
            &source,
            "--output",
            &human,
            "--max-input-bytes",
            "1",
        ]);
        assert!(output.status.success());
        assert!(output.stderr.is_empty());
        assert!(
            String::from_utf8(output.stdout)
                .unwrap()
                .starts_with("Compressed 1 bytes into ")
        );
    }
}

#[test]
fn bounds_and_existing_destinations_fail_without_named_partial_outputs() {
    for compression in ["gzip", "bgzf"] {
        let directory = Directory::new();
        let source = directory.path("source");
        let target = directory.path("target");
        fs::write(&source, b"ACGT").unwrap();
        let error = failure(
            &run(&[
                "compress-file",
                "--compression",
                compression,
                "--input",
                &source,
                "--output",
                &target,
                "--max-input-bytes",
                "3",
                "--json",
            ]),
            "limit_exceeded",
        );
        assert_eq!(error["details"]["resource"], "input_bytes");
        assert!(!PathBuf::from(&target).exists());
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 1);
        fs::write(&target, b"preserve").unwrap();
        failure(
            &run(&[
                "compress-file",
                "--compression",
                compression,
                "--input",
                &source,
                "--output",
                &target,
                "--json",
            ]),
            "output_error",
        );
        assert_eq!(fs::read(&target).unwrap(), b"preserve");
        assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 2);
    }
}

#[test]
fn malformed_repeated_and_inapplicable_flags_fail_before_publication() {
    let directory = Directory::new();
    let source = directory.path("source");
    let target = directory.path("target");
    fs::write(&source, b"ACGT").unwrap();
    for flags in [
        vec!["--level", "-1"],
        vec!["--level", "1.5"],
        vec!["--level", "1", "--level", "2"],
        vec!["--block-size", "20"],
        vec!["--max-blocks", "20"],
        vec!["--unknown", "1"],
        vec!["--max-input-bytes", "18446744073709551616"],
        vec!["--json"],
    ] {
        let mut args = vec![
            "compress-file",
            "--compression",
            "gzip",
            "--input",
            &source,
            "--output",
            &target,
            "--json",
        ];
        args.extend(flags);
        failure(&run(&args), "invalid_arguments");
        assert!(!PathBuf::from(&target).exists());
    }
    failure(
        &run(&[
            "compress-file",
            "--compression",
            "gzip",
            "--input",
            &source,
            "--output",
            &target,
            "--level",
            "10",
            "--json",
        ]),
        "invalid_gzip_options",
    );
    failure(
        &run(&[
            "compress-file",
            "--compression",
            "bgzf",
            "--input",
            &source,
            "--output",
            &target,
            "--block-size",
            "0",
            "--json",
        ]),
        "invalid_bgzf_options",
    );
    failure(
        &run(&[
            "compress-file",
            "--compression",
            "bgzf",
            "--input",
            &source,
            "--output",
            "-",
            "--json",
        ]),
        "invalid_arguments",
    );
}

#[test]
fn indexed_fastq_and_fasta_preserve_the_established_record_json_and_human_output() {
    for (format, input, human, sequence, qualities) in [
        (
            "fasta",
            b">r description\r\naCgT\r\n>empty\r\n".as_slice(),
            b">r description\naCgT\n".as_slice(),
            "ACGT",
            Value::Null,
        ),
        (
            "fastq-sanger",
            b"@r description\r\naCgT\r\n+\r\n!5I~\r\n",
            b"@r description\naCgT\n+\n!5I~\n",
            "ACGT",
            serde_json::json!([0, 20, 40, 93]),
        ),
    ] {
        let directory = Directory::new();
        let source = compress(&directory, "bgzf", input);
        let sidecar = directory.path("index.json");
        let result = success(&run(&[
            "index-bgzf",
            "dna",
            "--format",
            format,
            "--input",
            &source,
            "--output",
            &sidecar,
            "--json",
        ]));
        assert_eq!(result["compression"], "bgzf");
        assert_eq!(result["uncompressed_bytes"], input.len());
        assert_eq!(result["blocks"], 2);
        let args = [
            "get-bgzf-record",
            "dna",
            "r",
            "--input",
            &source,
            "--index",
            &sidecar,
            "--json",
        ];
        let result = success(&run(&args));
        assert_eq!(result["record"]["header"], "r description");
        assert_eq!(result["record"]["source_sequence"], "aCgT");
        assert_eq!(result["record"]["sequence"], sequence);
        assert_eq!(result["record"]["qualities"], qualities);
        let output = run(&args[..args.len() - 1]);
        assert!(output.status.success());
        assert!(output.stderr.is_empty());
        assert_eq!(output.stdout, human);
        failure(
            &run(&[
                "get-bgzf-record",
                "dna",
                "absent",
                "--input",
                &source,
                "--index",
                &sidecar,
                "--json",
            ]),
            "record_not_found",
        );
        failure(
            &run(&[
                "index-bgzf",
                "dna",
                "--format",
                format,
                "--input",
                &source,
                "--output",
                &sidecar,
                "--json",
            ]),
            "output_error",
        );
    }
}

#[test]
fn bgzf_indexes_reject_ordinary_gzip_inapplicable_options_and_forged_sidecars() {
    let directory = Directory::new();
    let gzip = compress(&directory, "gzip", b">r\nAC\n");
    let sidecar = directory.path("index.json");
    failure(
        &run(&[
            "index-bgzf",
            "dna",
            "--format",
            "fasta",
            "--input",
            &gzip,
            "--output",
            &sidecar,
            "--json",
        ]),
        "invalid_bgzf_header",
    );
    let source = compress(&directory, "bgzf", b">r\nAC\n");
    for flags in [
        vec!["--level", "6"],
        vec!["--block-size", "100"],
        vec!["--max-input-bytes", "100"],
        vec!["--compression", "bgzf"],
    ] {
        let mut args = vec![
            "index-bgzf",
            "dna",
            "--format",
            "fasta",
            "--input",
            &source,
            "--output",
            &sidecar,
            "--json",
        ];
        args.extend(flags);
        failure(&run(&args), "invalid_arguments");
    }
    for (flag, resource) in [
        ("--max-blocks", "bgzf_blocks"),
        ("--max-records", "index_records"),
        ("--max-id-bytes", "index_id_bytes"),
        ("--max-index-bytes", "index_bytes"),
        ("--max-compressed-bytes", "compressed_source_bytes"),
        ("--max-total-bytes", "source_bytes"),
    ] {
        let error = failure(
            &run(&[
                "index-bgzf",
                "dna",
                "--format",
                "fasta",
                "--input",
                &source,
                "--output",
                &sidecar,
                flag,
                "0",
                "--json",
            ]),
            "limit_exceeded",
        );
        assert_eq!(error["details"]["resource"], resource);
    }
    success(&run(&[
        "index-bgzf",
        "dna",
        "--format",
        "fasta",
        "--input",
        &source,
        "--output",
        &sidecar,
        "--json",
    ]));
    let mut value: Value = serde_json::from_slice(&fs::read(&sidecar).unwrap()).unwrap();
    value["entries"][0]["virtual_start"] = 1.into();
    fs::write(&sidecar, serde_json::to_vec(&value).unwrap()).unwrap();
    failure(
        &run(&[
            "get-bgzf-record",
            "dna",
            "r",
            "--input",
            &source,
            "--index",
            &sidecar,
            "--json",
        ]),
        "invalid_index",
    );
}

#[cfg(unix)]
#[test]
fn dangling_symlink_destinations_and_invalid_output_directories_are_preserved() {
    let directory = Directory::new();
    let source = directory.path("source");
    let target = directory.path("dangling");
    fs::write(&source, b"ACGT").unwrap();
    std::os::unix::fs::symlink("absent", &target).unwrap();
    for compression in ["gzip", "bgzf"] {
        failure(
            &run(&[
                "compress-file",
                "--compression",
                compression,
                "--input",
                &source,
                "--output",
                &target,
                "--json",
            ]),
            "output_error",
        );
        assert!(
            fs::symlink_metadata(&target)
                .unwrap()
                .file_type()
                .is_symlink()
        );
        failure(
            &run(&[
                "compress-file",
                "--compression",
                compression,
                "--input",
                &source,
                "--output",
                &directory.path("absent/target"),
                "--json",
            ]),
            "io_error",
        );
    }
}

#[test]
fn stdin_copy_exceeding_its_budget_aborts_publication() {
    let directory = Directory::new();
    let target = directory.path("target");
    let mut child = Command::new(binary())
        .args([
            "compress-file",
            "--compression",
            "bgzf",
            "--input",
            "-",
            "--output",
            &target,
            "--max-input-bytes",
            "2",
            "--json",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(b"ACG").unwrap();
    failure(&child.wait_with_output().unwrap(), "limit_exceeded");
    assert!(!PathBuf::from(&target).exists());
    assert_eq!(fs::read_dir(&directory.0).unwrap().count(), 0);
}
