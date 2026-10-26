//! Independent, writer-readback and archived-response checks of graph/result imports.

use std::io::{self, BufReader, Cursor, Read, Write};

use bioferris::io::graph::{parse_dot, parse_dot_reader, parse_gfa, parse_gfa_reader};
use bioferris::io::result_json::{
    AlignmentResultKind, CliResultKind, ImportedResult, parse_result_json, parse_result_json_reader,
};
use bioferris::io::{BioIoError, IoErrorKind, write_poa_dot_to, write_poa_gfa_to};
use bioferris::prelude::*;
use flate2::{Compression, write::GzEncoder};
use serde_json::{Value, json};

const GFA: &str = "H\tVN:Z:1.0\nS\tn0\tA\tCL:i:0\nS\tn1\tC\tCL:i:1\nL\tn0\t+\tn1\t+\t0M\tRC:i:1\nP\tread1\tn0+,n1+\t*\n";
const DOT: &str = "digraph POASTA {\n  rankdir=LR;\n  n0 [label=\"n0\\nA\\nclass 0\"];\n  n1 [label=\"n1\\nC\\nclass 1\"];\n  n0 -> n1 [label=\"1\"];\n}\n";
const DISTANCE: &str = r#"{"schema":"align-output-v1","kind":"myers_distance","algorithm":"Myers bit-vector","mode":"global unit edit distance","alphabet":"DNA","pattern":{"length":4,"sequence":"ACGT"},"text":{"length":3,"sequence":"AGT"},"distance":1}"#;

fn native(input: &str, kind: AlignmentResultKind) -> Value {
    let ImportedResult::Alignment(document) = parse_result_json(input).unwrap() else {
        panic!("expected a native result");
    };
    assert_eq!(document.kind(), kind);
    assert_eq!(
        document.value(),
        &serde_json::from_str::<Value>(input).unwrap()
    );
    document.value().clone()
}

fn write_json(write: impl FnOnce(&mut Vec<u8>) -> bioferris::align::Result<()>) -> String {
    let mut bytes = Vec::new();
    write(&mut bytes).unwrap();
    String::from_utf8(bytes).unwrap()
}

fn records() -> Vec<SequenceRecord<Dna>> {
    [
        ("reference one", "ACGT"),
        ("sample_2", "AGT"),
        ("deletion", "ACT"),
    ]
    .into_iter()
    .map(|(id, input)| SequenceRecord::new(id, dna(input).unwrap()).unwrap())
    .collect()
}

fn error_kind<T: std::fmt::Debug>(result: Result<T, BioIoError>, kind: IoErrorKind) {
    let error = result.unwrap_err();
    assert_eq!(error.kind(), kind, "{error}");
}

#[test]
fn independently_specified_gfa_preserves_paths_classes_and_support() {
    let graph = parse_gfa::<Dna>(GFA).unwrap();
    assert_eq!(graph.nodes().len(), 2);
    assert_eq!(graph.class_count(), 2);
    assert_eq!(graph.nodes()[1].id(), 1);
    assert_eq!(graph.nodes()[1].symbol().as_ascii(), b'C');
    assert_eq!(graph.nodes()[1].class(), 1);
    assert_eq!(graph.edges()[0].source(), 0);
    assert_eq!(graph.edges()[0].target(), 1);
    assert_eq!(graph.edges()[0].support(), 1);
    assert_eq!(graph.paths()[0].name(), "read1");
    assert_eq!(graph.paths()[0].nodes(), [0, 1]);
    assert_eq!(graph.path_sequence(0).as_deref(), Some("AC"));
    assert_eq!(graph.path_sequence(1), None);
    let crlf = format!("# comment\r\n{}\r\n", GFA.replace('\n', "\r\n"));
    assert_eq!(parse_gfa::<Dna>(&crlf).unwrap(), graph);
    let lowercase = GFA.replace("\tA\t", "\ta\t");
    assert_eq!(parse_gfa::<Dna>(&lowercase).unwrap(), graph);
    let empty = parse_gfa::<Dna>("H\tVN:Z:1.0\nP\tempty\t*\t*\n").unwrap();
    assert_eq!(empty.path_sequence(0).as_deref(), Some(""));
    assert!(empty.nodes().is_empty());
}

#[test]
fn own_graph_writers_read_back_without_reconstructing_missing_metadata() {
    let run = Graph::poasta().run(&records()).unwrap();
    let source = run.result().graph();
    let mut gfa = Vec::new();
    write_poa_gfa_to(&mut gfa, source).unwrap();
    let imported = parse_gfa_reader::<Dna>(gfa.as_slice()).unwrap();
    let mut dot = Vec::new();
    write_poa_dot_to(&mut dot, source).unwrap();
    let topology = parse_dot_reader::<Dna>(dot.as_slice()).unwrap();
    assert_eq!(imported.nodes(), topology.nodes());
    assert_eq!(imported.edges(), topology.edges());
    assert_eq!(imported.class_count(), source.class_count());
    assert_eq!(imported.nodes().len(), source.node_count());
    assert_eq!(imported.edges().len(), source.edge_count());
    assert_eq!(imported.paths().len(), source.path_count());
    for (imported, original) in imported.nodes().iter().zip(source.nodes()) {
        assert_eq!(imported.id(), original.id().index());
        assert_eq!(imported.class(), original.class().index());
        assert_eq!(imported.symbol(), original.symbol());
    }
    for (imported, original) in imported.edges().iter().zip(source.edges()) {
        assert_eq!(imported.source(), original.source().index());
        assert_eq!(imported.target(), original.target().index());
        assert_eq!(imported.support(), original.support());
    }
    for (index, record) in records().iter().enumerate() {
        assert_eq!(
            imported.path_sequence(index).unwrap(),
            record.sequence().to_string()
        );
    }
    assert_eq!(imported.paths()[0].name(), "reference%20one");
}

#[test]
fn gfa_rejects_unsupported_semantics_instead_of_discarding_them() {
    for input in [
        GFA.replace("VN:Z:1.0", "VN:Z:1.1"),
        GFA.replace("\tn0\tA\t", "\tn0\tAC\t"),
        GFA.replace("L\tn0\t+", "L\tn0\t-"),
        GFA.replace("0M", "1M"),
        GFA.replace("CL:i:0", "CL:i:0\tXX:Z:metadata"),
        GFA.replace("n0+,n1+", "n0+,n1-"),
        format!("{GFA}W\tsample\t0\tread\t0\t2\t>n0>n1\n"),
    ] {
        error_kind(parse_gfa::<Dna>(&input), IoErrorKind::Unsupported);
    }
}

#[test]
fn gfa_rejects_bad_names_topology_path_coverage_and_edge_support() {
    for input in [
        GFA.replace("RC:i:1", "RC:i:2"),
        GFA.replace("RC:i:1", "RC:i:0"),
        GFA.replace("CL:i:1", "CL:i:0"),
        GFA.replace("CL:i:1", "CL:i:2"),
        GFA.replace("n0+,n1+", "n0+,n9+"),
        GFA.replace("n0+,n1+", "n1+,n0+"),
        GFA.replace("n0+,n1+", "n0+"),
        GFA.replace("read1", "n0"),
        format!("{GFA}P\tread1\tn0+,n1+\t*\n"),
        GFA.replace("P\tread1", "L\tn1\t+\tn0\t+\t0M\tRC:i:1\nP\tread1"),
        GFA.replace("L\tn0\t+\tn1\t+\t0M\tRC:i:1\n", ""),
    ] {
        error_kind(parse_gfa::<Dna>(&input), IoErrorKind::InvalidStructure);
    }
    error_kind(parse_gfa::<Dna>(""), IoErrorKind::InvalidHeader);
    let error = parse_gfa::<Dna>(&GFA.replace("\tA\t", "\t!\t")).unwrap_err();
    assert_eq!(error.kind(), IoErrorKind::InvalidSequence);
    assert_eq!(error.record(), Some(1));
    assert_eq!(error.line(), Some(2));
}

#[test]
fn dot_preserves_only_the_topology_and_rejects_general_graphviz_syntax() {
    let topology = parse_dot::<Dna>(DOT).unwrap();
    let gfa = parse_gfa::<Dna>(GFA).unwrap();
    assert_eq!(topology.nodes(), gfa.nodes());
    assert_eq!(topology.edges(), gfa.edges());
    assert_eq!(topology.class_count(), 2);
    error_kind(
        parse_dot::<Dna>(DOT.trim_end_matches("}\n")),
        IoErrorKind::Truncated,
    );
    error_kind(
        parse_dot::<Dna>(&DOT.replace(" [label=", " [color=red, label=")),
        IoErrorKind::Unsupported,
    );
    error_kind(
        parse_dot::<Dna>(&DOT.replace("POASTA", "arbitrary")),
        IoErrorKind::Unsupported,
    );
    error_kind(
        parse_dot::<Dna>(&format!("{DOT}n9;\n")),
        IoErrorKind::InvalidStructure,
    );
    error_kind(
        parse_dot::<Dna>(&DOT.replace("n0 -> n1", "n0 -> n9")),
        IoErrorKind::InvalidStructure,
    );
    error_kind(
        parse_dot::<Dna>(&DOT.replace("n0 -> n1", "n0 -> n0")),
        IoErrorKind::InvalidStructure,
    );
}

#[test]
fn graph_and_json_streams_report_io_errors_and_document_limits() {
    struct Broken;
    impl Read for Broken {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::other("injected read failure"))
        }
    }
    error_kind(
        parse_gfa_reader::<Dna>(BufReader::new(Broken)),
        IoErrorKind::Io,
    );
    error_kind(
        parse_dot_reader::<Dna>(BufReader::new(Broken)),
        IoErrorKind::Io,
    );
    error_kind(
        parse_result_json_reader(BufReader::new(Broken)),
        IoErrorKind::Io,
    );
    let oversized = " ".repeat(64 * 1024 * 1024 + 1);
    error_kind(parse_gfa::<Dna>(&oversized), IoErrorKind::LimitExceeded);
    error_kind(parse_dot::<Dna>(&oversized), IoErrorKind::LimitExceeded);
    error_kind(parse_result_json(&oversized), IoErrorKind::LimitExceeded);
    error_kind(
        parse_gfa_reader::<Dna>(Cursor::new(oversized.as_bytes())),
        IoErrorKind::LimitExceeded,
    );
}

#[test]
fn graph_and_result_path_readers_accept_gzip_by_signature() {
    let directory =
        std::env::temp_dir().join(format!("seq-bio-graph-import-{}", std::process::id()));
    std::fs::create_dir(&directory).unwrap();
    let entries = [
        ("graph.data", GFA),
        ("topology.data", DOT),
        ("result.data", DISTANCE),
    ];
    for (name, text) in entries {
        let path = directory.join(name);
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(text.as_bytes()).unwrap();
        std::fs::write(&path, encoder.finish().unwrap()).unwrap();
    }
    assert_eq!(
        bioferris::io::graph::read_gfa::<Dna>(directory.join("graph.data")).unwrap(),
        parse_gfa::<Dna>(GFA).unwrap()
    );
    assert_eq!(
        bioferris::io::graph::read_dot::<Dna>(directory.join("topology.data")).unwrap(),
        parse_dot::<Dna>(DOT).unwrap()
    );
    assert_eq!(
        bioferris::io::result_json::read_result_json(directory.join("result.data")).unwrap(),
        parse_result_json(DISTANCE).unwrap()
    );
    for (name, _) in entries {
        std::fs::remove_file(directory.join(name)).unwrap();
    }
    std::fs::remove_dir(directory).unwrap();
}

#[test]
fn independent_json_fixture_checks_bounds_without_claiming_optimality() {
    assert_eq!(
        native(DISTANCE, AlignmentResultKind::MyersDistance)["distance"],
        1
    );
    let mut value: Value = serde_json::from_str(DISTANCE).unwrap();
    value["distance"] = json!(2);
    // A plausible number in a file is not proof of the optimal edit distance.
    assert!(parse_result_json(&value.to_string()).is_ok());
    for distance in [0, 5, -1] {
        value["distance"] = json!(distance);
        error_kind(
            parse_result_json(&value.to_string()),
            IoErrorKind::InvalidSchema,
        );
    }
    value["distance"] = json!(1);
    value["text"]["length"] = json!(4);
    error_kind(
        parse_result_json(&value.to_string()),
        IoErrorKind::InvalidSchema,
    );
}

#[test]
fn all_six_native_json_writer_shapes_are_accepted_and_preserved() {
    let first = dna("ACGT").unwrap();
    let second = dna("AGT").unwrap();
    let records = records();
    let graph = Graph::poasta().run(&records).unwrap();
    let cases = [
        (
            write_json(|writer| {
                Pairwise::gotoh()
                    .run(&first, &second)?
                    .write_json_to(writer)
            }),
            AlignmentResultKind::PairwiseAlignment,
        ),
        (
            write_json(|writer| {
                Pairwise::wfa()
                    .run_penalty(&first, &second)?
                    .write_json_to(writer)
            }),
            AlignmentResultKind::PairwiseValue,
        ),
        (
            write_json(|writer| {
                Multiple::mafft(MafftStrategy::FftNs1)
                    .run(&records)?
                    .write_json_to(writer)
            }),
            AlignmentResultKind::MultipleAlignment,
        ),
        (
            write_json(|writer| {
                Pairwise::myers()
                    .run_distance(&first, &second)?
                    .write_json_to(writer)
            }),
            AlignmentResultKind::MyersDistance,
        ),
        (
            write_json(|writer| {
                Pairwise::myers()
                    .run_find_all(&second, &first, 1)?
                    .write_json_to(writer)
            }),
            AlignmentResultKind::MyersSearch,
        ),
        (
            write_json(|writer| {
                Graph::poasta()
                    .run_query(graph.result().graph(), &first)?
                    .write_json_to(writer)
            }),
            AlignmentResultKind::PoastaGraphQuery,
        ),
    ];
    for (text, kind) in cases {
        native(&text, kind);
    }
}

#[test]
fn pairwise_readback_covers_local_affine_distance_and_penalty_objectives() {
    let first = dna("TTACCG").unwrap();
    let second = dna("ACG").unwrap();
    let cases = [
        write_json(|writer| {
            Pairwise::needleman_wunsch()
                .run(&first, &second)?
                .write_json_to(writer)
        }),
        write_json(|writer| {
            Pairwise::smith_waterman()
                .run(&first, &second)?
                .write_json_to(writer)
        }),
        write_json(|writer| {
            Pairwise::gotoh()
                .run(&first, &second)?
                .write_json_to(writer)
        }),
        write_json(|writer| {
            Pairwise::myers()
                .run(&first, &second)?
                .write_json_to(writer)
        }),
        write_json(|writer| Pairwise::wfa().run(&first, &second)?.write_json_to(writer)),
        write_json(|writer| {
            Pairwise::adaptive_affine()
                .run(&first, &second)?
                .write_json_to(writer)
        }),
        write_json(|writer| {
            Pairwise::xdrop()
                .run(&first, &second)?
                .write_json_to(writer)
        }),
    ];
    for text in cases {
        native(&text, AlignmentResultKind::PairwiseAlignment);
    }
}

#[test]
fn json_pairwise_shape_trace_coordinates_and_statistics_are_checked() {
    let first = dna("ACGT").unwrap();
    let second = dna("AGT").unwrap();
    let text = write_json(|writer| {
        Pairwise::gotoh()
            .run(&first, &second)?
            .write_json_to(writer)
    });
    let source: Value = serde_json::from_str(&text).unwrap();
    for (pointer, replacement) in [
        ("/statistics/columns", json!(99)),
        ("/statistics/identity", json!(0.9)),
        ("/statistics/deletion_runs", json!(99)),
        ("/cigar", json!("4=")),
        ("/alignment/second", json!("---")),
        ("/ranges/first/end", json!(99)),
        ("/native_alignment_score", json!(99)),
        ("/inputs/first/sequence", json!("ACG!")),
        ("/objective/direction", json!("minimize")),
    ] {
        let mut altered = source.clone();
        *altered.pointer_mut(pointer).unwrap() = replacement;
        error_kind(
            parse_result_json(&altered.to_string()),
            IoErrorKind::InvalidSchema,
        );
    }
}

#[test]
fn multiple_and_graph_json_reject_inconsistent_dimensions_and_traces() {
    let records = records();
    let graph = Graph::poasta().run(&records).unwrap();
    let multiple = write_json(|writer| graph.write_json_to(writer));
    let multiple: Value = serde_json::from_str(&multiple).unwrap();
    for (pointer, replacement) in [
        ("/statistics/sequences", json!(9)),
        ("/statistics/gap_characters", json!(999)),
        ("/alignment/0/sequence", json!("A")),
        ("/alignment/0/identifier", json!("sample_2")),
    ] {
        let mut altered = multiple.clone();
        *altered.pointer_mut(pointer).unwrap() = replacement;
        error_kind(
            parse_result_json(&altered.to_string()),
            IoErrorKind::InvalidSchema,
        );
    }
    let query = dna("ACGT").unwrap();
    let poasta = write_json(|writer| {
        Graph::poasta()
            .run_query(graph.result().graph(), &query)?
            .write_json_to(writer)
    });
    let poasta: Value = serde_json::from_str(&poasta).unwrap();
    for (pointer, replacement) in [
        ("/graph/classes", json!(999)),
        ("/costs/gap_opening", json!(0)),
        ("/costs/mismatch", json!(-1)),
        ("/operations/0/node", json!(999)),
        ("/operations/0/query_index", json!(3)),
        ("/trace", json!("999A")),
        ("/aligned_pairs/0/node", json!(999)),
        ("/search/discovered_states", json!(999_999)),
    ] {
        let mut altered = poasta.clone();
        *altered.pointer_mut(pointer).unwrap() = replacement;
        error_kind(
            parse_result_json(&altered.to_string()),
            IoErrorKind::InvalidSchema,
        );
    }
}

#[test]
fn json_rejects_duplicate_keys_requests_versions_and_unknown_fields() {
    for text in [
        DISTANCE.replace("\"distance\":1", "\"distance\":1,\"distance\":2"),
        DISTANCE.replace("\"distance\":1", r#""distance":1,"\u0064istance":2"#),
        DISTANCE.replace("\"length\":4", "\"length\":4,\"length\":3"),
        DISTANCE.replace("\"distance\":1", "\"distance\":1,\"extra\":true"),
        DISTANCE.replace("align-output-v1", "align-output-v2"),
        format!("{DISTANCE}{{}}"),
        "{\"schema_version\":1,\"operation\":\"complement\",\"sequence\":\"AC\"}".to_owned(),
        "{\"schema_version\":1,\"schema\":\"align-output-v1\"}".to_owned(),
    ] {
        error_kind(parse_result_json(&text), IoErrorKind::InvalidSchema);
    }
    error_kind(
        parse_result_json(&DISTANCE.replace("myers_distance", "future_result")),
        IoErrorKind::Unsupported,
    );
    error_kind(
        parse_result_json(
            r#"{"schema_version":1,"ok":true,"result":{"operation":"file-summary"}}"#,
        ),
        IoErrorKind::Unsupported,
    );
}

/// A response written by the removed CLI, kept as a test fixture.
fn archived(path: &str) -> String {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap();
    std::fs::read_to_string(root.join(path)).unwrap()
}

fn transform(
    operation: &str,
    molecule: &str,
    input: &str,
    output_molecule: &str,
    sequence: &str,
) -> String {
    json!({
        "schema_version": 1,
        "ok": true,
        "result": {
            "operation": operation,
            "molecule": molecule,
            "input": input,
            "input_length": input.len(),
            "output_molecule": output_molecule,
            "sequence": sequence,
            "length": sequence.len(),
        },
    })
    .to_string()
}

const ARCHIVED_TRANSLATION: &str =
    "seq-bio/tests/fixtures/removed-cli-responses/release-translation.json";
const ARCHIVED_PAIRWISE: &str =
    "seq-bio/tests/fixtures/removed-cli-responses/release-pairwise.json";
const ARCHIVED_ERROR: &str =
    "seq-bio/tests/fixtures/removed-cli-responses/protein-unsupported-release.json";

#[test]
fn removed_cli_responses_import_separately_from_requests() {
    for (text, kind) in [
        (archived(ARCHIVED_TRANSLATION), CliResultKind::Translation),
        (archived(ARCHIVED_PAIRWISE), CliResultKind::Pairwise),
        (archived(ARCHIVED_ERROR), CliResultKind::Error),
        (
            transform("complement", "dna", "ACGT", "dna", "TGCA"),
            CliResultKind::SequenceTransform,
        ),
        (
            transform("reverse-complement", "rna", "AUGC", "rna", "GCAU"),
            CliResultKind::SequenceTransform,
        ),
        (
            transform("transcribe", "dna", "ACGT", "rna", "ACGU"),
            CliResultKind::SequenceTransform,
        ),
        (
            transform("back-transcribe", "rna", "ACGU", "dna", "ACGT"),
            CliResultKind::SequenceTransform,
        ),
    ] {
        let ImportedResult::Cli(document) = parse_result_json_reader(text.as_bytes()).unwrap()
        else {
            panic!("expected a CLI response");
        };
        assert_eq!(document.kind(), kind);
        assert_eq!(
            document.value(),
            &serde_json::from_str::<Value>(&text).unwrap()
        );
    }
}

#[test]
fn cli_import_recomputes_transforms_and_translation_and_checks_nested_pairwise_inputs() {
    let pairwise = archived(ARCHIVED_PAIRWISE);
    for (text, pointer, replacement) in [
        (
            transform("complement", "dna", "ACG", "dna", "TGC"),
            "/result/sequence",
            json!("ACG"),
        ),
        (
            archived(ARCHIVED_TRANSLATION),
            "/result/peptide",
            json!("MX"),
        ),
        (pairwise.clone(), "/result/target", json!("TCGT")),
        (
            pairwise.clone(),
            "/result/alignment/mode",
            json!("local affine"),
        ),
        (
            pairwise,
            "/result/alignment/scoring",
            json!("match=99 mismatch=-1 gap_open=-5 gap_extend=-1"),
        ),
    ] {
        let mut value: Value = serde_json::from_str(&text).unwrap();
        *value.pointer_mut(pointer).unwrap() = replacement;
        error_kind(
            parse_result_json(&value.to_string()),
            IoErrorKind::InvalidSchema,
        );
    }
}

#[test]
fn documented_output_matches_imported_values_and_errors() {
    let guide = include_str!("../../docs/api/graph-result-input-v1.md");
    let marker = "<!-- graph-result-output -->\n```text\n";
    let expected = guide
        .split_once(marker)
        .unwrap()
        .1
        .split_once("\n```")
        .unwrap()
        .0;
    let graph = parse_gfa::<Dna>(GFA).unwrap();
    let topology = parse_dot::<Dna>(DOT).unwrap();
    let ImportedResult::Alignment(result) = parse_result_json(DISTANCE).unwrap() else {
        panic!("expected a native result");
    };
    let bad_gfa = parse_gfa::<Dna>(&GFA.replace("0M", "1M")).unwrap_err();
    let bad_json =
        parse_result_json(&DISTANCE.replace("\"distance\":1", "\"distance\":0")).unwrap_err();
    let observed = format!(
        "GFA: {} nodes, {} edge, {} classes\npath {}: {}\nDOT: {} nodes, {} edge; paths unavailable\nJSON: {:?}, distance={}\nbad GFA: {}\nbad JSON: {}",
        graph.nodes().len(),
        graph.edges().len(),
        graph.class_count(),
        graph.paths()[0].name(),
        graph.path_sequence(0).unwrap(),
        topology.nodes().len(),
        topology.edges().len(),
        result.kind(),
        result.value()["distance"],
        bad_gfa.code(),
        bad_json.code()
    );
    assert_eq!(observed, expected);
    println!("{observed}");
}
