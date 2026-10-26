//! Golden and boundary tests for POASTA graph presentation.

use align::prelude::*;

fn records() -> Vec<SequenceRecord<Dna>> {
    [
        ("reference", "ACGTACGT"),
        ("allele", "ACGTTCGT"),
        ("deletion", "ACGTCGT"),
    ]
    .into_iter()
    .map(|(identifier, sequence)| {
        SequenceRecord::new(identifier, DnaSequence::try_from(sequence).unwrap()).unwrap()
    })
    .collect()
}

#[test]
fn progressive_poasta_report_exposes_graph_and_alignment() {
    let records = records();
    let run = Graph::poasta().run(&records).unwrap();
    let output = format!("{}", run.report());
    assert!(!output.ends_with('\n'));
    assert!(output.contains("Algorithm"));
    assert!(output.contains(": POASTA"));
    assert!(output.contains("Graph nodes"));
    assert!(output.contains("Graph edges"));
    assert!(output.contains("Graph classes"));
    assert!(output.contains("Graph paths"));
    assert!(output.contains("deletion   1  ACGT-CGT  7"));
    assert_eq!(run.result(), &Graph::poasta().align(&records).unwrap());

    let detailed = format!("{}", run.report().verbosity(Verbosity::Detailed));
    assert!(detailed.contains("Construction steps"));
    assert!(detailed.contains("Graph storage bytes"));
    assert!(detailed.contains("Estimated peak bytes"));
}

#[test]
fn query_to_branched_graph_matches_the_golden_contract() {
    let records = records();
    let graph_run = Graph::poasta().run(&records).unwrap();
    assert!(graph_run.result().graph().node_count() > graph_run.result().graph().class_count());
    let query = DnaSequence::try_from("ACGTTACGT").unwrap();
    let run = Graph::poasta()
        .run_query(graph_run.result().graph(), &query)
        .unwrap();
    let actual = format!("{}", run.report());
    assert!(!actual.ends_with('\n'));
    assert_eq!(
        actual,
        include_str!("fixtures/reports/poasta-query-normal.txt").trim_end_matches('\n')
    );
    assert_eq!(
        run.result(),
        &Graph::poasta()
            .align_query(graph_run.result().graph(), &query)
            .unwrap()
    );
}

#[test]
fn graph_query_blocks_keep_node_and_query_coordinates() {
    let records = records();
    let graph = Graph::poasta().align(&records).unwrap();
    let query = DnaSequence::try_from("ACGTTACGT").unwrap();
    let run = Graph::poasta().run_query(graph.graph(), &query).unwrap();
    for width in [1, 3, 60, 100] {
        let output = format!("{}", run.report().block_width(width).unwrap());
        assert!(output.contains("Graph"));
        assert!(output.contains("Query"));
        assert!(output.contains("n0"));
    }
    assert!(run.report().block_width(0).is_err());
    let rust = format!(
        "{}",
        run.report()
            .coordinates(CoordinateStyle::Rust)
            .block_width(4)
            .unwrap()
    );
    assert!(rust.contains("ACGT"), "{rust}");
    assert!(rust.contains("Query"), "{rust}");
    assert!(rust.contains(" 0 "), "{rust}");
    assert!(rust.contains(" 4"), "{rust}");

    let colored = run.report().color(ColorMode::Always);
    assert_eq!(colored.options().color(), ColorMode::Always);
    assert_eq!(format!("{colored}"), format!("{}", run.report()));
}

#[test]
fn detailed_query_report_contains_exact_search_evidence() {
    let records = records();
    let graph = Graph::poasta().align(&records).unwrap();
    let query = DnaSequence::try_from("ACGTTACGT").unwrap();
    let run = Graph::poasta().run_query(graph.graph(), &query).unwrap();
    let detailed = format!("{}", run.report().verbosity(Verbosity::Detailed));
    for field in [
        "Possible states",
        "Discovered states",
        "Expanded states",
        "A* expansions",
        "DFS expansions",
        "Bubble prunes",
        "Modeled peak bytes",
        "Accelerations",
    ] {
        assert!(detailed.contains(field), "missing {field}");
    }
    assert!(!detailed.contains("\u{1b}["));
}
