//! Inert POASTA black-box comparisons against the documented SPOA 4.1.5 CLI.

use align_core::{Alphabet, Dna, Protein, Rna, Score, Sequence, Symbol};
use align_msa::{
    AlignedRow, AlignmentMetrics, MultipleAlignment, SequenceRecord, compare_alignments,
};
use align_poa::{
    LEE_2002_POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_ITERATIVE_V1, PoastaAffineCosts,
    PoastaPoa,
};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/spoa-4.1.5/global-affine-iterative-v1.json");

#[derive(Deserialize)]
struct Fixture {
    profile: String,
    scoring: Scoring,
    cases: Vec<ReferenceCase>,
}

#[derive(Clone, Copy, Deserialize)]
struct Scoring {
    mismatch: Score,
    gap_open_score: Score,
    gap_extension_score: Score,
    poasta_mismatch_cost: Score,
    poasta_gap_opening_cost: Score,
    poasta_gap_extension_cost: Score,
}

#[derive(Deserialize)]
struct ReferenceCase {
    id: String,
    molecule: String,
    records: Vec<ReferenceRecord>,
    reference_alignment: Vec<ReferenceRow>,
    pairwise_alignment_score: Option<Score>,
    reference_graph: ReferenceGraph,
}

#[derive(Deserialize)]
struct ReferenceRecord {
    identifier: String,
    sequence: String,
}

#[derive(Deserialize)]
struct ReferenceRow {
    identifier: String,
    aligned: String,
}

#[derive(Deserialize)]
struct ReferenceGraph {
    segments: usize,
    links: usize,
    paths: usize,
}

fn reference_alignment<A: Alphabet>(case: &ReferenceCase) -> MultipleAlignment<A> {
    let rows = case
        .reference_alignment
        .iter()
        .map(|row| {
            let cells = row
                .aligned
                .bytes()
                .map(|byte| {
                    if byte == b'-' {
                        None
                    } else {
                        Some(Symbol::<A>::from_ascii(byte).expect("fixture symbol is valid"))
                    }
                })
                .collect();
            AlignedRow::new(&row.identifier, cells).expect("fixture row is valid")
        })
        .collect();
    MultipleAlignment::try_new(rows).expect("fixture alignment is valid")
}

fn run_poasta_case<A: Alphabet>(
    case: &ReferenceCase,
    costs: PoastaAffineCosts,
) -> (AlignmentMetrics, bool) {
    let records = case
        .records
        .iter()
        .map(|record| {
            SequenceRecord::new(
                &record.identifier,
                Sequence::<A>::from_ascii(&record.sequence).expect("fixture sequence is valid"),
            )
            .expect("fixture record is valid")
        })
        .collect::<Vec<_>>();
    let result = PoastaPoa::new(costs)
        .align(&records)
        .expect("Rust accelerated POASTA case succeeds");

    result.graph().validate().expect("Rust graph remains valid");
    assert_eq!(result.graph().path_count(), records.len(), "{}", case.id);
    assert_eq!(
        result.graph().path_count(),
        case.reference_graph.paths,
        "{} graph path count",
        case.id
    );
    assert_eq!(
        result.graph().node_count(),
        case.reference_graph.segments,
        "{} graph node/segment count",
        case.id
    );
    assert_eq!(
        result.graph().edge_count(),
        case.reference_graph.links,
        "{} graph edge/link count",
        case.id
    );

    if let Some(reference_score) = case.pairwise_alignment_score {
        assert_eq!(records.len(), 2, "{} score fixture arity", case.id);
        let observed = result.steps()[1]
            .alignment()
            .alignment()
            .expect("second nonempty record is aligned")
            .score();
        assert_eq!(observed, reference_score, "{} affine score", case.id);
    }

    for (row, record) in result.alignment().rows().iter().zip(&records) {
        assert_eq!(
            row.identifier(),
            record.identifier(),
            "{} row identity",
            case.id
        );
        assert_eq!(
            row.to_ungapped_ascii(),
            record.sequence().to_string(),
            "{} source recovery",
            case.id
        );
    }

    let reference = reference_alignment(case);
    let metrics = compare_alignments(result.alignment(), &reference)
        .expect("candidate and reference recover the same sources");
    (metrics, result.alignment() == &reference)
}

#[test]
fn packaged_spoa_matches_poasta_scores_graph_sizes_and_homology() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).expect("SPOA fixture JSON is valid");
    assert_eq!(fixture.cases.len(), 10);
    assert_eq!(fixture.scoring.mismatch, -4);
    assert_eq!(fixture.scoring.gap_open_score, -8);
    assert_eq!(fixture.scoring.gap_extension_score, -2);
    assert_eq!(
        PoastaPoa::new(PoastaAffineCosts::new(4, 6, 2).unwrap()).profile_identifier(),
        LEE_2002_POASTA_2025_GLOBAL_AFFINE_ASTAR_DFS_SUPERBUBBLE_ITERATIVE_V1
    );
    assert_eq!(
        fixture.profile,
        "LEE-2002-POASTA-2025-global-affine-iterative-v1"
    );
    let costs = PoastaAffineCosts::new(
        fixture.scoring.poasta_mismatch_cost,
        fixture.scoring.poasta_gap_opening_cost,
        fixture.scoring.poasta_gap_extension_cost,
    )
    .expect("fixture POASTA costs are valid");

    let mut exact_count = 0_usize;
    for case in &fixture.cases {
        let (metrics, exact) = match case.molecule.as_str() {
            "dna" => run_poasta_case::<Dna>(case, costs),
            "rna" => run_poasta_case::<Rna>(case, costs),
            "protein" => run_poasta_case::<Protein>(case, costs),
            molecule => panic!("unsupported fixture molecule {molecule:?}"),
        };
        println!(
            "accelerated {}: SP={:.9}, TC={:.9}, exact={exact}",
            case.id,
            metrics.sp_score(),
            metrics.tc_score()
        );
        assert_eq!(
            metrics.correct_pairs(),
            metrics.reference_pairs(),
            "{} accelerated SP differs from SPOA 4.1.5",
            case.id
        );
        assert_eq!(
            metrics.correct_columns(),
            metrics.reference_columns(),
            "{} accelerated TC differs from SPOA 4.1.5",
            case.id
        );
        assert!(
            exact,
            "{} accelerated rectangular output differs from SPOA 4.1.5",
            case.id
        );
        exact_count += usize::from(exact);
    }
    assert_eq!(exact_count, fixture.cases.len());
}
