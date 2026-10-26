//! Opt-in iterative POASTA phase-timing contract tests.

use align_core::{Dna, DnaSequence};
use align_msa::SequenceRecord;
use align_poa::{PoastaAffineCosts, PoastaPoa};

fn record(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
    SequenceRecord::new(
        identifier,
        DnaSequence::from_ascii(sequence.as_bytes()).expect("DNA sequence"),
    )
    .expect("sequence record")
}

#[test]
fn profiled_iterative_poasta_preserves_the_exact_result_and_phase_contract() {
    let records = [
        record("a", "ACGTACGT"),
        record("b", "ACGTCGT"),
        record("c", "ACTTACGT"),
        record("d", "ACGTACGTT"),
    ];
    let costs = PoastaAffineCosts::new(4, 6, 2).expect("valid costs");
    let aligner = PoastaPoa::optimized(costs);

    let ordinary = aligner.align(&records).expect("ordinary alignment");
    let profiled = aligner
        .align_profiled(&records)
        .expect("profiled alignment");

    assert_eq!(profiled.result(), &ordinary);
    let timings = profiled.timings();
    for name in [
        "input_validation",
        "integration_targets",
        "graph_integration",
        "graph_index",
        "search_setup",
        "search",
        "traceback",
        "result_accounting",
        "verification",
        "alignment_serialization",
    ] {
        assert!(timings.phase(name).is_some(), "missing phase {name}");
    }
    assert!(timings.total() >= timings.attributed());

    let detailed = aligner
        .align_detailed_profiled(&records)
        .expect("detailed profiled alignment");
    assert_eq!(detailed.result(), &ordinary);
    assert!(detailed.search_timings().phase("frontier_search").is_some());
    let search = detailed
        .timings()
        .phase("search")
        .expect("coarse search phase");
    assert!(detailed.search_timings().attributed() <= search.elapsed());
    for name in [
        "exact_affine_band_layout",
        "exact_affine_table_prepare",
        "exact_affine_fill",
    ] {
        assert!(
            detailed.search_timings().phase(name).is_some(),
            "missing detailed phase {name}"
        );
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn detailed_exact_affine_row_counters_partition_the_filled_tables() {
    let first = "A".repeat(64);
    let second = "C".repeat(64);
    let records = [record("a", &first), record("b", &second)];
    let costs = PoastaAffineCosts::new(4, 6, 2).expect("valid costs");
    let detailed = PoastaPoa::optimized(costs)
        .align_detailed_profiled(&records)
        .expect("detailed hard alignment");
    let counters = detailed.search_work_counters();

    let reported_cells = detailed
        .result()
        .steps()
        .iter()
        .filter_map(|step| step.alignment().alignment())
        .map(|alignment| alignment.statistics().exact_affine_certificate_cells())
        .sum::<usize>();
    let reported_relaxed_cells = detailed
        .result()
        .steps()
        .iter()
        .filter_map(|step| step.alignment().alignment())
        .map(|alignment| alignment.statistics().linear_heuristic_cells())
        .sum::<usize>();
    let reported_exact_affine_certificates = detailed
        .result()
        .steps()
        .iter()
        .filter_map(|step| step.alignment().alignment())
        .filter(|alignment| alignment.statistics().used_exact_affine_certificate())
        .count();
    let partitioned_cells = counters
        .exact_affine_source_cells()
        .checked_add(counters.exact_affine_single_predecessor_cells())
        .and_then(|cells| cells.checked_add(counters.exact_affine_two_predecessor_cells()))
        .and_then(|cells| cells.checked_add(counters.exact_affine_many_predecessor_cells()))
        .expect("small counter partition");

    assert!(counters.exact_affine_attempts() > 0);
    assert!(counters.exact_affine_filled_attempts() > 0);
    assert!(counters.exact_affine_single_predecessor_rows() > 0);
    assert!(counters.relaxed_linear_rows() > 0);
    assert_eq!(
        counters.relaxed_linear_table_attempts(),
        counters
            .relaxed_linear_table_diagonal_certificates()
            .checked_add(counters.relaxed_linear_table_affine_opening_certificates())
            .and_then(|outcomes| outcomes.checked_add(counters.relaxed_linear_table_misses()))
            .expect("small relaxed-table outcome partition")
    );
    assert_eq!(
        counters
            .relaxed_linear_table_certificate_cells()
            .checked_add(counters.relaxed_linear_table_miss_cells())
            .expect("small relaxed-table cell partition"),
        reported_relaxed_cells
    );
    assert!(counters.relaxed_linear_insertion_cells() > 0);
    assert!(
        counters.relaxed_linear_insertion_avx2_cells() <= counters.relaxed_linear_insertion_cells()
    );
    assert_eq!(
        counters
            .relaxed_linear_zero_node_successor_rows()
            .checked_add(counters.relaxed_linear_single_node_successor_rows())
            .and_then(|rows| {
                rows.checked_add(counters.relaxed_linear_multiple_node_successor_rows())
            })
            .expect("small relaxed-row partition"),
        counters.relaxed_linear_rows()
    );
    assert_eq!(
        counters
            .relaxed_linear_two_node_successor_rows()
            .checked_add(counters.relaxed_linear_three_or_four_node_successor_rows())
            .and_then(|rows| {
                rows.checked_add(counters.relaxed_linear_five_plus_node_successor_rows())
            })
            .expect("small multiple-successor row partition"),
        counters.relaxed_linear_multiple_node_successor_rows()
    );
    assert!(
        counters.relaxed_linear_multiple_node_distinct_symbol_cell_visits()
            <= counters.relaxed_linear_multiple_node_successor_cell_visits()
    );
    assert!(
        counters.relaxed_linear_equivalent_multiple_transition_rows()
            <= counters.relaxed_linear_multiple_node_successor_rows()
    );
    assert!(
        counters.relaxed_linear_equivalent_multiple_transition_cells()
            <= counters.relaxed_linear_multiple_node_successor_cells()
    );
    assert!(
        counters.relaxed_linear_equivalent_multiple_transition_insertion_cells()
            <= counters.relaxed_linear_insertion_cells()
    );
    assert!(
        counters.relaxed_linear_equivalent_multiple_transition_successor_cell_visits()
            <= counters.relaxed_linear_multiple_node_successor_cell_visits()
    );
    assert!(counters.relaxed_linear_single_node_successor_rows() > 0);
    assert!(counters.relaxed_linear_successor_cell_visits() > 0);
    assert!(
        counters.relaxed_linear_equivalent_single_transition_rows()
            <= counters.relaxed_linear_rows()
    );
    assert!(
        counters.relaxed_linear_equivalent_single_transition_cells()
            <= counters.relaxed_linear_insertion_cells()
                + counters.relaxed_linear_equivalent_single_transition_rows()
    );
    assert!(
        counters.relaxed_linear_equivalent_single_transition_insertion_cells()
            <= counters.relaxed_linear_insertion_cells()
    );
    assert!(
        counters.relaxed_linear_equivalent_single_transition_successor_cell_visits()
            <= counters.relaxed_linear_successor_cell_visits()
    );
    assert!(counters.relaxed_linear_initialized_rows() > 0);
    assert!(counters.relaxed_linear_initialized_cells() > 0);
    assert!(
        counters
            .relaxed_linear_initialized_rows()
            .checked_add(counters.relaxed_linear_reused_rows())
            .is_some_and(|rows| rows >= counters.relaxed_linear_rows())
    );
    assert_eq!(partitioned_cells, reported_cells);
    assert_eq!(
        counters
            .exact_affine_bypass_attempts()
            .checked_add(counters.exact_affine_after_relaxed_attempts())
            .expect("small exact-affine route partition"),
        counters.exact_affine_attempts()
    );
    assert_eq!(
        counters
            .exact_affine_bypass_cells()
            .checked_add(counters.exact_affine_after_relaxed_cells())
            .expect("small exact-affine route cell partition"),
        partitioned_cells
    );
    assert_eq!(
        counters
            .exact_affine_bypass_certificates()
            .checked_add(counters.exact_affine_after_relaxed_certificates())
            .expect("small exact-affine certificate route partition"),
        reported_exact_affine_certificates
    );
    assert!(
        counters.exact_affine_single_predecessor_fast_cells()
            <= counters.exact_affine_single_predecessor_cells()
    );
    assert!(
        counters.exact_affine_single_predecessor_avx2_cells()
            <= counters.exact_affine_single_predecessor_fast_cells()
    );
    assert!(
        counters.exact_affine_two_predecessor_fast_cells()
            <= counters.exact_affine_two_predecessor_cells()
    );
    assert!(
        counters.exact_affine_two_predecessor_avx2_cells()
            <= counters.exact_affine_two_predecessor_fast_cells()
    );
    assert!(
        counters.exact_affine_many_predecessor_avx2_cells()
            <= counters.exact_affine_many_predecessor_cells()
    );
}
