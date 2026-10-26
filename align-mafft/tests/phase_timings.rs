//! Opt-in FFT-NS-i phase-timing contract tests.

use align_core::{Dna, DnaSequence};
use align_mafft::MafftFftNsi;
use align_msa::SequenceRecord;

fn record(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
    SequenceRecord::new(
        identifier,
        DnaSequence::from_ascii(sequence.as_bytes()).expect("DNA sequence"),
    )
    .expect("sequence record")
}

#[test]
fn profiled_fft_nsi_preserves_the_scientific_result_and_phase_contract() {
    let records = [
        record("a", "ACGTACGT"),
        record("b", "ACGTCGT"),
        record("c", "ACTTACGT"),
        record("d", "ACGTACGTT"),
    ];
    let aligner = MafftFftNsi::<Dna>::new();

    let ordinary = aligner.align(&records).expect("ordinary alignment");
    let profiled = aligner
        .align_profiled(&records)
        .expect("profiled alignment");

    assert_eq!(profiled.result().alignment(), ordinary.alignment());
    assert_eq!(
        profiled.result().initial().alignment(),
        ordinary.initial().alignment()
    );
    assert_eq!(
        profiled.result().accepted_cycle_count(),
        ordinary.accepted_cycle_count()
    );
    assert_eq!(profiled.result().termination(), ordinary.termination());
    assert_eq!(
        profiled.result().final_distance_matrix(),
        ordinary.final_distance_matrix()
    );
    assert_eq!(
        profiled.result().final_pair_weights(),
        ordinary.final_pair_weights()
    );
    assert_eq!(
        profiled.result().final_score().to_bits(),
        ordinary.final_score().to_bits()
    );
    assert_eq!(
        profiled.result().memory_estimate(),
        ordinary.memory_estimate()
    );
    let timings = profiled.timings();
    assert_eq!(timings.phases().len(), 7);
    for name in [
        "initial_alignment",
        "refinement_setup",
        "cycle_setup",
        "candidate_evaluation",
        "cycle_commit",
        "final_scoring",
        "result_accounting",
    ] {
        assert!(timings.phase(name).is_some(), "missing phase {name}");
    }
    assert!(timings.total() >= timings.attributed());
}

#[test]
fn detailed_candidate_timings_preserve_results_and_account_for_every_candidate() {
    let records = [
        record("a", "ACGTACGT"),
        record("b", "ACGTCGT"),
        record("c", "ACTTACGT"),
        record("d", "ACGTACGTT"),
    ];
    let aligner = MafftFftNsi::<Dna>::new();

    let ordinary = aligner.align(&records).expect("ordinary alignment");
    let detailed = aligner
        .align_detailed_profiled(&records)
        .expect("detailed profiled alignment");

    assert_eq!(detailed.result().alignment(), ordinary.alignment());
    assert_eq!(
        detailed.result().accepted_cycle_count(),
        ordinary.accepted_cycle_count()
    );
    assert_eq!(detailed.result().termination(), ordinary.termination());
    assert_eq!(
        detailed.result().final_score().to_bits(),
        ordinary.final_score().to_bits()
    );
    assert_eq!(
        detailed.result().memory_estimate(),
        ordinary.memory_estimate()
    );

    let candidate_timings = detailed.candidate_timings();
    let expected_candidates = detailed
        .result()
        .cycles()
        .iter()
        .map(align_mafft::FftNsiCycle::partition_count)
        .sum::<usize>();
    assert_eq!(candidate_timings.candidate_count(), expected_candidates);
    assert!(candidate_timings.total() >= candidate_timings.attributed());
    assert_eq!(
        candidate_timings.phases().map(|(name, _)| name),
        [
            "projection",
            "profile_preparation",
            "fft_anchor_search",
            "profile_dp",
            "profile_assembly",
            "alignment_materialization",
            "validation",
            "wsp_scoring",
            "accounting",
        ]
    );
    assert!(candidate_timings.profile_dp_block_count() > 0);
    assert!(candidate_timings.profile_dp_cell_count() > 0);
    assert_eq!(
        candidate_timings
            .profile_dp_homology_term_cell_counts()
            .into_iter()
            .sum::<usize>(),
        candidate_timings.profile_dp_cell_count()
    );
    assert!(candidate_timings.profile_dp_homology_term_evaluations() > 0);
    assert!(candidate_timings.profile_dp_total() >= candidate_timings.profile_dp_attributed());
    assert_eq!(
        candidate_timings.profile_dp_phases().map(|(name, _)| name),
        [
            "workspace_initialization",
            "recurrence_fill",
            "endpoint_traceback",
            "block_accounting",
            "score_replay",
        ]
    );
    assert!(detailed.timings().phase("candidate_evaluation").is_some());
    assert!(detailed.timings().total() >= detailed.timings().attributed());
}
