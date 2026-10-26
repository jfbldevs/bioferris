//! Opt-in Clustal Omega phase-timing contract tests.

use align_clustal::{ClustalOmega, ClustalOmegaParameters, OmegaEngineMode, OmegaGuideTreeMode};
use align_core::{Dna, DnaSequence};
use align_msa::SequenceRecord;

fn record(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
    SequenceRecord::new(
        identifier,
        DnaSequence::from_ascii(sequence.as_bytes()).expect("DNA sequence"),
    )
    .expect("sequence record")
}

#[test]
fn profiled_clustal_preserves_the_scientific_result_and_phase_contract() {
    let records = [
        record("a", "ACGTACGT"),
        record("b", "ACGTCGT"),
        record("c", "ACTTACGT"),
        record("d", "ACGTACGTT"),
    ];
    let parameters = ClustalOmegaParameters::nucleotide_v1_2_4()
        .with_engine_mode(OmegaEngineMode::Viterbi)
        .with_guide_tree_mode(OmegaGuideTreeMode::Full);
    let aligner = ClustalOmega::<Dna>::with_parameters(parameters).expect("valid parameters");

    let ordinary = aligner.align(&records).expect("ordinary alignment");
    let profiled = aligner
        .align_profiled(&records)
        .expect("profiled alignment");

    assert_eq!(profiled.result().alignment(), ordinary.alignment());
    assert_eq!(profiled.result().guide_tree(), ordinary.guide_tree());
    assert_eq!(
        profiled.result().progressive_steps(),
        ordinary.progressive_steps()
    );
    assert_eq!(
        profiled.result().memory_estimate(),
        ordinary.memory_estimate()
    );
    assert_eq!(
        profiled.result().profile_identifier(),
        ordinary.profile_identifier()
    );
    let timings = profiled.timings();
    assert_eq!(timings.phases().len(), 4);
    for name in [
        "input_validation",
        "guide_tree",
        "progressive_alignment",
        "finalization",
    ] {
        assert!(timings.phase(name).is_some(), "missing phase {name}");
    }
    assert!(timings.total() >= timings.attributed());
}
