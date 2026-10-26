//! Exact equivalence tests for the pairwise facade.

use align::PairwiseMetric;
use align::prelude::*;
use align_core::{AffineGap, LinearGap, MatchMismatch};
use align_fast::{WavefrontAffine, WavefrontAffinePenalties, XDrop};
use align_pairwise::{
    AffineEndGapAligner, AffineLocalAligner, Gotoh, NeedlemanWunsch, SmithWaterman,
};

fn sequences() -> (DnaSequence, DnaSequence) {
    (
        DnaSequence::try_from("GATTACA").unwrap(),
        DnaSequence::try_from("GCATGCT").unwrap(),
    )
}

#[test]
fn linear_facades_are_exact_kernel_adapters() {
    let (first, second) = sequences();
    let scoring = MatchMismatch::new(2, -1);
    let gap = LinearGap::new(-2).unwrap();

    let facade_global = Pairwise::needleman_wunsch()
        .scores(2, -1)
        .gap(-2)
        .align(&first, &second)
        .unwrap();
    let direct_global = NeedlemanWunsch::new(scoring, gap)
        .align(&first, &second)
        .unwrap();
    assert_eq!(facade_global, direct_global);

    let facade_local = Pairwise::smith_waterman()
        .scores(2, -1)
        .gap(-2)
        .align(&first, &second)
        .unwrap();
    let direct_local = SmithWaterman::new(scoring, gap)
        .align(&first, &second)
        .unwrap();
    assert_eq!(facade_local, direct_local);

    // `score()` uses the score-only kernels, which return the alignment's score.
    let global_score = Pairwise::needleman_wunsch()
        .scores(2, -1)
        .gap(-2)
        .score(&first, &second)
        .unwrap();
    assert_eq!(global_score, direct_global.score());
    assert_eq!(
        global_score,
        NeedlemanWunsch::new(scoring, gap)
            .score(&first, &second)
            .unwrap()
    );
    let local_score = Pairwise::smith_waterman()
        .scores(2, -1)
        .gap(-2)
        .score(&first, &second)
        .unwrap();
    assert_eq!(local_score, direct_local.score());
    assert_eq!(
        local_score,
        SmithWaterman::new(scoring, gap)
            .score(&first, &second)
            .unwrap()
    );
}

#[test]
fn every_gotoh_boundary_mode_delegates_exactly() {
    let (first, second) = sequences();
    let scoring = MatchMismatch::new(3, -4);
    let gap = AffineGap::new(-5, -1).unwrap();

    let facade_global = Pairwise::gotoh()
        .scores(3, -4)
        .gaps(-5, -1)
        .global()
        .align(&first, &second)
        .unwrap();
    let direct_global = Gotoh::new(scoring, gap).align(&first, &second).unwrap();
    assert_eq!(facade_global, direct_global);
    assert_eq!(
        Pairwise::gotoh()
            .scores(3, -4)
            .gaps(-5, -1)
            .score(&first, &second)
            .unwrap(),
        direct_global.score()
    );

    let facade_local = Pairwise::gotoh()
        .scores(3, -4)
        .gaps(-5, -1)
        .local()
        .align(&first, &second)
        .unwrap();
    let direct_local = AffineLocalAligner::new(scoring, gap)
        .align(&first, &second)
        .unwrap();
    assert_eq!(facade_local, direct_local);

    let configured = Pairwise::gotoh().scores(3, -4).gaps(-5, -1);
    for (facade_aligner, free_ends) in [
        (configured.semiglobal(), FreeEndGaps::semiglobal()),
        (
            configured.fit_first_in_second(),
            FreeEndGaps::fit_first_in_second(),
        ),
        (
            configured.fit_second_in_first(),
            FreeEndGaps::fit_second_in_first(),
        ),
        (
            configured.first_suffix_second_prefix(),
            FreeEndGaps::first_suffix_second_prefix(),
        ),
        (
            configured.second_suffix_first_prefix(),
            FreeEndGaps::second_suffix_first_prefix(),
        ),
        (
            configured.free_ends(FreeEndGaps::semiglobal()),
            FreeEndGaps::semiglobal(),
        ),
    ] {
        let facade = facade_aligner.align(&first, &second).unwrap();
        let direct = AffineEndGapAligner::new(scoring, gap, free_ends)
            .align(&first, &second)
            .unwrap();
        assert_eq!(facade, direct);
    }
}

#[test]
fn accelerated_facades_preserve_native_results() {
    let (first, second) = sequences();

    let myers = Pairwise::myers();
    let direct_pattern = align_fast::MyersPattern::new(&first).unwrap();
    assert_eq!(
        myers.distance(&first, &second).unwrap(),
        direct_pattern.global_distance(&second).unwrap()
    );
    assert_eq!(
        myers.align(&first, &second).unwrap(),
        direct_pattern.align_global(&second).unwrap()
    );
    assert_eq!(
        myers.find_all(&first, &second, 3).unwrap(),
        direct_pattern.find_all(&second, 3).unwrap()
    );

    let facade_wfa = Pairwise::wfa().penalties(4, 6, 2);
    let direct_wfa = WavefrontAffine::new(WavefrontAffinePenalties::new(4, 6, 2).unwrap());
    assert_eq!(
        facade_wfa.penalty(&first, &second).unwrap(),
        direct_wfa.penalty(&first, &second).unwrap()
    );
    assert_eq!(
        facade_wfa.align(&first, &second).unwrap(),
        direct_wfa.align(&first, &second).unwrap()
    );

    let facade_xdrop = Pairwise::xdrop().scores(2, -3).gaps(-6, -1).threshold(2);
    let direct_xdrop = XDrop::new(
        MatchMismatch::new(2, -3),
        AffineGap::new(-6, -1).unwrap(),
        2,
    )
    .unwrap();
    assert_eq!(
        facade_xdrop.score(&first, &second).unwrap(),
        direct_xdrop.score(&first, &second).unwrap()
    );
    assert_eq!(
        facade_xdrop.align(&first, &second).unwrap(),
        direct_xdrop.align(&first, &second).unwrap()
    );
    let seed = XDropSeed::new(0, 0);
    assert_eq!(
        facade_xdrop
            .score_from_seed(&first, &second, seed, XDropDirection::Both)
            .unwrap(),
        direct_xdrop
            .score_from_seed(&first, &second, seed, XDropDirection::Both)
            .unwrap()
    );
    assert_eq!(
        facade_xdrop
            .align_from_seed(&first, &second, seed, XDropDirection::Both)
            .unwrap(),
        direct_xdrop
            .align_from_seed(&first, &second, seed, XDropDirection::Both)
            .unwrap()
    );
}

#[test]
fn adaptive_affine_facade_is_exact_and_auditable_across_all_backends() {
    let adaptive = Pairwise::adaptive_affine().penalties(4, 6, 2);
    let identity = DnaSequence::from_ascii(vec![b'A'; 256]).unwrap();
    let short_first = DnaSequence::try_from("ACGTACGT").unwrap();
    let short_second = DnaSequence::try_from("ACGTTTGT").unwrap();
    let dense_first = DnaSequence::from_ascii(vec![b'A'; 256]).unwrap();
    let dense_second = DnaSequence::from_ascii(vec![b'C'; 256]).unwrap();

    let mut workspace = AdaptiveAffineWorkspace::new();
    for (first, second, expected) in [
        (&identity, &identity, AdaptiveAffineBackend::Identity),
        (
            &short_first,
            &short_second,
            AdaptiveAffineBackend::Wavefront,
        ),
        (&dense_first, &dense_second, AdaptiveAffineBackend::Gotoh),
    ] {
        let score = adaptive
            .penalty_with_workspace(first, second, &mut workspace)
            .unwrap();
        let alignment = adaptive
            .align_with_workspace(first, second, &mut workspace)
            .unwrap();

        assert_eq!(score.backend(), expected);
        assert_eq!(alignment.backend(), expected);
        assert_eq!(score.penalty(), alignment.penalty());
        assert_eq!(score.alignment_score(), alignment.alignment().score());
        assert_eq!(alignment.penalty(), -alignment.alignment().score());
        assert_eq!(
            alignment.decision(),
            adaptive.decision(first, second).unwrap()
        );
    }

    workspace.release();
    let replay = adaptive
        .align_with_workspace(&dense_first, &dense_second, &mut workspace)
        .unwrap();
    assert_eq!(replay.backend(), AdaptiveAffineBackend::Gotoh);
    assert_eq!(replay.penalty(), 1_024);
}

#[test]
fn adaptive_affine_reporting_and_batching_preserve_backend_provenance() {
    let adaptive = Pairwise::adaptive_affine().penalties(4, 6, 2);
    let identity = DnaSequence::from_ascii(vec![b'A'; 256]).unwrap();
    let dense = DnaSequence::from_ascii(vec![b'C'; 256]).unwrap();

    let run = adaptive.run(&identity, &dense).unwrap();
    assert_eq!(run.result().backend(), AdaptiveAffineBackend::Gotoh);
    assert_eq!(run.info().metric(), PairwiseMetric::Penalty(1_024));
    let report = format!("{}", run.report().verbosity(Verbosity::Detailed));
    assert!(report.contains("Selected backend"));
    assert!(report.contains("simd-gotoh"));
    assert!(report.contains("Dense cells"));
    assert!(report.contains("65536"));

    let penalty_run = adaptive.run_penalty(&identity, &identity).unwrap();
    assert_eq!(
        penalty_run.result().backend(),
        AdaptiveAffineBackend::Identity
    );
    let penalty_report = format!("{}", penalty_run.report().verbosity(Verbosity::Detailed));
    assert!(penalty_report.contains("Selected backend"));
    assert!(penalty_report.contains("identity"));

    let pairs = [(&identity, &identity), (&identity, &dense)];
    let results = adaptive.threads(2).align_many(pairs).unwrap();
    assert_eq!(results[0].backend(), AdaptiveAffineBackend::Identity);
    assert_eq!(results[1].backend(), AdaptiveAffineBackend::Gotoh);
}

#[test]
fn common_result_helpers_are_deterministic() {
    let alignment = PairwiseAlignment::try_new(
        3,
        0..4,
        0..5,
        vec![
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Substitution,
            AlignmentOperation::Insertion,
            AlignmentOperation::Insertion,
            AlignmentOperation::Deletion,
        ],
    )
    .unwrap();
    assert_eq!(alignment.cigar(), "2=1X2I1D");
    assert_eq!(alignment.matches(), 2);
    assert_eq!(alignment.identity(), Some(2.0 / 6.0));
}

#[test]
fn configuration_errors_remain_structured() {
    let (first, second) = sequences();
    assert!(matches!(
        Pairwise::needleman_wunsch().gap(1).align(&first, &second),
        Err(Error::Scoring(_))
    ));
    assert!(matches!(
        Pairwise::wfa().penalties(0, 6, 2).align(&first, &second),
        Err(Error::Wavefront(_))
    ));
    assert!(matches!(
        Pairwise::xdrop().threshold(-1).align(&first, &second),
        Err(Error::XDrop(_))
    ));
}

#[test]
fn reusable_aligners_accept_batches() {
    let (first, second) = sequences();
    let pairs = [(&first, &second), (&second, &first)];
    let results = Pairwise::gotoh().threads(2).align_many(pairs).unwrap();
    assert_eq!(results.len(), 2);
    assert_eq!(
        results[0],
        Pairwise::gotoh().align(&first, &second).unwrap()
    );

    let linear = Pairwise::needleman_wunsch()
        .threads(2)
        .align_many(pairs)
        .unwrap();
    assert_eq!(
        linear[1],
        Pairwise::needleman_wunsch().align(&second, &first).unwrap()
    );

    let local = Pairwise::smith_waterman()
        .threads(2)
        .align_many(pairs)
        .unwrap();
    assert_eq!(
        local[0],
        Pairwise::smith_waterman().align(&first, &second).unwrap()
    );

    let myers = Pairwise::myers().threads(2).align_many(pairs).unwrap();
    assert_eq!(myers[0], Pairwise::myers().align(&first, &second).unwrap());

    let wavefront = Pairwise::wfa().threads(2).align_many(pairs).unwrap();
    assert_eq!(
        wavefront[1],
        Pairwise::wfa().align(&second, &first).unwrap()
    );

    let xdrop = Pairwise::xdrop().threads(2).align_many(pairs).unwrap();
    assert_eq!(xdrop[0], Pairwise::xdrop().align(&first, &second).unwrap());

    assert!(matches!(
        Pairwise::gotoh().threads(0).align_many(pairs),
        Err(Error::InvalidThreadCount { value: 0 })
    ));
}
