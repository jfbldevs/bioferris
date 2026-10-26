//! Golden and boundary tests for human-readable pairwise output.

use align::prelude::*;
use proptest::prelude::*;

fn sequences() -> (DnaSequence, DnaSequence) {
    (
        DnaSequence::try_from("GATTACA").unwrap(),
        DnaSequence::try_from("GCATGCT").unwrap(),
    )
}

#[test]
fn gotoh_normal_report_matches_the_golden_contract() {
    let (first, second) = sequences();
    let run = Pairwise::gotoh()
        .scores(2, -1)
        .gaps(-5, -1)
        .run(&first, &second)
        .unwrap();
    let actual = format!("{}", run.report().labels("Target", "Query"));
    assert_eq!(
        actual,
        include_str!("fixtures/reports/gotoh-normal.txt").trim_end_matches('\n')
    );
    assert_eq!(
        run.result(),
        &Pairwise::gotoh().align(&first, &second).unwrap()
    );
}

#[test]
fn all_seven_pairwise_families_have_exact_unambiguous_compact_output() {
    let (first, second) = sequences();
    let reports = [
        format!(
            "{}",
            Pairwise::needleman_wunsch()
                .run(&first, &second)
                .unwrap()
                .report()
                .verbosity(Verbosity::Compact)
        ),
        format!(
            "{}",
            Pairwise::smith_waterman()
                .run(&first, &second)
                .unwrap()
                .report()
                .verbosity(Verbosity::Compact)
        ),
        format!(
            "{}",
            Pairwise::gotoh()
                .run(&first, &second)
                .unwrap()
                .report()
                .verbosity(Verbosity::Compact)
        ),
        format!(
            "{}",
            Pairwise::myers()
                .run(&first, &second)
                .unwrap()
                .report()
                .verbosity(Verbosity::Compact)
        ),
        format!(
            "{}",
            Pairwise::wfa()
                .run(&first, &second)
                .unwrap()
                .report()
                .verbosity(Verbosity::Compact)
        ),
        format!(
            "{}",
            Pairwise::adaptive_affine()
                .run(&first, &second)
                .unwrap()
                .report()
                .verbosity(Verbosity::Compact)
        ),
        format!(
            "{}",
            Pairwise::xdrop()
                .run(&first, &second)
                .unwrap()
                .report()
                .verbosity(Verbosity::Compact)
        ),
    ];
    let expected = [
        "Algorithm     : Needleman-Wunsch\nMode          : global linear\nScoring       : match=2 mismatch=-1 gap=-2\nSequences     : First=7 Second=7\nRange         : First 1..7, Second 1..7\nScore         : 2",
        "Algorithm     : Smith-Waterman\nMode          : local linear\nScoring       : match=2 mismatch=-1 gap=-2\nSequences     : First=7 Second=7\nRange         : First 2..3, Second 3..4\nScore         : 4",
        "Algorithm     : Gotoh\nMode          : global affine\nScoring       : match=2 mismatch=-1 gap_open=-5 gap_extend=-1\nSequences     : First=7 Second=7\nRange         : First 1..7, Second 1..7\nScore         : 2",
        "Algorithm     : Myers bit-vector\nMode          : global unit edit traceback\nScoring       : match=0 mismatch=1 insertion=1 deletion=1\nSequences     : First=7 Second=7\nRange         : First 1..7, Second 1..7\nDistance      : 4",
        "Algorithm     : WFA\nMode          : global affine wavefront\nScoring       : mismatch=4 gap_open=6 gap_extend=2 (positive penalties)\nSequences     : First=7 Second=7\nRange         : First 1..7, Second 1..7\nPenalty       : 16",
        "Algorithm     : Adaptive exact affine\nMode          : global identity/WFA/SIMD-Gotoh crossover\nScoring       : mismatch=4 gap_open=6 gap_extend=2 (positive WFA penalties)\nSequences     : First=7 Second=7\nRange         : First 1..7, Second 1..7\nPenalty       : 16",
        "Algorithm     : X-drop\nMode          : forward from origin\nScoring       : match=2 mismatch=-1 gap_open=-5 gap_extend=-1\nSequences     : First=7 Second=7\nRange         : First 1..6, Second 1..6\nScore         : 3",
    ];
    assert_eq!(reports, expected);
    assert!(reports[3].contains("Distance"));
    assert!(reports[4].contains("Penalty"));
    assert!(reports[5].contains("Penalty"));
    assert!(reports.iter().all(|report| !report.contains("\u{1b}[")));
    assert!(reports.iter().all(|report| !report.ends_with('\n')));
}

#[test]
fn every_gotoh_mode_is_named_and_custom_free_ends_are_explained() {
    let (first, second) = sequences();
    let cases = [
        (Pairwise::gotoh().global(), "global affine"),
        (Pairwise::gotoh().local(), "local affine"),
        (Pairwise::gotoh().semiglobal(), "semiglobal affine"),
        (
            Pairwise::gotoh().fit_first_in_second(),
            "fit first in second",
        ),
        (
            Pairwise::gotoh().fit_second_in_first(),
            "fit second in first",
        ),
        (
            Pairwise::gotoh().first_suffix_second_prefix(),
            "first suffix to second prefix",
        ),
        (
            Pairwise::gotoh().second_suffix_first_prefix(),
            "second suffix to first prefix",
        ),
        (
            Pairwise::gotoh().free_ends(FreeEndGaps::new(true, false, true, true)),
            "custom free-end affine",
        ),
    ];
    for (aligner, name) in cases {
        let run = aligner.run(&first, &second).unwrap();
        assert_eq!(run.info().mode(), name);
        let score_run = aligner.run_score(&first, &second).unwrap();
        assert_eq!(score_run.info().mode(), format!("{name} score-only"));
        let detailed = format!("{}", run.report().verbosity(Verbosity::Detailed));
        if matches!(aligner.mode(), GotohMode::EndGaps(_)) {
            assert!(detailed.contains("Free ends"));
        }
    }
}

#[test]
fn block_coordinates_cover_widths_gaps_and_rust_ranges() {
    let first = DnaSequence::try_from("ACGTAC").unwrap();
    let second = DnaSequence::try_from("ACGAC").unwrap();
    let run = Pairwise::needleman_wunsch().run(&first, &second).unwrap();
    for width in [1, 2, 5, 60, 100] {
        let report = run.report().block_width(width).unwrap();
        let output = format!("{report}");
        assert!(output.contains("First"));
        assert!(output.contains("Second"));
        assert!(output.contains("1D"));
    }
    assert!(matches!(
        run.report().block_width(0),
        Err(align::ReportError::InvalidBlockWidth { value: 0 })
    ));
    let rust = format!(
        "{}",
        run.report()
            .coordinates(CoordinateStyle::Rust)
            .verbosity(Verbosity::Compact)
    );
    assert!(rust.contains("First 0..6, Second 0..5"));
}

#[test]
fn empty_and_no_match_alignments_render_without_invented_identity() {
    let empty = DnaSequence::try_from("").unwrap();
    let empty_run = Pairwise::needleman_wunsch().run(&empty, &empty).unwrap();
    let output = format!("{}", empty_run.report());
    assert!(output.contains("Identity      : 0/0 (-)"));
    assert!(output.contains("CIGAR         : -"));
    assert!(output.contains("(empty alignment)"));

    let first = DnaSequence::try_from("AAAA").unwrap();
    let second = DnaSequence::try_from("CCCC").unwrap();
    let no_match = Pairwise::needleman_wunsch()
        .scores(1, -1)
        .gap(-2)
        .run(&first, &second)
        .unwrap();
    assert!(format!("{}", no_match.report()).contains("Identity      : 0/4 (0.00%)"));
}

#[test]
fn myers_distance_search_and_xdrop_evidence_have_dedicated_reports() {
    let (pattern, text) = sequences();
    let distance = Pairwise::myers().run_distance(&pattern, &text).unwrap();
    assert_eq!(distance.distance(), 4);
    assert_eq!(
        format!("{}", distance.report()),
        "Algorithm      : Myers bit-vector\nMode           : global unit edit distance\nPattern length : 7\nText length    : 7\nDistance       : 4"
    );
    assert_eq!(
        distance.report().color(ColorMode::Auto).options().color(),
        ColorMode::Auto
    );

    let short_text = DnaSequence::try_from("TTT").unwrap();
    let none = Pairwise::myers()
        .run_find_all(&pattern, &short_text, 0)
        .unwrap();
    assert!(format!("{}", none.report()).contains("Matches          : 0"));
    assert_eq!(
        none.report().color(ColorMode::Always).options().color(),
        ColorMode::Always
    );
    let short_pattern = DnaSequence::try_from("A").unwrap();
    let repeated_text = DnaSequence::try_from("AAAA").unwrap();
    let abundant = Pairwise::myers()
        .run_find_all(&short_pattern, &repeated_text, 0)
        .unwrap();
    assert_eq!(abundant.matches().len(), 4);
    assert!(format!("{}", abundant.report()).contains("4              4  0"));
    assert!(!format!("{}", distance.report()).ends_with('\n'));
    assert!(!format!("{}", none.report()).ends_with('\n'));

    let seeded = Pairwise::xdrop()
        .run_from_seed(&pattern, &text, XDropSeed::new(0, 0), XDropDirection::Both)
        .unwrap();
    let detailed = format!("{}", seeded.report().verbosity(Verbosity::Detailed));
    assert!(detailed.contains("Seed"));
    assert!(detailed.contains("Upstream termination"));
    assert!(detailed.contains("Downstream termination"));
}

#[test]
fn pairwise_human_output_is_independent_of_worker_configuration() {
    let (first, second) = sequences();
    let serial = Pairwise::gotoh().threads(1).run(&first, &second).unwrap();
    let parallel = Pairwise::gotoh().threads(8).run(&first, &second).unwrap();
    assert_eq!(
        format!("{}", serial.report().verbosity(Verbosity::Detailed)),
        format!("{}", parallel.report().verbosity(Verbosity::Detailed))
    );

    let colored = serial.report().color(ColorMode::Always);
    assert_eq!(colored.options().color(), ColorMode::Always);
    assert_eq!(format!("{colored}"), format!("{}", serial.report()));

    let escaped = format!(
        "{}",
        serial
            .report()
            .labels("\x1b[31mTarget\n", "Query\r")
            .verbosity(Verbosity::Compact)
    );
    assert!(!escaped.contains("\x1b["));
    assert!(escaped.contains(r"\u{1b}[31mTarget\n"));
    assert!(escaped.contains(r"Query\r"));
}

#[test]
fn score_only_paths_report_the_native_objective_without_faking_a_traceback() {
    let (first, second) = sequences();

    let linear = Pairwise::needleman_wunsch()
        .run_score(&first, &second)
        .unwrap();
    assert_eq!(
        *linear.result(),
        Pairwise::needleman_wunsch().score(&first, &second).unwrap()
    );
    assert_eq!(
        format!("{}", linear.report()),
        "Algorithm     : Needleman-Wunsch\nMode          : global linear score-only\nScoring       : match=2 mismatch=-1 gap=-2\nSequences     : First=7 Second=7\nScore         : 2"
    );

    let local = Pairwise::smith_waterman()
        .run_score(&first, &second)
        .unwrap();
    assert_eq!(
        *local.result(),
        Pairwise::smith_waterman().score(&first, &second).unwrap()
    );
    assert!(format!("{}", local.report()).contains("Mode          : local linear score-only"));

    let affine = Pairwise::gotoh().run_score(&first, &second).unwrap();
    assert_eq!(
        *affine.result(),
        Pairwise::gotoh().score(&first, &second).unwrap()
    );
    assert!(format!("{}", affine.report()).contains("Mode          : global affine score-only"));

    let wavefront = Pairwise::wfa().run_penalty(&first, &second).unwrap();
    assert_eq!(wavefront.result().penalty(), 16);
    let wavefront_detail = format!("{}", wavefront.report().verbosity(Verbosity::Detailed));
    assert!(wavefront_detail.contains("Penalty            : 16"));
    assert!(wavefront_detail.contains("Alignment score    : -16"));
    assert!(wavefront_detail.contains("Computed offsets"));
    assert_eq!(
        wavefront.report().color(ColorMode::Auto).options().color(),
        ColorMode::Auto
    );

    let xdrop = Pairwise::xdrop().run_score(&first, &second).unwrap();
    assert_eq!(
        xdrop.result(),
        &Pairwise::xdrop().score(&first, &second).unwrap()
    );
    let xdrop_detail = format!("{}", xdrop.report().verbosity(Verbosity::Detailed));
    assert!(xdrop_detail.contains(": forward score-only from origin"));
    assert!(xdrop_detail.contains("First native end"));
    assert!(xdrop_detail.contains("(exclusive)"));
    assert!(xdrop_detail.contains("Termination"));
    assert!(!xdrop_detail.ends_with('\n'));

    let seed = XDropSeed::new(0, 0);
    let seeded = Pairwise::xdrop()
        .run_score_from_seed(&first, &second, seed, XDropDirection::Both)
        .unwrap();
    assert_eq!(
        seeded.result(),
        &Pairwise::xdrop()
            .score_from_seed(&first, &second, seed, XDropDirection::Both)
            .unwrap()
    );
    let seeded_detail = format!("{}", seeded.report().verbosity(Verbosity::Detailed));
    assert!(seeded_detail.contains(": seeded two-way score-only"));
    assert!(seeded_detail.contains("First native range"));
    assert!(!seeded_detail.contains("CIGAR"));

    for (direction, expected_mode) in [
        (XDropDirection::Upstream, "seeded upstream score-only"),
        (XDropDirection::Downstream, "seeded downstream score-only"),
        (XDropDirection::Both, "seeded two-way score-only"),
    ] {
        let run = Pairwise::xdrop()
            .run_score_from_seed(&first, &second, XDropSeed::new(2, 2), direction)
            .unwrap();
        assert_eq!(run.info().mode(), expected_mode);
    }
}

proptest! {
    #[test]
    fn pairwise_reports_preserve_traceback_accounting_for_arbitrary_short_dna(
        first_symbols in prop::collection::vec(prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')], 0..20),
        second_symbols in prop::collection::vec(prop_oneof![Just('A'), Just('C'), Just('G'), Just('T')], 0..20),
        width in 1_usize..25,
    ) {
        let first_text = first_symbols.into_iter().collect::<String>();
        let second_text = second_symbols.into_iter().collect::<String>();
        let first = DnaSequence::try_from(first_text.as_str()).unwrap();
        let second = DnaSequence::try_from(second_text.as_str()).unwrap();
        let run = Pairwise::needleman_wunsch().run(&first, &second).unwrap();
        let statistics = run.statistics();

        prop_assert_eq!(
            statistics.columns(),
            statistics.matches()
                + statistics.substitutions()
                + statistics.insertion_residues()
                + statistics.deletion_residues()
        );
        let report = format!("{}", run.report().block_width(width).unwrap());
        prop_assert!(report.contains("Algorithm     : Needleman-Wunsch"));
        prop_assert!(!report.contains("\x1b["));
        prop_assert_eq!(report, format!("{}", run.report().block_width(width).unwrap()));
    }
}
