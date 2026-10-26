//! Golden and boundary tests for human-readable multiple-alignment output.

use align::prelude::*;

fn protein_records() -> Vec<SequenceRecord<Protein>> {
    [
        ("human", "MKTAYIAKQRQ"),
        ("mouse", "MKTAYIAKQKQ"),
        ("fish", "MKTAYIAQRQ"),
    ]
    .into_iter()
    .map(|(identifier, sequence)| {
        SequenceRecord::new(identifier, ProteinSequence::try_from(sequence).unwrap()).unwrap()
    })
    .collect()
}

#[test]
fn clustal_muscle_and_six_mafft_profiles_have_exact_compact_output() {
    let records = protein_records();
    let clustal = Multiple::clustal_omega().run(&records).unwrap();
    let muscle = Multiple::muscle5().run(&records).unwrap();
    let mut reports = vec![
        format!("{}", clustal.report().verbosity(Verbosity::Compact)),
        format!("{}", muscle.report().verbosity(Verbosity::Compact)),
    ];
    for strategy in [
        MafftStrategy::FftNs1,
        MafftStrategy::FftNs2,
        MafftStrategy::FftNsi,
        MafftStrategy::Ginsi,
        MafftStrategy::Linsi,
        MafftStrategy::Einsi,
    ] {
        reports.push(format!(
            "{}",
            Multiple::mafft(strategy)
                .run(&records)
                .unwrap()
                .report()
                .verbosity(Verbosity::Compact)
        ));
    }
    let expected = [
        "Algorithm      : Clustal Omega\nProfile        : CLUSTALO-1.0.2-default-protein-v1\nMode           : published progressive HMM\nSequences      : 3\nColumns        : 11",
        "Algorithm      : MUSCLE 5\nProfile        : MUSCLE-5.1.0-Nature-2022-default-v2\nMode           : Align\nSequences      : 3\nColumns        : 11",
        "Algorithm      : MAFFT\nProfile        : MAFFT-2002-FFT-NS-1-v1\nMode           : FFT-NS-1\nSequences      : 3\nColumns        : 11",
        "Algorithm      : MAFFT\nProfile        : MAFFT-2002-FFT-NS-2-v1\nMode           : FFT-NS-2\nSequences      : 3\nColumns        : 11",
        "Algorithm      : MAFFT\nProfile        : MAFFT-2002-FFT-NS-i-v2\nMode           : FFT-NS-i\nSequences      : 3\nColumns        : 11",
        "Algorithm      : MAFFT\nProfile        : MAFFT-2005-G-INS-i-v1\nMode           : G-INS-i\nSequences      : 3\nColumns        : 11",
        "Algorithm      : MAFFT\nProfile        : MAFFT-2005-L-INS-i-v1\nMode           : L-INS-i\nSequences      : 3\nColumns        : 11",
        "Algorithm      : MAFFT\nProfile        : MAFFT-2008-E-INS-i-v1\nMode           : E-INS-i\nSequences      : 3\nColumns        : 11",
    ];
    assert_eq!(reports, expected);
    assert!(reports.iter().all(|report| !report.contains("\u{1b}[")));
    assert!(reports.iter().all(|report| !report.ends_with('\n')));
    assert_eq!(
        clustal.result().alignment(),
        Multiple::clustal_omega()
            .align(&records)
            .unwrap()
            .alignment()
    );
}

#[test]
fn normal_msa_output_has_blocks_coordinates_conservation_and_gaps() {
    let records = protein_records();
    let run = Multiple::mafft(MafftStrategy::Linsi).run(&records).unwrap();
    let output = format!("{}", run.report());
    assert!(output.contains("Conserved      : 7"));
    assert!(output.contains("Gap characters : 1"));
    assert!(output.contains("human   1  MKTAYIAKQRQ  11"));
    assert!(output.contains("fish    1  MKTAYIAQRQ-  10"));
    assert!(output.contains("*******"));

    for width in [1, 3, 60, 100] {
        let split = format!("{}", run.report().block_width(width).unwrap());
        assert!(split.contains("human"));
        assert!(split.contains("Legend:"));
    }
    assert!(run.report().block_width(0).is_err());
}

#[test]
fn detailed_outputs_retain_family_specific_evidence() {
    let records = protein_records();
    let clustal = format!(
        "{}",
        Multiple::clustal_omega()
            .run(&records)
            .unwrap()
            .report()
            .verbosity(Verbosity::Detailed)
    );
    assert!(clustal.contains("Guide method"));
    assert!(clustal.contains("Progressive merges"));

    let muscle = format!(
        "{}",
        Multiple::muscle5_super5()
            .run(&records)
            .unwrap()
            .report()
            .verbosity(Verbosity::Detailed)
    );
    assert!(muscle.contains("Mode                 : Super5"));
    assert!(muscle.contains("Final objective"));

    for strategy in [
        MafftStrategy::FftNs1,
        MafftStrategy::FftNs2,
        MafftStrategy::FftNsi,
        MafftStrategy::Ginsi,
        MafftStrategy::Linsi,
        MafftStrategy::Einsi,
    ] {
        let output = format!(
            "{}",
            Multiple::mafft(strategy)
                .run(&records)
                .unwrap()
                .report()
                .verbosity(Verbosity::Detailed)
        );
        assert!(output.contains(strategy.profile_identifier()));
        assert!(output.contains("Estimated peak bytes"));
    }
}

#[test]
fn mafft_7526_reports_compatibility_provenance_and_workers() {
    let records = protein_records();
    let linsi = Multiple::mafft_7526_linsi()
        .threads(2)
        .run(&records)
        .unwrap();
    assert_eq!(linsi.info().algorithm(), "MAFFT 7.526");
    assert_eq!(linsi.info().profile(), "MAFFT-7.526-L-INS-i-localpair-v1");
    let output = format!("{}", linsi.report().verbosity(Verbosity::Detailed));
    assert!(output.contains("L-INS-i native-output-compatible"));
    assert!(output.contains("Native target"));
    assert!(output.contains("7.526 (2024/Apr/26)"));
    assert!(output.contains("Rust port"));
    assert!(output.contains("mafft 0.1.2"));
    assert!(output.contains("Workers"));
    assert!(
        output
            .lines()
            .any(|line| line.contains("Workers") && line.ends_with(": 2"))
    );

    let einsi = Multiple::mafft_7526_einsi().run(&records).unwrap();
    assert_eq!(
        einsi.info().profile(),
        "MAFFT-7.526-E-INS-i-genafpair-protein-v1"
    );
}

#[test]
fn identifiers_and_all_three_builtin_alphabets_render_without_truncation() {
    let protein = vec![
        SequenceRecord::new(
            "proteína humana completa",
            ProteinSequence::try_from("MKT").unwrap(),
        )
        .unwrap(),
    ];
    let protein_output = format!(
        "{}",
        Multiple::mafft(MafftStrategy::FftNs1)
            .run(&protein)
            .unwrap()
            .report()
    );
    assert!(protein_output.contains("proteína humana completa"));

    let dna = vec![
        SequenceRecord::new("dna-a", DnaSequence::try_from("ACGT").unwrap()).unwrap(),
        SequenceRecord::new("dna-b", DnaSequence::try_from("AGT").unwrap()).unwrap(),
    ];
    assert!(
        format!(
            "{}",
            Multiple::mafft(MafftStrategy::FftNs1)
                .run(&dna)
                .unwrap()
                .report()
        )
        .contains("dna-a")
    );

    let rna = vec![
        SequenceRecord::new("rna-a", RnaSequence::try_from("ACGU").unwrap()).unwrap(),
        SequenceRecord::new("rna-b", RnaSequence::try_from("AGU").unwrap()).unwrap(),
    ];
    assert!(
        format!(
            "{}",
            Multiple::mafft(MafftStrategy::FftNs1)
                .run(&rna)
                .unwrap()
                .report()
        )
        .contains("ACGU")
    );
}

#[test]
fn clustal_reports_the_published_profile_selected_for_each_molecule_kind() {
    let protein = protein_records();
    assert_eq!(
        Multiple::clustal_omega()
            .run(&protein)
            .unwrap()
            .info()
            .profile(),
        "CLUSTALO-1.0.2-default-protein-v1"
    );

    let dna = vec![
        SequenceRecord::new("dna-a", DnaSequence::try_from("ACGT").unwrap()).unwrap(),
        SequenceRecord::new("dna-b", DnaSequence::try_from("AGT").unwrap()).unwrap(),
    ];
    let rna = vec![
        SequenceRecord::new("rna-a", RnaSequence::try_from("ACGU").unwrap()).unwrap(),
        SequenceRecord::new("rna-b", RnaSequence::try_from("AGU").unwrap()).unwrap(),
    ];
    for profile in [
        Multiple::clustal_omega()
            .run(&dna)
            .unwrap()
            .info()
            .profile(),
        Multiple::clustal_omega()
            .run(&rna)
            .unwrap()
            .info()
            .profile(),
    ] {
        assert_eq!(profile, "CLUSTALO-1.2.4-default-nucleotide-v1");
    }
}

#[test]
fn row_coordinates_can_follow_native_rust_style() {
    let records = protein_records();
    let run = Multiple::mafft(MafftStrategy::FftNs1)
        .run(&records)
        .unwrap();
    let report = format!(
        "{}",
        run.report()
            .coordinates(CoordinateStyle::Rust)
            .block_width(5)
            .unwrap()
    );
    assert!(report.contains("0  MKTAY   5"), "{report}");
    assert!(report.contains("5  IAKQR  10"), "{report}");
}

#[test]
fn msa_human_output_is_independent_of_worker_configuration() {
    let records = protein_records();
    let serial = Multiple::clustal_omega().threads(1).run(&records).unwrap();
    let parallel = Multiple::clustal_omega().threads(2).run(&records).unwrap();
    assert_eq!(
        format!("{}", serial.report()),
        format!("{}", parallel.report())
    );

    let serial = Multiple::mafft(MafftStrategy::FftNsi)
        .threads(1)
        .run(&records)
        .unwrap();
    let parallel = Multiple::mafft(MafftStrategy::FftNsi)
        .threads(2)
        .run(&records)
        .unwrap();
    assert_eq!(
        format!("{}", serial.report()),
        format!("{}", parallel.report())
    );

    let colored = serial.report().color(ColorMode::Always);
    assert_eq!(colored.options().color(), ColorMode::Always);
    assert_eq!(format!("{colored}"), format!("{}", serial.report()));
}
