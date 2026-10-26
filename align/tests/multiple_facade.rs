//! Exact equivalence tests for the multiple-alignment facade.

use std::io::Cursor;

use align::prelude::*;
use align_clustal::{ClustalOmega, ClustalOmegaParameters};
use align_mafft::{MafftEinsi, MafftFftNs1, MafftFftNs2, MafftFftNsi, MafftGinsi, MafftLinsi};
use align_mafft_compat::{Mafft7526Einsi, Mafft7526Linsi};
use align_muscle::Muscle5;

fn dna_records() -> Vec<SequenceRecord<Dna>> {
    [
        ("first", "ACGTACGT"),
        ("second", "ACGTTCGT"),
        ("third", "ACGTCGT"),
    ]
    .into_iter()
    .map(|(identifier, sequence)| {
        SequenceRecord::new(identifier, DnaSequence::try_from(sequence).unwrap()).unwrap()
    })
    .collect()
}

fn protein_records() -> Vec<SequenceRecord<Protein>> {
    [
        ("first", "MKTAYIAKQRQ"),
        ("second", "MKTAYIAKQKQ"),
        ("third", "MKTAYIAQRQ"),
    ]
    .into_iter()
    .map(|(identifier, sequence)| {
        SequenceRecord::new(identifier, ProteinSequence::try_from(sequence).unwrap()).unwrap()
    })
    .collect()
}

#[test]
fn fasta_is_available_from_the_prelude() {
    let records = parse_fasta::<Dna>(Cursor::new(b">a\nacgt\n>b\nagt\n")).unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0].sequence().to_string(), "ACGT");
    let alignment = records
        .align_with(Multiple::mafft(MafftStrategy::FftNs1))
        .unwrap();
    assert_eq!(alignment.alignment().sequence_count(), 2);
}

#[test]
fn clustal_and_muscle_facades_match_published_kernels() {
    let records = dna_records();

    let facade_clustal = Multiple::clustal_omega()
        .threads(2)
        .align(&records)
        .unwrap();
    let parameters = ClustalOmegaParameters::nucleotide_v1_2_4()
        .with_thread_count(2)
        .unwrap();
    let direct_clustal = ClustalOmega::<Dna>::with_parameters(parameters)
        .unwrap()
        .align(&records)
        .unwrap();
    assert_eq!(facade_clustal.alignment(), direct_clustal.alignment());

    let facade_muscle = Multiple::muscle5().align(&records).unwrap();
    let direct_muscle = Muscle5::<Dna>::published_align().align(&records).unwrap();
    assert_eq!(facade_muscle.alignment(), direct_muscle.alignment());
    assert!(facade_muscle.alignment().to_string().contains("first"));

    let facade_super5 = Multiple::muscle5_super5().align(&records).unwrap();
    let direct_super5 = Muscle5::<Dna>::published_super5().align(&records).unwrap();
    assert_eq!(facade_super5.alignment(), direct_super5.alignment());
}

#[test]
fn all_six_mafft_strategies_match_their_frozen_entry_points() {
    let records = dna_records();

    let cases = [
        (
            MafftStrategy::FftNs1,
            MafftFftNs1::<Dna>::new()
                .align(&records)
                .unwrap()
                .into_alignment(),
        ),
        (
            MafftStrategy::FftNs2,
            MafftFftNs2::<Dna>::new()
                .align(&records)
                .unwrap()
                .into_alignment(),
        ),
        (
            MafftStrategy::FftNsi,
            MafftFftNsi::<Dna>::new()
                .align(&records)
                .unwrap()
                .into_alignment(),
        ),
        (
            MafftStrategy::Ginsi,
            MafftGinsi::<Dna>::new()
                .align(&records)
                .unwrap()
                .into_alignment(),
        ),
        (
            MafftStrategy::Linsi,
            MafftLinsi::<Dna>::new()
                .align(&records)
                .unwrap()
                .into_alignment(),
        ),
        (
            MafftStrategy::Einsi,
            MafftEinsi::<Dna>::new()
                .align(&records)
                .unwrap()
                .into_alignment(),
        ),
    ];

    for (strategy, direct) in cases {
        let facade = Multiple::mafft(strategy).align(&records).unwrap();
        assert_eq!(facade.strategy(), strategy);
        assert_eq!(facade.profile_identifier(), strategy.profile_identifier());
        assert_eq!(facade.alignment(), &direct);
    }
}

#[test]
fn mafft_threads_are_explicit_and_never_silently_ignored() {
    let records = dna_records();
    let facade = Multiple::mafft(MafftStrategy::FftNsi)
        .threads(2)
        .align(&records)
        .unwrap();
    let direct = MafftFftNsi::<Dna>::new()
        .with_thread_count(2)
        .unwrap()
        .align(&records)
        .unwrap();
    assert_eq!(facade.alignment(), direct.alignment());

    let facade = Multiple::mafft(MafftStrategy::Linsi)
        .threads(2)
        .align(&records)
        .unwrap();
    let direct = MafftLinsi::<Dna>::new()
        .with_thread_count(2)
        .unwrap()
        .align(&records)
        .unwrap();
    assert_eq!(facade.alignment(), direct.alignment());

    assert!(matches!(
        Multiple::mafft(MafftStrategy::Ginsi)
            .threads(2)
            .align(&records),
        Err(Error::UnsupportedOption {
            option: "threads",
            ..
        })
    ));
}

#[test]
fn mafft_7526_facades_match_the_pinned_compatibility_entry_points() {
    let dna = dna_records();
    {
        let facade = dna
            .align_with(Multiple::mafft_7526_linsi().threads(2))
            .unwrap();
        let direct = Mafft7526Linsi::<Dna>::new()
            .with_thread_count(2)
            .unwrap()
            .align(&dna)
            .unwrap();
        assert_eq!(facade.alignment(), direct.alignment());
        assert_eq!(facade.profile(), Mafft7526Profile::Linsi);
    }

    let protein = protein_records();
    {
        let facade = protein
            .align_with(Multiple::mafft_7526_einsi().threads(2))
            .unwrap();
        let direct = Mafft7526Einsi::new()
            .with_thread_count(2)
            .unwrap()
            .align(&protein)
            .unwrap();
        assert_eq!(facade.alignment(), direct.alignment());
        assert_eq!(facade.profile(), Mafft7526Profile::Einsi);
    }
}

#[test]
fn mafft_7526_rejects_zero_workers_through_the_facade() {
    assert!(matches!(
        Multiple::mafft_7526_linsi()
            .threads(0)
            .align(&dna_records()),
        Err(Error::Mafft(align_mafft::MafftError::InvalidThreadCount {
            value: 0
        }))
    ));
}
