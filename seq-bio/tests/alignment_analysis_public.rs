//! Checked analysis through public import, mapping and accepted-engine adapters.

use bioferris::io::alignment::{AlignmentFormat, parse_alignment};
use bioferris::prelude::*;

fn input() -> AlignmentAnalysis<Dna> {
    let document = parse_alignment::<Dna>(
        include_bytes!("../tests/fixtures/alignment-analysis.fa").as_slice(),
        AlignmentFormat::Fasta,
    )
    .unwrap();
    AlignmentAnalysis::from_document(&document).unwrap()
}

#[test]
fn imported_statistics_have_independently_counted_cells_pairs_and_gap_runs() {
    let alignment = input();
    let stats = alignment
        .statistics(AlignmentStatisticsOptions::default())
        .unwrap();
    assert_eq!((stats.rows, stats.columns), (3, 5));
    assert_eq!(
        stats.cells,
        AlignmentCellCounts {
            residues: 10,
            gaps: 5,
            missing: 0,
            ambiguous: 1,
        }
    );
    assert_eq!(
        stats.gap_runs,
        AlignmentGapRuns {
            leading: 0,
            internal: 3,
            trailing: 1,
        }
    );
    assert_eq!(
        stats.pairs,
        AlignmentPairCounts {
            total: 15,
            missing: 0,
            excluded_ambiguity: 0,
            double_gap: 3,
            single_gap: 4,
            residue_pairs: 8,
            identities: 7,
            mismatches: 1,
        }
    );
    assert_eq!(stats.denominator, 8);
    assert_eq!(stats.identity.unwrap().to_bits(), 0.875_f64.to_bits());
    let ignored = alignment
        .statistics(AlignmentStatisticsOptions {
            ambiguity: AlignmentAmbiguityPolicy::Ignore,
            denominator: IdentityDenominator::ResidueAndSingleGap,
        })
        .unwrap();
    assert_eq!(ignored.pairs.excluded_ambiguity, 2);
    assert_eq!(ignored.denominator, 10);
    assert_eq!(ignored.identity.unwrap().to_bits(), 0.7_f64.to_bits());
}

#[test]
fn reverse_origins_and_slices_retain_full_source_coordinates() {
    let alignment = input();
    let oriented = alignment
        .with_row_coordinates(0, RowCoordinates::new(20, 10..14, Strand::Reverse).unwrap())
        .unwrap();
    let positions = (0..5)
        .map(|column| oriented.column_to_sequence(0, column).unwrap())
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        vec![Some(13), Some(12), None, Some(11), Some(10)]
    );
    let slice = oriented.slice_columns(1..4).unwrap();
    assert_eq!(slice.rows()[0].aligned(), "C-G");
    let origin = slice.rows()[0].coordinates().unwrap();
    assert_eq!(origin.source_length(), 20);
    assert_eq!(origin.range(), 11..13);
    assert_eq!(origin.strand(), Strand::Reverse);
    assert_eq!(slice.column_to_sequence(0, 0).unwrap(), Some(12));
    assert_eq!(slice.sequence_to_column(0, 11).unwrap(), Some(2));
    assert_eq!(slice.sequence_to_column(0, 10).unwrap(), None);
    assert!(matches!(
        slice.sequence_to_column(0, 20),
        Err(AlignmentAnalysisError::InvalidSourcePosition { .. })
    ));
    assert_eq!(alignment.rows()[0].aligned(), "AC-GT");
    assert_eq!(alignment.rows()[0].coordinates().unwrap().range(), 0..4);
}

#[test]
fn missing_cells_never_gain_inferred_coordinates_after_slicing() {
    let document = parse_alignment::<Dna>(
        include_bytes!("../tests/fixtures/alignment-analysis.nex").as_slice(),
        AlignmentFormat::Nexus,
    )
    .unwrap();
    let alignment = AlignmentAnalysis::from_document(&document).unwrap();
    assert_eq!(alignment.rows()[0].aligned(), "A~!C");
    assert!(matches!(
        alignment.column_to_sequence(0, 0),
        Err(AlignmentAnalysisError::MissingCoordinates { row: 0 })
    ));
    let slice = alignment.slice_columns(0..2).unwrap();
    assert_eq!(slice.rows()[0].aligned(), "A~");
    assert!(slice.rows()[0].coordinates().is_none());
    assert!(matches!(
        slice.column_to_sequence(0, 0),
        Err(AlignmentAnalysisError::MissingCoordinates { row: 0 })
    ));
    let asserted = slice
        .with_row_coordinates(0, RowCoordinates::new(7, 4..5, Strand::Forward).unwrap())
        .unwrap();
    assert_eq!(asserted.column_to_sequence(0, 0).unwrap(), Some(4));
    assert_eq!(asserted.column_to_sequence(0, 1).unwrap(), None);
    assert_eq!(alignment.rows()[0].aligned(), "A~!C");
}

#[test]
fn annotated_slices_refuse_loss_or_return_the_complete_ordered_report() {
    let document = parse_alignment::<Dna>(
        include_bytes!("../tests/fixtures/alignment-analysis.sto").as_slice(),
        AlignmentFormat::Stockholm,
    )
    .unwrap();
    let alignment = AlignmentAnalysis::from_document(&document).unwrap();
    assert_eq!(alignment.annotations().len(), 5);
    assert!(matches!(
        alignment.slice_columns(1..3),
        Err(AlignmentAnalysisError::AnnotationLoss { .. })
    ));
    let identity = alignment
        .slice_with(0..2, 0..4, AlignmentAnnotationPolicy::Discard)
        .unwrap();
    assert_eq!(identity.alignment().annotations(), document.annotations());
    assert!(identity.discarded_annotations().is_empty());
    let slice = alignment
        .slice_with(0..2, 1..3, AlignmentAnnotationPolicy::Discard)
        .unwrap();
    assert_eq!(slice.discarded_annotations(), document.annotations());
    assert!(slice.alignment().annotations().is_empty());
    assert_eq!(slice.alignment().rows()[0].aligned(), "C-");
    assert_eq!(slice.alignment().rows()[1].aligned(), "CT");
}

#[test]
fn consensus_policies_distinguish_tie_resolution_thresholds_and_no_evidence() {
    let alignment = input();
    let consensus = alignment.consensus(ConsensusOptions::default()).unwrap();
    assert_eq!(consensus.sequence, "ANNGT");
    assert_eq!(consensus.columns[1].status, ConsensusStatus::Tie);
    assert_eq!(consensus.columns[2].status, ConsensusStatus::NoEvidence);
    let options = ConsensusOptions {
        gaps: ConsensusGapPolicy::Vote,
        ties: ConsensusTiePolicy::AsciiOrder,
        min_support: 0.3,
        ..ConsensusOptions::default()
    };
    let resolved = alignment.consensus(options).unwrap();
    assert_eq!(resolved.sequence, "A--GT");
    assert_eq!(resolved.columns[1].status, ConsensusStatus::TieResolved);
    assert_eq!(resolved.columns[1].winners, b"-CN");
    let rejected = alignment
        .consensus(ConsensusOptions {
            min_support: 0.5,
            ..options
        })
        .unwrap();
    assert_eq!(rejected.sequence, "AN-GT");
    assert_eq!(rejected.columns[1].status, ConsensusStatus::BelowThreshold);
    let empty = alignment.slice_rows(0..0).unwrap();
    assert_eq!(empty.sequence_count(), 0);
    assert_eq!(empty.width(), 5);
    let consensus = empty.consensus(ConsensusOptions::default()).unwrap();
    assert_eq!(consensus.sequence, "NNNNN");
    assert!(consensus.columns.iter().all(|column| {
        column.status == ConsensusStatus::NoEvidence && column.support.is_none()
    }));
    assert!(
        empty
            .statistics(AlignmentStatisticsOptions::default())
            .unwrap()
            .identity
            .is_none()
    );
}

#[test]
fn local_pairwise_adapter_retains_the_accepted_traceback_and_original_ranges() {
    let first = dna("TTACGTAA").unwrap();
    let second = dna("GGACGTCC").unwrap();
    let run = Pairwise::smith_waterman().run(&first, &second).unwrap();
    let analysis = AlignmentAnalysis::from_pairwise(run.first(), run.second(), &run).unwrap();
    assert_eq!(analysis.width(), 4);
    for (row, original) in [first, second].iter().enumerate() {
        assert_eq!(analysis.rows()[row].aligned(), "ACGT");
        let origin = analysis.rows()[row].coordinates().unwrap();
        assert_eq!(origin.source_length(), original.len());
        assert_eq!(origin.range(), 2..6);
        assert_eq!(analysis.column_to_sequence(row, 0).unwrap(), Some(2));
        assert_eq!(analysis.sequence_to_column(row, 0).unwrap(), None);
    }
}

#[test]
fn dimensions_and_strict_ranges_fail_without_changing_the_imported_document() {
    let document = parse_alignment::<Dna>(
        include_bytes!("../tests/fixtures/alignment-analysis.fa").as_slice(),
        AlignmentFormat::Fasta,
    )
    .unwrap();
    assert!(matches!(
        AlignmentAnalysis::from_document_with_limits(
            &document,
            AlignmentAnalysisLimits {
                max_cells: 14,
                ..AlignmentAnalysisLimits::default()
            }
        ),
        Err(AlignmentAnalysisError::LimitExceeded {
            resource: "cells",
            observed: 15,
            ..
        })
    ));
    let alignment = AlignmentAnalysis::from_document(&document).unwrap();
    assert!(matches!(
        alignment.slice_columns(0..6),
        Err(AlignmentAnalysisError::InvalidRange { .. })
    ));
    assert_eq!(document.rows()[0].aligned(), "AC-GT");
    assert_eq!(alignment.rows()[0].aligned(), "AC-GT");
}
