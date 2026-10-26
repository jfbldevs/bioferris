//! Checked alignment statistics, source coordinates, consensus and annotation loss.

use bioferris::io::alignment::{AlignmentFormat, parse_alignment};
use bioferris::prelude::*;

const FASTA: &str = ">reference\nAC-GT\n>sample\nA--GT\n>ambiguous\nAN-G-\n";
const STOCKHOLM: &str = "# STOCKHOLM 1.0\n#=GF ID original\n# comment with coordinates 1-4\n#=GS reference DE original reference row\nreference AC-G\nsample ACTG\n#=GC SS_cons <..>\n#=GR reference SS <..>\n//\n";
const NEXUS: &str = "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=4;\nFORMAT DATATYPE=DNA GAP=~ MISSING=!;\nMATRIX\na A~!C\nb CG!C\n;\nEND;\n";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let document = parse_alignment::<Dna>(FASTA.as_bytes(), AlignmentFormat::Fasta)?;
    let alignment = AlignmentAnalysis::from_document(&document)?;
    report_statistics(&alignment)?;
    report_coordinates_and_slice(&alignment)?;
    report_consensus(&alignment)?;
    report_missing_data()?;
    report_annotation_loss()?;
    report_local_pairwise_origin()?;
    Ok(())
}

fn report_statistics(alignment: &AlignmentAnalysis<Dna>) -> Result<(), Box<dyn std::error::Error>> {
    let stats = alignment.statistics(AlignmentStatisticsOptions::default())?;
    println!("Rectangle: {} rows x {} columns", stats.rows, stats.columns);
    println!("Cells: {:?}", stats.cells);
    println!("Gap runs: {:?}", stats.gap_runs);
    println!("Pairs: {:?}", stats.pairs);
    println!("Identity: {:?} / {}", stats.identity, stats.denominator);
    let ignored = alignment.statistics(AlignmentStatisticsOptions {
        ambiguity: AlignmentAmbiguityPolicy::Ignore,
        denominator: IdentityDenominator::ResidueAndSingleGap,
    })?;
    println!(
        "Ignore ambiguity, include single gaps: {:?} / {}",
        ignored.identity, ignored.denominator
    );
    Ok(())
}

fn report_coordinates_and_slice(
    alignment: &AlignmentAnalysis<Dna>,
) -> Result<(), Box<dyn std::error::Error>> {
    let oriented =
        alignment.with_row_coordinates(0, RowCoordinates::new(20, 10..14, Strand::Reverse)?)?;
    let positions: Result<Vec<_>, _> = (0..oriented.width())
        .map(|column| oriented.column_to_sequence(0, column))
        .collect();
    println!("Reverse source positions: {:?}", positions?);
    println!(
        "Source position 12 -> column {:?}",
        oriented.sequence_to_column(0, 12)?
    );
    println!(
        "Unaligned source position 0 -> {:?}",
        oriented.sequence_to_column(0, 0)?
    );
    let sliced = oriented.slice_with(0..2, 1..4, AlignmentAnnotationPolicy::Reject)?;
    for row in sliced.alignment().rows() {
        println!(
            "Slice {}: {} {:?}",
            row.identifier(),
            row.aligned(),
            row.coordinates()
        );
    }
    Ok(())
}

fn report_consensus(alignment: &AlignmentAnalysis<Dna>) -> Result<(), Box<dyn std::error::Error>> {
    println!(
        "Default consensus: {}",
        alignment.consensus(ConsensusOptions::default())?.sequence
    );
    let voted = alignment.consensus(ConsensusOptions {
        gaps: ConsensusGapPolicy::Vote,
        ..Default::default()
    })?;
    println!("Gap-voting consensus: {}", voted.sequence);
    for column in &voted.columns {
        println!(
            "Column {}: {:?}, votes {}/{}, winners {:?}",
            column.column,
            column.status,
            column.winning_count,
            column.denominator,
            String::from_utf8(column.winners.clone())?
        );
    }
    for min_support in [0.3, 0.5] {
        let resolved = alignment.consensus(ConsensusOptions {
            gaps: ConsensusGapPolicy::Vote,
            ties: ConsensusTiePolicy::AsciiOrder,
            min_support,
            ..Default::default()
        })?;
        println!(
            "ASCII ties at support {min_support}: {}, column 1 {:?}",
            resolved.sequence, resolved.columns[1].status
        );
    }
    println!(
        "Consensus ignoring ambiguity: {}",
        alignment
            .consensus(ConsensusOptions {
                ambiguity: AlignmentAmbiguityPolicy::Ignore,
                ..Default::default()
            })?
            .sequence
    );
    report_consensus_rejections(alignment);
    let empty = alignment.slice_rows(0..0)?;
    println!(
        "Zero-row slice: {} x {}, consensus {}",
        empty.sequence_count(),
        empty.width(),
        empty.consensus(ConsensusOptions::default())?.sequence
    );
    Ok(())
}

fn report_consensus_rejections(alignment: &AlignmentAnalysis<Dna>) {
    for (label, options) in [
        (
            "Reject gaps",
            ConsensusOptions {
                gaps: ConsensusGapPolicy::Reject,
                ..Default::default()
            },
        ),
        (
            "Reject ambiguity",
            ConsensusOptions {
                ambiguity: AlignmentAmbiguityPolicy::Reject,
                ..Default::default()
            },
        ),
        (
            "Reject ties",
            ConsensusOptions {
                ties: ConsensusTiePolicy::Error,
                ..Default::default()
            },
        ),
        (
            "Invalid support",
            ConsensusOptions {
                min_support: 1.1,
                ..Default::default()
            },
        ),
    ] {
        if let Err(error) = alignment.consensus(options) {
            println!("{label}: {}: {error}", error.code());
        }
    }
}

fn report_missing_data() -> Result<(), Box<dyn std::error::Error>> {
    let nexus = parse_alignment::<Dna>(NEXUS.as_bytes(), AlignmentFormat::Nexus)?;
    let missing = AlignmentAnalysis::from_document(&nexus)?;
    println!(
        "Custom glyphs: {} / {}",
        missing.rows()[0].aligned(),
        missing.rows()[1].aligned()
    );
    if let Err(error) = missing.column_to_sequence(0, 0) {
        println!("Missing origin: {}: {error}", error.code());
    }
    println!(
        "Consensus ignoring missing data: {}",
        missing.consensus(ConsensusOptions::default())?.sequence
    );
    if let Err(error) = missing.consensus(ConsensusOptions {
        missing: ConsensusMissingPolicy::Reject,
        ..Default::default()
    }) {
        println!("Reject missing data: {}: {error}", error.code());
    }
    Ok(())
}

fn report_annotation_loss() -> Result<(), Box<dyn std::error::Error>> {
    let annotated = parse_alignment::<Dna>(STOCKHOLM.as_bytes(), AlignmentFormat::Stockholm)?;
    let annotated = AlignmentAnalysis::from_document(&annotated)?;
    if let Err(error) = annotated.slice_columns(1..3) {
        println!("Strict annotation slice: {}: {error}", error.code());
    }
    let discarded = annotated.slice_with(0..2, 1..3, AlignmentAnnotationPolicy::Discard)?;
    for annotation in discarded.discarded_annotations() {
        println!(
            "Discard: {:?} {:?} {:?} line {} {:?}",
            annotation.scope(),
            annotation.sequence(),
            annotation.tag(),
            annotation.line(),
            annotation.value()
        );
    }
    Ok(())
}

fn report_local_pairwise_origin() -> Result<(), Box<dyn std::error::Error>> {
    let first = dna("TTACGTAA")?;
    let second = dna("GGACGTCC")?;
    let run = Pairwise::smith_waterman().run(&first, &second)?;
    let local = AlignmentAnalysis::from_pairwise(run.first(), run.second(), &run)?;
    println!(
        "Local first-row origin: {:?}",
        local.rows()[0].coordinates()
    );
    println!(
        "Local first column -> source {:?}",
        local.column_to_sequence(0, 0)?
    );
    Ok(())
}
