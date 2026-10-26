//! Structural tests for shared multiple-alignment primitives.

use align_core::{Dna, DnaSequence, Sequence, Symbol};
use align_msa::{
    AlignedRow, DistanceMatrix, GuideTree, GuideTreeEdge, GuideTreeNode, MsaError,
    MultipleAlignment, Profile, SequenceRecord,
};

fn row(identifier: &str, aligned: &str) -> AlignedRow<Dna> {
    let symbols = aligned
        .bytes()
        .map(|byte| {
            if byte == b'-' {
                None
            } else {
                Some(Symbol::<Dna>::from_ascii(byte).expect("test symbol is valid DNA"))
            }
        })
        .collect();
    AlignedRow::new(identifier, symbols).expect("test identifier is valid")
}

#[test]
fn sequence_records_validate_identifiers_without_rewriting_sequences() {
    let sequence = DnaSequence::from_ascii(b"ACGT").unwrap();
    let record = SequenceRecord::new("sample 1", sequence.clone()).unwrap();
    assert_eq!(record.identifier(), "sample 1");
    assert_eq!(record.sequence(), &sequence);
    assert_eq!(
        SequenceRecord::<Dna>::new("  ", sequence.clone()),
        Err(MsaError::EmptyIdentifier)
    );
    assert!(matches!(
        SequenceRecord::<Dna>::new("bad\nname", sequence),
        Err(MsaError::ControlCharacterInIdentifier { index: 3, .. })
    ));
}

#[test]
fn multiple_alignment_enforces_rectangle_identity_and_nonempty_columns() {
    let alignment = MultipleAlignment::try_new(vec![row("a", "A-C"), row("b", "ATC")])
        .expect("valid alignment");
    assert_eq!(alignment.sequence_count(), 2);
    assert_eq!(alignment.width(), 3);
    assert_eq!(alignment.row(0).unwrap().to_ungapped_ascii(), "AC");
    assert_eq!(
        alignment
            .column(1)
            .expect("column exists")
            .map(|cell| cell.map(Symbol::as_ascii))
            .collect::<Vec<_>>(),
        vec![None, Some(b'T')]
    );

    assert_eq!(
        MultipleAlignment::try_new(vec![row("a", "A-C"), row("b", "AT")]),
        Err(MsaError::RowWidthMismatch {
            row: 1,
            expected: 3,
            actual: 2,
        })
    );
    assert_eq!(
        MultipleAlignment::try_new(vec![row("same", "AC"), row("same", "AC")]),
        Err(MsaError::DuplicateIdentifier {
            identifier: "same".to_owned(),
        })
    );
    assert_eq!(
        MultipleAlignment::try_new(vec![row("a", "A-C"), row("b", "T-G")]),
        Err(MsaError::AllGapColumn { column: 1 })
    );
}

#[test]
fn zero_width_alignment_represents_empty_input_sequences() {
    let alignment = MultipleAlignment::try_new(vec![row("a", ""), row("b", "")]).unwrap();
    assert_eq!(alignment.width(), 0);
    assert!(alignment.column(0).is_none());
}

#[test]
fn profiles_require_one_positive_weight_per_row() {
    let alignment = MultipleAlignment::try_new(vec![row("a", "AC"), row("b", "AT")]).unwrap();
    let profile = Profile::try_new(alignment.clone(), vec![1_000_000, 500_000]).unwrap();
    assert_eq!(profile.weights(), &[1_000_000, 500_000]);
    assert_eq!(
        Profile::try_new(alignment.clone(), vec![1]),
        Err(MsaError::WeightCountMismatch {
            rows: 2,
            weights: 1,
        })
    );
    assert_eq!(
        Profile::try_new(alignment, vec![1, 0]),
        Err(MsaError::ZeroWeight { index: 1 })
    );
}

#[test]
fn distance_matrices_validate_shape_range_diagonal_and_symmetry() {
    let matrix = DistanceMatrix::try_from_flat(2, vec![0.0, 0.25, 0.25, 0.0]).unwrap();
    assert_eq!(matrix.size(), 2);
    assert_eq!(matrix.get(0, 1), Some(0.25));
    assert_eq!(matrix.get(2, 0), None);

    assert!(matches!(
        DistanceMatrix::try_from_flat(2, vec![0.0, 0.25, 0.5, 0.0]),
        Err(MsaError::AsymmetricDistance { .. })
    ));
    assert!(matches!(
        DistanceMatrix::try_from_flat(1, vec![f64::NAN]),
        Err(MsaError::InvalidDistance { .. })
    ));
}

#[test]
fn guide_tree_contains_each_leaf_exactly_once() {
    let left = GuideTreeEdge::new(0.1, GuideTreeNode::leaf(0)).unwrap();
    let right_child = GuideTreeNode::branch(
        GuideTreeEdge::new(0.2, GuideTreeNode::leaf(1)).unwrap(),
        GuideTreeEdge::new(0.3, GuideTreeNode::leaf(2)).unwrap(),
    );
    let right = GuideTreeEdge::new(0.4, right_child).unwrap();
    let tree = GuideTree::try_new(GuideTreeNode::branch(left, right), 3).unwrap();
    assert_eq!(tree.postorder_leaves(), vec![0, 1, 2]);
    assert_eq!(tree.root().minimum_leaf_index(), 0);

    let duplicate = GuideTreeNode::branch(
        GuideTreeEdge::new(0.0, GuideTreeNode::leaf(0)).unwrap(),
        GuideTreeEdge::new(0.0, GuideTreeNode::leaf(0)).unwrap(),
    );
    assert_eq!(
        GuideTree::try_new(duplicate, 2),
        Err(MsaError::DuplicateLeaf { index: 0 })
    );
}

#[test]
fn aligned_rows_round_trip_typed_cells() {
    let source = Sequence::<Dna>::from_ascii(b"ACGT").unwrap();
    let cells = source.iter().map(Some).collect::<Vec<_>>();
    let row = AlignedRow::new("dna", cells).unwrap();
    assert_eq!(row.to_ascii_with_gaps(), "ACGT");
    let (identifier, cells) = row.into_parts();
    assert_eq!(identifier, "dna");
    assert_eq!(cells.len(), 4);
}
