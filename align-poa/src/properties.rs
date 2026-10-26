use align_core::{AffineGap, Dna, DnaSequence, LinearGap, MatchMismatch, Score};
use align_msa::SequenceRecord;
use align_pairwise::Gotoh;
use proptest::prelude::*;

use crate::{
    AffinePoa, GlobalAffinePoa, LeePoa, LocalPoa, PoaGraph, PoaNodeId, Poasta, PoastaAffineCosts,
    PoastaPoa, PoastaSearchOptions, verify_global_affine_sequence_graph_alignment,
    verify_poasta_alignment, verify_sequence_graph_alignment,
};

const SYMBOLS: &[u8; 4] = b"ACGT";
const MATCH_SCORE: Score = 2;
const MISMATCH_SCORE: Score = -1;
const GAP_SCORE: Score = -2;

fn ascii(indices: &[u8]) -> Vec<u8> {
    indices
        .iter()
        .map(|index| SYMBOLS[usize::from(*index)])
        .collect()
}

fn record(identifier: &str, indices: &[u8]) -> SequenceRecord<Dna> {
    SequenceRecord::new(
        identifier,
        DnaSequence::from_ascii(ascii(indices)).expect("generated DNA is valid"),
    )
    .expect("generated identifier is valid")
}

fn branched_graph(
    paired_symbols: &[(u8, u8)],
    insertion_position: usize,
    insertion_symbol: u8,
) -> PoaGraph<Dna> {
    branched_graph_in_order(paired_symbols, insertion_position, insertion_symbol, false)
}

fn branched_graph_in_order(
    paired_symbols: &[(u8, u8)],
    insertion_position: usize,
    insertion_symbol: u8,
    insertion_first: bool,
) -> PoaGraph<Dna> {
    let backbone = paired_symbols
        .iter()
        .map(|(symbol, _)| *symbol)
        .collect::<Vec<_>>();
    let allele = paired_symbols
        .iter()
        .map(|(_, symbol)| *symbol)
        .collect::<Vec<_>>();
    let mut graph = PoaGraph::new();
    graph
        .append_record(&record("backbone", &backbone), &vec![None; backbone.len()])
        .expect("backbone integration succeeds");
    let backbone_targets = (0..backbone.len())
        .map(|index| Some(PoaNodeId::new(index)))
        .collect::<Vec<_>>();
    let insertion_position = insertion_position % (backbone.len() + 1);
    let mut inserted = Vec::with_capacity(backbone.len() + 1);
    let mut insertion_targets = Vec::with_capacity(backbone.len() + 1);
    for (index, symbol) in backbone.iter().copied().enumerate() {
        if index == insertion_position {
            inserted.push(insertion_symbol);
            insertion_targets.push(None);
        }
        inserted.push(symbol);
        insertion_targets.push(Some(PoaNodeId::new(index)));
    }
    if insertion_position == backbone.len() {
        inserted.push(insertion_symbol);
        insertion_targets.push(None);
    }
    let allele_record = record("allele", &allele);
    let insertion_record = record("insertion", &inserted);
    if insertion_first {
        graph
            .append_record(&insertion_record, &insertion_targets)
            .expect("insertion integration succeeds");
        graph
            .append_record(&allele_record, &backbone_targets)
            .expect("allele integration succeeds");
    } else {
        graph
            .append_record(&allele_record, &backbone_targets)
            .expect("allele integration succeeds");
        graph
            .append_record(&insertion_record, &insertion_targets)
            .expect("insertion integration succeeds");
    }
    graph.validate().expect("generated graph is valid");
    graph
}

fn all_source_to_sink_paths(graph: &PoaGraph<Dna>) -> Vec<Vec<u8>> {
    fn visit(
        graph: &PoaGraph<Dna>,
        node: PoaNodeId,
        prefix: &mut Vec<u8>,
        paths: &mut Vec<Vec<u8>>,
    ) {
        let graph_node = graph.node(node).expect("enumerated node exists");
        prefix.push(graph_node.symbol().as_ascii());
        if graph_node.successors().is_empty() {
            paths.push(prefix.clone());
        } else {
            for successor in graph_node.successors() {
                visit(graph, *successor, prefix, paths);
            }
        }
        prefix.pop();
    }

    let mut paths = Vec::new();
    for node in graph.nodes() {
        if node.predecessors().is_empty() {
            visit(graph, node.id(), &mut Vec::new(), &mut paths);
        }
    }
    paths
}

fn independent_smith_waterman(first: &[u8], second: &[u8]) -> Score {
    let mut previous = vec![0_i64; second.len() + 1];
    let mut current = vec![0_i64; second.len() + 1];
    let mut optimum = 0_i64;
    for first_symbol in first {
        current[0] = 0;
        for (second_index, second_symbol) in second.iter().enumerate() {
            let column = second_index + 1;
            let substitution = if first_symbol == second_symbol {
                MATCH_SCORE
            } else {
                MISMATCH_SCORE
            };
            current[column] = 0
                .max(previous[column - 1] + substitution)
                .max(previous[column] + GAP_SCORE)
                .max(current[column - 1] + GAP_SCORE);
            optimum = optimum.max(current[column]);
        }
        core::mem::swap(&mut previous, &mut current);
    }
    optimum
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum ExhaustiveGap {
    FirstOnly,
    SecondOnly,
}

fn independent_exhaustive_global_affine(
    first: &[u8],
    second: &[u8],
    match_score: Score,
    mismatch_score: Score,
    open_score: Score,
    extension_score: Score,
) -> Score {
    struct Search<'a> {
        first: &'a [u8],
        second: &'a [u8],
        match_score: Score,
        mismatch_score: Score,
        open_score: Score,
        extension_score: Score,
        optimum: i128,
    }

    impl Search<'_> {
        fn visit(
            &mut self,
            first_index: usize,
            second_index: usize,
            active_gap: Option<ExhaustiveGap>,
            score: i128,
        ) {
            if first_index == self.first.len() && second_index == self.second.len() {
                self.optimum = self.optimum.max(score);
                return;
            }
            if first_index < self.first.len() && second_index < self.second.len() {
                let substitution = if self.first[first_index] == self.second[second_index] {
                    self.match_score
                } else {
                    self.mismatch_score
                };
                self.visit(
                    first_index + 1,
                    second_index + 1,
                    None,
                    score + i128::from(substitution),
                );
            }
            if first_index < self.first.len() {
                let increment = if active_gap == Some(ExhaustiveGap::FirstOnly) {
                    self.extension_score
                } else {
                    self.open_score
                };
                self.visit(
                    first_index + 1,
                    second_index,
                    Some(ExhaustiveGap::FirstOnly),
                    score + i128::from(increment),
                );
            }
            if second_index < self.second.len() {
                let increment = if active_gap == Some(ExhaustiveGap::SecondOnly) {
                    self.extension_score
                } else {
                    self.open_score
                };
                self.visit(
                    first_index,
                    second_index + 1,
                    Some(ExhaustiveGap::SecondOnly),
                    score + i128::from(increment),
                );
            }
        }
    }

    let mut search = Search {
        first,
        second,
        match_score,
        mismatch_score,
        open_score,
        extension_score,
        optimum: i128::MIN,
    };
    search.visit(0, 0, None, 0);
    Score::try_from(search.optimum).expect("bounded exhaustive score fits")
}

fn linear_graph(indices: &[u8]) -> PoaGraph<Dna> {
    let mut graph = PoaGraph::new();
    graph
        .append_record(&record("linear", indices), &vec![None; indices.len()])
        .expect("linear graph integration succeeds");
    graph.validate().expect("linear graph is valid");
    graph
}

fn binary_sequences(maximum_len: usize) -> Vec<Vec<u8>> {
    let mut sequences = Vec::new();
    for length in 0..=maximum_len {
        for bits in 0..(1_usize << length) {
            sequences.push(
                (0..length)
                    .map(|offset| u8::from((bits >> offset) & 1 != 0))
                    .collect(),
            );
        }
    }
    sequences
}

#[test]
fn affine_oracle_matches_every_binary_linear_pair_through_length_four() {
    let sequences = binary_sequences(4);
    let scoring = MatchMismatch::new(0, -4);
    let gap = AffineGap::new(-8, -2).unwrap();
    let aligner = GlobalAffinePoa::new(scoring, gap);
    for first in &sequences {
        let graph = linear_graph(first);
        let first_ascii = ascii(first);
        for second in &sequences {
            let query = DnaSequence::from_ascii(ascii(second)).unwrap();
            let alignment = aligner.align(&graph, &query).unwrap();
            let oracle =
                independent_exhaustive_global_affine(&first_ascii, &ascii(second), 0, -4, -8, -2);
            assert_eq!(
                alignment.score(),
                oracle,
                "first={first:?}, second={second:?}"
            );
        }
    }
}

#[test]
fn poasta_modes_match_every_binary_linear_pair_through_length_four() {
    let sequences = binary_sequences(4);
    let costs = PoastaAffineCosts::new(4, 6, 2).unwrap();
    let oracle = GlobalAffinePoa::from_poasta_costs(costs);
    let modes = [
        PoastaSearchOptions::astar_baseline(),
        PoastaSearchOptions::astar_baseline().with_exact_match_dfs(true),
        PoastaSearchOptions::astar_baseline().with_superbubble_pruning(true),
        PoastaSearchOptions::published(),
    ];
    for first in &sequences {
        let graph = linear_graph(first);
        for second in &sequences {
            let query = DnaSequence::from_ascii(ascii(second)).unwrap();
            let expected = oracle.align(&graph, &query).unwrap().score();
            for options in modes {
                let actual = Poasta::with_options(costs, options)
                    .align(&graph, &query)
                    .unwrap();
                assert_eq!(
                    actual.score(),
                    expected,
                    "first={first:?}, second={second:?}, options={options:?}"
                );
            }
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(192))]

    #[test]
    fn dag_dp_matches_exhaustive_path_enumeration(
        paired_symbols in prop::collection::vec((0_u8..4, 0_u8..4), 1..7),
        insertion_position in 0_usize..8,
        insertion_symbol in 0_u8..4,
        query_symbols in prop::collection::vec(0_u8..4, 0..8),
    ) {
        let graph = branched_graph(
            &paired_symbols,
            insertion_position,
            insertion_symbol,
        );
        let query = DnaSequence::from_ascii(ascii(&query_symbols)).unwrap();
        let scoring = MatchMismatch::new(MATCH_SCORE, MISMATCH_SCORE);
        let gap = LinearGap::new(GAP_SCORE).unwrap();
        let alignment = LocalPoa::new(scoring, gap).align(&graph, &query).unwrap();
        let query_ascii = query
            .iter()
            .map(align_core::Symbol::as_ascii)
            .collect::<Vec<_>>();
        let oracle = all_source_to_sink_paths(&graph)
            .iter()
            .map(|path| independent_smith_waterman(path, &query_ascii))
            .max()
            .unwrap_or(0);
        prop_assert_eq!(alignment.score(), oracle);
        prop_assert_eq!(
            verify_sequence_graph_alignment(&alignment, &graph, &query, &scoring, gap)
                .unwrap()
                .recomputed_score(),
            oracle
        );
    }

    #[test]
    fn iterative_construction_is_deterministic_valid_and_source_preserving(
        generated in prop::collection::vec(
            prop::collection::vec(0_u8..4, 0..7),
            1..6,
        ),
    ) {
        let records = generated
            .iter()
            .enumerate()
            .map(|(index, symbols)| record(&format!("sequence-{index}"), symbols))
            .collect::<Vec<_>>();
        let aligner = LeePoa::new(
            MatchMismatch::new(MATCH_SCORE, MISMATCH_SCORE),
            LinearGap::new(GAP_SCORE).unwrap(),
        );
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        prop_assert_eq!(&first, &second);
        first.graph().validate().unwrap();
        prop_assert_eq!(first.graph().path_count(), records.len());
        for (row, source) in first.alignment().rows().iter().zip(&records) {
            prop_assert_eq!(row.identifier(), source.identifier());
            prop_assert_eq!(row.to_ungapped_ascii(), source.sequence().to_string());
        }

        let mut reversed = records.clone();
        reversed.reverse();
        let reordered = aligner.align(&reversed).unwrap();
        reordered.graph().validate().unwrap();
        for (row, source) in reordered.alignment().rows().iter().zip(&reversed) {
            prop_assert_eq!(row.identifier(), source.identifier());
            prop_assert_eq!(row.to_ungapped_ascii(), source.sequence().to_string());
        }
    }

    #[test]
    fn affine_dag_matches_exhaustive_graph_paths_and_alignment_paths(
        paired_symbols in prop::collection::vec((0_u8..4, 0_u8..4), 1..5),
        insertion_position in 0_usize..6,
        insertion_symbol in 0_u8..4,
        query_symbols in prop::collection::vec(0_u8..4, 0..6),
    ) {
        let graph = branched_graph(
            &paired_symbols,
            insertion_position,
            insertion_symbol,
        );
        let query = DnaSequence::from_ascii(ascii(&query_symbols)).unwrap();
        let costs = PoastaAffineCosts::new(4, 6, 2).unwrap();
        let aligner = GlobalAffinePoa::from_poasta_costs(costs);
        let first = aligner.align(&graph, &query).unwrap();
        let second = aligner.align(&graph, &query).unwrap();
        prop_assert_eq!(&first, &second);
        let query_ascii = ascii(&query_symbols);
        let oracle = all_source_to_sink_paths(&graph)
            .iter()
            .map(|path| {
                independent_exhaustive_global_affine(path, &query_ascii, 0, -4, -8, -2)
            })
            .max()
            .expect("generated graph has a source-to-sink path");
        prop_assert_eq!(first.score(), oracle);
        prop_assert_eq!(
            verify_global_affine_sequence_graph_alignment(
                &first,
                &graph,
                &query,
                &costs,
                costs.affine_gap(),
            )
            .unwrap()
            .recomputed_score(),
            oracle
        );
    }

    #[test]
    fn every_poasta_acceleration_mode_matches_the_scalar_oracle(
        paired_symbols in prop::collection::vec((0_u8..4, 0_u8..4), 1..7),
        insertion_position in 0_usize..8,
        insertion_symbol in 0_u8..4,
        query_symbols in prop::collection::vec(0_u8..4, 0..8),
        mismatch in 0_i64..8,
        extension in 0_i64..5,
        opening_delta in 1_i64..6,
    ) {
        let graph = branched_graph(
            &paired_symbols,
            insertion_position,
            insertion_symbol,
        );
        let query = DnaSequence::from_ascii(ascii(&query_symbols)).unwrap();
        let costs = PoastaAffineCosts::new(mismatch, extension + opening_delta, extension).unwrap();
        let expected = GlobalAffinePoa::from_poasta_costs(costs)
            .align(&graph, &query)
            .unwrap()
            .score();
        for options in [
            PoastaSearchOptions::astar_baseline(),
            PoastaSearchOptions::astar_baseline().with_exact_match_dfs(true),
            PoastaSearchOptions::astar_baseline().with_superbubble_pruning(true),
            PoastaSearchOptions::published(),
        ] {
            let aligner = Poasta::with_options(costs, options);
            let first = aligner.align(&graph, &query).unwrap();
            let second = aligner.align(&graph, &query).unwrap();
            prop_assert_eq!(&first, &second);
            prop_assert_eq!(first.score(), expected, "options={:?}", options);
            prop_assert_eq!(
                verify_poasta_alignment(&first, &graph, &query, costs)
                    .unwrap()
                    .recomputed_score(),
                expected
            );
        }
    }

    #[test]
    fn poasta_score_is_invariant_to_equivalent_branch_insertion_order(
        paired_symbols in prop::collection::vec((0_u8..4, 0_u8..4), 1..7),
        insertion_position in 0_usize..8,
        insertion_symbol in 0_u8..4,
        query_symbols in prop::collection::vec(0_u8..4, 0..8),
    ) {
        let first = branched_graph_in_order(
            &paired_symbols,
            insertion_position,
            insertion_symbol,
            false,
        );
        let reordered = branched_graph_in_order(
            &paired_symbols,
            insertion_position,
            insertion_symbol,
            true,
        );
        let query = DnaSequence::from_ascii(ascii(&query_symbols)).unwrap();
        let costs = PoastaAffineCosts::new(4, 6, 2).unwrap();
        let first_alignment = Poasta::new(costs).align(&first, &query).unwrap();
        let reordered_alignment = Poasta::new(costs).align(&reordered, &query).unwrap();
        prop_assert_eq!(first_alignment.score(), reordered_alignment.score());
        prop_assert_eq!(
            first_alignment.score(),
            GlobalAffinePoa::from_poasta_costs(costs)
                .align(&first, &query)
                .unwrap()
                .score()
        );
        prop_assert_eq!(
            reordered_alignment.score(),
            GlobalAffinePoa::from_poasta_costs(costs)
                .align(&reordered, &query)
                .unwrap()
                .score()
        );
    }

    #[test]
    fn affine_linear_graph_reduces_to_global_gotoh(
        first_symbols in prop::collection::vec(0_u8..4, 0..9),
        second_symbols in prop::collection::vec(0_u8..4, 0..9),
    ) {
        let first = DnaSequence::from_ascii(ascii(&first_symbols)).unwrap();
        let second = DnaSequence::from_ascii(ascii(&second_symbols)).unwrap();
        let graph = linear_graph(&first_symbols);
        let scoring = MatchMismatch::new(3, -4);
        let gap = AffineGap::new(-6, -2).unwrap();
        let graph_alignment = GlobalAffinePoa::new(scoring, gap)
            .align(&graph, &second)
            .unwrap();
        let pairwise = Gotoh::new(scoring, gap).align(&first, &second).unwrap();
        prop_assert_eq!(graph_alignment.score(), pairwise.score());
    }

    #[test]
    fn iterative_affine_poa_is_deterministic_valid_and_source_preserving(
        generated in prop::collection::vec(
            prop::collection::vec(0_u8..4, 0..7),
            1..6,
        ),
    ) {
        let records = generated
            .iter()
            .enumerate()
            .map(|(index, symbols)| record(&format!("affine-{index}"), symbols))
            .collect::<Vec<_>>();
        let costs = PoastaAffineCosts::new(4, 6, 2).unwrap();
        let aligner = AffinePoa::from_poasta_costs(costs);
        let first = aligner.align(&records).unwrap();
        let second = aligner.align(&records).unwrap();
        prop_assert_eq!(&first, &second);
        first.graph().validate().unwrap();
        prop_assert_eq!(first.graph().path_count(), records.len());
        prop_assert_eq!(first.steps().len(), records.len());
        for (index, ((row, source), step)) in first
            .alignment()
            .rows()
            .iter()
            .zip(&records)
            .zip(first.steps())
            .enumerate()
        {
            prop_assert_eq!(row.identifier(), source.identifier());
            prop_assert_eq!(row.to_ungapped_ascii(), source.sequence().to_string());
            prop_assert_eq!(step.record_index(), index);
            if let (Some(alignment), Some(verification)) = (
                step.alignment().alignment(),
                step.alignment().verification(),
            ) {
                prop_assert_eq!(alignment.score(), verification.recomputed_score());
                prop_assert_eq!(alignment.query_range(), 0..source.sequence().len());
            } else {
                prop_assert!(step.alignment().is_seed());
            }
        }

        let mut reversed = records.clone();
        reversed.reverse();
        let reordered = aligner.align(&reversed).unwrap();
        reordered.graph().validate().unwrap();
        for (row, source) in reordered.alignment().rows().iter().zip(&reversed) {
            prop_assert_eq!(row.identifier(), source.identifier());
            prop_assert_eq!(row.to_ungapped_ascii(), source.sequence().to_string());
        }
    }

    #[test]
    fn iterative_poasta_is_deterministic_and_matches_independent_scalar_scores(
        generated in prop::collection::vec(
            prop::collection::vec(0_u8..4, 0..7),
            1..6,
        ),
    ) {
        let records = generated
            .iter()
            .enumerate()
            .map(|(index, symbols)| record(&format!("poasta-{index}"), symbols))
            .collect::<Vec<_>>();
        let costs = PoastaAffineCosts::new(4, 6, 2).unwrap();
        let scalar = AffinePoa::from_poasta_costs(costs).align(&records).unwrap();
        let baseline = PoastaPoa::with_options(costs, PoastaSearchOptions::astar_baseline())
            .align(&records)
            .unwrap();
        prop_assert_eq!(baseline.graph(), scalar.graph());
        prop_assert_eq!(baseline.alignment(), scalar.alignment());

        let accelerated_aligner = PoastaPoa::new(costs);
        let accelerated = accelerated_aligner.align(&records).unwrap();
        let repeated = accelerated_aligner.align(&records).unwrap();
        prop_assert_eq!(&accelerated, &repeated);
        accelerated.graph().validate().unwrap();
        prop_assert_eq!(accelerated.graph().path_count(), records.len());
        for (row, source) in accelerated.alignment().rows().iter().zip(&records) {
            prop_assert_eq!(row.identifier(), source.identifier());
            prop_assert_eq!(row.to_ungapped_ascii(), source.sequence().to_string());
        }

        for record_count in 2..=records.len() {
            let step = &accelerated_aligner
                .align(&records[..record_count])
                .unwrap()
                .steps()[record_count - 1]
                .clone();
            let Some(observed) = step.alignment().alignment() else {
                prop_assert!(records[..record_count - 1]
                    .iter()
                    .all(|record| record.sequence().is_empty()));
                continue;
            };
            let before = accelerated_aligner
                .align(&records[..record_count - 1])
                .unwrap();
            let oracle = GlobalAffinePoa::from_poasta_costs(costs)
                .align(before.graph(), records[record_count - 1].sequence())
                .unwrap();
            prop_assert_eq!(observed.score(), oracle.score());
        }
    }

    #[test]
    fn every_iterative_affine_step_matches_independent_path_enumeration(
        generated in prop::collection::vec(
            prop::collection::vec(0_u8..4, 0..5),
            2..5,
        ),
    ) {
        let records = generated
            .iter()
            .enumerate()
            .map(|(index, symbols)| record(&format!("step-{index}"), symbols))
            .collect::<Vec<_>>();
        let costs = PoastaAffineCosts::new(4, 6, 2).unwrap();
        let aligner = AffinePoa::from_poasta_costs(costs);
        for record_count in 2..=records.len() {
            let result = aligner.align(&records[..record_count]).unwrap();
            let step = &result.steps()[record_count - 1];
            let Some(alignment) = step.alignment().alignment() else {
                prop_assert!(records[..record_count - 1]
                    .iter()
                    .all(|record| record.sequence().is_empty()));
                continue;
            };
            let before = aligner.align(&records[..record_count - 1]).unwrap();
            let query_ascii = records[record_count - 1]
                .sequence()
                .iter()
                .map(align_core::Symbol::as_ascii)
                .collect::<Vec<_>>();
            let oracle = all_source_to_sink_paths(before.graph())
                .iter()
                .map(|path| {
                    independent_exhaustive_global_affine(
                        path,
                        &query_ascii,
                        0,
                        -4,
                        -8,
                        -2,
                    )
                })
                .max()
                .expect("non-seed graph has a source-to-sink path");
            prop_assert_eq!(alignment.score(), oracle);
        }
    }
}
