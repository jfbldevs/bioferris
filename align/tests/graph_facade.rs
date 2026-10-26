//! Exact equivalence tests for the POASTA facade.

use align::prelude::*;
use align_poa::{Poasta, PoastaAffineCosts, PoastaPoa};

fn records() -> Vec<SequenceRecord<Dna>> {
    [
        ("reference", "ACGTACGT"),
        ("allele", "ACGTTCGT"),
        ("deletion", "ACGTCGT"),
    ]
    .into_iter()
    .map(|(identifier, sequence)| {
        SequenceRecord::new(identifier, DnaSequence::try_from(sequence).unwrap()).unwrap()
    })
    .collect()
}

#[test]
fn progressive_poasta_facade_matches_the_frozen_kernel() {
    let records = records();
    let costs = PoastaAffineCosts::new(4, 6, 2).unwrap();
    let facade = Graph::poasta().costs(4, 6, 2).align(&records).unwrap();
    let direct = PoastaPoa::new(costs).align(&records).unwrap();
    assert_eq!(facade, direct);
}

#[test]
fn sequence_to_graph_facade_matches_the_frozen_kernel() {
    let records = records();
    let costs = PoastaAffineCosts::new(4, 6, 2).unwrap();
    let progressive = PoastaPoa::new(costs).align(&records).unwrap();
    let query = DnaSequence::try_from("ACGTTACGT").unwrap();

    let facade = Graph::poasta()
        .costs(4, 6, 2)
        .align_query(progressive.graph(), &query)
        .unwrap();
    let direct = Poasta::new(costs)
        .align(progressive.graph(), &query)
        .unwrap();
    assert_eq!(facade, direct);
}
