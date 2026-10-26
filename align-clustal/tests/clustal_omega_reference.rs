//! Inert black-box comparison against an isolated Clustal Omega 1.2.4 executable.
//!
//! The executable and fixture generator live in a sibling validation laboratory.
//! This test only decodes reviewed JSON and never starts Python or Clustal Omega.

use std::collections::{BTreeMap, BTreeSet};

use align_clustal::{ClustalAlphabet, ClustalOmega, OmegaGuideTreeMode};
use align_core::{Alphabet, Dna, Protein, Rna, Sequence, Symbol};
use align_msa::{AlignedRow, GuideTreeNode, MultipleAlignment, SequenceRecord, compare_alignments};
use serde::Deserialize;

const FIXTURE: &str = include_str!("fixtures/clustalo-1.2.4/default.json");
const PRINTED_DISTANCE_TOLERANCE: f64 = 5.1e-7;
const PRINTED_BRANCH_TOLERANCE: f64 = 1.1e-6;

#[derive(Debug, Deserialize)]
struct Fixture {
    schema: u32,
    reference: Provenance,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Provenance {
    name: String,
    version: String,
    ubuntu_package: String,
    executable_sha256: String,
    generated_utc: String,
}

#[derive(Debug, Deserialize)]
struct Case {
    name: String,
    molecule: String,
    records: Vec<(String, String)>,
    alignment: Vec<(String, String)>,
    distance_matrix: Vec<Vec<f64>>,
    guide_tree_newick: String,
}

#[derive(Debug)]
struct ParsedNode {
    name: Option<String>,
    branch_length: Option<f64>,
    children: Vec<ParsedNode>,
}

#[derive(Clone, Copy, Debug)]
struct ObservedMetrics {
    sp: f64,
    tc: f64,
    exact: bool,
}

struct NewickParser<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> NewickParser<'a> {
    fn parse(text: &'a str) -> ParsedNode {
        let mut parser = Self {
            bytes: text.as_bytes(),
            position: 0,
        };
        let root = parser.node();
        parser.consume(b';');
        assert_eq!(parser.position, parser.bytes.len(), "trailing Newick input");
        root
    }

    fn node(&mut self) -> ParsedNode {
        let (name, children) = if self.peek() == Some(b'(') {
            self.position += 1;
            let mut children = vec![self.node()];
            while self.peek() == Some(b',') {
                self.position += 1;
                children.push(self.node());
            }
            self.consume(b')');
            (self.optional_name(), children)
        } else {
            (Some(self.required_name()), Vec::new())
        };
        let branch_length = if self.peek() == Some(b':') {
            self.position += 1;
            let start = self.position;
            while self
                .peek()
                .is_some_and(|byte| !matches!(byte, b',' | b')' | b';'))
            {
                self.position += 1;
            }
            Some(
                std::str::from_utf8(&self.bytes[start..self.position])
                    .unwrap()
                    .parse::<f64>()
                    .unwrap(),
            )
        } else {
            None
        };
        ParsedNode {
            name,
            branch_length,
            children,
        }
    }

    fn optional_name(&mut self) -> Option<String> {
        (!matches!(self.peek(), Some(b':' | b',' | b')' | b';') | None))
            .then(|| self.required_name())
    }

    fn required_name(&mut self) -> String {
        let start = self.position;
        while self
            .peek()
            .is_some_and(|byte| !matches!(byte, b':' | b',' | b')' | b';'))
        {
            self.position += 1;
        }
        assert!(self.position > start, "Newick leaf name is empty");
        std::str::from_utf8(&self.bytes[start..self.position])
            .unwrap()
            .to_owned()
    }

    fn consume(&mut self, expected: u8) {
        assert_eq!(self.peek(), Some(expected), "unexpected Newick token");
        self.position += 1;
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.position).copied()
    }
}

fn reference_alignment<A: Alphabet>(case: &Case) -> MultipleAlignment<A> {
    let rows = case
        .alignment
        .iter()
        .map(|(identifier, aligned)| {
            let symbols = aligned
                .bytes()
                .map(|byte| {
                    if byte == b'-' {
                        None
                    } else {
                        Some(Symbol::<A>::from_ascii(byte).unwrap_or_else(|| {
                            panic!(
                                "case {}: invalid {} symbol {byte:?} in {identifier}",
                                case.name,
                                A::NAME
                            )
                        }))
                    }
                })
                .collect();
            AlignedRow::new(identifier, symbols).unwrap()
        })
        .collect();
    MultipleAlignment::try_new(rows).unwrap()
}

fn records<A: Alphabet>(case: &Case) -> Vec<SequenceRecord<A>> {
    case.records
        .iter()
        .map(|(identifier, sequence)| {
            SequenceRecord::new(
                identifier,
                Sequence::<A>::from_ascii(sequence.as_bytes()).unwrap(),
            )
            .unwrap()
        })
        .collect()
}

fn clade_key(members: &BTreeSet<String>) -> String {
    members.iter().cloned().collect::<Vec<_>>().join("\u{1f}")
}

fn reference_edges(
    node: &ParsedNode,
    edges: &mut BTreeMap<String, f64>,
    is_root: bool,
) -> BTreeSet<String> {
    let members = if node.children.is_empty() {
        BTreeSet::from([node.name.clone().expect("Newick leaf has a name")])
    } else {
        node.children
            .iter()
            .flat_map(|child| reference_edges(child, edges, false))
            .collect()
    };
    if !is_root {
        edges.insert(
            clade_key(&members),
            node.branch_length
                .expect("non-root Newick node has a length"),
        );
    }
    members
}

fn observed_edges(
    node: &GuideTreeNode,
    identifiers: &[String],
    edges: &mut BTreeMap<String, f64>,
) -> BTreeSet<String> {
    match node {
        GuideTreeNode::Leaf { sequence_index } => {
            BTreeSet::from([identifiers[*sequence_index].clone()])
        }
        GuideTreeNode::Branch { left, right } => {
            let left_members = observed_edges(left.child(), identifiers, edges);
            let right_members = observed_edges(right.child(), identifiers, edges);
            edges.insert(clade_key(&left_members), left.length());
            edges.insert(clade_key(&right_members), right.length());
            left_members.union(&right_members).cloned().collect()
        }
    }
}

fn assert_guide_tree<A: ClustalAlphabet>(
    case: &Case,
    result: &align_clustal::ClustalOmegaAlignment<A>,
) {
    let parsed = NewickParser::parse(&case.guide_tree_newick);
    let mut expected = BTreeMap::new();
    let expected_members = reference_edges(&parsed, &mut expected, true);
    let identifiers = case
        .records
        .iter()
        .map(|(identifier, _)| identifier.clone())
        .collect::<Vec<_>>();
    let mut observed = BTreeMap::new();
    let observed_members = observed_edges(result.guide_tree().root(), &identifiers, &mut observed);
    assert_eq!(
        observed_members, expected_members,
        "case {} leaves",
        case.name
    );
    assert_eq!(
        observed.keys().collect::<Vec<_>>(),
        expected.keys().collect::<Vec<_>>(),
        "case {} clades",
        case.name
    );
    for (clade, expected_length) in expected {
        let observed_length = observed[&clade];
        assert!(
            (observed_length - expected_length).abs() <= PRINTED_BRANCH_TOLERANCE,
            "case {} clade {clade:?}: branch {observed_length} versus {expected_length}",
            case.name
        );
    }
}

fn run_case<A: ClustalAlphabet>(case: &Case) -> ObservedMetrics {
    let records = records::<A>(case);
    let reference = reference_alignment::<A>(case);
    let parameters = ClustalOmega::<A>::published()
        .parameters()
        .clone()
        .with_guide_tree_mode(OmegaGuideTreeMode::Full);
    let result = ClustalOmega::<A>::with_parameters(parameters)
        .unwrap()
        .align(&records)
        .unwrap_or_else(|error| panic!("case {} failed: {error}", case.name));

    let distances = result.guide_evidence().distances().unwrap();
    assert_eq!(case.distance_matrix.len(), distances.size());
    for (row, expected_row) in case.distance_matrix.iter().enumerate() {
        assert_eq!(expected_row.len(), distances.size());
        for (column, expected) in expected_row.iter().copied().enumerate() {
            let observed = distances.get(row, column).unwrap();
            assert!(
                (observed - expected).abs() <= PRINTED_DISTANCE_TOLERANCE,
                "case {} distance ({row},{column}): {observed} versus {expected}",
                case.name
            );
        }
    }
    assert_guide_tree(case, &result);

    let metrics = compare_alignments(result.alignment(), &reference).unwrap();
    if metrics.sp_score() < 0.8 {
        eprintln!("{} alignment disagreement:", case.name);
        eprintln!("steps: {:?}", result.progressive_steps());
        for row in result.alignment().rows() {
            eprintln!(
                "Rust      {:<12} {}",
                row.identifier(),
                row.to_ascii_with_gaps()
            );
        }
        for row in reference.rows() {
            eprintln!(
                "clustalo  {:<12} {}",
                row.identifier(),
                row.to_ascii_with_gaps()
            );
        }
    }
    ObservedMetrics {
        sp: metrics.sp_score(),
        tc: metrics.tc_score(),
        exact: result.alignment() == &reference,
    }
}

#[test]
fn published_pipeline_matches_static_clustalo_observations() {
    let fixture: Fixture = serde_json::from_str(FIXTURE).expect("fixture JSON must be valid");
    assert_eq!(fixture.schema, 1);
    assert_eq!(fixture.reference.name, "Clustal Omega");
    assert_eq!(fixture.reference.version, "1.2.4");
    assert_eq!(fixture.reference.ubuntu_package, "1.2.4-8build1");
    assert_eq!(fixture.reference.generated_utc, "2026-08-03");
    assert_eq!(
        fixture.reference.executable_sha256,
        "e695ae9d1fdc9b27540eaf041740b3dafc6e38670cbadee8afe3f6d9114062ea"
    );
    assert_eq!(fixture.cases.len(), 5);

    let mut observed = BTreeMap::new();
    for case in &fixture.cases {
        let metrics = match case.molecule.as_str() {
            "protein" => run_case::<Protein>(case),
            "dna" => run_case::<Dna>(case),
            "rna" => run_case::<Rna>(case),
            molecule => panic!("case {}: unknown molecule {molecule}", case.name),
        };
        observed.insert(case.name.as_str(), metrics);
    }

    for name in [
        "protein_internal_gap",
        "protein_terminal_and_ambiguity",
        "dna_internal_and_terminal_gaps",
        "rna_internal_gap",
    ] {
        let metrics = observed[name];
        assert!(metrics.exact, "case {name} must agree exactly");
        assert_eq!((metrics.sp, metrics.tc), (1.0, 1.0));
    }

    let repetitive = observed["dna_repetitive"];
    assert!(
        repetitive.sp >= 0.35,
        "repetitive-case SP {} < 0.35",
        repetitive.sp
    );

    let count = f64::from(u32::try_from(observed.len()).expect("case count fits u32"));
    let mean_sp = observed.values().map(|metrics| metrics.sp).sum::<f64>() / count;
    let mean_tc = observed.values().map(|metrics| metrics.tc).sum::<f64>() / count;
    assert!(mean_sp >= 0.87, "mean SP {mean_sp} is below 0.87");
    assert!(mean_tc >= 0.80, "mean TC {mean_tc} is below 0.80");
}
