//! Validated import of the residue-level GFA and DOT dialects emitted by `RustBio`.
//!
//! These immutable documents are not computational `PoaGraph` values. GFA keeps
//! paths, classes and support; DOT contains topology only and cannot recover paths.
//!
//! ```
//! use bioferris::{Dna, io::graph::parse_gfa};
//! let input = "H\tVN:Z:1.0\nS\tn0\tA\tCL:i:0\nS\tn1\tC\tCL:i:1\n\
//!              L\tn0\t+\tn1\t+\t0M\tRC:i:1\nP\tread1\tn0+,n1+\t*\n";
//! let graph = parse_gfa::<Dna>(input)?;
//! assert_eq!(graph.nodes().len(), 2);
//! assert_eq!(graph.path_sequence(0).as_deref(), Some("AC"));
//! println!("{} nodes; {} path: {}", graph.nodes().len(), graph.paths().len(),
//!          graph.path_sequence(0).unwrap_or_default());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::io::BufRead;
use std::path::Path;

use align_core::{Alphabet, Symbol};

use super::{BioIoError, IoErrorKind, open_input, read_text};

const MAX_DOCUMENT_BYTES: usize = 64 * 1024 * 1024;

/// One validated residue node, with its original stable numeric identifier.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GfaNode<A: Alphabet> {
    id: usize,
    symbol: Symbol<A>,
    class: usize,
}

impl<A: Alphabet> GfaNode<A> {
    /// Returns the integer in the segment name `n<integer>`.
    #[must_use]
    pub const fn id(&self) -> usize {
        self.id
    }

    /// Returns the normalized typed residue.
    #[must_use]
    pub const fn symbol(&self) -> Symbol<A> {
        self.symbol
    }

    /// Returns the original alignment-class identifier.
    #[must_use]
    pub const fn class(&self) -> usize {
        self.class
    }
}

/// One forward, zero-overlap link and its stated path support.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GfaEdge {
    source: usize,
    target: usize,
    support: usize,
}

impl GfaEdge {
    /// Returns the zero-based source node identifier.
    #[must_use]
    pub const fn source(&self) -> usize {
        self.source
    }

    /// Returns the zero-based target node identifier.
    #[must_use]
    pub const fn target(&self) -> usize {
        self.target
    }

    /// Returns the positive path count from the `RC:i` tag.
    #[must_use]
    pub const fn support(&self) -> usize {
        self.support
    }
}

/// One GFA path, retaining its file-level name and ordered node identifiers.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GfaPath {
    name: String,
    nodes: Vec<usize>,
}

impl GfaPath {
    /// Returns the exact GFA name; percent sequences are not implicitly decoded.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns node identifiers in source-residue order.
    #[must_use]
    pub fn nodes(&self) -> &[usize] {
        &self.nodes
    }
}

/// A complete, internally consistent import of the project's GFA subset.
///
/// Node and class identifiers, edge support and input path order are preserved.
/// There is deliberately no conversion to the computational graph: its current
/// public API cannot reconstruct these components without changing their identity.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GfaDocument<A: Alphabet> {
    topology: DotTopology<A>,
    paths: Vec<GfaPath>,
}

impl<A: Alphabet> GfaDocument<A> {
    /// Returns nodes in their original stable identifier order.
    #[must_use]
    pub fn nodes(&self) -> &[GfaNode<A>] {
        self.topology.nodes()
    }

    /// Returns links in source/target order.
    #[must_use]
    pub fn edges(&self) -> &[GfaEdge] {
        self.topology.edges()
    }

    /// Returns paths in their file order.
    #[must_use]
    pub fn paths(&self) -> &[GfaPath] {
        &self.paths
    }

    /// Returns the number of contiguous alignment classes.
    #[must_use]
    pub const fn class_count(&self) -> usize {
        self.topology.class_count
    }

    /// Spells a path's normalized sequence, or returns `None` for an absent path.
    #[must_use]
    pub fn path_sequence(&self, index: usize) -> Option<String> {
        self.paths.get(index).map(|path| {
            path.nodes
                .iter()
                .map(|node| char::from(self.topology.nodes[*node].symbol.as_ascii()))
                .collect()
        })
    }
}

/// Validated topology from the project's restricted DOT dialect.
///
/// Source paths, occurrences and sequence identifiers are absent from DOT. Edge
/// support is positive but cannot be checked against paths in this representation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DotTopology<A: Alphabet> {
    nodes: Vec<GfaNode<A>>,
    edges: Vec<GfaEdge>,
    class_count: usize,
}

impl<A: Alphabet> DotTopology<A> {
    /// Returns the residue nodes in stable identifier order.
    #[must_use]
    pub fn nodes(&self) -> &[GfaNode<A>] {
        &self.nodes
    }

    /// Returns the ordered directed edges and their stated support.
    #[must_use]
    pub fn edges(&self) -> &[GfaEdge] {
        &self.edges
    }

    /// Returns the number of alignment classes.
    #[must_use]
    pub const fn class_count(&self) -> usize {
        self.class_count
    }
}

/// Reads a GFA document from a path, including gzip input.
///
/// # Errors
/// Returns a bounded-input, I/O, unsupported-format or invariant error.
pub fn read_gfa<A: Alphabet>(path: impl AsRef<Path>) -> Result<GfaDocument<A>, BioIoError> {
    parse_gfa_reader(open_input(path)?)
}

/// Reads a GFA document from a buffered input stream.
///
/// # Errors
/// Returns an I/O, size, unsupported-format or invariant error.
pub fn parse_gfa_reader<A: Alphabet>(reader: impl BufRead) -> Result<GfaDocument<A>, BioIoError> {
    parse_gfa(&read_text(
        reader.take(MAX_DOCUMENT_BYTES as u64 + 1),
        "GFA",
    )?)
}

/// Parses the explicit GFA 1.0 `H/S/L/P` residue-graph subset.
///
/// # Errors
/// Rejects missing/extra tags, arbitrary segment names, multi-residue segments,
/// reverse orientations, nonzero overlaps, invalid paths/classes/support, cycles,
/// and documents larger than 64 MiB. Comments and blank lines are ignored.
#[allow(clippy::too_many_lines)]
pub fn parse_gfa<A: Alphabet>(input: &str) -> Result<GfaDocument<A>, BioIoError> {
    document_limit(input, "GFA")?;
    let mut nodes = Vec::new();
    let mut edges = Vec::new();
    let mut paths = Vec::new();
    let mut stage = 0;
    for (index, line) in input.lines().enumerate() {
        let line_number = index + 1;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parsed = (|| {
            if !line.is_ascii() {
                return Err(invalid("GFA", "GFA input must be ASCII"));
            }
            let fields: Vec<_> = line.split('\t').collect();
            match fields.first().copied().unwrap_or_default() {
                "H" => {
                    if stage != 0 {
                        return Err(invalid("GFA", "header must occur once before records"));
                    }
                    if fields != ["H", "VN:Z:1.0"] {
                        return Err(unsupported(
                            "GFA",
                            "only header H with VN:Z:1.0 is supported",
                        ));
                    }
                    stage = 1;
                }
                "S" => {
                    if stage != 1 {
                        return Err(invalid(
                            "GFA",
                            "segments must follow the header and precede links/paths",
                        ));
                    }
                    if fields.len() != 4 {
                        return Err(unsupported(
                            "GFA",
                            "segments require exactly name, residue and CL:i tag",
                        ));
                    }
                    let node = parse_node::<A>(
                        fields[1],
                        fields[2],
                        fields[3]
                            .strip_prefix("CL:i:")
                            .ok_or_else(|| unsupported("GFA", "segment requires CL:i class tag"))?,
                        "GFA",
                    )?;
                    if node.id != nodes.len() {
                        return Err(invalid(
                            "GFA",
                            "node names must be contiguous n0, n1, ... in order",
                        ));
                    }
                    nodes.push(node);
                }
                "L" => {
                    if !(1..=2).contains(&stage) {
                        return Err(invalid(
                            "GFA",
                            "links must follow segments and precede paths",
                        ));
                    }
                    if fields.len() != 7 {
                        return Err(unsupported(
                            "GFA",
                            "links require exactly two orientations, overlap and RC:i tag",
                        ));
                    }
                    if fields[2] != "+" || fields[4] != "+" || fields[5] != "0M" {
                        return Err(unsupported(
                            "GFA",
                            "only forward (+/+) zero-overlap (0M) links are supported",
                        ));
                    }
                    edges.push(GfaEdge {
                        source: node_id(fields[1], "GFA")?,
                        target: node_id(fields[3], "GFA")?,
                        support: natural(
                            fields[6].strip_prefix("RC:i:").ok_or_else(|| {
                                unsupported("GFA", "link requires RC:i support tag")
                            })?,
                            "GFA",
                        )?,
                    });
                    stage = 2;
                }
                "P" => {
                    if stage == 0 {
                        return Err(invalid("GFA", "a version header is required before paths"));
                    }
                    if fields.len() != 4 || fields[3] != "*" {
                        return Err(unsupported(
                            "GFA",
                            "paths require a name, forward node list, and * overlap field",
                        ));
                    }
                    let name = fields[1];
                    if name.is_empty()
                        || name.starts_with(['*', '='])
                        || name.bytes().any(|byte| !byte.is_ascii_graphic())
                        || name.contains("+,")
                        || name.contains("-,")
                    {
                        return Err(invalid("GFA", "invalid GFA path name"));
                    }
                    let path_nodes = if fields[2] == "*" {
                        Vec::new()
                    } else {
                        fields[2]
                            .split(',')
                            .map(|item| {
                                node_id(
                                    item.strip_suffix('+').ok_or_else(|| {
                                        unsupported("GFA", "path steps must use + orientation")
                                    })?,
                                    "GFA",
                                )
                            })
                            .collect::<Result<Vec<_>, _>>()?
                    };
                    paths.push(GfaPath {
                        name: name.to_owned(),
                        nodes: path_nodes,
                    });
                    stage = 3;
                }
                _ => {
                    return Err(unsupported(
                        "GFA",
                        format!("unsupported record type {:?}", fields[0]),
                    ));
                }
            }
            Ok(())
        })();
        parsed.map_err(|error| error.at(1, line_number))?;
    }
    if stage == 0 {
        return Err(BioIoError::new(
            IoErrorKind::InvalidHeader,
            "GFA",
            "missing GFA 1.0 header",
        ));
    }
    let topology = validate_topology(nodes, edges, "GFA")?;
    validate_paths(&topology, &paths)?;
    Ok(GfaDocument { topology, paths })
}

/// Reads the project's DOT topology dialect from a path, including gzip input.
///
/// # Errors
/// Returns I/O, size, unsupported-syntax or topology-invariant errors.
pub fn read_dot<A: Alphabet>(path: impl AsRef<Path>) -> Result<DotTopology<A>, BioIoError> {
    parse_dot_reader(open_input(path)?)
}

/// Reads the project's DOT topology dialect from a buffered stream.
///
/// # Errors
/// Returns I/O, size, unsupported-syntax or topology-invariant errors.
pub fn parse_dot_reader<A: Alphabet>(reader: impl BufRead) -> Result<DotTopology<A>, BioIoError> {
    parse_dot(&read_text(
        reader.take(MAX_DOCUMENT_BYTES as u64 + 1),
        "DOT",
    )?)
}

/// Parses only the line-oriented DOT syntax produced by `write_poa_dot_to`.
///
/// # Errors
/// Rejects arbitrary Graphviz attributes/statements, missing delimiters, invalid
/// residues/classes/edges, cycles, and documents larger than 64 MiB.
pub fn parse_dot<A: Alphabet>(input: &str) -> Result<DotTopology<A>, BioIoError> {
    document_limit(input, "DOT")?;
    let mut lines = input.lines().enumerate();
    if lines.next().map(|(_, line)| line.trim()) != Some("digraph POASTA {")
        || lines.next().map(|(_, line)| line.trim()) != Some("rankdir=LR;")
    {
        return Err(unsupported(
            "DOT",
            "expected the project's digraph POASTA and rankdir=LR header",
        ));
    }
    let mut nodes = Vec::new();
    let mut edges = Vec::new();
    let mut closed = false;
    for (index, line) in lines {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if closed {
            return Err(invalid("DOT", "unexpected content after closing brace").at(1, index + 1));
        }
        if line == "}" {
            closed = true;
            continue;
        }
        let parsed = (|| {
            let (head, label) = line
                .strip_suffix("\"];")
                .and_then(|line| line.split_once(" [label=\""))
                .ok_or_else(|| unsupported("DOT", "unsupported DOT statement or attributes"))?;
            if let Some((source, target)) = head.split_once(" -> ") {
                edges.push(GfaEdge {
                    source: node_id(source, "DOT")?,
                    target: node_id(target, "DOT")?,
                    support: natural(label, "DOT")?,
                });
            } else {
                if !edges.is_empty() {
                    return Err(invalid("DOT", "nodes must precede edges"));
                }
                let fields: Vec<_> = label.split("\\n").collect();
                if fields.len() != 3 || fields[0] != head {
                    return Err(invalid("DOT", "node label must repeat its node identifier"));
                }
                let class = fields[2]
                    .strip_prefix("class ")
                    .ok_or_else(|| unsupported("DOT", "node label requires class identifier"))?;
                let node = parse_node(head, fields[1], class, "DOT")?;
                if node.id != nodes.len() {
                    return Err(invalid("DOT", "nodes must occur as contiguous n0, n1, ..."));
                }
                nodes.push(node);
            }
            Ok(())
        })();
        parsed.map_err(|error| error.at(1, index + 1))?;
    }
    if !closed {
        return Err(BioIoError::new(
            IoErrorKind::Truncated,
            "DOT",
            "missing closing brace",
        ));
    }
    validate_topology(nodes, edges, "DOT")
}

fn document_limit(input: &str, format: &'static str) -> Result<(), BioIoError> {
    if input.len() > MAX_DOCUMENT_BYTES {
        return Err(BioIoError::new(
            IoErrorKind::LimitExceeded,
            format,
            "graph document exceeds 64 MiB",
        ));
    }
    Ok(())
}

fn natural(input: &str, format: &'static str) -> Result<usize, BioIoError> {
    if input.is_empty()
        || !input.bytes().all(|byte| byte.is_ascii_digit())
        || (input.len() > 1 && input.starts_with('0'))
    {
        return Err(invalid(
            format,
            "expected canonical nonnegative decimal integer",
        ));
    }
    input
        .parse()
        .map_err(|_| invalid(format, "integer exceeds platform range"))
}

fn node_id(input: &str, format: &'static str) -> Result<usize, BioIoError> {
    natural(
        input
            .strip_prefix('n')
            .ok_or_else(|| unsupported(format, "only node names n<integer> are supported"))?,
        format,
    )
}

fn parse_node<A: Alphabet>(
    id: &str,
    residue: &str,
    class: &str,
    format: &'static str,
) -> Result<GfaNode<A>, BioIoError> {
    if residue.len() != 1 {
        return Err(unsupported(
            format,
            "only one-residue segments are supported",
        ));
    }
    let symbol = Symbol::from_ascii(residue.as_bytes()[0]).ok_or_else(|| {
        BioIoError::new(
            IoErrorKind::InvalidSequence,
            format,
            format!("invalid {} node residue", A::NAME),
        )
    })?;
    Ok(GfaNode {
        id: node_id(id, format)?,
        symbol,
        class: natural(class, format)?,
    })
}

fn validate_topology<A: Alphabet>(
    nodes: Vec<GfaNode<A>>,
    edges: Vec<GfaEdge>,
    format: &'static str,
) -> Result<DotTopology<A>, BioIoError> {
    let mut class_symbols = BTreeSet::new();
    let mut classes = BTreeSet::new();
    for node in &nodes {
        if !class_symbols.insert((node.class, node.symbol.as_ascii())) {
            return Err(invalid(format, "duplicate residue in an alignment class"));
        }
        classes.insert(node.class);
    }
    if classes.iter().copied().ne(0..classes.len()) {
        return Err(invalid(
            format,
            "class identifiers must be contiguous from zero",
        ));
    }
    let class_count = classes.len();
    let mut previous = None;
    let mut class_edges = BTreeSet::new();
    for edge in &edges {
        if edge.source >= nodes.len() || edge.target >= nodes.len() {
            return Err(invalid(format, "edge refers to an absent node"));
        }
        if edge.support == 0 {
            return Err(invalid(format, "edge support must be positive"));
        }
        let key = (edge.source, edge.target);
        if previous.is_some_and(|previous| previous >= key) {
            return Err(invalid(
                format,
                "edges must be unique and sorted by source/target",
            ));
        }
        previous = Some(key);
        let source_class = nodes[edge.source].class;
        let target_class = nodes[edge.target].class;
        if source_class == target_class {
            return Err(invalid(
                format,
                "an edge cannot connect nodes in the same class",
            ));
        }
        class_edges.insert((source_class, target_class));
    }
    // Acyclicity of the class quotient also proves acyclicity of the node graph.
    validate_dag(class_count, class_edges, format)?;
    Ok(DotTopology {
        nodes,
        edges,
        class_count,
    })
}

fn validate_dag(
    count: usize,
    edges: BTreeSet<(usize, usize)>,
    format: &'static str,
) -> Result<(), BioIoError> {
    let mut degree = vec![0_usize; count];
    let mut successors = vec![Vec::new(); count];
    for (source, target) in edges {
        degree[target] += 1;
        successors[source].push(target);
    }
    let mut ready: VecDeque<_> = degree
        .iter()
        .enumerate()
        .filter_map(|(id, degree)| (*degree == 0).then_some(id))
        .collect();
    let mut visited = 0;
    while let Some(node) = ready.pop_front() {
        visited += 1;
        for target in &successors[node] {
            degree[*target] -= 1;
            if degree[*target] == 0 {
                ready.push_back(*target);
            }
        }
    }
    if visited != count {
        return Err(invalid(format, "alignment-class graph contains a cycle"));
    }
    Ok(())
}

fn validate_paths<A: Alphabet>(
    topology: &DotTopology<A>,
    paths: &[GfaPath],
) -> Result<(), BioIoError> {
    let mut names: BTreeSet<String> = topology
        .nodes
        .iter()
        .map(|node| format!("n{}", node.id))
        .collect();
    let mut covered = vec![false; topology.nodes.len()];
    let mut support: BTreeMap<_, usize> = topology
        .edges
        .iter()
        .map(|edge| ((edge.source, edge.target), 0))
        .collect();
    for path in paths {
        if !names.insert(path.name.clone()) {
            return Err(invalid(
                "GFA",
                "path name duplicates another path or segment name",
            ));
        }
        let mut classes = BTreeSet::new();
        for id in &path.nodes {
            let node = topology
                .nodes
                .get(*id)
                .ok_or_else(|| invalid("GFA", "path refers to an absent node"))?;
            if !classes.insert(node.class) {
                return Err(invalid(
                    "GFA",
                    "path visits an alignment class more than once",
                ));
            }
            covered[*id] = true;
        }
        for pair in path.nodes.windows(2) {
            let count = support
                .get_mut(&(pair[0], pair[1]))
                .ok_or_else(|| invalid("GFA", "path traverses an absent edge"))?;
            *count = count
                .checked_add(1)
                .ok_or_else(|| invalid("GFA", "path support overflows"))?;
        }
    }
    if covered.contains(&false) {
        return Err(invalid("GFA", "every node must occur in a source path"));
    }
    for edge in &topology.edges {
        if support.get(&(edge.source, edge.target)) != Some(&edge.support) {
            return Err(invalid(
                "GFA",
                "RC support differs from the number of stored paths traversing the edge",
            ));
        }
    }
    Ok(())
}

fn invalid(format: &'static str, message: impl Into<String>) -> BioIoError {
    BioIoError::new(IoErrorKind::InvalidStructure, format, message)
}

fn unsupported(format: &'static str, message: impl Into<String>) -> BioIoError {
    BioIoError::new(IoErrorKind::Unsupported, format, message)
}
