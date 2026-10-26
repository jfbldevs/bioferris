use core::fmt;

use align_core::{Alphabet, Sequence};
use align_poa::{PoaGraph, PoaNodeId, PoaOperation, PoastaAffineCosts, PoastaAlignment};

use super::{ColorMode, CoordinateStyle, ReportError, ReportOptions, Verbosity, sanitized_label};

/// One exact POASTA query-to-graph result plus its graph, query, and costs.
#[derive(Clone, Debug)]
pub struct PoastaQueryRun<'a, A: Alphabet> {
    graph: &'a PoaGraph<A>,
    query: &'a Sequence<A>,
    result: PoastaAlignment,
    costs: PoastaAffineCosts,
    profile: &'static str,
}

impl<'a, A: Alphabet> PoastaQueryRun<'a, A> {
    pub(crate) const fn new(
        graph: &'a PoaGraph<A>,
        query: &'a Sequence<A>,
        result: PoastaAlignment,
        costs: PoastaAffineCosts,
        profile: &'static str,
    ) -> Self {
        Self {
            graph,
            query,
            result,
            costs,
            profile,
        }
    }

    /// Returns the authoritative input graph.
    #[must_use]
    pub const fn graph(&self) -> &'a PoaGraph<A> {
        self.graph
    }

    /// Returns the query sequence.
    #[must_use]
    pub const fn query(&self) -> &'a Sequence<A> {
        self.query
    }

    /// Returns the untouched native POASTA result.
    #[must_use]
    pub const fn result(&self) -> &PoastaAlignment {
        &self.result
    }

    /// Returns the exact immutable POASTA profile identifier.
    #[must_use]
    pub const fn profile(&self) -> &'static str {
        self.profile
    }

    /// Returns the positive affine costs used by the query alignment.
    #[must_use]
    pub const fn costs(&self) -> PoastaAffineCosts {
        self.costs
    }

    /// Consumes this run and returns the native result.
    #[must_use]
    pub fn into_result(self) -> PoastaAlignment {
        self.result
    }

    /// Creates a normal query-to-graph report.
    #[must_use]
    pub fn report(&self) -> PoastaQueryReport<'_, 'a, A> {
        PoastaQueryReport {
            run: self,
            options: ReportOptions::new(),
            graph_label: String::from("Graph"),
            query_label: String::from("Query"),
        }
    }
}

/// Borrowed configurable POASTA sequence-to-graph report.
#[derive(Clone, Debug)]
pub struct PoastaQueryReport<'run, 'input, A: Alphabet> {
    run: &'run PoastaQueryRun<'input, A>,
    options: ReportOptions,
    graph_label: String,
    query_label: String,
}

impl<A: Alphabet> PoastaQueryReport<'_, '_, A> {
    /// Replaces graph and query row labels.
    #[must_use]
    pub fn labels(mut self, graph: impl Into<String>, query: impl Into<String>) -> Self {
        self.graph_label = sanitized_label(graph);
        self.query_label = sanitized_label(query);
        self
    }

    /// Replaces the positive traceback block width.
    ///
    /// # Errors
    ///
    /// Returns [`ReportError::InvalidBlockWidth`] for zero.
    pub fn block_width(mut self, width: usize) -> Result<Self, ReportError> {
        self.options = self.options.with_block_width(width)?;
        Ok(self)
    }

    /// Replaces report verbosity.
    #[must_use]
    pub const fn verbosity(mut self, verbosity: Verbosity) -> Self {
        self.options = self.options.with_verbosity(verbosity);
        self
    }

    /// Replaces query-coordinate representation.
    #[must_use]
    pub const fn coordinates(mut self, coordinates: CoordinateStyle) -> Self {
        self.options = self.options.with_coordinates(coordinates);
        self
    }

    /// Retains a color preference for terminal-aware writers.
    ///
    /// Stable [`fmt::Display`] output remains colorless.
    #[must_use]
    pub const fn color(mut self, color: ColorMode) -> Self {
        self.options = self.options.with_color(color);
        self
    }

    /// Returns the complete rendering options.
    #[must_use]
    pub const fn options(&self) -> ReportOptions {
        self.options
    }
}

impl<A: Alphabet> fmt::Display for PoastaQueryReport<'_, '_, A> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write_query_summary(formatter, self.run)?;

        if self.options.verbosity() == Verbosity::Compact {
            return Ok(());
        }
        if self.options.verbosity() == Verbosity::Detailed {
            writeln!(formatter)?;
            write_query_details(formatter, self.run.result())?;
        } else {
            writeln!(formatter)?;
        }
        formatter.write_str("\n")?;
        write_graph_blocks(
            formatter,
            self.run.graph,
            self.run.query,
            self.run.result(),
            &self.graph_label,
            &self.query_label,
            self.options,
        )?;
        formatter.write_str(
            "\nLegend: | exact match, . substitution, blank graph/query gap; nX is a stable graph node",
        )
    }
}

const QUERY_FIELD_WIDTH: usize = 20;

fn write_query_summary<A: Alphabet>(
    formatter: &mut fmt::Formatter<'_>,
    run: &PoastaQueryRun<'_, A>,
) -> fmt::Result {
    let result = run.result();
    let width = QUERY_FIELD_WIDTH;
    writeln!(formatter, "{:<width$} : POASTA", "Algorithm")?;
    writeln!(formatter, "{:<width$} : {}", "Profile", run.profile)?;
    writeln!(
        formatter,
        "{:<width$} : global affine sequence-to-graph",
        "Mode"
    )?;
    writeln!(
        formatter,
        "{:<width$} : mismatch={} gap_open={} gap_extend={} (positive costs)",
        "Scoring",
        run.costs.mismatch(),
        run.costs.gap_opening(),
        run.costs.gap_extension()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Score (negated cost)",
        result.score()
    )?;
    writeln!(formatter, "{:<width$} : {:?}", "Backend", result.backend())?;
    writeln!(
        formatter,
        "{:<width$} : nodes={} edges={} classes={} paths={}",
        "Graph",
        run.graph.node_count(),
        run.graph.edge_count(),
        run.graph.class_count(),
        run.graph.path_count()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Query length",
        run.query.len()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Trace columns",
        result.operations().len()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Aligned pairs",
        result.aligned_pairs().len()
    )?;
    write!(formatter, "{:<width$} : ", "Trace")?;
    write_trace(formatter, result.operations())
}

fn write_query_details(
    formatter: &mut fmt::Formatter<'_>,
    result: &PoastaAlignment,
) -> fmt::Result {
    let statistics = result.statistics();
    let memory = result.memory_estimate();
    let options = result.options();
    let width = QUERY_FIELD_WIDTH;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Possible states",
        statistics.total_possible_states()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Discovered states",
        statistics.discovered_states()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Expanded states",
        statistics.expanded_states()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "A* expansions",
        statistics.astar_expanded_states()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "DFS expansions",
        statistics.dfs_expanded_states()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Greedy matches",
        statistics.greedy_exact_matches()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Bubble prunes",
        statistics.bubble_pruned_states()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Peak queue entries",
        statistics.peak_queue_entries()
    )?;
    writeln!(
        formatter,
        "{:<width$} : {}",
        "Modeled peak bytes",
        memory.modeled_peak_bytes()
    )?;
    writeln!(
        formatter,
        "{:<width$} : exact_match_dfs={} superbubble_pruning={}",
        "Accelerations",
        options.exact_match_dfs(),
        options.superbubble_pruning()
    )
}

fn write_trace(formatter: &mut fmt::Formatter<'_>, operations: &[PoaOperation]) -> fmt::Result {
    let Some(first) = operations.first().copied() else {
        return formatter.write_str("-");
    };
    let mut selected = poa_operation_code(first);
    let mut count = 0_usize;
    for operation in operations.iter().copied() {
        let code = poa_operation_code(operation);
        if code == selected {
            count += 1;
        } else {
            write!(formatter, "{count}{selected}")?;
            selected = code;
            count = 1;
        }
    }
    write!(formatter, "{count}{selected}")
}

const fn poa_operation_code(operation: PoaOperation) -> char {
    match operation {
        PoaOperation::Aligned { .. } => 'A',
        PoaOperation::GraphOnly { .. } => 'G',
        PoaOperation::QueryOnly { .. } => 'Q',
    }
}

#[allow(clippy::too_many_arguments)]
fn write_graph_blocks<A: Alphabet>(
    formatter: &mut fmt::Formatter<'_>,
    graph: &PoaGraph<A>,
    query: &Sequence<A>,
    result: &PoastaAlignment,
    graph_label: &str,
    query_label: &str,
    options: ReportOptions,
) -> fmt::Result {
    if result.operations().is_empty() {
        return formatter.write_str("(empty alignment)");
    }
    let label_width = graph_label.chars().count().max(query_label.chars().count());
    let coordinate_width = format!("n{}", graph.node_count().max(query.len())).len();
    let mut query_position = 0_usize;

    for (block_index, operations) in result
        .operations()
        .chunks(options.block_width())
        .enumerate()
    {
        if block_index != 0 {
            formatter.write_str("\n\n")?;
        }
        let query_before = query_position;
        let mut first_node = None;
        let mut last_node = None;
        let mut graph_row = String::with_capacity(operations.len());
        let mut query_row = String::with_capacity(operations.len());
        let mut middle = String::with_capacity(operations.len());
        for operation in operations.iter().copied() {
            match operation {
                PoaOperation::Aligned { node, query_index } => {
                    let graph_symbol = graph
                        .node(node)
                        .expect("validated POASTA node is present")
                        .symbol();
                    let query_symbol = query
                        .get(query_index)
                        .expect("validated POASTA query index is present");
                    remember_node(&mut first_node, &mut last_node, node);
                    graph_row.push(char::from(graph_symbol.as_ascii()));
                    query_row.push(char::from(query_symbol.as_ascii()));
                    middle.push(if graph_symbol == query_symbol {
                        '|'
                    } else {
                        '.'
                    });
                    query_position += 1;
                }
                PoaOperation::GraphOnly { node } => {
                    let graph_symbol = graph
                        .node(node)
                        .expect("validated POASTA node is present")
                        .symbol();
                    remember_node(&mut first_node, &mut last_node, node);
                    graph_row.push(char::from(graph_symbol.as_ascii()));
                    query_row.push('-');
                    middle.push(' ');
                }
                PoaOperation::QueryOnly { query_index } => {
                    let query_symbol = query
                        .get(query_index)
                        .expect("validated POASTA query index is present");
                    graph_row.push('-');
                    query_row.push(char::from(query_symbol.as_ascii()));
                    middle.push(' ');
                    query_position += 1;
                }
            }
        }
        let graph_start = first_node.map_or_else(|| String::from("-"), node_label);
        let graph_end = last_node.map_or_else(|| String::from("-"), node_label);
        writeln!(
            formatter,
            "{graph_label:<label_width$}  {graph_start:>coordinate_width$}  {graph_row}  {graph_end:>coordinate_width$}"
        )?;
        writeln!(
            formatter,
            "{:<label_width$}  {:>coordinate_width$}  {middle}",
            "", ""
        )?;
        let (query_start, query_end) =
            query_coordinates(query_before, query_position, options.coordinates());
        writeln!(
            formatter,
            "{query_label:<label_width$}  {query_start:>coordinate_width$}  {query_row}  {query_end:>coordinate_width$}"
        )?;
    }
    Ok(())
}

fn remember_node(first: &mut Option<PoaNodeId>, last: &mut Option<PoaNodeId>, node: PoaNodeId) {
    first.get_or_insert(node);
    *last = Some(node);
}

fn node_label(node: PoaNodeId) -> String {
    format!("n{}", node.index())
}

fn query_coordinates(before: usize, after: usize, style: CoordinateStyle) -> (String, String) {
    if before == after {
        return (String::from("-"), String::from("-"));
    }
    match style {
        CoordinateStyle::Human => ((before + 1).to_string(), after.to_string()),
        CoordinateStyle::Rust => (before.to_string(), after.to_string()),
    }
}
