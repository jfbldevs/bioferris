use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use align_core::Alphabet;
use align_poa::PoaGraph;

use crate::FormatError;

/// Writes a POASTA graph in GFA 1.0 segment/link/path form.
///
/// # Errors
///
/// Returns [`FormatError::Io`] for failed output.
pub fn write_poa_gfa<A: Alphabet>(
    path: impl AsRef<Path>,
    graph: &PoaGraph<A>,
) -> Result<(), FormatError> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_poa_gfa_to(&mut writer, graph)?;
    writer.flush()?;
    Ok(())
}

/// Writes a POASTA graph in GFA 1.0 segment/link/path form to a stream.
///
/// # Errors
///
/// Returns [`FormatError::Io`] for failed output.
pub fn write_poa_gfa_to<A: Alphabet>(
    mut writer: impl Write,
    graph: &PoaGraph<A>,
) -> Result<(), FormatError> {
    writer.write_all(b"H\tVN:Z:1.0\n")?;
    for node in graph.nodes() {
        writeln!(
            writer,
            "S\tn{}\t{}\tCL:i:{}",
            node.id().index(),
            char::from(node.symbol().as_ascii()),
            node.class().index()
        )?;
    }
    for edge in graph.edges().iter().copied() {
        writeln!(
            writer,
            "L\tn{}\t+\tn{}\t+\t0M\tRC:i:{}",
            edge.source().index(),
            edge.target().index(),
            edge.support()
        )?;
    }
    for path in graph.paths() {
        write!(writer, "P\t{}\t", gfa_name(path.identifier()))?;
        if path.nodes().is_empty() {
            writer.write_all(b"*")?;
        } else {
            for (index, node) in path.nodes().iter().copied().enumerate() {
                if index != 0 {
                    writer.write_all(b",")?;
                }
                write!(writer, "n{}+", node.index())?;
            }
        }
        writer.write_all(b"\t*\n")?;
    }
    Ok(())
}

/// Writes a POASTA graph as deterministic Graphviz DOT.
///
/// # Errors
///
/// Returns [`FormatError::Io`] for failed output.
pub fn write_poa_dot<A: Alphabet>(
    path: impl AsRef<Path>,
    graph: &PoaGraph<A>,
) -> Result<(), FormatError> {
    let mut writer = BufWriter::new(File::create(path)?);
    write_poa_dot_to(&mut writer, graph)?;
    writer.flush()?;
    Ok(())
}

/// Writes a POASTA graph as deterministic Graphviz DOT to a stream.
///
/// # Errors
///
/// Returns [`FormatError::Io`] for failed output.
pub fn write_poa_dot_to<A: Alphabet>(
    mut writer: impl Write,
    graph: &PoaGraph<A>,
) -> Result<(), FormatError> {
    writer.write_all(b"digraph POASTA {\n  rankdir=LR;\n")?;
    for node in graph.nodes() {
        writeln!(
            writer,
            "  n{} [label=\"n{}\\n{}\\nclass {}\"];",
            node.id().index(),
            node.id().index(),
            char::from(node.symbol().as_ascii()),
            node.class().index()
        )?;
    }
    for edge in graph.edges().iter().copied() {
        writeln!(
            writer,
            "  n{} -> n{} [label=\"{}\"];",
            edge.source().index(),
            edge.target().index(),
            edge.support()
        )?;
    }
    writer.write_all(b"}\n")?;
    Ok(())
}

fn gfa_name(identifier: &str) -> String {
    let mut output = String::with_capacity(identifier.len());
    for byte in identifier.bytes() {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-') {
            output.push(char::from(byte));
        } else {
            use core::fmt::Write as _;
            write!(output, "%{byte:02X}").expect("writing to String cannot fail");
        }
    }
    output
}

/// Convenient graph output methods.
pub trait PoaGraphWrite {
    /// Writes this graph in GFA 1.0 format.
    ///
    /// # Errors
    ///
    /// Returns [`FormatError::Io`] for failed output.
    fn write_gfa(&self, path: impl AsRef<Path>) -> Result<(), FormatError>;

    /// Writes this graph in Graphviz DOT format.
    ///
    /// # Errors
    ///
    /// Returns [`FormatError::Io`] for failed output.
    fn write_dot(&self, path: impl AsRef<Path>) -> Result<(), FormatError>;
}

impl<A: Alphabet> PoaGraphWrite for PoaGraph<A> {
    fn write_gfa(&self, path: impl AsRef<Path>) -> Result<(), FormatError> {
        write_poa_gfa(path, self)
    }

    fn write_dot(&self, path: impl AsRef<Path>) -> Result<(), FormatError> {
        write_poa_dot(path, self)
    }
}
