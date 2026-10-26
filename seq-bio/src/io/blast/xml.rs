//! Namespace-aware XML event transport and complete-query streaming.

use std::collections::HashSet;
use std::io::{self, BufRead, Read};
use std::iter::FusedIterator;
use std::sync::Arc;

use quick_xml::events::{BytesDecl, BytesStart, Event};
use quick_xml::name::ResolveResult;
use quick_xml::reader::NsReader;

use super::{BlastComment, BlastQuery, BlastSourceField, BlastXmlFormat, BlastXmlOptions};
use crate::io::{BioIoError, IoErrorKind};

#[path = "xml_model.rs"]
mod model;

const FORMAT: &str = "BLAST XML";
const NAMESPACE: &str = "http://www.ncbi.nlm.nih.gov";
const MAX_NAMESPACE_DECLARATIONS_PER_ELEMENT: usize = 256;

fn error(kind: IoErrorKind, message: impl Into<String>) -> BioIoError {
    BioIoError::new(kind, FORMAT, message)
}

pub(super) fn relabel(failure: BioIoError) -> BioIoError {
    failure.in_format(FORMAT)
}

struct BoundedInput<R> {
    reader: R,
    total: usize,
    line: usize,
    line_bytes: usize,
    total_limit: usize,
    line_limit: usize,
    event_remaining: usize,
    capture_remaining: Option<usize>,
    limited: bool,
}

impl<R: BufRead> BufRead for BoundedInput<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        let available = self.reader.fill_buf()?;
        if available.is_empty() {
            return Ok(available);
        }
        let remaining = self
            .total_limit
            .saturating_sub(self.total)
            .min(self.event_remaining)
            .min(self.capture_remaining.unwrap_or(usize::MAX));
        if remaining == 0 || self.line_bytes >= self.line_limit {
            self.limited = true;
            return Err(io::Error::other(
                "BLAST XML source/event/capture byte limit exceeded",
            ));
        }
        let count = available
            .len()
            .min(remaining)
            .min(self.line_limit - self.line_bytes);
        let count = available[..count]
            .iter()
            .position(|byte| *byte == b'\n')
            .map_or(count, |index| index + 1);
        Ok(&available[..count])
    }

    fn consume(&mut self, amount: usize) {
        if let Ok(available) = self.reader.fill_buf() {
            for byte in &available[..amount.min(available.len())] {
                if *byte == b'\n' {
                    self.line += 1;
                    self.line_bytes = 0;
                } else {
                    self.line_bytes += 1;
                }
            }
        }
        self.total += amount;
        self.event_remaining = self.event_remaining.saturating_sub(amount);
        if let Some(remaining) = &mut self.capture_remaining {
            *remaining = remaining.saturating_sub(amount);
        }
        self.reader.consume(amount);
    }
}

impl<R: BufRead> Read for BoundedInput<R> {
    fn read(&mut self, destination: &mut [u8]) -> io::Result<usize> {
        if destination.is_empty() {
            return Ok(0);
        }
        let input = self.fill_buf()?;
        let count = input.len().min(destination.len());
        destination[..count].copy_from_slice(&input[..count]);
        self.consume(count);
        Ok(count)
    }
}

pub(super) struct Node {
    pub name: String,
    pub text: String,
    pub attributes: Vec<BlastSourceField>,
    pub children: Vec<Node>,
    pub line: usize,
}

impl Node {
    pub fn one(&self, name: &str) -> Result<Option<&Self>, BioIoError> {
        let mut matching = self.children.iter().filter(|node| node.name == name);
        let first = matching.next();
        if matching.next().is_some() {
            return Err(error(
                IoErrorKind::InvalidStructure,
                format!("repeated XML element {name}"),
            ));
        }
        Ok(first)
    }

    pub fn required(&self, name: &str) -> Result<&Self, BioIoError> {
        self.one(name)?.ok_or_else(|| {
            error(
                IoErrorKind::InvalidStructure,
                format!("missing required XML element {name}"),
            )
        })
    }

    pub fn value(&self) -> Result<&str, BioIoError> {
        if !self.children.is_empty() {
            return Err(error(
                IoErrorKind::InvalidStructure,
                format!("expected text-only XML element {}", self.name),
            ));
        }
        Ok(&self.text)
    }

    pub fn leaf(&self, name: &str) -> Result<Option<&str>, BioIoError> {
        self.one(name)?.map(Self::value).transpose()
    }

    pub fn required_value(&self, name: &str) -> Result<&str, BioIoError> {
        self.required(name)?.value()
    }

    pub fn check(&self, single: &[&str], repeated: &[&str]) -> Result<(), BioIoError> {
        if !self.text.trim().is_empty() {
            return Err(error(
                IoErrorKind::InvalidStructure,
                "mixed text in an XML result container",
            ));
        }
        let mut seen = HashSet::new();
        for child in &self.children {
            if !single.contains(&child.name.as_str()) && !repeated.contains(&child.name.as_str()) {
                return Err(error(
                    IoErrorKind::Unsupported,
                    format!("unsupported XML element {}/{}", self.name, child.name),
                ));
            }
            if single.contains(&child.name.as_str()) && !seen.insert(child.name.as_str()) {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    format!("repeated XML element {}", child.name),
                ));
            }
        }
        Ok(())
    }
}

enum XmlEvent {
    Start(Node, bool),
    End(String),
    Text(String),
    Reference(String),
    CData(String),
    Comment(String),
    Declaration(String),
    Doctype(String),
    Eof,
}

enum OuterAction {
    Container,
    Metadata,
    Statistics,
    SourceError,
}

#[derive(Default)]
struct Report {
    fields: Vec<BlastSourceField>,
    frozen: Option<Arc<[BlastSourceField]>>,
    seen: HashSet<String>,
    program: Option<Arc<str>>,
    version: Option<Arc<str>>,
    target: Option<Arc<str>>,
    root_query_id: Option<String>,
    root_query_description: Option<String>,
    root_query_length: Option<u64>,
    queries: usize,
}

/// Complete-query iterator for explicit legacy NCBI XML or single-file XML2.
///
/// XML is tokenized by quick-xml, with namespace checks and bounded input/event/
/// depth/query buffers. Unknown semantics fail explicitly. The parser never
/// fetches external DTDs, expands declared entities or processes `XInclude`.
// These independent lifecycle/prolog flags also track encoding and completion;
// they are not mutually exclusive alternatives of one state enum.
#[allow(clippy::struct_excessive_bools)]
pub struct BlastXmlReader<R: BufRead> {
    reader: NsReader<BoundedInput<R>>,
    buffer: Vec<u8>,
    options: BlastXmlOptions,
    outer: Vec<String>,
    root_seen: bool,
    root_closed: bool,
    prolog_started: bool,
    declaration_seen: bool,
    doctype_seen: bool,
    ascii_only: bool,
    report: Report,
    report_wrapper_seen: bool,
    queries: usize,
    completed: HashSet<String>,
    completed_id_bytes: usize,
    comments: Vec<BlastComment>,
    source_fields: Vec<BlastSourceField>,
    metadata_bytes: usize,
    done: bool,
}

impl<R: BufRead> BlastXmlReader<R> {
    /// Creates a reader without consuming source input.
    ///
    /// # Errors
    /// Configuration errors are reserved for future option constraints. Zero
    /// budgets are valid and reject the first event requiring that resource.
    pub fn new(reader: R, options: BlastXmlOptions) -> Result<Self, BioIoError> {
        let input = BoundedInput {
            reader,
            total: 0,
            line: 1,
            line_bytes: 0,
            total_limit: options.limits.read.max_total_bytes,
            line_limit: options.limits.read.max_line_bytes,
            event_remaining: options.max_event_bytes,
            capture_remaining: None,
            limited: false,
        };
        let mut reader = NsReader::from_reader(input);
        reader
            .resolver_mut()
            .set_max_declarations_per_element(MAX_NAMESPACE_DECLARATIONS_PER_ELEMENT);
        reader.config_mut().check_comments = true;
        reader.config_mut().check_end_names = true;
        Ok(Self {
            reader,
            buffer: Vec::new(),
            options,
            outer: Vec::new(),
            root_seen: false,
            root_closed: false,
            prolog_started: false,
            declaration_seen: false,
            doctype_seen: false,
            ascii_only: false,
            report: Report::default(),
            report_wrapper_seen: false,
            queries: 0,
            completed: HashSet::new(),
            completed_id_bytes: 0,
            comments: Vec::new(),
            source_fields: Vec::new(),
            metadata_bytes: 0,
            done: false,
        })
    }

    /// Returns document-level comments encountered so far; consume to EOF for all.
    #[must_use]
    pub fn comments(&self) -> &[BlastComment] {
        &self.comments
    }

    /// Returns declarations, outer attributes, global statistics and orphan report metadata.
    ///
    /// Values use decoded XML 1.0 text; references are resolved only for predefined
    /// entities and legal numeric character references, never external definitions.
    #[must_use]
    pub fn source_fields(&self) -> &[BlastSourceField] {
        &self.source_fields
    }

    fn charge_metadata(&mut self, bytes: usize) -> Result<(), BioIoError> {
        self.metadata_bytes = self
            .metadata_bytes
            .checked_add(bytes)
            .filter(|bytes| *bytes <= self.options.max_metadata_bytes)
            .ok_or_else(|| {
                error(
                    IoErrorKind::LimitExceeded,
                    "XML document/report metadata byte limit exceeded",
                )
            })?;
        Ok(())
    }

    fn validate_namespace(&self, namespace: Option<&str>) -> Result<(), BioIoError> {
        let expected = match self.options.format {
            BlastXmlFormat::Legacy => None,
            BlastXmlFormat::Xml2 => Some(NAMESPACE),
        };
        if namespace != expected {
            return Err(error(
                IoErrorKind::Unsupported,
                "XML namespace does not match the selected BLAST profile",
            ));
        }
        Ok(())
    }

    fn node(&self, start: &BytesStart<'_>, line: usize) -> Result<Node, BioIoError> {
        let name = std::str::from_utf8(start.local_name().as_ref())
            .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?
            .to_owned();
        let mut attributes = Vec::new();
        for attribute in start.attributes() {
            let attribute = attribute
                .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?;
            let key = std::str::from_utf8(attribute.key.as_ref())
                .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?;
            let namespace_attribute = key == "xmlns" || key.starts_with("xmlns:");
            let schema_location = if key.ends_with(":schemaLocation") {
                let (namespace, local) = self.reader.resolver().resolve_attribute(attribute.key);
                matches!(namespace, ResolveResult::Bound(namespace) if namespace.as_ref() == b"http://www.w3.org/2001/XMLSchema-instance")
                    && local.as_ref() == b"schemaLocation"
                    && self.outer.is_empty()
                    && self.options.format == BlastXmlFormat::Xml2
            } else {
                false
            };
            if !namespace_attribute && !schema_location {
                return Err(error(
                    IoErrorKind::Unsupported,
                    format!("unsupported XML attribute {key}"),
                ));
            }
            // Retain the profile's decoded source whitespace rather than applying
            // the new dependency API's additional XML attribute normalization.
            let decoded = self
                .reader
                .decoder()
                .decode(attribute.value.as_ref())
                .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?;
            let value = quick_xml::escape::unescape(&decoded)
                .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?
                .into_owned();
            validate_text(&value)?;
            attributes.push(BlastSourceField {
                name: key.to_owned(),
                value,
            });
        }
        Ok(Node {
            name,
            text: String::new(),
            attributes,
            children: Vec::new(),
            line,
        })
    }

    fn event(&mut self) -> Result<(XmlEvent, usize, usize), BioIoError> {
        self.reader.get_mut().event_remaining = self.options.max_event_bytes;
        self.buffer.clear();
        let before = self.reader.get_ref().total;
        let line = self.reader.get_ref().line;
        let event = self.reader.read_resolved_event_into(&mut self.buffer);
        let (namespace, event) = match event {
            Ok((namespace, event)) => {
                let namespace = match namespace {
                    ResolveResult::Bound(namespace) => Some(
                        std::str::from_utf8(namespace.as_ref())
                            .map_err(|failure| {
                                error(IoErrorKind::InvalidStructure, failure.to_string())
                            })?
                            .to_owned(),
                    ),
                    ResolveResult::Unbound => None,
                    ResolveResult::Unknown(_) => {
                        return Err(error(
                            IoErrorKind::Unsupported,
                            "unbound XML namespace prefix",
                        ));
                    }
                };
                (namespace, event.into_owned())
            }
            Err(failure) => {
                let kind = if self.reader.get_ref().limited {
                    IoErrorKind::LimitExceeded
                } else if matches!(failure, quick_xml::Error::Io(_)) {
                    IoErrorKind::Io
                } else if matches!(
                    failure,
                    quick_xml::Error::Syntax(_)
                        | quick_xml::Error::IllFormed(
                            quick_xml::errors::IllFormedError::MissingEndTag(_)
                        )
                ) {
                    IoErrorKind::Truncated
                } else {
                    IoErrorKind::InvalidStructure
                };
                return Err(error(kind, failure.to_string()));
            }
        };
        if self.ascii_only && !event.as_ref().is_ascii() {
            return Err(error(
                IoErrorKind::InvalidStructure,
                "non-ASCII source bytes contradict the US-ASCII XML declaration",
            ));
        }
        if matches!(event, Event::Start(_) | Event::Empty(_) | Event::End(_)) {
            self.validate_namespace(namespace.as_deref())?;
        }
        let event = match event {
            Event::Start(start) => XmlEvent::Start(self.node(&start, line)?, false),
            Event::Empty(start) => XmlEvent::Start(self.node(&start, line)?, true),
            Event::End(end) => XmlEvent::End(
                std::str::from_utf8(end.local_name().as_ref())
                    .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?
                    .to_owned(),
            ),
            Event::Text(text) => XmlEvent::Text(
                text.xml10_content()
                    .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?
                    .into_owned(),
            ),
            Event::CData(text) => XmlEvent::CData(
                text.xml10_content()
                    .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?
                    .into_owned(),
            ),
            Event::Comment(text) => XmlEvent::Comment(
                text.xml10_content()
                    .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?
                    .into_owned(),
            ),
            Event::GeneralRef(reference) => {
                XmlEvent::Reference(reference_value(&reference.decode().map_err(
                    |failure| error(IoErrorKind::InvalidStructure, failure.to_string()),
                )?)?)
            }
            Event::Decl(declaration) => self.declaration(&declaration)?,
            Event::DocType(doctype) => XmlEvent::Doctype(
                doctype
                    .xml10_content()
                    .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?
                    .into_owned(),
            ),
            Event::PI(_) => {
                return Err(error(
                    IoErrorKind::Unsupported,
                    "XML processing instructions are outside this BLAST profile",
                ));
            }
            Event::Eof => XmlEvent::Eof,
        };
        Ok((event, before, self.reader.get_ref().total - before))
    }

    fn declaration(&mut self, declaration: &BytesDecl<'_>) -> Result<XmlEvent, BioIoError> {
        let declaration_text = std::str::from_utf8(declaration.as_ref())
            .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?;
        let start = BytesStart::from_content(declaration_text, 3);
        let mut previous = None;
        for attribute in start.attributes() {
            let attribute = attribute
                .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?;
            let position = match attribute.key.as_ref() {
                b"version" => 0,
                b"encoding" => 1,
                b"standalone" => 2,
                _ => {
                    return Err(error(
                        IoErrorKind::InvalidStructure,
                        "unknown XML declaration attribute",
                    ));
                }
            };
            if previous.is_some_and(|previous| position <= previous) {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "XML declaration attributes are repeated or out of order",
                ));
            }
            previous = Some(position);
        }
        let version = declaration
            .version()
            .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?;
        if version.as_ref() != b"1.0" {
            return Err(error(IoErrorKind::Unsupported, "only XML 1.0 is supported"));
        }
        if let Some(encoding) = declaration.encoding() {
            let encoding = encoding
                .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?;
            if !encoding.eq_ignore_ascii_case(b"utf-8")
                && !encoding.eq_ignore_ascii_case(b"us-ascii")
            {
                return Err(error(
                    IoErrorKind::Unsupported,
                    "only UTF-8/US-ASCII XML encodings are supported",
                ));
            }
            self.ascii_only = encoding.eq_ignore_ascii_case(b"us-ascii");
        }
        if let Some(standalone) = declaration.standalone() {
            let standalone = standalone
                .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?;
            if !matches!(standalone.as_ref(), b"yes" | b"no") {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "XML standalone must be yes or no",
                ));
            }
        }
        Ok(XmlEvent::Declaration(
            std::str::from_utf8(declaration.as_ref())
                .map_err(|failure| error(IoErrorKind::InvalidStructure, failure.to_string()))?
                .to_owned(),
        ))
    }

    fn capture(
        &mut self,
        root: Node,
        empty: bool,
        first_bytes: usize,
        query: bool,
    ) -> Result<(Node, Vec<BlastComment>), BioIoError> {
        let budget = if query {
            self.options.limits.max_query_bytes
        } else {
            self.options
                .max_metadata_bytes
                .saturating_sub(self.metadata_bytes)
        };
        if first_bytes > budget {
            return Err(error(
                IoErrorKind::LimitExceeded,
                "XML capture byte limit exceeded",
            ));
        }
        if self.outer.len() + 1 > self.options.max_depth {
            return Err(error(
                IoErrorKind::LimitExceeded,
                "XML nesting depth limit exceeded",
            ));
        }
        if empty {
            if !query {
                self.charge_metadata(first_bytes)?;
            }
            return Ok((root, Vec::new()));
        }
        self.reader.get_mut().capture_remaining = Some(budget - first_bytes);
        let start_position = self.reader.get_ref().total - first_bytes;
        let mut stack = vec![root];
        let mut comments = Vec::new();
        let mut counts = (0usize, 0usize);
        loop {
            let (event, _, _) = self.event()?;
            match event {
                XmlEvent::Start(node, empty) => {
                    self.capture_limits(&node, stack.len(), query, &mut counts)?;
                    if empty {
                        stack
                            .last_mut()
                            .ok_or_else(|| {
                                error(IoErrorKind::InvalidStructure, "missing XML parent")
                            })?
                            .children
                            .push(node);
                    } else {
                        stack.push(node);
                    }
                }
                XmlEvent::End(name) => {
                    let node = stack.pop().ok_or_else(|| {
                        error(IoErrorKind::InvalidStructure, "unmatched XML end tag")
                    })?;
                    if node.name != name {
                        return Err(error(
                            IoErrorKind::InvalidStructure,
                            "mismatched XML end tag",
                        ));
                    }
                    if let Some(parent) = stack.last_mut() {
                        parent.children.push(node);
                    } else {
                        self.reader.get_mut().capture_remaining = None;
                        if !query {
                            self.charge_metadata(self.reader.get_ref().total - start_position)?;
                        }
                        return Ok((node, comments));
                    }
                }
                XmlEvent::Text(text) | XmlEvent::Reference(text) | XmlEvent::CData(text) => {
                    validate_text(&text)?;
                    stack
                        .last_mut()
                        .ok_or_else(|| {
                            error(IoErrorKind::InvalidStructure, "text outside an XML node")
                        })?
                        .text
                        .push_str(&text);
                }
                XmlEvent::Comment(text) => {
                    validate_text(&text)?;
                    comments.push(BlastComment {
                        text,
                        line: self.reader.get_ref().line,
                    });
                }
                XmlEvent::Eof => {
                    return Err(error(
                        IoErrorKind::Truncated,
                        "XML ended inside a query/metadata element",
                    ));
                }
                XmlEvent::Declaration(_) | XmlEvent::Doctype(_) => {
                    return Err(error(
                        IoErrorKind::InvalidStructure,
                        "declaration inside an XML result element",
                    ));
                }
            }
        }
    }

    fn capture_limits(
        &self,
        node: &Node,
        depth: usize,
        query: bool,
        counts: &mut (usize, usize),
    ) -> Result<(), BioIoError> {
        if self.outer.len() + depth + 1 > self.options.max_depth {
            return Err(error(
                IoErrorKind::LimitExceeded,
                "XML nesting depth limit exceeded",
            ));
        }
        if query && node.name == "Hit" {
            counts.0 += 1;
            if counts.0 > self.options.limits.max_hits_per_query {
                return Err(error(IoErrorKind::LimitExceeded, "XML hit limit exceeded"));
            }
        }
        if query && node.name == "Hsp" {
            counts.1 += 1;
            if counts.1 > self.options.limits.max_hsps_per_query {
                return Err(error(IoErrorKind::LimitExceeded, "XML HSP limit exceeded"));
            }
        }
        Ok(())
    }

    fn add_query(
        &mut self,
        node: &Node,
        comments: Vec<BlastComment>,
    ) -> Result<BlastQuery, BioIoError> {
        if self.queries == self.options.limits.max_queries {
            return Err(error(
                IoErrorKind::LimitExceeded,
                "XML query-count limit exceeded",
            ));
        }
        let report_fields = self
            .report
            .frozen
            .get_or_insert_with(|| std::mem::take(&mut self.report.fields).into())
            .clone();
        let mut query = model::query(node, self.options.format, &self.report)?;
        query.report_fields = report_fields;
        query.comments = comments;
        if self.completed.contains(&query.id) {
            return Err(error(
                IoErrorKind::Unsupported,
                "repeated XML query ID/iteration rounds are outside the noniterative profile",
            ));
        }
        self.completed_id_bytes = self
            .completed_id_bytes
            .checked_add(query.id.len())
            .filter(|bytes| *bytes <= self.options.limits.max_query_id_bytes)
            .ok_or_else(|| {
                error(
                    IoErrorKind::LimitExceeded,
                    "XML query identifier-byte limit exceeded",
                )
            })?;
        self.completed.insert(query.id.clone());
        self.queries += 1;
        self.report.queries += 1;
        Ok(query)
    }

    fn outer_attributes(&mut self, node: &Node) {
        let prefix = self.outer.join("/");
        self.source_fields
            .extend(node.attributes.iter().map(|attribute| BlastSourceField {
                name: format!("{prefix}/{}@{}", node.name, attribute.name),
                value: attribute.value.clone(),
            }));
    }

    fn unique(&mut self, name: &str) -> Result<(), BioIoError> {
        if !self.report.seen.insert(name.to_owned()) {
            return Err(error(
                IoErrorKind::InvalidStructure,
                format!("repeated XML container/metadata {name}"),
            ));
        }
        Ok(())
    }

    fn metadata(&mut self, node: &Node) -> Result<(), BioIoError> {
        self.unique(&node.name)?;
        model::report_metadata(node, self.options.format, &mut self.report)?;
        self.report.fields.extend(model::flatten(node, ""));
        Ok(())
    }

    fn begin_outer(
        &mut self,
        node: Node,
        empty: bool,
        bytes: usize,
    ) -> Result<Option<BlastQuery>, BioIoError> {
        let name = node.name.as_str();
        if self.outer.is_empty() {
            if self.options.max_depth == 0 {
                return Err(error(
                    IoErrorKind::LimitExceeded,
                    "XML nesting depth limit exceeded",
                ));
            }
            if self.root_seen {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "multiple XML document roots",
                ));
            }
            let expected = match self.options.format {
                BlastXmlFormat::Legacy => "BlastOutput",
                BlastXmlFormat::Xml2 => "BlastXML2",
            };
            if name != expected {
                return Err(error(
                    IoErrorKind::Unsupported,
                    "XML root does not match selected BLAST profile",
                ));
            }
            self.charge_metadata(bytes)?;
            self.root_seen = true;
            self.outer_attributes(&node);
            let name = node.name;
            self.outer.push(name.clone());
            if empty {
                self.end_outer(&name)?;
            }
            return Ok(None);
        }
        let parent = self.outer.last().cloned().unwrap_or_default();
        let query_node = matches!(
            (self.options.format, parent.as_str(), name),
            (
                BlastXmlFormat::Legacy,
                "BlastOutput_iterations",
                "Iteration"
            ) | (BlastXmlFormat::Xml2, "search", "Search")
        );
        if query_node {
            if self.queries == self.options.limits.max_queries {
                return Err(error(
                    IoErrorKind::LimitExceeded,
                    "XML query-count limit exceeded",
                ));
            }
            let (node, comments) = self.capture(node, empty, bytes, true)?;
            return self.add_query(&node, comments).map(Some);
        }
        let action = match self.options.format {
            BlastXmlFormat::Legacy => self.legacy_action(&parent, name)?,
            BlastXmlFormat::Xml2 => self.xml2_action(&parent, name)?,
        };
        if !matches!(action, OuterAction::Container) {
            let (node, comments) = self.capture(node, empty, bytes, false)?;
            self.comments.extend(comments);
            match action {
                OuterAction::Metadata => self.metadata(&node)?,
                OuterAction::Statistics => {
                    model::statistics_wrapper(&node, "Statistics", true)?;
                    self.source_fields.extend(model::flatten(&node, ""));
                }
                OuterAction::SourceError => return Err(model::source_error(&node)?),
                OuterAction::Container => {
                    return Err(error(
                        IoErrorKind::InvalidStructure,
                        "invalid XML parser capture action",
                    ));
                }
            }
            return Ok(None);
        }
        self.charge_metadata(bytes)?;
        if self.outer.len() + 1 > self.options.max_depth {
            return Err(error(
                IoErrorKind::LimitExceeded,
                "XML nesting depth limit exceeded",
            ));
        }
        self.outer_attributes(&node);
        let name = node.name;
        self.outer.push(name.clone());
        if empty {
            self.end_outer(&name)?;
        }
        Ok(None)
    }

    fn legacy_action(&mut self, parent: &str, name: &str) -> Result<OuterAction, BioIoError> {
        if parent != "BlastOutput" {
            return Err(error(
                IoErrorKind::Unsupported,
                format!("unsupported legacy XML structure {parent}/{name}"),
            ));
        }
        match name {
            "BlastOutput_iterations" => {
                self.unique(name)?;
                model::validate_report(&self.report, BlastXmlFormat::Legacy, false)?;
                Ok(OuterAction::Container)
            }
            "BlastOutput_mbstat" => {
                self.unique(name)?;
                if !self.report.seen.contains("BlastOutput_iterations") {
                    return Err(error(
                        IoErrorKind::InvalidStructure,
                        "mbstat precedes XML iterations",
                    ));
                }
                Ok(OuterAction::Statistics)
            }
            _ => {
                if self.report.seen.contains("BlastOutput_iterations") {
                    return Err(error(
                        IoErrorKind::InvalidStructure,
                        "report metadata follows XML query results",
                    ));
                }
                Ok(OuterAction::Metadata)
            }
        }
    }

    fn xml2_action(&mut self, parent: &str, name: &str) -> Result<OuterAction, BioIoError> {
        match (parent, name) {
            ("BlastXML2", "BlastOutput2") => {
                self.report_wrapper_seen = false;
                self.report = Report::default();
            }
            ("BlastOutput2", "report") => {
                if self.report_wrapper_seen {
                    return Err(error(IoErrorKind::InvalidStructure, "repeated XML2 report"));
                }
                self.report_wrapper_seen = true;
            }
            ("report", "Report") => self.unique("Report")?,
            ("BlastOutput2", "error") => return Ok(OuterAction::SourceError),
            ("Report", "results") => {
                self.unique("results")?;
                model::validate_report(&self.report, BlastXmlFormat::Xml2, false)?;
            }
            ("results", "Results") => self.unique("Results")?,
            ("Results", "search") => self.unique("search")?,
            ("Results", "iterations" | "bl2seq") => {
                return Err(error(
                    IoErrorKind::Unsupported,
                    "XML2 iterative/bl2seq result modes require a separately defined profile",
                ));
            }
            ("Report", _) => {
                if self.report.seen.contains("results") {
                    return Err(error(
                        IoErrorKind::InvalidStructure,
                        "XML2 metadata follows query results",
                    ));
                }
                return Ok(OuterAction::Metadata);
            }
            _ => {
                return Err(error(
                    IoErrorKind::Unsupported,
                    format!("unsupported XML2 structure {parent}/{name}"),
                ));
            }
        }
        Ok(OuterAction::Container)
    }

    fn end_outer(&mut self, name: &str) -> Result<(), BioIoError> {
        if self.outer.pop().as_deref() != Some(name) {
            return Err(error(
                IoErrorKind::InvalidStructure,
                "unmatched XML container end",
            ));
        }
        match (self.options.format, name) {
            (BlastXmlFormat::Legacy, "BlastOutput") => {
                model::validate_report(&self.report, BlastXmlFormat::Legacy, true)?;
                if self.report.queries == 0 {
                    self.source_fields.append(&mut self.report.fields);
                }
                self.root_closed = true;
            }
            (BlastXmlFormat::Xml2, "Search") => {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "query end escaped its capture",
                ));
            }
            (BlastXmlFormat::Xml2, "search") if self.report.queries != 1 => {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "XML2 search must contain exactly one Search",
                ));
            }
            (BlastXmlFormat::Xml2, "report") if !self.report.seen.contains("Report") => {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "XML2 report lacks Report",
                ));
            }
            (BlastXmlFormat::Xml2, "Report") => {
                model::validate_report(&self.report, BlastXmlFormat::Xml2, true)?;
            }
            (BlastXmlFormat::Xml2, "BlastOutput2") if !self.report_wrapper_seen => {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "XML2 result lacks a report",
                ));
            }
            (BlastXmlFormat::Xml2, "BlastXML2") => self.root_closed = true,
            _ => {}
        }
        Ok(())
    }

    fn read_query(&mut self) -> Result<Option<BlastQuery>, BioIoError> {
        loop {
            let (event, _, bytes) = self.event()?;
            match event {
                XmlEvent::Start(node, empty) => {
                    self.prolog_started = true;
                    if let Some(query) = self.begin_outer(node, empty, bytes)? {
                        return Ok(Some(query));
                    }
                }
                XmlEvent::End(name) => {
                    self.charge_metadata(bytes)?;
                    self.end_outer(&name)?;
                }
                XmlEvent::Text(text) => {
                    self.charge_metadata(bytes)?;
                    validate_text(&text)?;
                    if !text.trim().is_empty() {
                        return Err(error(
                            IoErrorKind::InvalidStructure,
                            "text outside XML data fields",
                        ));
                    }
                    self.prolog_started = true;
                }
                XmlEvent::Reference(_) | XmlEvent::CData(_) => {
                    return Err(error(
                        IoErrorKind::InvalidStructure,
                        "reference/CDATA outside an XML data field",
                    ));
                }
                XmlEvent::Comment(text) => {
                    self.charge_metadata(bytes)?;
                    validate_text(&text)?;
                    self.comments.push(BlastComment {
                        text,
                        line: self.reader.get_ref().line,
                    });
                    self.prolog_started = true;
                }
                XmlEvent::Declaration(value) => {
                    self.charge_metadata(bytes)?;
                    if self.declaration_seen || self.prolog_started {
                        return Err(error(
                            IoErrorKind::InvalidStructure,
                            "XML declaration must occur once at the beginning",
                        ));
                    }
                    self.declaration_seen = true;
                    self.source_fields.push(BlastSourceField {
                        name: "xml-declaration".to_owned(),
                        value,
                    });
                }
                XmlEvent::Doctype(value) => {
                    self.charge_metadata(bytes)?;
                    if self.root_seen
                        || self.doctype_seen
                        || self.options.format != BlastXmlFormat::Legacy
                        || !known_doctype(&value)
                    {
                        return Err(error(
                            IoErrorKind::Unsupported,
                            "only the known legacy BLAST external DOCTYPE is accepted; internal entities are unsupported",
                        ));
                    }
                    self.doctype_seen = true;
                    self.prolog_started = true;
                    self.source_fields.push(BlastSourceField {
                        name: "doctype".to_owned(),
                        value,
                    });
                }
                XmlEvent::Eof => {
                    if !self.outer.is_empty() || !self.root_closed {
                        return Err(error(
                            IoErrorKind::Truncated,
                            "incomplete XML result document",
                        ));
                    }
                    return Ok(None);
                }
            }
        }
    }
}

impl<R: BufRead> Iterator for BlastXmlReader<R> {
    type Item = Result<BlastQuery, BioIoError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        match self.read_query() {
            Ok(Some(query)) => Some(Ok(query)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(failure) => {
                self.done = true;
                Some(Err(
                    failure.at((self.queries + 1).max(1), self.reader.get_ref().line)
                ))
            }
        }
    }
}

impl<R: BufRead> FusedIterator for BlastXmlReader<R> {}

fn reference_value(reference: &str) -> Result<String, BioIoError> {
    let escaped = format!("&{reference};");
    quick_xml::escape::unescape(&escaped)
        .map(std::borrow::Cow::into_owned)
        .map_err(|failure| {
            error(
                IoErrorKind::Unsupported,
                format!("unsupported XML entity/reference: {failure}"),
            )
        })
}

fn validate_text(value: &str) -> Result<(), BioIoError> {
    if value.chars().any(|character| {
        let value = u32::from(character);
        !(matches!(value, 0x9 | 0xa | 0xd)
            || (0x20..=0xd7ff).contains(&value)
            || (0xe000..=0xfffd).contains(&value)
            || (0x0001_0000..=0x0010_ffff).contains(&value))
    }) {
        return Err(error(
            IoErrorKind::InvalidStructure,
            "illegal XML 1.0 character",
        ));
    }
    Ok(())
}

fn known_doctype(value: &str) -> bool {
    if value.contains(['[', ']', '<', '>']) {
        return false;
    }
    let normalized = value.split_whitespace().collect::<Vec<_>>().join(" ");
    normalized == "BlastOutput"
        || [
            "NCBI_BlastOutput.dtd",
            "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd",
            "https://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd",
        ]
        .iter()
        .any(|location| {
            normalized
                == format!("BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"{location}\"")
        })
}
