//! NCBI DTD/XSD model checks and source-preserving result construction.

use std::collections::HashMap;
use std::sync::Arc;

use super::{Node, Report, error, relabel};
use crate::io::blast::fields;
use crate::io::blast::{
    BlastHit, BlastHsp, BlastQuery, BlastSourceField, BlastStrand, BlastXmlFormat,
};
use crate::io::{BioIoError, IoErrorKind};

pub(super) fn flatten(node: &Node, prefix: &str) -> Vec<BlastSourceField> {
    let path = if prefix.is_empty() {
        node.name.clone()
    } else {
        format!("{prefix}/{}", node.name)
    };
    let mut output: Vec<_> = node
        .attributes
        .iter()
        .map(|attribute| BlastSourceField {
            name: format!("{path}@{}", attribute.name),
            value: attribute.value.clone(),
        })
        .collect();
    if node.children.is_empty() {
        output.push(BlastSourceField {
            name: path,
            value: node.text.clone(),
        });
    } else {
        let mut positions = HashMap::<&str, usize>::new();
        for child in &node.children {
            let repeated = matches!(child.name.as_str(), "HitDescr" | "subjects" | "Range");
            if repeated {
                let position = positions.entry(child.name.as_str()).or_default();
                let nested = flatten(child, "");
                let base = format!("{path}/{}[{}]", child.name, *position);
                output.extend(nested.into_iter().map(|field| BlastSourceField {
                    name: format!("{base}{}", &field.name[child.name.len()..]),
                    value: field.value,
                }));
                *position += 1;
            } else {
                output.extend(flatten(child, &path));
            }
        }
    }
    output
}

fn flat_children(node: &Node, excluded: &[&str]) -> Vec<BlastSourceField> {
    let mut fields = Vec::new();
    fields.extend(node.attributes.iter().map(|attribute| BlastSourceField {
        name: format!("@{}", attribute.name),
        value: attribute.value.clone(),
    }));
    for child in &node.children {
        if !excluded.contains(&child.name.as_str()) {
            fields.extend(flatten(child, ""));
        }
    }
    fields
}

fn uint(value: &str, name: &str, positive: bool) -> Result<u64, BioIoError> {
    fields::unsigned(value.trim(), name, positive).map_err(relabel)
}

fn number(value: &str, name: &str, nonnegative: bool) -> Result<f64, BioIoError> {
    fields::finite(value.trim(), name, nonnegative).map_err(relabel)
}

fn signed(value: &str, name: &str) -> Result<i64, BioIoError> {
    value.trim().parse().map_err(|_| {
        error(
            IoErrorKind::InvalidStructure,
            format!("{name} must be a signed 64-bit integer"),
        )
    })
}

fn identity(value: &str, name: &str) -> Result<String, BioIoError> {
    if value.is_empty() || value.chars().any(char::is_whitespace) {
        return Err(error(
            IoErrorKind::InvalidHeader,
            format!("{name} requires a nonempty whitespace-free source identity"),
        ));
    }
    Ok(value.to_owned())
}

fn program(value: &str) -> Result<(), BioIoError> {
    if !matches!(
        value.to_ascii_lowercase().as_str(),
        "blastn" | "blastp" | "blastx" | "tblastn" | "tblastx"
    ) {
        return Err(error(
            IoErrorKind::Unsupported,
            "XML profile supports blastn/blastp/blastx/tblastn/tblastx; iterative/other programs need separate semantics",
        ));
    }
    Ok(())
}

pub(super) fn report_metadata(
    node: &Node,
    format: BlastXmlFormat,
    report: &mut Report,
) -> Result<(), BioIoError> {
    match format {
        BlastXmlFormat::Legacy => match node.name.as_str() {
            "BlastOutput_program" => {
                let value = node.value()?;
                program(value)?;
                report.program = Some(value.into());
            }
            "BlastOutput_version" => {
                report.version = Some(node.value()?.into());
            }
            "BlastOutput_reference" | "BlastOutput_query-seq" => {
                node.value()?;
            }
            "BlastOutput_db" => {
                report.target = Some(node.value()?.into());
            }
            "BlastOutput_query-ID" => {
                report.root_query_id = Some(identity(node.value()?, &node.name)?);
            }
            "BlastOutput_query-def" => {
                report.root_query_description = Some(node.value()?.to_owned());
            }
            "BlastOutput_query-len" => {
                report.root_query_length = Some(uint(node.value()?, &node.name, true)?);
            }
            "BlastOutput_param" => {
                node.check(&["Parameters"], &[])?;
                parameters(node.required("Parameters")?, true)?;
            }
            _ => {
                return Err(error(
                    IoErrorKind::Unsupported,
                    format!("unsupported legacy XML report field {}", node.name),
                ));
            }
        },
        BlastXmlFormat::Xml2 => match node.name.as_str() {
            "program" => {
                let value = node.value()?;
                program(value)?;
                report.program = Some(value.into());
            }
            "version" => {
                report.version = Some(node.value()?.into());
            }
            "reference" => {
                node.value()?;
            }
            "search-target" => {
                node.check(&["Target"], &[])?;
                let target = node.required("Target")?;
                target.check(&["db"], &["subjects"])?;
                let db = target.leaf("db")?;
                let subjects: Vec<_> = target
                    .children
                    .iter()
                    .filter(|child| child.name == "subjects")
                    .collect();
                if db.is_some() != subjects.is_empty() {
                    return Err(error(
                        IoErrorKind::InvalidStructure,
                        "XML2 Target requires exactly one of db or subjects",
                    ));
                }
                for subject in subjects {
                    subject.value()?;
                }
                report.target = db.map(Into::into);
            }
            "params" => {
                node.check(&["Parameters"], &[])?;
                parameters(node.required("Parameters")?, false)?;
            }
            _ => {
                return Err(error(
                    IoErrorKind::Unsupported,
                    format!("unsupported XML2 report field {}", node.name),
                ));
            }
        },
    }
    Ok(())
}

pub(super) fn validate_report(
    report: &Report,
    format: BlastXmlFormat,
    final_check: bool,
) -> Result<(), BioIoError> {
    let required: &[&str] = match format {
        BlastXmlFormat::Legacy => &[
            "BlastOutput_program",
            "BlastOutput_version",
            "BlastOutput_reference",
            "BlastOutput_db",
            "BlastOutput_query-ID",
            "BlastOutput_query-def",
            "BlastOutput_query-len",
            "BlastOutput_param",
        ],
        BlastXmlFormat::Xml2 => &["program", "version", "reference", "search-target", "params"],
    };
    for name in required {
        if !report.seen.contains(*name) {
            return Err(error(
                IoErrorKind::InvalidStructure,
                format!("missing XML report field {name}"),
            ));
        }
    }
    if final_check {
        let required: &[&str] = match format {
            BlastXmlFormat::Legacy => &["BlastOutput_iterations"],
            BlastXmlFormat::Xml2 => &["results", "Results", "search"],
        };
        for name in required {
            if !report.seen.contains(*name) {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    format!("missing XML result container {name}"),
                ));
            }
        }
    }
    Ok(())
}

fn parameters(node: &Node, legacy: bool) -> Result<(), BioIoError> {
    let prefix = if legacy { "Parameters_" } else { "" };
    let common = [
        "matrix",
        "expect",
        "include",
        "sc-match",
        "sc-mismatch",
        "gap-open",
        "gap-extend",
        "filter",
        "pattern",
        "entrez-query",
    ];
    let extension = ["cbs", "query-gencode", "db-gencode", "bl2seq-mode"];
    let allowed: Vec<_> = common
        .iter()
        .chain(extension.iter().filter(|_| !legacy))
        .map(|name| format!("{prefix}{name}"))
        .collect();
    let allowed: Vec<_> = allowed.iter().map(String::as_str).collect();
    node.check(&allowed, &[])?;
    number(
        node.required_value(&format!("{prefix}expect"))?,
        "expect",
        true,
    )?;
    if legacy {
        uint(
            node.required_value("Parameters_gap-open")?,
            "gap-open",
            false,
        )?;
        uint(
            node.required_value("Parameters_gap-extend")?,
            "gap-extend",
            false,
        )?;
    }
    for child in &node.children {
        let name = child.name.strip_prefix(prefix).unwrap_or(&child.name);
        let value = child.value()?;
        match name {
            "expect" | "include" => {
                number(value, name, true)?;
            }
            "sc-match" | "sc-mismatch" | "cbs" => {
                signed(value, name)?;
            }
            "gap-open" | "gap-extend" => {
                uint(value, name, false)?;
            }
            "query-gencode" | "db-gencode" => {
                uint(value, name, true)?;
            }
            _ => {}
        }
    }
    Ok(())
}

pub(super) fn statistics_wrapper(node: &Node, child: &str, legacy: bool) -> Result<(), BioIoError> {
    node.check(&[child], &[])?;
    statistics(node.required(child)?, legacy)
}

fn statistics(node: &Node, legacy: bool) -> Result<(), BioIoError> {
    let prefix = if legacy { "Statistics_" } else { "" };
    let allowed = [
        "db-num",
        "db-len",
        "hsp-len",
        "eff-space",
        "kappa",
        "lambda",
        "entropy",
    ];
    let allowed: Vec<_> = allowed.map(|name| format!("{prefix}{name}")).to_vec();
    let names: Vec<_> = allowed.iter().map(String::as_str).collect();
    node.check(&names, &[])?;
    for name in &allowed {
        if legacy || !name.ends_with("db-num") && !name.ends_with("db-len") {
            node.required(name)?;
        }
    }
    for child in &node.children {
        let name = child.name.strip_prefix(prefix).unwrap_or(&child.name);
        let value = child.value()?;
        match name {
            "db-num" | "db-len" | "hsp-len" => {
                uint(value, name, false)?;
            }
            "eff-space" => {
                number(value, name, true)?;
            }
            _ => {
                number(value, name, false)?;
            }
        }
    }
    Ok(())
}

pub(super) fn query(
    node: &Node,
    format: BlastXmlFormat,
    report: &Report,
) -> Result<BlastQuery, BioIoError> {
    let legacy = format == BlastXmlFormat::Legacy;
    let QueryMetadata {
        id,
        identity_field,
        description,
        length,
        hits_name,
    } = if legacy {
        legacy_query_metadata(node, report)?
    } else {
        xml2_query_metadata(node)?
    };
    let mut hits = Vec::new();
    if let Some(container) = node.one(hits_name)? {
        container.check(&[], &["Hit"])?;
        for hit in &container.children {
            hits.push(parse_hit(hit, legacy, length)?);
        }
    }
    Ok(BlastQuery {
        reported_id: Some(id.clone()),
        id,
        identity_field: identity_field.to_owned(),
        description,
        length,
        program: report.program.clone(),
        version: report.version.clone(),
        target: report.target.clone(),
        fields: Vec::new(),
        source_fields: flat_children(node, &[hits_name]),
        report_fields: Arc::default(),
        comments: Vec::new(),
        declared_rows: None,
        hits,
        line: node.line,
    })
}

struct QueryMetadata {
    id: String,
    identity_field: &'static str,
    description: Option<String>,
    length: Option<u64>,
    hits_name: &'static str,
}

fn legacy_query_metadata(node: &Node, report: &Report) -> Result<QueryMetadata, BioIoError> {
    node.check(
        &[
            "Iteration_iter-num",
            "Iteration_query-ID",
            "Iteration_query-def",
            "Iteration_query-len",
            "Iteration_hits",
            "Iteration_stat",
            "Iteration_message",
        ],
        &[],
    )?;
    uint(
        node.required_value("Iteration_iter-num")?,
        "iteration number",
        true,
    )?;
    let explicit = node.leaf("Iteration_query-ID")?;
    let id = explicit
        .or(report.root_query_id.as_deref())
        .ok_or_else(|| {
            error(
                IoErrorKind::InvalidHeader,
                "legacy XML query has no source ID",
            )
        })?;
    let description = node
        .leaf("Iteration_query-def")?
        .map(str::to_owned)
        .or_else(|| report.root_query_description.clone());
    let length = node
        .leaf("Iteration_query-len")?
        .map(|value| uint(value, "query length", true))
        .transpose()?
        .or(report.root_query_length);
    if let Some(stat) = node.one("Iteration_stat")? {
        statistics_wrapper(stat, "Statistics", true)?;
    }
    if let Some(message) = node.one("Iteration_message")? {
        message.value()?;
    }
    Ok(QueryMetadata {
        id: identity(id, "query ID")?,
        identity_field: if explicit.is_some() {
            "Iteration_query-ID"
        } else {
            "BlastOutput_query-ID"
        },
        description,
        length,
        hits_name: "Iteration_hits",
    })
}

fn xml2_query_metadata(node: &Node) -> Result<QueryMetadata, BioIoError> {
    node.check(
        &[
            "query-id",
            "query-title",
            "query-len",
            "query-masking",
            "hits",
            "stat",
            "message",
        ],
        &[],
    )?;
    let id = node.leaf("query-id")?.ok_or_else(|| {
        error(
            IoErrorKind::Unsupported,
            "XML2 queries without explicit query-id require a different identity policy",
        )
    })?;
    let description = node.leaf("query-title")?.map(str::to_owned);
    let length = node
        .leaf("query-len")?
        .map(|value| uint(value, "query length", true))
        .transpose()?;
    if let Some(stat) = node.one("stat")? {
        statistics_wrapper(stat, "Statistics", false)?;
    }
    if let Some(message) = node.one("message")? {
        message.value()?;
    }
    if let Some(masking) = node.one("query-masking")? {
        masking.check(&[], &["Range"])?;
        for range in &masking.children {
            range.check(&["from", "to"], &[])?;
            let coordinates = fields::coordinates(
                Some(range.required_value("from")?.trim()),
                Some(range.required_value("to")?.trim()),
                length,
                "masking",
            )
            .map_err(relabel)?
            .ok_or_else(|| error(IoErrorKind::InvalidStructure, "missing masking range"))?;
            if coordinates.reported_start > coordinates.reported_end {
                return Err(error(
                    IoErrorKind::InvalidStructure,
                    "masking range endpoints must be ascending",
                ));
            }
        }
    }
    Ok(QueryMetadata {
        id: identity(id, "query ID")?,
        identity_field: "query-id",
        description,
        length,
        hits_name: "hits",
    })
}

fn parse_hit(node: &Node, legacy: bool, query_length: Option<u64>) -> Result<BlastHit, BioIoError> {
    let (id, identity_field, description, length, hsps_name) = if legacy {
        node.check(
            &[
                "Hit_num",
                "Hit_id",
                "Hit_def",
                "Hit_accession",
                "Hit_len",
                "Hit_hsps",
            ],
            &[],
        )?;
        uint(node.required_value("Hit_num")?, "hit number", true)?;
        node.required_value("Hit_accession")?;
        let id = identity(node.required_value("Hit_id")?, "Hit_id")?;
        let description = Some(node.required_value("Hit_def")?.to_owned());
        let length = uint(node.required_value("Hit_len")?, "Hit_len", true)?;
        (id, "Hit_id", description, length, "Hit_hsps")
    } else {
        node.check(&["num", "description", "len", "hsps"], &[])?;
        uint(node.required_value("num")?, "hit number", true)?;
        let descriptions = node.required("description")?;
        descriptions.check(&[], &["HitDescr"])?;
        for description in &descriptions.children {
            description.check(&["id", "accession", "title", "taxid", "sciname"], &[])?;
            identity(description.required_value("id")?, "HitDescr/id")?;
            for name in ["accession", "title", "sciname"] {
                description.leaf(name)?;
            }
            if let Some(taxid) = description.leaf("taxid")? {
                uint(taxid, "taxid", false)?;
            }
        }
        let first = descriptions.children.first().ok_or_else(|| {
            error(
                IoErrorKind::Unsupported,
                "XML2 hit without a HitDescr identity is outside this profile",
            )
        })?;
        let id = identity(first.required_value("id")?, "HitDescr/id")?;
        let description = first.leaf("title")?.map(str::to_owned);
        let length = uint(node.required_value("len")?, "subject length", true)?;
        (
            id,
            "description/HitDescr[0]/id",
            description,
            length,
            "hsps",
        )
    };
    let mut hsps = Vec::new();
    if let Some(container) = node.one(hsps_name)? {
        container.check(&[], &["Hsp"])?;
        for hsp in &container.children {
            hsps.push(parse_hsp(hsp, legacy, query_length, Some(length))?);
        }
    }
    Ok(BlastHit {
        id,
        identity_field: identity_field.to_owned(),
        description,
        length: Some(length),
        source_fields: flat_children(node, &[hsps_name]),
        hsps,
    })
}

fn parse_hsp(
    node: &Node,
    legacy: bool,
    query_length: Option<u64>,
    hit_length: Option<u64>,
) -> Result<BlastHsp, BioIoError> {
    let prefix = if legacy { "Hsp_" } else { "" };
    check_hsp_fields(node, legacy, prefix)?;
    let get = |name: &str| node.leaf(&format!("{prefix}{name}"));
    let required = |name: &str| node.required_value(&format!("{prefix}{name}"));
    uint(required("num")?, "HSP number", true)?;
    let query_coordinates = fields::coordinates(
        Some(required("query-from")?.trim()),
        Some(required("query-to")?.trim()),
        query_length,
        "query",
    )
    .map_err(relabel)?;
    let hit_coordinates = fields::coordinates(
        Some(required("hit-from")?.trim()),
        Some(required("hit-to")?.trim()),
        hit_length,
        "subject",
    )
    .map_err(relabel)?;
    let query_frame = get("query-frame")?
        .map(|value| fields::frame(value.trim(), "query-frame").map_err(relabel))
        .transpose()?;
    let hit_frame = get("hit-frame")?
        .map(|value| fields::frame(value.trim(), "hit-frame").map_err(relabel))
        .transpose()?;
    let query_strand = merged_strand(
        query_frame,
        if legacy { None } else { get("query-strand")? },
    )?;
    let hit_strand = merged_strand(hit_frame, if legacy { None } else { get("hit-strand")? })?;
    fields::validate_direction(query_coordinates, query_strand, "query").map_err(relabel)?;
    fields::validate_direction(hit_coordinates, hit_strand, "subject").map_err(relabel)?;
    let integer = |name: &str| get(name)?.map(|value| uint(value, name, false)).transpose();
    let alignment_length = get("align-len")?
        .map(|value| uint(value, "align-len", true))
        .transpose()?;
    let pattern_from = get("pattern-from")?;
    let pattern_to = get("pattern-to")?;
    fields::coordinates(
        pattern_from.map(str::trim),
        pattern_to.map(str::trim),
        query_length,
        "pattern",
    )
    .map_err(relabel)?;
    if let Some(density) = get("density")? {
        uint(density, "density", false)?;
    }
    let hsp = BlastHsp {
        line: node.line,
        source_fields: flat_children(node, &[]),
        query_coordinates,
        hit_coordinates,
        query_frame,
        hit_frame,
        query_strand,
        hit_strand,
        evalue: Some(number(required("evalue")?, "evalue", true)?),
        bit_score: Some(number(required("bit-score")?, "bit-score", true)?),
        raw_score: Some(number(required("score")?, "score", false)?),
        alignment_length,
        identical: integer("identity")?,
        identity_percent: None,
        mismatches: None,
        positives: integer("positive")?,
        positive_percent: None,
        gap_openings: None,
        gaps: integer("gaps")?,
        query_alignment: Some(required("qseq")?.to_owned()),
        hit_alignment: Some(required("hseq")?.to_owned()),
        midline: get("midline")?.map(str::to_owned),
    };
    fields::validate_hsp(&hsp).map_err(relabel)?;
    if let Some(midline) = &hsp.midline {
        if !midline.is_ascii()
            || midline
                .bytes()
                .any(|byte| !byte.is_ascii_graphic() && byte != b' ')
            || hsp
                .query_alignment
                .as_ref()
                .is_some_and(|query| query.len() != midline.len())
        {
            return Err(error(
                IoErrorKind::LengthMismatch,
                "XML midline must preserve one printable ASCII character per alignment column",
            ));
        }
    }
    Ok(hsp)
}

fn check_hsp_fields(node: &Node, legacy: bool, prefix: &str) -> Result<(), BioIoError> {
    let common = [
        "num",
        "bit-score",
        "score",
        "evalue",
        "query-from",
        "query-to",
        "hit-from",
        "hit-to",
        "pattern-from",
        "pattern-to",
        "query-frame",
        "hit-frame",
        "identity",
        "positive",
        "gaps",
        "align-len",
        "density",
        "qseq",
        "hseq",
        "midline",
    ];
    let additional = ["query-strand", "hit-strand"];
    let allowed: Vec<_> = common
        .iter()
        .chain(additional.iter().filter(|_| !legacy))
        .map(|name| format!("{prefix}{name}"))
        .collect();
    let allowed: Vec<_> = allowed.iter().map(String::as_str).collect();
    node.check(&allowed, &[])
}

fn merged_strand(
    frame: Option<i8>,
    explicit: Option<&str>,
) -> Result<Option<BlastStrand>, BioIoError> {
    let inferred = fields::frame_strand(frame);
    let explicit = explicit
        .map(|value| fields::strand(value.trim()).map_err(relabel))
        .transpose()?
        .flatten();
    if explicit.is_some() && inferred.is_some() && explicit != inferred {
        return Err(error(
            IoErrorKind::InvalidStructure,
            "XML strand contradicts its explicit frame",
        ));
    }
    Ok(explicit.or(inferred))
}

pub(super) fn source_error(node: &Node) -> Result<BioIoError, BioIoError> {
    node.check(&["Err"], &[])?;
    let source = node.required("Err")?;
    source.check(&["code", "message"], &[])?;
    let code = signed(source.required_value("code")?, "BLAST error code")?;
    let message = source.leaf("message")?.unwrap_or("");
    Ok(error(
        IoErrorKind::InvalidStructure,
        format!("BLAST XML2 reports error {code}: {message}"),
    ))
}
