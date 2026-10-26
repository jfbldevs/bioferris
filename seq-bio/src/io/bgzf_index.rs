//! Verified, source-ordered FASTA/Sanger FASTQ indexes over stable BGZF files.
//!
//! Sidecars contain no path and are never trusted as proof of block boundaries.
//! Opening rederives every member and record span. Subsequent retrieval checks the
//! retained source handle, named path and selected decompressed raw-span digest.
//! Keep sources stable; these checks are not a filesystem snapshot or lock.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::Record;
use crate::publication::{PublishError, publish_new};

use super::bgzf::{BgzfError, BgzfReader, BgzfVirtualOffset, read_block};
use super::index::{
    IndexAlphabet, IndexError, IndexOptions, IndexedFormat, SourceState, absolute_path,
    open_regular, verify_source,
};
use super::{BioIoError, FastaReader, FastqReader, ReadLimits, RecordSpan, one_record};

#[path = "bgzf_index/error.rs"]
mod error;
pub use error::BgzfIndexError;
#[path = "bgzf_index/manifest.rs"]
mod manifest;
use manifest::{BlockEntry, Entry, Manifest, read_manifest, serialized_size};

const FORMAT_ID: &str = "bioferris-bgzf-sequence-index";
const KEY_POLICY: &str = "record-id-first-token-v1";
const EMPTY_DIGEST: &str = "0000000000000000000000000000000000000000000000000000000000000000";

/// Explicit inclusive budgets for compressed bytes, block metadata and records.
#[derive(Clone, Copy, Debug)]
pub struct BgzfIndexOptions {
    /// Existing FASTA/FASTQ and index budgets; total bytes are decompressed bytes.
    pub index: IndexOptions,
    /// Compressed source-byte limit; defaults to 1 GiB.
    pub max_compressed_bytes: u64,
    /// Maximum retained members, including empty members; defaults to 1,000,000.
    pub max_blocks: usize,
}

impl Default for BgzfIndexOptions {
    fn default() -> Self {
        Self {
            index: IndexOptions::default(),
            max_compressed_bytes: 1024 * 1024 * 1024,
            max_blocks: 1_000_000,
        }
    }
}

/// An immutable verified BGZF index retaining its source for checked retrieval.
#[derive(Debug)]
pub struct BgzfSequenceIndex<A: IndexAlphabet> {
    source: File,
    source_path: PathBuf,
    source_state: SourceState,
    manifest: Manifest,
    lookup: Vec<usize>,
    options: BgzfIndexOptions,
    alphabet: PhantomData<A>,
}

impl<A: IndexAlphabet> BgzfSequenceIndex<A> {
    /// Validates every BGZF member and record, without creating a sidecar.
    ///
    /// # Errors
    /// Rejects nonregular sources, malformed/ordinary gzip, missing final empty
    /// members, invalid records, duplicate IDs, exceeded budgets and source changes.
    pub fn build(
        source: impl AsRef<Path>,
        format: IndexedFormat,
        options: BgzfIndexOptions,
    ) -> Result<Self, BgzfIndexError> {
        let source_path = absolute_path(source.as_ref())?;
        let (mut source, source_state) = open_regular(&source_path)?;
        let result = derive_manifest::<A>(&mut source, &source_state, format, options);
        verify_source(&source, &source_path, &source_state)?;
        Self::verified(source, source_path, source_state, result?, options)
    }

    /// Fully revalidates a bounded sidecar against its explicitly supplied source.
    ///
    /// Opening is linear in source size, not metadata-only. No serialized path is
    /// followed. Every block, raw span, digest and canonical offset is rederived.
    ///
    /// # Errors
    /// Returns schema/forgery, source-integrity, BGZF, record, budget and I/O errors.
    pub fn open(
        source: impl AsRef<Path>,
        index_path: impl AsRef<Path>,
        options: BgzfIndexOptions,
    ) -> Result<Self, BgzfIndexError> {
        let manifest = read_manifest(index_path.as_ref(), options)?;
        manifest.validate::<A>(options)?;
        let source_path = absolute_path(source.as_ref())?;
        let (mut source, source_state) = open_regular(&source_path)?;
        if source_state.len != manifest.source_len {
            return Err(stale("compressed source length differs from the sidecar"));
        }
        let result =
            derive_manifest::<A>(&mut source, &source_state, manifest.source_format, options);
        verify_source(&source, &source_path, &source_state)?;
        let derived = result?;
        if derived.source_sha256 != manifest.source_sha256 {
            return Err(stale("compressed source SHA-256 differs from the sidecar"));
        }
        if derived != manifest {
            return Err(invalid(
                "sidecar differs from fully rederived blocks or records",
            ));
        }
        Self::verified(source, source_path, source_state, manifest, options)
    }

    fn verified(
        source: File,
        source_path: PathBuf,
        source_state: SourceState,
        manifest: Manifest,
        options: BgzfIndexOptions,
    ) -> Result<Self, BgzfIndexError> {
        let mut lookup = Vec::new();
        lookup
            .try_reserve_exact(manifest.entries.len())
            .map_err(allocation_error)?;
        lookup.extend(0..manifest.entries.len());
        lookup.sort_unstable_by(|a, b| manifest.entries[*a].id.cmp(&manifest.entries[*b].id));
        Ok(Self {
            source,
            source_path,
            source_state,
            manifest,
            lookup,
            options,
            alphabet: PhantomData,
        })
    }

    /// Publishes compact deterministic JSON without replacing any existing path.
    ///
    /// # Errors
    /// Returns source-change, existing-path, temporary-file and publication errors.
    /// The destination directory must be trusted; failed temporary files are removed.
    pub fn save_new(&self, path: impl AsRef<Path>) -> Result<(), BgzfIndexError> {
        self.verify_source()?;
        publish_new(path.as_ref(), |writer| {
            serde_json::to_writer(writer, &self.manifest).map_err(manifest::json_error)?;
            self.verify_source()
        })
        .map_err(|error| match error {
            PublishError::Write(error) => error,
            PublishError::Io(error) => error.into(),
            PublishError::DestinationExists(path) => IndexError::DestinationExists(path).into(),
            PublishError::TemporaryUnavailable => IndexError::TemporaryUnavailable.into(),
        })
    }

    /// Retrieves a preserved record by case-sensitive first-token ID.
    ///
    /// Missing IDs return `None` only after checking both source handle and path.
    /// No partially verified record is returned.
    ///
    /// # Errors
    /// Returns stale-source, I/O, decompression, raw-span or parsed-identity errors.
    pub fn get(&mut self, id: &str) -> Result<Option<Record<A>>, BgzfIndexError> {
        self.verify_source()?;
        let found = self
            .lookup
            .binary_search_by(|index| self.manifest.entries[*index].id.as_str().cmp(id));
        let Ok(position) = found else {
            self.verify_source()?;
            return Ok(None);
        };
        let entry = &self.manifest.entries[self.lookup[position]];
        let result = read_entry::<A>(
            &mut self.source,
            entry,
            &self.manifest,
            self.options.index.read_limits,
        );
        self.verify_source()?;
        result.map(Some)
    }

    /// Iterates keys in source order without copying.
    #[must_use]
    pub fn ids(&self) -> impl ExactSizeIterator<Item = &str> + DoubleEndedIterator {
        self.manifest.entries.iter().map(|entry| entry.id.as_str())
    }

    /// Returns the number of indexed records.
    #[must_use]
    pub fn len(&self) -> usize {
        self.manifest.entries.len()
    }

    /// Returns whether there are no records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.manifest.entries.is_empty()
    }

    /// Returns the explicit biological input profile.
    #[must_use]
    pub const fn format(&self) -> IndexedFormat {
        self.manifest.source_format
    }

    /// Returns the compressed source size in bytes.
    #[must_use]
    pub const fn source_len(&self) -> u64 {
        self.manifest.source_len
    }

    /// Returns the decompressed source size in bytes.
    #[must_use]
    pub const fn uncompressed_len(&self) -> u64 {
        self.manifest.uncompressed_len
    }

    /// Returns the lowercase SHA-256 of every compressed source byte.
    #[must_use]
    pub fn source_sha256(&self) -> &str {
        &self.manifest.source_sha256
    }

    /// Returns the number of validated members, including empty members.
    #[must_use]
    pub fn block_count(&self) -> usize {
        self.manifest.blocks.len()
    }

    /// Converts a decompressed byte position to its canonical verified virtual position.
    ///
    /// At a shared empty-member boundary the earliest physical block wins.
    ///
    /// # Errors
    /// Returns `InvalidIndex` if the position exceeds the decompressed source.
    pub fn virtual_offset(&self, offset: u64) -> Result<BgzfVirtualOffset, BgzfIndexError> {
        self.manifest.virtual_offset(offset)
    }

    /// Resolves a virtual position against the fully verified block table.
    ///
    /// Equivalent short-block ends and physical EOF are accepted, not arithmetic
    /// on packed offsets. This method does not perform I/O or recheck source state.
    ///
    /// # Errors
    /// Returns `InvalidIndex` for non-boundaries or beyond-block data positions.
    pub fn uncompressed_offset(&self, offset: BgzfVirtualOffset) -> Result<u64, BgzfIndexError> {
        self.manifest.uncompressed_offset(offset)
    }

    fn verify_source(&self) -> Result<(), BgzfIndexError> {
        verify_source(&self.source, &self.source_path, &self.source_state).map_err(Into::into)
    }
}

enum RecordReader<R: BufRead, A: IndexAlphabet> {
    Fasta(FastaReader<R, A>),
    Fastq(FastqReader<R, A>),
}

impl<R: BufRead, A: IndexAlphabet> RecordReader<R, A> {
    fn new(reader: R, format: IndexedFormat, limits: ReadLimits) -> Self {
        match format {
            IndexedFormat::Fasta => Self::Fasta(FastaReader::with_limits(reader, limits)),
            IndexedFormat::FastqSanger => Self::Fastq(FastqReader::with_limits(reader, limits)),
        }
    }
    fn span(&self) -> Option<RecordSpan> {
        match self {
            Self::Fasta(reader) => reader.record_span(),
            Self::Fastq(reader) => reader.record_span(),
        }
    }
}

impl<R: BufRead, A: IndexAlphabet> Iterator for RecordReader<R, A> {
    type Item = Result<Record<A>, BioIoError>;
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Fasta(reader) => reader.next(),
            Self::Fastq(reader) => reader.next(),
        }
    }
}

fn derive_manifest<A: IndexAlphabet>(
    file: &mut File,
    state: &SourceState,
    format: IndexedFormat,
    options: BgzfIndexOptions,
) -> Result<Manifest, BgzfIndexError> {
    let mut manifest = scan_blocks::<A>(file, state, format, options)?;
    derive_entries::<A>(file, &mut manifest, options)?;
    fingerprint_spans(file, &mut manifest)?;
    manifest.validate::<A>(options)?;
    Ok(manifest)
}

fn fingerprint_spans(file: &mut File, manifest: &mut Manifest) -> Result<(), BgzfIndexError> {
    file.seek(SeekFrom::Start(0))?;
    let mut decoder = BgzfReader::new(file);
    let mut offset = 0_u64;
    // One streaming pass hashes exact raw spans. Do not repeatedly decompress
    // the same member for each of its potentially thousands of short records.
    for entry in &mut manifest.entries {
        let gap = entry
            .offset
            .checked_sub(offset)
            .ok_or_else(|| invalid("raw record spans overlap"))?;
        if io::copy(&mut (&mut decoder).take(gap), &mut io::sink())? != gap {
            return Err(stale("source ended before a raw record span"));
        }
        let mut span = HashingReader::new((&mut decoder).take(entry.byte_length));
        io::copy(&mut span, &mut io::sink())?;
        if span.reader.limit() != 0 {
            return Err(stale("record raw span ended early"));
        }
        entry.sha256 = format!("{:x}", span.digest.finalize());
        offset = entry
            .offset
            .checked_add(entry.byte_length)
            .ok_or(IndexError::CountOverflow)?;
    }
    let tail = manifest
        .uncompressed_len
        .checked_sub(offset)
        .ok_or_else(|| invalid("raw record span exceeds decompressed size"))?;
    if io::copy(&mut (&mut decoder).take(tail), &mut io::sink())? != tail {
        return Err(stale(
            "decompressed source ended before its verified length",
        ));
    }
    if decoder.read_bytes(&mut [0])? != 0 {
        return Err(stale("decompressed source grew during fingerprinting"));
    }
    decoder.finish()?;
    Ok(())
}

fn scan_blocks<A: IndexAlphabet>(
    file: &mut File,
    state: &SourceState,
    format: IndexedFormat,
    options: BgzfIndexOptions,
) -> Result<Manifest, BgzfIndexError> {
    check_limit(
        "compressed_source_bytes",
        options.max_compressed_bytes,
        state.len,
    )?;
    file.seek(SeekFrom::Start(0))?;
    let mut hashing = HashingReader::new(&mut *file);
    let mut compressed = 0;
    let mut uncompressed = 0_u64;
    let mut blocks = Vec::new();
    while let Some(block) = read_block(&mut hashing, &mut compressed)? {
        check_limit(
            "bgzf_blocks",
            as_u64(options.max_blocks)?,
            as_u64(blocks.len())?
                .checked_add(1)
                .ok_or(IndexError::CountOverflow)?,
        )?;
        let next = uncompressed
            .checked_add(as_u64(block.data.len())?)
            .ok_or(IndexError::CountOverflow)?;
        check_limit(
            "source_bytes",
            as_u64(options.index.read_limits.max_total_bytes)?,
            next,
        )?;
        blocks.try_reserve(1).map_err(allocation_error)?;
        blocks.push(BlockEntry {
            compressed_offset: block.info.compressed_offset,
            compressed_size: u32::try_from(block.info.compressed_size)
                .map_err(|_| IndexError::CountOverflow)?,
            uncompressed_offset: uncompressed,
            uncompressed_size: u32::try_from(block.data.len())
                .map_err(|_| IndexError::CountOverflow)?,
            canonical_eof_marker: block.info.canonical_eof_marker,
        });
        uncompressed = next;
    }
    if compressed != state.len {
        return Err(stale("source size changed during member validation"));
    }
    if blocks
        .last()
        .is_none_or(|block| block.uncompressed_size != 0)
    {
        return Err(BgzfError::MissingEof { offset: compressed }.into());
    }
    let source_sha256 = format!("{:x}", hashing.digest.finalize());
    Ok(Manifest {
        format_id: FORMAT_ID.into(),
        schema_version: 1,
        molecule: A::MOLECULE.into(),
        source_format: format,
        key_policy: KEY_POLICY.into(),
        source_len: state.len,
        source_sha256,
        uncompressed_len: uncompressed,
        blocks,
        entries: Vec::new(),
    })
}

fn derive_entries<A: IndexAlphabet>(
    file: &mut File,
    manifest: &mut Manifest,
    options: BgzfIndexOptions,
) -> Result<(), BgzfIndexError> {
    let mut size = serialized_size(manifest)?;
    check_limit("index_bytes", as_u64(options.index.max_index_bytes)?, size)?;
    file.seek(SeekFrom::Start(0))?;
    let mut reader = RecordReader::<_, A>::new(
        BgzfReader::new(&mut *file),
        manifest.source_format,
        options.index.read_limits,
    );
    let mut seen = BTreeMap::<String, usize>::new();
    let mut id_bytes = 0_u64;
    while let Some(record) = reader.next() {
        let record = record?;
        let ordinal = manifest
            .entries
            .len()
            .checked_add(1)
            .ok_or(IndexError::CountOverflow)?;
        check_limit(
            "index_records",
            as_u64(options.index.max_records)?,
            as_u64(ordinal)?,
        )?;
        id_bytes = id_bytes
            .checked_add(as_u64(record.id().len())?)
            .ok_or(IndexError::CountOverflow)?;
        check_limit(
            "index_id_bytes",
            as_u64(options.index.max_id_bytes)?,
            id_bytes,
        )?;
        if let Some(first) = seen.get(record.id()) {
            return Err(IndexError::DuplicateId {
                id: record.id().into(),
                first_record: *first,
                duplicate_record: ordinal,
            }
            .into());
        }
        let span = reader
            .span()
            .ok_or_else(|| invalid("reader returned a record without raw coordinates"))?;
        let byte_length = span
            .end
            .checked_sub(span.start)
            .ok_or(IndexError::CountOverflow)?;
        if byte_length == 0 || span.end > manifest.uncompressed_len {
            return Err(stale("record extends beyond decompressed source"));
        }
        let entry = Entry {
            id: record.id().into(),
            offset: span.start,
            byte_length,
            virtual_start: manifest.virtual_offset(span.start)?.raw(),
            virtual_end: manifest.virtual_offset(span.end)?.raw(),
            sequence_length: as_u64(record.sequence().len())?,
            record: ordinal,
            header_line: span.header_line,
            sha256: EMPTY_DIGEST.into(),
        };
        size = size
            .checked_add(serialized_size(&entry)?)
            .and_then(|size| size.checked_add(u64::from(!manifest.entries.is_empty())))
            .ok_or(IndexError::CountOverflow)?;
        check_limit("index_bytes", as_u64(options.index.max_index_bytes)?, size)?;
        manifest.entries.try_reserve(1).map_err(allocation_error)?;
        seen.insert(entry.id.clone(), ordinal);
        manifest.entries.push(entry);
    }
    Ok(())
}

fn read_entry<A: IndexAlphabet>(
    file: &mut File,
    entry: &Entry,
    manifest: &Manifest,
    limits: ReadLimits,
) -> Result<Record<A>, BgzfIndexError> {
    let mut decoder = BgzfReader::new(file);
    decoder.seek_virtual(BgzfVirtualOffset::from_raw(entry.virtual_start))?;
    let mut hashing = HashingReader::new((&mut decoder).take(entry.byte_length));
    let result = one_record(RecordReader::<_, A>::new(
        BufReader::new(&mut hashing),
        manifest.source_format,
        limits,
    ));
    io::copy(&mut hashing, &mut io::sink())?;
    if hashing.reader.limit() != 0 || format!("{:x}", hashing.digest.finalize()) != entry.sha256 {
        return Err(stale(
            "decompressed record span differs from its verified SHA-256",
        ));
    }
    let end = decoder.virtual_position()?;
    let logical_end = manifest.uncompressed_offset(end)?;
    if logical_end
        != entry
            .offset
            .checked_add(entry.byte_length)
            .ok_or(IndexError::CountOverflow)?
    {
        return Err(invalid(
            "retrieved span ended at a different logical position",
        ));
    }
    let record = result
        .map_err(|error| invalid(format!("verified raw span failed record parsing: {error}")))?;
    if record.id() != entry.id || as_u64(record.sequence().len())? != entry.sequence_length {
        return Err(invalid(
            "retrieved record identity differs from its index entry",
        ));
    }
    Ok(record)
}

struct HashingReader<R> {
    reader: R,
    digest: Sha256,
}
impl<R> HashingReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            digest: Sha256::new(),
        }
    }
}
impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        let count = self.reader.read(output)?;
        self.digest.update(&output[..count]);
        Ok(count)
    }
}

fn check_limit(resource: &'static str, limit: u64, observed: u64) -> Result<(), BgzfIndexError> {
    if observed > limit {
        Err(IndexError::LimitExceeded {
            resource,
            limit,
            observed,
        }
        .into())
    } else {
        Ok(())
    }
}
fn as_u64(value: usize) -> Result<u64, BgzfIndexError> {
    u64::try_from(value).map_err(|_| IndexError::CountOverflow.into())
}
fn allocation_error(error: std::collections::TryReserveError) -> BgzfIndexError {
    io::Error::other(error).into()
}
fn stale(message: impl Into<String>) -> BgzfIndexError {
    IndexError::StaleSource(message.into()).into()
}
fn invalid(message: impl Into<String>) -> BgzfIndexError {
    IndexError::InvalidIndex(message.into()).into()
}

#[cfg(test)]
#[path = "bgzf_index/tests.rs"]
mod tests;
