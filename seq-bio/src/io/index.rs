//! Verified, source-ordered indexes for stable plain FASTA and Sanger FASTQ files.
//!
//! Building validates the whole source, and so does opening by default: it also
//! rederives every record boundary instead of trusting saved offsets. An opt-in
//! quick check compares the source length and modification time and verifies
//! sampled records instead. Sidecars are JSON (schema 1) or binary (schema 2);
//! an index opened from a binary sidecar reads record locations from it on
//! demand. Lookups seek on a retained file and verify the selected raw span
//! before returning the existing [`Record`] representation. Source files must
//! remain unchanged; metadata and digest checks do not provide a filesystem
//! snapshot or protect against every hostile concurrent mutation.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::fmt;
use std::fs::{File, Metadata};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use align_core::{Alphabet, Dna, Protein, Rna};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::Record;
use crate::publication::{PublishError, publish_new};

use super::{
    BioIoError, FastaReader, FastqReader, IoErrorKind, ReadLimits, RecordSpan, one_record,
};

mod binary;

const FORMAT_ID: &str = "bioferris-sequence-index";
const SCHEMA_VERSION: u32 = 1;
const KEY_POLICY: &str = "record-id-first-token-v1";
const BUFFER_BYTES: usize = 16 * 1024;
const EMPTY_DIGEST: &str = "0000000000000000000000000000000000000000000000000000000000000000";

mod sealed {
    pub trait Sealed {}

    impl Sealed for align_core::Dna {}
    impl Sealed for align_core::Rna {}
    impl Sealed for align_core::Protein {}
}

/// The three explicit molecule types supported by persisted sequence indexes.
///
/// This trait is sealed so a custom alphabet cannot silently reuse a built-in
/// molecule identifier in a persisted index.
pub trait IndexAlphabet: Alphabet + sealed::Sealed {
    /// Stable lowercase molecule identifier in the index schema.
    const MOLECULE: &'static str;
}

impl IndexAlphabet for Dna {
    const MOLECULE: &'static str = "dna";
}
impl IndexAlphabet for Rna {
    const MOLECULE: &'static str = "rna";
}
impl IndexAlphabet for Protein {
    const MOLECULE: &'static str = "protein";
}

/// The explicit, uncompressed source profiles supported by a sequence index.
#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
pub enum IndexedFormat {
    /// Metadata-preserving FASTA, including empty records and wrapped lines.
    #[serde(rename = "fasta")]
    Fasta,
    /// Sanger FASTQ with Phred+33 qualities and length-delimited wrapped records.
    #[serde(rename = "fastq-sanger")]
    FastqSanger,
}

impl IndexedFormat {
    /// Returns the stable source-profile identifier.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Fasta => "fasta",
            Self::FastqSanger => "fastq-sanger",
        }
    }
}

/// The sidecar encodings [`SequenceIndex::save_new_as`] can write.
///
/// [`SequenceIndex::open`] recognizes either from the file's first bytes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SidecarFormat {
    /// The schema-1 JSON object, which opening reads whole into memory.
    Json,
    /// The schema-2 binary layout, whose record locations stay on disk after
    /// opening; it also records the source modification time.
    Binary,
}

/// How much of the source [`SequenceIndex::open_with`] checks against a sidecar.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SourceCheck {
    /// Hash the whole source and rederive every record, as [`SequenceIndex::open`] does.
    #[default]
    Full,
    /// Compare the source length and modification time with those saved in a
    /// binary sidecar and verify up to 32 sampled records. It does not read the
    /// rest of the source, so it cannot detect every change.
    Quick,
}

/// Reader and index budgets, checked before retaining new index entries.
///
/// Zero limits are valid. Limits are inclusive. ID bytes count UTF-8 bytes of
/// stored keys, not headers. Building counts the complete compact JSON sidecar
/// without a trailing newline, which is always larger than the binary one;
/// opening counts the sidecar file and saving the bytes written. These are not
/// exact peak-memory limits.
#[derive(Clone, Copy, Debug)]
pub struct IndexOptions {
    /// Existing physical-line, record-symbol and total source-byte limits.
    pub read_limits: ReadLimits,
    /// Maximum indexed records; defaults to 100,000.
    pub max_records: usize,
    /// Maximum aggregate UTF-8 key bytes; defaults to 16 MiB.
    pub max_id_bytes: usize,
    /// Maximum serialized, loaded or written sidecar bytes; defaults to 64 MiB.
    pub max_index_bytes: usize,
}

impl Default for IndexOptions {
    fn default() -> Self {
        Self {
            read_limits: ReadLimits::default(),
            max_records: 100_000,
            max_id_bytes: 16 * 1024 * 1024,
            max_index_bytes: 64 * 1024 * 1024,
        }
    }
}

/// Typed failures from index construction, verification, persistence or lookup.
#[derive(Debug)]
#[non_exhaustive]
pub enum IndexError {
    /// A filesystem or stream operation failed.
    Io(io::Error),
    /// An existing biological reader rejected input or exceeded a reader limit.
    Input(BioIoError),
    /// A source or sidecar is not a supported regular, uncompressed file.
    UnsupportedSource(String),
    /// The sidecar schema, metadata or rederived record boundaries are invalid.
    InvalidIndex(String),
    /// Distinct source records share the same case-sensitive first-token key.
    DuplicateId {
        /// The duplicated record identifier.
        id: String,
        /// One-based ordinal of its first occurrence.
        first_record: usize,
        /// One-based ordinal of the repeated occurrence.
        duplicate_record: usize,
    },
    /// Source metadata or observed bytes no longer match the verified source.
    StaleSource(String),
    /// An explicit inclusive resource budget would be exceeded.
    LimitExceeded {
        /// Stable resource name, such as `index_records` or `index_bytes`.
        resource: &'static str,
        /// Configured inclusive limit.
        limit: u64,
        /// Attempted size, including the entry or bytes that would exceed it.
        observed: u64,
    },
    /// A count or coordinate cannot be represented without integer overflow.
    CountOverflow,
    /// Create-only publication refused an existing destination.
    DestinationExists(PathBuf),
    /// A same-directory temporary publication file could not be reserved.
    TemporaryUnavailable,
}

impl IndexError {
    /// Returns the stable machine-readable error category.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            Self::Io(_) => "io_error",
            Self::Input(error) => error.code(),
            Self::UnsupportedSource(_) => "unsupported_index_source",
            Self::InvalidIndex(_) => "invalid_index",
            Self::DuplicateId { .. } => "duplicate_record_id",
            Self::StaleSource(_) => "stale_source",
            Self::LimitExceeded { .. } => "limit_exceeded",
            Self::CountOverflow => "count_overflow",
            Self::DestinationExists(_) | Self::TemporaryUnavailable => "output_error",
        }
    }
}

impl fmt::Display for IndexError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => error.fmt(formatter),
            Self::Input(error) => error.fmt(formatter),
            Self::UnsupportedSource(message) => {
                write!(formatter, "unsupported index source: {message}")
            }
            Self::InvalidIndex(message) => write!(formatter, "invalid sequence index: {message}"),
            Self::DuplicateId {
                id,
                first_record,
                duplicate_record,
            } => write!(
                formatter,
                "duplicate record ID {id:?} at record {duplicate_record}; first seen at record {first_record}"
            ),
            Self::StaleSource(message) => write!(formatter, "indexed source changed: {message}"),
            Self::LimitExceeded {
                resource,
                limit,
                observed,
            } => write!(
                formatter,
                "{resource} limit {limit} exceeded by attempted size {observed}"
            ),
            Self::CountOverflow => formatter.write_str("index count or coordinate overflow"),
            Self::DestinationExists(path) => write!(
                formatter,
                "destination already exists; output is create-only: {}",
                path.display()
            ),
            Self::TemporaryUnavailable => {
                formatter.write_str("could not reserve a temporary output file")
            }
        }
    }
}

impl std::error::Error for IndexError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            Self::Input(error) => Some(error),
            _ => None,
        }
    }
}

impl From<io::Error> for IndexError {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}

impl From<BioIoError> for IndexError {
    fn from(error: BioIoError) -> Self {
        Self::Input(error)
    }
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {
    format_id: String,
    schema_version: u32,
    molecule: String,
    source_format: IndexedFormat,
    key_policy: String,
    source_len: u64,
    source_sha256: String,
    #[serde(deserialize_with = "deserialize_entries")]
    entries: Vec<Entry>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Entry {
    id: String,
    offset: u64,
    byte_length: u64,
    sequence_length: u64,
    record: usize,
    header_line: usize,
    sha256: String,
}

// Serde's derived struct deserializer also accepts positional sequences. The
// persisted contract requires objects, including each entry, and must retain
// the derived visitor's duplicate/unknown-field checks rather than pass through
// a generic JSON value that would overwrite duplicate keys.
struct ObjectOnly<T>(T);

impl<'de, T: Deserialize<'de>> Deserialize<'de> for ObjectOnly<T> {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct ObjectVisitor<T>(PhantomData<T>);

        impl<'de, T: Deserialize<'de>> serde::de::Visitor<'de> for ObjectVisitor<T> {
            type Value = ObjectOnly<T>;

            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str("a sequence-index JSON object")
            }

            fn visit_map<M: serde::de::MapAccess<'de>>(
                self,
                map: M,
            ) -> Result<Self::Value, M::Error> {
                T::deserialize(serde::de::value::MapAccessDeserializer::new(map)).map(ObjectOnly)
            }
        }

        deserializer.deserialize_map(ObjectVisitor(PhantomData))
    }
}

fn deserialize_entries<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<Vec<Entry>, D::Error> {
    Vec::<ObjectOnly<Entry>>::deserialize(deserializer)
        .map(|entries| entries.into_iter().map(|entry| entry.0).collect())
}

/// A verified, single-source index with source-ordered identifiers.
///
/// Keys are [`Record::id`] values, not complete headers. A sidecar is either the
/// deterministic `BioFerris` JSON object (schema 1) or the binary layout (schema
/// 2), not a Biopython `SQLite` index; [`Self::open`] recognizes both. Paths are
/// supplied explicitly and never serialized; non-UTF-8 Rust paths are preserved.
///
/// Construction and verified opening are linear in source size. By default,
/// opening is not a fast metadata-only operation: it validates every source
/// record and compares all saved metadata to newly derived entries.
/// [`Self::open_with`] with [`SourceCheck::Quick`] instead compares the source
/// length and modification time and verifies sampled records, which cannot
/// detect every change. Later lookups verify only the selected span plus source
/// metadata, not the whole file. An index opened from a binary sidecar holds only
/// the identifiers in memory and reads each record's location from the sidecar,
/// which must then stay unchanged as well. Keep the source stable for the index's
/// lifetime; no filesystem snapshot or mutation lock is implied.
#[derive(Debug)]
pub struct SequenceIndex<A: IndexAlphabet> {
    source: File,
    source_path: PathBuf,
    source_state: SourceState,
    format: IndexedFormat,
    source_len: u64,
    source_sha256: String,
    keys: Keys,
    locations: Store,
    options: IndexOptions,
    alphabet: PhantomData<A>,
}

impl<A: IndexAlphabet> SequenceIndex<A> {
    /// Validates and fingerprints a regular, uncompressed source file.
    ///
    /// Every record is parsed with the existing bounded FASTA/FASTQ readers.
    /// No sidecar is created until [`Self::save_new`] or [`Self::save_new_as`] is
    /// called.
    ///
    /// # Errors
    /// Rejects compressed/special files, malformed input, duplicate IDs, exceeded
    /// budgets, I/O failures and detected changes during validation.
    pub fn build(
        source: impl AsRef<Path>,
        format: IndexedFormat,
        options: IndexOptions,
    ) -> Result<Self, IndexError> {
        let source_path = absolute_path(source.as_ref())?;
        let (mut source, source_state) = open_regular(&source_path)?;
        check_limit(
            "source_bytes",
            options.read_limits.max_total_bytes,
            source_state.len,
        )?;
        let result = derive_entries::<A>(&mut source, &source_state, format, options);
        verify_source(&source, &source_path, &source_state)?;
        let derived = result?;
        Self::from_derived(source, source_path, source_state, format, derived, options)
    }

    /// Opens a sidecar and fully verifies its explicitly supplied source.
    ///
    /// This is [`Self::open_with`] with [`SourceCheck::Full`]. The sidecar may be
    /// JSON or binary. It is strictly validated within its byte budget, then the
    /// source digest and every rederived record boundary are checked. This reads
    /// the whole source; only subsequent [`Self::get`] calls use direct record seeks.
    ///
    /// # Errors
    /// Rejects malformed/forged indexes, wrong molecules, duplicate keys, exceeded
    /// budgets, unsupported files, I/O failures and stale or changing sources.
    pub fn open(
        source: impl AsRef<Path>,
        index_path: impl AsRef<Path>,
        options: IndexOptions,
    ) -> Result<Self, IndexError> {
        Self::open_with(source, index_path, options, SourceCheck::Full)
    }

    /// Opens a sidecar and checks its explicitly supplied source as selected.
    ///
    /// [`SourceCheck::Full`] behaves as [`Self::open`]. [`SourceCheck::Quick`]
    /// needs a binary sidecar. It validates the whole sidecar as a full check
    /// does, then requires the source's length and modification time to equal the
    /// saved ones and verifies up to 32 sampled records as [`Self::get`] would,
    /// without reading the rest of the source. It cannot detect a change that
    /// keeps the length and keeps or restores the modification time outside the
    /// sampled records: lookups still refuse any changed record they read, but
    /// [`Self::len`], [`Self::ids`] and [`Self::source_sha256`] then describe the
    /// indexed bytes rather than the current file.
    ///
    /// # Errors
    /// As [`Self::open`]. A quick check also rejects a JSON sidecar as
    /// `invalid_index` and a different modification time as `stale_source`.
    pub fn open_with(
        source: impl AsRef<Path>,
        index_path: impl AsRef<Path>,
        options: IndexOptions,
        check: SourceCheck,
    ) -> Result<Self, IndexError> {
        let index_path = index_path.as_ref();
        if binary::is_binary(index_path)? {
            return Self::open_binary(source.as_ref(), index_path, options, check);
        }
        if check == SourceCheck::Quick {
            return Err(IndexError::InvalidIndex(
                "quick source checks need a binary sidecar, which records the source modification time"
                    .into(),
            ));
        }
        Self::open_json(source.as_ref(), index_path, options)
    }

    fn open_json(
        source: &Path,
        index_path: &Path,
        options: IndexOptions,
    ) -> Result<Self, IndexError> {
        let manifest = read_manifest(index_path, options)?;
        validate_manifest::<A>(&manifest, options)?;
        let source_path = absolute_path(source)?;
        let (mut source, source_state) = open_regular(&source_path)?;
        ensure_plain(&mut source)?;
        if source_state.len != manifest.source_len {
            return Err(IndexError::StaleSource(
                "source byte length differs from the index".into(),
            ));
        }
        check_limit(
            "source_bytes",
            options.read_limits.max_total_bytes,
            source_state.len,
        )?;
        let observed_hash = fingerprint(&mut source, source_state.len, &mut [])?;
        verify_source(&source, &source_path, &source_state)?;
        if observed_hash != manifest.source_sha256 {
            return Err(IndexError::StaleSource(
                "source SHA-256 differs from the index".into(),
            ));
        }
        let format = manifest.source_format;
        let derived = derive_entries::<A>(&mut source, &source_state, format, options);
        verify_source(&source, &source_path, &source_state)?;
        let derived = derived.map_err(profile_error)?;
        if derived.source_sha256 != manifest.source_sha256 {
            return Err(IndexError::StaleSource(
                "source bytes changed during index verification".into(),
            ));
        }
        if !derived.matches(&manifest) {
            return Err(IndexError::InvalidIndex(
                "saved entries differ from the fully rederived source records".into(),
            ));
        }
        drop(manifest);
        Self::from_derived(source, source_path, source_state, format, derived, options)
    }

    fn open_binary(
        source: &Path,
        index_path: &Path,
        options: IndexOptions,
        check: SourceCheck,
    ) -> Result<Self, IndexError> {
        let loaded = binary::load::<A>(index_path, options)?;
        let source_path = absolute_path(source)?;
        let (mut source, source_state) = open_regular(&source_path)?;
        ensure_plain(&mut source)?;
        if source_state.len != loaded.header.source_len {
            return Err(IndexError::StaleSource(
                "source byte length differs from the index".into(),
            ));
        }
        check_limit(
            "source_bytes",
            options.read_limits.max_total_bytes,
            source_state.len,
        )?;
        match check {
            SourceCheck::Full => binary::verify_full::<A>(
                &mut source,
                &source_path,
                &source_state,
                &loaded,
                options,
            )?,
            SourceCheck::Quick => {
                binary::verify_quick::<A>(&mut source, &source_state, &loaded, options)?;
            }
        }
        verify_source(&source, &source_path, &source_state)?;
        loaded.sidecar.verify()?;
        let binary::Loaded {
            header,
            keys,
            sidecar,
        } = loaded;
        Ok(Self {
            source,
            source_path,
            source_state,
            format: header.format,
            source_len: header.source_len,
            source_sha256: hex(&header.source_sha256),
            keys,
            locations: Store::Binary(sidecar),
            options,
            alphabet: PhantomData,
        })
    }

    fn from_derived(
        source: File,
        source_path: PathBuf,
        source_state: SourceState,
        format: IndexedFormat,
        derived: Derived,
        options: IndexOptions,
    ) -> Result<Self, IndexError> {
        let Derived {
            text,
            ends,
            locations,
            source_sha256,
        } = derived;
        Ok(Self {
            source,
            source_path,
            source_len: source_state.len,
            source_state,
            format,
            source_sha256,
            keys: Keys::new(text, ends)?,
            locations: Store::Memory(locations),
            options,
            alphabet: PhantomData,
        })
    }

    /// Publishes the compact JSON sidecar without replacing an existing path.
    ///
    /// This is [`Self::save_new_as`] with [`SidecarFormat::Json`].
    ///
    /// # Errors
    /// As [`Self::save_new_as`].
    pub fn save_new(&self, path: impl AsRef<Path>) -> Result<(), IndexError> {
        self.save_new_as(path, SidecarFormat::Json)
    }

    /// Publishes a JSON or binary sidecar without replacing an existing path.
    ///
    /// Uses a same-directory temporary file and create-only hard-link publication.
    /// The destination directory must be trusted. A failed write leaves no named
    /// partial sidecar; publication does not change or copy the source file. Both
    /// encodings hold the same entries and reopen to the same records; a binary
    /// sidecar also records the source modification time for
    /// [`SourceCheck::Quick`].
    ///
    /// # Errors
    /// Returns existing-destination, temporary-file or I/O failures, refuses a
    /// sidecar larger than `max_index_bytes`, and rejects a source whose metadata
    /// changed before or during serialization.
    pub fn save_new_as(
        &self,
        path: impl AsRef<Path>,
        format: SidecarFormat,
    ) -> Result<(), IndexError> {
        self.verify_source()?;
        publish_new(path.as_ref(), |writer| {
            let mut output = BudgetWriter {
                writer: BufWriter::new(writer),
                written: 0,
                limit: as_u64(self.options.max_index_bytes)?,
                failure: None,
            };
            let written = match format {
                SidecarFormat::Json => self.write_json(&mut output),
                SidecarFormat::Binary => self.write_binary(&mut output),
            }
            .and_then(|()| output.flush().map_err(IndexError::from));
            if let Some(error) = output.failure.take() {
                return Err(error);
            }
            written?;
            self.verify_source()
        })
        .map_err(|error| match error {
            PublishError::Io(error) => IndexError::Io(error),
            PublishError::DestinationExists(path) => IndexError::DestinationExists(path),
            PublishError::TemporaryUnavailable => IndexError::TemporaryUnavailable,
            PublishError::Write(error) => error,
        })
    }

    fn write_json(&self, writer: &mut dyn Write) -> Result<(), IndexError> {
        let entries = JsonEntries {
            index: self,
            failure: RefCell::new(None),
        };
        let manifest = JsonManifest {
            format_id: FORMAT_ID,
            schema_version: SCHEMA_VERSION,
            molecule: A::MOLECULE,
            source_format: self.format,
            key_policy: KEY_POLICY,
            source_len: self.source_len,
            source_sha256: &self.source_sha256,
            entries: &entries,
        };
        serde_json::to_writer(writer, &manifest)
            .map_err(|error| entries.failure.take().unwrap_or_else(|| json_error(error)))
    }

    fn write_binary(&self, writer: &mut dyn Write) -> Result<(), IndexError> {
        let sha256 = decode_digest(&self.source_sha256)?;
        let source = binary::Source {
            format: self.format,
            len: self.source_len,
            modified: self.source_state.modified,
            sha256: &sha256,
        };
        let mut locations = self.locations.cursor(self.keys.len())?;
        binary::write::<A>(writer, &source, &self.keys, || locations.next_location())
    }

    /// Retrieves one preserved record by its case-sensitive first-token ID.
    ///
    /// A missing ID returns `None`. Metadata is checked even for a missing ID.
    /// Present records are parsed through the existing reader after seeking to
    /// their bounded span; raw SHA-256, parsed ID and sequence length must match.
    /// An index opened from a binary sidecar first reads the record's location
    /// from the sidecar, whose metadata is checked the same way.
    ///
    /// # Errors
    /// Returns I/O failures or a stale-source error on detected changes. No
    /// partial record is returned after a digest, parse or metadata failure.
    pub fn get(&mut self, id: &str) -> Result<Option<Record<A>>, IndexError> {
        self.verify_source()?;
        let Some(position) = self.keys.find(id) else {
            self.verify_source()?;
            return Ok(None);
        };
        let result = self.locations.location(position).and_then(|location| {
            read_entry::<A>(
                &mut self.source,
                self.keys.id(position),
                &location,
                self.format,
                self.options.read_limits,
            )
        });
        self.verify_source()?;
        result.map(Some)
    }

    /// Iterates identifiers in original source order without copying them.
    #[must_use]
    pub fn ids(&self) -> impl ExactSizeIterator<Item = &str> + DoubleEndedIterator {
        (0..self.keys.len()).map(|position| self.keys.id(position))
    }

    /// Returns the number of indexed records.
    #[must_use]
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    /// Returns whether the source contains no records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.keys.len() == 0
    }

    /// Returns the explicitly selected source profile.
    #[must_use]
    pub const fn format(&self) -> IndexedFormat {
        self.format
    }

    /// Returns the verified raw source byte length, including line endings.
    #[must_use]
    pub const fn source_len(&self) -> u64 {
        self.source_len
    }

    /// Returns the lowercase hexadecimal SHA-256 of all source bytes.
    ///
    /// Building and full checks verify it; after a quick check it is the saved
    /// value.
    #[must_use]
    pub fn source_sha256(&self) -> &str {
        &self.source_sha256
    }

    fn verify_source(&self) -> Result<(), IndexError> {
        verify_source(&self.source, &self.source_path, &self.source_state)?;
        self.locations.verify()
    }
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) struct SourceState {
    pub(crate) len: u64,
    modified: SystemTime,
    created: Option<SystemTime>,
    #[cfg(unix)]
    device: u64,
    #[cfg(unix)]
    inode: u64,
    #[cfg(unix)]
    change_seconds: i64,
    #[cfg(unix)]
    change_nanoseconds: i64,
}

impl SourceState {
    fn from_metadata(metadata: &Metadata) -> Result<Self, IndexError> {
        #[cfg(unix)]
        use std::os::unix::fs::MetadataExt;

        Ok(Self {
            len: metadata.len(),
            modified: metadata.modified()?,
            created: metadata.created().ok(),
            #[cfg(unix)]
            device: metadata.dev(),
            #[cfg(unix)]
            inode: metadata.ino(),
            #[cfg(unix)]
            change_seconds: metadata.ctime(),
            #[cfg(unix)]
            change_nanoseconds: metadata.ctime_nsec(),
        })
    }
}

pub(crate) fn absolute_path(path: &Path) -> Result<PathBuf, IndexError> {
    if path.is_absolute() {
        Ok(path.to_path_buf())
    } else {
        Ok(std::env::current_dir()?.join(path))
    }
}

pub(crate) fn open_regular(path: &Path) -> Result<(File, SourceState), IndexError> {
    // Inspect the type before opening so an ordinary FIFO cannot block here.
    // A hostile path replacement between these operations is not a snapshot.
    let before = std::fs::metadata(path)?;
    if !before.is_file() {
        return Err(IndexError::UnsupportedSource(
            "a regular file is required; directories, devices and pipes are unsupported".into(),
        ));
    }
    let state = SourceState::from_metadata(&before)?;
    let file = File::open(path)?;
    let opened = file.metadata()?;
    if !opened.is_file() || SourceState::from_metadata(&opened)? != state {
        return Err(IndexError::StaleSource(
            "file identity changed while opening".into(),
        ));
    }
    Ok((file, state))
}

pub(crate) fn verify_source(
    file: &File,
    path: &Path,
    expected: &SourceState,
) -> Result<(), IndexError> {
    check_unchanged(file, path, expected, "source", IndexError::StaleSource)
}

/// Rejects a retained file whose handle or path no longer matches `expected`.
fn check_unchanged(
    file: &File,
    path: &Path,
    expected: &SourceState,
    role: &str,
    changed: fn(String) -> IndexError,
) -> Result<(), IndexError> {
    let handle = file
        .metadata()
        .map_err(|error| changed(format!("cannot inspect retained {role}: {error}")))?;
    let named = std::fs::metadata(path)
        .map_err(|error| changed(format!("cannot inspect {role} path: {error}")))?;
    if !handle.is_file()
        || !named.is_file()
        || SourceState::from_metadata(&handle)? != *expected
        || SourceState::from_metadata(&named)? != *expected
    {
        return Err(changed(format!(
            "{role} identity, size or timestamps differ from the verified file"
        )));
    }
    Ok(())
}

fn ensure_plain(file: &mut File) -> Result<(), IndexError> {
    file.seek(SeekFrom::Start(0))?;
    let mut signature = [0; 2];
    let mut read = 0usize;
    while read < signature.len() {
        match file.read(&mut signature[read..]) {
            Ok(0) => break,
            Ok(count) => read = read.checked_add(count).ok_or(IndexError::CountOverflow)?,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error.into()),
        }
    }
    file.seek(SeekFrom::Start(0))?;
    if read == 2 && signature == [0x1f, 0x8b] {
        return Err(IndexError::UnsupportedSource(
            "gzip and BGZF sources cannot use plain-file byte offsets".into(),
        ));
    }
    Ok(())
}

enum RecordReader<R: BufRead, A: Alphabet> {
    Fasta(FastaReader<R, A>),
    Fastq(FastqReader<R, A>),
}

impl<R: BufRead, A: Alphabet> RecordReader<R, A> {
    fn new(reader: R, format: IndexedFormat, limits: ReadLimits) -> Self {
        match format {
            IndexedFormat::Fasta => Self::Fasta(FastaReader::with_limits(reader, limits)),
            IndexedFormat::FastqSanger => Self::Fastq(FastqReader::with_limits(reader, limits)),
        }
    }

    fn record_span(&self) -> Option<RecordSpan> {
        match self {
            Self::Fasta(reader) => reader.record_span(),
            Self::Fastq(reader) => reader.record_span(),
        }
    }
}

impl<R: BufRead, A: Alphabet> Iterator for RecordReader<R, A> {
    type Item = Result<Record<A>, BioIoError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Fasta(reader) => reader.next(),
            Self::Fastq(reader) => reader.next(),
        }
    }
}

/// Identifiers and locations derived from every source record, with the digest
/// of the whole source.
struct Derived {
    text: String,
    ends: Vec<usize>,
    locations: Vec<Location>,
    source_sha256: String,
}

impl Derived {
    /// Whether a validated JSON manifest holds exactly these entries.
    fn matches(&self, manifest: &Manifest) -> bool {
        manifest.source_sha256 == self.source_sha256
            && manifest.entries.len() == self.locations.len()
            && manifest
                .entries
                .iter()
                .zip(&self.locations)
                .enumerate()
                .all(|(position, (entry, location))| {
                    entry.id == id_at(&self.text, &self.ends, position)
                        && entry.offset == location.offset
                        && entry.byte_length == location.byte_length
                        && entry.sequence_length == location.sequence_length
                        && Some(entry.record) == position.checked_add(1)
                        && as_u64(entry.header_line).ok() == Some(location.header_line)
                        && decode_digest(&entry.sha256).ok() == Some(location.sha256)
                })
    }
}

fn derive_entries<A: IndexAlphabet>(
    file: &mut File,
    state: &SourceState,
    format: IndexedFormat,
    options: IndexOptions,
) -> Result<Derived, IndexError> {
    ensure_plain(file)?;
    // Budget the JSON form, which is always larger than the binary one.
    let mut serialized_bytes = serialized_size(&JsonManifest {
        format_id: FORMAT_ID,
        schema_version: SCHEMA_VERSION,
        molecule: A::MOLECULE,
        source_format: format,
        key_policy: KEY_POLICY,
        source_len: state.len,
        source_sha256: EMPTY_DIGEST,
        entries: &[] as &[JsonEntry<'_>],
    })?;
    check_limit("index_bytes", options.max_index_bytes, serialized_bytes)?;
    let mut ids = BTreeMap::<String, usize>::new();
    let mut id_bytes = 0u64;
    let mut text = String::new();
    let mut ends = Vec::new();
    let mut locations = Vec::<Location>::new();
    {
        let mut reader =
            RecordReader::<_, A>::new(BufReader::new(&mut *file), format, options.read_limits);
        while let Some(record) = reader.next() {
            let record = record?;
            let ordinal = locations
                .len()
                .checked_add(1)
                .ok_or(IndexError::CountOverflow)?;
            check_limit("index_records", options.max_records, as_u64(ordinal)?)?;
            id_bytes = id_bytes
                .checked_add(as_u64(record.id().len())?)
                .ok_or(IndexError::CountOverflow)?;
            check_limit("index_id_bytes", options.max_id_bytes, id_bytes)?;
            if let Some(first_record) = ids.get(record.id()) {
                return Err(IndexError::DuplicateId {
                    id: record.id().into(),
                    first_record: *first_record,
                    duplicate_record: ordinal,
                });
            }
            let span = reader.record_span().ok_or_else(|| {
                IndexError::InvalidIndex("reader returned a record without a source span".into())
            })?;
            let byte_length = span
                .end
                .checked_sub(span.start)
                .ok_or(IndexError::CountOverflow)?;
            if byte_length == 0 || span.end > state.len {
                return Err(IndexError::StaleSource(
                    "parsed record extends outside the original source length".into(),
                ));
            }
            let location = Location {
                offset: span.start,
                byte_length,
                sequence_length: as_u64(record.sequence().len())?,
                header_line: as_u64(span.header_line)?,
                sha256: [0; 32],
            };
            let separator_bytes = u64::from(!locations.is_empty());
            serialized_bytes = serialized_bytes
                .checked_add(serialized_size(&JsonEntry {
                    id: record.id(),
                    offset: location.offset,
                    byte_length,
                    sequence_length: location.sequence_length,
                    record: ordinal,
                    header_line: location.header_line,
                    sha256: EMPTY_DIGEST,
                })?)
                .and_then(|count| count.checked_add(separator_bytes))
                .ok_or(IndexError::CountOverflow)?;
            check_limit("index_bytes", options.max_index_bytes, serialized_bytes)?;
            let terminated = record
                .id()
                .len()
                .checked_add(1)
                .ok_or(IndexError::CountOverflow)?;
            text.try_reserve(terminated).map_err(allocation_error)?;
            ends.try_reserve(1).map_err(allocation_error)?;
            locations.try_reserve(1).map_err(allocation_error)?;
            ids.insert(record.id().into(), ordinal);
            text.push_str(record.id());
            ends.push(text.len());
            text.push('\n');
            locations.push(location);
        }
    }
    // Temporary duplicate detection owns no source sequence data and is released
    // before fingerprinting. The final lookup table stores positions, not ID copies.
    drop(ids);
    let source_sha256 = fingerprint(file, state.len, &mut locations)?;
    Ok(Derived {
        text,
        ends,
        locations,
        source_sha256,
    })
}

fn fingerprint(
    file: &mut File,
    source_len: u64,
    locations: &mut [Location],
) -> Result<String, IndexError> {
    file.seek(SeekFrom::Start(0))?;
    let mut source_digest = Sha256::new();
    let mut offset = 0u64;
    let mut buffer = [0u8; BUFFER_BYTES];
    for location in locations {
        let gap = location.offset.checked_sub(offset).ok_or_else(|| {
            IndexError::InvalidIndex("record spans overlap or are out of order".into())
        })?;
        hash_section(file, gap, &mut source_digest, None, &mut buffer)?;
        let mut span_digest = Sha256::new();
        hash_section(
            file,
            location.byte_length,
            &mut source_digest,
            Some(&mut span_digest),
            &mut buffer,
        )?;
        location.sha256 = span_digest.finalize().into();
        offset = location
            .offset
            .checked_add(location.byte_length)
            .ok_or(IndexError::CountOverflow)?;
    }
    let tail = source_len
        .checked_sub(offset)
        .ok_or_else(|| IndexError::InvalidIndex("record span exceeds source length".into()))?;
    hash_section(file, tail, &mut source_digest, None, &mut buffer)?;
    if !at_eof(file)? {
        return Err(IndexError::StaleSource(
            "source grew during fingerprinting".into(),
        ));
    }
    Ok(format!("{:x}", source_digest.finalize()))
}

fn hash_section(
    file: &mut File,
    mut remaining: u64,
    source_digest: &mut Sha256,
    mut span_digest: Option<&mut Sha256>,
    buffer: &mut [u8],
) -> Result<(), IndexError> {
    while remaining > 0 {
        let count = usize::try_from(remaining.min(as_u64(buffer.len())?))
            .map_err(|_| IndexError::CountOverflow)?;
        let read = match file.read(&mut buffer[..count]) {
            Ok(0) => {
                return Err(IndexError::StaleSource(
                    "source ended before its verified byte length".into(),
                ));
            }
            Ok(count) => count,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error.into()),
        };
        source_digest.update(&buffer[..read]);
        if let Some(digest) = &mut span_digest {
            digest.update(&buffer[..read]);
        }
        remaining = remaining
            .checked_sub(as_u64(read)?)
            .ok_or(IndexError::CountOverflow)?;
    }
    Ok(())
}

fn read_manifest(path: &Path, options: IndexOptions) -> Result<Manifest, IndexError> {
    let (mut file, state) = open_regular(path)?;
    check_limit("index_bytes", options.max_index_bytes, state.len)?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(usize::try_from(state.len).map_err(|_| IndexError::CountOverflow)?)
        .map_err(allocation_error)?;
    (&mut file)
        .take(as_u64(options.max_index_bytes)?)
        .read_to_end(&mut bytes)?;
    if !at_eof(&mut file)? {
        return Err(IndexError::LimitExceeded {
            resource: "index_bytes",
            limit: as_u64(options.max_index_bytes)?,
            observed: as_u64(options.max_index_bytes)?
                .checked_add(1)
                .ok_or(IndexError::CountOverflow)?,
        });
    }
    verify_source(&file, path, &state).map_err(|error| {
        IndexError::InvalidIndex(format!("sidecar changed while reading: {error}"))
    })?;
    serde_json::from_slice::<ObjectOnly<Manifest>>(&bytes)
        .map(|manifest| manifest.0)
        .map_err(json_error)
}

fn at_eof(reader: &mut impl Read) -> io::Result<bool> {
    let mut byte = [0u8; 1];
    loop {
        match reader.read(&mut byte) {
            Ok(count) => return Ok(count == 0),
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(error),
        }
    }
}

fn validate_manifest<A: IndexAlphabet>(
    manifest: &Manifest,
    options: IndexOptions,
) -> Result<(), IndexError> {
    if manifest.format_id != FORMAT_ID
        || manifest.schema_version != SCHEMA_VERSION
        || manifest.key_policy != KEY_POLICY
    {
        return Err(IndexError::InvalidIndex(
            "unsupported format ID, schema version or key policy".into(),
        ));
    }
    if manifest.molecule != A::MOLECULE {
        return Err(IndexError::InvalidIndex(format!(
            "index molecule {:?} does not match {}",
            manifest.molecule,
            A::MOLECULE
        )));
    }
    validate_digest(&manifest.source_sha256)?;
    check_limit(
        "index_records",
        options.max_records,
        as_u64(manifest.entries.len())?,
    )?;
    let mut id_bytes = 0u64;
    let mut spans = SpanCheck::new(manifest.source_len, manifest.source_format);
    for (position, entry) in manifest.entries.iter().enumerate() {
        id_bytes = id_bytes
            .checked_add(as_u64(entry.id.len())?)
            .ok_or(IndexError::CountOverflow)?;
        check_limit("index_id_bytes", options.max_id_bytes, id_bytes)?;
        if !valid_id(&entry.id) {
            return Err(IndexError::InvalidIndex(
                "entry ID must be a nonempty first-token identifier".into(),
            ));
        }
        spans.check(
            entry.offset,
            entry.byte_length,
            entry.sequence_length,
            as_u64(entry.header_line)?,
        )?;
        let ordinal = position.checked_add(1).ok_or(IndexError::CountOverflow)?;
        if entry.record != ordinal {
            return Err(inconsistent_entry());
        }
        validate_digest(&entry.sha256)?;
    }
    // Entry/ID budgets have been validated before building any lookup structure.
    let mut seen = BTreeMap::new();
    for entry in &manifest.entries {
        if seen.insert(entry.id.as_str(), entry.record).is_some() {
            return Err(IndexError::InvalidIndex(
                "sidecar contains duplicate record IDs".into(),
            ));
        }
    }
    Ok(())
}

fn validate_digest(digest: &str) -> Result<(), IndexError> {
    if digest.len() != 64
        || !digest
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(IndexError::InvalidIndex(
            "SHA-256 must be exactly 64 lowercase hexadecimal characters".into(),
        ));
    }
    Ok(())
}

fn valid_id(id: &str) -> bool {
    !id.is_empty()
        && !id
            .chars()
            .any(|character| character.is_whitespace() || character.is_control())
}

/// The order and bounds every saved entry must respect, checked in source order.
struct SpanCheck {
    previous_end: u64,
    previous_line: u64,
    source_len: u64,
    format: IndexedFormat,
}

impl SpanCheck {
    const fn new(source_len: u64, format: IndexedFormat) -> Self {
        Self {
            previous_end: 0,
            previous_line: 0,
            source_len,
            format,
        }
    }

    fn check(
        &mut self,
        offset: u64,
        byte_length: u64,
        sequence_length: u64,
        header_line: u64,
    ) -> Result<(), IndexError> {
        let end = offset.checked_add(byte_length).ok_or_else(|| {
            IndexError::InvalidIndex("record span overflows its coordinate type".into())
        })?;
        if byte_length == 0
            || offset < self.previous_end
            || end > self.source_len
            || sequence_length > byte_length
            || header_line <= self.previous_line
            || (self.format == IndexedFormat::FastqSanger && sequence_length == 0)
        {
            return Err(inconsistent_entry());
        }
        self.previous_end = end;
        self.previous_line = header_line;
        Ok(())
    }
}

fn inconsistent_entry() -> IndexError {
    IndexError::InvalidIndex("entry coordinates, order or record metadata are inconsistent".into())
}

/// Reports a source that no longer parses under the saved profile as an invalid
/// index, keeping I/O and limit failures as they are.
fn profile_error(error: IndexError) -> IndexError {
    match error {
        IndexError::Input(input)
            if !matches!(input.kind(), IoErrorKind::Io | IoErrorKind::LimitExceeded) =>
        {
            IndexError::InvalidIndex(format!(
                "source does not match the declared profile: {input}"
            ))
        }
        other => other,
    }
}

struct HashingReader<R> {
    reader: R,
    digest: Sha256,
}

impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let count = self.reader.read(buffer)?;
        self.digest.update(&buffer[..count]);
        Ok(count)
    }
}

/// Identifiers in source order, each followed by a line feed, and their order.
#[derive(Debug)]
struct Keys {
    text: String,
    /// Byte offset of each identifier's terminating line feed.
    ends: Vec<usize>,
    /// Source positions in ascending identifier order.
    sorted: Vec<usize>,
}

impl Keys {
    /// Orders source positions by identifier for binary search.
    fn new(text: String, ends: Vec<usize>) -> Result<Self, IndexError> {
        let mut sorted = Vec::new();
        sorted
            .try_reserve_exact(ends.len())
            .map_err(allocation_error)?;
        sorted.extend(0..ends.len());
        sorted.sort_unstable_by(|first, second| {
            id_at(&text, &ends, *first).cmp(id_at(&text, &ends, *second))
        });
        Ok(Self { text, ends, sorted })
    }

    fn len(&self) -> usize {
        self.ends.len()
    }

    fn id(&self, position: usize) -> &str {
        id_at(&self.text, &self.ends, position)
    }

    fn find(&self, id: &str) -> Option<usize> {
        self.sorted
            .binary_search_by(|position| self.id(*position).cmp(id))
            .ok()
            .map(|found| self.sorted[found])
    }
}

/// The identifier at a source position in line-feed-terminated `text`.
fn id_at<'a>(text: &'a str, ends: &[usize], position: usize) -> &'a str {
    let start = position
        .checked_sub(1)
        .map_or(0, |previous| ends[previous] + 1);
    &text[start..ends[position]]
}

/// Where one record's raw span lies in the source, and its digest.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Location {
    offset: u64,
    byte_length: u64,
    sequence_length: u64,
    header_line: u64,
    sha256: [u8; 32],
}

/// Where an index keeps its record locations.
#[derive(Debug)]
enum Store {
    /// Built indexes and those opened from JSON sidecars.
    Memory(Vec<Location>),
    /// Indexes opened from binary sidecars read locations on demand.
    Binary(binary::Sidecar),
}

impl Store {
    fn location(&self, position: usize) -> Result<Location, IndexError> {
        match self {
            Self::Memory(locations) => Ok(locations[position].clone()),
            Self::Binary(sidecar) => sidecar.location(position),
        }
    }

    fn cursor(&self, records: usize) -> Result<Cursor<'_>, IndexError> {
        Ok(match self {
            Self::Memory(locations) => Cursor::Memory(locations.iter()),
            Self::Binary(sidecar) => Cursor::Binary(sidecar.locations(records)?),
        })
    }

    fn verify(&self) -> Result<(), IndexError> {
        match self {
            Self::Memory(_) => Ok(()),
            Self::Binary(sidecar) => sidecar.verify(),
        }
    }
}

/// Source-ordered locations from either store.
enum Cursor<'a> {
    Memory(std::slice::Iter<'a, Location>),
    Binary(binary::Locations<'a>),
}

impl Cursor<'_> {
    fn next_location(&mut self) -> Result<Option<Location>, IndexError> {
        match self {
            Self::Memory(locations) => Ok(locations.next().cloned()),
            Self::Binary(locations) => locations.next_location(),
        }
    }
}

/// The schema-1 JSON layout, with the same fields in the same order as [`Manifest`].
#[derive(Serialize)]
struct JsonManifest<'a, E> {
    format_id: &'static str,
    schema_version: u32,
    molecule: &'static str,
    source_format: IndexedFormat,
    key_policy: &'static str,
    source_len: u64,
    source_sha256: &'a str,
    entries: E,
}

/// One schema-1 entry, with the same fields in the same order as [`Entry`].
#[derive(Serialize)]
struct JsonEntry<'a> {
    id: &'a str,
    offset: u64,
    byte_length: u64,
    sequence_length: u64,
    record: usize,
    header_line: u64,
    sha256: &'a str,
}

/// Serializes an index's entries in source order, reading them as needed.
struct JsonEntries<'a, A: IndexAlphabet> {
    index: &'a SequenceIndex<A>,
    failure: RefCell<Option<IndexError>>,
}

impl<A: IndexAlphabet> Serialize for JsonEntries<'_, A> {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::{Error as _, SerializeSeq as _};

        let fail = |error: IndexError| {
            self.failure.replace(Some(error));
            S::Error::custom("could not read an index entry")
        };
        let keys = &self.index.keys;
        let mut locations = self.index.locations.cursor(keys.len()).map_err(fail)?;
        let mut sequence = serializer.serialize_seq(Some(keys.len()))?;
        for position in 0..keys.len() {
            let location = locations.next_location().map_err(fail)?.ok_or_else(|| {
                fail(IndexError::InvalidIndex(
                    "location count differs from the ID count".into(),
                ))
            })?;
            let digest = hex(&location.sha256);
            sequence.serialize_element(&JsonEntry {
                id: keys.id(position),
                offset: location.offset,
                byte_length: location.byte_length,
                sequence_length: location.sequence_length,
                record: position + 1,
                header_line: location.header_line,
                sha256: &digest,
            })?;
        }
        sequence.end()
    }
}

/// Counts written sidecar bytes and refuses to exceed the byte budget.
struct BudgetWriter<W> {
    writer: W,
    written: u64,
    limit: u64,
    failure: Option<IndexError>,
}

impl<W: Write> Write for BudgetWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let attempted = u64::try_from(bytes.len())
            .ok()
            .and_then(|length| self.written.checked_add(length));
        match attempted {
            Some(attempted) if attempted <= self.limit => {
                let count = self.writer.write(bytes)?;
                self.written = self
                    .written
                    .saturating_add(u64::try_from(count).unwrap_or(u64::MAX));
                Ok(count)
            }
            attempted => {
                self.failure = Some(attempted.map_or(IndexError::CountOverflow, |observed| {
                    IndexError::LimitExceeded {
                        resource: "index_bytes",
                        limit: self.limit,
                        observed,
                    }
                }));
                Err(io::Error::other("sidecar byte budget exceeded"))
            }
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

fn hex(digest: &[u8; 32]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut text = String::with_capacity(64);
    for byte in digest {
        text.push(char::from(DIGITS[usize::from(byte >> 4)]));
        text.push(char::from(DIGITS[usize::from(byte & 0x0f)]));
    }
    text
}

fn decode_digest(text: &str) -> Result<[u8; 32], IndexError> {
    validate_digest(text)?;
    let mut digest = [0; 32];
    for (byte, pair) in digest.iter_mut().zip(text.as_bytes().chunks_exact(2)) {
        *byte = (nibble(pair[0]) << 4) | nibble(pair[1]);
    }
    Ok(digest)
}

/// Decodes one lowercase hexadecimal digit already checked by [`validate_digest`].
const fn nibble(digit: u8) -> u8 {
    if digit.is_ascii_digit() {
        digit - b'0'
    } else {
        digit - b'a' + 10
    }
}

fn read_entry<A: IndexAlphabet>(
    file: &mut File,
    id: &str,
    location: &Location,
    format: IndexedFormat,
    limits: ReadLimits,
) -> Result<Record<A>, IndexError> {
    file.seek(SeekFrom::Start(location.offset))?;
    let mut hashing = HashingReader {
        reader: file.take(location.byte_length),
        digest: Sha256::new(),
    };
    let parsed = one_record(RecordReader::<_, A>::new(
        BufReader::new(&mut hashing),
        format,
        limits,
    ));
    // Even a parser failure must not stop digest validation at a valid prefix.
    io::copy(&mut hashing, &mut io::sink())?;
    if hashing.reader.limit() != 0 || hashing.digest.finalize().as_slice() != location.sha256 {
        return Err(IndexError::StaleSource(
            "record span bytes differ from the verified SHA-256".into(),
        ));
    }
    let record = parsed.map_err(|error| {
        IndexError::InvalidIndex(format!("verified record span could not be parsed: {error}"))
    })?;
    if record.id() != id || as_u64(record.sequence().len())? != location.sequence_length {
        return Err(IndexError::InvalidIndex(
            "retrieved record ID or sequence length differs from the index".into(),
        ));
    }
    Ok(record)
}

#[derive(Default)]
struct SizeWriter {
    bytes: u64,
    failure: Option<IndexError>,
}

impl Write for SizeWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let count = as_u64(bytes.len()).and_then(|length| {
            self.bytes
                .checked_add(length)
                .ok_or(IndexError::CountOverflow)
        });
        match count {
            Ok(count) => {
                self.bytes = count;
                Ok(bytes.len())
            }
            Err(error) => {
                self.failure = Some(error);
                Err(io::Error::other("index byte count overflow"))
            }
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn serialized_size(value: &impl Serialize) -> Result<u64, IndexError> {
    let mut writer = SizeWriter::default();
    if let Err(error) = serde_json::to_writer(&mut writer, value) {
        return Err(writer.failure.unwrap_or_else(|| json_error(error)));
    }
    Ok(writer.bytes)
}

fn check_limit(resource: &'static str, limit: usize, observed: u64) -> Result<(), IndexError> {
    let limit = as_u64(limit)?;
    if observed > limit {
        return Err(IndexError::LimitExceeded {
            resource,
            limit,
            observed,
        });
    }
    Ok(())
}

fn as_u64(value: usize) -> Result<u64, IndexError> {
    u64::try_from(value).map_err(|_| IndexError::CountOverflow)
}

fn allocation_error(error: std::collections::TryReserveError) -> IndexError {
    IndexError::Io(io::Error::other(error))
}

fn json_error(error: serde_json::Error) -> IndexError {
    if let Some(kind) = error.io_error_kind() {
        IndexError::Io(io::Error::new(kind, error))
    } else {
        IndexError::InvalidIndex(error.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn selected_span_digest_is_checked_independently_of_source_metadata() {
        struct Directory(PathBuf);
        impl Drop for Directory {
            fn drop(&mut self) {
                std::fs::remove_dir_all(&self.0).unwrap();
            }
        }

        const SOURCE: &[u8] = b">a\nACGT\n";
        let stamp = SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let directory = std::env::temp_dir().join(format!(
            "bioferris-index-digest-test-{}-{stamp}",
            std::process::id()
        ));
        std::fs::create_dir(&directory).unwrap();
        let directory = Directory(directory);
        let path = directory.0.join("source.fa");
        let mut source = File::options()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        source.write_all(SOURCE).unwrap();
        source.flush().unwrap();
        let state = SourceState::from_metadata(&source.metadata().unwrap()).unwrap();
        let mut location = Location {
            offset: 0,
            byte_length: u64::try_from(SOURCE.len()).unwrap(),
            sequence_length: 4,
            header_line: 1,
            sha256: Sha256::digest(SOURCE).into(),
        };
        let record = read_entry::<Dna>(
            &mut source,
            "a",
            &location,
            IndexedFormat::Fasta,
            ReadLimits::default(),
        )
        .unwrap();
        assert_eq!(record.id(), "a");
        assert_eq!(record.source_sequence(), "ACGT");
        verify_source(&source, &path, &state).unwrap();

        // Corrupt only the expected span digest. Source bytes and metadata stay
        // stable, so the independent metadata layer cannot mask this assertion.
        location.sha256 = [0; 32];
        let result = read_entry::<Dna>(
            &mut source,
            "a",
            &location,
            IndexedFormat::Fasta,
            ReadLimits::default(),
        );
        assert!(matches!(result, Err(IndexError::StaleSource(ref message))
            if message == "record span bytes differ from the verified SHA-256"));
        verify_source(&source, &path, &state).unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), SOURCE);
    }

    #[test]
    fn serialized_size_uses_exact_json_escaping_and_detects_overflow() {
        let escaped = "a\\b\"\n\t\u{0000}µ";
        assert_eq!(
            serialized_size(&escaped).unwrap(),
            u64::try_from(serde_json::to_vec(escaped).unwrap().len()).unwrap()
        );
        let mut writer = SizeWriter {
            bytes: u64::MAX,
            failure: None,
        };
        assert!(writer.write_all(b"a").is_err());
        assert!(matches!(writer.failure, Some(IndexError::CountOverflow)));
        assert_eq!(writer.bytes, u64::MAX);
    }

    #[test]
    fn exact_and_zero_limits_are_inclusive() {
        assert!(check_limit("test", 0, 0).is_ok());
        assert!(check_limit("test", 5, 5).is_ok());
        assert!(matches!(
            check_limit("test", 0, 1),
            Err(IndexError::LimitExceeded {
                resource: "test",
                limit: 0,
                observed: 1
            })
        ));
    }

    #[test]
    fn digest_syntax_and_error_source_contracts() {
        use std::error::Error;

        assert!(validate_digest(EMPTY_DIGEST).is_ok());
        assert!(validate_digest(&"F".repeat(64)).is_err());
        assert!(validate_digest(&"g".repeat(64)).is_err());
        assert!(validate_digest(&"0".repeat(63)).is_err());
        let input = IndexError::Input(BioIoError::new(
            IoErrorKind::InvalidSequence,
            "fasta",
            "invalid residue",
        ));
        assert_eq!(input.code(), "invalid_sequence");
        assert!(input.source().is_some());
        let output = IndexError::DestinationExists(PathBuf::from("existing.json"));
        assert_eq!(output.code(), "output_error");
        assert!(output.source().is_none());
        assert_eq!(IndexError::CountOverflow.code(), "count_overflow");
    }

    #[test]
    fn sidecar_requires_objects_and_rejects_duplicate_fields() {
        let entry = serde_json::json!({
            "id": "a", "offset": 0, "byte_length": 5, "sequence_length": 1,
            "record": 1, "header_line": 1, "sha256": EMPTY_DIGEST
        });
        let mut manifest = serde_json::json!({
            "format_id": FORMAT_ID, "schema_version": 1, "molecule": "dna",
            "source_format": "fasta", "key_policy": KEY_POLICY, "source_len": 5,
            "source_sha256": EMPTY_DIGEST, "entries": [entry]
        });
        let valid = serde_json::to_string(&manifest).unwrap();
        assert!(serde_json::from_str::<ObjectOnly<Manifest>>(&valid).is_ok());
        let positional = serde_json::json!([
            FORMAT_ID,
            1,
            "dna",
            "fasta",
            KEY_POLICY,
            5,
            EMPTY_DIGEST,
            manifest["entries"]
        ]);
        assert!(serde_json::from_str::<ObjectOnly<Manifest>>(&positional.to_string()).is_err());
        manifest["entries"][0] = serde_json::json!(["a", 0, 5, 1, 1, 1, EMPTY_DIGEST]);
        assert!(serde_json::from_str::<ObjectOnly<Manifest>>(&manifest.to_string()).is_err());
        let duplicate_root = valid.replacen(
            "\"format_id\":",
            "\"format_id\":\"duplicate\",\"format_id\":",
            1,
        );
        assert!(serde_json::from_str::<ObjectOnly<Manifest>>(&duplicate_root).is_err());
        let duplicate_entry = valid.replacen("\"id\":", "\"id\":\"duplicate\",\"id\":", 1);
        assert!(serde_json::from_str::<ObjectOnly<Manifest>>(&duplicate_entry).is_err());
    }

    #[test]
    fn eof_probe_retries_interruption_without_mistaking_bytes_for_eof() {
        struct InterruptedOnce {
            interrupt: bool,
            bytes: io::Cursor<Vec<u8>>,
        }

        impl Read for InterruptedOnce {
            fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
                if self.interrupt {
                    self.interrupt = false;
                    return Err(io::Error::from(io::ErrorKind::Interrupted));
                }
                self.bytes.read(buffer)
            }
        }

        for bytes in [Vec::new(), vec![1]] {
            let empty = bytes.is_empty();
            let mut reader = InterruptedOnce {
                interrupt: true,
                bytes: io::Cursor::new(bytes),
            };
            assert_eq!(at_eof(&mut reader).unwrap(), empty);
            assert!(at_eof(&mut reader).unwrap());
        }
    }
}
