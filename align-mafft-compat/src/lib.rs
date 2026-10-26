//! Native-output-compatible MAFFT 7.526 profiles backed by a pinned pure-Rust port.
//!
//! This crate is deliberately separate from [`align_mafft`]. The latter keeps
//! the independently reconstructed, publication-versioned profiles. This crate
//! adapts `mafft` 0.1.2, a source-guided Rust port of native MAFFT 7.526, into
//! the workspace's typed sequence and multiple-alignment contracts.
//!
//! No C code, subprocess, Python runtime, or external executable is used at
//! runtime. This is nevertheless a delegated compatibility implementation, not
//! an independent reconstruction from the scientific publications.
//!
//! Exact-output evidence is corpus-scoped, not a guarantee for arbitrary inputs.
//! Expanded diagnostics expose L-INS-i differences on divergent DNA and RNA;
//! see `docs/validation/msa-native-fidelity-v1.md` in the repository.
#![forbid(unsafe_code)]

use core::any::TypeId;
use core::marker::PhantomData;
use std::collections::BTreeSet;

use align_core::{Dna, Protein, Rna, Symbol};
use align_mafft::{MafftAlphabet, MafftError};
use align_msa::{AlignedRow, MsaError, MultipleAlignment, SequenceRecord};
use mafft_7526_engine::{AlignmentMode, MafftEngine, SeqType, Sequence, SequenceSet};

/// Stable identifier for native MAFFT 7.526 L-INS-i compatibility.
pub const MAFFT_7526_LINSI_PROFILE_IDENTIFIER: &str = "MAFFT-7.526-L-INS-i-localpair-v1";

/// Stable identifier for native MAFFT 7.526 protein E-INS-i compatibility.
pub const MAFFT_7526_EINSI_PROFILE_IDENTIFIER: &str = "MAFFT-7.526-E-INS-i-genafpair-protein-v1";

/// Version of the pinned source-guided Rust engine used by this adapter.
pub const MAFFT_7526_PORT_VERSION: &str = "0.1.2";

/// Native MAFFT release reproduced by the pinned compatibility engine.
pub const MAFFT_7526_NATIVE_VERSION: &str = "7.526 (2024/Apr/26)";

/// The selected native-compatible MAFFT profile.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Mafft7526Profile {
    /// L-INS-i using all-pairs local alignment information.
    Linsi,
    /// E-INS-i using generalized-affine pair alignment information.
    Einsi,
}

impl Mafft7526Profile {
    /// Returns the immutable compatibility-profile identifier.
    #[must_use]
    pub const fn identifier(self) -> &'static str {
        match self {
            Self::Linsi => MAFFT_7526_LINSI_PROFILE_IDENTIFIER,
            Self::Einsi => MAFFT_7526_EINSI_PROFILE_IDENTIFIER,
        }
    }
}

/// A typed MAFFT 7.526-compatible alignment and adapter evidence.
#[derive(Clone, Debug)]
pub struct Mafft7526Alignment<A: MafftAlphabet> {
    alignment: MultipleAlignment<A>,
    profile: Mafft7526Profile,
    engine_score: f64,
    progressive_step_count: usize,
    thread_count: usize,
}

impl<A: MafftAlphabet> Mafft7526Alignment<A> {
    /// Returns the canonical rectangular alignment in input order.
    #[must_use]
    pub const fn alignment(&self) -> &MultipleAlignment<A> {
        &self.alignment
    }

    /// Returns the selected native-compatible profile.
    #[must_use]
    pub const fn profile(&self) -> Mafft7526Profile {
        self.profile
    }

    /// Returns the immutable compatibility-profile identifier.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        self.profile.identifier()
    }

    /// Returns the score exposed by the pinned compatibility engine.
    #[must_use]
    pub const fn engine_score(&self) -> f64 {
        self.engine_score
    }

    /// Returns the number of recorded progressive merge steps.
    #[must_use]
    pub const fn progressive_step_count(&self) -> usize {
        self.progressive_step_count
    }

    /// Returns the deterministic Rayon worker count used for this call.
    #[must_use]
    pub const fn thread_count(&self) -> usize {
        self.thread_count
    }

    /// Consumes the evidence wrapper and returns the canonical alignment.
    #[must_use]
    pub fn into_alignment(self) -> MultipleAlignment<A> {
        self.alignment
    }
}

/// MAFFT 7.526 L-INS-i compatibility for protein, DNA, or RNA records.
#[derive(Clone, Copy, Debug)]
pub struct Mafft7526Linsi<A: MafftAlphabet> {
    thread_count: usize,
    marker: PhantomData<A>,
}

impl<A: MafftAlphabet> Mafft7526Linsi<A> {
    /// Creates the pinned native-compatible L-INS-i profile.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            thread_count: 1,
            marker: PhantomData,
        }
    }

    /// Returns the immutable compatibility-profile identifier.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        MAFFT_7526_LINSI_PROFILE_IDENTIFIER
    }

    /// Returns the deterministic Rayon worker count.
    #[must_use]
    pub const fn thread_count(self) -> usize {
        self.thread_count
    }

    /// Replaces the deterministic worker count.
    ///
    /// # Errors
    ///
    /// Returns [`MafftError::InvalidThreadCount`] when `count` is zero.
    pub const fn with_thread_count(mut self, count: usize) -> Result<Self, MafftError> {
        if count == 0 {
            return Err(MafftError::InvalidThreadCount { value: count });
        }
        self.thread_count = count;
        Ok(self)
    }

    /// Aligns records with the pinned MAFFT 7.526 L-INS-i engine.
    ///
    /// # Errors
    ///
    /// Returns [`MafftError`] for invalid records, worker-pool failure, or an
    /// invalid/source-changing result from the delegated compatibility engine.
    pub fn align(
        &self,
        records: &[SequenceRecord<A>],
    ) -> Result<Mafft7526Alignment<A>, MafftError> {
        run_engine(
            records,
            AlignmentMode::LInsi { iterations: 1_000 },
            Mafft7526Profile::Linsi,
            self.thread_count,
        )
    }
}

impl<A: MafftAlphabet> Default for Mafft7526Linsi<A> {
    fn default() -> Self {
        Self::new()
    }
}

/// MAFFT 7.526 E-INS-i compatibility for protein records.
///
/// The pinned upstream port is byte-compatible with the pinned native protein
/// reference corpus. Its generalized-affine DNA/RNA route still diverges on
/// adversarial inputs, so those alphabets are intentionally unavailable here.
#[derive(Clone, Copy, Debug)]
pub struct Mafft7526Einsi {
    thread_count: usize,
}

impl Mafft7526Einsi {
    /// Creates the pinned native-compatible protein E-INS-i profile.
    #[must_use]
    pub const fn new() -> Self {
        Self { thread_count: 1 }
    }

    /// Returns the immutable compatibility-profile identifier.
    #[must_use]
    pub const fn profile_identifier(&self) -> &'static str {
        MAFFT_7526_EINSI_PROFILE_IDENTIFIER
    }

    /// Returns the deterministic Rayon worker count.
    #[must_use]
    pub const fn thread_count(self) -> usize {
        self.thread_count
    }

    /// Replaces the deterministic worker count.
    ///
    /// # Errors
    ///
    /// Returns [`MafftError::InvalidThreadCount`] when `count` is zero.
    pub const fn with_thread_count(mut self, count: usize) -> Result<Self, MafftError> {
        if count == 0 {
            return Err(MafftError::InvalidThreadCount { value: count });
        }
        self.thread_count = count;
        Ok(self)
    }

    /// Aligns protein records with the pinned MAFFT 7.526 E-INS-i engine.
    ///
    /// # Errors
    ///
    /// Returns [`MafftError`] for invalid records, worker-pool failure, or an
    /// invalid/source-changing result from the delegated compatibility engine.
    pub fn align(
        &self,
        records: &[SequenceRecord<Protein>],
    ) -> Result<Mafft7526Alignment<Protein>, MafftError> {
        run_engine(
            records,
            AlignmentMode::EInsi { iterations: 1_000 },
            Mafft7526Profile::Einsi,
            self.thread_count,
        )
    }
}

impl Default for Mafft7526Einsi {
    fn default() -> Self {
        Self::new()
    }
}

fn run_engine<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    mode: AlignmentMode,
    profile: Mafft7526Profile,
    thread_count: usize,
) -> Result<Mafft7526Alignment<A>, MafftError> {
    run_engine_with_type(records, mode, profile, thread_count, sequence_type::<A>())
}

fn native_sequence_bytes<A: MafftAlphabet>(
    record: &SequenceRecord<A>,
    native_type: SeqType,
) -> Vec<u8> {
    record
        .sequence()
        .iter()
        .map(|symbol| {
            let ascii = symbol.as_ascii();
            if native_type.is_nucleotide() {
                ascii.to_ascii_lowercase()
            } else {
                ascii
            }
        })
        .collect()
}

fn run_engine_with_type<A: MafftAlphabet>(
    records: &[SequenceRecord<A>],
    mode: AlignmentMode,
    profile: Mafft7526Profile,
    thread_count: usize,
    native_type: SeqType,
) -> Result<Mafft7526Alignment<A>, MafftError> {
    validate_records(records)?;

    let active_records = records
        .iter()
        .filter(|record| !record.sequence().is_empty())
        .collect::<Vec<_>>();
    if active_records.is_empty() {
        let rows = records
            .iter()
            .map(|record| AlignedRow::new(record.identifier(), Vec::new()))
            .collect::<Result<Vec<_>, _>>()?;
        return Ok(Mafft7526Alignment {
            alignment: MultipleAlignment::try_new(rows)?,
            profile,
            engine_score: 0.0,
            progressive_step_count: 0,
            thread_count,
        });
    }

    let input = SequenceSet {
        sequences: active_records
            .iter()
            .map(|record| Sequence {
                name: record.identifier().to_owned(),
                data: native_sequence_bytes(record, native_type),
            })
            .collect(),
        seq_type: native_type,
    };
    let engine = MafftEngine::new(mode).with_reorder(false);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .map_err(|_| MafftError::InternalInvariant {
            message: "MAFFT 7.526 compatibility worker pool could not be constructed",
        })?;
    let delegated = pool.install(|| engine.align(&input));
    if delegated.sequences.len() != active_records.len()
        || delegated.names.len() != active_records.len()
    {
        return Err(MafftError::InternalInvariant {
            message: "MAFFT 7.526 compatibility engine changed the row count",
        });
    }
    if !delegated.score.is_finite() {
        return Err(MafftError::NonFiniteComputation {
            stage: "7.526 compatibility engine score",
        });
    }
    let width = delegated.width();
    let mut active_index = 0usize;
    let mut rows = Vec::with_capacity(records.len());
    for record in records {
        let cells = if record.sequence().is_empty() {
            vec![None; width]
        } else {
            let expected = active_records[active_index];
            let name = &delegated.names[active_index];
            let sequence = &delegated.sequences[active_index];
            active_index += 1;
            if name != expected.identifier() || sequence.len() != width {
                return Err(MafftError::InternalInvariant {
                    message: "MAFFT 7.526 compatibility engine changed row order or width",
                });
            }
            sequence
                .iter()
                .copied()
                .map(|byte| match byte {
                    b'-' | b'.' => Ok(None),
                    residue => Symbol::<A>::from_ascii(residue).map(Some).ok_or(
                        MafftError::InternalInvariant {
                            message: "MAFFT 7.526 compatibility engine emitted an invalid residue",
                        },
                    ),
                })
                .collect::<Result<Vec<_>, _>>()?
        };
        rows.push(AlignedRow::new(record.identifier(), cells)?);
    }
    let alignment = MultipleAlignment::try_new(rows)?;
    for (record, row) in records.iter().zip(alignment.rows()) {
        if row.to_ungapped_ascii() != record.sequence().to_ascii_string() {
            return Err(MafftError::InternalInvariant {
                message: "MAFFT 7.526 compatibility engine changed an input sequence",
            });
        }
    }
    Ok(Mafft7526Alignment {
        alignment,
        profile,
        engine_score: delegated.score,
        progressive_step_count: delegated.step_trace.len(),
        thread_count,
    })
}

fn validate_records<A: MafftAlphabet>(records: &[SequenceRecord<A>]) -> Result<(), MafftError> {
    if records.is_empty() {
        return Err(MsaError::EmptyCollection.into());
    }
    let mut identifiers = BTreeSet::new();
    for record in records {
        if !identifiers.insert(record.identifier()) {
            return Err(MsaError::DuplicateIdentifier {
                identifier: record.identifier().to_owned(),
            }
            .into());
        }
    }
    Ok(())
}

fn sequence_type<A: MafftAlphabet>() -> SeqType {
    if TypeId::of::<A>() == TypeId::of::<Protein>() {
        SeqType::Protein
    } else if TypeId::of::<A>() == TypeId::of::<Dna>() {
        SeqType::Dna
    } else if TypeId::of::<A>() == TypeId::of::<Rna>() {
        SeqType::Rna
    } else {
        unreachable!("MafftAlphabet is sealed to protein, DNA, and RNA")
    }
}

/// Opt-in research frontend for the pinned CLI's implicit molecule selection.
///
/// The existing typed L-INS-i/E-INS-i APIs above remain unchanged. This module
/// does not assert general CLI, error-ABI or nucleotide E-INS-i fidelity.
pub mod native_cli_research {
    use super::{
        AlignmentMode, Mafft7526Alignment, Mafft7526Profile, MafftAlphabet, MafftError, Rna,
        SeqType, SequenceRecord, TypeId, run_engine_with_type, validate_records,
    };
    use core::fmt;

    /// Actual native model category, independent of the input's syntax type.
    #[derive(Clone, Copy, Debug, Eq, PartialEq)]
    pub enum Molecule {
        /// Nucleotide model, including original RNA spelling.
        Nucleotide,
        /// Protein model.
        Protein,
    }

    /// Owned arithmetic evidence of the ordered original-record decision.
    #[derive(Clone, Copy, Debug, Eq, PartialEq)]
    pub struct Decision {
        molecule: Molecule,
        sites: usize,
        bases: usize,
    }

    impl Decision {
        /// Returns the selected native model category.
        #[must_use]
        pub const fn molecule(self) -> Molecule {
            self.molecule
        }
        /// Returns alphabetic sites in the complete sampled record prefix.
        #[must_use]
        pub const fn sites(self) -> usize {
            self.sites
        }
        /// Returns A/T/G/C/U/N sites in that same prefix.
        #[must_use]
        pub const fn bases(self) -> usize {
            self.bases
        }
    }

    /// Safe research rejection, not an emulated native process/error ABI.
    #[derive(Debug)]
    pub enum Error {
        /// The native default protein frontend rejects U without `--anysymbol`.
        NativeProteinU,
        /// Empty source is outside this controlled frontend boundary.
        EmptySource,
        /// Existing checked alignment or original-record validation failed.
        Alignment(MafftError),
    }

    impl fmt::Display for Error {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                Self::NativeProteinU => f.write_str("native automatic protein frontend rejects U"),
                Self::EmptySource => {
                    f.write_str("no alphabetic source for native molecule research")
                }
                Self::Alignment(error) => error.fmt(f),
            }
        }
    }
    impl std::error::Error for Error {}
    impl From<MafftError> for Error {
        fn from(value: MafftError) -> Self {
            Self::Alignment(value)
        }
    }

    /// Determines the implicit model from a whole-record prefix reaching 100k sites.
    ///
    /// The strict three-quarter decision is expressed with exact integer ratios,
    /// not a copied native scanner or a reference-provided decision.
    ///
    /// # Errors
    /// Rejects empty/duplicate identifiers, empty source or count overflow.
    pub fn decision<A: MafftAlphabet>(records: &[SequenceRecord<A>]) -> Result<Decision, Error> {
        validate_records(records)?;
        let (mut sites, mut bases) = (0usize, 0usize);
        for record in records {
            if sites >= 100_000 {
                break;
            }
            sites = sites.checked_add(record.sequence().len()).ok_or(
                MafftError::ArithmeticOverflow {
                    stage: "native molecule site count",
                },
            )?;
            bases = bases
                .checked_add(
                    record
                        .sequence()
                        .iter()
                        .filter(|symbol| {
                            matches!(symbol.as_ascii(), b'A' | b'T' | b'G' | b'C' | b'U' | b'N')
                        })
                        .count(),
                )
                .ok_or(MafftError::ArithmeticOverflow {
                    stage: "native molecule base count",
                })?;
        }
        if sites == 0 {
            return Err(Error::EmptySource);
        }
        Ok(Decision {
            molecule: if (bases as u128) * 4 > (sites as u128) * 3 {
                Molecule::Nucleotide
            } else {
                Molecule::Protein
            },
            sites,
            bases,
        })
    }

    /// Aligns original typed records with the separately observed implicit CLI model.
    ///
    /// The type validates original syntax/ownership; it no longer silently forces
    /// the statistical model in this explicitly separate research entrypoint.
    ///
    /// # Errors
    /// Returns safe frontend/alignment errors, not native stderr or OS signals.
    pub fn align<A: MafftAlphabet>(
        records: &[SequenceRecord<A>],
        profile: Mafft7526Profile,
        threads: usize,
    ) -> Result<(Mafft7526Alignment<A>, Decision), Error> {
        if threads == 0 {
            return Err(MafftError::InvalidThreadCount { value: threads }.into());
        }
        let selected = decision(records)?;
        let native_type = match selected.molecule {
            Molecule::Protein => {
                if records.iter().any(|record| {
                    record
                        .sequence()
                        .iter()
                        .any(|symbol| symbol.as_ascii() == b'U')
                }) {
                    return Err(Error::NativeProteinU);
                }
                SeqType::Protein
            }
            Molecule::Nucleotide if TypeId::of::<A>() == TypeId::of::<Rna>() => SeqType::Rna,
            Molecule::Nucleotide => SeqType::Dna,
        };
        let mode = match profile {
            Mafft7526Profile::Linsi => AlignmentMode::LInsi { iterations: 1_000 },
            Mafft7526Profile::Einsi => AlignmentMode::EInsi { iterations: 1_000 },
        };
        Ok((
            run_engine_with_type(records, mode, profile, threads, native_type)?,
            selected,
        ))
    }
}

#[cfg(test)]
mod tests {
    use align_core::{Dna, DnaSequence, ProteinSequence};
    use align_mafft::MafftError;
    use align_msa::{MsaError, SequenceRecord};

    use super::{
        MAFFT_7526_EINSI_PROFILE_IDENTIFIER, MAFFT_7526_LINSI_PROFILE_IDENTIFIER, Mafft7526Einsi,
        Mafft7526Linsi,
    };

    fn dna(identifier: &str, sequence: &str) -> SequenceRecord<Dna> {
        SequenceRecord::new(identifier, DnaSequence::from_ascii(sequence).unwrap()).unwrap()
    }

    #[test]
    fn original_native_prefix_scalars_match_the_owned_decision_on_large_source() {
        use super::native_cli_research::{Molecule, decision};
        use serde_json::Value;
        let fixture: Value = serde_json::from_str(include_str!(
            "../tests/fixtures/mafft-7.526/native-prefix-scalar-v1.json"
        ))
        .unwrap();
        assert_eq!(fixture["qualified"], false);
        assert_eq!(fixture["large_alignment_fidelity_claimed"], false);
        assert_eq!(fixture["cases"].as_array().unwrap().len(), 5);
        for case in fixture["cases"].as_array().unwrap() {
            let records = case["records"]
                .as_array()
                .unwrap()
                .iter()
                .map(|record| {
                    let source = record["segments"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|segment| {
                            segment["symbol"].as_str().unwrap().repeat(
                                usize::try_from(segment["length"].as_u64().unwrap()).unwrap(),
                            )
                        })
                        .collect::<String>();
                    dna(record["id"].as_str().unwrap(), &source)
                })
                .collect::<Vec<_>>();
            let result = decision(&records).unwrap();
            assert_eq!(
                result.sites() as u64,
                case["source_defined_sampled_sites"].as_u64().unwrap(),
                "{}",
                case["id"]
            );
            assert_eq!(
                result.bases() as u64,
                case["source_defined_sampled_bases"].as_u64().unwrap(),
                "{}",
                case["id"]
            );
            assert_eq!(
                result.molecule(),
                if case["native_type"] == "d" {
                    Molecule::Nucleotide
                } else {
                    Molecule::Protein
                },
                "{}",
                case["id"]
            );
        }
    }

    #[test]
    fn forced_iupac_pair_distances_and_progressive_final_rows_match_original_binary() {
        use super::{
            AlignmentMode, MafftAlphabet, MafftEngine, SeqType, Sequence, SequenceSet,
            native_sequence_bytes,
        };
        use align_core::{Rna, Sequence as TypedSequence};
        use serde_json::{Value, json};

        fn check<A: MafftAlphabet>(case: &Value, kind: SeqType) {
            let records = case["request"]["records"]
                .as_array()
                .unwrap()
                .iter()
                .map(|row| {
                    SequenceRecord::<A>::new(
                        row["id"].as_str().unwrap(),
                        TypedSequence::from_ascii(row["sequence"].as_str().unwrap()).unwrap(),
                    )
                    .unwrap()
                })
                .collect::<Vec<_>>();
            let typed = Mafft7526Linsi::<A>::new().align(&records).unwrap();
            let expected=case["native_final_rows"].as_array().unwrap().iter().map(|row|
                json!({"id":row["id"],"aligned":row["aligned"].as_str().unwrap().to_ascii_uppercase()})
            ).collect::<Vec<_>>();
            assert_eq!(
                typed
                    .alignment()
                    .rows()
                    .iter()
                    .map(|row| json!({"id":row.identifier(),"aligned":row.to_ascii_with_gaps()}))
                    .collect::<Vec<_>>(),
                expected,
                "{}",
                case["id"]
            );
            for iterations in [0, 1000] {
                let input = SequenceSet {
                    seq_type: kind,
                    sequences: records
                        .iter()
                        .map(|record| Sequence {
                            name: record.identifier().to_owned(),
                            data: native_sequence_bytes(record, kind),
                        })
                        .collect(),
                };
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(1)
                    .build()
                    .unwrap();
                let result = pool.install(|| {
                    MafftEngine::new(AlignmentMode::LInsi { iterations })
                        .with_reorder(false)
                        .align(&input)
                });
                let dm = result.distance_matrix.as_ref().unwrap();
                let bits = (0..dm.nseq)
                    .flat_map(|i| (i + 1..dm.nseq).map(move |j| dm.get(i, j).to_bits()))
                    .collect::<Vec<_>>();
                assert_eq!(
                    json!(bits),
                    case["native_distance_upper_bits"],
                    "{} / {iterations}",
                    case["id"]
                );
                let rows=result.names.iter().zip(&result.sequences).map(|(name,sequence)|
                    json!({"id":name,"aligned":std::str::from_utf8(sequence).unwrap()})
                ).collect::<Vec<_>>();
                assert_eq!(
                    json!(rows),
                    case[if iterations == 0 {
                        "native_progressive_rows"
                    } else {
                        "native_final_rows"
                    }],
                    "{} / {iterations}",
                    case["id"]
                );
            }
        }
        let fixture: Value = serde_json::from_str(include_str!(
            "../tests/fixtures/mafft-7.526/forced-iupac-pair-context-v1.json"
        ))
        .unwrap();
        assert_eq!(fixture["cases"].as_array().unwrap().len(), 6);
        assert_eq!(fixture["qualified"], false);
        for case in fixture["cases"].as_array().unwrap() {
            match case["request"]["alphabet"].as_str().unwrap() {
                "dna" => check::<Dna>(case, SeqType::Dna),
                "rna" => check::<Rna>(case, SeqType::Rna),
                _ => panic!("fixture molecule scope changed"),
            }
        }
    }

    #[test]
    fn native_cli_decision_has_a_strict_three_quarter_boundary_and_counts_n() {
        use super::native_cli_research::{Molecule, decision};
        for count in [74, 75, 76] {
            let source = "A".repeat(count) + &"R".repeat(100 - count);
            let result = decision(&[dna("source", &source)]).unwrap();
            assert_eq!((result.sites(), result.bases()), (100, count));
            assert_eq!(
                result.molecule(),
                if count > 75 {
                    Molecule::Nucleotide
                } else {
                    Molecule::Protein
                }
            );
        }
        assert_eq!(
            decision(&[dna("all-n", "NNNN")]).unwrap().molecule(),
            Molecule::Nucleotide
        );
    }

    #[test]
    fn native_cli_prefix_keeps_complete_crossing_records_and_stops_afterwards() {
        use super::native_cli_research::{Molecule, decision};
        let records = [
            dna("a", &"A".repeat(99_999)),
            dna("b", &"R".repeat(100_001)),
            dna("ignored", &"A".repeat(1_000_000)),
        ];
        let result = decision(&records).unwrap();
        assert_eq!(
            (result.sites(), result.bases(), result.molecule()),
            (200_000, 99_999, Molecule::Protein)
        );
        let result = decision(&[dna("a", &"A".repeat(100_000)), dna("ignored", "RRRR")]).unwrap();
        assert_eq!((result.sites(), result.bases()), (100_000, 100_000));
    }

    #[test]
    fn native_cli_model_does_not_mutate_the_existing_typed_api() {
        use super::{Mafft7526Profile, native_cli_research};
        let records = [dna("a", "ACGTRRRR"), dna("b", "ACGTRRRR")];
        let before = records.clone();
        let typed = Mafft7526Linsi::<Dna>::new().align(&records).unwrap();
        let (implicit, selected) =
            native_cli_research::align(&records, Mafft7526Profile::Linsi, 1).unwrap();
        assert_eq!(selected.molecule(), native_cli_research::Molecule::Protein);
        assert_eq!(
            typed.profile_identifier(),
            MAFFT_7526_LINSI_PROFILE_IDENTIFIER
        );
        assert_eq!(typed.alignment(), implicit.alignment());
        assert_eq!(records, before);
    }

    #[test]
    fn native_cli_research_safely_rejects_u_in_a_selected_protein_model() {
        use super::{Mafft7526Profile, native_cli_research};
        use align_core::{Rna, RnaSequence};
        let records =
            [
                SequenceRecord::<Rna>::new("rna", RnaSequence::from_ascii("ACGURRRR").unwrap())
                    .unwrap(),
            ];
        assert!(matches!(
            native_cli_research::align(&records, Mafft7526Profile::Linsi, 1),
            Err(native_cli_research::Error::NativeProteinU)
        ));
        assert!(matches!(
            native_cli_research::align(&records, Mafft7526Profile::Linsi, 0),
            Err(native_cli_research::Error::Alignment(
                MafftError::InvalidThreadCount { value: 0 }
            ))
        ));
        assert!(Mafft7526Linsi::<Rna>::new().align(&records).is_ok());
    }

    #[test]
    fn linsi_is_bit_exact_across_worker_counts() {
        let records = [
            dna("a", "TTACGTACGTAA"),
            dna("b", "GGACGTACTTCC"),
            dna("c", "CCACGTACGTGG"),
        ];
        let serial = Mafft7526Linsi::<Dna>::new().align(&records).unwrap();
        let parallel = Mafft7526Linsi::<Dna>::new()
            .with_thread_count(4)
            .unwrap()
            .align(&records)
            .unwrap();
        assert_eq!(serial.alignment(), parallel.alignment());
        assert_eq!(
            serial.profile_identifier(),
            MAFFT_7526_LINSI_PROFILE_IDENTIFIER
        );
        assert_eq!(serial.thread_count(), 1);
        assert_eq!(parallel.thread_count(), 4);
    }

    #[test]
    fn empty_rows_are_preserved_without_delegating_the_upstream_empty_row_bug() {
        let records = [dna("empty", ""), dna("full", "ACGT")];
        let result = Mafft7526Linsi::<Dna>::new().align(&records).unwrap();
        assert_eq!(result.alignment().rows()[0].to_ascii_with_gaps(), "----");
        assert_eq!(result.alignment().rows()[1].to_ascii_with_gaps(), "ACGT");
    }

    #[test]
    fn all_empty_rows_and_invalid_collections_are_totalized() {
        let empty_rows = [dna("first", ""), dna("second", "")];
        let result = Mafft7526Linsi::<Dna>::new().align(&empty_rows).unwrap();
        assert_eq!(result.alignment().sequence_count(), 2);
        assert_eq!(result.alignment().width(), 0);
        assert_eq!(result.progressive_step_count(), 0);

        assert!(matches!(
            Mafft7526Linsi::<Dna>::new().align(&[]),
            Err(MafftError::Msa(MsaError::EmptyCollection))
        ));
        let duplicates = [dna("same", "ACGT"), dna("same", "AGT")];
        assert!(matches!(
            Mafft7526Linsi::<Dna>::new().align(&duplicates),
            Err(MafftError::Msa(MsaError::DuplicateIdentifier { .. }))
        ));
    }

    #[test]
    fn einsi_is_explicitly_protein_typed() {
        let records = [
            SequenceRecord::new("a", ProteinSequence::from_ascii("ACDEFG").unwrap()).unwrap(),
            SequenceRecord::new("b", ProteinSequence::from_ascii("ACDFG").unwrap()).unwrap(),
        ];
        let result = Mafft7526Einsi::new().align(&records).unwrap();
        let parallel = Mafft7526Einsi::new()
            .with_thread_count(4)
            .unwrap()
            .align(&records)
            .unwrap();
        assert_eq!(
            result.profile_identifier(),
            MAFFT_7526_EINSI_PROFILE_IDENTIFIER
        );
        assert_eq!(result.alignment().sequence_count(), 2);
        assert_eq!(result.alignment(), parallel.alignment());
    }
}
