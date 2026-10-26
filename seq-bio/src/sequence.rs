use align_core::{DnaSequence, RnaSequence, SequenceError};

use crate::translation::{
    Translation, TranslationConfig, TranslationError, TranslationOptions, translate_configured,
    translate_standard,
};

/// Translation with an explicit supported NCBI code and forward reading frame.
///
/// This separate trait keeps existing [`NucleotideSequenceExt`] implementations
/// and standard-code convenience calls unchanged. DNA and RNA use equivalent
/// concrete assignments; the whole input is validated before frame selection.
///
/// ```
/// use bioferris::prelude::*;
/// let translation = dna("CATATGAAGA")?.translate_with(TranslationConfig {
///     genetic_code: GeneticCode::VertebrateMitochondrial,
///     reading_frame: ReadingFrame::One,
///     ..TranslationConfig::default()
/// })?;
/// assert_eq!(translation.as_str(), "MW*");
/// assert_eq!(translation.skipped_leading_bases(), 1);
/// assert_eq!(translation.ignored_trailing_bases(), 0);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait ConfiguredTranslationExt {
    /// Translates the suffix selected by the code, frame and explicit policies.
    ///
    /// Offsets equal to the input length select an empty suffix; larger offsets
    /// are rejected. Partial-codon and CDS-length errors describe the suffix,
    /// while internal-stop offsets always refer to the original input.
    ///
    /// # Errors
    /// Returns [`TranslationError`] for an out-of-bounds frame, rejected partial
    /// codons or a CDS failing the selected table's start/stop/interior checks.
    fn translate_with(&self, config: TranslationConfig) -> Result<Translation, TranslationError>;
}

/// Immutable biological operations for the existing typed DNA and RNA sequences.
pub trait NucleotideSequenceExt: Sized {
    /// Translates with standard code 1, retaining stops and rejecting partial codons.
    ///
    /// This is shorthand for [`Self::translate_standard`] with default options.
    ///
    /// # Errors
    ///
    /// Returns [`TranslationError::PartialCodon`] for an incomplete trailing codon.
    fn translate(&self) -> Result<Translation, TranslationError> {
        self.translate_standard(TranslationOptions::default())
    }

    /// Translates a complete standard-code CDS, recoding its start and omitting its stop.
    ///
    /// # Errors
    ///
    /// Returns [`TranslationError`] for incomplete codons, an invalid start,
    /// missing terminal stop, too-short CDS, or internal definite stops.
    fn translate_cds(&self) -> Result<Translation, TranslationError> {
        self.translate_standard(TranslationOptions {
            mode: crate::TranslationMode::Cds,
            ..TranslationOptions::default()
        })
    }

    /// Translates up to the first definite stop, rejecting incomplete trailing codons.
    ///
    /// # Errors
    ///
    /// Returns [`TranslationError::PartialCodon`] for an incomplete trailing codon,
    /// even when a stop appears earlier in the input.
    fn translate_to_stop(&self) -> Result<Translation, TranslationError> {
        self.translate_standard(TranslationOptions {
            mode: crate::TranslationMode::ToStop,
            ..TranslationOptions::default()
        })
    }

    /// Returns the IUPAC complement in the same molecule's alphabet.
    ///
    /// # Errors
    ///
    /// Propagates validation errors from constructing the result. Every valid
    /// input DNA/RNA symbol has a valid complement, including ambiguity codes.
    fn complement(&self) -> Result<Self, SequenceError>;

    /// Returns the reversed IUPAC complement without modifying the input.
    ///
    /// # Errors
    ///
    /// Propagates result validation as described by [`Self::complement`].
    fn reverse_complement(&self) -> Result<Self, SequenceError>;

    /// Translates reading frame zero using only NCBI standard genetic code 1.
    ///
    /// See [`TranslationOptions`] for stop handling and the explicit policy for
    /// trailing incomplete codons. DNA and corresponding RNA give equal results.
    ///
    /// # Errors
    ///
    /// Returns [`TranslationError`] for incomplete codons under the reject policy
    /// or invalid complete CDS boundaries/interior stops. Gaps and invalid symbols
    /// are rejected earlier by the typed sequence parser.
    fn translate_standard(
        &self,
        options: TranslationOptions,
    ) -> Result<Translation, TranslationError>;
}

/// Coding-strand transcription on typed DNA, preserving IUPAC ambiguity.
pub trait DnaSequenceExt {
    /// Returns RNA by replacing T with U, not by complementing the sequence.
    ///
    /// # Errors
    ///
    /// Propagates validation of the RNA result; every typed DNA symbol maps to RNA.
    fn transcribe(&self) -> Result<RnaSequence, SequenceError>;
}

/// Back-transcription on typed RNA, preserving IUPAC ambiguity.
pub trait RnaSequenceExt {
    /// Returns DNA by replacing U with T.
    ///
    /// # Errors
    ///
    /// Propagates validation of the DNA result; every typed RNA symbol maps to DNA.
    fn back_transcribe(&self) -> Result<DnaSequence, SequenceError>;
}

pub(crate) fn complement_byte(byte: u8, adenine_complement: u8) -> u8 {
    match byte {
        b'A' => adenine_complement,
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        // S, W and N are self-complementary. The caller's validated sequence
        // type cannot contain any other byte.
        _ => byte,
    }
}

macro_rules! nucleotide_operations {
    ($sequence:ty, $adenine_complement:expr) => {
        impl ConfiguredTranslationExt for $sequence {
            fn translate_with(
                &self,
                config: TranslationConfig,
            ) -> Result<Translation, TranslationError> {
                let bytes: Vec<_> = self.iter().map(|symbol| symbol.as_ascii()).collect();
                translate_configured(&bytes, config)
            }
        }

        impl NucleotideSequenceExt for $sequence {
            fn complement(&self) -> Result<Self, SequenceError> {
                let bytes: Vec<_> = self
                    .iter()
                    .map(|symbol| complement_byte(symbol.as_ascii(), $adenine_complement))
                    .collect();
                Self::from_ascii(bytes)
            }

            fn reverse_complement(&self) -> Result<Self, SequenceError> {
                let bytes: Vec<_> = self
                    .as_slice()
                    .iter()
                    .rev()
                    .map(|symbol| complement_byte(symbol.as_ascii(), $adenine_complement))
                    .collect();
                Self::from_ascii(bytes)
            }

            fn translate_standard(
                &self,
                options: TranslationOptions,
            ) -> Result<Translation, TranslationError> {
                let bytes: Vec<_> = self.iter().map(|symbol| symbol.as_ascii()).collect();
                translate_standard(&bytes, options)
            }
        }
    };
}

nucleotide_operations!(DnaSequence, b'T');
nucleotide_operations!(RnaSequence, b'U');

impl DnaSequenceExt for DnaSequence {
    fn transcribe(&self) -> Result<RnaSequence, SequenceError> {
        let bytes: Vec<_> = self
            .iter()
            .map(|symbol| match symbol.as_ascii() {
                b'T' => b'U',
                byte => byte,
            })
            .collect();
        RnaSequence::from_ascii(bytes)
    }
}

impl RnaSequenceExt for RnaSequence {
    fn back_transcribe(&self) -> Result<DnaSequence, SequenceError> {
        let bytes: Vec<_> = self
            .iter()
            .map(|symbol| match symbol.as_ascii() {
                b'U' => b'T',
                byte => byte,
            })
            .collect();
        DnaSequence::from_ascii(bytes)
    }
}
