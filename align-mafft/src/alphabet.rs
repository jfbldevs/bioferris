use align_core::{Alphabet, Dna, Protein, Rna};

mod sealed {
    pub trait Sealed {}
}

impl sealed::Sealed for Dna {}
impl sealed::Sealed for Rna {}
impl sealed::Sealed for Protein {}

/// An alphabet supported by the MAFFT reconstruction.
///
/// This sealed marker is implemented for [`Protein`], [`Dna`], and [`Rna`].
pub trait MafftAlphabet: Alphabet + sealed::Sealed {}

impl MafftAlphabet for Dna {}
impl MafftAlphabet for Rna {}
impl MafftAlphabet for Protein {}
