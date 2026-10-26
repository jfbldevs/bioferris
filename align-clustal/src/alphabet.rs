use align_core::{Alphabet, Dna, Protein, Rna};

mod sealed {
    pub trait Sealed {}
}

impl sealed::Sealed for Dna {}
impl sealed::Sealed for Rna {}
impl sealed::Sealed for Protein {}

/// An alphabet supported by the versioned Clustal Omega profiles.
///
/// This sealed marker is implemented for [`Protein`], [`Dna`], and [`Rna`].
pub trait ClustalAlphabet: Alphabet + sealed::Sealed {}

impl ClustalAlphabet for Dna {}
impl ClustalAlphabet for Rna {}
impl ClustalAlphabet for Protein {}
