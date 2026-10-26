use core::mem::size_of;

use align_core::{Alphabet, Symbol};
use align_msa::{AlignedRow, MultipleAlignment, SequenceRecord};

use crate::PoaError;

pub(crate) fn verify_source_recovery<A: Alphabet>(
    alignment: &MultipleAlignment<A>,
    records: &[SequenceRecord<A>],
) -> Result<(), PoaError> {
    if alignment.sequence_count() != records.len() {
        return verification("rectangular serialization changed the sequence count");
    }
    for (row, record) in alignment.rows().iter().zip(records) {
        if row.identifier() != record.identifier() {
            return verification("rectangular serialization changed input order or identifiers");
        }
        if !row
            .symbols()
            .iter()
            .filter_map(|cell| *cell)
            .eq(record.sequence().iter())
        {
            return verification(format!(
                "rectangular serialization does not recover source {:?}",
                record.identifier()
            ));
        }
    }
    Ok(())
}

pub(crate) fn modeled_alignment_bytes<A: Alphabet>(
    alignment: &MultipleAlignment<A>,
) -> Result<usize, PoaError> {
    let cells = alignment
        .sequence_count()
        .checked_mul(alignment.width())
        .ok_or(PoaError::ArithmeticOverflow {
            context: "rectangular POA cell count",
        })?;
    let mut bytes =
        cells
            .checked_mul(size_of::<Option<Symbol<A>>>())
            .ok_or(PoaError::ArithmeticOverflow {
                context: "rectangular POA cell bytes",
            })?;
    bytes = bytes
        .checked_add(
            alignment
                .sequence_count()
                .checked_mul(size_of::<AlignedRow<A>>())
                .ok_or(PoaError::ArithmeticOverflow {
                    context: "rectangular POA row bytes",
                })?,
        )
        .ok_or(PoaError::ArithmeticOverflow {
            context: "rectangular POA storage",
        })?;
    for row in alignment.rows() {
        bytes = bytes
            .checked_add(row.identifier().len())
            .ok_or(PoaError::ArithmeticOverflow {
                context: "rectangular POA identifiers",
            })?;
    }
    Ok(bytes)
}

fn verification<T>(message: impl Into<String>) -> Result<T, PoaError> {
    Err(PoaError::VerificationFailed {
        message: message.into(),
    })
}
