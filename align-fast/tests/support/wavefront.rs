#![allow(dead_code)]

use align_core::AlignmentOperation;

const INFINITY: i64 = i64::MAX / 4;

/// Independent dense Gotoh minimization using the WFA `o + k*e` convention.
pub fn dense_affine(
    first: &[u8],
    second: &[u8],
    mismatch: i64,
    opening: i64,
    extension: i64,
) -> i64 {
    let columns = second.len() + 1;
    let cells = (first.len() + 1) * columns;
    let mut main = vec![INFINITY; cells];
    let mut insertion = vec![INFINITY; cells];
    let mut deletion = vec![INFINITY; cells];
    main[0] = 0;
    let opening_transition = opening + extension;

    for first_prefix in 0..=first.len() {
        for second_prefix in 0..=second.len() {
            if first_prefix == 0 && second_prefix == 0 {
                continue;
            }
            let index = first_prefix * columns + second_prefix;
            if second_prefix > 0 {
                let left = index - 1;
                insertion[index] =
                    add(main[left], opening_transition).min(add(insertion[left], extension));
            }
            if first_prefix > 0 {
                let upper = index - columns;
                deletion[index] =
                    add(main[upper], opening_transition).min(add(deletion[upper], extension));
            }
            let diagonal = if first_prefix > 0 && second_prefix > 0 {
                add(
                    main[index - columns - 1],
                    if first[first_prefix - 1] == second[second_prefix - 1] {
                        0
                    } else {
                        mismatch
                    },
                )
            } else {
                INFINITY
            };
            main[index] = diagonal.min(deletion[index]).min(insertion[index]);
        }
    }
    main[first.len() * columns + second.len()]
}

fn add(value: i64, increment: i64) -> i64 {
    if value == INFINITY {
        INFINITY
    } else {
        value + increment
    }
}

/// Exhaustively enumerates all affine edit paths for tiny inputs.
pub fn exhaustive_affine(
    first: &[u8],
    second: &[u8],
    mismatch: i64,
    opening: i64,
    extension: i64,
) -> i64 {
    #[allow(clippy::too_many_arguments)]
    fn visit(
        first: &[u8],
        second: &[u8],
        first_index: usize,
        second_index: usize,
        previous_gap: Option<AlignmentOperation>,
        mismatch: i64,
        opening: i64,
        extension: i64,
    ) -> i64 {
        if first_index == first.len() && second_index == second.len() {
            return 0;
        }
        let mut best = INFINITY;
        if first_index < first.len() && second_index < second.len() {
            let diagonal = if first[first_index] == second[second_index] {
                0
            } else {
                mismatch
            } + visit(
                first,
                second,
                first_index + 1,
                second_index + 1,
                None,
                mismatch,
                opening,
                extension,
            );
            best = best.min(diagonal);
        }
        if first_index < first.len() {
            let operation = AlignmentOperation::Deletion;
            let gap = if previous_gap == Some(operation) {
                extension
            } else {
                opening + extension
            };
            best = best.min(
                gap + visit(
                    first,
                    second,
                    first_index + 1,
                    second_index,
                    Some(operation),
                    mismatch,
                    opening,
                    extension,
                ),
            );
        }
        if second_index < second.len() {
            let operation = AlignmentOperation::Insertion;
            let gap = if previous_gap == Some(operation) {
                extension
            } else {
                opening + extension
            };
            best = best.min(
                gap + visit(
                    first,
                    second,
                    first_index,
                    second_index + 1,
                    Some(operation),
                    mismatch,
                    opening,
                    extension,
                ),
            );
        }
        best
    }

    visit(first, second, 0, 0, None, mismatch, opening, extension)
}

/// Replays a complete path and returns its WFA affine penalty.
pub fn replay_affine(
    first: &[u8],
    second: &[u8],
    operations: &[AlignmentOperation],
    mismatch: i64,
    opening: i64,
    extension: i64,
) -> i64 {
    let (penalty, first_end, second_end) =
        replay(first, second, operations, mismatch, opening, extension);
    assert_eq!(first_end, first.len(), "path must consume the first input");
    assert_eq!(
        second_end,
        second.len(),
        "path must consume the second input"
    );
    penalty
}

fn replay(
    first: &[u8],
    second: &[u8],
    operations: &[AlignmentOperation],
    mismatch: i64,
    opening: i64,
    extension: i64,
) -> (i64, usize, usize) {
    let mut first_index = 0_usize;
    let mut second_index = 0_usize;
    let mut penalty = 0_i64;
    let mut previous_gap = None;
    for &operation in operations {
        match operation {
            AlignmentOperation::Match => {
                assert_eq!(first[first_index], second[second_index]);
                first_index += 1;
                second_index += 1;
                previous_gap = None;
            }
            AlignmentOperation::Substitution => {
                assert_ne!(first[first_index], second[second_index]);
                first_index += 1;
                second_index += 1;
                penalty += mismatch;
                previous_gap = None;
            }
            AlignmentOperation::Deletion => {
                first_index += 1;
                penalty += if previous_gap == Some(operation) {
                    extension
                } else {
                    opening + extension
                };
                previous_gap = Some(operation);
            }
            AlignmentOperation::Insertion => {
                second_index += 1;
                penalty += if previous_gap == Some(operation) {
                    extension
                } else {
                    opening + extension
                };
                previous_gap = Some(operation);
            }
        }
    }
    (penalty, first_index, second_index)
}

/// Enumerates all binary strings with lengths through `maximum_len`.
pub fn binary_strings(maximum_len: usize) -> Vec<Vec<u8>> {
    let mut strings = Vec::new();
    for length in 0..=maximum_len {
        for bits in 0..(1_usize << length) {
            strings.push(
                (0..length)
                    .map(|index| if bits & (1 << index) == 0 { b'A' } else { b'C' })
                    .collect(),
            );
        }
    }
    strings
}
