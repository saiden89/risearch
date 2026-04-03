use crate::types::Base;

use super::sa_char;

const LINEAR_PARTITION_CUTOFF: usize = 1024;

/// Base discriminants that mark partition boundaries (sorted by SA order).
const BOUNDARIES: [u8; 5] = [
    Base::A as u8,
    Base::C as u8,
    Base::G as u8,
    Base::N as u8,
    Base::U as u8,
];

/// Partition an SA interval by base at given offset.
///
/// Returns boundaries: `[A_start, C_start, G_start, N_start, U_start, end]`.
#[inline(always)]
pub(super) fn partition_interval_into(
    sa: &[u64],
    seq: &[Base],
    start: usize,
    end: usize,
    offset: usize,
    out: &mut [usize; 6],
) {
    if start >= end {
        *out = [start; 6];
        return;
    }

    if end - start <= LINEAR_PARTITION_CUTOFF {
        let mut i = start;
        for (slot, &target) in BOUNDARIES.iter().enumerate() {
            while i < end && sa_char(sa, seq, i, offset) < target {
                i += 1;
            }
            out[slot] = i;
        }
    } else {
        for (slot, &target) in BOUNDARIES.iter().enumerate() {
            out[slot] = sa_search_left(sa, seq, start, end, offset, target);
        }
    }
    out[5] = end;
}

#[cfg(test)]
pub(super) fn partition_interval(
    sa: &[u64],
    seq: &[Base],
    start: usize,
    end: usize,
    offset: usize,
) -> [usize; 6] {
    let mut out = [0usize; 6];
    partition_interval_into(sa, seq, start, end, offset, &mut out);
    out
}

/// Binary search for leftmost position where character at `offset` >= `target`.
///
/// Mirrors C's sa_search_left exactly: pre-halving pattern where `half` is
/// halved each iteration instead of recomputed from `(end - start) / 2`.
/// SA must have SA_CHAR_PADDING sentinel entries beyond the real data.
#[inline(always)]
fn sa_search_left(
    sa: &[u64],
    seq: &[Base],
    mut start: usize,
    mut end: usize,
    offset: usize,
    target: u8,
) -> usize {
    let mut half = (end - start) >> 1;
    while start < end {
        let mid = start + half;
        let char_val = sa_char(sa, seq, mid, offset);

        if char_val >= target {
            end = start + half;
        } else {
            start += if half != 0 { half } else { 1 };
        }
        half >>= 1;
    }
    start
}
