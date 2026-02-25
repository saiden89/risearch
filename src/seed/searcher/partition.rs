use crate::types::{Base, PackedSaEntry};

use super::LINEAR_PARTITION_CUTOFF;

/// C-style SA class rank at `offset` for the suffix ranked at `sa_idx`.
///
/// Rank order: Gap(0) < A(1) < C(2) < G(3) < N(4) < U(5)
/// This mirrors `XSTRM` + `"acgnu"` comparisons in RIsearch2.
#[inline(always)]
fn sa_char_rank(sa: &[u64], seq: &[Base], sa_idx: usize, offset: usize) -> u8 {
    const BASE_TO_C_ORDER: [u8; 8] = [
        0, // Gap
        1, // A
        3, // G
        2, // C
        5, // U
        4, // N
        0, 0,
    ];
    // SAFETY: `sa_idx` is bounded by caller; `suffix_pos + offset` is valid
    // because the backing sequence slices are padded with Gap sentinels.
    let suffix_pos = (unsafe { *sa.get_unchecked(sa_idx) } & PackedSaEntry::POS_MASK) as usize;
    let base = unsafe { *seq.get_unchecked(suffix_pos + offset) } as u8;
    BASE_TO_C_ORDER[base as usize]
}

/// Partition an SA interval by base at given offset.
///
/// Returns C-style boundaries:
/// `[A_start, C_start, G_start, N_start, U_start, end]`.
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

    // Tiny intervals dominate deeper recursion; linear partition wins there.
    if end - start <= LINEAR_PARTITION_CUTOFF {
        let mut i = start;

        while i < end && sa_char_rank(sa, seq, i, offset) < 1 {
            i += 1;
        }
        out[0] = i;
        while i < end && sa_char_rank(sa, seq, i, offset) < 2 {
            i += 1;
        }
        out[1] = i;
        while i < end && sa_char_rank(sa, seq, i, offset) < 3 {
            i += 1;
        }
        out[2] = i;
        while i < end && sa_char_rank(sa, seq, i, offset) < 4 {
            i += 1;
        }
        out[3] = i;
        while i < end && sa_char_rank(sa, seq, i, offset) < 5 {
            i += 1;
        }
        out[4] = i;
        out[5] = end;
        return;
    }

    out[0] = sa_search_left(sa, seq, start, end, offset, 1); // >= A
    out[1] = sa_search_left(sa, seq, start, end, offset, 2); // >= C
    out[2] = sa_search_left(sa, seq, start, end, offset, 3); // >= G
    out[3] = sa_search_left(sa, seq, start, end, offset, 4); // >= N
    out[4] = sa_search_left(sa, seq, start, end, offset, 5); // >= U
    out[5] = end;
}

#[cfg(test)]
#[inline(always)]
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
/// Mirrors C's sa_search_left exactly: halving pattern, no bounds checks.
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
        // SAFETY: mid is within [start, end) which is within the real SA range.
        let suffix_pos = (unsafe { *sa.get_unchecked(mid) } & PackedSaEntry::POS_MASK) as usize;
        // SAFETY: suffix_pos + offset is within the padded sequence slice.
        let base = unsafe { *seq.get_unchecked(suffix_pos + offset) } as u8;
        let char_val = match base {
            1 => 1, // A
            3 => 2, // C
            2 => 3, // G
            5 => 4, // N
            4 => 5, // U
            _ => 0, // Gap/sentinel
        };

        if char_val >= target {
            end = start + half;
        } else {
            start += if half != 0 { half } else { 1 };
        }
        half >>= 1;
    }
    start
}
