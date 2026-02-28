use crate::types::Base;

use super::{sa_char_rank, RANK_A, RANK_C, RANK_G, RANK_N, RANK_U};

const LINEAR_PARTITION_CUTOFF: usize = 1024;

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

        while i < end && sa_char_rank(sa, seq, i, offset) < RANK_A {
            i += 1;
        }
        out[0] = i;
        while i < end && sa_char_rank(sa, seq, i, offset) < RANK_C {
            i += 1;
        }
        out[1] = i;
        while i < end && sa_char_rank(sa, seq, i, offset) < RANK_G {
            i += 1;
        }
        out[2] = i;
        while i < end && sa_char_rank(sa, seq, i, offset) < RANK_N {
            i += 1;
        }
        out[3] = i;
        while i < end && sa_char_rank(sa, seq, i, offset) < RANK_U {
            i += 1;
        }
        out[4] = i;
        out[5] = end;
        return;
    }

    out[0] = sa_search_left(sa, seq, start, end, offset, RANK_A);
    out[1] = sa_search_left(sa, seq, start, end, offset, RANK_C);
    out[2] = sa_search_left(sa, seq, start, end, offset, RANK_G);
    out[3] = sa_search_left(sa, seq, start, end, offset, RANK_N);
    out[4] = sa_search_left(sa, seq, start, end, offset, RANK_U);
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
        let char_val = sa_char_rank(sa, seq, mid, offset);

        if char_val >= target {
            end = start + half;
        } else {
            start += if half != 0 { half } else { 1 };
        }
        half >>= 1;
    }
    start
}
