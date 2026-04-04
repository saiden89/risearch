use crate::types::Base;

use super::SeedSaView;

const LINEAR_PARTITION_CUTOFF: usize = 1024;

/// Base discriminants that mark partition boundaries (sorted by SA order).
const BOUNDARIES: [u8; 5] = [
    Base::A as u8,
    Base::C as u8,
    Base::G as u8,
    Base::N as u8,
    Base::U as u8,
];

/// Partition a sorted SA interval by base character at `depth`.
///
/// Returns 6 boundary positions `[A, C, G, N, U, end]` in suffix-array order.
/// Sub-interval for slot `k` is `bounds[k]..bounds[k+1]`.
/// Callers that search only matchable RNA bases should skip the `N` bucket.
#[inline(always)]
pub(super) fn partition<V: SeedSaView>(
    view: V,
    start: usize,
    end: usize,
    depth: usize,
) -> [usize; 6] {
    if start >= end {
        return [start; 6];
    }

    let mut out = [0usize; 6];
    if end - start <= LINEAR_PARTITION_CUTOFF {
        let mut i = start;
        for (slot, &target) in BOUNDARIES.iter().enumerate() {
            while i < end && (view.sa_base(i, depth) as u8) < target {
                i += 1;
            }
            out[slot] = i;
        }
    } else {
        for (slot, &target) in BOUNDARIES.iter().enumerate() {
            out[slot] = binary_search(view, start, end, depth, target);
        }
    }
    out[5] = end;
    out
}

/// Binary search for leftmost position where character at `depth` >= `target`.
#[inline(always)]
fn binary_search<V: SeedSaView>(
    view: V,
    mut start: usize,
    mut end: usize,
    depth: usize,
    target: u8,
) -> usize {
    let mut half = (end - start) >> 1;
    while start < end {
        let mid = start + half;
        if (view.sa_base(mid, depth) as u8) >= target {
            end = start + half;
        } else {
            start += if half != 0 { half } else { 1 };
        }
        half >>= 1;
    }
    start
}
