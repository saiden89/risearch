use crate::alignment::Pairing;
use crate::seq::Sequence;
/// Build alignment for the seed region.
///
/// Creates a vector of Pairing entries representing the base-pair interactions
/// in the seed region. Used by both extension and seed-only paths.
///
/// # Arguments
/// * `query` - Query sequence wrapper
/// * `target` - Target sequence wrapper
/// * `q_pos` - Starting position in query (0-based)
/// * `t_match_end` - Ending position in target (0-based, antiparallel)
/// * `len` - Length of seed
pub fn build_seed_alignment(
    query: &Sequence,
    target: &Sequence,
    q_pos: usize,
    t_match_end: usize,
    len: usize,
) -> smallvec::SmallVec<[Pairing; 64]> {
    let mut seed_alignment = smallvec::SmallVec::new();
    for n in 0..len {
        let q_idx = q_pos + n;
        let t_idx = t_match_end.saturating_sub(n);
        let q_b = query[q_idx];
        let t_b = target[t_idx];
        seed_alignment.push(Pairing::from_bases(q_b, t_b));
    }
    seed_alignment
}
