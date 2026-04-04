use super::*;
use crate::config::{MismatchSpec, SeedConfig, SeedSpec};
use crate::index::sa::SuffixArray;
use crate::seq::Sequence;
use crate::types::Base;

fn build_padded_sa(seq: &Sequence) -> (Vec<u64>, Vec<Base>, usize) {
    let sa = SuffixArray::try_from(&seq[..]).expect("SA construction failed");
    let real_len = sa.len();
    let mut padded: Vec<u64> = sa.into_inner();
    padded.resize(padded.len() + crate::index::store::SA_CHAR_PADDING, 0u64);
    let mut padded_seq: Vec<Base> = seq.iter().copied().collect();
    padded_seq.resize(
        padded_seq.len() + crate::index::store::SA_CHAR_PADDING,
        Base::Gap,
    );
    (padded, padded_seq, real_len)
}

#[test]
fn partition_splits_by_base() {
    let seq_bases = vec![Base::A, Base::G, Base::C, Base::U];
    let seq = Sequence::from(seq_bases);
    let (padded_sa, padded_seq, real_len) = build_padded_sa(&seq);
    let parts = partition((padded_sa.as_slice(), padded_seq.as_slice(), real_len), 0, real_len, 0);

    assert_eq!(parts[1] - parts[0], 1); // A
    assert_eq!(parts[2] - parts[1], 1); // C
    assert_eq!(parts[3] - parts[2], 1); // G
    assert_eq!(parts[4] - parts[3], 0); // N
    assert_eq!(parts[5] - parts[4], 1); // U
}

#[test]
fn singleton_handoff_does_not_double_emit() {
    // Distinct first symbols force singleton intervals at depth=1 in both
    // query and target branches, exercising the singleton fast-path handoff.
    let q_seq = Sequence::from(vec![Base::A, Base::C]);
    let t_seq = Sequence::from(vec![Base::U, Base::G]);
    let (q_sa, q_seq_padded, q_sa_len) = build_padded_sa(&q_seq);
    let (t_sa, t_seq_padded, t_sa_len) = build_padded_sa(&t_seq);

    let cfg = SeedConfig::with_wobble(SeedSpec::LengthOnly(1), MismatchSpec::exact(), false);
    let searcher = SeedSearcher::new(
        (q_sa.as_slice(), q_seq_padded.as_slice(), q_sa_len),
        0,
        (t_sa.as_slice(), t_seq_padded.as_slice(), t_sa_len),
        &cfg,
    );

    let mut results = Vec::new();
    searcher.search_length_range(1, 2, &mut results);

    let mut seen = std::collections::HashSet::new();
    for m in &results {
        assert!(
            seen.insert((
                m.query_interval.start,
                m.query_interval.end,
                m.target_interval.start,
                m.target_interval.end,
                m.seed_len
            )),
            "duplicate seed match emitted: {:?}",
            m
        );
    }
}
