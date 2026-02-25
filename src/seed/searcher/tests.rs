use super::*;
use crate::index::sa::SuffixArray;
use crate::seq::Sequence;
use crate::types::Base;

fn build_padded_sa(seq: &Sequence) -> (Vec<u64>, Vec<Base>, usize) {
    let sa = SuffixArray::try_from(seq).expect("SA construction failed");
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
fn test_partition_basic() {
    let seq_bases = vec![Base::A, Base::G, Base::C, Base::U];
    let seq = Sequence::from(seq_bases);
    let (padded_sa, padded_seq, real_len) = build_padded_sa(&seq);
    let parts = partition_interval(&padded_sa, &padded_seq, 0, real_len, 0);

    assert_eq!(parts[1] - parts[0], 1); // A
    assert_eq!(parts[2] - parts[1], 1); // C
    assert_eq!(parts[3] - parts[2], 1); // G
    assert_eq!(parts[4] - parts[3], 0); // N
    assert_eq!(parts[5] - parts[4], 1); // U
}
