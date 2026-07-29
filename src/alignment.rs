use std::cmp::Ordering;
use std::ops::Range;

use smallvec::SmallVec;

use crate::types::{Base, PairType};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairClass {
    Canonical,
    Wobble,
    Mismatch,
    TargetBulge,
    QueryBulge,
}

impl PairClass {
    /// Classify two physical bases from the same duplex column.
    #[inline]
    pub const fn from_bases(query: Base, target: Base) -> Self {
        match (query, target) {
            (Base::Gap, _) => Self::TargetBulge,
            (_, Base::Gap) => Self::QueryBulge,
            _ => match query.pair_type(target) {
                PairType::Canonical => Self::Canonical,
                PairType::Wobble => Self::Wobble,
                PairType::Mismatch => Self::Mismatch,
            },
        }
    }

    pub const fn symbol(self) -> char {
        match self {
            Self::Canonical => 'P',
            Self::Wobble => 'W',
            Self::Mismatch => 'U',
            Self::TargetBulge => 'T',
            Self::QueryBulge => 'Q',
        }
    }

    /// Canonical rank for ordering pairings. Derived from [`symbol`](Self::symbol)
    /// so a printed fingerprint and a sort of the same columns can never disagree;
    /// note this is deliberately not the variant declaration order.
    #[inline]
    pub const fn rank(self) -> u8 {
        self.symbol() as u8
    }

    pub const fn alignment_symbol(self) -> char {
        match self {
            Self::Canonical => '|',
            Self::Wobble => ':',
            _ => ' ',
        }
    }

    #[inline]
    pub const fn consumes_query(self) -> bool {
        !matches!(self, Self::TargetBulge)
    }

    #[inline]
    pub const fn consumes_target(self) -> bool {
        !matches!(self, Self::QueryBulge)
    }
}

/// One rendered column: what the pair is, and the bases on either side.
///
/// `query` and `target` are real bases, `Base::Gap` on the side a bulge skips,
/// so a formatter can emit a column without consulting a sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignColumn {
    pub class: PairClass,
    pub query: Base,
    pub target: Base,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    columns: SmallVec<[AlignColumn; 128]>,
    seed: Option<Range<usize>>,
}

impl Alignment {
    /// Resolve flank classes and a gap-free core of `core_len` into columns.
    ///
    /// Every column's class must agree with its own bases, which is what pins the
    /// flanks' order: a flank handed over reversed lands its classes on the wrong
    /// bases. Checked under `debug_assertions`.
    pub fn from_parts(
        prefix: &[PairClass],
        core_len: usize,
        suffix: &[PairClass],
        query: &[Base],
        target: &[Base],
    ) -> Self {
        let classes = prefix
            .iter()
            .copied()
            .map(Some)
            .chain(std::iter::repeat_n(None, core_len))
            .chain(suffix.iter().copied().map(Some));
        let out = Self::resolve(
            classes,
            Some(prefix.len()..prefix.len() + core_len),
            query,
            target,
        );
        debug_assert!(
            out.columns
                .iter()
                .all(|c| c.class == PairClass::from_bases(c.query, c.target)),
            "alignment column disagrees with its own bases: {:?}",
            out.columns
        );
        out
    }

    /// Resolve an alignment produced elsewhere, keeping its classes as given.
    ///
    /// Unlike [`Self::from_parts`] the classes are not re-derived: a foreign
    /// aligner's chemistry is data to compare against, not an invariant to hold.
    ///
    /// Exists so the parity harness can populate [`SearchHit::alignment`] from
    /// another implementation's output. Not part of the library's own pipeline.
    ///
    /// [`SearchHit::alignment`]: crate::SearchHit::alignment
    #[doc(hidden)]
    pub fn from_classes(
        classes: &[PairClass],
        seed: Option<Range<usize>>,
        query: &[Base],
        target: &[Base],
    ) -> Self {
        Self::resolve(classes.iter().copied().map(Some), seed, query, target)
    }

    /// Walk physical duplex columns in their canonical order: 5'->3' along the
    /// query and 3'->5' along the target. A `None` class is derived, which only
    /// works where the column is known to be gap-free.
    fn resolve(
        classes: impl Iterator<Item = Option<PairClass>>,
        seed: Option<Range<usize>>,
        query: &[Base],
        target: &[Base],
    ) -> Self {
        let mut columns = SmallVec::with_capacity(classes.size_hint().0);
        let (mut q, mut t) = (0usize, 0usize);
        for class in classes {
            let query_base = if class.is_none_or(PairClass::consumes_query) {
                let base = query.get(q).copied().unwrap_or(Base::Gap);
                q += 1;
                base
            } else {
                Base::Gap
            };
            let target_base = if class.is_none_or(PairClass::consumes_target) {
                let base = target.get(t).copied().unwrap_or(Base::Gap);
                t += 1;
                base
            } else {
                Base::Gap
            };
            columns.push(AlignColumn {
                class: class.unwrap_or(PairClass::from_bases(query_base, target_base)),
                query: query_base,
                target: target_base,
            });
        }
        Self { columns, seed }
    }

    #[inline]
    pub fn columns(&self) -> &[AlignColumn] {
        &self.columns
    }

    /// Span of the seed core within [`columns`](Self::columns).
    #[inline]
    pub fn seed(&self) -> Option<Range<usize>> {
        self.seed.clone()
    }

    pub fn fingerprint(&self) -> String {
        self.columns.iter().map(|c| c.class.symbol()).collect()
    }

    /// Column pairings as ranks, for lexicographic comparison without allocating.
    pub fn pairing_ranks(&self) -> impl Iterator<Item = u8> + '_ {
        self.columns.iter().map(|c| c.class.rank())
    }

    /// Compare pairing fingerprints lexicographically without allocating.
    pub(crate) fn fingerprint_cmp(&self, other: &Self) -> Ordering {
        self.pairing_ranks().cmp(other.pairing_ranks())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn columns_keep_physical_order_and_advance_only_the_consumed_strand() {
        let query = [Base::A, Base::C];
        let target = [Base::U, Base::G];
        let prefix = [PairClass::TargetBulge, PairClass::QueryBulge];

        let alignment = Alignment::from_parts(&prefix, 1, &[], &query, &target);

        assert_eq!(alignment.seed(), Some(2..3));
        let columns = alignment
            .columns()
            .iter()
            .map(|c| (c.class, c.query, c.target))
            .collect::<Vec<_>>();
        assert_eq!(
            columns,
            [
                (PairClass::TargetBulge, Base::Gap, Base::U),
                (PairClass::QueryBulge, Base::A, Base::Gap),
                (PairClass::Canonical, Base::C, Base::G),
            ]
        );
    }
}
