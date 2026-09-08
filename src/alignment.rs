//! Duplex alignment columns and how each column is classified.
//!
//! A traceback yields a run of [`AlignColumn`]s in physical duplex order, each
//! carrying its [`PairClass`].

use crate::types::{Base, PairType};

/// How one duplex column pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairClass {
    /// Watson-Crick pair: A↔U or C↔G.
    Canonical,
    /// G↔U wobble pair.
    Wobble,
    /// Both sides consumed, neither canonical nor wobble.
    Mismatch,
    /// Target-side insertion: the query side is a gap.
    TargetBulge,
    /// Query-side insertion: the target side is a gap.
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

    /// Single-character class code used by the detailed output format.
    pub const fn symbol(self) -> char {
        match self {
            Self::Canonical => 'P',
            Self::Wobble => 'W',
            Self::Mismatch => 'U',
            Self::TargetBulge => 'T',
            Self::QueryBulge => 'Q',
        }
    }

    /// Connector drawn between the two strands: `|` for canonical, `:` for
    /// wobble, blank otherwise.
    pub const fn alignment_symbol(self) -> char {
        match self {
            Self::Canonical => '|',
            Self::Wobble => ':',
            _ => ' ',
        }
    }
}

/// One resolved alignment column in physical duplex order.
///
/// [`Base::Gap`] marks the side a traceback step did not consume. Normalized
/// query and target sequences never contain gap-ranked input bases, so every
/// other value is an actual consumed input.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignColumn {
    class: PairClass,
    query: Base,
    target: Base,
}

impl AlignColumn {
    /// How this column pairs.
    #[inline]
    pub const fn class(&self) -> PairClass {
        self.class
    }

    /// Query-side base, or [`Base::Gap`] on a target bulge.
    #[inline]
    pub const fn query(&self) -> Base {
        self.query
    }

    /// Target-side base, or [`Base::Gap`] on a query bulge.
    #[inline]
    pub const fn target(&self) -> Base {
        self.target
    }

    /// A diagonal traceback step consuming both inputs.
    #[inline]
    pub const fn paired(query: Base, target: Base) -> Self {
        debug_assert!(
            !matches!(query, Base::Gap) && !matches!(target, Base::Gap),
            "a paired column must consume two non-gap inputs"
        );
        Self {
            class: PairClass::from_bases(query, target),
            query,
            target,
        }
    }

    /// A traceback step consuming only the query input.
    #[inline]
    pub const fn query_only(query: Base) -> Self {
        debug_assert!(
            !matches!(query, Base::Gap),
            "a query-only column must consume a non-gap input"
        );
        Self {
            class: PairClass::QueryBulge,
            query,
            target: Base::Gap,
        }
    }

    /// A traceback step consuming only the target input.
    #[inline]
    pub const fn target_only(target: Base) -> Self {
        debug_assert!(
            !matches!(target, Base::Gap),
            "a target-only column must consume a non-gap input"
        );
        Self {
            class: PairClass::TargetBulge,
            query: Base::Gap,
            target,
        }
    }

    /// Record a paired column classified by another aligner.
    ///
    /// The reported chemistry is comparison data, so it is not re-derived from
    /// the bases. Consumption remains valid because both inputs are required.
    #[doc(hidden)]
    #[inline]
    pub const fn reported_pair(class: PairClass, query: Base, target: Base) -> Self {
        assert!(
            matches!(
                class,
                PairClass::Canonical | PairClass::Wobble | PairClass::Mismatch
            ),
            "a reported pair must consume both inputs"
        );
        assert!(
            !matches!(query, Base::Gap) && !matches!(target, Base::Gap),
            "a reported pair must consume two non-gap inputs"
        );
        Self {
            class,
            query,
            target,
        }
    }
}

/// The canonical fingerprint spelling and traversal order for an alignment.
///
/// Formatting, parity comparison, and deduplication all consume this iterator,
/// so the visible fingerprint and its lexicographic order cannot drift apart.
#[inline]
pub fn fingerprint_symbols(columns: &[AlignColumn]) -> impl Iterator<Item = char> + '_ {
    columns.iter().map(|column| column.class().symbol())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn columns_mark_the_unconsumed_input_with_gap() {
        let target_only = AlignColumn::target_only(Base::U);

        assert_eq!(
            target_only,
            AlignColumn {
                class: PairClass::TargetBulge,
                query: Base::Gap,
                target: Base::U,
            }
        );
    }
}
