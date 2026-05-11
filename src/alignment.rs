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

    /// Classify pairing between a query base and a transformed target-view base.
    #[inline]
    pub const fn from_view_bases(query: Base, target_view: Base) -> Self {
        Self::from_bases(query, target_view.complement())
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    steps: SmallVec<[PairClass; 128]>,
}

impl Alignment {
    pub fn from_steps(steps: SmallVec<[PairClass; 128]>) -> Self {
        Self { steps }
    }

    pub fn from_parts(prefix: &[PairClass], core: &[PairClass], suffix: &[PairClass]) -> Self {
        let mut steps = SmallVec::with_capacity(prefix.len() + core.len() + suffix.len());
        steps.extend_from_slice(prefix);
        steps.extend_from_slice(core);
        steps.extend_from_slice(suffix);
        Self { steps }
    }

    #[inline]
    pub fn steps(&self) -> &[PairClass] {
        &self.steps
    }

    pub fn fingerprint(&self) -> String {
        self.steps.iter().map(|&p| p.symbol()).collect()
    }

    pub fn from_c_output(interaction: &str, _target_seq: &str) -> Self {
        let steps = interaction
            .chars()
            .map(|c| match c {
                'P' => PairClass::Canonical,
                'W' => PairClass::Wobble,
                'U' => PairClass::Mismatch,
                'T' => PairClass::TargetBulge,
                'Q' => PairClass::QueryBulge,
                _ => PairClass::Mismatch,
            })
            .collect();
        Self::from_steps(steps)
    }
}
