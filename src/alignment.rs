use smallvec::SmallVec;

use crate::types::Base;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairClass {
    Match,
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
            (Base::A, Base::U) | (Base::U, Base::A) | (Base::G, Base::C) | (Base::C, Base::G) => {
                Self::Match
            }
            (Base::G, Base::U) | (Base::U, Base::G) => Self::Wobble,
            _ => Self::Mismatch,
        }
    }

    pub const fn symbol(self) -> char {
        match self {
            Self::Match => 'P',
            Self::Wobble => 'W',
            Self::Mismatch => 'U',
            Self::TargetBulge => 'T',
            Self::QueryBulge => 'Q',
        }
    }

    pub const fn alignment_symbol(self) -> char {
        match self {
            Self::Match => '|',
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
    pub fn new(left: &[PairClass], seed: &[PairClass], right: &[PairClass]) -> Self {
        let mut steps = SmallVec::with_capacity(left.len() + seed.len() + right.len());
        steps.extend_from_slice(left);
        steps.extend_from_slice(seed);
        steps.extend_from_slice(right);
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
                'P' => PairClass::Match,
                'W' => PairClass::Wobble,
                'U' => PairClass::Mismatch,
                'T' => PairClass::TargetBulge,
                'Q' => PairClass::QueryBulge,
                _ => PairClass::Mismatch,
            })
            .collect::<Vec<_>>();
        Self::new(&[], &steps, &[])
    }
}
