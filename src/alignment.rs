use smallvec::SmallVec;
use std::fmt;
use std::ops::Range;

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
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Pairing {
    query: Base,
    target: Base,
}

impl Pairing {
    #[inline]
    pub const fn new(query: Base, target: Base) -> Self {
        Self { query, target }
    }

    #[inline]
    pub const fn from_bases(query: Base, target: Base) -> Self {
        Self::new(query, target)
    }

    #[inline]
    pub const fn query_bulge(query_base: Base) -> Self {
        Self::new(query_base, Base::Gap)
    }

    #[inline]
    pub const fn target_bulge(target_base: Base) -> Self {
        Self::new(Base::Gap, target_base)
    }

    pub fn class(self) -> PairClass {
        match (self.query, self.target) {
            (Base::Gap, _) => PairClass::TargetBulge,
            (_, Base::Gap) => PairClass::QueryBulge,
            (q, t) => match q.pairing_class(t) {
                'P' => PairClass::Match,
                'W' => PairClass::Wobble,
                _ => PairClass::Mismatch,
            },
        }
    }

    pub fn query_char(self) -> char {
        self.query.as_char().to_ascii_lowercase()
    }

    pub fn target_char(self) -> char {
        self.target.as_char().to_ascii_lowercase()
    }

    fn from_fingerprint_char(c: char, target: Base) -> Self {
        match c {
            'P' => Pairing::from_bases(Self::canonical_mate(target), target),
            'W' => Pairing::from_bases(Self::wobble_mate(target), target),
            'U' => Pairing::from_bases(Self::mismatch_base(target), target),
            'T' => Pairing::target_bulge(target),
            'Q' => Pairing::query_bulge(Base::N),
            _ => Pairing::from_bases(Base::N, target),
        }
    }

    const fn canonical_mate(base: Base) -> Base {
        match base {
            Base::A => Base::U,
            Base::U => Base::A,
            Base::G => Base::C,
            Base::C => Base::G,
            _ => Base::N,
        }
    }

    const fn wobble_mate(base: Base) -> Base {
        match base {
            Base::G => Base::U,
            Base::U => Base::G,
            _ => Base::N,
        }
    }

    const fn mismatch_base(base: Base) -> Base {
        match base {
            Base::A | Base::U => Base::C,
            Base::G | Base::C => Base::A,
            _ => Base::N,
        }
    }
}

impl fmt::Display for Pairing {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.class().symbol())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    steps: SmallVec<[Pairing; 128]>,
    seed_range: Range<usize>,
}

impl Alignment {
    pub fn new(left: &[Pairing], seed: &[Pairing], right: &[Pairing]) -> Self {
        let left_len = left.len();
        let seed_len = seed.len();

        let mut steps = SmallVec::with_capacity(left_len + seed_len + right.len());
        steps.extend_from_slice(left);
        steps.extend_from_slice(seed);
        steps.extend_from_slice(right);

        Self {
            steps,
            seed_range: left_len..(left_len + seed_len),
        }
    }

    #[inline]
    pub fn steps(&self) -> &[Pairing] {
        &self.steps
    }

    #[inline]
    fn write_mapped(&self, buf: &mut Vec<u8>, map: fn(Pairing) -> u8) {
        for &p in self.steps.iter() {
            buf.push(map(p));
        }
    }

    pub fn write_query_seq(&self, buf: &mut Vec<u8>) {
        self.write_mapped(buf, |p| p.query_char() as u8);
    }

    pub fn write_target_seq(&self, buf: &mut Vec<u8>) {
        self.write_mapped(buf, |p| p.target_char() as u8);
    }

    pub fn write_alignment_line(&self, buf: &mut Vec<u8>) {
        self.write_mapped(buf, |p| p.class().alignment_symbol() as u8);
    }

    pub fn write_pairing_string(&self, buf: &mut Vec<u8>) {
        self.write_mapped(buf, |p| p.class().symbol() as u8);
    }

    pub fn seed(&self) -> &[Pairing] {
        &self.steps[self.seed_range.start..self.seed_range.end]
    }

    pub fn left_extension(&self) -> &[Pairing] {
        &self.steps[..self.seed_range.start]
    }

    pub fn fingerprint(&self) -> String {
        self.steps.iter().map(|&p| p.class().symbol()).collect()
    }

    pub fn target_sequence(&self) -> String {
        self.steps.iter().map(|&p| p.target_char()).collect()
    }

    pub fn query_sequence(&self) -> String {
        self.steps.iter().map(|&p| p.query_char()).collect()
    }

    pub fn from_c_output(
        interaction: &str,
        target_seq: &str,
        seed_start: Option<usize>,
        seed_end: Option<usize>,
    ) -> Self {
        let mut target_iter = target_seq.chars();
        let mut steps = Vec::with_capacity(interaction.len());
        for fp in interaction.chars() {
            let target = Base::from_byte(target_iter.next().unwrap_or('-') as u8);
            steps.push(Pairing::from_fingerprint_char(fp, target));
        }

        match (seed_start, seed_end) {
            (Some(s), Some(e)) => {
                let s = s.min(steps.len());
                let e = e.min(steps.len()).max(s);
                Self::new(&steps[..s], &steps[s..e], &steps[e..])
            }
            _ => Self::new(&[], &steps, &[]),
        }
    }
}
