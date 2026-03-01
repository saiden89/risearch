//! Borrowed view over a normalized base sequence.
//!
//! This mirrors the `String`/`str` split at the API level:
//! - `Sequence` owns `Vec<Base>`
//! - `SeqView` borrows `&[Base]`
//!
//! Keep hot-path internals on raw slices for now; `SeqView` is primarily for
//! boundary type safety and API clarity.

use std::ops::Deref;

use crate::types::Base;

use super::Sequence;

/// Borrowed sequence view.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeqView<'a>(&'a [Base]);

impl<'a> SeqView<'a> {
    #[inline]
    pub const fn new(bases: &'a [Base]) -> Self {
        Self(bases)
    }

    #[inline]
    pub const fn as_slice(&self) -> &'a [Base] {
        self.0
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = &Base> {
        self.0.iter()
    }
}

impl<'a> Deref for SeqView<'a> {
    type Target = [Base];

    #[inline]
    fn deref(&self) -> &[Base] {
        self.0
    }
}

impl<'a> AsRef<[Base]> for SeqView<'a> {
    #[inline]
    fn as_ref(&self) -> &[Base] {
        self.0
    }
}

impl<'a> From<&'a [Base]> for SeqView<'a> {
    #[inline]
    fn from(value: &'a [Base]) -> Self {
        Self(value)
    }
}

impl<'a> From<&'a Sequence> for SeqView<'a> {
    #[inline]
    fn from(value: &'a Sequence) -> Self {
        Self(&value[..])
    }
}
