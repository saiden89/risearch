//! Output formatting and destinations.

#[doc(hidden)]
pub mod format;

pub mod sink;
pub use sink::TextSink;

pub(crate) mod writer;
