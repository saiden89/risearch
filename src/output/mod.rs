//! Output formatting and destinations.

pub mod format;
pub use format::format_hit_into;

pub mod sink;
pub use sink::TextSink;

pub mod writer;
pub use writer::OutputWriter;
