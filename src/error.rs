//! The error type returned by this crate's library API.
//!
//! Variants exist to be matched on, and to map onto distinct exception classes
//! in the Python bindings. The binary layer keeps `anyhow`, which absorbs these
//! through `?`.

/// Anything a library call can fail with.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// Filesystem, mapping, or write failure. `fs-err` carries the path.
    #[error(transparent)]
    Io(#[from] std::io::Error),

    /// A suffix-array index that could not be built, written, or read.
    #[error("{0}")]
    Index(String),

    /// A bundled energy model that could not be selected or parsed.
    #[error("{0}")]
    Dsm(String),

    /// Query or target input that cannot be searched.
    #[error("{0}")]
    Input(String),

    /// Parameters that do not describe a runnable search.
    #[error("{0}")]
    Config(String),

    /// An output destination or codec that cannot be used as requested.
    #[error("{0}")]
    Output(String),
}

/// [`Result`](std::result::Result) with this crate's [`Error`].
pub type Result<T> = std::result::Result<T, Error>;
