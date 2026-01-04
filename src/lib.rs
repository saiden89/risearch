pub mod args;
pub mod dsm;
#[cfg(feature = "fm-index")]
pub mod fm;
pub mod io;
pub mod sa;
pub mod search;
pub mod seed;
pub mod types; // Core domain types - must be first

// Re-exports for convenience
pub use types::{BASE_COUNT, Base, SeedPairing, Strand};
