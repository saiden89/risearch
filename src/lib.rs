pub mod args;
pub mod dp;
pub mod dsm;
#[cfg(feature = "fm-index")]
pub mod fm;
pub mod io;
pub mod sa;
pub mod search;
pub mod seed;
pub mod seq;
pub mod types; // Core domain types - must be first

// Re-exports for convenience
pub use dsm::{Dsm, StackPair};
pub use seq::Seq;
pub use types::{BASE_COUNT, Base, Energy, QueryId, SeedPairing, Strand, TargetId};

// Search API re-exports for library usage
pub use search::{Alignment, SearchHit, run_search_collect};
