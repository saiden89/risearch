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
pub use seq::{Seq, reverse_complement_dna, reverse_complement_rna};
pub use types::{BASE_COUNT, Base, Energy, QueryId, SeedPairing, Span, Strand, TargetId};

// Search API re-exports for library usage
pub use search::{SearchHit, run_search_collect};
pub use seed::{SeedCandidate, SeedSpec};
pub use types::{Alignment, Pairing};
