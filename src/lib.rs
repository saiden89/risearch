pub mod args;
pub mod dp;
pub mod dsm;
pub mod io;
pub mod parallel_sa;
pub mod sa;
pub mod search;
pub mod seed;
pub mod seq;
pub mod types; // Core domain types - must be first

// Re-exports for convenience
pub use dsm::{Dsm, StackPair};
pub use seq::{AlignedSeq, Seq, reverse_complement_dna, reverse_complement_rna};
pub use types::{BASE_COUNT, Base, Energy, QueryId, SeedPairing, Span, Strand, TargetId};

// Search API re-exports for library usage
pub use search::{SearchHit, run_search, write_results};
pub use seed::{SeedCandidate, SeedSpec};
pub use types::{Alignment, Pairing};
