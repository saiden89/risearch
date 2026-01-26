pub mod cli_args;
pub mod config;
pub mod dp;
pub mod dsm;
pub mod fastx;
pub mod index;
pub mod output;
pub mod sa;
pub mod search;
pub mod seed;
pub mod seq;
pub mod types; // Core domain types - must be first

// Re-exports for convenience
pub use dsm::{Dsm, StackPair};
pub use seq::{AlignedSeq, Sequence, reverse_complement_dna, reverse_complement_rna};
pub use types::{BASE_COUNT, Base, Energy, QueryId, SeedPairingMode, Span, Strand, TargetId};

// Search API re-exports for library usage
pub use output::{write_results, write_results_with_format};
pub use search::SearchHit;
pub use seed::{SeedCandidate, SeedSpec};
pub use types::{Alignment, Pairing};
