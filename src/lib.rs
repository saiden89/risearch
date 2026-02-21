pub mod alignment;
pub mod cli {
    #[path = "args/mod.rs"]
    pub mod args;
}
pub mod config;
pub mod dp;
pub mod dsm;
mod dsm_extend;
pub mod fastx;
pub mod index;
pub mod output;
pub mod registry;
pub mod search;
pub mod seed;
pub mod seq;
pub mod types; // Core domain types

// Re-exports for convenience
pub use alignment::{Alignment, PairClass};
pub use index::store::TargetStore;
pub use registry::{QueryRegistry, TargetRegistry};
pub use seq::{AlignedSeq, Sequence};
pub use types::{Base, Energy, SeedPairingMode, Strand, BASE_COUNT};

// Search API re-exports for library usage
pub use config::{MismatchSpec, SeedSpec};
pub use search::SearchHit;
pub use seed::SeedHit;
