pub mod adapter;
pub mod alignment;
pub mod cli {
    #[path = "args/mod.rs"]
    pub mod args;
}
pub mod config;
pub mod dp;
pub mod dsm;
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
pub use index::store::TargetRegistry;
pub use registry::QueryRegistry;
pub use seq::{AlignedSeq, Sequence};
pub use types::{Base, Energy, PairType, Strand, BASE_COUNT, GAP};

// Search API re-exports for library usage
pub use config::{
    ExtendConfig, FilterConfig, OutputCompression, OutputConfig, OutputFormat, ScoreConfig,
    SearchConfig, SeedConfig,
};
pub use search::{run_search, run_search_in_memory, SearchHit};
pub use seed::SeedHit;
pub use types::DsmId;
