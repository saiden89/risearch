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
pub use registry::QueryRegistry;
pub use seq::{AlignedSeq, SeqView, Sequence};
pub use types::{Base, Energy, SeedPairingMode, Strand, BASE_COUNT, GAP};

// Search API re-exports for library usage
pub use config::{
    ExtendConfig, FilterConfig, Matrix, MismatchSpec, OutputCompression, OutputConfig,
    OutputFormat, ScoreConfig, SearchConfig, SeedConfig, SeedSpec,
};
pub use search::{run_search, run_search_in_memory, SearchHit};
pub use seed::SeedHit;
