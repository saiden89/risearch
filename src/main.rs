//! Command-line entry point for `risearch`.

#[cfg(not(risearch_sanitize))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> anyhow::Result<()> {
    risearch::cli_main()
}
