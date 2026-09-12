//! Command-line entry point for `risearch`.

// Sanitizer builds use the instrumented system allocator.
#[cfg(not(risearch_sanitize))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() -> anyhow::Result<()> {
    risearch::cli_main()
}
