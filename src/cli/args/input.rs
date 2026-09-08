use std::path::{Path, PathBuf};

#[derive(clap::Args, Debug, Clone)]
pub(crate) struct InputArgs {
    /// FASTA file for query sequence(s) (.fa or .fa.gz) -- use '-' for stdin
    #[arg(short = 'q', long = "query", value_name = "FILE")]
    pub(crate) query: PathBuf,

    /// Target index file (created by `index` command)
    #[arg(
        short = 't',
        long = "target",
        value_name = "TARGET",
        required_unless_present = "legacy_target",
        conflicts_with = "legacy_target"
    )]
    pub(crate) target: Option<PathBuf>,

    /// DEPRECATED: legacy alias for -t/--target
    #[arg(
        short = 'i',
        value_name = "TARGET",
        hide = true,
        required_unless_present = "target",
        conflicts_with = "target"
    )]
    pub(crate) legacy_target: Option<PathBuf>,
}

impl InputArgs {
    pub(crate) fn target_path(&self) -> &Path {
        match (&self.target, &self.legacy_target) {
            (Some(target), None) | (None, Some(target)) => target.as_path(),
            _ => unreachable!("clap should enforce exactly one target argument"),
        }
    }

    pub(crate) fn uses_legacy_target(&self) -> bool {
        self.legacy_target.is_some()
    }
}
