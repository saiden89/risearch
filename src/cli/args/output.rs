use crate::config::{self, OutputCompression, OutputFormat};
use anyhow::{bail, Result};
use std::path::PathBuf;

/// Boundary CLI arguments for output destination, formatting, and compression.
#[derive(clap::Args, Debug, Clone)]
pub struct OutputArgs {
    /// Output file for search results (use '-' for stdout)
    #[arg(short = 'o', long = "output", value_name = "FILE", default_value = "-")]
    pub path: PathBuf,

    /// Output format
    #[arg(
        short = 'f',
        long = "format",
        value_name = "FORMAT",
        default_missing_value = "minimal",
        value_enum
    )]
    pub report_format: Option<OutputFormat>,

    /// DEPRECATED: Legacy argument for output format (1=detailed, 2=cigar, 3=binding_site, 4=minimal)
    #[arg(
        short = 'p',
        long = "report-alignment",
        value_name = "MODE",
        num_args = 0..=1,
        default_missing_value = "1",
        help_heading = "Deprecated"
    )]
    pub report_legacy: Option<u8>,

    /// Output compression codec (overrides file extension inference; gzip/gz, zstd/zst accepted)
    #[arg(long = "compress", value_enum)]
    pub output_compress: Option<OutputCompression>,

    /// Output compression level (codec-specific)
    #[arg(long = "compress-level", value_name = "LEVEL")]
    pub output_level: Option<i32>,

    /// Write one output file per query into the directory given by -o
    #[arg(long = "multifile", action = clap::ArgAction::SetTrue)]
    pub output_multifile: bool,
}

impl OutputArgs {
    pub fn validate(&self) -> Result<()> {
        if self.output_multifile && self.path.as_os_str() == "-" {
            bail!("--multifile requires -o/--output to be a directory path; '-' (stdout) is not allowed.");
        }

        Ok(())
    }

    #[inline]
    fn resolved_format(&self) -> OutputFormat {
        if let Some(f) = self.report_format {
            f
        } else if let Some(legacy_mode) = self.report_legacy {
            match legacy_mode {
                1 => config::OutputFormat::Detailed,
                2 => config::OutputFormat::Cigar,
                3 => config::OutputFormat::BindingSite,
                4 => config::OutputFormat::Minimal,
                _ => config::OutputFormat::Detailed,
            }
        } else {
            config::OutputFormat::Minimal
        }
    }
}

impl From<OutputArgs> for config::OutputConfig {
    fn from(value: OutputArgs) -> Self {
        config::OutputConfig {
            format: value.resolved_format(),
            compress: value.output_compress,
            level: value.output_level,
            multifile: value.output_multifile,
        }
    }
}
