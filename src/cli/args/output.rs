use crate::config::{self, OutputCodec, OutputCompression, OutputFormat};
use config::OutputConfig;

use anyhow::{bail, Error};
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
        value_parser = parse_legacy_format,
        help_heading = "Deprecated"
    )]
    pub report_legacy: Option<OutputFormat>,

    /// Output compression codec (overrides file extension inference; gzip/gz, zstd/zst accepted)
    #[arg(long = "compress", value_enum)]
    pub output_compress: Option<OutputCodec>,

    /// Output compression level (codec-specific: gzip 0–9, zstd -7..22)
    #[arg(long = "compress-level", value_name = "LEVEL")]
    pub output_level: Option<i32>,

    /// Write one output file per query into the directory given by -o
    #[arg(long = "multifile", action = clap::ArgAction::SetTrue)]
    pub output_multifile: bool,
}

fn parse_legacy_format(s: &str) -> Result<OutputFormat, String> {
    match s.parse::<u8>().map_err(|e| e.to_string())? {
        1 => Ok(OutputFormat::Detailed),
        2 => Ok(OutputFormat::Cigar),
        3 => Ok(OutputFormat::BindingSite),
        4 => Ok(OutputFormat::Minimal),
        n => Err(format!("unknown format mode {n}, expected 1–4")),
    }
}

impl TryFrom<OutputArgs> for OutputConfig {
    type Error = Error;

    fn try_from(value: OutputArgs) -> Result<Self, Self::Error> {
        if value.output_multifile && value.path.as_os_str() == "-" {
            bail!("--multifile requires -o/--output to be a directory path; '-' (stdout) is not allowed.");
        }

        if value.path.as_os_str() != "-" {
            if let Some(parent) = value.path.parent() {
                if !parent.as_os_str().is_empty() && !parent.exists() {
                    bail!("output directory '{}' does not exist", parent.display());
                }
            }
        }

        let codec = value
            .output_compress
            .unwrap_or_else(|| OutputCodec::from(value.path.as_path()));

        if let Some(level) = value.output_level {
            match codec {
                OutputCodec::None => bail!("--compress-level requires compressed output"),
                OutputCodec::Gzip if !(0..=9).contains(&level) => {
                    bail!("gzip level must be 0-9 (got {})", level)
                }
                OutputCodec::Zstd if !(-7..=22).contains(&level) => {
                    bail!("zstd level must be -7..22 (got {})", level)
                }
                _ => {}
            }
        }

        let compress = match (codec, value.output_level) {
            (OutputCodec::None, _) => OutputCompression::None,
            (OutputCodec::Gzip, lvl) => OutputCompression::Gzip(lvl.unwrap_or(6) as u8),
            (OutputCodec::Zstd, lvl) => OutputCompression::Zstd(lvl.unwrap_or(3)),
        };

        Ok(OutputConfig {
            format: value
                .report_format
                .or(value.report_legacy)
                .unwrap_or_default(),
            compress,
            multifile: value.output_multifile,
        })
    }
}
