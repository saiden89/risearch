use crate::config::{OutputCompression, OutputConfig, OutputFormat};

use anyhow::{bail, Context, Error, Result};
use clap::ValueEnum;
use log::warn;
use std::path::{Path, PathBuf};

/// CLI-facing codec selector — what the `--compress` flag parses into.
#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
#[clap(rename_all = "lowercase")]
pub(crate) enum OutputCodec {
    None,
    #[value(alias = "gz")]
    Gzip,
    #[value(alias = "zst")]
    Zstd,
}

impl From<&std::path::Path> for OutputCodec {
    /// Infer codec from file extension. Unrecognised or absent → `None`.
    fn from(path: &std::path::Path) -> Self {
        match path.extension().and_then(|e| e.to_str()) {
            Some(ext) => match ext.to_ascii_lowercase().as_str() {
                "gz" | "gzip" => Self::Gzip,
                "zst" | "zstd" => Self::Zstd,
                _ => Self::None,
            },
            None => Self::None,
        }
    }
}

/// Reject an output path whose parent directory is missing or is not a directory.
///
/// A CLI-boundary check: it exists so a bad `-o` fails before the work that
/// would fill it, not to guard the write itself.
pub(crate) fn validate_output_parent(path: &Path) -> Result<()> {
    let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) else {
        return Ok(());
    };

    let md = fs_err::metadata(parent)
        .with_context(|| format!("output directory '{}' does not exist", parent.display()))?;
    if !md.is_dir() {
        bail!("output path '{}' is not a directory", parent.display());
    }
    Ok(())
}

/// Boundary CLI arguments for output destination, formatting, and compression.
#[derive(clap::Args, Debug, Clone)]
pub(crate) struct OutputArgs {
    /// Output file for search results (use '-' for stdout)
    #[arg(short = 'o', long = "output", value_name = "FILE", default_value = "-")]
    pub(crate) path: PathBuf,

    /// Output format
    #[arg(
        short = 'f',
        long = "format",
        value_name = "FORMAT",
        default_missing_value = "minimal",
        value_enum
    )]
    pub(crate) report_format: Option<OutputFormat>,

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
    pub(crate) report_legacy: Option<OutputFormat>,

    /// Output compression codec (overrides file extension inference; gzip/gz, zstd/zst accepted)
    #[arg(long = "compress", value_enum)]
    pub(crate) output_compress: Option<OutputCodec>,

    /// Output compression level (codec-specific: gzip 0–9, zstd -7..22)
    #[arg(long = "compress-level", value_name = "LEVEL")]
    pub(crate) output_level: Option<i32>,

    /// Write one output file per query into the directory given by -o
    #[arg(long = "multifile", action = clap::ArgAction::SetTrue)]
    pub(crate) output_multifile: bool,
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
            validate_output_parent(&value.path)?;
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

        let format = match (value.report_format, value.report_legacy) {
            (Some(fmt), Some(_)) => {
                warn!("Both legacy -p/--report-alignment and --format were provided; --format takes precedence.");
                fmt
            }
            (Some(fmt), None) => fmt,
            (None, Some(fmt)) => {
                warn!(
                    "Legacy -p/--report-alignment is deprecated; use --format {}.",
                    match fmt {
                        OutputFormat::Cigar => "cigar",
                        OutputFormat::BindingSite => "bindingsite",
                        OutputFormat::Minimal => "minimal",
                        OutputFormat::Detailed => "detailed",
                    }
                );
                fmt
            }
            (None, None) => OutputFormat::default(),
        };

        Ok(OutputConfig {
            format,
            compress,
            multifile: value.output_multifile,
        })
    }
}
