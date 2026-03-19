use anyhow::Result;
use log::warn;

use risearch::cli::args::SearchArgs;
use risearch::config::SeedSpec;

fn warn_if(used: bool, flag: &str, replacement: &str) {
    if used {
        warn!("'{}' is deprecated; use {} instead.", flag, replacement);
    }
}

pub(crate) fn emit_legacy_warnings(args: &SearchArgs, legacy_target: bool) -> Result<()> {
    // Simple flag deprecations
    warn_if(args.seed.no_guseed_legacy, "--noGUseed", "--no-seed-wobble");
    warn_if(legacy_target, "-i", "-t/--target");

    // -p/--report-alignment: warn with smart suggestion
    if let Some(fmt) = &args.output.report_legacy {
        if args.output.report_format.is_some() {
            warn!(
                "Both legacy -p/--report-alignment and --format were provided; --format takes precedence."
            );
        } else {
            use risearch::config::OutputFormat;
            let replacement = match fmt {
                OutputFormat::Cigar => "cigar",
                OutputFormat::BindingSite => "bindingsite",
                OutputFormat::Minimal => "minimal",
                OutputFormat::Detailed => "detailed",
            };
            warn!("Legacy -p/--report-alignment is deprecated; use --format {}.", replacement);
        }
    }

    // -m/--mismatch: warn with smart suggestion
    if let Some(spec) = &args.seed.mismatch_legacy {
        let s = spec.0;
        warn!(
            "Legacy -m/--mismatch is deprecated; use --mismatch-max {} --mismatch-prefix {} --mismatch-suffix {}.",
            s.max_mismatches, s.min_prefix_matches, s.min_suffix_matches
        );
    }

    // -s/--seed: warn with smart suggestion
    if let Some(spec) = &args.seed.seed_legacy {
        let suggestion = match spec.0 {
            SeedSpec::LengthOnly(len) => format!("--seed-length {}", len),
            SeedSpec::Interval { start, end, length: None } => {
                format!("--seed-start {} --seed-end {}", start, end)
            }
            SeedSpec::Interval { start, end, length: Some(length) } => {
                format!("--seed-start {} --seed-end {} --seed-length {}", start, end, length)
            }
        };
        warn!("Legacy -s/--seed is deprecated; use {}.", suggestion);
    }

    Ok(())
}
