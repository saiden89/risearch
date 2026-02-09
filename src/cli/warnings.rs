use anyhow::{bail, Result};
use log::warn;

use risearch::config::{MismatchSpec, SearchArgs, SeedSpec};

fn iter_user_args(args: &[String]) -> impl Iterator<Item = &str> {
    args.iter().skip(1).map(String::as_str)
}

fn matches_long_option(arg: &str, long: &str) -> bool {
    arg == long || arg.strip_prefix(long).is_some_and(|suffix| suffix.starts_with('='))
}

fn matches_short_with_attached_value(arg: &str, short: &str) -> bool {
    arg.starts_with(short) && arg.len() > short.len() && !arg.starts_with("--")
}

fn has_long_option(args: &[String], long: &str) -> bool {
    iter_user_args(args).any(|arg| matches_long_option(arg, long))
}

fn has_any_long_option(args: &[String], longs: &[&str]) -> bool {
    iter_user_args(args).any(|arg| longs.iter().any(|long| matches_long_option(arg, long)))
}

fn has_any_exact_flag(args: &[String], flags: &[&str]) -> bool {
    iter_user_args(args).any(|arg| flags.iter().any(|flag| arg == *flag))
}

fn extract_value_option(args: &[String], short: &str, long: &str) -> Option<String> {
    let mut iter = args.iter().skip(1).map(String::as_str);
    while let Some(arg) = iter.next() {
        if arg == short || arg == long {
            if let Some(val) = iter.next() {
                return Some(val.to_string());
            }
        } else if let Some(val) = arg
            .strip_prefix(long)
            .and_then(|suffix| suffix.strip_prefix('='))
        {
            return Some(val.to_string());
        } else if matches_short_with_attached_value(arg, short) {
            return Some(arg[short.len()..].to_string());
        }
    }
    None
}

fn extract_legacy_mismatch_arg(args: &[String]) -> Option<String> {
    extract_value_option(args, "-m", "--mismatch")
}

fn extract_legacy_seed_arg(args: &[String]) -> Option<String> {
    extract_value_option(args, "-s", "--seed")
}

fn has_seed_pairing_arg(args: &[String]) -> bool {
    let mut iter = args.iter().skip(1).map(String::as_str).peekable();
    while let Some(arg) = iter.next() {
        if arg == "--seed-pairing" {
            if iter.peek().is_some() {
                return true;
            }
        } else if matches_long_option(arg, "--seed-pairing") {
            return true;
        }
    }
    false
}

fn has_no_guseed_arg(args: &[String]) -> bool {
    has_any_exact_flag(args, &["-U", "--no-guseed", "--noGUseed"])
}

fn has_seed_override_args(args: &[String]) -> bool {
    has_any_long_option(args, &["--seed-start", "--seed-end", "--seed-length"])
}

fn has_mismatch_override_args(args: &[String]) -> bool {
    has_any_long_option(
        args,
        &["--mismatch-max", "--mismatch-prefix", "--mismatch-suffix"],
    )
}

fn extract_legacy_report_arg(args: &[String]) -> Option<String> {
    for arg in iter_user_args(args) {
        if arg == "-p" || arg == "--report-alignment" {
            return Some("1".to_string());
        } else if let Some(val) = arg
            .strip_prefix("--report-alignment")
            .and_then(|suffix| suffix.strip_prefix('='))
        {
            return Some(val.to_string());
        } else if matches_short_with_attached_value(arg, "-p") {
            return Some(arg[2..].to_string());
        }
    }
    None
}

fn has_format_arg(args: &[String]) -> bool {
    iter_user_args(args).any(|arg| arg == "-f" || matches_long_option(arg, "--format"))
}

fn extract_legacy_target_flag(args: &[String]) -> Option<&'static str> {
    for arg in iter_user_args(args) {
        if arg == "-i" || matches_short_with_attached_value(arg, "-i") {
            return Some("-i");
        }
        if matches_long_option(arg, "--index") {
            return Some("--index");
        }
    }
    None
}

pub(crate) fn emit_legacy_warnings(raw_args: &[String], opts: &mut SearchArgs) -> Result<()> {
    let legacy_mismatch = extract_legacy_mismatch_arg(raw_args);
    let legacy_seed = extract_legacy_seed_arg(raw_args);
    let legacy_no_guseed = has_no_guseed_arg(raw_args);
    let explicit_pairing = has_seed_pairing_arg(raw_args);
    let has_seed_overrides = has_seed_override_args(raw_args);
    let legacy_report = extract_legacy_report_arg(raw_args);
    let legacy_target_flag = extract_legacy_target_flag(raw_args);
    let has_format = has_format_arg(raw_args);
    let has_mismatch_overrides = has_mismatch_override_args(raw_args);

    if legacy_seed.is_some() && has_seed_overrides {
        bail!(
            "Conflicting seed specification: legacy -s cannot be combined with --seed-start/--seed-end/--seed-length."
        );
    }

    if legacy_mismatch.is_some() && has_mismatch_overrides {
        bail!(
            "Conflicting mismatch specification: legacy -m cannot be combined with --mismatch-max/--mismatch-prefix/--mismatch-suffix."
        );
    }

    let has_seed_start = has_long_option(raw_args, "--seed-start");
    let has_seed_end = has_long_option(raw_args, "--seed-end");
    let has_seed_length = has_long_option(raw_args, "--seed-length");
    let has_any_seed_flag = has_seed_start || has_seed_end || has_seed_length;
    let valid_seed_flags = (!has_any_seed_flag)
        || (!has_seed_start && !has_seed_end && has_seed_length)
        || (has_seed_start && has_seed_end);
    if has_any_seed_flag && !valid_seed_flags {
        bail!(
            "Invalid seed flags: use --seed-length alone, or --seed-start + --seed-end (optionally with --seed-length)."
        );
    }

    if let Some(raw) = legacy_report {
        if has_format {
            warn!(
                "Both legacy -p/--report-alignment ({}) and --format flags were provided; --format takes precedence.",
                raw
            );
        } else {
            let replacement = match raw.as_str() {
                "2" => "cigar",
                "3" => "bindingsite",
                "4" => "minimal",
                _ => "detailed",
            };
            warn!(
                "Legacy report syntax '-p{}' is deprecated; use --format {}.",
                if raw == "1" { "" } else { &raw },
                replacement
            );
        }
    }

    if let Some(flag) = legacy_target_flag {
        if flag == "-i" {
            warn!("Legacy target flag '-i' is deprecated; use -t/--target.");
        } else {
            warn!("Legacy target flag '--index' is deprecated; use --target.");
        }
    }

    if let Some(raw) = legacy_mismatch.as_deref() {
        if opts.seed.has_named_mismatch() {
            warn!(
                "Both legacy -m/--mismatch ({}) and named --mismatch-* flags were provided; named flags take precedence.",
                raw
            );
        }
        let suggestion = raw.parse::<MismatchSpec>().ok().map(|spec| {
            format!(
                "--mismatch-max {} --mismatch-prefix {} --mismatch-suffix {}",
                spec.max_mismatches, spec.min_prefix_matches, spec.min_suffix_matches
            )
        });
        if let Some(s) = suggestion {
            warn!(
                "Legacy mismatch syntax '-m c[:ps[:pe]]' is deprecated; use {}.",
                s
            );
        } else {
            warn!(
                "Legacy mismatch syntax '-m c[:ps[:pe]]' is deprecated; use --mismatch-max/--mismatch-prefix/--mismatch-suffix."
            );
        }
    }

    if let Some(raw) = legacy_seed.as_deref() {
        if has_seed_overrides {
            warn!(
                "Both legacy -s/--seed ({}) and --seed-* flags were provided; --seed-* flags take precedence.",
                raw
            );
        }
        let suggestion = raw.parse::<SeedSpec>().ok().map(|spec| match spec {
            SeedSpec::LengthOnly(len) => format!("--seed-length {}", len),
            SeedSpec::Interval { start, end } => {
                format!("--seed-start {} --seed-end {}", start, end)
            }
            SeedSpec::IntervalWithLength { start, end, length } => {
                format!(
                    "--seed-start {} --seed-end {} --seed-length {}",
                    start, end, length
                )
            }
        });
        if let Some(s) = suggestion {
            warn!("Legacy seed syntax '-s ...' is deprecated; use {}.", s);
        } else {
            warn!(
                "Legacy seed syntax '-s ...' is deprecated; use --seed-start/--seed-end/--seed-length."
            );
        }
    }

    if legacy_no_guseed {
        warn!(
            "Legacy --no-guseed is deprecated; wobble is disabled by default. Use --seed-pairing allow_wobble to enable."
        );
        if explicit_pairing {
            warn!(
                "Both --no-guseed and --seed-pairing were provided; --seed-pairing takes precedence."
            );
        }
    }

    opts.seed.apply_mismatch_overrides();
    opts.seed.apply_pairing_overrides(explicit_pairing);

    Ok(())
}
