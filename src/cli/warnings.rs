use log::warn;

use risearch::config::SearchArgs;
use risearch::seed::{MismatchSpec, SeedSpec};

fn extract_legacy_mismatch_arg(args: &[String]) -> Option<String> {
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        if arg == "-m" || arg == "--mismatch" {
            if let Some(val) = iter.next() {
                return Some(val.clone());
            }
        } else if let Some(val) = arg.strip_prefix("--mismatch=") {
            return Some(val.to_string());
        } else if arg.starts_with("-m") && arg.len() > 2 {
            return Some(arg[2..].to_string());
        }
    }
    None
}

fn extract_legacy_seed_arg(args: &[String]) -> Option<String> {
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        if arg == "-s" || arg == "--seed" {
            if let Some(val) = iter.next() {
                return Some(val.clone());
            }
        } else if let Some(val) = arg.strip_prefix("--seed=") {
            return Some(val.to_string());
        } else if arg.starts_with("-s") && arg.len() > 2 {
            return Some(arg[2..].to_string());
        }
    }
    None
}

fn has_seed_pairing_arg(args: &[String]) -> bool {
    let mut iter = args.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        if arg == "--seed-pairing" {
            if iter.peek().is_some() {
                return true;
            }
        } else if arg.starts_with("--seed-pairing=") {
            return true;
        }
    }
    false
}

fn has_no_guseed_arg(args: &[String]) -> bool {
    args.iter()
        .skip(1)
        .any(|arg| arg == "-U" || arg == "--no-guseed" || arg == "--noGUseed")
}

fn has_seed_override_args(args: &[String]) -> bool {
    args.iter().skip(1).any(|arg| {
        arg == "--seed-start"
            || arg.starts_with("--seed-start=")
            || arg == "--seed-end"
            || arg.starts_with("--seed-end=")
            || arg == "--seed-length"
            || arg.starts_with("--seed-length=")
    })
}

fn extract_legacy_report_arg(args: &[String]) -> Option<String> {
    for arg in args.iter().skip(1) {
        if arg == "-p" || arg == "--report-alignment" {
            return Some("1".to_string());
        } else if let Some(val) = arg.strip_prefix("--report-alignment=") {
            return Some(val.to_string());
        } else if arg.starts_with("-p") && arg.len() > 2 && !arg.starts_with("--") {
            return Some(arg[2..].to_string());
        }
    }
    None
}

fn has_format_arg(args: &[String]) -> bool {
    args.iter()
        .skip(1)
        .any(|arg| arg == "-f" || arg == "--format" || arg.starts_with("--format="))
}

pub fn emit_legacy_warnings(raw_args: &[String], opts: &mut SearchArgs) {
    let legacy_mismatch = extract_legacy_mismatch_arg(raw_args);
    let legacy_seed = extract_legacy_seed_arg(raw_args);
    let legacy_no_guseed = has_no_guseed_arg(raw_args);
    let explicit_pairing = has_seed_pairing_arg(raw_args);
    let has_seed_overrides = has_seed_override_args(raw_args);
    let legacy_report = extract_legacy_report_arg(raw_args);
    let has_format = has_format_arg(raw_args);

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
            SeedSpec::Length(len) => format!("--seed-length {}", len),
            SeedSpec::Interval { start, end, length } => {
                let mut s = format!("--seed-start {} --seed-end {}", start, end);
                if let Some(l) = length {
                    s.push_str(&format!(" --seed-length {}", l));
                }
                s
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
}
