use log::warn;

use risearch::args::SearchArgs;
use risearch::seed::MismatchSpec;

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

fn has_wobble_arg(args: &[String]) -> bool {
    args.iter().skip(1).any(|arg| arg == "-w" || arg == "--wobble")
}

fn has_no_guseed_arg(args: &[String]) -> bool {
    args.iter()
        .skip(1)
        .any(|arg| arg == "-U" || arg == "--no-guseed" || arg == "--noGUseed")
}

pub fn emit_legacy_warnings(raw_args: &[String], opts: &mut SearchArgs) {
    let legacy_mismatch = extract_legacy_mismatch_arg(raw_args);
    let wobble_arg = has_wobble_arg(raw_args);
    let legacy_no_guseed = has_no_guseed_arg(raw_args);
    let explicit_pairing = has_seed_pairing_arg(raw_args);

    if let Some(raw) = legacy_mismatch.as_deref() {
        if opts.seed.has_named_mismatch() {
            warn!(
                "Both legacy -m/--mismatch ({}) and named --mismatch-* flags were provided; named flags take precedence.",
                raw
            );
        }
        let suggestion = raw
            .parse::<MismatchSpec>()
            .ok()
            .map(|spec| {
                format!(
                    "--mismatch-max {} --mismatch-prefix {} --mismatch-suffix {}",
                    spec.max_mismatches, spec.min_position, spec.min_matches_after
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

    if legacy_no_guseed {
        warn!("Legacy --no-guseed is deprecated; use --seed-pairing strict instead.");
        if wobble_arg {
            warn!(
                "Both -w/--wobble and --no-guseed were provided; --no-guseed takes precedence."
            );
        } else if explicit_pairing {
            warn!(
                "Both --no-guseed and --seed-pairing were provided; --seed-pairing takes precedence."
            );
        }
    }

    opts.seed.apply_mismatch_overrides();
    opts.seed.apply_pairing_overrides(explicit_pairing);
}
