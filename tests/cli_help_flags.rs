//! Regression test for the auto-generated `--help`/`-h` flags.
//!
//! These were silently broken when the `clap` dependency was built with
//! `default-features = false, features = ["std"]`: dropping clap's default
//! `help` feature means the `--help`/`-h` flags are never registered, so they
//! fall through to the "unexpected argument" path and exit with code 2 instead
//! of printing help. The fix re-enables the `help` (+ `usage`) features.
//!
//! Asserting on both the top level and each subcommand guards against a
//! regression at any single level of the command tree.

use assert_cmd::cargo::cargo_bin_cmd;
use predicates::prelude::*;

/// Top-level `--help` / `-h` must exit 0 and print usage + the command list.
#[test]
fn top_level_help_flags_print_usage() {
    for flag in ["--help", "-h"] {
        cargo_bin_cmd!("risearch")
            .arg(flag)
            .assert()
            .success()
            .stdout(predicate::str::contains("Usage: risearch"))
            .stdout(predicate::str::contains("Commands:"))
            .stdout(predicate::str::contains("search"));
    }
}

/// `search --help` / `search -h` must exit 0 and print the subcommand usage.
/// `search` flattens several arg groups, so this also proves none of them
/// swallow the help flag.
#[test]
fn search_subcommand_help_flags_print_usage() {
    for flag in ["--help", "-h"] {
        cargo_bin_cmd!("risearch")
            .arg("search")
            .arg(flag)
            .assert()
            .success()
            .stdout(predicate::str::contains("Usage: risearch search"))
            .stdout(predicate::str::contains("--query"));
    }
}

/// `index --help` must exit 0 and print the subcommand usage. `index` has only
/// plain positional args, so it isolates the failure to help-flag generation
/// rather than any per-arg parsing config.
#[test]
fn index_subcommand_help_flag_prints_usage() {
    cargo_bin_cmd!("risearch")
        .arg("index")
        .arg("--help")
        .assert()
        .success()
        .stdout(predicate::str::contains("Usage: risearch index"))
        .stdout(predicate::str::contains("<INPUT>"));
}
