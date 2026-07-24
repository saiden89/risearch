# CI operator guide

`ci.yml` runs on pushes and pull requests and gates merges.

## ci.yml — build and test

Runs the workspace build and the test suite on every push and PR.
Uses `dtolnay/rust-toolchain@stable`, `Swatinem/rust-cache@v2`, and `taiki-e/install-action@nextest` (cargo-nextest as the test runner).
A concurrency group cancels superseded runs on the same ref.
Default job permissions are `contents: read`, `actions/checkout` runs with `persist-credentials: false`, and every job sets a `timeout-minutes`.

Every dependency-resolving cargo and `uv` invocation passes `--locked` so a manifest change that is not accompanied by its lockfile update fails CI instead of silently resolving.
`Cargo.lock` and `risearch-python/uv.lock` are both committed for this reason.

The `lint` job runs rustfmt and clippy (`-D warnings`), then two extra guards.
It builds the `openmp` feature (`cargo build --locked -p risearch --features openmp`): it is a documented, user-facing feature that the default-feature test jobs never exercise, so without this it could rot unnoticed.
`build`, not `check`, is deliberate; only a full build links `libgomp`, which is where the native OpenMP failures this lane targets actually surface.
This runs only on Ubuntu because OpenMP under Apple clang is fragile; that is also why the C oracle is Linux-only.
It then builds the docs with `RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps -p risearch --all-features`, so a refactor that breaks an intra-doc link fails here rather than silently degrading the published docs.rs page.

The `actions-lint` job lints the workflow files with `actionlint` (shellcheck-backed), pinned to its `rhysd/actionlint` Docker image by digest.

The `msrv` job pins the declared minimum supported Rust version (`rust-version = "1.88"` in the root `Cargo.toml`, inherited by both crates) and runs `cargo check --locked --workspace --all-targets --all-features` on it.
The floor is 1.88 because `libsais 0.2` uses `let`-chains (stabilised in Rust 1.88) and declares no `rust-version` of its own, so `cargo msrv find` is the source of truth.
Bump both the `rust-version` field and this job's pinned toolchain together when the floor moves.

The `minimal-versions` job proves the declared dependency floors of the published `risearch` crate are real: it resolves every direct dependency down to its declared minimum and builds against that (`cargo minimal-versions check --direct --no-dev-deps -p risearch`, via `taiki-e/install-action`).
`--locked` is intentionally absent because the tool rewrites the lock down to those minimums.
Where `msrv` fixes the compiler (1.88) and uses the latest deps, this fixes the deps to their floors and uses the current compiler; together they bound the support envelope.

The `test-wheel` job (matrix over OS) builds a real abi3 wheel, installs it into a clean venv, and runs the suite against the installed package.
This catches packaging defects `maturin develop` cannot (module-name, `__init__.py`/`.pyi` inclusion, runtime dependency resolution).
Because the wheel is abi3 a single binary per OS works across every Python version, so testing one Python per OS is sufficient.
pytest runs from the repository root, not `risearch-python/`, so `import risearch` resolves to the installed wheel rather than the source package (which also carries a stale committed `_risearch.abi3.so`).

### OS matrix

- `ubuntu-latest` (x86_64)
- `macos-latest` (arm64)

### Why the C oracle is Linux-only

Eight integration test binaries depend on the legacy C oracle: they `Command::new` it (directly or via the shared parity runner) and panic if it is absent.
Six are the core differential parity binaries.
Two more are regression binaries that also drive the oracle: `parity_mismatch_regression` hard-asserts `risearch2.x` exists, and `parity_seed_prune_regression` runs the oracle through the shared `SingleSeqRunner`.
All eight must therefore be filtered out where the oracle is not built.
The oracle is built from committed source under `legacy_c/RIsearch2/` with GCC + OpenMP, which is clean on Ubuntu.
On macOS the toolchain is Apple clang (no OpenMP) and the brew GCC / PCRE stack drifts, so the oracle build is fragile and is deliberately not attempted there.

The Linux job builds both oracle binaries to `legacy_c/RIsearch2/bin/` (`risearch2.x` and the debug `risearch2.dbg.x`, both auto-detected; no env var needed) and runs the FULL suite including parity.
The macOS job runs everything EXCEPT the eight oracle-dependent binaries, filtered out via nextest with exact `binary(=...)` matchers.
The skip is logged with its reason (macOS has no C oracle) so it is never silent.

Excluded on macOS:

- `parity_seed`
- `parity_seed_prune_regression`
- `parity_extension`
- `parity_dsm`
- `parity_mismatch`
- `parity_mismatch_regression`
- `parity_wobble`
- `parity_debug`

Ubuntu system packages needed for the oracle build: `build-essential`, `cmake`, `libpcre3-dev`, `zlib1g-dev`.

## Action pinning and Dependabot

Official actions from well-governed orgs (GitHub `actions/*`) use a major tag (e.g. `@v7`); every other third-party action (`dtolnay/rust-toolchain`, `Swatinem/rust-cache`, `taiki-e/install-action`, `astral-sh/setup-uv`) is pinned to a full commit SHA with a `# version` comment, and `rhysd/actionlint` to its image digest.
SHA/digest pins are immutable (a moved tag cannot inject code); Dependabot bumps them and keeps the comment current.
`astral-sh/setup-uv` additionally pins the installed `uv` version via its `version:` input, since an unpinned `setup-uv` installs the latest `uv` at run time.

Dependabot tracks three ecosystems: `cargo` (root workspace), `github-actions` (workflows), and `uv` (`risearch-python`, so `uv.lock` and the PEP 735 `[dependency-groups]` stay current).
