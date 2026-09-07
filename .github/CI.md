# CI operator guide

`ci.yml` runs on pushes and pull requests and gates merges.

## ci.yml — build and test

Runs the workspace build and the test suite on every push and PR.
Uses `dtolnay/rust-toolchain@stable`, `Swatinem/rust-cache@v2`, and `taiki-e/install-action@nextest` (cargo-nextest as the test runner).
A concurrency group cancels superseded runs on the same ref.
Default job permissions are `contents: read`, `actions/checkout` runs with `persist-credentials: false`, and every job sets a `timeout-minutes`.

Cargo validation uses `--locked`. The Python package builds through
`uv run --locked`, and Maturin is configured to use the committed Cargo
lockfile. A manifest change that is not accompanied by its lockfile update
therefore fails the build instead of silently resolving. `Cargo.lock` and
`bindings/python/uv.lock` are both committed for this reason.

The `lint` job runs rustfmt and clippy (`-D warnings`), then two extra guards.
It builds the `openmp` feature (`cargo build --locked -p risearch --features openmp`): it is a documented, user-facing feature that the default-feature test jobs never exercise, so without this it could rot unnoticed.
`build`, not `check`, is deliberate; only a full build links `libgomp`, which is where the native OpenMP failures this lane targets actually surface.
This runs only on Ubuntu because Apple clang ships no OpenMP, so the feature cannot link there with the default toolchain.
It then builds the docs with `cargo doc --locked --no-deps -p risearch --all-features` under `RUSTDOCFLAGS: -D warnings`, so a refactor that breaks an intra-doc link fails here rather than silently degrading the published docs.rs page.

The `actions-lint` job lints the workflow files with `actionlint` (shellcheck-backed), pinned to its `rhysd/actionlint` Docker image by digest.

The `msrv` job pins the declared minimum supported Rust version (`rust-version = "1.88"` in the root `Cargo.toml`, inherited by both crates) and runs `cargo check --locked --workspace --all-targets --all-features` on it.
The floor is 1.88 because `libsais 0.2` uses `let`-chains (stabilised in Rust 1.88) and declares no `rust-version` of its own, so `cargo msrv find` is the source of truth.
Bump both the `rust-version` field and this job's pinned toolchain together when the floor moves.

The `minimal-versions` job proves the declared dependency floors of the published `risearch` crate are real: it resolves every direct dependency down to its declared minimum and builds against that.
`--locked` is intentionally absent because the tool rewrites the lock down to those minimums.
`--ignore-private` drops `risearch-python` (`publish = false`) from the graph: its `pyo3-log` dependency requires `log ~0.4.21`, which cannot unify with the `log = "0.4.8"` floor this job exists to verify, so without it the resolution fails outright and the published crate's floors are never tested.
Where `msrv` fixes the compiler (1.88) and uses the latest deps, this fixes the deps to their floors and uses the current compiler; together they bound the support envelope.

The `test-python` job builds the extension with `maturin develop` into the
locked `uv` environment and runs the Python suite from `bindings/python/`.
`uv sync` prunes packages it does not track, so it must run before
`maturin develop`; `uv run` does not prune, so the suite step leaves the
installed extension in place. The job pins CPython 3.10, the declared floor,
because that is where the Python layer's syntax and typing assumptions break
first — the `cp310-abi3` extension itself is identical on every supported
interpreter. The Rust `test` job excludes `risearch-python`, so this is the
only lane that exercises the binding.

There is no wheel or source-distribution build in CI. Packaging defects that
`maturin develop` cannot surface — module placement, missing Python package or
workspace files, runtime dependency metadata — are therefore unguarded until a
release lane exists.

### OS matrix

The `test` job is a matrix over both runners, reported as `Test ubuntu` and `Test macOS`:

- `ubuntu-latest` (x86_64)
- `macos-latest` (arm64)

Every step is shared except installing the C-oracle system dependencies, which branches on `runner.os`.
`fail-fast: false` keeps one OS from cancelling the other, since an arch-specific parity break is precisely what the second runner is there to catch.

### The C oracle on both runners

Nine integration test binaries depend on the legacy C oracle: they `Command::new` it (directly or via the shared parity runner) and panic if it is absent.
Both test jobs therefore build it from the committed source under `legacy_c/RIsearch2/` via `.github/scripts/build-c-oracle.sh`, then run the FULL suite with no nextest filter.
The script emits `risearch2.x` and the debug `risearch2.dbg.x` into `legacy_c/RIsearch2/bin/`, where the tests auto-detect them; no env var is needed (`PARITY_C_BIN` overrides the choice for local work).
Those binaries are gitignored and built fresh in every run, and the script's own final check exits non-zero unless both are present and executable, so an oracle that silently failed to build cannot masquerade as a passing run.

The oracle needs GCC with OpenMP plus pcre and zlib.
On Ubuntu that is stock: `build-essential`, `cmake`, `libpcre3-dev`, `zlib1g-dev`.
On macOS it is `brew install gcc pcre`, because Apple clang has no OpenMP; the script detects the newest Homebrew `gcc-<N>` and injects the Homebrew pcre include/lib paths, while `-lz` resolves against the SDK.
Running parity on both arches is the point: the Rust and C sides are compared on the same machine, so arm64 and x86_64 exercise independent floating-point paths through the energy model.

The cmake call passes `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` because libdivsufsort's `CMakeLists.txt` declares a pre-3.5 minimum that cmake 4 rejects outright.
Older cmake ignores the variable, so the flag is safe on both runners.

## Action pinning and Dependabot

Official actions from well-governed orgs (GitHub `actions/*`) use a major tag (e.g. `@v7`); every other third-party action (`dtolnay/rust-toolchain`, `Swatinem/rust-cache`, `taiki-e/install-action`, `astral-sh/setup-uv`) is pinned to a full commit SHA with a `# version` comment, and `rhysd/actionlint` to its image digest.
SHA/digest pins are immutable (a moved tag cannot inject code); Dependabot bumps them and keeps the comment current.
`astral-sh/setup-uv` additionally pins the installed `uv` version via its `version:` input, since an unpinned `setup-uv` installs the latest `uv` at run time.

Dependabot tracks three ecosystems: `cargo` (root workspace), `github-actions`
(workflows), and `uv` (`bindings/python`, so `uv.lock` and the PEP 735
`[dependency-groups]` stay current).
