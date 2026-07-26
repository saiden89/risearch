#!/usr/bin/env bash
#
# build-c-oracle.sh — build the legacy RIsearch2 C oracle for parity tests.
#
# Produces both parity-test binaries under legacy_c/RIsearch2/bin/:
#   - risearch2.x       (optimized)
#   - risearch2.dbg.x   (debug; the one parity tests prefer)
#
# The Rust parity integration tests (tests/parity_*.rs) auto-detect these
# binaries at legacy_c/RIsearch2/bin/ and PANIC when they are absent, so this
# script must leave both present and executable.
#
# Build steps mirror legacy_c/RIsearch2/rebuild.sh (RIsearch2 half only):
#   1. build libdivsufsort64 as a Release static lib via cmake, and
#   2. compile the RIsearch2 sources against it.
#
# Notes on the toolchain:
#   - src/Makefile hardcodes `CC = gcc-15`; we override it with `make CC=...`.
#     Honor a caller-supplied $CC, else default to plain `gcc`.
#   - The Makefile links -fopenmp / -lz / -lpcre / -lm, so the host needs a
#     GCC with OpenMP plus pcre and zlib. On Linux: build-essential, cmake,
#     libpcre3-dev, zlib1g-dev. On macOS: `brew install gcc pcre` — Apple clang
#     has no OpenMP, so this picks the newest Homebrew `gcc-<N>` and adds the
#     Homebrew pcre paths; -lz resolves against the SDK.
#
# Safe to re-run: the libdivsufsort build dir is recreated from scratch and the
# RIsearch2 objects are `make clean`ed before each build.

set -euo pipefail

# --- helpers ---------------------------------------------------------------

section() {
	printf '\n==> %s\n' "$*"
}

die() {
	printf 'ERROR: %s\n' "$*" >&2
	exit 1
}

# --- toolchain --------------------------------------------------------------
# Exported so cmake configures libdivsufsort with the same compiler that links
# the RIsearch2 binaries. CFLAGS/CLIBS are appended to by the Makefile's `+=`.

if [ "$(uname -s)" = "Darwin" ]; then
	command -v brew >/dev/null 2>&1 || die "Homebrew not found; needed for gcc and pcre"
	BREW_PREFIX="$(brew --prefix)"

	if [ -z "${CC:-}" ]; then
		CC="$(find "$BREW_PREFIX/bin" -maxdepth 1 -name 'gcc-[0-9]*' 2>/dev/null |
			grep -E '/gcc-[0-9]+$' | sort -V | tail -1)"
		[ -n "$CC" ] || die "no Homebrew gcc-<N> on PATH; run: brew install gcc"
	fi

	PCRE_PREFIX="$(brew --prefix pcre 2>/dev/null)" ||
		die "Homebrew pcre not installed; run: brew install pcre"
	export CFLAGS="-I$PCRE_PREFIX/include ${CFLAGS:-}"
	export CLIBS="-L$PCRE_PREFIX/lib ${CLIBS:-}"
fi

export CC="${CC:-gcc}"

# --- resolve paths ----------------------------------------------------------
# Resolve everything relative to this script so it works from any CWD.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

ORACLE_DIR="$REPO_ROOT/legacy_c/RIsearch2"
DSS_DIR="$ORACLE_DIR/libdivsufsort-2.0.1"
DSS_BUILD_DIR="$DSS_DIR/build"
SRC_DIR="$ORACLE_DIR/src"
BIN_DIR="$ORACLE_DIR/bin"

BIN_OPT="$BIN_DIR/risearch2.x"
BIN_DBG="$BIN_DIR/risearch2.dbg.x"

# Pick a sensible parallelism level, portable across Linux and macOS.
JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)"

# --- preflight --------------------------------------------------------------

section "RIsearch2 C oracle build"
echo "repo root      : $REPO_ROOT"
echo "oracle dir     : $ORACLE_DIR"
echo "C compiler (CC): $CC"
echo "parallel jobs  : $JOBS"

[ -d "$ORACLE_DIR" ] || die "oracle directory not found: $ORACLE_DIR"
[ -d "$DSS_DIR" ]    || die "libdivsufsort directory not found: $DSS_DIR"
[ -d "$SRC_DIR" ]    || die "RIsearch2 src directory not found: $SRC_DIR"

command -v cmake >/dev/null 2>&1 || die "cmake not found on PATH"
command -v make  >/dev/null 2>&1 || die "make not found on PATH"
command -v "$CC" >/dev/null 2>&1 || die "C compiler '$CC' not found on PATH"

# --- 1. build libdivsufsort64 (Release static lib) --------------------------

section "Building libdivsufsort64 (fresh build dir)"
rm -rf "$DSS_BUILD_DIR"
mkdir -p "$DSS_BUILD_DIR"

cmake \
	-S "$DSS_DIR" \
	-B "$DSS_BUILD_DIR" \
	-DCMAKE_BUILD_TYPE=Release \
	-DBUILD_DIVSUFSORT64:BOOL=ON \
	-DUSE_OPENMP:BOOL=ON \
	-DBUILD_SHARED_LIBS:BOOL=OFF \
	-DCMAKE_POLICY_VERSION_MINIMUM=3.5

cmake --build "$DSS_BUILD_DIR" --parallel "$JOBS"

DSS_LIB="$DSS_BUILD_DIR/lib/libdivsufsort64.a"
[ -f "$DSS_LIB" ] || die "libdivsufsort64 static lib not produced: $DSS_LIB"
echo "built static lib: $DSS_LIB"

# --- 2. build RIsearch2 binaries -------------------------------------------

section "Building RIsearch2 binaries (target: all)"
mkdir -p "$BIN_DIR"

# `make clean` removes any stale binaries so re-runs relink cleanly.
# The default `all` target builds BOTH risearch2.x and risearch2.dbg.x.
# CC=... on the command line overrides the Makefile's hardcoded gcc-15.
make -C "$SRC_DIR" CC="$CC" clean
make -C "$SRC_DIR" CC="$CC" -j"$JOBS" all

# --- 3. verify outputs ------------------------------------------------------

section "Verifying oracle binaries"

missing=0
for bin in "$BIN_OPT" "$BIN_DBG"; do
	if [ -x "$bin" ]; then
		echo "OK  $bin"
	else
		echo "MISSING/NOT-EXECUTABLE  $bin" >&2
		missing=1
	fi
done

[ "$missing" -eq 0 ] || die "one or more oracle binaries were not built"

section "C oracle build complete"
echo "$BIN_OPT"
echo "$BIN_DBG"
