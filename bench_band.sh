#!/bin/bash
# Experimental banding benchmark: Rust banded DP vs C baseline
# Usage: BAND=8 MODE=hard ./bench_band.sh

set -e

MIRNAS="legacy_c/RIsearch2/test_suite/mirnas.fa"
RUST_IDX="chr22.idx"
C_IDX="chr22_c.idx"
RUST_BIN="./target/release/risearch"
C_BIN="legacy_c/RIsearch2/bin/risearch2.x"
OUT="/tmp/bench_out"

BAND="${BAND:-8}"
MODE="${MODE:-hard}" # hard|soft
BUILD="${BUILD:-1}"

if [[ "$BUILD" == "1" ]]; then
  echo "Building Rust release..."
  cargo build --release 2>&1 | tail -1
fi

if [[ ! -f "$RUST_IDX" ]]; then
  cargo run -- index chr22.fa chr22.idx
fi

echo ""
echo "=============================================="
echo "BANDED DP BENCHMARKS (band=$BAND, mode=$MODE)"
echo "=============================================="

for SEED in 6 8 10; do
    echo ""
    echo "--- Seed length $SEED, max_ext 20 ---"
    echo "[Rust banded]"
    time $RUST_BIN search -q $MIRNAS -i $RUST_IDX -s $SEED -l 20 -e 100 -o $OUT -t 1 \
        --experimental --dp-band $BAND --dp-band-mode $MODE 2>&1
    echo "[C] (1 thread)"
    time $C_BIN -q $MIRNAS -i $C_IDX -s $SEED -l 20 -e 100 -t 1 > $OUT 2>&1
done

echo ""
echo "=============================================="
echo "MISMATCH BENCHMARKS"
echo "=============================================="

for MISMATCH in "1:0" "1:3" "2:2"; do
    for SEED in 6 8; do
        echo ""
        echo "--- Mismatch $MISMATCH, Seed $SEED, max_ext 10 ---"
        echo "[Rust banded]"
        time $RUST_BIN search -q $MIRNAS -i $RUST_IDX -s $SEED -l 10 -m $MISMATCH -e 100 -o $OUT -t 1 \
            --experimental --dp-band $BAND --dp-band-mode $MODE 2>&1
        echo "[C] (1 thread)"
        time $C_BIN -q $MIRNAS -i $C_IDX -s $SEED -l 10 -m $MISMATCH -e 100 -t 1 > $OUT 2>&1
    done
done

echo ""
echo "Done!"
