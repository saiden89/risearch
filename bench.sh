#!/bin/bash
# Quick performance comparison: Rust vs C RIsearch
# Usage: ./bench.sh

set -e

MIRNAS="legacy_c/RIsearch2/test_suite/mirnas.fa"
RUST_IDX="chr22.idx"
C_IDX="chr22_c.idx"
RUST_BIN="./target/release/risearch"
C_BIN="legacy_c/RIsearch2/bin/risearch2.x"
OUT="/tmp/bench_out"

# Build release
echo "Building Rust release..."
cargo build --release 2>&1 | tail -1

echo ""
echo "=============================================="
echo "EXACT MATCH BENCHMARKS"
echo "=============================================="

for SEED in 6 8 10; do
    echo ""
    echo "--- Seed length $SEED, max_ext 20 ---"
    echo "[Rust]"
    time $RUST_BIN search -q $MIRNAS -i $RUST_IDX -s $SEED -l 20 -e 100 -o $OUT 2>&1
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
        echo "[Rust]"
        time $RUST_BIN search -q $MIRNAS -i $RUST_IDX -s $SEED -l 10 -m $MISMATCH -e 100 -o $OUT 2>&1
        echo "[C] (1 thread)"
        time $C_BIN -q $MIRNAS -i $C_IDX -s $SEED -l 10 -m $MISMATCH -e 100 -t 1 > $OUT 2>&1
    done
done

echo ""
echo "=============================================="
echo "HEAVY MISMATCH (3:0)"
echo "=============================================="

echo ""
echo "--- Mismatch 3:0, Seed 6, max_ext 10 ---"
echo "[Rust]"
time $RUST_BIN search -q $MIRNAS -i $RUST_IDX -s 6 -l 10 -m 3:0 -e 100 -o $OUT 2>&1
echo "[C] (1 thread)"
time $C_BIN -q $MIRNAS -i $C_IDX -s 6 -l 10 -m 3:0 -e 100 -t 1 > $OUT 2>&1

echo ""
echo "Done!"
