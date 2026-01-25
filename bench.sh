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
RUSTFLAGS="-C target-cpu=native" cargo build --release 2>&1

echo ""
echo "=============================================="
echo "EXACT MATCH BENCHMARKS"
echo "=============================================="

# for SEED in 6 8 10; do
#     echo "--- Seed length $SEED, max_ext 20 ---"
#     echo "[Rust]"
#     time $RUST_BIN search -q $MIRNAS -i $RUST_IDX -s $SEED -l 20 -e 100 -o $OUT --output-compress gzip -t 1 --format=bindingsite 2>&1
#     echo "[C]"
#     time $C_BIN -q $MIRNAS -i $C_IDX -s $SEED -l 20 -e 100 -t 1 -p3 > $OUT 2>&1
#     echo ""
# done

echo ""
echo "=============================================="
echo "MISMATCH BENCHMARKS"
echo "=============================================="

for MISMATCH in "1:0" "1:3" "2:2"; do
    for SEED in 8; do
        echo ""
        echo "--- Mismatch $MISMATCH, Seed $SEED, max_ext 10 ---"
        echo "[Rust]"
        time $RUST_BIN search -q $MIRNAS -i $RUST_IDX -s $SEED -l 10 -m $MISMATCH -e 100 -o $OUT -t 2 --output-compress gzip --format=bindingsite  2>&1
        echo "[C]" 2>&1
        echo "[C]"
        time $C_BIN -q $MIRNAS -i $C_IDX -s $SEED -l 10 -m $MISMATCH -e 100 -t 2 -p3 > $OUT 2>&1
    done
done


echo ""
echo "Done!"
