#!/bin/bash
# Profile script for comparing C and Rust RIsearch implementations
# Usage: ./profile_comparison.sh

set -e

# Test data paths (adjust as needed)
QUERY="test_data/hsa-miR-24-3p.fa"
INDEX="test_data/test.idx"
OUTPUT_DIR="profiling_results"

mkdir -p "$OUTPUT_DIR"

echo "=== Profiling C implementation ==="
# Use macOS sample for quick profiling
sample legacy_c/RIsearch2/bin/risearch2.prof.x 5 -wait -file "$OUTPUT_DIR/c_sample.txt" &
C_PID=$!

# Run C binary (adjust arguments as needed)
legacy_c/RIsearch2/bin/risearch2.prof.x -q "$QUERY" -i "$INDEX" -o /dev/null -s 6 2>&1 | head -5

wait $C_PID 2>/dev/null || true

echo ""
echo "=== Profiling Rust implementation ==="
# Build Rust with profiling symbols
cargo build --release

sample target/release/risearch 5 -wait -file "$OUTPUT_DIR/rust_sample.txt" &
R_PID=$!

# Run Rust binary
target/release/risearch -q "$QUERY" -i "$INDEX" -o /dev/null -s 6 2>&1 | head -5

wait $R_PID 2>/dev/null || true

echo ""
echo "=== Results ==="
echo "C profile: $OUTPUT_DIR/c_sample.txt"
echo "Rust profile: $OUTPUT_DIR/rust_sample.txt"
echo ""
echo "To view, open in a text editor or use:"
echo "  open $OUTPUT_DIR/c_sample.txt"
echo "  open $OUTPUT_DIR/rust_sample.txt"
echo ""
echo "For interactive flamegraph, use Instruments.app:"
echo "  instruments -t 'Time Profiler' legacy_c/RIsearch2/bin/risearch2.prof.x [args]"
