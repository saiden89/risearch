#!/bin/bash
# Build C risearch2 binaries for parity testing
# Usage: ./scripts/build_c.sh [release|debug|all]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
C_SRC="$PROJECT_ROOT/legacy_c/RIsearch2"
LIBDIV_BUILD="$C_SRC/libdivsufsort-2.0.1/build"

BUILD_TYPE="${1:-all}"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${YELLOW}Building C risearch2...${NC}"

# Find real GCC (not Apple clang pretending to be gcc)
if [[ "$(uname)" == "Darwin" ]]; then
    # macOS - use homebrew gcc
    for v in 15 14 13 12; do
        if [ -x "/opt/homebrew/bin/gcc-$v" ]; then
            export CC="/opt/homebrew/bin/gcc-$v"
            break
        elif [ -x "/usr/local/bin/gcc-$v" ]; then
            export CC="/usr/local/bin/gcc-$v"
            break
        fi
    done
    if [ -z "$CC" ]; then
        echo "Error: No homebrew GCC found. Install with: brew install gcc"
        exit 1
    fi
else
    export CC="gcc"
fi

echo "Using compiler: $CC"

# Check if libdivsufsort is built
if [ ! -f "$LIBDIV_BUILD/lib/libdivsufsort64.a" ]; then
    echo -e "${YELLOW}Building libdivsufsort-2.0.1...${NC}"
    mkdir -p "$LIBDIV_BUILD"
    cd "$LIBDIV_BUILD"
    cmake -DCMAKE_BUILD_TYPE="Release" \
          -DBUILD_DIVSUFSORT64:BOOL=ON \
          -DUSE_OPENMP:BOOL=ON \
          -DBUILD_SHARED_LIBS:BOOL=OFF \
          ..
    make -j$(sysctl -n hw.ncpu 2>/dev/null || nproc)
    cd "$PROJECT_ROOT"
fi

# Build risearch2
cd "$C_SRC/src"
mkdir -p ../bin

case "$BUILD_TYPE" in
    release)
        echo -e "${YELLOW}Building release binary...${NC}"
        CC=$CC make ../bin/risearch2.x
        echo -e "${GREEN}Built: $C_SRC/bin/risearch2.x${NC}"
        ;;
    debug)
        echo -e "${YELLOW}Building debug binary (with x/y seed markers)...${NC}"
        CC=$CC make ../bin/risearch2.dbg.x
        echo -e "${GREEN}Built: $C_SRC/bin/risearch2.dbg.x${NC}"
        ;;
    all|*)
        echo -e "${YELLOW}Building all binaries...${NC}"
        CC=$CC make
        echo -e "${GREEN}Built:${NC}"
        echo -e "  - risearch2.x (release)"
        echo -e "  - risearch2.dbg.x (debug with x/y markers)"
        ;;
esac

echo -e "${GREEN}Done!${NC}"
