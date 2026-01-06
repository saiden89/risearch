
# Optimizations from Legacy C (RIsearch2)

Analysis of `legacy_c/RIsearch2/src/search.c` and `sa.h` reveals several optimization techniques used in the C implementation that differ from the current Rust approach.

## 1. Suffix Array "Packing" (Cache Locality)

**Concept:**
Standard Suffix Arrays store indices $p$ into the text. To check a character during search, you must read `Text[p]`. This causes random memory access (cache misses) across the huge `Text` array.

**C Implementation (`sa.h`):**
RIsearch2 "packs" the character at `Text[p]` *directly* into the 64-bit Suffix Array entry.

```c
// sa.h
/* SA is stored in bits 34-0 */
#define XSA(_x) ((_x)&0x00000003ffffffff)
/* Character is stored in bits 38-34 */
#define XRIS(_x) (((_x)&0x0000001c00000000)>>34)
```

* **Bit 0-33:** The Suffix Array Index (positions up to ~17 billion).
* **Bit 34-38:** The character code (A, C, G, U, N).
* **Result:** When performing binary search or traversal, the algorithm reads the 64-bit integer and immediately knows the character *without* dereferencing the text pointer.

**Rust Implementation:**
Currently likely uses `Vec<u64>` or `Vec<usize>` for SA values and looks up characters in `Vec<u8>` text.

* **Optimization Opportunity:** Implement a custom Struct or bit-packed `u64` for the SA to include the first character (or first 2 characters) of the suffix.

## 2. Sentinel Values vs Option Types (Memory Density)

**Concept:**
Rust's `Option<i32>` is safe but consumes 8 bytes (4 bytes data + 4 bytes tag/alignment). For a DP matrix, this doubles memory bandwidth usage.

**C Implementation (`search.c`):**
Uses a sentinel value (`INT_MIN / 2`) to represent "None"/"Not Reachable".

```c
#define NA INT_MIN / 2
// ...
int *M = (int *)calloc(dp_size, sizeof(int)); // 4 bytes per cell
M[TI(0, 0)] = NA;
```

* **Result:** 2x better cache density for DP matrices.

**Rust Implementation:**
Uses `Grid<Option<i32>>`.

* **Optimization Opportunity:** Switch to `Grid<i32>` and use a constant like `i32::MIN / 2` to represent "Empty".

## 3. Subseed Filtering Heuristic

**Concept:**
Sometimes a long seed has a poor overall energy, but contains a short "subseed" with very good energy. Conversely, a seed might pass the threshold but only because of one good region.

**C Implementation (`sa_evaluate_interval`):**
If `seed_threshold_flag` is set, it performs a sliding window check over every seed match:

```c
// search.c:1877
for (j = seed_len - 1; j >= (min_seed_length - 1); j--) {
  // ... efficient O(N) scan to find best sub-seed energy ...
}
if (best_perlength_score < threshold) continue;
```

It effectively maximizes the "Energy Per Nucleotide" metric to filter out weak seeds early.

## 4. Flattened 1D Array Access

**Concept:**
Multi-dimensional array access `M[i][j]` can involve pointer indirection.

**C Implementation:**
Explicitly manages 1D flat arrays with macros:

```c
#define TI(q, t) ((q) * (lt) + (t))
M[TI(i, j)] = ...
```

**Rust Implementation:**
Already does this in `Grid<T>`, so parity is good here.

## 5. "Double Walking" Suffix Trees

**Concept:**
`sa_parallel_match_neg` recursively steps through *both* the Query Suffix Array and Target Suffix Array simultaneously.
It uses small arrays `qint[6]` and `sint[6]` to store the intervals for A,C,G,U,N at the current depth.
It then recurses only on matching or valid wobble pairs.

```c
// e.g. Match G with C
if (qint[3] - qint[2]) { // Query has G
  if (sint[2] - sint[1]) // Target has C
     recurse(...)
}
```

This is likely similar to the Rust `search_sa_simple` stack approach but ensures very explicit handling of mismatches/wobbles at the tree level.
