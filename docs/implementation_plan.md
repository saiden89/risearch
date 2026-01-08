# RIsearch Implementation Plan

## Table of Contents

- [Overview](#overview)
- [Completed Work](#completed-work)
- [TODO: High Priority](#todo-high-priority)
- [TODO: Medium Priority](#todo-medium-priority)
- [TODO: Low Priority / Research](#todo-low-priority--research)
- [Code Structure Improvements](#code-structure-improvements)
- [Genome-Scale Optimizations](#genome-scale-optimizations)

---

## Overview

This document tracks implementation priorities for the `risearch` Rust port. Items are organized by impact and effort.

**Current Status:** Core functionality complete with C parity (588/588 tests passing). Focus now shifts to performance optimization and code cleanup.

---

## Completed Work

- [x] **Unified `DpView` abstraction** - Merged `rustify` branch, eliminating ~500 lines of duplication in `dp.rs`
- [x] **Pre-allocated trace vectors** - `Vec::with_capacity(best_i + best_j)` optimization
- [x] **Grid resize optimization** - Capacity check before reallocation
- [x] **`BYTE_TO_BASE` lookup table** - O(1) branchless base conversion
- [x] **Pre-computed reverse complement** in index (`sequence_rc` field)
- [x] **C parity testing** - Full test suite with 588 passing tests
- [x] **`SeedSpec` parsing** - Clean implementation with negative index support

---

## TODO: High Priority

### P0-1: Cache Query Suffix Arrays

**Location:** `src/search.rs` - `find_seeds_parallel()`

**Problem:** Query SA rebuilt for every query. For M queries, this is M × O(n log n) wasted work.

**Fix:**

```rust
// Before: inside loop
let query_sa = build_suffix_array(&q_norm);

// After: outside loop or cached
struct QueryCache {
    normalized: Vec<u8>,
    suffix_array: Vec<i64>,
}
```

**Effort:** 2-3 hours  
**Impact:** High for multi-query searches

---

### P0-2: Replace HashMap with Vec in ParallelSearchCache

**Location:** `src/search.rs:752-791`

**Problem:** HashMap indirection + heap allocations per target.

**Fix:**

```rust
struct ParallelSearchCache {
    cached: Vec<Option<CachedTarget>>,  // Index = target_idx
}
```

**Effort:** 1-2 hours  
**Impact:** Medium (reduces GC pressure)

---

### P0-3: Memory-Mapped Index Files

**Location:** `src/sa.rs:213-216`

**Problem:** `std::fs::read()` loads entire index into RAM. Human genome = 29+ GB.

**Fix:**

```rust
use memmap2::Mmap;

pub fn load_index_mmap(path: &Path) -> Result<MappedIndex> {
    let file = File::open(path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    // Deserialize header only, lazy-load sequences
}
```

**Effort:** 4-6 hours  
**Impact:** Critical for genome-scale

---

### P0-4: Compact Suffix Array (u32)

**Location:** `src/sa.rs`

**Problem:** `Vec<i64>` uses 8 bytes per position. Human genome: 25.6 GB just for SA.

**Fix:**

```rust
// For genomes < 4.3B bases
pub forward_sa: Vec<u32>,  // 4 bytes → 50% reduction
```

**Effort:** 2-3 hours  
**Impact:** 50% SA memory reduction

---

## TODO: Medium Priority

### P1-1: Implement LCP Array

**Problem:** No O(1) lookup for "how far can this seed extend?"

**Options:**

- Use `bio::data_structures::suffix_array::lcp()`
- Implement Kasai's algorithm directly

**Effort:** 4-6 hours  
**Impact:** Faster seed extension decisions

---

### P1-2: Add X-drop Heuristic to DP

**Location:** `src/dp.rs`

**Problem:** DP explores unpromising paths too long.

**Fix:**

```rust
const X_DROP: i32 = 100;  // Configurable

if current_score < best_score - X_DROP {
    break;  // Terminate early
}
```

**Effort:** 2-3 hours  
**Impact:** Faster DP for long sequences

---

### P1-3: Stratified Seed Search

**Problem:** All seed lengths searched uniformly.

**Fix:**

1. Start with longest seeds (high specificity)
2. Only fall back to shorter if needed
3. Priority queue by estimated hit quality

**Effort:** 4-8 hours  
**Impact:** Fewer candidates to extend

---

## TODO: Low Priority / Research

### P2-1: SIMD Interval Partitioning

**Location:** `src/parallel_sa.rs:493-500`

**Concept:** Use `std::simd` to accelerate base counting in partition search.

**Effort:** 8-16 hours (research needed)  
**Impact:** Marginal for small alphabets (σ=5)

---

### P2-2: FM-Index Hybrid

**Concept:** For very large genomes, FM-Index is more space-efficient than full SA.

**Trade-offs:**

- Pro: 0.5-2 bytes/base vs 8-12 bytes/base
- Con: Slower position retrieval

**When:** Only if memory-mapped full SA is still too large.

**Effort:** 16-24 hours  
**Impact:** Enables very large genome searches

---

### P2-3: Bit-Packed Sequences

**Concept:** 2-bit encoding (A=00, C=01, G=10, U=11) for 4× memory reduction.

**Challenge:** `Base::Gap` and `Base::N` need 3 bits or sparse bitvector.

**Effort:** 8-12 hours  
**Impact:** 75% sequence memory reduction

---

## Code Structure Improvements

### CS-1: Delete Orphan `lists.rs`

**Location:** `src/lists.rs` (7 lines)

**Problem:** Contains unused `SeedInterval` struct.

**Action:** Delete file, remove from `lib.rs`.

**Effort:** 5 minutes

---

### CS-2: Move Types from `search.rs`

**Problem:** `search.rs` is 1800+ lines mixing concerns.

**Candidates to extract:**

- `Pairing` enum → `alignment.rs` or `types.rs`
- `Alignment` struct → `alignment.rs`
- `SearchHit` → `types.rs` or `hits.rs`
- `Sequence` trait → `seq.rs`

**Effort:** 2-4 hours  
**Impact:** Improved maintainability

---

### CS-3: Export `SeedSpec` from `lib.rs`

**Location:** `src/seed.rs`

**Problem:** `SeedSpec` is well-designed but not exported.

**Action:** Add to `lib.rs` for library consumers.

**Effort:** 5 minutes

---

## Genome-Scale Optimizations

### Summary Table

| Optimization | Memory Impact | Implementation Effort |
|-------------|---------------|----------------------|
| `u32` SA positions | -50% SA | Low |
| Memory-mapped indexes | Enables genome-scale | Medium |
| Bit-packed sequences | -75% seq | Medium |
| FM-Index hybrid | -90% index | High |

### Memory Projection (Human Genome: 3.2B bases)

| Configuration | Total RAM |
|--------------|-----------|
| Current (`i64` SA, `Vec<u8>` seq) | ~29 GB |
| + u32 SA | ~16 GB |
| + Bit-packed seq | ~8 GB |
| + Memory-mapped | ~4 GB active |

---

## Next Actions

1. [ ] Implement P0-1 (Query SA caching)
2. [ ] Implement P0-2 (Vec-based cache)
3. [ ] Implement P0-4 (u32 SA)
4. [ ] Implement P0-3 (mmap indexes)
5. [ ] Delete `lists.rs` (CS-1)
6. [ ] Run benchmarks to measure impact
