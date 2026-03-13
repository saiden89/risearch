# Optimizations from Legacy C (RIsearch2)

This note tracks which ideas from the legacy C implementation still look
relevant against the current Rust codebase, and which ones are already matched
or no longer look attractive.

The old version of this file had become stale. In particular, it overstated the
value of SA character packing and described a DP representation that Rust no
longer uses.

## Still Worth Considering

### 1. O(1) Target-ID Recovery for SA Hits

**What C does**

C stores the target/transcript id in each packed SA entry via `XIDX(...)` in
`legacy_c/RIsearch2/src/sa.h`, written during `sa_merge(...)` in
`legacy_c/RIsearch2/src/sa.c`.

That means hit emission can do:

```c
t_start = XSA(sa[j]);
idx = XIDX(sa[j]);
```

with no extra search.

**What Rust does now**

Rust enumerates target suffix positions from the global SA, then remaps each
global position back to a target id with a binary search over target offsets in
`src/seed/search.rs`:

```rust
let target_idx = match offsets.partition_point(|&o| o <= t_pos as u64) {
    0 => continue,
    i => i - 1,
};
```

**Why this still looks attractive**

This is the clearest remaining constant-factor win from the C data model. It is
not about recovering the base. Rust already gets the base in O(1) from
`combined_seq[p]`. The missing piece is target ownership metadata.

**Likely Rust-friendly version**

Do not copy the full packed C layout blindly. The cleaner Rust version is
probably:

- add a sidecar `sa_target_idx: &[u32]`, or
- pack only target id into spare bits if the index format can support it safely.

The sidecar is the safer first step because it avoids overloading one array with
two index spaces.

### 2. Interval-Level Seed Energy Reuse and Early Pruning

**What C does**

In `sa_evaluate_interval(...)` in `legacy_c/RIsearch2/src/search.c`, C computes
the seed energy once for the whole `(query_interval, target_interval, seed_len)`
match interval before expanding the nested `for (j)` / `for (k)` loops.

When seed-threshold filtering is enabled, it also runs the subseed scan once and
can `continue` before materializing any hits from that interval.

**What Rust does now**

Rust expands the interval into individual `SeedHit`s in `src/seed/search.rs`,
and later computes `seed_energy(...)` per hit inside
`compute_seed_extension(...)` in `src/search/mod.rs`.

The relevant call is:

```rust
let seed_e = crate::dsm::seed_energy(model, query_bases, target_trans, q_start, t_start, len);
```

Because all suffixes inside a finished SA match interval share the same seed
string on both sides, this seed energy is interval-invariant. Recomputing it per
hit is unnecessary work.

**Why this still looks attractive**

This is a real C idea that still applies to the Rust pipeline. It can help in
two ways:

- compute the seed energy once and carry it with the interval or `SeedHit`
- apply seed-energy rejection before exploding the query×target cross product

There is also a parity angle here: Rust exposes `filter.seed_energy` in config,
but the current search path only applies `delta_g` filtering in
`src/search/mod.rs`. That suggests interval-level seed filtering may still be
missing, not just unoptimized.

### 3. C-Style Subseed Thresholding

**What C does**

C does not only score the full seed. With `seed_threshold_flag`, it scans all
subseeds down to the minimum seed length and keeps the best energy-per-length
score for the interval.

**Why this still matters**

If Rust wants true RIsearch2 seed-threshold parity, the relevant inspiration is
the whole interval-level subseed heuristic, not just a cheap whole-seed cutoff.

This is more of a feature/parity item than a pure micro-optimization, but it
could also save work by rejecting weak intervals early.

## Already Matched or Better in Rust

### 4. DP Sentinel Storage

This used to be a gap. It no longer is.

Rust does not use `Option<i32>` DP cells in the hot path. The current DP grid
stores plain `i32` values in `DpCell` and reuses the same allocation across
extensions:

- `src/dp/mod.rs`
- `src/dp/core.rs`
- `src/dp/init.rs`

This area is already in better shape than the stale version of this note
claimed.

### 5. Flat DP Arrays

Also already matched.

The legacy C code uses flat 1D indexing macros. Rust already uses flat storage
and contiguous traversal in the current DP grid implementation, so there is no
obvious C-specific gain left here.

### 6. Parallel SA Walking

Also already matched.

The useful part of C's `sa_parallel_match_neg(...)` is the parallel walk over
query and target SA intervals with explicit canonical/wobble branching. Rust's
`SeedSearcher` already does this and adds extra optimizations on top:

- singleton fast paths
- small-interval linear partitioning
- mismatch pruning

This is not an area where Rust is obviously missing a C trick.

## Probably Not Worth Copying Directly

### 7. Packing the Base into the SA Entry

**What C does**

C can recover the base at text position `p` with `XSTRM(sa[p])` because the base
code is packed into each entry.

**Why this is lower priority now**

Rust already gets the base in O(1) from `combined_seq[p]`, where `Base` is
`#[repr(u8)]`. Copying the C layout would:

- complicate the index format
- make the data model harder to reason about
- trade a 1-byte base load for an 8-byte packed-entry load plus bit extraction

That might still be worth measuring one day, but it is not the first place to
look for gains anymore. The target-id lookup path is a much stronger remaining
candidate.

## Current Priority Order

If we want to borrow from C today, the priority should be:

1. Add O(1) target-id recovery for emitted SA hits.
2. Move seed-energy computation and seed-threshold pruning to the interval level.
3. Decide whether full RIsearch2 subseed-threshold parity is required.
4. Leave full SA/base packing alone unless profiles say otherwise.
