# Type Safety Refactor: Sequence Abstraction

## Summary

Replace raw `Vec<u8>` / `&[u8]` sequence handling with a `Sequence` newtype wrapping `Vec<Base>`. Normalize once at input boundary, work with `Base` throughout, convert to bytes only for SA construction (once at index build time).

## Core Design

### Base (existing, minor changes)

```rust
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub enum Base {
    #[default]
    Gap = 0,
    A = 1,
    G = 2,
    C = 3,
    U = 4,
    N = 5,
}

impl Base {
    // Existing
    pub fn from_byte(b: u8) -> Self;     // LUT: ASCII → Base
    pub fn idx(self) -> usize;           // Free: self as usize
    pub fn complement(self) -> Self;     // LUT: Base → Base

    // Add
    pub fn to_byte(self) -> u8;          // LUT: Base → ASCII (for display)
}
```

**Key insight**: `idx()` is free because discriminant values (0-5) ARE the indices. No LUT needed.

### Sequence (new)

```rust
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequence(Vec<Base>);

impl Sequence {
    /// Normalize raw bytes into a Sequence. Only place ASCII→Base conversion happens.
    pub fn normalize(id: &str, raw: &[u8]) -> Result<(Self, NormalizationStats)>;

    /// Length
    pub fn len(&self) -> usize { self.0.len() }
    pub fn is_empty(&self) -> bool { self.0.is_empty() }

    /// Reverse complement
    pub fn reverse_complement(&self) -> Sequence;

    /// For SA construction only (O(n), done once at index build)
    pub fn to_bytes(&self) -> Vec<u8> {
        self.0.iter().map(|b| *b as u8).collect()
    }
}

impl Deref for Sequence {
    type Target = [Base];
    fn deref(&self) -> &[Base] { &self.0 }
}

impl Index<usize> for Sequence {
    type Output = Base;
    fn index(&self, i: usize) -> &Base { &self.0[i] }
}
```

## Data Flow

```
FASTA input (&[u8])
    ↓
Sequence::normalize()  ← Only ASCII→Base conversion
    ↓
Sequence (Vec<Base>)
    ↓
├── SA construction: sequence.to_bytes() → libsais (once at index build)
├── SA search: compare Base directly (no conversion)
├── DP: pass &[Base], use base.idx() for DSM
└── Output: base.to_byte() for display
```

## SA Ordering Change

**Old**: SA sorted by ASCII bytes (a=97 < c=99 < g=103 < n=110 < t=116)

**New**: SA sorted by discriminant values (A=1 < G=2 < C=3 < U=4 < N=5)

This is a **breaking change** to index format. Old indices must be rebuilt.

### Partition logic change

```rust
// Old (ASCII boundaries)
partition_point(|&idx| seq[idx + offset] < b'c')
partition_point(|&idx| seq[idx + offset] < b'g')
partition_point(|&idx| seq[idx + offset] < b'n')
partition_point(|&idx| seq[idx + offset] < b't')

// New (Base boundaries)
partition_point(|&idx| seq[idx + offset] < Base::C)
partition_point(|&idx| seq[idx + offset] < Base::U)
partition_point(|&idx| seq[idx + offset] < Base::N)
// Note: order changes because A(1) < G(2) < C(3) < U(4) < N(5)
```

## Storage Changes

### SequenceIndex

```rust
// Old
pub struct SequenceIndex {
    pub name: String,
    pub sequence: Vec<u8>,
    pub sequence_rc: Vec<u8>,
    pub forward_sa: Vec<u32>,
    pub reverse_sa: Vec<u32>,
}

// New
pub struct SequenceIndex {
    pub name: String,
    pub sequence: Sequence,
    pub sequence_rc: Sequence,
    pub forward_sa: Vec<u32>,
    pub reverse_sa: Vec<u32>,
}
```

## DP Interface Change

```rust
// Old
pub fn dp_left(query: &[u8], target: &[u8], ...) -> ...
pub fn dp_right(query: &[u8], target: &[u8], ...) -> ...

// New
pub fn dp_left(query: &[Base], target: &[Base], ...) -> ...
pub fn dp_right(query: &[Base], target: &[Base], ...) -> ...

// Inside DP, change:
// Old: Base::from_byte(seq[i]).idx()
// New: seq[i].idx()
```

## Files to Modify

### Phase 1: Core types
| File | Change |
|------|--------|
| `src/types.rs` | Add `Ord` derive to Base, add `to_byte()` method |
| `src/seq/sequence.rs` | NEW - Create Sequence type |
| `src/seq/mod.rs` | Export Sequence, delete `Seq<'a>` |

### Phase 2: Normalization & operations
| File | Change |
|------|--------|
| `src/seq/normalize.rs` | Return `Sequence` instead of `Vec<u8>` |
| `src/seq/rc.rs` | Operate on `&Sequence` → `Sequence` |

### Phase 3: Storage
| File | Change |
|------|--------|
| `src/index/sa.rs` | `SequenceIndex` uses `Sequence` |

### Phase 4: SA construction & search
| File | Change |
|------|--------|
| `src/seed/sa.rs` | Build SA from `Sequence::to_bytes()` |
| `src/seed/searcher.rs` | Take `&[Base]`, update partition boundaries |
| `src/seed/search.rs` | Use `Sequence` |

### Phase 5: DP interface
| File | Change |
|------|--------|
| `src/dp/core.rs` | Change interface to `&[Base]` |
| `src/dp/init.rs` | Change interface to `&[Base]` |
| `src/dp/traceback.rs` | Change interface to `&[Base]` |
| `src/dp/mod.rs` | Update exports if needed |

### Phase 6: Search & output
| File | Change |
|------|--------|
| `src/search/core.rs` | Use `Sequence`, pass `&[Base]` to DP |
| `src/search/stream.rs` | Use `Sequence` |
| `src/search/output.rs` | Use `Base::to_byte()` for display |

### Phase 7: App entry points
| File | Change |
|------|--------|
| `src/app/mod.rs` | Use `Sequence` |

## Execution Plan

Each phase can be a separate subagent task:

1. **Phase 1**: Core types (Base changes, create Sequence)
2. **Phase 2**: Normalization & reverse complement
3. **Phase 3**: SequenceIndex storage
4. **Phase 4**: SA construction & searcher
5. **Phase 5**: DP interface
6. **Phase 6**: Search & output
7. **Phase 7**: App integration, fix compilation errors
8. **Phase 8**: Run tests, fix any issues

## Performance Summary

| Operation | Cost |
|-----------|------|
| `base.idx()` | Free (discriminant cast) |
| `base.complement()` | 1 LUT lookup |
| `base.to_byte()` | 1 LUT lookup |
| `Base::from_byte(b)` | 1 LUT lookup |
| `sequence.to_bytes()` | O(n), once at index build |
| SA partition compare | Free (Base comparison) |
| DSM access | Free (idx is free) |

## Breaking Changes

1. **Index format**: Old indices incompatible, must rebuild
2. **SA ordering**: Results may come in different order (semantically equivalent)
3. **DP interface**: `&[u8]` → `&[Base]` (internal change, not user-facing)
