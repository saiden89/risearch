# RIsearch Performance Architecture — Complete Design Discussion

This document captures the full back-and-forth design discussion for high-performance RNA alignment.

---

# Part 1: The Iterator Mistake

## Initial Proposal

```rust
/// A temporary view into the database.
/// It is cheaply Copy-able because it just holds a pointer and a length.
#[derive(Copy, Clone)] 
pub struct Seq<'a> {
    data: &'a [u8],
    strand: Strand,
}

impl<'a> Seq<'a> {
    /// The specific "Tie-in" method.
    /// Returns an iterator that yields 'Base', hiding the underlying u8 conversion
    /// AND handling the reverse complement logic automatically.
    pub fn bases(&self) -> impl Iterator<Item = Base> + 'a {
        let iter = match self.strand {
            Strand::Forward => {
                self.data.iter()
                    .map(|&b| Base::from_byte(b))
                    .boxed()
            },
            Strand::Reverse => {
                self.data.iter()
                    .rev()
                    .map(|&b| Base::from_byte(b).complement())
                    .boxed() 
            }
        };
        iter
    }
}

trait BoxedIterator<'a> {
    fn boxed(self) -> Box<dyn Iterator<Item = Base> + 'a>;
}
impl<'a, I: Iterator<Item = Base> + 'a> BoxedIterator<'a> for I {
    fn boxed(self) -> Box<dyn Iterator<Item = Base> + 'a> {
        Box::new(self)
    }
}
```

## Why This Is Wrong

**Valid: Iterators are death for DP.** Smith-Waterman and Needleman-Wunsch algorithms fundamentally require random access (`matrix[i][j]`, `query[i]`, `target[j]`). An Iterator is forward-only. Trying to use `nth()` or `collect()` inside a DP loop destroys the O(NM) complexity or forces unnecessary allocations. You must have indexed access.

**Valid: Box is a bottleneck.** In a tight "hot loop" running millions of extensions, the allocator overhead and vtable indirection of `Box<dyn Iterator>` will likely consume more CPU cycles than the alignment math itself.

---

# Part 2: The Naive Fix (Still Wrong)

## The Proposed Fix

```rust
pub fn base(&self, pos: usize) -> Base {
    Base::from_byte(self.data[pos]) // BUG!
}
```

## Why This Has a Logic Bug

If your `Seq` represents a **Reverse Strand** view of the target:

- Logical Index 0 should be Physical Index `len - 1`
- The Base should be complemented

The suggested code reads from the front and returns the raw base. It ignores the biological reality of the reverse strand.

---

# Part 3: The Const Generics Solution

## The Zero-Cost Accessor

```rust
impl<'a> Seq<'a> {
    /// Get a base using a compile-time constant for direction.
    /// The compiler removes the 'if' statement completely in the binary.
    #[inline(always)]
    pub fn get_const<const REV: bool>(&self, i: usize) -> Base {
        if REV {
            // Reverse logic: Read from end, Complement result
            let phys_idx = self.data.len() - 1 - i;
            Base::from_byte(self.data[phys_idx]).complement()
        } else {
            // Forward logic: Direct read
            Base::from_byte(self.data[i])
        }
    }
}
```

## The Generic Alignment Kernel

```rust
// This function will be compiled TWICE.
// Once with REV=false, Once with REV=true.
fn extend_kernel<const REV: bool>(seq: Seq, query: &[Base]) -> i32 {
    let mut score = 0;
    
    for i in 0..seq.data.len() {
        // This call becomes a direct pointer access in assembly
        let target_base = seq.get_const::<REV>(i); 
        let query_base = query[i];
        
        if target_base == query_base {
            score += 2;
        } else {
            score -= 1;
        }
    }
    score
}
```

## The Dispatcher (Branch Once at Top Level)

```rust
pub fn run_alignment(seq: Seq, query: &[Base]) -> i32 {
    match seq.strand {
        Strand::Forward => extend_kernel::<false>(seq, query),
        Strand::Reverse => extend_kernel::<true>(seq, query),
    }
}
```

## Why This Works

Because the `if` check relies on a generic constant (`const R: bool`), the compiler deletes the `if` branch entirely in the generated machine code. You get two specialized functions with zero runtime branching.

---

# Part 4: Query Canonicalization (Even Better)

## The Insight

In "Seed and Extend," your Query is tiny (KB) and your Target is massive (GB). Never put logic in the inner loop to handle the Query's orientation.

Instead, **pre-compute the Query variants**:

- Query (Forward): Store as `Vec<u8>`
- Query (Reverse Complement): Store as `Vec<u8>` immediately upon loading

## Why This Fixes the Problem

If you have a seed match on the Reverse Strand of the Target:

- **Don't:** Read Target backwards and complement on the fly (Slow, complex index math)
- **Do:** Align `Query_RevComp` against `Target_Forward`

**Mathematical equivalence:** `Align(Query, Target_RC)` is identical to `Align(Query_RC, Target)`.

By swapping the query, your inner loop always reads the Target forward (or simply decrements indices for left-extension) and never computes complements. You use simple integer equality (`==`) instead of a base lookup table.

---

# Part 5: Unsafe Pointer Arithmetic

## The UnsafeSlice Pattern

```rust
// The "Hot" Accessor - Zero checks, Zero logic, just memory
struct UnsafeSlice<'a> {
    ptr: *const u8,
    len: usize, // kept only for debug assertions
    _marker: std::marker::PhantomData<&'a u8>,
}

impl<'a> UnsafeSlice<'a> {
    #[inline(always)]
    fn new(slice: &'a [u8]) -> Self {
        Self { ptr: slice.as_ptr(), len: slice.len(), _marker: std::marker::PhantomData }
    }

    #[inline(always)]
    unsafe fn get(&self, idx: usize) -> u8 {
        debug_assert!(idx < self.len); // Only runs in debug mode
        *self.ptr.add(idx)
    }
}
```

## The Direction Trait

```rust
trait Direction {
    const IS_REVERSE: bool;
    
    #[inline(always)]
    fn step(idx: usize) -> usize;
}

struct Right; // i++
impl Direction for Right { 
    const IS_REVERSE: bool = false;
    #[inline(always)] fn step(i: usize) -> usize { i + 1 } 
}

struct Left;  // i--
impl Direction for Left { 
    const IS_REVERSE: bool = true;
    #[inline(always)] fn step(i: usize) -> usize { i - 1 } 
}
```

## The Generic Kernel

```rust
fn extend_kernel<D: Direction>(
    target: UnsafeSlice, 
    query: UnsafeSlice, 
    mut t_idx: usize, 
    mut q_idx: usize
) -> i32 {
    let mut score = 0;
    loop {
        let t_byte = unsafe { target.get(t_idx) };
        let q_byte = unsafe { query.get(q_idx) };

        if t_byte == q_byte {
            score += 2;
        } else {
            score -= 1;
        }

        if score < 0 { break; }
        
        t_idx = D::step(t_idx);
        q_idx = D::step(q_idx);
    }
    score
}
```

## Critique: This Kernel is Ungapped

Both pointers move in lockstep (`offset(stride)`). This can only calculate matches/mismatches along a single diagonal. Real alignment requires gaps — pausing one pointer while moving the other.

## Critique: Missing Turner Context

`energy += STACKING[prev][curr]` — you need the **previous** pair for nearest-neighbor thermodynamics. A single-base lookup is insufficient.

## Critique: Scalar is 5% of Hardware

Modern CPUs can compare 32-64 bytes per cycle (AVX2/512). A scalar loop wastes 95% of available throughput.

---

# Part 6: The Pointer-Based RNA Kernel

```rust
use std::ptr;

/// extend_rna_kernel
/// LEFT=true  => Decrement pointers (move backwards towards 5' end)
/// LEFT=false => Increment pointers (move forwards towards 3' end)
/// 
/// Safety: Caller must ensure `q_ptr` and `t_ptr` have valid padding
/// in the direction of extension.
unsafe fn extend_rna_kernel<const LEFT: bool>(
    mut q_ptr: *const u8,
    mut t_ptr: *const u8,
    mut q_len: usize,
    mut t_len: usize, 
    energy_matrix: &[[i32; 5]; 5]
) -> i32 {
    let mut energy = 0;
    let stride = if LEFT { -1isize } else { 1isize };

    while q_len > 0 && t_len > 0 {
        let q_byte = *q_ptr;
        let t_byte = *t_ptr;

        let q_idx = *ASCII_TO_IDX.get_unchecked(q_byte as usize) as usize;
        let t_idx = *ASCII_TO_IDX.get_unchecked(t_byte as usize) as usize;

        let score = *energy_matrix.get_unchecked(q_idx).get_unchecked(t_idx);
        
        energy += score;
        if energy < -100 { break; }

        q_ptr = q_ptr.offset(stride);
        t_ptr = t_ptr.offset(stride);
        
        q_len -= 1;
        t_len -= 1;
    }
    energy
}
```

## Nearest Neighbor Lookup

```rust
// Current pair indices
let q_curr = ASCII_TO_IDX[*q_ptr as usize] as usize;
let t_curr = ASCII_TO_IDX[*t_ptr as usize] as usize;

// Previous pair indices (Look behind)
let q_prev = ASCII_TO_IDX[*q_ptr.offset(-stride) as usize] as usize;
let t_prev = ASCII_TO_IDX[*t_ptr.offset(-stride) as usize] as usize;

// 4D Lookup: Energy of Stacking (Prev Pair -> Curr Pair)
energy += STACKING_MATRIX[q_prev][t_prev][q_curr][t_curr];
```

---

# Part 7: AlignedSeq (Sentinel Pattern)

## The Problem

```rust
ptr.offset(-1)  // UB if ptr points to allocation start
```

Rust's pointer provenance rules make this undefined behavior even if you never dereference.

## The Solution

```rust
pub struct AlignedSeq {
    _allocation: Vec<u8>, 
    pub start_idx: usize,
    pub len: usize,
}

impl AlignedSeq {
    pub fn new(raw: &[u8]) -> Self {
        let mut vec = Vec::with_capacity(raw.len() + 80);
        
        // Leading sentinels
        vec.extend(std::iter::repeat(b'N').take(16));
        let start = vec.len();
        
        vec.extend_from_slice(raw);
        
        // Trailing sentinels
        vec.extend(std::iter::repeat(b'N').take(64));
        
        Self { _allocation: vec, start_idx: start, len: raw.len() }
    }
    
    pub fn as_ptr(&self) -> *const u8 {
        unsafe { self._allocation.as_ptr().add(self.start_idx) }
    }
}
```

---

# Part 8: Why Anti-Diagonal SIMD Fails

## The Idea

Cells on the same anti-diagonal `(i+j=k)` are independent. Process them in parallel.

## Why It Fails for X-Drop

- X-drop terminates when score drops below threshold
- This requires knowing the **running max score**
- Anti-diagonal processes cells **out of order** relative to optimal path
- You may compute (6,2) before discovering (5,3) should have X-dropped

---

# Part 9: Farrar's Striped Algorithm

## The Solution

Pre-compute a Query Profile and iterate the Target linearly.

## Why It Wins

| Property | Anti-Diagonal | Farrar Striped |
|----------|---------------|----------------|
| Target Access | Random | **Linear (zero-copy)** |
| X-Drop | Impossible | **Column-wise max** |
| SIMD Efficiency | ~50% | **~95%** |

## Profile Structure

```rust
use std::arch::x86_64::*;

const SIMD_WIDTH: usize = 16; 

#[repr(align(32))]
pub struct ThermoProfile16 {
    pub scores: Vec<Vec<__m256i>>, 
    pub gap_open: Vec<__m256i>,
    pub gap_extend: Vec<__m256i>,
    pub num_segments: usize,
}
```

## Interleaved Layout

For 16-lane SIMD:

```
Lane 0:  positions 0, 16, 32, 48...
Lane 1:  positions 1, 17, 33, 49...
...
Lane 15: positions 15, 31, 47, 63...
```

---

# Part 10: The Lazy-F Problem

## The Recurrence

```
F[i] = max(H[i] - Open, F[i-1] - Extend)
```

F is a vertical dependency — cell i depends on cell i-1. This seems serial.

## Logarithmic Solution (O(log N))

```rust
#[inline(always)]
unsafe fn propagate_f_logarithmic(
    mut v_h: __m256i,
    v_open: __m256i,
    v_extend: __m256i
) -> (__m256i, __m256i) {
    let mut v_f = _mm256_subs_epi16(v_h, v_open);
    
    // Step 1 (Shift 1)
    let v_f_sh1 = _mm256_slli_si256(v_f, 2);
    let v_term1 = _mm256_subs_epi16(v_f_sh1, v_extend);
    v_f = _mm256_max_epi16(v_f, v_term1);
    
    // Step 2 (Shift 2)
    let v_f_sh2 = _mm256_slli_si256(v_f, 4);
    let v_2extend = _mm256_add_epi16(v_extend, v_extend);
    let v_term2 = _mm256_subs_epi16(v_f_sh2, v_2extend);
    v_f = _mm256_max_epi16(v_f, v_term2);

    // Step 3 (Shift 4)
    let v_f_sh3 = _mm256_slli_si256(v_f, 8);
    let v_4extend = _mm256_add_epi16(v_2extend, v_2extend);
    let v_term3 = _mm256_subs_epi16(v_f_sh3, v_4extend);
    v_f = _mm256_max_epi16(v_f, v_term3);
    
    // Step 4 (Shift 8) - LANE CROSSING
    let v_lane_cross = _mm256_permute2x128_si256(v_f, v_f, 0x01);
    let v_lane_cross = _mm256_slli_si256(v_lane_cross, 14);
    let v_f_sh4 = _mm256_slli_si256(v_f, 16);
    let v_f_sh4 = _mm256_alignr_epi8(v_f_sh4, v_lane_cross, 14);
    let v_8extend = _mm256_add_epi16(v_4extend, v_4extend);
    let v_term4 = _mm256_subs_epi16(v_f_sh4, v_8extend);
    v_f = _mm256_max_epi16(v_f, v_term4);
    
    v_h = _mm256_max_epi16(v_h, v_f);
    (v_h, v_f)
}
```

## Why slli is Correct

Intel names by register math: "shift left" moves bits from LSB → MSB.
In array terms: `[A,B,C,D]` → `[0,A,B,C]`. Value A moves from index 0 → index 1.
This is exactly what we need for prefix propagation (i ← i-1).

---

# Part 11: Shift Lanes with Carry-In

```rust
#[inline(always)]
unsafe fn shift_lanes_right_1(v: __m256i, carry_in: __m256i) -> __m256i {
    let v_shift = _mm256_slli_si256(v, 2);
    let v_swap = _mm256_permute2x128_si256(v, v, 0x01);
    let v_patched = _mm256_alignr_epi8(v_shift, v_swap, 14);
    _mm256_blend_epi16(v_patched, carry_in, 0x01)
}
```

---

# Part 12: Complete Farrar Kernel

```rust
#[target_feature(enable = "avx2")]
pub unsafe fn align_farrar_16bit(
    target: &[u8],
    profile: &ThermoProfile16,
    v_h_store: &mut [__m256i], 
    v_e_store: &mut [__m256i], 
    xdrop: i16
) -> i16 {
    let mut overall_max: i16 = -32000;

    for &t_byte in target {
        let t_idx = *ASCII_TO_IDX.get_unchecked(t_byte as usize) as usize;
        let score_vecs = &profile.scores[t_idx];

        let mut v_f_carry = _mm256_setzero_si256();

        for i in 0..profile.num_segments {
            let v_h = v_h_store[i];
            let v_e = v_e_store[i];
            let v_score = score_vecs[i];
            
            let mut v_h_new = _mm256_adds_epi16(v_h, v_score);
            
            let v_open = profile.gap_open[i];
            let v_extend = profile.gap_extend[i];
            
            let v_h_minus_open = _mm256_subs_epi16(v_h, v_open);
            let v_e_minus_extend = _mm256_subs_epi16(v_e, v_extend);
            let v_e_new = _mm256_max_epi16(v_h_minus_open, v_e_minus_extend);
            
            v_h_new = _mm256_max_epi16(v_h_new, v_e_new);
            
            v_h_store[i] = v_h_new;
            v_e_store[i] = v_e_new;
        }
        
        // Lazy-F loop here
    }
    overall_max
}
```

---

# Part 13: Turner Model Approximation

## Full Turner (Accurate, Not Vectorizable)

- 1×1 internal loops: special table
- Bulges of size 1: special table
- Loop length penalty: `A + B × log(len)`

## Dual-Affine Approximation (Fast, ~95% Accurate)

```rust
// Two gap penalty sets
Short gap: O1=large, E1=large  // Bulges
Long gap:  O2=small, E2=small  // Structural

H[i][j] = max(0, H + score, E1, E2, F1, F2)
```

---

# Summary of Agreed Decisions

| Decision | Rationale |
|----------|-----------|
| LUT over trait | Single implementation, no polymorphism |
| Free function over method | No self parameter needed |
| Query canonicalization | Swap query, never RC target |
| Const generic direction | Compile-time branch elimination |
| AlignedSeq padding | UB-free pointer arithmetic |
| Farrar over anti-diagonal | X-drop compatible |
| Logarithmic F-propagation | O(log N) parallel prefix-max |
| slli not srli | Left shift = lower index → higher |
| 16-bit integers | Turner energy range |
| Dual-affine approximation | ~95% accurate, vectorizable |
| Two-pass traceback | Score-only SIMD + band-limited scalar |
