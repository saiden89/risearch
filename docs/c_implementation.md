# RIsearch2 C Implementation Technical Reference

This document catalogs the implementation details, quirks, and algorithmic tricks in the C implementation of RIsearch2, structured as a pipeline from input to output.

---

## 0. Data Types and Representations

### Base Encoding

C uses integer indices to represent nucleotides throughout:

| Index | Base | Notes |
|-------|------|-------|
| 0 | Gap | `-` or `.` in alignments |
| 1 | A | Adenine |
| 2 | G | Guanine |
| 3 | C | Cytosine |
| 4 | U/T | Uracil (T treated as U) |
| 5 | N | Ambiguous/unknown |

### XRIS Macro - Byte to Index Conversion

```c
#define XRIS(x) ((x) & 0x07)  // Extract lower 3 bits
```

**The Packed Byte Format**

C stores sequences as packed bytestrings where each byte encodes a base in its lower 3 bits:

```
Byte layout: [unused bits][base index (0-5)]
             ─────────────┬───────────────
             bits 7-3     │  bits 2-0
```

**Example:**

```
ASCII 'a' = 0x61 = 0b01100001
XRIS('a') = 0x61 & 0x07 = 0b00000001 = 1 (A)

ASCII 'g' = 0x67 = 0b01100111  
XRIS('g') = 0x67 & 0x07 = 0b00000111 = 7  ← Wait, that's wrong!
```

**The Trick:** C doesn't store ASCII - it stores **pre-converted** indices:

```c
// During index building, sequences are converted:
text[i] = base_to_index(fasta_char);  // Stores 0-5, not 'a'-'u'
```

So the suffix array text contains bytes like `0x01`, `0x02`, `0x03`, `0x04` (not `'a'`, `'g'`, `'c'`, `'u'`).

**Why This Matters for Rust Parity:**

| Aspect | C | Rust |
|--------|---|------|
| Sequence storage | Pre-converted indices (0-5) | ASCII bytes (`b'a'`, `b'g'`, etc.) |
| DSM lookup | Direct: `DSM[XRIS(seq[i])]` | Needs conversion: `DSM[Base::from_byte(seq[i]).idx()]` |
| Memory per base | Uses only 3 bits conceptually | Uses full byte as ASCII |

Rust's `DsmAccessor::dsm_at()` performs the equivalent conversion at lookup time rather than at index build time.

### Suffix Array Types

```c
typedef int64_t saidx64_t;    // SA positions (64-bit for large genomes)
saidx64_t *sa;                // Target suffix array
saidx64_t *qsa;               // Query suffix array
```

- **sa[]**: Suffix array over concatenated target sequences
- **qsa[]**: Suffix array over query sequence (allows efficient seed lookup)
- **XSA macro**: Extracts position from SA entry

### Text/Sequence Storage

```c
unsigned char *text;          // Concatenated sequences as bytestring
size_t *idx_lengths;          // Length of each sequence
char **name;                  // Sequence identifiers
```

- Sequences stored as **packed bytestrings** (not ASCII)
- Each byte encodes one base via lower 3 bits
- Sequences delimited by special separator in concatenated text

### DSM Table Type

```c
typedef int dsm_t[6][6][6][6];  // 4D array: [Q1][Q2][T1][T2]
const dsm_t dsm_t04_pos;        // Turner 2004 stacking energies
```

**Visual: Double-Strand Stacking**

The DSM table scores the energy of two adjacent base pairs (a "stack"). For RNA-RNA interaction with antiparallel binding:

```
        5' ──────────────────→ 3'   (Query strand, e.g., miRNA)
           Q1       Q2
            │       │
            ▼       ▼
           ═════════════         ← Stacking interaction
            ▲       ▲
            │       │
           T1       T2
        3' ←────────────────── 5'   (Target strand, e.g., mRNA)
```

**DSM Index Mapping:**

```
DSM[Q1][Q2][T1][T2] = stacking energy in centiunits

Where (for a 5'→3' query direction):
  Q1 = query base at position i     (5' side of stack)
  Q2 = query base at position i+1   (3' side of stack)
  T1 = target base pairing with Q1  (3' side on target)
  T2 = target base pairing with Q2  (5' side on target)
```

**Example: AU/UA Stack**

```
Query:   5' - A - U - 3'         Q1=A(1), Q2=U(4)
              |   |
Target:  3' - U - A - 5'         T1=U(4), T2=A(1)

Energy = DSM[1][4][4][1] = stacking energy for AU/UA
```

**Key Points:**

- **Dimensions**: `[query_base_5'][query_base_3'][target_base_5'][target_base_3']`
- **Values**: Centiunits (1/100 kcal/mol), negative = favorable
- **Size**: 6×6×6×6 = 1296 integers
- **Antiparallel**: Target indices run opposite to query direction

### comp[] Complement Table

```c
static const int comp[6] = {0, 4, 3, 2, 1, 5};
// Maps: Gap→Gap, A→U, G→C, C→G, U→A, N→N
```

Used to get Watson-Crick complement index for antiparallel binding semantics.

### Result Structure

```c
typedef struct {
  int score;              // Raw score in centiunits
  float energy;           // Final energy in kcal/mol
  saidx64_t best_left_i;  // Left extension query length
  saidx64_t best_left_j;  // Left extension target length
  saidx64_t best_right_i; // Right extension query length
  saidx64_t best_right_j; // Right extension target length
  // ... alignment strings, coordinates, etc.
} aln_result_t;
```

---

## 1. Index Building

### Input

- FASTA file with RNA/DNA sequences

### Processing

```
For each sequence S:
  1. Normalize: lowercase, T→U, remove gaps
  2. Create RC = reverse_complement(S)
  3. Append to unified text: [S₀][RC₀][S₁][RC₁]...
  4. Build single suffix array over concatenated text
```

> [!WARNING]
> **Key Quirks:**
>
> - **Odd/Even indexing**: Even indices = original sequence, Odd indices = RC
> - **Memory trick**: Single SA covers both strands, halving memory overhead
> - **Sequence delimiter**: Uses special character between sequences

### Output

- Single suffix array `sa[]`
- Sequence length table `idx_lengths[]`
- Base encoding via suffix array positions

---

## 2. Seed Search

### Input

- Query sequence (miRNA)
- Suffix array + target text
- Seed length parameter

### Processing

```
For query seed S:
  1. Convert to search alphabet (lowercase DNA)
  2. Binary search SA for interval matching S
  3. All positions in [start, end) are candidate hits
```

> [!WARNING]
> **Key Quirks:**
>
> - **No explicit strand flag during search** - strand determined by which sequence index (odd/even) the hit lands in
> - **Wobble handling**: Separate SA queries for G-U pairs

### Output

- List of SA positions where seed matches

---

## 3. Position Resolution

### Input

- SA position from seed search

### Processing (search.c:1570-1583)

```c
if (idx % 2 == 1) {           // Hit on RC sequence
  lpos = idx_lengths[idx] - rpos;   // Transform position!
  rpos = idx_lengths[idx] - tpos;
  strand = '+';               // Report as FORWARD strand
} else {                      // Hit on original sequence
  strand = '-';               // Report as REVERSE strand
  lpos += 1;
  rpos += 1;
}
idx = (idx / 2) * 2;          // Map to original sequence index
```

> [!WARNING]
> **Key Quirks:**
>
> - **Strand inversion**: `+` means hit on RC (odd), `-` means hit on original (even)
> - **Position transformation**: For RC hits, positions are mirrored via `length - pos`
> - **1-based adjustment**: Original strand hits get `+1` for 1-based output

### Output

- Strand character (`+` or `-`)
- Transformed target positions (1-based)

---

## 4. Seed Energy Calculation

### Input

- Query positions
- Target positions (from SA)
- Seed length

### Processing (search.c:2118-2134)

```c
#define Q(ix) (XRIS(qsa[q_start+(ix)]))
#define T(ix) (comp[XRIS(sa[t_start+(ix)])])  // <-- comp[] applied!

for (j = 0; j < seed_len - 1; j++) {
  seed_score += (*S)[Q(j)][Q(j+1)][T(j)][T(j+1)];
}
```

> [!WARNING]
> **Key Quirks:**
>
> - **comp[] on target**: Target base indices are complemented via `comp[]` table
> - **Query stays raw**: Query base indices used directly
> - **Stacking model**: Energy for adjacent base pairs (doublets)

### comp[] Table

```
comp[0] = 0  // Gap → Gap
comp[1] = 4  // A → U
comp[2] = 3  // G → C
comp[3] = 2  // C → G
comp[4] = 1  // U → A
comp[5] = 5  // N → N
```

### Output

- Seed energy score (centiunits, negative = favorable)

---

## 5. DP Extension

### Input

- Seed boundaries
- Maximum extension length
- DSM scoring matrix

### Processing - DP_left (5' extension)

```c
#define Q(ix) (XRIS(qsa[q_start-(ix)]))        // Query decreasing
#define T(ix) (comp[XRIS(tsa[t_start-(ix)])])  // Target decreasing + comp

M[i,j] = max(M[i-1,j-1], Bq[i-1,j-1], Bt[i-1,j-1]) + S[Q(i-1),Q(i),T(j-1),T(j)]
Bq[i,j] = max(Bq[i-1,j] + gap_extend, M[i-1,j] + gap_open)
Bt[i,j] = max(Bt[i,j-1] + gap_extend, M[i,j-1] + gap_open)
```

### Processing - DP_right (3' extension)

```c
#define Q(ix) (XRIS(qsa[q_start+(ix)]))        // Query increasing
#define T(ix) (comp[XRIS(tsa[t_start+(ix)])])  // Target decreasing + comp

// Same DP recurrence, different index directions
```

> [!WARNING]
> **Key Quirks:**
>
> - **comp[] applied in BOTH DP functions**
> - **Index directions differ**: DP_left decreases, DP_right increases
> - **Antiparallel binding**: Query 5'→3' pairs with Target 3'→5'
> - **Tie-breaking order**: When scores are equal, preference is `M > Bq > Bt` (match before gaps, query gap before target gap). C's `max(a, b)` returns first argument when tied.

### Output

- Best extension score
- Extension lengths (query and target)
- Traceback for alignment reconstruction

---

## 6. Energy Finalization

### Input

- Seed energy
- Left extension score
- Right extension score
- Nucleotide count

### Processing (search.c:1868-1870)

```c
result->score = total_score + nt_count * extPen;
result->energy = (result->score - 559.0f) / -100.0f;
```

> [!WARNING]
> **Key Quirks:**
>
> - **559 correction term**: Initiation free energy adjustment
> - **Division by -100**: Convert centiunits to kcal/mol with sign flip
> - **Extension penalty**: Per-nucleotide penalty (`extPen`, typically 3.0)

### Output

- Final energy in kcal/mol (negative = favorable)

---

## 7. Alignment String Construction

### Input

- Traceback from DP
- Query and target sequences

### Processing

- Reconstruct alignment by following traceback
- Build interaction string with pairing symbols (P, W, U, Q)

> [!WARNING]
> **Key Quirks:**
>
> - **Multiple print formats**: `-p1`, `-p2`, `-p3` with different column layouts
> - **Seed markers**: `[` and `]` delimit seed in interaction string

---

## Summary of C Quirks for Rust Parity

| Aspect | C Behavior | Parity Impact |
|--------|-----------|---------------|
| SA structure | Interleaved original+RC | Rust uses separate forward/reverse SA |
| Strand convention | `+` = RC hit, `-` = original | Rust has explicit Strand::Forward/Reverse |
| Position transform | `length - pos` for RC hits | Rust uses direct positions on RC'd sequence |
| comp[] in DSM | Always applied to target | **Unknown if Rust needs this** |
| Tie-breaking | First argument wins | Rust must match with `>=` |
| 1-based output | Applied differently by strand | Rust applies uniformly |

---

## Open Questions

1. **Why does applying comp[] in Rust break forward strand hits?**
   - C applies comp[] everywhere, but Rust's sequence representation may already encode this differently

2. **Is the 31 centiunit discrepancy in seed energy or DP extension?**
   - Needs targeted trace logging to isolate

3. **Does Rust's on-the-fly RC interact poorly with DSM indexing?**
   - C never generates RC on-the-fly; it's pre-computed in the SA
