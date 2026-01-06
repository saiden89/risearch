
# Flaws and Suboptimal Patterns in Legacy C (RIsearch2)

While the legacy C code contains clever optimizations, it also demonstrates several patterns that are unsafe, unidiomatic, or suboptimal compared to modern Rust practices.

## 1. Quadratic String Construction (`strcat`)

**The Flaw:**
In `print_alignment` (and similar functions), the code builds result strings character-by-character using `strcat` inside a loop.

```c
// search.c:840
while (i > 0 || j > 0) {
    sprintf(buf, "%c", QS(i));
    strcat(temp_query_str, buf); 
    // ...
}
```

**Why it's bad:**
`strcat` must walk the entire length of the destination string to find the null terminator *every single time*.

* **Complexity:** $O(N^2)$ instead of $O(N)$.
* **Impact:** For short strings (30 chars), this is negligible. For longer alignments, this becomes a severe performance bottleneck.
* **Rust Solution:** `String::push()` maintains the length/capacity and appends in $O(1)$.

## 2. Hard exits inside Library Code

**The Flaw:**
The search logic calls `exit()` directly when it encounters mostly recoverable errors or unexpected states.

```c
// search.c:1785
if (!M_left || !Bq_left ...) {
    fprintf(stderr, "Failed allocating DP matrices\n");
    exit(EXIT_FAILURE);
}
// search.c:1344
if (idx == -1) {
    fprintf(stderr, "binary search failed\n");
    exit(1);
}
```

**Why it's bad:**

* **Usability:** This makes the code impossible to use as a library (e.g., from Python or R) because a single error kills the entire host process.
* **Rust Solution:** Return `Result<T, E>` and let the caller decide how to handle panic/failure (propagate up).

## 3. Unsafe Macro Complexity

**The Flaw:**
The codebase relies heavily on complex macros for data structure access, often hiding bitwise operations.

```c
#define Q(ix) (XRIS(qsa[q_start - (ix)]))
#define T(ix) (comp[XRIS(tsa[t_start - (ix)])])
#define TI(q, t) ((q) * (lt) + (t))
```

**Why it's bad:**

* **Debuggability:** You cannot step into a macro or inspect `Q` easily in a debugger.
* **Safety:** No bounds checking. If `q_start - ix` underflows, it accesses random memory (segfault or security vulnerability).
* **Rust Solution:** Helper functions with `#[inline(always)]`. They provide the same speed but with type safety and namespace scoping.

## 4. Manual Memory Management Risks

**The Flaw:**
Extensive use of `calloc` and `free`.

```c
int *M_left = (int *)calloc(dp_size, sizeof(int));
// ... 200 lines of code ...
free(M_left);
```

**Why it's bad:**

* **Leaks:** If the function returns early (e.g., via one of the `if (print_debug) return;` checks seen in debugging variants), the memory is leaked.
* **Correctness:** It requires manual verification that every allocations is paired with a free in all control paths.
* **Rust Solution:** RAII (`Vec`, `Box`). Memory is automatically freed when the variable goes out of scope, regardless of how the function exits.

## 5. Global State and Threading

**The Flaw:**
The code mixes OpenMP `#pragma omp` with global/external variables (`extern int show_alignment;`).
**Why it's bad:**

* While read-only globals are generally thread-safe, this pattern discourages reentrancy and makes it hard to run two concurrent searches with *different* configurations in the same process.

## Summary

The C code prioritizes raw manual control (pointers, macros, custom memory layout) at the cost of safety, readability, and API flexibility. The Rust rewrite correctly addresses the $O(N^2)$ string building and resource safety (RAII), though it currently lacks some of the bit-packing optimizations.
