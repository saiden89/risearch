# RIsearch Implementation Strategies: A Deep Dive

This document explains the algorithms and strategies used in `risearch`. It is written to be accessible to readers who may not be experts in bioinformatics algorithms, breaking down the *how* and *why* of each step.

## Table of Contents

- [1. The Core Problem](#1-the-core-problem)
- [2. Seeding Strategy](#2-seeding-strategy)
  - [2.1 The Suffix Array Index](#21-the-suffix-array-index-srcsars)
  - [2.2 Generating Seeds](#22-generating-seeds-srcseedrs)
  - [2.3 Searching the Index](#23-searching-the-index-srcsearchrs)
- [3. Extension Strategy (Dynamic Programming)](#3-extension-strategy-dynamic-programming)
  - [3.1 The Energy Model](#31-the-energy-model-srcdsmrs)
  - [3.2 Dynamic Programming (DP)](#32-dynamic-programming-dp)
  - [3.3 The X-Drop Heuristic](#33-the-x-drop-heuristic)
- [4. Pruning and Optimization](#4-pruning-and-optimization)
  - [4.1 The Maximality Check](#41-the-maximality-check)
  - [4.2 Deduplication](#42-deduplication)
- [5. Future Work & TODOs](#5-future-work--todos)

## 1. The Core Problem

The goal of `risearch` is to find where a short RNA molecule (the **Query**, e.g., a microRNA) binds to a long RNA molecule (the **Target**, e.g., a messenger RNA).

- **Binding Rule**: RNA binds in an **antiparallel** fashion. The 5' end of the query binds to the 3' end of the target.
- **Base Pairing**:
  - **G** pairs with **C** (Strong)
  - **A** pairs with **U** (Strong)
  - **G** can also pair with **U** (Wobble, weaker but valid)

Because full comparisons are slow, we use a **Seed-and-Extend** strategy:

1. **Seed**: Quickly find exact (or near-exact) short matches (e.g., 7 letters long).
2. **Extend**: From those anchor points, carefully expand the alignment to see if it forms a stable, longer interaction.

---

## 2. Seeding Strategy

Seeding is the "fast search" phase. Instead of scanning the entire target sequence byte-by-byte (which is slow), we use a pre-built index called a **Suffix Array**.

### 2.1 The Suffix Array Index (`src/sa.rs`)

Imagine a book's index. Instead of reading the whole book to find "apple", you look it up in the sorted index. A Suffix Array is similar but effectively indexes *every possible substring* of the text.

- **Construction**: We take the target sequence and list every suffix (starting from index 0, index 1, etc.). We then sort these suffixes alphabetically.
- **Why?**: Once sorted, all occurrences of any substring (like `ACGT`) appear next to each other in the list.
- **Strand Logic**: Since RNA is single-stranded but can bind in two ways, we build two indexes for every target:
    1. **Forward Index**: Represents the target sequence as written. Used to find if the query binds to the target's "sense" strand.
    2. **Reverse Index**: Represents the *Reverse Complement* of the target. Used to find if the query binds to the target's "antisense" strand.

### 2.2 Generating Seeds (`src/seed.rs`)

A seed is just a small chunk of the query. If the user asks for a seed length of 7 (`-s 7`), we chop the query into every overlapping 7-mer:

- Query: `ACGTACGT`
- Seeds: `ACGTACG`, `CGTACGT`...

We also support **Intervals** (e.g., "only use seeds from the first 8 bases of the query") because biologically, the beginning of a microRNA (the "seed region") is the most critical part for binding.

### 2.3 Searching the Index (`src/search.rs`)

For each seed, we look it up in the Suffix Array. This gives us a list of **Candidates** (locations in the target) where this short pattern appears.

- **Wobble Support**: If enabled, the search becomes a "branching" search. As we walk down the index, if we see a `G` in the query, we look for both `C` *and* `U` in the target. This finds more matches but is slower.

---

## 3. Extension Strategy (Dynamic Programming)

Once we have a "Candidate" (a seed match), we have a starting point. Now we need to determine the *full* quality of the interaction.

### 3.1 The Energy Model (`src/dsm.rs`)

Biology isn't just about "matching letters"; it's about **Thermodynamics**.

- We don't count "+1 for match, -1 for mismatch".
- We use **Free Energy ($\Delta G$)** values measured in labs (Turner 2004 parameters).
- **Negative Energy is Good**: Lower energy means a more stable bond. A score of `-20.0` is better than `-5.0`.
- The energies are for **stacked pairs**. A `GC` pair followed by an `AU` pair (`GC/AU`) has a specific energy value based on how those atoms stack in 3D space.

### 3.2 Dynamic Programming (DP)

We use an algorithm similar to Smith-Waterman to extend the alignment. Imagine a grid where the Query is on one axis and the Target on the other.
From our seed (the center), we extend in two directions:

1. **Left**: Query moves 5' $\leftarrow$, Target moves 3' $\rightarrow$.
2. **Right**: Query moves 3' $\rightarrow$, Target moves 5' $\leftarrow$.

At each step, we calculate the energy of three possibilities/states:

- **Match (`M`)**: The bases pair (or mismatch). We add the stacking energy.
- **Gap in Target (`Bq`)**: The query has a base that "bulges" out because the target has skipped a spot. This costs a penalty (Gap Open + Gap Extend).
- **Gap in Query (`Bt`)**: The target has a bulge. This also costs a penalty.

### 3.3 The X-Drop Heuristic

To save time, we don't align forever. We keep track of the *best score seen so far*.

- If the current path's score drops too far below the best score (e.g., by more than a set threshold), we assume this alignment has "gone bad" and **stop extending**.
- This "drop-off" ensures we focus only on biologically relevant, high-affinity regions.

---

## 4. Pruning and Optimization

We find *many* seeds. Many are redundant or useless. We need strict rules to filter them.

### 4.1 The Maximality Check

**Concept**: If we find a match for the word **"banana"**, we will also find matches for **"ana"**, **"nan"**, and **"nana"** inside it. We don't want to report all of them separately; we only want the longest one, **"banana"**.

**Implementation (`extend_seed`)**:
Before running the expensive DP extension, we check the letters *immediately outside* the seed.

- If the seed is `ACGT`...
- We look at `Query[left]` and `Target[right]` (and vice versa).
- **Rule**: If we could have extended this seed by one more base *without* using a gap/mismatch (i.e., it forms a perfect pair), then this current seed is **not maximal**. It is just a substring of a longer perfect match that we will find (or have found) elsewhere.
- **Action**: We discard this seed immediately.

### 4.2 Deduplication

Even with the maximality check, we might find overlapping alignments.

- **Example**: A seed at position 1 yields an alignment scoring -20.0. A seed at position 3 yields a shorter part of that same alignment, scoring -15.0.
- The second one is redundant information.

**Strategy**:

1. Collect all successful alignments.
2. **Sort** them by Energy (Best/Lowest first).
3. Walk through the list. Keep a "Verified" list.
4. For each new hit, check if it is **"Shadowed"** by any hit in the Verified list.
    - **Shadowed means**: The new hit is physically *inside* the Verified hit (start/end coordinates are contained) AND it effectively offers no new information (worse score).
    - **Special Case**: If a hit starts at the *exact same position* as a better hit, we might keep it (it could be a shorter alternative nested structure). But if it starts *later* (is a sub-segment), we discard it.

---

## Summary of the Pipeline

1. **Index**: Pre-calculate where every string exists in the Target.
2. **Seed**: Chop Query into small words.
3. **Search**: Find all locations of these words.
4. **Filter (Maximality)**: Ignore words that are just substrings of longer words.
5. **Extend**: Use physics-based energies to grow the alignment from the center out.
6. **Filter (Deduplication)**: Remove weaker alignments that are already covered by stronger ones.
7. **Output**: Print the result.
