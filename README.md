# riordan-poset-enumeration

**Construction, enumeration and analysis of Riordan posets** — partially ordered
sets whose order relation is a binary Riordan matrix (Cheon, Curtis, Kwon and
Mesinga Mwafise, *Linear Algebra Appl.* 632 (2022) 308–331).

[![CI](https://github.com/arnauldmwafise/riordan-poset-enumeration/actions/workflows/ci.yml/badge.svg)](https://github.com/arnauldmwafise/riordan-poset-enumeration/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23188305.svg)](https://doi.org/10.5281/zenodo.23188305)

The repository contains a header-only C++20 library, command-line tools, a large
test suite, independent verification programs, and the data it produced: the
number of Riordan posets for every size up to **n = 28**, each value confirmed by
two separately written programs (previously known: n ≤ 15, OEIS A379608).


## Contents

1. [Background for first-time readers](#1-background-for-first-time-readers)
2. [What m(n) and r(n) count — a worked example](#2-what-mn-and-rn-count--a-worked-example)
3. [Why computing m(n) and r(n) is not trivial](#3-why-computing-mn-and-rn-is-not-trivial)
4. [Results](#4-results)
5. [Quick start](#5-quick-start)
6. [Tools](#6-tools)
7. [Repository layout](#7-repository-layout)
8. [Documentation](#8-documentation)
9. [Verification and reproducibility](#9-verification-and-reproducibility)
10. [Citing](#citing)
11. [License and acknowledgements](#11-license-and-acknowledgements)

---

## 1. Background for first-time readers

### 1.1 Riordan matrices

Take two formal power series

```
g(t) = 1 + g1 t + g2 t^2 + ...        (constant term 1)
f(t) =     t + f2 t^2 + ...            (starts with t)
```

The **Riordan matrix** `(g, f)` is the infinite lower-triangular matrix whose
column `j` (counting from 0) holds the coefficients of `g(t) f(t)^j`:

```
a(i, j) = [t^i] g(t) f(t)^j       (the coefficient of t^i in g f^j)
```

Since `f^j` starts at `t^j`, column `j` starts on the diagonal, so the matrix is
lower-triangular with ones on the diagonal. The classic example is Pascal's
triangle, `g = 1/(1−t)`, `f = t/(1−t)`:

```
(1/(1-t), t/(1-t)) =  1
                      1 1
                      1 2 1
                      1 3 3 1
                      1 4 6 4 1     ...  a(i,j) = C(i, j)
```

Riordan matrices form a group under matrix multiplication and unify many
counting arguments (Shapiro et al., 1991).

### 1.2 Binary Riordan matrices

Reduce everything modulo 2: `g` and `f` get coefficients in {0, 1} and
arithmetic is done in GF(2). The result `{g, f}` is a 0/1 matrix. Pascal's
triangle modulo 2 is the Sierpiński pattern:

```
{1/(1-t), t/(1-t)}  (first 5 rows)   1 0 0 0 0
                                     1 1 0 0 0
                                     1 0 1 0 0
                                     1 1 1 1 0
                                     1 0 0 0 1
```

The `n × n` top-left block is written `{g, f}_n`. It depends only on
`g mod t^n` and `f mod t^n`, i.e. on the bits `g1 … g(n−1)` and `f2 … f(n−1)`,
so there are exactly `2^(2n−3)` of them (for n ≥ 2).

### 1.3 Posets and poset matrices

A **partially ordered set** (poset) is a set with a relation `≤` that is
reflexive (`x ≤ x`), antisymmetric (`x ≤ y ≤ x` ⇒ `x = y`) and transitive
(`x ≤ y ≤ z` ⇒ `x ≤ z`). Label the elements `x_0, …, x_(n−1)` so that a smaller
element always has a smaller label (a *linear extension*). Its **poset matrix** is

```
a(i, j) = 1  if  x_j ≤ x_i,   and 0 otherwise.
```

It is lower-triangular with a unit diagonal (reflexivity and the labelling) and
**transitive**: `a(i,j) = a(j,k) = 1` implies `a(i,k) = 1`. Conversely every
unit lower-triangular transitive 0/1 matrix is the poset matrix of a poset.

### 1.4 The connection: Riordan posets

Every binary Riordan matrix is lower-triangular with ones on the diagonal, so it
is already reflexive and antisymmetric. **When it is also transitive, it is a
poset matrix**, and the poset is called a **Riordan poset**. For example the
matrix above, `{1/(1−t), t/(1−t)}_5`, is transitive; its poset `P_5` (a
"Pascal poset") has the cover relations `0<1, 0<2, 1<3, 2<3, 0<4`:

```
        3
       / \
      1   2   4
       \  |  /
          0
```

Not every binary Riordan matrix is transitive. With `g = 1 + t`, `f = t`:

```
{1+t, t}_3 =  1 0 0
              1 1 0        a(2,1) = 1 and a(1,0) = 1, but a(2,0) = 0:
              0 1 1        x_0 ≤ x_1 ≤ x_2 without x_0 ≤ x_2. Not a poset.
```

So Riordan posets connect two worlds: the algebra of generating functions
(`g`, `f`) and the combinatorics of partial orders. Generating functions give a
very compact description of a whole infinite family of posets, and structural
facts about Riordan matrices turn into facts about posets.

---

## 2. What m(n) and r(n) count — a worked example

* **m(n)** = the number of `n × n` binary Riordan matrices `{g, f}_n` that are
  poset matrices. Each one is a *labelled* poset together with its unique pair
  `(g mod t^n, f mod t^n)`.
* **r(n)** = the number of **non-isomorphic** posets among them, i.e. the number
  of different posets on `n` elements that have at least one Riordan labelling.

  > **Terminology.** *Riordan labelling* is a term coined for this project; it does
> not appear in Cheon et al. (2022). A **labelling** of a poset P on n elements
> numbers its elements `x_0, …, x_(n−1)`; it is a **Riordan labelling** if the
> resulting poset matrix A_P (`a(i,j) = 1` iff `x_j ≤ x_i`) is a binary Riordan
> matrix `{g, f}_n`. Thus P is a Riordan poset exactly when it has at least one
> Riordan labelling, and each Riordan labelling gives one Riordan poset matrix
> (one labeled copy of P). A Riordan poset can have several Riordan labellings,
> each giving a different poset matrix: m(n) counts all of these labeled Riordan
> posets, while r(n) counts only the non-isomorphic (unlabeled) Riordan posets.

A poset can have several labellings that are Riordan matrices, so
`r(n) ≤ m(n)`. Below is the complete list for **n = 4** (generated by
`bin/sample_listing 4`): **m(4) = 22 matrices**, falling into **r(4) = 11
isomorphism classes**. Rows are written `row0 / row1 / row2 / row3`; for
instance `1000 / 0100 / 0010 / 0101` is the matrix

```
1 0 0 0
0 1 0 0
0 0 1 0
0 1 0 1       (x_1 ≤ x_3: the poset "one 2-chain plus two isolated points")
```

| class (poset) | # | matrix rows | g mod t⁴ | f mod t⁴ |
|---|---|---|---|---|
| **1.** antichain `I_4` (no relations) | 1 | `1000 / 0100 / 0010 / 0001` | 1 | t |
| **2.** `C_2 + I_2` (cover 1<3) | 2 | `1000 / 0100 / 0010 / 0101` | 1 | t + t³ |
| | 3 | `1000 / 0100 / 0010 / 1001` | 1 + t³ | t |
| | 4 | `1000 / 0100 / 0110 / 0001` | 1 | t + t² |
| | 5 | `1000 / 0100 / 1010 / 0001` | 1 + t² | t + t³ |
| **3.** two elements below a third, plus a point: `B_{2,1} + 1` (0<3, 1<3) | 6 | `1000 / 0100 / 0010 / 1101` | 1 + t³ | t + t³ |
| | 7 | `1000 / 0100 / 1110 / 0001` | 1 + t² | t + t² + t³ |
| **4.** one element below two others, plus a point: `B_{1,2} + 1` (1<2, 1<3) | 8 | `1000 / 0100 / 0110 / 0101` | 1 | t + t² + t³ |
| | 9 | `1000 / 0100 / 1010 / 1001` | 1 + t² + t³ | t + t³ |
| **5.** two 2-chains, `C_2 + C_2` (1<2, 0<3) | 10 | `1000 / 0100 / 0110 / 1001` | 1 + t³ | t + t² |
| | 11 | `1000 / 0100 / 1010 / 0101` | 1 + t² | t |
| | 12 | `1000 / 1100 / 0010 / 0011` | 1 + t | t + t² + t³ |
| **6.** the "N" poset (1<2, 0<3, 1<3) | 13 | `1000 / 0100 / 0110 / 1101` | 1 + t³ | t + t² + t³ |
| | 14 | `1000 / 0100 / 1010 / 1101` | 1 + t² + t³ | t |
| | 15 | `1000 / 0100 / 1110 / 0101` | 1 + t² | t + t² |
| | 16 | `1000 / 0100 / 1110 / 1001` | 1 + t² + t³ | t + t² + t³ |
| | 17 | `1000 / 1100 / 0010 / 1011` | 1 + t + t³ | t + t² + t³ |
| **7.** complete bipartite `B_{2,2}` (0,1 < 2,3) | 18 | `1000 / 0100 / 1110 / 1101` | 1 + t² + t³ | t + t² |
| **8.** `(C_2 + 1) ⊕ 1` (0<1, 1<3, 2<3) | 19 | `1000 / 1100 / 0010 / 1111` | 1 + t + t³ | t + t² |
| **9.** `1 ⊕ (1 + C_2)` (0<1, 0<2, 2<3) | 20 | `1000 / 1100 / 1010 / 1011` | 1 + t + t² + t³ | t + t² |
| **10.** the diamond / Boolean lattice `B_2` (0<1, 0<2, 1<3, 2<3) | 21 | `1000 / 1100 / 1010 / 1111` | 1 + t + t² + t³ | t + t² + t³ |
| **11.** the chain `C_4` (0<1<2<3) | 22 | `1000 / 1100 / 1110 / 1111` | 1 + t + t² + t³ | t |

(`+` is disjoint union, `⊕` puts the second poset entirely above the first.)
There are 16 posets on 4 elements altogether; the other **5 are not Riordan**
(`B_{3,1}`, `B_{1,3}`, `C_2 ⊕ I_2`, `I_2 ⊕ C_2` and `C_3 + 1`), exactly as the
paper states. Matrix #21 is the Pascal matrix `{1/(1−t), t/(1−t)}_4`
(`1 + t + t² + t³` is `1/(1−t) mod t⁴`), and #22 is `{1/(1−t), t}_4`.

For a smaller warm-up, `bin/sample_listing 3` shows m(3) = 7 and r(3) = 5
(every poset on 3 elements is Riordan; the eighth binary Riordan matrix,
`{1+t, t}_3` from §1.4, is not transitive).

---

## 3. Why computing m(n) and r(n) is not trivial

* **The search space explodes.** There are `2^(2n−3)` binary Riordan matrices of
  order n — about 9 × 10¹⁵ at n = 28 — and each must be tested for
  transitivity (cubic time naively). Only a vanishing fraction are posets
  (5 × 10⁻⁶ at n = 22), so brute force is hopeless beyond n ≈ 12.
* **m(n) needs mathematical structure, not just speed.** The library exploits
  that a new row of `{g,f}` has only two free entries (the rest is forced by
  `g f^j = (g f^{j−1}) f`) and that only the new row can break transitivity, so
  the matrices form a pruned tree with ≈ 2.05ⁿ nodes per level instead of 4ⁿ.
* **r(n) is an isomorphism problem.** Counting posets "up to isomorphism" means
  deciding, for hundreds of millions of matrices, which ones describe the same
  poset under some relabelling. This is graph isomorphism in disguise; it needs
  canonical forms (we use nauty and an independent canonical form of our own)
  and dominates the running time.
* **Scale.** At n = 28 there are 799 million matrices and 573 million classes.
  Storing even a 16-byte fingerprint per matrix needs ~13 GB, so the counter
  partitions the work by an isomorphism invariant, deduplicates in memory pass
  by pass, checkpoints, and can run on many cores.
* **Trust.** A single bug in an enumerator or canonical form silently changes a
  count. Every published value here was therefore obtained by two programs with
  separately written enumerators (and, up to n = 23, two different canonical
  forms), and checked against brute force where brute force is possible.

---

## 4. Results

Counts of Riordan poset matrices m(n) and of Riordan posets up to isomorphism
r(n), for 1 ≤ n ≤ 28. Relative to the published record these counts reach two
milestones:

1. **Unlabelled Riordan posets r(n), extended from n = 15 to n = 28.** Cheon,
   Curtis, Kwon and Mesinga Mwafise (2022) gave r(n) for n ≤ 8, and OEIS
   [A379608](https://oeis.org/A379608) currently lists r(n) for n ≤ 15
   (a(9) added June 2025, a(10)–a(15) January 2026). This software reproduces
   all 15 known terms exactly and adds the **13 new terms r(16), …, r(28)**.
2. **Labeled Riordan posets m(n): a new sequence, n = 1 to 28.** m(n) is the
   number of n × n binary Riordan matrices that are poset matrices. Each such
   matrix is a labeled copy of a Riordan poset, and a single Riordan poset can
   admit several different labeled copies, that is, several distinct Riordan poset
   matrices describing the same poset under different labelings (for example, the
   "N" poset on 4 elements has 5). Hence m(n) ≥ r(n): m(n) counts the labeled
   Riordan posets (all labeled copies), while r(n) counts only the non-isomorphic
   (unlabeled) Riordan posets. This sequence has not been published
   before, and we found no OEIS entry for it (October 2026). **All 28 terms are new.**

Every value was confirmed by two programs built on separately written
enumerators (see the "checked by" column). Both counters are multi-threaded and
resumable, so larger n (29, 30, …) can be reached with the same code on a
multi-core machine (see [§6](#6-tools) and [`verify/README.md`](verify/README.md)).

| n | m(n) | r(n) | checked by |
|---|---|---|---|
| 1 | 1 | 1 | A, B, C |
| 2 | 2 | 2 | A, B, C |
| 3 | 7 | 5 | A, B, C |
| 4 | 22 | 11 | A, B, C |
| 5 | 55 | 33 | A, B, C |
| 6 | 121 | 74 | A, B, C |
| 7 | 214 | 144 | A, B, C |
| 8 | 475 | 232 | A, B, C |
| 9 | 970 | 639 | A, B, C |
| 10 | 2122 | 1406 | A, B, C |
| 11 | 4316 | 3164 | A, B, C |
| 12 | 8876 | 4992 | A, B, C |
| 13 | 16443 | 12501 | A, B, C |
| 14 | 36281 | 26973 | A, B, C |
| 15 | 72909 | 55937 | A, B, C |
| 16 | 162536 | 104169 | A, B, C |
| 17 | 321419 | 266428 | A, B, C |
| 18 | 656225 | 517316 | A, B, C |
| 19 | 1198005 | 941970 | A, B, C |
| 20 | 2634534 | 1742608 | A, B, C |
| 21 | 5243230 | 4424434 | A, C |
| 22 | 11497714 | 9366627 | A, C |
| 23 | 22571793 | 18579721 | A, C |
| 24 | 47110953 | 32869555 | C, D |
| 25 | 86046455 | 73299324 | C, D |
| 26 | 186702111 | 153092515 | C, D |
| 27 | 364554584 | 301274742 | C, D |
| 28 | 799401460 | 572885293 | C, D |

Methods: **A** library enumerator + library canonical form; **B** independent
enumerator (`verify/independent`) + nauty via Python; **C** `bin/rcount`
(fast enumerator + nauty from C++); **D** `verify/independent_rcount.cpp`
(independent enumerator + nauty from C++, different hash; multi-threaded). Every row has two
methods with different enumerators. Machine-readable copies, including OEIS
b-files, are in [`data/`](data/).

Other results of the project (details and proofs in
[`docs/RESULTS.md`](docs/RESULTS.md)):

* **Open question answered.** The paper asks whether every non-Riordan poset can
  be expressed in terms of Riordan posets. It cannot: four 5-element posets are
  neither Riordan nor a direct/ordinal sum of smaller posets. Full
  classification for n ≤ 9.
* **Deciding infiniteness.** For rational `g, f` the question "is `{g,f}_n` a poset
  matrix for every n?" is decidable via 2-automata; implemented in `riordan decide`.
* **A new infinite family**: `{(1+t)^k, t/(1+t)} = P*_{k+1} + P` for all `k ≥ 0`.
* **Errata** to Theorems 5.3 and 5.6 of the paper.
* **Sequence behaviour** ([`docs/SEQUENCE_ANALYSIS.md`](docs/SEQUENCE_ANALYSIS.md)):
  no low-order recurrence or algebraic generating function; m(n)^(1/n) → 2
  (conjecturally) with a slowly growing correction factor.

---

## 5. Quick start

Requirements: a C++20 compiler (GCC ≥ 10 or Clang ≥ 12) and either `make` or
CMake ≥ 3.16. nauty is bundled. Python 3 with `pynauty`, `numpy` and `sympy` is
needed only for the optional verification and analysis scripts.

```sh
git clone https://github.com/arnauldmwafise/riordan-poset-enumeration.git
cd riordan-poset-enumeration

make            # builds bin/riordan, bin/rcount, examples and test suites
make test       # runs all test suites (~950,000 checks, a few seconds)

# or with CMake
cmake -B build && cmake --build build && ctest --test-dir build
```

Five minutes with the tools:

```sh
bin/riordan build "1/(1-t)" "t/(1-t)" 8     # the Pascal poset matrix P_8, poset test, (g,f)
bin/sample_listing 4                        # all m(4) = 22 matrices in r(4) = 11 classes
bin/riordan census 12                       # m(n) and r(n) for n <= 12
bin/rcount 20                               # r(20) = 1742608 in about 10 s
bin/riordan decide "1+t^2" "t/(1+t)"        # proves {1+t^2, t/(1+t)} is an infinite Riordan poset
bin/riordan express 6                       # which posets are Riordan / expressible
```

Using the library (header-only, `#include "riordan/riordan.hpp"`):

```cpp
#include "riordan/riordan.hpp"
#include <iostream>
using namespace riordan;

int main() {
    const std::size_t n = 8;
    Series g = Series::parse("1/(1-t)", n), f = Series::parse("t/(1-t)", n);
    BoolMatrix A = riordan_matrix(g, f, n);          // {g, f}_n
    std::cout << A.str() << "poset matrix: " << A.is_poset_matrix() << "\n";

    // every Riordan poset matrix of order 5, with its generating pair
    for_each_riordan_poset(5, [](const RiordanBuilder& b) {
        if (b.order() == 5) std::cout << b.pair().g.str(false) << " | " << b.pair().f.str(false) << "\n";
        return true;                                  // keep descending
    });
}
```

---

## 6. Tools

`bin/riordan` (general-purpose command-line tool):

| command | purpose |
|---|---|
| `build <g> <f> <n> [--latex\|--json\|--dot]` | the matrix `{g,f}_n`, poset test, output formats |
| `info <g> <f> <n>` | structure: height, width, series-parallel decomposition, forest/binomial tests, Hasse edges, linear extensions, chain counts |
| `horizon <g> <f> [cap]` | largest n with `{g,f}_n` a poset matrix |
| `dual <g> <f> <n>` | dual poset and its pair (Theorem 2.8) |
| `census <N> [--labelled-only] [--list <n>]` | m(n), r(n) for n ≤ N |
| `recognize <file>` | is a given poset Riordan? find all Riordan labellings and their (g, f) |
| `decide <g> <f>` | exact test: infinite Riordan poset? (rational g, f) |
| `toeplitz <n>` / `toeplitz --gens a b ...` | Toeplitz posets, numerical semigroups, Frobenius numbers |
| `exp <family> <param> <n>` | binary exponential Riordan matrices (Stirling, Laguerre, Hermite, Jacobi) |
| `express <N>` | all posets up to N elements: Riordan / expressible / indecomposable |

Series are written with `t`, `+ - * / ^` and parentheses and are reduced mod 2,
e.g. `"1/(1-t)"`, `"t/(1-t^2)"`, `"(1+t)^4"`, `"1+t^3+t^4"`.

`bin/rcount N [options]` (large-scale counter for r(N) and m(N)):

| option | meaning |
|---|---|
| `--threads T` | worker threads (default: all cores) |
| `--split D` | split the search tree into jobs at depth D (default 12; use 16 for > 32 threads) |
| `--passes P` | in-memory mode: P passes over an isomorphism-invariant partition, no disk |
| `--checkpoint FILE` | with `--passes`: record finished passes, resume after interruption |
| `--canon nauty\|library` | canonical form (nauty is ~5× faster; library is independent) |
| `--tmp DIR`, `--buckets B` | disk mode (default without `--passes`) |

Measured on one core: n = 20 in ~10 s, n = 24 in ~4 min, n = 27 in 40 min,
n = 28 in ~2.5 h. Time roughly doubles per step; with T threads divide by
about T (see [`docs/ALGORITHMS.md`](docs/ALGORITHMS.md) for load balancing).

---

## 7. Repository layout

```
include/riordan/      header-only library (one module per file)
  bits.hpp              dynamic bitset with word-level operations
  gf2_series.hpp        truncated power series over GF(2): product, inverse, composition, reversion, parser
  bool_matrix.hpp       0/1 matrices, Boolean products, poset-matrix test
  binary_riordan.hpp    {g,f}_n, the incremental builder (two free bits per row), Lucas test
  poset.hpp             Hasse diagram, Möbius function, height/width, duals, sums, series-parallel
  canon.hpp             canonical form (individualisation-refinement) for posets
  enumerate.hpp         tree enumeration and census of Riordan posets
  recognize.hpp         recognition: Riordan labellings of a given poset
  toeplitz.hpp          numerical semigroups, Frobenius numbers, Toeplitz posets
  exponential.hpp       binary exponential Riordan matrices (exact mod-2 recurrences)
  infinite.hpp          2-automata and the decision procedure for infinite Riordan posets
  decompose.hpp         direct/ordinal decompositions, catalogue of all posets
  io.hpp                LaTeX, Graphviz DOT and JSON output
  riordan.hpp           umbrella header
tools/                riordan_cli.cpp (bin/riordan), rcount.cpp (bin/rcount)
examples/             quickstart.cpp, sample_listing.cpp
tests/                test_common.hpp + one test suite per module (CTest)
verify/               independent re-implementations and cross-check scripts
scripts/              check_rcount.sh, analyze_sequences.py
data/                 counts (CSV) and OEIS-style b-files
docs/                 theory, algorithms, results, sequence analysis
third_party/nauty/    nauty 2.8.8 (Apache-2.0)
```

## 8. Documentation

| document | contents |
|---|---|
| [`docs/THEORY.md`](docs/THEORY.md) | the mathematics behind every routine, with references to the paper |
| [`docs/ALGORITHMS.md`](docs/ALGORITHMS.md) | the enumeration, canonical forms, partitioning, parallelisation and their complexity |
| [`docs/RESULTS.md`](docs/RESULTS.md) | all results: survey of Riordan posets, open question, infinite posets, Toeplitz and exponential posets, errata |
| [`docs/SEQUENCE_ANALYSIS.md`](docs/SEQUENCE_ANALYSIS.md) | recurrence and generating-function searches, growth-rate estimates, conjectures |
| [`verify/README.md`](verify/README.md) | how each published number was cross-checked |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | development workflow, coding and testing conventions |
| [`CHANGELOG.md`](CHANGELOG.md) | release history |

## 9. Verification and reproducibility

* `make test` / `ctest`: ~950,000 checks in 8 suites — every computable statement
  of the paper, fast routines against naive reference implementations, brute
  force over all `2^(2n−3)` matrices for n ≤ 11, and the corrected theorems.
* `scripts/check_rcount.sh`: `bin/rcount` against r(1..16) in all its modes.
* `verify/`: independent implementations (separately written enumerator, Python
  brute force, exact rational arithmetic, nauty) used to confirm every published
  number; see [`verify/README.md`](verify/README.md) for the commands.
* `scripts/analyze_sequences.py`: reproduces the analysis of
  [`docs/SEQUENCE_ANALYSIS.md`](docs/SEQUENCE_ANALYSIS.md).

## Citing

If you use this software or its data, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff) (GitHub shows a "Cite this repository" button),
and cite the paper that introduced Riordan posets:

> G.-S. Cheon, B. Curtis, G. Kwon, A. Mesinga Mwafise, *Riordan posets and
> associated incidence matrices*, Linear Algebra and its Applications 632 (2022)
> 308–331. https://doi.org/10.1016/j.laa.2021.10.002

To obtain a DOI for the software: connect the repository to
[Zenodo](https://zenodo.org/), create a GitHub release, and put the DOI Zenodo
assigns into `CITATION.cff` and the badge above. `.zenodo.json` supplies the
metadata for the Zenodo record.

## 11. License and acknowledgements

The code is released under the [MIT License](LICENSE). The bundled
[nauty](https://pallini.di.uniroma1.it/) by Brendan McKay and Adolfo Piperno is
distributed under the Apache License 2.0 (`third_party/nauty/`). The
mathematical framework is that of Cheon, Curtis, Kwon and Mesinga Mwafise (2022).
