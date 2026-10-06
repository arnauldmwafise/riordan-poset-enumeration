# Algorithms

How m(n) and r(n) are computed, which mathematical shortcuts make it possible,
how the work is partitioned and parallelised, and what it costs. Notation:
`a(i,j) = [t^i] g f^j mod 2`, `{g,f}_n` the leading n × n block, rows and
columns indexed from 0, `a(i,j) = 1` iff `x_j ≤ x_i`.

## 1. Not brute force: the mathematical shortcuts

Brute force would build all `2^(2n−3)` binary Riordan matrices of order n and
test each for transitivity — 2.2 × 10¹² matrices at n = 22, 9 × 10¹⁵ at n = 28.
The implementation instead relies on four facts.

**F1 (bijection).** `{g,f}_n` determines `g mod t^n` (column 0) and
`f mod t^n` (column 1 is `g f`, and g is invertible). So matrices and truncated
pairs correspond one-to-one, and the generating functions of any matrix found are
read off from its first two columns (`recover_pair`).

**F2 (two free entries per row).** Since `g f^j = (g f^{j−1}) f`,

```
a(n, j) = Σ_k a(k, j−1) f_{n−k}   (mod 2),   j ≥ 2,
```

and the right-hand side involves only rows < n and `f_1 … f_{n−1}`. Hence when a
row is added, entries 2 … n−1 are *forced*, the diagonal is 1, and only
`a(n,0) = g_n` and `a(n,1) = f_n + Σ_{k=1}^{n−1} g_k f_{n−k}` are free. Each
extension has at most four children, in bijection with `(g_n, f_n)`.

**F3 (local transitivity).** In a lower-triangular matrix a violated chain
`i ≥ j ≥ k` must have the new element on top, so only the new row has to be
checked: `row_j ⊆ row_n` for every j with `a(n,j) = 1`.

**F4 (heredity).** Leading principal submatrices of poset matrices are poset
matrices (Theorem 3.6 of the paper), so a branch that fails can be discarded
without losing anything below it.

Together: the Riordan poset matrices form a rooted tree (root `[1]`, children =
valid extensions); level n holds exactly the m(n) matrices; every matrix is
reached exactly once (by F1 its parent and its choice bits are unique).

A further cheap pruning step: if the *forced* entries (columns ≥ 2) of the new row
are already not closed, no choice of the two free entries can repair it (they
only pull in rows 0 and 1, whose entries lie in columns 0 and 1), so the node has
no children.

Empirically the tree has about 2.05ⁿ nodes on level n instead of 4ⁿ/8: in most
rows transitivity forces one of the two free bits (if any element below the new
one is above x_0, then `a(n,0)` must be 1; similarly for column 1).

## 2. Enumeration engines

* **Library** (`include/riordan/binary_riordan.hpp`, `RiordanBuilder`): rows
  and columns as dynamic bitsets; a forced entry is the parity of
  `column_{j−1} AND reversed(f)` — one AND and one popcount per 64-bit word.
  Works for any n.
* **`rcount`** (`tools/rcount.cpp`, `struct Node`): the same recurrence
  specialised to n ≤ 64, every row, column, g and f in a single machine word,
  no heap allocation, candidates in a fixed array. About twice as fast per node.
* **Independent enumerator** (`verify/independent/`): written separately with its
  own bitsets, series arithmetic and extension rule; used only for verification.

## 3. Counting up to isomorphism

r(n) needs, for every matrix, a **canonical form**: a certificate equal for two
matrices exactly when their posets are isomorphic.

* **Library canonical form** (`canon.hpp`): individualisation–refinement
  specialised to posets. Colour refinement uses the multisets of colours of the
  strict down-set and up-set; when refinement stalls, the first non-singleton
  cell is individualised; *twins* (equal strict down- and up-sets) are
  interchangeable, so only one per cell is tried; the certificate is the
  lexicographically smallest relabelled relation matrix over the leaves.
* **nauty** (bundled, called from C++): the digraph of the strict order, with an
  initial ordered partition by `(down-set size, up-set size)`. Because this
  colouring is isomorphism-invariant, the canonical graph together with the
  sorted colour sequence is a valid certificate. The colouring gives large
  speed-ups on posets with many interchangeable elements.

Each certificate is reduced to a 128-bit hash; with k ≤ 10⁹ certificates the
probability of any collision is below k²/2¹²⁹ ≈ 10⁻²¹.

## 4. Scaling: partitioning, memory and checkpoints

**Disk mode (default).** Hashes are written into B bucket files chosen by hash;
each bucket is sorted and deduplicated separately. Memory O(m(n)/B), disk 16
bytes per matrix.

**Pass mode (`--passes P`).** Choose an isomorphism *invariant* — here the sorted
multiset of `(down-set size, up-set size)` over all elements — and hash it into
P classes. Pass p enumerates the whole tree but canonicalises only matrices whose
invariant falls in class p, deduplicating in memory. Since isomorphic matrices
share the invariant, every isomorphism class is counted in exactly one pass, and
every matrix is canonicalised exactly once overall; only the (cheap) enumeration
is repeated. Memory ≈ 16 bytes × r(n)/P × ~2.5, no disk.

Two refinements make pass mode efficient:

* **last-level shortcut**: the invariant of each leaf is computed from the
  parent's row/column counts plus the candidate row, so only the 1-in-P leaves
  of the current pass are built;
* **bounded compaction**: the in-memory buffer is sorted and deduplicated, then
  regrown by a bounded step (exact `reserve`), avoiding the 2× peak of vector
  doubling.

**Checkpoints (`--checkpoint FILE`).** After each pass the line
`N P pass r_pass m` is appended; a rerun skips recorded passes. Every pass
re-enumerates all m(n) matrices, so the m column must agree across passes — a
built-in consistency check (the tools abort otherwise). The n = 27 and n = 28
results were obtained across many interrupted sessions this way.

These are standard techniques (orderly generation by canonical augmentation is
the general method; partitioning by an invariant and external-memory
deduplication are common in large isomorphism censuses). What is specific here
is the combination with facts F1–F4, which turn the Riordan structure into a
2-way branching tree with constant-time-per-word steps.

## 5. Parallelisation

The tree is split at depth D (`--split`, default 12): every node of order D is a
job, recorded as its sequence of choices `(a0, a1)` from the root. Worker threads
take jobs from a shared atomic counter, replay the choices and finish the subtree
independently, writing to private buffers (no locks on the hot path). nauty is
compiled with thread-local storage.

Load balance: the tree is very unbalanced — split at depth 6, one job holds 97% of
the leaves; at depth 12 no job holds more than 2.5%, and the heaviest jobs come
first in depth-first order, which suits dynamic scheduling. Hence the default
scales to about 30–40 threads; use `--split 16` beyond that. (Development and all
runs reported here used a single core; the multi-thread code paths are exercised by
`scripts/check_rcount.sh` and `scripts/check_independent.sh`, which compare runs with
2-4 threads against the known values, but large multi-core speed-ups have not been
measured.)

The independent counter (`verify/independent_rcount.cpp`) is parallelised without
touching its enumerator: every thread runs its own copy and walks the top of the tree,
and at the split depth each node is claimed with an atomic exchange by the first thread
to reach it; the others prune it through the enumerator's ordinary "skip this subtree"
callback. Since the depth-first order is deterministic, node k is the same node in every
thread, and threads that finish early simply claim more nodes.

## 6. Verification design

| method | enumerator | canonical form | used for |
|---|---|---|---|
| A | library | library | n ≤ 23 |
| B | independent | nauty (Python) | n ≤ 20 |
| C | `rcount` Node | nauty (C++) | n ≤ 28 |
| D | independent | nauty (C++), other hash, other invariant split | n ≤ 28 (24–28 essential) |

Every published value has two methods with different enumerators; up to n = 23
also two different canonical forms. In addition: brute force over all
`2^(2n−3)` matrices for n ≤ 11; p(n) (all posets, OEIS A000112) reproduced for
n ≤ 9 by the same canonical form; Python brute force for n ≤ 6.

## 7. Complexity

Let w = 64 and n ≤ 64.

* **Per tree node**: forced row O(n) word operations (n parities), closure test
  O(n) words, push/pop O(n). Library builder: O(n²/w) generally.
* **Tree size**: Σ_{k≤n} m(k) ≈ 2 m(n) nodes, with m(n) ≈ 2.03ⁿ … 2.06ⁿ
  empirically (no proven bound better than 4ⁿ).
* **Canonical forms**: individualisation–refinement is exponential in the
  worst case but cheap here; measured ≈ 4–8 µs per matrix with nauty, 20 µs with
  the library canonical form, at n = 20–28.
* **Total**: Θ(m(n) · c(n)) + P × enumeration, where c(n) is the canonical-form
  cost. Measured single-core times with `rcount`: n = 20: 10 s; 24: 4 min;
  26: 27 min; 27: 40 min (P = 5); 28: ≈ 2.5 h (P = 8). Each step roughly
  doubles the time.
* **Memory**: builder O(n²) bits; deduplication 16 bytes per class per pass
  (pass mode) or per matrix on disk (disk mode).

## 8. Exact decision procedure for infinite Riordan posets

For rational `g = a/b`, `f = c/d` over GF(2), `G(t,u) = g/(1−uf) = P/Q` with
`P = ad`, `Q = bd + u·bc`. Sections `Λ_{r,s}` (residues of i, j mod 2) satisfy
`Λ_{r,s}(R/Q) = Λ_{r,s}(RQ)/Q` because `Q² = Q(t², u²)` in characteristic 2, and
numerator degrees stay bounded. So the reachable numerators form a finite
automaton reading the binary digits of (i, j) and outputting `a(i,j)`.
Transitivity fails iff the product of three copies (reading (i,j), (j,k), (i,k))
reaches a state with outputs (1,1,0): a finite search that returns a proof or a
witness. See `infinite.hpp` and `docs/RESULTS.md`, Section 5.
