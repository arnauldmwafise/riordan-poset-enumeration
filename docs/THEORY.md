# Theory notes behind the implementation

Conventions follow Cheon, Curtis, Kwon, Mesinga Mwafise, *Riordan posets and
associated incidence matrices*, Linear Algebra Appl. 632 (2022) 308–331.
A poset matrix has `a_ij = 1` iff `x_j <= x_i`; binary Riordan matrices are
`{g, f} = [a_ij]`, `a_ij = [t^i] g f^j mod 2`, with `g = 1 + g_1 t + ...` and
`f = t + f_2 t^2 + ...`. `{g, f}_n` is the leading `n x n` block.

The three facts below are not stated in the paper; they are elementary
consequences of its definitions and are what make automatic construction
cheap. Each is also verified by the test suite.

## F1. Matrix ↔ generating functions is a bijection

`a_{i,0} = g_i`, so column 0 is `g mod t^n`. Column 1 is `g f mod t^n`, and
`g` is invertible, so `f mod t^n = (column 1) · g^{-1}`. Conversely `a_ij` for
`i < n` only involves coefficients of index `<= i`. Hence
`(g mod t^n, f mod t^n) ↔ {g, f}_n` is a bijection. The free coefficients are
`g_1..g_{n-1}` and `f_2..f_{n-1}`, so there are exactly `2^(2n-3)` binary
Riordan matrices of order `n >= 2`.

*Consequence:* the generating functions of any constructed Riordan poset are
read off its matrix (`recover_pair`); they never have to be guessed.

## F2. Each new row has exactly two free entries

From `g f^j = (g f^{j-1}) f`,

    a_{n,j} = sum_{k=0}^{n} a_{k,j-1} f_{n-k}.

Since `f_0 = 0` the `k = n` term vanishes, and for `j >= 2` the term with
`f_n` would need `k = 0 < j-1`, where `a_{0,j-1} = 0`. So for `j >= 2` the
entries of row `n` depend only on rows `0..n-1` and `f_2..f_{n-1}`:
they are **forced**. The remaining two entries are

    a_{n,0} = g_n,     a_{n,1} = f_n + sum_{k=1}^{n-1} g_k f_{n-k},

and each toggles freely with the new coefficient (`g_n` resp. `f_n`).
For `n = 1` column 1 is the diagonal, so only `g_1` is free.

`RiordanBuilder::forced_next_row` evaluates the forced entries as parities of
`column_{j-1} AND reverse(f)`, `O(n^2/64)` per row.

## F3. Transitivity only has to be checked on the new row

If `a_ij = a_jk = 1` in a lower triangular matrix then `i >= j >= k`. Adding
row `n` (the largest label) creates new triples only with `i = n`, so the
extended matrix is transitive iff the old one is and `row_j ⊆ row_n` for every
`j` with `a_{n,j} = 1`.

Rows 0 and 1 are supported on columns `{0, 1}`, so the choice of `a_{n,0}`,
`a_{n,1}` cannot repair a violation among columns `>= 2`; such nodes are
pruned before trying the four choices.

## Enumeration tree

Leading principal submatrices of a transitive lower triangular matrix are
transitive (Theorem 3.6 states this for Riordan poset matrices). With F1–F3,
the Riordan poset matrices of order `<= N` form a rooted tree in which each
node has at most 4 children `(g_n, f_n)`, every node is a Riordan poset
matrix, and every Riordan poset matrix appears exactly once.
`for_each_riordan_poset` walks it depth-first.

## Recognition

A Riordan labelling is lower triangular, hence a linear extension. The search
assigns labels `0, 1, 2, ...`; placing element `x` at label `d` fixes row `d`,
whose columns `2..d-1` must equal the forced bits from F2 (a necessary
condition), after which columns 0 and 1 determine `(g_d, f_d)`. Two elements
with the same strict down-set and up-set give identical matrices when
swapped, so one representative per such twin class is tried. The search is
otherwise exhaustive, so a negative answer is a proof that the poset is not
Riordan.

## Canonical forms

Individualization–refinement on the order relation, with colour refinement on
the multisets of colours of strict down- and up-sets, branching on the first
non-singleton cell, and pruning of twins (a twin transposition is an
automorphism fixing every previously individualised vertex). The minimum leaf
certificate is canonical.

## Notes on the paper found while implementing

* **Proof of Theorem 5.1.** The paper writes `P4 = (1 ⊕ (1+1)) + 1`.
  The binary Pascal matrix gives `P4 = B_2 = 1 ⊕ (1+1) ⊕ 1`, which is also
  what the paper's own `P5 = 1 ⊕ (((1+1) ⊕ 1) + 1)` requires. The theorem's
  conclusion (P_n is series-parallel iff n <= 5) is confirmed for n <= 20.
* **Theorem 2.8.** The formula for the dual pair is garbled in some PDF
  extractions; the library obtains the dual by flip-transpose plus F1, and the
  tests confirm that the flip-transpose of every Riordan poset matrix of order
  <= 9 is again a Riordan poset matrix of the dual poset.
  For `P5` it returns `((1+t)^4, t/(1-t))`, as in the paper.
* **Theorem 5.3 (resolved in Phase 6: the statement is wrong, the proof is right).** The proof says `g_1 = α + 1 mod 2`,
  which with Theorem 5.2 makes *odd* `α` give the antichain, while the
  statement assigns the antichain to *even* `α`.
