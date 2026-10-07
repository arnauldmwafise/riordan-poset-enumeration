# Computational results on Riordan posets

This note records what the library established about the Riordan posets of
G.-S. Cheon, B. Curtis, G. Kwon and A. Mesinga Mwafise, *Riordan posets and
associated incidence matrices*, Linear Algebra Appl. 632 (2022) 308–331
(below: "the paper"). Every number is reproducible with the commands in the
last section. The verification column says how each claim was checked; a
claim marked *two independent methods* was obtained by two programs that share
no code.

Convention throughout: `a_ij = 1` iff `x_j <= x_i`, indices from 0,
`{g,f}_n = [[t^i] g f^j mod 2]`, `P_n` the Pascal poset, `P` the infinite one,
`I_n` the antichain, `+` direct sum, `⊕` ordinal sum, `Q*` the dual of `Q`.

## 1. Survey of Riordan posets

`r(n)` is the number of non-isomorphic Riordan posets on `n` elements; `m(n)`
the number of (labelled) Riordan poset matrices of order `n`.

| n | m(n) | r(n) | methods |
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

*  `n <= 8` agrees with the paper and `n <= 15` with OEIS A379608 (terms added in 2025-2026);
  r(n) for `16 <= n <= 28` is new. m(n) (labelled) does not appear in the paper or, as far as
  we could determine, in the OEIS (October 2026); all its terms are new.
* **Methods.** A: the library's tree enumerator (Fact F2 in `THEORY.md`) with
  its own individualization–refinement canonical form (`bin/riordan census`,
  `verify/census_lean.cpp`, or `bin/rcount --canon library`). B: a separately
  written enumerator (different code base, `verify/independent`) dumping every
  matrix, canonicalised by **nauty** through pynauty 2.8.8.1. C: `bin/rcount`,
  the library enumerator with nauty 2.8.8 called from C++ and disk-bucketed
  deduplication of 128-bit hashes of the canonical forms. A and C use different
  canonical forms; B shares no code with A or C. In B and C the nauty input is
  the digraph of the strict order with the isomorphism-invariant initial
  colouring by (down-set size, up-set size). D: `verify/independent_rcount.cpp`,
  the independent enumerator of B with nauty called from C++ and a different
  128-bit hash; it shares no code with A or C (also run at n = 5, 8, 12, 16, 20,
  where it matches).
* **Agreement.** Wherever two or three methods were run they agree exactly.
  `m(n)` was computed by both enumerators for `n <= 22`, and for `n <= 11` also
  by brute force over all `2^(2n-3)` binary Riordan matrices. Every row up to
  n = 28 was obtained by at least two methods built on different enumerators (A/C use the
  library's, B/D the independent one). For n <= 23 two different canonical
  forms also agree (library and nauty); n = 24-28 rely on nauty for the
  isomorphism step. The 128-bit hashes give a collision probability below
  `10^-20` even at `n = 25`.
* **Cost.** On one core `bin/rcount` takes 12.5 s for `n = 20`, 59 s for
  `n = 22`, 4 min for `n = 24` and 8 min for `n = 25`; time roughly doubles
  per step, and the work splits into independent jobs for parallel runs.
* As further checks of the isomorphism code: the same canonical form yields
  `p(n) = 1, 2, 5, 16, 63, 318, 2045, 16999, 183231` posets for `n <= 9`
  (OEIS A000112) and the 785 Toeplitz posets of order 17 (Section 3).

## 2. The open question at the end of Section 2 of the paper

The paper observes that the five non-Riordan posets on four elements are each a
direct or ordinal sum of two Riordan posets, and asks whether every non-Riordan
poset with `n >= 4` elements can be expressed in terms of Riordan posets.

**Answer: no. The first counterexamples occur at n = 5.**

Call a poset *indecomposable* if it is neither a direct sum nor an ordinal sum
of two nonempty posets. An indecomposable poset can be produced from smaller
posets by `+` and `⊕` in no way at all, so if it is not Riordan it is not
expressible in terms of Riordan posets under any reading of the question.

**Theorem 1 (computer-assisted).** There are exactly four indecomposable
non-Riordan posets on five elements. With Hasse diagrams given by cover
relations `a<b`:

| | cover relations | description |
|---|---|---|
| X1 | 0<2, 1<2, 2<3, 1<4 | the N poset `{0<2>1<4}` with an element above its upper-left vertex |
| X2 | 0<2, 1<2, 1<3, 1<4 | `1` has three upper covers, and `0<2` |
| X3 | 0<1, 1<3, 2<3, 1<4 | the dual of X1 |
| X4 | 2<3, 0<4, 1<4, 2<4 | the dual of X2: `4` covers three elements |

X1/X3 and X2/X4 are dual pairs (checked by canonical forms). Each contains N, as every indecomposable poset
on at least two elements must (series-parallel posets are exactly the N-free
ones, and they are decomposable).

*Proof of non-Riordan-ness.* For `n = 5` there are `2^(2·5−3) = 128` binary
Riordan matrices; by exhaustion none is a poset matrix isomorphic to X1–X4. This
was checked three ways: the census, the library's separate recognizer
(search over linear extensions), and an independent short Python script that
recomputes all 128 matrices from the definition and tests isomorphism by trying
all 120 permutations (`verify/decomp_check.py`).

Full classification for `n <= 9` (all posets up to isomorphism):

| n | posets | Riordan | non-Riordan | `A+B` or `A⊕B`, A, B Riordan | in the `+`/`⊕` closure of Riordan posets | indecomposable non-Riordan |
|---|---|---|---|---|---|---|
| 4 | 16 | 11 | 5 | 5 | 5 | 0 |
| 5 | 63 | 33 | 30 | 19 | 26 | 4 |
| 6 | 318 | 74 | 244 | 118 | 166 | 68 |
| 7 | 2045 | 144 | 1901 | 420 | 785 | 888 |
| 8 | 16999 | 232 | 16767 | 1252 | 3356 | 9897 |
| 9 | 183231 | 639 | 182592 | 3146 | 13625 | 126201 |

Here "closure" means the smallest class containing the Riordan posets and closed
under `+` and `⊕` (all binary splits tried, memoised over smaller posets). For
`n <= 6` every entry was confirmed by the independent Python brute force. So
under the strict reading (two Riordan summands) the property already fails for
11 of the 30 non-Riordan 5-element posets, under the closure reading for 4, and
the expressible fraction falls quickly (7.5% of non-Riordan posets at `n = 9`).

## 3. Toeplitz posets (Section 4) — Phase 5

* The number of Toeplitz poset matrices of order `n` equals the paper's list
  1, 2, 3, 5, 7, 12, 16, 27, 37, 58, 80, 131, 171, 277, 380, 580, 785 for
  `n <= 17`.
* **Theorem 4.3 holds up to isomorphism**: distinct supports never give
  isomorphic posets (checked by canonical forms for `n <= 17`), so the number of
  non-isomorphic Toeplitz posets of order `n` is exactly `|S(n)|`.
* Corollary 4.2 checked exhaustively for `n = 12` (all 2048 `g`); Corollary 4.4
  (self-duality) for every Toeplitz poset of order `<= 17`.
* Frobenius numbers by Nijenhuis' algorithm, checked against Sylvester's
  `ab − a − b`, brute force on 200 random triples, and known values
  `g(3,4) = 5`, `g(4,6,9) = 11`, `g(6,9,20) = 43`.

## 4. Exponential Riordan posets (Section 5) — Phase 6

The binary exponential matrix `[g,f]` is computed exactly in `Z/2` from the
integer recurrences `B_{m,k} = Σ C(m−1,i−1) F_i B_{m−i,k−1}` and
`b_{n,k} = Σ C(n,m) G_{n−m} B_{m,k}` (EGF coefficients `G_n`, `F_n`), and was
checked against exact rational arithmetic (`b_{n,k} = n!/k! [t^n] g f^k` with
Python `Fraction`) for 19 families and `n < 26`.

| Statement | Result |
|---|---|
| Thm 5.2 | Confirmed (`n <= 30`). |
| Stirling example | `S_5` confirmed, poset but not binary Riordan; `[e^t, e^t−1]_n` is a poset matrix exactly for `n <= 5`. |
| **Thm 5.3 (Laguerre)** | **Statement has the parities swapped.** `α` odd gives the antichain `I_n`; `α` even gives `(1⊕1)+(1⊕1)+…`. This is what the paper's own proof (`g_1 = α+1`) and Theorem 5.2 give; confirmed for `1 <= α <= 9`, `n <= 40`. |
| Thm 5.5 (Hermite) | Confirmed for `0 <= v <= 5`, `n <= 30`, including "series-parallel iff `n <= 10`" for odd `v`. |
| **Thm 5.6 (Jacobi cn, sn)** | Odd `m`: `P_⌈n/2⌉ + P_⌊n/2⌋` confirmed. **Even `m`: the printed `I_⌈n/2⌉ + P_⌊n/2⌋` is wrong for every odd `n`; the correct form is `I_⌊n/2⌋ + P_⌈n/2⌉`.** Smallest case by hand: `[cos t, sin t]_5 mod 2` has rows 10000, 01000, 10100, 00010, 10001, which is `I_2 + P_3`, not `I_3 + P_2`. Confirmed for `m = 0..4`, `n <= 60`; series-parallel exactly for `n <= 10`. |
| Thm 5.7 | Confirmed with the reading `D_m(i,j) = b_{i+1,j+1}` of `[cos, sin]` (the `t`-derivative shifts rows up; removing the first column and last row leaves this matrix): `D_m = I_m` for `m <= 5`, `I_{m−3} + P_3*` for `6 <= m <= 9`, not series-parallel for `m >= 10`. |

## 5. Infinite Riordan posets — Phase 7

### 5.1 Deciding infiniteness exactly

**Theorem 2.** For `g = a/b`, `f = c/d` rational over GF(2), it is decidable
whether `{g,f}_n` is a poset matrix for every `n`, and the algorithm
(`infinite.hpp`, `riordan decide`) returns either a proof or an explicit triple
`(i,j,k)` with `a_ij = a_jk = 1`, `a_ik = 0`.

*Proof.* `G(t,u) = Σ a_ij t^i u^j = g/(1 − uf) = P/Q` with `P = ad`,
`Q = bd + u·bc`. Let `Λ_{r,s}` extract the coefficients with `i ≡ r`,
`j ≡ s (mod 2)` and halve the exponents. Since `Q^2 = Q(t^2,u^2)` in
characteristic 2, `Λ_{r,s}(R/Q) = Λ_{r,s}(RQ)/Q`, and `deg_t R <= max(deg P, deg Q)`,
`deg_u R <= 1` are preserved. So the numerators `R` reachable from `P` form a
finite automaton that reads the binary digits of `(i,j)` (least significant
first) and outputs `R(0,0) = a_ij` (Christol/Salon). Transitivity is a
universal statement over triples; run three copies on `(i,j)`, `(j,k)`, `(i,k)`
in lockstep over digit triples. The finite product graph has a reachable state
with outputs `(1,1,0)` iff transitivity fails somewhere, and the path to it
spells out the witness. ∎

Validation: on 3000 random rational pairs (degrees <= 4) the automaton's
entries match the matrix (300 pairs × 160²), every "infinite" verdict (62)
holds to order 2048, and every one of the 2938 witnesses is a genuine violation.
It proves the paper's infinite families automatically (Thm 3.7, Cor 3.8 for
`k <= 12`, chains, `C_n + C_n`, Toeplitz `Γ(S)`), with automata of 2–7 states.

### 5.2 A new infinite family

**Theorem 3.** For every integer `k >= 0`,
`{(1+t)^k, t/(1+t)} = P*_{k+1} + P`,
the disjoint union of the dual of the finite Pascal poset `P_{k+1}` and the
infinite Pascal poset. In particular it is an infinite Riordan poset.

*Proof.* `a_ij = [t^{i−j}](1+t)^{k−j} mod 2`. For `i, j <= k` this is
`C(k−j, i−j)`, the dual `P*_{k+1}` (it is the flip-transpose of `P_{k+1}`). For
`i > k >= j`, `(1+t)^{k−j}` has degree `k − j < i − j`, so `a_ij = 0`. For
`i >= j > k`, `C(−(j−k), i−j) = ± C(i−k−1, i−j) = ± C(i−k−1, j−k−1)`, the Pascal
matrix shifted by `k+1`. ∎

Together with Corollary 3.8 (`k < 0`), `{(1+t)^k, t/(1+t)}` is infinite for
**every integer** `k`. For `k = n−1` the leading `n × n` block recovers the
paper's dual Pascal matrix (Theorem 2.8, `P_5* = ((1+t)^4, t/(1−t))_5`).

### 5.3 Classification for small rational degree

**Theorem 4 (computer-assisted).** Let `g = a/b`, `f = c/d` be reduced
fractions over GF(2) with `deg a, b, c, d <= 5`. Then `{g,f}` is an infinite
Riordan poset matrix if and only if

* `f = t` and `supp(g)` is a numerical semigroup (13 such `g`), or
* `f = t/(1+t)` and `g = (1+t)^k`, `−5 <= k <= 5`.

| max degree D | distinct pairs | infinite | largest finite horizon (attained by) |
|---|---|---|---|
| 2 | 66 | 9 | 7 (`g = 1`, `f = t/(1+t^2)`) |
| 3 | 946 | 13 | 12 (`g = 1`, `f = t/(1+t^3)`) |
| 4 | 14706 | 19 | 17 (`g = 1/(1+t^3)`, `f = t+t^4`) |
| 5 | 233586 | 24 | 30 (`g = 1`, `f = t/(1+t^5)`) |

Every pair was decided by Theorem 2 (each finite verdict comes with a verified
witness). The data suggest, but do not prove, that for rational `g, f` of any
degree the only infinite Riordan posets are the Toeplitz ones and
`{(1+t)^k, t/(1+t)}`. This is stated here as a conjecture.

## 6. Statistics of Riordan posets — Phase 7

Over the isomorphism classes counted by `r(n)`:

| n | r(n) | connected | series-parallel | self-dual | bounded (0̂ and 1̂) | mean height | mean width |
|---|---|---|---|---|---|---|---|
| 4 | 11 | 6 | 10 | 7 | 2 | 2.36 | 2.36 |
| 5 | 33 | 19 | 22 | 11 | 1 | 2.64 | 2.79 |
| 6 | 74 | 45 | 26 | 28 | 1 | 2.73 | 3.23 |
| 7 | 144 | 79 | 26 | 24 | 1 | 2.83 | 3.72 |
| 8 | 232 | 150 | 28 | 70 | 4 | 2.78 | 4.12 |
| 9 | 639 | 390 | 39 | 81 | 1 | 2.87 | 4.61 |
| 10 | 1406 | 862 | 52 | 332 | 1 | 2.78 | 5.15 |
| 11 | 3164 | 2093 | 56 | 200 | 1 | 2.82 | 5.52 |
| 12 | 4992 | 3509 | 56 | 894 | 1 | 2.73 | 6.03 |
| 13 | 12501 | 8340 | 72 | 679 | 1 | 2.87 | 6.56 |
| 14 | 26973 | 18115 | 84 | 3331 | 1 | 2.72 | 7.13 |

Observations: the mean height stays below 3 while the mean width grows roughly
like `n/2`; series-parallel Riordan posets become rare; and apart from the chain,
bounded Riordan posets occur only at `n = 4` and `n = 8` in this range (powers of
2, where `P_n` is the Boolean lattice `B_{log2 n}`).

## 7. Errata to the paper found by these computations

1. **Theorem 5.3**: parities of `α` swapped in the statement (proof is right).
2. **Theorem 5.6**: for even modulus the decomposition is `I_⌊n/2⌋ + P_⌈n/2⌉`.
3. **Proof of Theorem 5.1**: `P_4 = 1 ⊕ (1+1) ⊕ 1 = B_2`, not `(1 ⊕ (1+1)) + 1`
   (see `THEORY.md`); the theorem itself is correct.
4. **Theorem 5.7** does not define the "derivative" precisely; the reading above
   is the one that reproduces the stated result.

## Reproducing

```sh
make && make test                       # ~950k checks in 8 suites, a few seconds: every claim above in test range
./bin/riordan census 18                  # Section 1 (method A)
./bin/riordan express 9                  # Section 2
./bin/riordan toeplitz 17                # Section 3
./bin/riordan exp jacobi 0 5             # Section 4 example
./bin/riordan decide "1+t^2" "t/(1+t)"   # Section 5
python3 verify/nauty_count.py 1 20       # Section 1 (method B; needs pynauty and verify/dump)
./verify/census_lean 20                  # Section 1, r(20) by the library's method with low memory (build: see verify/README.md)
python3 verify/decomp_check.py           # Section 2, independent check for n <= 6
python3 verify/expcheck.py               # Section 4, exact rational check
```
