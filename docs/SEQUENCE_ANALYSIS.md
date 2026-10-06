# Analysis of the sequences m(n) and r(n)

Data: `data/riordan_poset_counts.csv`, n = 1..28 (every term confirmed by two
independent programs). Everything below is reproduced by
`python3 scripts/analyze_sequences.py` (needs numpy and sympy).

* m(n) = number of n × n binary Riordan matrices that are poset matrices;
* r(n) = number of non-isomorphic posets among them.

```
m: 1, 2, 7, 22, 55, 121, 214, 475, 970, 2122, 4316, 8876, 16443, 36281, 72909, 162536, 321419,
   656225, 1198005, 2634534, 5243230, 11497714, 22571793, 47110953, 86046455, 186702111,
   364554584, 799401460
r: 1, 2, 5, 11, 33, 74, 144, 232, 639, 1406, 3164, 4992, 12501, 26973, 55937, 104169, 266428,
   517316, 941970, 1742608, 4424434, 9366627, 18579721, 32869555, 73299324, 153092515,
   301274742, 572885293
```

## 1. Searches for exact structure (all negative)

A relation is only reported if it leaves at least three equations unused, so that
it predicts terms it was not fitted to. All arithmetic is exact.

| test | m(n) | r(n) |
|---|---|---|
| linear recurrence with constant coefficients (Berlekamp–Massey over Q) | minimal order 14 = N/2 | 14 = N/2 |
| P-recursive recurrence, order ≤ 6, polynomial degree ≤ 4 | none | none |
| algebraic generating function, degree ≤ 3, coefficient degree ≤ 7 | none | none |

Berlekamp–Massey returning exactly N/2 is what a sequence *without* a
linear recurrence of order < N/2 produces; the data therefore rule out any
rational generating function of denominator degree ≤ 13, and the other searches
rule out small holonomic and algebraic descriptions. This is evidence, not proof:
a relation of higher order or degree cannot be excluded with 28 terms.

There are also structural reasons to expect no simple formula: the sequences
show arithmetic effects tied to the binary expansion and to n mod 6 and n mod 4
(Section 3), as one would expect from objects defined by GF(2) arithmetic
(Lucas' theorem governs the Pascal poset, and binary Riordan arrays with rational
g, f are 2-automatic). Generating functions of such quantities are typically
not D-finite.

## 2. Growth rate

Fitting `log a(n) = n log α + β log n + c(n mod 6)` by least squares (residue-class
offsets absorb the period-6 oscillation described below):

| sequence, range | α (pure exponential) | α, β (with n^β) |
|---|---|---|
| m, n ≥ 8 | 2.042 | 2.006, +0.30 |
| m, n ≥ 12 | 2.037 | 1.996, +0.39 |
| m, n ≥ 16 | 2.033 | 2.014, +0.20 |
| r, n ≥ 16 | 2.047 | 1.89, +1.7 (unstable) |

The six-step rate `(m(n)/m(n−6))^(1/6)` decreases steadily, from 2.062 at n = 13
to 2.028 at n = 28. The ratio `m(n)/(2^n log₂ n)` stays between 0.54 and 0.62
for 8 ≤ n ≤ 28 (oscillating with period 6), while `m(n)/2^n` grows slowly from
1.9 to 3.0.

**Conjecture 1.** `lim m(n)^(1/n) = 2`. More precisely `m(n) = 2^n φ(n)` with a
subexponential, slowly growing φ; both `φ(n) ≈ 0.6 log₂ n` and `φ(n) ≈ c n^β`,
β ≈ 0.2–0.4, fit the data, and 28 terms cannot separate them.

For r(n) the fits are noisier because of the strong mod-4 effect below, but
`r(n)/m(n)` increases within every residue class mod 4 and lies in [0.49, 0.86]
for n ≥ 8, so r(n) and m(n) have the same exponential growth if the ratio stays
bounded away from 0, which the data strongly suggest.

**Conjecture 2.** `lim r(n)^(1/n) = 2`, and `r(n)/m(n)` converges (separately
for n ≡ 0 mod 4 and n ≢ 0 mod 4) to a positive limit, plausibly about 0.85 for
n ≢ 0 mod 4.

Proven bounds are much weaker: `r(n) ≤ m(n) ≤ 2^(2n−3)` and `r(n) ≥ m(n)/n!`.
Closing the gap between 4ⁿ and the observed 2ⁿ is an open problem.

## 3. Where the 2ⁿ comes from (an empirical decomposition)

Group the m(n) matrices by their series f. Very few f occur at all, and each
admits about equally many g:

| n | m(n) | distinct f | possible f (2^(n−2)) | mean g per f | most g for one f |
|---|---|---|---|---|---|
| 12 | 8876 | 144 | 1024 | 61.6 | 131 (f = t) |
| 16 | 162536 | 791 | 16384 | 205 | 580 (f = t) |
| 20 | 2634534 | 3952 | 262144 | 667 | 2616 (f = t) |
| 22 | 11497714 | 11276 | 1048576 | 1020 | 5344 (f = t) |

From n = 16 to 22 the number of admissible f grows like 1.56ⁿ and the mean number
of g per f like 1.31ⁿ; the product, 2.04ⁿ, is the growth of m(n). The most
productive f is always f = t, whose g are exactly the Toeplitz posets of
Section 4 of the paper (supports closed under addition). This suggests studying
the two factors separately, e.g. bounding the number of admissible f.

## 4. Oscillations

* **Period 6 in m(n).** The step ratio `m(n)/m(n−1)` is lowest at
  n = 7, 13, 19, 25 (1.77–1.85) and highest at n = 8, 14, 20, 26 and
  n = 10, 16, 22, 28 (≈ 2.2).
* **Period 4 in r(n)/m(n).** At multiples of 4 the ratio is much lower
  (0.49, 0.56, 0.64, 0.66, 0.70, 0.72 at n = 8, …, 28) than in the other residue
  classes (0.82–0.85 at n ≈ 25). At n ≡ 0 mod 4 many posets therefore have
  several Riordan labellings; Pascal posets are Boolean lattices exactly at powers
  of 2, which is one source of extra symmetry, but we have no full explanation.

These effects are why no low-order recurrence exists and why growth fits need
residue-class terms.

## 5. Summary of what the data do and do not support

| statement | status |
|---|---|
| no linear recurrence of order ≤ 13 | proved by the data (exact) |
| no small P-recurrence / algebraic GF | shown for the ranges tested (exact) |
| growth rate 2 for m(n) and r(n) | conjecture, strong numerical evidence |
| polynomial vs logarithmic correction factor | undetermined |
| periodic modulation (6 and 4) | observed; unexplained |
| closed form | none found; unlikely given the above |
