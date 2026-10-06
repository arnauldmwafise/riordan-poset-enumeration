# Changelog

All notable changes are recorded here. The project follows [Semantic Versioning](https://semver.org/).

## [1.1.0] - 2026-10-04

### Added
- `verify/independent_rcount.cpp` is multi-threaded (`--threads`, `--split`): threads claim
  subtrees at the split depth atomically, without modifying the independent enumerator.
- Build targets `bin/independent_rcount` (Makefile) and `independent_rcount` (CMake), compiled
  only against `verify/independent`.
- `scripts/check_independent.sh` regression check (n <= 14, several threads, both modes),
  registered with CTest and run by `make check-rcount`.

## [1.0.0] - 2026-10-04

First public release.

### Library (`include/riordan/`)
- GF(2) power series (product, Newton inverse, composition, reversion, expression parser).
- Binary Riordan matrices `{g,f}_n`; incremental builder with two free entries per row; Lucas test.
- Poset toolkit: transitivity witnesses, Hasse diagrams, Möbius function, height/width, duals,
  direct and ordinal sums, series-parallel decomposition, linear extensions.
- Canonical forms for posets (individualisation-refinement with twin pruning).
- Tree enumeration and census of Riordan posets; recognition of Riordan labellings.
- Toeplitz posets and numerical semigroups (closure, minimal generators, Frobenius numbers).
- Binary exponential Riordan matrices by exact mod-2 recurrences (Stirling, Laguerre, Hermite, Jacobi).
- 2-automata and an exact decision procedure for infinite Riordan posets with rational g, f.
- Direct/ordinal decompositions and a catalogue of all posets (expressibility question).

### Tools
- `bin/riordan`: build, info, horizon, dual, census, recognize, decide, toeplitz, exp, express.
- `bin/rcount`: large-scale counter (one-word-per-row enumerator, nauty, threads, disk buckets,
  invariant-partitioned passes, checkpoints).

### Data and results
- Verified counts m(n), r(n) for n <= 28 (`data/`), each confirmed by two independent programs.
- Answer to the open question of Section 2 of Cheon et al. (2022); errata to Theorems 5.3 and 5.6;
  new infinite family; sequence analysis.

### Quality
- Eight test suites (~950,000 checks) registered with CTest, `rcount` regression script,
  independent verification programs in `verify/`, CI workflow.
