# Independent verification

These scripts re-derive the published numbers by methods that share no code
with the library's census and canonical forms (except `census_lean.cpp`, which
re-runs the library's own method with little memory).

| Script | Checks | How |
|---|---|---|
| `dump.cpp` + `nauty_count.py` | r(n), m(n) for n <= 20 | `independent/` is a separately written enumerator (own bitsets, GF(2) series, transitivity test); canonical certificates from nauty via pynauty |
| `../tools/rcount.cpp` (`bin/rcount`) | r(n) for n <= 25 | library enumerator + nauty from C++ (default) or the library canonical form (`--canon library`); threads and disk buckets |
| `independent_rcount.cpp` | m(n), r(n) for n <= 28 | the independent enumerator (`independent/`) + nauty from C++, with its own hash; shares no code with the library or `rcount`; multi-threaded (`--threads T`, `--split D`) |
| `census_lean.cpp` | r(n) by the library's method | library enumerator and canonical form, keeping only a 128-bit hash of each certificate |
| `decomp_check.py` | expressibility table, n <= 6 | pure Python: all 2^(2n-3) binary Riordan matrices from the definition, all posets by brute force, isomorphism by trying every permutation |
| `expdump.cpp` + `expcheck.py` | exponential module | exact rational arithmetic (`fractions.Fraction`), b_{n,k} = n!/k! [t^n] g f^k |

```sh
pip install pynauty
g++ -std=c++20 -O2 -Iindependent/include dump.cpp -o dump && ./dump 20       # writes m1.bin .. m20.bin (~400 MB)
python3 nauty_count.py 1 20        # prints n, m(n), r(n), seconds; streams, n = 20 takes ~25 min
g++ -std=c++20 -O2 -I../include census_lean.cpp -o census_lean && ./census_lean 20   # -> 20 2634534 1742608
(cd .. && make bin/rcount && ./bin/rcount 22 && ./bin/rcount 22 --canon library)   # nauty vs library canonical form
(cd .. && make bin/rcount)   # also builds the nauty objects in ../bin/nauty
(cd .. && make bin/independent_rcount)   # or the CMake target independent_rcount
../bin/independent_rcount 25          # -> 25 86046455 73299324   (~10 min on one core, disk buckets)
../bin/independent_rcount 27 --passes 5 --checkpoint d27.ckpt   # in memory, resumable (its own invariant split)
python3 decomp_check.py
g++ -std=c++20 -O2 -I../include expdump.cpp -o expdump && ./expdump > exp_cpp.txt && python3 expcheck.py
```

Results obtained (all methods that were run agree; see the README for which method checked which n):

| n | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| m(n) | 72909 | 162536 | 321419 | 656225 | 1198005 | 2634534 | 5243230 | 11497714 | 22571793 | 47110953 | 86046455 | 186702111 | 364554584 | 799401460 |
| r(n) | 55937 | 104169 | 266428 | 517316 | 941970 | 1742608 | 4424434 | 9366627 | 18579721 | 32869555 | 73299324 | 153092515 | 301274742 | 572885293 |

## Running the independent counter on a multi-core machine

`independent_rcount` uses all cores by default (`--threads T` to choose). The
independent enumerator itself is unchanged: each thread runs its own copy, walks the
top of the tree, and at the split depth (`--split D`, default 12) atomically claims
nodes; only the claiming thread descends into a node. Work is therefore divided
dynamically without modifying the code being used for verification.

```sh
make bin/independent_rcount bin/rcount
make check-rcount                                   # both counters against known values
../bin/independent_rcount 28 --passes 8 --threads 32 --checkpoint d28.ckpt   # expect 799401460 572885293
../bin/independent_rcount 29 --passes 8 --threads 32 --checkpoint d29.ckpt
../bin/independent_rcount 30 --passes 16 --threads 32 --split 14 --checkpoint d30.ckpt
```

A new value is confirmed when `bin/rcount` and `bin/independent_rcount` report the same
m(n) and r(n) (the per-pass numbers differ by design: each tool splits by its own invariant).
