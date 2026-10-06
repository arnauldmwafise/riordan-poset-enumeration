# Contributing

Thank you for helping improve riordan-poset-enumeration. This page describes how the code is
organised and what a change needs before it is merged.

## Development workflow

```sh
make              # or: cmake -B build && cmake --build build
make test         # or: ctest --test-dir build
make check-rcount # regression check of the large-scale counter
```

All of these must pass, with no compiler warnings (`-Wall -Wextra -Wpedantic`), before
a pull request is merged; the CI workflow runs them with GCC and Clang.

## Code organisation

* `include/riordan/` is a header-only library, one module per file, all in namespace
  `riordan`. Modules depend only on modules listed before them in `riordan.hpp`.
* `tools/` contains the command-line programs, `examples/` small runnable examples,
  `verify/` independent re-implementations used to cross-check published numbers.
  Code in `verify/` must not include the library: its value is independence.

## Tests

* Every public function needs a test in the matching `tests/test_<module>.cpp`.
  Tests use the tiny framework in `tests/test_common.hpp`:
  `TEST(name) { CHECK(cond); CHECK_EQ(a, b); }`.
* Prefer tests against a naive reference implementation or against a statement of the
  paper (cite theorem numbers in the test name), and exhaustive checks for small n.
* A new test suite is picked up automatically by both the Makefile and CMake if its
  file name matches `tests/test_*.cpp`.

## Numbers that go into `data/` or the documentation

A new count is only recorded after two methods with *different enumerators* agree
(see `verify/README.md`), and the table states which methods were used. Do not round,
extrapolate or estimate in `data/`.

## Style

C++20, 4-space indentation, `snake_case` functions, `PascalCase` types, see
`.clang-format`. Comments explain *why*; mathematical routines cite the fact or theorem
they implement (`docs/THEORY.md`).
