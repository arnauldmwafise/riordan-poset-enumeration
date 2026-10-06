// poset.hpp - poset-matrix tests and the matrix <-> (g,f) correspondence (Phase 2).
#pragma once
#include <optional>

#include "binary_riordan.hpp"

namespace riordan {

/// Witness that transitivity fails: a(i,j) = a(j,k) = 1 but a(i,k) = 0,
/// i.e. x_k <= x_j <= x_i while x_k <= x_i is missing.
struct Violation {
  std::size_t i, j, k;
};

/// Scans rows in order and returns the first transitivity failure. Only the row
/// of the top element can fail for a chain i >= j >= k, so each row is tested
/// against the rows below it: row j must be a subset of row i whenever a(i,j) = 1.
inline std::optional<Violation> first_violation(const LowerBitMatrix& A) {
  for (std::size_t i = 0; i < A.order(); ++i) {
    const Bits& ri = A.row(i);
    for (auto j = ri.next(0); j < i; j = ri.next(j + 1)) {
      const auto k = first_and_not(A.row(j), ri);
      if (k < A.row(j).size()) return Violation{i, j, k};
    }
  }
  return std::nullopt;
}

inline bool has_unit_diagonal(const LowerBitMatrix& A) {
  for (std::size_t i = 0; i < A.order(); ++i)
    if (!A(i, i)) return false;
  return true;
}

/// Theorem 2.1: a unit lower-triangular 0/1 matrix is a poset matrix iff it is transitive.
inline bool is_poset_matrix(const LowerBitMatrix& A) {
  return has_unit_diagonal(A) && !first_violation(A);
}

/// True if the candidate last row r (size n+1) is transitively closed over rows 0..n-1.
inline bool row_is_closed(const std::vector<Bits>& rows, const Bits& r) {
  const std::size_t n = r.size() - 1;
  for (auto j = r.next(0); j < n; j = r.next(j + 1))
    if (!is_subset(rows[j], r)) return false;
  return true;
}

/// Poset horizon N*(g,f): the largest n such that {g,f}_n is a poset matrix.
/// Leading submatrices of poset matrices are poset matrices (Theorem 3.6), so if
/// row i is the first to fail, the horizon is exactly i.
struct Horizon {
  std::size_t order;                   // largest poset order found (== cap if none failed)
  bool reached_cap;                    // true: no failure up to cap (infinite not proven!)
  std::optional<Violation> witness;    // first failure, if any
};
inline Horizon poset_horizon(const Gf2Series& g, const Gf2Series& f, std::size_t cap) {
  const auto v = first_violation(binary_riordan(g, f, cap));
  if (!v) return {cap, true, std::nullopt};
  return {v->i, false, v};
}

/// The pair (g,f), truncated to the order of the matrix, that generates a binary Riordan matrix.
struct RiordanPair {
  Gf2Series g, f;
};

/// Inverse map: g = column 0 and f = (column 1) * g^{-1}. The order-n matrix fixes
/// g and f modulo t^n, so the answer is unique. Returns nullopt if A is not
/// {g,f}_n for any pair (checked by rebuilding).
inline std::optional<RiordanPair> recover_gf(const LowerBitMatrix& A) {
  const auto n = A.order();
  if (n == 0) return std::nullopt;
  Gf2Series g(A.column(0));
  if (!g[0]) return std::nullopt;
  Gf2Series f(n);
  if (n >= 2) {
    f = Gf2Series(A.column(1)) * g.inverse();
    if (f[0] || !f[1]) return std::nullopt;
  }
  if (!(binary_riordan(g, f, n) == A)) return std::nullopt;
  return RiordanPair{g, f};
}

/// Dual poset (Theorem 2.8): A* = E A^T E, i.e. A*(a,b) = A(n-1-b, n-1-a).
inline LowerBitMatrix flip_transpose(const LowerBitMatrix& A) {
  const auto n = A.order();
  LowerBitMatrix B(n);
  for (std::size_t i = 0; i < n; ++i)
    A.row(i).for_each([&](std::size_t j) { B.set(n - 1 - j, n - 1 - i); });
  return B;
}

/// Direct sum P + Q (block diagonal) and ordinal sum P (+) Q (all-ones block below), eq. (3).
inline LowerBitMatrix direct_sum(const LowerBitMatrix& P, const LowerBitMatrix& Q, bool ordinal = false) {
  const auto m = P.order(), n = Q.order();
  LowerBitMatrix S(m + n);
  for (std::size_t i = 0; i < m; ++i) P.row(i).for_each([&](std::size_t j) { S.set(i, j); });
  for (std::size_t i = 0; i < n; ++i) {
    Q.row(i).for_each([&](std::size_t j) { S.set(m + i, m + j); });
    if (ordinal)
      for (std::size_t j = 0; j < m; ++j) S.set(m + i, j);
  }
  return S;
}
inline LowerBitMatrix ordinal_sum(const LowerBitMatrix& P, const LowerBitMatrix& Q) {
  return direct_sum(P, Q, true);
}

/// Theorem 3.6: deleting the first n rows and columns of {g,f} gives {g f^n / t^n, f}.
inline RiordanPair shifted_pair(const Gf2Series& g, const Gf2Series& f, std::size_t n) {
  const Gf2Series gfn = g * f.pow(static_cast<long long>(n));
  const Gf2Series g2 = gfn.shift_down(n);
  return {g2, f.truncated(g2.prec())};
}

}  // namespace riordan
