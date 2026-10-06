// binary_riordan.hpp - lower-triangular bit matrices and {g,f}_n (Phase 2 of the plan).
//
// Convention (as in Cheon et al., LAA 632 (2022)): rows/columns start at 0 and
// a(i,j) = 1  <=>  x_j <= x_i.  The matrix is lower-triangular; its transpose is
// the zeta matrix of the poset, and every labeling is a linear extension.
#pragma once
#include <string>
#include <vector>

#include "gf2_series.hpp"

namespace riordan {

class LowerBitMatrix {
 public:
  LowerBitMatrix() = default;
  explicit LowerBitMatrix(std::size_t n) {
    rows_.reserve(n);
    for (std::size_t i = 0; i < n; ++i) rows_.emplace_back(i + 1);
  }
  /// Rows written as 0/1 strings of length n, e.g. {"100","110","101"}.
  static LowerBitMatrix from_rows(const std::vector<std::string>& rows) {
    const auto n = rows.size();
    LowerBitMatrix A(n);
    for (std::size_t i = 0; i < n; ++i) {
      if (rows[i].size() != n) throw std::invalid_argument("from_rows: matrix must be square");
      for (std::size_t j = 0; j < n; ++j) {
        if (rows[i][j] != '1') continue;
        if (j > i) throw std::invalid_argument("from_rows: matrix is not lower-triangular");
        A.set(i, j);
      }
    }
    return A;
  }

  std::size_t order() const { return rows_.size(); }
  bool operator()(std::size_t i, std::size_t j) const { return j <= i && rows_[i].test(j); }
  void set(std::size_t i, std::size_t j, bool v = true) {
    if (j > i) throw std::out_of_range("LowerBitMatrix::set above the diagonal");
    rows_[i].set(j, v);
  }
  const Bits& row(std::size_t i) const { return rows_[i]; }
  Bits column(std::size_t j) const {
    Bits c(order());
    for (std::size_t i = j; i < order(); ++i)
      if (rows_[i].test(j)) c.set(i);
    return c;
  }
  LowerBitMatrix leading(std::size_t m) const {
    LowerBitMatrix B;
    B.rows_.assign(rows_.begin(), rows_.begin() + static_cast<std::ptrdiff_t>(std::min(m, order())));
    return B;
  }
  bool operator==(const LowerBitMatrix& o) const { return rows_ == o.rows_; }

  std::string str() const {
    std::string s;
    for (std::size_t i = 0; i < order(); ++i) {
      for (std::size_t j = 0; j < order(); ++j) s += (*this)(i, j) ? '1' : '0';
      s += '\n';
    }
    return s;
  }

 private:
  std::vector<Bits> rows_;
};

/// The truncated binary Riordan matrix {g,f}_n with a(i,j) = [t^i] g f^j mod 2.
/// Built column by column: column 0 is g and column j+1 is (column j) * f.
inline LowerBitMatrix binary_riordan(const Gf2Series& g, const Gf2Series& f, std::size_t n) {
  if (g.prec() < n || f.prec() < n)
    throw std::invalid_argument("binary_riordan: g and f need precision >= n");
  if (n >= 1 && (!g[0] || f[0]))
    throw std::invalid_argument("binary_riordan: need g(0) = 1 and f(0) = 0");
  if (n >= 2 && !f[1]) throw std::invalid_argument("binary_riordan: need f = t + O(t^2)");
  LowerBitMatrix A(n);
  Gf2Series col = g.truncated(n);
  const Gf2Series ft = f.truncated(n);
  for (std::size_t j = 0; j < n; ++j) {
    col.bits().for_each([&](std::size_t i) { A.set(i, j); });  // valuation(g f^j) = j
    if (j + 1 < n) col = col * ft;
  }
  return A;
}

/// Lucas: C(i,j) is odd iff the binary digits of j are a subset of those of i.
inline bool pascal_mod2(std::size_t i, std::size_t j) { return (j & ~i) == 0; }

}  // namespace riordan
