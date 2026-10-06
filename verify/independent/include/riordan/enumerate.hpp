// enumerate.hpp - automatic construction of all Riordan poset matrices.
//
// Extending {g,f}_n to {g,f}_{n+1} adds row n. Because g f^j = (g f^{j-1}) f,
//     a(n,j) = sum_k a(k,j-1) f_{n-k}  (mod 2)   for j >= 2,
// so columns 2..n of the new row are forced by data already present. Only
//     a(n,0) = g_n   and   a(n,1) = f_n + sum_{k=1}^{n-1} g_k f_{n-k}
// are free, giving at most 4 children per node (2 when n = 1, where column 1
// is the diagonal). The new element has the largest label, so transitivity only
// has to be checked on the new row. Poset matrices are closed under taking
// leading submatrices, so pruning a failed child loses nothing.
#pragma once
#include <vector>

#include "poset.hpp"

namespace riordan {

class RiordanPosetEnumerator {
 public:
  explicit RiordanPosetEnumerator(std::size_t max_order) : N_(max_order) {}

  /// visit(order, g_bits, f_bits, rows) is called once for every Riordan poset matrix
  /// of order 1..max_order. g_bits/f_bits hold g_0..g_{order-1} and f_0..f_{order-1}.
  /// Return false from visit to skip that node's subtree.
  template <class Visit>
  void run(Visit&& visit) {
    if (N_ == 0) return;
    rows_.clear();
    cols_.assign(N_, Bits(N_));
    g_ = Bits(N_);
    f_ = Bits(N_);
    rows_.emplace_back(1);
    rows_[0].set(0);
    cols_[0].set(0);
    g_.set(0);
    dfs(visit);
  }

 private:
  template <class Visit>
  void dfs(Visit& visit) {
    const std::size_t n = rows_.size();
    if (!visit(n, g_, f_, rows_) || n == N_) return;

    Bits R(N_);  // R_k = f_{n-k} for k = 1..n-1 (R_0 stays 0: f_n is still unknown)
    for (std::size_t k = 1; k < n; ++k)
      if (f_.test(n - k)) R.set(k);
    Bits forced(n + 1);
    forced.set(n);
    for (std::size_t j = 2; j < n; ++j)
      if (parity_and(cols_[j - 1], R)) forced.set(j);
    const bool conv = parity_and(cols_[0], R);

    const int col1_choices = n >= 2 ? 2 : 1;
    for (int gb = 0; gb < 2; ++gb) {
      for (int c = 0; c < col1_choices; ++c) {
        Bits r = forced;
        r.set(0, gb != 0);
        if (n >= 2) r.set(1, c != 0);
        if (!row_is_closed(rows_, r)) continue;
        const bool fn = n >= 2 ? ((c != 0) != conv) : true;
        push(r, n, gb != 0, fn);
        dfs(visit);
        pop(r, n);
      }
    }
  }
  void push(const Bits& r, std::size_t n, bool gn, bool fn) {
    rows_.push_back(r);
    r.for_each([&](std::size_t j) { cols_[j].set(n); });
    g_.set(n, gn);
    f_.set(n, fn);
  }
  void pop(const Bits& r, std::size_t n) {
    rows_.pop_back();
    r.for_each([&](std::size_t j) { cols_[j].set(n, false); });
    g_.set(n, false);
    f_.set(n, false);
  }

  std::size_t N_;
  std::vector<Bits> rows_, cols_;
  Bits g_, f_;
};

/// Convenience: the first n coefficients of a coefficient bitset as a series.
inline Gf2Series to_series(const Bits& b, std::size_t n) {
  Gf2Series s(n);
  b.for_each([&](std::size_t i) { if (i < n) s.set(i); });
  return s;
}

inline LowerBitMatrix to_matrix(const std::vector<Bits>& rows) {
  LowerBitMatrix A(rows.size());
  for (std::size_t i = 0; i < rows.size(); ++i) rows[i].for_each([&](std::size_t j) { A.set(i, j); });
  return A;
}

}  // namespace riordan
