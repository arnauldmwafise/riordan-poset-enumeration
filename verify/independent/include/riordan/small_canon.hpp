// small_canon.hpp - isomorphism-invariant canonical form for SMALL posets (n <= ~9).
// Used to reproduce the paper's counts r(n). For larger n, replace with nauty/bliss.
#pragma once
#include <algorithm>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

#include "enumerate.hpp"

namespace riordan {

/// Canonical form = lexicographically smallest relation code over all relabelings that
/// sort elements by an isomorphism invariant (down-set size, up-set size). Since the
/// invariant is preserved by isomorphisms, isomorphic posets get identical codes.
inline std::string small_canonical_form(const LowerBitMatrix& A) {
  const auto n = A.order();
  auto rel = [&](std::size_t u, std::size_t v) { return u >= v ? A(u, v) : false; };  // x_v <= x_u
  std::vector<std::pair<std::size_t, std::size_t>> key(n);
  for (std::size_t u = 0; u < n; ++u)
    for (std::size_t v = 0; v < n; ++v) {
      if (rel(u, v)) ++key[u].first;   // elements below u
      if (rel(v, u)) ++key[u].second;  // elements above u
    }
  std::vector<std::size_t> perm(n);
  std::iota(perm.begin(), perm.end(), 0);
  std::sort(perm.begin(), perm.end(), [&](auto a, auto b) { return std::tie(key[a], a) < std::tie(key[b], b); });

  std::vector<std::pair<std::size_t, std::size_t>> seg;  // groups with equal keys
  for (std::size_t s = 0; s < n;) {
    std::size_t e = s + 1;
    while (e < n && key[perm[e]] == key[perm[s]]) ++e;
    seg.emplace_back(s, e);
    s = e;
  }
  std::string best, code(n * n, '0');
  for (;;) {
    for (std::size_t a = 0; a < n; ++a)
      for (std::size_t b = 0; b < n; ++b) code[a * n + b] = rel(perm[a], perm[b]) ? '1' : '0';
    if (best.empty() || code < best) best = code;
    std::size_t g = seg.size();  // odometer over permutations within each group
    for (; g-- > 0;) {
      auto first = perm.begin() + static_cast<std::ptrdiff_t>(seg[g].first);
      auto last = perm.begin() + static_cast<std::ptrdiff_t>(seg[g].second);
      if (std::next_permutation(first, last)) break;
    }
    if (g == static_cast<std::size_t>(-1)) break;
  }
  return best;
}

/// All naturally labelled posets of order n (every unit lower-triangular transitive
/// 0/1 matrix), for brute-force cross-checks.
template <class Visit>
void for_each_poset_matrix(std::size_t n, Visit&& visit) {
  std::vector<Bits> rows;
  auto rec = [&](auto&& self) -> void {
    const auto m = rows.size();
    if (m == n) { visit(to_matrix(rows)); return; }
    for (std::size_t mask = 0; mask < (std::size_t{1} << m); ++mask) {
      Bits r(m + 1);
      r.set(m);
      for (std::size_t j = 0; j < m; ++j) if (mask >> j & 1) r.set(j);
      if (!row_is_closed(rows, r)) continue;
      rows.push_back(r);
      self(self);
      rows.pop_back();
    }
  };
  if (n) rec(rec);
}

}  // namespace riordan
