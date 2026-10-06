// riordan/binary_riordan.hpp -- binary Riordan matrices {g, f}_n over GF(2).
//
// {g, f} = [a_ij] with a_ij = [t^i] g f^j  (mod 2),  g = 1 + g_1 t + ...,
// f = t + f_2 t^2 + ... .  The truncation {g, f}_n is its leading n x n block.
//
// Facts used here (proved in docs/THEORY.md):
//  (F1) {g, f}_n determines g mod t^n (column 0) and f mod t^n
//       (column 1 = g f), so  (g mod t^n, f mod t^n)  <->  {g, f}_n  is a
//       bijection and there are 2^(2n-3) binary Riordan matrices of order n>=2.
//  (F2) Row n of {g, f}_{n+1}:  a_{n,j} = sum_{k<n} a_{k,j-1} f_{n-k}  for j >= 2,
//       so entries in columns 2..n are forced by rows 0..n-1 and f_2..f_{n-1};
//       a_{n,0} = g_n and a_{n,1} = f_n + sum_{k=1}^{n-1} g_k f_{n-k} are free.
//  (F3) With the new element carrying the largest label, the extended matrix is
//       transitive iff the old one is and row_j ⊆ row_n for every j with a_{n,j}=1.
// Together these turn "find all Riordan poset matrices" into a depth-first
// search in a tree with branching factor <= 4 (see enumerate.hpp).
#pragma once
#include "bool_matrix.hpp"
#include "gf2_series.hpp"

#include <array>
#include <optional>
#include <utility>

namespace riordan {

struct RiordanPair {
    Series g, f;
    std::string str() const { return "g = " + g.str() + ",  f = " + f.str(); }
    bool operator==(const RiordanPair& o) const { return g == o.g && f == o.f; }
};

inline void check_riordan_pair(const Series& g, const Series& f, std::size_t n) {
    if (n == 0) return;
    if (g.precision() < n) throw std::invalid_argument("binary Riordan: g known to fewer than n terms");
    if (!g[0]) throw std::invalid_argument("binary Riordan: g(0) must be 1");
    if (n >= 2) {
        if (f.precision() < n) throw std::invalid_argument("binary Riordan: f known to fewer than n terms");
        if (f[0] || !f[1]) throw std::invalid_argument("binary Riordan: f must be t + O(t^2)");
    }
}

// {g, f}_n, built column by column: column j = g f^j mod t^n.
inline BoolMatrix riordan_matrix(const Series& g, const Series& f, std::size_t n) {
    check_riordan_pair(g, f, n);
    BoolMatrix a(n);
    Series col = g.truncated(n);
    Series fn = n >= 2 ? f.truncated(n) : Series::t(std::max<std::size_t>(n, 2)).truncated(n);
    for (std::size_t j = 0; j < n; ++j) {
        col.bits().for_each([&](std::size_t i) { a.set(i, j); });
        if (j + 1 < n) col = Series::mul(col, fn, n);
    }
    return a;
}
inline BoolMatrix riordan_matrix(const RiordanPair& p, std::size_t n) { return riordan_matrix(p.g, p.f, n); }

// Inverse of riordan_matrix: recover (g mod t^n, f mod t^n) from a candidate
// matrix, or nullopt if the matrix is not a binary Riordan matrix.
inline std::optional<RiordanPair> recover_pair(const BoolMatrix& a) {
    const std::size_t n = a.size();
    if (n == 0 || !a.is_unit_lower_triangular()) return std::nullopt;
    Series g(a.col(0));
    Series f = Series::t(std::max<std::size_t>(n, 2)).truncated(n);
    if (n >= 2) f = Series(a.col(1)) * g.inverse();
    if (riordan_matrix(g, f, n) != a) return std::nullopt;
    return RiordanPair{g, f};
}
inline bool is_binary_riordan(const BoolMatrix& a) { return recover_pair(a).has_value(); }
inline bool is_riordan_poset_matrix(const BoolMatrix& a) { return a.is_poset_matrix() && is_binary_riordan(a); }

// Number of binary Riordan matrices of order n (fact F1).
inline std::uint64_t count_binary_riordan(std::size_t n) {
    if (n <= 1) return 1;
    if (2 * n - 3 >= 64) throw std::overflow_error("count_binary_riordan: too large");
    return std::uint64_t(1) << (2 * n - 3);
}

// ------------------------------------------------------------------ builder --
// Incremental construction of {g, f}_n one row at a time (facts F2, F3).
// Capacity is fixed at construction; push/pop are O(n^2 / 64).
class RiordanBuilder {
public:
    explicit RiordanBuilder(std::size_t capacity)
        : cap_(std::max<std::size_t>(capacity, 1)), rows_(cap_, Bits(cap_)), cols_(cap_, Bits(cap_)),
          g_(cap_), f_(std::max<std::size_t>(cap_, 2)), frev_(cap_) {
        rows_[0].set(0); cols_[0].set(0);
        g_.set(0); f_.set(1);
        n_ = 1;
    }

    std::size_t order() const { return n_; }
    std::size_t capacity() const { return cap_; }
    bool full() const { return n_ >= cap_; }
    bool get(std::size_t i, std::size_t j) const { return rows_[i].test(j); }
    const Bits& row(std::size_t i) const { return rows_[i]; }
    const Bits& col(std::size_t j) const { return cols_[j]; }

    // Forced part of the next row (index n = order()): columns 2..n-1 and the
    // diagonal.  Columns 0 and 1 are left 0.  Also returns
    // conv = sum_{k=1}^{n-1} g_k f_{n-k}, so that a_{n,1} = f_n + conv.
    std::pair<Bits, bool> forced_next_row() {
        const std::size_t n = n_;
        if (n >= cap_) throw std::length_error("RiordanBuilder: capacity reached");
        Bits r(cap_);
        r.set(n);
        frev_.clear();
        for (std::size_t k = 1; k < n; ++k) if (f_.test(n - k)) frev_.set(k);   // frev[k] = f_{n-k}
        for (std::size_t j = 2; j < n; ++j)
            if (Bits::and_parity(cols_[j - 1], frev_)) r.set(j);
        bool conv = (n >= 2) && Bits::and_parity(g_, frev_);
        return {r, conv};
    }

    // Would a row with these entries keep the matrix transitive?  (fact F3)
    bool row_is_closed(const Bits& r) const {
        bool ok = true;
        r.for_each([&](std::size_t j) {
            if (ok && j < n_ && !Bits::subset(rows_[j], r)) ok = false;
        });
        return ok;
    }

    // Append row n with a_{n,0} = a0 and a_{n,1} = a1 (a1 is ignored for n = 1,
    // where column 1 is the diagonal).  Returns the coefficients (g_n, f_n).
    std::pair<bool, bool> push_entries(bool a0, bool a1) {
        auto [r, conv] = forced_next_row();
        return push_row_impl(r, conv, a0, a1);
    }
    // Append row n using the coefficients g_n and f_n.
    void push_coeffs(bool gn, bool fn) {
        auto [r, conv] = forced_next_row();
        bool a1 = fn ^ conv;
        push_row_impl(r, conv, gn, a1);
    }

    // The valid next rows: every (a0, a1) whose row keeps transitivity.
    struct Candidate { bool a0, a1; Bits row; };
    std::vector<Candidate> candidates() {
        std::vector<Candidate> out;
        auto [r, conv] = forced_next_row();
        (void)conv;
        if (!forced_part_closed(r)) return out;   // no choice of columns 0,1 can repair it
        const int a1_choices = (n_ >= 2) ? 2 : 1;
        for (int a0 = 0; a0 < 2; ++a0)
            for (int a1 = 0; a1 < a1_choices; ++a1) {
                Bits c = r;
                c.set(0, a0);
                if (n_ >= 2) c.set(1, a1);
                if (row_is_closed(c)) out.push_back({bool(a0), bool(a1), std::move(c)});
            }
        return out;
    }

    void pop() {
        if (n_ <= 1) throw std::logic_error("RiordanBuilder::pop on order 1");
        --n_;
        rows_[n_].for_each([&](std::size_t j) { cols_[j].reset(n_); });
        rows_[n_].clear();
        g_.reset(n_);
        if (n_ >= 2) f_.reset(n_);   // f_1 = 1 is permanent
    }

    BoolMatrix matrix() const {
        BoolMatrix a(n_);
        for (std::size_t i = 0; i < n_; ++i)
            rows_[i].for_each([&](std::size_t j) { a.set(i, j); });
        return a;
    }
    // (g mod t^n, f mod t^n) of the current matrix.
    RiordanPair pair() const {
        Series g(n_), f(std::max<std::size_t>(n_, 2));
        for (std::size_t k = 0; k < n_; ++k) { if (g_.test(k)) g.set(k); if (f_.test(k)) f.set(k); }
        f.set(1);
        return {g, f};
    }
    bool g_coeff(std::size_t k) const { return g_.test(k); }
    bool f_coeff(std::size_t k) const { return f_.test(k); }

private:
    // Entries in columns >= 2 of r must already be closed among themselves,
    // since columns 0 and 1 only ever add rows 0 and 1 (which live in {0,1}).
    bool forced_part_closed(const Bits& r) const {
        bool ok = true;
        r.for_each([&](std::size_t j) {
            if (!ok || j < 2 || j >= n_) return;
            const Bits& rj = rows_[j];
            rj.for_each([&](std::size_t k) { if (k >= 2 && !r.test(k)) ok = false; });
        });
        return ok;
    }
    std::pair<bool, bool> push_row_impl(Bits r, bool conv, bool a0, bool a1) {
        const std::size_t n = n_;
        r.set(0, a0);
        if (n >= 2) r.set(1, a1);
        const bool gn = a0;
        const bool fn = (n >= 2) ? (a1 ^ conv) : false;
        rows_[n] = r;
        r.for_each([&](std::size_t j) { cols_[j].set(n); });
        g_.set(n, gn);
        if (n >= 2) f_.set(n, fn);
        ++n_;
        return {gn, fn};
    }

    std::size_t cap_;
    std::vector<Bits> rows_, cols_;
    Bits g_, f_, frev_;
    std::size_t n_ = 0;
};

// Largest m <= cap such that {g, f}_m is a poset matrix (well defined because
// leading principal submatrices of poset matrices are poset matrices).
// Returns cap when no failure is found up to cap ("horizon >= cap").
inline std::size_t poset_horizon(const Series& g, const Series& f, std::size_t cap) {
    check_riordan_pair(g, f, std::min({cap, g.precision(), f.precision()}));
    cap = std::min({cap, g.precision(), f.precision()});
    if (cap <= 1) return cap;
    RiordanBuilder b(cap);
    while (!b.full()) {
        const std::size_t n = b.order();
        auto [r, conv] = b.forced_next_row();
        r.set(0, g[n]);
        if (n >= 2) r.set(1, f[n] ^ conv);
        if (!b.row_is_closed(r)) return n;
        b.push_coeffs(g[n], n >= 2 ? f[n] : false);
    }
    return cap;
}

// ------------------------------------------------------- derived pairs --
// Dual poset (Theorem 2.8): the flip-transpose of {g, f}_n is again a binary
// Riordan matrix; we compute it from the matrix and read its pair off (F1),
// which avoids depending on the closed formula.
inline RiordanPair dual_pair(const Series& g, const Series& f, std::size_t n) {
    auto p = recover_pair(riordan_matrix(g, f, n).flip_transpose());
    if (!p) throw std::logic_error("dual_pair: flip-transpose not Riordan (contradicts Theorem 2.8)");
    return *p;
}

// Theorem 3.6: deleting the first m rows and columns of {g, f} gives
// {g f^m / t^m, f}.  Result precision = prec - m.
inline RiordanPair shifted_pair(const Series& g, const Series& f, std::size_t m) {
    std::size_t p = std::min(g.precision(), f.precision());
    if (m >= p) throw std::invalid_argument("shifted_pair: m too large for precision");
    Series gm = Series::mul(g.truncated(p), f.truncated(p).pow(m), p).shift_down(m);
    return {gm, f.truncated(p - m)};
}

// Lemma 3.1: subdiagonal a_{k+1,k} = g_1 + k f_2 (mod 2).
inline bool subdiagonal_entry(const Series& g, const Series& f, std::size_t k) {
    return g[1] ^ ((k & 1) && f[2]);
}
// Theorem 3.2: for a Riordan poset matrix, the order is total iff g_1 = 1, f_2 = 0.
inline bool is_chain_pair(const Series& g, const Series& f) { return g[1] && !f[2]; }

}  // namespace riordan
