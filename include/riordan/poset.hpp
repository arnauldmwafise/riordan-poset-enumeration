// riordan/poset.hpp -- order-theoretic toolkit on poset matrices.
//
// Every function takes an order-relation matrix R with R(i,j) = 1 iff x_j <= x_i
// (reflexive, antisymmetric, transitive).  The labelling need not be a linear
// extension unless a function says so; Riordan poset matrices always are.
#pragma once
#include "bool_matrix.hpp"

#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace riordan {

// ------------------------------------------------------------ constructors --
inline BoolMatrix antichain(std::size_t n) { return BoolMatrix::identity(n); }
inline BoolMatrix chain(std::size_t n) {
    BoolMatrix m(n);
    for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j <= i; ++j) m.set(i, j);
    return m;
}
// Equation (3): direct sum P + Q (disjoint union) and ordinal sum P ⊕ Q (Q on top).
inline BoolMatrix direct_sum(const BoolMatrix& a, const BoolMatrix& b) {
    const std::size_t p = a.size(), n = p + b.size();
    BoolMatrix m(n);
    for (std::size_t i = 0; i < p; ++i) a.row(i).for_each([&](std::size_t j) { m.set(i, j); });
    for (std::size_t i = 0; i < b.size(); ++i) b.row(i).for_each([&](std::size_t j) { m.set(p + i, p + j); });
    return m;
}
inline BoolMatrix ordinal_sum(const BoolMatrix& a, const BoolMatrix& b) {
    BoolMatrix m = direct_sum(a, b);
    for (std::size_t i = a.size(); i < m.size(); ++i)
        for (std::size_t j = 0; j < a.size(); ++j) m.set(i, j);
    return m;
}
// B_{m,n}: m minimal elements, each below all n maximal elements (= I_m ⊕ I_n).
inline BoolMatrix complete_bipartite(std::size_t m, std::size_t n) { return ordinal_sum(antichain(m), antichain(n)); }
// Dual poset: reverse the order.  For a lower triangular poset matrix use
// BoolMatrix::flip_transpose() to stay lower triangular (Section 2).
inline BoolMatrix dual(const BoolMatrix& r) { return r.transpose(); }

// --------------------------------------------------------------- basics ----
inline bool is_poset(const BoolMatrix& r) { return r.is_order_relation(); }

inline bool comparable(const BoolMatrix& r, std::size_t a, std::size_t b) { return r.get(a, b) || r.get(b, a); }

inline bool is_chain(const BoolMatrix& r) {
    for (std::size_t i = 0; i < r.size(); ++i)
        for (std::size_t j = 0; j < i; ++j) if (!comparable(r, i, j)) return false;
    return true;
}
inline bool is_antichain(const BoolMatrix& r) {
    for (std::size_t i = 0; i < r.size(); ++i) if (r.row(i).count() != 1) return false;
    return true;
}

// A linear extension: elements sorted by size of their down-set (x < y implies
// |down(x)| < |down(y)|).  order[k] = element placed at label k.
inline std::vector<std::size_t> linear_extension(const BoolMatrix& r) {
    std::vector<std::size_t> order(r.size());
    for (std::size_t v = 0; v < r.size(); ++v) order[v] = v;
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return r.row(a).count() < r.row(b).count(); });
    return order;
}
// Relabel along a linear extension so that the matrix is unit lower triangular.
inline BoolMatrix naturally_labelled(const BoolMatrix& r, std::vector<std::size_t>* order_out = nullptr) {
    auto order = linear_extension(r);
    std::vector<std::size_t> pos(r.size());
    for (std::size_t k = 0; k < order.size(); ++k) pos[order[k]] = k;
    if (order_out) *order_out = order;
    return r.permuted(pos);
}

// Cover matrix: C(i,j) = 1 iff x_i covers x_j.
inline BoolMatrix cover_matrix(const BoolMatrix& r) {
    const std::size_t n = r.size();
    BoolMatrix c(n);
    for (std::size_t i = 0; i < n; ++i) {
        Bits strict = r.row(i);
        strict.reset(i);
        Bits below2(n);   // elements strictly below something strictly below i
        strict.for_each([&](std::size_t k) {
            Bits rk = r.row(k);
            rk.reset(k);
            below2 |= rk;
        });
        strict.for_each([&](std::size_t j) { if (!below2.test(j)) c.set(i, j); });
    }
    return c;
}
// Hasse diagram edges (lower, upper).
inline std::vector<std::pair<std::size_t, std::size_t>> hasse_edges(const BoolMatrix& r) {
    std::vector<std::pair<std::size_t, std::size_t>> e;
    BoolMatrix c = cover_matrix(r);
    for (std::size_t i = 0; i < r.size(); ++i) c.row(i).for_each([&](std::size_t j) { e.emplace_back(j, i); });
    return e;
}

// ------------------------------------------------- incidence algebra (Thm 2.7) --
// mu[i][j] = μ(x_j, x_i) (lower element j, upper element i), 0 unless x_j <= x_i.
// Equivalently the integer inverse of the poset matrix.  Overflow throws.
inline std::vector<std::vector<long long>> mobius(const BoolMatrix& r) {
    const std::size_t n = r.size();
    auto order = linear_extension(r);
    std::vector<std::vector<long long>> mu(n, std::vector<long long>(n, 0));
    for (std::size_t j = 0; j < n; ++j) {
        mu[j][j] = 1;
        for (std::size_t i : order) {
            if (i == j || !r.get(i, j)) continue;
            long long s = 0;   // sum over x_j <= z < x_i of μ(x_j, z)
            r.row(i).for_each([&](std::size_t z) {
                if (z != i && r.get(z, j))
                    if (__builtin_add_overflow(s, mu[z][j], &s)) throw std::overflow_error("mobius overflow");
            });
            mu[i][j] = -s;
        }
    }
    return mu;
}

// Elements of the matrix algebra A(g,f)_n of Theorem 2.7: lower triangular
// matrices X with X(i,j) = 0 whenever A(i,j) = 0.  T is any ring type.
template <class T>
class IncidenceAlgebra {
public:
    using Element = std::vector<std::vector<T>>;
    explicit IncidenceAlgebra(BoolMatrix pattern) : a_(std::move(pattern)) {
        if (!a_.is_poset_matrix()) throw std::invalid_argument("IncidenceAlgebra: pattern must be a poset matrix");
    }
    std::size_t size() const { return a_.size(); }
    const BoolMatrix& pattern() const { return a_; }

    Element zero() const { return Element(size(), std::vector<T>(size(), T(0))); }
    Element identity() const { Element x = zero(); for (std::size_t i = 0; i < size(); ++i) x[i][i] = T(1); return x; }
    Element zeta() const {   // the pattern itself; its transpose is ζ
        Element x = zero();
        for (std::size_t i = 0; i < size(); ++i) a_.row(i).for_each([&](std::size_t j) { x[i][j] = T(1); });
        return x;
    }
    bool contains(const Element& x) const {
        for (std::size_t i = 0; i < size(); ++i)
            for (std::size_t j = 0; j < size(); ++j)
                if (!a_.get(i, j) && x[i][j] != T(0)) return false;
        return true;
    }
    // Closed in the algebra because the pattern is transitive (convolution (1)).
    Element multiply(const Element& x, const Element& y) const {
        Element z = zero();
        for (std::size_t i = 0; i < size(); ++i)
            a_.row(i).for_each([&](std::size_t k) {
                if (x[i][k] == T(0)) return;
                a_.row(k).for_each([&](std::size_t j) { z[i][j] += x[i][k] * y[k][j]; });
            });
        return z;
    }
    // Inverse by forward substitution; requires invertible diagonal in T.
    Element inverse(const Element& x) const {
        Element y = zero();
        for (std::size_t j = 0; j < size(); ++j) {
            y[j][j] = T(1) / x[j][j];
            for (std::size_t i = j + 1; i < size(); ++i) {
                if (!a_.get(i, j)) continue;
                T s = T(0);
                a_.row(i).for_each([&](std::size_t k) { if (k >= j && k < i && a_.get(k, j)) s += x[i][k] * y[k][j]; });
                y[i][j] = -s / x[i][i];
            }
        }
        return y;
    }

private:
    BoolMatrix a_;
};

// ------------------------------------------------------------ statistics --
// c[k] = number of chains with exactly k elements (c[0] = 1).  Overflow throws.
inline std::vector<std::uint64_t> chain_counts(const BoolMatrix& r) {
    const std::size_t n = r.size();
    auto order = linear_extension(r);
    // top[v][k] = chains of k elements whose maximum is v
    std::vector<std::vector<std::uint64_t>> top(n);
    std::vector<std::uint64_t> c{1};
    for (std::size_t v : order) {
        std::vector<std::uint64_t> t{0, 1};
        r.row(v).for_each([&](std::size_t u) {
            if (u == v) return;
            if (top[u].size() + 1 > t.size()) t.resize(top[u].size() + 1, 0);
            for (std::size_t k = 1; k < top[u].size(); ++k)
                if (__builtin_add_overflow(t[k + 1], top[u][k], &t[k + 1])) throw std::overflow_error("chain_counts overflow");
        });
        top[v] = t;
        if (t.size() > c.size()) c.resize(t.size(), 0);
        for (std::size_t k = 1; k < t.size(); ++k)
            if (__builtin_add_overflow(c[k], t[k], &c[k])) throw std::overflow_error("chain_counts overflow");
    }
    while (c.size() > 1 && c.back() == 0) c.pop_back();
    return c;
}
// Number of elements in a longest chain.
inline std::size_t height(const BoolMatrix& r) { return r.size() ? chain_counts(r).size() - 1 : 0; }

// Size of a largest antichain = n - maximum matching in the comparability
// bipartite graph (Dilworth / König).
inline std::size_t width(const BoolMatrix& r) {
    const std::size_t n = r.size();
    std::vector<long> match(n, -1);
    std::function<bool(std::size_t, std::vector<char>&)> aug = [&](std::size_t u, std::vector<char>& seen) {
        bool found = false;
        r.row(u).for_each([&](std::size_t v) {
            if (found || v == u || seen[v]) return;
            seen[v] = 1;
            if (match[v] < 0 || aug(std::size_t(match[v]), seen)) { match[v] = long(u); found = true; }
        });
        return found;
    };
    std::size_t m = 0;
    for (std::size_t u = 0; u < n; ++u) {
        std::vector<char> seen(n, 0);
        if (aug(u, seen)) ++m;
    }
    return n - m;
}

// Number of linear extensions (DP over order ideals), n <= 26.
inline std::uint64_t linear_extension_count(const BoolMatrix& r) {
    const std::size_t n = r.size();
    if (n > 26) throw std::length_error("linear_extension_count: n > 26");
    std::vector<std::uint32_t> pred(n, 0);
    for (std::size_t v = 0; v < n; ++v) r.row(v).for_each([&](std::size_t u) { if (u != v) pred[v] |= 1u << u; });
    std::vector<std::uint64_t> f(std::size_t(1) << n, 0);
    f[0] = 1;
    for (std::uint32_t mask = 0; mask < (1u << n); ++mask) {
        if (!f[mask]) continue;
        for (std::size_t v = 0; v < n; ++v)
            if (!(mask >> v & 1) && (pred[v] & ~mask) == 0)
                if (__builtin_add_overflow(f[mask | (1u << v)], f[mask], &f[mask | (1u << v)]))
                    throw std::overflow_error("linear_extension_count overflow");
    }
    return f[(std::size_t(1) << n) - 1];
}

// ------------------------------------------------------ structural classes --
// Forest poset (Section 5): no element covers two different elements.
inline bool is_forest(const BoolMatrix& r) {
    BoolMatrix c = cover_matrix(r);
    for (std::size_t i = 0; i < r.size(); ++i) if (c.row(i).count() > 1) return false;
    return true;
}

// Contains an induced N (a<b, c<b, c<d; a||c, a||d, b||d)?  O(n^4), for testing.
inline bool contains_N(const BoolMatrix& r) {
    const std::size_t n = r.size();
    auto lt = [&](std::size_t x, std::size_t y) { return x != y && r.get(y, x); };
    auto inc = [&](std::size_t x, std::size_t y) { return x != y && !comparable(r, x, y); };
    for (std::size_t b = 0; b < n; ++b)
        for (std::size_t a = 0; a < n; ++a) if (lt(a, b))
            for (std::size_t c = 0; c < n; ++c) if (c != a && lt(c, b) && inc(a, c))
                for (std::size_t d = 0; d < n; ++d)
                    if (d != b && lt(c, d) && inc(a, d) && inc(b, d)) return true;
    return false;
}

// Series-parallel decomposition: returns an expression in "1", "+" (direct sum)
// and "⊕" (ordinal sum, lower part first), or nullopt if the poset is not
// series-parallel.  Recursion: split by connected components of the
// comparability graph (direct sum) or of the incomparability graph (ordinal sum).
inline std::optional<std::string> series_parallel_decomposition(const BoolMatrix& r, const std::string& plus = " + ",
                                                                const std::string& oplus = " ⊕ ") {
    std::function<std::optional<std::string>(const std::vector<std::size_t>&)> rec =
        [&](const std::vector<std::size_t>& s) -> std::optional<std::string> {
        if (s.size() == 1) return std::string("1");
        auto components = [&](bool comparability) {
            std::vector<long> comp(s.size(), -1);
            long nc = 0;
            for (std::size_t a = 0; a < s.size(); ++a) {
                if (comp[a] >= 0) continue;
                std::vector<std::size_t> stack{a};
                comp[a] = nc;
                while (!stack.empty()) {
                    std::size_t x = stack.back(); stack.pop_back();
                    for (std::size_t y = 0; y < s.size(); ++y)
                        if (comp[y] < 0 && y != x && comparable(r, s[x], s[y]) == comparability) { comp[y] = nc; stack.push_back(y); }
                }
                ++nc;
            }
            std::vector<std::vector<std::size_t>> parts(nc);
            for (std::size_t a = 0; a < s.size(); ++a) parts[comp[a]].push_back(s[a]);
            return parts;
        };
        auto wrap = [](const std::string& e) { return e == "1" ? e : "(" + e + ")"; };
        auto parts = components(true);
        std::string sep = plus;
        if (parts.size() == 1) {
            parts = components(false);
            if (parts.size() == 1) return std::nullopt;   // both graphs connected: not SP
            // Order the blocks: a block lower in the poset has elements below the others.
            std::sort(parts.begin(), parts.end(), [&](const auto& p, const auto& q) { return r.get(q[0], p[0]); });
            sep = oplus;
        }
        std::string out;
        for (std::size_t k = 0; k < parts.size(); ++k) {
            auto e = rec(parts[k]);
            if (!e) return std::nullopt;
            if (k) out += sep;
            out += wrap(*e);
        }
        return out;
    };
    if (r.size() == 0) return std::string("");
    std::vector<std::size_t> all(r.size());
    for (std::size_t v = 0; v < r.size(); ++v) all[v] = v;
    return rec(all);
}
inline bool is_series_parallel(const BoolMatrix& r) { return series_parallel_decomposition(r).has_value(); }

// Finite binomial-type test: every interval [x,y] is graded (all maximal chains
// have equal length) and the number of maximal chains depends only on the length.
inline bool is_binomial(const BoolMatrix& r) {
    const std::size_t n = r.size();
    BoolMatrix c = cover_matrix(r);
    auto order = linear_extension(r);
    std::map<std::size_t, std::uint64_t> chains_of_length;
    for (std::size_t x = 0; x < n; ++x) {
        std::vector<std::size_t> mn(n, 0), mx(n, 0);
        std::vector<std::uint64_t> cnt(n, 0);
        std::vector<char> in(n, 0);
        in[x] = 1; cnt[x] = 1;
        for (std::size_t y : order) {
            if (y == x || !r.get(y, x)) continue;
            bool first = true;
            c.row(y).for_each([&](std::size_t z) {
                if (!in[z]) return;
                if (first) { mn[y] = mn[z] + 1; mx[y] = mx[z] + 1; first = false; }
                else { mn[y] = std::min(mn[y], mn[z] + 1); mx[y] = std::max(mx[y], mx[z] + 1); }
                cnt[y] += cnt[z];
            });
            in[y] = 1;
            if (mn[y] != mx[y]) return false;
            auto [it, fresh] = chains_of_length.emplace(mn[y], cnt[y]);
            if (!fresh && it->second != cnt[y]) return false;
        }
    }
    return true;
}

}  // namespace riordan
