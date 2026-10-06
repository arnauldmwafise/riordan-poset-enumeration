// riordan/toeplitz.hpp -- Phase 5: Toeplitz posets {g,t} and numerical semigroups
// (Section 4 of Cheon-Curtis-Kwon-Mesinga Mwafise), plus certified infinite families.
//
// {g,t} is a poset matrix iff supp(g) is closed under addition (Corollary 4.2); for an
// infinite matrix that means supp(g) is a submonoid Gamma(S) of N (Theorem 4.1). For a
// truncation {g,t}_n the condition is closure for sums < n.
#pragma once
#include "binary_riordan.hpp"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <numeric>
#include <vector>

namespace riordan {

/// Gamma(S) = { sum x_i s_i : x_i in N } restricted to [0, n), as a bitset.
inline Bits monoid_closure(const std::vector<std::size_t>& gens, std::size_t n) {
    Bits b(n);
    if (n == 0) return b;
    b.set(0);
    for (std::size_t v = 1; v < n; ++v)
        for (auto s : gens)
            if (s && s <= v && b.test(v - s)) { b.set(v); break; }
    return b;
}

/// True if the support of the bitset (containing 0) is closed under sums < size().
inline bool is_additively_closed(const Bits& supp) {
    const std::size_t n = supp.size();
    if (n == 0 || !supp.test(0)) return false;
    for (std::size_t a = 1; a < n; ++a) {
        if (!supp.test(a)) continue;
        for (std::size_t b = a; a + b < n; ++b)
            if (supp.test(b) && !supp.test(a + b)) return false;
    }
    return true;
}

/// Minimal generating set of a closed support: elements s > 0 of supp that are not a
/// sum of two smaller nonzero elements of supp (this is the minimal spanning set S).
inline std::vector<std::size_t> minimal_generators(const Bits& supp) {
    std::vector<std::size_t> gens;
    for (std::size_t s = 1; s < supp.size(); ++s) {
        if (!supp.test(s)) continue;
        bool decomposable = false;
        for (std::size_t a = 1; 2 * a <= s && !decomposable; ++a)
            decomposable = supp.test(a) && supp.test(s - a);
        if (!decomposable) gens.push_back(s);
    }
    return gens;
}

/// Frobenius number g(a_1..a_k): the largest integer not in Gamma(S). Returns -1 when
/// Gamma(S) is cofinite with no gaps (S contains 1), and -2 when gcd(S) > 1 (infinitely
/// many gaps). Exact, via shortest paths on residues mod the smallest generator
/// (Nijenhuis' algorithm).
inline long long frobenius_number(std::vector<std::size_t> gens) {
    gens.erase(std::remove(gens.begin(), gens.end(), std::size_t{0}), gens.end());
    if (gens.empty()) return -2;
    std::size_t g = 0;
    for (auto s : gens) g = std::gcd(g, s);
    if (g != 1) return -2;
    const std::size_t a = *std::min_element(gens.begin(), gens.end());
    if (a == 1) return -1;
    const std::uint64_t inf = ~std::uint64_t{0};
    std::vector<std::uint64_t> dist(a, inf);
    std::vector<char> done(a, 0);
    dist[0] = 0;
    for (std::size_t it = 0; it < a; ++it) {  // Dijkstra with a dense scan (a is small)
        std::size_t u = a;
        for (std::size_t r = 0; r < a; ++r)
            if (!done[r] && dist[r] != inf && (u == a || dist[r] < dist[u])) u = r;
        if (u == a) break;
        done[u] = 1;
        for (auto s : gens) {
            const std::size_t v = (u + s) % a;
            if (dist[u] + s < dist[v]) dist[v] = dist[u] + s;
        }
    }
    return static_cast<long long>(*std::max_element(dist.begin(), dist.end())) - static_cast<long long>(a);
}

/// The Toeplitz pair (indicator series of Gamma(S), t), to precision p.
inline RiordanPair toeplitz_pair(const std::vector<std::size_t>& gens, std::size_t p) {
    return {Series(monoid_closure(gens, p)), Series::t(p)};
}

/// Number of Toeplitz poset matrices {g,t}_n = number of subsets of {1..n-1} closed under
/// sums < n (Theorem 4.3: these correspond to minimal spanning sets, |S(n)|).
inline std::uint64_t count_toeplitz(std::size_t n) {
    if (n <= 1) return 1;
    std::uint64_t count = 0;
    std::vector<char> in(n, 0);
    in[0] = 1;
    // Decide k = 1..n-1 in order; adding k is allowed iff no a + b = k with a, b chosen is
    // left out and ... we simply branch and check closure for the newly decided element.
    std::function<void(std::size_t)> rec = [&](std::size_t k) {
        if (k == n) { ++count; return; }
        bool forced = false;  // k = a + b with a, b in the set forces k in
        for (std::size_t a = 1; 2 * a <= k && !forced; ++a) forced = in[a] && in[k - a];
        in[k] = 1;
        rec(k + 1);
        in[k] = 0;
        if (!forced) rec(k + 1);
    };
    rec(1);
    return count;
}

/// Enumerate the supports (as bitsets of length n) of all Toeplitz poset matrices of order n.
template <class Visit>
void for_each_toeplitz_support(std::size_t n, Visit&& visit) {
    Bits in(n);
    if (n == 0) return;
    in.set(0);
    std::function<void(std::size_t)> rec = [&](std::size_t k) {
        if (k == n) { visit(static_cast<const Bits&>(in)); return; }
        bool forced = false;
        for (std::size_t a = 1; 2 * a <= k && !forced; ++a) forced = in.test(a) && in.test(k - a);
        in.set(k);
        rec(k + 1);
        in.reset(k);
        if (!forced) rec(k + 1);
    };
    rec(1);
}

// ---------------------------------------------------------------- certified infinite families
// Each pair below is an infinite Riordan poset matrix by a theorem of the paper; the
// precision p only controls how many coefficients are materialised.

/// Corollary 3.8: {1/(1-t)^k, t/(1-t)} (k = 1 is the Pascal poset, Theorem 3.7).
inline RiordanPair pascal_family_pair(std::size_t k, std::size_t p) {
    Series g = Series::geometric(1, p).pow(k);
    Series f = Series::t(p) * Series::geometric(1, p);
    return {g.truncated(p), f.truncated(p)};
}
/// Theorem 4.1: {indicator(Gamma(S)), t}.
inline RiordanPair toeplitz_family_pair(const std::vector<std::size_t>& gens, std::size_t p) {
    return toeplitz_pair(gens, p);
}

}  // namespace riordan
