// riordan/infinite.hpp -- Phase 7: deciding whether {g,f} is an INFINITE Riordan poset
// matrix when g and f are rational functions over GF(2).
//
// Theory. Let g = a/b, f = c/d with a,b,c,d in GF(2)[t], b(0) = d(0) = 1. The bivariate
// generating function of the whole matrix is
//     G(t,u) = sum_{i,j} a(i,j) t^i u^j = g / (1 - u f) = P/Q,
//     P = a d,   Q = b d + u b c          (signs vanish mod 2).
// For r,s in {0,1} let Lambda_{r,s} F = sum_{i,j} [t^{2i+r} u^{2j+s}]F t^i u^j. In
// characteristic 2, Q^2 = Q(t^2,u^2), hence
//     Lambda_{r,s}(R/Q) = Lambda_{r,s}(R Q) / Q.
// So every section of P/Q is R/Q with R a polynomial; deg_t R <= max(deg P, deg Q) and
// deg_u R <= 1 are preserved, so finitely many R occur. Reading the binary digits of
// (i,j), least significant first, from state P and outputting R(0,0) (= the constant
// term of R/Q, as Q(0,0) = 1) computes a(i,j): the array is 2-automatic (Christol/Salon).
//
// Transitivity "for all i,j,k: a(i,j) & a(j,k) -> a(i,k)" is decided by running three
// copies of the automaton on the digit triples of (i,j,k) -- copies reading (i,j), (j,k),
// (i,k) -- and searching the finite product graph for a state with outputs (1,1,0).
// None reachable: {g,f}_n is a poset matrix for EVERY n (a proof). Otherwise the path
// spells out an explicit violating triple (i,j,k).
#pragma once
#include "gf2_series.hpp"

#include <array>
#include <cstdint>
#include <optional>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace riordan {

using Poly = std::uint64_t;  // GF(2)[t], bit i = coefficient of t^i (degree <= 63)

inline int poly_degree(Poly p) { return p ? 63 - __builtin_clzll(p) : -1; }

/// Carry-less product. The degree limits enforced by RiordanAutomaton (deg P <= 31,
/// deg Q <= 30, numerators stay of degree <= max(deg P, deg Q)) keep every product
/// below degree 62, so 64-bit words suffice; overflow is still checked.
inline Poly poly_mul(Poly a, Poly b) {
    if (a && b && poly_degree(a) + poly_degree(b) > 63) throw std::overflow_error("poly_mul: degree exceeds 63");
    Poly r = 0;
    for (int i = 0; a; ++i, a >>= 1)
        if (a & 1) r ^= b << i;
    return r;
}
/// Lambda_r on a univariate polynomial: sum_i p_{2i+r} t^i.
inline Poly poly_section(Poly p, int r) {
    Poly out = 0;
    p >>= r;
    for (int i = 0; p; ++i, p >>= 2)
        if (p & 1) out |= Poly{1} << i;
    return out;
}

struct RationalGF2 {
    Poly num = 1, den = 1;  // num / den, den(0) = 1
};

/// Parse "P/Q" or "P" where P, Q are polynomial expressions accepted by Series::parse.
inline RationalGF2 parse_rational(const std::string& s) {
    int depth = 0;
    std::size_t cut = std::string::npos;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '(') ++depth;
        else if (s[i] == ')') --depth;
        else if (s[i] == '/' && depth == 0) { cut = i; break; }
    }
    auto poly = [](const std::string& e) {
        const Series x = Series::parse(e, 64);
        Poly p = 0;
        for (std::size_t i = 0; i < 64; ++i)
            if (x[i]) p |= Poly{1} << i;
        // Guard against a power series that is not a polynomial of degree < 48.
        for (std::size_t i = 48; i < 64; ++i)
            if (x[i]) throw std::invalid_argument("parse_rational: '" + e + "' is not a polynomial of degree < 48");
        return p;
    };
    if (cut == std::string::npos) return {poly(s), 1};
    return {poly(s.substr(0, cut)), poly(s.substr(cut + 1))};
}

/// Power series of a rational function to precision p (for cross-checks).
inline Series rational_series(const RationalGF2& r, std::size_t p) {
    Series n(p), d(p);
    for (std::size_t i = 0; i < 64 && i < p; ++i) {
        if (r.num >> i & 1) n.set(i);
        if (r.den >> i & 1) d.set(i);
    }
    return n * d.inverse();
}

struct InfinitenessResult {
    bool infinite = false;                  // true: poset matrix for every order (proved)
    std::size_t i = 0, j = 0, k = 0;        // otherwise a violating triple a(i,j)=a(j,k)=1, a(i,k)=0
    std::size_t automaton_states = 0;       // size of the 2-automaton for a(i,j)
    std::size_t product_states = 0;         // product states explored
};

class RiordanAutomaton {
public:
    RiordanAutomaton(const RationalGF2& g, const RationalGF2& f) {
        if (!(g.num & 1) || !(g.den & 1) || !(f.den & 1))
            throw std::invalid_argument("RiordanAutomaton: need g(0) = 1 and denominators with constant term 1");
        if ((f.num & 1) || !(f.num & 2))
            throw std::invalid_argument("RiordanAutomaton: need f = t + O(t^2)");
        Q0_ = poly_mul(g.den, f.den);
        Q1_ = poly_mul(g.den, f.num);
        if (poly_degree(Q0_) > 30 || poly_degree(Q1_) > 30)
            throw std::invalid_argument("RiordanAutomaton: degrees too large (deg b + deg d <= 30)");
        const Poly P = poly_mul(g.num, f.den);
        if (poly_degree(P) > 31) throw std::invalid_argument("RiordanAutomaton: numerator degree too large");
        build({P, 0});
    }
    std::size_t size() const { return out_.size(); }
    bool output(std::size_t s) const { return out_[s]; }
    std::size_t next(std::size_t s, int r, int c) const { return delta_[s][2 * r + c]; }
    /// a(i,j) computed by running the automaton (for testing).
    bool entry(std::size_t i, std::size_t j) const {
        std::size_t s = 0;
        while (i || j) { s = next(s, static_cast<int>(i & 1), static_cast<int>(j & 1)); i >>= 1; j >>= 1; }
        return out_[s];
    }

    InfinitenessResult decide() const {
        InfinitenessResult res;
        res.automaton_states = size();
        const std::size_t S = size();
        auto key = [S](std::size_t a, std::size_t b, std::size_t c) { return (a * S + b) * S + c; };
        struct Parent { std::uint64_t prev; int digit; };
        std::unordered_map<std::uint64_t, Parent> seen;
        std::queue<std::uint64_t> q;
        const std::uint64_t start = key(0, 0, 0);
        seen[start] = {start, -1};
        q.push(start);
        while (!q.empty()) {
            const std::uint64_t cur = q.front();
            q.pop();
            const std::size_t A = cur / (S * S), B = (cur / S) % S, C = cur % S;
            if (out_[A] && out_[B] && !out_[C]) {  // violation reachable: reconstruct digits
                std::vector<int> digits;
                for (std::uint64_t x = cur; seen[x].digit >= 0; x = seen[x].prev) digits.push_back(seen[x].digit);
                std::size_t i = 0, j = 0, k = 0;
                for (std::size_t t = digits.size(); t-- > 0;) {  // digits were collected last-first
                    const std::size_t pos = digits.size() - 1 - t;
                    const int d = digits[t];
                    i |= std::size_t(d >> 2 & 1) << pos;
                    j |= std::size_t(d >> 1 & 1) << pos;
                    k |= std::size_t(d & 1) << pos;
                }
                res.i = i; res.j = j; res.k = k;
                res.product_states = seen.size();
                return res;
            }
            for (int d = 0; d < 8; ++d) {
                const int x = d >> 2 & 1, y = d >> 1 & 1, z = d & 1;
                const std::uint64_t nk = key(next(A, x, y), next(B, y, z), next(C, x, z));
                if (seen.emplace(nk, Parent{cur, d}).second) q.push(nk);
            }
        }
        res.infinite = true;
        res.product_states = seen.size();
        return res;
    }

private:
    using State = std::pair<Poly, Poly>;  // R = R0 + u R1
    struct Hash {
        std::size_t operator()(const State& s) const { return std::hash<Poly>()(s.first * 0x9E3779B97F4A7C15ULL ^ s.second); }
    };
    State step(const State& R, int r, int c) const {
        const auto A0 = poly_mul(R.first, Q0_);
        const auto A1 = poly_mul(R.first, Q1_) ^ poly_mul(R.second, Q0_);
        const auto A2 = poly_mul(R.second, Q1_);
        if (c == 0) return {poly_section(A0, r), poly_section(A2, r)};
        return {poly_section(A1, r), 0};
    }
    void build(const State& init) {
        std::unordered_map<State, std::size_t, Hash> id;
        std::vector<State> states{init};
        id[init] = 0;
        for (std::size_t s = 0; s < states.size(); ++s) {
            std::array<std::size_t, 4> tr{};
            for (int r = 0; r < 2; ++r)
                for (int c = 0; c < 2; ++c) {
                    const State nx = step(states[s], r, c);
                    auto [it, fresh] = id.emplace(nx, states.size());
                    if (fresh) states.push_back(nx);
                    tr[2 * r + c] = it->second;
                }
            delta_.push_back(tr);
        }
        for (const auto& s : states) out_.push_back(s.first & 1);
    }
    Poly Q0_ = 1, Q1_ = 0;
    std::vector<std::array<std::size_t, 4>> delta_;
    std::vector<bool> out_;
};

/// Decide whether {g,f} is an infinite Riordan poset matrix (g, f rational over GF(2)).
inline InfinitenessResult decide_infinite(const RationalGF2& g, const RationalGF2& f) {
    return RiordanAutomaton(g, f).decide();
}

}  // namespace riordan
