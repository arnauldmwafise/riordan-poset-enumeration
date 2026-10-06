// riordan/exponential.hpp -- Phase 6: binary exponential Riordan matrices [g,f] (Section 5).
//
// <g,f> = [b_{n,k}] has k-th column with exponential generating function g f^k / k!.
// Write g = sum G_n t^n/n!, f = sum F_n t^n/n!. When the EGF coefficients G_n, F_n are
// integers (true for every family in the paper: Stirling, Laguerre, Hermite, Jacobi
// cn/sn), every b_{n,k} is an integer given by integer recurrences:
//
//   partial Bell polynomials  B_{0,0} = 1,
//       B_{m,k} = sum_{i=1}^{m-k+1} C(m-1,i-1) F_i B_{m-i,k-1}         (= m! [t^m] f^k/k!)
//   matrix entries            b_{n,k} = sum_{m=k}^{n} C(n,m) G_{n-m} B_{m,k}.
//
// Everything is a ring operation in Z, so it can be evaluated exactly in Z/2: binomial
// coefficients mod 2 by Lucas' theorem, G and F as bitsets. No rationals, no big integers.
#pragma once
#include "bool_matrix.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace riordan {

/// EGF coefficients mod 2: bit n of G is G_n mod 2, bit n of F is F_n mod 2.
struct EgfPair {
    Bits G, F;
};

inline bool binom_mod2(std::size_t n, std::size_t k) { return k <= n && (k & ~n) == 0; }

/// The binary exponential Riordan matrix [g,f]_n = <g,f>_n mod 2.
inline BoolMatrix exponential_riordan_mod2(const EgfPair& p, std::size_t n) {
    if (p.G.size() < n || p.F.size() < n) throw std::invalid_argument("exponential_riordan_mod2: too few coefficients");
    if (n >= 1 && !p.G.test(0)) throw std::invalid_argument("exponential_riordan_mod2: need G_0 odd");
    if (n >= 1 && p.F.test(0)) throw std::invalid_argument("exponential_riordan_mod2: need F_0 even (f(0) = 0)");
    if (n >= 2 && !p.F.test(1)) throw std::invalid_argument("exponential_riordan_mod2: need F_1 odd");
    // bell[m] = bitset over k of B_{m,k} mod 2
    std::vector<Bits> bell(n, Bits(n));
    if (n) bell[0].set(0);
    for (std::size_t m = 1; m < n; ++m)
        for (std::size_t k = 1; k <= m; ++k) {
            bool s = false;
            for (std::size_t i = 1; i <= m - k + 1; ++i)
                if (p.F.test(i) && binom_mod2(m - 1, i - 1) && bell[m - i].test(k - 1)) s = !s;
            if (s) bell[m].set(k);
        }
    BoolMatrix B(n);
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t k = 0; k <= r; ++k) {
            bool s = false;
            for (std::size_t m = k; m <= r; ++m)
                if (p.G.test(r - m) && binom_mod2(r, m) && bell[m].test(k)) s = !s;
            if (s) B.set(r, k);
        }
    return B;
}

// ------------------------------------------------------------------ families of Section 5

/// Ordinary integer coefficients g_n, f_n: G_n = n! g_n, so G_n, F_n are even for n >= 2
/// (Theorem 5.2's setting). Only g_0, g_1, f_1 matter mod 2.
inline EgfPair egf_from_integer_ogf(bool g1, std::size_t p) {
    EgfPair e{Bits(p), Bits(p)};
    if (p) e.G.set(0);
    if (p > 1) { e.G.set(1, g1); e.F.set(1); }
    return e;
}

/// Stirling numbers of the second kind: <e^t, e^t - 1>.
inline EgfPair stirling2_egf(std::size_t p) {
    EgfPair e{Bits(p), Bits(p)};
    for (std::size_t i = 0; i < p; ++i) { e.G.set(i); if (i) e.F.set(i); }
    return e;
}

/// Generalized Laguerre: <(1-t)^{-alpha-1}, t/(t-1)>.
/// G_n = (alpha+1)(alpha+2)...(alpha+n) (rising factorial), F_n = -n!.
inline EgfPair laguerre_egf(long long alpha, std::size_t p) {
    EgfPair e{Bits(p), Bits(p)};
    bool prod = true;  // empty product
    for (std::size_t i = 0; i < p; ++i) {
        if (i) prod = prod && (((alpha + static_cast<long long>(i)) % 2 + 2) % 2 == 1);
        e.G.set(i, prod);
    }
    if (p > 1) e.F.set(1);  // n! is even for n >= 2
    return e;
}

/// Generalized Hermite: <e^{-v t^2/2}, t>. G_{2m} = (-v)^m (2m-1)!!, odd parts zero;
/// (2m-1)!! is odd, so G_{2m} = v^m mod 2.
inline EgfPair hermite_egf(long long v, std::size_t p) {
    EgfPair e{Bits(p), Bits(p)};
    const bool vodd = (v % 2 + 2) % 2 == 1;
    for (std::size_t i = 0; i < p; i += 2) e.G.set(i, i == 0 || vodd);
    if (p > 1) e.F.set(1);
    return e;
}

/// Jacobi elliptic functions sn, cn, dn with modulus m, EGF coefficients mod 2, from
///   sn' = cn dn,  cn' = -sn dn,  dn' = -m sn cn,   sn(0) = 0, cn(0) = dn(0) = 1
/// (an EGF derivative is a shift and an EGF product a binomial convolution).
struct JacobiEgf { Bits sn, cn, dn; };
inline JacobiEgf jacobi_egf(long long m, std::size_t p) {
    JacobiEgf J{Bits(p), Bits(p), Bits(p)};
    const bool modd = (m % 2 + 2) % 2 == 1;
    if (p) { J.cn.set(0); J.dn.set(0); }
    auto conv = [&](const Bits& a, const Bits& b, std::size_t n) {
        bool s = false;
        for (std::size_t i = 0; i <= n; ++i)
            if (binom_mod2(n, i) && a.test(i) && b.test(n - i)) s = !s;
        return s;
    };
    for (std::size_t n = 0; n + 1 < p; ++n) {
        J.sn.set(n + 1, conv(J.cn, J.dn, n));
        J.cn.set(n + 1, conv(J.sn, J.dn, n));
        J.dn.set(n + 1, modd && conv(J.sn, J.cn, n));
    }
    return J;
}
/// [cn(t,m), sn(t,m)].
inline EgfPair jacobi_cn_sn_egf(long long m, std::size_t p) {
    const auto J = jacobi_egf(m, p);
    return {J.cn, J.sn};
}

}  // namespace riordan
