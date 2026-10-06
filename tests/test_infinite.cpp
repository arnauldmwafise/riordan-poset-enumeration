// tests/test_infinite.cpp -- decision procedure for infinite Riordan posets and the new family.
// Reference: Cheon, Curtis, Kwon, Mesinga Mwafise, Linear Algebra Appl. 632 (2022) 308-331.
#include "test_common.hpp"

TEST(phase7_automaton_agrees_with_matrix) {
    auto rnd = [](int deg) { return Poly(rng()) & ((Poly{1} << (deg + 1)) - 1); };
    for (int it = 0; it < 400; ++it) {
        RationalGF2 g{rnd(4) | 1, rnd(4) | 1}, f{(rnd(4) & ~Poly{1}) | 2, rnd(4) | 1};
        RiordanAutomaton A(g, f);
        const std::size_t n = 96;
        auto M = riordan_matrix(rational_series(g, n), rational_series(f, n), n);
        bool same = true;
        for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j) same &= M(i, j) == A.entry(i, j);
        CHECK(same);
        auto r = A.decide();
        if (r.infinite) {
            CHECK_EQ(poset_horizon(rational_series(g, 1024), rational_series(f, 1024), 1024), std::size_t(1024));
        } else {
            const std::size_t m = r.i + 1;
            auto W = riordan_matrix(rational_series(g, m), rational_series(f, m), m);
            CHECK(W(r.i, r.j) && W(r.j, r.k) && !W(r.i, r.k));
        }
    }
}

TEST(phase7_infinite_families_proved) {
    auto inf = [](const char* g, const char* f) { return decide_infinite(parse_rational(g), parse_rational(f)).infinite; };
    CHECK(inf("1/(1+t)", "t/(1+t)"));                    // Theorem 3.7
    for (int k = 1; k <= 12; ++k)                        // Corollary 3.8
        CHECK(inf(("1/(1+t)^" + std::to_string(k)).c_str(), "t/(1+t)"));
    CHECK(inf("1/(1+t)", "t") && inf("1/(1+t^2)", "t") && inf("1", "t"));
    CHECK(inf("(1+t+t^3)/(1+t)", "t"));                  // Gamma({3,4,5})
    CHECK(!inf("1+t^3+t^5", "t"));
    CHECK(!inf("1", "t/(1+t^2)"));
}

TEST(phase7_new_family_dual_pascal_plus_pascal) {
    // New: {(1+t)^k, t/(1+t)} = P*_{k+1} + P (block diagonal) for every k >= 0.
    for (std::size_t k = 0; k <= 12; ++k) {
        const std::size_t N = 48;
        auto g = Series::parse("(1+t)^" + std::to_string(k), N);
        auto M = riordan_matrix(g, Series::parse("t/(1+t)", N), N);
        CHECK(M == direct_sum(pascal_poset(k + 1).flip_transpose(), pascal_poset(N - k - 1)));
        CHECK(decide_infinite(parse_rational("(1+t)^" + std::to_string(k)), parse_rational("t/(1+t)")).infinite);
    }
}
