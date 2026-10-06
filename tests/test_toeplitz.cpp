// tests/test_toeplitz.cpp -- Toeplitz posets and numerical semigroups (Section 4).
// Reference: Cheon, Curtis, Kwon, Mesinga Mwafise, Linear Algebra Appl. 632 (2022) 308-331.
#include "test_common.hpp"

TEST(phase5_toeplitz_counts_and_theorem_4_3) {
    const std::uint64_t paper[] = {0, 1, 2, 3, 5, 7, 12, 16, 27, 37, 58, 80, 131, 171, 277, 380, 580, 785};
    for (std::size_t n = 1; n <= 17; ++n) CHECK_EQ(count_toeplitz(n), paper[n]);
    for (std::size_t n = 1; n <= 12; ++n) {   // Theorem 4.3 holds up to isomorphism; Corollary 4.4
        std::set<std::string> iso;
        std::uint64_t k = 0;
        for_each_toeplitz_support(n, [&](const Bits& s) {
            ++k;
            auto m = riordan_matrix(Series(s), Series::t(n), n);
            CHECK(m.is_poset_matrix());
            auto cf = canonical_form(m);
            CHECK(canonical_form(dual(m)) == cf);
            iso.insert(cf.matrix().compact());
        });
        CHECK_EQ(iso.size(), k);
    }
}

TEST(phase5_corollary_4_2_exhaustive) {
    const std::size_t n = 12;
    for (std::uint32_t mask = 0; mask < (1u << (n - 1)); ++mask) {
        Bits b(n); b.set(0);
        for (std::size_t i = 1; i < n; ++i) if (mask >> (i - 1) & 1) b.set(i);
        CHECK_EQ(riordan_matrix(Series(b), Series::t(n), n).is_poset_matrix(), is_additively_closed(b));
    }
}

TEST(phase5_frobenius_and_generators) {
    CHECK_EQ(frobenius_number({3, 4}), 5LL);
    CHECK_EQ(frobenius_number({6, 9, 20}), 43LL);
    CHECK_EQ(frobenius_number({4, 6, 9}), 11LL);
    CHECK_EQ(frobenius_number({2, 4}), -2LL);
    CHECK_EQ(frobenius_number({1, 5}), -1LL);
    for (std::size_t a = 2; a < 15; ++a) for (std::size_t b = a + 1; b < 20; ++b) {
        std::size_t g = a, h = b; while (h) { auto r = g % h; g = h; h = r; }
        if (g == 1) CHECK_EQ(frobenius_number({a, b}), (long long)(a * b - a - b));   // Sylvester
    }
    for (int it = 0; it < 200; ++it) {   // against brute force on Gamma(S)
        std::vector<std::size_t> S{2 + rng() % 9, 2 + rng() % 13, 2 + rng() % 17};
        const Bits gam = monoid_closure(S, 400);
        long long brute = -2;
        std::size_t g = 0; for (auto s : S) { std::size_t a = g, b = s; while (b) { auto r = a % b; a = b; b = r; } g = a; }
        if (g == 1) { brute = -1; for (std::size_t v = 0; v < 400; ++v) if (!gam.test(v)) brute = (long long)v; }
        CHECK_EQ(frobenius_number(S), brute);
        CHECK(is_additively_closed(gam));
        CHECK(monoid_closure(minimal_generators(gam), 400) == gam);
    }
    auto p = toeplitz_pair({3, 4}, 300);   // example after Theorem 4.1
    CHECK_EQ(poset_horizon(p.g, p.f, 300), std::size_t(300));
}
