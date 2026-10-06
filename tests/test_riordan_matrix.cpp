// tests/test_riordan_matrix.cpp -- binary Riordan matrices, the tree construction and basic poset facts.
// Reference: Cheon, Curtis, Kwon, Mesinga Mwafise, Linear Algebra Appl. 632 (2022) 308-331.
#include "test_common.hpp"

TEST(builder_matches_definition) {   // facts F1 and F2
    for (int it = 0; it < 300; ++it) {
        std::size_t n = 1 + rng() % 70;
        Series g = random_g(n + 1), f = random_f(n + 2);
        BoolMatrix a = riordan_matrix(g, f, n);
        CHECK(a == naive_riordan(g, f, n));
        RiordanBuilder b(n);
        while (!b.full()) { std::size_t m = b.order(); b.push_coeffs(g[m], m >= 2 && f[m]); }
        CHECK(b.matrix() == a);
        auto p = recover_pair(a);
        CHECK(p.has_value());
        if (p) {
            CHECK(p->g == g.truncated(n));
            if (n >= 2) CHECK(p->f == f.truncated(n));
        }
    }
}

TEST(bijection_count_2_pow_2n_minus_3) {   // fact F1
    for (std::size_t n = 2; n <= 7; ++n) {
        std::set<std::string> distinct;
        std::size_t total = 0;
        for (std::uint64_t gb = 0; gb < (1ull << (n - 1)); ++gb)
            for (std::uint64_t fb = 0; fb < (1ull << (n - 2)); ++fb) {
                Series g(n), f(n);
                g.set(0); f.set(1);
                for (std::size_t k = 1; k < n; ++k) if (gb >> (k - 1) & 1) g.set(k);
                for (std::size_t k = 2; k < n; ++k) if (fb >> (k - 2) & 1) f.set(k);
                distinct.insert(riordan_matrix(g, f, n).compact());
                ++total;
            }
        CHECK_EQ(distinct.size(), total);
        CHECK_EQ(std::uint64_t(total), count_binary_riordan(n));
    }
}

TEST(theorem_2_1_poset_iff_idempotent) {
    for (int it = 0; it < 400; ++it) {
        std::size_t n = 1 + rng() % 14;
        BoolMatrix a = riordan_matrix(random_g(n), random_f(std::max<std::size_t>(n, 2)), n);
        CHECK_EQ(a.is_poset_matrix(), naive_bool_square(a) == a);
        CHECK(bool_product(a, a) == naive_bool_square(a));
    }
}

TEST(example_2_3_pascal_P5) {
    BoolMatrix p5 = riordan_matrix(S("1/(1-t)", 5), S("t/(1-t)", 5), 5);
    CHECK(p5 == BoolMatrix::from_rows({"10000", "11000", "10100", "11110", "10001"}));
    CHECK(p5.is_poset_matrix());
    auto sp = series_parallel_decomposition(p5);
    CHECK(sp.has_value());
}

TEST(theorem_3_7_pascal_lucas_and_horizon) {
    const std::size_t N = 512;
    BoolMatrix p = riordan_matrix(S("1/(1-t)", N), S("t/(1-t)", N), N);
    bool lucas = true;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j <= i; ++j) if (p.get(i, j) != ((i & j) == j)) lucas = false;
    CHECK(lucas);
    CHECK(p.is_poset_matrix());
    // Corollary 3.8: {1/(1-t)^k, t/(1-t)} are Riordan poset matrices.
    for (int k = 1; k <= 6; ++k) {
        std::string g = "1/(1-t)^" + std::to_string(k);
        CHECK_EQ(poset_horizon(S(g.c_str(), 600), S("t/(1-t)", 600), 600), std::size_t(600));
    }
    // P_{2^m} is the Boolean lattice: 2^m elements, height m+1, width C(m, m/2).
    BoolMatrix b3 = p.leading(8);
    CHECK_EQ(height(b3), std::size_t(4));
    CHECK_EQ(width(b3), std::size_t(3));
    CHECK_EQ(linear_extension_count(b3), std::uint64_t(48));
}

TEST(chains_antichains_and_two_chains) {
    for (std::size_t n = 1; n <= 12; ++n) {
        BoolMatrix in = riordan_matrix(S("1", n), S("t", std::max<std::size_t>(n, 2)), n);
        BoolMatrix cn = riordan_matrix(S("1/(1-t)", n), S("t", std::max<std::size_t>(n, 2)), n);
        BoolMatrix cc = riordan_matrix(S("1/(1-t^2)", n), S("t", std::max<std::size_t>(n, 2)), n);
        CHECK(in == antichain(n));
        CHECK(cn == chain(n));
        CHECK(cc.is_poset_matrix());
        CHECK(isomorphic(cc, direct_sum(chain((n + 1) / 2), chain(n / 2))));
        CHECK(is_series_parallel(cc));
    }
}
