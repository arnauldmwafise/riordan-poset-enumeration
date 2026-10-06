// tests/test_paper_sections_2_3.cpp -- statements of Sections 2-3 and 5.1 of the paper.
// Reference: Cheon, Curtis, Kwon, Mesinga Mwafise, Linear Algebra Appl. 632 (2022) 308-331.
#include "test_common.hpp"

TEST(theorem_2_4_complete_bipartite) {
    struct Case { std::size_t m, n; const char* g; const char* f; };
    for (auto c : {Case{1, 1, "1+t", "t"}, Case{1, 2, "1+t+t^2", "t+t^2"}, Case{2, 1, "1+t^2", "t+t^2"},
                   Case{2, 2, "1+t^2+t^3", "t+t^2"}}) {
        std::size_t N = c.m + c.n;
        BoolMatrix a = riordan_matrix(S(c.g, N), S(c.f, N), N);
        CHECK(a.is_poset_matrix());
        CHECK(isomorphic(a, complete_bipartite(c.m, c.n)));
    }
    for (std::size_t m = 1; m <= 4; ++m)
        for (std::size_t n = 1; n <= 4; ++n)
            CHECK_EQ(is_riordan_poset(complete_bipartite(m, n)), (m <= 2 && n <= 2));
}

TEST(theorem_2_5_and_four_element_posets) {
    // Every poset on <= 3 elements is Riordan; exactly 5 of the 16 on 4 are not.
    for (std::size_t n = 1; n <= 4; ++n) {
        std::unordered_map<CanonicalForm, BoolMatrix, CanonicalFormHash> classes;
        for_each_naturally_labelled_poset(n, [&](const BoolMatrix& m) { classes.emplace(canonical_form(m), m); });
        std::vector<BoolMatrix> non;
        for (auto& [cf, m] : classes) if (!is_riordan_poset(m)) non.push_back(m);
        if (n <= 3) CHECK(non.empty());
        if (n == 4) {
            CHECK_EQ(classes.size(), std::size_t(16));
            CHECK_EQ(non.size(), std::size_t(5));
            BoolMatrix one = chain(1);
            std::vector<BoolMatrix> expected = {
                complete_bipartite(3, 1), complete_bipartite(1, 3), ordinal_sum(chain(2), antichain(2)),
                ordinal_sum(antichain(2), chain(2)), direct_sum(chain(3), one)};
            for (auto& e : expected) {
                bool found = false;
                for (auto& m : non) if (isomorphic(m, e)) found = true;
                CHECK(found);
            }
        }
    }
}

TEST(theorem_2_8_duals) {
    BoolMatrix p5 = riordan_matrix(S("1/(1-t)", 5), S("t/(1-t)", 5), 5);
    BoolMatrix d5 = p5.flip_transpose();
    CHECK(d5 == BoolMatrix::from_rows({"10000", "01000", "01100", "01010", "11111"}));
    RiordanPair dp = dual_pair(S("1/(1-t)", 5), S("t/(1-t)", 5), 5);
    CHECK(dp.g == S("(1+t)^4", 5));
    CHECK(dp.f == S("t/(1-t)", 5));
    // B_{1,2}* = B_{2,1}
    CHECK(isomorphic(riordan_matrix(S("1+t+t^2", 3), S("t+t^2", 3), 3).flip_transpose(), complete_bipartite(2, 1)));
    // Flip-transposes of all Riordan poset matrices up to order 9 are Riordan
    // poset matrices of the dual poset.
    for_each_riordan_poset(9, [&](const RiordanBuilder& b) {
        BoolMatrix a = b.matrix(), d = a.flip_transpose();
        CHECK(is_riordan_poset_matrix(d));
        CHECK(isomorphic(d, dual(a)));
        return true;
    });
}

TEST(theorem_3_6_shift) {
    for (int it = 0; it < 100; ++it) {
        std::size_t N = 8 + rng() % 40, m = rng() % 6;
        Series g = random_g(N), f = random_f(N);
        BoolMatrix big = riordan_matrix(g, f, N);
        RiordanPair sh = shifted_pair(g, f, m);
        CHECK(riordan_matrix(sh, N - m) == big.submatrix_range(m, N));
    }
    // Riordan poset matrices stay poset matrices after the shift (Pascal).
    RiordanPair sh = shifted_pair(S("1/(1-t)", 300), S("t/(1-t)", 300), 5);
    CHECK_EQ(poset_horizon(sh.g, sh.f, 295), std::size_t(295));
}

TEST(lemma_3_1_theorem_3_2_theorem_3_4) {
    for_each_riordan_poset(14, [&](const RiordanBuilder& b) {
        const std::size_t n = b.order();
        RiordanPair p = b.pair();
        BoolMatrix a = b.matrix();
        for (std::size_t k = 0; k + 1 < n; ++k) CHECK_EQ(a.get(k + 1, k), subdiagonal_entry(p.g, p.f, k));
        if (n >= 3) CHECK_EQ(is_chain(a), is_chain_pair(p.g, p.f));
        const bool g1 = p.g[1], f2 = p.f[2];
        if (g1 && f2)
            for (std::size_t k = 1; 2 * k + 1 < n; ++k) if (p.g[2 * k]) CHECK(p.g[2 * k + 1]);
        if (!g1 && f2)
            for (std::size_t k = 2; 2 * k < n; ++k) if (p.g[2 * k - 1]) CHECK(p.g[2 * k]);
        if (!g1 && !f2 && p.g[2])
            for (std::size_t k = 2; 2 * k < n; ++k) CHECK(p.g[2 * k]);
        return true;
    });
}

TEST(theorem_3_5_block_form) {
    // B = [I_k O; C D_m] unit binary Riordan with m < k  =>  poset matrix, D_m = I_m.
    for (std::size_t k = 1; k <= 6; ++k)
        for (std::size_t m = 0; m < k; ++m) {
            std::size_t n = k + m;
            for (int it = 0; it < 200; ++it) {
                // g, f with g_1..g_{k-1} = 0 and f_2..f_{k-1} = 0 give the I_k block.
                Series g(n), f(std::max<std::size_t>(n, 2));
                g.set(0); f.set(1);
                for (std::size_t t = k; t < n; ++t) { if (rng() & 1) g.set(t); if (t >= 2 && (rng() & 1)) f.set(t); }
                BoolMatrix b = riordan_matrix(g, f, n);
                CHECK(b.is_poset_matrix());
                CHECK(b.submatrix_range(k, n) == BoolMatrix::identity(m));
            }
        }
}

TEST(mobius_and_incidence_algebra) {
    // Boolean lattice: μ(S, T) = (-1)^{|T|-|S|}.
    const std::size_t N = 32;
    BoolMatrix p = riordan_matrix(S("1/(1-t)", N), S("t/(1-t)", N), N);
    auto mu = mobius(p);
    bool ok = true;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) {
            long long expect = ((i & j) == j) ? (((__builtin_popcountll(i) - __builtin_popcountll(j)) & 1) ? -1 : 1) : 0;
            if (mu[i][j] != expect) ok = false;
        }
    CHECK(ok);
    // Chain: μ = 1 on the diagonal, -1 on covers, 0 otherwise.
    auto mc = mobius(chain(6));
    CHECK(mc[3][3] == 1 && mc[3][2] == -1 && mc[3][1] == 0);
    // Theorem 2.7: the pattern algebra is closed and μ is the inverse of ζ.
    for_each_riordan_poset(8, [&](const RiordanBuilder& b) {
        if (b.order() != 8) return true;
        IncidenceAlgebra<long long> alg(b.matrix());
        auto x = alg.zero(), y = alg.zero();
        for (std::size_t i = 0; i < 8; ++i)
            for (std::size_t j = 0; j <= i; ++j)
                if (b.get(i, j)) { x[i][j] = long(rng() % 7) - 3; y[i][j] = long(rng() % 7) - 3; }
        CHECK(alg.contains(alg.multiply(x, y)));
        auto m = mobius(b.matrix());
        CHECK(alg.multiply(alg.zeta(), m) == alg.identity());
        CHECK(alg.inverse(alg.zeta()) == m);
        return true;
    });
}

TEST(theorem_5_1_pascal_series_parallel) {
    const std::size_t N = 20;
    BoolMatrix p = riordan_matrix(S("1/(1-t)", N), S("t/(1-t)", N), N);
    for (std::size_t n = 1; n <= N; ++n) CHECK_EQ(is_series_parallel(p.leading(n)), n <= 5);
    // SP decomposition agrees with the N-free characterisation on all posets n <= 6.
    for (std::size_t n = 1; n <= 6; ++n)
        for_each_naturally_labelled_poset(n, [&](const BoolMatrix& m) { CHECK_EQ(is_series_parallel(m), !contains_N(m)); });
    // Decompositions in the proof of Theorem 5.1.  (The paper prints
    // P4 = (1 ⊕ (1+1)) + 1, but P4 is the Boolean lattice B2 = 1 ⊕ (1+1) ⊕ 1,
    // consistent with its own P5 = 1 ⊕ (((1+1) ⊕ 1) + 1).)
    BoolMatrix one = chain(1), two = antichain(2);
    CHECK(isomorphic(p.leading(3), ordinal_sum(one, two)));
    CHECK(isomorphic(p.leading(4), ordinal_sum(ordinal_sum(one, two), one)));
    CHECK(!isomorphic(p.leading(4), direct_sum(ordinal_sum(one, two), one)));
    CHECK(isomorphic(p.leading(5), ordinal_sum(one, direct_sum(ordinal_sum(two, one), one))));
}

TEST(structure_helpers) {
    BoolMatrix p5 = riordan_matrix(S("1/(1-t)", 5), S("t/(1-t)", 5), 5);
    auto e = hasse_edges(p5);
    CHECK_EQ(e.size(), std::size_t(5));                 // 0-1, 0-2, 1-3, 2-3, 0-4
    CHECK(is_binomial(p5));                             // remark after eq. (3)
    CHECK(!is_forest(p5));                              // 3 covers both 1 and 2
    CHECK(is_forest(riordan_matrix(S("1/(1-t^2)", 9), S("t", 9), 9)));
    auto cc = chain_counts(chain(5));                    // C(5,k)
    CHECK(cc == (std::vector<std::uint64_t>{1, 5, 10, 10, 5, 1}));
    CHECK_EQ(linear_extension_count(antichain(6)), std::uint64_t(720));
    CHECK_EQ(width(antichain(6)), std::size_t(6));
    CHECK_EQ(height(chain(6)), std::size_t(6));
}
