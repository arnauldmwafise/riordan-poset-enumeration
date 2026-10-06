// tests/test_exponential.cpp -- binary exponential Riordan matrices (Section 5), including the errata.
// Reference: Cheon, Curtis, Kwon, Mesinga Mwafise, Linear Algebra Appl. 632 (2022) 308-331.
#include "test_common.hpp"

TEST(phase6_theorem_5_2_and_stirling) {
    for (bool g1 : {false, true}) {
        auto B = exponential_riordan_mod2(egf_from_integer_ogf(g1, 30), 30);
        CHECK(B == (g1 ? two_chains_pairs(30) : antichain(30)));
    }
    auto S5 = exponential_riordan_mod2(stirling2_egf(5), 5);   // example after Theorem 5.2
    CHECK(S5 == BoolMatrix::from_rows({"10000", "11000", "11100", "11010", "11101"}));
    CHECK(S5.is_poset_matrix() && !is_binary_riordan(S5));
    CHECK(!exponential_riordan_mod2(stirling2_egf(6), 6).is_poset_matrix());   // poset horizon 5
}

TEST(phase6_theorem_5_3_laguerre_parity_corrected) {
    // The proof of Theorem 5.3 (g1 = alpha + 1) and Theorem 5.2 give: alpha ODD -> antichain,
    // alpha EVEN -> pairs of 2-chains. The printed statement has the parities swapped.
    for (long long a = 1; a <= 9; ++a) {
        auto L = exponential_riordan_mod2(laguerre_egf(a, 40), 40);
        CHECK(L == (a % 2 ? antichain(40) : two_chains_pairs(40)));
    }
}

TEST(phase6_theorem_5_5_hermite) {
    for (long long v = 0; v <= 5; ++v)
        for (std::size_t n = 1; n <= 30; ++n) {
            auto H = exponential_riordan_mod2(hermite_egf(v, n), n);
            CHECK(H.is_poset_matrix());
            if (v % 2 == 0) { CHECK(H == antichain(n)); continue; }
            CHECK(isomorphic(H, direct_sum(pascal_poset((n + 1) / 2), pascal_poset(n / 2))));
            CHECK_EQ(is_series_parallel(H), n <= 10);
        }
}

TEST(phase6_theorem_5_6_jacobi_corrected) {
    // m odd: P_ceil(n/2) + P_floor(n/2) as printed. m even: I_floor(n/2) + P_ceil(n/2); the
    // printed I_ceil(n/2) + P_floor(n/2) is wrong for every odd n (n = 5 by hand: I_2 + P_3).
    CHECK(exponential_riordan_mod2(jacobi_cn_sn_egf(0, 5), 5) ==
          BoolMatrix::from_rows({"10000", "01000", "10100", "00010", "10001"}));
    for (long long m = 0; m <= 3; ++m)
        for (std::size_t n = 1; n <= 40; ++n) {
            auto J = exponential_riordan_mod2(jacobi_cn_sn_egf(m, n), n);
            CHECK(J.is_poset_matrix());
            const BoolMatrix want = m % 2 ? direct_sum(pascal_poset((n + 1) / 2), pascal_poset(n / 2))
                                          : direct_sum(antichain(n / 2), pascal_poset((n + 1) / 2));
            CHECK(isomorphic(J, want));
            CHECK_EQ(is_series_parallel(J), n <= 10);
        }
}

TEST(phase6_theorem_5_7_derivative) {
    // D_m = [cos t, sin t]_{m+1} with its first row and column removed (the t-derivative
    // shifts rows up; removing the first column and last row leaves b_{i+1,j+1}).
    const BoolMatrix P3s = pascal_poset(3).flip_transpose();
    for (std::size_t m = 1; m <= 12; ++m) {
        auto D = exponential_riordan_mod2(jacobi_cn_sn_egf(0, m + 1), m + 1).submatrix_range(1, m + 1);
        CHECK(D.is_poset_matrix());
        if (m <= 5) CHECK(D == antichain(m));
        else if (m <= 9) CHECK(isomorphic(D, direct_sum(antichain(m - 3), P3s)));
        CHECK_EQ(is_series_parallel(D), m <= 9);
    }
}
