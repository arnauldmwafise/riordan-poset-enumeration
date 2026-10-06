// tests/test_expressibility.cpp -- the open question of Section 2 (expressibility in Riordan posets).
// Reference: Cheon, Curtis, Kwon, Mesinga Mwafise, Linear Algebra Appl. 632 (2022) 308-331.
#include "test_common.hpp"

TEST(phase7_expressibility_answers_open_question) {
    PosetCatalogue cat(7);
    const std::size_t expect[][7] = {{1, 1, 1, 0, 0, 0, 0}, {2, 2, 2, 0, 0, 0, 0}, {3, 5, 5, 0, 0, 0, 0},
                                     {4, 16, 11, 5, 5, 5, 0}, {5, 63, 33, 30, 19, 26, 4},
                                     {6, 318, 74, 244, 118, 166, 68}, {7, 2045, 144, 1901, 420, 785, 888}};
    for (auto& e : expect) {
        auto s = cat.stats(e[0]);
        CHECK_EQ(s.posets, e[1]); CHECK_EQ(s.riordan, e[2]); CHECK_EQ(s.non_riordan, e[3]);
        CHECK_EQ(s.two_summand, e[4]); CHECK_EQ(s.closure, e[5]); CHECK_EQ(s.indecomposable_non_riordan, e[6]);
    }
    // The smallest counterexample: N with an extra element above its "Lambda" vertex.
    auto X = BoolMatrix::from_rows({"10000", "01000", "11100", "11110", "01001"});
    CHECK(X.is_poset_matrix() && is_indecomposable(X) && !is_riordan_poset(X));
}
