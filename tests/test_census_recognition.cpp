// tests/test_census_recognition.cpp -- census r(n) and recognition of Riordan posets.
// Reference: Cheon, Curtis, Kwon, Mesinga Mwafise, Linear Algebra Appl. 632 (2022) 308-331.
#include "test_common.hpp"

TEST(r_n_matches_paper_and_recognizer) {
    // r(n) for n = 1..8 from the paper, then new values.
    const std::vector<std::size_t> paper = {1, 2, 5, 11, 33, 74, 144, 232};
    auto census = riordan_census(10);
    for (std::size_t n = 1; n <= 8; ++n) CHECK_EQ(census[n - 1].classes.size(), paper[n - 1]);
    CHECK_EQ(census[8].classes.size(), std::size_t(639));
    CHECK_EQ(census[9].classes.size(), std::size_t(1406));
    // Every enumerated matrix is a binary Riordan poset matrix with the stored pair.
    for (auto& lvl : census)
        for (auto& c : lvl.classes) {
            CHECK(is_riordan_poset_matrix(c.matrix));
            for (auto& p : c.pairs) CHECK(isomorphic(riordan_matrix(p, lvl.order), c.matrix));
        }
    // Independent route: recognise every poset on n <= 7 elements.
    const std::vector<std::size_t> p_n = {1, 2, 5, 16, 63, 318, 2045};   // OEIS A000112
    const std::vector<std::size_t> nl = {1, 2, 7, 40, 357, 4824, 96428};  // OEIS A006455
    for (std::size_t n = 1; n <= 7; ++n) {
        std::unordered_map<CanonicalForm, BoolMatrix, CanonicalFormHash> classes;
        std::size_t labelled = 0;
        for_each_naturally_labelled_poset(n, [&](const BoolMatrix& m) { ++labelled; classes.emplace(canonical_form(m), m); });
        CHECK_EQ(labelled, nl[n - 1]);
        CHECK_EQ(classes.size(), p_n[n - 1]);
        std::size_t riordan = 0;
        for (auto& [cf, m] : classes) {
            auto rec = recognize_riordan(m);
            if (!rec) continue;
            ++riordan;
            CHECK(is_riordan_poset_matrix(rec->matrix));
            CHECK(riordan_matrix(rec->pair, n) == rec->matrix);
            CHECK(isomorphic(rec->matrix, m));
        }
        CHECK_EQ(riordan, paper[n - 1]);
    }
}

TEST(recognizer_finds_all_labellings) {
    // Every Riordan poset matrix of order 7 is found by the recogniser from its class.
    auto census = riordan_census(7);
    for (auto& c : census[6].classes) {
        auto all = riordan_labellings(c.matrix, {.find_all = true});
        std::set<std::string> found;
        for (auto& l : all) found.insert(l.matrix.compact());
        for (auto& p : c.pairs) CHECK(found.count(riordan_matrix(p, 7).compact()) == 1);
        CHECK_EQ(found.size(), c.pairs.size());
    }
}
