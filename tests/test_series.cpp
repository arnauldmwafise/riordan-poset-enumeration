// tests/test_series.cpp -- GF(2) power series: parsing and arithmetic against naive references.
// Reference: Cheon, Curtis, Kwon, Mesinga Mwafise, Linear Algebra Appl. 632 (2022) 308-331.
#include "test_common.hpp"

TEST(series_parser) {
    CHECK(S("1/(1-t)", 10).bits().str() == "1111111111");
    CHECK(S("t/(1-t)", 6).bits().str() == "011111");
    CHECK(S("(1+t)^4", 8).bits().str() == "10001000");          // Frobenius
    CHECK(S("t/(1-t)^2", 8).bits().str() == "01010101");        // 1/(1+t)^2 = 1/(1+t^2)
    CHECK(S("1+t^3+t^4+t^6/(1-t)", 10).bits().str() == "1001101111");
    CHECK(S("3t + 2t^2 - t^3", 5).bits().str() == "01010");
    CHECK(S("1/(1-t^2)", 7).bits().str() == "1010101");
}

TEST(series_arithmetic_random) {
    for (int it = 0; it < 150; ++it) {
        std::size_t p = 1 + rng() % 400;
        Series g = random_g(p), f = random_f(std::max<std::size_t>(p, 2)).truncated(p), h = random_g(p);
        CHECK(g * g.inverse() == Series::one(p));
        CHECK((g * h).square() == g.square() * h.square());
        CHECK(g.pow(5) == g * g * g * g * g);
        if (p >= 2) {
            Series fb = f.reversion();
            CHECK(Series::compose(f, fb, p) == Series::t(p));
            CHECK(Series::compose(fb, f, p) == Series::t(p));
            Series acc(p), pw = Series::one(p);           // Horner-free reference for g(f)
            for (std::size_t k = 0; k < p; ++k) { if (g[k]) acc = acc + pw; pw = pw * f; }
            CHECK(Series::compose(g, f, p) == acc);
        }
    }
}
