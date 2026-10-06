// tests/test_common.hpp -- tiny self-contained test framework and shared helpers.
//
// Each tests/test_*.cpp file is its own executable (registered with CTest). A file
// defines tests with TEST(name) { ... CHECK(...); CHECK_EQ(a, b); ... } and includes
// this header once; the header supplies main(), which runs every registered test and
// returns non-zero if any check failed.
#pragma once
#include "riordan/riordan.hpp"

#include <chrono>
#include <cstdio>
#include <functional>
#include <iostream>
#include <random>
#include <set>
#include <unordered_map>

using namespace riordan;

// ---------------------------------------------------------- mini framework --
static int g_failed = 0, g_checks = 0;
#define CHECK(cond)                                                                       \
    do {                                                                                  \
        ++g_checks;                                                                       \
        if (!(cond)) { ++g_failed; std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); } \
    } while (0)
#define CHECK_EQ(a, b)                                                                    \
    do {                                                                                  \
        ++g_checks;                                                                       \
        auto _va = (a); auto _vb = (b);                                                   \
        if (!(_va == _vb)) { ++g_failed; std::cout << "  FAIL " << __FILE__ << ":" << __LINE__ \
                              << "  " #a " == " #b "  (" << _va << " vs " << _vb << ")\n"; } \
    } while (0)

struct Test { const char* name; std::function<void()> fn; };
static std::vector<Test>& registry() { static std::vector<Test> t; return t; }
struct Reg { Reg(const char* n, std::function<void()> f) { registry().push_back({n, std::move(f)}); } };
#define TEST(name) static void name(); static Reg reg_##name(#name, name); static void name()

[[maybe_unused]] static Series S(const char* e, std::size_t p) { return Series::parse(e, p); }
[[maybe_unused]] static std::mt19937_64 rng(20220101);

[[maybe_unused]] static Series random_g(std::size_t p) { Series g(p); g.set(0); for (std::size_t k = 1; k < p; ++k) if (rng() & 1) g.set(k); return g; }
[[maybe_unused]] static Series random_f(std::size_t p) { Series f(p); if (p > 1) f.set(1); for (std::size_t k = 2; k < p; ++k) if (rng() & 1) f.set(k); return f; }

// Naive reference: a_ij = [t^i] g f^j by schoolbook integer convolution mod 2.
[[maybe_unused]] static BoolMatrix naive_riordan(const Series& g, const Series& f, std::size_t n) {
    std::vector<int> col(n), fv(n);
    for (std::size_t i = 0; i < n; ++i) { col[i] = g[i]; fv[i] = f[i]; }
    BoolMatrix a(n);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) if (col[i] & 1) a.set(i, j);
        std::vector<int> nx(n, 0);
        for (std::size_t i = 0; i < n; ++i) for (std::size_t k = 0; k + i < n; ++k) nx[i + k] ^= (col[i] & fv[k]) & 1;
        col = nx;
    }
    return a;
}
[[maybe_unused]] static BoolMatrix naive_bool_square(const BoolMatrix& a) {
    BoolMatrix c(a.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < a.size(); ++j)
            for (std::size_t k = 0; k < a.size(); ++k) if (a.get(i, k) && a.get(k, j)) { c.set(i, j); break; }
    return c;
}

// ============================================================ Phase 1 =====

[[maybe_unused]] static BoolMatrix pascal_poset(std::size_t n) { auto p = pascal_family_pair(1, n); return riordan_matrix(p.g, p.f, n); }
[[maybe_unused]] static BoolMatrix two_chains_pairs(std::size_t n) {
    BoolMatrix m(n);
    for (std::size_t i = 0; i < n; ++i) { m.set(i, i); if (i % 2) m.set(i, i - 1); }
    return m;
}

int main() {
    auto t0 = std::chrono::steady_clock::now();
    for (auto& t : registry()) {
        int before = g_failed;
        auto s = std::chrono::steady_clock::now();
        t.fn();
        double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - s).count();
        std::printf("%-48s %s  (%.2fs)\n", t.name, g_failed == before ? "ok" : "FAILED", dt);
    }
    double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("\n%d checks, %d failed, %.2fs\n", g_checks, g_failed, dt);
    return g_failed ? 1 : 0;
}
