// examples/sample_listing.cpp -- list every Riordan poset matrix of order n together with
// its generating pair (g, f), grouped into isomorphism classes. This is what m(n) and r(n)
// count: m(n) = number of matrices printed, r(n) = number of groups.
//
//   bin/sample_listing 4             # human-readable listing (used for the README)
//   bin/sample_listing 4 --markdown  # the same as Markdown tables
#include "riordan/riordan.hpp"

#include <cstdio>
#include <map>
#include <string>
#include <vector>

using namespace riordan;

namespace {

struct Entry {
    BoolMatrix matrix;
    RiordanPair pair;
};

std::string covers(const BoolMatrix& m) {
    std::string s;
    for (auto [lo, hi] : hasse_edges(m)) {
        if (!s.empty()) s += ", ";
        s += std::to_string(lo) + "<" + std::to_string(hi);
    }
    return s.empty() ? "(no relations)" : s;
}

std::string rows(const BoolMatrix& m) {
    std::string s;
    for (std::size_t i = 0; i < m.size(); ++i) {
        if (i) s += " / ";
        for (std::size_t j = 0; j < m.size(); ++j) s += m.get(i, j) ? '1' : '0';
    }
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    const std::size_t n = argc > 1 ? std::stoul(argv[1]) : 4;
    const bool markdown = argc > 2 && std::string(argv[2]) == "--markdown";

    // Group matrices by the canonical form of the poset they define.
    std::map<std::string, std::vector<Entry>> classes;
    std::vector<std::string> order;  // classes in order of first appearance
    std::size_t m = 0;
    for_each_riordan_poset(n, [&](const RiordanBuilder& b) {
        if (b.order() != n) return true;
        ++m;
        const BoolMatrix A = b.matrix();
        const std::string key = canonical_form(A).matrix().compact();
        if (!classes.count(key)) order.push_back(key);
        classes[key].push_back({A, b.pair()});
        return true;
    });

    std::printf(markdown ? "m(%zu) = %zu matrices in r(%zu) = %zu isomorphism classes\n\n"
                         : "m(%zu) = %zu Riordan poset matrices, r(%zu) = %zu isomorphism classes\n\n",
                n, m, n, classes.size());
    std::size_t k = 0, idx = 0;
    for (const auto& key : order) {
        const auto& cls = classes[key];
        ++k;
        if (markdown) {
            std::printf("**Class %zu** (%zu labelling%s), cover relations of the first: %s\n\n", k, cls.size(),
                        cls.size() == 1 ? "" : "s", covers(cls.front().matrix).c_str());
            std::printf("| # | rows of the matrix | g(t) mod t^%zu | f(t) mod t^%zu |\n|---|---|---|---|\n", n, n);
            for (const auto& e : cls)
                std::printf("| %zu | `%s` | %s | %s |\n", ++idx, rows(e.matrix).c_str(), e.pair.g.str(false).c_str(),
                            e.pair.f.str(false).c_str());
            std::printf("\n");
        } else {
            std::printf("class %zu: %zu labelling(s); cover relations: %s\n", k, cls.size(),
                        covers(cls.front().matrix).c_str());
            for (const auto& e : cls)
                std::printf("  #%-3zu %s   g = %s   f = %s\n", ++idx, rows(e.matrix).c_str(),
                            e.pair.g.str(false).c_str(), e.pair.f.str(false).c_str());
        }
    }
    return 0;
}
