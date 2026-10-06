// riordan/decompose.hpp -- Phase 7: direct/ordinal decompositions and the question of
// whether non-Riordan posets can be "expressed in terms of Riordan posets" (end of
// Section 2 of the paper).
//
// Two readings are computed:
//   two-summand: P = A + B or P = A (+) B with A and B Riordan posets (the paper's n = 4 remark);
//   closure:     P lies in the smallest class containing the Riordan posets and closed
//                under + and (+); equivalently P is Riordan, or P = A + B / A (+) B with A, B
//                in the class (all splits are tried).
// A poset that is neither a direct sum nor an ordinal sum ("indecomposable") is in the
// closure only if it is itself Riordan.
#pragma once
#include "enumerate.hpp"

#include <algorithm>
#include <functional>
#include <unordered_map>
#include <vector>

namespace riordan {

namespace detail {
inline std::vector<std::vector<std::size_t>> graph_components(const BoolMatrix& r, bool incomparability) {
    const std::size_t n = r.size();
    std::vector<int> comp(n, -1);
    std::vector<std::vector<std::size_t>> out;
    for (std::size_t s = 0; s < n; ++s) {
        if (comp[s] >= 0) continue;
        out.emplace_back();
        std::vector<std::size_t> stack{s};
        comp[s] = static_cast<int>(out.size()) - 1;
        while (!stack.empty()) {
            const std::size_t v = stack.back();
            stack.pop_back();
            out.back().push_back(v);
            for (std::size_t w = 0; w < n; ++w) {
                if (comp[w] >= 0 || w == v) continue;
                const bool comparable = r(v, w) || r(w, v);
                if (comparable != incomparability) { comp[w] = comp[v]; stack.push_back(w); }
            }
        }
    }
    for (auto& c : out) std::sort(c.begin(), c.end());
    return out;
}
}  // namespace detail

/// Ordinal-sum factors bottom to top: components of the incomparability graph.
inline std::vector<std::vector<std::size_t>> ordinal_factors(const BoolMatrix& r) {
    auto f = detail::graph_components(r, true);
    auto below = [&](std::size_t v) { std::size_t c = 0; for (std::size_t w = 0; w < r.size(); ++w) c += r(v, w); return c; };
    std::sort(f.begin(), f.end(), [&](const auto& a, const auto& b) { return below(a[0]) < below(b[0]); });
    return f;
}
inline std::vector<std::vector<std::size_t>> components(const BoolMatrix& r) { return detail::graph_components(r, false); }
/// Neither a direct sum nor an ordinal sum of smaller posets.
inline bool is_indecomposable(const BoolMatrix& r) {
    return components(r).size() == 1 && ordinal_factors(r).size() == 1;
}

/// Every way to write P as A + B or A (+) B (unordered for +, by cut for (+)).
inline std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> binary_splits(const BoolMatrix& r) {
    std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> out;
    const auto co = components(r);
    for (std::size_t mask = 1; co.size() > 1 && mask < (std::size_t{1} << (co.size() - 1)); ++mask) {
        std::vector<std::size_t> A(co[0]), B;
        for (std::size_t c = 1; c < co.size(); ++c)
            for (auto v : co[c]) ((mask >> (c - 1) & 1) ? B : A).push_back(v);
        out.emplace_back(A, B);
    }
    const auto of = ordinal_factors(r);
    for (std::size_t cut = 1; cut < of.size(); ++cut) {
        std::vector<std::size_t> A, B;
        for (std::size_t c = 0; c < of.size(); ++c) for (auto v : of[c]) (c < cut ? A : B).push_back(v);
        out.emplace_back(A, B);
    }
    return out;
}

struct PosetRecord {
    BoolMatrix rep;               // naturally labelled representative
    bool riordan = false;
    bool two_summand = false;     // A + B or A (+) B with A, B Riordan
    bool closure = false;         // in the +/(+)-closure of the Riordan posets
    bool indecomposable = false;
};
struct ExpressibilityLevel {
    std::size_t order = 0, posets = 0, riordan = 0, non_riordan = 0;
    std::size_t two_summand = 0, closure = 0, indecomposable_non_riordan = 0;
};

/// All posets on n <= N elements up to isomorphism (grown by adding a maximal element over
/// every order ideal), classified as above. Feasible up to N = 9 (183231 posets).
class PosetCatalogue {
public:
    explicit PosetCatalogue(std::size_t N) : lv_(N + 1) {
        if (N == 0) return;
        std::vector<std::unordered_map<CanonicalForm, char, CanonicalFormHash>> riordan(N + 1);
        for (auto& l : riordan_census(N))
            for (auto& c : l.classes) riordan[l.order][c.canon] = 1;
        lv_[1].emplace(canonical_form(antichain(1)), PosetRecord{antichain(1)});
        for (std::size_t n = 2; n <= N; ++n)
            for (auto& [cf, rec] : lv_[n - 1]) grow(rec.rep, n);
        for (std::size_t n = 1; n <= N; ++n)
            for (auto& [cf, x] : lv_[n]) {
                x.riordan = riordan[n].count(cf) > 0;
                x.indecomposable = is_indecomposable(x.rep);
                x.closure = x.riordan;
                for (const auto& [A, B] : binary_splits(x.rep)) {
                    const auto& a = find(x.rep.principal(A));
                    const auto& b = find(x.rep.principal(B));
                    x.two_summand = x.two_summand || (a.riordan && b.riordan);
                    x.closure = x.closure || (a.closure && b.closure);
                }
            }
    }
    const PosetRecord& find(const BoolMatrix& m) const { return lv_[m.size()].at(canonical_form(m)); }
    const std::unordered_map<CanonicalForm, PosetRecord, CanonicalFormHash>& level(std::size_t n) const { return lv_[n]; }
    ExpressibilityLevel stats(std::size_t n) const {
        ExpressibilityLevel s;
        s.order = n;
        for (const auto& [cf, x] : lv_[n]) {
            ++s.posets;
            if (x.riordan) { ++s.riordan; continue; }
            ++s.non_riordan;
            s.two_summand += x.two_summand;
            s.closure += x.closure;
            s.indecomposable_non_riordan += x.indecomposable;
        }
        return s;
    }

private:
    void grow(const BoolMatrix& p, std::size_t n) {  // p lower triangular on n-1 elements
        const std::size_t m = p.size();
        Bits ideal(m);
        std::function<void(std::size_t)> rec = [&](std::size_t k) {
            if (k == m) {
                BoolMatrix q(n);
                for (std::size_t i = 0; i < m; ++i) p.row(i).for_each([&](std::size_t j) { q.set(i, j); });
                ideal.for_each([&](std::size_t j) { q.set(n - 1, j); });
                q.set(n - 1, n - 1);
                lv_[n].try_emplace(canonical_form(q), PosetRecord{q});
                return;
            }
            rec(k + 1);
            bool closed = true;
            for (std::size_t j = 0; j < k && closed; ++j) closed = !p(k, j) || ideal.test(j);
            if (closed) { ideal.set(k); rec(k + 1); ideal.reset(k); }
        };
        rec(0);
    }
    std::vector<std::unordered_map<CanonicalForm, PosetRecord, CanonicalFormHash>> lv_;
};

}  // namespace riordan
