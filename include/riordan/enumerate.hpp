// riordan/enumerate.hpp -- exhaustive generation of Riordan poset matrices.
//
// The Riordan poset matrices of order <= N form a prefix-closed tree: the
// children of {g,f}_n are the valid choices (g_n, f_n), at most 4 of them,
// and validity only needs the new row to be closed (binary_riordan.hpp, F1-F3).
// Every node is visited exactly once and carries its own generating functions.
#pragma once
#include "binary_riordan.hpp"
#include "canon.hpp"

#include <functional>
#include <unordered_map>

namespace riordan {

// Calls visit(builder) at every Riordan poset matrix of order 1..max_order.
// If visit returns false the subtree below that node is skipped.
template <class Visit>
void for_each_riordan_poset(std::size_t max_order, Visit&& visit) {
    if (max_order == 0) return;
    RiordanBuilder b(max_order);
    std::function<void()> dfs = [&]() {
        if (!visit(static_cast<const RiordanBuilder&>(b))) return;
        if (b.full()) return;
        for (const auto& c : b.candidates()) {
            b.push_entries(c.a0, c.a1);
            dfs();
            b.pop();
        }
    };
    dfs();
}

// All Riordan poset matrices of exactly order n, with their pairs.
inline std::vector<std::pair<BoolMatrix, RiordanPair>> riordan_poset_matrices(std::size_t n) {
    std::vector<std::pair<BoolMatrix, RiordanPair>> out;
    for_each_riordan_poset(n, [&](const RiordanBuilder& b) {
        if (b.order() == n) out.emplace_back(b.matrix(), b.pair());
        return true;
    });
    return out;
}

struct RiordanPosetClass {
    CanonicalForm canon;
    BoolMatrix matrix;                 // first Riordan poset matrix found in the class
    std::vector<RiordanPair> pairs;    // every (g mod t^n, f mod t^n) realising the class
};

struct CensusLevel {
    std::size_t order = 0;
    std::uint64_t labelled = 0;        // Riordan poset matrices of this order
    std::vector<RiordanPosetClass> classes;   // isomorphism classes (empty if not requested)
};

// Census of Riordan posets for orders 1..max_order.  With classify = false only
// the labelled counts are computed (much cheaper).  classes are sorted by
// canonical certificate, so the output is deterministic.
inline std::vector<CensusLevel> riordan_census(std::size_t max_order, bool classify = true) {
    std::vector<CensusLevel> lv(max_order + 1);
    std::vector<std::unordered_map<CanonicalForm, std::size_t, CanonicalFormHash>> idx(max_order + 1);
    for (std::size_t n = 0; n <= max_order; ++n) lv[n].order = n;
    for_each_riordan_poset(max_order, [&](const RiordanBuilder& b) {
        const std::size_t n = b.order();
        ++lv[n].labelled;
        if (classify) {
            BoolMatrix m = b.matrix();
            CanonicalForm cf = canonical_form(m);
            auto [it, fresh] = idx[n].emplace(cf, lv[n].classes.size());
            if (fresh) lv[n].classes.push_back({cf, m, {}});
            lv[n].classes[it->second].pairs.push_back(b.pair());
        }
        return true;
    });
    for (auto& l : lv)
        std::sort(l.classes.begin(), l.classes.end(),
                  [](const auto& a, const auto& b) { return a.canon.cert < b.canon.cert; });
    lv.erase(lv.begin());
    return lv;
}

// -------------------------------------------------------- all posets ----
// Calls visit(m) for every unit lower triangular poset matrix of order n
// (= every naturally labelled poset on n elements; OEIS A006455).
// Row k is {k} ∪ D for an order ideal D of the poset on 0..k-1.
template <class Visit>
void for_each_naturally_labelled_poset(std::size_t n, Visit&& visit) {
    if (n == 0) return;
    std::vector<Bits> rows(n, Bits(n));
    rows[0].set(0);
    std::function<void(std::size_t)> place = [&](std::size_t k) {
        if (k == n) {
            BoolMatrix m(n);
            for (std::size_t i = 0; i < n; ++i) rows[i].for_each([&](std::size_t j) { m.set(i, j); });
            visit(m);
            return;
        }
        // Enumerate order ideals of {0..k-1} by deciding elements from the top down:
        // including x forces its whole down-set.
        Bits cur(n);
        std::function<void(long)> ideal = [&](long x) {
            if (x < 0) {
                rows[k] = cur;
                rows[k].set(k);
                place(k + 1);
                return;
            }
            if (cur.test(std::size_t(x))) { ideal(x - 1); return; }   // already forced
            ideal(x - 1);                                              // exclude x
            Bits saved = cur;
            cur |= rows[std::size_t(x)];                               // include x and its down-set
            ideal(x - 1);
            cur = saved;
        };
        ideal(long(k) - 1);
    };
    place(1);
}

}  // namespace riordan
