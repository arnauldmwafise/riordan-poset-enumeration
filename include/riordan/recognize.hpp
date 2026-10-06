// riordan/recognize.hpp -- is a given finite poset a Riordan poset?
//
// A poset P is Riordan iff some labelling makes its poset matrix a binary
// Riordan matrix.  Such a labelling is a linear extension, so we build it one
// label at a time.  Placing element x at label d fixes row d of the matrix;
// by fact F2 its entries in columns 2..d-1 are already forced by the earlier
// rows and f, so x is admissible only if its relations to the elements at
// labels 2..d-1 match the forced bits exactly.  Columns 0 and 1 then set
// g_d and f_d.  Elements with identical strict down- and up-sets (twins) are
// interchangeable, so one representative per twin class is tried.
#pragma once
#include "binary_riordan.hpp"
#include "poset.hpp"

#include <functional>
#include <map>
#include <set>

namespace riordan {

struct RiordanLabelling {
    std::vector<std::size_t> labelling;   // labelling[k] = element that receives label k
    BoolMatrix matrix;                    // the resulting Riordan poset matrix
    RiordanPair pair;                     // its generating functions (mod t^n)
};

struct RecognizeOptions {
    bool find_all = false;          // collect every distinct Riordan poset matrix of P
    std::size_t max_results = 0;    // 0 = unlimited (only with find_all)
};

inline std::vector<RiordanLabelling> riordan_labellings(const BoolMatrix& r, RecognizeOptions opt = {}) {
    if (!is_poset(r)) throw std::invalid_argument("riordan_labellings: input is not an order relation");
    const std::size_t n = r.size();
    std::vector<RiordanLabelling> out;
    if (n == 0) return out;

    // Twin classes.
    std::vector<std::size_t> twin(n);
    {
        std::map<std::pair<std::string, std::string>, std::size_t> cls;
        for (std::size_t v = 0; v < n; ++v) {
            Bits d = r.row(v), u = r.col(v);
            d.reset(v); u.reset(v);
            twin[v] = cls.emplace(std::make_pair(d.str(), u.str()), cls.size()).first->second;
        }
    }
    std::vector<Bits> strict_down(n, Bits(n));
    for (std::size_t v = 0; v < n; ++v) { strict_down[v] = r.row(v); strict_down[v].reset(v); }

    RiordanBuilder b(n);
    std::vector<std::size_t> lab;    // lab[k] = element at label k
    Bits used(n);
    std::set<std::string> seen;      // distinct matrices found (compact form)
    bool stop = false;

    std::function<void()> dfs = [&]() {
        if (stop) return;
        const std::size_t d = lab.size();
        if (d == n) {
            BoolMatrix m = b.matrix();
            if (seen.insert(m.compact()).second) {
                out.push_back({lab, m, b.pair()});
                if (!opt.find_all || (opt.max_results && out.size() >= opt.max_results)) stop = true;
            }
            return;
        }
        std::vector<char> tried(n, 0);
        auto forced = (d >= 1) ? b.forced_next_row() : std::pair<Bits, bool>{Bits(n), false};
        for (std::size_t x = 0; x < n && !stop; ++x) {
            if (used.test(x) || tried[twin[x]]) continue;
            if (!Bits::subset(strict_down[x], used)) continue;   // must be minimal among the rest
            tried[twin[x]] = 1;
            if (d == 0) {
                lab.push_back(x); used.set(x);
                dfs();
                lab.pop_back(); used.reset(x);
                continue;
            }
            // Target row d: relations of x to the elements already labelled.
            Bits row(n);
            for (std::size_t k = 0; k < d; ++k) if (r.get(x, lab[k])) row.set(k);
            row.set(d);
            Bits diff = row ^ forced.first;
            diff.reset(0);
            if (d >= 2) diff.reset(1);
            if (diff.any()) continue;            // contradicts the forced entries
            b.push_entries(row.test(0), d >= 2 ? row.test(1) : true);
            lab.push_back(x); used.set(x);
            dfs();
            lab.pop_back(); used.reset(x);
            b.pop();
        }
    };
    dfs();
    return out;
}

// First Riordan labelling of the poset, if any.
inline std::optional<RiordanLabelling> recognize_riordan(const BoolMatrix& r) {
    auto v = riordan_labellings(r);
    if (v.empty()) return std::nullopt;
    return v.front();
}
inline bool is_riordan_poset(const BoolMatrix& r) { return recognize_riordan(r).has_value(); }

}  // namespace riordan
