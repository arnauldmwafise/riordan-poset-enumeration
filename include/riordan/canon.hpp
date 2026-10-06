// riordan/canon.hpp -- canonical forms (isomorphism testing) for finite posets.
//
// Self-contained individualization-refinement, in the spirit of nauty but
// specialised to order relations:
//   * refinement: colour refinement using the multisets of colours of the
//     strict down-set and strict up-set of each element;
//   * individualization: branch on the vertices of the first non-singleton cell;
//   * pruning: two elements with identical strict down- and up-sets ("twins")
//     are swapped by an automorphism fixing everything else, so only one
//     representative per twin class is tried in a cell.
// The certificate is the n*n relation matrix in the canonical order, and the
// lexicographically smallest leaf certificate is kept.  Two posets are
// isomorphic iff their certificates are equal.
//
// This is exact for every input; it is fast for the poset sizes reachable by
// exhaustive enumeration (n <= ~20).  For large highly symmetric posets a
// nauty/bliss backend can be plugged in behind the same interface.
#pragma once
#include "bool_matrix.hpp"

#include <algorithm>
#include <map>
#include <optional>
#include <unordered_map>

namespace riordan {

struct CanonicalForm {
    std::size_t n = 0;
    Bits cert;                         // cert[p*n + q] = 1 iff y_q <= y_p in canonical order
    std::vector<std::size_t> perm;     // perm[v] = canonical position of element v

    bool operator==(const CanonicalForm& o) const { return n == o.n && cert == o.cert; }
    bool operator!=(const CanonicalForm& o) const { return !(*this == o); }
    // Canonical relabelling of the poset, made lower triangular when possible.
    BoolMatrix matrix() const {
        BoolMatrix m(n);
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q) if (cert.test(p * n + q)) m.set(p, q);
        return m;
    }
};
struct CanonicalFormHash { std::size_t operator()(const CanonicalForm& c) const { return c.cert.hash(); } };

namespace detail {

class PosetCanonizer {
public:
    explicit PosetCanonizer(const BoolMatrix& r) : r_(r), n_(r.size()) {
        down_.resize(n_);
        up_.resize(n_);
        for (std::size_t v = 0; v < n_; ++v) {
            r_.row(v).for_each([&](std::size_t u) { if (u != v) down_[v].push_back(u); });
            r_.col(v).for_each([&](std::size_t u) { if (u != v) up_[v].push_back(u); });
        }
        // Twin classes: equal strict down-set and strict up-set.
        std::map<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>, std::size_t> cls;
        twin_.resize(n_);
        for (std::size_t v = 0; v < n_; ++v)
            twin_[v] = cls.emplace(std::make_pair(down_[v], up_[v]), cls.size()).first->second;
    }

    CanonicalForm run() {
        CanonicalForm best;
        best.n = n_;
        if (n_ == 0) return best;
        std::vector<std::size_t> col(n_, 0);
        refine(col);
        search(col, best);
        return best;
    }

private:
    // Colour refinement to the coarsest equitable partition; colours are
    // renumbered canonically (sorted signatures), preserving the old order.
    void refine(std::vector<std::size_t>& col) const {
        std::size_t classes = count(col);
        while (true) {
            std::vector<std::pair<std::vector<std::size_t>, std::size_t>> sig(n_);
            for (std::size_t v = 0; v < n_; ++v) {
                std::vector<std::size_t> d, u;
                for (auto x : down_[v]) d.push_back(col[x]);
                for (auto x : up_[v]) u.push_back(col[x]);
                std::sort(d.begin(), d.end());
                std::sort(u.begin(), u.end());
                std::vector<std::size_t> s;
                s.reserve(d.size() + u.size() + 3);
                s.push_back(col[v]);
                s.push_back(d.size());
                s.insert(s.end(), d.begin(), d.end());
                s.push_back(u.size());
                s.insert(s.end(), u.begin(), u.end());
                sig[v] = {std::move(s), v};
            }
            std::vector<std::size_t> order(n_);
            for (std::size_t v = 0; v < n_; ++v) order[v] = v;
            std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return sig[a].first < sig[b].first; });
            std::vector<std::size_t> nc(n_);
            std::size_t id = 0;
            for (std::size_t k = 0; k < n_; ++k) {
                if (k && sig[order[k]].first != sig[order[k - 1]].first) ++id;
                nc[order[k]] = id;
            }
            // Use "number of elements with smaller colour" as the colour, so a
            // discrete partition's colours are directly the positions 0..n-1.
            std::vector<std::size_t> start(id + 1, 0);
            for (std::size_t k = 0; k < n_; ++k) if (k == 0 || nc[order[k]] != nc[order[k - 1]]) start[nc[order[k]]] = k;
            for (std::size_t v = 0; v < n_; ++v) nc[v] = start[nc[v]];
            col.swap(nc);
            std::size_t c = count(col);
            if (c == classes) return;
            classes = c;
        }
    }
    std::size_t count(const std::vector<std::size_t>& col) const {
        std::vector<char> seen(n_, 0);
        std::size_t c = 0;
        for (auto x : col) if (!seen[x]) { seen[x] = 1; ++c; }
        return c;
    }

    void search(const std::vector<std::size_t>& col, CanonicalForm& best) const {
        // Find the first (smallest colour) non-singleton cell.
        std::vector<std::size_t> size(n_, 0);
        for (auto x : col) ++size[x];
        std::size_t target = n_;
        for (std::size_t c = 0; c < n_; ++c) if (size[c] > 1) { target = c; break; }
        if (target == n_) { leaf(col, best); return; }

        std::vector<char> tried;   // twin classes already tried in this cell
        tried.assign(n_, 0);
        for (std::size_t v = 0; v < n_; ++v) {
            if (col[v] != target || tried[twin_[v]]) continue;
            tried[twin_[v]] = 1;
            // Individualize v: it keeps colour `target`, the rest of the cell moves up by one.
            std::vector<std::size_t> c2 = col;
            for (std::size_t u = 0; u < n_; ++u) if (col[u] == target && u != v) c2[u] = target + 1;
            refine(c2);
            search(c2, best);
        }
    }

    void leaf(const std::vector<std::size_t>& pos, CanonicalForm& best) const {
        Bits cert(n_ * n_);
        for (std::size_t v = 0; v < n_; ++v)
            r_.row(v).for_each([&](std::size_t u) { cert.set(pos[v] * n_ + pos[u]); });
        if (best.perm.empty() || cert < best.cert) {
            best.cert = std::move(cert);
            best.perm = pos;
        }
    }

    const BoolMatrix& r_;
    std::size_t n_;
    std::vector<std::vector<std::size_t>> down_, up_;
    std::vector<std::size_t> twin_;
};

}  // namespace detail

// Canonical form of a poset given by its order-relation matrix (any labelling;
// r(i,j) = 1 iff x_j <= x_i).
inline CanonicalForm canonical_form(const BoolMatrix& r) { return detail::PosetCanonizer(r).run(); }
inline bool isomorphic(const BoolMatrix& a, const BoolMatrix& b) {
    return a.size() == b.size() && canonical_form(a) == canonical_form(b);
}

}  // namespace riordan
