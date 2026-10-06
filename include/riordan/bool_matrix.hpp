// riordan/bool_matrix.hpp -- square (0,1)-matrices with row and column views.
//
// Convention used throughout the library (as in Cheon et al., LAA 2022):
// a poset matrix A = [a_ij] of a labelled poset {x_0, x_1, ...} has
//        a_ij = 1   iff   x_j <= x_i,
// so A is unit lower triangular whenever the labelling is a linear extension,
// and A^T is the zeta (incidence) matrix of the poset.
#pragma once
#include "bits.hpp"

#include <sstream>
#include <string>
#include <vector>

namespace riordan {

class BoolMatrix {
public:
    BoolMatrix() = default;
    explicit BoolMatrix(std::size_t n) : n_(n), rows_(n, Bits(n)), cols_(n, Bits(n)) {}

    static BoolMatrix identity(std::size_t n) {
        BoolMatrix m(n);
        for (std::size_t i = 0; i < n; ++i) m.set(i, i);
        return m;
    }
    // Rows given as strings of '0'/'1', e.g. {"1", "11", "101"} (short rows are zero-padded).
    static BoolMatrix from_rows(const std::vector<std::string>& rows) {
        BoolMatrix m(rows.size());
        for (std::size_t i = 0; i < rows.size(); ++i)
            for (std::size_t j = 0; j < rows[i].size(); ++j) {
                char c = rows[i][j];
                if (c == '1') m.set(i, j);
                else if (c != '0') throw std::invalid_argument("BoolMatrix::from_rows: expected 0/1");
            }
        return m;
    }

    std::size_t size() const { return n_; }
    bool get(std::size_t i, std::size_t j) const { return rows_[i].test(j); }
    bool operator()(std::size_t i, std::size_t j) const { return get(i, j); }
    void set(std::size_t i, std::size_t j, bool v = true) { rows_[i].set(j, v); cols_[j].set(i, v); }
    const Bits& row(std::size_t i) const { return rows_[i]; }
    const Bits& col(std::size_t j) const { return cols_[j]; }

    bool operator==(const BoolMatrix& o) const { return n_ == o.n_ && rows_ == o.rows_; }
    bool operator!=(const BoolMatrix& o) const { return !(*this == o); }

    BoolMatrix transpose() const {
        BoolMatrix t(n_);
        t.rows_ = cols_;
        t.cols_ = rows_;
        return t;
    }
    // E A^T E with E the backward identity: relabel i -> n-1-i and reverse order.
    // For a poset matrix this is the poset matrix of the dual poset.
    BoolMatrix flip_transpose() const {
        BoolMatrix t(n_);
        for (std::size_t i = 0; i < n_; ++i)
            rows_[i].for_each([&](std::size_t j) { t.set(n_ - 1 - j, n_ - 1 - i); });
        return t;
    }
    BoolMatrix leading(std::size_t m) const { return submatrix_range(0, m); }
    // Principal submatrix on rows/cols [lo, hi).
    BoolMatrix submatrix_range(std::size_t lo, std::size_t hi) const {
        BoolMatrix s(hi - lo);
        for (std::size_t i = lo; i < hi; ++i)
            for (std::size_t j = lo; j < hi; ++j) if (get(i, j)) s.set(i - lo, j - lo);
        return s;
    }
    BoolMatrix principal(const std::vector<std::size_t>& idx) const {
        BoolMatrix s(idx.size());
        for (std::size_t a = 0; a < idx.size(); ++a)
            for (std::size_t b = 0; b < idx.size(); ++b) if (get(idx[a], idx[b])) s.set(a, b);
        return s;
    }
    // Relabel: result(p[i], p[j]) = this(i, j).
    BoolMatrix permuted(const std::vector<std::size_t>& p) const {
        BoolMatrix s(n_);
        for (std::size_t i = 0; i < n_; ++i)
            rows_[i].for_each([&](std::size_t j) { s.set(p[i], p[j]); });
        return s;
    }

    // Boolean product over ({0,1}, OR, AND).
    friend BoolMatrix bool_product(const BoolMatrix& a, const BoolMatrix& b) {
        BoolMatrix c(a.n_);
        for (std::size_t i = 0; i < a.n_; ++i) {
            Bits r(a.n_);
            a.rows_[i].for_each([&](std::size_t k) { r |= b.rows_[k]; });
            r.for_each([&](std::size_t j) { c.set(i, j); });
        }
        return c;
    }

    bool is_unit_lower_triangular() const {
        for (std::size_t i = 0; i < n_; ++i) {
            if (!get(i, i)) return false;
            std::size_t l = rows_[i].last();
            if (l != i) return false;
        }
        return true;
    }
    bool is_reflexive() const {
        for (std::size_t i = 0; i < n_; ++i) if (!get(i, i)) return false;
        return true;
    }
    bool is_antisymmetric() const {
        for (std::size_t i = 0; i < n_; ++i)
            for (std::size_t j = i + 1; j < n_; ++j) if (get(i, j) && get(j, i)) return false;
        return true;
    }
    // a_ij = a_jk = 1  =>  a_ik = 1, checked as: row_j ⊆ row_i whenever a_ij = 1.
    bool is_transitive() const {
        for (std::size_t i = 0; i < n_; ++i) {
            bool ok = true;
            rows_[i].for_each([&](std::size_t j) { if (ok && !Bits::subset(rows_[j], rows_[i])) ok = false; });
            if (!ok) return false;
        }
        return true;
    }
    // Is this the matrix of a partial order (any labelling)?
    bool is_order_relation() const { return is_reflexive() && is_antisymmetric() && is_transitive(); }
    // Theorem 2.1: a unit lower triangular binary matrix is a poset matrix
    // iff it is transitive iff A = A^2 over the Boolean semiring.
    bool is_poset_matrix() const { return is_unit_lower_triangular() && is_transitive(); }

    std::string str() const {
        std::ostringstream os;
        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t j = 0; j < n_; ++j) os << (get(i, j) ? '1' : '0');
            os << '\n';
        }
        return os.str();
    }
    // Rows truncated at the diagonal: "1/11/101" -- compact for lower triangular matrices.
    std::string compact() const {
        std::string s;
        for (std::size_t i = 0; i < n_; ++i) {
            if (i) s += '/';
            for (std::size_t j = 0; j <= i; ++j) s += get(i, j) ? '1' : '0';
        }
        return s;
    }

private:
    std::size_t n_ = 0;
    std::vector<Bits> rows_, cols_;
};

}  // namespace riordan
