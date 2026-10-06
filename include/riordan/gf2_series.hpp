// riordan/gf2_series.hpp -- truncated formal power series over GF(2).
//
// A Series stores the coefficients c_0 .. c_{p-1} of a power series known
// modulo t^p (p = precision()).  Results of binary operations carry the
// smaller of the two precisions, so precision is never silently invented.
//
// Algorithms
//   multiply      carry-less shift/XOR over 64-bit words, O(p * popcount / 64)
//   square        Frobenius: (sum c_k t^k)^2 = sum c_k t^{2k}, linear time
//   inverse       Newton iteration h <- g h^2 (char-2 form of h(2 - g h))
//   compose g(h)  even/odd split g(t) = E(t^2) + t O(t^2) gives
//                 g(h) = E(h)^2 + h * O(h)^2  -> O(M(p) log p)
//   revert f      Newton h <- h + (f(h) + t) / f'(h)
//
// The parser accepts expressions such as  "1/(1-t)",  "t/(1-t)^2",
// "1+t^3+t^4+t^6/(1-t)",  "(1+t)^4".  Integer constants are reduced mod 2 and
// '-' is the same as '+'.
#pragma once
#include "bits.hpp"

#include <cctype>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

namespace riordan {

class Series {
public:
    Series() = default;
    explicit Series(std::size_t precision) : c_(precision) {}
    Series(Bits coeffs) : c_(std::move(coeffs)) {}

    static Series zero(std::size_t p) { return Series(p); }
    static Series one(std::size_t p) { Series s(p); if (p) s.c_.set(0); return s; }
    static Series monomial(std::size_t k, std::size_t p) { Series s(p); if (k < p) s.c_.set(k); return s; }
    static Series t(std::size_t p) { return monomial(1, p); }
    // 1/(1 - t^k)
    static Series geometric(std::size_t k, std::size_t p) {
        Series s(p);
        for (std::size_t i = 0; i < p; i += k) s.c_.set(i);
        return s;
    }
    static Series from_exponents(const std::vector<std::size_t>& e, std::size_t p) {
        Series s(p);
        for (auto k : e) if (k < p) s.c_.flip(k);
        return s;
    }
    // Reduce an integer coefficient sequence mod 2 (e.g. pasted from OEIS).
    template <class Int>
    static Series from_integers(const std::vector<Int>& a) {
        Series s(a.size());
        for (std::size_t i = 0; i < a.size(); ++i) if ((a[i] % 2 + 2) % 2) s.c_.set(i);
        return s;
    }
    static Series parse(const std::string& expr, std::size_t p);

    std::size_t precision() const { return c_.size(); }
    bool operator[](std::size_t k) const { return k < c_.size() && c_.test(k); }
    void set(std::size_t k, bool v = true) { c_.set(k, v); }
    const Bits& bits() const { return c_; }
    Bits& bits() { return c_; }

    std::size_t valuation() const { return c_.first(); }   // Bits::npos if zero
    std::size_t degree() const { return c_.last(); }       // highest known nonzero
    bool is_zero() const { return c_.none(); }

    Series truncated(std::size_t p) const {
        Series s = *this;
        if (p < s.precision()) s.c_.resize(p);
        return s;
    }
    // Pads with zeros: only valid when the series is known to be a polynomial.
    Series padded(std::size_t p) const { Series s = *this; s.c_.resize(p); return s; }

    bool operator==(const Series& o) const { return c_ == o.c_; }
    bool operator!=(const Series& o) const { return !(*this == o); }

    friend Series operator+(const Series& a, const Series& b) {
        std::size_t p = std::min(a.precision(), b.precision());
        Series r = a.truncated(p);
        r.c_ ^= b.truncated(p).c_;
        return r;
    }
    friend Series operator*(const Series& a, const Series& b) { return mul(a, b, std::min(a.precision(), b.precision())); }

    static Series mul(const Series& a, const Series& b, std::size_t p) {
        p = std::min({p, a.precision(), b.precision()});
        const Series* x = &a; const Series* y = &b;
        if (x->c_.count() > y->c_.count()) std::swap(x, y);   // iterate the sparser one
        Series r(p);
        x->c_.for_each([&](std::size_t k) { if (k < p) r.c_.xor_shifted(y->c_, k); });
        return r;
    }

    Series square() const {
        Series r(precision());
        c_.for_each([&](std::size_t k) { if (2 * k < r.precision()) r.c_.set(2 * k); });
        return r;
    }

    Series pow(std::uint64_t e) const {
        Series result = one(precision()), base = *this;
        while (e) {
            if (e & 1) result = result * base;
            e >>= 1;
            if (e) base = base.square();
        }
        return result;
    }

    // Multiplicative inverse; requires c_0 = 1.
    Series inverse() const {
        if (!(*this)[0]) throw std::domain_error("Series::inverse: constant term is 0");
        const std::size_t p = precision();
        Series h = one(std::min<std::size_t>(1, p));
        std::size_t k = 1;
        while (k < p) {
            k = std::min(2 * k, p);
            Series hk = h.padded(k);              // h is exact mod t^{old k}
            h = mul(truncated(k), hk.square(), k);  // g h^2
        }
        return h.truncated(p);
    }
    friend Series operator/(const Series& a, const Series& b) { return a * b.inverse(); }

    // Multiply by t^m (precision grows by m) / divide by t^m (low terms must be 0).
    Series shift_up(std::size_t m) const {
        Series r(precision() + m);
        r.c_.xor_shifted(c_, m);
        return r;
    }
    Series shift_down(std::size_t m) const {
        if (m > precision()) throw std::domain_error("shift_down beyond precision");
        for (std::size_t i = 0; i < m; ++i)
            if ((*this)[i]) throw std::domain_error("shift_down: series not divisible by t^m");
        Series r(precision() - m);
        for (std::size_t i = m; i < precision(); ++i) if ((*this)[i]) r.c_.set(i - m);
        return r;
    }

    // Formal derivative: coefficient of t^{k-1} is k*c_k mod 2.
    Series derivative() const {
        Series r(precision() ? precision() - 1 : 0);
        for (std::size_t k = 1; k < precision(); k += 2) if ((*this)[k]) r.c_.set(k - 1);
        return r;
    }

    // g(h) with h(0) = 0, result mod t^p.
    static Series compose(const Series& g, const Series& h, std::size_t p) {
        if (p == 0) return Series(0);
        if (h[0]) throw std::domain_error("compose: inner series must have zero constant term");
        // h has valuation >= 1, so g(h) mod t^p depends only on g mod t^p and
        // h mod t^p; the result is known to min(p, prec g, prec h).
        p = std::min({p, g.precision(), h.precision()});
        return compose_rec(g.truncated(p), h.truncated(p), p);
    }
    Series operator()(const Series& h) const { return compose(*this, h, std::min(precision(), h.precision())); }

    // Compositional inverse of f = t + f_2 t^2 + ...
    Series reversion() const {
        const std::size_t p = precision();
        if (p < 2 || (*this)[0] || !(*this)[1])
            throw std::domain_error("reversion: need f = t + O(t^2)");
        Series h = t(2);
        std::size_t k = 2;
        Series fp = derivative();
        while (k < p) {
            k = std::min(2 * k, p);
            Series hk = h.padded(k);
            Series err = compose(truncated(k), hk, k) + t(k);       // f(h) - t
            // err has valuation >= k_old, so 1/f'(h) is only needed mod t^{k-k_old};
            // padding the (possibly one-shorter) derivative series is therefore exact.
            Series d = compose(fp.truncated(k), hk, k).inverse().padded(k);
            h = hk + mul(err, d, k);
        }
        return h.truncated(p);
    }

    // Human-readable: "1 + t + t^3 + O(t^8)".
    std::string str(bool show_order = true) const {
        std::ostringstream os;
        bool first = true;
        c_.for_each([&](std::size_t k) {
            if (!first) os << " + ";
            first = false;
            if (k == 0) os << "1";
            else if (k == 1) os << "t";
            else os << "t^" << k;
        });
        if (first) os << "0";
        if (show_order) os << " + O(t^" << precision() << ")";
        return os.str();
    }

private:
    static Series horner(const Series& g, const Series& h, std::size_t p) {
        Series r(p);
        for (std::size_t k = g.precision(); k-- > 0;) {
            r = mul(r, h, p);
            if (g[k]) r.c_.flip(0);
        }
        return r;
    }
    static Series compose_rec(const Series& g, const Series& h, std::size_t p) {
        if (g.precision() <= 32 || p <= 64) return horner(g, h, p);
        const std::size_t q = (p + 1) / 2;
        Series E((g.precision() + 1) / 2), O(g.precision() / 2);
        g.c_.for_each([&](std::size_t k) {
            if (k % 2 == 0) E.c_.set(k / 2); else O.c_.set(k / 2);
        });
        Series hq = h.truncated(q);
        Series e = compose_rec(E.truncated(std::min(E.precision(), q)), hq, q).padded(p).square();
        Series o = compose_rec(O.truncated(std::min(O.precision(), q)), hq, q).padded(p).square();
        return e.truncated(p) + mul(h, o, p);
    }

    Bits c_;
};

// ---------------------------------------------------------------- parser --
namespace detail {
class SeriesParser {
public:
    SeriesParser(const std::string& s, std::size_t p) : s_(s), p_(p) {}
    Series run() {
        Series r = expr();
        skip();
        if (i_ != s_.size()) fail("unexpected character");
        return r;
    }

private:
    // expr   := term (('+'|'-') term)*
    // term   := factor (('*'|'/'|<implicit>) factor)*
    // factor := atom ('^' integer)?
    // atom   := integer | 't' | 'x' | 'z' | '(' expr ')'
    Series expr() {
        skip();
        if (peek() == '-' || peek() == '+') ++i_;   // leading sign is irrelevant mod 2
        Series r = term();
        while (true) {
            skip();
            if (peek() == '+' || peek() == '-') { ++i_; r = r + term(); }
            else return r;
        }
    }
    Series term() {
        Series r = factor();
        while (true) {
            skip();
            char c = peek();
            if (c == '*') { ++i_; r = r * factor(); }
            else if (c == '/') { ++i_; r = r / factor(); }
            else if (c == '(' || c == 't' || c == 'x' || c == 'z' || std::isdigit((unsigned char)c)) r = r * factor();
            else return r;
        }
    }
    Series factor() {
        Series a = atom();
        skip();
        if (peek() == '^' || (peek() == '*' && i_ + 1 < s_.size() && s_[i_ + 1] == '*')) {
            i_ += (peek() == '^') ? 1 : 2;
            skip();
            std::uint64_t e = number();
            a = a.pow(e);
        }
        return a;
    }
    Series atom() {
        skip();
        char c = peek();
        if (c == '(') {
            ++i_;
            Series r = expr();
            skip();
            if (peek() != ')') fail("expected ')'");
            ++i_;
            return r;
        }
        if (c == 't' || c == 'x' || c == 'z') { ++i_; return Series::t(p_); }
        if (std::isdigit((unsigned char)c)) {
            std::uint64_t v = number();
            return (v & 1) ? Series::one(p_) : Series::zero(p_);
        }
        fail("expected a number, t or '('");
        return Series();
    }
    std::uint64_t number() {
        skip();
        if (!std::isdigit((unsigned char)peek())) fail("expected integer");
        std::uint64_t v = 0;
        while (std::isdigit((unsigned char)peek())) v = v * 10 + (s_[i_++] - '0');
        return v;
    }
    void skip() { while (i_ < s_.size() && std::isspace((unsigned char)s_[i_])) ++i_; }
    char peek() const { return i_ < s_.size() ? s_[i_] : '\0'; }
    [[noreturn]] void fail(const char* msg) {
        throw std::invalid_argument(std::string("Series::parse: ") + msg + " at position " +
                                    std::to_string(i_) + " in \"" + s_ + "\"");
    }

    const std::string& s_;
    std::size_t p_, i_ = 0;
};
}  // namespace detail

inline Series Series::parse(const std::string& expr, std::size_t p) {
    return detail::SeriesParser(expr, p).run();
}

}  // namespace riordan
