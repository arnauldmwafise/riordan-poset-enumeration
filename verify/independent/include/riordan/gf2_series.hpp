// gf2_series.hpp - truncated formal power series over GF(2) (Phase 1 of the plan).
#pragma once
#include <string>
#include <vector>

#include "bits.hpp"

namespace riordan {

/// Power series c_0 + c_1 t + ... + c_{p-1} t^{p-1} + O(t^p) over GF(2); p = prec().
/// Binary operations return a result whose precision is the smaller of the inputs'.
class Gf2Series {
 public:
  Gf2Series() = default;
  explicit Gf2Series(std::size_t prec) : c_(prec) {}
  explicit Gf2Series(Bits b) : c_(std::move(b)) {}

  static Gf2Series zero(std::size_t p) { return Gf2Series(p); }
  static Gf2Series one(std::size_t p) { return monomial(0, p); }
  static Gf2Series t(std::size_t p) { return monomial(1, p); }
  static Gf2Series monomial(std::size_t k, std::size_t p) {
    Gf2Series s(p);
    if (k < p) s.c_.set(k);
    return s;
  }
  /// Integer coefficients (lowest degree first) reduced mod 2; missing ones are 0.
  static Gf2Series from_integers(const std::vector<long long>& a, std::size_t p) {
    Gf2Series s(p);
    for (std::size_t i = 0; i < a.size() && i < p; ++i) s.c_.set(i, (a[i] % 2) != 0);
    return s;
  }
  static Gf2Series from_support(const std::vector<std::size_t>& support, std::size_t p) {
    Gf2Series s(p);
    for (auto i : support) if (i < p) s.c_.set(i);
    return s;
  }
  /// "1101" -> 1 + t + t^3 + O(t^4).
  static Gf2Series from_bitstring(const std::string& bits) {
    Gf2Series s(bits.size());
    for (std::size_t i = 0; i < bits.size(); ++i) s.c_.set(i, bits[i] == '1');
    return s;
  }

  std::size_t prec() const { return c_.size(); }
  bool operator[](std::size_t i) const { return c_.test(i); }
  void set(std::size_t i, bool v = true) { c_.set(i, v); }
  const Bits& bits() const { return c_; }
  bool operator==(const Gf2Series& o) const { return c_ == o.c_; }
  /// Index of the lowest nonzero coefficient (prec() for the zero series).
  std::size_t valuation() const { return c_.next(0); }

  Gf2Series truncated(std::size_t p) const {
    if (p > prec()) throw std::invalid_argument("truncated: cannot raise precision");
    Gf2Series r(p);
    r.c_.xor_shifted(c_, 0);
    return r;
  }
  /// Raises the precision to p by declaring the new coefficients zero (internal use).
  Gf2Series padded(std::size_t p) const {
    Gf2Series r(*this);
    if (p > prec()) r.c_.resize(p);
    return r;
  }

  friend Gf2Series operator+(const Gf2Series& a, const Gf2Series& b) {
    Gf2Series r = a.truncated(std::min(a.prec(), b.prec()));
    r.c_.xor_shifted(b.c_, 0);
    return r;
  }
  /// Carry-less product; loops over the set bits of the sparser factor.
  friend Gf2Series operator*(const Gf2Series& a, const Gf2Series& b) {
    const auto p = std::min(a.prec(), b.prec());
    const bool a_sparse = a.c_.count() <= b.c_.count();
    const Gf2Series& sp = a_sparse ? a : b;
    const Gf2Series& dn = a_sparse ? b : a;
    Gf2Series r(p);
    sp.c_.for_each([&](std::size_t i) {
      if (i < p) r.c_.xor_shifted(dn.c_, i);
    });
    return r;
  }

  /// Frobenius: f(t)^2 = f(t^2) in characteristic 2, so squaring spreads the bits.
  Gf2Series square() const {
    Gf2Series r(prec());
    c_.for_each([&](std::size_t i) {
      if (2 * i < prec()) r.c_.set(2 * i);
    });
    return r;
  }

  Gf2Series pow(long long e) const {
    if (e < 0) return inverse().pow(-e);
    Gf2Series r = one(prec()), b = *this;
    while (e) {
      if (e & 1) r = r * b;
      e >>= 1;
      if (e) b = b.square();
    }
    return r;
  }

  /// Multiplicative inverse. Newton over GF(2): h <- g * h^2 doubles the precision.
  Gf2Series inverse() const {
    const auto p = prec();
    if (p == 0) return *this;
    if (!(*this)[0]) throw std::domain_error("inverse: series has zero constant term");
    Gf2Series h = one(1);
    for (std::size_t k = 1; k < p;) {
      k = std::min(2 * k, p);
      h = truncated(k) * h.padded(k).square();
    }
    return h;
  }

  /// Composition g(f(t)); requires f(0) = 0.
  Gf2Series compose(const Gf2Series& f) const {
    if (f.prec() && f[0]) throw std::domain_error("compose: inner series must have f(0) = 0");
    const auto p = std::min(prec(), f.prec());
    Gf2Series r(p);
    for (std::size_t i = p; i-- > 0;) {  // Horner
      r = r * f;
      if ((*this)[i]) r.c_.flip(0);
    }
    return r;
  }

  /// Formal derivative (only odd-degree terms survive); precision drops by one.
  Gf2Series derivative() const {
    const auto p = prec();
    if (p == 0) return *this;
    Gf2Series r(p - 1);
    for (std::size_t i = 1; i < p; i += 2)
      if ((*this)[i]) r.c_.set(i - 1);
    return r;
  }

  /// Compositional inverse fbar with f(fbar(t)) = t; requires f = t + O(t^2).
  /// Newton iteration h <- h + (f(h) - t) / f'(h), valid in characteristic 2.
  Gf2Series reversion() const {
    const auto p = prec();
    if ((p >= 1 && (*this)[0]) || (p >= 2 && !(*this)[1]))
      throw std::domain_error("reversion: series must be t + O(t^2)");
    if (p <= 2) return t(p);
    Gf2Series h = t(2);
    for (std::size_t k = 2; k < p;) {
      k = std::min(2 * k, p);
      const Gf2Series fk = truncated(k), hk = h.padded(k);
      const Gf2Series resid = fk.compose(hk) + t(k);
      const Gf2Series dinv = fk.derivative().compose(hk.truncated(k - 1)).inverse().padded(k);
      h = hk + resid * dinv;
    }
    return h;
  }

  /// Division by t^k; the k lowest coefficients must vanish.
  Gf2Series shift_down(std::size_t k) const {
    if (k > prec() || valuation() < k) throw std::domain_error("shift_down: not divisible by t^k");
    Gf2Series r(prec() - k);
    c_.for_each([&](std::size_t i) { r.c_.set(i - k); });
    return r;
  }

  std::string to_string(std::size_t max_terms = 40) const {
    std::string s;
    std::size_t shown = 0;
    bool more = false;
    c_.for_each([&](std::size_t i) {
      if (shown == max_terms) { more = true; return; }
      if (!s.empty()) s += " + ";
      s += i == 0 ? "1" : i == 1 ? "t" : "t^" + std::to_string(i);
      ++shown;
    });
    if (s.empty()) s = "0";
    if (more) s += " + ...";
    return s + " + O(t^" + std::to_string(prec()) + ")";
  }

 private:
  Bits c_;
};

}  // namespace riordan
