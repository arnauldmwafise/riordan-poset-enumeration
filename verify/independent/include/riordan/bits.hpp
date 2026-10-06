// bits.hpp - dynamic bitset with word-level operations used throughout the library.
#pragma once
#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace riordan {

/// Dynamic bitset. Invariant: bits at positions >= size() are always zero,
/// so word-level operations between bitsets of different sizes are safe.
class Bits {
 public:
  Bits() = default;
  explicit Bits(std::size_t n) : n_(n), w_((n + 63) / 64, 0) {}

  std::size_t size() const { return n_; }
  void resize(std::size_t n) {
    n_ = n;
    w_.resize((n + 63) / 64, 0);
    clear_tail();
  }

  bool test(std::size_t i) const { return i < n_ && ((w_[i >> 6] >> (i & 63)) & 1u); }
  void set(std::size_t i, bool v = true) {
    check(i);
    const std::uint64_t m = std::uint64_t{1} << (i & 63);
    if (v) w_[i >> 6] |= m; else w_[i >> 6] &= ~m;
  }
  void flip(std::size_t i) {
    check(i);
    w_[i >> 6] ^= std::uint64_t{1} << (i & 63);
  }
  void reset() { std::fill(w_.begin(), w_.end(), 0); }

  std::size_t count() const {
    std::size_t c = 0;
    for (auto x : w_) c += static_cast<std::size_t>(std::popcount(x));
    return c;
  }
  bool none() const {
    for (auto x : w_) if (x) return false;
    return true;
  }

  /// Lowest set bit at position >= from, or size() if there is none.
  std::size_t next(std::size_t from) const {
    if (from >= n_) return n_;
    std::size_t wi = from >> 6;
    std::uint64_t x = w_[wi] & (~std::uint64_t{0} << (from & 63));
    for (;;) {
      if (x) return std::min(n_, (wi << 6) + static_cast<std::size_t>(std::countr_zero(x)));
      if (++wi >= w_.size()) return n_;
      x = w_[wi];
    }
  }
  template <class F>
  void for_each(F&& fn) const {
    for (auto i = next(0); i < n_; i = next(i + 1)) fn(i);
  }

  /// this ^= (src << shift), truncated to size(). Core of GF(2) multiplication.
  void xor_shifted(const Bits& src, std::size_t shift) {
    const std::size_t ws = shift >> 6, bs = shift & 63;
    for (std::size_t i = 0; i < src.w_.size(); ++i) {
      const std::size_t d = i + ws;
      if (d >= w_.size()) break;
      const std::uint64_t v = src.w_[i];
      if (!v) continue;
      w_[d] ^= v << bs;
      if (bs && d + 1 < w_.size()) w_[d + 1] ^= v >> (64 - bs);
    }
    clear_tail();
  }

  bool operator==(const Bits& o) const { return n_ == o.n_ && w_ == o.w_; }

  /// Parity of |a AND b|, i.e. the GF(2) dot product.
  friend bool parity_and(const Bits& a, const Bits& b) {
    std::uint64_t acc = 0;
    const auto m = std::min(a.w_.size(), b.w_.size());
    for (std::size_t i = 0; i < m; ++i) acc ^= a.w_[i] & b.w_[i];
    return std::popcount(acc) & 1;
  }
  /// First position set in a but not in b, or a.size() if a is a subset of b.
  friend std::size_t first_and_not(const Bits& a, const Bits& b) {
    for (std::size_t i = 0; i < a.w_.size(); ++i) {
      const std::uint64_t bw = i < b.w_.size() ? b.w_[i] : 0;
      if (const std::uint64_t x = a.w_[i] & ~bw)
        return (i << 6) + static_cast<std::size_t>(std::countr_zero(x));
    }
    return a.n_;
  }
  friend bool is_subset(const Bits& a, const Bits& b) { return first_and_not(a, b) == a.n_; }

  std::string str() const {
    std::string s(n_, '0');
    for_each([&](std::size_t i) { s[i] = '1'; });
    return s;
  }

 private:
  void check(std::size_t i) const {
    if (i >= n_) throw std::out_of_range("Bits: index out of range");
  }
  void clear_tail() {
    if (n_ & 63) w_.back() &= (std::uint64_t{1} << (n_ & 63)) - 1;
  }
  std::size_t n_ = 0;
  std::vector<std::uint64_t> w_;
};

}  // namespace riordan
