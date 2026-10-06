// riordan/bits.hpp -- dynamic bitset tuned for GF(2) work.
//
// Bits is a fixed-length vector of bits packed into 64-bit words.  It is the
// storage type for rows and columns of binary matrices and for the
// coefficients of GF(2) power series.  Binary operators require equal lengths
// unless documented otherwise.
#pragma once
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace riordan {

class Bits {
public:
    using word = std::uint64_t;
    static constexpr std::size_t W = 64;
    static constexpr std::size_t npos = static_cast<std::size_t>(-1);

    Bits() = default;
    explicit Bits(std::size_t n) : n_(n), w_((n + W - 1) / W, 0) {}

    std::size_t size() const { return n_; }
    std::size_t words() const { return w_.size(); }
    word word_at(std::size_t k) const { return w_[k]; }   // raw 64-bit word k

    bool test(std::size_t i) const { assert(i < n_); return (w_[i / W] >> (i % W)) & 1u; }
    bool operator[](std::size_t i) const { return test(i); }
    void set(std::size_t i, bool v = true) {
        assert(i < n_);
        if (v) w_[i / W] |= word(1) << (i % W);
        else   w_[i / W] &= ~(word(1) << (i % W));
    }
    void reset(std::size_t i) { set(i, false); }
    void flip(std::size_t i) { assert(i < n_); w_[i / W] ^= word(1) << (i % W); }
    void clear() { std::fill(w_.begin(), w_.end(), 0); }

    // Grow or shrink; new bits are zero, dropped bits are discarded.
    void resize(std::size_t n) {
        n_ = n;
        w_.resize((n + W - 1) / W, 0);
        trim();
    }

    std::size_t count() const {
        std::size_t c = 0;
        for (word x : w_) c += std::popcount(x);
        return c;
    }
    bool any() const { for (word x : w_) if (x) return true; return false; }
    bool none() const { return !any(); }

    std::size_t first() const {
        for (std::size_t k = 0; k < w_.size(); ++k)
            if (w_[k]) return k * W + std::countr_zero(w_[k]);
        return npos;
    }
    std::size_t last() const {
        for (std::size_t k = w_.size(); k-- > 0;)
            if (w_[k]) return k * W + (W - 1 - std::countl_zero(w_[k]));
        return npos;
    }
    template <class F> void for_each(F&& f) const {
        for (std::size_t k = 0; k < w_.size(); ++k) {
            word x = w_[k];
            while (x) { f(k * W + std::countr_zero(x)); x &= x - 1; }
        }
    }

    Bits& operator^=(const Bits& o) { chk(o); for (std::size_t k = 0; k < w_.size(); ++k) w_[k] ^= o.w_[k]; return *this; }
    Bits& operator&=(const Bits& o) { chk(o); for (std::size_t k = 0; k < w_.size(); ++k) w_[k] &= o.w_[k]; return *this; }
    Bits& operator|=(const Bits& o) { chk(o); for (std::size_t k = 0; k < w_.size(); ++k) w_[k] |= o.w_[k]; return *this; }
    friend Bits operator^(Bits a, const Bits& b) { return a ^= b; }
    friend Bits operator&(Bits a, const Bits& b) { return a &= b; }
    friend Bits operator|(Bits a, const Bits& b) { return a |= b; }
    bool operator==(const Bits& o) const { return n_ == o.n_ && w_ == o.w_; }
    bool operator!=(const Bits& o) const { return !(*this == o); }
    bool operator<(const Bits& o) const { return n_ != o.n_ ? n_ < o.n_ : w_ < o.w_; }

    // popcount(a & b) mod 2, over the common prefix of words.
    static bool and_parity(const Bits& a, const Bits& b) {
        std::size_t m = std::min(a.w_.size(), b.w_.size());
        word x = 0;
        for (std::size_t k = 0; k < m; ++k) x ^= a.w_[k] & b.w_[k];
        return std::popcount(x) & 1;
    }
    // true iff a ⊆ b (bits of a beyond b's storage must be zero).
    static bool subset(const Bits& a, const Bits& b) {
        std::size_t m = std::min(a.w_.size(), b.w_.size());
        for (std::size_t k = 0; k < m; ++k) if (a.w_[k] & ~b.w_[k]) return false;
        for (std::size_t k = m; k < a.w_.size(); ++k) if (a.w_[k]) return false;
        return true;
    }
    // OR `src` into this (src may be shorter).
    void or_prefix(const Bits& src) {
        std::size_t m = std::min(src.w_.size(), w_.size());
        for (std::size_t k = 0; k < m; ++k) w_[k] |= src.w_[k];
        trim();
    }

    // this ^= (src << shift), truncated to size().  src may have any length.
    // This is the inner loop of carry-less (GF(2)) polynomial multiplication.
    void xor_shifted(const Bits& src, std::size_t shift) {
        if (shift >= n_) return;
        const std::size_t ws = shift / W, bs = shift % W, nw = w_.size();
        for (std::size_t k = 0; k < src.w_.size(); ++k) {
            std::size_t d = k + ws;
            if (d >= nw) break;
            word v = src.w_[k];
            if (!v) continue;
            if (bs == 0) w_[d] ^= v;
            else {
                w_[d] ^= v << bs;
                if (d + 1 < nw) w_[d + 1] ^= v >> (W - bs);
            }
        }
        trim();
    }

    // "0110..." with bit 0 first.
    std::string str() const {
        std::string s(n_, '0');
        for (std::size_t i = 0; i < n_; ++i) if (test(i)) s[i] = '1';
        return s;
    }
    static Bits from_string(const std::string& s) {
        Bits b(s.size());
        for (std::size_t i = 0; i < s.size(); ++i) {
            if (s[i] == '1') b.set(i);
            else if (s[i] != '0') throw std::invalid_argument("Bits::from_string: expected 0/1");
        }
        return b;
    }

    std::size_t hash() const {
        std::uint64_t h = 1469598103934665603ull ^ n_;
        for (word x : w_) { h ^= x; h *= 1099511628211ull; h ^= h >> 31; }
        return static_cast<std::size_t>(h);
    }

private:
    void chk(const Bits& o) const { if (o.n_ != n_) throw std::invalid_argument("Bits: size mismatch"); }
    void trim() { if (n_ % W && !w_.empty()) w_.back() &= (word(1) << (n_ % W)) - 1; }

    std::size_t n_ = 0;
    std::vector<word> w_;
};

struct BitsHash { std::size_t operator()(const Bits& b) const { return b.hash(); } };

}  // namespace riordan
