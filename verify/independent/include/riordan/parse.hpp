// parse.hpp - parse expressions such as "1/(1-t)", "t/(1-t^2)", "(1+t)^4", "1+t^3+t^4"
// into GF(2) power series. Integers are reduced mod 2 and '-' equals '+'.
//
// Grammar:  expr  := term (('+'|'-') term)*
//           term  := power (('*'|'/')? power)*      (juxtaposition = multiplication)
//           power := atom ('^' '-'? integer)?
//           atom  := '(' expr ')' | 't' | integer | '-' atom
#pragma once
#include <cctype>
#include <string>

#include "gf2_series.hpp"

namespace riordan {

class SeriesParser {
 public:
  SeriesParser(const std::string& text, std::size_t prec) : p_(prec) {
    for (char c : text)
      if (!std::isspace(static_cast<unsigned char>(c))) s_ += c;
  }
  Gf2Series parse() {
    Gf2Series r = expr();
    if (pos_ != s_.size()) fail("unexpected character");
    return r;
  }

 private:
  char peek() const { return pos_ < s_.size() ? s_[pos_] : '\0'; }
  [[noreturn]] void fail(const std::string& why) const {
    throw std::invalid_argument("parse error at position " + std::to_string(pos_) + " in '" +
                                s_ + "': " + why);
  }
  long long integer() {
    if (!std::isdigit(static_cast<unsigned char>(peek()))) fail("expected an integer");
    long long v = 0;
    while (std::isdigit(static_cast<unsigned char>(peek()))) v = v * 10 + (s_[pos_++] - '0');
    return v;
  }
  Gf2Series expr() {
    Gf2Series r = term();
    while (peek() == '+' || peek() == '-') {
      ++pos_;
      r = r + term();
    }
    return r;
  }
  Gf2Series term() {
    Gf2Series r = power();
    for (;;) {
      const char c = peek();
      if (c == '*') { ++pos_; r = r * power(); }
      else if (c == '/') {
        ++pos_;
        const Gf2Series d = power();
        if (!d[0]) fail("division by a series with zero constant term");
        r = r * d.inverse();
      } else if (c == '(' || c == 't' || std::isdigit(static_cast<unsigned char>(c))) {
        r = r * power();
      } else {
        return r;
      }
    }
  }
  Gf2Series power() {
    Gf2Series b = atom();
    if (peek() == '^') {
      ++pos_;
      bool neg = false;
      if (peek() == '-') { neg = true; ++pos_; }
      const long long e = integer();
      if (neg && !b[0]) fail("negative power of a non-invertible series");
      b = b.pow(neg ? -e : e);
    }
    return b;
  }
  Gf2Series atom() {
    const char c = peek();
    if (c == '(') {
      ++pos_;
      Gf2Series r = expr();
      if (peek() != ')') fail("expected ')'");
      ++pos_;
      return r;
    }
    if (c == 't') { ++pos_; return Gf2Series::t(p_); }
    if (c == '-') { ++pos_; return atom(); }  // -x == x over GF(2)
    if (std::isdigit(static_cast<unsigned char>(c)))
      return (integer() & 1) ? Gf2Series::one(p_) : Gf2Series::zero(p_);
    fail("expected '(', 't' or an integer");
  }

  std::string s_;
  std::size_t pos_ = 0;
  std::size_t p_;
};

inline Gf2Series parse_series(const std::string& text, std::size_t prec) {
  return SeriesParser(text, prec).parse();
}

}  // namespace riordan
