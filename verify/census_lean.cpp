// Memory-light count of r(N) with the library's enumerator and canonical form:
// keeps only a 128-bit hash of each canonical certificate at order N.
#include <cstdio>
#include <string>
#include <unordered_set>
#include "riordan/riordan.hpp"
using namespace riordan;
struct H { std::size_t operator()(const std::pair<std::uint64_t, std::uint64_t>& p) const { return p.first ^ (p.second * 0x9E3779B97F4A7C15ULL); } };
int main(int argc, char** argv) {
  const std::size_t N = std::stoul(argv[1]);
  std::unordered_set<std::pair<std::uint64_t, std::uint64_t>, H> seen;
  std::uint64_t labelled = 0;
  for_each_riordan_poset(N, [&](const RiordanBuilder& b) {
    if (b.order() < N) return true;
    ++labelled;
    const Bits& c = canonical_form(b.matrix()).cert;
    std::uint64_t h1 = 1469598103934665603ULL, h2 = 0x9E3779B97F4A7C15ULL;  // two independent FNV/mix hashes
    for (std::size_t i = 0; i < c.size(); ++i) if (c.test(i)) {
      h1 = (h1 ^ i) * 1099511628211ULL;
      h2 = (h2 ^ (i + 0x5bd1e995)) * 0xff51afd7ed558ccdULL; h2 ^= h2 >> 33;
    }
    seen.insert({h1, h2});
    return true;
  });
  std::printf("%zu %llu %zu\n", N, (unsigned long long)labelled, seen.size());
}
