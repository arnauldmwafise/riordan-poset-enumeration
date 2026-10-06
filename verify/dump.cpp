// Writes every Riordan poset matrix of order n (one file per n) as n uint32 row masks.
#include <cstdio>
#include <string>
#include "riordan/riordan.hpp"
using namespace riordan;
int main(int argc, char** argv) {
  const std::size_t N = std::stoul(argv[1]);
  std::vector<FILE*> out(N + 1);
  for (std::size_t n = 1; n <= N; ++n) out[n] = std::fopen(("m" + std::to_string(n) + ".bin").c_str(), "wb");
  RiordanPosetEnumerator(N).run([&](std::size_t n, const Bits&, const Bits&, const std::vector<Bits>& rows) {
    for (auto& r : rows) { std::uint32_t m = 0; r.for_each([&](std::size_t j) { m |= 1u << j; }); std::fwrite(&m, 4, 1, out[n]); }
    return true;
  });
  for (std::size_t n = 1; n <= N; ++n) std::fclose(out[n]);
}
