#include <cstdio>
#include "riordan/riordan.hpp"
#include "riordan/exponential.hpp"
using namespace riordan;
static void out(const char* name, const BoolMatrix& B) {
  std::printf("%s", name);
  for (std::size_t i = 0; i < B.size(); ++i) { std::printf(" "); for (std::size_t j = 0; j <= i; ++j) std::printf("%d", (int)B(i, j)); }
  std::printf("\n");
}
int main() {
  const std::size_t N = 26; char buf[64];
  out("stirling", exponential_riordan_mod2(stirling2_egf(N), N));
  for (int a = 0; a <= 5; ++a) { std::snprintf(buf, 64, "laguerre%d", a); out(buf, exponential_riordan_mod2(laguerre_egf(a, N), N)); }
  for (int v = 0; v <= 5; ++v) { std::snprintf(buf, 64, "hermite%d", v); out(buf, exponential_riordan_mod2(hermite_egf(v, N), N)); }
  for (int m = 0; m <= 5; ++m) { std::snprintf(buf, 64, "jacobi%d", m); out(buf, exponential_riordan_mod2(jacobi_cn_sn_egf(m, N), N)); }
}
