// Quickstart: build a Riordan poset, read off its generating functions,
// enumerate small Riordan posets and recognise a poset given abstractly.
#include "riordan/riordan.hpp"
#include <iostream>

using namespace riordan;

int main() {
    // 1. The Pascal poset P_8 = {1/(1-t), t/(1-t)}_8 mod 2 (Boolean lattice B_3).
    const std::size_t n = 8;
    Series g = Series::parse("1/(1-t)", n), f = Series::parse("t/(1-t)", n);
    BoolMatrix p8 = riordan_matrix(g, f, n);
    std::cout << p8.str() << "poset matrix: " << p8.is_poset_matrix() << "\n";

    // 2. Generating functions are recovered from any binary Riordan matrix.
    auto pair = recover_pair(p8);
    std::cout << pair->str() << "\n";

    // 3. How far is {g,f}_n a poset matrix?
    std::cout << "horizon of {1+t^2, t+t^2}: "
              << poset_horizon(Series::parse("1+t^2", 100), Series::parse("t+t^2", 100), 100) << "\n";

    // 4. Every Riordan poset matrix of order 5, each with its (g, f).
    for (auto& [m, p] : riordan_poset_matrices(5))
        std::cout << m.compact() << "   " << p.str() << "\n";

    // 5. Is the "N" poset Riordan?  a<b, c<b, c<d on elements a=0,b=1,c=2,d=3.
    BoolMatrix N = BoolMatrix::from_rows({"1000", "1110", "0010", "0011"});
    if (auto r = recognize_riordan(N)) std::cout << "N is Riordan: " << r->pair.str() << "\n" << r->matrix.str();

    // 6. Counts of non-isomorphic Riordan posets r(n).
    for (auto& level : riordan_census(9)) std::cout << "r(" << level.order << ") = " << level.classes.size() << "\n";
}
