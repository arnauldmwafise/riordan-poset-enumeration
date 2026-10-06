// riordan -- command-line front end.
//
//   riordan build   <g> <f> <n> [--latex|--json|--dot]
//   riordan info    <g> <f> <n>
//   riordan horizon <g> <f> [cap]
//   riordan dual    <g> <f> <n>
//   riordan census  <N> [--labelled-only] [--list <n>]
//   riordan recognize <file|-> [--all]
//   riordan decide  <g> <f>            (g, f rational: "P/Q")
//   riordan toeplitz <n> | --gens a b ...
//   riordan exp     <family> <param> <n>   (stirling | laguerre | hermite | jacobi)
//   riordan express <N>
//
// <g>, <f> are GF(2) series expressions, e.g. "1/(1-t)", "t/(1-t)", "1+t^2+t^3".
// recognize reads a 0/1 matrix (one row per line; a_ij = 1 iff x_j <= x_i, any
// labelling, full square or lower-triangular rows).
#include "riordan/riordan.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace riordan;

static void usage() {
    std::cerr << "usage:\n"
                 "  riordan build   <g> <f> <n> [--latex|--json|--dot]\n"
                 "  riordan info    <g> <f> <n>\n"
                 "  riordan horizon <g> <f> [cap=1024]\n"
                 "  riordan dual    <g> <f> <n>\n"
                 "  riordan census  <N> [--labelled-only] [--list <n>]\n"
                 "  riordan recognize <file|-> [--all]\n"
                 "  riordan decide  <g> <f>                 infinite Riordan poset? (g, f rational P/Q; exact)\n"
                 "  riordan toeplitz <n>                    count/list Toeplitz posets of order n\n"
                 "  riordan toeplitz --gens a b ...         Gamma(S), Frobenius number, minimal generators\n"
                 "  riordan exp <stirling|laguerre|hermite|jacobi> <param> <n>   binary exponential Riordan matrix\n"
                 "  riordan express <N>                     all posets n <= N: which are Riordan / expressible\n";
}

static bool has(const std::vector<std::string>& a, const std::string& flag) {
    for (auto& s : a) if (s == flag) return true;
    return false;
}

static void print_info(const BoolMatrix& a) {
    std::cout << "poset matrix      : " << (a.is_poset_matrix() ? "yes" : "no") << "\n";
    if (!a.is_poset_matrix()) return;
    std::cout << "height / width    : " << height(a) << " / " << width(a) << "\n";
    auto sp = series_parallel_decomposition(a);
    std::cout << "series-parallel   : " << (sp ? *sp : std::string("no")) << "\n";
    std::cout << "forest / binomial : " << (is_forest(a) ? "yes" : "no") << " / " << (is_binomial(a) ? "yes" : "no") << "\n";
    std::cout << "Hasse edges       : " << hasse_edges(a).size() << "\n";
    if (a.size() <= 22) std::cout << "linear extensions : " << linear_extension_count(a) << "\n";
    try {
        auto cc = chain_counts(a);
        std::cout << "chains by size    :";
        for (std::size_t k = 1; k < cc.size(); ++k) std::cout << ' ' << cc[k];
        std::cout << "\n";
    } catch (const std::overflow_error&) {}
}

static BoolMatrix read_matrix(std::istream& in) {
    std::vector<std::string> rows;
    std::string line;
    while (std::getline(in, line)) {
        std::string r;
        for (char c : line) if (c == '0' || c == '1') r += c;
        if (!r.empty()) rows.push_back(r);
    }
    return BoolMatrix::from_rows(rows);
}

int main(int argc, char** argv) {
    std::vector<std::string> args(argv + 1, argv + argc);
    if (args.empty()) { usage(); return 2; }
    const std::string cmd = args[0];
    try {
        if (cmd == "build" || cmd == "info" || cmd == "dual") {
            if (args.size() < 4) { usage(); return 2; }
            std::size_t n = std::stoul(args[3]);
            std::size_t p = std::max<std::size_t>(n, 2);
            Series g = Series::parse(args[1], p), f = Series::parse(args[2], p);
            if (cmd == "dual") {
                RiordanPair d = dual_pair(g, f, n);
                std::cout << "dual: " << d.str() << "\n" << riordan_matrix(d, n).str();
                return 0;
            }
            BoolMatrix a = riordan_matrix(g, f, n);
            RiordanPair pr{g.truncated(n), f.truncated(p)};
            if (has(args, "--latex")) std::cout << to_latex(a) << "\n";
            else if (has(args, "--json")) std::cout << to_json(a, &pr) << "\n";
            else if (has(args, "--dot")) std::cout << to_dot(a);
            else {
                std::cout << pr.str() << "\n" << a.str();
                if (cmd == "build") std::cout << "poset matrix: " << (a.is_poset_matrix() ? "yes" : "no") << "\n";
            }
            if (cmd == "info") print_info(a);
            return 0;
        }
        if (cmd == "horizon") {
            if (args.size() < 3) { usage(); return 2; }
            std::size_t cap = args.size() > 3 ? std::stoul(args[3]) : 1024;
            Series g = Series::parse(args[1], cap), f = Series::parse(args[2], std::max<std::size_t>(cap, 2));
            std::size_t h = poset_horizon(g, f, cap);
            if (h == cap) std::cout << "{g,f}_n is a poset matrix for every n <= " << cap << "\n";
            else std::cout << "horizon " << h << ": {g,f}_" << h << " is a poset matrix, {g,f}_" << h + 1 << " is not\n";
            return 0;
        }
        if (cmd == "census") {
            if (args.size() < 2) { usage(); return 2; }
            std::size_t N = std::stoul(args[1]);
            bool classify = !has(args, "--labelled-only");
            std::size_t list = 0;
            for (std::size_t k = 0; k + 1 < args.size(); ++k) if (args[k] == "--list") list = std::stoul(args[k + 1]);
            auto census = riordan_census(N, classify);
            std::cout << "n\tRiordan poset matrices\tnon-isomorphic Riordan posets\n";
            for (auto& l : census)
                std::cout << l.order << '\t' << l.labelled << '\t' << (classify ? std::to_string(l.classes.size()) : "-") << "\n";
            if (list && list <= N && classify) {
                std::cout << "\nclasses of order " << list << ":\n";
                std::size_t k = 0;
                for (auto& c : census[list - 1].classes) {
                    auto sp = series_parallel_decomposition(c.matrix);
                    std::cout << "#" << ++k << "  " << c.matrix.compact() << "   SP: " << (sp ? *sp : "no")
                              << "   realisations: " << c.pairs.size() << "\n";
                    for (auto& p : c.pairs) std::cout << "      " << p.str() << "\n";
                }
            }
            return 0;
        }
        if (cmd == "recognize") {
            if (args.size() < 2) { usage(); return 2; }
            BoolMatrix r;
            if (args[1] == "-") r = read_matrix(std::cin);
            else {
                std::ifstream in(args[1]);
                if (!in) { std::cerr << "cannot open " << args[1] << "\n"; return 2; }
                r = read_matrix(in);
            }
            if (!is_poset(r)) { std::cerr << "input is not an order relation\n"; return 1; }
            auto res = riordan_labellings(r, {.find_all = has(args, "--all")});
            if (res.empty()) { std::cout << "not a Riordan poset\n"; return 0; }
            std::cout << "Riordan poset; " << res.size() << " distinct Riordan matrix(es) shown\n";
            for (auto& l : res) {
                std::cout << "labelling (label -> element):";
                for (std::size_t k = 0; k < l.labelling.size(); ++k) std::cout << ' ' << k << "->" << l.labelling[k];
                std::cout << "\n" << l.pair.str() << "\n" << l.matrix.str() << "\n";
            }
            return 0;
        }
        if (cmd == "decide") {
            if (args.size() < 3) { usage(); return 2; }
            const auto g = parse_rational(args[1]), f = parse_rational(args[2]);
            const auto r = decide_infinite(g, f);
            std::cout << "automaton states: " << r.automaton_states << ", product states explored: " << r.product_states << "\n";
            if (r.infinite) {
                std::cout << "INFINITE: {g,f}_n is a Riordan poset matrix for every n (proved)\n";
            } else {
                std::cout << "not infinite: a(" << r.i << "," << r.j << ") = a(" << r.j << "," << r.k << ") = 1 but a("
                          << r.i << "," << r.k << ") = 0\n";
                const std::size_t cap = r.i + 1;
                std::cout << "poset horizon N* = " << poset_horizon(rational_series(g, cap), rational_series(f, cap), cap) << "\n";
            }
            return 0;
        }
        if (cmd == "toeplitz") {
            if (args.size() < 2) { usage(); return 2; }
            if (args[1] == "--gens") {
                std::vector<std::size_t> S;
                for (std::size_t i = 2; i < args.size(); ++i) S.push_back(std::stoul(args[i]));
                const long long F = frobenius_number(S);
                const Bits gam = monoid_closure(S, 64);
                std::cout << "Gamma(S) up to 63:";
                gam.for_each([](std::size_t v) { std::cout << " " << v; });
                std::cout << "\nminimal generators:";
                for (auto x : minimal_generators(monoid_closure(S, 512))) std::cout << " " << x;
                std::cout << "\nFrobenius number: ";
                if (F == -2) std::cout << "none (gcd > 1: infinitely many gaps)\n";
                else if (F == -1) std::cout << "none (no gaps)\n";
                else std::cout << F << "\n";
                std::cout << "{g,t} with supp(g) = Gamma(S) is an infinite Toeplitz poset (Theorem 4.1)\n";
                return 0;
            }
            const std::size_t n = std::stoul(args[1]);
            std::size_t k = 0;
            for_each_toeplitz_support(n, [&](const Bits& b) {
                std::cout << "#" << ++k << "  supp(g) =";
                b.for_each([](std::size_t v) { std::cout << " " << v; });
                std::cout << "   minimal generators:";
                for (auto x : minimal_generators(b)) std::cout << " " << x;
                std::cout << "\n";
            });
            std::cout << k << " Toeplitz poset matrices of order " << n << "\n";
            return 0;
        }
        if (cmd == "exp") {
            if (args.size() < 4) { usage(); return 2; }
            const std::string fam = args[1];
            const long long par = std::stoll(args[2]);
            const std::size_t n = std::stoul(args[3]);
            EgfPair e;
            if (fam == "stirling") e = stirling2_egf(n);
            else if (fam == "laguerre") e = laguerre_egf(par, n);
            else if (fam == "hermite") e = hermite_egf(par, n);
            else if (fam == "jacobi") e = jacobi_cn_sn_egf(par, n);
            else { usage(); return 2; }
            const BoolMatrix B = exponential_riordan_mod2(e, n);
            std::cout << B.str();
            std::cout << "poset matrix: " << (B.is_poset_matrix() ? "yes" : "no") << "\n";
            if (B.is_poset_matrix()) {
                auto d = series_parallel_decomposition(B);
                std::cout << "series-parallel: " << (d ? *d : std::string("no")) << "\n";
            }
            return 0;
        }
        if (cmd == "express") {
            if (args.size() < 2) { usage(); return 2; }
            const std::size_t N = std::stoul(args[1]);
            PosetCatalogue cat(N);
            std::cout << "n\tposets\tRiordan\tnon-Riordan\tA+B or A(+)B, A,B Riordan\tin +/(+)-closure\tindecomposable non-Riordan\n";
            for (std::size_t n = 1; n <= N; ++n) {
                auto s = cat.stats(n);
                std::cout << n << "\t" << s.posets << "\t" << s.riordan << "\t" << s.non_riordan << "\t" << s.two_summand
                          << "\t" << s.closure << "\t" << s.indecomposable_non_riordan << "\n";
            }
            return 0;
        }
        usage();
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
