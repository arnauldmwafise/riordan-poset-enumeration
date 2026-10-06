// riordan/io.hpp -- LaTeX, Graphviz and JSON output.
#pragma once
#include "binary_riordan.hpp"
#include "poset.hpp"

#include <sstream>
#include <string>

namespace riordan {

inline std::string to_latex(const BoolMatrix& a) {
    std::ostringstream os;
    os << "\\begin{bmatrix}\n";
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < a.size(); ++j) os << (j ? " & " : "") << (a.get(i, j) ? 1 : 0);
        os << (i + 1 < a.size() ? " \\\\\n" : "\n");
    }
    os << "\\end{bmatrix}";
    return os.str();
}

// Polynomial part of a series in LaTeX, e.g. "1 + t + t^{3}".
inline std::string to_latex(const Series& s) {
    std::ostringstream os;
    bool first = true;
    s.bits().for_each([&](std::size_t k) {
        os << (first ? "" : " + ");
        first = false;
        if (k == 0) os << "1"; else if (k == 1) os << "t"; else os << "t^{" << k << "}";
    });
    if (first) os << "0";
    return os.str();
}

// Hasse diagram in Graphviz DOT (edges drawn upward).
inline std::string to_dot(const BoolMatrix& r, const std::string& name = "P") {
    std::ostringstream os;
    os << "digraph " << name << " {\n  rankdir=BT;\n  node [shape=circle];\n";
    for (std::size_t v = 0; v < r.size(); ++v) os << "  " << v << ";\n";
    for (auto [lo, hi] : hasse_edges(r)) os << "  " << lo << " -> " << hi << " [arrowhead=none];\n";
    os << "}\n";
    return os.str();
}

inline std::string json_bits(const Bits& b) {
    std::ostringstream os;
    os << '[';
    bool first = true;
    b.for_each([&](std::size_t k) { os << (first ? "" : ",") << k; first = false; });
    os << ']';
    return os.str();
}

// {"order":n,"g":{"support":[...],"precision":n},"f":{...},"rows":["1","11",...]}
inline std::string to_json(const BoolMatrix& a, const RiordanPair* p = nullptr) {
    std::ostringstream os;
    os << "{\"order\":" << a.size();
    if (p) {
        os << ",\"g\":{\"support\":" << json_bits(p->g.bits()) << ",\"precision\":" << p->g.precision() << "}";
        os << ",\"f\":{\"support\":" << json_bits(p->f.bits()) << ",\"precision\":" << p->f.precision() << "}";
    }
    os << ",\"rows\":[";
    for (std::size_t i = 0; i < a.size(); ++i) {
        os << (i ? "," : "") << '"';
        for (std::size_t j = 0; j <= i && j < a.size(); ++j) os << (a.get(i, j) ? '1' : '0');
        os << '"';
    }
    os << "]}";
    return os.str();
}

}  // namespace riordan
