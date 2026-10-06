#!/usr/bin/env python3
"""Reproduce docs/SEQUENCE_ANALYSIS.md from data/riordan_poset_counts.csv.

Tests run on m(n) (Riordan poset matrices) and r(n) (Riordan posets up to isomorphism):
  1. Berlekamp-Massey over Q: minimal linear recurrence with constant coefficients.
  2. P-recursive (holonomic) recurrences sum_i P_i(n) a(n+i) = 0, order <= 6, deg P_i <= 4.
  3. Algebraic generating functions sum_j Q_j(x) A(x)^j = 0 mod x^(N+1), degree <= 3.
  4. Growth-rate fits with residue-class offsets, and local growth rates.
Exact arithmetic (fractions / sympy) is used for 1-3. A relation is reported only if it
leaves at least 3 equations unused, so that it is a genuine prediction, not a fit.

Requires: numpy, sympy.   Usage: python3 scripts/analyze_sequences.py
"""
import csv
import os
from fractions import Fraction as F

import numpy as np
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, "..", "data", "riordan_poset_counts.csv"))))
m = [int(r["m_n"]) for r in rows]
r = [int(r["r_n"]) for r in rows]
N = len(m)


def berlekamp_massey(s):
    s = [F(x) for x in s]
    C, B, L, shift, b = [F(1)], [F(1)], 0, 1, F(1)
    for n in range(len(s)):
        d = s[n] + sum(C[i] * s[n - i] for i in range(1, L + 1))
        if d == 0:
            shift += 1
            continue
        T, coef = C[:], d / b
        C = C + [F(0)] * (len(B) + shift - len(C))
        for i in range(len(B)):
            C[i + shift] -= coef * B[i]
        if 2 * L <= n:
            L, B, b, shift = n + 1 - L, T, d, 1
        else:
            shift += 1
    return L


def holonomic_relations(s, order, deg):
    rows_ = [[F(n) ** k * s[n + i] for i in range(order + 1) for k in range(deg + 1)]
             for n in range(len(s) - order)]
    return len(rows_), (order + 1) * (deg + 1), len(sympy.Matrix(rows_).nullspace())


def algebraic_relations(s, D, k):
    x = sympy.symbols("x")
    G = sum(sympy.Integer(s[i]) * x ** (i + 1) for i in range(len(s)))
    powers = [sympy.Integer(1)]
    for _ in range(D):
        powers.append(sympy.expand(powers[-1] * G))
    cols = []
    for pw in powers:
        poly = sympy.Poly(pw, x)
        for kk in range(k + 1):
            c = [0] * (len(s) + 1)
            for (e,), v in poly.terms():
                if e + kk <= len(s):
                    c[e + kk] += v
            cols.append(c)
    M = sympy.Matrix(len(s) + 1, len(cols), lambda i, j: cols[j][i])
    return len(s) + 1, len(cols), len(M.nullspace())


def growth_fit(s, n0, with_beta, period=6):
    ns = np.arange(1, len(s) + 1)
    y = np.log(np.array(s, float))
    sel = ns >= n0
    X = [ns[sel]] + ([np.log(ns[sel])] if with_beta else [])
    X += [(ns[sel] % period == k).astype(float) for k in range(period)]
    X = np.array(X).T
    coef = np.linalg.lstsq(X, y[sel], rcond=None)[0]
    rms = float(np.sqrt(np.mean((X @ coef - y[sel]) ** 2)))
    return float(np.exp(coef[0])), (float(coef[1]) if with_beta else None), rms


def main():
    print(f"data: n = 1..{N}\n")
    print("1. Berlekamp-Massey minimal linear recurrence order (max possible = N/2 = %d):" % (N // 2))
    print(f"   m: {berlekamp_massey(m)}   r: {berlekamp_massey(r)}\n")

    print("2. P-recursive relations (order <= 6, degree <= 4, >= 3 spare equations):")
    for name, s in (("m", m), ("r", r)):
        found = [(o, d) for o in range(1, 7) for d in range(5)
                 for eq, unk, nd in [holonomic_relations(s, o, d)] if eq >= unk + 3 and nd]
        print(f"   {name}: {found if found else 'none'}")

    print("\n3. Algebraic generating functions (degree <= 3, coefficient degree <= 7, >= 3 spare):")
    for name, s in (("m", m), ("r", r)):
        found = [(D, k) for D in (1, 2, 3) for k in range(8)
                 for eq, unk, nd in [algebraic_relations(s, D, k)] if eq >= unk + 3 and nd]
        print(f"   {name}: {found if found else 'none'}")

    print("\n4. Growth (log a(n) = n log(alpha) [+ beta log n] + offset(n mod 6)):")
    for name, s in (("m", m), ("r", r)):
        for n0 in (8, 12, 16):
            a1, _, e1 = growth_fit(s, n0, False)
            a2, b2, e2 = growth_fit(s, n0, True)
            print(f"   {name}, n >= {n0:2d}: alpha = {a1:.4f} (rms {e1:.4f});"
                  f"  with n^beta: alpha = {a2:.4f}, beta = {b2:+.3f} (rms {e2:.4f})")
    print("\n   local rate (a(n)/a(n-6))^(1/6):")
    for name, s in (("m", m), ("r", r)):
        print(f"   {name}: " + " ".join(f"{n}:{(s[n-1]/s[n-7])**(1/6):.4f}" for n in range(13, N + 1)))
    print("\n   m(n) / (2^n log2 n): " + " ".join(f"{n}:{m[n-1]/(2**n*np.log2(n)):.3f}" for n in range(8, N + 1)))
    print("\n   r(n) / m(n) by n mod 4:")
    for k in range(4):
        print(f"   n = {k} mod 4: " + " ".join(f"{n}:{r[n-1]/m[n-1]:.3f}" for n in range(5, N + 1) if n % 4 == k))


if __name__ == "__main__":
    main()
