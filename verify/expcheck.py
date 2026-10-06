# Independent exact check: b_{n,k} = n!/k! [t^n] g f^k with Fraction arithmetic.
from fractions import Fraction as Fr
from math import factorial
N = 26
def mul(a, b):
    c = [Fr(0)] * N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i): c[i + j] += x * b[j]
    return c
def matrix(g, f):
    rows = [[0] * (n + 1) for n in range(N)]
    col = g[:]
    for k in range(N):
        for n in range(k, N):
            v = col[n] * Fr(factorial(n), factorial(k))
            assert v.denominator == 1, "non-integer entry"
            rows[n][k] = int(v) % 2
        col = mul(col, f)
    return " ".join("".join(map(str, r)) for r in rows)
def binom_series(e):  # (1-t)^e, e integer (may be negative)
    c, coef = [], Fr(1)
    for n in range(N):
        c.append(coef); coef = coef * (e - n) / (n + 1) * -1
    return c
def jacobi(m):  # ordinary coefficients of sn, cn, dn from the ODEs
    sn, cn, dn = [Fr(0)] * N, [Fr(0)] * N, [Fr(0)] * N
    cn[0] = dn[0] = Fr(1)
    for n in range(N - 1):  # (n+1) x_{n+1} = [t^n] rhs
        a = sum(cn[i] * dn[n - i] for i in range(n + 1))
        b = sum(sn[i] * dn[n - i] for i in range(n + 1))
        c = sum(sn[i] * cn[n - i] for i in range(n + 1))
        sn[n + 1], cn[n + 1], dn[n + 1] = a / (n + 1), -b / (n + 1), -m * c / (n + 1)
    return sn, cn
exp = [Fr(1, factorial(n)) for n in range(N)]
t = [Fr(0)] * N; t[1] = Fr(1)
ref = {"stirling": matrix(exp, [Fr(0)] + exp[1:])}
for a in range(6):
    ref[f"laguerre{a}"] = matrix(binom_series(-a - 1), [Fr(0)] + [Fr(-1)] * (N - 1))  # t/(t-1) = -t - t^2 - ...
for v in range(6):
    h = [Fr(0)] * N
    for m in range(0, N, 2): h[m] = Fr(-v, 2) ** (m // 2) / factorial(m // 2)
    ref[f"hermite{v}"] = matrix(h, t)
for m in range(6):
    sn, cn = jacobi(m); ref[f"jacobi{m}"] = matrix(cn, sn)
ok = 0
for line in open("exp_cpp.txt"):
    name, rest = line.split(" ", 1)
    assert ref[name] == rest.strip(), name
    ok += 1
print("exponential module matches exact rational computation for", ok, "families, n <", N)
