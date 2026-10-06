# Independent brute-force check (shares no code with the C++ library) of the
# "expressible in terms of Riordan posets" classification for n <= 6.
from itertools import permutations, product
def riordan_matrix(g, f, n):  # a_ij = [t^i] g f^j mod 2, polys as int bitmasks
    A, col = [[0]*n for _ in range(n)], g
    for j in range(n):
        for i in range(n): A[i][j] = col >> i & 1
        nxt = 0
        for k in range(n):
            if f >> k & 1: nxt ^= col << k
        col = nxt & ((1 << n) - 1)
    return A
def is_poset(rel, n):
    return all(rel[i][i] for i in range(n)) and all(not (rel[i][j] and rel[j][k]) or rel[i][k]
               for i in range(n) for j in range(n) for k in range(n))
def canon(rel, n):  # rel[i][j] = 1 iff x_j <= x_i
    return min(tuple(rel[p[a]][p[b]] for a in range(n) for b in range(n)) for p in permutations(range(n)))
def sub(rel, idx): return [[rel[a][b] for b in idx] for a in idx]
def comps(rel, n, inc):
    seen, out = [False]*n, []
    for s in range(n):
        if seen[s]: continue
        st, cur = [s], []; seen[s] = True
        while st:
            v = st.pop(); cur.append(v)
            for w in range(n):
                if not seen[w] and w != v and bool(rel[v][w] or rel[w][v]) != inc: seen[w] = True; st.append(w)
        out.append(cur)
    return out
N = 6
riordan, allp = {}, {}
for n in range(1, N + 1):
    R = set()
    for gm, fm in product(range(1 << (n - 1)), range(1 << max(n - 2, 0))):
        g = 1 | gm << 1; f = (2 | fm << 2) if n >= 2 else 0
        A = riordan_matrix(g, f, n)
        if is_poset(A, n): R.add(canon(A, n))
    riordan[n] = R
    P = {}
    pairs = [(i, j) for i in range(n) for j in range(i)]
    for mask in range(1 << len(pairs)):
        rel = [[int(i == j) for j in range(n)] for i in range(n)]
        for b, (i, j) in enumerate(pairs):
            if mask >> b & 1: rel[i][j] = 1
        if is_poset(rel, n): P.setdefault(canon(rel, n), rel)
    allp[n] = P
closure = {}
for n in range(1, N + 1):
    stats = [0, 0, 0, 0]
    for c, rel in allp[n].items():
        isR = c in riordan[n]
        co, ic = comps(rel, n, False), comps(rel, n, True)
        two, clo = False, isR
        splits = []
        if len(co) > 1:
            for mask in range(1, 1 << (len(co) - 1)):
                A = co[0] + [v for k in range(1, len(co)) if not mask >> (k - 1) & 1 for v in co[k]]
                B = [v for k in range(1, len(co)) if mask >> (k - 1) & 1 for v in co[k]]
                splits.append((A, B))
        if len(ic) > 1:
            ic.sort(key=lambda K: sum(rel[K[0]][w] for w in range(n)))  # down-set size orders the chain of blocks
            for cut in range(1, len(ic)):
                splits.append(([v for K in ic[:cut] for v in K], [v for K in ic[cut:] for v in K]))
        for A, B in splits:
            ca, cb = canon(sub(rel, A), len(A)), canon(sub(rel, B), len(B))
            two |= ca in riordan[len(A)] and cb in riordan[len(B)]
            clo |= closure[ca] and closure[cb]
        closure[c] = clo
        if not isR:
            stats[0] += 1; stats[1] += two; stats[2] += clo; stats[3] += len(co) == 1 and len(ic) == 1
    print(n, len(allp[n]), len(riordan[n]), *stats)
