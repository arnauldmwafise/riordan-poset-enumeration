# Count non-isomorphic posets among the dumped matrices with nauty canonical certificates.
# Streams each file matrix by matrix (constant memory apart from the certificate set) and
# stores a 128-bit BLAKE2 digest of each certificate. An isomorphism-invariant initial
# colouring (down-set size, up-set size) speeds nauty up; certificates of isomorphic posets
# still coincide because the ordered colour partition is itself canonical.
import sys, struct, time, hashlib, pynauty
lo, hi = int(sys.argv[1]), int(sys.argv[2])
for n in range(lo, hi + 1):
    t0 = time.time()
    certs, total = set(), 0
    with open(f"m{n}.bin", "rb") as fh:
        while True:
            buf = fh.read(4 * n)
            if len(buf) < 4 * n: break
            R = struct.unpack(f"<{n}I", buf)
            total += 1
            down = [bin(r).count("1") for r in R]
            up = [sum(R[i] >> j & 1 for i in range(n)) for j in range(n)]
            keys = sorted(set(zip(down, up)))
            cells = [set(v for v in range(n) if (down[v], up[v]) == k) for k in keys]
            adj = {i: [j for j in range(i) if R[i] >> j & 1] for i in range(n)}
            g = pynauty.Graph(n, directed=True, adjacency_dict=adj, vertex_coloring=cells)
            h = hashlib.blake2b(repr((keys, [len(c) for c in cells])).encode() + pynauty.certificate(g), digest_size=16)
            certs.add(h.digest())
    print(n, total, len(certs), f"{time.time()-t0:.1f}s", flush=True)
