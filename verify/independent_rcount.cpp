// independent_rcount -- second, code-independent route to m(N) and r(N) for large N.
//
// Enumeration uses verify/independent (written separately from the library: its own
// bitsets, GF(2) series and extension rule); isomorphism classes use nauty from C++.
// Shares no code with include/riordan or tools/rcount.cpp, and uses a different
// 128-bit hash (FNV-1a / murmur-style chains) and a different invariant split.
//
//   g++ -std=c++20 -O3 -DMAXN=WORDSIZE -DWORDSIZE=64 -Iindependent/include -I../third_party/nauty
//       independent_rcount.cpp ../bin/nauty/*.o -o independent_rcount -pthread   (see verify/README.md)
//
//   ./independent_rcount N [tmpdir] [--threads T] [--split D]                   (disk buckets)
//   ./independent_rcount N --passes P [--checkpoint FILE] [--threads T] [--split D]
//
// Pass mode: pass p handles the matrices whose invariant (sorted multiset of (down-set
// size, up-set size), hashed with this file's own function) is p mod P; the line
// "N P pass r_pass m" is appended to the checkpoint after each pass and reused on rerun.
//
// Threads (default: all cores). The independent enumerator is NOT modified: every thread
// runs its own copy and walks the top of the tree (cheap); at the split depth D each node
// is claimed atomically by the first thread to reach it, and only the claiming thread
// descends into it (the others prune it via the enumerator's "return false"). The
// depth-first order is deterministic, so node k is the same node for every thread.
// Fast threads simply claim more nodes, which balances the load dynamically.
#include "riordan/riordan.hpp"   // verify/independent/include, NOT the library
#include "nauty.h"
#include "nautinv.h"

#include <algorithm>
#include <atomic>
#include <bit>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <string>
#include <thread>
#include <vector>

namespace fs = std::filesystem;
using u64 = std::uint64_t;

struct H128 {
    u64 a, b;
    bool operator<(const H128& o) const { return a != o.a ? a < o.a : b < o.b; }
    bool operator==(const H128& o) const { return a == o.a && b == o.b; }
};
static inline u64 fmix(u64 k) {  // murmur3 finaliser
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL;
    return k ^ (k >> 33);
}

static DEFAULTOPTIONS_DIGRAPH(base_opts);

namespace {

struct Config {
    std::size_t N = 0, P = 0, threads = 1, split = 12;
    std::string dir = "irc_tmp", ckpt;
};

// Per-thread state: no sharing on the hot path.
struct ThreadState {
    u64 m = 0;                               // order-N matrices visited
    std::vector<H128> mem;                   // pass mode: class hashes (deduplicated in steps)
    std::vector<std::vector<H128>> buf;      // disk mode: bucket buffers
};

constexpr std::size_t kBuckets = 256;

void mem_push(std::vector<H128>& mem, const H128& h) {
    if (mem.size() == mem.capacity()) {  // compact, then grow by a bounded step
        std::sort(mem.begin(), mem.end());
        mem.erase(std::unique(mem.begin(), mem.end()), mem.end());
        std::vector<H128> next;
        next.reserve(mem.size() + std::max<std::size_t>(std::size_t{1} << 22, mem.size() / 4));
        next.assign(mem.begin(), mem.end());
        mem.swap(next);
    }
    mem.push_back(h);
}

void flush(const Config& c, std::size_t tid, ThreadState& s, std::size_t k) {
    auto& v = s.buf[k];
    if (v.empty()) return;
    const std::string name = c.dir + "/t" + std::to_string(tid) + "_b" + std::to_string(k);
    FILE* f = std::fopen(name.c_str(), "ab");
    if (!f || std::fwrite(v.data(), sizeof(H128), v.size(), f) != v.size()) {
        std::fprintf(stderr, "cannot write %s\n", name.c_str());
        std::exit(1);
    }
    std::fclose(f);
    v.clear();
}

// Canonicalise one order-N matrix (rows given by the independent enumerator) and record it.
void process_leaf(const Config& c, std::size_t pass, std::size_t tid, ThreadState& s,
                  const std::vector<riordan::Bits>& rows) {
    const std::size_t n = c.N;
    ++s.m;
    u64 row[64] = {}, col[64] = {};
    for (std::size_t i = 0; i < n; ++i)
        rows[i].for_each([&](std::size_t j) { row[i] |= u64{1} << j; col[j] |= u64{1} << i; });
    u64 key[64];  // initial ordered partition by (down-set size, up-set size)
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    for (std::size_t i = 0; i < n; ++i) {
        key[i] = u64(std::popcount(row[i])) * 64 + std::popcount(col[i]);
        lab[i] = int(i);
    }
    if (c.P) {  // pass filter on an isomorphism invariant
        u64 sk[64];
        std::copy(key, key + n, sk);
        std::sort(sk, sk + n);
        u64 iv = 0x6a09e667f3bcc908ULL;
        for (std::size_t i = 0; i < n; ++i) iv = fmix(iv + sk[i] * 0x9e3779b97f4a7c15ULL + i);
        if (iv % c.P != pass) return;
    }
    std::stable_sort(lab, lab + n, [&](int x, int y) { return key[x] < key[y]; });
    for (std::size_t k = 0; k < n; ++k) ptn[k] = (k + 1 < n && key[lab[k + 1]] == key[lab[k]]) ? 1 : 0;
    graph g[MAXN], cg[MAXN];
    EMPTYGRAPH(g, 1, int(n));
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < i; ++j)
            if (row[i] >> j & 1) ADDONEARC(g, int(i), int(j), 1);
    optionblk opt = base_opts;
    opt.getcanon = TRUE;
    opt.defaultptn = FALSE;
    statsblk st;
    densenauty(g, lab, ptn, orbits, &opt, &st, 1, int(n), cg);
    u64 a = 0xcbf29ce484222325ULL, b = 0x84222325cbf29ce4ULL;  // FNV-1a style chains
    auto eat = [&](u64 w) {
        for (int sh = 0; sh < 64; sh += 8) a = (a ^ ((w >> sh) & 0xff)) * 0x100000001b3ULL;
        b = fmix(b ^ w) + 0x9e3779b97f4a7c15ULL;
    };
    eat(n);
    for (std::size_t k = 0; k < n; ++k) eat(key[lab[k]]);
    for (std::size_t k = 0; k < n; ++k) eat(u64(cg[k]));
    const H128 h{fmix(a), fmix(b ^ a)};
    if (c.P) { mem_push(s.mem, h); return; }
    auto& v = s.buf[h.b % kBuckets];
    v.push_back(h);
    if (v.size() >= 32768) flush(c, tid, s, h.b % kBuckets);
}

// Number of nodes of order D in the tree (all threads number them in the same DFS order).
std::size_t count_split_nodes(std::size_t N, std::size_t D) {
    std::size_t k = 0;
    riordan::RiordanPosetEnumerator(N).run(
        [&](std::size_t n, const riordan::Bits&, const riordan::Bits&, const std::vector<riordan::Bits>&) {
            if (n == D) { ++k; return false; }
            return true;
        });
    return k;
}

// One full traversal by T threads with atomic claiming at depth D; returns per-thread states.
void run_traversal(const Config& c, std::size_t pass, std::vector<ThreadState>& st) {
    const std::size_t D = std::min(c.split, c.N);
    const std::size_t K = count_split_nodes(c.N, D);
    std::unique_ptr<std::atomic<bool>[]> claimed(new std::atomic<bool>[K]);
    for (std::size_t k = 0; k < K; ++k) claimed[k].store(false, std::memory_order_relaxed);
    auto work = [&](std::size_t tid) {
        ThreadState& s = st[tid];
        std::size_t idx = 0;  // position among order-D nodes in DFS order
        riordan::RiordanPosetEnumerator(c.N).run(
            [&](std::size_t n, const riordan::Bits&, const riordan::Bits&, const std::vector<riordan::Bits>& rows) {
                if (n == D) {
                    const std::size_t k = idx++;
                    if (claimed[k].exchange(true, std::memory_order_acq_rel)) return false;  // someone else's
                }
                if (n == c.N) { process_leaf(c, pass, tid, s, rows); return false; }
                return true;
            });
        if (!c.P) for (std::size_t b = 0; b < kBuckets; ++b) flush(c, tid, s, b);
    };
    std::vector<std::thread> pool;
    for (std::size_t t = 1; t < c.threads; ++t) pool.emplace_back(work, t);
    work(0);
    for (auto& th : pool) th.join();
}

u64 count_unique_in_memory(std::vector<ThreadState>& st) {
    std::size_t total = 0;
    for (auto& s : st) total += s.mem.size();
    std::vector<H128> all;
    all.reserve(total);
    for (auto& s : st) {
        all.insert(all.end(), s.mem.begin(), s.mem.end());
        std::vector<H128>().swap(s.mem);
    }
    std::sort(all.begin(), all.end());
    return u64(std::unique(all.begin(), all.end()) - all.begin());
}

u64 count_unique_on_disk(const Config& c) {
    u64 r = 0;
    for (std::size_t k = 0; k < kBuckets; ++k) {
        std::vector<H128> v;
        for (std::size_t t = 0; t < c.threads; ++t) {
            const std::string name = c.dir + "/t" + std::to_string(t) + "_b" + std::to_string(k);
            if (!fs::exists(name)) continue;
            const std::size_t old = v.size(), add = fs::file_size(name) / sizeof(H128);
            v.resize(old + add);
            FILE* f = std::fopen(name.c_str(), "rb");
            if (!f || std::fread(v.data() + old, sizeof(H128), add, f) != add) {
                std::fprintf(stderr, "cannot read %s\n", name.c_str());
                std::exit(1);
            }
            std::fclose(f);
        }
        std::sort(v.begin(), v.end());
        r += u64(std::unique(v.begin(), v.end()) - v.begin());
    }
    return r;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
                     "usage: independent_rcount N [tmpdir] [--threads T] [--split D]\n"
                     "       independent_rcount N --passes P [--checkpoint FILE] [--threads T] [--split D]\n");
        return 2;
    }
    Config c;
    c.N = std::stoul(argv[1]);
    c.threads = std::max(1u, std::thread::hardware_concurrency());
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--passes" && i + 1 < argc) c.P = std::stoul(argv[++i]);
        else if (a == "--checkpoint" && i + 1 < argc) c.ckpt = argv[++i];
        else if (a == "--threads" && i + 1 < argc) c.threads = std::max<std::size_t>(1, std::stoul(argv[++i]));
        else if (a == "--split" && i + 1 < argc) c.split = std::max<std::size_t>(1, std::stoul(argv[++i]));
        else c.dir = a;
    }
    if (c.N < 1 || c.N > 64) return 2;
    nauty_check(WORDSIZE, 1, static_cast<int>(c.N), NAUTYVERSIONID);
    std::fprintf(stderr, "independent_rcount N=%zu, %zu thread(s), split depth %zu, %s\n", c.N, c.threads,
                 std::min(c.split, c.N), c.P ? "pass mode" : "disk mode");

    if (c.P) {  // in-memory, resumable
        std::vector<char> done(c.P, 0);
        u64 r = 0, mm = 0;
        if (!c.ckpt.empty() && fs::exists(c.ckpt)) {
            FILE* f = std::fopen(c.ckpt.c_str(), "r");
            unsigned long long cn, cp, cq, cr, cm;
            while (f && std::fscanf(f, "%llu %llu %llu %llu %llu", &cn, &cp, &cq, &cr, &cm) == 5)
                if (cn == c.N && cp == c.P && cq < c.P && !done[cq]) {
                    if (mm && mm != cm) { std::fprintf(stderr, "checkpoint disagrees on m\n"); return 1; }
                    done[cq] = 1; r += cr; mm = cm;
                }
            if (f) std::fclose(f);
        }
        for (std::size_t pass = 0; pass < c.P; ++pass) {
            if (done[pass]) continue;
            std::vector<ThreadState> st(c.threads);
            run_traversal(c, pass, st);
            u64 m = 0;
            for (auto& s : st) m += s.m;
            const u64 rp = count_unique_in_memory(st);
            if (mm && mm != m) { std::fprintf(stderr, "passes disagree on m\n"); return 1; }
            mm = m; r += rp;
            std::fprintf(stderr, " pass %zu/%zu: %llu classes\n", pass + 1, c.P, (unsigned long long)rp);
            if (!c.ckpt.empty()) {
                FILE* f = std::fopen(c.ckpt.c_str(), "a");
                if (!f) { std::fprintf(stderr, "cannot write checkpoint\n"); return 1; }
                std::fprintf(f, "%zu %zu %zu %llu %llu\n", c.N, c.P, pass, (unsigned long long)rp,
                             (unsigned long long)m);
                std::fclose(f);
            }
        }
        std::printf("%zu %llu %llu\n", c.N, (unsigned long long)mm, (unsigned long long)r);
        return 0;
    }

    fs::remove_all(c.dir);
    fs::create_directories(c.dir);
    std::vector<ThreadState> st(c.threads);
    for (auto& s : st) s.buf.assign(kBuckets, {});
    run_traversal(c, 0, st);
    u64 m = 0;
    for (auto& s : st) m += s.m;
    const u64 r = count_unique_on_disk(c);
    fs::remove_all(c.dir);
    std::printf("%zu %llu %llu\n", c.N, (unsigned long long)m, (unsigned long long)r);
}
