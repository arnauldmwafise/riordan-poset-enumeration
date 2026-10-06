// rcount -- count r(N), the number of non-isomorphic Riordan posets on N elements,
// at scale: parallel tree enumeration, nauty canonical forms called from C++, and
// disk-backed deduplication of 128-bit class hashes.
//
//   rcount N [--threads T] [--split D] [--canon nauty|library] [--buckets B] [--tmp DIR]
//            [--passes P] [--checkpoint FILE]
//
// Pipeline
//   1. Enumerate the tree of Riordan poset matrices (struct Node below, the one-word-per-row
//      form of RiordanBuilder) down to order D and
//      record the (a0,a1) choices leading to every node of order D: these are the jobs.
//   2. T worker threads take jobs from an atomic counter, replay the choices, and finish
//      the subtree to order N. Every order-N matrix is canonicalised (nauty by default) and
//      reduced to a 128-bit hash, which goes into one of B per-thread buffers chosen by the
//      hash; full buffers are appended to files DIR/t<thread>_b<bucket>.bin.
//   3. For each bucket: load its files, sort, count distinct hashes. r(N) is the sum.
// Memory is O(m(N)/B) hashes at a time, so N is limited by time, not RAM.
//
// --passes P (memory mode, no disk): the enumeration is repeated P times. Pass p
// canonicalises only the matrices whose isomorphism INVARIANT (the sorted multiset of
// (down-set size, up-set size) pairs) hashes to p mod P, and deduplicates them in RAM.
// Isomorphic matrices share the invariant, so every class is counted in exactly one pass;
// canonicalisation is still done once per matrix overall, only the cheap enumeration
// repeats. Peak memory is about 16 bytes x 2.5 x r(N)/P per thread (compaction in bounded steps).
//
// --checkpoint FILE (with --passes): after each pass the line "N P pass r_pass m" is appended
// to FILE; a rerun with the same N and P skips the passes already recorded. Every pass
// enumerates all m(N) matrices, so the m column of every line must agree (checked).
//
// nauty must be compiled with MAXN=WORDSIZE=64 and thread-local storage (see Makefile).
#include "riordan/riordan.hpp"
#include "nauty.h"
#include "nautinv.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using namespace riordan;
namespace fs = std::filesystem;

struct Hash128 {
    std::uint64_t a = 0, b = 0;
    bool operator<(const Hash128& o) const { return a != o.a ? a < o.a : b < o.b; }
    bool operator==(const Hash128& o) const { return a == o.a && b == o.b; }
};

static inline std::uint64_t mix64(std::uint64_t x) {  // splitmix64 finaliser
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
// Two independently seeded position-dependent chains -> 128 bits.
struct Hasher {
    std::uint64_t a = 0x243F6A8885A308D3ULL, b = 0x13198A2E03707344ULL, i = 0;
    void add(std::uint64_t w) {
        ++i;
        a = mix64(a ^ mix64(w + i * 0x9E3779B97F4A7C15ULL));
        b = mix64(b + mix64(w ^ (i * 0xD6E8FEB86659FD93ULL + 0x632BE59BD9B4E019ULL)));
    }
    Hash128 get() const { return {mix64(a ^ i), mix64(b + i)}; }
};


// ------------------------------------------------------------- fast enumerator --
// The same tree as RiordanBuilder (riordan/binary_riordan.hpp), specialised to N <= 64 so
// that every row, column, g and f fits in one 64-bit word and nothing is heap-allocated.
//   forced entries:  a(n,j) = parity(col[j-1] & frev), j = 2..n-1, frev_k = f_{n-k} (k >= 1)
//   free entries:    a(n,0) = g_n,  a(n,1) = f_n + parity(g & frev)
//   transitivity:    row[j] must be a subset of the new row for every j in it
struct Node {
    struct Child { std::uint64_t r; bool a0, a1, fn; };
    struct Children {
        int k = 0;
        Child c[4];
        const Child* begin() const { return c; }
        const Child* end() const { return c + k; }
    };
    std::size_t n = 1;
    std::uint64_t row[64] = {1}, col[64] = {1}, g = 1, f = 0;

    std::size_t order() const { return n; }
    std::uint64_t row_word(std::size_t i) const { return row[i]; }
    std::uint64_t col_word(std::size_t i) const { return col[i]; }

    Children candidates() const {
        Children out;
        std::uint64_t frev = 0;
        for (std::size_t k = 1; k < n; ++k) frev |= ((f >> (n - k)) & 1) << k;
        const std::uint64_t diag = std::uint64_t{1} << n;
        std::uint64_t r = diag;
        for (std::size_t j = 2; j < n; ++j)
            if (std::popcount(col[j - 1] & frev) & 1) r |= std::uint64_t{1} << j;
        const bool conv = std::popcount(g & frev) & 1;
        for (std::uint64_t x = r & ~diag & ~std::uint64_t{3}; x; x &= x - 1)   // forced part closed?
            if (row[std::countr_zero(x)] & ~r & ~std::uint64_t{3}) return out;
        for (int a0 = 0; a0 < 2; ++a0)
            for (int a1 = (n >= 2 ? 0 : 1); a1 < 2; ++a1) {
                const std::uint64_t rr = r | std::uint64_t(a0) | (n >= 2 ? std::uint64_t(a1) << 1 : 0);
                bool closed = true;
                for (std::uint64_t x = rr & ~diag; x && closed; x &= x - 1)
                    closed = !(row[std::countr_zero(x)] & ~rr);
                if (closed) out.c[out.k++] = {rr, a0 != 0, a1 != 0, n >= 2 ? ((a1 != 0) != conv) : true};
            }
        return out;
    }
    void push(const Child& c) {
        const std::uint64_t bit = std::uint64_t{1} << n;
        row[n] = c.r;
        col[n] = 0;
        for (std::uint64_t x = c.r; x; x &= x - 1) col[std::countr_zero(x)] |= bit;
        if (c.a0) g |= bit;
        if (c.fn) f |= bit;
        ++n;
    }
    void pop() {
        --n;
        const std::uint64_t bit = std::uint64_t{1} << n;
        for (std::uint64_t x = row[n]; x; x &= x - 1) col[std::countr_zero(x)] &= ~bit;
        g &= ~bit;
        f &= ~bit;
    }
    void push_choice(bool a0, bool a1) {  // replay a recorded (a0, a1)
        for (const auto& c : candidates())
            if (c.a0 == a0 && c.a1 == a1) { push(c); return; }
        throw std::logic_error("rcount: invalid recorded path");
    }
    BoolMatrix matrix() const {
        BoolMatrix m(n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::uint64_t x = row[i]; x; x &= x - 1) m.set(i, std::countr_zero(x));
        return m;
    }
};

// ------------------------------------------------------------ canonical hashes --
static DEFAULTOPTIONS_DIGRAPH(g_nauty_options);

// nauty: digraph of the strict order (arc i -> j iff x_j < x_i), initial ordered partition by
// the invariant (down-set size, up-set size). The canonical graph together with the sorted
// invariant sequence identifies the isomorphism class.
static Hash128 class_hash_nauty(const Node& bld, std::size_t n) {
    graph g[MAXN], cg[MAXN];
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    std::uint64_t key[MAXN];
    for (std::size_t i = 0; i < n; ++i)
        key[i] = (std::uint64_t(std::popcount(bld.row_word(i))) << 8) | std::popcount(bld.col_word(i));
    for (std::size_t i = 0; i < n; ++i) lab[i] = static_cast<int>(i);
    std::sort(lab, lab + n, [&](int x, int y) { return key[x] != key[y] ? key[x] < key[y] : x < y; });
    for (std::size_t k = 0; k < n; ++k) ptn[k] = (k + 1 < n && key[lab[k + 1]] == key[lab[k]]) ? 1 : 0;
    EMPTYGRAPH(g, 1, static_cast<int>(n));
    for (std::size_t i = 0; i < n; ++i) {
        std::uint64_t r = bld.row_word(i) & ~(std::uint64_t{1} << i);
        while (r) {
            const int j = std::countr_zero(r);
            r &= r - 1;
            ADDONEARC(g, static_cast<int>(i), j, 1);
        }
    }
    optionblk opt = g_nauty_options;
    opt.getcanon = TRUE;
    opt.defaultptn = FALSE;
    statsblk stats;
    densenauty(g, lab, ptn, orbits, &opt, &stats, 1, static_cast<int>(n), cg);
    Hasher h;
    h.add(n);
    for (std::size_t k = 0; k < n; ++k) h.add(key[lab[k]]);
    for (std::size_t k = 0; k < n; ++k) h.add(static_cast<std::uint64_t>(cg[k]));
    return h.get();
}

// The library's own canonical form (independent of nauty), for cross-checks.
static Hash128 class_hash_library(const Node& bld, std::size_t) {
    const CanonicalForm cf = canonical_form(bld.matrix());
    Hasher h;
    h.add(cf.n);
    for (std::size_t k = 0; k < cf.cert.words(); ++k) h.add(cf.cert.word_at(k));
    return h.get();
}

// Isomorphism invariant used to split the work between passes: hash of the sorted
// multiset of (down-set size, up-set size).
static std::uint64_t hash_sorted_keys(std::uint64_t* key, std::size_t n) {
    std::sort(key, key + n);
    std::uint64_t h = 0x452821E638D01377ULL;
    for (std::size_t i = 0; i < n; ++i) h = mix64(h ^ (key[i] + 0x9E3779B97F4A7C15ULL * (i + 1)));
    return h;
}
static std::uint64_t invariant_hash(const Node& bld, std::size_t n) {
    std::uint64_t key[64];
    for (std::size_t i = 0; i < n; ++i)
        key[i] = (std::uint64_t(std::popcount(bld.row_word(i))) << 8) | std::popcount(bld.col_word(i));
    return hash_sorted_keys(key, n);
}

// ------------------------------------------------------------------ the count --
struct Options {
    std::size_t N = 0, threads = 1, split = 12, buckets = 256, passes = 0;
    bool use_nauty = true;
    std::string tmp = "rcount_tmp";
    std::string checkpoint;
};

using Path = std::vector<std::pair<bool, bool>>;

static void collect_jobs(Node& b, std::size_t depth, Path& path, std::vector<Path>& jobs) {
    if (b.order() == depth) { jobs.push_back(path); return; }
    for (const auto& c : b.candidates()) {
        b.push(c);
        path.emplace_back(c.a0, c.a1);
        collect_jobs(b, depth, path, jobs);
        path.pop_back();
        b.pop();
    }
}

class Worker {
public:
    Worker(const Options& o, std::size_t id) : o_(o), id_(id), buf_(o.buckets) {}
    template <class Canon>
    void run_job(const Path& p, Canon canon) {
        Node b;
        for (auto [a0, a1] : p) b.push_choice(a0, a1);
        dfs(b, canon);
    }
    void flush_all() { for (std::size_t k = 0; k < buf_.size(); ++k) flush(k); }
    std::uint64_t leaves() const { return leaves_; }
    void set_pass(std::size_t p) { pass_ = p; leaves_ = 0; }
    std::vector<Hash128>& memory() { return mem_; }
    static void compact(std::vector<Hash128>& v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

private:
    template <class Canon>
    void dfs(Node& b, Canon canon) {
        if (o_.passes && b.order() + 1 == o_.N) {
            // Last level in pass mode: compute each child's invariant from the parent's row and
            // column counts plus the candidate row, and build (push) only the children of this pass.
            const std::size_t n = o_.N, m = n - 1;
            std::uint64_t down[64], up[64], key[64];
            for (std::size_t i = 0; i < m; ++i) {
                down[i] = std::popcount(b.row_word(i));
                up[i] = std::popcount(b.col_word(i));
            }
            for (const auto& c : b.candidates()) {
                ++leaves_;
                const std::uint64_t r = c.r;
                for (std::size_t i = 0; i < m; ++i) key[i] = (down[i] << 8) | (up[i] + (r >> i & 1));
                key[m] = (std::uint64_t(std::popcount(r)) << 8) | 1;  // new element: only itself above
                if (hash_sorted_keys(key, n) % o_.passes != pass_) continue;
                b.push(c);
                mem_push(canon(b, o_.N));
                b.pop();
            }
            return;
        }
        if (b.order() == o_.N) {
            if (o_.passes) {  // memory mode
                ++leaves_;
                if (invariant_hash(b, o_.N) % o_.passes != pass_) return;
                mem_push(canon(b, o_.N));
                return;
            }
            ++leaves_;
            const Hash128 h = canon(b, o_.N);
            auto& v = buf_[h.a % o_.buckets];
            v.push_back(h);
            if (v.size() >= (1u << 15)) flush(h.a % o_.buckets);
            return;
        }
        for (const auto& c : b.candidates()) {
            b.push(c);
            dfs(b, canon);
            b.pop();
        }
    }
    void mem_push(const Hash128& h) {
        if (mem_.size() == mem_.capacity()) {  // compact, then grow by a bounded step
            compact(mem_);
            const std::size_t step = std::max<std::size_t>(std::size_t{1} << 22, mem_.size() / 4);
            std::vector<Hash128> next;
            next.reserve(mem_.size() + step);
            next.assign(mem_.begin(), mem_.end());
            mem_.swap(next);
        }
        mem_.push_back(h);
    }
    void flush(std::size_t k) {
        auto& v = buf_[k];
        if (v.empty()) return;
        const std::string name = o_.tmp + "/t" + std::to_string(id_) + "_b" + std::to_string(k) + ".bin";
        FILE* f = std::fopen(name.c_str(), "ab");
        if (!f || std::fwrite(v.data(), sizeof(Hash128), v.size(), f) != v.size())
            throw std::runtime_error("rcount: cannot write " + name);
        std::fclose(f);
        v.clear();
    }
    const Options& o_;
    std::size_t id_;
    std::vector<std::vector<Hash128>> buf_;
    std::uint64_t leaves_ = 0;
    std::size_t pass_ = 0;
    std::vector<Hash128> mem_;
};

static std::uint64_t count_bucket(const Options& o, std::size_t k) {
    std::vector<Hash128> all;
    for (std::size_t t = 0; t < o.threads; ++t) {
        const std::string name = o.tmp + "/t" + std::to_string(t) + "_b" + std::to_string(k) + ".bin";
        if (!fs::exists(name)) continue;
        const auto bytes = fs::file_size(name);
        const std::size_t old = all.size();
        all.resize(old + bytes / sizeof(Hash128));
        FILE* f = std::fopen(name.c_str(), "rb");
        if (!f || std::fread(all.data() + old, sizeof(Hash128), bytes / sizeof(Hash128), f) != bytes / sizeof(Hash128))
            throw std::runtime_error("rcount: cannot read " + name);
        std::fclose(f);
        fs::remove(name);
    }
    std::sort(all.begin(), all.end());
    return static_cast<std::uint64_t>(std::unique(all.begin(), all.end()) - all.begin());
}

int main(int argc, char** argv) {
    Options o;
    std::vector<std::string> a(argv + 1, argv + argc);
    if (a.empty()) {
        std::cerr << "usage: rcount N [--threads T] [--split D] [--canon nauty|library] [--buckets B] [--tmp DIR]\n";
        return 2;
    }
    o.N = std::stoul(a[0]);
    o.threads = std::max(1u, std::thread::hardware_concurrency());
    for (std::size_t i = 1; i + 1 < a.size(); i += 2) {
        if (a[i] == "--threads") o.threads = std::stoul(a[i + 1]);
        else if (a[i] == "--split") o.split = std::stoul(a[i + 1]);
        else if (a[i] == "--buckets") o.buckets = std::stoul(a[i + 1]);
        else if (a[i] == "--tmp") o.tmp = a[i + 1];
        else if (a[i] == "--canon") o.use_nauty = a[i + 1] != "library";
        else if (a[i] == "--passes") o.passes = std::stoul(a[i + 1]);
        else if (a[i] == "--checkpoint") o.checkpoint = a[i + 1];
        else { std::cerr << "unknown option " << a[i] << "\n"; return 2; }
    }
    if (o.N < 1 || o.N > 64) { std::cerr << "rcount: need 1 <= N <= 64\n"; return 2; }
    nauty_check(WORDSIZE, 1, static_cast<int>(o.N), NAUTYVERSIONID);
    o.split = std::min(o.split, o.N);
    fs::remove_all(o.tmp);
    fs::create_directories(o.tmp);

    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };

    std::vector<Path> jobs;
    {
        Node b;
        Path p;
        collect_jobs(b, o.split, p, jobs);
    }
    std::cerr << "rcount N=" << o.N << ": " << jobs.size() << " jobs at depth " << o.split << ", " << o.threads
              << " thread(s), canonical form: " << (o.use_nauty ? "nauty" : "library") << "\n";

    std::vector<Worker> workers;
    for (std::size_t t = 0; t < o.threads; ++t) workers.emplace_back(o, t);
    std::atomic<std::size_t> next{0}, done{0};
    std::uint64_t r_mem = 0, m_mem = 0;
    const std::size_t npass = o.passes ? o.passes : 1;
    std::vector<char> pass_done(npass, 0);
    if (o.passes && !o.checkpoint.empty() && fs::exists(o.checkpoint)) {
        FILE* f = std::fopen(o.checkpoint.c_str(), "r");
        unsigned long long cn, cp, cq, cr, cm;
        while (f && std::fscanf(f, "%llu %llu %llu %llu %llu", &cn, &cp, &cq, &cr, &cm) == 5) {
            if (cn != o.N || cp != o.passes || cq >= npass || pass_done[cq]) continue;
            if (m_mem && m_mem != cm) throw std::runtime_error("rcount: checkpoint lines disagree on m(N)");
            pass_done[cq] = 1; r_mem += cr; m_mem = cm;
            std::cerr << " resuming: pass " << cq + 1 << " already done (" << cr << " classes)\n";
        }
        if (f) std::fclose(f);
    }
    for (std::size_t pass = 0; pass < npass; ++pass) {
    if (pass_done[pass]) continue;
    next = 0;
    done = 0;
    for (auto& w : workers) w.set_pass(pass);
    if (o.passes) std::cerr << " pass " << pass + 1 << "/" << o.passes << ", " << elapsed() << " s\n";
    std::mutex io;
    auto work = [&](std::size_t t) {
        for (std::size_t j; (j = next.fetch_add(1)) < jobs.size();) {
            if (o.use_nauty) workers[t].run_job(jobs[j], class_hash_nauty);
            else workers[t].run_job(jobs[j], class_hash_library);
            const std::size_t d = ++done;
            if (d % std::max<std::size_t>(1, jobs.size() / 20) == 0) {
                std::lock_guard<std::mutex> lk(io);
                std::cerr << "  " << d << "/" << jobs.size() << " jobs, " << elapsed() << " s\n";
            }
        }
        workers[t].flush_all();
    };
    std::vector<std::thread> pool;
    for (std::size_t t = 1; t < o.threads; ++t) pool.emplace_back(work, t);
    work(0);
    for (auto& th : pool) th.join();
    if (o.passes) {  // merge this pass's classes and count them
        auto& all = workers[0].memory();
        for (std::size_t t = 1; t < workers.size(); ++t) {
            auto& v = workers[t].memory();
            all.insert(all.end(), v.begin(), v.end());
            std::vector<Hash128>().swap(v);
        }
        Worker::compact(all);
        std::uint64_t mp = 0;
        for (auto& w : workers) mp += w.leaves();
        if (m_mem && m_mem != mp) throw std::runtime_error("rcount: passes disagree on m(N)");
        m_mem = mp;
        r_mem += all.size();
        if (!o.checkpoint.empty()) {
            FILE* f = std::fopen(o.checkpoint.c_str(), "a");
            if (!f) throw std::runtime_error("rcount: cannot write checkpoint");
            std::fprintf(f, "%zu %zu %zu %llu %llu\n", o.N, o.passes, pass, (unsigned long long)all.size(),
                         (unsigned long long)mp);
            std::fclose(f);
        }
        std::vector<Hash128>().swap(all);
    }
    }

    std::uint64_t m = 0, r = 0;
    if (o.passes) m = m_mem;
    else for (auto& w : workers) m += w.leaves();
    const double t_enum = elapsed();
    if (o.passes) r = r_mem;
    else for (std::size_t k = 0; k < o.buckets; ++k) r += count_bucket(o, k);
    fs::remove_all(o.tmp);
    std::printf("%zu %llu %llu  (enumerate+canonicalise %.1f s, deduplicate %.1f s)\n", o.N, (unsigned long long)m,
                (unsigned long long)r, t_enum, elapsed() - t_enum);
    return 0;
}
