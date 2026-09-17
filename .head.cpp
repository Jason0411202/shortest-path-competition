// solver.cpp
//
// Derived from dijkstra_foundation.cpp for the Shortest-Path Competition.
// Same command line, same input formats, same output as the foundation:
//
//     ./solver <graph_file> <query_file> <output_file>
//
// The foundation's file-format documentation still applies verbatim; see
// dijkstra_foundation.cpp / README.md.  What changed:
//
//   1. I/O.  Whole files are read into one buffer and parsed by hand; answers
//      are formatted into one buffer and written once.
//
//   2. Contraction Hierarchies (Geisberger, Sanders, Schultes, Delling 2008).
//      Vertices are contracted one by one in order of a lazily updated
//      priority (edge quotient + hop quotient + level, as in RoutingKit),
//      adding shortcuts unless a bounded witness search proves them redundant.
//      Queries are bidirectional upward Dijkstra searches with
//      stall-on-demand, over a CSR renumbered by contraction rank.
//
// Everything is C++17 and the standard library only, in this one file.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <string>
#include <vector>

using std::vector;

static const int32_t FLAG_COORDS   = 1;
static const int32_t FLAG_DIRECTED = 2;

typedef uint64_t W;
static const W INF = (~0ull) >> 2;

static bool g_debug = false;
static std::chrono::steady_clock::time_point g_t0;
static void tlog(const char* what) {
    if (!g_debug) return;
    double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - g_t0).count();
    std::fprintf(stderr, "[%8.3f s] %s\n", s, what);
}

// ---------------------------------------------------------------- input ----

struct Scanner {
    vector<char> buf;
    const char* p = nullptr;
    const char* end = nullptr;

    void load(const char* path) {
        std::FILE* f = std::fopen(path, "rb");
        if (!f) { std::fprintf(stderr, "cannot open file: %s\n", path); std::exit(1); }
        std::fseek(f, 0, SEEK_END);
        long sz = std::ftell(f);
        std::fseek(f, 0, SEEK_SET);
        if (sz < 0) sz = 0;
        buf.resize((size_t)sz + 1);
        size_t got = std::fread(buf.data(), 1, (size_t)sz, f);
        buf[got] = '\0';
        std::fclose(f);
        p = buf.data();
        end = buf.data() + got;
    }

    void release() { vector<char>().swap(buf); p = end = nullptr; }

    inline uint64_t uint_() {
        while (p < end && (*p < '0' || *p > '9')) ++p;
        uint64_t x = 0;
        while (p < end && *p >= '0' && *p <= '9') { x = x * 10 + (uint64_t)(*p - '0'); ++p; }
        return x;
    }
};

// ---------------------------------------------------------------- graph ----

static int32_t V = 0, E = 0, FLAGS = 0;
static bool directed = false;

// Dynamic graph used during contraction.  Each vertex keeps one list of
// neighbours; fw is the weight of v -> to, bw the weight of to -> v (INF when
// that direction does not exist).  Undirected graphs have fw == bw.
struct DArc { int32_t to; int32_t hops; W fw; W bw; };
static vector<vector<DArc>> g;

static void add_arc(int32_t u, int32_t v, W fw, W bw, int32_t hops) {
    for (DArc& a : g[u]) {
        if (a.to == v) {
            if (fw < a.fw) a.fw = fw;
            if (bw < a.bw) a.bw = bw;
            if (hops > a.hops) a.hops = hops;
            return;
        }
    }
    g[u].push_back({v, hops, fw, bw});
}

static void read_graph(const char* path) {
    Scanner sc;
    sc.load(path);
    V = (int32_t)sc.uint_();
    E = (int32_t)sc.uint_();
    FLAGS = (int32_t)sc.uint_();
    directed = (FLAGS & FLAG_DIRECTED) != 0;
    g.assign((size_t)V, {});
    for (int32_t i = 0; i < E; ++i) {
        int32_t u = (int32_t)sc.uint_();
        int32_t v = (int32_t)sc.uint_();
        W w = (W)sc.uint_();
        if (u == v) continue;
        if (directed) { add_arc(u, v, w, INF, 1); add_arc(v, u, INF, w, 1); }
        else          { add_arc(u, v, w, w, 1);   add_arc(v, u, w, w, 1); }
    }
    // Coordinates (if any) are not needed.
    sc.release();
}

// -------------------------------------------------------------- witness ----

struct MinHeap {
    struct Item { W d; int32_t v; };
    vector<Item> a;
    inline void clear() { a.clear(); }
    inline bool empty() const { return a.empty(); }
    inline W top_key() const { return a[0].d; }
    inline void push(W d, int32_t v) {
        size_t i = a.size();
        a.push_back({d, v});
        while (i) {
            size_t par = (i - 1) >> 2;
            if (a[par].d <= d) break;
            a[i] = a[par];
            i = par;
        }
        a[i] = {d, v};
    }
    inline Item pop() {
        Item best = a[0];
        Item last = a.back();
        a.pop_back();
        if (a.empty()) return best;
        size_t n = a.size(), i = 0;
        for (;;) {
            size_t c = (i << 2) + 1;
            if (c >= n) break;
            size_t lim = c + 4 < n ? c + 4 : n;
            size_t m = c;
            for (size_t j = c + 1; j < lim; ++j) if (a[j].d < a[m].d) m = j;
            if (a[m].d >= last.d) break;
            a[i] = a[m];
            i = m;
        }
        a[i] = last;
        return best;
    }
};

static int SIM_SETTLE = 60, CON_SETTLE = 1000;
static uint64_t g_settles = 0, g_searches = 0;
static vector<W> wdist;
static vector<uint32_t> wstamp;
static uint32_t wcur = 0;
static MinHeap wheap;

static inline W wget(int32_t v) { return wstamp[v] == wcur ? wdist[v] : INF; }

static vector<uint32_t> tstamp;   // target marks for the current search
static uint32_t tcur = 0;

// Forward Dijkstra from src in the remaining graph, skipping `avoid`, until the
// key exceeds `limit`, `max_settle` vertices are settled, or all `ntargets`
// vertices marked with tcur are settled.
static void witness_search(int32_t src, int32_t avoid, W limit, int max_settle, int ntargets) {
    ++g_searches;
    ++wcur;
    if (wcur == 0) { std::fill(wstamp.begin(), wstamp.end(), 0); wcur = 1; }
    wheap.clear();
    wdist[src] = 0; wstamp[src] = wcur;
    wheap.push(0, src);
    int settled = 0;
    while (!wheap.empty()) {
        MinHeap::Item it = wheap.pop();
        if (it.d > wdist[it.v]) continue;
        if (it.d > limit) break;
        if (tstamp[it.v] == tcur && --ntargets <= 0) break;
        if (++settled > max_settle) break;
        ++g_settles;

        for (const DArc& a : g[it.v]) {
            if (a.fw >= INF || a.to == avoid) continue;
            W nd = it.d + a.fw;
            if (nd > limit) continue;
            if (wstamp[a.to] != wcur || nd < wdist[a.to]) {
                wstamp[a.to] = wcur;
                wdist[a.to] = nd;
                wheap.push(nd, a.to);
            }
        }
    }
}

static inline void new_targets() {
    ++tcur;
    if (tcur == 0) { std::fill(tstamp.begin(), tstamp.end(), 0); tcur = 1; }
}

