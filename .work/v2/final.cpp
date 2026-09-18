// solver.cpp
//
// Derived from dijkstra_foundation.cpp for the Shortest-Path Competition.
// Same command line, same input formats, same output as the foundation:
//
//     ./solver <graph_file> <query_file> <output_file>
//
// The foundation's file-format documentation still applies verbatim; see
// dijkstra_foundation.cpp / README.md.  The foundation answers every query
// with its own Dijkstra search.  This file keeps its read_graph /
// run_queries / main structure, reads both input files in one go, and picks a
// preprocessing scheme from the shape of the graph it was given:
//
//   * undirected torus lattices (and any other undirected graph): contraction
//     hierarchies (Geisberger et al. 2008) with landmark-bounded, goal-directed
//     witness searches; the last few thousand vertices form a core whose
//     all-pairs distances are tabulated; queries are bidirectional upward
//     searches that meet either below the core or through the table;
//   * directed graphs with coordinates (road networks): see namespace road;
//   * directed graphs without coordinates (power-law graphs): see namespace sf.
//
// Every scheme is exact on any input of its kind; structure detection only
// selects faster code paths.  Single-threaded, C++17, standard library only.

// ===================================================================== common
// Shared I/O layer.  Pasted verbatim at the top of every family prototype and
// of the final solver.cpp (no #include of own headers in the final file).
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <vector>
#include <type_traits>
#include <cmath>

using std::vector;

static const int32_t FLAG_COORDS   = 1;
static const int32_t FLAG_DIRECTED = 2;

static bool g_debug = false;
static std::chrono::steady_clock::time_point g_t0;
static void tlog(const char* what) {
    if (!g_debug) return;
    double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - g_t0).count();
    std::fprintf(stderr, "[%8.3f s] %s\n", s, what);
}

// Whole file in one buffer, NUL-terminated (the parser relies on a non-digit
// sentinel after the last number).
static char* slurp(const char* path, size_t& n) {
    std::FILE* f = std::fopen(path, "rb");
    if (!f) { std::fprintf(stderr, "cannot open file: %s\n", path); std::exit(1); }
    std::fseek(f, 0, SEEK_END);
    long sz = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    if (sz < 0) sz = 0;
    char* b = (char*)std::malloc((size_t)sz + 16);
    n = std::fread(b, 1, (size_t)sz, f);
    std::memset(b + n, 0, 16);
    std::fclose(f);
    return b;
}

// Reads the next non-negative integer; skips any non-digit characters first.
static inline uint64_t rd_u(const char*& p) {
    while ((unsigned)(*p - '0') > 9u) { if (!*p) return 0; ++p; }
    uint64_t x = (uint64_t)(*p++ - '0');
    unsigned c;
    while ((c = (unsigned)(*p - '0')) <= 9u) { x = x * 10 + c; ++p; }
    return x;
}

static int32_t V = 0, E = 0, FLAGS = 0;
static bool directed = false, has_coords = false;
static vector<int32_t> eu, ev;     // edge endpoints as listed in the file
static vector<uint32_t> ew;        // weights (positive, fit in int32)
static vector<int64_t> cx, cy;     // coordinates (only if FLAGS bit 0)
static int32_t Q = 0;
static vector<int32_t> qs, qt;
static vector<int64_t> qans;       // -1 = unreachable

static void read_graph(const char* path) {
    size_t n;
    char* buf = slurp(path, n);
    const char* p = buf;
    V = (int32_t)rd_u(p); E = (int32_t)rd_u(p); FLAGS = (int32_t)rd_u(p);
    directed = (FLAGS & FLAG_DIRECTED) != 0;
    has_coords = (FLAGS & FLAG_COORDS) != 0;
    eu.resize((size_t)E); ev.resize((size_t)E); ew.resize((size_t)E);
    for (int32_t i = 0; i < E; ++i) {
        eu[i] = (int32_t)rd_u(p); ev[i] = (int32_t)rd_u(p); ew[i] = (uint32_t)rd_u(p);
    }
    if (has_coords) {
        cx.resize((size_t)V); cy.resize((size_t)V);
        for (int32_t i = 0; i < V; ++i) { cx[i] = (int64_t)rd_u(p); cy[i] = (int64_t)rd_u(p); }
    }
    std::free(buf);
}

static void read_queries(const char* path) {
    size_t n;
    char* buf = slurp(path, n);
    const char* p = buf;
    Q = (int32_t)rd_u(p);
    qs.resize((size_t)Q); qt.resize((size_t)Q); qans.assign((size_t)Q, -1);
    for (int32_t i = 0; i < Q; ++i) { qs[i] = (int32_t)rd_u(p); qt[i] = (int32_t)rd_u(p); }
    std::free(buf);
}

static void write_answers(const char* path) {
    char* out = (char*)std::malloc((size_t)Q * 21 + 64);
    char* o = out;
    for (int32_t i = 0; i < Q; ++i) {
        int64_t d = qans[i];
        if (d < 0) { *o++ = '-'; *o++ = '1'; }
        else {
            char tmp[24]; int k = 0; uint64_t x = (uint64_t)d;
            do { tmp[k++] = (char)('0' + x % 10); x /= 10; } while (x);
            while (k) *o++ = tmp[--k];
        }
        *o++ = '\n';
    }
    std::FILE* f = std::fopen(path, "wb");
    if (!f) { std::fprintf(stderr, "cannot open output file: %s\n", path); std::exit(1); }
    std::fwrite(out, 1, (size_t)(o - out), f);
    std::fclose(f);
}

// 4-ary min-heap of (key, vertex) with lazy deletion (callers skip stale pops).
template <class K>
struct Heap4 {
    struct Item { K d; int32_t v; };
    vector<Item> a;
    inline void clear() { a.clear(); }
    inline bool empty() const { return a.empty(); }
    inline K top_key() const { return a[0].d; }
    inline void push(K d, int32_t v) {
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
        size_t n = a.size();
        if (!n) return best;
        size_t i = 0;
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
// =================================================================== /common


#pragma GCC push_options
#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")
#include <type_traits>
#pragma GCC push_options
#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")

// ============================================================ wide64
// Undirected graphs with 64-bit distances (the 2-D torus lattice family).
//
// 1. Contraction hierarchy (greedy priorities, contracted in rounds of
//    independent sets in Morton order) down to a core of C vertices.
// 2. Queries: bidirectional upward searches in the hierarchy (complete, so
//    the core entry distances are exact), then a bidirectional ALT search
//    inside the core (landmarks computed on the core graph), seeded with the
//    entries.  s-t paths that never reach the core are found by the upward
//    searches meeting.
namespace w64 {

typedef uint64_t D;
static const D INF = (~0ull) >> 2;
static int64_t envi(const char* k, int64_t def) { const char* e = std::getenv(k); return e ? std::atoll(e) : def; }

// ------------------------------------------------------------ lattice detect
static int dim = 0; static int64_t side = 0; static int64_t strd[3];
static bool detect() {
    if (directed || has_coords) return false;
    for (int dd = 2; dd <= 3; ++dd) {
        int64_t s = (int64_t)(std::pow((double)V, 1.0 / dd) + 0.5);
        for (int64_t c = std::max<int64_t>(3, s - 1); c <= s + 1; ++c) {
            int64_t p = 1; for (int k = 0; k < dd; ++k) p *= c;
            if (p == V && (int64_t)E == dd * (int64_t)V) { dim = dd; side = c; }
        }
    }
    if (!dim) return false;
    strd[0] = 1; strd[1] = side; strd[2] = side * side;
    for (int64_t i = 0; i < V; ++i) for (int k = 0; k < dim; ++k) {
        int64_t e = i * dim + k; int64_t c = (i / strd[k]) % side;
        int64_t j = c + 1 < side ? i + strd[k] : i - c * strd[k];
        if (eu[e] != i || ev[e] != j) { dim = 0; return false; }
    }
    return true;
}

// Dispatch test: undirected 2-D torus lattice (edges in generator order)
// whose distances may exceed 32 bits.  solve() itself is correct for any
// undirected graph.
bool applicable() {
    if (directed || has_coords) return false;
    if (!detect() || dim != 2) return false;
    uint64_t maxw = 0;
    for (int32_t e = 0; e < E; ++e) if (ew[e] > maxw) maxw = ew[e];
    return maxw * (uint64_t)(dim * (side / 2) + 1) >= (1ull << 32);
}

// ------------------------------------------------------------ arc formats
typedef double v4d __attribute__((vector_size(32)));
struct A16 { int32_t t; uint32_t pad; D w; };
static inline int32_t ato(const A16& a) { return a.t; }
static inline D aw(const A16& a) { return a.w; }
static inline A16 amk16(int32_t t, D w) { return {t, 0, w}; }
// packed: weight << 22 | target  (needs V < 2^22 and every useful weight < 2^42)
struct A8 { uint64_t x; };
static inline int32_t ato(const A8& a) { return (int32_t)(a.x & 0x3FFFFFu); }
static inline D aw(const A8& a) { return a.x >> 22; }
template <class A> static inline A amk(int32_t t, D w);
template <> inline A16 amk<A16>(int32_t t, D w) { return {t, 0, w}; }
template <> inline A8 amk<A8>(int32_t t, D w) { return {(w << 22) | (uint64_t)(uint32_t)t}; }

// 4-ary min-heap of packed (key << 22 | vertex), padded with ~0 so the
// four children of a node are compared without bounds checks.  Keys < 2^42.
struct PH {
    struct Item { D d; int32_t v; };
    vector<uint64_t> a; uint32_t n = 0;
    PH() { a.assign(256, ~0ull); }
    inline void clear() { for (uint32_t i = 0; i < n; ++i) a[i] = ~0ull; n = 0; }
    inline bool empty() const { return n == 0; }
    inline D top_key() const { return a[0] >> 22; }
    inline void push(D d, int32_t v) {
        if (n + 8 >= a.size()) a.resize(a.size() * 2, ~0ull);
        uint64_t x = (d << 22) | (uint64_t)(uint32_t)v;
        uint32_t i = n++;
        while (i) { uint32_t p = (i - 1) >> 2; if (a[p] <= x) break; a[i] = a[p]; i = p; }
        a[i] = x;
    }
    inline Item pop() {
        uint64_t best = a[0];
        uint64_t x = a[--n];
        a[n] = ~0ull;
        uint32_t i = 0;
        for (;;) {
            uint32_t c = 4 * i + 1;
            if (c >= n) break;
            uint64_t m0 = a[c], m1 = a[c + 1], m2 = a[c + 2], m3 = a[c + 3];
            uint32_t j0 = m1 < m0 ? c + 1 : c;       uint64_t v0 = m1 < m0 ? m1 : m0;
            uint32_t j1 = m3 < m2 ? c + 3 : c + 2;   uint64_t v1 = m3 < m2 ? m3 : m2;
            uint32_t j = v1 < v0 ? j1 : j0;          uint64_t vv = v1 < v0 ? v1 : v0;
            if (vv >= x) break;
            a[i] = vv; i = j;
        }
        if (n) a[i] = x;
        return {best >> 22, (int32_t)(best & 0x3FFFFFu)};
    }
};

// ------------------------------------------------------------ CH engine
template <class A>
struct Engine {
    int32_t m = 0, alive = 0;
    vector<A> pool;
    struct VR { uint32_t off, deg, cap, ms; D md; };        // adjacency + 2-hop marks
    struct WS { D d, tcost, tlim; uint32_t stamp, tm; };   // witness search state
    vector<VR> vr; vector<WS> ws;
    vector<int32_t> orig, level; vector<int64_t> pr; vector<uint8_t> dirty, dead;
    uint32_t wc = 0, tc = 0, mc = 0;
    typename std::conditional<std::is_same<A, A8>::value, PH, Heap4<D>>::type heap;
    D UB = INF;                        // arcs longer than this are never needed
    // output hierarchy (vertex ids)
    vector<int32_t> rorder;
    vector<A16> uarc; vector<uint32_t> uoff;
    vector<int32_t> vlevel;
    // knobs / stats
    int CON_SET = 200; int64_t LC = 250; int32_t STOPC = 120000;
    uint64_t searches = 0, settles = 0, sims = 0, lbneed = 0;
    // landmark distances by vertex id (computed on the current graph; the
    // contraction preserves distances between remaining vertices)
    int LK = 0; vector<D> LDv;
    inline bool lb_needed(int32_t u, int32_t x, D c) const {
        const D* a = &LDv[(size_t)orig[u] * LK]; const D* b = &LDv[(size_t)orig[x] * LK];
        for (int l = 0; l < LK; ++l) { D d = a[l] > b[l] ? a[l] - b[l] : b[l] - a[l]; if (d >= c) return true; }
        return false;
    }
    void landmarks_now(int k) {
        LK = k; if (!LK) return;
        int32_t nV = (int32_t)vlevel.size();
        LDv.assign((size_t)nV * LK, 0);
        vector<D> d(m), mind(m, INF);
        PH h;
        int32_t src = 0;
        for (int l = 0; l < LK; ++l) {
            std::fill(d.begin(), d.end(), INF);
            d[src] = 0; h.clear(); h.push(0, src);
            while (!h.empty()) {
                auto it = h.pop(); if (it.d > d[it.v]) continue;
                const A* a = nb(it.v); uint32_t kk = vr[it.v].deg;
                for (uint32_t q = 0; q < kk; ++q) { int32_t y = ato(a[q]); D nd = it.d + aw(a[q]); if (nd < d[y]) { d[y] = nd; h.push(nd, y); } }
            }
            int32_t far = src; D fd = 0;
            for (int32_t v = 0; v < m; ++v) {
                if (dead[v]) continue;
                if (d[v] >= INF) { LK = 0; LDv.clear(); return; }
                LDv[(size_t)orig[v] * LK + l] = d[v];
                if (d[v] < mind[v]) mind[v] = d[v];
                if (mind[v] > fd) { fd = mind[v]; far = v; }
            }
            src = far;
        }
    }

    inline A* nb(int32_t v) { return &pool[vr[v].off]; }
    void alloc_slots(int32_t n) {
        m = n; vr.assign(n, VR{0, 0, 0, 0, INF}); ws.assign(n, WS{INF, 0, 0, 0, 0});
        wc = tc = mc = 0;
    }
    inline void push_arc(int32_t u, int32_t v, D w) {
        VR& r = vr[u];
        if (r.deg == r.cap) {
            uint32_t nc = r.cap * 2 + 4;
            uint32_t no = (uint32_t)pool.size();
            pool.resize(pool.size() + nc);
            if (r.deg) std::memcpy(&pool[no], &pool[r.off], sizeof(A) * r.deg);
            r.off = no; r.cap = nc;
        }
        pool[r.off + r.deg++] = amk<A>(v, w);
    }
    void add_edge(int32_t u, int32_t v, D w) {
        if (w > UB) return;
        A* a = nb(u); uint32_t k = vr[u].deg;
        for (uint32_t i = 0; i < k; ++i) if (ato(a[i]) == v) {
            if (w < aw(a[i])) {
                a[i] = amk<A>(v, w);
                A* b = nb(v); uint32_t kb = vr[v].deg;
                for (uint32_t j = 0; j < kb; ++j) if (ato(b[j]) == u) { b[j] = amk<A>(u, w); break; }
            }
            return;
        }
        push_arc(u, v, w); push_arc(v, u, w);
    }
    inline void remove_arc(int32_t u, int32_t v) {
        A* a = nb(u); uint32_t k = vr[u].deg;
        for (uint32_t i = 0; i < k; ++i) if (ato(a[i]) == v) { a[i] = a[k - 1]; --vr[u].deg; return; }
    }

    // --- witness search from src avoiding `avoid`, targets tl (tm == tc)
    vector<int32_t> tl;
    __attribute__((noinline)) void witness(int32_t src, int32_t avoid, int pending) {
        ++searches;
        if (++wc == 0) { for (WS& s : ws) s.stamp = 0; wc = 1; }
        D limit = 0, tmax = 0;
        for (int32_t x : tl) if (ws[x].tm == tc) { limit = std::max(limit, ws[x].tlim); tmax = std::max(tmax, ws[x].tcost); }
        heap.clear();
        ws[src].d = 0; ws[src].stamp = wc; heap.push(0, src);
        int settled = 0;
        while (!heap.empty()) {
            auto it = heap.pop();
            if (it.d > ws[it.v].d) continue;
            if (it.d > limit) break;
            if (++settled > CON_SET) break;
            ++settles;
            bool relim = false;
            const A* a = nb(it.v); uint32_t k = vr[it.v].deg;
            for (uint32_t i = 0; i < k; ++i) {
                int32_t y = ato(a[i]);
                if (y == avoid) continue;
                D nd = it.d + aw(a[i]);
                if (nd > tmax) continue;           // neither a target hit nor within the limit
                WS& sy = ws[y];
                if (sy.tm == tc && nd <= sy.tcost) {
                    sy.tm = 0;
                    if (--pending <= 0) return;
                    if (sy.tlim >= limit) relim = true;
                }
                if (nd > limit) continue;
                if (sy.stamp != wc || nd < sy.d) { sy.stamp = wc; sy.d = nd; heap.push(nd, y); }
            }
            if (relim) { limit = 0; for (int32_t x : tl) if (ws[x].tm == tc) limit = std::max(limit, ws[x].tlim); }
        }
    }
    // mark u's neighbours (except v) for 2-hop checks
    inline void mark(int32_t u, int32_t v) {
        if (++mc == 0) { for (VR& r : vr) r.ms = 0; mc = 1; }
        const A* au = nb(u); uint32_t ku = vr[u].deg;
        for (uint32_t q = 0; q < ku; ++q) { int32_t b = ato(au[q]); if (b != v) { vr[b].ms = mc; vr[b].md = aw(au[q]); } }
    }
    // is there a 1- or 2-arc u-x path avoiding v of length <= c (u marked)?
    inline bool hop2(int32_t x, int32_t v, D c) const {
        const VR& rx = vr[x];
        if (rx.ms == mc && rx.md <= c) return true;
        const A* ax = &pool[rx.off]; uint32_t kx = rx.deg;
        for (uint32_t q = 0; q < kx; ++q) {
            int32_t b = ato(ax[q]);
            if (b == v) continue;
            const VR& rb = vr[b];
            if (rb.ms == mc && rb.md + aw(ax[q]) <= c) return true;
        }
        return false;
    }
    struct SC { int32_t a, b; D w; };
    vector<D> mwv; vector<SC> scs;
    vector<uint8_t> hw; int32_t simv = -1;   // 2-hop result of the last simulate()
    __attribute__((noinline)) void shortcuts(int32_t v) {
        scs.clear();
        const A* nv = nb(v); uint32_t k = vr[v].deg;
        mwv.resize(k);
        for (uint32_t j = 0; j < k; ++j) {
            int32_t x = ato(nv[j]);
            D mn = INF; const A* a = nb(x); uint32_t kk = vr[x].deg;
            for (uint32_t q = 0; q < kk; ++q) if (ato(a[q]) != v && aw(a[q]) < mn) mn = aw(a[q]);
            mwv[j] = mn;
        }
        for (uint32_t i = 0; i + 1 < k; ++i) {
            if (++tc == 0) { for (WS& s : ws) s.tm = 0; tc = 1; }
            tl.clear();
            int pending = 0;
            int32_t u = ato(nv[i]); D wu = aw(nv[i]);
            const bool reuse = simv == v;         // simulate(v) just ran on this graph
            if (!reuse) mark(u, v);
            for (uint32_t j = i + 1; j < k; ++j) {
                int32_t x = ato(nv[j]); D c = wu + aw(nv[j]);
                WS& sx = ws[x];
                sx.tm = tc;
                if (reuse ? hw[i * k + j] != 0 : (c > UB || hop2(x, v, c))) { sx.tm = 0; continue; }
                if (mwv[j] > c || mwv[i] > c) { sx.tcost = 0; continue; }   // no witness possible
                if (LK && lb_needed(u, x, c)) { sx.tcost = 0; ++lbneed; continue; }
                sx.tcost = c; sx.tlim = c - mwv[j]; tl.push_back(x); ++pending;
            }
            if (pending) witness(u, v, pending);
            for (uint32_t j = i + 1; j < k; ++j) {
                int32_t x = ato(nv[j]);
                if (ws[x].tm == tc) scs.push_back({u, x, wu + aw(nv[j])});
            }
        }
    }
    // number of shortcuts contracting v would add (1- and 2-hop witnesses only)
    __attribute__((noinline)) int simulate(int32_t v) {
        ++sims;
        const A* nv = nb(v); uint32_t k = vr[v].deg;
        int added = 0;
        hw.resize((size_t)k * k);
        for (uint32_t i = 0; i + 1 < k; ++i) {
            int32_t u = ato(nv[i]); D wu = aw(nv[i]);
            mark(u, v);
            for (uint32_t j = i + 1; j < k; ++j) {
                D c = wu + aw(nv[j]);
                bool w = c > UB || hop2(ato(nv[j]), v, c);
                hw[i * k + j] = w;
                added += !w;
            }
        }
        simv = v;
        return added;
    }
    int LMK_N = 0, LMK_AT = 0;
    int LAZY = 0, LAZYT = 1000000; vector<int32_t> plevel, sdeg;    // level / degree used in pr (lazy mode)
    inline int64_t prio(int32_t v) {
        int deg = (int)vr[v].deg;
        if (!plevel.empty()) { plevel[v] = level[v]; sdeg[v] = vr[v].deg; }
        int add = simulate(v);
        return (int64_t)level[v] * LC + (1000 * (int64_t)add) / std::max(1, deg);
    }
    void contract_one(int32_t v) {
        shortcuts(v);
        simv = -1;
        rorder.push_back(orig[v]);
        vlevel[orig[v]] = level[v];
        const A* nv = nb(v); uint32_t k = vr[v].deg;
        for (uint32_t i = 0; i < k; ++i) { int32_t u = ato(nv[i]); uarc.push_back(amk16(orig[u], aw(nv[i]))); remove_arc(u, v); }
        uoff.push_back((uint32_t)uarc.size());
        vr[v].deg = 0; dead[v] = 1; --alive;
        for (const SC& s : scs) add_edge(s.a, s.b, s.w);
    }
    // renumber alive slots 0..alive-1 (order kept) and rebuild storage
    void compact() {
        simv = -1;
        vector<int32_t> ns(m, -1);
        int32_t n2 = 0;
        for (int32_t v = 0; v < m; ++v) if (!dead[v]) ns[v] = n2++;
        size_t tot = 0;
        for (int32_t v = 0; v < m; ++v) if (!dead[v]) tot += vr[v].deg + vr[v].deg / 2 + 4;
        vector<A> np(tot);
        vector<VR> vr2(n2, VR{0, 0, 0, 0, INF});
        vector<int32_t> or2(n2), lv2(n2), pl2(n2), sd2(n2); vector<int64_t> pr2(n2); vector<uint8_t> di2(n2);
        size_t pos = 0;
        for (int32_t v = 0; v < m; ++v) if (!dead[v]) {
            int32_t w = ns[v];
            uint32_t d = vr[v].deg;
            vr2[w].off = (uint32_t)pos; vr2[w].deg = d; vr2[w].cap = d + d / 2 + 4;
            const A* a = nb(v);
            for (uint32_t i = 0; i < d; ++i) np[pos + i] = amk<A>(ns[ato(a[i])], aw(a[i]));
            pos += vr2[w].cap;
            or2[w] = orig[v]; lv2[w] = level[v]; pr2[w] = pr[v]; di2[w] = dirty[v]; pl2[w] = plevel.empty() ? 0 : plevel[v]; sd2[w] = sdeg.empty() ? 0 : sdeg[v];
        }
        pool.swap(np);
        alloc_slots(n2);
        vr.swap(vr2);
        orig.swap(or2); level.swap(lv2); pr.swap(pr2); dirty.swap(di2);
        if (!plevel.empty()) { plevel.swap(pl2); sdeg.swap(sd2); }
        dead.assign(n2, 0);
        alive = n2;
    }
    void debug_state(const char* what) {
        if (!g_debug) return;
        size_t arcs = 0; for (int32_t w = 0; w < m; ++w) if (!dead[w]) arcs += vr[w].deg;
        std::fprintf(stderr, "  %s: alive=%d lbneed=%llu avgdeg=%.2f searches=%llu settles=%llu sims=%llu\n", what, alive,
                     (unsigned long long)lbneed, (double)arcs / std::max(1, alive), (unsigned long long)searches, (unsigned long long)settles, (unsigned long long)sims);
        tlog(what);
    }
    // Rounds of independent sets: every vertex whose priority is a strict
    // local minimum (ties by slot) is contracted, in slot (= spatial) order.
    void contract_rounds(int32_t until) {
        vector<int32_t> cand, touched;
        // bitmap of vertices whose local-minimum status may have changed
        vector<uint64_t> chk(((size_t)m + 63) / 64, 0);
        vector<int32_t> tstamp(m, 0);
        int32_t next_compact = alive / 2;
        int rounds = 0;
        for (int32_t v = 0; v < m; ++v) if (!dead[v]) chk[v >> 6] |= 1ull << (v & 63);
        auto is_min = [&](int32_t v) {
            const A* a = nb(v); uint32_t k = vr[v].deg;
            int64_t p = pr[v];
            for (uint32_t i = 0; i < k; ++i) { int32_t u = ato(a[i]); if (pr[u] < p || (pr[u] == p && u < v)) return false; }
            return true;
        };
        auto add_check = [&](int32_t v) { chk[v >> 6] |= 1ull << (v & 63); };
        bool any = true;
        while (alive > until && any) {
            ++rounds;
            cand.clear();
            for (size_t w = 0; w < chk.size(); ++w) {
                uint64_t b = chk[w]; chk[w] = 0;
                while (b) {
                    int32_t v = (int32_t)(w * 64 + (size_t)__builtin_ctzll(b)); b &= b - 1;
                    if (!dead[v] && is_min(v)) cand.push_back(v);
                }
            }
            any = false;
            touched.clear();
            for (int32_t v : cand) {
                if (alive <= until) break;
                if (LAZY && dirty[v]) {
                    // stale quotient: refresh and re-check against the neighbours
                    dirty[v] = 0;
                    int64_t old = pr[v];
                    pr[v] = prio(v);
                    if (!is_min(v)) {
                        add_check(v); any = true;
                        if (pr[v] > old) { const A* a = nb(v); for (uint32_t i = 0; i < vr[v].deg; ++i) add_check(ato(a[i])); }
                        continue;
                    }
                }
                const A* nv = nb(v); uint32_t k = vr[v].deg;
                for (uint32_t i = 0; i < k; ++i) { int32_t u = ato(nv[i]); level[u] = std::max(level[u], level[v] + 1); dirty[u] = 1; if (tstamp[u] != rounds) { tstamp[u] = rounds; touched.push_back(u); } }
                contract_one(v);
            }
            if (LAZY) {
                // keep the (stale) quotient, update the level term only
                for (int32_t u : touched) if (!dead[u]) {
                    int dd = (int)vr[u].deg - (int)sdeg[u];
                    if (dd >= LAZYT || -dd >= LAZYT) { dirty[u] = 0; pr[u] = prio(u); }
                    else pr[u] = (int64_t)level[u] * LC + (pr[u] - (int64_t)plevel[u] * LC), plevel[u] = level[u];
                }
            } else {
                std::sort(touched.begin(), touched.end());
                for (int32_t u : touched) { dirty[u] = 0; if (!dead[u]) pr[u] = prio(u); }
            }
            for (int32_t u : touched) if (!dead[u]) { any = true; add_check(u); const A* a = nb(u); for (uint32_t i = 0; i < vr[u].deg; ++i) add_check(ato(a[i])); }
            if (alive > until && !any) {   // safety: full rescan
                for (int32_t v = 0; v < m; ++v) if (!dead[v]) add_check(v);
                any = !cand.empty();
            }
            if (alive <= next_compact && alive > until) {
                vector<int32_t> ns(m, -1); int32_t c2 = 0;
                for (int32_t v = 0; v < m; ++v) if (!dead[v]) ns[v] = c2++;
                vector<int32_t> keep;
                for (size_t w = 0; w < chk.size(); ++w) for (uint64_t b = chk[w]; b; b &= b - 1) { int32_t v = (int32_t)(w * 64 + (size_t)__builtin_ctzll(b)); if (ns[v] >= 0) keep.push_back(ns[v]); }
                compact();
                next_compact = alive / 2;
                chk.assign(((size_t)m + 63) / 64, 0);
                for (int32_t v : keep) add_check(v);
                tstamp.assign(m, 0);
                if (!LK && LMK_N > 0 && LMK_AT >= alive) { landmarks_now(LMK_N); tlog("  engine landmarks"); }
                debug_state("compact");
            }
        }
        if (g_debug) std::fprintf(stderr, "  rounds=%d\n", rounds);
    }
    void contract_greedy(int32_t until) {
        typedef std::pair<int64_t, int32_t> P;
        vector<P> init; init.reserve(alive);
        for (int32_t v = 0; v < m; ++v) if (!dead[v]) init.push_back({pr[v], v});
        std::priority_queue<P, vector<P>, std::greater<P>> pq(std::greater<P>(), std::move(init));
        int32_t next_compact = alive / 2;
        while (!pq.empty() && alive > until) {
            P t = pq.top(); pq.pop();
            int32_t v = t.second;
            if (dead[v] || t.first != pr[v]) continue;
            if (dirty[v]) {
                dirty[v] = 0;
                int64_t np = prio(v);
                pr[v] = np;
                if (!pq.empty() && np > pq.top().first) { pq.push({np, v}); continue; }
            }
            const A* nv = nb(v); uint32_t k = vr[v].deg;
            for (uint32_t i = 0; i < k; ++i) { int32_t u = ato(nv[i]); level[u] = std::max(level[u], level[v] + 1); dirty[u] = 1; }
            contract_one(v);
            if (alive <= next_compact && alive > until && alive >= 1000) {
                compact();
                vector<P> re; re.reserve(alive);
                for (int32_t w = 0; w < m; ++w) re.push_back({pr[w], w});
                pq = std::priority_queue<P, vector<P>, std::greater<P>>(std::greater<P>(), std::move(re));
                next_compact = alive / 2;
                debug_state("compact");
            }
        }
    }
    void run(int32_t rounds_until) {
        uoff.assign(1, 0);
        level.assign(m, 0); pr.assign(m, 0); dirty.assign(m, 0); dead.assign(m, 0);
        plevel.assign(m, 0); sdeg.assign(m, 0);
        alive = m;
        for (int32_t v = 0; v < m; ++v) pr[v] = prio(v);
        int32_t ru = std::max(rounds_until, STOPC);
        contract_rounds(ru);
        if (alive > STOPC) contract_greedy(STOPC);
        if (m != alive) compact();
        debug_state("core");
    }
};

// ------------------------------------------------------------ query structure
template <class QH>
struct Hier {
    int32_t n = 0, nc = 0, C = 0;          // vertices, contracted, core
    vector<int32_t> qid;                    // vertex id -> query index (core last)
    vector<uint32_t> H; vector<A16> Au;     // upward arcs by query index (contracted only)
    vector<uint32_t> cH; vector<A16> cA;    // core graph (core indices)
    vector<D> dF, dB; vector<int32_t> tF, tB; QH hF, hB;
    vector<int32_t> eF, eB;
    uint64_t qsettle = 0, qrel = 0, qes = 0, qcore = 0;
    // landmarks on the core
    int K = 0;
    static constexpr int64_t LINF = INT64_MAX / 4;
    vector<double> LD;                      // LD[c*K + l] = d(core c, landmark l) (exact: < 2^53)
    vector<int64_t> pot, potT, potS; vector<uint32_t> potst; uint32_t potc = 0;
    alignas(32) double Us[64], Ut[64];
    bool packedOK = false;
    D UBq = INF;

    template <class A>
    void build(Engine<A>& g, int32_t nv) {
        n = nv; nc = (int32_t)g.rorder.size();
        packedOK = std::is_same<A, A8>::value;
        C = g.m;                            // after the final compaction only core slots remain
        vector<int32_t> ordv(g.rorder);
        std::sort(ordv.begin(), ordv.end(), [&](int32_t a, int32_t b) { return g.vlevel[a] != g.vlevel[b] ? g.vlevel[a] < g.vlevel[b] : a < b; });
        qid.assign(n, -1);
        for (int32_t i = 0; i < nc; ++i) qid[ordv[i]] = i;
        for (int32_t c = 0; c < C; ++c) qid[g.orig[c]] = nc + c;
        {
            vector<int32_t> rank(n, -1);
            for (int32_t r = 0; r < nc; ++r) rank[g.rorder[r]] = r;
            H.assign(n + 1, 0);
            for (int32_t i = 0; i < nc; ++i) { int32_t r = rank[ordv[i]]; H[i + 1] = H[i] + (g.uoff[r + 1] - g.uoff[r]); }
            for (int32_t i = nc; i < n; ++i) H[i + 1] = H[i];
            Au.resize(H[nc]);
            for (int32_t i = 0; i < nc; ++i) {
                int32_t r = rank[ordv[i]]; uint32_t k = H[i];
                for (uint32_t j = g.uoff[r]; j < g.uoff[r + 1]; ++j) Au[k++] = amk16(qid[g.uarc[j].t], g.uarc[j].w);
            }
        }
        vector<A16>().swap(g.uarc);
        cH.assign(C + 1, 0);
        for (int32_t i = 0; i < C; ++i) cH[i + 1] = cH[i] + g.vr[i].deg;
        cA.resize(cH[C]);
        for (int32_t i = 0; i < C; ++i) {
            const A* a = g.nb(i);
            for (uint32_t q = 0; q < g.vr[i].deg; ++q) cA[cH[i] + q] = amk16(ato(a[q]), aw(a[q]));
        }
        dF.assign(n, INF); dB.assign(n, INF);
        tlog("query graph");
        landmarks((int)envi("W_LMK", 64));
        if (g_debug) std::fprintf(stderr, "core C=%d arcs=%u upward arcs=%u landmarks=%d\n", C, cH[C], H[n], K);
        tlog("landmarks");
    }
    void landmarks(int k) {
        K = std::min(k, 64) & ~3; if (!K || !C) { K = 0; return; }
        LD.assign((size_t)C * K, 0.0);
        vector<int64_t> d(C), mind(C, LINF);
        Heap4<D> h; PH ph;
        // packed core arcs for the landmark Dijkstras (when weights allow)
        bool packed = packedOK && C < (1 << 22);
        vector<uint64_t> pa;
        if (packed) { pa.resize(cA.size()); for (size_t q = 0; q < cA.size(); ++q) pa[q] = (cA[q].w << 22) | (uint32_t)cA[q].t; }
        int32_t src = 0;
        for (int l = 0; l < K; ++l) {
            std::fill(d.begin(), d.end(), LINF);
            d[src] = 0;
            if (packed) {
                ph.clear(); ph.push(0, src);
                while (!ph.empty()) {
                    auto it = ph.pop(); if ((int64_t)it.d > d[it.v]) continue;
                    for (uint32_t q = cH[it.v]; q < cH[it.v + 1]; ++q) {
                        int32_t y = (int32_t)(pa[q] & 0x3FFFFFu);
                        int64_t nd = (int64_t)it.d + (int64_t)(pa[q] >> 22);
                        if (nd < d[y]) { d[y] = nd; ph.push((D)nd, y); }
                    }
                }
            } else {
                h.clear(); h.push(0, src);
                while (!h.empty()) {
                    auto it = h.pop(); if ((int64_t)it.d > d[it.v]) continue;
                    for (uint32_t q = cH[it.v]; q < cH[it.v + 1]; ++q) {
                        int64_t nd = (int64_t)it.d + (int64_t)cA[q].w;
                        if (nd < d[cA[q].t]) { d[cA[q].t] = nd; h.push((D)nd, cA[q].t); }
                    }
                }
            }
            int32_t far = src; int64_t fd = -1;
            for (int32_t v = 0; v < C; ++v) {
                if (d[v] == LINF || d[v] >= (1ll << 52)) { K = 0; LD.clear(); return; }   // disconnected core: no ALT
                LD[(size_t)v * K + l] = (double)d[v];
                if (d[v] < mind[v]) mind[v] = d[v];
                if (mind[v] > fd) { fd = mind[v]; far = v; }
            }
            src = far;
        }
        pot.assign(C, 0); potT.assign(C, 0); potS.assign(C, 0); potst.assign(C, 0); potc = 0;
    }
    // x2-scaled potential of core vertex c: pi_t(c) - pi_s(c), with
    // pi_x(c) = max_l |d(c,l) - d(x,l)| (exact in doubles)
    inline int64_t potential(int32_t c) {
        if (potst[c] == potc) return pot[c];
        const double* ld = &LD[(size_t)c * K];
        v4d mt = {0, 0, 0, 0}, ms = {0, 0, 0, 0};
        for (int l = 0; l < K; l += 4) {
            v4d a; std::memcpy(&a, ld + l, sizeof(a));
            v4d x = a - *(const v4d*)(Ut + l); x = x < 0 ? -x : x;
            v4d y = a - *(const v4d*)(Us + l); y = y < 0 ? -y : y;
            mt = x > mt ? x : mt; ms = y > ms ? y : ms;
        }
        double pt = std::max(std::max(mt[0], mt[1]), std::max(mt[2], mt[3]));
        double ps = std::max(std::max(ms[0], ms[1]), std::max(ms[2], ms[3]));
        potst[c] = potc; potT[c] = (int64_t)pt; potS[c] = (int64_t)ps; pot[c] = potT[c] - potS[c];
        return pot[c];
    }
    inline void step(vector<D>& d1, vector<int32_t>& tt, QH& h, const vector<D>& d2, D& mu, vector<int32_t>& ent) {
        auto it = h.pop();
        int32_t u = it.v; D d = it.d;
        if (d > d1[u]) return;
        ++qsettle;
        if (d2[u] < INF && d + d2[u] < mu) mu = d + d2[u];
        if (u >= nc) { ent.push_back(u); return; }
        const uint32_t e0 = H[u], e1 = H[u + 1];
        for (uint32_t k = e0; k < e1; ++k) { const A16& a = Au[k]; if (d1[a.t] + a.w < d) return; }   // stall
        qrel += e1 - e0;
        for (uint32_t k = e0; k < e1; ++k) {
            const A16& a = Au[k]; D nd = d + a.w;
            if (nd > UBq) continue;            // longer than any shortest path
            if (nd < d1[a.t]) { if (d1[a.t] == INF) tt.push_back(a.t); d1[a.t] = nd; h.push(nd, a.t); }
        }
    }
    int64_t query(int32_t s, int32_t t) {
        if (s == t) return 0;
        s = qid[s]; t = qid[t];
        D mu = INF;
        eF.clear(); eB.clear();
        hF.clear(); hB.clear();
        dF[s] = 0; tF.push_back(s); hF.push(0, s);
        dB[t] = 0; tB.push_back(t); hB.push(0, t);
        // complete upward searches (exact entry distances needed by the landmarks)
        while (!hF.empty() || !hB.empty()) {
            if (!hF.empty() && (hB.empty() || hF.top_key() <= hB.top_key())) step(dF, tF, hF, dB, mu, eF);
            else step(dB, tB, hB, dF, mu, eB);
        }
        qes += eF.size() + eB.size();
        if (!eF.empty() && !eB.empty()) {
            if (K > 0) {
                for (int l = 0; l < K; ++l) { Us[l] = 1e300; Ut[l] = 1e300; }
                for (int32_t a : eF) { const double* ld = &LD[(size_t)(a - nc) * K]; double da = (double)dF[a];
                    for (int l = 0; l < K; ++l) Us[l] = std::min(Us[l], da + ld[l]); }
                for (int32_t b : eB) { const double* ld = &LD[(size_t)(b - nc) * K]; double db = (double)dB[b];
                    for (int l = 0; l < K; ++l) Ut[l] = std::min(Ut[l], db + ld[l]); }
                if (++potc == 0) { std::fill(potst.begin(), potst.end(), 0); potc = 1; }
            }
            auto P = [&](int32_t c) -> int64_t { return K > 0 ? potential(c) : 0; };
            hF.clear(); hB.clear();
            for (int32_t a : eF) hF.push((D)(2 * (int64_t)dF[a] + P(a - nc)), a);
            for (int32_t b : eB) hB.push((D)(2 * (int64_t)dB[b] - P(b - nc)), b);
            for (;;) {
                if (hF.empty() || hB.empty()) break;
                D kf = hF.top_key(), kb = hB.top_key();
                if (mu < INF && kf + kb >= 2 * mu) break;
                bool fw = kf <= kb;
                QH& h = fw ? hF : hB;
                vector<D>& d1 = fw ? dF : dB; const vector<D>& d2 = fw ? dB : dF;
                vector<int32_t>& tt = fw ? tF : tB;
                auto it = h.pop();
                int32_t cu = it.v - nc;
                int64_t pu = P(cu);
                D du = d1[it.v];
                if ((int64_t)it.d > 2 * (int64_t)du + (fw ? pu : -pu)) continue;   // stale
                ++qcore;
                for (uint32_t k = cH[cu]; k < cH[cu + 1]; ++k) {
                    int32_t x = cA[k].t + nc; D nd = du + cA[k].w;
                    if (nd > UBq) continue;
                    if (nd < d1[x]) {
                        if (d1[x] == INF) tt.push_back(x);
                        d1[x] = nd;
                        if (d2[x] < INF && nd + d2[x] < mu) mu = nd + d2[x];
                        int64_t px = P(cA[k].t);
                        // A* pruning: nd + (lower bound towards the other side) >= mu
                        if (K > 0 && (int64_t)nd + (fw ? potT[cA[k].t] : potS[cA[k].t]) >= (int64_t)mu) continue;
                        h.push((D)(2 * (int64_t)nd + (fw ? px : -px)), x);
                    }
                }
            }
        }
        for (int32_t v : tF) dF[v] = INF;
        for (int32_t v : tB) dB[v] = INF;
        tF.clear(); tB.clear();
        return mu >= INF ? -1 : (int64_t)mu;
    }
};

// ------------------------------------------------------------ driver
template <class HQ, class A>
static void answer(Engine<A>& g, const vector<int32_t>& perm) {
    HQ hq;
    hq.UBq = g.UB;
    hq.build(g, V);
    {
        vector<std::pair<int64_t, int32_t>> qo(Q);
        for (int32_t i = 0; i < Q; ++i) qo[i] = {((int64_t)perm[qs[i]] << 32) | (uint32_t)perm[qt[i]], i};
        std::sort(qo.begin(), qo.end());
        for (int32_t j = 0; j < Q; ++j) { int32_t i = qo[j].second; qans[i] = hq.query(perm[qs[i]], perm[qt[i]]); }
    }
    if (g_debug) std::fprintf(stderr, "per query: settles %.1f relaxed %.1f entries %.1f core %.1f\n", (double)hq.qsettle / Q,
                              (double)hq.qrel / Q, (double)hq.qes / Q, (double)hq.qcore / Q);
}

template <class A>
static void run(const vector<int32_t>& perm, const vector<uint8_t>& deadE, D ub) {
    Engine<A> g;
    g.UB = ub;
    g.CON_SET = (int)envi("W_CONSET", 200); g.LC = envi("W_LC", 250); g.STOPC = (int32_t)envi("W_STOPC", 120000); g.LAZY = (int)envi("W_LAZY", 1); g.LMK_N = (int)envi("W_ELMK", 0); g.LMK_AT = (int)envi("W_ELMK_AT", 1100000); g.LAZYT = (int)envi("W_LAZYT", 1000000);
    {
        vector<uint32_t> deg0(V, 0);
        for (int32_t e = 0; e < E; ++e) if (!deadE[e]) { ++deg0[perm[eu[e]]]; ++deg0[perm[ev[e]]]; }
        g.alloc_slots(V);
        size_t tot = 0;
        for (int32_t v = 0; v < V; ++v) { g.vr[v].off = (uint32_t)tot; g.vr[v].cap = deg0[v] + 4; tot += g.vr[v].cap; }
        g.pool.reserve(tot * 2); g.pool.resize(tot);
        g.orig.resize(V); for (int32_t v = 0; v < V; ++v) g.orig[v] = v;
        for (int32_t e = 0; e < E; ++e) if (!deadE[e]) g.add_edge(perm[eu[e]], perm[ev[e]], (D)ew[e]);
    }
    g.vlevel.assign(V, 0);
    tlog("build");
    g.run((int32_t)envi("W_ROUNDS_UNTIL", 0));
    tlog("contract");
    if (envi("W_EXIT", 0)) std::_Exit(0);
    // packed query heaps when every key (< 3*UB for the ALT keys) fits 42 bits
    if (ub < (1ull << 42) / 3 && V < (1 << 22)) answer<Hier<PH>>(g, perm);
    else answer<Hier<Heap4<D>>>(g, perm);
    tlog("queries");
}

void solve() {
    detect();
    // vertex numbering: Morton order of lattice coordinates, else identity
    vector<int32_t> perm(V);
    for (int32_t v = 0; v < V; ++v) perm[v] = v;
    if (dim) {
        // rank of each vertex's Morton code: mark present codes, then scan
        int bits = 0; while ((1ll << bits) < side) ++bits;
        const uint64_t ncode = 1ull << (bits * dim);
        vector<int32_t> rk(ncode, -1);
        vector<uint64_t> code(V);
        for (int64_t v = 0; v < V; ++v) {
            uint64_t mo = 0;
            for (int k = 0; k < dim; ++k) {
                uint64_t c = (uint64_t)((v / strd[k]) % side);
                for (int b = 0; b < bits; ++b) mo |= ((c >> b) & 1ull) << (b * dim + k);
            }
            code[v] = mo; rk[mo] = 0;
        }
        int32_t r = 0;
        for (uint64_t c = 0; c < ncode; ++c) if (rk[c] == 0) rk[c] = r++;
        for (int64_t v = 0; v < V; ++v) perm[v] = rk[code[v]];
    }
    tlog("order");
    // edge pruning: an edge (u,v,w) is dropped if some path of up to three
    // edges is strictly shorter; parallel edges and self-loops go too
    vector<uint8_t> deadE(E, 0);
    D maxw = 0;
    for (int32_t e = 0; e < E; ++e) if (ew[e] > maxw) maxw = ew[e];
    if (dim == 2) {
        // 2-D torus: an edge is dominated by one of its two square detours
        const int64_t n = side;
        int64_t cut = 0;
        auto R = [&](int64_t i) { int64_t x = i % n; return x + 1 < n ? i + 1 : i - x; };
        auto L = [&](int64_t i) { int64_t x = i % n; return x > 0 ? i - 1 : i + n - 1; };
        auto U = [&](int64_t i) { int64_t y = i / n; return y + 1 < n ? i + n : i - y * n; };
        auto Dn = [&](int64_t i) { int64_t y = i / n; return y > 0 ? i - n : i + (n - 1) * n; };
        for (int64_t i = 0; i < V; ++i) {
            {   // horizontal edge i -> R(i)
                D w = ew[2 * i]; int64_t r = R(i), u = U(i), d = Dn(i), dr = R(d);
                D up = (D)ew[2 * i + 1] + ew[2 * u] + ew[2 * r + 1];
                D dn = (D)ew[2 * d + 1] + ew[2 * d] + ew[2 * dr + 1];
                if (up < w || dn < w) { deadE[2 * i] = 1; ++cut; }
            }
            {   // vertical edge i -> U(i)
                D w = ew[2 * i + 1]; int64_t u = U(i), r = R(i), l = L(i), lu = U(l);
                D rt = (D)ew[2 * i] + ew[2 * r + 1] + ew[2 * u];
                D lt = (D)ew[2 * l] + ew[2 * l + 1] + ew[2 * lu];
                if (rt < w || lt < w) { deadE[2 * i + 1] = 1; ++cut; }
            }
        }
        if (g_debug) std::fprintf(stderr, "prune (lattice): %lld of %d edges\n", (long long)cut, E);
    } else {

        vector<uint32_t> h(V + 1, 0);
        for (int32_t e = 0; e < E; ++e) if (eu[e] != ev[e]) { ++h[eu[e] + 1]; ++h[ev[e] + 1]; }
        for (int32_t v = 0; v < V; ++v) h[v + 1] += h[v];
        vector<std::pair<int32_t, uint32_t>> adj(h[V]);
        { vector<uint32_t> c(h.begin(), h.end() - 1);
          for (int32_t e = 0; e < E; ++e) if (eu[e] != ev[e]) { adj[c[eu[e]]++] = {ev[e], ew[e]}; adj[c[ev[e]]++] = {eu[e], ew[e]}; } }
        int64_t cut = 0;
        for (int32_t e = 0; e < E; ++e) {
            int32_t u = eu[e], v = ev[e]; D w = ew[e];
            if (u == v) { deadE[e] = 1; continue; }
            bool dom = false;
            for (uint32_t i = h[u]; i < h[u + 1] && !dom; ++i) {
                int32_t a = adj[i].first; D wa = adj[i].second;
                if (a == v) { if (wa < w) dom = true; continue; }
                if (wa >= w) continue;
                for (uint32_t j = h[a]; j < h[a + 1] && !dom; ++j) {
                    int32_t b = adj[j].first; D wb = wa + adj[j].second;
                    if (b == u || wb >= w) continue;
                    if (b == v) { dom = true; break; }
                    for (uint32_t q = h[b]; q < h[b + 1]; ++q) if (adj[q].first == v && wb + adj[q].second < w) { dom = true; break; }
                }
            }
            if (dom) { deadE[e] = 1; ++cut; }
        }
        if (g_debug) std::fprintf(stderr, "prune: %lld of %d edges\n", (long long)cut, E);
    }
    tlog("prune");
    // Upper bound on any shortest-path length (lattice: max weight times the
    // torus diameter in hops); arcs longer than it are dropped.  Packed 8-byte
    // arcs need it below 2^42 and V below 2^22.
    D ub = INF;
    if (dim) ub = maxw * (D)(dim * (side / 2) + 1);
    if (ub < (1ull << 42) && V < (1 << 22) && !envi("W_A16", 0)) run<A8>(perm, deadE, ub);
    else run<A16>(perm, deadE, ub);
}

}  // namespace w64
#pragma GCC pop_options



namespace lat {

// Monotone radix heap (keys never decrease below the last popped key).
template <class K>
struct RadixHeap {
    static const int NB = (int)sizeof(K) * 8 + 1;
    struct Item { K d; int32_t v; };
    vector<Item> b[NB]; K last = 0; size_t sz = 0; int hi = 0;
    static inline int bucket(K x, K l) {
        if (x == l) return 0;
        if (sizeof(K) == 8) return 64 - __builtin_clzll((unsigned long long)(x ^ l));
        return 32 - __builtin_clz((unsigned)(x ^ l));
    }
    inline void clear() { for (int i = 0; i <= hi; ++i) b[i].clear(); last = 0; sz = 0; hi = 0; }
    inline bool empty() const { return sz == 0; }
    inline void push(K d, int32_t v) { int i = bucket(d, last); b[i].push_back({d, v}); if (i > hi) hi = i; ++sz; }
    inline void refill() {
        if (!b[0].empty()) return;
        int i = 1; while (b[i].empty()) ++i;
        K m = b[i][0].d; for (const Item& it : b[i]) if (it.d < m) m = it.d;
        last = m;
        for (const Item& it : b[i]) b[bucket(it.d, last)].push_back(it);
        b[i].clear();
    }
    inline K top_key() { refill(); return last; }
    inline Item pop() { refill(); Item it = b[0].back(); b[0].pop_back(); --sz; return it; }
};

static int64_t envi(const char* k, int64_t def) { const char* e = std::getenv(k); return e ? std::atoll(e) : def; }
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// ------------------------------------------------------------ lattice detect
static int dim = 0; static int64_t side = 0; static int64_t strd[3];
static bool detect() {
    if (directed) return false;
    for (int dd = 2; dd <= 3; ++dd) {
        int64_t s = llround(std::pow((double)V, 1.0 / dd));
        for (int64_t c = std::max<int64_t>(3, s - 1); c <= s + 1; ++c) {
            int64_t p = 1; for (int k = 0; k < dd; ++k) p *= c;
            if (p == V && (int64_t)E == dd * (int64_t)V) { dim = dd; side = c; }
        }
    }
    if (!dim) return false;
    strd[0] = 1; strd[1] = side; strd[2] = side * side;
    for (int64_t i = 0; i < V; ++i) for (int k = 0; k < dim; ++k) {
        int64_t e = i * dim + k; int64_t c = (i / strd[k]) % side;
        int64_t j = c + 1 < side ? i + strd[k] : i - c * strd[k];
        if (eu[e] != i || ev[e] != j) { dim = 0; return false; }
    }
    return true;
}

// ------------------------------------------------- nested dissection order
// Geometric ND on the torus: every periodic axis is first opened by removing
// one hyperplane, then boxes are bisected along their longest axis.  A
// separator is itself ordered by ND (it is a lower-dimensional box).  Output:
// vertices in contraction order (least important first).
struct Box { int64_t lo[3], len[3]; bool per[3]; };
static vector<int32_t> nd_out;
static int64_t vid(const int64_t* c) { int64_t id = 0; for (int k = 0; k < dim; ++k) id += (c[k] % side) * strd[k]; return id; }
static void nd_rec(Box b) {
    int64_t cnt = 1; for (int k = 0; k < dim; ++k) cnt *= b.len[k];
    if (cnt <= 0) return;
    if (cnt == 1) { int64_t c[3] = {b.lo[0], b.lo[1], b.lo[2]}; nd_out.push_back((int32_t)vid(c)); return; }
    // choose axis: a periodic axis first, else the longest
    int ax = -1;
    for (int k = 0; k < dim; ++k) if (b.per[k] && b.len[k] > 1) { ax = k; break; }
    if (ax < 0) { ax = 0; for (int k = 1; k < dim; ++k) if (b.len[k] > b.len[ax]) ax = k; }
    if (b.len[ax] == 1) { // should not happen unless cnt==1
        return;
    }
    Box sep = b, l = b, r = b;
    if (b.per[ax]) {
        // remove hyperplane at lo; remaining is a non-periodic run of len-1
        sep.len[ax] = 1; sep.per[ax] = false;
        l.lo[ax] = b.lo[ax] + 1; l.len[ax] = b.len[ax] - 1; l.per[ax] = false;
        nd_rec(l);
        nd_rec(sep);
        return;
    }
    int64_t mid = b.len[ax] / 2;
    l.len[ax] = mid;
    sep.lo[ax] = b.lo[ax] + mid; sep.len[ax] = 1;
    r.lo[ax] = b.lo[ax] + mid + 1; r.len[ax] = b.len[ax] - mid - 1;
    nd_rec(l); nd_rec(r); nd_rec(sep);
}
static void nd_order() {
    nd_out.clear(); nd_out.reserve(V);
    Box b; for (int k = 0; k < 3; ++k) { b.lo[k] = 0; b.len[k] = k < dim ? side : 1; b.per[k] = k < dim; }
    nd_rec(b);
}

// ------------------------------------------------------------------- CH
template <class D>
struct CH {
    struct Arc { int32_t to; D w; };
    static constexpr D INF = (D)(~(D)0) >> 2;
    int32_t n;
    // flat adjacency arena
    vector<Arc> pool; vector<uint32_t> aoff, adeg, acap;
    struct St { D d; uint32_t stamp; uint32_t tm; uint32_t ms; };   // search dist/stamp, target mark, sim stamp
    vector<St> st; vector<D> md; uint32_t wc = 0, tc = 0, mc = 0;
    vector<D> tcost, tlim;
    vector<int32_t> rank_of;
    // upward arcs in contraction order
    vector<Arc> uarc; vector<uint32_t> uoff;   // uoff indexed by rank
    Heap4<D> heap;
    uint64_t settles = 0, searches = 0, capped = 0, capset = 0; uint64_t caphist[64] = {0}; uint64_t allhist[64] = {0};
    int CON_SET = 200;
    int64_t LC = 1000; int32_t STOPC = 0;
    double tprio = 0, tcon = 0;
    D UB = INF;
    int K = 0; vector<D> lmd; uint64_t lm_skip = 0, lm_need = 0;
    // K landmarks by farthest-point selection; Dijkstra over the current graph.
    void landmarks(int k) {
        K = k; if (!K) return;
        lmd.assign((size_t)n * K, INF);
        vector<D> mind(n, INF), d(n);
        RadixHeap<D> hp;
        int32_t src = 0;
        for (int q = 0; q < K; ++q) {
            std::fill(d.begin(), d.end(), INF);
            d[src] = 0; hp.clear(); hp.push(0, src);
            while (!hp.empty()) {
                auto it = hp.pop(); if (it.d > d[it.v]) continue;
                const Arc* a = nb(it.v); uint32_t kk = adeg[it.v];
                for (uint32_t i = 0; i < kk; ++i) { D nd = it.d + a[i].w; if (nd < d[a[i].to]) { d[a[i].to] = nd; hp.push(nd, a[i].to); } }
            }
            int32_t far = 0; D fd = 0;
            for (int32_t v = 0; v < n; ++v) { lmd[(size_t)v * K + q] = d[v]; if (d[v] < mind[v]) mind[v] = d[v]; if (mind[v] != INF && mind[v] > fd) { fd = mind[v]; far = v; } }
            src = far;
        }
    }

    inline Arc* nb(int32_t v) { return &pool[aoff[v]]; }
    void init(int32_t n_, const vector<int32_t>& deg0) {
        n = n_;
        aoff.resize(n); adeg.assign(n, 0); acap.resize(n);
        size_t tot = 0;
        for (int32_t v = 0; v < n; ++v) { aoff[v] = (uint32_t)tot; acap[v] = (uint32_t)deg0[v] + 2; tot += acap[v]; }
        pool.resize(tot + tot / 2 + 16); pool.resize(tot);
        st.assign(n, {INF, 0, 0, 0}); md.assign(n, INF); tcost.assign(n, 0); tlim.assign(n, 0);
        rank_of.assign(n, -1);
    }
    inline void push_arc(int32_t u, int32_t v, D w) {
        if (adeg[u] == acap[u]) {
            uint32_t nc = acap[u] * 2;
            uint32_t no = (uint32_t)pool.size();
            pool.resize(pool.size() + nc);
            std::memcpy(&pool[no], &pool[aoff[u]], sizeof(Arc) * adeg[u]);
            aoff[u] = no; acap[u] = nc;
        }
        pool[aoff[u] + adeg[u]++] = {v, w};
    }
    void add_edge(int32_t u, int32_t v, D w) {
        Arc* a = nb(u); uint32_t k = adeg[u];
        for (uint32_t i = 0; i < k; ++i) if (a[i].to == v) {
            if (w < a[i].w) { a[i].w = w; Arc* b = nb(v); for (uint32_t j = 0; j < adeg[v]; ++j) if (b[j].to == u) { b[j].w = w; break; } }
            return;
        }
        push_arc(u, v, w); push_arc(v, u, w);
    }
    inline void remove_arc(int32_t u, int32_t v) {
        Arc* a = nb(u); uint32_t k = adeg[u];
        for (uint32_t i = 0; i < k; ++i) if (a[i].to == v) { a[i] = a[k - 1]; --adeg[u]; return; }
    }

    // Witness search from src avoiding `avoid` for targets tl (tm == tc).
    vector<int32_t> tl;
    // Goal-directed variant: key = d + pi(y), pi(y) = min over pending targets of
    // the landmark lower bound LB(y, x).  Stops when key > max pending tcost.
    int ASTAR = 0;
    inline D lbk(int32_t a, int32_t b) const {
        const D* x = &lmd[(size_t)a * K]; const D* y = &lmd[(size_t)b * K];
        D m = 0; for (int q = 0; q < K; ++q) { D df = x[q] > y[q] ? x[q] - y[q] : y[q] - x[q]; if (df > m) m = df; }
        return m;
    }
    vector<int32_t> pend;
    inline D pot(int32_t y) {
        D m = INF; for (int32_t x : pend) if (st[x].tm == tc) { D l = lbk(y, x); if (l < m) m = l; }
        return m;
    }
    vector<D> potc;   // potential cache, valid while st[v].stamp == wc
    void witness_astar(int32_t src, int32_t avoid, int budget, int pending) {
        ++searches;
        if (++wc == 0) { for (St& s : st) s.stamp = 0; wc = 1; }
        if (potc.size() != (size_t)n) potc.assign(n, 0);
        pend.clear(); D limit = 0;
        for (int32_t x : tl) if (st[x].tm == tc) { pend.push_back(x); limit = std::max(limit, tcost[x]); }
        heap.clear();
        st[src].d = 0; st[src].stamp = wc; potc[src] = pot(src); heap.push(potc[src], src);
        int settled = 0;
        while (!heap.empty()) {
            auto it = heap.pop();
            if (it.d > limit) break;
            D d = st[it.v].d;
            if (it.d != d + potc[it.v]) continue;               // stale entry
            if (++settled > budget) { ++capped; capset += settled; break; }
            ++settles;
            bool relim = false;
            const Arc* a = nb(it.v); uint32_t k = adeg[it.v];
            for (uint32_t i = 0; i < k; ++i) {
                int32_t y = a[i].to;
                if (y == avoid) continue;
                D nd = d + a[i].w;
                St& sy = st[y];
                if (sy.tm == tc && nd <= tcost[y]) {
                    sy.tm = 0;
                    if (--pending <= 0) return;
                    relim = true;
                }
                if (nd > limit) continue;
                if (sy.stamp != wc) { sy.stamp = wc; sy.d = nd; potc[y] = pot(y); if (nd + potc[y] <= limit) heap.push(nd + potc[y], y); }
                else if (nd < sy.d) { sy.d = nd; if (nd + potc[y] <= limit) heap.push(nd + potc[y], y); }
            }
            if (relim) { limit = 0; for (int32_t x : pend) if (st[x].tm == tc) limit = std::max(limit, tcost[x]); }
        }
    }
    void witness(int32_t src, int32_t avoid, int budget, int pending) {
        ++searches;
        if (++wc == 0) { for (St& s : st) s.stamp = 0; wc = 1; }
        D limit = 0;
        for (int32_t x : tl) if (st[x].tm == tc) limit = std::max(limit, tlim[x]);
        heap.clear();
        st[src].d = 0; st[src].stamp = wc; heap.push(0, src);
        int settled = 0;
        while (!heap.empty()) {
            auto it = heap.pop();
            if (it.d > st[it.v].d) continue;
            if (it.d > limit) break;
            if (++settled > budget) { ++capped; capset += settled; break; }
            ++settles;
            bool relim = false;
            const Arc* a = nb(it.v); uint32_t k = adeg[it.v];
            for (uint32_t i = 0; i < k; ++i) {
                int32_t y = a[i].to;
                if (y == avoid) continue;
                D nd = it.d + a[i].w;
                St& sy = st[y];
                if (sy.tm == tc && nd <= tcost[y]) {
                    sy.tm = 0;
                    if (--pending <= 0) return;
                    if (tlim[y] >= limit) relim = true;
                }
                if (nd > limit) continue;
                if (sy.stamp != wc || nd < sy.d) { sy.stamp = wc; sy.d = nd; heap.push(nd, y); }
            }
            if (relim) { limit = 0; for (int32_t x : tl) if (st[x].tm == tc) limit = std::max(limit, tlim[x]); }
        }
    }
    struct SC { int32_t a, b; D w; };
    vector<D> mwv; vector<SC> scs;
    // Needed shortcuts of v into scs.
    void shortcuts(int32_t v, int budget) {
        scs.clear();
        const Arc* nv = nb(v); uint32_t k = adeg[v];
        mwv.resize(k);
        for (uint32_t j = 0; j < k; ++j) {
            D m = INF; const Arc* a = nb(nv[j].to); uint32_t kk = adeg[nv[j].to];
            for (uint32_t q = 0; q < kk; ++q) if (a[q].to != v && a[q].w < m) m = a[q].w;
            mwv[j] = m;
        }
        for (uint32_t i = 0; i + 1 < k; ++i) {
            if (++tc == 0) { for (St& s : st) s.tm = 0; tc = 1; }
            tl.clear();
            int pending = 0;
            for (uint32_t j = i + 1; j < k; ++j) {
                int32_t x = nv[j].to; D c = nv[i].w + nv[j].w;
                st[x].tm = tc;
                if (c > UB) { st[x].tm = 0; continue; }                 // never a shortest path
                if (K) {
                    const D* lu = &lmd[(size_t)nv[i].to * K]; const D* lx = &lmd[(size_t)x * K];
                    D ub = INF, lb = 0;
                    for (int q = 0; q < K; ++q) { D a = lu[q], b = lx[q]; ub = std::min<D>(ub, a + b); D df = a > b ? a - b : b - a; lb = std::max(lb, df); }
                    if (ub < c) { st[x].tm = 0; ++lm_skip; continue; }     // a shorter u-x path exists
                    if (lb >= c) { tcost[x] = 0; ++lm_need; continue; }     // u-v-x is a shortest path
                }
                if (mwv[j] > c || mwv[i] > c) { tcost[x] = 0; continue; }  // no witness possible
                tcost[x] = c; tlim[x] = c - mwv[j]; tl.push_back(x); ++pending;
            }
            if (pending) { if (ASTAR && K) witness_astar(nv[i].to, v, budget, pending); else witness(nv[i].to, v, budget, pending); }
            for (uint32_t j = i + 1; j < k; ++j)
                if (st[nv[j].to].tm == tc) scs.push_back({nv[i].to, nv[j].to, nv[i].w + nv[j].w});
        }
    }
    // 2-hop simulation for the priority
    int simulate(int32_t v) {
        const Arc* nv = nb(v); uint32_t k = adeg[v];
        int added = 0;
        for (uint32_t i = 0; i + 1 < k; ++i) {
            int32_t u = nv[i].to;
            if (++mc == 0) { for (St& s : st) s.ms = 0; mc = 1; }
            const Arc* au = nb(u); uint32_t ku = adeg[u];
            for (uint32_t q = 0; q < ku; ++q) if (au[q].to != v) { st[au[q].to].ms = mc; md[au[q].to] = au[q].w; }
            for (uint32_t j = i + 1; j < k; ++j) {
                int32_t x = nv[j].to; D c = nv[i].w + nv[j].w;
                bool found = st[x].ms == mc && md[x] <= c;
                if (!found && k <= SIM2) { const Arc* ax = nb(x); uint32_t kx = adeg[x];
                    for (uint32_t q = 0; q < kx; ++q) { int32_t b = ax[q].to; if (b != v && st[b].ms == mc && md[b] + ax[q].w <= c) { found = true; break; } } }
                if (!found) ++added;
            }
        }
        return added;
    }
    vector<int32_t> level;
    int PMODE = 0, QORD = 1, STALL = 1, TIMEC = 0; double tcomb = 0; int64_t PED = 0; int ASORT = 1; uint32_t SIM2 = 1000000;
    int64_t prio(int32_t v) {
        int deg = (int)adeg[v];
        if (PMODE == 1) return (int64_t)level[v] * LC + 250 * (int64_t)deg;
        int add = simulate(v);
        return (int64_t)level[v] * LC + (1000 * (int64_t)add) / std::max(1, deg) + PED * ((int64_t)add - deg);
    }
    int32_t ncontracted = 0;
    void contract_one(int32_t v) {
        { double t0 = now(); shortcuts(v, CON_SET); tcon += now() - t0; }
        int32_t r = ncontracted++;
        rank_of[v] = r;
        const Arc* nv = nb(v); uint32_t k = adeg[v];
        for (uint32_t i = 0; i < k; ++i) { uarc.push_back(nv[i]); remove_arc(nv[i].to, v); }
        uoff.push_back((uint32_t)uarc.size());
        adeg[v] = 0;
        for (const SC& s : scs) add_edge(s.a, s.b, s.w);
    }
    void progress() {
        int32_t r = ncontracted;
        if (g_debug && (r % (n / 10) == 0 || n - r == 3000 || n - r == 1000)) {
            size_t arcs = 0, mx = 0; for (int32_t x = 0; x < n; ++x) if (rank_of[x] < 0) { arcs += adeg[x]; mx = std::max<size_t>(mx, adeg[x]); }
            std::fprintf(stderr, "  %d/%d searches=%llu settles=%llu capped=%llu (%llu) avgdeg=%.2f maxdeg=%zu tprio=%.2f tcon=%.2f pool=%zu\n", r, n,
                (unsigned long long)searches, (unsigned long long)settles, (unsigned long long)capped, (unsigned long long)capset,
                (double)arcs / std::max(1, n - r), mx, tprio, tcon, pool.size());
            tlog("progress");
        }
    }
    void contract_greedy() {
        level.assign(n, 0);
        uoff.assign(1, 0);
        vector<int64_t> pr(n);
        vector<uint8_t> dirty(n, 0);
        typedef std::pair<int64_t, int32_t> P;
        vector<P> init(n);
        double t0 = now();
        for (int32_t v = 0; v < n; ++v) { pr[v] = prio(v); init[v] = {pr[v], v}; }
        tprio += now() - t0;
        std::priority_queue<P, vector<P>, std::greater<P>> pq(std::greater<P>(), std::move(init));
        tlog("init prio");
        while (!pq.empty()) {
            P t = pq.top(); pq.pop();
            int32_t v = t.second;
            if (rank_of[v] >= 0 || t.first != pr[v]) continue;
            if (n - ncontracted <= STOPC) break;
            if (dirty[v]) {
                dirty[v] = 0;
                double t1 = now();
                int64_t np = prio(v);
                tprio += now() - t1;
                pr[v] = np;
                if (!pq.empty() && np > pq.top().first) { pq.push({np, v}); continue; }
            }
            const Arc* nv = nb(v); uint32_t k = adeg[v];
            for (uint32_t i = 0; i < k; ++i) { level[nv[i].to] = std::max(level[nv[i].to], level[v] + 1); dirty[nv[i].to] = 1; }
            contract_one(v);
            progress();
        }
    }
    // Rebuilds the arena so that adjacency lists are stored in vertex order.
    void compact_pool() {
        vector<Arc> np; size_t tot = 0;
        for (int32_t v = 0; v < n; ++v) if (rank_of[v] < 0) tot += adeg[v] + adeg[v] / 2 + 2;
        np.reserve(tot + tot / 2);
        for (int32_t v = 0; v < n; ++v) {
            if (rank_of[v] >= 0) { aoff[v] = 0; acap[v] = 0; adeg[v] = 0; continue; }
            uint32_t o = (uint32_t)np.size(); uint32_t c = adeg[v] + adeg[v] / 2 + 2;
            np.insert(np.end(), pool.begin() + aoff[v], pool.begin() + aoff[v] + adeg[v]);
            np.resize(o + c);
            aoff[v] = o; acap[v] = c;
        }
        pool.swap(np);
    }
    // Rounds of independent sets: every vertex whose priority is a strict
    // local minimum (ties broken by id) is contracted, in id order.
    void contract_rounds() {
        level.assign(n, 0);
        uoff.assign(1, 0);
        vector<int64_t> pr(n);
        vector<uint8_t> dirty(n, 0);
        double t0 = now();
        for (int32_t v = 0; v < n; ++v) pr[v] = prio(v);
        tprio += now() - t0;
        vector<int32_t> alive(n); for (int32_t v = 0; v < n; ++v) alive[v] = v;
        vector<int32_t> cand, touched;
        int rounds = 0;
        while ((int32_t)alive.size() > STOPC) {
            ++rounds;
            cand.clear();
            for (int32_t v : alive) {
                const Arc* a = nb(v); uint32_t k = adeg[v]; bool mn = true;
                int64_t p = pr[v];
                for (uint32_t i = 0; i < k; ++i) { int32_t u = a[i].to; if (pr[u] < p || (pr[u] == p && u < v)) { mn = false; break; } }
                if (mn) cand.push_back(v);
            }
            int64_t room = (int64_t)alive.size() - STOPC;
            if ((int64_t)cand.size() > room) cand.resize(room);
            touched.clear();
            for (int32_t v : cand) {
                const Arc* nv = nb(v); uint32_t k = adeg[v];
                for (uint32_t i = 0; i < k; ++i) { int32_t u = nv[i].to; level[u] = std::max(level[u], level[v] + 1); if (!dirty[u]) { dirty[u] = 1; touched.push_back(u); } }
                contract_one(v);
                progress();
            }
            double t1 = now();
            std::sort(touched.begin(), touched.end());
            for (int32_t u : touched) { dirty[u] = 0; if (rank_of[u] < 0) pr[u] = prio(u); }
            tprio += now() - t1;
            size_t w = 0; for (int32_t v : alive) if (rank_of[v] < 0) alive[w++] = v; alive.resize(w);
            if (rounds % 4 == 0) compact_pool();
        }
        if (g_debug) std::fprintf(stderr, "rounds=%d\n", rounds);
    }
    void contract_order(const vector<int32_t>& ord) {
        uoff.assign(1, 0);
        for (int32_t v : ord) {
            if (n - ncontracted <= STOPC) break;
            contract_one(v);
            progress();
        }
    }

    // ---------------------------------------------------------------- core
    int32_t C = 0;
    vector<int32_t> core_ids;
    vector<D> T;
    void build_core() {
        for (int32_t v = 0; v < n; ++v) if (rank_of[v] < 0) { rank_of[v] = ncontracted + (int32_t)core_ids.size(); core_ids.push_back(v); }
        C = (int32_t)core_ids.size();
        if (!C) return;
        vector<uint32_t> ch(C + 1, 0); vector<Arc> ca;
        for (int32_t i = 0; i < C; ++i) {
            const Arc* a = nb(core_ids[i]);
            for (uint32_t q = 0; q < adeg[core_ids[i]]; ++q) ca.push_back({rank_of[a[q].to] - ncontracted, a[q].w});
            ch[i + 1] = (uint32_t)ca.size();
        }
        T.assign((size_t)C * C, INF);
        vector<int32_t> cto(ca.size()); vector<D> cw(ca.size());
        for (size_t k = 0; k < ca.size(); ++k) { cto[k] = ca[k].to; cw[k] = ca[k].w; }
        RadixHeap<D> h;
        for (int32_t s = 0; s < C; ++s) {
            D* row = &T[(size_t)s * C];
            row[s] = 0; h.clear(); h.push(0, s);
            while (!h.empty()) {
                auto it = h.pop();
                if (it.d > row[it.v]) continue;
                const uint32_t e1 = ch[it.v + 1];
                for (uint32_t k = ch[it.v]; k < e1; ++k) {
                    D nd = it.d + cw[k]; int32_t y = cto[k];
                    if (nd < row[y]) { row[y] = nd; h.push(nd, y); }
                }
            }
        }
        if (g_debug) std::fprintf(stderr, "core C=%d arcs=%zu table=%.1f MB\n", C, ca.size(), (double)T.size() * sizeof(D) / 1e6);
    }
    // query structure: CSR by rank (contracted part only)
    vector<uint32_t> H; vector<Arc> A;
    struct DD { D d[2]; };            // forward / backward tentative distance
    vector<DD> dd; vector<int32_t> tF, tB; RadixHeap<D> hF, hB;
    // Query ids: contracted vertices ordered by (CH level, id) so that one
    // query touches nearby memory; core vertices last, in core-table order.
    vector<int32_t> qid;
    void build_query() {
        int32_t m = ncontracted;
        vector<int32_t> ord; ord.reserve(m);
        for (int32_t v = 0; v < n; ++v) if (rank_of[v] < m) ord.push_back(v);
        if (!level.empty() && QORD) std::sort(ord.begin(), ord.end(), [&](int32_t a, int32_t b) { return level[a] != level[b] ? level[a] < level[b] : a < b; });
        else std::sort(ord.begin(), ord.end(), [&](int32_t a, int32_t b) { return rank_of[a] < rank_of[b]; });
        qid.assign(n, -1);
        for (int32_t i = 0; i < m; ++i) qid[ord[i]] = i;
        for (int32_t c = 0; c < C; ++c) qid[core_ids[c]] = m + c;
        H.assign(n + 1, 0);
        for (int32_t i = 0; i < m; ++i) { int32_t r = rank_of[ord[i]]; H[i + 1] = H[i] + (uoff[r + 1] - uoff[r]); }
        for (int32_t i = m; i < n; ++i) H[i + 1] = H[i];
        A.resize(H[m]);
        for (int32_t i = 0; i < m; ++i) {
            int32_t r = rank_of[ord[i]]; uint32_t k = H[i];
            for (uint32_t j = uoff[r]; j < uoff[r + 1]; ++j) A[k++] = {qid[uarc[j].to], uarc[j].w};
            // highest targets first: the stall test then tends to fire on the first arc
            if (ASORT) std::sort(A.begin() + H[i], A.begin() + H[i + 1], [](const Arc& x, const Arc& y) { return x.to > y.to; });
        }
        vector<Arc>().swap(uarc);
        dd.assign(n, DD{{INF, INF}});
    }
    uint64_t qsettle = 0, qrel = 0, qes = 0, qlook = 0;
    void upsearch(int32_t s, vector<D>& d1, vector<int32_t>& tt, Heap4<D>& h, const vector<D>& d2, D& mu, vector<int32_t>& ent) {
        const int32_t base = ncontracted;
        h.clear(); d1[s] = 0; tt.push_back(s); h.push(0, s);
        while (!h.empty()) {
            auto it = h.pop();
            int32_t u = it.v; D d = it.d;
            if (d > d1[u]) continue;
            if (d >= mu) break;
            ++qsettle;
            if (d2[u] < INF && d + d2[u] < mu) mu = d + d2[u];
            if (u >= base) { ent.push_back(u); continue; }
            bool stalled = false;
            for (uint32_t k = H[u]; k < H[u + 1]; ++k) { const Arc& a = A[k]; if (d1[a.to] < INF && d1[a.to] + a.w < d) { stalled = true; break; } }
            if (stalled) continue;
            qrel += H[u + 1] - H[u];
            for (uint32_t k = H[u]; k < H[u + 1]; ++k) {
                const Arc& a = A[k]; D nd = d + a.w;
                if (nd < d1[a.to]) { if (d1[a.to] == INF) tt.push_back(a.to); d1[a.to] = nd; h.push(nd, a.to); }
            }
        }
    }
    vector<int32_t> eF, eB;
    struct CE { D d; int32_t i; };
    vector<CE> cf, cb;
    // one step of an upward search; returns false when that side is exhausted
    // One settle of the upward search on side `sd` (0 forward, 1 backward).
    template <int sd, class HP>
    inline void step(vector<int32_t>& tt, HP& h, D& mu, vector<int32_t>& ent) {
        auto it = h.pop();
        int32_t u = it.v; D d = it.d;
        DD* X = dd.data();
        if (d > X[u].d[sd]) return;
        ++qsettle;
        { D o = X[u].d[1 - sd]; if (o < INF && d + o < mu) mu = d + o; }
        if (u >= ncontracted) { ent.push_back(u); return; }
        const uint32_t e0 = H[u], e1 = H[u + 1];
        if (STALL) for (uint32_t k = e0; k < e1; ++k) { const Arc& a = A[k]; if (X[a.to].d[sd] + a.w < d) return; }   // stall (INF sums stay large)
        qrel += e1 - e0;
        const int32_t base = ncontracted;
        for (uint32_t k = e0; k < e1; ++k) {
            const Arc& a = A[k]; D nd = d + a.w;
            D& t = X[a.to].d[sd];
            if (nd < t) {
                if (t == INF) { tt.push_back(a.to); if (a.to >= base) ent.push_back(a.to); }
                t = nd;
                if (a.to < base) h.push(nd, a.to);      // core vertices are never expanded
            }
        }
    }
    __attribute__((noinline)) int64_t query(int32_t s, int32_t t) {
        if (s == t) return 0;
        s = qid[s]; t = qid[t];
        D mu = INF;
        eF.clear(); eB.clear();
        hF.clear(); hB.clear();
        dd[s].d[0] = 0; tF.push_back(s); hF.push(0, s);
        dd[t].d[1] = 0; tB.push_back(t); hB.push(0, t);
        for (;;) {
            bool f = !hF.empty() && hF.top_key() < mu, b = !hB.empty() && hB.top_key() < mu;
            if (!f && !b) break;
            if (f && (!b || hF.top_key() <= hB.top_key())) step<0>(tF, hF, mu, eF);
            else step<1>(tB, hB, mu, eB);
        }
        qes += eF.size() + eB.size();
        std::chrono::steady_clock::time_point tc0;
        if (TIMEC) tc0 = std::chrono::steady_clock::now();
        if (!eF.empty() && !eB.empty()) {
            const int32_t base = ncontracted;
            cb.clear();
            for (int32_t b : eB) if (dd[b].d[1] < mu) cb.push_back({dd[b].d[1], b - base});
            if (!cb.empty()) {
                std::sort(cb.begin(), cb.end(), [](const CE& x, const CE& y) { return x.d < y.d; });
                cf.clear();
                for (int32_t a : eF) if (dd[a].d[0] + cb[0].d < mu) cf.push_back({dd[a].d[0], a - base});
                std::sort(cf.begin(), cf.end(), [](const CE& x, const CE& y) { return x.d < y.d; });
                for (const CE& fa : cf) {
                    if (fa.d + cb[0].d >= mu) break;
                    const D* row = &T[(size_t)fa.i * C];
                    for (const CE& fb : cb) {
                        D base2 = fa.d + fb.d; if (base2 >= mu) break;
                        D c = base2 + row[fb.i]; if (c < mu) mu = c;
                        ++qlook;
                    }
                }
            }
        }
        if (TIMEC) tcomb += std::chrono::duration<double>(std::chrono::steady_clock::now() - tc0).count();
        for (int32_t v : tF) dd[v].d[0] = INF;
        for (int32_t v : tB) dd[v].d[1] = INF;
        tF.clear(); tB.clear();
        return mu >= INF ? -1 : (int64_t)mu;
    }
};

// Plain bidirectional Dijkstra on the (pruned) input graph with a settle budget:
// cheap for the many very short queries of a lattice; returns false when the
// budget runs out (the caller then uses the hierarchy).
template <class D>
struct LocalBidir {
    static constexpr D INF = (D)(~(D)0) >> 2;
    vector<uint32_t> H; vector<int32_t> to; vector<D> w;
    vector<D> g[2]; vector<int32_t> tt[2]; Heap4<D> h[2];
    void build(int32_t n, const vector<int32_t>& a, const vector<int32_t>& b, const vector<D>& ww) {
        H.assign(n + 1, 0);
        for (size_t i = 0; i < a.size(); ++i) { ++H[a[i] + 1]; ++H[b[i] + 1]; }
        for (int32_t v = 0; v < n; ++v) H[v + 1] += H[v];
        to.resize(H[n]); w.resize(H[n]);
        vector<uint32_t> c(H.begin(), H.end() - 1);
        for (size_t i = 0; i < a.size(); ++i) { to[c[a[i]]] = b[i]; w[c[a[i]]++] = ww[i]; to[c[b[i]]] = a[i]; w[c[b[i]]++] = ww[i]; }
        g[0].assign(n, INF); g[1].assign(n, INF);
    }
    bool run(int32_t s, int32_t t, int budget, int64_t& ans) {
        if (s == t) { ans = 0; return true; }
        D mu = INF; int settled = 0; bool ok = true;
        for (int k = 0; k < 2; ++k) h[k].clear();
        g[0][s] = 0; tt[0].push_back(s); h[0].push(0, s);
        g[1][t] = 0; tt[1].push_back(t); h[1].push(0, t);
        while (!h[0].empty() && !h[1].empty()) {
            if (h[0].top_key() + h[1].top_key() >= mu) break;
            if (++settled > budget) { ok = false; break; }
            int k = h[0].top_key() <= h[1].top_key() ? 0 : 1;
            auto it = h[k].pop();
            D d = it.d; int32_t u = it.v;
            if (d > g[k][u]) continue;
            vector<D>& gk = g[k]; const vector<D>& go = g[1 - k];
            for (uint32_t e = H[u]; e < H[u + 1]; ++e) {
                int32_t x = to[e]; D nd = d + w[e];
                if (nd < gk[x]) {
                    if (gk[x] == INF) tt[k].push_back(x);
                    gk[x] = nd; h[k].push(nd, x);
                    if (go[x] != INF && nd + go[x] < mu) mu = nd + go[x];
                }
            }
        }
        for (int k = 0; k < 2; ++k) { for (int32_t v : tt[k]) g[k][v] = INF; tt[k].clear(); }
        if (ok) ans = mu >= INF ? -1 : (int64_t)mu;
        return ok;
    }
};

// Defaults chosen by solve() from the graph's shape.
static int64_t DEF_STOPC = 3000, DEF_LMK = 4, DEF_ASTAR = 1, DEF_LOCALR = 0;

template <class D>
static void run(uint64_t ub) {
    CH<D> ch;
    ch.UB = (D)ub;
    ch.LC = envi("LC", 250); ch.STOPC = (int32_t)std::min<int64_t>(V, envi("STOPC", DEF_STOPC));
    ch.CON_SET = (int)envi("CON_SET", 200); ch.ASTAR = (int)envi("ASTAR", DEF_ASTAR); ch.PMODE = (int)envi("PMODE", 0); ch.QORD = (int)envi("QORD", 1); ch.STALL = (int)envi("STALL", 1); ch.TIMEC = (int)envi("TIMEC", 0); ch.PED = envi("PED", 0); ch.ASORT = (int)envi("ASORT", 1); ch.SIM2 = (uint32_t)envi("SIM2", 1000000);
    vector<uint8_t> dead(E, 0);
    if (dim && envi("PRUNE3", 1)) {
        auto nbi = [&](int64_t i, int k, int dir) -> int64_t { int64_t c = (i / strd[k]) % side; if (dir > 0) return c + 1 < side ? i + strd[k] : i - c * strd[k]; return c > 0 ? i - strd[k] : i + (side - 1) * strd[k]; };
        auto wd = [&](int64_t i, int k, int dir) -> uint64_t { return dir > 0 ? ew[i * dim + k] : ew[nbi(i, k, -1) * dim + k]; };
        int64_t cut = 0;
        for (int64_t i = 0; i < V; ++i) for (int k = 0; k < dim; ++k) {
            int64_t j = nbi(i, k, 1); uint64_t w = ew[i * dim + k];
            bool dom = false;
            for (int k2 = 0; k2 < dim && !dom; ++k2) if (k2 != k) for (int s2 = -1; s2 <= 1 && !dom; s2 += 2) {
                int64_t a = nbi(i, k2, s2);
                if (wd(i, k2, s2) + wd(a, k, 1) + wd(j, k2, s2) < w) dom = true;
            }
            if (dom) { dead[i * dim + k] = 1; ++cut; }
        }
        if (g_debug) std::fprintf(stderr, "prune3: %lld\n", (long long)cut);
    }
    // vertex renumbering: Morton order of lattice coordinates (locality)
    vector<int32_t> perm(V);
    for (int32_t v = 0; v < V; ++v) perm[v] = v;
    if (dim && envi("MORTON", 1)) {
        vector<std::pair<uint64_t, int32_t>> key(V);
        for (int64_t v = 0; v < V; ++v) {
            uint64_t c[3] = {0, 0, 0};
            for (int k = 0; k < dim; ++k) c[k] = (uint64_t)((v / strd[k]) % side);
            uint64_t m = 0;
            for (int b = 0; b < 21; ++b) for (int k = 0; k < dim; ++k) m |= ((c[k] >> b) & 1ull) << (b * dim + k);
            key[v] = {m, (int32_t)v};
        }
        std::sort(key.begin(), key.end());
        for (int32_t i = 0; i < V; ++i) perm[key[i].second] = i;
    }
    vector<int32_t> deg0(V, 0);
    for (int32_t e = 0; e < E; ++e) if (!dead[e] && eu[e] != ev[e]) { ++deg0[perm[eu[e]]]; ++deg0[perm[ev[e]]]; }
    ch.init(V, deg0);
    const int64_t LOCALR = dim == 2 ? envi("LOCALR", DEF_LOCALR) : 0;
    LocalBidir<D> lb;
    if (LOCALR > 0) {
        vector<int32_t> la, lbv; vector<D> lw;
        for (int32_t e = 0; e < E; ++e) if (!dead[e] && eu[e] != ev[e]) { la.push_back(perm[eu[e]]); lbv.push_back(perm[ev[e]]); lw.push_back((D)ew[e]); }
        lb.build(V, la, lbv, lw);
    }
    const int LOCALB = (int)envi("LOCALB", 300);
    for (int32_t e = 0; e < E; ++e) if (!dead[e] && eu[e] != ev[e]) ch.add_edge(perm[eu[e]], perm[ev[e]], (D)ew[e]);
    tlog("build");
    ch.landmarks((int)std::min<int64_t>(V, envi("LMK", DEF_LMK)));
    if (ch.K && envi("LMPRUNE", 1)) {
        // an edge longer than some landmark detour lies on no shortest path
        int64_t cut = 0;
        for (int32_t u = 0; u < V; ++u) {
            typename CH<D>::Arc* a = ch.nb(u);
            for (uint32_t i = 0; i < ch.adeg[u]; ) {
                int32_t x = a[i].to;
                if (x > u) {
                    const D* lu = &ch.lmd[(size_t)u * ch.K]; const D* lx = &ch.lmd[(size_t)x * ch.K];
                    D ubd = CH<D>::INF; for (int q = 0; q < ch.K; ++q) ubd = std::min<D>(ubd, lu[q] + lx[q]);
                    if (ubd < a[i].w) { ch.remove_arc(x, u); a[i] = a[ch.adeg[u] - 1]; --ch.adeg[u]; ++cut; continue; }
                }
                ++i;
            }
        }
        if (g_debug) std::fprintf(stderr, "landmark prune: %lld\n", (long long)cut);
    }
    tlog("landmarks");
    int ord = (int)envi("ORDER", 0);
    if (ord == 1 && dim) { nd_order(); tlog("nd order"); ch.contract_order(nd_out); }
    else if (ord == 2) ch.contract_rounds();
    else ch.contract_greedy();
    tlog("contract");
    if (g_debug) std::fprintf(stderr, "searches=%llu settles=%llu capped=%llu (%llu) tprio=%.2f tcon=%.2f\n",
        (unsigned long long)ch.searches, (unsigned long long)ch.settles, (unsigned long long)ch.capped, (unsigned long long)ch.capset, ch.tprio, ch.tcon);
    if (g_debug) std::fprintf(stderr, "landmark skip %llu need %llu\n", (unsigned long long)ch.lm_skip, (unsigned long long)ch.lm_need);
    if (g_debug) { for (int b = 0; b < 64; ++b) if (ch.allhist[b]) std::fprintf(stderr, "  limit 2^%d: searches %llu capped %llu\n", b, (unsigned long long)ch.allhist[b], (unsigned long long)ch.caphist[b]); }
    ch.build_core();
    tlog("core table");
    ch.build_query();
    if (g_debug) std::fprintf(stderr, "upward arcs: %u\n", ch.H[V]);
    {
        // process queries in source order (Morton ids): consecutive queries share
        // most of their upward search spaces, which then stay in cache
        vector<std::pair<int64_t, int32_t>> qo(Q);
        for (int32_t i = 0; i < Q; ++i) qo[i] = {((int64_t)perm[qs[i]] << 32) | (uint32_t)perm[qt[i]], i};
        if (envi("QSORT", 1)) std::sort(qo.begin(), qo.end());
        double tl_ = 0, tg_ = 0; int64_t nl = 0, ng = 0, nfast = 0;
        int localonly = (int)envi("LOCALONLY", 0);
        for (int32_t j = 0; j < Q; ++j) {
            int32_t i = qo[j].second;
            if (localonly && dim == 2) { int64_t a = qs[i], b = qt[i]; int64_t dx = std::llabs(a % side - b % side), dy = std::llabs(a / side - b / side); dx = std::min(dx, side - dx); dy = std::min(dy, side - dy); if ((dx + dy <= 80) != (localonly == 1)) continue; }
            if (g_debug) {
                auto t0 = std::chrono::steady_clock::now();
                qans[i] = ch.query(perm[qs[i]], perm[qt[i]]);
                double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                bool loc = false;
                if (dim == 2) { int64_t a = qs[i], b = qt[i]; int64_t dx = std::llabs(a % side - b % side), dy = std::llabs(a / side - b / side); dx = std::min(dx, side - dx); dy = std::min(dy, side - dy); loc = dx + dy <= 80; }
                if (loc) { tl_ += dt; ++nl; } else { tg_ += dt; ++ng; }
            } else {
                if (LOCALR > 0) {
                    int64_t a = qs[i], b = qt[i];
                    int64_t dx = std::llabs(a % side - b % side), dy = std::llabs(a / side - b / side);
                    dx = std::min(dx, side - dx); dy = std::min(dy, side - dy);
                    if (dx + dy <= LOCALR && lb.run(perm[a], perm[b], LOCALB, qans[i])) { ++nfast; continue; }
                }
                qans[i] = ch.query(perm[qs[i]], perm[qt[i]]);
            }
        }
        if (g_debug) std::fprintf(stderr, "no-core queries %lld: %.3f s (%.2f us)  core queries %lld: %.3f s (%.2f us)\n", (long long)nl, tl_, 1e6 * tl_ / std::max<int64_t>(1, nl), (long long)ng, tg_, 1e6 * tg_ / std::max<int64_t>(1, ng));
    }
    if (g_debug) std::fprintf(stderr, "combination time %.3f s ", ch.tcomb);
    if (g_debug) std::fprintf(stderr, "per query: settles %.1f relaxed %.1f entries %.1f lookups %.1f\n", (double)ch.qsettle / Q, (double)ch.qrel / Q, (double)ch.qes / Q, (double)ch.qlook / Q);
    tlog("queries");
}

static bool applicable() { return !directed; }
static uint64_t distance_bound() {
    vector<uint32_t> h(V + 1, 0);
    for (int32_t e = 0; e < E; ++e) { ++h[eu[e] + 1]; ++h[ev[e] + 1]; }
    for (int32_t v = 0; v < V; ++v) h[v + 1] += h[v];
    vector<std::pair<int32_t, uint32_t>> adj(h[V]);
    { vector<uint32_t> c(h.begin(), h.end() - 1);
      for (int32_t e = 0; e < E; ++e) { adj[c[eu[e]]++] = {ev[e], ew[e]}; adj[c[ev[e]]++] = {eu[e], ew[e]}; } }
    vector<uint64_t> d(V, ~0ull);
    Heap4<uint64_t> hp; d[0] = 0; hp.push(0, 0);
    uint64_t ecc = 0; int32_t cnt = 0;
    while (!hp.empty()) {
        auto it = hp.pop(); if (it.d > d[it.v]) continue;
        ecc = it.d; ++cnt;
        for (uint32_t k = h[it.v]; k < h[it.v + 1]; ++k) { uint64_t nd = it.d + adj[k].second; if (nd < d[adj[k].first]) { d[adj[k].first] = nd; hp.push(nd, adj[k].first); } }
    }
    return cnt == V ? 2 * ecc : ~0ull;
}
static void solve() {
    detect();
    uint64_t ub = distance_bound();
    if (g_debug) std::fprintf(stderr, "lattice dim=%d side=%lld UB=%llu\n", dim, (long long)side, (unsigned long long)ub);
    tlog("bound");
    bool small = ub < (1ull << 29);
    if (dim == 3) { DEF_STOPC = 3000; DEF_LMK = 8; DEF_ASTAR = 1; }         // 3-D lattice
    else if (dim == 2 && small) { DEF_STOPC = 2000; DEF_LMK = 4; DEF_ASTAR = 1; DEF_LOCALR = 6; }  // 2-D, 32-bit
    else if (dim == 2) { DEF_STOPC = 3000; DEF_LMK = 0; DEF_ASTAR = 0; }   // 2-D, 64-bit
    else { DEF_STOPC = std::max<int64_t>(1, std::min<int64_t>(3000, V / 8)); DEF_LMK = 4; DEF_ASTAR = 1; }
    if (small) run<uint32_t>(ub); else run<uint64_t>(ub == ~0ull ? (~0ull >> 3) : ub);
}
}  // namespace lat


#include <cmath>
// Pragmas go after all standard headers: a target pragma in front of them
// breaks a build without -march flags (always_inline target mismatch).
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3")
#endif
#if defined(__GNUC__) && !defined(__clang__) && (defined(__x86_64__) || defined(__i386__))
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")
#endif

// ======================================================================= road
// Directed graphs with coordinates (road2d / hugeq families).
//
// On a row-major S x S lattice whose rows/columns carry road classes, the
// fast lines (arterials, highways) cut the lattice into small rectangular
// "cells" of local streets.  Every path leaving a cell crosses its boundary.
//   1. Cells: exact dense min-plus elimination of the interior gives, for
//      every interior vertex, its first-exit distances to the boundary
//      (and first-entry distances from it) and the boundary-to-boundary
//      distances through the interior (added as shortcuts when shorter than
//      any path along the boundary).
//   2. The remaining "network" graph (boundary lines + shortcuts) gets a
//      contraction hierarchy (lazy priority, bounded witness searches) and
//      hub labels built top-down.
//   3. Interior labels = merge of (exit distance + label of exit vertex).
//   4. Queries: label intersection by scatter/gather; a pair inside one cell
//      additionally gets a Dijkstra restricted to that cell's interior.
// Without lattice structure everything is "network": plain CH + labels.
namespace road {

static int env_int(const char* name, int def) {
    const char* e = std::getenv(name);
    return e ? std::atoi(e) : def;
}
static double now_s() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - g_t0).count();
}

// Detects a row-major S x S lattice: V = S*S and (almost) all arcs join
// lattice neighbours (dx, dy in {-1,0,1}).  Returns S or 0.
static int detect_lattice() {
    int64_t S = 1;
    while (S * S < V) ++S;
    if (S * S != V || S < 8) return 0;
    int64_t bad = 0;
    for (int32_t i = 0; i < E; ++i) {
        int32_t a = eu[i], b = ev[i];
        int32_t ax = a % S, ay = a / S, bx = b % S, by = b / S;
        if (std::abs(ax - bx) > 1 || std::abs(ay - by) > 1) ++bad;
    }
    return bad * 50 <= E ? (int)S : 0;
}

// Road class per lattice line from the travel-time / length ratio of the
// arcs running along it: 0 local, 1 arterial, 2 highway.
static void line_classes(int S, vector<int>& rowc, vector<int>& colc) {
    vector<double> rs((size_t)S, 0), cs((size_t)S, 0);
    vector<int> rn((size_t)S, 0), cn((size_t)S, 0);
    for (int32_t i = 0; i < E; ++i) {
        int32_t a = eu[i], b = ev[i];
        double dx = (double)(cx[a] - cx[b]), dy = (double)(cy[a] - cy[b]);
        double len = std::sqrt(dx * dx + dy * dy);
        if (!(len > 0)) continue;
        double r = (double)ew[i] / len;
        if ((b == a + 1 || a == b + 1) && a / S == b / S) {
            rs[a / S] += r; rn[a / S]++;
        } else if (b == a + S || a == b + S) {
            cs[a % S] += r; cn[a % S]++;
        }
    }
    auto classify = [&](const vector<double>& s, const vector<int>& n, vector<int>& out) {
        vector<double> m;
        for (int i = 0; i < S; ++i) if (n[i]) m.push_back(s[i] / n[i]);
        out.assign((size_t)S, 1);
        if (m.empty()) return;
        std::nth_element(m.begin(), m.begin() + m.size() / 2, m.end());
        double med = m[m.size() / 2];
        for (int i = 0; i < S; ++i) {
            if (!n[i]) continue;
            double q = (s[i] / n[i]) / med;
            out[i] = q < 0.2 ? 2 : q < 0.55 ? 1 : 0;
        }
    };
    classify(rs, rn, rowc);
    classify(cs, cn, colc);
}

struct Plan {
    int S = 0;
    vector<int32_t> cellv;    // cell interiors, each in elimination order
    vector<uint32_t> cellst;  // cell start offsets (size = #cells + 1)
    vector<int32_t> cellx;    // 4 extra boundary vertices per cell (rectangle corners; -1 = none)
    vector<int32_t> cell_of;  // vertex -> cell index, -1 for network vertices
};

// 1-D nested dissection order of [lo, hi] (middle last)
static void nd1(int lo, int hi, vector<int>& out) {
    if (lo > hi) return;
    int m = (lo + hi) / 2;
    nd1(lo, m - 1, out); nd1(m + 1, hi, out);
    out.push_back(m);
}
// 2-D nested dissection order of the lattice rectangle [x0,x1] x [y0,y1]
static void nd2(int S, int x0, int x1, int y0, int y1, vector<int32_t>& out) {
    if (x0 > x1 || y0 > y1) return;
    int w = x1 - x0 + 1, h = y1 - y0 + 1;
    if (w * h <= 1) { out.push_back((int32_t)(y0 * S + x0)); return; }
    vector<int> line;
    if (w >= h) {
        int m = (x0 + x1) / 2;
        nd2(S, x0, m - 1, y0, y1, out); nd2(S, m + 1, x1, y0, y1, out);
        nd1(y0, y1, line);
        for (int y : line) out.push_back((int32_t)(y * S + m));
    } else {
        int m = (y0 + y1) / 2;
        nd2(S, x0, x1, y0, m - 1, out); nd2(S, x0, x1, m + 1, y1, out);
        nd1(x0, x1, line);
        for (int x : line) out.push_back((int32_t)(m * S + x));
    }
}

static Plan make_plan() {
    Plan p;
    if (env_int("NOCELLS", 0)) return p;
    int S = detect_lattice();
    if (!S) return p;
    p.S = S;
    vector<int> rowc, colc;
    line_classes(S, rowc, colc);
    if (g_debug) {
        int na = 0, nh = 0;
        for (int i = 0; i < S; ++i) { na += rowc[i] == 1; nh += rowc[i] == 2; }
        std::fprintf(stderr, "lattice S=%d: arterial rows %d, highway rows %d\n", S, na, nh);
    }
    const int max_interior = env_int("CELL_MAXI", 400);
    auto runs = [&](const vector<int>& c) {
        vector<std::pair<int, int>> r;
        for (int i = 0; i < S;) {
            if (c[i] != 0) { ++i; continue; }
            int j = i;
            while (j + 1 < S && c[j + 1] == 0) ++j;
            r.push_back({i, j});
            i = j + 1;
        }
        return r;
    };
    auto xr = runs(colc), yr = runs(rowc);
    vector<int32_t> cellv;
    vector<uint32_t> cellst;
    vector<int32_t> cellx;
    for (auto& yy : yr)
        for (auto& xx : xr) {
            int w = xx.second - xx.first + 1, h = yy.second - yy.first + 1;
            if (w * h > max_interior) continue;
            cellst.push_back((uint32_t)cellv.size());
            nd2(S, xx.first, xx.second, yy.first, yy.second, cellv);
            int cxs[2] = {xx.first - 1, xx.second + 1}, cys[2] = {yy.first - 1, yy.second + 1};
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b) {
                    int X = cxs[a], Y = cys[b];
                    cellx.push_back(X >= 0 && X < S && Y >= 0 && Y < S ? (int32_t)(Y * S + X) : -1);
                }
        }
    cellst.push_back((uint32_t)cellv.size());
    size_t nc = cellst.size() - 1;
    vector<int32_t> cell_of((size_t)V, -1);
    for (size_t c = 0; c < nc; ++c)
        for (uint32_t k = cellst[c]; k < cellst[c + 1]; ++k) cell_of[cellv[k]] = (int32_t)c;
    // Cells must be separated by network vertices: drop cells joined by an arc.
    vector<uint8_t> drop(nc, 0);
    for (int32_t i = 0; i < E; ++i) {
        int32_t a = cell_of[eu[i]], b = cell_of[ev[i]];
        if (a >= 0 && b >= 0 && a != b) { drop[a] = 1; drop[b] = 1; }
    }
    // Boundary size guard (a few stray arcs could make a boundary large).
    {
        vector<uint32_t> cnt(nc, 0);
        for (int32_t i = 0; i < E; ++i) {
            int32_t a = cell_of[eu[i]], b = cell_of[ev[i]];
            if (a >= 0 && b < 0) ++cnt[a];
            if (b >= 0 && a < 0) ++cnt[b];
        }
        for (size_t c = 0; c < nc; ++c) if (cnt[c] > 1200) drop[c] = 1;
    }
    for (size_t c = 0; c < nc; ++c) {
        if (drop[c]) {
            for (uint32_t k = cellst[c]; k < cellst[c + 1]; ++k) cell_of[cellv[k]] = -1;
            continue;
        }
        int32_t id = (int32_t)(p.cellst.size());
        p.cellst.push_back((uint32_t)p.cellv.size());
        for (uint32_t k = cellst[c]; k < cellst[c + 1]; ++k) { p.cellv.push_back(cellv[k]); cell_of[cellv[k]] = id; }
        for (int k = 0; k < 4; ++k) p.cellx.push_back(cellx[4 * c + k]);
    }
    p.cellst.push_back((uint32_t)p.cellv.size());
    if (p.cellst.size() == 1) p.cellst.clear();
    p.cell_of.swap(cell_of);
    return p;
}

template <class DT>
struct Engine {
    static constexpr DT INF = ~DT(0);
    static constexpr DT CAP = (DT)(~DT(0) >> 1);
    bool overflow = false;

    static inline DT sat_add(DT a, DT b) {   // a <= CAP; b may be INF
        DT c = a + b;
        return c < b ? INF : c;
    }

    struct Arc { int32_t to; int32_t hops; DT fw, bw; };   // fw: v->to, bw: to->v
    vector<Arc> pool;                                      // pooled adjacency
    vector<uint32_t> aoff, asz, acap;

    // --- witness search state
    vector<DT> wd;
    vector<uint32_t> wst;
    uint32_t wcur = 0;
    vector<uint32_t> tst;
    uint32_t tcur = 0;
    vector<DT> md;
    vector<uint32_t> mst;
    uint32_t mcur = 0;
    Heap4<DT> hp;
    int con_settle = 300;
    uint64_t n_settle = 0, n_search = 0;
    uint64_t work = 0, work_budget = ~0ull;   // CH effort guard (deterministic)
    bool aborted = false;
    vector<int32_t> lvl;
    int64_t level_coef = 1000;

    // --- network CH (ranks 0..NR-1 over network vertices)
    int32_t NR = 0;
    vector<int32_t> rank_of, order;
    vector<uint32_t> uhead;
    struct UArc { int32_t to; DT fw, bw; };
    vector<UArc> ua;
    vector<UArc> tmp_up;
    vector<uint32_t> tmp_head;

    // --- cells: exits (interior -> boundary) and entries (boundary -> interior)
    vector<uint32_t> exs, ens;            // per vertex start (interior only)
    vector<uint32_t> exn, enn;            // per vertex count
    struct XE { int32_t b; DT d; };
    vector<XE> exl, enl;
    size_t exl_n = 0, enl_n = 0;
    vector<uint8_t> done;
    vector<int32_t> iord;                 // interior vertices, cell by cell
    vector<uint32_t> iost;                // cell starts in iord
    const vector<int32_t>* cell_of = nullptr;

    inline DT addc(DT a, DT b) {
        uint64_t s = (uint64_t)a + (uint64_t)b;
        if (s > (uint64_t)CAP) { overflow = true; return CAP; }
        return (DT)s;
    }
    inline Arc* adj(int32_t v) { return pool.data() + aoff[v]; }

    void push_arc(int32_t u, const Arc& a) {
        if (asz[u] == acap[u]) {
            uint32_t nc = acap[u] ? acap[u] * 2 : 4;
            size_t no = pool.size();
            pool.resize(pool.size() + nc);
            std::memcpy(pool.data() + no, pool.data() + aoff[u], sizeof(Arc) * asz[u]);
            aoff[u] = (uint32_t)no; acap[u] = nc;
        }
        pool[aoff[u] + asz[u]++] = a;
    }
    void add_arc(int32_t u, int32_t x, DT fw, DT bw, int32_t hops) {
        Arc* a = adj(u);
        uint32_t n = asz[u];
        for (uint32_t k = 0; k < n; ++k) {
            if (a[k].to == x) {
                if (fw < a[k].fw) a[k].fw = fw;
                if (bw < a[k].bw) a[k].bw = bw;
                if (hops > a[k].hops) a[k].hops = hops;
                return;
            }
        }
        push_arc(u, {x, hops, fw, bw});
    }

    void build_graph() {
        vector<uint32_t> deg((size_t)V, 0);
        for (int32_t i = 0; i < E; ++i) if (eu[i] != ev[i]) { ++deg[eu[i]]; ++deg[ev[i]]; }
        aoff.assign((size_t)V, 0); asz.assign((size_t)V, 0); acap.assign((size_t)V, 0);
        uint64_t tot = 0;
        for (int32_t v = 0; v < V; ++v) { aoff[v] = (uint32_t)tot; acap[v] = deg[v] + 2; tot += acap[v]; }
        pool.reserve(tot * 2);
        pool.resize(tot);
        for (int32_t i = 0; i < E; ++i) {
            int32_t a = eu[i], b = ev[i];
            if (a == b) continue;
            DT w = (DT)ew[i];
            if (w > CAP) { overflow = true; w = CAP; }
            add_arc(a, b, w, INF, 1);
            add_arc(b, a, INF, w, 1);
        }
        wd.assign((size_t)V, 0); wst.assign((size_t)V, 0);
        tst.assign((size_t)V, 0);
        md.assign((size_t)V, 0); mst.assign((size_t)V, 0);
        lvl.assign((size_t)V, 0);
        done.assign((size_t)V, 0);
    }

    inline void new_search() { if (++wcur == 0) { std::fill(wst.begin(), wst.end(), 0); wcur = 1; } }
    inline void new_targets() { if (++tcur == 0) { std::fill(tst.begin(), tst.end(), 0); tcur = 1; } }
    inline void new_marks() { if (++mcur == 0) { std::fill(mst.begin(), mst.end(), 0); mcur = 1; } }

    // Dijkstra from src over the remaining graph, bounded by limit; stops once
    // `pending` targets (tst == tcur) are settled or con_settle vertices were.
    void wsearch(int32_t src, DT limit, int pending) {
        ++n_search;
        new_search();
        hp.clear();
        wd[src] = 0; wst[src] = wcur;
        hp.push(0, src);
        int settled = 0;
        while (!hp.empty()) {
            auto it = hp.pop();
            int32_t v = it.v; DT d = it.d;
            if (d > wd[v]) continue;
            if (d > limit) break;
            if (tst[v] == tcur) { if (--pending <= 0) break; }
            if (++settled > con_settle) break;
            const Arc* a = adj(v);
            uint32_t n = asz[v];
            for (uint32_t k = 0; k < n; ++k) {
                if (a[k].fw == INF) continue;
                DT nd = d + a[k].fw;
                if (nd > limit) continue;
                int32_t x = a[k].to;
                if (wst[x] != wcur || nd < wd[x]) {
                    wst[x] = wcur; wd[x] = nd;
                    hp.push(nd, x);
                }
            }
        }
        n_settle += (uint64_t)settled;
    }

    inline void mark_out(int32_t u, int32_t skip) {
        new_marks();
        const Arc* b = adj(u);
        uint32_t n = asz[u];
        work += n;
        for (uint32_t k = 0; k < n; ++k)
            if (b[k].fw != INF && b[k].to != skip) { mst[b[k].to] = mcur; md[b[k].to] = b[k].fw; }
    }
    // u->x path of <= 2 arcs (u's out-arcs marked), avoiding skip, cost <= c?
    inline bool two_hop(int32_t x, DT c, int32_t skip) {
        if (mst[x] == mcur && md[x] <= c) return true;
        const Arc* b = adj(x);
        uint32_t n = asz[x];
        work += n;
        for (uint32_t k = 0; k < n; ++k) {
            int32_t y = b[k].to;
            if (b[k].bw != INF && y != skip && mst[y] == mcur && (uint64_t)md[y] + b[k].bw <= (uint64_t)c) return true;
        }
        return false;
    }

    // Priority: level + edge quotient + hop quotient (2-hop witness estimate).
    int64_t priority(int32_t v) {
        const Arc* nb = adj(v);
        uint32_t n = asz[v];
        int64_t added = 0, added_hops = 0, removed = 0, removed_hops = 0;
        for (uint32_t i = 0; i < n; ++i) {
            int c = (nb[i].fw != INF) + (nb[i].bw != INF);
            removed += c; removed_hops += (int64_t)nb[i].hops * c;
        }
        for (uint32_t i = 0; i < n; ++i) {
            if (nb[i].bw == INF) continue;
            mark_out(nb[i].to, v);
            for (uint32_t j = 0; j < n; ++j) {
                if (j == i || nb[j].fw == INF) continue;
                DT c = (DT)std::min<uint64_t>((uint64_t)nb[i].bw + nb[j].fw, (uint64_t)CAP);
                if (!two_hop(nb[j].to, c, v)) { ++added; added_hops += nb[i].hops + nb[j].hops; }
            }
        }
        if (removed == 0) removed = 1;
        if (removed_hops == 0) removed_hops = 1;
        return (int64_t)lvl[v] * level_coef + (1000 * added) / removed + (1000 * added_hops) / removed_hops;
    }

    struct SC { int32_t u, x; DT w; int32_t hops; };
    vector<SC> sc;
    vector<uint32_t> pend;

    void contract(int32_t v) {
        uint32_t n = asz[v];
        Arc* nb = adj(v);
        for (uint32_t i = 0; i < n; ++i) {
            int32_t w = nb[i].to;
            Arc* l = adj(w);
            uint32_t m = asz[w];
            for (uint32_t k = 0; k < m; ++k)
                if (l[k].to == v) { l[k] = l[m - 1]; asz[w] = m - 1; break; }
        }
        sc.clear();
        for (uint32_t i = 0; i < n; ++i) {
            if (nb[i].bw == INF) continue;
            int32_t u = nb[i].to;
            DT wu = nb[i].bw;
            mark_out(u, v);
            pend.clear();
            new_targets();
            DT lim = 0;
            for (uint32_t j = 0; j < n; ++j) {
                if (j == i || nb[j].fw == INF) continue;
                DT c = addc(wu, nb[j].fw);
                if (two_hop(nb[j].to, c, v)) continue;
                pend.push_back(j);
                tst[nb[j].to] = tcur;
                if (c > lim) lim = c;
            }
            if (pend.empty()) continue;
            wsearch(u, lim, (int)pend.size());
            for (uint32_t j : pend) {
                int32_t x = nb[j].to;
                DT c = addc(wu, nb[j].fw);
                if (!(wst[x] == wcur && wd[x] <= c)) sc.push_back({u, x, c, nb[i].hops + nb[j].hops});
            }
        }
        tmp_head.push_back((uint32_t)tmp_up.size());
        for (uint32_t i = 0; i < n; ++i) {
            tmp_up.push_back({nb[i].to, nb[i].fw, nb[i].bw});
            int32_t w = nb[i].to;
            if (lvl[w] < lvl[v] + 1) lvl[w] = lvl[v] + 1;
        }
        asz[v] = 0;
        for (const SC& s : sc) {
            add_arc(s.u, s.x, s.w, INF, s.hops);
            add_arc(s.x, s.u, INF, s.w, s.hops);
        }
        order.push_back(v);
    }

    // ------------------------------------------------ exact cell elimination
    vector<int32_t> loc;
    vector<DT> W, dring, Df, Db;
    vector<int32_t> cl;
    vector<uint32_t> uhd;
    vector<int32_t> ula;
    vector<DT> ulf, ulb;
    vector<int> inl, outl, nbl, rnf, rnb;
    vector<uint32_t> rnfh, rnbh;
    vector<DT> rnfw, rnbw;
    vector<uint64_t> nbm;
    vector<uint32_t> rh;
    struct RA { int j; DT w; };
    vector<RA> ra;
    vector<DT> dringT, Mx;
    vector<DT> nbw[4], Tpad;
    struct XA { int i, j; DT w; };
    vector<XA> xarc;
    vector<int64_t> rkey;
    int S = 0;
    int dense_k = 6;
    uint64_t st_exits = 0, st_entries = 0, st_ring = 0, st_cells = 0;

    // same-cell queries: interior-only distance via the cell's elimination
    // structure (up-down paths whose top vertex is interior)
    vector<uint32_t> cqh, cql;
    struct SCQ { uint32_t qi; DT d; };
    vector<SCQ> scq;                      // interior-only distances of same-cell queries
    vector<DT> upS, upT;
    void same_cell_queries(size_t c, int ni) {
        if (cqh.empty()) return;
        if (upS.size() < (size_t)ni) { upS.assign((size_t)ni, INF); upT.assign((size_t)ni, INF); }
        for (uint32_t q = cqh[c]; q < cqh[c + 1]; ++q) {
            uint32_t qi = cql[q];
            int s0 = loc[qs[qi]], t0 = loc[qt[qi]];
            upS[s0] = 0; upT[t0] = 0;
            int lo = s0 < t0 ? s0 : t0;
            DT best = INF;
            for (int k = lo; k < ni; ++k) {
                DT us = upS[k], ut = upT[k];
                if (us == INF && ut == INF) continue;
                if (us != INF && ut != INF) { DT c2 = us + ut; if (c2 < best) best = c2; }
                for (uint32_t i = uhd[k]; i < uhd[k + 1]; ++i) {
                    int a = ula[i];
                    if (us != INF && ulf[i] != INF) { DT c2 = us + ulf[i]; if (c2 < upS[a]) upS[a] = c2; }
                    if (ut != INF && ulb[i] != INF) { DT c2 = ut + ulb[i]; if (c2 < upT[a]) upT[a] = c2; }
                }
                upS[k] = INF; upT[k] = INF;
            }
            if (best != INF) scq.push_back({qi, best});
        }
    }

    // Processes one cell: I = interior (elimination order), X = extra
    // boundary vertices (may be -1).
    void eliminate_cell(size_t cidx, const int32_t* I, int ni, const int32_t* X, int nx) {
        cl.assign(I, I + ni);
        for (int k = 0; k < ni; ++k) loc[I[k]] = k;
        for (int k = 0; k < ni; ++k) {
            const Arc* a = adj(I[k]);
            for (uint32_t j = 0; j < asz[I[k]]; ++j)
                if (loc[a[j].to] < 0) { loc[a[j].to] = (int32_t)cl.size(); cl.push_back(a[j].to); }
        }
        for (int k = 0; k < nx; ++k)
            if (X[k] >= 0 && !done[X[k]] && loc[X[k]] < 0) { loc[X[k]] = (int32_t)cl.size(); cl.push_back(X[k]); }
        const int n = (int)cl.size(), nb = n - ni;
        if (S > 0 && nb > 1) {
            int x0 = S, x1 = -1, y0 = S, y1 = -1;
            for (int k = 0; k < ni; ++k) {
                int x = I[k] % S, y = I[k] / S;
                x0 = std::min(x0, x); x1 = std::max(x1, x); y0 = std::min(y0, y); y1 = std::max(y1, y);
            }
            const int X0 = x0 - 1, X1 = x1 + 1, Y0 = y0 - 1, Y1 = y1 + 1, Wd = X1 - X0, Ht = Y1 - Y0;
            rkey.clear();
            for (int i = ni; i < n; ++i) {
                int g = cl[i], x = g % S, y = g / S;
                int64_t key;
                if (y == Y0 && x >= X0 && x <= X1) key = x - X0;
                else if (x == X1 && y >= Y0 && y <= Y1) key = Wd + (y - Y0);
                else if (y == Y1 && x >= X0 && x <= X1) key = Wd + Ht + (X1 - x);
                else if (x == X0 && y >= Y0 && y <= Y1) key = 2 * Wd + Ht + (Y1 - y);
                else key = 2 * (int64_t)(Wd + Ht) + 1 + g;
                rkey.push_back((key << 32) | (uint32_t)g);
            }
            std::sort(rkey.begin(), rkey.end());
            for (int i = ni; i < n; ++i) { cl[i] = (int32_t)(rkey[i - ni] & 0xFFFFFFFFu); loc[cl[i]] = i; }
        }
        ++st_cells;
        W.assign((size_t)n * n, INF);
        DT maxw = 0;
        const int words = (n + 63) >> 6;
        nbm.assign((size_t)n * words, 0);
        for (int i = 0; i < n; ++i) W[(size_t)i * n + i] = 0;
        for (int i = 0; i < n; ++i) {
            const Arc* a = adj(cl[i]);
            for (uint32_t j = 0; j < asz[cl[i]]; ++j) {
                int32_t t = loc[a[j].to];
                if (t < 0 || t == i) continue;
                DT* f = &W[(size_t)i * n + t];
                DT* b = &W[(size_t)t * n + i];
                if (a[j].fw < *f) *f = a[j].fw;
                if (a[j].bw < *b) *b = a[j].bw;
                if (a[j].fw != INF && a[j].fw > maxw) maxw = a[j].fw;
                if (a[j].bw != INF && a[j].bw > maxw) maxw = a[j].bw;
                nbm[(size_t)i * words + (t >> 6)] |= 1ull << (t & 63);
                nbm[(size_t)t * words + (i >> 6)] |= 1ull << (i & 63);
            }
        }
        // boundary-only all-pairs distances (Dijkstra over the sparse boundary graph)
        rh.assign((size_t)nb + 1, 0); ra.clear();
        for (int i = 0; i < nb; ++i) {
            const DT* row = &W[(size_t)(ni + i) * n + ni];
            for (int j = 0; j < nb; ++j) if (j != i && row[j] != INF) ra.push_back({j, row[j]});
            rh[i + 1] = (uint32_t)ra.size();
        }
        {
            // boundary arcs between cyclic neighbours in perimeter order, plus extras
            for (int q = 0; q < 4; ++q) nbw[q].assign((size_t)nb, INF);
            xarc.clear();
            for (int i = 0; i < nb; ++i)
                for (uint32_t e = rh[i]; e < rh[i + 1]; ++e) {
                    int j = ra[e].j;
                    DT w = ra[e].w;
                    if (j == (i + 1) % nb) { nbw[0][j] = w; nbw[2][i] = w; }        // i = j-1 -> j ; i -> i+1
                    else if (i == (j + 1) % nb) { nbw[1][j] = w; nbw[3][i] = w; }   // i = j+1 -> j ; i -> i-1
                    else xarc.push_back({i, j, w});
                }
        }
        // boundary-only all-pairs distances: Bellman-Ford rounds of one sweep
        // each way around the perimeter cycle plus the extra arcs
        dring.assign((size_t)nb * nb, INF);
        {
            const DT* wp = nbw[0].data();   // arc (j-1 -> j)
            const DT* wn = nbw[1].data();   // arc (j+1 -> j)
            for (int s0 = 0; s0 < nb; ++s0) {
                DT* D = &dring[(size_t)s0 * nb];
                D[s0] = 0;
                for (int round = 0; round < nb + 1; ++round) {
                    bool ch = false;
                    for (int step = 1; step <= nb; ++step) {
                        int j = s0 + step; if (j >= nb) j -= nb;
                        int i = j == 0 ? nb - 1 : j - 1;
                        DT c = sat_add(wp[j], D[i]);
                        if (c < D[j]) { D[j] = c; ch = true; }
                    }
                    for (int step = 1; step <= nb; ++step) {
                        int j = s0 - step; if (j < 0) j += nb;
                        int i = j + 1 == nb ? 0 : j + 1;
                        DT c = sat_add(wn[j], D[i]);
                        if (c < D[j]) { D[j] = c; ch = true; }
                    }
                    bool chx = false;
                    for (const XA& e : xarc) {
                        DT c = sat_add(e.w, D[e.i]);
                        if (c < D[e.j]) { D[e.j] = c; chx = true; }
                    }
                    if (xarc.empty() || (!ch && !chx)) break;   // a cycle alone is exact after one round
                }
            }
        }
        // bottom-up elimination of the interior (structure kept as bitsets)
        uhd.assign((size_t)ni + 1, 0);
        rnfh.assign((size_t)ni + 1, 0); rnbh.assign((size_t)ni + 1, 0);
        {
            const size_t cu = (size_t)ni * (size_t)ni / 2 + (size_t)ni, cr = (size_t)ni * (size_t)nb + 1;
            if (ula.size() < cu) { ula.resize(cu); ulf.resize(cu); ulb.resize(cu); }
            if (rnf.size() < cr) { rnf.resize(cr); rnfw.resize(cr); rnb.resize(cr); rnbw.resize(cr); }
            if (nbl.size() < (size_t)n) { nbl.resize(n); inl.resize(n); outl.resize(n); }
        }
        size_t nu = 0, nrf = 0, nrb = 0;
        int* NB = nbl.data();
        int* IN = inl.data();
        int* OUT = outl.data();
        for (int k = 0; k < ni; ++k) {
            int nn = 0, nin = 0, nout = 0;
            uhd[k] = (uint32_t)nu;
            rnfh[k] = (uint32_t)nrf; rnbh[k] = (uint32_t)nrb;
            const DT* rowk = &W[(size_t)k * n];
            uint64_t* bk = &nbm[(size_t)k * words];
            const int w0 = (k + 1) >> 6;
            for (int wi = w0; wi < words; ++wi) {
                uint64_t m = bk[wi];
                if (wi == w0) m &= ~0ull << ((k + 1) & 63);
                while (m) {
                    int x = (wi << 6) + __builtin_ctzll(m);
                    m &= m - 1;
                    DT f = rowk[x], b = W[(size_t)x * n + k];
                    if (f == INF && b == INF) continue;
                    NB[nn++] = x;
                    if (x < ni) { ula[nu] = x; ulf[nu] = f; ulb[nu] = b; ++nu; }
                    else {
                        if (f != INF) { rnf[nrf] = x - ni; rnfw[nrf] = f; ++nrf; }
                        if (b != INF) { rnb[nrb] = x - ni; rnbw[nrb] = b; ++nrb; }
                    }
                    if (f != INF) OUT[nout++] = x;
                    if (b != INF) IN[nin++] = x;
                }
            }
            // structural fill: the neighbours of k become pairwise adjacent
            for (int i = 0; i < nn; ++i) {
                uint64_t* bx = &nbm[(size_t)NB[i] * words];
                for (int wi = w0; wi < words; ++wi) bx[wi] |= bk[wi];
            }
            if (!nout || !nin) continue;
            // out-neighbours split into an interior range and a boundary range;
            // each range is updated densely (vectorised) when it is well filled
            int oi = 0;
            while (oi < nout && OUT[oi] < ni) ++oi;
            int rl[2] = {oi ? OUT[0] : 0, oi < nout ? OUT[oi] : 0};
            int rr[2] = {oi ? OUT[oi - 1] + 1 : 0, oi < nout ? OUT[nout - 1] + 1 : 0};
            int rc[2] = {oi, nout - oi};
            int rs[2] = {0, oi};
            bool dn[2];
            for (int q = 0; q < 2; ++q) dn[q] = rc[q] * dense_k >= rr[q] - rl[q];
            for (int ii = 0; ii < nin; ++ii) {
                const int u = IN[ii];
                const DT wu = W[(size_t)u * n + k];
                DT* rowu = &W[(size_t)u * n];
                const DT keep = rowu[u];
                for (int q = 0; q < 2; ++q) {
                    if (!rc[q]) continue;
                    if (dn[q]) {
                        for (int x = rl[q]; x < rr[q]; ++x) {
                            DT c = sat_add(wu, rowk[x]);
                            rowu[x] = c < rowu[x] ? c : rowu[x];
                        }
                    } else {
                        const int* ol = OUT + rs[q];
                        for (int i = 0; i < rc[q]; ++i) {
                            int x = ol[i];
                            DT c = wu + rowk[x];
                            if (c < rowu[x]) rowu[x] = c;
                        }
                    }
                }
                rowu[u] = keep;
            }
        }
        uhd[ni] = (uint32_t)nu;
        rnfh[ni] = (uint32_t)nrf; rnbh[ni] = (uint32_t)nrb;
        // every value computed in this cell is a path of < 2n arcs of this cell
        if ((uint64_t)maxw * (uint64_t)(2 * n + 2) > (uint64_t)CAP) {
            for (size_t i = 0; i < W.size(); ++i) if (W[i] > CAP / 2 && W[i] != INF) { overflow = true; break; }
            for (size_t i = 0; i < dring.size(); ++i) if (dring[i] > CAP / 2 && dring[i] != INF) { overflow = true; break; }
        }
        same_cell_queries(cidx, ni);
        // top-down closures: Tf[v][b] = min over boundary b' of (first-exit
        // distance v->b') + (boundary path b'->b); Tb symmetric for entries.
        // b is a useful exit of v iff Tf[v][b] is not matched by arriving at b
        // along a boundary arc (then Tf[v][b] is the first-exit distance).
        for (int i = 0; i < nb; ++i) { dring[(size_t)i * nb + i] = 0; }
        dringT.resize((size_t)nb * nb);
        for (int i = 0; i < nb; ++i)
            for (int j = 0; j < nb; ++j) dringT[(size_t)j * nb + i] = dring[(size_t)i * nb + j];
        Df.assign((size_t)ni * nb, INF);
        Db.assign((size_t)ni * nb, INF);
        Mx.resize((size_t)nb);
        Tpad.resize((size_t)nb + 2);
        for (int k = ni - 1; k >= 0; --k) {
            const int32_t g = cl[k];
            for (int dir = 0; dir < 2; ++dir) {
                DT* T = dir == 0 ? &Df[(size_t)k * nb] : &Db[(size_t)k * nb];
                const DT* R = dir == 0 ? dring.data() : dringT.data();
                const vector<int>& rl_ = dir == 0 ? rnf : rnb;
                const vector<DT>& rw_ = dir == 0 ? rnfw : rnbw;
                const vector<uint32_t>& rh_ = dir == 0 ? rnfh : rnbh;
                const vector<DT>& uw = dir == 0 ? ulf : ulb;
                const vector<DT>& TT = dir == 0 ? Df : Db;
                for (uint32_t i = uhd[k]; i < uhd[k + 1]; ++i) {
                    const DT w = uw[i];
                    if (w == INF) continue;
                    const DT* ta = &TT[(size_t)ula[i] * nb];
                    for (int b = 0; b < nb; ++b) { DT c = sat_add(w, ta[b]); T[b] = c < T[b] ? c : T[b]; }
                }
                // T is closed under boundary paths; a direct boundary arc k->b0
                // only matters if it improves T[b0]
                for (uint32_t i = rh_[k]; i < rh_[k + 1]; ++i) {
                    const DT w = rw_[i];
                    if (w >= T[rl_[i]]) continue;
                    const DT* r = R + (size_t)rl_[i] * nb;
                    for (int b = 0; b < nb; ++b) { DT c = sat_add(w, r[b]); T[b] = c < T[b] ? c : T[b]; }
                }
                // best arrival at each b over one boundary arc
                DT* M = Mx.data();
                {
                    DT* Tp = Tpad.data();       // Tp[i + 1] = T[i], cyclic ends
                    Tp[0] = T[nb - 1]; Tp[nb + 1] = T[0];
                    for (int b = 0; b < nb; ++b) Tp[b + 1] = T[b];
                    const DT* wa = nbw[dir == 0 ? 0 : 3].data();   // weight of the arc to/from b-1
                    const DT* wb = nbw[dir == 0 ? 1 : 2].data();   // weight of the arc to/from b+1
                    for (int b = 0; b < nb; ++b) {
                        DT c1 = sat_add(wa[b], Tp[b]), c2 = sat_add(wb[b], Tp[b + 2]);
                        M[b] = c1 < c2 ? c1 : c2;
                    }
                    for (const XA& e : xarc) {
                        if (dir == 0) { DT c = sat_add(e.w, T[e.i]); if (c < M[e.j]) M[e.j] = c; }
                        else { DT c = sat_add(e.w, T[e.j]); if (c < M[e.i]) M[e.i] = c; }
                    }
                }
                vector<XE>& L = dir == 0 ? exl : enl;
                size_t& Ln = dir == 0 ? exl_n : enl_n;
                if (L.size() < Ln + (size_t)nb) L.resize(std::max(L.size() * 2, Ln + (size_t)nb + 1024));
                XE* o = L.data() + Ln;
                uint32_t cnt = 0;
                const int32_t* rg = cl.data() + ni;
                for (int b = 0; b < nb; ++b) {
                    o[cnt] = {rg[b], T[b]};
                    cnt += T[b] < M[b];
                }
                uint32_t st = (uint32_t)Ln;
                Ln += cnt;
                if (dir == 0) { exs[g] = st; exn[g] = cnt; st_exits += cnt; }
                else { ens[g] = st; enn[g] = cnt; st_entries += cnt; }
            }
            iord.push_back(g);
            done[g] = 1;
        }
        // detach the interior from the boundary; boundary shortcuts through it
        for (int i = ni; i < n; ++i) {
            int32_t g = cl[i];
            Arc* a = adj(g);
            uint32_t m = asz[g];
            for (uint32_t j = 0; j < m;) {
                int32_t t = loc[a[j].to];
                if (t >= 0 && t < ni) a[j] = a[--m]; else ++j;
            }
            asz[g] = m;
        }
        for (int i = 0; i < nb; ++i) {
            const DT* row = &W[(size_t)(ni + i) * n + ni];
            const DT* dr = &dring[(size_t)i * nb];
            for (int j = 0; j < nb; ++j) {
                if (i == j || row[j] >= dr[j]) continue;
                add_arc(cl[ni + i], cl[ni + j], row[j], INF, 2);
                add_arc(cl[ni + j], cl[ni + i], INF, row[j], 2);
                ++st_ring;
            }
        }
        for (int32_t g : cl) loc[g] = -1;
    }

    void run_cells(const Plan& p) {
        exs.assign((size_t)V, 0); ens.assign((size_t)V, 0);
        exn.assign((size_t)V, 0); enn.assign((size_t)V, 0);
        if (p.cellst.empty()) return;
        loc.assign((size_t)V, -1);
        {
            // bucket same-cell queries by cell
            size_t nc = p.cellst.size() - 1;
            cqh.assign(nc + 1, 0);
            for (int32_t i = 0; i < Q; ++i) {
                int32_t a = p.cell_of[qs[i]];
                if (a >= 0 && a == p.cell_of[qt[i]] && qs[i] != qt[i]) ++cqh[a + 1];
            }
            for (size_t c = 0; c < nc; ++c) cqh[c + 1] += cqh[c];
            cql.resize(cqh[nc]);
            vector<uint32_t> cur(cqh.begin(), cqh.end() - 1);
            for (int32_t i = 0; i < Q; ++i) {
                int32_t a = p.cell_of[qs[i]];
                if (a >= 0 && a == p.cell_of[qt[i]] && qs[i] != qt[i]) cql[cur[a]++] = (uint32_t)i;
            }
        }
        for (size_t c = 0; c + 1 < p.cellst.size(); ++c) {
            iost.push_back((uint32_t)iord.size());
            eliminate_cell(c, p.cellv.data() + p.cellst[c], (int)(p.cellst[c + 1] - p.cellst[c]),
                           p.cellx.data() + 4 * c, 4);
        }
        iost.push_back((uint32_t)iord.size());
        vector<int32_t>().swap(loc);
        if (g_debug) std::fprintf(stderr, "  [%7.3f] cells %llu: exits %llu entries %llu boundary shortcuts %llu\n",
                                  now_s(), (unsigned long long)st_cells, (unsigned long long)st_exits,
                                  (unsigned long long)st_entries, (unsigned long long)st_ring);
    }

    // ------------------------------------------------------ network CH
    void contract_network() {
        order.clear();
        tmp_up.clear(); tmp_head.clear();
        vector<uint8_t> dirty((size_t)V, 0);
        vector<int64_t> pri((size_t)V, 0);
        typedef std::pair<int64_t, int32_t> PQ;
        vector<PQ> init;
        for (int32_t v = 0; v < V; ++v) if (!done[v]) { pri[v] = priority(v); init.push_back({pri[v], v}); }
        std::priority_queue<PQ, vector<PQ>, std::greater<PQ>> pq(std::greater<PQ>(), std::move(init));
        if (g_debug) std::fprintf(stderr, "  [%7.3f] network vertices %zu\n", now_s(), pq.size());
        while (!pq.empty()) {
            if (work + n_settle * 8 > work_budget) { aborted = true; return; }
            auto [p, v] = pq.top(); pq.pop();
            if (done[v] || p != pri[v]) continue;
            if (dirty[v]) {
                dirty[v] = 0;
                int64_t np = priority(v);
                pri[v] = np;
                if (!pq.empty() && np > pq.top().first) { pq.push({np, v}); continue; }
            }
            const Arc* nb = adj(v);
            for (uint32_t i = 0; i < asz[v]; ++i) dirty[nb[i].to] = 1;
            contract(v); done[v] = 1;
        }
        tmp_head.push_back((uint32_t)tmp_up.size());
        vector<Arc>().swap(pool);
        NR = (int32_t)order.size();
        rank_of.assign((size_t)V, -1);
        for (int32_t r = 0; r < NR; ++r) rank_of[order[r]] = r;
        uhead.assign((size_t)NR + 1, 0);
        ua.resize(tmp_up.size());
        for (int32_t r = 0; r < NR; ++r) {
            uhead[r] = tmp_head[r];
            for (uint32_t k = tmp_head[r]; k < tmp_head[r + 1]; ++k) {
                UArc a = tmp_up[k];
                a.to = rank_of[a.to];
                ua[k] = a;
            }
        }
        uhead[NR] = tmp_head[NR];
        vector<UArc>().swap(tmp_up);
        vector<uint32_t>().swap(tmp_head);
        if (g_debug) std::fprintf(stderr, "  [%7.3f] network CH: %d vertices, %u up arcs, searches %llu settles %llu work %llu (budget %llu)\n",
                                  now_s(), NR, uhead[NR], (unsigned long long)n_search, (unsigned long long)n_settle,
                                  (unsigned long long)(work + n_settle * 8), (unsigned long long)work_budget);
    }

    // ------------------------------------------------------------ labels
    struct LE { uint32_t h; DT d; };      // hub = network rank
    struct LBuf {                         // growable POD array
        LE* p = nullptr;
        size_t n = 0, cap = 0;
        ~LBuf() { std::free(p); }
        inline LE* need(size_t k) {
            if (n + k > cap) {
                size_t nc = std::max(cap * 2, n + k + 4096);
                LE* q = (LE*)std::realloc(p, nc * sizeof(LE));
                if (!q) { std::fprintf(stderr, "out of memory\n"); std::exit(1); }
                p = q; cap = nc;
            }
            return p + n;
        }
        inline LE* data() { return p; }
        inline const LE* data() const { return p; }
        inline size_t size() const { return n; }
    };
    vector<uint64_t> fst, bst;            // per vertex id
    vector<uint32_t> flen, blen;
    LBuf fl, bl;
    vector<DT> tmp;
    vector<uint32_t> hcand;
    size_t nhc = 0;

    // appends the pruned candidate set (tmp/hcand) as the label of v
    int stall_int = 1;
    inline void emit(int dir, int32_t v, bool stall = true) {
        LBuf& L = dir == 0 ? fl : bl;
        LE* o = L.need(nhc);
        const uint32_t* H = uhead.data();
        const UArc* A = ua.data();
        uint32_t cnt = 0;
        for (size_t i = 0; i < nhc; ++i) {
            const uint32_t h = hcand[i];
            const DT d = tmp[h];
            bool keep = true;
            if (stall) for (uint32_t k = H[h]; k < H[h + 1]; ++k) {
                DT w = dir == 0 ? A[k].bw : A[k].fw;
                if (w == INF) continue;
                DT tx = tmp[A[k].to];
                if (tx != INF && (uint64_t)tx + w < (uint64_t)d) { keep = false; break; }
            }
            if (keep) {
                if (d > CAP) overflow = true;
                o[cnt++] = {h, d};
            }
        }
        (dir == 0 ? fst : bst)[v] = L.n;
        (dir == 0 ? flen : blen)[v] = cnt;
        L.n += cnt;
        for (size_t i = 0; i < nhc; ++i) tmp[hcand[i]] = INF;
        nhc = 0;
    }
    inline void merge_label(int dir, int32_t u, DT w) {
        const LBuf& L = dir == 0 ? fl : bl;
        const LE* e = L.data() + (dir == 0 ? fst : bst)[u];
        uint32_t m = (dir == 0 ? flen : blen)[u];
        uint32_t* hc = hcand.data();
        for (uint32_t i = 0; i < m; ++i) {
            uint32_t h = e[i].h;
            DT nd = w + e[i].d;
            if (tmp[h] == INF) { hc[nhc++] = h; tmp[h] = nd; }
            else if (nd < tmp[h]) tmp[h] = nd;
        }
    }

    void build_labels() {
        fst.assign((size_t)V, 0); bst.assign((size_t)V, 0);
        flen.assign((size_t)V, 0); blen.assign((size_t)V, 0);
        fl.need((size_t)V * 44); bl.need((size_t)V * 44);
        tmp.assign((size_t)NR, INF);
        hcand.assign((size_t)NR + 1, 0);
        nhc = 0;
        const uint32_t* H = uhead.data();
        const UArc* A = ua.data();
        for (int32_t r = NR - 1; r >= 0; --r) {
            int32_t v = order[r];
            for (int dir = 0; dir < 2; ++dir) {
                tmp[r] = 0; hcand[nhc++] = (uint32_t)r;
                for (uint32_t k = H[r]; k < H[r + 1]; ++k) {
                    DT w = dir == 0 ? A[k].fw : A[k].bw;
                    if (w == INF) continue;
                    merge_label(dir, order[A[k].to], w);
                }
                emit(dir, v);
            }
        }
        tlog("network labels");
        if (g_debug) std::fprintf(stderr, "network label entries: fwd %zu bwd %zu (avg %.1f)\n",
                                  fl.size(), bl.size(), (double)(fl.size() + bl.size()) / (2.0 * std::max(1, NR)));
        // Interior labels, cell by cell, over a cell-local hub index space:
        // the labels of the cell's exit/entry vertices become dense rows over
        // the local hubs; a label is the min over its exits of (d + row).
        vector<int32_t> hloc((size_t)NR, -1), lh;
        vector<int32_t> bloc((size_t)V, -1), bl_used;
        vector<uint8_t> need((size_t)V, 0);      // bit 0: forward label, bit 1: backward
        for (int32_t i = 0; i < Q; ++i) { need[qs[i]] |= 1; need[qt[i]] |= 2; }
        vector<DT> lt;
        vector<DT> dm[2];
        for (size_t c = 0; c + 1 < iost.size(); ++c) {
            lh.clear(); bl_used.clear();
            for (uint32_t i = iost[c]; i < iost[c + 1]; ++i) {
                int32_t v = iord[i];
                for (int dir = 0; dir < 2; ++dir) {
                    if (!(need[v] >> dir & 1)) continue;
                    const XE* x = dir == 0 ? &exl[exs[v]] : &enl[ens[v]];
                    uint32_t m = dir == 0 ? exn[v] : enn[v];
                    for (uint32_t j = 0; j < m; ++j) {
                        int32_t b = x[j].b;
                        if (bloc[b] >= 0) continue;
                        bloc[b] = (int32_t)bl_used.size();
                        bl_used.push_back(b);
                    }
                }
            }
            for (int32_t b : bl_used)
                for (int dir = 0; dir < 2; ++dir) {
                    const LE* e = (dir == 0 ? fl.data() + fst[b] : bl.data() + bst[b]);
                    uint32_t m = dir == 0 ? flen[b] : blen[b];
                    for (uint32_t j = 0; j < m; ++j)
                        if (hloc[e[j].h] < 0) { hloc[e[j].h] = (int32_t)lh.size(); lh.push_back((int32_t)e[j].h); }
                }
            const size_t Hn = lh.size(), Hp = (Hn + 7) & ~(size_t)7;
            for (int dir = 0; dir < 2; ++dir) {
                dm[dir].assign(bl_used.size() * Hp, INF);
                for (size_t bi = 0; bi < bl_used.size(); ++bi) {
                    int32_t b = bl_used[bi];
                    const LE* e = (dir == 0 ? fl.data() + fst[b] : bl.data() + bst[b]);
                    uint32_t m = dir == 0 ? flen[b] : blen[b];
                    DT* row = &dm[dir][bi * Hp];
                    for (uint32_t j = 0; j < m; ++j) row[hloc[e[j].h]] = e[j].d;
                }
            }
            lt.assign(Hp, INF);
            DT mx = 0;
            for (uint32_t i = iost[c]; i < iost[c + 1]; ++i) {
                int32_t v = iord[i];
                for (int dir = 0; dir < 2; ++dir) {
                    if (!(need[v] >> dir & 1)) continue;
                    const XE* x = dir == 0 ? &exl[exs[v]] : &enl[ens[v]];
                    uint32_t m = dir == 0 ? exn[v] : enn[v];
                    DT* T = lt.data();
                    if (m) {
                        const DT w0 = x[0].d;
                        const DT* row0 = &dm[dir][(size_t)bloc[x[0].b] * Hp];
                        for (size_t h = 0; h < Hp; ++h) T[h] = sat_add(w0, row0[h]);
                    } else {
                        for (size_t h = 0; h < Hp; ++h) T[h] = INF;
                    }
                    for (uint32_t j = 1; j < m; ++j) {
                        const DT w = x[j].d;
                        const DT* row = &dm[dir][(size_t)bloc[x[j].b] * Hp];
                        for (size_t h = 0; h < Hp; ++h) { DT c2 = sat_add(w, row[h]); T[h] = c2 < T[h] ? c2 : T[h]; }
                    }
                    LBuf& L = dir == 0 ? fl : bl;
                    LE* o = L.need(Hn);
                    uint32_t cnt = 0;
                    for (size_t h = 0; h < Hn; ++h) {
                        const DT t = T[h];
                        o[cnt] = {(uint32_t)lh[h], t};
                        const bool fin = t != INF;
                        cnt += fin;
                        const DT tf = fin ? t : 0;
                        mx = tf > mx ? tf : mx;
                    }
                    (dir == 0 ? fst : bst)[v] = L.n;
                    (dir == 0 ? flen : blen)[v] = cnt;
                    L.n += cnt;
                }
            }
            if (mx > CAP) overflow = true;
            for (int32_t h : lh) hloc[h] = -1;
            for (int32_t b : bl_used) bloc[b] = -1;
        }
        if (g_debug) std::fprintf(stderr, "label entries: fwd %zu bwd %zu (avg %.1f)\n",
                                  fl.size(), bl.size(), (double)(fl.size() + bl.size()) / (2.0 * V));
    }

    // ------------------------------------------------------------ queries
    void answer() {
        vector<DT> D((size_t)NR, INF);
        // queries in source order: target and original index, sequential
        vector<uint32_t> head((size_t)V + 1, 0);
        for (int32_t i = 0; i < Q; ++i) ++head[qs[i] + 1];
        for (int32_t v = 0; v < V; ++v) head[v + 1] += head[v];
        vector<int32_t> st((size_t)Q + 16, 0);
        vector<uint32_t> si((size_t)Q);
        {
            vector<uint32_t> cur(head.begin(), head.end() - 1);
            for (int32_t i = 0; i < Q; ++i) { uint32_t k = cur[qs[i]]++; st[k] = qt[i]; si[k] = (uint32_t)i; }
        }
        for (int k = 0; k < 16; ++k) st[(size_t)Q + k] = st.empty() || Q == 0 ? 0 : st[(size_t)Q - 1];
        vector<int64_t> sa((size_t)Q);
        const int PD = 8;   // prefetch distance (queries)
        const LE* FL = fl.data();
        const LE* BL = bl.data();
        int32_t nexts = 0;
        for (int32_t s = 0; s < V; ++s) {
            uint32_t a = head[s], b = head[s + 1];
            if (a == b) continue;
            // prefetch the forward label of the next source
            nexts = s + 1;
            while (nexts < V && head[nexts] == head[nexts + 1]) ++nexts;
            if (nexts < V) {
                const char* pf = (const char*)(FL + fst[nexts]);
                __builtin_prefetch(pf); __builtin_prefetch(pf + 64); __builtin_prefetch(pf + 128);
                __builtin_prefetch(pf + 192); __builtin_prefetch(pf + 256);
            }
            const LE* F = FL + fst[s];
            uint32_t nf = flen[s];
            for (uint32_t i = 0; i < nf; ++i) D[F[i].h] = F[i].d;
            for (uint32_t k = a; k < b; ++k) {
                {
                    const char* pb = (const char*)(BL + bst[st[k + PD]]);
                    __builtin_prefetch(pb); __builtin_prefetch(pb + 64); __builtin_prefetch(pb + 128);
                    __builtin_prefetch(pb + 192); __builtin_prefetch(pb + 256);
                }
                int32_t t = st[k];
                if (t == s) { sa[k] = 0; continue; }
                const LE* B = BL + bst[t];
                uint32_t nb = blen[t];
                uint64_t best = ~0ull;
                if (sizeof(DT) == 4) {
                    for (uint32_t i = 0; i < nb; ++i) {
                        uint64_t c = (uint64_t)D[B[i].h] + (uint64_t)B[i].d;
                        best = c < best ? c : best;
                    }
                    if (best >= (uint64_t)INF) best = ~0ull;
                } else {
                    for (uint32_t i = 0; i < nb; ++i) {
                        DT x = D[B[i].h];
                        if (x != INF) {
                            uint64_t c = (uint64_t)x + B[i].d;
                            if (c < best) best = c;
                        }
                    }
                }
                sa[k] = best == ~0ull ? -1 : (int64_t)best;
            }
            for (uint32_t i = 0; i < nf; ++i) D[F[i].h] = INF;
        }
        for (int32_t k = 0; k < Q; ++k) qans[si[k]] = sa[k];
        for (const SCQ& c : scq)
            if (qans[c.qi] < 0 || (uint64_t)c.d < (uint64_t)qans[c.qi]) qans[c.qi] = (int64_t)c.d;
    }
};

// Fallback for graphs a hierarchy does not suit: one Dijkstra per distinct
// source, stopped once all of that source's targets are settled.
static void dijkstra_queries() {
    vector<uint32_t> h((size_t)V + 1, 0);
    for (int32_t i = 0; i < E; ++i) ++h[eu[i] + 1];
    for (int32_t v = 0; v < V; ++v) h[v + 1] += h[v];
    vector<int32_t> to((size_t)E);
    vector<uint32_t> w((size_t)E);
    {
        vector<uint32_t> cur(h.begin(), h.end() - 1);
        for (int32_t i = 0; i < E; ++i) { uint32_t k = cur[eu[i]]++; to[k] = ev[i]; w[k] = ew[i]; }
    }
    vector<uint32_t> qh((size_t)V + 1, 0), idx((size_t)Q);
    for (int32_t i = 0; i < Q; ++i) ++qh[qs[i] + 1];
    for (int32_t v = 0; v < V; ++v) qh[v + 1] += qh[v];
    {
        vector<uint32_t> cur(qh.begin(), qh.end() - 1);
        for (int32_t i = 0; i < Q; ++i) idx[cur[qs[i]]++] = (uint32_t)i;
    }
    const uint64_t INF64 = ~0ull;
    vector<uint64_t> dist((size_t)V, INF64);
    vector<uint32_t> want((size_t)V, 0);
    vector<int32_t> touched;
    Heap4<uint64_t> heap;
    for (int32_t s = 0; s < V; ++s) {
        if (qh[s] == qh[s + 1]) continue;
        int pending = 0;
        for (uint32_t k = qh[s]; k < qh[s + 1]; ++k) { int32_t t = qt[idx[k]]; if (!want[t]++) ++pending; }
        heap.clear();
        dist[s] = 0; touched.push_back(s); heap.push(0, s);
        while (!heap.empty() && pending > 0) {
            auto it = heap.pop();
            if (it.d > dist[it.v]) continue;
            if (want[it.v]) --pending;
            for (uint32_t k = h[it.v]; k < h[it.v + 1]; ++k) {
                uint64_t nd = it.d + w[k];
                if (nd < dist[to[k]]) {
                    if (dist[to[k]] == INF64) touched.push_back(to[k]);
                    dist[to[k]] = nd; heap.push(nd, to[k]);
                }
            }
        }
        for (uint32_t k = qh[s]; k < qh[s + 1]; ++k) {
            int32_t t = qt[idx[k]];
            qans[idx[k]] = dist[t] == INF64 ? -1 : (int64_t)dist[t];
            want[t] = 0;
        }
        for (int32_t v : touched) dist[v] = INF64;
        touched.clear();
    }
}

template <class DT>
static bool run(const Plan& plan, bool& aborted) {
    Engine<DT> en;
    en.con_settle = env_int("CON_SETTLE", 300);
    en.level_coef = env_int("LEVEL_COEF", 1000);
    en.stall_int = env_int("STALL_INT", 1);
    en.cell_of = &plan.cell_of;
    en.S = plan.S;
    en.dense_k = env_int("DENSE_K", 6);
    en.build_graph();
    tlog("graph built");
    en.run_cells(plan);
    tlog("cells");
    if (en.overflow) return false;
    en.work_budget = (uint64_t)env_int("CH_BUDGET", 60) * ((uint64_t)V + (uint64_t)E) + 20000000ull;
    en.contract_network();
    tlog("network CH");
    if (en.aborted) { tlog("CH effort budget exceeded"); aborted = true; return false; }
    if (en.overflow) return false;
    en.build_labels();
    tlog("labels");
    if (en.overflow) return false;
    en.answer();
    tlog("queries");
    return true;
}

bool applicable() { return directed && has_coords; }

void solve() {
    Plan plan = make_plan();
    tlog("plan");
    bool aborted = false;
    if (env_int("FORCE64", 0) || !run<uint32_t>(plan, aborted)) {
        if (!aborted) {
            tlog("32-bit overflow, redo with 64-bit");
            run<uint64_t>(plan, aborted);
        }
        if (aborted) { dijkstra_queries(); tlog("dijkstra queries"); }
    }
}

}  // namespace road


// After the standard headers (a target pragma in front of them breaks GCC's
// always_inline allocator helpers when compiled without -march), and scoped
// to this namespace with push/pop so it does not leak into other code.
#pragma GCC push_options
#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")

// ============================================================ scalefree
// Directed graphs without coordinates (the power-law family: out-degree 3,
// power-law in-degree, a few percent sinks).
//
// Hubs H = the K vertices of highest in-degree.  Forward searches reach them
// within a few hops, so every s-t path either
//   (a) has no hub in its interior, or
//   (b) has a first hub h, reached by a hub-free prefix.
// A table DIN[x][k] = d(hub k, x) (exact, all x) closes case (b): the forward
// search treats hubs as leaves and a settled hub h gives the candidate
// f(h) + d(h, t).  Case (a) is an ordinary bidirectional meeting, with the
// backward search never entering a hub.  The same table gives landmark lower
// bounds d(x,t) >= d(h,t) - d(h,x), used to prune the forward search.
//
// The table is filled by K-lane Bellman-Ford sweeps over the vertices in
// Dijkstra order from the top hub (~10 sweeps, vectorised).  It is kept in
// 16 bits with saturation; if some finite distance does not fit (or the hubs
// are not mutually reachable) an exact 32-bit table is built as well.
// Graphs whose distances might not fit 31 bits use a plain 64-bit
// bidirectional Dijkstra instead (solve_wide).
namespace sf {

typedef uint32_t D;
static const D INF = 0x7FFFFFFFu;   // every finite distance is < INF (checked); INF+INF fits in 32 bits

// 32-byte vertex record: the three lightest arcs inline, the rest in `ovf`.
struct Rec { uint32_t off; uint32_t deg; int32_t to[3]; D w[3]; };
struct Arc { int32_t to; D w; };
struct Graph { vector<Rec> r; vector<Arc> ovf; };

static Graph GF, GB;                 // out-arcs / in-arcs
// In-trees of in-degree-1 vertices: a "chain" vertex has exactly one in-arc
// (from its parent) and is not a hub; every path into it passes through its
// root, the nearest non-chain ancestor, at chain distance cdelta.  The
// backward search runs on GQ, where each in-arc (x -> y) from a chain vertex
// x is replaced by (root(x) -> y) with weight cdelta(x) + w: chain vertices
// are never visited backwards.  A chain source s pre-expands its own chain
// subtree forwards (those arcs are not represented from s in GQ).
static Graph GQ;
static vector<uint64_t> chainbits;
static int CHAIN = 1;
static vector<int32_t> CROOT;          // chain root (valid for chain vertices)
static vector<D> CDL;                 // chain distance from the root
static vector<uint64_t> sinkbits;    // out-degree 0
static vector<uint64_t> hubbits;
static inline bool bit(const vector<uint64_t>& b, int32_t v) { return (b[(size_t)v >> 6] >> (v & 63)) & 1; }

static int K = 64;                   // hubs
static int LBL = 64;                 // landmark lanes used for forward pruning
static const int SWPD = 16;          // sweep prefetch distance (vertices)
static vector<int32_t> hubs, hid;
static vector<D> DIN;                // DIN[x*K + k] = d(hub k, x)
static vector<uint16_t> T16;         // T16[x*K + k] = min(d(hub k, x), 65535)
static bool EXACT16 = false;         // T16 is exact (65535 = unreachable); DIN unused
static int32_t NREACH = 0;           // vertices reachable from hub 0 (ids 0..NREACH-1)
static D ROWT[512];                  // d(hub k, t) for the current query

bool applicable() { return directed && !has_coords; }

static void build_graph(const vector<int32_t>& a, const vector<int32_t>& b, Graph& g, const vector<uint32_t>& wt = ew) {
    const int32_t m = (int32_t)a.size();
    vector<uint32_t> off((size_t)V + 1, 0);
    for (int32_t i = 0; i < m; ++i) if (a[i] != b[i]) ++off[a[i] + 1];
    for (int32_t v = 0; v < V; ++v) off[v + 1] += off[v];
    vector<Arc> tmp(off[V]);
    {
        vector<uint32_t> cur(off.begin(), off.end() - 1);
        for (int32_t i = 0; i < m; ++i) if (a[i] != b[i]) tmp[cur[a[i]]++] = {b[i], (D)wt[i]};
    }
    g.r.assign((size_t)V, Rec{});
    size_t novf = 0;
    for (int32_t v = 0; v < V; ++v) { uint32_t d = off[v + 1] - off[v]; if (d > 3) novf += d - 3; }
    g.ovf.assign(novf, Arc{});
    size_t o = 0;
    for (int32_t v = 0; v < V; ++v) {
        Arc* s = tmp.data() + off[v];
        uint32_t d = off[v + 1] - off[v];
        std::sort(s, s + d, [](const Arc& x, const Arc& y) { return x.w != y.w ? x.w < y.w : x.to < y.to; });
        Rec& R = g.r[v];
        R.deg = d; R.off = (uint32_t)o;
        for (uint32_t j = 0; j < 3; ++j) { R.to[j] = j < d ? s[j].to : 0; R.w[j] = j < d ? s[j].w : INF; }
        for (uint32_t j = 3; j < d; ++j) g.ovf[o++] = s[j];
    }
}

static inline void arc_at(const Rec& R, const Arc* OV, uint32_t j, int32_t& x, D& w) {
    if (j < 3) { x = R.to[j]; w = R.w[j]; } else { x = OV[R.off + j - 3].to; w = OV[R.off + j - 3].w; }
}

// ------------------------------------------------------------ search state
// One open-addressing table holds both directions' tentative distances; heap
// items carry the table slot, so a pop needs no hash probe.
struct Slot { int32_t v; D d[2]; uint32_t pad; };
struct Table {
    vector<Slot> s;
    vector<uint32_t> used;
    uint32_t mask = 0;
    void init(uint32_t cap) { s.assign(cap, Slot{-1, {INF, INF}, 0}); mask = cap - 1; used.clear(); used.reserve(cap / 2 + 16); }
    inline uint32_t find(int32_t v) {          // inserts if absent
        uint32_t i = ((uint32_t)v * 0x9E3779B1u) >> 8 & mask;
        for (;;) {
            int32_t k = s[i].v;
            if (k == v) return i;
            if (k == -1) { s[i].v = v; used.push_back(i); return i; }
            i = (i + 1) & mask;
        }
    }
    inline D get(int32_t v, int side) const {  // no insertion
        uint32_t i = ((uint32_t)v * 0x9E3779B1u) >> 8 & mask;
        for (;;) {
            int32_t k = s[i].v;
            if (k == v) return s[i].d[side];
            if (k == -1) return INF;
            i = (i + 1) & mask;
        }
    }
    void clear() { for (uint32_t i : used) s[i] = Slot{-1, {INF, INF}, 0}; used.clear(); }
    bool full() const { return used.size() * 2 > mask; }
};

// 4-ary min-heap of packed (key << 32 | payload); the array is padded with ~0
// so all four children of a node are compared without bounds checks.
struct PHeap {
    vector<uint64_t> a;
    uint32_t n = 0;
    PHeap() { a.assign(1024, ~0ull); }
    void clear() { for (uint32_t i = 0; i < n; ++i) a[i] = ~0ull; n = 0; }
    inline bool empty() const { return n == 0; }
    inline D topkey() const { return (D)(a[0] >> 32); }
    inline uint32_t toppay() const { return (uint32_t)a[0]; }
    inline void push(D k, uint32_t pay) {
        if (n + 8 >= a.size()) a.resize(a.size() * 2, ~0ull);
        uint64_t x = (uint64_t)k << 32 | pay;
        uint32_t i = n++;
        while (i) {
            uint32_t p = (i - 1) >> 2;
            if (a[p] <= x) break;
            a[i] = a[p];
            i = p;
        }
        a[i] = x;
    }
    inline uint64_t pop() {
        uint64_t best = a[0];
        uint64_t x = a[--n];
        a[n] = ~0ull;
        uint32_t i = 0;
        for (;;) {
            uint32_t c = 4 * i + 1;
            if (c >= n) break;
            uint64_t m0 = a[c], m1 = a[c + 1], m2 = a[c + 2], m3 = a[c + 3];
            uint32_t j0 = m1 < m0 ? c + 1 : c;       uint64_t v0 = m1 < m0 ? m1 : m0;
            uint32_t j1 = m3 < m2 ? c + 3 : c + 2;   uint64_t v1 = m3 < m2 ? m3 : m2;
            uint32_t j = v1 < v0 ? j1 : j0;          uint64_t v = v1 < v0 ? v1 : v0;
            if (v >= x) break;
            a[i] = v;
            i = j;
        }
        if (n) a[i] = x;
        return best;
    }
    void remap(const vector<uint32_t>& r) { for (uint32_t i = 0; i < n; ++i) a[i] = (a[i] & 0xFFFFFFFF00000000ull) | r[(uint32_t)a[i]]; }
};

static Table TB;
static uint32_t TB_DEFAULT = 1u << 13;
static PHeap hF, hB;
static uint64_t st_settleF = 0, st_settleB = 0, st_search = 0, st_prune = 0;

static void rehash() {
    vector<Slot> old = TB.s;
    vector<uint32_t> oldused = TB.used;
    TB.init((TB.mask + 1) * 2);
    vector<uint32_t> remap(old.size(), 0);
    for (uint32_t i : oldused) { uint32_t j = TB.find(old[i].v); TB.s[j].d[0] = old[i].d[0]; TB.s[j].d[1] = old[i].d[1]; remap[i] = j; }
    hF.remap(remap); hB.remap(remap);
}

// ------------------------------------------------------------ query
// HUBS = false: plain bidirectional Dijkstra (no table).
// Pruning: an arc is skipped when its tentative value plus the other side's
// smallest key reaches mu (backward always; forward only without hubs, since
// with hubs the forward search also closes paths through the table).
template <bool HUBS>
static D query(int32_t s, int32_t t) {
    ++st_search;
    hF.clear(); hB.clear();
    { uint32_t a = TB.find(t); TB.s[a].d[1] = 0; hB.push(0, a); }
    D mu = INF;
    uint64_t wF = 0, wB = 0;
    const Rec* RF = GF.r.data();
    const Rec* RB = GQ.r.data();
    const Arc* OF = GF.ovf.data();
    const Arc* OB = GQ.ovf.data();
    if (CHAIN && bit(chainbits, s)) {
        // settle s's chain subtree (distances along the tree are exact)
        { uint32_t a = TB.find(s); TB.s[a].d[0] = 0; }
        static vector<int32_t> stk;
        stk.clear(); stk.push_back(s);
        while (!stk.empty()) {
            int32_t c = stk.back(); stk.pop_back();
            D dc = TB.s[TB.find(c)].d[0];
            const Rec& R = RF[c];
            wF += R.deg + 1;
            for (uint32_t j = 0; j < R.deg; ++j) {
                int32_t x; D w;
                arc_at(R, OF, j, x, w);
                D nd = dc + w;
                if (bit(sinkbits, x)) { if (x == t && nd < mu) mu = nd; continue; }
                uint32_t sx = TB.find(x);
                Slot& S = TB.s[sx];
                if (nd < S.d[0]) {
                    S.d[0] = nd;
                    D o = S.d[1];
                    if (o != INF && nd + o < mu) mu = nd + o;
                    if (bit(chainbits, x)) stk.push_back(x);   // child in the subtree
                    else hF.push(nd, sx);
                    if (TB.full()) rehash();
                }
            }
        }
    } else {
        uint32_t a = TB.find(s); TB.s[a].d[0] = 0; hF.push(0, a);
    }
    const D* rowt = ROWT;
    D minin = INF;                  // min over hubs h != t of d(h, t)
    if (HUBS) {
        if (EXACT16) {
            const uint16_t* r = &T16[(size_t)t * K];
            for (int k = 0; k < K; ++k) ROWT[k] = r[k] == 0xFFFF ? INF : (D)r[k];
        } else {
            const D* r = &DIN[(size_t)t * K];
            for (int k = 0; k < K; ++k) ROWT[k] = r[k];
        }
        for (int k = 0; k < K; ++k) { D x = hubs[k] != t ? rowt[k] : INF; minin = x < minin ? x : minin; }
    }
    for (;;) {
        if (hF.empty()) break;
        D kf = hF.topkey();
        bool meet = false;
        D kb = INF;
        if (!hB.empty()) { kb = hB.topkey(); meet = kf + kb < mu; }
        bool fwd;
        if (meet) {
            // balance the scanned arcs, counting the arcs the next settle would scan
            int32_t uf = TB.s[hF.toppay()].v, ub = TB.s[hB.toppay()].v;
            fwd = wF + RF[uf].deg <= wB + RB[ub].deg;
        } else {
            if (!HUBS || kf + minin >= mu) break;   // only first-hub candidates remain
            fwd = true;
        }
        if (fwd) {
            uint64_t it = hF.pop();
            D d = (D)(it >> 32);
            Slot& U = TB.s[(uint32_t)it];
            if (d > U.d[0]) continue;
            int32_t u = U.v;
            ++st_settleF;
            if (HUBS) {
                if (bit(hubbits, u)) {                     // leaf: close through the table
                    D c = d + rowt[hid[u]];
                    mu = c < mu ? c : mu;
                    continue;
                }
                if (mu != INF) {                           // landmark bound on d(u, t)
                    int32_t lb = 0;
                    const uint16_t* ru = &T16[(size_t)u * K];
                    for (int k = 0; k < LBL; ++k) {
                        int32_t v = ru[k] != 0xFFFF ? (int32_t)(rowt[k] - (D)ru[k]) : 0;
                        lb = v > lb ? v : lb;
                    }
                    if (d + (D)lb >= mu) { ++st_prune; continue; }
                }
            }
            const Rec& R = RF[u];
            uint32_t deg = R.deg;
            wF += deg + 1;
            for (uint32_t j = 0; j < deg; ++j) {
                int32_t x; D w;
                arc_at(R, OF, j, x, w);
                D nd = d + w;
                if (nd >= mu) break;
                if (bit(sinkbits, x)) {                  // dead end: only useful if it is t
                    if (x == t && nd < mu) mu = nd;
                    continue;
                }
                if (HUBS && bit(hubbits, x)) {           // hub leaf: close at once, never queued
                    D c = nd + rowt[hid[x]];
                    mu = c < mu ? c : mu;
                    continue;
                }
                if (!HUBS && nd + kb >= mu && !bit(chainbits, x)) {
                    // x cannot lie on a better path unless it is already met
                    D o = TB.get(x, 1);
                    if (o != INF && nd + o < mu) mu = nd + o;
                    continue;
                }
                uint32_t sx = TB.find(x);
                Slot& S = TB.s[sx];
                if (nd < S.d[0]) {
                    S.d[0] = nd;
                    D o = S.d[1];
                    D c = nd + o;
                    if (o != INF && c < mu) mu = c;
                    __builtin_prefetch(&RF[x]);
                    if (HUBS) {                          // its landmark row, read when settled
                        const char* q = (const char*)&T16[(size_t)x * K];
                        for (int b = 0; b < LBL * 2; b += 64) __builtin_prefetch(q + b);
                    }
                    hF.push(nd, sx);
                    if (TB.full()) rehash();
                }
            }
        } else {
            uint64_t it = hB.pop();
            D d = (D)(it >> 32);
            Slot& U = TB.s[(uint32_t)it];
            if (d > U.d[1]) continue;
            int32_t u = U.v;
            ++st_settleB;
            const Rec& R = RB[u];
            uint32_t deg = R.deg;
            wB += deg + 1;
            for (uint32_t j = 0; j < deg; ++j) {
                int32_t x; D w;
                arc_at(R, OB, j, x, w);
                D nd = d + w;
                if (nd >= mu) break;
                if (HUBS && bit(hubbits, x)) continue;       // paths through hubs are case (b)
                if (nd + kf >= mu) {                         // no push, but it may meet
                    D o = TB.get(x, 0);
                    if (o != INF && nd + o < mu) mu = nd + o;
                    continue;
                }
                uint32_t sx = TB.find(x);
                Slot& S = TB.s[sx];
                if (nd < S.d[1]) {
                    S.d[1] = nd;
                    D o = S.d[0];
                    D c = nd + o;
                    if (o != INF && c < mu) mu = c;
                    __builtin_prefetch(&RB[x]);
                    hB.push(nd, sx);
                    if (TB.full()) rehash();
                }
            }
        }
    }
    TB.clear();
    if (TB.mask + 1 > TB_DEFAULT) TB.init(TB_DEFAULT);
    return mu;
}

// ------------------------------------------------------------ hub table
static void dijkstra_lane(int k) {
    vector<D> dist((size_t)V, INF);
    PHeap heap;
    int32_t src = hubs[k];
    dist[src] = 0; heap.push(0, (uint32_t)src);
    while (!heap.empty()) {
        uint64_t it = heap.pop();
        D d = (D)(it >> 32);
        int32_t u = (int32_t)(uint32_t)it;
        if (d > dist[u]) continue;
        const Rec& R = GF.r[u];
        for (uint32_t j = 0; j < R.deg; ++j) {
            int32_t x; D w;
            arc_at(R, GF.ovf.data(), j, x, w);
            D nd = d + w;
            if (nd < dist[x]) { dist[x] = nd; heap.push(nd, (uint32_t)x); }
        }
    }
    for (int32_t v = 0; v < V; ++v) DIN[(size_t)v * K + k] = dist[v];
}

// Pull-form Bellman-Ford, all K lanes at once, vertices in id order (= settle
// order from the top hub).  A vertex is re-evaluated only when one of its
// in-neighbours changed after its previous evaluation.  Returns false if it
// did not converge within the sweep budget.
// 16-bit table: T16[x*K + k] = min(d(hub k, x), 65535), same sweeps with
// saturating arithmetic (min and saturating add commute, so every value
// below 65535 is exact).  It is used for the landmark bounds (65535 =
// "unknown", lane skipped) and, when verified (all finite distances from the
// hubs below 65535, hubs mutually reachable), also as the exact table.
static bool sweep16(int max_sweeps) {
    T16.assign((size_t)V * K, 0xFFFF);
    for (int k = 0; k < K; ++k) T16[(size_t)hubs[k] * K + k] = 0;
    // prv: row changed in the previous pass, cur: changed in this pass.  In
    // id order, x must be re-evaluated if an in-neighbour u < x changed in
    // this pass or an in-neighbour u > x changed in the previous one.
    const size_t NW = ((size_t)V + 63) / 64;
    vector<uint64_t> prv(NW, ~0ull), cur(NW, 0);
    auto chgd = [&](int32_t u, int32_t x) {
        const vector<uint64_t>& b = u < x ? cur : prv;
        return (b[(size_t)u >> 6] >> (u & 63)) & 1;
    };
    // Chain vertices are skipped: their rows follow from their roots at the
    // end, and GQ already routes arcs out of chains from the roots.
    const Rec* RB = GQ.r.data();
    const Arc* OV = GQ.ovf.data();
    vector<uint16_t> nr((size_t)K);
    for (int pass = 0; pass < max_sweeps; ++pass) {
        uint64_t nch = 0;
        std::fill(cur.begin(), cur.end(), 0);
        for (int32_t x = 0; x < V; ++x) {
            if (x + SWPD < V) {                    // prefetch the rows a later vertex will pull
                const Rec& P = RB[x + SWPD];
                uint32_t pd = bit(chainbits, x + SWPD) ? 0 : P.deg < 3 ? P.deg : 3;
                for (uint32_t j = 0; j < pd; ++j) {
                    int32_t u = P.to[j];
                    if (chgd(u, x + SWPD)) {
                        const char* q = (const char*)&T16[(size_t)u * K];
                        for (int b = 0; b < K * 2; b += 64) __builtin_prefetch(q + b);
                    }
                }
            }
            const Rec& R = RB[x];
            uint32_t deg = R.deg;
            if (deg == 0 || bit(chainbits, x)) continue;
            bool need = false;
            for (uint32_t j = 0; j < deg && !need; ++j) {
                int32_t u = j < 3 ? R.to[j] : OV[R.off + j - 3].to;
                need = chgd(u, x);
            }
            if (!need) continue;
            uint16_t* row = &T16[(size_t)x * K];
            uint16_t* __restrict q = nr.data();
            for (int k = 0; k < K; ++k) q[k] = row[k];
            for (uint32_t j = 0; j < deg; ++j) {
                int32_t u; D w;
                arc_at(R, OV, j, u, w);
                if (!chgd(u, x)) continue;         // rows only decrease: unchanged ones are already in
                const uint16_t ww = w < 0xFFFF ? (uint16_t)w : (uint16_t)0xFFFF;
                const uint16_t* __restrict ru = &T16[(size_t)u * K];
                for (int k = 0; k < K; ++k) {
                    uint16_t c = (uint16_t)(ru[k] + ww);
                    c = c < ru[k] ? (uint16_t)0xFFFF : c;          // saturate
                    q[k] = c < q[k] ? c : q[k];
                }
            }
            uint32_t diff = 0;
            for (int k = 0; k < K; ++k) diff |= (uint32_t)(q[k] ^ row[k]);
            if (diff) { for (int k = 0; k < K; ++k) row[k] = q[k]; cur[(size_t)x >> 6] |= 1ull << (x & 63); ++nch; }
        }
        prv.swap(cur);
        if (g_debug) std::fprintf(stderr, "  sweep16 %d: changed %llu\n", pass, (unsigned long long)nch);
        if (nch == 0) {
            for (int32_t x = 0; x < V; ++x) {
                if (!bit(chainbits, x)) continue;
                const uint16_t* rr = &T16[(size_t)CROOT[x] * K];
                uint16_t* row = &T16[(size_t)x * K];
                const uint16_t ww = CDL[x] < 0xFFFF ? (uint16_t)CDL[x] : (uint16_t)0xFFFF;
                for (int k = 0; k < K; ++k) {
                    uint16_t c = (uint16_t)(rr[k] + ww);
                    row[k] = c < rr[k] ? (uint16_t)0xFFFF : c;
                }
            }
            return true;
        }
    }
    return false;
}

// Exactness of T16: vertices [0, nreach) are those reachable from hub 0 (ids
// are its settle order).  If every hub reaches hub 0 and is reached from it,
// all hubs reach exactly that set, so 65535 below nreach would mean a
// saturated distance (not allowed) and above it means unreachable.
static bool verify16(int32_t nreach) {
    const int32_t h0 = hubs[0];
    for (int k = 0; k < K; ++k) {
        if (T16[(size_t)hubs[k] * K + 0] == 0xFFFF) return false;
        if (T16[(size_t)h0 * K + k] == 0xFFFF) return false;
    }
    const uint16_t* p = T16.data();
    const size_t n = (size_t)nreach * K;
    uint32_t bad = 0;
    for (size_t i = 0; i < n; ++i) bad |= p[i] == 0xFFFF;
    return bad == 0;
}

static bool sweep_tables(int max_sweeps) {
    std::fill(DIN.begin(), DIN.end(), INF);
    for (int k = 0; k < K; ++k) DIN[(size_t)hubs[k] * K + k] = 0;
    vector<uint32_t> chg((size_t)V, 1), evl((size_t)V, 0);
    uint32_t now = 2;
    const Rec* RB = GB.r.data();
    const Arc* OV = GB.ovf.data();
    vector<D> nr((size_t)K);
    for (int pass = 0; pass < max_sweeps; ++pass) {
        uint64_t nch = 0;
        for (int32_t x = 0; x < V; ++x) {
            if (x + SWPD < V) {                    // prefetch the rows a later vertex will pull
                const Rec& P = RB[x + SWPD];
                uint32_t e2 = evl[x + SWPD];
                uint32_t pd = P.deg < 3 ? P.deg : 3;
                for (uint32_t j = 0; j < pd; ++j) {
                    int32_t u = P.to[j];
                    if (chg[u] > e2) {
                        const char* q = (const char*)&DIN[(size_t)u * K];
                        for (int b = 0; b < K * 4; b += 64) __builtin_prefetch(q + b);
                    }
                }
            }
            const Rec& R = RB[x];
            uint32_t deg = R.deg;
            if (deg == 0) continue;
            uint32_t ex = evl[x];
            bool need = false;
            for (uint32_t j = 0; j < deg && !need; ++j) {
                int32_t u = j < 3 ? R.to[j] : OV[R.off + j - 3].to;
                need = chg[u] > ex;
            }
            evl[x] = ++now;
            if (!need) continue;
            D* row = &DIN[(size_t)x * K];
            D* __restrict q = nr.data();
            for (int k = 0; k < K; ++k) q[k] = row[k];
            for (uint32_t j = 0; j < deg; ++j) {
                int32_t u; D w;
                arc_at(R, OV, j, u, w);
                const D* __restrict ru = &DIN[(size_t)u * K];
                for (int k = 0; k < K; ++k) { D c = ru[k] + w; q[k] = c < q[k] ? c : q[k]; }
            }
            uint32_t diff = 0;
            for (int k = 0; k < K; ++k) diff |= q[k] ^ row[k];
            if (diff) { for (int k = 0; k < K; ++k) row[k] = q[k]; chg[x] = now; ++nch; }
        }
        if (g_debug) std::fprintf(stderr, "  sweep %d: changed %llu\n", pass, (unsigned long long)nch);
        if (nch == 0) return true;
    }
    return false;
}

// Dijkstra order from `src` (unreached vertices appended in id order).
static vector<int32_t> settle_order(int32_t src) {
    // light CSR (unsorted) + monotone radix heap
    vector<uint32_t> off((size_t)V + 1, 0);
    for (int32_t i = 0; i < E; ++i) if (eu[i] != ev[i]) ++off[eu[i] + 1];
    for (int32_t v = 0; v < V; ++v) off[v + 1] += off[v];
    vector<uint64_t> arc(off[V]);                    // (w << 32) | to
    {
        vector<uint32_t> cur(off.begin(), off.end() - 1);
        for (int32_t i = 0; i < E; ++i) if (eu[i] != ev[i]) arc[cur[eu[i]]++] = (uint64_t)ew[i] << 32 | (uint32_t)ev[i];
    }
    vector<D> dist((size_t)V, INF);
    vector<char> done((size_t)V, 0);
    vector<int32_t> ord;
    ord.reserve(V);
    vector<uint64_t> bk[33];
    D last = 0;
    size_t n = 0;
    auto bucket = [&](D k) { return k == last ? 0 : 32 - __builtin_clz(k ^ last); };
    dist[src] = 0; bk[0].push_back((uint64_t)src); ++n;
    while (n) {
        if (bk[0].empty()) {
            int i = 1;
            while (bk[i].empty()) ++i;
            D m = INF;
            for (uint64_t it : bk[i]) { D k = (D)(it >> 32); m = k < m ? k : m; }
            last = m;
            for (uint64_t it : bk[i]) bk[bucket((D)(it >> 32))].push_back(it);
            bk[i].clear();
        }
        uint64_t it = bk[0].back(); bk[0].pop_back(); --n;
        D d = (D)(it >> 32);
        int32_t u = (int32_t)(uint32_t)it;
        if (d > dist[u] || done[u]) continue;
        done[u] = 1; ord.push_back(u);
        for (uint32_t a = off[u]; a < off[u + 1]; ++a) {
            int32_t x = (int32_t)(uint32_t)arc[a];
            D nd = d + (D)(arc[a] >> 32);
            if (nd < dist[x]) { dist[x] = nd; bk[bucket(nd)].push_back((uint64_t)nd << 32 | (uint32_t)x); ++n; }
        }
    }
    NREACH = (int32_t)ord.size();
    for (int32_t v = 0; v < V; ++v) if (!done[v]) ord.push_back(v);
    return ord;
}

// ------------------------------------------------------------ 64-bit fallback
// Plain bidirectional Dijkstra for graphs whose distances may not fit the
// 32-bit fast path.
static void solve_wide() {
    typedef uint64_t W;
    const W WINF = ~0ull;
    vector<uint32_t> fo((size_t)V + 1, 0), bo((size_t)V + 1, 0);
    for (int32_t i = 0; i < E; ++i) if (eu[i] != ev[i]) { ++fo[eu[i] + 1]; ++bo[ev[i] + 1]; }
    for (int32_t v = 0; v < V; ++v) { fo[v + 1] += fo[v]; bo[v + 1] += bo[v]; }
    vector<int32_t> ft(fo[V]), bt(bo[V]);
    vector<W> fw(fo[V]), bw(bo[V]);
    {
        vector<uint32_t> cf(fo.begin(), fo.end() - 1), cb(bo.begin(), bo.end() - 1);
        for (int32_t i = 0; i < E; ++i) if (eu[i] != ev[i]) {
            uint32_t a = cf[eu[i]]++; ft[a] = ev[i]; fw[a] = ew[i];
            uint32_t b = cb[ev[i]]++; bt[b] = eu[i]; bw[b] = ew[i];
        }
    }
    vector<W> dF((size_t)V, WINF), dB((size_t)V, WINF);
    vector<int32_t> tF, tB;
    Heap4<W> HF, HB;
    for (int32_t q = 0; q < Q; ++q) {
        int32_t s = qs[q], t = qt[q];
        if (s == t) { qans[q] = 0; continue; }
        W mu = WINF;
        HF.clear(); HB.clear();
        dF[s] = 0; tF.push_back(s); HF.push(0, s);
        dB[t] = 0; tB.push_back(t); HB.push(0, t);
        while (!HF.empty() && !HB.empty()) {
            W kf = HF.top_key(), kb = HB.top_key();
            if (kf >= mu || kb >= mu || kf + kb >= mu) break;
            bool fwd = tF.size() <= tB.size();
            Heap4<W>& H = fwd ? HF : HB;
            vector<W>& d1 = fwd ? dF : dB;
            const vector<W>& d2 = fwd ? dB : dF;
            vector<int32_t>& tc = fwd ? tF : tB;
            const vector<uint32_t>& O = fwd ? fo : bo;
            const vector<int32_t>& TO = fwd ? ft : bt;
            const vector<W>& WT = fwd ? fw : bw;
            auto it = H.pop();
            if (it.d > d1[it.v]) continue;
            for (uint32_t a = O[it.v]; a < O[it.v + 1]; ++a) {
                int32_t x = TO[a];
                W nd = it.d + WT[a];
                if (nd < d1[x]) {
                    if (d1[x] == WINF) tc.push_back(x);
                    d1[x] = nd;
                    if (d2[x] != WINF && nd + d2[x] < mu) mu = nd + d2[x];
                    H.push(nd, x);
                }
            }
        }
        for (int32_t v : tF) dF[v] = WINF;
        for (int32_t v : tB) dB[v] = WINF;
        tF.clear(); tB.clear();
        qans[q] = mu == WINF ? -1 : (int64_t)mu;
    }
}

// ------------------------------------------------------------ driver
void solve() {
    if (const char* e = std::getenv("SF_K")) K = std::atoi(e);
    if (K > 512) K = 512;                        // ROWT capacity
    if (const char* e = std::getenv("SF_LBL")) LBL = std::atoi(e);
    vector<int32_t> indeg((size_t)V, 0), outdeg((size_t)V, 0);
    vector<D> maxout((size_t)V, 0);
    uint64_t maxw = 0;
    for (int32_t i = 0; i < E; ++i) if (eu[i] != ev[i]) {
        ++outdeg[eu[i]]; ++indeg[ev[i]];
        if (ew[i] > maxout[eu[i]]) maxout[eu[i]] = ew[i];
        if (ew[i] > maxw) maxw = ew[i];
    }
    // A shortest path uses each vertex's out-arc at most once, so its length
    // is below sum_v maxout(v); tentative values add one more arc.
    uint64_t bound = 2 * maxw;
    for (int32_t v = 0; v < V; ++v) bound += maxout[v];
    if (bound >= INF) { solve_wide(); tlog("queries (64-bit)"); return; }

    // hubs: highest in-degree among vertices that have both in- and out-arcs
    vector<int32_t> byin;
    for (int32_t v = 0; v < V; ++v) if (outdeg[v] > 0 && indeg[v] > 0) byin.push_back(v);
    std::stable_sort(byin.begin(), byin.end(), [&](int32_t a, int32_t b) { return indeg[a] > indeg[b]; });
    if ((int)byin.size() < K) K = (int)byin.size();
    if (K < 0) K = 0;
    if (LBL > K) LBL = K;

    // Renumber vertices in Dijkstra order from the top hub: the table sweeps
    // then converge in a few passes.
    vector<int32_t> nid((size_t)V);
    for (int32_t v = 0; v < V; ++v) nid[v] = v;
    if (K > 0) {
        vector<int32_t> order = settle_order(byin[0]);
        for (int32_t k = 0; k < V; ++k) nid[order[k]] = k;
        for (int32_t i = 0; i < E; ++i) { eu[i] = nid[eu[i]]; ev[i] = nid[ev[i]]; }
        vector<int32_t> a((size_t)V), b((size_t)V);
        for (int32_t v = 0; v < V; ++v) { a[nid[v]] = indeg[v]; b[nid[v]] = outdeg[v]; }
        indeg.swap(a); outdeg.swap(b);
        for (int32_t& v : byin) v = nid[v];
    }
    build_graph(eu, ev, GF);
    build_graph(ev, eu, GB);
    sinkbits.assign(((size_t)V + 63) / 64, 0);
    hubbits.assign(((size_t)V + 63) / 64, 0);
    for (int32_t v = 0; v < V; ++v) if (outdeg[v] == 0) sinkbits[v >> 6] |= 1ull << (v & 63);
    tlog("graph built");

    if (K > 0) {
        hubs.assign(byin.begin(), byin.begin() + K);
        hid.assign((size_t)V, -1);
        for (int k = 0; k < K; ++k) { hubbits[hubs[k] >> 6] |= 1ull << (hubs[k] & 63); hid[hubs[k]] = k; }
    }

    // chain (in-degree-1) trees and the compressed backward graph GQ
    if (const char* e = std::getenv("SF_CHAIN")) CHAIN = std::atoi(e);
    chainbits.assign(((size_t)V + 63) / 64, 0);
    if (CHAIN) {
        vector<int32_t> par((size_t)V, -1), root((size_t)V);
        vector<D> dl((size_t)V, 0);
        vector<char> st((size_t)V, 0);   // 0 new, 1 on stack, 2 done
        for (int32_t v = 0; v < V; ++v) {
            root[v] = v;
            if (indeg[v] == 1 && !bit(hubbits, v) && GB.r[v].deg == 1) {
                par[v] = GB.r[v].to[0];
                dl[v] = GB.r[v].w[0];      // weight of the parent arc (for now)
            }
        }
        vector<int32_t> path;
        for (int32_t v = 0; v < V; ++v) {
            if (st[v]) continue;
            path.clear();
            int32_t u = v;
            while (par[u] >= 0 && st[u] == 0) { st[u] = 1; path.push_back(u); u = par[u]; }
            if (par[u] >= 0 && st[u] == 1) {
                // cycle of chain vertices: make its members roots
                size_t k = path.size();
                while (path[k - 1] != u) --k;
                for (size_t i = k - 1; i < path.size(); ++i) { par[path[i]] = -1; dl[path[i]] = 0; root[path[i]] = path[i]; st[path[i]] = 2; }
                path.resize(k - 1);
            }
            st[u] = 2;
            for (size_t i = path.size(); i-- > 0;) {
                int32_t w = path[i];
                int32_t p = par[w];
                root[w] = root[p];
                dl[w] += dl[p];
                st[w] = 2;
            }
        }
        for (int32_t v = 0; v < V; ++v) if (par[v] >= 0) chainbits[v >> 6] |= 1ull << (v & 63);
        CROOT = root; CDL = dl;
        vector<int32_t> qa, qb;
        vector<uint32_t> qw;
        qa.reserve(E); qb.reserve(E); qw.reserve(E);
        for (int32_t i = 0; i < E; ++i) {
            int32_t a = eu[i], b = ev[i];
            if (a == b) continue;
            D w = ew[i];
            if (par[a] >= 0) { w += dl[a]; a = root[a]; }
            if (a == b) continue;
            qa.push_back(a); qb.push_back(b); qw.push_back(w);
        }
        build_graph(qb, qa, GQ, qw);
        tlog("chains");
    } else {
        GQ = GB;
    }

    if (K > 0) {
        bool conv = sweep16(40);
        EXACT16 = conv && verify16(NREACH);
        if (!conv) {                               // bounds must still be valid: recompute exactly
            DIN.assign((size_t)V * K, INF);
            if (!sweep_tables(40)) for (int k = 0; k < K; ++k) dijkstra_lane(k);
            for (size_t i = 0; i < DIN.size(); ++i) T16[i] = DIN[i] < 0xFFFF ? (uint16_t)DIN[i] : (uint16_t)0xFFFF;
        } else if (!EXACT16) {
            DIN.assign((size_t)V * K, INF);
            if (!sweep_tables(40)) for (int k = 0; k < K; ++k) dijkstra_lane(k);
        }
        if (g_debug) std::fprintf(stderr, "  table: 16-bit %s\n", EXACT16 ? "exact" : "bounds only");
        tlog("hub table");
    }

    TB.init(TB_DEFAULT);
    // map queries first so the next query's data can be prefetched
    for (int32_t i = 0; i < Q; ++i) { qs[i] = nid[qs[i]]; qt[i] = nid[qt[i]]; }
    for (int32_t i = 0; i < Q; ++i) {
        int32_t s = qs[i], t = qt[i];
        if (i + 1 < Q) {                       // warm the next query's first reads
            int32_t s2 = qs[i + 1], t2 = qt[i + 1];
            __builtin_prefetch(&GF.r[s2]);
            __builtin_prefetch(&GQ.r[t2]);
            if (K > 0) { const char* q = (const char*)&T16[(size_t)t2 * K]; for (int b = 0; b < K * 2; b += 64) __builtin_prefetch(q + b); }
        }
        if (s == t) { qans[i] = 0; continue; }
        if (outdeg[s] == 0 || indeg[t] == 0) { qans[i] = -1; continue; }
        D r = K > 0 ? query<true>(s, t) : query<false>(s, t);
        qans[i] = r >= INF ? -1 : (int64_t)r;
    }
    if (g_debug && st_search) std::fprintf(stderr, "K=%d searches %llu settled fwd %.1f bwd %.1f pruned %.1f\n", K,
                                           (unsigned long long)st_search, (double)st_settleF / st_search,
                                           (double)st_settleB / st_search, (double)st_prune / st_search);
    tlog("queries");
}

}  // namespace sf
#pragma GCC pop_options


#pragma GCC pop_options

// Picks the scheme for the graph at hand (first module that applies).
static void run_queries_all() {
    if (w64::applicable()) { w64::solve(); return; }
    if (lat::applicable()) { lat::solve(); return; }
    if (road::applicable()) { road::solve(); return; }
    if (sf::applicable()) { sf::solve(); return; }
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: %s <graph_file> <query_file> <output_file>\n", argv[0]);
        return 1;
    }
    g_debug = std::getenv("SPC_DEBUG") != nullptr;
    g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]);
    read_queries(argv[2]);
    tlog("read");
    run_queries_all();
    write_answers(argv[3]);
    tlog("write");
    // Nothing left to do: skip the destructors of the large containers.
    std::fflush(nullptr);
    std::_Exit(0);
}
