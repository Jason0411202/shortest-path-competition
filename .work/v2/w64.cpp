#include "common.inc"
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

int main(int argc, char** argv) {
    if (argc != 4) { std::fprintf(stderr, "usage: %s graph queries out\n", argv[0]); return 1; }
    g_debug = std::getenv("SPC_DEBUG") != nullptr;
    g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]); tlog("graph");
    read_queries(argv[2]); tlog("queries read");
    if (g_debug) std::fprintf(stderr, "w64 applicable=%d\n", (int)w64::applicable());
    w64::solve();
    write_answers(argv[3]); tlog("write");
    std::_Exit(0);
}
