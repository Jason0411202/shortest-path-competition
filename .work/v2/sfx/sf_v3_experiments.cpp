#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")
#include "common.inc"

// Directed graphs without coordinates (the power-law family).
namespace sf {

typedef uint32_t D;
static const D INF = 0x7FFFFFFFu;   // all finite distances stay below this (checked)

// 32-byte vertex record: the three lightest arcs inline, the rest in `ovf`.
struct Rec { uint32_t off; uint32_t deg; int32_t to[3]; D w[3]; };
struct Arc { int32_t to; D w; };

struct Graph {
    vector<Rec> r;
    vector<Arc> ovf;
};
static Graph GF, GB;                 // forward / backward
static vector<uint64_t> sinkbits;    // out-degree 0
static vector<uint64_t> hubbits;     // vertex is a hub (skipped by the hub-free search)
static inline bool bit(const vector<uint64_t>& b, int32_t v) { return (b[(size_t)v >> 6] >> (v & 63)) & 1; }

static int K = 0;
static vector<int32_t> hubs;
static vector<D> DOUT, DIN;          // DOUT[v*K+k] = d(v, hub k), DIN[v*K+k] = d(hub k, v)

bool applicable() { return directed && !has_coords; }

static void build_graph(const vector<int32_t>& a, const vector<int32_t>& b, Graph& g) {
    vector<uint32_t> off((size_t)V + 1, 0);
    for (int32_t i = 0; i < E; ++i) if (a[i] != b[i]) ++off[a[i] + 1];
    for (int32_t v = 0; v < V; ++v) off[v + 1] += off[v];
    vector<Arc> tmp(off[V]);
    {
        vector<uint32_t> cur(off.begin(), off.end() - 1);
        for (int32_t i = 0; i < E; ++i) if (a[i] != b[i]) tmp[cur[a[i]]++] = {b[i], (D)ew[i]};
    }
    g.r.assign((size_t)V, Rec{});
    size_t novf = 0;
    for (int32_t v = 0; v < V; ++v) { uint32_t d = off[v + 1] - off[v]; if (d > 3) novf += d - 3; }
    g.ovf.resize(novf);
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

// ------------------------------------------------------------ per-query state
// One open-addressing table holds both search directions' distances; heap
// items refer to table slots, so a pop needs no hash probe.
struct Slot { int32_t v; D d[2]; uint32_t pad; };
struct Table {
    vector<Slot> s;
    vector<uint32_t> used;
    uint32_t mask = 0;
    void init(uint32_t cap) { s.assign(cap, Slot{-1, {INF, INF}, 0}); mask = cap - 1; used.clear(); used.reserve(cap / 2 + 16); }
    inline uint32_t find(int32_t v) {          // insert if absent
        uint32_t i = ((uint32_t)v * 0x9E3779B1u) >> 8 & mask;
        for (;;) {
            int32_t k = s[i].v;
            if (k == v) return i;
            if (k == -1) { s[i].v = v; used.push_back(i); return i; }
            i = (i + 1) & mask;
        }
    }
    void clear() { for (uint32_t i : used) s[i] = Slot{-1, {INF, INF}, 0}; used.clear(); }
    bool full() const { return used.size() * 2 > mask; }
};

// direct-indexed alternative (slot = vertex id)
struct DTable {
    vector<Slot> s;
    vector<uint32_t> used;
    uint32_t mask = 0;
    void init(uint32_t) { if (s.empty()) s.assign((size_t)V, Slot{-1, {INF, INF}, 0}); mask = 0xFFFFFFFF; used.clear(); }
    inline uint32_t find(int32_t v) { if (s[v].v == -1) { s[v].v = v; used.push_back((uint32_t)v); } return (uint32_t)v; }
    void clear() { for (uint32_t i : used) s[i] = Slot{-1, {INF, INF}, 0}; used.clear(); }
    bool full() const { return false; }
};
#ifdef USE_DIRECT
#define Table DTable
#endif
struct HItem { D k; uint32_t slot; };
// Monotone radix heap on 32-bit keys.
struct RHeap {
    vector<HItem> b[33];
    D last = 0;
    uint32_t n = 0;
    static inline int bk(D k, D last) { return k == last ? 0 : 32 - __builtin_clz(k ^ last); }
    void clear() { for (auto& x : b) x.clear(); last = 0; n = 0; }
    inline void push(D k, uint32_t slot) { b[bk(k, last)].push_back({k, slot}); ++n; }
    inline bool empty() const { return n == 0; }
    // make b[0] non-empty (requires n > 0)
    inline void prep() {
        if (!b[0].empty()) return;
        int i = 1;
        while (b[i].empty()) ++i;
        D m = b[i][0].k;
        for (const HItem& it : b[i]) m = it.k < m ? it.k : m;
        last = m;
        for (const HItem& it : b[i]) b[bk(it.k, m)].push_back(it);
        b[i].clear();
    }
    inline HItem top() const { return b[0].back(); }
    inline HItem pop() { HItem it = b[0].back(); b[0].pop_back(); --n; return it; }
    void remap(const vector<uint32_t>& r) { for (auto& bb : b) for (HItem& it : bb) it.slot = r[it.slot]; }
};


// Branch-light 4-ary heap on packed (key << 32 | slot) items; the array is
// padded with ~0 so the four children of any node can be compared blindly.
struct PHeap {
    vector<uint64_t> a;
    uint32_t n = 0;
    D last = 0;
    PHeap() { a.assign(1024, ~0ull); }
    void clear() { for (uint32_t i = 0; i < n; ++i) a[i] = ~0ull; n = 0; }
    inline bool empty() const { return n == 0; }
    inline void push(D k, uint32_t slot) {
        if (n + 8 >= a.size()) a.resize(a.size() * 2, ~0ull);
        uint64_t x = (uint64_t)k << 32 | slot;
        uint32_t i = n++;
        while (i) {
            uint32_t p = (i - 1) >> 2;
            if (a[p] <= x) break;
            a[i] = a[p];
            i = p;
        }
        a[i] = x;
    }
    inline void prep() { last = (D)(a[0] >> 32); }
    inline HItem top() const { return HItem{(D)(a[0] >> 32), (uint32_t)a[0]}; }
    inline HItem pop() {
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
        if (n) a[i] = x; else a[0] = ~0ull;
        return HItem{(D)(best >> 32), (uint32_t)best};
    }
    void remap(const vector<uint32_t>& r) { for (uint32_t i = 0; i < n; ++i) a[i] = (a[i] & 0xFFFFFFFF00000000ull) | r[(uint32_t)a[i]]; }
};
#ifdef USE_RADIX
typedef RHeap QHeap;
#else
typedef PHeap QHeap;
#endif
static Table TB;
static uint32_t TB_DEFAULT = 1u << 13;
static QHeap hF, hB;
static uint64_t st_settle = 0, st_relax = 0, st_search = 0, st_onlyF = 0;

static uint64_t st_rehash = 0;
static vector<uint32_t> st_hist;
static void rehash() {
    ++st_rehash;
    // grow the table; heap items are remapped through the old table
    vector<Slot> old = TB.s;
    uint32_t cap = (TB.mask + 1) * 2;
    vector<uint32_t> oldused = TB.used;
    TB.init(cap);
    vector<uint32_t> remap(old.size(), 0);
    for (uint32_t i : oldused) { uint32_t j = TB.find(old[i].v); TB.s[j].d[0] = old[i].d[0]; TB.s[j].d[1] = old[i].d[1]; remap[i] = j; }
    hF.remap(remap); hB.remap(remap);
}

static vector<int32_t> hid;          // hub index (valid for hubs only)
static int LBMODE = 0, LBL = 1024, NOPEEK = 0;
static uint64_t BAL_NUM = 1, BAL_DEN = 1;
static uint64_t st_lbprune = 0, st_sf = 0;
static D rs[1024]; static char rs_ok[1024];

// MODE 0: plain bidirectional Dijkstra.
// MODE 1: forward search treats hubs as leaves; a settled hub h yields the
//         candidate f(h) + DIN[t][h] (exact hub -> t distance).  The backward
//         search never enters hubs.  Every path either has no interior hub
//         (found where the two searches meet) or has a first hub h, whose
//         prefix is hub-free and found by the forward search.
// MODE 2: exact UB = min_h d(s,h) + d(h,t) from both tables; both searches
//         avoid hubs entirely.
template <int MODE>
static D bidir(int32_t s, int32_t t, D mu) {
    ++st_search;
    hF.clear(); hB.clear();
    { uint32_t a = TB.find(s); TB.s[a].d[0] = 0; hF.push(0, a); }
    { uint32_t a = TB.find(t); TB.s[a].d[1] = 0; hB.push(0, a); }
    uint64_t wF = 0, wB = 0;
    const Rec* RF = GF.r.data();
    const Rec* RB = GB.r.data();
    const D* rowt = MODE == 1 ? &DIN[(size_t)t * K] : nullptr;
    D minin = INF;   // min over hubs h != t of d(h, t)
    if (MODE == 1) {
        for (int k = 0; k < K; ++k) if (hubs[k] != t && rowt[k] < minin) minin = rowt[k];
        if (LBMODE) { const D* r0 = &DIN[(size_t)s * K]; for (int k = 0; k < K; ++k) { rs[k] = r0[k]; rs_ok[k] = r0[k] != INF; } }
    }
    for (;;) {
        bool fe = hF.empty(), be = hB.empty();
        if (fe) break;
        hF.prep();
        D kf = hF.last;
        bool needF = MODE == 1 && kf + minin < mu;      // first-hub candidates still possible
        D kb = INF;
        bool meet = false;
        if (!be) { hB.prep(); kb = hB.last; meet = kf + kb < mu; }
        if (!meet && !needF) break;
        if (!meet) ++st_onlyF;
        bool fwd;
        if (!meet) fwd = true;
        else {
            if (NOPEEK) fwd = wF <= wB;
            else {
                int32_t uf = TB.s[hF.top().slot].v, ub = TB.s[hB.top().slot].v;
                fwd = (wF + RF[uf].deg) * BAL_DEN <= (wB + RB[ub].deg) * BAL_NUM;
            }
        }
        const int side = fwd ? 0 : 1;
        QHeap& hp = fwd ? hF : hB;
        const Arc* OV = (fwd ? GF : GB).ovf.data();
        const Rec* RX = fwd ? RF : RB;
        D other = fwd ? kb : kf;
        HItem it = hp.pop();
        D d = it.k;
        if (d > TB.s[it.slot].d[side]) continue;
        int32_t u = TB.s[it.slot].v;
        ++st_settle; if (fwd) ++st_sf;
        if (MODE == 1 && fwd && bit(hubbits, u)) {       // leaf: jump through the table
            D c = d + rowt[hid[u]];
            if (c < mu) mu = c;
            continue;
        }
        if (MODE == 1 && LBMODE && mu != INF && (fwd ? (LBMODE & 2) : (LBMODE & 1))) {
            // hub-landmark lower bound on the remaining distance
            const D* ry = &DIN[(size_t)u * K];
            int32_t lb = 0;
            if (!fwd) { for (int k = 0; k < K; ++k) { int32_t v = rs_ok[k] ? (int32_t)(ry[k] - rs[k]) : 0; lb = v > lb ? v : lb; } }
            else { for (int k = 0; k < LBL; ++k) { int32_t v = ry[k] != INF ? (int32_t)(rowt[k] - ry[k]) : 0; lb = v > lb ? v : lb; } }
            if ((uint64_t)d + (uint32_t)lb >= mu) { ++st_lbprune; continue; }
        }
        const Rec& R = RX[u];
        uint32_t deg = R.deg;
        (fwd ? wF : wB) += deg + 1;
        for (uint32_t j = 0; j < deg; ++j) {
            int32_t x; D w;
            if (j < 3) { x = R.to[j]; w = R.w[j]; } else { x = OV[R.off + j - 3].to; w = OV[R.off + j - 3].w; }
            D nd = d + w;
            if (MODE == 1 && fwd) { if (nd >= mu) break; }
            else if (nd + other >= mu) break;
            if (MODE == 2 && bit(hubbits, x)) continue;
            if (MODE == 1 && !fwd && bit(hubbits, x)) continue;
            ++st_relax;
            if (fwd && bit(sinkbits, x)) {   // dead end: only useful if it is t
                if (x == t && nd < mu) mu = nd;
                continue;
            }
            uint32_t sx = TB.find(x);
            Slot& S = TB.s[sx];
            if (nd < S.d[side]) {
                S.d[side] = nd;
                D o = S.d[side ^ 1];
                if (o != INF && nd + o < mu) mu = nd + o;
                __builtin_prefetch(&RX[x]);
                hp.push(nd, sx);
                if (TB.full()) rehash();
            }
        }
    }
    TB.clear();
    if (TB.mask + 1 > TB_DEFAULT) TB.init(TB_DEFAULT);
    return mu;
}

// MODE 3: unidirectional A* from s.  Potential pi(x) = max_k d(h_k,t) - d(h_k,x)
// (hub landmarks, consistent).  Hubs are leaves: a settled hub h closes the
// path with the exact table distance d(h,t).  Any s-t path is either hub-free
// (found when t is reached) or has a first hub h with a hub-free prefix.
static uint64_t st_pi = 0;
static inline D potential(const D* rowt, int32_t x) {
    const D* rx = &DIN[(size_t)x * K];
    int32_t lb = 0;
    bool dead = false;
    for (int k = 0; k < K; ++k) {
        D a = rx[k], b = rowt[k];
        if (a == INF) continue;
        dead |= b == INF;                       // h_k reaches x but not t: x cannot reach t
        int32_t v = (int32_t)(b - a);
        lb = v > lb ? v : lb;
    }
    return dead ? INF : (D)lb;
}

static D astar(int32_t s, int32_t t) {
    ++st_search;
    hF.clear();
    const D* rowt = &DIN[(size_t)t * K];
    D mu = INF;
    const Rec* RF = GF.r.data();
    const Arc* OV = GF.ovf.data();
    {
        D p = potential(rowt, s);
        if (p == INF) return INF;
        uint32_t a = TB.find(s); TB.s[a].d[0] = 0; TB.s[a].d[1] = p; hF.push(p, a);
        ++st_pi;
    }
    while (!hF.empty()) {
        hF.prep();
        if (hF.last >= mu) break;
        HItem it = hF.pop();
        Slot& U = TB.s[it.slot];
        D f = it.k - U.d[1];
        if (f > U.d[0]) continue;
        int32_t u = U.v;
        ++st_settle;
        if (u == t) { if (f < mu) mu = f; continue; }
        if (bit(hubbits, u)) { D c = f + rowt[hid[u]]; if (c < mu) mu = c; continue; }
        const Rec& R = RF[u];
        for (uint32_t j = 0; j < R.deg; ++j) {
            int32_t x; D w;
            if (j < 3) { x = R.to[j]; w = R.w[j]; } else { x = OV[R.off + j - 3].to; w = OV[R.off + j - 3].w; }
            D nd = f + w;
            if (nd >= mu) break;
            ++st_relax;
            if (x == t) { if (nd < mu) mu = nd; continue; }
            if (bit(sinkbits, x)) continue;
            uint32_t sx = TB.find(x);
            Slot& S = TB.s[sx];
            if (S.d[0] == INF && S.d[1] == INF) { S.d[1] = potential(rowt, x); ++st_pi; if (S.d[1] == INF) { S.d[0] = 0; continue; } }
            if (nd < S.d[0]) {
                S.d[0] = nd;
                D key = nd + S.d[1];
                if (key < mu) {
                    __builtin_prefetch(&RF[x]);
                    hF.push(key, sx);
                    if (TB.full()) rehash();
                }
            }
        }
    }
    TB.clear();
    if (TB.mask + 1 > TB_DEFAULT) TB.init(TB_DEFAULT);
    return mu;
}

// ------------------------------------------------------------ hub tables
static void full_dijkstra(const Graph& g, int32_t src, vector<D>& dist, RHeap& heap) {
    std::fill(dist.begin(), dist.end(), INF);
    heap.clear();
    dist[src] = 0; heap.push(0, (uint32_t)src);
    const Rec* RR = g.r.data();
    const Arc* OV = g.ovf.data();
    while (!heap.empty()) {
        heap.prep();
        HItem it = heap.pop();
        int32_t u = (int32_t)it.slot; D d = it.k;
        if (d > dist[u]) continue;
        const Rec& R = RR[u];
        for (uint32_t j = 0; j < R.deg; ++j) {
            int32_t x; D w;
            if (j < 3) { x = R.to[j]; w = R.w[j]; } else { x = OV[R.off + j - 3].to; w = OV[R.off + j - 3].w; }
            D nd = d + w;
            if (nd < dist[x]) { dist[x] = nd; __builtin_prefetch(&RR[x]); heap.push(nd, (uint32_t)x); }
        }
    }
}

static int MODE_SEL = 0;
static int RENUM = 0;

static vector<int32_t> settle_order(const Graph& g, const vector<int32_t>& srcs) {
    vector<D> dist((size_t)V, INF);
    vector<char> done((size_t)V, 0);
    vector<int32_t> ord;
    ord.reserve(V);
    RHeap heap;
    for (int32_t src : srcs) { dist[src] = 0; heap.push(0, (uint32_t)src); }
    while (!heap.empty()) {
        heap.prep();
        HItem it = heap.pop();
        int32_t u = (int32_t)it.slot; D d = it.k;
        if (d > dist[u] || done[u]) continue;
        done[u] = 1; ord.push_back(u);
        const Rec& R = g.r[u];
        for (uint32_t j = 0; j < R.deg; ++j) {
            int32_t x = j < 3 ? R.to[j] : g.ovf[R.off + j - 3].to;
            D w = j < 3 ? R.w[j] : g.ovf[R.off + j - 3].w;
            if (d + w < dist[x]) { dist[x] = d + w; heap.push(d + w, (uint32_t)x); }
        }
    }
    for (int32_t v = 0; v < V; ++v) if (!done[v]) ord.push_back(v);
    return ord;
}

// K-lane Bellman-Ford sweeps (pull form) in vertex-id order: row[x] =
// min(row[x], row[u] + w(u,x)) over in-arcs.  With vertices numbered in
// settle order from the top hub most shortest-path arcs point forward in the
// order, so few sweeps are needed; a vertex is re-evaluated only if an
// in-neighbour changed since its last evaluation.
static int SWEEPALT = 0;
static void build_hub_tables_vec() {
    if (std::getenv("SF_ALT")) SWEEPALT = 1;
    DIN.assign((size_t)V * K, INF);
    for (int k = 0; k < K; ++k) DIN[(size_t)hubs[k] * K + k] = 0;
    vector<uint32_t> chg((size_t)V, 1), evl((size_t)V, 0);   // change / evaluation times
    uint32_t now = 2;
    const Rec* RB = GB.r.data();
    const Arc* OV = GB.ovf.data();
    vector<D> nr((size_t)K);
    for (int pass = 0; ; ++pass) {
        uint64_t nev = 0, nch = 0;
        for (int32_t xi = 0; xi < V; ++xi) {
            int32_t x = (SWEEPALT && (pass & 1)) ? V - 1 - xi : xi;
            const Rec& R = RB[x];
            uint32_t deg = R.deg;
            if (deg == 0) continue;
            uint32_t ex = evl[x];
            bool need = false;
            for (uint32_t j = 0; j < deg; ++j) {
                int32_t u = j < 3 ? R.to[j] : OV[R.off + j - 3].to;
                if (chg[u] > ex) { need = true; break; }
            }
            ++now;
            evl[x] = now;
            if (!need) continue;
            ++nev;
            D* row = &DIN[(size_t)x * K];
            for (int k = 0; k < K; ++k) nr[k] = row[k];
            for (uint32_t j = 0; j < deg; ++j) {
                int32_t u; D w;
                if (j < 3) { u = R.to[j]; w = R.w[j]; } else { u = OV[R.off + j - 3].to; w = OV[R.off + j - 3].w; }
                const D* ru = &DIN[(size_t)u * K];
                D* __restrict q = nr.data();
                for (int k = 0; k < K; ++k) { D c = ru[k] + w; q[k] = c < q[k] ? c : q[k]; }
            }
            bool diff = false;
            for (int k = 0; k < K; ++k) diff |= nr[k] != row[k];
            if (diff) { for (int k = 0; k < K; ++k) row[k] = nr[k]; chg[x] = now; ++nch; }
        }
        if (g_debug) std::fprintf(stderr, "  sweep %d: evaluated %llu changed %llu\n", pass, (unsigned long long)nev, (unsigned long long)nch);
        if (nch == 0) break;
    }
}

static void build_hub_tables() {
    DIN.assign((size_t)V * K, INF);
    if (MODE_SEL == 2) DOUT.assign((size_t)V * K, INF);
    vector<D> dist((size_t)V);
    RHeap heap;
    for (int k = 0; k < K; ++k) {
        full_dijkstra(GF, hubs[k], dist, heap);
        for (int32_t v = 0; v < V; ++v) DIN[(size_t)v * K + k] = dist[v];
        if (MODE_SEL == 2) {
            full_dijkstra(GB, hubs[k], dist, heap);
            for (int32_t v = 0; v < V; ++v) DOUT[(size_t)v * K + k] = dist[v];
        }
    }
}

static inline D hub_ub(int32_t s, int32_t t) {
    const D* a = &DOUT[(size_t)s * K];
    const D* b = &DIN[(size_t)t * K];
    D best = INF;
    for (int k = 0; k < K; ++k) { D c = a[k] + b[k]; best = c < best ? c : best; }
    return best;
}

// ------------------------------------------------------------ interleaved queries
// Several independent queries advance in round-robin, one settle each, so
// the memory latency of one search overlaps with the work of the others.
// Each step ends by prefetching what the next step of that query will read.
struct BItem { D k; int32_t v; uint32_t slot; };
struct BHeap {
    vector<BItem> b[33];
    D last = 0;
    uint32_t n = 0;
    static inline int bk(D k, D last) { return k == last ? 0 : 32 - __builtin_clz(k ^ last); }
    void clear() { for (auto& x : b) x.clear(); last = 0; n = 0; }
    inline void push(D k, int32_t v, uint32_t slot) { b[bk(k, last)].push_back({k, v, slot}); ++n; }
    inline bool empty() const { return n == 0; }
    inline void prep() {
        if (!b[0].empty()) return;
        int i = 1;
        while (b[i].empty()) ++i;
        D m = b[i][0].k;
        for (const BItem& it : b[i]) m = it.k < m ? it.k : m;
        last = m;
        for (const BItem& it : b[i]) b[bk(it.k, m)].push_back(it);
        b[i].clear();
    }
    inline const BItem& top() const { return b[0].back(); }
    inline BItem pop() { BItem it = b[0].back(); b[0].pop_back(); --n; return it; }
};

struct Ctx {
    int32_t qi = -1, s = 0, t = 0;
    D mu = INF, minin = INF;
    uint64_t wF = 0, wB = 0;
    const D* rowt = nullptr;
    BHeap hF, hB;
    Table T;
    void rehash_ctx() {
        vector<Slot> old = T.s;
        vector<uint32_t> oldused = T.used;
        T.init((T.mask + 1) * 2);
        vector<uint32_t> remap(old.size(), 0);
        for (uint32_t i : oldused) { uint32_t j = T.find(old[i].v); T.s[j].d[0] = old[i].d[0]; T.s[j].d[1] = old[i].d[1]; remap[i] = j; }
        for (BHeap* h : {&hF, &hB}) for (auto& bb : h->b) for (BItem& it : bb) it.slot = remap[it.slot];
    }
};

static inline void ctx_start(Ctx& c, int32_t qi, int32_t s, int32_t t) {
    c.qi = qi; c.s = s; c.t = t; c.mu = INF; c.wF = c.wB = 0;
    c.hF.clear(); c.hB.clear();
    c.rowt = &DIN[(size_t)t * K];
    D m = INF;
    for (int k = 0; k < K; ++k) if (hubs[k] != t && c.rowt[k] < m) m = c.rowt[k];
    c.minin = m;
    { uint32_t a = c.T.find(s); c.T.s[a].d[0] = 0; c.hF.push(0, s, a); }
    { uint32_t a = c.T.find(t); c.T.s[a].d[1] = 0; c.hB.push(0, t, a); }
    __builtin_prefetch(&GF.r[s]);
    __builtin_prefetch(&GB.r[t]);
}

static inline void ctx_finish(Ctx& c) {
    c.T.clear();
    if (c.T.mask + 1 > TB_DEFAULT) c.T.init(TB_DEFAULT);
}

static inline void pf_row(const D* p) {
    const char* q = (const char*)p;
    for (int i = 0; i < K * 4; i += 64) __builtin_prefetch(q + i);
}

// one settle; returns false when the query is finished
static inline bool ctx_step(Ctx& c) {
    const Rec* RF = GF.r.data();
    const Rec* RB = GB.r.data();
    for (;;) {
        if (c.hF.empty()) return false;
        c.hF.prep();
        D kf = c.hF.last;
        bool needF = kf + c.minin < c.mu;
        D kb = INF;
        bool meet = false;
        if (!c.hB.empty()) { c.hB.prep(); kb = c.hB.last; meet = kf + kb < c.mu; }
        if (!meet && !needF) return false;
        bool fwd = !meet || c.wF + RF[c.hF.top().v].deg <= c.wB + RB[c.hB.top().v].deg;
        const int side = fwd ? 0 : 1;
        BHeap& hp = fwd ? c.hF : c.hB;
        const Arc* OV = (fwd ? GF : GB).ovf.data();
        const Rec* RX = fwd ? RF : RB;
        D other = fwd ? kb : kf;
        BItem it = hp.pop();
        D d = it.k;
        if (d > c.T.s[it.slot].d[side]) continue;          // stale: try again
        int32_t u = it.v;
        ++st_settle;
        if (fwd) {
            if (bit(hubbits, u)) {                          // leaf: jump through the table
                D cc = d + c.rowt[hid[u]];
                if (cc < c.mu) c.mu = cc;
                break;
            }
            if (LBMODE && c.mu != INF) {
                const D* ry = &DIN[(size_t)u * K];
                int32_t lb = 0;
                for (int k = 0; k < LBL; ++k) { int32_t v = ry[k] != INF ? (int32_t)(c.rowt[k] - ry[k]) : 0; lb = v > lb ? v : lb; }
                if ((uint64_t)d + (uint32_t)lb >= c.mu) { ++st_lbprune; break; }
            }
        }
        const Rec& R = RX[u];
        uint32_t deg = R.deg;
        (fwd ? c.wF : c.wB) += deg + 1;
        for (uint32_t j = 0; j < deg; ++j) {
            int32_t x; D w;
            if (j < 3) { x = R.to[j]; w = R.w[j]; } else { x = OV[R.off + j - 3].to; w = OV[R.off + j - 3].w; }
            D nd = d + w;
            if (fwd) { if (nd >= c.mu) break; }
            else {
                if (nd + other >= c.mu) break;
                if (bit(hubbits, x)) continue;
            }
            ++st_relax;
            if (fwd && bit(sinkbits, x)) {
                if (x == c.t && nd < c.mu) c.mu = nd;
                continue;
            }
            uint32_t sx = c.T.find(x);
            Slot& S = c.T.s[sx];
            if (nd < S.d[side]) {
                S.d[side] = nd;
                D o = S.d[side ^ 1];
                if (o != INF && nd + o < c.mu) c.mu = nd + o;
                __builtin_prefetch(&RX[x]);
                hp.push(nd, x, sx);
                if (c.T.full()) c.rehash_ctx();
            }
        }
        break;
    }
    // prefetch for the next step of this query
    if (!c.hF.empty() && !c.hF.b[0].empty()) {
        int32_t v = c.hF.b[0].back().v;
        __builtin_prefetch(&RF[v]);
        if (LBMODE) pf_row(&DIN[(size_t)v * K]);
    }
    if (!c.hB.empty() && !c.hB.b[0].empty()) __builtin_prefetch(&RB[c.hB.b[0].back().v]);
    return true;
}

static int NCTX = 4;

// answers queries whose index is in `list` (all need a search)
static void run_batch(const vector<int32_t>& list, const vector<int32_t>& S, const vector<int32_t>& Tt) {
    vector<Ctx> cx((size_t)NCTX);
    for (Ctx& c : cx) c.T.init(TB_DEFAULT);
    size_t next = 0;
    int active = 0;
    for (Ctx& c : cx) {
        if (next < list.size()) { int32_t q = list[next++]; ctx_start(c, q, S[q], Tt[q]); ++st_search; ++active; }
        else c.qi = -1;
    }
    while (active > 0) {
        for (Ctx& c : cx) {
            if (c.qi < 0) continue;
            if (!ctx_step(c)) {
                qans[c.qi] = c.mu >= INF ? -1 : (int64_t)c.mu;
                ctx_finish(c);
                if (next < list.size()) { int32_t q = list[next++]; ctx_start(c, q, S[q], Tt[q]); ++st_search; }
                else { c.qi = -1; --active; }
            }
        }
    }
}

void solve() {
    if (const char* e = std::getenv("SF_K")) K = std::atoi(e);
    if (const char* e = std::getenv("SF_MODE")) MODE_SEL = std::atoi(e);
    if (const char* e = std::getenv("SF_RENUM")) RENUM = std::atoi(e);
    if (const char* e = std::getenv("SF_LB")) LBMODE = std::atoi(e);
    if (const char* e = std::getenv("SF_LBL")) LBL = std::atoi(e);
    if (const char* e = std::getenv("SF_NOPEEK")) NOPEEK = std::atoi(e);
    if (const char* e = std::getenv("SF_BALN")) BAL_NUM = std::atoi(e);
    if (const char* e = std::getenv("SF_BALD")) BAL_DEN = std::atoi(e);
    if (MODE_SEL == 0) K = 0;
    if (K == 0) MODE_SEL = 0;
    if (LBL > K) LBL = K;
    vector<int32_t> indeg((size_t)V, 0), outdeg((size_t)V, 0);
    vector<D> maxout((size_t)V, 0);
    D maxw = 0;
    for (int32_t i = 0; i < E; ++i) if (eu[i] != ev[i]) {
        ++outdeg[eu[i]]; ++indeg[ev[i]];
        if (ew[i] > maxout[eu[i]]) maxout[eu[i]] = ew[i];
        if (ew[i] > maxw) maxw = ew[i];
    }
    uint64_t bound = maxw;
    for (int32_t v = 0; v < V; ++v) bound += maxout[v];
    if (bound >= INF / 2) { std::fprintf(stderr, "distance bound too large\n"); std::exit(2); }  // TODO fallback
    vector<int32_t> nid((size_t)V);
    for (int32_t v = 0; v < V; ++v) nid[v] = v;
    if (RENUM) {
        vector<int32_t> order;
        vector<int32_t> byin((size_t)V);
        for (int32_t v = 0; v < V; ++v) byin[v] = v;
        std::stable_sort(byin.begin(), byin.end(), [&](int32_t a, int32_t b) { return indeg[a] > indeg[b]; });
        if (RENUM == 1) { build_graph(eu, ev, GF); order = settle_order(GF, vector<int32_t>{byin[0]}); }
        else if (RENUM == 3) {
            vector<int32_t> src;
            for (int32_t k = 0; k < V && (int)src.size() < std::max(K, 1); ++k) if (outdeg[byin[k]] > 0) src.push_back(byin[k]);
            build_graph(eu, ev, GF); order = settle_order(GF, src);
        }
        else order = byin;
        for (int32_t k = 0; k < V; ++k) nid[order[k]] = k;
        for (int32_t i = 0; i < E; ++i) { eu[i] = nid[eu[i]]; ev[i] = nid[ev[i]]; }
        vector<int32_t> a((size_t)V), b((size_t)V);
        for (int32_t v = 0; v < V; ++v) { a[nid[v]] = indeg[v]; b[nid[v]] = outdeg[v]; }
        indeg.swap(a); outdeg.swap(b);
    }
    build_graph(eu, ev, GF);
    build_graph(ev, eu, GB);
    sinkbits.assign(((size_t)V + 63) / 64, 0);
    hubbits.assign(((size_t)V + 63) / 64, 0);
    for (int32_t v = 0; v < V; ++v) if (outdeg[v] == 0) sinkbits[v >> 6] |= 1ull << (v & 63);
    tlog("graph built");
    if (K > 0) {
        vector<int32_t> byin;
        for (int32_t v = 0; v < V; ++v) if (outdeg[v] > 0 && indeg[v] > 0) byin.push_back(v);
        std::stable_sort(byin.begin(), byin.end(), [&](int32_t a, int32_t b) { return indeg[a] > indeg[b]; });
        if ((int)byin.size() < K) K = (int)byin.size();
        hubs.assign(byin.begin(), byin.begin() + K);
        hid.assign((size_t)V, -1);
        for (int k = 0; k < K; ++k) { hubbits[hubs[k] >> 6] |= 1ull << (hubs[k] & 63); hid[hubs[k]] = k; }
        if (std::getenv("SF_VEC")) build_hub_tables_vec(); else build_hub_tables();
        tlog("hub tables");
    }
    if (const char* e = std::getenv("SF_TB")) TB_DEFAULT = 1u << std::atoi(e);
    TB.init(TB_DEFAULT);
    vector<char> bighub((size_t)V, 0);
    {
        vector<int32_t> o((size_t)V);
        for (int32_t v = 0; v < V; ++v) o[v] = v;
        std::stable_sort(o.begin(), o.end(), [&](int32_t a, int32_t b) { return indeg[a] > indeg[b]; });
        for (int k = 0; k < V / 1000; ++k) bighub[o[k]] = 1;
    }
    double tcat[2] = {0, 0}; uint64_t scat[2] = {0, 0}, ncat[2] = {0, 0};
    if (MODE_SEL == 5) {
        if (const char* e = std::getenv("SF_NCTX")) NCTX = std::atoi(e);
        vector<int32_t> S((size_t)Q), Tt((size_t)Q), list;
        for (int32_t i = 0; i < Q; ++i) {
            int32_t s = nid[qs[i]], t = nid[qt[i]];
            S[i] = s; Tt[i] = t;
            if (s == t) qans[i] = 0;
            else if (outdeg[s] == 0 || indeg[t] == 0) qans[i] = -1;
            else list.push_back(i);
        }
        run_batch(list, S, Tt);
    } else
    for (int32_t i = 0; i < Q; ++i) {
        auto T0 = std::chrono::steady_clock::now(); uint64_t S0 = st_settle; int cat = bighub[nid[qt[i]]];
        int32_t s = nid[qs[i]], t = nid[qt[i]];
        if (s == t) { qans[i] = 0; continue; }
        if (outdeg[s] == 0 || indeg[t] == 0) { qans[i] = -1; continue; }
        D r;
        if (MODE_SEL == 2) {
            D ub = hub_ub(s, t);
            if (bit(hubbits, s) || bit(hubbits, t)) r = ub;
            else r = bidir<2>(s, t, ub);
        } else if (MODE_SEL == 3) {
            r = astar(s, t);
        } else if (MODE_SEL == 1) {
            uint64_t s0 = st_settle;
            r = bidir<1>(s, t, INF);
            st_hist.push_back((uint32_t)(st_settle - s0));
        } else {
            { static int rep = std::getenv("SF_REP") ? std::atoi(std::getenv("SF_REP")) : 1; for (int q = 0; q < rep; ++q) r = bidir<0>(s, t, INF); }
        }
        qans[i] = r >= INF ? -1 : (int64_t)r;
        tcat[cat] += std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count(); scat[cat] += st_settle - S0; ncat[cat]++;
    }
    if (g_debug && !st_hist.empty()) {
        std::sort(st_hist.begin(), st_hist.end());
        size_t n = st_hist.size(); uint64_t tot = 0, top = 0;
        for (size_t i = 0; i < n; ++i) { tot += st_hist[i]; if (i >= n - n / 100) top += st_hist[i]; }
        std::fprintf(stderr, "settles p50 %u p90 %u p99 %u p999 %u max %u; top1%% share %.3f rehash %llu\n", st_hist[n / 2], st_hist[n * 9 / 10], st_hist[n * 99 / 100], st_hist[n * 999 / 1000], st_hist[n - 1], (double)top / tot, (unsigned long long)st_rehash);
    }
    if (g_debug) for (int c = 0; c < 2; ++c) std::fprintf(stderr, "cat%d n=%llu time %.3f settles/q %.1f\n", c, (unsigned long long)ncat[c], tcat[c], (double)scat[c] / ncat[c]);
    if (g_debug) std::fprintf(stderr, "pi/search %.1f fwd settles %.1f lbprune/search %.1f  onlyF settles/search %.1f\n", (double)st_pi / st_search, (double)st_sf / st_search, (double)st_lbprune / st_search, (double)st_onlyF / st_search);
    if (g_debug) std::fprintf(stderr, "searches %llu settle/search %.1f relax/search %.1f\n",
                              (unsigned long long)st_search, (double)st_settle / st_search, (double)st_relax / st_search);
    tlog("queries");
}

}  // namespace sf

int main(int argc, char** argv) {
    if (argc != 4) { std::fprintf(stderr, "usage: %s graph queries out\n", argv[0]); return 1; }
    g_debug = std::getenv("SPC_DEBUG") != nullptr;
    g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]); tlog("graph");
    read_queries(argv[2]); tlog("queries read");
    sf::solve();
    write_answers(argv[3]); tlog("write");
    std::_Exit(0);
}
