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
};

static Table TB;
static RHeap hF, hB;
static uint64_t st_settle = 0, st_relax = 0, st_search = 0;

static void rehash() {
    // grow the table; heap items are remapped through the old table
    vector<Slot> old = TB.s;
    uint32_t cap = (TB.mask + 1) * 2;
    vector<uint32_t> oldused = TB.used;
    TB.init(cap);
    vector<uint32_t> remap(old.size(), 0);
    for (uint32_t i : oldused) { uint32_t j = TB.find(old[i].v); TB.s[j].d[0] = old[i].d[0]; TB.s[j].d[1] = old[i].d[1]; remap[i] = j; }
    for (RHeap* h : {&hF, &hB}) for (auto& bb : h->b) for (HItem& it : bb) it.slot = remap[it.slot];
}

template <bool HUBFREE>
static D bidir(int32_t s, int32_t t, D mu) {
    ++st_search;
    hF.clear(); hB.clear();
    { uint32_t a = TB.find(s); TB.s[a].d[0] = 0; hF.push(0, a); }
    { uint32_t a = TB.find(t); TB.s[a].d[1] = 0; hB.push(0, t == s ? a : a); }
    uint64_t wF = 0, wB = 0;
    const Rec* RF = GF.r.data();
    const Rec* RB = GB.r.data();
    while (!hF.empty() && !hB.empty()) {
        hF.prep(); hB.prep();
        D kf = hF.last, kb = hB.last;
        if (kf + kb >= mu) break;
        HItem tf = hF.top(), tb = hB.top();
        int32_t uf = TB.s[tf.slot].v, ub = TB.s[tb.slot].v;
        const bool fwd = wF + RF[uf].deg <= wB + RB[ub].deg;
        const int side = fwd ? 0 : 1;
        RHeap& hp = fwd ? hF : hB;
        const Rec& R = fwd ? RF[uf] : RB[ub];
        const Arc* OV = (fwd ? GF : GB).ovf.data();
        const Rec* RX = fwd ? RF : RB;
        D other = fwd ? kb : kf;
        HItem it = hp.pop();
        D d = it.k;
        if (d > TB.s[it.slot].d[side]) continue;
        ++st_settle;
        uint32_t deg = R.deg;
        (fwd ? wF : wB) += deg + 1;
        for (uint32_t j = 0; j < deg; ++j) {
            int32_t x; D w;
            if (j < 3) { x = R.to[j]; w = R.w[j]; } else { x = OV[R.off + j - 3].to; w = OV[R.off + j - 3].w; }
            D nd = d + w;
            if (nd + other >= mu) break;
            if (HUBFREE && bit(hubbits, x)) continue;
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
    return mu;
}
// ------------------------------------------------------------ hub tables
static void full_dijkstra(const Graph& g, int32_t src, vector<D>& dist, Heap4<D>& heap) {
    std::fill(dist.begin(), dist.end(), INF);
    heap.clear();
    dist[src] = 0; heap.push(0, src);
    const Rec* RR = g.r.data();
    const Arc* OV = g.ovf.data();
    while (!heap.empty()) {
        auto it = heap.pop();
        int32_t u = it.v; D d = it.d;
        if (d > dist[u]) continue;
        const Rec& R = RR[u];
        for (uint32_t j = 0; j < R.deg; ++j) {
            int32_t x; D w;
            if (j < 3) { x = R.to[j]; w = R.w[j]; } else { x = OV[R.off + j - 3].to; w = OV[R.off + j - 3].w; }
            D nd = d + w;
            if (nd < dist[x]) { dist[x] = nd; __builtin_prefetch(&RR[x]); heap.push(nd, x); }
        }
    }
}

static void build_hub_tables() {
    DOUT.assign((size_t)V * K, INF);
    DIN.assign((size_t)V * K, INF);
    vector<D> dist((size_t)V);
    Heap4<D> heap;
    for (int k = 0; k < K; ++k) {
        full_dijkstra(GF, hubs[k], dist, heap);
        for (int32_t v = 0; v < V; ++v) DIN[(size_t)v * K + k] = dist[v];
        full_dijkstra(GB, hubs[k], dist, heap);
        for (int32_t v = 0; v < V; ++v) DOUT[(size_t)v * K + k] = dist[v];
    }
}

static inline D hub_ub(int32_t s, int32_t t) {
    const D* a = &DOUT[(size_t)s * K];
    const D* b = &DIN[(size_t)t * K];
    D best = INF;
    for (int k = 0; k < K; ++k) { D c = a[k] + b[k]; best = c < best ? c : best; }
    return best;
}

static int RENUM = 0;

static vector<int32_t> settle_order(const Graph& g, int32_t src) {
    vector<D> dist((size_t)V, INF);
    vector<char> done((size_t)V, 0);
    vector<int32_t> ord;
    ord.reserve(V);
    Heap4<D> heap;
    dist[src] = 0; heap.push(0, src);
    while (!heap.empty()) {
        auto it = heap.pop();
        int32_t u = it.v; D d = it.d;
        if (d > dist[u] || done[u]) continue;
        done[u] = 1; ord.push_back(u);
        const Rec& R = g.r[u];
        for (uint32_t j = 0; j < R.deg; ++j) {
            int32_t x = j < 3 ? R.to[j] : g.ovf[R.off + j - 3].to;
            D w = j < 3 ? R.w[j] : g.ovf[R.off + j - 3].w;
            if (d + w < dist[x]) { dist[x] = d + w; heap.push(d + w, x); }
        }
    }
    for (int32_t v = 0; v < V; ++v) if (!done[v]) ord.push_back(v);
    return ord;
}

void solve() {
    if (const char* e = std::getenv("SF_K")) K = std::atoi(e);
    if (const char* e = std::getenv("SF_RENUM")) RENUM = std::atoi(e);
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
    // vertex order: new id = position
    vector<int32_t> byin((size_t)V);
    for (int32_t v = 0; v < V; ++v) byin[v] = v;
    std::stable_sort(byin.begin(), byin.end(), [&](int32_t a, int32_t b) { return indeg[a] > indeg[b]; });
    vector<int32_t> order;
    if (RENUM == 1) { build_graph(eu, ev, GF); order = settle_order(GF, byin[0]); }
    else if (RENUM == 2) order = byin;
    else { order.resize(V); for (int32_t v = 0; v < V; ++v) order[v] = v; }
    vector<int32_t> nid((size_t)V);
    for (int32_t k = 0; k < V; ++k) nid[order[k]] = k;
    {
        vector<int32_t> a((size_t)E), b((size_t)E);
        for (int32_t i = 0; i < E; ++i) { a[i] = nid[eu[i]]; b[i] = nid[ev[i]]; }
        build_graph(a, b, GF);
        build_graph(b, a, GB);
    }
    vector<int32_t> in2((size_t)V), out2((size_t)V);
    for (int32_t v = 0; v < V; ++v) { in2[nid[v]] = indeg[v]; out2[nid[v]] = outdeg[v]; }
    indeg.swap(in2); outdeg.swap(out2);
    sinkbits.assign(((size_t)V + 63) / 64, 0);
    hubbits.assign(((size_t)V + 63) / 64, 0);
    for (int32_t v = 0; v < V; ++v) if (outdeg[v] == 0) sinkbits[v >> 6] |= 1ull << (v & 63);
    tlog("graph built");
    if (K > 0) {
        hubs.clear();
        for (int32_t k = 0; k < V && (int)hubs.size() < K; ++k) {
            int32_t v = nid[byin[k]];
            if (outdeg[v] > 0 && indeg[v] > 0) hubs.push_back(v);
        }
        K = (int)hubs.size();
        for (int32_t h : hubs) hubbits[h >> 6] |= 1ull << (h & 63);
        build_hub_tables();
        tlog("hub tables");
    }
    TB.init(1 << 14);
    for (int32_t i = 0; i < Q; ++i) {
        int32_t s = nid[qs[i]], t = nid[qt[i]];
        if (s == t) { qans[i] = 0; continue; }
        if (outdeg[s] == 0 || indeg[t] == 0) { qans[i] = -1; continue; }
        D r;
        if (K > 0) {
            D ub = hub_ub(s, t);
            if (bit(hubbits, s) || bit(hubbits, t)) r = ub;
            else r = bidir<true>(s, t, ub);
        } else {
            { static int rep = std::getenv("SF_REP") ? std::atoi(std::getenv("SF_REP")) : 1; for (int q = 0; q < rep; ++q) r = bidir<false>(s, t, INF); }
        }
        qans[i] = r >= INF ? -1 : (int64_t)r;
    }
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
