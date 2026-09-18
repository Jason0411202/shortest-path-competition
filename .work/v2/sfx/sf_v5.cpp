#include "common.inc"
// After the standard headers: a target pragma in front of them breaks
// GCC's always_inline allocator helpers when compiled without -march.
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
// Dijkstra order from the top hub (few sweeps, vectorised, streaming).
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
static vector<uint64_t> sinkbits;    // out-degree 0
static vector<uint64_t> hubbits;
static inline bool bit(const vector<uint64_t>& b, int32_t v) { return (b[(size_t)v >> 6] >> (v & 63)) & 1; }

static int K = 64;                   // hubs
static int LBL = 64;                 // landmark lanes used for forward pruning
static int PFROW = 1, LBGATE = 0, SWPD = 8;
static uint64_t BALF = 1, BALB = 1;
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
            int32_t uf = TB.s[hF.toppay()].v, ub = TB.s[hB.toppay()].v;
            fwd = (wF + RF[uf].deg) * BALF <= (wB + RB[ub].deg) * BALB;
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
                    if (HUBS && PFROW) {
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
    vector<uint32_t> chg((size_t)V, 1), evl((size_t)V, 0);
    uint32_t now = 2;
    const Rec* RB = GB.r.data();
    const Arc* OV = GB.ovf.data();
    vector<uint16_t> nr((size_t)K);
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
                        const char* q = (const char*)&T16[(size_t)u * K];
                        for (int b = 0; b < K * 2; b += 64) __builtin_prefetch(q + b);
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
            uint16_t* row = &T16[(size_t)x * K];
            uint16_t* __restrict q = nr.data();
            for (int k = 0; k < K; ++k) q[k] = row[k];
            for (uint32_t j = 0; j < deg; ++j) {
                int32_t u; D w;
                arc_at(R, OV, j, u, w);
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
            if (diff) { for (int k = 0; k < K; ++k) row[k] = q[k]; chg[x] = now; ++nch; }
        }
        if (g_debug) std::fprintf(stderr, "  sweep16 %d: changed %llu\n", pass, (unsigned long long)nch);
        if (nch == 0) return true;
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
static vector<int32_t> settle_order(const Graph& g, int32_t src) {
    vector<D> dist((size_t)V, INF);
    vector<char> done((size_t)V, 0);
    vector<int32_t> ord;
    ord.reserve(V);
    PHeap heap;
    dist[src] = 0; heap.push(0, (uint32_t)src);
    while (!heap.empty()) {
        uint64_t it = heap.pop();
        D d = (D)(it >> 32);
        int32_t u = (int32_t)(uint32_t)it;
        if (d > dist[u] || done[u]) continue;
        done[u] = 1; ord.push_back(u);
        const Rec& R = g.r[u];
        for (uint32_t j = 0; j < R.deg; ++j) {
            int32_t x; D w;
            arc_at(R, g.ovf.data(), j, x, w);
            if (d + w < dist[x]) { dist[x] = d + w; heap.push(d + w, (uint32_t)x); }
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
    if (const char* e = std::getenv("SF_LBL")) LBL = std::atoi(e);
    if (const char* e = std::getenv("SF_PFROW")) PFROW = std::atoi(e);
    if (const char* e = std::getenv("SF_SWPD")) SWPD = std::atoi(e);
    if (const char* e = std::getenv("SF_BALF")) BALF = std::atoi(e);
    if (const char* e = std::getenv("SF_BALB")) BALB = std::atoi(e);
    if (const char* e = std::getenv("SF_LBGATE")) LBGATE = std::atoi(e);
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
        build_graph(eu, ev, GF);
        vector<int32_t> order = settle_order(GF, byin[0]);
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
        bool conv = sweep16(40);
        EXACT16 = conv && verify16(NREACH) && !std::getenv("SF_NO16");
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

    TB.init(TB_DEFAULT);
    for (int32_t i = 0; i < Q; ++i) {
        int32_t s = nid[qs[i]], t = nid[qt[i]];
        if (s == t) { qans[i] = 0; continue; }
        if (outdeg[s] == 0 || indeg[t] == 0) { qans[i] = -1; continue; }
        auto T0 = std::chrono::steady_clock::now(); uint64_t s0 = st_settleF + st_settleB;
        D r = K > 0 ? query<true>(s, t) : query<false>(s, t);
        { int c = indeg[t] >= 130 ? (bit(hubbits, t) ? 2 : 1) : 0; static double tt[3]; static uint64_t ss[3], nn[3];
          tt[c] += std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count(); ss[c] += st_settleF + st_settleB - s0; nn[c]++;
          if (i == Q - 1 || (i == Q - 2)) for (int z = 0; z < 3; ++z) std::fprintf(stderr, "cat%d n=%llu time %.3f settles %.1f us/q %.2f\n", z, (unsigned long long)nn[z], tt[z], (double)ss[z] / (nn[z] + 1e-9), tt[z] / (nn[z] + 1e-9) * 1e6); }
        qans[i] = r >= INF ? -1 : (int64_t)r;
    }
    if (g_debug && st_search) std::fprintf(stderr, "K=%d searches %llu settled fwd %.1f bwd %.1f pruned %.1f\n", K,
                                           (unsigned long long)st_search, (double)st_settleF / st_search,
                                           (double)st_settleB / st_search, (double)st_prune / st_search);
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
