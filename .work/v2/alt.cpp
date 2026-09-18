// Prototype: bidirectional ALT (A*, landmarks, triangle inequality) for undirected graphs.
#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")
#include "common.inc"

static int64_t envi(const char* k, int64_t def) { const char* e = std::getenv(k); return e ? std::atoll(e) : def; }
typedef uint64_t D;
static const D INF = ~0ull >> 2;
static vector<uint32_t> H; static vector<int32_t> TO; static vector<uint32_t> WT;
static int K;
static vector<D> LM;   // LM[v*K+i]

static void dijkstra_all(int32_t s, vector<D>& d) {
    std::fill(d.begin(), d.end(), INF);
    Heap4<D> h; d[s] = 0; h.push(0, s);
    while (!h.empty()) {
        auto it = h.pop(); if (it.d > d[it.v]) continue;
        for (uint32_t k = H[it.v]; k < H[it.v + 1]; ++k) { D nd = it.d + WT[k]; if (nd < d[TO[k]]) { d[TO[k]] = nd; h.push(nd, TO[k]); } }
    }
}
static inline D lb(int32_t a, int32_t b) {
    const D* x = &LM[(size_t)a * K]; const D* y = &LM[(size_t)b * K];
    D m = 0;
    for (int i = 0; i < K; ++i) { D df = x[i] > y[i] ? x[i] - y[i] : y[i] - x[i]; if (df > m) m = df; }
    return m;
}
static vector<D> dF, dB; static vector<int32_t> tF, tB; static vector<D> potc; static vector<uint32_t> pstamp; static uint32_t pc = 0;
static uint64_t settl = 0;
// consistent average potential: p(v) = (lb(v,t) - lb(v,s)) / 2 ; forward reduced key = d + p(v)
static int32_t S_, T_;
static inline int64_t pot(int32_t v) {
    if (pstamp[v] == pc) return (int64_t)potc[v];
    int64_t p = ((int64_t)lb(v, T_) - (int64_t)lb(v, S_)) / 2;
    pstamp[v] = pc; potc[v] = (D)p;
    return p;
}
struct HeapS { struct It { int64_t k; int32_t v; }; vector<It> a;
    void clear() { a.clear(); } bool empty() const { return a.empty(); } int64_t top() const { return a[0].k; }
    void push(int64_t k, int32_t v) { size_t i = a.size(); a.push_back({k, v}); while (i) { size_t p = (i - 1) >> 2; if (a[p].k <= k) break; a[i] = a[p]; i = p; } a[i] = {k, v}; }
    It pop() { It b = a[0], l = a.back(); a.pop_back(); size_t n = a.size(); if (!n) return b; size_t i = 0; for (;;) { size_t c = 4 * i + 1; if (c >= n) break; size_t m = c, lim = std::min(n, c + 4); for (size_t j = c + 1; j < lim; ++j) if (a[j].k < a[m].k) m = j; if (a[m].k >= l.k) break; a[i] = a[m]; i = m; } a[i] = l; return b; }
};
static HeapS hf, hb;
static int64_t query(int32_t s, int32_t t) {
    if (s == t) return 0;
    S_ = s; T_ = t; if (++pc == 0) { std::fill(pstamp.begin(), pstamp.end(), 0); pc = 1; }
    D mu = INF;
    hf.clear(); hb.clear();
    dF[s] = 0; tF.push_back(s); hf.push(pot(s), s);
    dB[t] = 0; tB.push_back(t); hb.push(-pot(t), t);
    // forward key = d + p(v); backward key = d - p(v). stop when topF + topB >= mu
    while (!hf.empty() && !hb.empty()) {
        if ((D)(hf.top() + hb.top()) >= mu && hf.top() + hb.top() >= (int64_t)mu) break;
        bool fw = hf.top() <= hb.top();
        if (fw) {
            auto it = hf.pop(); int32_t u = it.v; D d = (D)(it.k - pot(u)); if (d > dF[u]) continue;
            ++settl;
            if (dB[u] != INF && d + dB[u] < mu) mu = d + dB[u];
            for (uint32_t k = H[u]; k < H[u + 1]; ++k) { int32_t x = TO[k]; D nd = d + WT[k]; if (nd < dF[x]) { if (dF[x] == INF) tF.push_back(x); dF[x] = nd; hf.push((int64_t)nd + pot(x), x); } }
        } else {
            auto it = hb.pop(); int32_t u = it.v; D d = (D)(it.k + pot(u)); if (d > dB[u]) continue;
            ++settl;
            if (dF[u] != INF && d + dF[u] < mu) mu = d + dF[u];
            for (uint32_t k = H[u]; k < H[u + 1]; ++k) { int32_t x = TO[k]; D nd = d + WT[k]; if (nd < dB[x]) { if (dB[x] == INF) tB.push_back(x); dB[x] = nd; hb.push((int64_t)nd - pot(x), x); } }
        }
    }
    for (int32_t v : tF) dF[v] = INF; for (int32_t v : tB) dB[v] = INF; tF.clear(); tB.clear();
    return mu >= INF ? -1 : (int64_t)mu;
}

int main(int argc, char** argv) {
    g_debug = std::getenv("SPC_DEBUG") != nullptr; g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]); read_queries(argv[2]); tlog("read");
    H.assign(V + 1, 0);
    for (int32_t e = 0; e < E; ++e) { ++H[eu[e] + 1]; ++H[ev[e] + 1]; }
    for (int32_t v = 0; v < V; ++v) H[v + 1] += H[v];
    TO.resize(H[V]); WT.resize(H[V]);
    { vector<uint32_t> c(H.begin(), H.end() - 1); for (int32_t e = 0; e < E; ++e) { TO[c[eu[e]]] = ev[e]; WT[c[eu[e]]++] = ew[e]; TO[c[ev[e]]] = eu[e]; WT[c[ev[e]]++] = ew[e]; } }
    tlog("csr");
    K = (int)envi("LMK", 8);
    LM.assign((size_t)V * K, 0);
    vector<D> d(V), mind(V, INF); int32_t src = 0;
    for (int q = 0; q < K; ++q) {
        dijkstra_all(src, d);
        int32_t far = 0; D fd = 0;
        for (int32_t v = 0; v < V; ++v) { LM[(size_t)v * K + q] = d[v]; mind[v] = std::min(mind[v], d[v]); if (mind[v] > fd && mind[v] != INF) { fd = mind[v]; far = v; } }
        src = far;
    }
    tlog("landmarks");
    dF.assign(V, INF); dB.assign(V, INF); potc.assign(V, 0); pstamp.assign(V, 0);
    int lim = (int)envi("QLIM", Q);
    for (int32_t i = 0; i < Q; ++i) {
        if (i < lim) { uint64_t s0 = settl; qans[i] = query(qs[i], qt[i]); if (g_debug && envi("PER", 0)) std::fprintf(stderr, "q %d settles %llu\n", i, (unsigned long long)(settl - s0)); }
    }
    if (g_debug) std::fprintf(stderr, "avg settles %.1f\n", (double)settl / std::min(lim, Q));
    tlog("queries");
    write_answers(argv[3]);
    std::_Exit(0);
}
