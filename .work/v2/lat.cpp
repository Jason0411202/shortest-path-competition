#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")
#include "common.inc"

namespace lat {

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

// ------------------------------------------------------------------- CH
template <class D>
struct CH {
    struct Arc { int32_t to; D w; };
    static constexpr D INF = (D)(~(D)0) >> 2;
    int32_t n;
    vector<vector<Arc>> g;           // remaining graph (undirected, symmetric)
    vector<int32_t> rank_of;
    vector<vector<Arc>> up;          // upward arcs (by original id)
    vector<D> wd; vector<uint32_t> ws; uint32_t wc = 0;
    vector<uint32_t> tmark; uint32_t tc = 0;
    Heap4<D> heap;
    uint64_t settles = 0, searches = 0, capped = 0, capset = 0;
    int SIM_SET = 0, CON_SET = 200;
    int64_t LC = 1000; int32_t STOPC = 0;
    double tprio = 0, tcon = 0;
    D UB = INF;

    void init(int32_t n_) {
        n = n_; g.assign(n, {}); wd.assign(n, INF); ws.assign(n, 0); tmark.assign(n, 0);
        rank_of.assign(n, -1); up.assign(n, {});
    }
    void add_edge(int32_t u, int32_t v, D w) {
        for (Arc& a : g[u]) if (a.to == v) { if (w < a.w) { a.w = w; for (Arc& b : g[v]) if (b.to == u) b.w = w; } return; }
        g[u].push_back({v, w}); g[v].push_back({u, w});
    }
    inline D getd(int32_t v) const { return ws[v] == wc ? wd[v] : INF; }
    // targets carry their own cost bound in tcost; limit shrinks as they settle
    vector<D> tcost; vector<int32_t> tl;
    // Witness search from src (avoiding `avoid`) for the pending targets in tl.
    // Target x needs a path of length <= tcost[x]; its last arc costs at least
    // mw(x), so only keys <= tlim[x] = tcost[x] - mw(x) can lead to it.  A
    // target is resolved as soon as a relaxation reaches it within tcost[x].
    vector<D> tlim;
    void witness(int32_t src, int32_t avoid, int budget, int pending) {
        ++searches;
        if (++wc == 0) { std::fill(ws.begin(), ws.end(), 0); wc = 1; }
        D limit = 0;
        for (int32_t x : tl) if (tmark[x] == tc) limit = std::max(limit, tlim[x]);
        heap.clear();
        wd[src] = 0; ws[src] = wc; heap.push(0, src);
        int settled = 0;
        while (!heap.empty()) {
            auto it = heap.pop();
            if (it.d > wd[it.v]) continue;
            if (it.d > limit) break;
            if (++settled > budget) { ++capped; capset += settled; break; }
            ++settles;
            bool relim = false;
            for (const Arc& a : g[it.v]) {
                if (a.to == avoid) continue;
                D nd = it.d + a.w;
                if (tmark[a.to] == tc && nd <= tcost[a.to]) {
                    tmark[a.to] = 0;
                    if (--pending <= 0) return;
                    if (tlim[a.to] >= limit) relim = true;
                }
                if (nd > limit) continue;
                if (ws[a.to] != wc || nd < wd[a.to]) { ws[a.to] = wc; wd[a.to] = nd; heap.push(nd, a.to); }
            }
            if (relim) { limit = 0; for (int32_t x : tl) if (tmark[x] == tc) limit = std::max(limit, tlim[x]); }
        }
    }
    struct SC { int32_t a, b; D w; };
    vector<D> mwv;
    int count_shortcuts(int32_t v, int budget, vector<SC>* out) {
        const vector<Arc>& nb = g[v];
        size_t k = nb.size();
        int cnt = 0;
        mwv.resize(k);
        for (size_t j = 0; j < k; ++j) {
            D m = INF;
            for (const Arc& a : g[nb[j].to]) if (a.to != v && a.w < m) m = a.w;
            mwv[j] = m;
        }
        for (size_t i = 0; i + 1 < k; ++i) {
            if (++tc == 0) { std::fill(tmark.begin(), tmark.end(), 0); tc = 1; }
            tl.clear();
            int pending = 0;
            for (size_t j = i + 1; j < k; ++j) {
                int32_t x = nb[j].to; D c = nb[i].w + nb[j].w;
                tl.push_back(x);
                if (mwv[j] > c || mwv[i] > c) { tmark[x] = tc; tcost[x] = 0; tlim[x] = 0; continue; }  // no witness possible
                tmark[x] = tc; tcost[x] = c; tlim[x] = c - mwv[j]; ++pending;
            }
            if (pending) {
                // targets that cannot have a witness must not block the search
                for (size_t j = i + 1; j < k; ++j) if (tcost[nb[j].to] == 0) tmark[nb[j].to] = tc + 0;
                witness_mask(nb[i].to, v, budget, pending, i, k, nb);
            }
            for (size_t j = i + 1; j < k; ++j) {
                if (tmark[nb[j].to] == tc && nb[i].w + nb[j].w <= UB) { ++cnt; if (out) out->push_back({nb[i].to, nb[j].to, nb[i].w + nb[j].w}); }
            }
        }
        return cnt;
    }
    // helper: run witness only over targets with a possible witness
    vector<int32_t> tl2;
    void witness_mask(int32_t src, int32_t avoid, int budget, int pending, size_t i, size_t k, const vector<Arc>& nb) {
        (void)i; (void)k; (void)nb;
        tl2.clear();
        for (int32_t x : tl) if (tcost[x] != 0) tl2.push_back(x);
        tl.swap(tl2);
        witness(src, avoid, budget, pending);
        tl.swap(tl2);
    }
    vector<int32_t> level;
    vector<D> md; vector<uint32_t> mst; uint32_t mc = 0;
    int simulate(int32_t v) {
        const vector<Arc>& nb = g[v];
        size_t k = nb.size();
        int added = 0;
        for (size_t i = 0; i + 1 < k; ++i) {
            int32_t u = nb[i].to;
            if (++mc == 0) { std::fill(mst.begin(), mst.end(), 0); mc = 1; }
            for (const Arc& a : g[u]) if (a.to != v) { mst[a.to] = mc; md[a.to] = a.w; }
            for (size_t j = i + 1; j < k; ++j) {
                int32_t x = nb[j].to; D c = nb[i].w + nb[j].w;
                bool found = mst[x] == mc && md[x] <= c;
                if (!found) for (const Arc& b : g[x]) if (b.to != v && mst[b.to] == mc && md[b.to] + b.w <= c) { found = true; break; }
                if (!found) ++added;
            }
        }
        return added;
    }
    int64_t prio(int32_t v) {
        double t0 = now();
        int deg = (int)g[v].size();
        int add = SIM_SET > 0 ? count_shortcuts(v, SIM_SET, nullptr) : simulate(v);
        int64_t r = (int64_t)level[v] * LC + (1000 * (int64_t)add) / std::max(1, deg);
        tprio += now() - t0;
        return r;
    }
    int32_t ncontracted = 0;
    void contract() {
        level.assign(n, 0); md.assign(n, INF); mst.assign(n, 0); tcost.assign(n, 0); tlim.assign(n, 0);
        vector<int64_t> pr(n);
        vector<uint8_t> dirty(n, 0);
        typedef std::pair<int64_t, int32_t> P;
        vector<P> init(n);
        for (int32_t v = 0; v < n; ++v) { pr[v] = prio(v); init[v] = {pr[v], v}; }
        std::priority_queue<P, vector<P>, std::greater<P>> pq(std::greater<P>(), std::move(init));
        tlog("init prio");
        vector<SC> scs;
        int32_t r = 0;
        while (!pq.empty()) {
            P t = pq.top(); pq.pop();
            int32_t v = t.second;
            if (rank_of[v] >= 0 || t.first != pr[v]) continue;
            if (n - r <= STOPC) break;
            if (dirty[v]) {
                dirty[v] = 0;
                int64_t np = prio(v);
                pr[v] = np;
                if (!pq.empty() && np > pq.top().first) { pq.push({np, v}); continue; }
            }
            scs.clear();
            { double t0 = now(); count_shortcuts(v, CON_SET, &scs); tcon += now() - t0; }
            rank_of[v] = r++;
            up[v] = g[v];
            for (const Arc& a : g[v]) {
                vector<Arc>& l = g[a.to];
                for (size_t q = 0; q < l.size(); ++q) if (l[q].to == v) { l[q] = l.back(); l.pop_back(); break; }
                level[a.to] = std::max(level[a.to], level[v] + 1);
                dirty[a.to] = 1;
            }
            vector<Arc>().swap(g[v]);
            for (const SC& s : scs) add_edge(s.a, s.b, s.w);
            if (g_debug && (r % (n / 10) == 0 || n - r == 3000 || n - r == 1000 || n - r == 300)) {
                size_t arcs = 0, mx = 0; for (int32_t x = 0; x < n; ++x) if (rank_of[x] < 0) { arcs += g[x].size(); mx = std::max(mx, g[x].size()); }
                std::fprintf(stderr, "  %d/%d searches=%llu settles=%llu capped=%llu (%llu) avgdeg=%.2f maxdeg=%zu tprio=%.2f tcon=%.2f\n", r, n,
                    (unsigned long long)searches, (unsigned long long)settles, (unsigned long long)capped, (unsigned long long)capset,
                    (double)arcs / std::max(1, n - r), mx, tprio, tcon);
                tlog("progress");
            }
        }
        ncontracted = r;
    }

    // ---------------------------------------------------------------- core
    // Remaining (uncontracted) vertices form the core; their all-pairs
    // distances go into a table, and they take the top ranks.
    int32_t C = 0;
    vector<int32_t> core_ids;        // core index -> original id
    vector<D> T;                     // C x C distances
    void build_core() {
        for (int32_t v = 0; v < n; ++v) if (rank_of[v] < 0) { rank_of[v] = ncontracted + (int32_t)core_ids.size(); core_ids.push_back(v); }
        C = (int32_t)core_ids.size();
        if (!C) return;
        // compact core graph
        vector<uint32_t> ch(C + 1, 0); vector<Arc> ca;
        for (int32_t i = 0; i < C; ++i) {
            for (const Arc& a : g[core_ids[i]]) ca.push_back({rank_of[a.to] - ncontracted, a.w});
            ch[i + 1] = (uint32_t)ca.size();
        }
        T.assign((size_t)C * C, INF);
        Heap4<D> h;
        for (int32_t s = 0; s < C; ++s) {
            D* row = &T[(size_t)s * C];
            row[s] = 0; h.clear(); h.push(0, s);
            while (!h.empty()) {
                auto it = h.pop();
                if (it.d > row[it.v]) continue;
                for (uint32_t k = ch[it.v]; k < ch[it.v + 1]; ++k) {
                    D nd = it.d + ca[k].w;
                    if (nd < row[ca[k].to]) { row[ca[k].to] = nd; h.push(nd, ca[k].to); }
                }
            }
        }
        if (g_debug) std::fprintf(stderr, "core C=%d arcs=%zu table=%.1f MB\n", C, ca.size(), (double)T.size() * sizeof(D) / 1e6);
    }

    // query structure: CSR by rank over contracted vertices only
    vector<uint32_t> H; vector<Arc> A;
    vector<D> dF, dB; vector<int32_t> tF, tB; Heap4<D> hF, hB;
    void build_query() {
        H.assign(n + 1, 0);
        for (int32_t v = 0; v < n; ++v) H[rank_of[v] + 1] = (uint32_t)up[v].size();
        for (int32_t r = 0; r < n; ++r) H[r + 1] += H[r];
        A.resize(H[n]);
        for (int32_t v = 0; v < n; ++v) { uint32_t k = H[rank_of[v]]; for (const Arc& a : up[v]) A[k++] = {rank_of[a.to], a.w}; }
        vector<vector<Arc>>().swap(up);
        dF.assign(n, INF); dB.assign(n, INF);
    }
    uint64_t qsettle = 0, qrel = 0, qes = 0, qlook = 0;
    // one upward search; core vertices are recorded, not expanded
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
    int64_t query(int32_t s, int32_t t) {
        if (s == t) return 0;
        s = rank_of[s]; t = rank_of[t];
        D mu = INF;
        eF.clear(); eB.clear();
        // forward fully first (mu unknown), then backward with mu pruning
        D inf = INF;
        upsearch(s, dF, tF, hF, dB, inf, eF);
        upsearch(t, dB, tB, hB, dF, mu, eB);
        // core combination
        qes += eF.size() + eB.size();
        const int32_t base = ncontracted;
        for (int32_t a : eF) {
            D da = dF[a]; if (da >= mu) continue;
            const D* row = &T[(size_t)(a - base) * C];
            for (int32_t b : eB) {
                D c = da + row[b - base] + dB[b];
                if (c < mu) mu = c;
            }
            qlook += eB.size();
        }
        for (int32_t v : tF) dF[v] = INF;
        for (int32_t v : tB) dB[v] = INF;
        tF.clear(); tB.clear();
        return mu >= INF ? -1 : (int64_t)mu;
    }
};

template <class D>
static void run(uint64_t ub) {
    CH<D> ch;
    ch.UB = (D)ub;
    ch.SIM_SET = (int)envi("SIM_SET", 0); ch.LC = envi("LC", 1000); ch.STOPC = (int32_t)envi("STOPC", 3000);
    ch.CON_SET = (int)envi("CON_SET", 200);
    ch.init(V);
    vector<uint8_t> dead(E, 0);
    if (dim && envi("PRUNE3", 1)) {
        // 3-hop square detours: drop edges that are strictly dominated
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
    for (int32_t e = 0; e < E; ++e) if (!dead[e] && eu[e] != ev[e]) ch.add_edge(eu[e], ev[e], (D)ew[e]);
    tlog("build");
    ch.contract();
    tlog("contract");
    if (g_debug) std::fprintf(stderr, "searches=%llu settles=%llu capped=%llu (%llu) tprio=%.2f tcon=%.2f\n",
        (unsigned long long)ch.searches, (unsigned long long)ch.settles, (unsigned long long)ch.capped, (unsigned long long)ch.capset, ch.tprio, ch.tcon);
    ch.build_core();
    tlog("core table");
    ch.build_query();
    if (g_debug) std::fprintf(stderr, "upward arcs: %u\n", ch.H[V]);
    for (int32_t i = 0; i < Q; ++i) qans[i] = ch.query(qs[i], qt[i]);
    if (g_debug) std::fprintf(stderr, "per query: settles %.1f relaxed %.1f entries %.1f lookups %.1f\n", (double)ch.qsettle / Q, (double)ch.qrel / Q, (double)ch.qes / Q, (double)ch.qlook / Q);
    tlog("queries");
}

static bool applicable() { return !directed; }
// Upper bound on every finite distance: twice the eccentricity of vertex 0
// (one Dijkstra).  Returns ~0 if the graph is disconnected.
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
    if (ub < (1ull << 29)) run<uint32_t>(ub); else run<uint64_t>(ub == ~0ull ? (~0ull >> 3) : ub);
}
}  // namespace lat

int main(int argc, char** argv) {
    if (argc != 4) return 1;
    g_debug = std::getenv("SPC_DEBUG") != nullptr;
    g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]); read_queries(argv[2]); tlog("read");
    (void)lat::applicable;
    lat::solve();
    write_answers(argv[3]); tlog("write");
    std::_Exit(0);
}
