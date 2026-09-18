#include "common.inc"
#pragma GCC push_options
#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")

// ============================================================ wide64
// Undirected graphs with 64-bit distances (the 2-D torus lattice family).
// Contraction hierarchy down to a core of C vertices, all-pairs table on the
// core, bidirectional upward queries combined through the table.
namespace w64 {

typedef uint64_t D;
static const D INF = (~0ull) >> 2;
static int64_t envi(const char* k, int64_t def) { const char* e = std::getenv(k); return e ? std::atoll(e) : def; }
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

bool applicable() { return !directed; }

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

// ------------------------------------------------------------ CH engine
struct Arc { int32_t to; uint32_t pad; D w; };
struct HI { D d; int32_t v; };

struct Engine {
    // dynamic graph over slots 0..m-1 (compacted periodically)
    int32_t m = 0, alive = 0;
    vector<Arc> pool; vector<uint32_t> aoff, adeg, acap;
    vector<int32_t> orig;          // slot -> vertex id
    vector<int64_t> pr; vector<int32_t> level; vector<uint8_t> dirty, dead;
    struct St { D d; uint32_t stamp, tm; };
    vector<St> st; vector<D> tcost, tlim, md; vector<uint32_t> ms;
    uint32_t wc = 0, tc = 0, mc = 0;
    Heap4<D> heap;
    // output hierarchy (vertex ids)
    vector<int32_t> rorder;        // contraction order
    vector<Arc> uarc; vector<uint32_t> uoff;
    vector<int32_t> vlevel;        // by vertex id
    // knobs
    int CON_SET = 200; int64_t LC = 1000; int32_t STOPC = 3000;
    int TWOHOP = 1;
    uint64_t st_2hop = 0, st_sim = 0;
    uint64_t searches = 0, settles = 0, st_pairs = 0, st_pairs_srch = 0, st_pairs_need = 0;
    double tprio = 0, tcon = 0;

    inline Arc* nb(int32_t v) { return &pool[aoff[v]]; }
    void alloc_slots(int32_t n) {
        m = n; aoff.assign(n, 0); adeg.assign(n, 0); acap.assign(n, 0);
        st.assign(n, {INF, 0, 0}); th.assign(n, {0, 0, 0, 0}); thc = 0; tcost.assign(n, 0); tlim.assign(n, 0); md.assign(n, INF); ms.assign(n, 0);
        wc = tc = mc = 0;
    }
    inline void push_arc(int32_t u, int32_t v, D w) {
        if (adeg[u] == acap[u]) {
            uint32_t nc = acap[u] * 2 + 4;
            uint32_t no = (uint32_t)pool.size();
            pool.resize(pool.size() + nc);
            if (adeg[u]) std::memcpy(&pool[no], &pool[aoff[u]], sizeof(Arc) * adeg[u]);
            aoff[u] = no; acap[u] = nc;
        }
        pool[aoff[u] + adeg[u]++] = {v, 0, w};
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

    // --- witness search from src avoiding `avoid`, targets in tl (tm == tc)
    vector<int32_t> tl;
    __attribute__((noinline)) void witness(int32_t src, int32_t avoid, int budget, int pending) {
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
            if (++settled > budget) break;
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
    // --- 1- and 2-hop witnesses for all neighbour pairs of v at once:
    // wit[i*k+j] (i<j) = 1 if u_i and u_j are joined by an arc or a 2-arc
    // path avoiding v of length <= w(u_i,v)+w(v,u_j).  Cost O(k*deg + #2-paths)
    // instead of O(k^2*deg).
    struct TH { uint32_t ns; int32_t ni; uint32_t hs; int32_t hd; };
    struct TE { int32_t i, next; D w; };
    vector<TH> th; vector<TE> te; vector<uint8_t> wit; uint32_t thc = 0; int THMODE = 0;
    void twohop(int32_t v) {
        const Arc* nv = nb(v); uint32_t k = adeg[v];
        if (++thc == 0) { for (TH& t : th) { t.ns = 0; t.hs = 0; } thc = 1; }
        wit.assign((size_t)k * k, 0);
        te.clear();
        for (uint32_t i = 0; i < k; ++i) { TH& t = th[nv[i].to]; t.ns = thc; t.ni = (int32_t)i; }
        for (uint32_t i = 0; i + 1 < k; ++i) {
            const Arc* au = nb(nv[i].to); uint32_t ku = adeg[nv[i].to];
            for (uint32_t q = 0; q < ku; ++q) {
                int32_t b = au[q].to; if (b == v) continue;
                TH& t = th[b];
                int32_t nx = t.hs == thc ? t.hd : -1;
                t.hs = thc; t.hd = (int32_t)te.size();
                te.push_back({(int32_t)i, nx, au[q].w});
            }
        }
        for (uint32_t j = 1; j < k; ++j) {
            const Arc* ax = nb(nv[j].to); uint32_t kx = adeg[nv[j].to];
            const D wj = nv[j].w;
            uint8_t* wcol = &wit[j];
            for (uint32_t q = 0; q < kx; ++q) {
                int32_t b = ax[q].to; if (b == v) continue;
                const D w2 = ax[q].w;
                const TH& t = th[b];
                if (t.ns == thc) { uint32_t i = (uint32_t)t.ni; if (i < j && w2 <= nv[i].w + wj) wcol[i * k] = 1; }
                if (t.hs == thc)
                    for (int32_t e = t.hd; e >= 0; e = te[e].next) {
                        uint32_t i = (uint32_t)te[e].i;
                        if (i < j && te[e].w + w2 <= nv[i].w + wj) wcol[i * k] = 1;
                    }
            }
        }
    }
    void shortcuts(int32_t v) {
        scs.clear();
        const Arc* nv = nb(v); uint32_t k = adeg[v];
        mwv.resize(k);
        for (uint32_t j = 0; j < k; ++j) {
            D mn = INF; const Arc* a = nb(nv[j].to); uint32_t kk = adeg[nv[j].to];
            for (uint32_t q = 0; q < kk; ++q) if (a[q].to != v && a[q].w < mn) mn = a[q].w;
            mwv[j] = mn;
        }
        if (THMODE) twohop(v);
        for (uint32_t i = 0; i + 1 < k; ++i) {
            if (++tc == 0) { for (St& s : st) s.tm = 0; tc = 1; }
            tl.clear();
            int pending = 0;
            int32_t u = nv[i].to;
            if (!THMODE) {
                if (++mc == 0) { std::fill(ms.begin(), ms.end(), 0); mc = 1; }
                const Arc* au = nb(u); uint32_t ku = adeg[u];
                for (uint32_t q = 0; q < ku; ++q) if (au[q].to != v) { ms[au[q].to] = mc; md[au[q].to] = au[q].w; }
            }
            for (uint32_t j = i + 1; j < k; ++j) {
                int32_t x = nv[j].to; D c = nv[i].w + nv[j].w;
                st[x].tm = tc;
                bool found;
                if (THMODE) found = wit[i * k + j];
                else {
                    found = ms[x] == mc && md[x] <= c;
                    if (!found) {
                        const Arc* ax = nb(x); uint32_t kx = adeg[x];
                        for (uint32_t q = 0; q < kx; ++q) { int32_t b = ax[q].to; if (b != v && ms[b] == mc && md[b] + ax[q].w <= c) { found = true; break; } }
                    }
                }
                if (found) { st[x].tm = 0; ++st_2hop; continue; }
                if (mwv[j] > c || mwv[i] > c) { tcost[x] = 0; continue; }   // no witness possible
                tcost[x] = c; tlim[x] = c - mwv[j]; tl.push_back(x); ++pending;
            }
            st_pairs_srch += pending;
            if (pending) witness(nv[i].to, v, CON_SET, pending);
            for (uint32_t j = i + 1; j < k; ++j)
                if (st[nv[j].to].tm == tc) { scs.push_back({nv[i].to, nv[j].to, nv[i].w + nv[j].w}); if (tcost[nv[j].to]) ++st_pairs_need; }
            st_pairs += k - 1 - i;
        }
    }
    __attribute__((noinline)) int simulate_old(int32_t v) {
        const Arc* nv = nb(v); uint32_t k = adeg[v];
        int added = 0;
        for (uint32_t i = 0; i + 1 < k; ++i) {
            int32_t u = nv[i].to;
            if (++mc == 0) { std::fill(ms.begin(), ms.end(), 0); mc = 1; }
            const Arc* au = nb(u); uint32_t ku = adeg[u];
            for (uint32_t q = 0; q < ku; ++q) if (au[q].to != v) { ms[au[q].to] = mc; md[au[q].to] = au[q].w; }
            for (uint32_t j = i + 1; j < k; ++j) {
                int32_t x = nv[j].to; D c = nv[i].w + nv[j].w;
                bool found = ms[x] == mc && md[x] <= c;
                if (!found) {
                    const Arc* ax = nb(x); uint32_t kx = adeg[x];
                    for (uint32_t q = 0; q < kx; ++q) { int32_t b = ax[q].to; if (b != v && ms[b] == mc && md[b] + ax[q].w <= c) { found = true; break; } }
                }
                if (!found) ++added;
            }
        }
        return added;
    }
    // 2-hop simulation: number of shortcuts that contracting v would add
    int simulate(int32_t v) {
        ++st_sim;
        uint32_t k = adeg[v];
        if (k < 2) return 0;
        if (!THMODE) return simulate_old(v);
        twohop(v);
        int added = 0;
        for (uint32_t i = 0; i + 1 < k; ++i) for (uint32_t j = i + 1; j < k; ++j) added += !wit[i * k + j];
        return added;
    }
    int PMODE = 0, SIM2 = 1, SIMEXACT = 0;
    int64_t prio(int32_t v) {
        int deg = (int)adeg[v];
        if (PMODE == 1) return (int64_t)level[v] * LC + 250 * (int64_t)deg;
        if (PMODE == 2) return (int64_t)level[v] * LC + 1000 * (int64_t)deg * (deg - 1) / 2 / 8 / std::max(1, deg);
        int add;
        if (SIMEXACT) { shortcuts(v); add = (int)scs.size(); } else add = simulate(v);
        return (int64_t)level[v] * LC + (1000 * (int64_t)add) / std::max(1, deg);
    }
    void contract_one(int32_t v) {
        double t0 = now();
        shortcuts(v);
        tcon += now() - t0;
        rorder.push_back(orig[v]);
        vlevel[orig[v]] = level[v];
        const Arc* nv = nb(v); uint32_t k = adeg[v];
        for (uint32_t i = 0; i < k; ++i) { uarc.push_back({orig[nv[i].to], 0, nv[i].w}); remove_arc(nv[i].to, v); }
        uoff.push_back((uint32_t)uarc.size());
        adeg[v] = 0; dead[v] = 1; --alive;
        for (const SC& s : scs) add_edge(s.a, s.b, s.w);
    }
    // renumber the alive vertices 0..alive-1 (slot order kept) and rebuild storage
    void compact() {
        vector<int32_t> ns(m, -1);
        int32_t n2 = 0;
        for (int32_t v = 0; v < m; ++v) if (!dead[v]) ns[v] = n2++;
        size_t tot = 0;
        for (int32_t v = 0; v < m; ++v) if (!dead[v]) tot += adeg[v] + adeg[v] / 2 + 4;
        vector<Arc> np(tot);
        vector<uint32_t> o2(n2), d2(n2), c2(n2);
        vector<int32_t> or2(n2), lv2(n2); vector<int64_t> pr2(n2); vector<uint8_t> di2(n2);
        size_t pos = 0;
        for (int32_t v = 0; v < m; ++v) if (!dead[v]) {
            int32_t w = ns[v];
            o2[w] = (uint32_t)pos; d2[w] = adeg[v]; c2[w] = adeg[v] + adeg[v] / 2 + 4;
            const Arc* a = nb(v);
            for (uint32_t i = 0; i < adeg[v]; ++i) np[pos + i] = {ns[a[i].to], 0, a[i].w};
            pos += c2[w];
            or2[w] = orig[v]; lv2[w] = level[v]; pr2[w] = pr[v]; di2[w] = dirty[v];
        }
        pool.swap(np);
        alloc_slots(n2);
        aoff.swap(o2); adeg.swap(d2); acap.swap(c2);
        orig.swap(or2); level.swap(lv2); pr.swap(pr2); dirty.swap(di2);
        dead.assign(n2, 0);
        alive = n2;
    }
    void debug_compact() {
        if (!g_debug) return;
        size_t arcs = 0; for (int32_t w = 0; w < m; ++w) arcs += adeg[w];
        std::fprintf(stderr, "  compact: alive=%d avgdeg=%.2f searches=%llu settles=%llu tprio=%.2f tcon=%.2f pairs=%llu searched=%llu needed=%llu twohop=%llu\n", alive,
                     (double)arcs / alive, (unsigned long long)searches, (unsigned long long)settles, tprio, tcon,
                     (unsigned long long)st_pairs, (unsigned long long)st_pairs_srch, (unsigned long long)st_pairs_need, (unsigned long long)st_2hop); std::fprintf(stderr, "  sims=%llu contracted=%d", (unsigned long long)st_sim, (int)rorder.size());
        tlog("  compact");
    }
    // Rounds of independent sets: every vertex whose priority is a strict
    // local minimum (ties by slot) is contracted, in slot (= spatial) order.
    void contract_rounds(int32_t until) {
        vector<int32_t> cand, touched;
        int32_t next_compact = alive / 2;
        int rounds = 0;
        while (alive > until && alive > STOPC) {
            ++rounds;
            cand.clear();
            for (int32_t v = 0; v < m; ++v) {
                if (dead[v]) continue;
                const Arc* a = nb(v); uint32_t k = adeg[v]; bool mn = true;
                int64_t p = pr[v];
                for (uint32_t i = 0; i < k; ++i) { int32_t u = a[i].to; if (pr[u] < p || (pr[u] == p && u < v)) { mn = false; break; } }
                if (mn) cand.push_back(v);
            }
            touched.clear();
            for (int32_t v : cand) {
                if (alive <= STOPC) break;
                const Arc* nv = nb(v); uint32_t k = adeg[v];
                for (uint32_t i = 0; i < k; ++i) { int32_t u = nv[i].to; level[u] = std::max(level[u], level[v] + 1); if (!dirty[u]) { dirty[u] = 1; touched.push_back(u); } }
                contract_one(v);
            }
            double t1 = now();
            std::sort(touched.begin(), touched.end());
            for (int32_t u : touched) { dirty[u] = 0; if (!dead[u]) pr[u] = prio(u); }
            tprio += now() - t1;
            if (alive <= next_compact) { compact(); next_compact = alive / 2; debug_compact(); }
        }
        if (g_debug) std::fprintf(stderr, "  rounds=%d\n", rounds);
    }
    void contract_greedy() {
        typedef std::pair<int64_t, int32_t> P;
        uoff.assign(1, 0);
        level.assign(m, 0); pr.assign(m, 0); dirty.assign(m, 0); dead.assign(m, 0);
        alive = m;
        double t0 = now();
        for (int32_t v = 0; v < m; ++v) pr[v] = prio(v);
        tprio += now() - t0;
        int cs = CON_SET; CON_SET = (int)envi("W_CONSET_BOT", CON_SET);
        contract_rounds((int32_t)envi("W_ROUNDS_UNTIL", 250000));
        CON_SET = cs;
        if (m != alive) { compact(); debug_compact(); }
        PMODE = (int)envi("W_PMODE_TOP", 0);
        SIM2 = (int)envi("W_SIM2_TOP", 1);
        SIMEXACT = (int)envi("W_SIMEXACT_TOP", 0);
        if (SIMEXACT) for (int32_t v = 0; v < m; ++v) pr[v] = prio(v);
        if (!SIM2) for (int32_t v = 0; v < m; ++v) pr[v] = prio(v);
        if (PMODE) for (int32_t v = 0; v < m; ++v) pr[v] = prio(v);
        vector<P> init(m);
        for (int32_t v = 0; v < m; ++v) init[v] = {pr[v], v};
        std::priority_queue<P, vector<P>, std::greater<P>> pq(std::greater<P>(), std::move(init));
        int32_t next_compact = alive / 2;
        const int32_t min_compact = (int32_t)envi("W_MINCOMPACT", 1000);
        while (!pq.empty() && alive > STOPC) {
            P t = pq.top(); pq.pop();
            int32_t v = t.second;
            if (dead[v] || t.first != pr[v]) continue;
            if (dirty[v]) {
                dirty[v] = 0;
                double t1 = now();
                int64_t np = prio(v);
                tprio += now() - t1;
                pr[v] = np;
                if (!pq.empty() && np > pq.top().first) { pq.push({np, v}); continue; }
            }
            const Arc* nv = nb(v); uint32_t k = adeg[v];
            for (uint32_t i = 0; i < k; ++i) { int32_t u = nv[i].to; level[u] = std::max(level[u], level[v] + 1); dirty[u] = 1; }
            contract_one(v);
            if (alive <= next_compact && alive >= min_compact) {
                compact();
                vector<P> re; re.reserve(alive);
                for (int32_t w = 0; w < m; ++w) re.push_back({pr[w], w});
                pq = std::priority_queue<P, vector<P>, std::greater<P>>(std::greater<P>(), std::move(re));
                next_compact = alive / 2;
                debug_compact();
            }
        }
    }
};

// ------------------------------------------------------------ core + query
struct Hier {
    int32_t n = 0, nc = 0, C = 0;          // vertices, contracted, core
    vector<int32_t> qid;                    // vertex id -> query index
    vector<uint32_t> H; vector<Arc> A;      // upward arcs by query index (contracted only)
    vector<D> T;                            // core table C x C
    vector<uint32_t> cH; vector<Arc> cA;    // core graph (core indices)
    bool useTable = true;
    uint64_t qcore = 0;
    // landmarks on the core (ALT for the core search)
    int K = 0;
    vector<int64_t> LD;                     // LD[c*K + l] = d(core c, landmark l), or LINF
    static constexpr int64_t LINF = INT64_MAX / 4;
    vector<int64_t> pot; vector<uint32_t> potst; uint32_t potc = 0;
    int64_t Us[64], Ut[64];
    void landmarks(int k) {
        K = k; if (!K || !C) { K = 0; return; }
        LD.assign((size_t)C * K, LINF);
        vector<int64_t> d(C), mind(C, LINF);
        Heap4<D> h;
        int32_t src = 0;
        for (int l = 0; l < K; ++l) {
            std::fill(d.begin(), d.end(), LINF);
            d[src] = 0; h.clear(); h.push(0, src);
            while (!h.empty()) {
                auto it = h.pop(); if ((int64_t)it.d > d[it.v]) continue;
                for (uint32_t q = cH[it.v]; q < cH[it.v + 1]; ++q) {
                    int64_t nd = (int64_t)it.d + (int64_t)cA[q].w;
                    if (nd < d[cA[q].to]) { d[cA[q].to] = nd; h.push((D)nd, cA[q].to); }
                }
            }
            int32_t far = src; int64_t fd = -1;
            for (int32_t v = 0; v < C; ++v) {
                LD[(size_t)v * K + l] = d[v];
                if (d[v] < mind[v]) mind[v] = d[v];
                if (mind[v] != LINF && mind[v] > fd) { fd = mind[v]; far = v; }
            }
            src = far;
        }
        pot.assign(C, 0); potst.assign(C, 0); potc = 0;
    }
    // potential (x2 scale): pi_t(v) - pi_s(v)
    inline int64_t potential(int32_t c) {
        if (potst[c] == potc) return pot[c];
        const int64_t* ld = &LD[(size_t)c * K];
        int64_t pt = 0, ps = 0;
        for (int l = 0; l < K; ++l) {
            int64_t a = ld[l];
            if (a == LINF) continue;
            int64_t x = Ut[l] == LINF ? 0 : (a > Ut[l] ? a - Ut[l] : Ut[l] - a);
            int64_t y = Us[l] == LINF ? 0 : (a > Us[l] ? a - Us[l] : Us[l] - a);
            pt = x > pt ? x : pt; ps = y > ps ? y : ps;
        }
        potst[c] = potc; pot[c] = pt - ps;
        return pot[c];
    }
    vector<D> dF, dB; vector<int32_t> tF, tB; Heap4<D> hF, hB;
    vector<int32_t> eF, eB;
    uint64_t qsettle = 0, qrel = 0, qes = 0, qlook = 0;

    void build(Engine& g, int32_t nv) {
        n = nv; nc = (int32_t)g.rorder.size();
        // core vertices in slot order (= vertex-id order after compaction)
        vector<int32_t> core;
        for (int32_t s = 0; s < g.m; ++s) if (!g.dead[s]) core.push_back(s);
        C = (int32_t)core.size();
        // query ids: contracted vertices ordered by (level, id), then core
        vector<int32_t> ordv(g.rorder);
        std::sort(ordv.begin(), ordv.end(), [&](int32_t a, int32_t b) { return g.vlevel[a] != g.vlevel[b] ? g.vlevel[a] < g.vlevel[b] : a < b; });
        qid.assign(n, -1);
        for (int32_t i = 0; i < nc; ++i) qid[ordv[i]] = i;
        for (int32_t c = 0; c < C; ++c) qid[g.orig[core[c]]] = nc + c;
        vector<int32_t> rank(n, -1);
        for (int32_t r = 0; r < nc; ++r) rank[g.rorder[r]] = r;
        H.assign(n + 1, 0);
        for (int32_t i = 0; i < nc; ++i) { int32_t r = rank[ordv[i]]; H[i + 1] = H[i] + (g.uoff[r + 1] - g.uoff[r]); }
        for (int32_t i = nc; i < n; ++i) H[i + 1] = H[i];
        A.resize(H[nc]);
        for (int32_t i = 0; i < nc; ++i) {
            int32_t r = rank[ordv[i]]; uint32_t k = H[i];
            for (uint32_t j = g.uoff[r]; j < g.uoff[r + 1]; ++j) A[k++] = {qid[g.uarc[j].to], 0, g.uarc[j].w};
        }
        vector<Arc>().swap(g.uarc);
        tlog("query graph");
        // core table
        vector<uint32_t> ch(C + 1, 0); vector<Arc> ca;
        for (int32_t i = 0; i < C; ++i) {
            const Arc* a = g.nb(core[i]);
            for (uint32_t q = 0; q < g.adeg[core[i]]; ++q) ca.push_back({qid[g.orig[a[q].to]] - nc, 0, a[q].w});
            ch[i + 1] = (uint32_t)ca.size();
        }
        cH.swap(ch); cA.swap(ca);
        useTable = C <= (int32_t)envi("W_TABLEMAX", 4000);
        dF.assign(n, INF); dB.assign(n, INF);
        if (!useTable) {
            landmarks((int)envi("W_LMK", 16));
            if (g_debug) std::fprintf(stderr, "core C=%d arcs=%zu upward arcs=%u landmarks=%d\n", C, cA.size(), H[n], K);
            tlog("core graph + landmarks");
            return;
        }
        const vector<uint32_t>& ch2 = cH; const vector<Arc>& ca2 = cA;
        T.assign((size_t)C * C, INF);
        Heap4<D> h;
        for (int32_t s = 0; s < C; ++s) {
            D* row = &T[(size_t)s * C];
            row[s] = 0; h.clear(); h.push(0, s);
            while (!h.empty()) {
                auto it = h.pop();
                if (it.d > row[it.v]) continue;
                for (uint32_t k = ch2[it.v]; k < ch2[it.v + 1]; ++k) {
                    D nd = it.d + ca2[k].w;
                    if (nd < row[ca2[k].to]) { row[ca2[k].to] = nd; h.push(nd, ca2[k].to); }
                }
            }
        }
        if (g_debug) std::fprintf(stderr, "core C=%d arcs=%zu upward arcs=%u\n", C, cA.size(), H[n]);
        tlog("core table");
        dF.assign(n, INF); dB.assign(n, INF);
    }
    inline void step(vector<D>& d1, vector<int32_t>& tt, Heap4<D>& h, const vector<D>& d2, D& mu, vector<int32_t>& ent) {
        auto it = h.pop();
        int32_t u = it.v; D d = it.d;
        if (d > d1[u]) return;
        ++qsettle;
        if (d2[u] < INF && d + d2[u] < mu) mu = d + d2[u];
        if (u >= nc) { ent.push_back(u); return; }
        const uint32_t e0 = H[u], e1 = H[u + 1];
        for (uint32_t k = e0; k < e1; ++k) { const Arc& a = A[k]; if (d1[a.to] + a.w < d) return; }   // stall
        qrel += e1 - e0;
        for (uint32_t k = e0; k < e1; ++k) {
            const Arc& a = A[k]; D nd = d + a.w;
            if (nd < d1[a.to]) { if (d1[a.to] == INF) tt.push_back(a.to); d1[a.to] = nd; h.push(nd, a.to); }
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
        const bool full = !useTable && K > 0;
        for (;;) {
            bool f = !hF.empty() && (full || hF.top_key() < mu), b = !hB.empty() && (full || hB.top_key() < mu);
            if (!f && !b) break;
            if (f && (!b || hF.top_key() <= hB.top_key())) step(dF, tF, hF, dB, mu, eF);
            else step(dB, tB, hB, dF, mu, eB);
        }
        qes += eF.size() + eB.size();
        if (!useTable && K > 0 && !eF.empty() && !eB.empty()) {
            // exact landmark distances of s and t through their core entries
            for (int l = 0; l < K; ++l) { Us[l] = LINF; Ut[l] = LINF; }
            for (int32_t a : eF) { const int64_t* ld = &LD[(size_t)(a - nc) * K]; int64_t da = (int64_t)dF[a];
                for (int l = 0; l < K; ++l) if (ld[l] != LINF && da + ld[l] < Us[l]) Us[l] = da + ld[l]; }
            for (int32_t b : eB) { const int64_t* ld = &LD[(size_t)(b - nc) * K]; int64_t db = (int64_t)dB[b];
                for (int l = 0; l < K; ++l) if (ld[l] != LINF && db + ld[l] < Ut[l]) Ut[l] = db + ld[l]; }
            if (++potc == 0) { std::fill(potst.begin(), potst.end(), 0); potc = 1; }
            hF.clear(); hB.clear();
            for (int32_t a : eF) hF.push((D)(2 * (int64_t)dF[a] + potential(a - nc)), a);
            for (int32_t b : eB) hB.push((D)(2 * (int64_t)dB[b] - potential(b - nc)), b);
            for (;;) {
                if (hF.empty() || hB.empty()) break;
                D kf = hF.top_key(), kb = hB.top_key();
                if (mu < INF && kf + kb >= 2 * mu) break;
                bool fw = kf <= kb;
                Heap4<D>& h = fw ? hF : hB;
                vector<D>& d1 = fw ? dF : dB; const vector<D>& d2 = fw ? dB : dF;
                vector<int32_t>& tt = fw ? tF : tB;
                auto it = h.pop();
                int32_t cu = it.v - nc;
                int64_t pu = potential(cu);
                D du = d1[it.v];
                if ((int64_t)it.d > 2 * (int64_t)du + (fw ? pu : -pu)) continue;   // stale
                ++qcore;
                for (uint32_t k = cH[cu]; k < cH[cu + 1]; ++k) {
                    int32_t x = cA[k].to + nc; D nd = du + cA[k].w;
                    if (nd < d1[x]) {
                        if (d1[x] == INF) tt.push_back(x);
                        d1[x] = nd;
                        if (d2[x] < INF && nd + d2[x] < mu) mu = nd + d2[x];
                        int64_t px = potential(cA[k].to);
                        h.push((D)(2 * (int64_t)nd + (fw ? px : -px)), x);
                    }
                }
            }
        }
        if (!useTable && K == 0 && !eF.empty() && !eB.empty()) {
            // bidirectional Dijkstra inside the core, seeded with the entries
            hF.clear(); hB.clear();
            for (int32_t a : eF) hF.push(dF[a], a);
            for (int32_t b : eB) hB.push(dB[b], b);
            for (;;) {
                if (hF.empty() || hB.empty()) break;
                D kf = hF.top_key(), kb = hB.top_key();
                if (kf + kb >= mu) break;
                bool fw = kf <= kb;
                Heap4<D>& h = fw ? hF : hB;
                vector<D>& d1 = fw ? dF : dB; const vector<D>& d2 = fw ? dB : dF;
                vector<int32_t>& tt = fw ? tF : tB;
                auto it = h.pop();
                if (it.d > d1[it.v]) continue;
                ++qcore;
                int32_t cu = it.v - nc;
                for (uint32_t k = cH[cu]; k < cH[cu + 1]; ++k) {
                    int32_t x = cA[k].to + nc; D nd = it.d + cA[k].w;
                    if (nd < d1[x]) {
                        if (d1[x] == INF) tt.push_back(x);
                        d1[x] = nd;
                        if (d2[x] < INF && nd + d2[x] < mu) mu = nd + d2[x];
                        h.push(nd, x);
                    }
                }
            }
        }
        if (useTable && !eF.empty() && !eB.empty()) {
            for (int32_t a : eF) {
                D da = dF[a]; if (da >= mu) continue;
                const D* row = &T[(size_t)(a - nc) * C];
                for (int32_t b : eB) { D c = da + row[b - nc] + dB[b]; if (c < mu) mu = c; }
                qlook += eB.size();
            }
        }
        for (int32_t v : tF) dF[v] = INF;
        for (int32_t v : tB) dB[v] = INF;
        tF.clear(); tB.clear();
        return mu >= INF ? -1 : (int64_t)mu;
    }
};

void solve() {
    detect();
    tlog("detect");
    // vertex numbering: Morton order of lattice coordinates, else identity
    vector<int32_t> perm(V);
    for (int32_t v = 0; v < V; ++v) perm[v] = v;
    if (dim) {
        vector<std::pair<uint64_t, int32_t>> key(V);
        for (int64_t v = 0; v < V; ++v) {
            uint64_t c[3] = {0, 0, 0};
            for (int k = 0; k < dim; ++k) c[k] = (uint64_t)((v / strd[k]) % side);
            uint64_t mo = 0;
            for (int b = 0; b < 21; ++b) for (int k = 0; k < dim; ++k) mo |= ((c[k] >> b) & 1ull) << (b * dim + k);
            key[v] = {mo, (int32_t)v};
        }
        std::sort(key.begin(), key.end());
        for (int32_t i = 0; i < V; ++i) perm[key[i].second] = i;
    }
    // edge pruning: an edge (u,v,w) is dropped if some path u-a-b-v of up to
    // three edges is strictly shorter (generic: over the adjacency)
    vector<uint8_t> deadE(E, 0);
    {
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
    Engine g;
    g.CON_SET = (int)envi("W_CONSET", 200); g.LC = envi("W_LC", 1000); g.STOPC = (int32_t)envi("W_STOPC", 3000); g.TWOHOP = (int)envi("W_TWOHOP", 1); g.THMODE = (int)envi("W_TH", 0);
    {
        vector<uint32_t> deg0(V, 0);
        for (int32_t e = 0; e < E; ++e) if (!deadE[e]) { ++deg0[perm[eu[e]]]; ++deg0[perm[ev[e]]]; }
        g.alloc_slots(V);
        size_t tot = 0;
        for (int32_t v = 0; v < V; ++v) { g.aoff[v] = (uint32_t)tot; g.acap[v] = deg0[v] + 2; tot += g.acap[v]; }
        g.pool.reserve(tot * 2); g.pool.resize(tot);
        g.orig.resize(V); for (int32_t v = 0; v < V; ++v) g.orig[v] = v;
        for (int32_t e = 0; e < E; ++e) if (!deadE[e]) g.add_edge(perm[eu[e]], perm[ev[e]], (D)ew[e]);
    }
    g.vlevel.assign(V, 0);
    tlog("build");
    g.contract_greedy();
    if (g_debug) std::fprintf(stderr, "contract: searches=%llu settles=%llu tprio=%.2f tcon=%.2f\n",
                              (unsigned long long)g.searches, (unsigned long long)g.settles, g.tprio, g.tcon);
    tlog("contract");
    if (envi("W_EXIT", 0)) std::_Exit(0);
    Hier hq;
    hq.build(g, V);
    {
        vector<std::pair<int64_t, int32_t>> qo(Q);
        for (int32_t i = 0; i < Q; ++i) qo[i] = {((int64_t)perm[qs[i]] << 32) | (uint32_t)perm[qt[i]], i};
        std::sort(qo.begin(), qo.end());
        for (int32_t j = 0; j < Q; ++j) { int32_t i = qo[j].second; qans[i] = hq.query(perm[qs[i]], perm[qt[i]]); }
    }
    if (g_debug) std::fprintf(stderr, "per query: settles %.1f relaxed %.1f entries %.1f lookups %.1f core %.1f\n", (double)hq.qsettle / Q,
                              (double)hq.qrel / Q, (double)hq.qes / Q, (double)hq.qlook / Q, (double)hq.qcore / Q);
    tlog("queries");
}

}  // namespace w64
#pragma GCC pop_options

int main(int argc, char** argv) {
    if (argc != 4) { std::fprintf(stderr, "usage: %s graph queries out\n", argv[0]); return 1; }
    g_debug = std::getenv("SPC_DEBUG") != nullptr;
    g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]); tlog("graph");
    read_queries(argv[2]); tlog("queries read");
    w64::solve();
    write_answers(argv[3]); tlog("write");
    std::_Exit(0);
}
