#pragma GCC optimize("O3")
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")
#include "common.inc"

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
        Heap4<D> hp;
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
    void witness(int32_t src, int32_t avoid, int budget, int pending) {
        ++searches;
        if (++wc == 0) { for (St& s : st) s.stamp = 0; wc = 1; }
        D limit = 0;
        for (int32_t x : tl) if (st[x].tm == tc) limit = std::max(limit, tlim[x]);
        { int b = 0; D L = limit; while (L > 1) { L >>= 1; ++b; } ++allhist[b]; }
        heap.clear();
        st[src].d = 0; st[src].stamp = wc; heap.push(0, src);
        int settled = 0;
        while (!heap.empty()) {
            auto it = heap.pop();
            if (it.d > st[it.v].d) continue;
            if (it.d > limit) break;
            if (++settled > budget) { ++capped; capset += settled; { int b = 0; D L = limit; while (L > 1) { L >>= 1; ++b; } ++caphist[b]; } break; }
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
            if (pending) witness(nv[i].to, v, budget, pending);
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
    int PMODE = 0, QORD = 1; uint32_t SIM2 = 1000000;
    int64_t prio(int32_t v) {
        int deg = (int)adeg[v];
        if (PMODE == 1) return (int64_t)level[v] * LC + 250 * (int64_t)deg;
        int add = simulate(v);
        return (int64_t)level[v] * LC + (1000 * (int64_t)add) / std::max(1, deg);
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
        if (g_debug && (r % (n / 10) == 0 || n - r == 3000 || n - r == 1000 || n - r == 100000 || n - r == 50000 || n - r == 25000 || n - r == 12000 || n - r == 6000)) {
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
    vector<D> dF, dB; vector<int32_t> tF, tB; RadixHeap<D> hF, hB;
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
        }
        vector<Arc>().swap(uarc);
        dF.assign(n, INF); dB.assign(n, INF);
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
    // one step of an upward search; returns false when that side is exhausted
    template <class HP>
    __attribute__((noinline)) void step(vector<D>& d1, vector<int32_t>& tt, HP& h, const vector<D>& d2, D& mu, vector<int32_t>& ent) {
        auto it = h.pop();
        int32_t u = it.v; D d = it.d;
        if (d > d1[u]) return;
        ++qsettle;
        if (d2[u] < INF && d + d2[u] < mu) mu = d + d2[u];
        if (u >= ncontracted) { ent.push_back(u); return; }
        const uint32_t e0 = H[u], e1 = H[u + 1];
        for (uint32_t k = e0; k < e1; ++k) { const Arc& a = A[k]; if (d1[a.to] + a.w < d) return; }   // stall (INF sums stay large)
        qrel += e1 - e0;
        for (uint32_t k = e0; k < e1; ++k) {
            const Arc& a = A[k]; D nd = d + a.w;
            if (nd < d1[a.to]) { if (d1[a.to] == INF) tt.push_back(a.to); d1[a.to] = nd; h.push(nd, a.to); }
        }
    }
    __attribute__((noinline)) int64_t query(int32_t s, int32_t t) {
        if (s == t) return 0;
        s = qid[s]; t = qid[t];
        D mu = INF;
        eF.clear(); eB.clear();
        hF.clear(); hB.clear();
        dF[s] = 0; tF.push_back(s); hF.push(0, s);
        dB[t] = 0; tB.push_back(t); hB.push(0, t);
        for (;;) {
            bool f = !hF.empty() && hF.top_key() < mu, b = !hB.empty() && hB.top_key() < mu;
            if (!f && !b) break;
            if (f && (!b || hF.top_key() <= hB.top_key())) step(dF, tF, hF, dB, mu, eF);
            else step(dB, tB, hB, dF, mu, eB);
        }
        qes += eF.size() + eB.size();
        if (!eF.empty() && !eB.empty()) {
            const int32_t base = ncontracted;
            for (int32_t a : eF) {
                D da = dF[a]; if (da >= mu) continue;
                const D* row = &T[(size_t)(a - base) * C];
                for (int32_t b : eB) { D c = da + row[b - base] + dB[b]; if (c < mu) mu = c; }
                qlook += eB.size();
            }
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
    ch.LC = envi("LC", 1000); ch.STOPC = (int32_t)envi("STOPC", 3000);
    ch.CON_SET = (int)envi("CON_SET", 200); ch.PMODE = (int)envi("PMODE", 0); ch.QORD = (int)envi("QORD", 1); ch.SIM2 = (uint32_t)envi("SIM2", 1000000);
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
    for (int32_t e = 0; e < E; ++e) if (!dead[e] && eu[e] != ev[e]) ch.add_edge(perm[eu[e]], perm[ev[e]], (D)ew[e]);
    tlog("build");
    ch.landmarks((int)envi("LMK", 0));
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
    if (envi("LBLSTAT", 0)) {
        uint64_t tot = 0, totc = 0, nstall = 0; int ns = 2000;
        vector<int32_t> ent; typename CH<D>::Arc dummy; (void)dummy;
        for (int k = 0; k < ns; ++k) {
            int32_t s = ch.qid[perm[(int32_t)((uint64_t)k * 2654435761ull % V)]];
            D inf = CH<D>::INF; ent.clear();
            uint64_t before = ch.qsettle;
            RadixHeap<D> h; h.clear();
            ch.dF[s] = 0; ch.tF.push_back(s); h.push(0, s);
            while (!h.empty()) ch.step(ch.dF, ch.tF, h, ch.dB, inf, ent);
            tot += ch.qsettle - before; totc += ent.size(); (void)nstall;
            for (int32_t v : ch.tF) ch.dF[v] = CH<D>::INF; ch.tF.clear();
        }
        std::fprintf(stderr, "full upward search: settled %.1f core entries %.1f\n", (double)tot / ns, (double)totc / ns);
    }
    {
        // process queries in source order (Morton ids): consecutive queries share
        // most of their upward search spaces, which then stay in cache
        vector<std::pair<int64_t, int32_t>> qo(Q);
        for (int32_t i = 0; i < Q; ++i) qo[i] = {((int64_t)perm[qs[i]] << 32) | (uint32_t)perm[qt[i]], i};
        if (envi("QSORT", 1)) std::sort(qo.begin(), qo.end());
        for (int32_t j = 0; j < Q; ++j) { int32_t i = qo[j].second; qans[i] = ch.query(perm[qs[i]], perm[qt[i]]); }
    }
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
