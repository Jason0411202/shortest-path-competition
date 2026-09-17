s = open('ch.cpp', 'rb').read().decode()
def rep(old, new, cnt=1):
    global s
    assert s.count(old) == cnt, (s.count(old), old[:80])
    s = s.replace(old, new)

rep('static int USE_PLL = 0;', 'static int USE_PLL = 0, USE_BIDIR = 0;')
rep('        if (USE_PLL) { raw_u.push_back(u);', '        if (USE_PLL || USE_BIDIR) { raw_u.push_back(u);')

BIDIR = r'''
// ------------------------------------------------ bidirectional dijkstra ----

// For directed graphs without coordinates (the power-law family): out-degree
// is tiny, and a bidirectional search balanced by scanned arcs meets long
// before either side reaches the hubs' huge in-lists.  Adjacency is sorted by
// weight so a scan stops at the first arc that cannot improve the best
// meeting distance.
static vector<uint32_t> fh, bh;
static vector<int32_t> fto, bto;
static vector<W> fwt, bwt;
static vector<W> bdF, bdB;
static vector<int32_t> bdTF, bdTB;
static MinHeap bhF, bhB;

static void build_csr(const vector<int32_t>& eu, const vector<int32_t>& ev, const vector<W>& ew,
                      vector<uint32_t>& h, vector<int32_t>& to, vector<W>& w) {
    size_t m = eu.size();
    h.assign((size_t)V + 1, 0);
    for (size_t i = 0; i < m; ++i) ++h[eu[i] + 1];
    for (int32_t v = 0; v < V; ++v) h[v + 1] += h[v];
    vector<std::pair<W, int32_t>> tmp(m);
    vector<uint32_t> cur(h.begin(), h.end() - 1);
    for (size_t i = 0; i < m; ++i) tmp[cur[eu[i]]++] = {ew[i], ev[i]};
    for (int32_t v = 0; v < V; ++v) std::sort(tmp.begin() + h[v], tmp.begin() + h[v + 1]);
    to.resize(m); w.resize(m);
    for (size_t i = 0; i < m; ++i) { w[i] = tmp[i].first; to[i] = tmp[i].second; }
}

static void build_bidir() {
    if (!directed) {
        size_t m = raw_u.size();
        for (size_t i = 0; i < m; ++i) { raw_u.push_back(raw_v[i]); raw_v.push_back(raw_u[i]); raw_w.push_back(raw_w[i]); }
    }
    build_csr(raw_u, raw_v, raw_w, fh, fto, fwt);
    build_csr(raw_v, raw_u, raw_w, bh, bto, bwt);
    vector<int32_t>().swap(raw_u); vector<int32_t>().swap(raw_v); vector<W>().swap(raw_w);
    bdF.assign((size_t)V, INF);
    bdB.assign((size_t)V, INF);
}

static int64_t bidir_query(int32_t s, int32_t t) {
    if (s == t) return 0;
    if (fh[s] == fh[s + 1] || bh[t] == bh[t + 1]) return -1;
    W mu = INF;
    bhF.clear(); bhB.clear();
    bdF[s] = 0; bdTF.push_back(s); bhF.push(0, s);
    bdB[t] = 0; bdTB.push_back(t); bhB.push(0, t);
    uint64_t workF = 0, workB = 0;
    while (!bhF.empty() && !bhB.empty()) {
        W kf = bhF.top_key(), kb = bhB.top_key();
        if (kf + kb >= mu) break;
        int32_t tf = bhF.a[0].v, tb = bhB.a[0].v;
        bool fwd = workF + (fh[tf + 1] - fh[tf]) <= workB + (bh[tb + 1] - bh[tb]);
        MinHeap& hp = fwd ? bhF : bhB;
        vector<W>& d1 = fwd ? bdF : bdB;
        const vector<W>& d2 = fwd ? bdB : bdF;
        vector<int32_t>& tch = fwd ? bdTF : bdTB;
        const uint32_t* H = fwd ? fh.data() : bh.data();
        const int32_t* TO = fwd ? fto.data() : bto.data();
        const W* WT = fwd ? fwt.data() : bwt.data();
        W other = fwd ? kb : kf;
        MinHeap::Item it = hp.pop();
        int32_t u = it.v; W d = it.d;
        if (d > d1[u]) continue;
        uint32_t k = H[u], e = H[u + 1];
        for (; k < e; ++k) {
            W nd = d + WT[k];
            if (nd >= mu) break;
            int32_t x = TO[k];
            if (nd < d1[x]) {
                if (d1[x] == INF) tch.push_back(x);
                d1[x] = nd;
                if (d2[x] < INF && nd + d2[x] < mu) mu = nd + d2[x];
                if (nd + other < mu && x != (fwd ? t : s) &&
                    (fwd ? fh[x] != fh[x + 1] : bh[x] != bh[x + 1])) hp.push(nd, x);
            }
        }
        (fwd ? workF : workB) += (k - H[u]) + 1;
    }
    for (int32_t v : bdTF) bdF[v] = INF;
    for (int32_t v : bdTB) bdB[v] = INF;
    bdTF.clear(); bdTB.clear();
    return mu >= INF ? -1 : (int64_t)mu;
}

// --------------------------------------------------------------- output ----
'''
rep('\n// --------------------------------------------------------------- output ----\n', BIDIR)
rep('    if (!USE_PLL) { dF.assign', '    if (!USE_PLL && !USE_BIDIR) { dF.assign')
rep('int64_t d = USE_PLL ? pll_query(s, t) :', 'int64_t d = USE_BIDIR ? bidir_query(s, t) : USE_PLL ? pll_query(s, t) :')
open('comb.cpp', 'wb').write(s.encode())

s = open('comb.cpp', 'rb').read().decode()
rep('''    directed = (FLAGS & FLAG_DIRECTED) != 0;
    g.assign((size_t)V, {});''', '''    directed = (FLAGS & FLAG_DIRECTED) != 0;
    // Method choice from the graph's shape (env overrides for experiments):
    // directed without coordinates -> bidirectional Dijkstra; otherwise a
    // contraction hierarchy, with hub labels when they fit in memory.
    USE_BIDIR = directed && !(FLAGS & FLAG_COORDS);
    USE_HL = V <= 1000000;
    if (const char* e = std::getenv("USE_BIDIR")) USE_BIDIR = std::atoi(e);
    if (const char* e = std::getenv("USE_HL")) USE_HL = std::atoi(e);
    if (!USE_PLL && !USE_BIDIR) g.assign((size_t)V, {});''')
rep('''    if (USE_PLL) {
        build_pll''', '''    if (USE_BIDIR) {
        build_bidir();
        tlog("bidir");
        run_queries(argv[2], argv[3]);
        tlog("queries");
        return 0;
    }
    if (USE_PLL) {
        build_pll''')
rep('''    if (const char* e = std::getenv("USE_HL")) USE_HL = std::atoi(e);
    if (const char* e = std::getenv("HL_PRUNE"))''', '''    if (const char* e = std::getenv("HL_PRUNE"))''')
open('comb.cpp', 'wb').write(s.encode())
