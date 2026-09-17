s = open('sp.cpp', 'rb').read().decode()
code = r'''
// ------------------------------------------------ bidirectional dijkstra ----

static vector<uint32_t> bo_h, bi_h;
static vector<int32_t> bo_to, bi_to;
static vector<W> bo_w, bi_w;

static void build_bidir(const vector<int32_t>& eu, const vector<int32_t>& ev, const vector<W>& ew) {
    size_t m = eu.size();
    bo_h.assign((size_t)V + 1, 0); bi_h.assign((size_t)V + 1, 0);
    for (size_t i = 0; i < m; ++i) {
        ++bo_h[eu[i] + 1]; ++bi_h[ev[i] + 1];
        if (!directed) { ++bo_h[ev[i] + 1]; ++bi_h[eu[i] + 1]; }
    }
    for (int32_t v = 0; v < V; ++v) { bo_h[v + 1] += bo_h[v]; bi_h[v + 1] += bi_h[v]; }
    bo_to.resize(bo_h[V]); bo_w.resize(bo_h[V]); bi_to.resize(bi_h[V]); bi_w.resize(bi_h[V]);
    vector<uint32_t> co(bo_h.begin(), bo_h.end() - 1), ci(bi_h.begin(), bi_h.end() - 1);
    for (size_t i = 0; i < m; ++i) {
        uint32_t k = co[eu[i]]++; bo_to[k] = ev[i]; bo_w[k] = ew[i];
        k = ci[ev[i]]++; bi_to[k] = eu[i]; bi_w[k] = ew[i];
        if (!directed) {
            k = co[ev[i]]++; bo_to[k] = eu[i]; bo_w[k] = ew[i];
            k = ci[eu[i]]++; bi_to[k] = ev[i]; bi_w[k] = ew[i];
        }
    }
}

static uint64_t g_bd_settled = 0;

static int64_t bidir_query(int32_t s, int32_t t) {
    if (s == t) return 0;
    W mu = INF;
    hF.clear(); hB.clear();
    dF[s] = 0; touchedF.push_back(s); hF.push(0, s);
    dB[t] = 0; touchedB.push_back(t); hB.push(0, t);
    while (!hF.empty() && !hB.empty()) {
        if (hF.top_key() + hB.top_key() >= mu) break;
        bool fwd = hF.a.size() <= hB.a.size();
        MinHeap& h = fwd ? hF : hB;
        vector<W>& d1 = fwd ? dF : dB;
        vector<W>& d2 = fwd ? dB : dF;
        vector<int32_t>& tch = fwd ? touchedF : touchedB;
        const uint32_t* H = fwd ? bo_h.data() : bi_h.data();
        const int32_t* TO = fwd ? bo_to.data() : bi_to.data();
        const W* WT = fwd ? bo_w.data() : bi_w.data();
        MinHeap::Item it = h.pop();
        int32_t u = it.v; W d = it.d;
        if (d > d1[u]) continue;
        ++g_bd_settled;
        for (uint32_t k = H[u]; k < H[u + 1]; ++k) {
            int32_t x = TO[k];
            W nd = d + WT[k];
            if (nd < d1[x]) {
                if (d1[x] == INF) tch.push_back(x);
                d1[x] = nd;
                h.push(nd, x);
                if (d2[x] < INF && nd + d2[x] < mu) mu = nd + d2[x];
            }
        }
    }
    for (int32_t v : touchedF) dF[v] = INF;
    for (int32_t v : touchedB) dB[v] = INF;
    touchedF.clear(); touchedB.clear();
    return mu >= INF ? -1 : (int64_t)mu;
}

'''
anchor = '// --------------------------------------------------------------- output ----'
s = s.replace(anchor, code + anchor)
s = s.replace('static int USE_PLL = 0;', 'static int USE_PLL = 0, USE_BIDIR = 0;')
s = s.replace('        if (USE_PLL) { raw_u.push_back(u);', '        if (USE_PLL || USE_BIDIR) { raw_u.push_back(u);')
s = s.replace('        int64_t d = USE_PLL ? pll_query(s, t) :', '        int64_t d = USE_BIDIR ? bidir_query(s, t) : USE_PLL ? pll_query(s, t) :')
s = s.replace('    if (!USE_PLL) { dF.assign', '    if (!USE_PLL || USE_BIDIR) { dF.assign')
old = '''    if (USE_PLL) {
        build_pll(raw_u, raw_v, raw_w);'''
assert old in s
s = s.replace(old, '''    if (USE_BIDIR) {
        build_bidir(raw_u, raw_v, raw_w);
        run_queries(argv[2], argv[3]);
        tlog("queries");
        if (g_debug) std::fprintf(stderr, "bidir settled %llu\\n", (unsigned long long)g_bd_settled);
        return 0;
    }
''' + old)
s = s.replace('    if (const char* e = std::getenv("USE_PLL")) USE_PLL = std::atoi(e);', '    if (const char* e = std::getenv("USE_PLL")) USE_PLL = std::atoi(e);\n    if (const char* e = std::getenv("USE_BIDIR")) USE_BIDIR = std::atoi(e);')
open('sp.cpp', 'wb').write(s.encode())
