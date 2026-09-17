// ---------------------------------------------------------------- query ----

// CSR over ranks: vertex r's arcs lead to higher ranks.
static vector<uint32_t> uhead;
static vector<UArc> uarcs;

static void build_query_graph() {
    uhead.assign((size_t)V + 1, 0);
    for (int32_t v = 0; v < V; ++v) uhead[rank_of[v] + 1] = (uint32_t)upw[v].size();
    for (int32_t r = 0; r < V; ++r) uhead[r + 1] += uhead[r];
    uarcs.resize(uhead[V]);
    for (int32_t v = 0; v < V; ++v) {
        uint32_t k = uhead[rank_of[v]];
        for (const UArc& a : upw[v]) uarcs[k++] = {rank_of[a.to], a.fw, a.bw};
    }
    vector<vector<UArc>>().swap(upw);
}

static vector<W> dF, dB;
static vector<int32_t> touchedF, touchedB;
static MinHeap hF, hB;

static int64_t ch_query(int32_t s, int32_t t) {
    if (s == t) return 0;
    W mu = INF;
    hF.clear(); hB.clear();
    dF[s] = 0; touchedF.push_back(s); hF.push(0, s);
    dB[t] = 0; touchedB.push_back(t); hB.push(0, t);
    const uint32_t* H = uhead.data();
    const UArc* A = uarcs.data();

    while (!hF.empty() || !hB.empty()) {
        if (!hF.empty() && hF.top_key() >= mu) hF.clear();
        if (!hB.empty() && hB.top_key() >= mu) hB.clear();
        bool fwd;
        if (hF.empty()) { if (hB.empty()) break; fwd = false; }
        else if (hB.empty()) fwd = true;
        else fwd = hF.top_key() <= hB.top_key();

        if (fwd) {
            MinHeap::Item it = hF.pop();
            int32_t u = it.v; W d = it.d;
            if (d > dF[u]) continue;
            if (dB[u] < INF && d + dB[u] < mu) mu = d + dB[u];
            bool stalled = false;
            for (uint32_t k = H[u]; k < H[u + 1]; ++k) {
                const UArc& a = A[k];
                if (a.bw < INF && dF[a.to] < INF && dF[a.to] + a.bw < d) { stalled = true; break; }
            }
            if (stalled) continue;
            for (uint32_t k = H[u]; k < H[u + 1]; ++k) {
                const UArc& a = A[k];
                if (a.fw >= INF) continue;
                W nd = d + a.fw;
                if (nd < dF[a.to]) {
                    if (dF[a.to] == INF) touchedF.push_back(a.to);
                    dF[a.to] = nd;
                    hF.push(nd, a.to);
                }
            }
        } else {
            MinHeap::Item it = hB.pop();
            int32_t u = it.v; W d = it.d;
            if (d > dB[u]) continue;
            if (dF[u] < INF && d + dF[u] < mu) mu = d + dF[u];
            bool stalled = false;
            for (uint32_t k = H[u]; k < H[u + 1]; ++k) {
                const UArc& a = A[k];
                if (a.fw < INF && dB[a.to] < INF && dB[a.to] + a.fw < d) { stalled = true; break; }
            }
            if (stalled) continue;
            for (uint32_t k = H[u]; k < H[u + 1]; ++k) {
                const UArc& a = A[k];
                if (a.bw >= INF) continue;
                W nd = d + a.bw;
                if (nd < dB[a.to]) {
                    if (dB[a.to] == INF) touchedB.push_back(a.to);
                    dB[a.to] = nd;
                    hB.push(nd, a.to);
                }
            }
        }
    }
    for (int32_t v : touchedF) dF[v] = INF;
    for (int32_t v : touchedB) dB[v] = INF;
    touchedF.clear(); touchedB.clear();
    return mu >= INF ? -1 : (int64_t)mu;
}

// --------------------------------------------------------------- output ----

static inline char* put_u64(char* o, uint64_t x) {
    char tmp[24];
    int n = 0;
    do { tmp[n++] = (char)('0' + (x % 10)); x /= 10; } while (x);
    while (n) *o++ = tmp[--n];
    return o;
}

static void run_queries(const char* qpath, const char* opath) {
    Scanner sc;
    sc.load(qpath);
    int32_t Q = (int32_t)sc.uint_();

    dF.assign((size_t)V, INF);
    dB.assign((size_t)V, INF);

    vector<char> out((size_t)Q * 21 + 64);
    char* o = out.data();

    for (int32_t i = 0; i < Q; ++i) {
        int32_t s = (int32_t)sc.uint_();
        int32_t t = (int32_t)sc.uint_();
        int64_t d = ch_query(rank_of[s], rank_of[t]);
        if (d < 0) { *o++ = '-'; *o++ = '1'; }
        else o = put_u64(o, (uint64_t)d);
        *o++ = '\n';
    }

    std::FILE* f = std::fopen(opath, "wb");
    if (!f) { std::fprintf(stderr, "cannot open output file: %s\n", opath); std::exit(1); }
    std::fwrite(out.data(), 1, (size_t)(o - out.data()), f);
    std::fclose(f);
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: %s <graph_file> <query_file> <output_file>\n", argv[0]);
        return 1;
    }
    g_debug = std::getenv("SPC_DEBUG") != nullptr;
    g_t0 = std::chrono::steady_clock::now();
    if (const char* e = std::getenv("SIM_SETTLE")) SIM_SETTLE = std::atoi(e);
    if (const char* e = std::getenv("CON_SETTLE")) CON_SETTLE = std::atoi(e);
    if (const char* e = std::getenv("LEVEL_COEF")) LEVEL_COEF = std::atoi(e);
    read_graph(argv[1]);
    tlog("read");
    contract_all();
    tlog("contract");
    build_query_graph();
    if (g_debug) std::fprintf(stderr, "upward arcs: %u\n", uhead[V]);
    run_queries(argv[2], argv[3]);
    tlog("queries");
    return 0;
}
