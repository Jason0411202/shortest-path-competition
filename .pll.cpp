// --------------------------------------------------- pruned landmark labels ----

// Pruned Landmark Labeling (Akiba, Iwata, Yoshida 2013), directed and
// weighted.  Vertices are processed by decreasing (in+1)*(out+1) degree; from
// each root a forward and a backward Dijkstra run, and a vertex whose
// distance is already answered by the labels built so far is neither labelled
// nor expanded.  Power-law graphs give tiny labels because nearly every
// shortest path touches one of the first few hubs.
struct PEntry { uint32_t hub; W d; };
static vector<vector<PEntry>> PIn, POut;   // PIn[u]: hub -> u, POut[u]: u -> hub
static vector<int32_t> pll_rank;           // original id -> processing index

static void build_pll(const vector<int32_t>& eu, const vector<int32_t>& ev, const vector<W>& ew) {
    vector<uint32_t> oh((size_t)V + 1, 0), ih((size_t)V + 1, 0);
    size_t m = eu.size();
    for (size_t i = 0; i < m; ++i) { ++oh[eu[i] + 1]; ++ih[ev[i] + 1]; }
    for (int32_t v = 0; v < V; ++v) { oh[v + 1] += oh[v]; ih[v + 1] += ih[v]; }
    vector<int32_t> oto(m), ito(m);
    vector<W> ow(m), iw(m);
    {
        vector<uint32_t> co(oh.begin(), oh.end() - 1), ci(ih.begin(), ih.end() - 1);
        for (size_t i = 0; i < m; ++i) {
            uint32_t k = co[eu[i]]++; oto[k] = ev[i]; ow[k] = ew[i];
            uint32_t j = ci[ev[i]]++; ito[j] = eu[i]; iw[j] = ew[i];
        }
    }
    vector<int32_t> order((size_t)V);
    for (int32_t v = 0; v < V; ++v) order[v] = v;
    std::sort(order.begin(), order.end(), [&](int32_t a, int32_t b) {
        uint64_t da = (uint64_t)(oh[a + 1] - oh[a] + 1) * (ih[a + 1] - ih[a] + 1);
        uint64_t db = (uint64_t)(oh[b + 1] - oh[b] + 1) * (ih[b + 1] - ih[b] + 1);
        return da != db ? da > db : a < b;
    });
    pll_rank.assign((size_t)V, 0);
    for (int32_t k = 0; k < V; ++k) pll_rank[order[k]] = k;

    PIn.assign((size_t)V, {});
    POut.assign((size_t)V, {});
    vector<W> T((size_t)V, INF);
    vector<W> dist((size_t)V, INF);
    vector<int32_t> touched;
    MinHeap heap;

    for (int32_t k = 0; k < V; ++k) {
        int32_t r = order[k];
        for (int pass = 0; pass < 2; ++pass) {
            const bool fwd = pass == 0;
            // T[hub] = dist(r -> hub) for the forward pass, dist(hub -> r) backward.
            const vector<PEntry>& own = fwd ? POut[r] : PIn[r];
            for (const PEntry& e : own) T[e.hub] = e.d;
            vector<vector<PEntry>>& target = fwd ? PIn : POut;
            const uint32_t* H = fwd ? oh.data() : ih.data();
            const int32_t* TO = fwd ? oto.data() : ito.data();
            const W* WT = fwd ? ow.data() : iw.data();
            heap.clear();
            dist[r] = 0; touched.push_back(r);
            heap.push(0, r);
            while (!heap.empty()) {
                MinHeap::Item it = heap.pop();
                int32_t u = it.v; W d = it.d;
                if (d > dist[u]) continue;
                bool pruned = false;
                for (const PEntry& e : target[u]) {
                    W t = T[e.hub];
                    if (t != INF && t + e.d <= d) { pruned = true; break; }
                }
                if (pruned) continue;
                target[u].push_back({(uint32_t)k, d});
                for (uint32_t a = H[u]; a < H[u + 1]; ++a) {
                    int32_t x = TO[a];
                    W nd = d + WT[a];
                    if (nd < dist[x]) {
                        if (dist[x] == INF) touched.push_back(x);
                        dist[x] = nd;
                        heap.push(nd, x);
                    }
                }
            }
            for (int32_t x : touched) dist[x] = INF;
            touched.clear();
            for (const PEntry& e : own) T[e.hub] = INF;
        }
    }
    if (g_debug) {
        size_t a = 0, b = 0;
        for (int32_t v = 0; v < V; ++v) { a += PIn[v].size(); b += POut[v].size(); }
        std::fprintf(stderr, "pll entries: in %zu out %zu (avg %.1f)\n", a, b, (double)(a + b) / (2.0 * V));
    }
}

static inline int64_t pll_query(int32_t s, int32_t t) {
    if (s == t) return 0;
    const vector<PEntry>& A = POut[s];
    const vector<PEntry>& B = PIn[t];
    size_t i = 0, j = 0;
    W best = INF;
    while (i < A.size() && j < B.size()) {
        if (A[i].hub < B[j].hub) ++i;
        else if (A[i].hub > B[j].hub) ++j;
        else { W c = A[i].d + B[j].d; if (c < best) best = c; ++i; ++j; }
    }
    return best >= INF ? -1 : (int64_t)best;
}

