// Kernel statistics for undirected lattices: full useless-edge pruning, then
// iterated removal of degree-1 vertices and degree-2 chain contraction.
#include "common.inc"
static int64_t envi(const char* k, int64_t def) { const char* e = std::getenv(k); return e ? std::atoll(e) : def; }
struct Arc { int32_t to; uint64_t w; };
int main(int argc, char** argv) {
    g_debug = 1; g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]); tlog("read");
    vector<vector<Arc>> g(V);
    for (int32_t e = 0; e < E; ++e) { g[eu[e]].push_back({ev[e], ew[e]}); g[ev[e]].push_back({eu[e], ew[e]}); }
    int budget = (int)envi("BUD", 1000);
    // prune: edge (u,v,w) removed if a strictly shorter path exists (sequential, current graph)
    vector<uint64_t> dist(V, ~0ull); vector<int32_t> touched; Heap4<uint64_t> h;
    int64_t removed = 0, capped = 0; uint64_t settles = 0;
    for (int32_t u = 0; u < V; ++u) {
        for (size_t ai = 0; ai < g[u].size(); ++ai) {
            int32_t v = g[u][ai].to; uint64_t w = g[u][ai].w;
            if (v < u) continue;
            // dijkstra from u, avoid direct arc, limit w-1
            h.clear(); dist[u] = 0; touched.push_back(u); h.push(0, u);
            bool found = false; int st = 0;
            while (!h.empty()) {
                auto it = h.pop(); if (it.d > dist[it.v]) continue;
                if (it.v == v) { found = true; break; }
                if (++st > budget) { ++capped; break; }
                ++settles;
                for (const Arc& a : g[it.v]) {
                    if (it.v == u && a.to == v) continue;
                    if (it.v == v && a.to == u) continue;
                    uint64_t nd = it.d + a.w; if (nd >= w) continue;
                    if (nd < dist[a.to]) { if (dist[a.to] == ~0ull) touched.push_back(a.to); dist[a.to] = nd; h.push(nd, a.to); }
                }
            }
            for (int32_t x : touched) dist[x] = ~0ull; touched.clear();
            if (found) {
                g[u][ai] = g[u].back(); g[u].pop_back(); --ai;
                for (size_t q = 0; q < g[v].size(); ++q) if (g[v][q].to == u) { g[v][q] = g[v].back(); g[v].pop_back(); break; }
                ++removed;
            }
        }
    }
    std::fprintf(stderr, "removed %lld capped %lld settles %llu\n", (long long)removed, (long long)capped, (unsigned long long)settles);
    tlog("prune");
    // degree histogram
    auto hist = [&](const char* tag) {
        int64_t c[8] = {0}; int64_t alive = 0, arcs = 0;
        for (int32_t v = 0; v < V; ++v) { if (g[v].empty()) continue; ++alive; arcs += g[v].size(); c[std::min<size_t>(7, g[v].size())]++; }
        std::fprintf(stderr, "%s alive=%lld edges=%lld deg1=%lld deg2=%lld deg3=%lld deg4=%lld deg5=%lld deg6=%lld deg7+=%lld\n", tag, (long long)alive, (long long)arcs / 2, (long long)c[1], (long long)c[2], (long long)c[3], (long long)c[4], (long long)c[5], (long long)c[6], (long long)c[7]);
    };
    hist("pruned");
    // iterated deg-1 removal and deg-2 contraction (with multi-edge min, no witness)
    vector<int32_t> stk; for (int32_t v = 0; v < V; ++v) if (g[v].size() <= 2) stk.push_back(v);
    vector<uint8_t> gone(V, 0);
    auto rm_arc = [&](int32_t a, int32_t b) { for (size_t q = 0; q < g[a].size(); ++q) if (g[a][q].to == b) { g[a][q] = g[a].back(); g[a].pop_back(); return; } };
    while (!stk.empty()) {
        int32_t v = stk.back(); stk.pop_back();
        if (gone[v]) continue;
        if (g[v].size() == 1) {
            int32_t a = g[v][0].to; rm_arc(a, v); g[v].clear(); gone[v] = 1; if (g[a].size() <= 2) stk.push_back(a);
        } else if (g[v].size() == 2) {
            Arc A = g[v][0], B = g[v][1];
            if (A.to == B.to) continue;
            rm_arc(A.to, v); rm_arc(B.to, v); g[v].clear(); gone[v] = 1;
            uint64_t w = A.w + B.w; bool ex = false;
            for (Arc& x : g[A.to]) if (x.to == B.to) { ex = true; if (w < x.w) { x.w = w; for (Arc& y : g[B.to]) if (y.to == A.to) y.w = w; } }
            if (!ex) { g[A.to].push_back({B.to, w}); g[B.to].push_back({A.to, w}); }
            if (g[A.to].size() <= 2) stk.push_back(A.to);
            if (g[B.to].size() <= 2) stk.push_back(B.to);
        }
    }
    hist("kernel");
    tlog("kernel");
}
