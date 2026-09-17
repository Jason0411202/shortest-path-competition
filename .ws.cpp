// Forward Dijkstra from src in the remaining graph, skipping `avoid`.  The
// targets are tlist, each with a threshold tcost; a target is resolved once
// it is settled.  The search stops when every target is resolved, when the
// key exceeds the largest threshold still pending, or after max_settle
// settled vertices.
static vector<int32_t> tlist;
static vector<W> tcost;

static void witness_search(int32_t src, int32_t avoid, int max_settle) {
    ++g_searches;
    ++wcur;
    if (wcur == 0) { std::fill(wstamp.begin(), wstamp.end(), 0); wcur = 1; }
    W limit = 0;
    for (int32_t x : tlist) limit = std::max(limit, tcost[x]);
    int pending = (int)tlist.size();
    wheap.clear();
    wdist[src] = 0; wstamp[src] = wcur;
    wheap.push(0, src);
    int settled = 0;
    while (!wheap.empty()) {
        MinHeap::Item it = wheap.pop();
        if (it.d > wdist[it.v]) continue;
        if (it.d > limit) break;
        if (tstamp[it.v] == tcur) {
            tstamp[it.v] = 0;
            if (--pending <= 0) break;
            if (tcost[it.v] >= limit) {
                limit = 0;
                for (int32_t x : tlist) if (tstamp[x] == tcur) limit = std::max(limit, tcost[x]);
                if (it.d > limit) break;
            }
        }
        if (++settled > max_settle) break;
        ++g_settles;
        for (const DArc& a : g[it.v]) {
            if (a.fw >= INF || a.to == avoid) continue;
            W nd = it.d + a.fw;
            if (nd > limit) continue;
            if (wstamp[a.to] != wcur || nd < wdist[a.to]) {
                wstamp[a.to] = wcur;
                wdist[a.to] = nd;
                wheap.push(nd, a.to);
            }
        }
    }
}

static inline void new_targets() {
    ++tcur;
    if (tcur == 0) { std::fill(tstamp.begin(), tstamp.end(), 0); tcur = 1; }
    tlist.clear();
}

// Deletes every arc u -> x for which a strictly shorter u -> x path exists
// (found by a bounded search).  Such arcs lie on no shortest path, and
// removing them one source at a time keeps all distances intact.
static void prune_redundant_arcs(int max_settle) {
    size_t removed = 0;
    for (int32_t u = 0; u < V; ++u) {
        new_targets();
        for (const DArc& a : g[u])
            if (a.fw < INF) {
                // tcost must be <= the arc weight minus one to prove strictness
                if (a.fw == 0) continue;
                tstamp[a.to] = tcur; tcost[a.to] = a.fw - 1; tlist.push_back(a.to);
            }
        if (tlist.size() < 2) continue;
        vector<int32_t> tg = tlist;
        witness_search(u, -1, max_settle);
        for (int32_t x : tg) {
            if (wget(x) >= tcost[x] + 1) continue;
            // strictly shorter path exists: drop u -> x
            for (size_t k = 0; k < g[u].size(); ++k) {
                DArc& a = g[u][k];
                if (a.to != x) continue;
                a.fw = INF;
                if (!directed) a.bw = INF;
                if (a.fw >= INF && a.bw >= INF) { a = g[u].back(); g[u].pop_back(); }
                break;
            }
            for (size_t k = 0; k < g[x].size(); ++k) {
                DArc& a = g[x][k];
                if (a.to != u) continue;
                a.bw = INF;
                if (!directed) a.fw = INF;
                if (a.fw >= INF && a.bw >= INF) { a = g[x].back(); g[x].pop_back(); }
                break;
            }
            ++removed;
        }
    }
    if (g_debug) std::fprintf(stderr, "pruned arcs: %zu\n", removed);
}
