// ---------------------------------------------------------- contraction ----

struct Shortcut { int32_t u, x; W w; int32_t hops; };
static vector<Shortcut> scs;

// Exact-as-bounded shortcut computation used when v is actually contracted:
// one witness Dijkstra per in-neighbour, targeting all out-neighbours.
static void find_shortcuts(int32_t v, int max_settle) {
    scs.clear();
    const vector<DArc>& nb = g[v];
    size_t n = nb.size();
    for (size_t i = 0; i < n; ++i) {
        if (nb[i].bw >= INF) continue;
        W wu = nb[i].bw;
        size_t j0 = directed ? 0 : i + 1;
        int cnt = 0;
        W mx = 0;
        new_targets();
        for (size_t j = j0; j < n; ++j)
            if (j != i && nb[j].fw < INF) { ++cnt; mx = std::max(mx, nb[j].fw); tstamp[nb[j].to] = tcur; }
        if (!cnt) continue;
        witness_search(nb[i].to, v, wu + mx, max_settle, cnt);
        for (size_t j = j0; j < n; ++j) {
            if (j == i || nb[j].fw >= INF) continue;
            W w = wu + nb[j].fw;
            if (wget(nb[j].to) > w)
                scs.push_back({nb[i].to, nb[j].to, w, nb[i].hops + nb[j].hops});
        }
    }
}

// Cheap estimate for the priority: a shortcut u -> x is considered redundant
// if a path of at most two arcs avoiding v is no longer.  No heap involved.
static vector<W> mdist;
static vector<uint32_t> mstamp;
static uint32_t mcur = 0;

static void simulate(int32_t v, int64_t& added, int64_t& added_hops) {
    added = 0; added_hops = 0;
    const vector<DArc>& nb = g[v];
    size_t n = nb.size();
    for (size_t i = 0; i < n; ++i) {
        if (nb[i].bw >= INF) continue;
        int32_t u = nb[i].to;
        W wu = nb[i].bw;
        ++mcur;
        if (mcur == 0) { std::fill(mstamp.begin(), mstamp.end(), 0); mcur = 1; }
        for (const DArc& a : g[u])
            if (a.fw < INF && a.to != v) { mstamp[a.to] = mcur; mdist[a.to] = a.fw; }
        size_t j0 = directed ? 0 : i + 1;
        for (size_t j = j0; j < n; ++j) {
            if (j == i || nb[j].fw >= INF) continue;
            int32_t x = nb[j].to;
            W c = wu + nb[j].fw;
            bool found = mstamp[x] == mcur && mdist[x] <= c;
            if (!found) {
                for (const DArc& b : g[x]) {
                    if (b.bw < INF && b.to != v && mstamp[b.to] == mcur && mdist[b.to] + b.bw <= c) {
                        found = true; break;
                    }
                }
            }
            if (!found) { ++added; added_hops += nb[i].hops + nb[j].hops; }
        }
    }
}

static vector<int32_t> level_;
static int64_t LEVEL_COEF = 1000;

static int64_t priority(int32_t v) {
    int64_t added, added_hops;
    simulate(v, added, added_hops);
    const vector<DArc>& nb = g[v];
    int64_t removed = 0, removed_hops = 0;
    for (const DArc& a : nb) {
        int c = directed ? (a.fw < INF) + (a.bw < INF) : 1;
        removed += c;
        removed_hops += (int64_t)a.hops * c;
    }
    if (removed == 0) removed = 1;
    if (removed_hops == 0) removed_hops = 1;
    return (int64_t)level_[v] * LEVEL_COEF + (1000 * added) / removed + (1000 * added_hops) / removed_hops;
}

// Upward graph, indexed by original id during construction.
struct UArc { int32_t to; W fw; W bw; };
static vector<vector<UArc>> upw;
static vector<int32_t> rank_of;   // original id -> rank

static void contract_all() {
    level_.assign((size_t)V, 0);
    wdist.assign((size_t)V, INF);
    wstamp.assign((size_t)V, 0);
    tstamp.assign((size_t)V, 0);
    mdist.assign((size_t)V, INF);
    mstamp.assign((size_t)V, 0);
    upw.assign((size_t)V, {});
    rank_of.assign((size_t)V, -1);
    vector<uint8_t> dirty((size_t)V, 0);

    vector<int64_t> pri((size_t)V);
    typedef std::pair<int64_t, int32_t> PQ;
    vector<PQ> init((size_t)V);
    for (int32_t v = 0; v < V; ++v) { pri[v] = priority(v); init[v] = {pri[v], v}; }
    std::priority_queue<PQ, vector<PQ>, std::greater<PQ>> pq(std::greater<PQ>(), std::move(init));
    tlog("initial priorities");

    int32_t next_rank = 0;
    while (!pq.empty()) {
        auto [p, v] = pq.top(); pq.pop();
        if (rank_of[v] >= 0 || p != pri[v]) continue;
        if (dirty[v]) {
            dirty[v] = 0;
            int64_t np = priority(v);
            pri[v] = np;
            if (!pq.empty() && np > pq.top().first) { pq.push({np, v}); continue; }
        }
        find_shortcuts(v, CON_SETTLE);
        rank_of[v] = next_rank++;

        vector<DArc>& nb = g[v];
        upw[v].reserve(nb.size());
        for (const DArc& a : nb) {
            upw[v].push_back({a.to, a.fw, a.bw});
            vector<DArc>& l = g[a.to];
            for (size_t k = 0; k < l.size(); ++k)
                if (l[k].to == v) { l[k] = l.back(); l.pop_back(); break; }
            level_[a.to] = std::max(level_[a.to], level_[v] + 1);
            dirty[a.to] = 1;
        }
        vector<DArc>().swap(nb);
        for (const Shortcut& s : scs) {
            if (!directed) {
                add_arc(s.u, s.x, s.w, s.w, s.hops);
                add_arc(s.x, s.u, s.w, s.w, s.hops);
            } else {
                add_arc(s.u, s.x, s.w, INF, s.hops);
                add_arc(s.x, s.u, INF, s.w, s.hops);
            }
        }
        if (g_debug && (next_rank % (V / 10) == 0)) {
            std::fprintf(stderr, "  %d/%d searches=%llu settles=%llu\n", next_rank, V,
                (unsigned long long)g_searches, (unsigned long long)g_settles);
            tlog("progress");
        }
    }
    vector<vector<DArc>>().swap(g);
    vector<W>().swap(wdist);
    vector<uint32_t>().swap(wstamp);
}

