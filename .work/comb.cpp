// solver.cpp
//
// Derived from dijkstra_foundation.cpp for the Shortest-Path Competition.
// Same command line, same input formats, same output as the foundation:
//
//     ./solver <graph_file> <query_file> <output_file>
//
// The foundation's file-format documentation still applies verbatim; see
// dijkstra_foundation.cpp / README.md.  What changed:
//
//   1. I/O.  Whole files are read into one buffer and parsed by hand; answers
//      are formatted into one buffer and written once.
//
//   2. Contraction Hierarchies (Geisberger, Sanders, Schultes, Delling 2008).
//      Vertices are contracted one by one in order of a lazily updated
//      priority (edge quotient + hop quotient + level, as in RoutingKit),
//      adding shortcuts unless a bounded witness search proves them redundant.
//      Queries are bidirectional upward Dijkstra searches with
//      stall-on-demand, over a CSR renumbered by contraction rank.
//
// Everything is C++17 and the standard library only, in this one file.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <string>
#include <vector>

using std::vector;

static const int32_t FLAG_COORDS   = 1;
static const int32_t FLAG_DIRECTED = 2;

typedef uint64_t W;
static const W INF = (~0ull) >> 2;

static bool g_debug = false;
static std::chrono::steady_clock::time_point g_t0;
static void tlog(const char* what) {
    if (!g_debug) return;
    double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - g_t0).count();
    std::fprintf(stderr, "[%8.3f s] %s\n", s, what);
}

// ---------------------------------------------------------------- input ----

struct Scanner {
    vector<char> buf;
    const char* p = nullptr;
    const char* end = nullptr;

    void load(const char* path) {
        std::FILE* f = std::fopen(path, "rb");
        if (!f) { std::fprintf(stderr, "cannot open file: %s\n", path); std::exit(1); }
        std::fseek(f, 0, SEEK_END);
        long sz = std::ftell(f);
        std::fseek(f, 0, SEEK_SET);
        if (sz < 0) sz = 0;
        buf.resize((size_t)sz + 1);
        size_t got = std::fread(buf.data(), 1, (size_t)sz, f);
        buf[got] = '\0';
        std::fclose(f);
        p = buf.data();
        end = buf.data() + got;
    }

    void release() { vector<char>().swap(buf); p = end = nullptr; }

    inline uint64_t uint_() {
        while (p < end && (*p < '0' || *p > '9')) ++p;
        uint64_t x = 0;
        while (p < end && *p >= '0' && *p <= '9') { x = x * 10 + (uint64_t)(*p - '0'); ++p; }
        return x;
    }
};

// ---------------------------------------------------------------- graph ----

static int32_t V = 0, E = 0, FLAGS = 0;
static bool directed = false;

// Dynamic graph used during contraction.  Each vertex keeps one list of
// neighbours; fw is the weight of v -> to, bw the weight of to -> v (INF when
// that direction does not exist).  Undirected graphs have fw == bw.
struct DArc { int32_t to; int32_t hops; W fw; W bw; };
static vector<vector<DArc>> g;

static void add_arc(int32_t u, int32_t v, W fw, W bw, int32_t hops) {
    for (DArc& a : g[u]) {
        if (a.to == v) {
            if (fw < a.fw) a.fw = fw;
            if (bw < a.bw) a.bw = bw;
            if (hops > a.hops) a.hops = hops;
            return;
        }
    }
    g[u].push_back({v, hops, fw, bw});
}

static int USE_PLL = 0, USE_BIDIR = 0, USE_HL = 0;
static vector<int32_t> raw_u, raw_v;
static vector<W> raw_w;

static void read_graph(const char* path) {
    Scanner sc;
    sc.load(path);
    V = (int32_t)sc.uint_();
    E = (int32_t)sc.uint_();
    FLAGS = (int32_t)sc.uint_();
    directed = (FLAGS & FLAG_DIRECTED) != 0;
    // Method choice from the graph's shape (env overrides for experiments):
    // directed without coordinates -> bidirectional Dijkstra; otherwise a
    // contraction hierarchy, answered with hub labels on road networks and
    // with bidirectional upward searches on undirected lattices (whose labels
    // are large and slow to build).
    USE_BIDIR = directed && !(FLAGS & FLAG_COORDS);
    USE_HL = directed && V <= 1000000;
    if (const char* e = std::getenv("USE_BIDIR")) USE_BIDIR = std::atoi(e);
    if (const char* e = std::getenv("USE_HL")) USE_HL = std::atoi(e);
    if (!USE_PLL && !USE_BIDIR) g.assign((size_t)V, {});
    for (int32_t i = 0; i < E; ++i) {
        int32_t u = (int32_t)sc.uint_();
        int32_t v = (int32_t)sc.uint_();
        W w = (W)sc.uint_();
        if (u == v) continue;
        if (USE_PLL || USE_BIDIR) { raw_u.push_back(u); raw_v.push_back(v); raw_w.push_back(w); continue; }
        if (directed) { add_arc(u, v, w, INF, 1); add_arc(v, u, INF, w, 1); }
        else          { add_arc(u, v, w, w, 1);   add_arc(v, u, w, w, 1); }
    }
    // Coordinates (if any) are not needed.
    sc.release();
}

// -------------------------------------------------------------- witness ----

struct MinHeap {
    struct Item { W d; int32_t v; };
    vector<Item> a;
    inline void clear() { a.clear(); }
    inline bool empty() const { return a.empty(); }
    inline W top_key() const { return a[0].d; }
    inline void push(W d, int32_t v) {
        size_t i = a.size();
        a.push_back({d, v});
        while (i) {
            size_t par = (i - 1) >> 2;
            if (a[par].d <= d) break;
            a[i] = a[par];
            i = par;
        }
        a[i] = {d, v};
    }
    inline Item pop() {
        Item best = a[0];
        Item last = a.back();
        a.pop_back();
        if (a.empty()) return best;
        size_t n = a.size(), i = 0;
        for (;;) {
            size_t c = (i << 2) + 1;
            if (c >= n) break;
            size_t lim = c + 4 < n ? c + 4 : n;
            size_t m = c;
            for (size_t j = c + 1; j < lim; ++j) if (a[j].d < a[m].d) m = j;
            if (a[m].d >= last.d) break;
            a[i] = a[m];
            i = m;
        }
        a[i] = last;
        return best;
    }
};

static int HL_PRUNE = 1;
static int SIM_SETTLE = 60, CON_SETTLE = 300, PRUNE_SETTLE = 300;
static uint64_t g_settles = 0, g_searches = 0;
static vector<W> wdist;
static vector<uint32_t> wstamp;
static uint32_t wcur = 0;
static MinHeap wheap;

static inline W wget(int32_t v) { return wstamp[v] == wcur ? wdist[v] : INF; }

static vector<uint32_t> tstamp;   // target marks for the current search
static uint32_t tcur = 0;

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
        new_targets();
        for (size_t j = j0; j < n; ++j)
            if (j != i && nb[j].fw < INF) {
                int32_t x = nb[j].to;
                tstamp[x] = tcur; tcost[x] = wu + nb[j].fw; tlist.push_back(x);
            }
        if (tlist.empty()) continue;
        witness_search(nb[i].to, v, max_settle);
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
    tcost.assign((size_t)V, 0);
    if (PRUNE_SETTLE > 0) { prune_redundant_arcs(PRUNE_SETTLE); tlog("prune"); }
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
            size_t arcs = 0, mx = 0;
            for (int32_t x = 0; x < V; ++x) if (rank_of[x] < 0) { arcs += g[x].size(); mx = std::max(mx, g[x].size()); }
            std::fprintf(stderr, "  %d/%d searches=%llu settles=%llu avgdeg=%.1f maxdeg=%zu\n", next_rank, V,
                (unsigned long long)g_searches, (unsigned long long)g_settles,
                (double)arcs / std::max(1, V - next_rank), mx);
            tlog("progress");
        }
    }
    vector<vector<DArc>>().swap(g);
    vector<W>().swap(wdist);
    vector<uint32_t>().swap(wstamp);
}

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

// ----------------------------------------------------------- hub labels ----

// Hierarchical hub labels derived from the CH (Abraham, Delling, Goldberg,
// Werneck 2012).  Labels are built from the most important vertex down: the
// forward label of r is {(r,0)} merged with fw(r->u) + L_f(u) over upward
// arcs, and an entry (h,d) is dropped when the other hubs already prove a
// strictly shorter r -> h path (checked against the finished backward label
// of h).  Entries are sorted by hub rank, so a query is a linear merge.
struct Labels {
    vector<uint64_t> start, len;   // per rank, into hub/dist (build order)
    vector<uint32_t> hub;
    vector<W> dist;
};
static Labels LF, LB;

static void build_one(int32_t r, bool forward, Labels& L, const Labels& other,
                      vector<W>& tmp, vector<uint32_t>& cand) {
    const uint32_t* H = uhead.data();
    const UArc* A = uarcs.data();
    cand.clear();
    tmp[r] = 0; cand.push_back((uint32_t)r);
    for (uint32_t k = H[r]; k < H[r + 1]; ++k) {
        W w = forward ? A[k].fw : A[k].bw;
        if (w >= INF) continue;
        int32_t u = A[k].to;
        const uint32_t* hb = L.hub.data() + L.start[u];
        const W* ds = L.dist.data() + L.start[u];
        uint64_t n = L.len[u];
        for (uint64_t i = 0; i < n; ++i) {
            uint32_t h = hb[i];
            W nd = w + ds[i];
            if (tmp[h] == INF) { cand.push_back(h); tmp[h] = nd; }
            else if (nd < tmp[h]) tmp[h] = nd;
        }
    }
    std::sort(cand.begin(), cand.end());
    L.start[r] = L.hub.size();
    for (uint32_t h : cand) {
        W d = tmp[h];
        bool keep = true;
        if (HL_PRUNE && h != (uint32_t)r) {
            const uint32_t* hb = other.hub.data() + other.start[h];
            const W* ds = other.dist.data() + other.start[h];
            uint64_t n = other.len[h];
            for (uint64_t i = 0; i < n; ++i) {
                uint32_t k2 = hb[i];
                if (k2 != h && tmp[k2] != INF && tmp[k2] + ds[i] < d) { keep = false; break; }
            }
        }
        if (keep) { L.hub.push_back(h); L.dist.push_back(d); }
    }
    L.len[r] = L.hub.size() - L.start[r];
    for (uint32_t h : cand) tmp[h] = INF;
}

// Reorders the storage so that labels are laid out by rank.
static void compact(Labels& L) {
    vector<uint32_t> hub2(L.hub.size());
    vector<W> dist2(L.dist.size());
    uint64_t pos = 0;
    for (int32_t r = 0; r < V; ++r) {
        std::memcpy(hub2.data() + pos, L.hub.data() + L.start[r], L.len[r] * sizeof(uint32_t));
        std::memcpy(dist2.data() + pos, L.dist.data() + L.start[r], L.len[r] * sizeof(W));
        L.start[r] = pos;
        pos += L.len[r];
    }
    L.hub.swap(hub2);
    L.dist.swap(dist2);
}

static void build_hub_labels() {
    for (Labels* L : {&LF, &LB}) {
        L->start.assign((size_t)V, 0);
        L->len.assign((size_t)V, 0);
        L->hub.clear(); L->dist.clear();
    }
    vector<W> tmp((size_t)V, INF);
    vector<uint32_t> cand;
    for (int32_t r = V - 1; r >= 0; --r) {
        build_one(r, true, LF, LB, tmp, cand);
        build_one(r, false, LB, LF, tmp, cand);
    }
    if (g_debug) std::fprintf(stderr, "label entries: fwd %zu bwd %zu (avg %.1f)\n",
        LF.hub.size(), LB.hub.size(), (double)(LF.hub.size() + LB.hub.size()) / (2.0 * V));
    compact(LF);
    compact(LB);
}

static inline int64_t hl_query(int32_t s, int32_t t) {
    if (s == t) return 0;
    const uint32_t* ha = LF.hub.data() + LF.start[s];
    const uint32_t* ea = ha + LF.len[s];
    const W* da = LF.dist.data() + LF.start[s];
    const uint32_t* hb = LB.hub.data() + LB.start[t];
    const uint32_t* eb = hb + LB.len[t];
    const W* db = LB.dist.data() + LB.start[t];
    W best = INF;
    while (ha < ea && hb < eb) {
        if (*ha < *hb) { ++ha; ++da; }
        else if (*ha > *hb) { ++hb; ++db; }
        else {
            W c = *da + *db;
            if (c < best) best = c;
            ++ha; ++da; ++hb; ++db;
        }
    }
    return best >= INF ? -1 : (int64_t)best;
}

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

    if (!USE_PLL && !USE_BIDIR) { dF.assign((size_t)V, INF); dB.assign((size_t)V, INF); }

    vector<char> out((size_t)Q * 21 + 64);
    char* o = out.data();

    for (int32_t i = 0; i < Q; ++i) {
        int32_t s = (int32_t)sc.uint_();
        int32_t t = (int32_t)sc.uint_();
        int64_t d = USE_BIDIR ? bidir_query(s, t) : USE_PLL ? pll_query(s, t) : USE_HL ? hl_query(rank_of[s], rank_of[t]) : ch_query(rank_of[s], rank_of[t]);
        if (d < 0) { *o++ = '-'; *o++ = '1'; }
        else o = put_u64(o, (uint64_t)d);
        *o++ = '\n';
    }

    std::FILE* f = std::fopen(opath, "wb");
    if (!f) { std::fprintf(stderr, "cannot open output file: %s\n", opath); std::exit(1); }
    std::fwrite(out.data(), 1, (size_t)(o - out.data()), f);
    std::fclose(f);
    if (g_debug) tlog("queries");
    // Nothing left to do: skip the destructors of the large containers.
    std::fflush(nullptr);
    std::_Exit(0);
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
    if (const char* e = std::getenv("PRUNE_SETTLE")) PRUNE_SETTLE = std::atoi(e);
    if (const char* e = std::getenv("LEVEL_COEF")) LEVEL_COEF = std::atoi(e);
    if (const char* e = std::getenv("USE_PLL")) USE_PLL = std::atoi(e);
    read_graph(argv[1]);
    tlog("read");
    if (USE_BIDIR) {
        build_bidir();
        tlog("bidir");
        run_queries(argv[2], argv[3]);
        tlog("queries");
        return 0;
    }
    if (USE_PLL) {
        build_pll(raw_u, raw_v, raw_w);
        tlog("pll");
        run_queries(argv[2], argv[3]);
        tlog("queries");
        return 0;
    }
    contract_all();
    tlog("contract");
    build_query_graph();
    if (g_debug) std::fprintf(stderr, "upward arcs: %u\n", uhead[V]);
    if (const char* e = std::getenv("HL_PRUNE")) HL_PRUNE = std::atoi(e);
    if (USE_HL) { build_hub_labels(); tlog("labels"); }
    run_queries(argv[2], argv[3]);
    tlog("queries");
    return 0;
}
