// solver.cpp
//
// Derived from dijkstra_foundation.cpp for the Shortest-Path Competition.
// Same command line, same input formats, same output as the foundation:
//
//     ./solver <graph_file> <query_file> <output_file>
//
// The foundation's file-format documentation still applies verbatim; see
// dijkstra_foundation.cpp / README.md.
//
// Method: hub labels over a contraction hierarchy with a PLL-labelled core.
//
//   1. Contraction (Geisberger et al. 2008).  Vertices are contracted in order
//      of a lazily updated priority, adding shortcuts that a bounded witness
//      search cannot rule out.  Contraction stops once the cheapest remaining
//      vertex is too dense; what is left is the core, which (by the usual CH
//      argument) preserves all distances between its vertices.
//   2. Core labels: Pruned Landmark Labeling (Akiba, Iwata, Yoshida 2013) on
//      the core graph, directed and weighted.
//   3. Label extension (as in hierarchical hub labels, Abraham et al. 2012):
//      walking the contraction order backwards, a vertex's label is itself
//      plus its upward neighbours' labels shifted by the arc weight, pruned
//      against the finished labels of the hubs.
//   4. A query is a merge of two hub-sorted labels.
//
// Everything is C++17 and the standard library only, in this one file.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <vector>

using std::vector;

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
static long env_int(const char* name, long def) {
    const char* e = std::getenv(name);
    return e ? std::atol(e) : def;
}

// Tunables (environment overrides exist only for experiments).
static int  CON_SETTLE   = 300;    // witness search settle limit
static int  PRUNE_SETTLE = 300;    // redundant-arc pre-pass settle limit
static long CORE_PROD    = 1L << 40; // stop contracting when in*out exceeds this
static long CORE_SIZE    = 0;      // stop contracting when this many remain
static int  LABEL_PRUNE  = 1;
static int  SPT_SAMPLES  = 0;
static long SPT_WEIGHT   = 1;

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

struct Arc { int32_t to; int32_t hops; W w; };

// Dynamic graph.  For undirected graphs only `out` is used and serves both
// directions.  Arcs to contracted vertices are deleted lazily.
static vector<vector<Arc>> out_, in_;
static vector<uint8_t> gone;       // contracted (not core)
static vector<int32_t> dead_;      // lazily deleted entries per vertex (upper bound)

static inline vector<Arc>& OUT(int32_t v) { return out_[v]; }
static inline vector<Arc>& IN(int32_t v) { return directed ? in_[v] : out_[v]; }

static void compact_list(vector<Arc>& l) {
    size_t k = 0;
    for (size_t i = 0; i < l.size(); ++i) if (!gone[l[i].to]) l[k++] = l[i];
    l.resize(k);
}

// Adds (or improves) u -> x.
static void add_arc(int32_t u, int32_t x, W w, int32_t hops) {
    vector<Arc>& o = OUT(u);
    if (o.size() <= 64) {
        for (Arc& a : o) {
            if (a.to == x && !gone[x]) {
                if (w < a.w) {
                    a.w = w; a.hops = hops;
                    vector<Arc>& i2 = IN(x);
                    bool found = false;
                    if (i2.size() <= 64 || !directed) {
                        for (Arc& b : i2) if (b.to == u) { if (w < b.w) { b.w = w; b.hops = hops; } found = true; break; }
                    }
                    if (!found && directed) i2.push_back({u, hops, w});
                }
                return;
            }
        }
    }
    o.push_back({x, hops, w});
    if (directed) IN(x).push_back({u, hops, w});
}

static void read_graph(const char* path) {
    Scanner sc;
    sc.load(path);
    V = (int32_t)sc.uint_();
    E = (int32_t)sc.uint_();
    FLAGS = (int32_t)sc.uint_();
    directed = (FLAGS & FLAG_DIRECTED) != 0;
    out_.assign((size_t)V, {});
    if (directed) in_.assign((size_t)V, {});
    gone.assign((size_t)V, 0);
    dead_.assign((size_t)V, 0);
    for (int32_t i = 0; i < E; ++i) {
        int32_t u = (int32_t)sc.uint_();
        int32_t v = (int32_t)sc.uint_();
        W w = (W)sc.uint_();
        if (u == v) continue;
        add_arc(u, v, w, 1);
        if (!directed) add_arc(v, u, w, 1);
    }
    sc.release();
}

// ----------------------------------------------------------------- heap ----

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

// -------------------------------------------------------------- witness ----

struct WNode { W dist; W tcost; uint32_t stamp; uint32_t tstamp; };
static vector<WNode> wn;
static uint32_t wcur = 0, tcur = 0;
static MinHeap wheap;
static vector<int32_t> tlist;
static uint64_t g_settles = 0, g_searches = 0;

static inline W wget(int32_t v) { return wn[v].stamp == wcur ? wn[v].dist : INF; }

static inline void new_targets() {
    if (++tcur == 0) { for (WNode& x : wn) x.tstamp = 0; tcur = 1; }
    tlist.clear();
}
static inline void add_target(int32_t x, W c) {
    if (wn[x].tstamp == tcur) { if (c > wn[x].tcost) wn[x].tcost = c; return; }
    wn[x].tstamp = tcur; wn[x].tcost = c; tlist.push_back(x);
}

// Dijkstra from src over live out-arcs avoiding `avoid`; stops when all
// targets are settled, the key passes the largest pending threshold, or
// max_settle vertices have been settled.
static void witness_search(int32_t src, int32_t avoid, int max_settle) {
    ++g_searches;
    if (++wcur == 0) { for (WNode& x : wn) x.stamp = 0; wcur = 1; }
    W limit = 0;
    for (int32_t x : tlist) limit = std::max(limit, wn[x].tcost);
    int pending = (int)tlist.size();
    wheap.clear();
    wn[src].dist = 0; wn[src].stamp = wcur;
    wheap.push(0, src);
    int settled = 0;
    while (!wheap.empty()) {
        MinHeap::Item it = wheap.pop();
        if (it.d > wn[it.v].dist) continue;
        if (it.d > limit) break;
        if (wn[it.v].tstamp == tcur) {
            wn[it.v].tstamp = 0;
            if (--pending <= 0) break;
            if (wn[it.v].tcost >= limit) {
                limit = 0;
                for (int32_t x : tlist) if (wn[x].tstamp == tcur) limit = std::max(limit, wn[x].tcost);
                if (it.d > limit) break;
            }
        }
        if (++settled > max_settle) break;
        ++g_settles;
        for (const Arc& a : OUT(it.v)) {
            if (a.to == avoid || gone[a.to]) continue;
            W nd = it.d + a.w;
            if (nd > limit) continue;
            WNode& n = wn[a.to];
            if (n.stamp != wcur || nd < n.dist) {
                n.stamp = wcur; n.dist = nd;
                wheap.push(nd, a.to);
            }
        }
    }
}

// Removes arcs u -> x for which a strictly shorter u -> x path exists.
static void prune_redundant_arcs() {
    size_t removed = 0;
    vector<int32_t> tg;
    for (int32_t u = 0; u < V; ++u) {
        vector<Arc>& o = OUT(u);
        if (o.size() < 2) continue;
        new_targets();
        for (const Arc& a : o) if (a.w > 0) add_target(a.to, a.w - 1);
        tg = tlist;
        witness_search(u, -1, PRUNE_SETTLE);
        for (int32_t x : tg) {
            W best = wget(x);
            // find arc weight
            for (size_t k = 0; k < o.size(); ++k) {
                if (o[k].to != x) continue;
                if (best < o[k].w) {
                    W wv = o[k].w;
                    o[k] = o.back(); o.pop_back();
                    vector<Arc>& i2 = IN(x);
                    for (size_t j = 0; j < i2.size(); ++j)
                        if (i2[j].to == u && i2[j].w == wv) { i2[j] = i2.back(); i2.pop_back(); break; }
                    ++removed;
                }
                break;
            }
        }
    }
    if (g_debug) std::fprintf(stderr, "pruned arcs: %zu\n", removed);
}

// ---------------------------------------------------------- contraction ----

struct Shortcut { int32_t u, x; W w; int32_t hops; };
static vector<Shortcut> scs;

static inline void live_compact(int32_t v) {
    if (dead_[v] > 8 && (size_t)dead_[v] * 2 > OUT(v).size()) {
        compact_list(OUT(v));
        if (directed) compact_list(IN(v));
        dead_[v] = 0;
    }
}

static void find_shortcuts(int32_t v) {
    scs.clear();
    vector<Arc>& in = IN(v);
    vector<Arc>& out = OUT(v);
    for (size_t i = 0; i < in.size(); ++i) {
        const Arc& ai = in[i];
        if (gone[ai.to]) continue;
        new_targets();
        size_t j0 = directed ? 0 : i + 1;
        for (size_t j = j0; j < out.size(); ++j) {
            const Arc& aj = out[j];
            if (gone[aj.to] || aj.to == ai.to) continue;
            add_target(aj.to, ai.w + aj.w);
        }
        if (tlist.empty()) continue;
        witness_search(ai.to, v, CON_SETTLE);
        for (size_t j = j0; j < out.size(); ++j) {
            const Arc& aj = out[j];
            if (gone[aj.to] || aj.to == ai.to) continue;
            W w = ai.w + aj.w;
            if (wget(aj.to) > w) scs.push_back({ai.to, aj.to, w, ai.hops + aj.hops});
        }
    }
}

// Two-hop estimate of the number of shortcuts, no heap.
static vector<W> mdist;
static vector<uint32_t> mstamp;
static uint32_t mcur = 0;

static void simulate(int32_t v, int64_t& added, int64_t& added_hops, int64_t& removed, int64_t& removed_hops) {
    added = added_hops = removed = removed_hops = 0;
    vector<Arc>& in = IN(v);
    vector<Arc>& out = OUT(v);
    for (const Arc& a : out) if (!gone[a.to]) { ++removed; removed_hops += a.hops; }
    if (directed) for (const Arc& a : in) if (!gone[a.to]) { ++removed; removed_hops += a.hops; }
    for (size_t i = 0; i < in.size(); ++i) {
        const Arc& ai = in[i];
        if (gone[ai.to]) continue;
        int32_t u = ai.to;
        if (++mcur == 0) { std::fill(mstamp.begin(), mstamp.end(), 0); mcur = 1; }
        // one and two hops from u, avoiding v
        vector<Arc>& ou = OUT(u);
        size_t budget = 256;
        for (const Arc& a : ou) {
            if (a.to == v || gone[a.to]) continue;
            if (mstamp[a.to] != mcur || a.w < mdist[a.to]) { mstamp[a.to] = mcur; mdist[a.to] = a.w; }
        }
        for (const Arc& a : ou) {
            if (a.to == v || gone[a.to]) continue;
            vector<Arc>& ow = OUT(a.to);
            if (ow.size() > budget) continue;
            budget -= ow.size();
            for (const Arc& b : ow) {
                if (b.to == v || gone[b.to]) continue;
                W d2 = a.w + b.w;
                if (mstamp[b.to] != mcur || d2 < mdist[b.to]) { mstamp[b.to] = mcur; mdist[b.to] = d2; }
            }
        }
        size_t j0 = directed ? 0 : i + 1;
        for (size_t j = j0; j < out.size(); ++j) {
            const Arc& aj = out[j];
            if (gone[aj.to] || aj.to == u) continue;
            W c = ai.w + aj.w;
            if (!(mstamp[aj.to] == mcur && mdist[aj.to] <= c)) { ++added; added_hops += ai.hops + aj.hops; }
        }
    }
}

static vector<int32_t> level_;

static inline int64_t live_prod(int32_t v) {
    int64_t o = 0, i = 0;
    for (const Arc& a : OUT(v)) o += !gone[a.to];
    if (!directed) return o * o;
    for (const Arc& a : IN(v)) i += !gone[a.to];
    return o * i;
}

static int64_t priority(int32_t v) {
    int64_t prod = live_prod(v);
    if (prod > 4096) return ((int64_t)1 << 40) + prod;
    int64_t added, ah, removed, rh;
    simulate(v, added, ah, removed, rh);
    if (removed == 0) removed = 1;
    if (rh == 0) rh = 1;
    return (int64_t)level_[v] * 1000 + (1000 * added) / removed + (1000 * ah) / rh;
}

// Upward arcs of contracted vertices (by original id).
static vector<vector<Arc>> up_out, up_in;
static vector<int32_t> con_order;   // contraction order
static vector<int32_t> core;        // remaining vertices

static void contract() {
    level_.assign((size_t)V, 0);
    wn.assign((size_t)V, WNode{INF, 0, 0, 0});
    mdist.assign((size_t)V, INF);
    mstamp.assign((size_t)V, 0);
    if (PRUNE_SETTLE > 0) { prune_redundant_arcs(); tlog("prune"); }

    up_out.assign((size_t)V, {});
    if (directed) up_in.assign((size_t)V, {});
    vector<uint8_t> dirty((size_t)V, 0);
    vector<int64_t> pri((size_t)V);
    typedef std::pair<int64_t, int32_t> PQ;
    vector<PQ> init((size_t)V);
    for (int32_t v = 0; v < V; ++v) { pri[v] = priority(v); init[v] = {pri[v], v}; }
    std::priority_queue<PQ, vector<PQ>, std::greater<PQ>> pq(std::greater<PQ>(), std::move(init));
    tlog("initial priorities");

    int32_t remaining = V;
    while (!pq.empty()) {
        auto [p, v] = pq.top();
        if (gone[v] || p != pri[v]) { pq.pop(); continue; }
        if (dirty[v]) {
            pq.pop();
            dirty[v] = 0;
            live_compact(v);
            int64_t np = priority(v);
            pri[v] = np;
            pq.push({np, v});
            continue;
        }
        if (remaining <= CORE_SIZE) break;
        if (live_prod(v) > CORE_PROD) break;
        pq.pop();

        live_compact(v);
        find_shortcuts(v);
        gone[v] = 1;
        --remaining;
        con_order.push_back(v);
        for (const Arc& a : OUT(v)) if (!gone[a.to]) up_out[v].push_back(a);
        if (directed) for (const Arc& a : IN(v)) if (!gone[a.to]) up_in[v].push_back(a);
        for (const Arc& a : OUT(v)) if (!gone[a.to]) {
            ++dead_[a.to]; dirty[a.to] = 1;
            level_[a.to] = std::max(level_[a.to], level_[v] + 1);
        }
        if (directed) for (const Arc& a : IN(v)) if (!gone[a.to]) {
            ++dead_[a.to]; dirty[a.to] = 1;
            level_[a.to] = std::max(level_[a.to], level_[v] + 1);
        }
        vector<Arc>().swap(OUT(v));
        if (directed) vector<Arc>().swap(IN(v));
        for (const Shortcut& s : scs) {
            add_arc(s.u, s.x, s.w, s.hops);
            if (!directed) add_arc(s.x, s.u, s.w, s.hops);
        }
        if (g_debug && (int32_t)con_order.size() % (V / 10 > 0 ? V / 10 : 1) == 0) {
            std::fprintf(stderr, "  %zu contracted, searches=%llu settles=%llu\n", con_order.size(),
                (unsigned long long)g_searches, (unsigned long long)g_settles);
            tlog("progress");
        }
    }
    for (int32_t v = 0; v < V; ++v) if (!gone[v]) core.push_back(v);
    if (g_debug) std::fprintf(stderr, "contracted %zu, core %zu\n", con_order.size(), core.size());
    vector<WNode>().swap(wn);
    vector<W>().swap(mdist);
    vector<uint32_t>().swap(mstamp);
}

// ---------------------------------------------------------------- labels ----

// Hub ids: core vertices by PLL order get 0..C-1, contracted vertices get
// C + (position counted from the end of the contraction order).
static vector<uint32_t> hub_id;     // original id -> hub id

struct Labels {
    vector<uint64_t> start;
    vector<uint32_t> len;
    vector<uint32_t> hub;
    vector<W> dist;
};
static Labels LO, LI;   // LO[v]: v -> hub, LI[v]: hub -> v   (indexed by hub id)

// PLL on the core.  Label storage is appended per vertex into vectors first.
static void core_pll(vector<vector<std::pair<uint32_t, W>>>& co, vector<vector<std::pair<uint32_t, W>>>& ci) {
    size_t C = core.size();
    // local ids
    vector<int32_t> lid((size_t)V, -1);
    for (size_t i = 0; i < C; ++i) lid[core[i]] = (int32_t)i;
    vector<uint32_t> oh(C + 1, 0), ih(C + 1, 0);
    for (size_t i = 0; i < C; ++i) {
        for (const Arc& a : OUT(core[i])) if (!gone[a.to]) { ++oh[i + 1]; if (directed) ++ih[lid[a.to] + 1]; }
    }
    for (size_t i = 0; i < C; ++i) { oh[i + 1] += oh[i]; ih[i + 1] += ih[i]; }
    vector<uint32_t> oto(oh[C]), ito(directed ? ih[C] : 0);
    vector<W> ow(oh[C]), iw(directed ? ih[C] : 0);
    {
        vector<uint32_t> c1(oh.begin(), oh.end() - 1), c2(ih.begin(), ih.end() - 1);
        for (size_t i = 0; i < C; ++i)
            for (const Arc& a : OUT(core[i])) if (!gone[a.to]) {
                uint32_t k = c1[i]++; oto[k] = (uint32_t)lid[a.to]; ow[k] = a.w;
                if (directed) { uint32_t j = c2[lid[a.to]]++; ito[j] = (uint32_t)i; iw[j] = a.w; }
            }
    }
    const vector<uint32_t>& IH = directed ? ih : oh;
    const vector<uint32_t>& ITO = directed ? ito : oto;
    const vector<W>& IW = directed ? iw : ow;

    vector<uint32_t> order(C);
    for (size_t i = 0; i < C; ++i) order[i] = (uint32_t)i;
    // Score: degree product, plus descendant counts in sampled shortest-path
    // trees (both directions) -- vertices covering many shortest paths first.
    vector<uint64_t> score(C);
    for (size_t i = 0; i < C; ++i)
        score[i] = (uint64_t)(oh[i + 1] - oh[i] + 1) * (IH[i + 1] - IH[i] + 1);
    if (SPT_SAMPLES > 0 && C > 0) {
        vector<W> dd(C, INF);
        vector<uint32_t> par(C), seq, sub(C);
        vector<uint32_t> tch;
        MinHeap hp;
        uint64_t rng = 0x9E3779B97F4A7C15ull;
        for (int smp = 0; smp < SPT_SAMPLES; ++smp) {
            rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
            uint32_t root = (uint32_t)(rng % C);
            for (int pass = 0; pass < (directed ? 2 : 1); ++pass) {
                const uint32_t* H = pass == 0 ? oh.data() : IH.data();
                const uint32_t* TO = pass == 0 ? oto.data() : ITO.data();
                const W* WT = pass == 0 ? ow.data() : IW.data();
                seq.clear();
                hp.clear();
                dd[root] = 0; par[root] = root; tch.push_back(root);
                hp.push(0, (int32_t)root);
                while (!hp.empty()) {
                    MinHeap::Item it = hp.pop();
                    uint32_t u = (uint32_t)it.v;
                    if (it.d > dd[u]) continue;
                    seq.push_back(u);
                    for (uint32_t a = H[u]; a < H[u + 1]; ++a) {
                        uint32_t x = TO[a];
                        W nd = it.d + WT[a];
                        if (nd < dd[x]) {
                            if (dd[x] == INF) tch.push_back(x);
                            dd[x] = nd; par[x] = u;
                            hp.push(nd, (int32_t)x);
                        }
                    }
                }
                for (uint32_t u : seq) sub[u] = 1;
                for (size_t i = seq.size(); i-- > 1;) sub[par[seq[i]]] += sub[seq[i]];
                for (uint32_t u : seq) score[u] += (uint64_t)sub[u] * SPT_WEIGHT;
                for (uint32_t x : tch) dd[x] = INF;
                tch.clear();
            }
        }
    }
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        return score[a] != score[b] ? score[a] > score[b] : a < b;
    });
    for (size_t k = 0; k < C; ++k) hub_id[core[order[k]]] = (uint32_t)k;

    // labels in local ids: lo[i] = (hub k, dist i -> hub), li[i] = (k, hub -> i)
    co.assign(C, {}); ci.assign(C, {});
    vector<W> T(C, INF), dist(C, INF);
    vector<uint32_t> touched;
    MinHeap heap;
    for (size_t k = 0; k < C; ++k) {
        uint32_t r = order[k];
        for (int pass = 0; pass < 2; ++pass) {
            const bool fwd = pass == 0;
            auto& own = fwd ? co[r] : ci[r];
            for (auto& e : own) T[e.first] = e.second;
            auto& target = (fwd && directed) ? ci : co;
            const uint32_t* H = fwd ? oh.data() : IH.data();
            const uint32_t* TO = fwd ? oto.data() : ITO.data();
            const W* WT = fwd ? ow.data() : IW.data();
            heap.clear();
            dist[r] = 0; touched.push_back(r);
            heap.push(0, (int32_t)r);
            while (!heap.empty()) {
                MinHeap::Item it = heap.pop();
                uint32_t u = (uint32_t)it.v; W d = it.d;
                if (d > dist[u]) continue;
                bool pruned = false;
                for (auto& e : target[u]) {
                    W t = T[e.first];
                    if (t != INF && t + e.second <= d) { pruned = true; break; }
                }
                if (pruned) continue;
                target[u].push_back({(uint32_t)k, d});
                for (uint32_t a = H[u]; a < H[u + 1]; ++a) {
                    uint32_t x = TO[a];
                    W nd = d + WT[a];
                    if (nd < dist[x]) {
                        if (dist[x] == INF) touched.push_back(x);
                        dist[x] = nd;
                        heap.push(nd, (int32_t)x);
                    }
                }
            }
            for (uint32_t x : touched) dist[x] = INF;
            touched.clear();
            for (auto& e : own) T[e.first] = INF;
            if (!directed) {
                // symmetric: backward label equals forward label
                break;
            }
        }
    }
}

static void build_labels() {
    size_t C = core.size();
    size_t N = con_order.size();
    hub_id.assign((size_t)V, 0);
    for (size_t i = 0; i < N; ++i) hub_id[con_order[i]] = (uint32_t)(C + (N - 1 - i));

    vector<vector<std::pair<uint32_t, W>>> co, ci;
    core_pll(co, ci);
    tlog("core pll");
    if (g_debug) {
        size_t a = 0, b = 0;
        for (auto& l : co) a += l.size();
        for (auto& l : ci) b += l.size();
        std::fprintf(stderr, "core label entries: out %zu in %zu (avg %.1f)\n", a, b, C ? (double)(a + b) / (2.0 * C) : 0.0);
    }

    for (Labels* L : {&LO, &LI}) {
        L->start.assign((size_t)V, 0);
        L->len.assign((size_t)V, 0);
        L->hub.clear(); L->dist.clear();
    }
    // core labels, stored under hub id
    for (size_t i = 0; i < C; ++i) {
        uint32_t h = hub_id[core[i]];
        LO.start[h] = LO.hub.size();
        for (auto& e : co[i]) { LO.hub.push_back(e.first); LO.dist.push_back(e.second); }
        LO.len[h] = (uint32_t)co[i].size();
        if (directed) {
            LI.start[h] = LI.hub.size();
            for (auto& e : ci[i]) { LI.hub.push_back(e.first); LI.dist.push_back(e.second); }
            LI.len[h] = (uint32_t)ci[i].size();
        }
    }
    vector<vector<std::pair<uint32_t, W>>>().swap(co);
    vector<vector<std::pair<uint32_t, W>>>().swap(ci);

    vector<W> tmp((size_t)V, INF);   // indexed by hub id
    vector<uint32_t> cand;
    for (size_t idx = N; idx-- > 0;) {
        int32_t v = con_order[idx];
        uint32_t hv = hub_id[v];
        for (int pass = 0; pass < 2; ++pass) {
            const bool fwd = pass == 0;
            if (!directed && !fwd) break;
            Labels& L = fwd ? LO : LI;
            const Labels& other = (fwd && directed) ? LI : LO;
            const vector<Arc>& ups = fwd ? up_out[v] : (directed ? up_in[v] : up_out[v]);
            cand.clear();
            tmp[hv] = 0; cand.push_back(hv);
            for (const Arc& a : ups) {
                uint32_t hu = hub_id[a.to];
                const uint32_t* hb = L.hub.data() + L.start[hu];
                const W* ds = L.dist.data() + L.start[hu];
                uint32_t n = L.len[hu];
                for (uint32_t i = 0; i < n; ++i) {
                    uint32_t h = hb[i];
                    W nd = a.w + ds[i];
                    if (tmp[h] == INF) { cand.push_back(h); tmp[h] = nd; }
                    else if (nd < tmp[h]) tmp[h] = nd;
                }
            }
            std::sort(cand.begin(), cand.end());
            L.start[hv] = L.hub.size();
            for (uint32_t h : cand) {
                W d = tmp[h];
                bool keep = true;
                if (LABEL_PRUNE && h != hv) {
                    const uint32_t* hb = other.hub.data() + other.start[h];
                    const W* ds = other.dist.data() + other.start[h];
                    uint32_t n = other.len[h];
                    for (uint32_t i = 0; i < n; ++i) {
                        uint32_t k2 = hb[i];
                        if (k2 != h && tmp[k2] != INF && tmp[k2] + ds[i] < d) { keep = false; break; }
                    }
                }
                if (keep) { L.hub.push_back(h); L.dist.push_back(d); }
            }
            L.len[hv] = (uint32_t)(L.hub.size() - L.start[hv]);
            for (uint32_t h : cand) tmp[h] = INF;
        }
    }
    if (g_debug) std::fprintf(stderr, "label entries: out %zu in %zu (avg %.1f)\n",
        LO.hub.size(), LI.hub.size(), (double)(LO.hub.size() + (directed ? LI.hub.size() : LO.hub.size())) / (2.0 * V));
}

static inline int64_t label_query(int32_t s, int32_t t) {
    if (s == t) return 0;
    uint32_t a = hub_id[s], b = hub_id[t];
    const uint32_t* ha = LO.hub.data() + LO.start[a];
    const uint32_t* ea = ha + LO.len[a];
    const W* da = LO.dist.data() + LO.start[a];
    const Labels& R = directed ? LI : LO;
    const uint32_t* hb = R.hub.data() + R.start[b];
    const uint32_t* eb = hb + R.len[b];
    const W* db = R.dist.data() + R.start[b];
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

// ------------------------------------------------ bidirectional dijkstra ----

// For the directed power-law family: out-degree is tiny, and a balanced
// bidirectional search meets long before either side reaches the hubs'
// huge in-lists.  Adjacency is sorted by weight so a scan can stop at the
// first arc that cannot lead below the current best meeting distance.
static vector<uint32_t> fh, bh;
static vector<int32_t> fto, bto;
static vector<W> fw, bw;
static vector<W> bdF, bdB;
static vector<int32_t> bdTF, bdTB;
static MinHeap bhF, bhB;
static int BD_BAL = 0;
static inline W alt_h(int32_t v, int32_t t);
static int NLM = 16;
static int LM_HUBS = 0, BD_UB = 0, BD_LB = 0;
static vector<W> lmF, lmB;
static bool g_bidir = false, g_alt = false, g_alt_bidir = false;
static uint64_t g_bsF = 0, g_bsB = 0, g_bwF = 0, g_bwB = 0;

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

static void bidir_prepare(const vector<int32_t>& eu, const vector<int32_t>& ev, const vector<W>& ew) {
    build_csr(eu, ev, ew, fh, fto, fw);
    build_csr(ev, eu, ew, bh, bto, bw);
    bdF.assign((size_t)V, INF);
    bdB.assign((size_t)V, INF);
}

static int64_t bidir_query(int32_t s, int32_t t) {
    if (s == t) return 0;
    if (fh[s] == fh[s + 1] || bh[t] == bh[t + 1]) return -1;
    W mu = INF;
    if (BD_UB) {
        const W* bs = &lmB[(size_t)s * NLM];
        const W* ft = &lmF[(size_t)t * NLM];
        for (int i = 0; i < NLM; ++i) if (bs[i] < INF && ft[i] < INF && bs[i] + ft[i] < mu) mu = bs[i] + ft[i];
        if (BD_LB) {
            W lb = alt_h(s, t);
            if (lb >= INF) return -1;
            if (lb >= mu) return (int64_t)mu;
        }
    }
    bhF.clear(); bhB.clear();
    bdF[s] = 0; bdTF.push_back(s); bhF.push(0, s);
    bdB[t] = 0; bdTB.push_back(t); bhB.push(0, t);
    uint64_t workF = 0, workB = 0;
    while (!bhF.empty() && !bhB.empty()) {
        W kf = bhF.top_key(), kb = bhB.top_key();
        if (kf + kb >= mu) break;
        bool fwd;
        if (BD_BAL == 3) {
            int32_t tf = bhF.a[0].v, tb = bhB.a[0].v;
            fwd = workF + (fh[tf + 1] - fh[tf]) <= workB + (bh[tb + 1] - bh[tb]);
        } else fwd = BD_BAL == 0 ? workF <= workB : (BD_BAL == 1 ? bhF.a.size() <= bhB.a.size() : kf <= kb);
        MinHeap& hp = fwd ? bhF : bhB;
        vector<W>& d1 = fwd ? bdF : bdB;
        const vector<W>& d2 = fwd ? bdB : bdF;
        vector<int32_t>& tch = fwd ? bdTF : bdTB;
        const uint32_t* H = fwd ? fh.data() : bh.data();
        const int32_t* TO = fwd ? fto.data() : bto.data();
        const W* WT = fwd ? fw.data() : bw.data();
        W other = fwd ? kb : kf;
        MinHeap::Item it = hp.pop();
        int32_t u = it.v; W d = it.d;
        if (d > d1[u]) continue;
        uint32_t k = H[u], e = H[u + 1];
        for (; k < e; ++k) {
            W nd = d + WT[k];
            if (nd >= mu) break;               // sorted: no later arc can help
            int32_t x = TO[k];
            if (nd < d1[x]) {
                if (d1[x] == INF) tch.push_back(x);
                d1[x] = nd;
                if (d2[x] < INF && nd + d2[x] < mu) mu = nd + d2[x];
                // A path through x is at least nd + (other side's key) unless
                // the other side already settled x, which mu now covers.
                if (nd + other < mu && x != (fwd ? t : s) &&
                    (fwd ? fh[x] != fh[x + 1] : bh[x] != bh[x + 1])) hp.push(nd, x);
            }
        }
        (fwd ? workF : workB) += (k - H[u]) + 1;
        if (g_debug) { (fwd ? g_bsF : g_bsB) += 1; }
    }
    if (g_debug) { g_bwF += workF; g_bwB += workB; }
    for (int32_t v : bdTF) bdF[v] = INF;
    for (int32_t v : bdTB) bdB[v] = INF;
    bdTF.clear(); bdTB.clear();
    return mu >= INF ? -1 : (int64_t)mu;
}

// ------------------------------------------------------------------ ALT ----

// Forward A* with landmark lower bounds (Goldberg & Harrelson 2005):
//   dist(v,t) >= d(L,t) - d(L,v)   and   dist(v,t) >= d(v,L) - d(t,L).
// Landmark distances also certify unreachability: if L reaches v but not t,
// or t reaches L but v does not, no v -> t path exists.
// lmF[v*NLM+i] = d(L_i, v), lmB[v*NLM+i] = d(v, L_i)

static void full_dijkstra(int32_t src, const vector<uint32_t>& h, const vector<int32_t>& to,
                          const vector<W>& w, W* outd, int stride, int idx) {
    MinHeap hp;
    for (int32_t v = 0; v < V; ++v) outd[(size_t)v * stride + idx] = INF;
    outd[(size_t)src * stride + idx] = 0;
    hp.push(0, src);
    while (!hp.empty()) {
        MinHeap::Item it = hp.pop();
        if (it.d > outd[(size_t)it.v * stride + idx]) continue;
        for (uint32_t k = h[it.v]; k < h[it.v + 1]; ++k) {
            W nd = it.d + w[k];
            W& dx = outd[(size_t)to[k] * stride + idx];
            if (nd < dx) { dx = nd; hp.push(nd, to[k]); }
        }
    }
}

static void alt_prepare() {
    lmF.assign((size_t)V * NLM, INF);
    lmB.assign((size_t)V * NLM, INF);
    // Farthest-first selection on the sum of forward and backward distances.
    if (LM_HUBS) {
        vector<int32_t> ids((size_t)V);
        for (int32_t v = 0; v < V; ++v) ids[v] = v;
        std::partial_sort(ids.begin(), ids.begin() + NLM, ids.end(), [&](int32_t a, int32_t b) {
            uint32_t da = bh[a + 1] - bh[a], db = bh[b + 1] - bh[b];
            return da != db ? da > db : a < b;
        });
        for (int i = 0; i < NLM; ++i) {
            full_dijkstra(ids[i], fh, fto, fw, lmF.data(), NLM, i);
            full_dijkstra(ids[i], bh, bto, bw, lmB.data(), NLM, i);
        }
        return;
    }
    vector<W> score((size_t)V, 0);
    int32_t cur = 0;
    for (int32_t v = 0; v < V; ++v)
        if ((fh[v + 1] - fh[v]) * (bh[v + 1] - bh[v]) > (fh[cur + 1] - fh[cur]) * (bh[cur + 1] - bh[cur])) cur = v;
    for (int i = 0; i < NLM; ++i) {
        full_dijkstra(cur, fh, fto, fw, lmF.data(), NLM, i);
        full_dijkstra(cur, bh, bto, bw, lmB.data(), NLM, i);
        int32_t best = -1; W bs = 0;
        for (int32_t v = 0; v < V; ++v) {
            W a = lmF[(size_t)v * NLM + i], b = lmB[(size_t)v * NLM + i];
            if (a >= INF || b >= INF) { score[v] = 0; continue; }
            if (i == 0) score[v] = a + b; else score[v] = std::min(score[v], a + b);
            if (score[v] > bs) { bs = score[v]; best = v; }
        }
        if (best < 0) break;
        cur = best;
    }
}

static vector<W> aD;
static vector<W> aH;
static vector<int32_t> aT;
static MinHeap aHeap;

static inline W alt_h(int32_t v, int32_t t) {
    const W* fv = &lmF[(size_t)v * NLM];
    const W* ft = &lmF[(size_t)t * NLM];
    const W* bv = &lmB[(size_t)v * NLM];
    const W* bt = &lmB[(size_t)t * NLM];
    W h = 0;
    for (int i = 0; i < NLM; ++i) {
        if (ft[i] < INF) {
            if (fv[i] >= INF) { /* no info */ }
            else if (ft[i] > fv[i] && ft[i] - fv[i] > h) h = ft[i] - fv[i];
        } else if (fv[i] < INF) return INF;       // L reaches v but not t
        if (bt[i] < INF) {
            if (bv[i] >= INF) return INF;          // t reaches L, v does not
            if (bv[i] > bt[i] && bv[i] - bt[i] > h) h = bv[i] - bt[i];
        }
    }
    return h;
}

static uint64_t g_alt_settled = 0;

static int64_t alt_query(int32_t s, int32_t t) {
    if (s == t) return 0;
    if (fh[s] == fh[s + 1] || bh[t] == bh[t + 1]) return -1;
    W hs = alt_h(s, t);
    if (hs >= INF) return -1;
    aHeap.clear();
    aD[s] = 0; aH[s] = hs; aT.push_back(s);
    aHeap.push(hs, s);
    W ans = INF;
    while (!aHeap.empty()) {
        MinHeap::Item it = aHeap.pop();
        int32_t u = it.v;
        W du = aD[u];
        if (it.d != du + aH[u]) continue;
        if (u == t) { ans = du; break; }
        ++g_alt_settled;
        for (uint32_t k = fh[u]; k < fh[u + 1]; ++k) {
            int32_t x = fto[k];
            W nd = du + fw[k];
            if (aD[x] == INF) {
                W hx = alt_h(x, t);
                aT.push_back(x);
                aH[x] = hx;
                aD[x] = nd;
                if (hx < INF) aHeap.push(nd + hx, x);
            } else if (nd < aD[x]) {
                aD[x] = nd;
                if (aH[x] < INF) aHeap.push(nd + aH[x], x);
            }
        }
    }
    for (int32_t v : aT) aD[v] = INF;
    aT.clear();
    return ans >= INF ? -1 : (int64_t)ans;
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
    vector<char> out((size_t)Q * 21 + 64);
    char* o = out.data();
    vector<int32_t> qs((size_t)Q), qt((size_t)Q);
    for (int32_t i = 0; i < Q; ++i) { qs[i] = (int32_t)sc.uint_(); qt[i] = (int32_t)sc.uint_(); }
    vector<uint32_t> tcnt;
    double tb[2] = {0, 0}; uint64_t nb[2] = {0, 0};
    if (g_debug) { tcnt.assign((size_t)V, 0); for (int32_t i = 0; i < Q; ++i) ++tcnt[qt[i]]; }
    for (int32_t i = 0; i < Q; ++i) {
        int32_t s = qs[i];
        int32_t t = qt[i];
        auto q0 = g_debug ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
        int64_t d = (g_alt && !g_alt_bidir) ? alt_query(s, t) : g_bidir ? bidir_query(s, t) : label_query(s, t);
        if (g_debug) {
            int c = tcnt[t] >= 8;
            tb[c] += std::chrono::duration<double>(std::chrono::steady_clock::now() - q0).count(); ++nb[c];
        }
        if (d < 0) { *o++ = '-'; *o++ = '1'; }
        else o = put_u64(o, (uint64_t)d);
        *o++ = '\n';
    }
    if (g_debug) std::fprintf(stderr, "bidir settles F %llu B %llu, arcs F %llu B %llu\n",
        (unsigned long long)g_bsF, (unsigned long long)g_bsB, (unsigned long long)g_bwF, (unsigned long long)g_bwB);
    if (g_debug) std::fprintf(stderr, "queries: rare targets %llu in %.3f s, frequent targets %llu in %.3f s\n",
        (unsigned long long)nb[0], tb[0], (unsigned long long)nb[1], tb[1]);
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
    CON_SETTLE = (int)env_int("CON_SETTLE", CON_SETTLE);
    PRUNE_SETTLE = (int)env_int("PRUNE_SETTLE", PRUNE_SETTLE);
    CORE_PROD = env_int("CORE_PROD", CORE_PROD);
    CORE_SIZE = env_int("CORE_SIZE", CORE_SIZE);
    LABEL_PRUNE = (int)env_int("LABEL_PRUNE", LABEL_PRUNE);
    SPT_SAMPLES = (int)env_int("SPT_SAMPLES", SPT_SAMPLES);
    SPT_WEIGHT = env_int("SPT_WEIGHT", SPT_WEIGHT);

    BD_BAL = (int)env_int("BD_BAL", BD_BAL);
    if (env_int("USE_BIDIR", 0)) {
        Scanner sc;
        sc.load(argv[1]);
        V = (int32_t)sc.uint_(); E = (int32_t)sc.uint_(); FLAGS = (int32_t)sc.uint_();
        directed = (FLAGS & FLAG_DIRECTED) != 0;
        vector<int32_t> eu, ev; vector<W> ew;
        eu.reserve((size_t)E * (directed ? 1 : 2)); ev.reserve(eu.capacity()); ew.reserve(eu.capacity());
        for (int32_t i = 0; i < E; ++i) {
            int32_t u = (int32_t)sc.uint_(), v = (int32_t)sc.uint_(); W w = (W)sc.uint_();
            eu.push_back(u); ev.push_back(v); ew.push_back(w);
            if (!directed) { eu.push_back(v); ev.push_back(u); ew.push_back(w); }
        }
        sc.release();
        bidir_prepare(eu, ev, ew);
        tlog("csr");
        g_bidir = true;
        if (env_int("USE_ALT", 0)) {
            NLM = (int)env_int("NLM", NLM);
            LM_HUBS = (int)env_int("LM_HUBS", 0); BD_UB = (int)env_int("BD_UB", 0); BD_LB = (int)env_int("BD_LB", 0);
            if (BD_UB) g_alt_bidir = true;
            alt_prepare();
            aD.assign((size_t)V, INF); aH.assign((size_t)V, 0);
            g_alt = true;
            tlog("landmarks");
        }
        run_queries(argv[2], argv[3]);
        tlog("queries");
        return 0;
    }
    read_graph(argv[1]);
    tlog("read");
    contract();
    tlog("contract");
    if (env_int("STOP_AFTER_CONTRACT", 0)) return 0;
    build_labels();
    tlog("labels");
    run_queries(argv[2], argv[3]);
    tlog("queries");
    return 0;
}
