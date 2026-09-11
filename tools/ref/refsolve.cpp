// refsolve.cpp -- INSTRUCTOR TOOL. Not a submission, not subject to the
// competition rules: it uses every core it can find, because the single-thread
// rule binds students and not the person producing the answer key.
//
// Modes
//   plan    <graph> <qplan> <out.queries> <out.answers>
//       Resolve a query plan into concrete (s,t) pairs and their distances.
//       An `R s rank` line is answered by a Dijkstra from s stopped at the
//       rank-th settled node: that node is t, and d(s,t) is known the moment it
//       settles, so query generation and answer generation are one pass.
//
//   answers <graph> <queries> <out.answers>
//       Distances for an existing query file. Output is byte-identical to what
//       dijkstra_foundation would produce, which is what makes it usable as the
//       shipped answer key.
//
//   verify  <graph> <queries> <answers> [samples]
//       Invariants that must hold before an instance may be released.
//
//   probe   <graph> <queries> [samples]
//       Dijkstra vs Euclidean-A* search space, plus the reachable fraction.
//       This is the regression guard for the v1 defect where edge weights were
//       exactly proportional to Euclidean length.
//
// Build: make -C tools/ref     (falls back to single-threaded without OpenMP)

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using std::vector;

static const int32_t FLAG_COORDS = 1;
static const int32_t FLAG_DIRECTED = 2;

// ---------------------------------------------------------------- input ----

// Whole-file slurp plus a hand-rolled integer scanner. The instructor runs this
// over ~600 MB of text; fscanf would dominate the wall clock for no reason.
struct Reader {
    vector<char> buf;
    size_t p = 0;

    explicit Reader(const char* path) {
        std::FILE* f = std::fopen(path, "rb");
        if (!f) { std::fprintf(stderr, "cannot open %s: %s\n", path, std::strerror(errno)); std::exit(1); }
        std::fseek(f, 0, SEEK_END);
        long n = std::ftell(f);
        std::fseek(f, 0, SEEK_SET);
        buf.resize((size_t)n + 1);
        if (n > 0 && std::fread(buf.data(), 1, (size_t)n, f) != (size_t)n) {
            std::fprintf(stderr, "short read on %s\n", path); std::exit(1);
        }
        buf[(size_t)n] = '\0';
        std::fclose(f);
    }

    bool eof() {
        while (p < buf.size() - 1 && (unsigned char)buf[p] <= ' ') ++p;
        return p >= buf.size() - 1;
    }

    int64_t integer() {
        while (p < buf.size() && (unsigned char)buf[p] <= ' ') ++p;
        bool neg = false;
        if (buf[p] == '-') { neg = true; ++p; }
        int64_t v = 0;
        while (p < buf.size() && buf[p] >= '0' && buf[p] <= '9') {
            v = v * 10 + (buf[p] - '0');
            ++p;
        }
        return neg ? -v : v;
    }

    char token() {
        while (p < buf.size() && (unsigned char)buf[p] <= ' ') ++p;
        return buf[p++];
    }

    // v2 instances write integer coordinates; the v1 files this tool is also
    // pointed at for comparison wrote "%.6f", so accept a fractional part.
    double number() {
        while (p < buf.size() && (unsigned char)buf[p] <= ' ') ++p;
        bool neg = false;
        if (buf[p] == '-') { neg = true; ++p; }
        double v = 0;
        while (p < buf.size() && buf[p] >= '0' && buf[p] <= '9')
            v = v * 10 + (buf[p++] - '0');
        if (p < buf.size() && buf[p] == '.') {
            ++p;
            double scale = 0.1;
            while (p < buf.size() && buf[p] >= '0' && buf[p] <= '9') {
                v += (buf[p++] - '0') * scale;
                scale *= 0.1;
            }
        }
        return neg ? -v : v;
    }
};

struct Graph {
    int32_t V = 0, E = 0, flags = 0;
    bool directed = false, has_coords = false;
    vector<int64_t> off;   // CSR over arcs in the search direction
    vector<int32_t> to;
    vector<int64_t> w;
    vector<double> cx, cy;

    void load(const char* path) {
        Reader r(path);
        V = (int32_t)r.integer();
        E = (int32_t)r.integer();
        flags = (int32_t)r.integer();
        directed = (flags & FLAG_DIRECTED) != 0;
        has_coords = (flags & FLAG_COORDS) != 0;

        vector<int32_t> eu((size_t)E), ev((size_t)E);
        vector<int64_t> ew((size_t)E);
        vector<int64_t> deg((size_t)V + 1, 0);
        for (int32_t i = 0; i < E; ++i) {
            eu[i] = (int32_t)r.integer();
            ev[i] = (int32_t)r.integer();
            ew[i] = r.integer();
            deg[eu[i]]++;
            if (!directed) deg[ev[i]]++;
        }
        off.assign((size_t)V + 1, 0);
        int64_t s = 0;
        for (int32_t i = 0; i < V; ++i) { off[i] = s; s += deg[i]; }
        off[V] = s;
        to.resize((size_t)s);
        w.resize((size_t)s);
        vector<int64_t> pos(off.begin(), off.end() - 1);
        for (int32_t i = 0; i < E; ++i) {
            to[pos[eu[i]]] = ev[i]; w[pos[eu[i]]] = ew[i]; pos[eu[i]]++;
            if (!directed) { to[pos[ev[i]]] = eu[i]; w[pos[ev[i]]] = ew[i]; pos[ev[i]]++; }
        }
        if (has_coords) {
            cx.resize((size_t)V); cy.resize((size_t)V);
            for (int32_t i = 0; i < V; ++i) {
                cx[i] = r.number();
                cy[i] = r.number();
            }
        }
    }
};

// ------------------------------------------------------------- searching ----

// Per-thread scratch, sized so that answering millions of queries in parallel
// stays within memory: threads x V x 12 bytes.
struct Scratch {
    vector<int64_t> dist;
    vector<uint32_t> stamp;
    uint32_t gen = 0;
    vector<std::pair<int64_t, int32_t>> heap;

    void init(int32_t V) {
        dist.assign((size_t)V, 0);
        stamp.assign((size_t)V, 0);
        gen = 0;
    }
    void begin() {
        if (++gen == 0) { std::fill(stamp.begin(), stamp.end(), 0); gen = 1; }
        heap.clear();
    }
    inline bool seen(int32_t v) const { return stamp[v] == gen; }
    inline void set(int32_t v, int64_t d) { stamp[v] = gen; dist[v] = d; }
};

struct HeapGreater {
    bool operator()(const std::pair<int64_t, int32_t>& a,
                    const std::pair<int64_t, int32_t>& b) const {
        return a.first > b.first;
    }
};

// Dijkstra from s. Stops when `target` settles (if >= 0) or when `stop_rank`
// nodes have settled (if > 0). Returns the settled count; `out_node` / `out_d`
// carry the node that triggered the stop.
static int64_t dijkstra(const Graph& g, Scratch& sc, int32_t s, int32_t target,
                        int64_t stop_rank, int32_t* out_node, int64_t* out_d) {
    sc.begin();
    sc.set(s, 0);
    sc.heap.assign(1, {0, s});
    int64_t settled = 0;
    int32_t last = s;
    int64_t last_d = 0;
    HeapGreater cmp;

    while (!sc.heap.empty()) {
        std::pop_heap(sc.heap.begin(), sc.heap.end(), cmp);
        auto [d, u] = sc.heap.back();
        sc.heap.pop_back();
        if (!sc.seen(u) || d != sc.dist[u]) continue;   // stale entry
        ++settled;
        last = u; last_d = d;
        if (u == target) { *out_node = u; *out_d = d; return settled; }
        if (stop_rank > 0 && settled >= stop_rank) { *out_node = u; *out_d = d; return settled; }
        for (int64_t j = g.off[u]; j < g.off[u + 1]; ++j) {
            int32_t v = g.to[j];
            int64_t nd = d + g.w[j];
            if (!sc.seen(v) || nd < sc.dist[v]) {
                sc.set(v, nd);
                sc.heap.push_back({nd, v});
                std::push_heap(sc.heap.begin(), sc.heap.end(), cmp);
            }
        }
    }
    *out_node = last;
    *out_d = (target >= 0) ? -1 : last_d;   // exhausted without reaching target
    return settled;
}

static int64_t astar(const Graph& g, Scratch& sc, int32_t s, int32_t t,
                     double lambda) {
    sc.begin();
    sc.set(s, 0);
    HeapGreater cmp;
    auto h = [&](int32_t u) {
        double dx = g.cx[u] - g.cx[t], dy = g.cy[u] - g.cy[t];
        return (int64_t)(lambda * std::sqrt(dx * dx + dy * dy));
    };
    sc.heap.assign(1, {h(s), s});
    int64_t settled = 0;
    while (!sc.heap.empty()) {
        std::pop_heap(sc.heap.begin(), sc.heap.end(), cmp);
        auto [f, u] = sc.heap.back();
        sc.heap.pop_back();
        int64_t d = sc.dist[u];
        if (!sc.seen(u) || f != d + h(u)) continue;
        ++settled;
        if (u == t) return settled;
        for (int64_t j = g.off[u]; j < g.off[u + 1]; ++j) {
            int32_t v = g.to[j];
            int64_t nd = d + g.w[j];
            if (!sc.seen(v) || nd < sc.dist[v]) {
                sc.set(v, nd);
                sc.heap.push_back({nd + h(v), v});
                std::push_heap(sc.heap.begin(), sc.heap.end(), cmp);
            }
        }
    }
    return settled;
}

// ----------------------------------------------------------------- modes ----

static int n_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

static void write_answers(const char* path, const vector<int64_t>& ans) {
    std::FILE* f = std::fopen(path, "w");
    if (!f) { std::fprintf(stderr, "cannot write %s\n", path); std::exit(1); }
    std::string out;
    out.reserve(ans.size() * 8);
    char tmp[32];
    for (int64_t d : ans) {
        int n = std::snprintf(tmp, sizeof tmp, "%lld\n", (long long)d);
        out.append(tmp, (size_t)n);
    }
    std::fwrite(out.data(), 1, out.size(), f);
    std::fclose(f);
}

static int mode_plan(int argc, char** argv) {
    if (argc != 6) { std::fprintf(stderr, "plan <graph> <qplan> <out.queries> <out.answers>\n"); return 2; }
    Graph g; g.load(argv[2]);
    Reader r(argv[3]);
    int64_t Q = r.integer();
    vector<char> kind((size_t)Q);
    vector<int32_t> src((size_t)Q);
    vector<int64_t> arg((size_t)Q);
    for (int64_t i = 0; i < Q; ++i) {
        kind[i] = r.token();
        src[i] = (int32_t)r.integer();
        arg[i] = r.integer();
    }

    vector<int32_t> tgt((size_t)Q);
    vector<int64_t> ans((size_t)Q);
    int64_t short_ranks = 0, unreachable = 0;

    int T = n_threads();
    vector<Scratch> scratch((size_t)T);
    for (auto& s : scratch) s.init(g.V);
    std::fprintf(stderr, "plan: V=%d arcs=%lld Q=%lld threads=%d\n",
                 g.V, (long long)g.off[g.V], (long long)Q, T);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64) reduction(+ : short_ranks, unreachable)
#endif
    for (int64_t i = 0; i < Q; ++i) {
#ifdef _OPENMP
        Scratch& sc = scratch[(size_t)omp_get_thread_num()];
#else
        Scratch& sc = scratch[0];
#endif
        int32_t node; int64_t d;
        if (kind[i] == 'R') {
            int64_t settled = dijkstra(g, sc, src[i], -1, arg[i], &node, &d);
            if (settled < arg[i]) short_ranks += 1;   // reachable set smaller than the rank
            tgt[i] = node; ans[i] = d;
        } else {
            int32_t t = (int32_t)arg[i];
            if (t == src[i]) { tgt[i] = t; ans[i] = 0; continue; }
            dijkstra(g, sc, src[i], t, 0, &node, &d);
            tgt[i] = t; ans[i] = d;
            if (d < 0) unreachable += 1;
        }
    }

    {
        std::FILE* f = std::fopen(argv[4], "w");
        if (!f) { std::fprintf(stderr, "cannot write %s\n", argv[4]); return 1; }
        std::string out;
        out.reserve((size_t)Q * 16);
        char tmp[48];
        int n = std::snprintf(tmp, sizeof tmp, "%lld\n", (long long)Q);
        out.append(tmp, (size_t)n);
        for (int64_t i = 0; i < Q; ++i) {
            n = std::snprintf(tmp, sizeof tmp, "%d %d\n", src[i], tgt[i]);
            out.append(tmp, (size_t)n);
        }
        std::fwrite(out.data(), 1, out.size(), f);
        std::fclose(f);
    }
    write_answers(argv[5], ans);

    int64_t maxd = 0;
    for (int64_t d : ans) maxd = std::max(maxd, d);
    std::printf("{\"Q\": %lld, \"unreachable\": %lld, \"rank_truncated\": %lld, "
                "\"max_finite_distance\": %lld, \"fits_uint32\": %s}\n",
                (long long)Q, (long long)unreachable, (long long)short_ranks,
                (long long)maxd, maxd <= 4294967295LL ? "true" : "false");
    return 0;
}

static int mode_answers(int argc, char** argv) {
    if (argc != 5) { std::fprintf(stderr, "answers <graph> <queries> <out.answers>\n"); return 2; }
    Graph g; g.load(argv[2]);
    Reader r(argv[3]);
    int64_t Q = r.integer();
    vector<int32_t> s((size_t)Q), t((size_t)Q);
    for (int64_t i = 0; i < Q; ++i) { s[i] = (int32_t)r.integer(); t[i] = (int32_t)r.integer(); }

    vector<int64_t> ans((size_t)Q);
    int T = n_threads();
    vector<Scratch> scratch((size_t)T);
    for (auto& x : scratch) x.init(g.V);
    std::fprintf(stderr, "answers: V=%d Q=%lld threads=%d\n", g.V, (long long)Q, T);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (int64_t i = 0; i < Q; ++i) {
#ifdef _OPENMP
        Scratch& sc = scratch[(size_t)omp_get_thread_num()];
#else
        Scratch& sc = scratch[0];
#endif
        if (s[i] == t[i]) { ans[i] = 0; continue; }
        int32_t node; int64_t d;
        dijkstra(g, sc, s[i], t[i], 0, &node, &d);
        ans[i] = d;
    }
    write_answers(argv[4], ans);
    return 0;
}

static int mode_verify(int argc, char** argv) {
    if (argc < 5) { std::fprintf(stderr, "verify <graph> <queries> <answers> [samples]\n"); return 2; }
    int samples = argc > 5 ? std::atoi(argv[5]) : 10000;
    Graph g; g.load(argv[2]);
    Reader rq(argv[3]);
    int64_t Q = rq.integer();
    vector<int32_t> s((size_t)Q), t((size_t)Q);
    for (int64_t i = 0; i < Q; ++i) { s[i] = (int32_t)rq.integer(); t[i] = (int32_t)rq.integer(); }
    Reader ra(argv[4]);
    vector<int64_t> ans((size_t)Q);
    for (int64_t i = 0; i < Q; ++i) ans[i] = ra.integer();

    int failures = 0;
    int T = n_threads();
    vector<Scratch> scratch((size_t)T);
    for (auto& x : scratch) x.init(g.V);

    // 1. recompute a sample of the answers independently
    int64_t step = std::max<int64_t>(1, Q / std::max(1, samples));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : failures)
#endif
    for (int64_t i = 0; i < Q; i += step) {
#ifdef _OPENMP
        Scratch& sc = scratch[(size_t)omp_get_thread_num()];
#else
        Scratch& sc = scratch[0];
#endif
        int64_t d;
        int32_t node;
        if (s[i] == t[i]) d = 0;
        else dijkstra(g, sc, s[i], t[i], 0, &node, &d);
        if (d != ans[i]) {
            std::fprintf(stderr, "MISMATCH q%lld (%d->%d): answers=%lld recomputed=%lld\n",
                         (long long)i, s[i], t[i], (long long)ans[i], (long long)d);
            failures += 1;
        }
    }

    // 2. symmetry, on undirected instances only
    if (!g.directed) {
        Scratch& sc = scratch[0];
        for (int64_t i = 0; i < Q; i += step * 4) {
            if (s[i] == t[i] || ans[i] < 0) continue;
            int32_t node; int64_t d;
            dijkstra(g, sc, t[i], s[i], 0, &node, &d);
            if (d != ans[i]) {
                std::fprintf(stderr, "ASYMMETRY %d<->%d: %lld vs %lld\n",
                             s[i], t[i], (long long)ans[i], (long long)d);
                failures++;
            }
        }
    }

    // 3. every reported -1 really is unreachable, by an independent full sweep
    int64_t checked_unreach = 0;
    {
        Scratch& sc = scratch[0];
        for (int64_t i = 0; i < Q && checked_unreach < 200; ++i) {
            if (ans[i] != -1) continue;
            checked_unreach++;
            sc.begin();
            sc.set(s[i], 0);
            vector<int32_t> stack{s[i]};
            bool found = false;
            while (!stack.empty() && !found) {
                int32_t u = stack.back(); stack.pop_back();
                for (int64_t j = g.off[u]; j < g.off[u + 1]; ++j) {
                    int32_t v = g.to[j];
                    if (!sc.seen(v)) { sc.set(v, 0); stack.push_back(v); if (v == t[i]) { found = true; break; } }
                }
            }
            if (found) {
                std::fprintf(stderr, "BAD -1: %d -> %d is actually reachable\n", s[i], t[i]);
                failures++;
            }
        }
    }

    int64_t maxd = 0, n_unreach = 0;
    for (int64_t d : ans) { if (d < 0) n_unreach++; else maxd = std::max(maxd, d); }
    std::printf("{\"sampled\": %lld, \"failures\": %d, \"unreachable\": %lld, "
                "\"unreachable_frac\": %.6f, \"unreachable_confirmed\": %lld, "
                "\"max_finite_distance\": %lld, \"fits_uint32\": %s}\n",
                (long long)((Q + step - 1) / step), failures, (long long)n_unreach,
                (double)n_unreach / (double)Q, (long long)checked_unreach,
                (long long)maxd, maxd <= 4294967295LL ? "true" : "false");
    return failures ? 1 : 0;
}

static int mode_probe(int argc, char** argv) {
    if (argc < 4) { std::fprintf(stderr, "probe <graph> <queries> [samples]\n"); return 2; }
    int samples = argc > 4 ? std::atoi(argv[4]) : 200;
    Graph g; g.load(argv[2]);
    if (!g.has_coords) { std::printf("{\"astar_reduction\": null, \"reason\": \"no coords\"}\n"); return 0; }

    // the tightest admissible scale for this instance
    double lambda = 1e300;
    for (int32_t u = 0; u < g.V; ++u) {
        for (int64_t j = g.off[u]; j < g.off[u + 1]; ++j) {
            int32_t v = g.to[j];
            double dx = g.cx[u] - g.cx[v], dy = g.cy[u] - g.cy[v];
            double d = std::sqrt(dx * dx + dy * dy);
            if (d > 0) lambda = std::min(lambda, (double)g.w[j] / d);
        }
    }

    Reader rq(argv[3]);
    int64_t Q = rq.integer();
    vector<int32_t> s((size_t)Q), t((size_t)Q);
    for (int64_t i = 0; i < Q; ++i) { s[i] = (int32_t)rq.integer(); t[i] = (int32_t)rq.integer(); }

    int64_t step = std::max<int64_t>(1, Q / samples);
    double sum_d = 0, sum_a = 0;
    int64_t n = 0;
    Scratch sc; sc.init(g.V);
    for (int64_t i = 0; i < Q; i += step) {
        if (s[i] == t[i]) continue;
        int32_t node; int64_t d;
        int64_t sd = dijkstra(g, sc, s[i], t[i], 0, &node, &d);
        if (d < 0) continue;
        int64_t sa = astar(g, sc, s[i], t[i], lambda);
        sum_d += (double)sd; sum_a += (double)sa; n++;
    }
    std::printf("{\"lambda_admissible\": %.6f, \"samples\": %lld, "
                "\"mean_settled_dijkstra\": %.1f, \"mean_settled_astar\": %.1f, "
                "\"astar_reduction\": %.3f}\n",
                lambda, (long long)n, sum_d / (double)n, sum_a / (double)n,
                sum_d / sum_a);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
            "usage: %s <plan|answers|verify|probe> ...\n"
            "  plan    <graph> <qplan> <out.queries> <out.answers>\n"
            "  answers <graph> <queries> <out.answers>\n"
            "  verify  <graph> <queries> <answers> [samples]\n"
            "  probe   <graph> <queries> [samples]\n", argv[0]);
        return 2;
    }
    std::string m = argv[1];
    if (m == "plan")    return mode_plan(argc, argv);
    if (m == "answers") return mode_answers(argc, argv);
    if (m == "verify")  return mode_verify(argc, argv);
    if (m == "probe")   return mode_probe(argc, argv);
    std::fprintf(stderr, "unknown mode: %s\n", argv[1]);
    return 2;
}
