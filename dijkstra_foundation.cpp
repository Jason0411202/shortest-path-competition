// dijkstra_foundation.cpp
//
// Foundation Dijkstra implementation for the Shortest-Path Competition.
//
// Students must build their solver by modifying and extending this file.
// The authoritative description of the formats and the rules is README.md.
//
// Usage:
//     ./solver <graph_file> <query_file> <output_file>
//
// File formats (0-based vertex indices):
//
//   <graph_file>
//     V E FLAGS
//     u_1 v_1 w_1
//     ...
//     u_E v_E w_E
//     [ if FLAGS has bit 0 set:
//       x_0 y_0
//       ...
//       x_{V-1} y_{V-1}
//     ]
//
//   FLAGS is a bitmask:
//     bit 0 (value 1) — the coordinate block is present
//     bit 1 (value 2) — the edge list is DIRECTED: "u v w" is the arc u -> v
//                       only. When the bit is clear the graph is undirected and
//                       each line contributes both u -> v and v -> u.
//   A FLAGS value of 0 or 1 therefore means exactly what HAS_COORDS used to.
//
//   Coordinates are integers in the v2 datasets (they are still read as
//   floating point here, so both integer and decimal text parse correctly).
//
//   <query_file>
//     Q
//     s_1 t_1
//     ...
//     s_Q t_Q
//
//   <output_file>
//     d_1
//     ...
//     d_Q                 (-1 if t_i unreachable from s_i)

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <vector>

using std::vector;

struct Edge {
    int32_t to;
    int32_t w;
};

static const int32_t FLAG_COORDS   = 1;
static const int32_t FLAG_DIRECTED = 2;

static int32_t V = 0;
static int32_t E = 0;
static int32_t FLAGS = 0;
static vector<vector<Edge>> adj;
static vector<double> coord_x, coord_y;

static void read_graph(const char* path) {
    std::FILE* f = std::fopen(path, "r");
    if (!f) { std::fprintf(stderr, "cannot open graph file: %s\n", path); std::exit(1); }

    if (std::fscanf(f, "%d %d %d", &V, &E, &FLAGS) != 3) {
        std::fprintf(stderr, "bad graph header\n"); std::exit(1);
    }
    const bool directed = (FLAGS & FLAG_DIRECTED) != 0;
    adj.assign(V, {});
    for (int32_t i = 0; i < E; ++i) {
        int32_t u, v, w;
        if (std::fscanf(f, "%d %d %d", &u, &v, &w) != 3) {
            std::fprintf(stderr, "bad edge at line %d\n", i + 2); std::exit(1);
        }
        adj[u].push_back({v, w});
        if (!directed) adj[v].push_back({u, w});
    }
    if (FLAGS & FLAG_COORDS) {
        coord_x.resize(V);
        coord_y.resize(V);
        for (int32_t i = 0; i < V; ++i) {
            if (std::fscanf(f, "%lf %lf", &coord_x[i], &coord_y[i]) != 2) {
                std::fprintf(stderr, "bad coords at vertex %d\n", i); std::exit(1);
            }
        }
    }
    std::fclose(f);
}

// Standard Dijkstra with binary heap.
// Returns shortest-path distance from s to t, or -1 if unreachable.
// Early-terminates once t is settled.
static int64_t dijkstra(int32_t s, int32_t t) {
    if (s == t) return 0;

    static vector<int64_t> dist;
    if ((int32_t)dist.size() != V) dist.assign(V, -1);
    else std::fill(dist.begin(), dist.end(), (int64_t)-1);

    using PQItem = std::pair<int64_t, int32_t>; // (dist, node)
    std::priority_queue<PQItem, vector<PQItem>, std::greater<PQItem>> pq;

    dist[s] = 0;
    pq.push({0, s});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (d != dist[u]) continue;   // stale entry
        if (u == t) return d;
        for (const Edge& e : adj[u]) {
            int64_t nd = d + e.w;
            if (dist[e.to] == -1 || nd < dist[e.to]) {
                dist[e.to] = nd;
                pq.push({nd, e.to});
            }
        }
    }
    return -1;
}

static void run_queries(const char* qpath, const char* opath) {
    std::FILE* qf = std::fopen(qpath, "r");
    if (!qf) { std::fprintf(stderr, "cannot open query file: %s\n", qpath); std::exit(1); }
    std::FILE* of = std::fopen(opath, "w");
    if (!of) { std::fprintf(stderr, "cannot open output file: %s\n", opath); std::exit(1); }

    int32_t Q;
    if (std::fscanf(qf, "%d", &Q) != 1) {
        std::fprintf(stderr, "bad query header\n"); std::exit(1);
    }
    for (int32_t i = 0; i < Q; ++i) {
        int32_t s, t;
        if (std::fscanf(qf, "%d %d", &s, &t) != 2) {
            std::fprintf(stderr, "bad query %d\n", i); std::exit(1);
        }
        int64_t d = dijkstra(s, t);
        std::fprintf(of, "%lld\n", (long long)d);
    }
    std::fclose(qf);
    std::fclose(of);
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: %s <graph_file> <query_file> <output_file>\n", argv[0]);
        return 1;
    }
    read_graph(argv[1]);
    run_queries(argv[2], argv[3]);
    return 0;
}
