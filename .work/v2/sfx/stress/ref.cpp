// reference: plain Dijkstra per query (64-bit)
#include "../../common.inc"
int main(int argc, char** argv) {
    read_graph(argv[1]); read_queries(argv[2]);
    vector<vector<std::pair<int,uint64_t>>> g(V);
    for (int i = 0; i < E; ++i) g[eu[i]].push_back({ev[i], ew[i]});
    vector<uint64_t> d(V);
    for (int q = 0; q < Q; ++q) {
        std::fill(d.begin(), d.end(), ~0ull);
        Heap4<uint64_t> h; d[qs[q]] = 0; h.push(0, qs[q]);
        while (!h.empty()) { auto it = h.pop(); if (it.d > d[it.v]) continue; for (auto& p : g[it.v]) if (it.d + p.second < d[p.first]) { d[p.first] = it.d + p.second; h.push(d[p.first], p.first); } }
        qans[q] = d[qt[q]] == ~0ull ? -1 : (int64_t)d[qt[q]];
    }
    write_answers(argv[3]);
}
