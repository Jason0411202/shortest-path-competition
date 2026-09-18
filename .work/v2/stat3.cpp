// lattice3d structure stats: how many edges are useless (not shortest u-v path) and
// how many are caught by 3-hop square detours.
#include "common.inc"
int main(int argc, char** argv) {
    g_debug = 1; g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]); tlog("read");
    // detect torus lattice
    int d = 0; int64_t n = 0;
    for (int dd = 2; dd <= 3; ++dd) { int64_t s = llround(pow((double)V, 1.0 / dd)); for (int64_t c = s - 1; c <= s + 1; ++c) { int64_t p = 1; for (int k = 0; k < dd; ++k) p *= c; if (p == V && (int64_t)E == dd * (int64_t)V) { d = dd; n = c; } } }
    std::fprintf(stderr, "d=%d n=%lld\n", d, (long long)n);
    int64_t st[3] = {1, n, n * n};
    bool ok = true;
    for (int64_t i = 0; i < V && ok; ++i) for (int k = 0; k < d; ++k) {
        int64_t e = i * d + k; int64_t c = (i / st[k]) % n;
        int64_t j = c + 1 < n ? i + st[k] : i - c * st[k];
        if (eu[e] != i || ev[e] != j) { ok = false; break; }
    }
    std::fprintf(stderr, "lattice order ok=%d\n", ok);
    tlog("detect");
    // weights per axis
    auto W = [&](int64_t i, int k) -> uint32_t { return ew[i * d + k]; };
    auto nb = [&](int64_t i, int k, int dir) -> int64_t { int64_t c = (i / st[k]) % n; if (dir > 0) return c + 1 < n ? i + st[k] : i - c * st[k]; return c > 0 ? i - st[k] : i + (n - 1) * st[k]; };
    // weight of edge from i in direction (k,dir)
    auto wd = [&](int64_t i, int k, int dir) -> uint32_t { return dir > 0 ? W(i, k) : W(nb(i, k, -1), k); };
    int64_t cut3 = 0;
    for (int64_t i = 0; i < V; ++i) for (int k = 0; k < d; ++k) {
        int64_t j = nb(i, k, 1); uint64_t w = W(i, k);
        bool dom = false;
        for (int k2 = 0; k2 < d && !dom; ++k2) if (k2 != k) for (int s2 = -1; s2 <= 1 && !dom; s2 += 2) {
            int64_t a = nb(i, k2, s2), b = nb(j, k2, s2);
            uint64_t alt = (uint64_t)wd(i, k2, s2) + wd(a, k, 1) + wd(j, k2, s2);
            (void)b;
            if (alt < w) dom = true;
        }
        cut3 += dom;
    }
    std::fprintf(stderr, "3-hop dominated edges: %lld of %d\n", (long long)cut3, E);
    tlog("3hop");
}
