#include "common.inc"
#include <cmath>
// Pragmas go after all standard headers: a target pragma in front of them
// breaks a build without -march flags (always_inline target mismatch).
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3")
#endif
#if defined(__GNUC__) && !defined(__clang__) && (defined(__x86_64__) || defined(__i386__))
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")
#endif

// ======================================================================= road
// Directed graphs with coordinates (road2d / hugeq families).
//
// On a row-major S x S lattice whose rows/columns carry road classes, the
// fast lines (arterials, highways) cut the lattice into small rectangular
// "cells" of local streets.  Every path leaving a cell crosses its boundary.
//   1. Cells: exact dense min-plus elimination of the interior gives, for
//      every interior vertex, its first-exit distances to the boundary
//      (and first-entry distances from it) and the boundary-to-boundary
//      distances through the interior (added as shortcuts when shorter than
//      any path along the boundary).
//   2. The remaining "network" graph (boundary lines + shortcuts) gets a
//      contraction hierarchy (lazy priority, bounded witness searches) and
//      hub labels built top-down.
//   3. Interior labels = merge of (exit distance + label of exit vertex).
//   4. Queries: label intersection by scatter/gather; a pair inside one cell
//      additionally gets a Dijkstra restricted to that cell's interior.
// Without lattice structure everything is "network": plain CH + labels.
namespace road {

static int env_int(const char* name, int def) {
    const char* e = std::getenv(name);
    return e ? std::atoi(e) : def;
}
static double now_s() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - g_t0).count();
}

// Detects a row-major S x S lattice: V = S*S and (almost) all arcs join
// lattice neighbours (dx, dy in {-1,0,1}).  Returns S or 0.
static int detect_lattice() {
    int64_t S = 1;
    while (S * S < V) ++S;
    if (S * S != V || S < 8) return 0;
    int64_t bad = 0;
    for (int32_t i = 0; i < E; ++i) {
        int32_t a = eu[i], b = ev[i];
        int32_t ax = a % S, ay = a / S, bx = b % S, by = b / S;
        if (std::abs(ax - bx) > 1 || std::abs(ay - by) > 1) ++bad;
    }
    return bad * 50 <= E ? (int)S : 0;
}

// Road class per lattice line from the travel-time / length ratio of the
// arcs running along it: 0 local, 1 arterial, 2 highway.
static void line_classes(int S, vector<int>& rowc, vector<int>& colc) {
    vector<double> rs((size_t)S, 0), cs((size_t)S, 0);
    vector<int> rn((size_t)S, 0), cn((size_t)S, 0);
    for (int32_t i = 0; i < E; ++i) {
        int32_t a = eu[i], b = ev[i];
        double dx = (double)(cx[a] - cx[b]), dy = (double)(cy[a] - cy[b]);
        double len = std::sqrt(dx * dx + dy * dy);
        if (!(len > 0)) continue;
        double r = (double)ew[i] / len;
        if ((b == a + 1 || a == b + 1) && a / S == b / S) {
            rs[a / S] += r; rn[a / S]++;
        } else if (b == a + S || a == b + S) {
            cs[a % S] += r; cn[a % S]++;
        }
    }
    auto classify = [&](const vector<double>& s, const vector<int>& n, vector<int>& out) {
        vector<double> m;
        for (int i = 0; i < S; ++i) if (n[i]) m.push_back(s[i] / n[i]);
        out.assign((size_t)S, 1);
        if (m.empty()) return;
        std::nth_element(m.begin(), m.begin() + m.size() / 2, m.end());
        double med = m[m.size() / 2];
        for (int i = 0; i < S; ++i) {
            if (!n[i]) continue;
            double q = (s[i] / n[i]) / med;
            out[i] = q < 0.2 ? 2 : q < 0.55 ? 1 : 0;
        }
    };
    classify(rs, rn, rowc);
    classify(cs, cn, colc);
}

struct Plan {
    int S = 0;
    vector<int32_t> cellv;    // cell interiors, each in elimination order
    vector<uint32_t> cellst;  // cell start offsets (size = #cells + 1)
    vector<int32_t> cellx;    // 4 extra boundary vertices per cell (rectangle corners; -1 = none)
    vector<int32_t> cell_of;  // vertex -> cell index, -1 for network vertices
};

// 1-D nested dissection order of [lo, hi] (middle last)
static void nd1(int lo, int hi, vector<int>& out) {
    if (lo > hi) return;
    int m = (lo + hi) / 2;
    nd1(lo, m - 1, out); nd1(m + 1, hi, out);
    out.push_back(m);
}
// 2-D nested dissection order of the lattice rectangle [x0,x1] x [y0,y1]
static void nd2(int S, int x0, int x1, int y0, int y1, vector<int32_t>& out) {
    if (x0 > x1 || y0 > y1) return;
    int w = x1 - x0 + 1, h = y1 - y0 + 1;
    if (w * h <= 1) { out.push_back((int32_t)(y0 * S + x0)); return; }
    vector<int> line;
    if (w >= h) {
        int m = (x0 + x1) / 2;
        nd2(S, x0, m - 1, y0, y1, out); nd2(S, m + 1, x1, y0, y1, out);
        nd1(y0, y1, line);
        for (int y : line) out.push_back((int32_t)(y * S + m));
    } else {
        int m = (y0 + y1) / 2;
        nd2(S, x0, x1, y0, m - 1, out); nd2(S, x0, x1, m + 1, y1, out);
        nd1(x0, x1, line);
        for (int x : line) out.push_back((int32_t)(m * S + x));
    }
}

static Plan make_plan() {
    Plan p;
    if (env_int("NOCELLS", 0)) return p;
    int S = detect_lattice();
    if (!S) return p;
    p.S = S;
    vector<int> rowc, colc;
    line_classes(S, rowc, colc);
    if (g_debug) {
        int na = 0, nh = 0;
        for (int i = 0; i < S; ++i) { na += rowc[i] == 1; nh += rowc[i] == 2; }
        std::fprintf(stderr, "lattice S=%d: arterial rows %d, highway rows %d\n", S, na, nh);
    }
    const int max_interior = env_int("CELL_MAXI", 400);
    auto runs = [&](const vector<int>& c) {
        vector<std::pair<int, int>> r;
        for (int i = 0; i < S;) {
            if (c[i] != 0) { ++i; continue; }
            int j = i;
            while (j + 1 < S && c[j + 1] == 0) ++j;
            r.push_back({i, j});
            i = j + 1;
        }
        return r;
    };
    auto xr = runs(colc), yr = runs(rowc);
    vector<int32_t> cellv;
    vector<uint32_t> cellst;
    vector<int32_t> cellx;
    for (auto& yy : yr)
        for (auto& xx : xr) {
            int w = xx.second - xx.first + 1, h = yy.second - yy.first + 1;
            if (w * h > max_interior) continue;
            cellst.push_back((uint32_t)cellv.size());
            nd2(S, xx.first, xx.second, yy.first, yy.second, cellv);
            int cxs[2] = {xx.first - 1, xx.second + 1}, cys[2] = {yy.first - 1, yy.second + 1};
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b) {
                    int X = cxs[a], Y = cys[b];
                    cellx.push_back(X >= 0 && X < S && Y >= 0 && Y < S ? (int32_t)(Y * S + X) : -1);
                }
        }
    cellst.push_back((uint32_t)cellv.size());
    size_t nc = cellst.size() - 1;
    vector<int32_t> cell_of((size_t)V, -1);
    for (size_t c = 0; c < nc; ++c)
        for (uint32_t k = cellst[c]; k < cellst[c + 1]; ++k) cell_of[cellv[k]] = (int32_t)c;
    // Cells must be separated by network vertices: drop cells joined by an arc.
    vector<uint8_t> drop(nc, 0);
    for (int32_t i = 0; i < E; ++i) {
        int32_t a = cell_of[eu[i]], b = cell_of[ev[i]];
        if (a >= 0 && b >= 0 && a != b) { drop[a] = 1; drop[b] = 1; }
    }
    // Boundary size guard (a few stray arcs could make a boundary large).
    {
        vector<uint32_t> cnt(nc, 0);
        for (int32_t i = 0; i < E; ++i) {
            int32_t a = cell_of[eu[i]], b = cell_of[ev[i]];
            if (a >= 0 && b < 0) ++cnt[a];
            if (b >= 0 && a < 0) ++cnt[b];
        }
        for (size_t c = 0; c < nc; ++c) if (cnt[c] > 1200) drop[c] = 1;
    }
    for (size_t c = 0; c < nc; ++c) {
        if (drop[c]) {
            for (uint32_t k = cellst[c]; k < cellst[c + 1]; ++k) cell_of[cellv[k]] = -1;
            continue;
        }
        int32_t id = (int32_t)(p.cellst.size());
        p.cellst.push_back((uint32_t)p.cellv.size());
        for (uint32_t k = cellst[c]; k < cellst[c + 1]; ++k) { p.cellv.push_back(cellv[k]); cell_of[cellv[k]] = id; }
        for (int k = 0; k < 4; ++k) p.cellx.push_back(cellx[4 * c + k]);
    }
    p.cellst.push_back((uint32_t)p.cellv.size());
    if (p.cellst.size() == 1) p.cellst.clear();
    p.cell_of.swap(cell_of);
    return p;
}

template <class DT>
struct Engine {
    static constexpr DT INF = ~DT(0);
    static constexpr DT CAP = (DT)(~DT(0) >> 1);
    bool overflow = false;

    static inline DT sat_add(DT a, DT b) {   // a <= CAP; b may be INF
        DT c = a + b;
        return c < b ? INF : c;
    }

    struct Arc { int32_t to; int32_t hops; DT fw, bw; };   // fw: v->to, bw: to->v
    vector<Arc> pool;                                      // pooled adjacency
    vector<uint32_t> aoff, asz, acap;

    // --- witness search state
    vector<DT> wd;
    vector<uint32_t> wst;
    uint32_t wcur = 0;
    vector<uint32_t> tst;
    uint32_t tcur = 0;
    vector<DT> md;
    vector<uint32_t> mst;
    uint32_t mcur = 0;
    Heap4<DT> hp;
    int con_settle = 300;
    uint64_t n_settle = 0, n_search = 0;
    uint64_t work = 0, work_budget = ~0ull;   // CH effort guard (deterministic)
    bool aborted = false;
    vector<int32_t> lvl;
    int64_t level_coef = 1000;

    // --- network CH (ranks 0..NR-1 over network vertices)
    int32_t NR = 0;
    vector<int32_t> rank_of, order;
    vector<uint32_t> uhead;
    struct UArc { int32_t to; DT fw, bw; };
    vector<UArc> ua;
    vector<UArc> tmp_up;
    vector<uint32_t> tmp_head;

    // --- cells: exits (interior -> boundary) and entries (boundary -> interior)
    vector<uint32_t> exs, ens;            // per vertex start (interior only)
    vector<uint32_t> exn, enn;            // per vertex count
    struct XE { int32_t b; DT d; };
    vector<XE> exl, enl;
    size_t exl_n = 0, enl_n = 0;
    vector<uint8_t> done;
    vector<int32_t> iord;                 // interior vertices, cell by cell
    vector<uint32_t> iost;                // cell starts in iord
    const vector<int32_t>* cell_of = nullptr;

    inline DT addc(DT a, DT b) {
        uint64_t s = (uint64_t)a + (uint64_t)b;
        if (s > (uint64_t)CAP) { overflow = true; return CAP; }
        return (DT)s;
    }
    inline Arc* adj(int32_t v) { return pool.data() + aoff[v]; }

    void push_arc(int32_t u, const Arc& a) {
        if (asz[u] == acap[u]) {
            uint32_t nc = acap[u] ? acap[u] * 2 : 4;
            size_t no = pool.size();
            pool.resize(pool.size() + nc);
            std::memcpy(pool.data() + no, pool.data() + aoff[u], sizeof(Arc) * asz[u]);
            aoff[u] = (uint32_t)no; acap[u] = nc;
        }
        pool[aoff[u] + asz[u]++] = a;
    }
    void add_arc(int32_t u, int32_t x, DT fw, DT bw, int32_t hops) {
        Arc* a = adj(u);
        uint32_t n = asz[u];
        for (uint32_t k = 0; k < n; ++k) {
            if (a[k].to == x) {
                if (fw < a[k].fw) a[k].fw = fw;
                if (bw < a[k].bw) a[k].bw = bw;
                if (hops > a[k].hops) a[k].hops = hops;
                return;
            }
        }
        push_arc(u, {x, hops, fw, bw});
    }

    void build_graph() {
        vector<uint32_t> deg((size_t)V, 0);
        for (int32_t i = 0; i < E; ++i) if (eu[i] != ev[i]) { ++deg[eu[i]]; ++deg[ev[i]]; }
        aoff.assign((size_t)V, 0); asz.assign((size_t)V, 0); acap.assign((size_t)V, 0);
        uint64_t tot = 0;
        for (int32_t v = 0; v < V; ++v) { aoff[v] = (uint32_t)tot; acap[v] = deg[v] + 2; tot += acap[v]; }
        pool.reserve(tot * 2);
        pool.resize(tot);
        for (int32_t i = 0; i < E; ++i) {
            int32_t a = eu[i], b = ev[i];
            if (a == b) continue;
            DT w = (DT)ew[i];
            if (w > CAP) { overflow = true; w = CAP; }
            add_arc(a, b, w, INF, 1);
            add_arc(b, a, INF, w, 1);
        }
        wd.assign((size_t)V, 0); wst.assign((size_t)V, 0);
        tst.assign((size_t)V, 0);
        md.assign((size_t)V, 0); mst.assign((size_t)V, 0);
        lvl.assign((size_t)V, 0);
        done.assign((size_t)V, 0);
    }

    inline void new_search() { if (++wcur == 0) { std::fill(wst.begin(), wst.end(), 0); wcur = 1; } }
    inline void new_targets() { if (++tcur == 0) { std::fill(tst.begin(), tst.end(), 0); tcur = 1; } }
    inline void new_marks() { if (++mcur == 0) { std::fill(mst.begin(), mst.end(), 0); mcur = 1; } }

    // Dijkstra from src over the remaining graph, bounded by limit; stops once
    // `pending` targets (tst == tcur) are settled or con_settle vertices were.
    void wsearch(int32_t src, DT limit, int pending) {
        ++n_search;
        new_search();
        hp.clear();
        wd[src] = 0; wst[src] = wcur;
        hp.push(0, src);
        int settled = 0;
        while (!hp.empty()) {
            auto it = hp.pop();
            int32_t v = it.v; DT d = it.d;
            if (d > wd[v]) continue;
            if (d > limit) break;
            if (tst[v] == tcur) { if (--pending <= 0) break; }
            if (++settled > con_settle) break;
            const Arc* a = adj(v);
            uint32_t n = asz[v];
            for (uint32_t k = 0; k < n; ++k) {
                if (a[k].fw == INF) continue;
                DT nd = d + a[k].fw;
                if (nd > limit) continue;
                int32_t x = a[k].to;
                if (wst[x] != wcur || nd < wd[x]) {
                    wst[x] = wcur; wd[x] = nd;
                    hp.push(nd, x);
                }
            }
        }
        n_settle += (uint64_t)settled;
    }

    inline void mark_out(int32_t u, int32_t skip) {
        new_marks();
        const Arc* b = adj(u);
        uint32_t n = asz[u];
        work += n;
        for (uint32_t k = 0; k < n; ++k)
            if (b[k].fw != INF && b[k].to != skip) { mst[b[k].to] = mcur; md[b[k].to] = b[k].fw; }
    }
    // u->x path of <= 2 arcs (u's out-arcs marked), avoiding skip, cost <= c?
    inline bool two_hop(int32_t x, DT c, int32_t skip) {
        if (mst[x] == mcur && md[x] <= c) return true;
        const Arc* b = adj(x);
        uint32_t n = asz[x];
        work += n;
        for (uint32_t k = 0; k < n; ++k) {
            int32_t y = b[k].to;
            if (b[k].bw != INF && y != skip && mst[y] == mcur && (uint64_t)md[y] + b[k].bw <= (uint64_t)c) return true;
        }
        return false;
    }

    // Priority: level + edge quotient + hop quotient (2-hop witness estimate).
    int64_t priority(int32_t v) {
        const Arc* nb = adj(v);
        uint32_t n = asz[v];
        int64_t added = 0, added_hops = 0, removed = 0, removed_hops = 0;
        for (uint32_t i = 0; i < n; ++i) {
            int c = (nb[i].fw != INF) + (nb[i].bw != INF);
            removed += c; removed_hops += (int64_t)nb[i].hops * c;
        }
        for (uint32_t i = 0; i < n; ++i) {
            if (nb[i].bw == INF) continue;
            mark_out(nb[i].to, v);
            for (uint32_t j = 0; j < n; ++j) {
                if (j == i || nb[j].fw == INF) continue;
                DT c = (DT)std::min<uint64_t>((uint64_t)nb[i].bw + nb[j].fw, (uint64_t)CAP);
                if (!two_hop(nb[j].to, c, v)) { ++added; added_hops += nb[i].hops + nb[j].hops; }
            }
        }
        if (removed == 0) removed = 1;
        if (removed_hops == 0) removed_hops = 1;
        return (int64_t)lvl[v] * level_coef + (1000 * added) / removed + (1000 * added_hops) / removed_hops;
    }

    struct SC { int32_t u, x; DT w; int32_t hops; };
    vector<SC> sc;
    vector<uint32_t> pend;

    void contract(int32_t v) {
        uint32_t n = asz[v];
        Arc* nb = adj(v);
        for (uint32_t i = 0; i < n; ++i) {
            int32_t w = nb[i].to;
            Arc* l = adj(w);
            uint32_t m = asz[w];
            for (uint32_t k = 0; k < m; ++k)
                if (l[k].to == v) { l[k] = l[m - 1]; asz[w] = m - 1; break; }
        }
        sc.clear();
        for (uint32_t i = 0; i < n; ++i) {
            if (nb[i].bw == INF) continue;
            int32_t u = nb[i].to;
            DT wu = nb[i].bw;
            mark_out(u, v);
            pend.clear();
            new_targets();
            DT lim = 0;
            for (uint32_t j = 0; j < n; ++j) {
                if (j == i || nb[j].fw == INF) continue;
                DT c = addc(wu, nb[j].fw);
                if (two_hop(nb[j].to, c, v)) continue;
                pend.push_back(j);
                tst[nb[j].to] = tcur;
                if (c > lim) lim = c;
            }
            if (pend.empty()) continue;
            wsearch(u, lim, (int)pend.size());
            for (uint32_t j : pend) {
                int32_t x = nb[j].to;
                DT c = addc(wu, nb[j].fw);
                if (!(wst[x] == wcur && wd[x] <= c)) sc.push_back({u, x, c, nb[i].hops + nb[j].hops});
            }
        }
        tmp_head.push_back((uint32_t)tmp_up.size());
        for (uint32_t i = 0; i < n; ++i) {
            tmp_up.push_back({nb[i].to, nb[i].fw, nb[i].bw});
            int32_t w = nb[i].to;
            if (lvl[w] < lvl[v] + 1) lvl[w] = lvl[v] + 1;
        }
        asz[v] = 0;
        for (const SC& s : sc) {
            add_arc(s.u, s.x, s.w, INF, s.hops);
            add_arc(s.x, s.u, INF, s.w, s.hops);
        }
        order.push_back(v);
    }

    // ------------------------------------------------ exact cell elimination
    vector<int32_t> loc;
    vector<DT> W, dring, Df, Db;
    vector<int32_t> cl;
    vector<uint32_t> uhd;
    vector<int32_t> ula;
    vector<DT> ulf, ulb;
    vector<int> inl, outl, nbl, rnf, rnb;
    vector<uint32_t> rnfh, rnbh;
    vector<DT> rnfw, rnbw;
    vector<uint64_t> nbm;
    vector<uint32_t> rh;
    struct RA { int j; DT w; };
    vector<RA> ra;
    vector<DT> dringT, Mx;
    vector<DT> nbw[4], Tpad;
    struct XA { int i, j; DT w; };
    vector<XA> xarc;
    vector<int64_t> rkey;
    int S = 0;
    int dense_k = 6;
    uint64_t st_exits = 0, st_entries = 0, st_ring = 0, st_cells = 0;

    // same-cell queries: interior-only distance via the cell's elimination
    // structure (up-down paths whose top vertex is interior)
    vector<uint32_t> cqh, cql;
    struct SCQ { uint32_t qi; DT d; };
    vector<SCQ> scq;                      // interior-only distances of same-cell queries
    vector<DT> upS, upT;
    void same_cell_queries(size_t c, int ni) {
        if (cqh.empty()) return;
        if (upS.size() < (size_t)ni) { upS.assign((size_t)ni, INF); upT.assign((size_t)ni, INF); }
        for (uint32_t q = cqh[c]; q < cqh[c + 1]; ++q) {
            uint32_t qi = cql[q];
            int s0 = loc[qs[qi]], t0 = loc[qt[qi]];
            upS[s0] = 0; upT[t0] = 0;
            int lo = s0 < t0 ? s0 : t0;
            DT best = INF;
            for (int k = lo; k < ni; ++k) {
                DT us = upS[k], ut = upT[k];
                if (us == INF && ut == INF) continue;
                if (us != INF && ut != INF) { DT c2 = us + ut; if (c2 < best) best = c2; }
                for (uint32_t i = uhd[k]; i < uhd[k + 1]; ++i) {
                    int a = ula[i];
                    if (us != INF && ulf[i] != INF) { DT c2 = us + ulf[i]; if (c2 < upS[a]) upS[a] = c2; }
                    if (ut != INF && ulb[i] != INF) { DT c2 = ut + ulb[i]; if (c2 < upT[a]) upT[a] = c2; }
                }
                upS[k] = INF; upT[k] = INF;
            }
            if (best != INF) scq.push_back({qi, best});
        }
    }

    // Processes one cell: I = interior (elimination order), X = extra
    // boundary vertices (may be -1).
    void eliminate_cell(size_t cidx, const int32_t* I, int ni, const int32_t* X, int nx) {
        cl.assign(I, I + ni);
        for (int k = 0; k < ni; ++k) loc[I[k]] = k;
        for (int k = 0; k < ni; ++k) {
            const Arc* a = adj(I[k]);
            for (uint32_t j = 0; j < asz[I[k]]; ++j)
                if (loc[a[j].to] < 0) { loc[a[j].to] = (int32_t)cl.size(); cl.push_back(a[j].to); }
        }
        for (int k = 0; k < nx; ++k)
            if (X[k] >= 0 && !done[X[k]] && loc[X[k]] < 0) { loc[X[k]] = (int32_t)cl.size(); cl.push_back(X[k]); }
        const int n = (int)cl.size(), nb = n - ni;
        if (S > 0 && nb > 1) {
            int x0 = S, x1 = -1, y0 = S, y1 = -1;
            for (int k = 0; k < ni; ++k) {
                int x = I[k] % S, y = I[k] / S;
                x0 = std::min(x0, x); x1 = std::max(x1, x); y0 = std::min(y0, y); y1 = std::max(y1, y);
            }
            const int X0 = x0 - 1, X1 = x1 + 1, Y0 = y0 - 1, Y1 = y1 + 1, Wd = X1 - X0, Ht = Y1 - Y0;
            rkey.clear();
            for (int i = ni; i < n; ++i) {
                int g = cl[i], x = g % S, y = g / S;
                int64_t key;
                if (y == Y0 && x >= X0 && x <= X1) key = x - X0;
                else if (x == X1 && y >= Y0 && y <= Y1) key = Wd + (y - Y0);
                else if (y == Y1 && x >= X0 && x <= X1) key = Wd + Ht + (X1 - x);
                else if (x == X0 && y >= Y0 && y <= Y1) key = 2 * Wd + Ht + (Y1 - y);
                else key = 2 * (int64_t)(Wd + Ht) + 1 + g;
                rkey.push_back((key << 32) | (uint32_t)g);
            }
            std::sort(rkey.begin(), rkey.end());
            for (int i = ni; i < n; ++i) { cl[i] = (int32_t)(rkey[i - ni] & 0xFFFFFFFFu); loc[cl[i]] = i; }
        }
        ++st_cells;
        W.assign((size_t)n * n, INF);
        DT maxw = 0;
        const int words = (n + 63) >> 6;
        nbm.assign((size_t)n * words, 0);
        for (int i = 0; i < n; ++i) W[(size_t)i * n + i] = 0;
        for (int i = 0; i < n; ++i) {
            const Arc* a = adj(cl[i]);
            for (uint32_t j = 0; j < asz[cl[i]]; ++j) {
                int32_t t = loc[a[j].to];
                if (t < 0 || t == i) continue;
                DT* f = &W[(size_t)i * n + t];
                DT* b = &W[(size_t)t * n + i];
                if (a[j].fw < *f) *f = a[j].fw;
                if (a[j].bw < *b) *b = a[j].bw;
                if (a[j].fw != INF && a[j].fw > maxw) maxw = a[j].fw;
                if (a[j].bw != INF && a[j].bw > maxw) maxw = a[j].bw;
                nbm[(size_t)i * words + (t >> 6)] |= 1ull << (t & 63);
                nbm[(size_t)t * words + (i >> 6)] |= 1ull << (i & 63);
            }
        }
        // boundary-only all-pairs distances (Dijkstra over the sparse boundary graph)
        rh.assign((size_t)nb + 1, 0); ra.clear();
        for (int i = 0; i < nb; ++i) {
            const DT* row = &W[(size_t)(ni + i) * n + ni];
            for (int j = 0; j < nb; ++j) if (j != i && row[j] != INF) ra.push_back({j, row[j]});
            rh[i + 1] = (uint32_t)ra.size();
        }
        {
            // boundary arcs between cyclic neighbours in perimeter order, plus extras
            for (int q = 0; q < 4; ++q) nbw[q].assign((size_t)nb, INF);
            xarc.clear();
            for (int i = 0; i < nb; ++i)
                for (uint32_t e = rh[i]; e < rh[i + 1]; ++e) {
                    int j = ra[e].j;
                    DT w = ra[e].w;
                    if (j == (i + 1) % nb) { nbw[0][j] = w; nbw[2][i] = w; }        // i = j-1 -> j ; i -> i+1
                    else if (i == (j + 1) % nb) { nbw[1][j] = w; nbw[3][i] = w; }   // i = j+1 -> j ; i -> i-1
                    else xarc.push_back({i, j, w});
                }
        }
        // boundary-only all-pairs distances: Bellman-Ford rounds of one sweep
        // each way around the perimeter cycle plus the extra arcs
        dring.assign((size_t)nb * nb, INF);
        {
            const DT* wp = nbw[0].data();   // arc (j-1 -> j)
            const DT* wn = nbw[1].data();   // arc (j+1 -> j)
            for (int s0 = 0; s0 < nb; ++s0) {
                DT* D = &dring[(size_t)s0 * nb];
                D[s0] = 0;
                for (int round = 0; round < nb + 1; ++round) {
                    bool ch = false;
                    for (int step = 1; step <= nb; ++step) {
                        int j = s0 + step; if (j >= nb) j -= nb;
                        int i = j == 0 ? nb - 1 : j - 1;
                        DT c = sat_add(wp[j], D[i]);
                        if (c < D[j]) { D[j] = c; ch = true; }
                    }
                    for (int step = 1; step <= nb; ++step) {
                        int j = s0 - step; if (j < 0) j += nb;
                        int i = j + 1 == nb ? 0 : j + 1;
                        DT c = sat_add(wn[j], D[i]);
                        if (c < D[j]) { D[j] = c; ch = true; }
                    }
                    bool chx = false;
                    for (const XA& e : xarc) {
                        DT c = sat_add(e.w, D[e.i]);
                        if (c < D[e.j]) { D[e.j] = c; chx = true; }
                    }
                    if (xarc.empty() || (!ch && !chx)) break;   // a cycle alone is exact after one round
                }
            }
        }
        // bottom-up elimination of the interior (structure kept as bitsets)
        uhd.assign((size_t)ni + 1, 0);
        rnfh.assign((size_t)ni + 1, 0); rnbh.assign((size_t)ni + 1, 0);
        {
            const size_t cu = (size_t)ni * (size_t)ni / 2 + (size_t)ni, cr = (size_t)ni * (size_t)nb + 1;
            if (ula.size() < cu) { ula.resize(cu); ulf.resize(cu); ulb.resize(cu); }
            if (rnf.size() < cr) { rnf.resize(cr); rnfw.resize(cr); rnb.resize(cr); rnbw.resize(cr); }
            if (nbl.size() < (size_t)n) { nbl.resize(n); inl.resize(n); outl.resize(n); }
        }
        size_t nu = 0, nrf = 0, nrb = 0;
        int* NB = nbl.data();
        int* IN = inl.data();
        int* OUT = outl.data();
        for (int k = 0; k < ni; ++k) {
            int nn = 0, nin = 0, nout = 0;
            uhd[k] = (uint32_t)nu;
            rnfh[k] = (uint32_t)nrf; rnbh[k] = (uint32_t)nrb;
            const DT* rowk = &W[(size_t)k * n];
            uint64_t* bk = &nbm[(size_t)k * words];
            const int w0 = (k + 1) >> 6;
            for (int wi = w0; wi < words; ++wi) {
                uint64_t m = bk[wi];
                if (wi == w0) m &= ~0ull << ((k + 1) & 63);
                while (m) {
                    int x = (wi << 6) + __builtin_ctzll(m);
                    m &= m - 1;
                    DT f = rowk[x], b = W[(size_t)x * n + k];
                    if (f == INF && b == INF) continue;
                    NB[nn++] = x;
                    if (x < ni) { ula[nu] = x; ulf[nu] = f; ulb[nu] = b; ++nu; }
                    else {
                        if (f != INF) { rnf[nrf] = x - ni; rnfw[nrf] = f; ++nrf; }
                        if (b != INF) { rnb[nrb] = x - ni; rnbw[nrb] = b; ++nrb; }
                    }
                    if (f != INF) OUT[nout++] = x;
                    if (b != INF) IN[nin++] = x;
                }
            }
            // structural fill: the neighbours of k become pairwise adjacent
            for (int i = 0; i < nn; ++i) {
                uint64_t* bx = &nbm[(size_t)NB[i] * words];
                for (int wi = w0; wi < words; ++wi) bx[wi] |= bk[wi];
            }
            if (!nout || !nin) continue;
            // out-neighbours split into an interior range and a boundary range;
            // each range is updated densely (vectorised) when it is well filled
            int oi = 0;
            while (oi < nout && OUT[oi] < ni) ++oi;
            int rl[2] = {oi ? OUT[0] : 0, oi < nout ? OUT[oi] : 0};
            int rr[2] = {oi ? OUT[oi - 1] + 1 : 0, oi < nout ? OUT[nout - 1] + 1 : 0};
            int rc[2] = {oi, nout - oi};
            int rs[2] = {0, oi};
            bool dn[2];
            for (int q = 0; q < 2; ++q) dn[q] = rc[q] * dense_k >= rr[q] - rl[q];
            for (int ii = 0; ii < nin; ++ii) {
                const int u = IN[ii];
                const DT wu = W[(size_t)u * n + k];
                DT* rowu = &W[(size_t)u * n];
                const DT keep = rowu[u];
                for (int q = 0; q < 2; ++q) {
                    if (!rc[q]) continue;
                    if (dn[q]) {
                        for (int x = rl[q]; x < rr[q]; ++x) {
                            DT c = sat_add(wu, rowk[x]);
                            rowu[x] = c < rowu[x] ? c : rowu[x];
                        }
                    } else {
                        const int* ol = OUT + rs[q];
                        for (int i = 0; i < rc[q]; ++i) {
                            int x = ol[i];
                            DT c = wu + rowk[x];
                            if (c < rowu[x]) rowu[x] = c;
                        }
                    }
                }
                rowu[u] = keep;
            }
        }
        uhd[ni] = (uint32_t)nu;
        rnfh[ni] = (uint32_t)nrf; rnbh[ni] = (uint32_t)nrb;
        // every value computed in this cell is a path of < 2n arcs of this cell
        if ((uint64_t)maxw * (uint64_t)(2 * n + 2) > (uint64_t)CAP) {
            for (size_t i = 0; i < W.size(); ++i) if (W[i] > CAP / 2 && W[i] != INF) { overflow = true; break; }
            for (size_t i = 0; i < dring.size(); ++i) if (dring[i] > CAP / 2 && dring[i] != INF) { overflow = true; break; }
        }
        same_cell_queries(cidx, ni);
        // top-down closures: Tf[v][b] = min over boundary b' of (first-exit
        // distance v->b') + (boundary path b'->b); Tb symmetric for entries.
        // b is a useful exit of v iff Tf[v][b] is not matched by arriving at b
        // along a boundary arc (then Tf[v][b] is the first-exit distance).
        for (int i = 0; i < nb; ++i) { dring[(size_t)i * nb + i] = 0; }
        dringT.resize((size_t)nb * nb);
        for (int i = 0; i < nb; ++i)
            for (int j = 0; j < nb; ++j) dringT[(size_t)j * nb + i] = dring[(size_t)i * nb + j];
        Df.assign((size_t)ni * nb, INF);
        Db.assign((size_t)ni * nb, INF);
        Mx.resize((size_t)nb);
        Tpad.resize((size_t)nb + 2);
        for (int k = ni - 1; k >= 0; --k) {
            const int32_t g = cl[k];
            for (int dir = 0; dir < 2; ++dir) {
                DT* T = dir == 0 ? &Df[(size_t)k * nb] : &Db[(size_t)k * nb];
                const DT* R = dir == 0 ? dring.data() : dringT.data();
                const vector<int>& rl_ = dir == 0 ? rnf : rnb;
                const vector<DT>& rw_ = dir == 0 ? rnfw : rnbw;
                const vector<uint32_t>& rh_ = dir == 0 ? rnfh : rnbh;
                const vector<DT>& uw = dir == 0 ? ulf : ulb;
                const vector<DT>& TT = dir == 0 ? Df : Db;
                for (uint32_t i = uhd[k]; i < uhd[k + 1]; ++i) {
                    const DT w = uw[i];
                    if (w == INF) continue;
                    const DT* ta = &TT[(size_t)ula[i] * nb];
                    for (int b = 0; b < nb; ++b) { DT c = sat_add(w, ta[b]); T[b] = c < T[b] ? c : T[b]; }
                }
                // T is closed under boundary paths; a direct boundary arc k->b0
                // only matters if it improves T[b0]
                for (uint32_t i = rh_[k]; i < rh_[k + 1]; ++i) {
                    const DT w = rw_[i];
                    if (w >= T[rl_[i]]) continue;
                    const DT* r = R + (size_t)rl_[i] * nb;
                    for (int b = 0; b < nb; ++b) { DT c = sat_add(w, r[b]); T[b] = c < T[b] ? c : T[b]; }
                }
                // best arrival at each b over one boundary arc
                DT* M = Mx.data();
                {
                    DT* Tp = Tpad.data();       // Tp[i + 1] = T[i], cyclic ends
                    Tp[0] = T[nb - 1]; Tp[nb + 1] = T[0];
                    for (int b = 0; b < nb; ++b) Tp[b + 1] = T[b];
                    const DT* wa = nbw[dir == 0 ? 0 : 3].data();   // weight of the arc to/from b-1
                    const DT* wb = nbw[dir == 0 ? 1 : 2].data();   // weight of the arc to/from b+1
                    for (int b = 0; b < nb; ++b) {
                        DT c1 = sat_add(wa[b], Tp[b]), c2 = sat_add(wb[b], Tp[b + 2]);
                        M[b] = c1 < c2 ? c1 : c2;
                    }
                    for (const XA& e : xarc) {
                        if (dir == 0) { DT c = sat_add(e.w, T[e.i]); if (c < M[e.j]) M[e.j] = c; }
                        else { DT c = sat_add(e.w, T[e.j]); if (c < M[e.i]) M[e.i] = c; }
                    }
                }
                vector<XE>& L = dir == 0 ? exl : enl;
                size_t& Ln = dir == 0 ? exl_n : enl_n;
                if (L.size() < Ln + (size_t)nb) L.resize(std::max(L.size() * 2, Ln + (size_t)nb + 1024));
                XE* o = L.data() + Ln;
                uint32_t cnt = 0;
                const int32_t* rg = cl.data() + ni;
                for (int b = 0; b < nb; ++b) {
                    o[cnt] = {rg[b], T[b]};
                    cnt += T[b] < M[b];
                }
                uint32_t st = (uint32_t)Ln;
                Ln += cnt;
                if (dir == 0) { exs[g] = st; exn[g] = cnt; st_exits += cnt; }
                else { ens[g] = st; enn[g] = cnt; st_entries += cnt; }
            }
            iord.push_back(g);
            done[g] = 1;
        }
        // detach the interior from the boundary; boundary shortcuts through it
        for (int i = ni; i < n; ++i) {
            int32_t g = cl[i];
            Arc* a = adj(g);
            uint32_t m = asz[g];
            for (uint32_t j = 0; j < m;) {
                int32_t t = loc[a[j].to];
                if (t >= 0 && t < ni) a[j] = a[--m]; else ++j;
            }
            asz[g] = m;
        }
        for (int i = 0; i < nb; ++i) {
            const DT* row = &W[(size_t)(ni + i) * n + ni];
            const DT* dr = &dring[(size_t)i * nb];
            for (int j = 0; j < nb; ++j) {
                if (i == j || row[j] >= dr[j]) continue;
                add_arc(cl[ni + i], cl[ni + j], row[j], INF, 2);
                add_arc(cl[ni + j], cl[ni + i], INF, row[j], 2);
                ++st_ring;
            }
        }
        for (int32_t g : cl) loc[g] = -1;
    }

    void run_cells(const Plan& p) {
        exs.assign((size_t)V, 0); ens.assign((size_t)V, 0);
        exn.assign((size_t)V, 0); enn.assign((size_t)V, 0);
        if (p.cellst.empty()) return;
        loc.assign((size_t)V, -1);
        {
            // bucket same-cell queries by cell
            size_t nc = p.cellst.size() - 1;
            cqh.assign(nc + 1, 0);
            for (int32_t i = 0; i < Q; ++i) {
                int32_t a = p.cell_of[qs[i]];
                if (a >= 0 && a == p.cell_of[qt[i]] && qs[i] != qt[i]) ++cqh[a + 1];
            }
            for (size_t c = 0; c < nc; ++c) cqh[c + 1] += cqh[c];
            cql.resize(cqh[nc]);
            vector<uint32_t> cur(cqh.begin(), cqh.end() - 1);
            for (int32_t i = 0; i < Q; ++i) {
                int32_t a = p.cell_of[qs[i]];
                if (a >= 0 && a == p.cell_of[qt[i]] && qs[i] != qt[i]) cql[cur[a]++] = (uint32_t)i;
            }
        }
        for (size_t c = 0; c + 1 < p.cellst.size(); ++c) {
            iost.push_back((uint32_t)iord.size());
            eliminate_cell(c, p.cellv.data() + p.cellst[c], (int)(p.cellst[c + 1] - p.cellst[c]),
                           p.cellx.data() + 4 * c, 4);
        }
        iost.push_back((uint32_t)iord.size());
        vector<int32_t>().swap(loc);
        if (g_debug) std::fprintf(stderr, "  [%7.3f] cells %llu: exits %llu entries %llu boundary shortcuts %llu\n",
                                  now_s(), (unsigned long long)st_cells, (unsigned long long)st_exits,
                                  (unsigned long long)st_entries, (unsigned long long)st_ring);
    }

    // ------------------------------------------------------ network CH
    void contract_network() {
        order.clear();
        tmp_up.clear(); tmp_head.clear();
        vector<uint8_t> dirty((size_t)V, 0);
        vector<int64_t> pri((size_t)V, 0);
        typedef std::pair<int64_t, int32_t> PQ;
        vector<PQ> init;
        for (int32_t v = 0; v < V; ++v) if (!done[v]) { pri[v] = priority(v); init.push_back({pri[v], v}); }
        std::priority_queue<PQ, vector<PQ>, std::greater<PQ>> pq(std::greater<PQ>(), std::move(init));
        if (g_debug) std::fprintf(stderr, "  [%7.3f] network vertices %zu\n", now_s(), pq.size());
        while (!pq.empty()) {
            if (work + n_settle * 8 > work_budget) { aborted = true; return; }
            auto [p, v] = pq.top(); pq.pop();
            if (done[v] || p != pri[v]) continue;
            if (dirty[v]) {
                dirty[v] = 0;
                int64_t np = priority(v);
                pri[v] = np;
                if (!pq.empty() && np > pq.top().first) { pq.push({np, v}); continue; }
            }
            const Arc* nb = adj(v);
            for (uint32_t i = 0; i < asz[v]; ++i) dirty[nb[i].to] = 1;
            contract(v); done[v] = 1;
        }
        tmp_head.push_back((uint32_t)tmp_up.size());
        vector<Arc>().swap(pool);
        NR = (int32_t)order.size();
        rank_of.assign((size_t)V, -1);
        for (int32_t r = 0; r < NR; ++r) rank_of[order[r]] = r;
        uhead.assign((size_t)NR + 1, 0);
        ua.resize(tmp_up.size());
        for (int32_t r = 0; r < NR; ++r) {
            uhead[r] = tmp_head[r];
            for (uint32_t k = tmp_head[r]; k < tmp_head[r + 1]; ++k) {
                UArc a = tmp_up[k];
                a.to = rank_of[a.to];
                ua[k] = a;
            }
        }
        uhead[NR] = tmp_head[NR];
        vector<UArc>().swap(tmp_up);
        vector<uint32_t>().swap(tmp_head);
        if (g_debug) std::fprintf(stderr, "  [%7.3f] network CH: %d vertices, %u up arcs, searches %llu settles %llu work %llu (budget %llu)\n",
                                  now_s(), NR, uhead[NR], (unsigned long long)n_search, (unsigned long long)n_settle,
                                  (unsigned long long)(work + n_settle * 8), (unsigned long long)work_budget);
    }

    // ------------------------------------------------------------ labels
    struct LE { uint32_t h; DT d; };      // hub = network rank
    struct LBuf {                         // growable POD array
        LE* p = nullptr;
        size_t n = 0, cap = 0;
        ~LBuf() { std::free(p); }
        inline LE* need(size_t k) {
            if (n + k > cap) {
                size_t nc = std::max(cap * 2, n + k + 4096);
                LE* q = (LE*)std::realloc(p, nc * sizeof(LE));
                if (!q) { std::fprintf(stderr, "out of memory\n"); std::exit(1); }
                p = q; cap = nc;
            }
            return p + n;
        }
        inline LE* data() { return p; }
        inline const LE* data() const { return p; }
        inline size_t size() const { return n; }
    };
    vector<uint64_t> fst, bst;            // per vertex id
    vector<uint32_t> flen, blen;
    LBuf fl, bl;
    vector<DT> tmp;
    vector<uint32_t> hcand;
    size_t nhc = 0;

    // appends the pruned candidate set (tmp/hcand) as the label of v
    int stall_int = 1;
    inline void emit(int dir, int32_t v, bool stall = true) {
        LBuf& L = dir == 0 ? fl : bl;
        LE* o = L.need(nhc);
        const uint32_t* H = uhead.data();
        const UArc* A = ua.data();
        uint32_t cnt = 0;
        for (size_t i = 0; i < nhc; ++i) {
            const uint32_t h = hcand[i];
            const DT d = tmp[h];
            bool keep = true;
            if (stall) for (uint32_t k = H[h]; k < H[h + 1]; ++k) {
                DT w = dir == 0 ? A[k].bw : A[k].fw;
                if (w == INF) continue;
                DT tx = tmp[A[k].to];
                if (tx != INF && (uint64_t)tx + w < (uint64_t)d) { keep = false; break; }
            }
            if (keep) {
                if (d > CAP) overflow = true;
                o[cnt++] = {h, d};
            }
        }
        (dir == 0 ? fst : bst)[v] = L.n;
        (dir == 0 ? flen : blen)[v] = cnt;
        L.n += cnt;
        for (size_t i = 0; i < nhc; ++i) tmp[hcand[i]] = INF;
        nhc = 0;
    }
    inline void merge_label(int dir, int32_t u, DT w) {
        const LBuf& L = dir == 0 ? fl : bl;
        const LE* e = L.data() + (dir == 0 ? fst : bst)[u];
        uint32_t m = (dir == 0 ? flen : blen)[u];
        uint32_t* hc = hcand.data();
        for (uint32_t i = 0; i < m; ++i) {
            uint32_t h = e[i].h;
            DT nd = w + e[i].d;
            if (tmp[h] == INF) { hc[nhc++] = h; tmp[h] = nd; }
            else if (nd < tmp[h]) tmp[h] = nd;
        }
    }

    void build_labels() {
        fst.assign((size_t)V, 0); bst.assign((size_t)V, 0);
        flen.assign((size_t)V, 0); blen.assign((size_t)V, 0);
        fl.need((size_t)V * 44); bl.need((size_t)V * 44);
        tmp.assign((size_t)NR, INF);
        hcand.assign((size_t)NR + 1, 0);
        nhc = 0;
        const uint32_t* H = uhead.data();
        const UArc* A = ua.data();
        for (int32_t r = NR - 1; r >= 0; --r) {
            int32_t v = order[r];
            for (int dir = 0; dir < 2; ++dir) {
                tmp[r] = 0; hcand[nhc++] = (uint32_t)r;
                for (uint32_t k = H[r]; k < H[r + 1]; ++k) {
                    DT w = dir == 0 ? A[k].fw : A[k].bw;
                    if (w == INF) continue;
                    merge_label(dir, order[A[k].to], w);
                }
                emit(dir, v);
            }
        }
        tlog("network labels");
        if (g_debug) std::fprintf(stderr, "network label entries: fwd %zu bwd %zu (avg %.1f)\n",
                                  fl.size(), bl.size(), (double)(fl.size() + bl.size()) / (2.0 * std::max(1, NR)));
        // Interior labels, cell by cell, over a cell-local hub index space:
        // the labels of the cell's exit/entry vertices become dense rows over
        // the local hubs; a label is the min over its exits of (d + row).
        vector<int32_t> hloc((size_t)NR, -1), lh;
        vector<int32_t> bloc((size_t)V, -1), bl_used;
        vector<uint8_t> need((size_t)V, 0);      // bit 0: forward label, bit 1: backward
        for (int32_t i = 0; i < Q; ++i) { need[qs[i]] |= 1; need[qt[i]] |= 2; }
        vector<DT> lt;
        vector<DT> dm[2];
        for (size_t c = 0; c + 1 < iost.size(); ++c) {
            lh.clear(); bl_used.clear();
            for (uint32_t i = iost[c]; i < iost[c + 1]; ++i) {
                int32_t v = iord[i];
                for (int dir = 0; dir < 2; ++dir) {
                    if (!(need[v] >> dir & 1)) continue;
                    const XE* x = dir == 0 ? &exl[exs[v]] : &enl[ens[v]];
                    uint32_t m = dir == 0 ? exn[v] : enn[v];
                    for (uint32_t j = 0; j < m; ++j) {
                        int32_t b = x[j].b;
                        if (bloc[b] >= 0) continue;
                        bloc[b] = (int32_t)bl_used.size();
                        bl_used.push_back(b);
                    }
                }
            }
            for (int32_t b : bl_used)
                for (int dir = 0; dir < 2; ++dir) {
                    const LE* e = (dir == 0 ? fl.data() + fst[b] : bl.data() + bst[b]);
                    uint32_t m = dir == 0 ? flen[b] : blen[b];
                    for (uint32_t j = 0; j < m; ++j)
                        if (hloc[e[j].h] < 0) { hloc[e[j].h] = (int32_t)lh.size(); lh.push_back((int32_t)e[j].h); }
                }
            const size_t Hn = lh.size(), Hp = (Hn + 7) & ~(size_t)7;
            for (int dir = 0; dir < 2; ++dir) {
                dm[dir].assign(bl_used.size() * Hp, INF);
                for (size_t bi = 0; bi < bl_used.size(); ++bi) {
                    int32_t b = bl_used[bi];
                    const LE* e = (dir == 0 ? fl.data() + fst[b] : bl.data() + bst[b]);
                    uint32_t m = dir == 0 ? flen[b] : blen[b];
                    DT* row = &dm[dir][bi * Hp];
                    for (uint32_t j = 0; j < m; ++j) row[hloc[e[j].h]] = e[j].d;
                }
            }
            lt.assign(Hp, INF);
            DT mx = 0;
            for (uint32_t i = iost[c]; i < iost[c + 1]; ++i) {
                int32_t v = iord[i];
                for (int dir = 0; dir < 2; ++dir) {
                    if (!(need[v] >> dir & 1)) continue;
                    const XE* x = dir == 0 ? &exl[exs[v]] : &enl[ens[v]];
                    uint32_t m = dir == 0 ? exn[v] : enn[v];
                    DT* T = lt.data();
                    if (m) {
                        const DT w0 = x[0].d;
                        const DT* row0 = &dm[dir][(size_t)bloc[x[0].b] * Hp];
                        for (size_t h = 0; h < Hp; ++h) T[h] = sat_add(w0, row0[h]);
                    } else {
                        for (size_t h = 0; h < Hp; ++h) T[h] = INF;
                    }
                    for (uint32_t j = 1; j < m; ++j) {
                        const DT w = x[j].d;
                        const DT* row = &dm[dir][(size_t)bloc[x[j].b] * Hp];
                        for (size_t h = 0; h < Hp; ++h) { DT c2 = sat_add(w, row[h]); T[h] = c2 < T[h] ? c2 : T[h]; }
                    }
                    LBuf& L = dir == 0 ? fl : bl;
                    LE* o = L.need(Hn);
                    uint32_t cnt = 0;
                    for (size_t h = 0; h < Hn; ++h) {
                        const DT t = T[h];
                        o[cnt] = {(uint32_t)lh[h], t};
                        const bool fin = t != INF;
                        cnt += fin;
                        const DT tf = fin ? t : 0;
                        mx = tf > mx ? tf : mx;
                    }
                    (dir == 0 ? fst : bst)[v] = L.n;
                    (dir == 0 ? flen : blen)[v] = cnt;
                    L.n += cnt;
                }
            }
            if (mx > CAP) overflow = true;
            for (int32_t h : lh) hloc[h] = -1;
            for (int32_t b : bl_used) bloc[b] = -1;
        }
        if (g_debug) std::fprintf(stderr, "label entries: fwd %zu bwd %zu (avg %.1f)\n",
                                  fl.size(), bl.size(), (double)(fl.size() + bl.size()) / (2.0 * V));
    }

    // ------------------------------------------------------------ queries
    void answer() {
        vector<DT> D((size_t)NR, INF);
        // queries in source order: target and original index, sequential
        vector<uint32_t> head((size_t)V + 1, 0);
        for (int32_t i = 0; i < Q; ++i) ++head[qs[i] + 1];
        for (int32_t v = 0; v < V; ++v) head[v + 1] += head[v];
        vector<int32_t> st((size_t)Q + 16, 0);
        vector<uint32_t> si((size_t)Q);
        {
            vector<uint32_t> cur(head.begin(), head.end() - 1);
            for (int32_t i = 0; i < Q; ++i) { uint32_t k = cur[qs[i]]++; st[k] = qt[i]; si[k] = (uint32_t)i; }
        }
        for (int k = 0; k < 16; ++k) st[(size_t)Q + k] = st.empty() || Q == 0 ? 0 : st[(size_t)Q - 1];
        vector<int64_t> sa((size_t)Q);
        const int PD = 8;   // prefetch distance (queries)
        const LE* FL = fl.data();
        const LE* BL = bl.data();
        int32_t nexts = 0;
        for (int32_t s = 0; s < V; ++s) {
            uint32_t a = head[s], b = head[s + 1];
            if (a == b) continue;
            // prefetch the forward label of the next source
            nexts = s + 1;
            while (nexts < V && head[nexts] == head[nexts + 1]) ++nexts;
            if (nexts < V) {
                const char* pf = (const char*)(FL + fst[nexts]);
                __builtin_prefetch(pf); __builtin_prefetch(pf + 64); __builtin_prefetch(pf + 128);
                __builtin_prefetch(pf + 192); __builtin_prefetch(pf + 256);
            }
            const LE* F = FL + fst[s];
            uint32_t nf = flen[s];
            for (uint32_t i = 0; i < nf; ++i) D[F[i].h] = F[i].d;
            for (uint32_t k = a; k < b; ++k) {
                {
                    const char* pb = (const char*)(BL + bst[st[k + PD]]);
                    __builtin_prefetch(pb); __builtin_prefetch(pb + 64); __builtin_prefetch(pb + 128);
                    __builtin_prefetch(pb + 192); __builtin_prefetch(pb + 256);
                }
                int32_t t = st[k];
                if (t == s) { sa[k] = 0; continue; }
                const LE* B = BL + bst[t];
                uint32_t nb = blen[t];
                uint64_t best = ~0ull;
                if (sizeof(DT) == 4) {
                    for (uint32_t i = 0; i < nb; ++i) {
                        uint64_t c = (uint64_t)D[B[i].h] + (uint64_t)B[i].d;
                        best = c < best ? c : best;
                    }
                    if (best >= (uint64_t)INF) best = ~0ull;
                } else {
                    for (uint32_t i = 0; i < nb; ++i) {
                        DT x = D[B[i].h];
                        if (x != INF) {
                            uint64_t c = (uint64_t)x + B[i].d;
                            if (c < best) best = c;
                        }
                    }
                }
                sa[k] = best == ~0ull ? -1 : (int64_t)best;
            }
            for (uint32_t i = 0; i < nf; ++i) D[F[i].h] = INF;
        }
        for (int32_t k = 0; k < Q; ++k) qans[si[k]] = sa[k];
        for (const SCQ& c : scq)
            if (qans[c.qi] < 0 || (uint64_t)c.d < (uint64_t)qans[c.qi]) qans[c.qi] = (int64_t)c.d;
    }
};

// Fallback for graphs a hierarchy does not suit: one Dijkstra per distinct
// source, stopped once all of that source's targets are settled.
static void dijkstra_queries() {
    vector<uint32_t> h((size_t)V + 1, 0);
    for (int32_t i = 0; i < E; ++i) ++h[eu[i] + 1];
    for (int32_t v = 0; v < V; ++v) h[v + 1] += h[v];
    vector<int32_t> to((size_t)E);
    vector<uint32_t> w((size_t)E);
    {
        vector<uint32_t> cur(h.begin(), h.end() - 1);
        for (int32_t i = 0; i < E; ++i) { uint32_t k = cur[eu[i]]++; to[k] = ev[i]; w[k] = ew[i]; }
    }
    vector<uint32_t> qh((size_t)V + 1, 0), idx((size_t)Q);
    for (int32_t i = 0; i < Q; ++i) ++qh[qs[i] + 1];
    for (int32_t v = 0; v < V; ++v) qh[v + 1] += qh[v];
    {
        vector<uint32_t> cur(qh.begin(), qh.end() - 1);
        for (int32_t i = 0; i < Q; ++i) idx[cur[qs[i]]++] = (uint32_t)i;
    }
    const uint64_t INF64 = ~0ull;
    vector<uint64_t> dist((size_t)V, INF64);
    vector<uint32_t> want((size_t)V, 0);
    vector<int32_t> touched;
    Heap4<uint64_t> heap;
    for (int32_t s = 0; s < V; ++s) {
        if (qh[s] == qh[s + 1]) continue;
        int pending = 0;
        for (uint32_t k = qh[s]; k < qh[s + 1]; ++k) { int32_t t = qt[idx[k]]; if (!want[t]++) ++pending; }
        heap.clear();
        dist[s] = 0; touched.push_back(s); heap.push(0, s);
        while (!heap.empty() && pending > 0) {
            auto it = heap.pop();
            if (it.d > dist[it.v]) continue;
            if (want[it.v]) --pending;
            for (uint32_t k = h[it.v]; k < h[it.v + 1]; ++k) {
                uint64_t nd = it.d + w[k];
                if (nd < dist[to[k]]) {
                    if (dist[to[k]] == INF64) touched.push_back(to[k]);
                    dist[to[k]] = nd; heap.push(nd, to[k]);
                }
            }
        }
        for (uint32_t k = qh[s]; k < qh[s + 1]; ++k) {
            int32_t t = qt[idx[k]];
            qans[idx[k]] = dist[t] == INF64 ? -1 : (int64_t)dist[t];
            want[t] = 0;
        }
        for (int32_t v : touched) dist[v] = INF64;
        touched.clear();
    }
}

template <class DT>
static bool run(const Plan& plan, bool& aborted) {
    Engine<DT> en;
    en.con_settle = env_int("CON_SETTLE", 300);
    en.level_coef = env_int("LEVEL_COEF", 1000);
    en.stall_int = env_int("STALL_INT", 1);
    en.cell_of = &plan.cell_of;
    en.S = plan.S;
    en.dense_k = env_int("DENSE_K", 6);
    en.build_graph();
    tlog("graph built");
    en.run_cells(plan);
    tlog("cells");
    if (en.overflow) return false;
    en.work_budget = (uint64_t)env_int("CH_BUDGET", 60) * ((uint64_t)V + (uint64_t)E) + 20000000ull;
    en.contract_network();
    tlog("network CH");
    if (en.aborted) { tlog("CH effort budget exceeded"); aborted = true; return false; }
    if (en.overflow) return false;
    en.build_labels();
    tlog("labels");
    if (en.overflow) return false;
    en.answer();
    tlog("queries");
    return true;
}

bool applicable() { return directed && has_coords; }

void solve() {
    Plan plan = make_plan();
    tlog("plan");
    bool aborted = false;
    if (env_int("FORCE64", 0) || !run<uint32_t>(plan, aborted)) {
        if (!aborted) {
            tlog("32-bit overflow, redo with 64-bit");
            run<uint64_t>(plan, aborted);
        }
        if (aborted) { dijkstra_queries(); tlog("dijkstra queries"); }
    }
}

}  // namespace road

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: %s <graph_file> <query_file> <output_file>\n", argv[0]);
        return 1;
    }
    g_debug = std::getenv("SPC_DEBUG") != nullptr;
    g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]); tlog("graph");
    read_queries(argv[2]); tlog("queries read");
    road::solve();
    write_answers(argv[3]); tlog("write");
    std::_Exit(0);
}
