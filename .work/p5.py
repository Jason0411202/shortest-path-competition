s=open('sp.cpp','rb').read().decode()
def rep(old,new):
    global s
    assert s.count(old)==1, old[:80]
    s=s.replace(old,new)
rep('''    W limit = 0;
    for (int32_t x : tlist) limit = std::max(limit, tcost[x]);
    int pending = (int)tlist.size();''','''    W limit = 0;
    for (int32_t x : tlist) limit = std::max(limit, tlim[x]);
    int pending = (int)tlist.size();
    if (pending == 0) return;''')
rep('''        if (tstamp[it.v] == tcur) {
            tstamp[it.v] = 0;
            if (--pending <= 0) break;
            if (tcost[it.v] >= limit) {
                limit = 0;
                for (int32_t x : tlist) if (tstamp[x] == tcur) limit = std::max(limit, tcost[x]);
                if (it.d > limit) break;
            }
        }
        if (++settled > max_settle) break;''','''        if (++settled > max_settle) break;''')
rep('''            W nd = it.d + a.fw;
            if (nd > limit) continue;
            if (wstamp[a.to] != wcur || nd < wdist[a.to]) {
                wstamp[a.to] = wcur;
                wdist[a.to] = nd;
                wheap.push(nd, a.to);
            }
        }
    }
}''','''            W nd = it.d + a.fw;
            int32_t x = a.to;
            if (tstamp[x] == tcur && nd <= tcost[x]) {
                // witness for x found (x is reached through its last arc)
                tstamp[x] = 0;
                wstamp[x] = wcur; wdist[x] = nd;
                if (--pending <= 0) return;
                if (tlim[x] >= limit) {
                    limit = 0;
                    for (int32_t y : tlist) if (tstamp[y] == tcur) limit = std::max(limit, tlim[y]);
                }
                continue;
            }
            if (nd > limit) continue;
            if (wstamp[x] != wcur || nd < wdist[x]) {
                wstamp[x] = wcur;
                wdist[x] = nd;
                wheap.push(nd, x);
            }
        }
    }
}

// Smallest weight of an arc entering x other than from `avoid`.
static inline W min_in(int32_t x, int32_t avoid) {
    W m = INF;
    for (const DArc& b : g[x]) if (b.bw < m && b.to != avoid) m = b.bw;
    return m;
}''')
rep('static vector<W> tcost;','static vector<W> tcost, tlim;')
rep('''                if (a.fw == 0) continue;
                tstamp[a.to] = tcur; tcost[a.to] = a.fw - 1; tlist.push_back(a.to);''','''                if (a.fw == 0) continue;
                W mi = min_in(a.to, -1);
                if (mi > a.fw - 1) continue;   // every path ends with an arc >= a.fw
                tstamp[a.to] = tcur; tcost[a.to] = a.fw - 1; tlim[a.to] = a.fw - 1 - mi; tlist.push_back(a.to);''')
rep('''        new_targets();
        for (size_t j = j0; j < n; ++j)
            if (j != i && nb[j].fw < INF) {
                int32_t x = nb[j].to;
                tstamp[x] = tcur; tcost[x] = wu + nb[j].fw; tlist.push_back(x);
            }''','''        new_targets();
        for (size_t j = j0; j < n; ++j)
            if (j != i && nb[j].fw < INF) {
                int32_t x = nb[j].to;
                W c = wu + nb[j].fw;
                if (minin[j] > c) continue;
                tstamp[x] = tcur; tcost[x] = c; tlim[x] = c - minin[j]; tlist.push_back(x);
            }''')
rep('''    scs.clear();
    const vector<DArc>& nb = g[v];
    size_t n = nb.size();
    for (size_t i = 0; i < n; ++i) {''','''    scs.clear();
    const vector<DArc>& nb = g[v];
    size_t n = nb.size();
    static vector<W> minin;
    minin.resize(n);
    for (size_t j = 0; j < n; ++j) minin[j] = nb[j].fw < INF ? min_in(nb[j].to, v) : INF;
    for (size_t i = 0; i < n; ++i) {''')
rep('    tcost.assign((size_t)V, 0);','    tcost.assign((size_t)V, 0);\n    tlim.assign((size_t)V, 0);')
open('sp.cpp','wb').write(s.encode())
