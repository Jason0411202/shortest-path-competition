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
        if (h != (uint32_t)r) {
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

