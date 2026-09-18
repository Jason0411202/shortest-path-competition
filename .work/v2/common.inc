// ===================================================================== common
// Shared I/O layer.  Pasted verbatim at the top of every family prototype and
// of the final solver.cpp (no #include of own headers in the final file).
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <vector>
#include <type_traits>
#include <cmath>

using std::vector;

static const int32_t FLAG_COORDS   = 1;
static const int32_t FLAG_DIRECTED = 2;

static bool g_debug = false;
static std::chrono::steady_clock::time_point g_t0;
static void tlog(const char* what) {
    if (!g_debug) return;
    double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - g_t0).count();
    std::fprintf(stderr, "[%8.3f s] %s\n", s, what);
}

// Whole file in one buffer, NUL-terminated (the parser relies on a non-digit
// sentinel after the last number).
static char* slurp(const char* path, size_t& n) {
    std::FILE* f = std::fopen(path, "rb");
    if (!f) { std::fprintf(stderr, "cannot open file: %s\n", path); std::exit(1); }
    std::fseek(f, 0, SEEK_END);
    long sz = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    if (sz < 0) sz = 0;
    char* b = (char*)std::malloc((size_t)sz + 16);
    n = std::fread(b, 1, (size_t)sz, f);
    std::memset(b + n, 0, 16);
    std::fclose(f);
    return b;
}

// Reads the next non-negative integer; skips any non-digit characters first.
static inline uint64_t rd_u(const char*& p) {
    while ((unsigned)(*p - '0') > 9u) { if (!*p) return 0; ++p; }
    uint64_t x = (uint64_t)(*p++ - '0');
    unsigned c;
    while ((c = (unsigned)(*p - '0')) <= 9u) { x = x * 10 + c; ++p; }
    return x;
}

static int32_t V = 0, E = 0, FLAGS = 0;
static bool directed = false, has_coords = false;
static vector<int32_t> eu, ev;     // edge endpoints as listed in the file
static vector<uint32_t> ew;        // weights (positive, fit in int32)
static vector<int64_t> cx, cy;     // coordinates (only if FLAGS bit 0)
static int32_t Q = 0;
static vector<int32_t> qs, qt;
static vector<int64_t> qans;       // -1 = unreachable

static void read_graph(const char* path) {
    size_t n;
    char* buf = slurp(path, n);
    const char* p = buf;
    V = (int32_t)rd_u(p); E = (int32_t)rd_u(p); FLAGS = (int32_t)rd_u(p);
    directed = (FLAGS & FLAG_DIRECTED) != 0;
    has_coords = (FLAGS & FLAG_COORDS) != 0;
    eu.resize((size_t)E); ev.resize((size_t)E); ew.resize((size_t)E);
    for (int32_t i = 0; i < E; ++i) {
        eu[i] = (int32_t)rd_u(p); ev[i] = (int32_t)rd_u(p); ew[i] = (uint32_t)rd_u(p);
    }
    if (has_coords) {
        cx.resize((size_t)V); cy.resize((size_t)V);
        for (int32_t i = 0; i < V; ++i) { cx[i] = (int64_t)rd_u(p); cy[i] = (int64_t)rd_u(p); }
    }
    std::free(buf);
}

static void read_queries(const char* path) {
    size_t n;
    char* buf = slurp(path, n);
    const char* p = buf;
    Q = (int32_t)rd_u(p);
    qs.resize((size_t)Q); qt.resize((size_t)Q); qans.assign((size_t)Q, -1);
    for (int32_t i = 0; i < Q; ++i) { qs[i] = (int32_t)rd_u(p); qt[i] = (int32_t)rd_u(p); }
    std::free(buf);
}

static void write_answers(const char* path) {
    char* out = (char*)std::malloc((size_t)Q * 21 + 64);
    char* o = out;
    for (int32_t i = 0; i < Q; ++i) {
        int64_t d = qans[i];
        if (d < 0) { *o++ = '-'; *o++ = '1'; }
        else {
            char tmp[24]; int k = 0; uint64_t x = (uint64_t)d;
            do { tmp[k++] = (char)('0' + x % 10); x /= 10; } while (x);
            while (k) *o++ = tmp[--k];
        }
        *o++ = '\n';
    }
    std::FILE* f = std::fopen(path, "wb");
    if (!f) { std::fprintf(stderr, "cannot open output file: %s\n", path); std::exit(1); }
    std::fwrite(out, 1, (size_t)(o - out), f);
    std::fclose(f);
}

// 4-ary min-heap of (key, vertex) with lazy deletion (callers skip stale pops).
template <class K>
struct Heap4 {
    struct Item { K d; int32_t v; };
    vector<Item> a;
    inline void clear() { a.clear(); }
    inline bool empty() const { return a.empty(); }
    inline K top_key() const { return a[0].d; }
    inline void push(K d, int32_t v) {
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
        size_t n = a.size();
        if (!n) return best;
        size_t i = 0;
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
// =================================================================== /common
