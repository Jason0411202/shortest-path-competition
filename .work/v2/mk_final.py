"""Builds the single-file solver from common.inc and the family modules.

usage: python mk_final.py out.cpp module1.cpp [module2.cpp ...]
Each module: optional leading #pragma lines, `#include "common.inc"`, code in a
namespace exposing applicable()/solve(), and a trailing `int main(...)` that is
dropped here.
"""
import re
import sys

HERE = __file__.rsplit('/', 1)[0] if '/' in __file__ else '.'
HERE = HERE.replace('\\', '/')

HEADER = r'''// solver.cpp
//
// Derived from dijkstra_foundation.cpp for the Shortest-Path Competition.
// Same command line, same input formats, same output as the foundation:
//
//     ./solver <graph_file> <query_file> <output_file>
//
// The foundation's file-format documentation still applies verbatim; see
// dijkstra_foundation.cpp / README.md.  The foundation answers every query
// with its own Dijkstra search.  This file keeps its read_graph /
// run_queries / main structure, reads both input files in one go, and picks a
// preprocessing scheme from the shape of the graph it was given:
//
//   * undirected torus lattices (and any other undirected graph): contraction
//     hierarchies (Geisberger et al. 2008) with landmark-bounded, goal-directed
//     witness searches; the last few thousand vertices form a core whose
//     all-pairs distances are tabulated; queries are bidirectional upward
//     searches that meet either below the core or through the table;
//   * directed graphs with coordinates (road networks): see namespace road;
//   * directed graphs without coordinates (power-law graphs): see namespace sf.
//
// Every scheme is exact on any input of its kind; structure detection only
// selects faster code paths.  Single-threaded, C++17, standard library only.
'''


def module_body(path):
    src = open(path, encoding='utf-8').read()
    lines = src.split('\n')
    out = []
    started = False
    for ln in lines:
        if not started:
            if ln.startswith('#pragma GCC optimize') or ln.startswith('#pragma GCC target'):
                continue
            if ln.strip() == '#include "common.inc"':
                started = True
                continue
            if ln.strip() == '' or ln.startswith('//'):
                out.append(ln)
                continue
            started = True
        out.append(ln)
    body = '\n'.join(out)
    m = re.search(r'\nint main\s*\(', body)
    if m:
        body = body[:m.start()] + '\n'
    return body


def main():
    outp = sys.argv[1]
    mods = sys.argv[2:]
    common = open(HERE + '/common.inc', encoding='utf-8').read()
    parts = [HEADER, common, '']
    parts.append('#pragma GCC push_options')
    parts.append('#pragma GCC optimize("O3")')
    parts.append('#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt,fma")')
    for m in mods:
        b = module_body(m)
        # modules that manage their own push/pop keep them (nesting is fine)
        parts.append(b)
    parts.append('#pragma GCC pop_options')
    names = [re.search(r'\nnamespace (\w+) \{', module_body(m)).group(1) for m in mods]
    disp = []
    for n in names:
        disp.append(f'    if ({n}::applicable()) {{ {n}::solve(); return; }}')
    parts.append('''
// Picks the scheme for the graph at hand (first module that applies).
static void run_queries_all() {
''' + '\n'.join(disp) + '''
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: %s <graph_file> <query_file> <output_file>\\n", argv[0]);
        return 1;
    }
    g_debug = std::getenv("SPC_DEBUG") != nullptr;
    g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]);
    read_queries(argv[2]);
    tlog("read");
    run_queries_all();
    write_answers(argv[3]);
    tlog("write");
    // Nothing left to do: skip the destructors of the large containers.
    std::fflush(nullptr);
    std::_Exit(0);
}
''')
    open(outp, 'w', encoding='utf-8', newline='\n').write('\n'.join(parts))


if __name__ == '__main__':
    main()
