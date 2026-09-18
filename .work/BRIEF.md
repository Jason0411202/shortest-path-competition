# Shared brief for solver work packages (2026-09-19)

## The task
`solver.cpp` (single C++17 file, std library only, single thread, derived from
`dijkstra_foundation.cpp`) answers shortest-path queries. Score per instance =
T_base / T_solver (wall clock, best of 3), overall = geometric mean of six.
Rules (README.md): same CLI `./solver <graph> <queries> <output>`, exact same
output as the foundation (`-1` for unreachable), no threads, < 4 GB, no
precomputed answers, only standard headers. Final grading re-runs solver.cpp on
fresh instances with the same generator parameters, so exploit graph
*structure* freely but never the exact bytes of the released files.

## Instances (large tier, scored) and current honest timings of solver.cpp
| instance | V | E | Q | type | query mix | T_base | now | phases now |
|---|---|---|---|---|---|---|---|---|
| road2d    | 300k | 1.63M | 200k | directed, coords | Dijkstra-rank 2^6..2^15 | 1013 s | 15.5 s | contract 9.6, hub labels 6.3, queries 0.15 |
| lattice3d | 300k | 902k | 30k | undirected 3-D torus 67^3, loguniform w 1..1e6 | uniform pairs | 1526 s | 30–45 s | contract 35, CH queries 8 |
| local2d   | 300k | 601k | 300k | undirected 2-D torus 548^2, uniform w 1..1000 | 90% rank 2^4..2^12, 10% uniform | 596 s | 10 s | contract 4.9, CH queries 9.8 |
| scalefree | 300k | 873k | 50k | directed power-law (gamma 2.1, m 3), loguniform w 1..1e4, 3% sinks | 70% uniform, 30% target = hub | 1378 s | 15.3 s | bidir Dijkstra queries 21 |
| hugeq     | 50k  | 271k | 600k | directed road, coords | uniform pairs | 1764 s | 1.3 s | contract 1.28, labels 0.58, queries 0.42 |
| wide64    | 2.0M | 4.0M | 30k  | undirected 2-D torus 1414^2, loguniform w 1e8..1e9, 64-bit dists | rank 2^8..2^15 | 882 s | 41.5 s | contract ~35+, queries rest |

"rank 2^k" = the target is the 2^k-th vertex settled by Dijkstra from the
source, i.e. these queries are local. Sources are (almost) all distinct.
Lattice vertex ids are row-major grid coordinates (edges (i,i+1) and (i,i+side)
etc., torus wrap-around), so grid coordinates can be inferred from the ids.
Dev-tier instances (`*_dev`, ~50k vertices, 10k queries) are in `instances/` for
correctness checks; the family/params match the large ones.

## Code map of solver.cpp (906 lines)
- 51–79 Scanner (fread whole file, hand parser)
- 108–136 read_graph: picks method. directed && !coords -> USE_BIDIR;
  directed -> CH + hub labels (USE_HL); undirected -> CH with bidirectional
  upward search (ch_query). Env overrides: USE_BIDIR, USE_HL, USE_PLL,
  SIM_SETTLE, CON_SETTLE, PRUNE_SETTLE, LEVEL_COEF, HL_PRUNE, SPC_DEBUG=1
  prints phase times to stderr.
- 140–176 MinHeap (4-ary); 199–283 witness search + prune_redundant_arcs
- 287–441 contraction: lazy priority queue, priority() = level*1000 +
  edge quotient + hop quotient (simulate() estimates shortcuts), find_shortcuts
  with bounded witness Dijkstra; dynamic graph g = vector<vector<DArc>>
- 449–530 build_query_graph (CSR by rank) and ch_query (bidir, stall-on-demand)
- 540–639 hub labels from CH; 649–743 pruned landmark labeling (PLL, unused)
- 752–842 bidirectional Dijkstra on CSR (scalefree path)
- 848–906 run_queries (formats output, `std::_Exit(0)` at the end) and main

## Environment and how to build / measure (all in WSL)
Run WSL commands from Windows with:  `wsl -e bash -lc '<command>'`
Data lives in `/root/spc/instances/` (large + dev, with `.answers` keys).
Each work package uses its OWN scratch dir `/root/spc_<TAG>` so agents do not
collide. Helper (from the Windows checkout):

    SRC=<file in .work> TAG=<tag> bash /mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/bench.sh "lattice3d:" "local2d:ENV=1"

It copies `.work/<SRC>` to `/root/spc_<TAG>/solver.cpp`, builds with
`g++ -O3 -march=native -std=c++17`, runs each listed large instance with
SPC_DEBUG=1, prints phase times + wall + RSS, and OK/WRONG against the key.
Env vars after the colon are passed to the run. Dev-tier check of all six:

    bash /mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/devcheck.sh <file in .work>

Machine: i9-13900H under WSL2, g++ 13. Other agents benchmark concurrently, so
expect ±15% noise; compare before/after within your own runs and repeat when
close. Do NOT run grade.py on the large tier (it takes ~1 h and needs a quiet
machine); the integrator does the final scoring. Never edit result.json,
`/root/spc/solver.cpp`, or the top-level `solver.cpp`; never submit anything.

## Deliverable per work package
1. Your variant as `.work/<name>.cpp` (a modified copy of `solver.cpp`; keep
   the other code paths working — devcheck.sh must print OK for all six).
2. A final report: what you changed, phase timings before/after on your
   instances (large tier), and what you tried that did not help.
Keep the code single-file, std-only, single-threaded, deterministic. Avoid
`#include <immintrin.h>` / `<sys/mman.h>` (not standard). `#pragma GCC
optimize`/`target` and GCC vector extensions are acceptable.
