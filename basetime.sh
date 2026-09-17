#!/usr/bin/env bash
set -euo pipefail
cd /root/spc
g++ -O2 -std=c++17 -o foundation dijkstra_foundation.cpp
tier="${1:-dev}"
for n in road2d lattice3d local2d scalefree hugeq wide64; do
  /usr/bin/time -f "  foundation ${n}_${tier}  %e s  %M KB" \
    ./foundation "instances/${n}_${tier}.graph" "instances/${n}_${tier}.queries" /tmp/f.out
done
