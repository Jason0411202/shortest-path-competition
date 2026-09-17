#!/usr/bin/env bash
# Helper: sync solver.cpp from the Windows checkout into the WSL working copy,
# build it, and run whatever was asked for.
#
#   run.sh build            -- just compile
#   run.sh check            -- compile + verify dev answers (fast, no timing)
#   run.sh dev              -- compile + grade.py on the dev tier
#   run.sh large            -- compile + grade.py on the scored tier (--json)
#   run.sh time <name>      -- compile + time one large instance directly
set -euo pipefail

WIN=/mnt/c/Home/programmingFiles/Github/shortest-path-competition
WORK=/root/spc
CXXFLAGS="-O3 -march=native -std=c++17 -Wall -Wextra"

cd "$WORK"
cp "$WIN/solver.cpp" "$WORK/solver.cpp"

cmd="${1:-check}"

build() {
  echo ">>> building (${CXXFLAGS})"
  /usr/bin/time -f "    compile: %e s" g++ $CXXFLAGS -o solver solver.cpp
}

case "$cmd" in
  build) build ;;
  check)
    build
    for n in road2d lattice3d local2d scalefree hugeq wide64; do
      g="instances/${n}_dev.graph"; q="instances/${n}_dev.queries"; a="instances/${n}_dev.answers"
      /usr/bin/time -f "%e" -o /tmp/t.$$ ./solver "$g" "$q" /tmp/out.$$
      tt=$(cat /tmp/t.$$)
      if diff -q <(sed 's/[[:space:]]*$//' /tmp/out.$$) <(sed 's/[[:space:]]*$//' "$a") >/dev/null; then
        printf "  ok    %-12s %8s s\n" "$n" "$tt"
      else
        printf "  WRONG %-12s %8s s\n" "$n" "$tt"
        diff <(sed 's/[[:space:]]*$//' /tmp/out.$$) <(sed 's/[[:space:]]*$//' "$a") | head -5
      fi
    done
    rm -f /tmp/t.$$ /tmp/out.$$
    ;;
  checklarge)
    build
    for n in road2d lattice3d local2d scalefree hugeq wide64; do
      g="instances/${n}_large.graph"; q="instances/${n}_large.queries"; a="instances/${n}_large.answers"
      /usr/bin/time -f "%e %M" -o /tmp/t.$$ ./solver "$g" "$q" /tmp/out.$$
      tt=$(cat /tmp/t.$$)
      if diff -q <(sed 's/[[:space:]]*$//' /tmp/out.$$) <(sed 's/[[:space:]]*$//' "$a") >/dev/null; then
        printf "  ok    %-12s %s (s, KB)\n" "$n" "$tt"
      else
        printf "  WRONG %-12s %s (s, KB)\n" "$n" "$tt"
        diff <(sed 's/[[:space:]]*$//' /tmp/out.$$) <(sed 's/[[:space:]]*$//' "$a") | head -5
      fi
    done
    rm -f /tmp/t.$$ /tmp/out.$$
    ;;
  one)
    build
    n="$2"; tier="${3:-large}"
    g="instances/${n}_${tier}.graph"; q="instances/${n}_${tier}.queries"; a="instances/${n}_${tier}.answers"
    /usr/bin/time -f "    time: %e s   maxrss: %M KB" ./solver "$g" "$q" /tmp/out.$$
    if diff -q <(sed 's/[[:space:]]*$//' /tmp/out.$$) <(sed 's/[[:space:]]*$//' "$a") >/dev/null; then
      echo "    ok"
    else
      echo "    WRONG"; diff <(sed 's/[[:space:]]*$//' /tmp/out.$$) <(sed 's/[[:space:]]*$//' "$a") | head -10
    fi
    rm -f /tmp/out.$$
    ;;
  dev)
    build
    python3 grade.py --solver ./solver --instances instances_dev.txt
    ;;
  large)
    build
    python3 grade.py --solver ./solver --instances instances.txt --json result.json
    cp result.json "$WIN/result.json"
    ;;
  *) echo "unknown: $cmd" >&2; exit 2 ;;
esac
