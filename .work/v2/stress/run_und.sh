# usage: bash run_und.sh <module.cpp> <count>
set -u
W=/mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/v2
D=/root/stress_und; mkdir -p $D; cd $D
cp $W/$1 m.cpp; cp $W/common.inc .; cp $W/stress/gen_und.py .
g++ -O2 -march=native -std=c++17 -o m m.cpp 2>&1 | grep -E "error" | head
g++ -O2 -std=c++17 -o ref /root/spc/dijkstra_foundation.cpp
fail=0
for i in $(seq 1 $2); do
  python3 gen_und.py $i t
  ./ref t.graph t.queries r.out
  for cfg in "" "STOPC=0" "STOPC=5,LMK=2" "LMK=0,ASTAR=0,MORTON=0" "STOPC=100000"; do
    env $(echo $cfg | tr ',' ' ') ./m t.graph t.queries m.out 2>/dev/null
    if ! cmp -s r.out m.out; then echo "FAIL seed=$i cfg=$cfg"; fail=$((fail+1)); fi
  done
done
echo "failures: $fail"
