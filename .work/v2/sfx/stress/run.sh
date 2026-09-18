# usage: run.sh <n>  -- random graphs vs reference
D=/root/v2_sfx/stress; mkdir -p $D; cd $D
W=/mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/v2
cp $W/sf.cpp $W/common.inc . ; cp $W/sfx/stress/ref.cpp . ; sed -i 's#../../common.inc#common.inc#' ref.cpp
g++ -O2 -march=native -std=c++17 -o sfb sf.cpp && g++ -O2 -std=c++17 -o ref ref.cpp || exit 1
fail=0
for i in $(seq $1); do
  info=$(python3 $W/sfx/stress/gen.py $i g)
  ./ref g.graph g.queries r.out
  for k in 0 1 3 8 64; do
    for c in 0 1; do
      SF_K=$k SF_CHAIN=$c ./sfb g.graph g.queries s.out 2>/dev/null
      if ! cmp -s r.out s.out; then echo "FAIL seed=$i k=$k chain=$c ($info)"; fail=$((fail+1)); fi
    done
  done
done
echo "done, failures=$fail"
