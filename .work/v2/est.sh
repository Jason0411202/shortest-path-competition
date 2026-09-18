# usage: bash est.sh final.cpp  -- build with -O2, time each large instance (best of N), estimate score
set -u
W=/mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/v2
D=/root/est; mkdir -p $D; cd $D
cp $W/$1 f.cpp
g++ -O2 -std=c++17 -o f f.cpp 2>&1 | grep error | head
declare -A TB=( [road2d]=1012.84 [lattice3d]=1526.03 [local2d]=596.09 [scalefree]=1378.19 [hugeq]=1764 [wide64]=882 )
N=${N:-2}
logsum=0
for n in road2d lattice3d local2d scalefree hugeq wide64; do
  best=999
  for r in $(seq 1 $N); do
    /usr/bin/time -f "%e" -o t.txt ./f /root/spc/instances/${n}_large.graph /root/spc/instances/${n}_large.queries o.txt
    t=$(cat t.txt); best=$(python3 -c "print(min($best,$t))")
  done
  cmp -s o.txt /root/spc/instances/${n}_large.answers && ok=OK || ok=WRONG
  sp=$(python3 -c "print(round(${TB[$n]}/max($best,0.01),1))")
  echo "$n $ok best=$best speedup=$sp"
  logsum=$(python3 -c "import math;print($logsum+math.log(${TB[$n]}/max($best,0.01)))")
done
python3 -c "import math;print('geomean', round(math.exp($logsum/6),1))"
