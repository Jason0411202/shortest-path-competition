# usage: bash b.sh <file.cpp in .work/v2> <tier: dev|large> inst[:ENV=1,ENV2=2] ...
# builds /root/v2_<file>/s from .work/v2 (common.inc included), runs each instance, checks answers
set -u
W=/mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/v2
F=$1; T=$2; shift 2
D=/root/v2_${F%.cpp}; mkdir -p $D
cp $W/$F $W/common.inc $D/
cd $D && g++ -O3 -march=native -std=c++17 -Wall -o s $F 2>&1 | grep -E "error|warning" | head -20
for spec in "$@"; do
  n=${spec%%:*}; cfg=${spec#*:}; [ "$cfg" = "$spec" ] && cfg=""
  r=$(env $(echo $cfg | tr ',' ' ') SPC_DEBUG=1 /usr/bin/time -f "wall=%e rss=%MKB" timeout ${TO:-600} ./s /root/spc/instances/${n}_$T.graph /root/spc/instances/${n}_$T.queries $D/out_$n 2>&1 | tr '\n' ' ')
  cmp -s $D/out_$n /root/spc/instances/${n}_$T.answers && ok=OK || ok=WRONG
  echo "$n [$cfg] $ok :: $r"
done
