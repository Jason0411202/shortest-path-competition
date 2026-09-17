# usage: SRC=file.cpp TAG=x bash bench.sh "inst:ENV=1,ENV2=2" ...
set -u
D=/root/spc_${TAG}; mkdir -p $D && cd $D && ln -sfn /root/spc/instances instances
cp /mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/$SRC solver.cpp
g++ -O3 -march=native -std=c++17 -o solver solver.cpp 2>&1 | grep -E "error" | head -5
for spec in "$@"; do
  n=${spec%%:*}; cfg=${spec#*:}
  r=$(env $(echo $cfg | tr ',' ' ') SPC_DEBUG=1 /usr/bin/time -f "wall=%e rss=%MKB" timeout ${TO:-600} ./solver instances/${n}_large.graph instances/${n}_large.queries /tmp/o_$TAG 2>&1 | grep -vE "^  |progress" | tr '\n' ' ')
  cmp -s /tmp/o_$TAG instances/${n}_large.answers && ok=OK || ok=WRONG
  echo "$n [$cfg] $ok :: $r"
done
