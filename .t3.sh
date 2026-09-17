# usage: SRC=hy.cpp INST=... ENVS="A=1,B=2 ..." bash .t3.sh
D=/root/spc_${SRC%.cpp}; mkdir -p $D && cd $D && ln -sfn /root/spc/instances instances
cp /mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/$SRC solver.cpp && g++ -O3 -march=native -std=c++17 -Wall -o solver solver.cpp 2>&1 | grep -E "error|warning" | head -5
for n in $INST; do for cfg in ${ENVS:-NONE=1}; do
  r=$(env $(echo $cfg | tr ',' ' ') SPC_DEBUG=1 /usr/bin/time -f "wall=%e rss=%MKB" timeout ${TO:-300} ./solver instances/${n}.graph instances/${n}.queries /tmp/o_$SRC 2>&1 | grep -vE "^  " | tr '\n' ' ')
  cmp -s /tmp/o_$SRC instances/${n}.answers && ok=OK || ok=WRONG
  echo "$n [$cfg] $ok :: $r"
done; done
