# usage: INST=... ENVS="A=1,B=2 ..." bash .t2.sh   (binary built in /root/spc2)
mkdir -p /root/spc2 && cd /root/spc2 && ln -sfn /root/spc/instances instances
cp /mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/sp.cpp solver.cpp && g++ -O3 -march=native -std=c++17 -o solver solver.cpp || exit 1
for n in $INST; do for cfg in ${ENVS:-NONE=1}; do
  r=$(env $(echo $cfg | tr ',' ' ') SPC_DEBUG=1 /usr/bin/time -f "wall=%e rss=%MKB" ./solver instances/${n}.graph instances/${n}.queries /tmp/o2 2>&1 | grep -E "\] (read|contract|queries|labels)|upward|label entries|pll|bidir|wall=" | tr '\n' ' ')
  cmp -s /tmp/o2 instances/${n}.answers && ok=OK || ok=WRONG
  echo "$n [$cfg] $ok :: $r"
done; done
