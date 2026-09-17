cd /root/spc && cp /mnt/c/Home/programmingFiles/Github/shortest-path-competition/solver.cpp . && g++ -O3 -march=native -std=c++17 -o solver solver.cpp || exit 1
for n in $INST; do for cfg in $CFGS; do
  export SIM_SETTLE=${cfg%/*} CON_SETTLE=${cfg#*/}
  r=$(SPC_DEBUG=1 ./solver instances/${n}.graph instances/${n}.queries /tmp/o 2>&1 | grep -E "\] (read|contract|queries|labels)|upward|label entries|10[0-9]*/[0-9]* deg" | tr '\n' ' ')
  cmp -s /tmp/o instances/${n}.answers && ok=OK || ok=WRONG
  echo "$n sim=$SIM_SETTLE con=$CON_SETTLE $ok :: $r"
done; done
