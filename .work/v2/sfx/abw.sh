# usage: abw.sh <reps> <core> "cfgA" "cfgB" ... : contraction-only timing of /root/v2_w64b/s (large), alternating
R=$1; CORE=$2; shift 2
declare -A B
for i in $(seq $R); do
  n=0
  for cfg in "$@"; do
    t=$(env $cfg W_EXIT=1 SPC_DEBUG=1 taskset -c $CORE /root/v2_w64b/s /root/spc/instances/wide64_large.graph /root/spc/instances/wide64_large.queries /tmp/abw 2>&1 | grep -E 'contract$' | sed -E 's/\[ *([0-9.]+) s\].*/\1/')
    B[$n]=$(echo "$t ${B[$n]:-999}" | awk '{print ($1<$2)?$1:$2}')
    n=$((n+1))
  done
done
n=0; for cfg in "$@"; do echo "[$cfg] contract-end(min)=${B[$n]}"; n=$((n+1)); done
