# usage: abf.sh <reps> <core> <bin> "cfgA" "cfgB" ... : full-run phase timing, alternating, min over reps
R=$1; CORE=$2; BIN=$3; shift 3
declare -A BC BT
for i in $(seq $R); do
  n=0
  for cfg in "$@"; do
    out=$(env $cfg SPC_DEBUG=1 taskset -c $CORE $BIN /root/spc/instances/wide64_large.graph /root/spc/instances/wide64_large.queries /tmp/abf 2>&1)
    ok=OK; cmp -s /tmp/abf /root/spc/instances/wide64_large.answers || ok=WRONG
    c=$(echo "$out" | grep -E 'contract$' | sed -E 's/\[ *([0-9.]+) s\].*/\1/')
    t=$(echo "$out" | grep -E 'queries$' | tail -1 | sed -E 's/\[ *([0-9.]+) s\].*/\1/')
    BC[$n]=$(echo "$c ${BC[$n]:-999}" | awk '{print ($1<$2)?$1:$2}')
    BT[$n]=$(echo "$t ${BT[$n]:-999}" | awk '{print ($1<$2)?$1:$2}')
    [ $ok = WRONG ] && BT[$n]=WRONG
    n=$((n+1))
  done
done
n=0; for cfg in "$@"; do echo "[$cfg] contract-end=${BC[$n]} total=${BT[$n]}"; n=$((n+1)); done
