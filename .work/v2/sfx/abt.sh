# usage: abt.sh <reps> <graphbase> "binA|cfgA" "binB|cfgB" ... : min contraction time (build->contract) over reps, alternating
R=$1; G=$2; shift 2
declare -A B
for i in $(seq $R); do
  n=0
  for spec in "$@"; do
    bin=${spec%%|*}; cfg=${spec#*|}; [ "$cfg" = "$spec" ] && cfg=""
    out=$(env $cfg W_EXIT=1 SPC_DEBUG=1 $bin $G.graph $G.queries /tmp/abt 2>&1)
    b=$(echo "$out" | grep -E 'build$' | sed -E 's/\[ *([0-9.]+) s\].*/\1/')
    c=$(echo "$out" | grep -E 'contract$' | sed -E 's/\[ *([0-9.]+) s\].*/\1/')
    t=$(echo "$c $b" | awk '{printf "%.3f", $1-$2}')
    B[$n]=$(echo "$t ${B[$n]:-999}" | awk '{print ($1<$2)?$1:$2}')
    n=$((n+1))
  done
done
n=0; for spec in "$@"; do echo "[$spec] contract(min)=${B[$n]}"; n=$((n+1)); done
