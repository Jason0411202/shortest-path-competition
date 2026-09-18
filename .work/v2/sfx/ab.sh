# usage: ab.sh <reps> <queries-file|large> "cfgA" "cfgB" ...  alternates configs each rep; prints min phases
R=$1; QF=$2; shift 2
[ "$QF" = large ] && QF=/root/spc/instances/scalefree_large.queries
declare -A BT BQ BTOT
for i in $(seq $R); do
  n=0
  for cfg in "$@"; do
    out=$(env $cfg SPC_DEBUG=1 taskset -c 8 /root/v2_sf/s1 /root/spc/instances/scalefree_large.graph $QF /tmp/abx 2>&1)
    tq=$(echo "$out" | grep -E 'queries$' | tail -1 | sed -E 's/\[ *([0-9.]+) s\].*/\1/')
    tb=$(echo "$out" | grep -E '(hub tables|graph built)$' | tail -1 | sed -E 's/\[ *([0-9.]+) s\].*/\1/')
    q=$(echo "$tq $tb" | awk '{printf "%.3f", $1-$2}')
    k="$n"
    BT[$k]=$(echo "$tb ${BT[$k]:-999}" | awk '{print ($1<$2)?$1:$2}')
    BQ[$k]=$(echo "$q ${BQ[$k]:-999}" | awk '{print ($1<$2)?$1:$2}')
    BTOT[$k]=$(echo "$tq ${BTOT[$k]:-999}" | awk '{print ($1<$2)?$1:$2}')
    n=$((n+1))
  done
done
n=0
for cfg in "$@"; do echo "[$cfg] pre=${BT[$n]} query=${BQ[$n]} total=${BTOT[$n]}"; n=$((n+1)); done
