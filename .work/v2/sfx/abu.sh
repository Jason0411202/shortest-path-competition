# usage: abu.sh <reps> <bin> "cfgA" "cfgB" ... : user+sys CPU seconds (min over reps), alternating
R=$1; BIN=$2; shift 2
declare -A B
for i in $(seq $R); do
  n=0
  for cfg in "$@"; do
    t=$( { env $cfg /usr/bin/time -f "%U %S" $BIN /root/spc/instances/wide64_large.graph /root/spc/instances/wide64_large.queries /tmp/abu > /dev/null 2>/tmp/abu.err; } ; tail -1 /tmp/abu.err | awk '{print $1+$2}')
    cmp -s /tmp/abu /root/spc/instances/wide64_large.answers || t=999
    B[$n]=$(echo "$t ${B[$n]:-999}" | awk '{print ($1<$2)?$1:$2}')
    n=$((n+1))
  done
done
n=0; for cfg in "$@"; do echo "[$cfg] cpu=${B[$n]}"; n=$((n+1)); done
