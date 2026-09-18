# usage: rb.sh <reps> "ENV=.. ENV2=.." ...   (binary /root/v2_sf/s1, large instance)
R=$1; shift
for cfg in "$@"; do
  best=999; bt=999; ok=OK
  for i in $(seq $R); do
    out=$(env $cfg SPC_DEBUG=1 taskset -c 8 /root/v2_sf/s1 /root/spc/instances/scalefree_large.graph /root/spc/instances/scalefree_large.queries /tmp/rbx 2>&1)
    cmp -s /tmp/rbx /root/spc/instances/scalefree_large.answers || ok=WRONG
    tot=$(echo "$out" | grep -E 'queries$' | tail -1 | sed -E 's/\[ *([0-9.]+) s\].*/\1/')
    tb=$(echo "$out" | grep -E 'hub tables$' | sed -E 's/\[ *([0-9.]+) s\].*/\1/'); [ -z "$tb" ] && tb=0
    best=$(echo "$tot $best" | awk '{print ($1<$2)?$1:$2}'); bt=$(echo "$tb $bt" | awk '{print ($1<$2)?$1:$2}')
  done
  st=$(echo "$out" | grep -oE 'settle/search [0-9.]+')
  echo "$cfg :: $ok total(min)=$best tables_end(min)=$bt $st"
done
