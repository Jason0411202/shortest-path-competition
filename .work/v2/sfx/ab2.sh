# usage: ab2.sh <reps> <graph-base> binA binB ...   (alternating; prints min total wall & phase)
R=$1; G=$2; shift 2
declare -A BEST
for i in $(seq $R); do
  for b in "$@"; do
    s=$(date +%s.%N)
    taskset -c 8 $b $G.graph $G.queries /tmp/ab2x 2>/dev/null
    e=$(date +%s.%N)
    t=$(echo "$e - $s" | bc)
    cmp -s /tmp/ab2x $G.answers || t="WRONG"
    BEST[$b]=$(echo "$t ${BEST[$b]:-999}" | awk '{print ($1<$2)?$1:$2}')
  done
done
for b in "$@"; do echo "$b wall(min of $R)=${BEST[$b]}"; done
