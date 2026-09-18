D=/root/v2_sfx/stressw; mkdir -p $D; cd $D
W=/mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/v2
cp $W/w64.cpp $W/common.inc $W/sfx/stressw/ref.cpp .
g++ -O2 -std=c++17 -o sw w64.cpp && g++ -O2 -std=c++17 -o ref ref.cpp || exit 1
fail=0
for i in $(seq $1 $2); do
  info=$(python3 $W/sfx/stressw/gen.py $i g)
  ./ref g.graph g.queries r.out
  for cfg in "W_STOPC=0" "W_STOPC=5" "W_STOPC=50" "W_STOPC=1000000" "W_STOPC=5 W_LMK=0" "W_STOPC=20 W_A16=1" "W_STOPC=3 W_LAZY=1" "W_STOPC=30 W_LAZY=1 W_LMK=8" "W_STOPC=10 W_ELMK=4 W_ELMK_AT=1000000" "W_STOPC=10 W_CONSET=3"; do
    env $cfg ./sw g.graph g.queries s.out 2>/dev/null
    if ! cmp -s r.out s.out; then echo "FAIL seed=$i [$cfg] ($info)"; fail=$((fail+1)); fi
  done
done
echo "done, failures=$fail"
