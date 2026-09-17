# usage: bash devcheck.sh file.cpp  -- compile and check all dev instances
set -u
D=/root/devchk; mkdir -p $D && cd $D && ln -sfn /root/spc/instances instances
cp /mnt/c/Home/programmingFiles/Github/shortest-path-competition/.work/$1 s.cpp
g++ -O3 -march=native -std=c++17 -Wall -o s s.cpp 2>&1 | head -20
for n in road2d lattice3d local2d scalefree hugeq wide64; do
  /usr/bin/time -f "%e s %M KB" -o /tmp/dct ./s instances/${n}_dev.graph instances/${n}_dev.queries /tmp/dco
  cmp -s /tmp/dco instances/${n}_dev.answers && echo "$n OK $(cat /tmp/dct)" || echo "$n WRONG $(cat /tmp/dct)"
done
