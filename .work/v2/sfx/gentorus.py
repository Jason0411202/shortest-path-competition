import random, math, sys
n = int(sys.argv[1]); out = sys.argv[2]; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
r = random.Random(seed)
V = n * n
lines = [f"{V} {2*V} 0\n"]
la, lb = math.log(1e8), math.log(1e9)
for i in range(V):
    for st in (1, n):
        c = (i // st) % n
        j = i + st if c + 1 < n else i - c * st
        w = int(math.exp(r.uniform(la, lb)))
        lines.append(f"{i} {j} {w}\n")
open(out + ".graph", "w").writelines(lines)
Q = 2000
with open(out + ".queries", "w") as f:
    f.write(f"{Q}\n")
    for _ in range(Q): f.write(f"{r.randrange(V)} {r.randrange(V)}\n")
