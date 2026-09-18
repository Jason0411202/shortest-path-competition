import random, sys
seed = int(sys.argv[1]); out = sys.argv[2]
r = random.Random(seed)
kind = r.choice(['torus', 'torus', 'rand', 'grid_plus', 'tiny', 'disc'])
lo = r.choice([0, 1, 1, 100000000])
hi = r.choice([1, 1000, 10**9, 4294967295])
if lo > hi: lo, hi = hi, lo
arcs = []
if kind == 'torus':
    n = r.randint(3, 40); V = n * n
    for i in range(V):
        for k, st in enumerate([1, n]):
            c = (i // st) % n
            j = i + st if c + 1 < n else i - c * st
            arcs.append((i, j, r.randint(lo, hi)))
elif kind == 'grid_plus':
    n = r.randint(3, 30); V = n * n
    for i in range(V):
        for k, st in enumerate([1, n]):
            c = (i // st) % n
            j = i + st if c + 1 < n else i - c * st
            arcs.append((i, j, r.randint(lo, hi)))
    for _ in range(r.randint(1, V)):
        a = r.randrange(V); b = r.randrange(V); arcs.append((a, b, r.randint(lo, hi)))
elif kind == 'rand':
    V = r.randint(2, 800)
    for _ in range(r.randint(0, 4 * V)): arcs.append((r.randrange(V), r.randrange(V), r.randint(lo, hi)))
elif kind == 'disc':
    V = r.randint(10, 600)
    for _ in range(r.randint(V // 2, 2 * V)):
        a = r.randrange(V // 2); b = r.randrange(V // 2)
        if r.random() < 0.5: a += V // 2; b += V // 2
        arcs.append((a, b, r.randint(lo, hi)))
else:
    V = r.randint(1, 6)
    for _ in range(r.randint(0, 10)): arcs.append((r.randrange(V), r.randrange(V), r.randint(lo, hi)))
with open(out + '.graph', 'w') as f:
    f.write(f"{V} {len(arcs)} 0\n")
    for a in arcs: f.write(f"{a[0]} {a[1]} {a[2]}\n")
Q = r.randint(1, 300)
with open(out + '.queries', 'w') as f:
    f.write(f"{Q}\n")
    for _ in range(Q):
        s = r.randrange(V); t = s if r.random() < 0.05 else r.randrange(V)
        f.write(f"{s} {t}\n")
print(kind, V, len(arcs), lo, hi)
