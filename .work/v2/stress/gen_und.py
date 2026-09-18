import random, sys
seed = int(sys.argv[1]); out = sys.argv[2]
r = random.Random(seed)
kind = seed % 4
V = r.choice([1, 2, 3, 5, 10, 40, 200, 1500])
edges = []
if kind == 0:   # random sparse, maybe disconnected
    for _ in range(r.randint(0, 3 * V)):
        edges.append((r.randrange(V), r.randrange(V), r.randint(1, 1000)))
elif kind == 1:  # torus lattice exact generator order
    n = r.choice([3, 4, 5, 7, 12]); d = r.choice([2, 3]); V = n ** d
    st = [n ** k for k in range(d)]
    for i in range(V):
        for k in range(d):
            c = (i // st[k]) % n
            j = i + st[k] if c + 1 < n else i - c * st[k]
            edges.append((i, j, int(round(10 ** r.uniform(0, 6)))))
elif kind == 2:  # big weights (64-bit distances), path-like
    for i in range(V - 1):
        edges.append((i, i + 1, r.randint(10 ** 8, 2 * 10 ** 9)))
    for _ in range(V):
        edges.append((r.randrange(V), r.randrange(V), r.randint(10 ** 8, 2 * 10 ** 9)))
else:            # dense-ish with duplicates and self loops
    for _ in range(r.randint(0, 8 * V)):
        u = r.randrange(V); v = u if r.random() < 0.05 else r.randrange(V)
        edges.append((u, v, r.randint(1, 20)))
edges = [e for e in edges if e[0] != e[1]] if r.random() < 0.5 else edges
with open(out + '.graph', 'w') as f:
    f.write(f"{V} {len(edges)} 0\n")
    for u, v, w in edges:
        f.write(f"{u} {v} {w}\n")
Q = r.randint(1, 300)
with open(out + '.queries', 'w') as f:
    f.write(f"{Q}\n")
    for _ in range(Q):
        s = r.randrange(V); t = s if r.random() < 0.1 else r.randrange(V)
        f.write(f"{s} {t}\n")
