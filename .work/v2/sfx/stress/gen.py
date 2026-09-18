import random, sys
seed = int(sys.argv[1]); out = sys.argv[2]
r = random.Random(seed)
kind = r.choice(['sf', 'rand', 'chain', 'cycle', 'tiny', 'big'])
V = {'tiny': r.randint(1, 8), 'big': r.randint(500, 3000)}.get(kind, r.randint(20, 400))
lo = r.choice([0, 1, 1]); maxw = r.choice([1, 3, 100, 10000, 60000, 2**20, 2**31 - 1])
arcs = []
if kind in ('sf', 'big'):
    wts = [r.paretovariate(1.1) for _ in range(V)]
    tot = sum(wts)
    for u in range(V):
        if r.random() < 0.05: continue
        for _ in range(3):
            v = r.choices(range(V), weights=wts)[0] if V < 2000 else int(min(V - 1, r.paretovariate(1.1)) ) if r.random() < 0.5 else r.randrange(V)
            arcs.append((u, v, r.randint(lo, maxw)))
elif kind == 'chain':
    for u in range(V):
        if r.random() < 0.3: arcs.append((u, (u + 1) % V, r.randint(lo, maxw)))
        for _ in range(r.randint(0, 2)): arcs.append((u, r.randrange(V), r.randint(lo, maxw)))
elif kind == 'cycle':
    for u in range(V): arcs.append((u, (u + 1) % V, r.randint(lo, maxw)))
    for _ in range(r.randint(0, V)): arcs.append((r.randrange(V), r.randrange(V), r.randint(lo, maxw)))
else:
    m = r.randint(0, 4 * V)
    for _ in range(m): arcs.append((r.randrange(V), r.randrange(V), r.randint(lo, maxw)))
with open(out + '.graph', 'w') as f:
    f.write(f"{V} {len(arcs)} 2\n")
    for a in arcs: f.write(f"{a[0]} {a[1]} {a[2]}\n")
Q = r.randint(1, 300)
with open(out + '.queries', 'w') as f:
    f.write(f"{Q}\n")
    for _ in range(Q):
        s = r.randrange(V); t = s if r.random() < 0.05 else r.randrange(V)
        f.write(f"{s} {t}\n")
print(kind, V, len(arcs), maxw)
