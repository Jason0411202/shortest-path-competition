import os, collections
os.chdir('/root/spc/instances')
for n in ['road2d','lattice3d','local2d','scalefree','hugeq','wide64']:
    g = n + '_large.graph'
    q = n + '_large.queries'
    with open(g) as f:
        hdr = f.readline().split()
    with open(q) as f:
        Q = int(f.readline())
        srcs = collections.Counter()
        dsts = collections.Counter()
        self_pairs = 0
        for _ in range(Q):
            s,t = f.readline().split()
            srcs[s]+=1; dsts[t]+=1
            if s==t: self_pairs+=1
    print(f"{n:11s} V={hdr[0]:>8} E={hdr[1]:>8} FLAGS={hdr[2]} | Q={Q:>7} distinct_src={len(srcs):>7} maxdup={srcs.most_common(1)[0][1]:>4} distinct_dst={len(dsts):>7} self={self_pairs}")
