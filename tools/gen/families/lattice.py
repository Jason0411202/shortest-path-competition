"""d-dimensional torus lattices.

One family covers three instances; what differs between them is the dimension,
the weight range and the query mix, not the construction.

``lattice3d`` (d=3, n=126) has degree 6 and, at equal V, far larger vertex
separators than a 2D lattice of the same size: O(n^2) rather than O(n), which
at V = 2M is a treewidth of roughly 15,876 against 1,414.  ``local2d`` and
``wide64`` (d=2, n=1414) are planar tori with degree 4.

``chord_frac`` adds uniformly random long-range edges on top of the lattice.
It is a continuous version of the same dial: measured on a 45x45 grid, chords
at 0 / 5 / 15 / 30 percent of V move a hierarchy's shortcut count from 1.41 to
1.80, 2.68 and 3.57 shortcuts per edge.  Note that small rates do essentially
nothing -- 0.5% measured 1.40, indistinguishable from none.

``wide64`` carries weights in [1e8, 1e9], which puts shortest-path distances
6.6x past 2^32 even on the dev tier.  The one-decade range is deliberate: a
wider [1e6, 1e9] range reads as more interesting but lets shortest paths route
almost entirely through the cheap decade, dragging typical path cost down by a
factor of ~25 so that the dev tier tops out at 1.75e9 -- *below* 2^32.  An
overflow that lands marginally is worse than none at all: it would show up on a
handful of queries out of half a million, which ranks luck and cannot be
debugged against a hidden set.
"""

from ..graph import Graph
from ..weights import make as make_weight


def build(rng, params):
    d = int(params.get("dim", 2))
    n = int(params["n"])
    torus = bool(params.get("torus", True))
    chord_frac = float(params.get("chord_frac", 0.0))
    wspec = params.get("weights", "uniform:1:1000")

    if torus and n < 3:
        raise ValueError("a torus needs n >= 3, otherwise +1 and its wraparound "
                         "are the same edge")
    V = n ** d
    g = Graph(V, directed=False)
    weight = make_weight(wspec, None)
    rw = rng["weights"]

    # strides[k] is the index delta for a step of +1 along axis k
    strides = [n ** k for k in range(d)]

    def coord(i, k):
        return (i // strides[k]) % n

    for i in range(V):
        for k in range(d):
            c = coord(i, k)
            if c + 1 < n:
                j = i + strides[k]
            elif torus:
                j = i - c * strides[k]
            else:
                continue
            g.add(i, j, weight(rw))

    if chord_frac > 0.0:
        rt = rng["topology"]
        for _ in range(int(chord_frac * V)):
            a = rt.randrange(V)
            b = rt.randrange(V)
            if a != b:
                g.add(a, b, weight(rw))

    meta = {"dim": d, "side": n, "torus": torus, "weights": wspec,
            "chord_frac": chord_frac}
    return g, meta
