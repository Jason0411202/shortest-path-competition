"""The in-memory graph a family builds, plus the invariants every family must
satisfy before its instance is allowed to be written out.

Edges are stored as three parallel ``array``s rather than a list of tuples: at
V = 2,000,000 the tuple representation costs about 2 GB, the arrays about 60 MB.
"""

from array import array
from collections import deque

FLAG_COORDS = 1
FLAG_DIRECTED = 2


class Graph:
    def __init__(self, V, directed=False):
        self.V = V
        self.directed = directed
        self.eu = array("i")
        self.ev = array("i")
        self.ew = array("q")
        self.x = None  # array("q") of integer coordinates, or None
        self.y = None

    # -- construction ------------------------------------------------------

    def add(self, u, v, w):
        self.eu.append(u)
        self.ev.append(v)
        self.ew.append(w)

    @property
    def E(self):
        return len(self.eu)

    @property
    def flags(self):
        f = 0
        if self.x is not None:
            f |= FLAG_COORDS
        if self.directed:
            f |= FLAG_DIRECTED
        return f

    # -- derived views -----------------------------------------------------

    def _undirected_csr(self):
        V = self.V
        deg = array("i", bytes(4 * V))
        for i in range(self.E):
            deg[self.eu[i]] += 1
            deg[self.ev[i]] += 1
        off = array("q", bytes(8 * (V + 1)))
        s = 0
        for i in range(V):
            off[i] = s
            s += deg[i]
        off[V] = s
        pos = array("q", off[:V])
        tgt = array("i", bytes(4 * s))
        for i in range(self.E):
            u, v = self.eu[i], self.ev[i]
            tgt[pos[u]] = v
            pos[u] += 1
            tgt[pos[v]] = u
            pos[v] += 1
        return off, tgt

    def component_labels(self):
        """(labels, representatives) under the *undirected* reading of the arcs.

        For a directed graph this is weak connectivity -- an upper bound on
        reachability, which is all the connectivity postcondition needs.
        """
        V = self.V
        off, tgt = self._undirected_csr()
        labels = array("i", [-1]) * V
        reps = []
        for root in range(V):
            if labels[root] >= 0:
                continue
            lab = len(reps)
            reps.append(root)
            labels[root] = lab
            dq = deque((root,))
            while dq:
                a = dq.popleft()
                for j in range(off[a], off[a + 1]):
                    b = tgt[j]
                    if labels[b] < 0:
                        labels[b] = lab
                        dq.append(b)
        return labels, reps

    def undirected_components(self):
        """Component sizes, largest first."""
        labels, reps = self.component_labels()
        sizes = [0] * len(reps)
        for lab in labels:
            sizes[lab] += 1
        sizes.sort(reverse=True)
        return sizes

    # -- postconditions ----------------------------------------------------

    def check(self, expect_components=1):
        """Postconditions every generated graph must satisfy.

        A lattice or road generator can quietly leave the graph disconnected;
        the README promises one component, so check it instead of assuming.
        """
        V, E = self.V, self.E
        if V <= 0 or E <= 0:
            raise AssertionError("empty graph")
        for name, arr in (("u", self.eu), ("v", self.ev)):
            lo, hi = min(arr), max(arr)
            if lo < 0 or hi >= V:
                raise AssertionError(f"endpoint {name} out of range: [{lo},{hi}] vs V={V}")
        if min(self.ew) < 1:
            raise AssertionError(f"non-positive edge weight: {min(self.ew)}")
        for i in range(E):
            if self.eu[i] == self.ev[i]:
                raise AssertionError(f"self-loop at edge {i}")
        if self.x is not None and (len(self.x) != V or len(self.y) != V):
            raise AssertionError("coordinate block length != V")

        sizes = self.undirected_components()
        if expect_components is not None and len(sizes) != expect_components:
            raise AssertionError(
                f"expected {expect_components} weakly-connected component(s), "
                f"got {len(sizes)} (largest {sizes[0]}, smallest {sizes[-1]})"
            )
        return sizes


def ensure_weakly_connected(g, rng, weight):
    """Merge every stray component into the largest one.

    Returns the number of arcs added.  Random constructions -- the directed
    configuration model in particular -- reliably leave a handful of tiny
    components behind, while the README promises one.  Repairing here means the
    postcondition in ``Graph.check`` can be an assertion rather than a hope.
    """
    labels, reps = g.component_labels()
    if len(reps) <= 1:
        return 0
    sizes = [0] * len(reps)
    for lab in labels:
        sizes[lab] += 1
    main = max(range(len(reps)), key=lambda k: sizes[k])
    added = 0
    for lab, rep in enumerate(reps):
        if lab == main:
            continue
        anchor = reps[main]
        # both directions, so the repair cannot create a one-way trap
        g.add(anchor, rep, weight())
        if g.directed:
            g.add(rep, anchor, weight())
        added += 1
    return added
