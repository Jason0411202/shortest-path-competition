"""ROAD: a perturbed lattice carrying a *travel-time* metric.

Topology is a jittered S x S lattice plus random diagonals.  The weight model
is the whole point of the family.

Weights are deliberately *not* ``w = round(1e6 * ||u-v||)``: that would make
the graph metric a scaled copy of the coordinate embedding, which hands out a
large advantage for nothing.

Here each edge belongs to a road class with a speed factor, and
``w = round(TIME_SCALE * ||u-v|| / speed * noise)``.  Shortest paths prefer the
fast classes while ``lambda`` collapses to the slowest one, so the Euclidean
potential goes slack by roughly ``max_speed * exp(clamp*sigma)``.  Measured
effect on the released instances: 7.30x down to 1.63x.

The road classes are laid out on lattice lines at fixed strides, so the graph
has a genuine hierarchy (local streets -> arterials -> highways) rather than
i.i.d. speeds.  Real road networks are hierarchical and it matters here, so it
is built in rather than sprinkled on.

Directedness: every lattice edge emits both arcs, with an *independent* noise
draw per direction, so the graph is genuinely asymmetric without any risk to
strong connectivity.  True one-way arcs are applied only to diagonals, which are
redundant chords over a fully bidirectional grid and therefore cannot disconnect
anything.
"""

import math

from ..graph import Graph
from ..weights import bounded_lognormal

COORD_SCALE = 10 ** 7  # integer coordinates over the unit square
TIME_SCALE = 16        # coordinate units -> time units, keeps distances < 2^31


def build(rng, params):
    V_target = int(params["V"])
    S = int(round(math.isqrt(V_target)))
    while S * S < V_target:
        S += 1
    V = S * S

    diag_p = float(params.get("diag_p", 0.50))
    anti_p = float(params.get("anti_p", 0.25))
    sigma = float(params.get("lognormal_sigma", 0.35))
    clamp = float(params.get("clamp_sigmas", 2.0))
    oneway = float(params.get("oneway_frac", 0.08))
    directed = bool(params.get("directed", True))
    arterial_stride = int(params.get("arterial_stride", 32))
    highway_stride = int(params.get("highway_stride", 256))
    speed_local = float(params.get("speed_local", 1.0))
    speed_arterial = float(params.get("speed_arterial", 3.0))
    speed_highway = float(params.get("speed_highway", 8.0))

    g = Graph(V, directed=directed)

    # --- coordinates: one node jittered inside each lattice cell -----------
    rc = rng["points"]
    xs = [0] * V
    ys = [0] * V
    cell = COORD_SCALE / S
    for r in range(S):
        base = r * S
        for c in range(S):
            i = base + c
            xs[i] = int((c + rc.random()) * cell)
            ys[i] = int((r + rc.random()) * cell)
    g.x, g.y = xs, ys

    def speed_of_line(line):
        if line % highway_stride == 0:
            return speed_highway
        if line % arterial_stride == 0:
            return speed_arterial
        return speed_local

    rw = rng["weights"]
    ro = rng["oneway"]
    rt = rng["topology"]

    def euclid(a, b):
        dx = xs[a] - xs[b]
        dy = ys[a] - ys[b]
        return math.sqrt(dx * dx + dy * dy)

    def emit(a, b, speed, both):
        d = euclid(a, b)
        base = TIME_SCALE * d / speed
        w1 = max(1, int(round(base * bounded_lognormal(rw, sigma, clamp))))
        g.add(a, b, w1)
        if both:
            if directed:
                w2 = max(1, int(round(base * bounded_lognormal(rw, sigma, clamp))))
                g.add(b, a, w2)
            # undirected: the single line already means both directions

    # --- lattice edges: always bidirectional -------------------------------
    for r in range(S):
        base = r * S
        s_row = speed_of_line(r)
        for c in range(S - 1):
            emit(base + c, base + c + 1, s_row, True)
    for c in range(S):
        s_col = speed_of_line(c)
        for r in range(S - 1):
            emit(r * S + c, (r + 1) * S + c, s_col, True)

    # --- diagonals: redundant chords, safe to make one-way -----------------
    for r in range(S - 1):
        base = r * S
        for c in range(S - 1):
            if rt.random() < diag_p:
                one = directed and ro.random() < oneway
                emit(base + c, base + S + c + 1, speed_local, not one)
        for c in range(1, S):
            if rt.random() < anti_p:
                one = directed and ro.random() < oneway
                emit(base + c, base + S + c - 1, speed_local, not one)

    meta = {
        "lattice_side": S,
        "time_scale": TIME_SCALE,
        "coord_scale": COORD_SCALE,
        "speed_classes": {
            "local": speed_local,
            "arterial": speed_arterial,
            "highway": speed_highway,
        },
        "arterial_stride": arterial_stride,
        "highway_stride": highway_stride,
        "lognormal_sigma": sigma,
        "oneway_frac": oneway,
    }
    return g, meta
