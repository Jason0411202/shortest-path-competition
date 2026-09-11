"""Query-plan construction.

The generator does not resolve rank queries itself -- resolving a rank-2^k query
means running a bounded Dijkstra, and for the flagship instance the whole set
costs ~1.9e11 settles.  That belongs in ``tools/ref/refsolve`` (C++/OpenMP),
which resolves the plan into the ``.queries`` file *and* falls out with the
reference answer for free: ``d(s,t)`` is known at the moment ``t`` is settled.

Plan file format, one line per query, already shuffled:

    Q
    P s t        explicit pair (uniform, hub or self)
    R s rank     rank query: t is the node Dijkstra from s settles at position
                 `rank`, resolved downstream

Why the mixes matter
--------------------
Uniform-random pairs on a two-million-node graph are almost all at Dijkstra
rank ~V.  Sampling nothing else would measure a single point on the
query-difficulty curve.  Rank stratification (Sanders & Schultes) sweeps
the curve instead: query cost is reported as a function of the rank at which
the target settles, and different ranks stress very different things.
"""

import math


class QueryPlan:
    def __init__(self):
        self.lines = []
        self.n_rank = 0
        self.n_pair = 0
        self.n_self = 0
        self.n_hub = 0

    def pair(self, s, t):
        self.lines.append(("P", s, t))
        self.n_pair += 1

    def rank(self, s, r):
        self.lines.append(("R", s, r))
        self.n_rank += 1


def _sources(rng, V, Q, policy):
    """Sources for the whole query set.

    Reusing a source lets one search answer several queries, which measures
    something other than what the instance is for -- so sources are distinct
    wherever Q allows it, rather than left to birthday collisions.
    """
    if policy == "strict":
        if Q > V:
            raise ValueError(
                f"distinct-sources=strict needs Q <= V, got Q={Q} V={V}; "
                f"use 'balanced' (each source reused at most ceil(Q/V) times)"
            )
        return rng.sample(range(V), Q)
    if policy == "balanced":
        reps = -(-Q // V)
        pool = list(range(V)) * reps
        rng.shuffle(pool)
        return pool[:Q]
    if policy == "free":
        return [rng.randrange(V) for _ in range(Q)]
    raise ValueError(f"unknown distinct-sources policy: {policy}")


def _parse_mix(spec):
    """``"rank:10-20"`` / ``"local:0.9@4096,uniform:0.1"`` / ``"uniform:0.7,hub:0.3"``"""
    out = []
    for part in spec.split(","):
        part = part.strip()
        kind, _, rest = part.partition(":")
        cfg = {"kind": kind}
        if kind == "rank":
            lo, _, hi = rest.partition("-")
            cfg["lo"], cfg["hi"] = int(lo), int(hi)
            cfg["frac"] = 1.0
        elif kind == "local":
            frac, _, cap = rest.partition("@")
            cfg["frac"] = float(frac)
            cfg["cap"] = int(cap) if cap else 4096
        else:
            cfg["frac"] = float(rest) if rest else 1.0
        out.append(cfg)
    return out


def build(rng, V, Q, mix_spec, hub_nodes=None, selfpair_frac=0.001,
          distinct_sources="strict"):
    rs = rng["queries.sources"]
    rk = rng["queries.ranks"]
    rh = rng["queries.hubs"]
    rz = rng["shuffle"]

    srcs = _sources(rs, V, Q, distinct_sources)
    mix = _parse_mix(mix_spec)
    plan = QueryPlan()

    n_self = int(selfpair_frac * Q)
    counts = []
    remaining = Q - n_self
    for i, cfg in enumerate(mix):
        n = remaining - sum(counts) if i == len(mix) - 1 else int(cfg["frac"] * remaining)
        counts.append(max(0, n))

    idx = 0
    for cfg, n in zip(mix, counts):
        kind = cfg["kind"]
        if kind == "rank":
            lo, hi = cfg["lo"], cfg["hi"]
            buckets = list(range(lo, hi + 1))
            for j in range(n):
                s = srcs[idx]
                idx += 1
                plan.rank(s, 1 << buckets[j % len(buckets)])
        elif kind == "local":
            top = max(1, int(math.log2(cfg["cap"])))
            for j in range(n):
                s = srcs[idx]
                idx += 1
                plan.rank(s, 1 << rk.randint(4, top))
        elif kind == "hub":
            if not hub_nodes:
                raise ValueError("hub mix requested but the family exposed no hubs")
            for _ in range(n):
                s = srcs[idx]
                idx += 1
                t = hub_nodes[rh.randrange(len(hub_nodes))]
                if t == s:
                    t = (t + 1) % V
                plan.pair(s, t)
                plan.n_hub += 1
        elif kind == "uniform":
            for _ in range(n):
                s = srcs[idx]
                idx += 1
                t = rh.randrange(V)
                if t == s:
                    t = (t + 1) % V
                plan.pair(s, t)
        else:
            raise ValueError(f"unknown query kind: {kind}")

    # s == t is free to include, the foundation defines the answer (0), and it
    # catches solvers that special-case the empty search incorrectly.
    for _ in range(n_self):
        s = srcs[idx]
        idx += 1
        plan.pair(s, s)
        plan.n_self += 1

    # Stratify, then shuffle, so the file leaks no block structure.  The bucket
    # labels survive in the metadata for per-rank timing reports.
    rz.shuffle(plan.lines)
    return plan


def write(plan, path):
    with open(path, "w") as f:
        f.write(f"{len(plan.lines)}\n")
        f.writelines(f"{k} {a} {b}\n" for k, a, b in plan.lines)
