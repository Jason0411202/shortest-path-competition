"""SCALEFREE: fixed out-degree, power-law in-degree, genuinely directed.

Every non-sink node has exactly ``m`` out-arcs, while in-degree follows a power
law with exponent ``gamma``, so the top nodes collect tens of thousands of
in-arcs.  Out-degree and in-degree are therefore wildly asymmetric, and the
graph is directed in earnest: there are no reverse arcs, and the query mix
deliberately aims a share of its targets at the highest-in-degree nodes.

A fraction of nodes are sinks with no out-arcs, which is where the instance's
unreachable pairs come from.  On a directed graph unreachability is a real
question; on an undirected one it would be a component lookup.

Construction is a directed configuration model rather than index-ordered
preferential attachment.  That matters: in plain BA-style attachment every arc
points at an *earlier* node, so the reachable set from node i is contained in
[0, i) and roughly half of all uniform-random pairs are unreachable by
construction.  Drawing in-degrees first and then matching stubs removes the
index bias, so unreachability comes only from sinks and from nodes that drew
zero in-degree -- a few percent, which is what the query design wants.
"""

from array import array

from ..graph import Graph, ensure_weakly_connected
from ..weights import make as make_weight


def build(rng, params):
    V = int(params["V"])
    m = int(params.get("m", 3))
    gamma = float(params.get("gamma", 2.1))
    sink_frac = float(params.get("sink_frac", 0.03))
    wspec = params.get("weights", "loguniform:1:10000")

    g = Graph(V, directed=True)
    weight = make_weight(wspec, None)
    rw = rng["weights"]
    rt = rng["topology"]
    rs = rng["sinks"]

    sink = bytearray(V)
    n_sinks = 0
    for i in range(V):
        if rs.random() < sink_frac:
            sink[i] = 1
            n_sinks += 1

    n_arcs = m * (V - n_sinks)

    # --- draw a power-law in-degree sequence summing to n_arcs ------------
    # Pareto weights p_i = u^(-1/(gamma-1)), normalised.  gamma = 2.1 gives the
    # heavy tail we want (max in-degree in the tens of thousands at V = 2M).
    #
    # One stub per node is reserved up front so that *every* node has in-degree
    # at least 1.  Without that floor a gamma = 2.1 draw leaves ~33% of nodes
    # with in-degree 0 (measured), and since those nodes can never be reached,
    # ~44% of uniform-random queries come back -1.  Unreachability is supposed
    # to come from the sinks -- a few percent, enough to catch a bad `-1` path
    # -- not to dominate the instance: an unreachable query costs plain Dijkstra
    # a full sweep and costs CH nothing, so a large fraction inflates the gap
    # between them without testing anything.
    reserved = min(V, n_arcs)
    tail = n_arcs - reserved
    expo = -1.0 / (gamma - 1.0)
    rp = rng["indegree"]
    wts = [0.0] * V
    total = 0.0
    for i in range(V):
        p = (1.0 - rp.random()) ** expo
        wts[i] = p
        total += p
    scale = tail / total if total else 0.0

    stubs = array("i", range(reserved))
    for i in range(V):
        k = int(wts[i] * scale)
        if k:
            stubs.extend([i] * k)
    # top up / trim to exactly n_arcs
    while len(stubs) < n_arcs:
        stubs.append(rp.randrange(V))
    if len(stubs) > n_arcs:
        del stubs[n_arcs:]

    # Fisher-Yates over the stub list, in place.
    for i in range(len(stubs) - 1, 0, -1):
        j = rt.randrange(i + 1)
        stubs[i], stubs[j] = stubs[j], stubs[i]

    # --- match stubs to sources -------------------------------------------
    indeg = array("i", bytes(4 * V))
    cursor = 0
    nstub = len(stubs)
    for i in range(V):
        if sink[i]:
            continue
        chosen = []
        for _ in range(m):
            # skip past self-loops and duplicates by swapping the offending
            # stub with a random unconsumed one
            for _retry in range(32):
                t = stubs[cursor]
                if t != i and t not in chosen:
                    break
                j = cursor + 1 + rt.randrange(max(1, nstub - cursor - 1))
                if j < nstub:
                    stubs[cursor], stubs[j] = stubs[j], stubs[cursor]
                else:
                    break
            t = stubs[cursor]
            cursor += 1
            if t == i or t in chosen:
                continue
            chosen.append(t)
        for t in chosen:
            g.add(i, t, weight(rw))
            indeg[t] += 1

    merged = ensure_weakly_connected(g, rt, lambda: weight(rw))

    meta = {
        "m": m,
        "gamma": gamma,
        "sink_frac": sink_frac,
        "n_sinks": n_sinks,
        "weights": wspec,
        "components_merged": merged,
        "max_indegree": max(indeg),
        "zero_indegree_nodes": sum(1 for i in range(V) if indeg[i] == 0),
    }
    return g, meta
