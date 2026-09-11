"""CLI and orchestration for the instance generator.

Usage
-----
    python3 -m tools.gen --manifest tools/manifest/public.json \
                         --instance road2d_dev --out instances

    python3 -m tools.gen --manifest tools/manifest/public.json --all --out instances

    python3 -m tools.gen --compare-manifests tools/manifest/public.json \
                                             tools/manifest/hidden.json

Emits ``<name>.graph``, ``<name>.qplan`` and a partial ``<name>.meta.json``.
Turning the plan into ``<name>.queries`` + ``<name>.answers`` is the job of
``tools/ref/refsolve``, which runs the bounded Dijkstras in parallel.

Manifests are JSON rather than YAML so the generator stays stdlib-only: the
students get this code, and a dependency they have to install is a dependency
that will differ between their machine and yours.
"""

import argparse
import json
import os
import secrets
import sys
from array import array

from . import queries as qmod
from . import writer
from .families import lattice, road2d, scalefree
from .rng import stream

FAMILIES = {
    "road2d": road2d.build,
    "lattice": lattice.build,
    "scalefree": scalefree.build,
}

STREAMS = ("points", "topology", "weights", "oneway", "sinks", "indegree",
           "queries.sources", "queries.ranks", "queries.hubs", "shuffle")


def rng_bundle(seed, salt):
    """A dict of independent named streams.

    Independent per name is the point: adding a query type must not perturb the
    graph, so regenerating queries against an existing graph reproduces it
    exactly.
    """
    return {name: stream(seed, salt, name) for name in STREAMS}


def admissible_lambda(g):
    """min_e w(e)/||u-v||, the largest scale an admissible Euclidean A*
    potential may use.

    With ``w = round(1e6 * ||u-v||)`` the ratio would be constant and the graph
    metric a scaled copy of the coordinate embedding.  A travel-time metric
    drags lambda down to the slowest road class, so coordinates only loosely
    predict distances.  Recorded in the metadata so
    a regression is visible rather than silent.
    """
    if g.x is None:
        return None
    import math
    x, y, eu, ev, ew = g.x, g.y, g.eu, g.ev, g.ew
    best = None
    for i in range(g.E):
        a, b = eu[i], ev[i]
        d = math.hypot(x[a] - x[b], y[a] - y[b])
        if d <= 0:
            continue
        r = ew[i] / d
        if best is None or r < best:
            best = r
    return best


def hub_nodes(g, frac=0.001):
    """The top `frac` of nodes by in-degree, for hub-target query mixes."""
    indeg = array("i", bytes(4 * g.V))
    for i in range(g.E):
        indeg[g.ev[i]] += 1
        if not g.directed:
            indeg[g.eu[i]] += 1
    k = max(1, int(frac * g.V))
    order = sorted(range(g.V), key=lambda v: indeg[v], reverse=True)
    return order[:k]


def generate_one(name, spec, seed, salt, outdir, argv):
    family = spec["family"]
    if family not in FAMILIES:
        raise SystemExit(f"unknown family: {family}")
    rng = rng_bundle(seed, salt)

    print(f"[{name}] building {family} ...", flush=True)
    g, fam_meta = FAMILIES[family](rng, spec["params"])

    print(f"[{name}] V={g.V:,} E={g.E:,} flags={g.flags} -- checking invariants",
          flush=True)
    sizes = g.check(expect_components=spec.get("expect_components", 1))

    hubs = None
    if "hub" in spec.get("mix", ""):
        hubs = hub_nodes(g, spec.get("hub_frac", 0.001))

    print(f"[{name}] planning {spec['Q']:,} queries ({spec['mix']})", flush=True)
    plan = qmod.build(
        rng, g.V, int(spec["Q"]), spec["mix"],
        hub_nodes=hubs,
        selfpair_frac=spec.get("selfpair_frac", 0.001),
        distinct_sources=spec.get("distinct_sources", "strict"),
    )

    os.makedirs(outdir, exist_ok=True)
    gpath = os.path.join(outdir, f"{name}.graph")
    writer.write_graph(g, gpath)
    qmod.write(plan, os.path.join(outdir, f"{name}.qplan"))

    meta = writer.base_meta(name, g, family, spec["params"], seed, salt, argv, sizes)
    meta["family_detail"] = fam_meta
    meta["lambda_admissible"] = admissible_lambda(g)
    meta["Q"] = len(plan.lines)
    meta["query_mix"] = spec["mix"]
    meta["distinct_sources_policy"] = spec.get("distinct_sources", "strict")
    meta["query_counts"] = {
        "rank": plan.n_rank,
        "pair": plan.n_pair,
        "hub": plan.n_hub,
        "selfpair": plan.n_self,
    }
    meta["distinct_sources"] = len({a for _, a, _ in plan.lines})
    writer.write_meta(meta, os.path.join(outdir, f"{name}.meta.json"))
    print(f"[{name}] wrote {gpath} ({os.path.getsize(gpath) / 1e6:.1f} MB)",
          flush=True)
    return meta


def compare_manifests(a_path, b_path):
    """Public and hidden manifests must differ only in seed and salt.

    Changing parameters between the released set and the graded set is a moving
    target, and students would be right to object to it.  This is the check that
    makes "same code, same parameters, different seed" enforceable rather than a
    promise.
    """
    a = json.load(open(a_path))
    b = json.load(open(b_path))
    if a["instances"] != b["instances"]:
        only_a = set(a["instances"]) - set(b["instances"])
        only_b = set(b["instances"]) - set(a["instances"])
        diffs = [k for k in set(a["instances"]) & set(b["instances"])
                 if a["instances"][k] != b["instances"][k]]
        print("MANIFESTS DIFFER")
        if only_a:
            print(f"  only in {a_path}: {sorted(only_a)}")
        if only_b:
            print(f"  only in {b_path}: {sorted(only_b)}")
        for k in sorted(diffs):
            print(f"  parameters differ for {k}")
        return 1
    if (a["seed"], a["salt"]) == (b["seed"], b["salt"]):
        print("MANIFESTS SHARE A SEED -- the hidden set would be the public set")
        return 1
    print(f"OK: {len(a['instances'])} instances, parameter-identical, distinct seeds")
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    ap = argparse.ArgumentParser(prog="python3 -m tools.gen")
    ap.add_argument("--manifest")
    ap.add_argument("--instance", action="append", default=[])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--out", default="instances")
    ap.add_argument("--seed", type=int, help="override the manifest seed")
    ap.add_argument("--salt", help="override the manifest salt")
    ap.add_argument("--compare-manifests", nargs=2, metavar=("A", "B"))
    ap.add_argument("--make-hidden", nargs=2, metavar=("PUBLIC", "OUT"),
                    help="copy a manifest with a fresh random seed and salt")
    args = ap.parse_args(argv)

    if args.compare_manifests:
        return compare_manifests(*args.compare_manifests)

    if args.make_hidden:
        src, dst = args.make_hidden
        man = json.load(open(src))
        man["seed"] = secrets.randbelow(2 ** 63)
        man["salt"] = "hidden-" + secrets.token_hex(8)
        man["_comment"] = [
            "INSTRUCTOR ONLY -- never publish this file or the instances it",
            "generates. Parameter-identical to the public manifest by",
            "construction; only seed and salt differ. Verify with:",
            f"  python3 -m tools.gen --compare-manifests {src} {dst}",
        ]
        with open(dst, "w") as f:
            json.dump(man, f, indent=2)
            f.write("\n")
        print(f"wrote {dst} (seed hidden -- keep this file out of the student repo)")
        return 0

    if not args.manifest:
        ap.error("--manifest is required (or use --compare-manifests)")
    man = json.load(open(args.manifest))
    seed = args.seed if args.seed is not None else man["seed"]
    salt = args.salt if args.salt is not None else man["salt"]

    names = list(man["instances"]) if args.all else args.instance
    if not names:
        ap.error("pass --instance NAME (repeatable) or --all")
    for name in names:
        if name not in man["instances"]:
            raise SystemExit(f"no such instance in manifest: {name}")
        generate_one(name, man["instances"][name], seed, salt, args.out,
                     " ".join(argv))
    return 0
