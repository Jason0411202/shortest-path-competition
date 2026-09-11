"""Instance serialisation.

Text format is kept deliberately -- see README.  Once Q is sized correctly,
parsing is a couple of percent of runtime, so a binary format would buy nothing
and would force every student to rewrite the parser, weakening the "derivative
of the foundation" rule.

The one change is that coordinates are written as *integers*.  ``strtod`` is
5-10x slower than integer parsing, and ranking a student on whether they
remembered to avoid ``%lf`` measures nothing anyone cares about.  It also takes
float determinism out of any A* potential.
"""

import hashlib
import json
import os


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def write_graph(g, path):
    with open(path, "w") as f:
        f.write(f"{g.V} {g.E} {g.flags}\n")
        eu, ev, ew = g.eu, g.ev, g.ew
        f.writelines(f"{eu[i]} {ev[i]} {ew[i]}\n" for i in range(g.E))
        if g.x is not None:
            x, y = g.x, g.y
            f.writelines(f"{x[i]} {y[i]}\n" for i in range(g.V))


def write_meta(meta, path):
    with open(path, "w") as f:
        json.dump(meta, f, indent=2, sort_keys=True)
        f.write("\n")


def base_meta(name, g, family, params, seed, salt, argv, component_sizes):
    return {
        "name": name,
        "family": family,
        "generator_version": 2,
        "argv": argv,
        "seed": seed,
        "salt": salt,
        "V": g.V,
        "E": g.E,
        "flags": g.flags,
        "directed": g.directed,
        "has_coords": g.x is not None,
        "n_components": len(component_sizes),
        "component_sizes_top10": component_sizes[:10],
        "w_min": min(g.ew),
        "w_max": max(g.ew),
        "params": params,
    }


def finalise_meta(meta, outdir, name):
    """Attach file digests last, once every artefact exists."""
    digests = {}
    for ext in ("graph", "queries", "answers"):
        p = os.path.join(outdir, f"{name}.{ext}")
        if os.path.exists(p):
            digests[ext] = sha256_file(p)
    meta["sha256"] = digests
    return meta
