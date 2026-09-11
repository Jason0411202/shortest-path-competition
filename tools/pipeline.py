#!/usr/bin/env python3
"""Build instances end to end and refuse to ship a bad one.

    python3 tools/pipeline.py --manifest tools/manifest/public.json --all \
                              --out instances

For each instance:
  1. tools.gen            -> .graph, .qplan, partial .meta.json
  2. refsolve plan        -> .queries, .answers   (one pass: a rank-2^k query
                             yields d(s,t) the moment t settles)
  3. refsolve verify      -> independent recomputation of a sample, symmetry on
                             undirected instances, every -1 confirmed by a
                             separate reachability sweep
  4. refsolve probe       -> Dijkstra vs A* search space
  5. thresholds           -> the checks below, which FAIL THE BUILD

The thresholds are the point of this script: a road family whose edge weights
are exactly proportional to Euclidean length would hand out a large unintended
advantage, so the build checks for it instead of trusting the generator.
"""

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REFSOLVE = os.path.join(HERE, "ref", "refsolve")

# family -> (check name, predicate on the merged metadata, message)
THRESHOLDS = {
    "road2d": [
        ("astar_not_a_giveaway",
         lambda m: (m.get("probe", {}).get("astar_reduction") or 0) <= 2.0,
         "a plain Euclidean A* should be worth <= 2x here; if it is more, the "
         "weights have drifted back towards being proportional to distance"),
    ],
    "lattice": [],
    "scalefree": [
        ("unreachable_fraction_sane",
         lambda m: 0.0 <= m.get("resolve", {}).get("unreachable", 0) / max(1, m["Q"]) <= 0.20,
         "unreachable pairs should be a few percent: high fractions inflate the "
         "CH-vs-Dijkstra gap without testing anything"),
    ],
}

# name substring -> extra checks
NAME_THRESHOLDS = {
    "wide64": [
        ("distances_overflow_uint32",
         lambda m: not m.get("resolve", {}).get("fits_uint32", True),
         "wide64 exists to break 32-bit distance packing; max distance must "
         "exceed 2^32 by a wide margin, not marginally"),
    ],
}


def run(cmd, capture_json=False):
    print("  $ " + " ".join(str(c) for c in cmd), flush=True)
    r = subprocess.run([str(c) for c in cmd], capture_output=capture_json, text=True)
    if r.returncode != 0:
        if capture_json:
            sys.stderr.write(r.stdout or "")
            sys.stderr.write(r.stderr or "")
        raise SystemExit(f"FAILED: {' '.join(str(c) for c in cmd)}")
    if capture_json:
        sys.stderr.write(r.stderr or "")
        return json.loads(r.stdout.strip().splitlines()[-1])
    return None


def build_one(name, outdir, manifest, keep_plan):
    print(f"\n=== {name} ===", flush=True)
    run([sys.executable, "-m", "tools.gen", "--manifest", manifest,
         "--instance", name, "--out", outdir])

    graph = os.path.join(outdir, f"{name}.graph")
    qplan = os.path.join(outdir, f"{name}.qplan")
    queries = os.path.join(outdir, f"{name}.queries")
    answers = os.path.join(outdir, f"{name}.answers")
    metap = os.path.join(outdir, f"{name}.meta.json")

    resolve = run([REFSOLVE, "plan", graph, qplan, queries, answers], capture_json=True)
    verify = run([REFSOLVE, "verify", graph, queries, answers, 5000], capture_json=True)
    probe = run([REFSOLVE, "probe", graph, queries, 200], capture_json=True)

    meta = json.load(open(metap))
    meta["resolve"] = resolve
    meta["verify"] = verify
    meta["probe"] = probe

    from tools.gen import writer
    writer.finalise_meta(meta, outdir, name)
    writer.write_meta(meta, metap)

    if not keep_plan:
        os.remove(qplan)

    problems = []
    if verify.get("failures", 1) != 0:
        problems.append(("answers_verified",
                         f"{verify['failures']} answers disagreed with an "
                         f"independent recomputation"))
    checks = list(THRESHOLDS.get(meta["family"], []))
    for frag, extra in NAME_THRESHOLDS.items():
        if frag in name:
            checks += extra
    for check_name, pred, msg in checks:
        try:
            ok = pred(meta)
        except Exception as exc:                      # a missing field is a failure
            ok, msg = False, f"{msg} (check raised {exc!r})"
        if not ok:
            problems.append((check_name, msg))

    if problems:
        print(f"  [{name}] REJECTED:", flush=True)
        for cn, msg in problems:
            print(f"    - {cn}: {msg}", flush=True)
        return meta, False

    print(f"  [{name}] ok  V={meta['V']:,} E={meta['E']:,} Q={meta['Q']:,} "
          f"unreachable={resolve['unreachable']:,} "
          f"astar={probe.get('astar_reduction')}", flush=True)
    return meta, True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--instance", action="append", default=[])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--out", default="instances")
    ap.add_argument("--keep-plan", action="store_true",
                    help="keep the .qplan files (they are an intermediate)")
    ap.add_argument("--checksums", help="write a sha256sum-format file here")
    args = ap.parse_args()

    if not os.path.exists(REFSOLVE):
        raise SystemExit(f"{REFSOLVE} not built -- run: make -C tools/ref")

    man = json.load(open(args.manifest))
    names = list(man["instances"]) if args.all else args.instance
    if not names:
        raise SystemExit("pass --instance NAME (repeatable) or --all")

    metas, failed = [], []
    for name in names:
        meta, ok = build_one(name, args.out, args.manifest, args.keep_plan)
        metas.append(meta)
        if not ok:
            failed.append(name)

    print("\n" + "=" * 78)
    print(f"{'instance':<20}{'V':>10}{'E':>12}{'Q':>10}{'unreach':>9}"
          f"{'maxdist':>14}{'A*':>7}")
    print("-" * 78)
    for m in metas:
        pr = m.get("probe", {}).get("astar_reduction")
        print(f"{m['name']:<20}{m['V']:>10,}{m['E']:>12,}{m['Q']:>10,}"
              f"{m['resolve']['unreachable']:>9,}"
              f"{m['resolve']['max_finite_distance']:>14,}"
              f"{('-' if pr is None else f'{pr:.2f}'):>7}")

    if args.checksums:
        with open(args.checksums, "w") as f:
            for m in sorted(metas, key=lambda m: m["name"]):
                for ext, digest in sorted(m["sha256"].items()):
                    f.write(f"{digest}  {m['name']}.{ext}\n")
        print(f"\nchecksums -> {args.checksums}")

    if failed:
        raise SystemExit(f"\n{len(failed)} instance(s) rejected: {', '.join(failed)}")
    print("\nall instances passed")
    return 0


if __name__ == "__main__":
    if ROOT not in sys.path:
        sys.path.insert(0, ROOT)
    sys.exit(main())
