#!/usr/bin/env python3
"""
grade_class.py -- map a class of competition scores onto 70-100 grades.

Input is the set of final `result.json` files, one per student, produced by
`grade.py --json` on the fresh instances.  Output is a grade per student.

Usage:
  python3 grade_class.py results/                    # one <student_id>.json each
  python3 grade_class.py results/ --csv grades.csv
  python3 grade_class.py a.json b.json c.json

A student ID is taken from the file stem (`B11012345.json` -> `B11012345`), or
from the containing directory when the file is named `result.json`
(`B11012345/result.json` -> `B11012345`).

--- How the mapping works ---------------------------------------------------

`overall_geomean` is a geometric mean of six speedups, so its logarithm is a
sum of six per-instance log-speedups and is close to normally distributed
across a class.  Grading therefore works on ln G rather than G:

    z     = (ln G - mean(ln G)) / sd(ln G)
    grade = clip(MEAN + SPREAD * z, FLOOR, CEIL)

The class mean is MEAN by construction, and the shape is normal.  Working in
the log is what keeps the top tail short: a student who is 10x faster than the
median is about 2.3 natural-log units out, not 10 grade points times ten, so
climbing from 3rd place to 1st is worth very little.  That is deliberate --
the marginal grade per hour (or per API dollar) spent at the very top of the
board is small, while the marginal grade for a first real algorithmic
improvement, down where the class bunches, is large.

SPREAD is the only dial.  Simulated over 300 synthetic classes of 100:

    spread   mean   >95   >90   at floor   top grade
         5   83.0   0.6   7.5        0.5        95.1
         6   83.0   1.8  11.7        1.9        97.3
         7   83.1   3.7  15.7        4.0        99.0
         8   83.2   6.1  19.2        6.3        99.7

7 is chosen so that a handful of students can realistically reach the high
90s, at the cost of a few landing on the 70 floor.

A student who fails an instance needs no special case: the 0.1 floor that
`grade.py` gives a WRONG / TIMEOUT / MEMORY / THREADS instance pulls their
geometric mean down far enough to put them at or near the grade floor.

Note this is norm-referenced -- a grade depends on the rest of the class, so
it can only be computed once every submission is in.
"""

import argparse
import csv
import json
import math
import signal
import statistics
import sys
from pathlib import Path

MEAN = 83.0      # class mean, by construction
SPREAD = 7.0     # grade points per standard deviation of ln G
FLOOR = 70.0
CEIL = 100.0


def student_id(path):
    """`B11012345.json` -> B11012345;  `B11012345/result.json` -> B11012345."""
    return path.parent.name if path.stem == "result" else path.stem


def collect(paths):
    """Expand directories, then read `overall_geomean` out of each file."""
    files = []
    for p in paths:
        if p.is_dir():
            files.extend(sorted(p.rglob("*.json")))
        else:
            files.append(p)

    entries, bad = [], []
    for f in files:
        try:
            payload = json.loads(f.read_text())
            g = float(payload["overall_geomean"])
        except (OSError, ValueError, KeyError, TypeError) as exc:
            bad.append((f, exc))
            continue
        if not math.isfinite(g) or g <= 0:
            bad.append((f, f"overall_geomean = {g!r}"))
            continue
        entries.append({"id": student_id(f), "path": f, "geomean": g})
    return entries, bad


def grade_all(entries, mean, spread):
    logs = [math.log(e["geomean"]) for e in entries]
    mu = statistics.mean(logs)
    # Population, not sample: this is the whole class, not a draw from it.
    sd = statistics.pstdev(logs)
    for e, x in zip(entries, logs):
        # A class where everyone scored the same has nothing to spread out.
        e["z"] = 0.0 if sd == 0 else (x - mu) / sd
        e["grade"] = min(CEIL, max(FLOOR, mean + spread * e["z"]))
    return mu, sd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", type=Path,
                    help="result.json files, or directories holding them")
    ap.add_argument("--mean", type=float, default=MEAN, help=f"class mean (default {MEAN:g})")
    ap.add_argument("--spread", type=float, default=SPREAD,
                    help=f"grade points per sd of ln G (default {SPREAD:g})")
    ap.add_argument("--csv", type=Path, help="also write the table as CSV")
    args = ap.parse_args()

    # A 100-line table invites `| head`, which otherwise ends in a traceback.
    if hasattr(signal, "SIGPIPE"):
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    entries, bad = collect(args.paths)
    for f, why in bad:
        print(f"skipped {f}: {why}", file=sys.stderr)
    if not entries:
        print("no usable result.json files found", file=sys.stderr)
        return 1

    mu, sd = grade_all(entries, args.mean, args.spread)
    entries.sort(key=lambda e: -e["geomean"])

    print(f"{'rank':>4} {'student':<16} {'geomean':>10} {'z':>7} {'grade':>7}")
    print("-" * 48)
    for rank, e in enumerate(entries, 1):
        print(f"{rank:>4} {e['id']:<16} {e['geomean']:>10.3f} "
              f"{e['z']:>7.2f} {e['grade']:>7.1f}")

    grades = [e["grade"] for e in entries]
    print()
    print(f"n = {len(grades)}   ln G: mean {mu:.3f}, sd {sd:.3f}")
    print(f"grades: mean {statistics.mean(grades):.1f}  "
          f"median {statistics.median(grades):.1f}  "
          f"min {min(grades):.1f}  max {max(grades):.1f}")
    print(f"        >95: {sum(g > 95 for g in grades)}   "
          f">90: {sum(g > 90 for g in grades)}   "
          f"at floor: {sum(g <= FLOOR + 0.01 for g in grades)}")

    if args.csv:
        with args.csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["rank", "student_id", "overall_geomean", "z", "grade"])
            for rank, e in enumerate(entries, 1):
                w.writerow([rank, e["id"], f"{e['geomean']:.6f}",
                            f"{e['z']:.4f}", f"{e['grade']:.1f}"])
        print(f"\nwrote {args.csv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
