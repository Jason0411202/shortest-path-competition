#!/usr/bin/env python3
"""Regression tests for grade.py.

    python3 tools/test_grade.py

Covers the defects the v2 grader was written to fix, so that they cannot come
back quietly:

  1. foundation-as-solver scores ~1.0 on every instance -- which is also an
     end-to-end check that baseline extrapolation agrees with reality.
  2. a wrong answer scores the floor AND still counts towards the geometric
     mean, so failing your weakest instance can never raise your score.
  3. the timeout, the 4 GB cap and the single-thread rule are enforced rather
     than merely written down, and one bad instance does not abort the run.
  4. `--baseline-probe` agrees with `--baseline-full`.

Timing thresholds here are deliberately loose: on a shared machine the dev-tier
instances run for only a few seconds, and a few seconds cannot be timed
tightly.  The tests check that the machinery works, not that the clock is
quiet.
"""

import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GRADE = ROOT / "grade.py"
FOUNDATION = ROOT / "foundation"
# The smallest committed instance, so the suite stays quick.
INST = ("LOCAL", ROOT / "instances/local2d_dev.graph",
        ROOT / "instances/local2d_dev.queries")

PASS, FAIL = [], []


def check(name, ok, detail=""):
    (PASS if ok else FAIL).append(name)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   {detail}" if detail else ""))


def run_grade(tmp, solver, extra=(), instances=None):
    ilist = Path(tmp) / "i.txt"
    rows = instances or [INST]
    ilist.write_text("".join(f"{c} {g} {q}\n" for c, g, q in rows))
    out = Path(tmp) / "r.json"
    out.unlink(missing_ok=True)
    r = subprocess.run(
        [sys.executable, str(GRADE), "--solver", str(solver),
         "--instances", str(ilist), "--json", str(out), *extra],
        cwd=ROOT, capture_output=True, text=True)
    if not out.exists():
        print(r.stdout, r.stderr)
        raise SystemExit(f"grade.py produced no result for {solver}")
    return json.loads(out.read_text()), r


def stub(tmp, name, body):
    p = Path(tmp) / name
    p.write_text("#!/usr/bin/env bash\n" + body)
    p.chmod(0o755)
    return p


def main():
    if not FOUNDATION.exists():
        raise SystemExit("build the foundation first: make foundation")

    with tempfile.TemporaryDirectory() as tmp:
        print("\n[1] foundation as its own solver")
        res, _ = run_grade(tmp, FOUNDATION)
        inst = res["instances"][0]
        check("scores ok", inst["note"] == "ok", f"note={inst['note']}")
        sp = inst["speedup"]
        # Wide on purpose. Measured run-to-run sd on this instance is ~7% and
        # block drift another ~10%, so a tight band here would flake without
        # telling anyone anything. The systematic part is checked in [5].
        check("speedup within 45% of 1.0", 0.7 <= sp <= 1.45, f"speedup={sp:.3f}")
        check("baseline was extrapolated, not run in full",
              inst.get("baseline", {}).get("mode") == "probe")

        print("\n[1b] the source digest binds solver.cpp to the timings")
        src = Path(tmp) / "mine.cpp"
        src.write_text("int main(){return 0;}\n")
        res, _ = run_grade(tmp, FOUNDATION, extra=["--source", str(src)])
        check("records the source digest",
              res["solver_source_sha256"] == hashlib.sha256(
                  src.read_bytes()).hexdigest(),
              f"got {str(res['solver_source_sha256'])[:16]}...")
        check("records which file it hashed",
              res["solver_source"] == str(src))
        res, _ = run_grade(tmp, FOUNDATION,
                           extra=["--source", str(Path(tmp) / "absent.cpp")])
        check("absent source gives null, not a crash",
              res["solver_source_sha256"] is None)

        print("\n[2] wrong answers")
        # copy the foundation's output, then corrupt one line
        bad = stub(tmp, "wrong.sh",
                   f'"{FOUNDATION}" "$1" "$2" "$3"\n'
                   'python3 - "$3" <<\'PY\'\n'
                   'import sys\n'
                   'p=sys.argv[1]; L=open(p).read().split("\\n")\n'
                   'L[0]=str(int(L[0] or 0)+1)\n'
                   'open(p,"w").write("\\n".join(L))\n'
                   'PY\n')
        res, _ = run_grade(tmp, bad)
        inst = res["instances"][0]
        check("flagged WRONG", inst["note"] == "WRONG", f"note={inst['note']}")
        check("scores the floor", inst["speedup"] == 0.1, f"speedup={inst['speedup']}")
        # a failed instance must drag the mean down, not vanish from it
        check("counts towards the geomean",
              abs(res["overall_geomean"] - 0.1) < 1e-9,
              f"overall={res['overall_geomean']}")

        print("\n[3] limits")
        res, _ = run_grade(tmp, stub(tmp, "slow.sh", "sleep 30\n"), ["--timeout", "3"])
        check("timeout is enforced", res["instances"][0]["note"] == "TIMEOUT",
              f"note={res['instances'][0]['note']}")

        # Fill 5 GB with incompressible data.  Linux kills this at the
        # allocation via RLIMIT_AS; macOS cannot lower RLIMIT_AS, so only peak
        # the sampled footprint can catch it, and `exec` keeps the hog in the
        # pid grade.py is watching.
        hog = stub(tmp, "hog.sh",
                   'exec python3 -c "import os; n=5*1024*1024*1024; b=bytearray(n); '
                   'm=memoryview(b); c=os.urandom(1<<20)\n'
                   'for i in range(0, n, 1<<20): m[i:i+(1<<20)]=c\n'
                   'print(n)"\n')
        res, _ = run_grade(tmp, hog)
        check("4 GB cap is enforced",
              res["instances"][0]["note"] in ("MEMORY", "CRASH"),
              f"note={res['instances'][0]['note']}")

        threads = stub(tmp, "threads.sh",
                       f'"{FOUNDATION}" "$1" "$2" "$3" &\n'
                       f'"{FOUNDATION}" "$1" "$2" "$3.b" &\n'
                       f'"{FOUNDATION}" "$1" "$2" "$3.c" &\n'
                       'wait\n')
        res, _ = run_grade(tmp, threads)
        check("parallel work is detected", res["instances"][0]["note"] == "THREADS",
              f"note={res['instances'][0]['note']}")

        print("\n[4] one bad instance does not abort the run")
        rows = [INST, ("GLOBAL", ROOT / "instances/lattice3d_dev.graph",
                       ROOT / "instances/lattice3d_dev.queries")]
        res, _ = run_grade(tmp, stub(tmp, "crash.sh", "exit 3\n"), instances=rows)
        check("both instances reported", len(res["instances"]) == 2,
              f"got {len(res['instances'])}")

        print("\n[5] extrapolated baseline vs a full baseline run")
        a, _ = run_grade(tmp, FOUNDATION)
        b, _ = run_grade(tmp, FOUNDATION, ["--baseline-full"])
        ta = a["instances"][0]["t_base"]
        tb = b["instances"][0]["t_base"]
        if ta is None or tb is None or tb == 0:
            check("agree within 30%", False,
                  f"no baseline measured: probe={ta} full={tb} "
                  f"note={a['instances'][0]['note']}/{b['instances'][0]['note']}")
            print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
            print("failed: " + ", ".join(FAIL))
            return 1
        err = abs(ta - tb) / tb
        # Single-shot comparison, so this is mostly a machine-noise check.
        # The systematic sampling bias, measured over 5 repetitions, is under
        # ~3%; anything beyond 30% here means the fit itself has broken.
        check("agree within 30%", err < 0.30,
              f"probe={ta:.3f}s full={tb:.3f}s err={err:.1%}")

    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        print("failed: " + ", ".join(FAIL))
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
