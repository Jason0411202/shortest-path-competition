#!/usr/bin/env python3
"""
grade.py -- score a student submission for the Shortest-Path competition.

Per instance:
  1. Check the solver's output against the shipped answer key (`<name>.answers`,
     pinned by SHA-256).  Any mismatch -> the instance scores CLIP_LO.
  2. Estimate T_base, the unmodified foundation's time on the full query set.
  3. Time the solver.
  4. speedup = T_base / T_solver, clipped to [CLIP_LO, CLIP_HI].

Aggregate: geometric mean per category, per track, and overall.

Usage:
  python3 grade.py --solver ./solver --instances instances.txt
  python3 grade.py --solver ./solver --instances instances.txt --json out.json

where instances.txt has one line per instance:
   <category> <graph_path> <query_path>

Categories: GLOBAL (throughput -- long-range queries dominate) or
LOCAL (latency -- short-range queries and per-query fixed costs dominate).

--- Why T_base is estimated rather than measured in full ---------------------

The scored instances carry up to 600,000 queries.  Running the foundation
over all of them takes about 1.7 hours, and the baseline has to be re-timed
on every grading machine, so it is not measured in full.

Instead the foundation is timed on two short prefixes of the query file and
T(q) = a + b*q is fitted: `a` absorbs graph parsing, `b` the per-query cost
(including the foundation's per-query std::fill over dist[], which is a large
constant at these sizes).  The form is exact, so the only error is in the
sample.

Measured on the dev tier, 5 repetitions per configuration: the systematic
sampling bias of a 2,000-query probe is +3.2% on lattice3d_dev and -0.1% on
local2d_dev.  What actually dominates is machine noise -- run-to-run standard
deviation was 3.9% on a 75-second instance and 7.2% on a 5-second one, with
block-to-block drift of about 10% on top.  So the extrapolation is not the
weak link; short runs are.  That is also why the dev tier is not scored.

Note that the sampling part of the error is a per-instance constant: it is the
same queries for everybody, so it scales every student's speedup on that
instance by the same factor and cancels out of the ranking.  Only the machine
noise moves people relative to each other, which is what the repeats are for.

Ratio scoring is kept because times are self-reported from student hardware --
the ratio is what normalises for the machine.  Pass --baseline-full on a small
instance to check the extrapolation against a real full run.
"""

import argparse
import hashlib
import json
import math
import os
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from collections import defaultdict
from pathlib import Path

GLOBAL_CATEGORIES = {"GLOBAL"}
LOCAL_CATEGORIES = {"LOCAL"}
ALL_CATEGORIES = GLOBAL_CATEGORIES | LOCAL_CATEGORIES

REPEATS = 3
# The baseline probes are short by construction, so extra repeats are nearly
# free and they buy back precision on T_base -- which is a per-instance
# constant shared by every student, so noise in it is noise in the ranking.
BASELINE_REPEATS = 5
CLIP_LO, CLIP_HI = 0.1, 1_000_000.0   # effectively uncapped; see below
MEM_CAP_BYTES = 4 << 30               # the 4 GB rule, now actually enforced
DEFAULT_TIMEOUT = 3600.0

# The score is deliberately uncapped: on the scored instances a good
# submission is worth several hundred x, and a cap would tie the top of the
# field.  CLIP_HI survives only as a guard against a division blowing up on a
# near-zero measurement.

FOUNDATION_SRC = "dijkstra_foundation.cpp"
FOUNDATION_FLAGS = ["-O2", "-std=c++17"]


def platform_cxxflags():
    """Header-search fix for broken macOS Command Line Tools (see Makefile).

    Some CLT installs ship an almost-empty /usr/include/c++/v1, so <algorithm>
    is not found.  Fall back to the SDK's copy.  Adds no codegen flags, so the
    pinned -O2 baseline is unaffected.  No-op elsewhere.
    """
    if sys.platform != "darwin":
        return []
    clt = Path("/Library/Developer/CommandLineTools/usr/include/c++/v1/cstdint")
    if clt.exists():
        return []
    r = subprocess.run(["xcrun", "--sdk", "macosx", "--show-sdk-path"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return []
    inc = Path(r.stdout.strip()) / "usr/include/c++/v1"
    return ["-isystem", str(inc)] if (inc / "cstdint").exists() else []


# --------------------------------------------------------------- running ----

class RunResult:
    def __init__(self, wall, cpu, maxrss, status, stderr=""):
        self.wall = wall
        self.cpu = cpu
        self.maxrss = maxrss
        self.status = status        # "ok" | "timeout" | "crash" | "memory" | "threads" | "nooutput"
        self.stderr = stderr

    @property
    def ok(self):
        return self.status == "ok"


_RSS_UNIT = 1 if sys.platform == "darwin" else 1024   # ru_maxrss: bytes vs KiB


class _Footprint:
    """Live memory footprint of a child on macOS, via proc_pid_rusage.

    macOS refuses to lower RLIMIT_AS, and ru_maxrss is useless as a cap there:
    the memory compressor evicts pages as fast as a hog dirties them, so a
    5 GB allocation can show a 3 GB peak RSS.  phys_footprint is what the
    kernel itself uses for memory limits -- resident plus compressed pages --
    so grade.py samples it while polling and kills a child that crosses the
    cap.  Elsewhere this is a no-op and RLIMIT_AS does the job.
    """
    _V2, _SIZE, _OFF = 2, 160, 72          # RUSAGE_INFO_V2, sizeof, ri_phys_footprint

    def __init__(self):
        self.lib = None
        if sys.platform == "darwin":
            try:
                import ctypes, ctypes.util
                self.lib = ctypes.CDLL(ctypes.util.find_library("proc") or "libproc.dylib")
                self.buf = ctypes.create_string_buffer(self._SIZE)
            except OSError:
                self.lib = None

    def sample(self, pid):
        if self.lib is None:
            return 0
        if self.lib.proc_pid_rusage(pid, self._V2, self.buf) != 0:
            return 0
        return int.from_bytes(self.buf.raw[self._OFF:self._OFF + 8], "little")


_FOOTPRINT = _Footprint()


def time_run(binary, graph, queries, out, timeout, enforce_limits=True):
    """Run one binary once, measuring wall time, CPU time and peak RSS.

    fork + execv + wait4 rather than subprocess.run, because wait4 reports the
    rusage of *this* child.  `getrusage(RUSAGE_CHILDREN)` looks like it would
    do, but it is a session-wide high-water mark: after one large run every
    later run appears to have used just as much memory, so a delta measures
    nothing and the absolute value belongs to whichever child peaked first.
    """
    errfd_r, errfd_w = os.pipe()
    t0 = time.perf_counter()
    pid = os.fork()
    if pid == 0:                                   # child
        try:
            os.close(errfd_r)
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, 1)
            os.dup2(errfd_w, 2)
            os.close(errfd_w)
            if enforce_limits:
                try:
                    resource.setrlimit(resource.RLIMIT_AS,
                                       (MEM_CAP_BYTES, MEM_CAP_BYTES))
                except (OSError, ValueError):
                    # macOS rejects lowering RLIMIT_AS.  The peak-RSS check
                    # after the run still enforces the cap; a failure here
                    # must not turn every run into a "crash".
                    pass
            os.execv(str(binary), [str(binary), str(graph), str(queries), str(out)])
        except BaseException:
            os._exit(127)
    os.close(errfd_w)

    # Drain stderr on a helper thread so a chatty solver cannot fill the pipe
    # and deadlock against our own wait.
    chunks = []

    def drain():
        with os.fdopen(errfd_r, "rb") as f:
            chunks.append(f.read())

    t = threading.Thread(target=drain, daemon=True)
    t.start()

    deadline = t0 + timeout
    status = rusage = None
    killed = over_cap = False
    peak_footprint = 0
    while True:
        done, st, ru = os.wait4(pid, os.WNOHANG)
        if done:
            status, rusage = st, ru
            break
        now = time.perf_counter()
        if enforce_limits:
            peak_footprint = max(peak_footprint, _FOOTPRINT.sample(pid))
        if now > deadline or (enforce_limits and peak_footprint > MEM_CAP_BYTES):
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            _, _, rusage = os.wait4(pid, 0)
            killed = now > deadline
            over_cap = not killed
            break
        # Poll tightly at first so short probe runs are not rounded up, then
        # back off once the run is clearly a long one.
        time.sleep(0.0005 if now - t0 < 2.0 else 0.005)
    t1 = time.perf_counter()
    t.join(timeout=5)
    err = (chunks[0] if chunks else b"").decode(errors="replace")

    if killed:
        return RunResult(timeout, 0.0, 0, "timeout", err)

    wall = t1 - t0
    cpu = rusage.ru_utime + rusage.ru_stime
    rss = max(rusage.ru_maxrss * _RSS_UNIT, peak_footprint)
    if over_cap:
        return RunResult(wall, cpu, rss, "memory", err)
    exited_ok = os.WIFEXITED(status) and os.WEXITSTATUS(status) == 0

    if not exited_ok:
        # A memory-limited process usually dies on a failed allocation, so say
        # that rather than reporting a bare non-zero exit.
        near_cap = enforce_limits and rss >= MEM_CAP_BYTES * 0.9
        return RunResult(wall, cpu, rss, "memory" if near_cap else "crash", err)
    if enforce_limits and rss > MEM_CAP_BYTES:
        return RunResult(wall, cpu, rss, "memory", err)
    # A single-threaded program cannot burn appreciably more CPU than wall
    # time.  1.4x leaves room for measurement slop.
    if enforce_limits and wall > 1.0 and cpu > 1.4 * wall:
        return RunResult(wall, cpu, rss, "threads", err)
    return RunResult(wall, cpu, rss, "ok", err)


def best_of(binary, graph, queries, out, n, timeout, enforce_limits=True):
    """Best-of-n wall clock, but any failing repeat fails the instance.

    Output files are kept per repeat so nondeterminism across runs is visible.
    """
    best = None
    digests = set()
    for i in range(n):
        r = time_run(binary, graph, queries, out, timeout, enforce_limits)
        if not r.ok:
            return r, digests
        if not Path(out).exists():
            # Exit 0 but no output file.  Seen with a memory hog on macOS,
            # where the missing RLIMIT_AS lets it finish "successfully".
            return RunResult(r.wall, r.cpu, r.maxrss, "nooutput", r.stderr), digests
        digests.add(sha256_file(out))
        if best is None or r.wall < best.wall:
            best = r
    return best, digests


# ------------------------------------------------------------ correctness ----

def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def normalised_digest(path):
    """Digest of the answer sequence, insensitive to trailing whitespace.

    Streamed rather than built into two Python lists, which at hundreds of
    thousands of queries would materialise tens of MB per side.
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for line in f:
            line = line.strip()
            if line:
                h.update(line)
                h.update(b"\n")
    return h.hexdigest()


def answers_path(query_path):
    p = Path(str(query_path).replace(".queries", ".answers"))
    return p if p.exists() else None


# --------------------------------------------------------------- baseline ----

def build_foundation(workdir, src=FOUNDATION_SRC):
    """Compile the baseline ourselves, with pinned flags.

    The Makefile uses `CXXFLAGS ?=` and students submit their own Makefile, so
    a baseline built at -O0 would inflate every ratio.  With self-reported times
    that is the largest integrity hole in the whole scheme, and it costs one
    subprocess call to close.
    """
    src = Path(src).resolve()
    if not src.exists():
        sys.exit(f"cannot find {src} -- run grade.py from the repository root")
    out = Path(workdir) / "foundation_pinned"
    cxx = os.environ.get("CXX") or shutil.which("c++") or "g++"
    r = subprocess.run([cxx, *FOUNDATION_FLAGS, *platform_cxxflags(),
                        "-o", str(out), str(src)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"failed to build the foundation:\n{r.stderr}")
    return out, sha256_file(src)


def prefix_query_file(src, n, dst):
    """First n queries of a query file, as a standalone query file."""
    with open(src) as f, open(dst, "w") as g:
        total = int(f.readline().split()[0])
        n = min(n, total)
        g.write(f"{n}\n")
        for _ in range(n):
            g.write(f.readline())
    return n


def estimate_baseline(foundation, graph, queries, workdir, probe, timeout,
                      full=False):
    """(T_base, detail).  Fits T(q) = a + b*q from two prefixes unless --baseline-full."""
    with open(queries) as f:
        Q = int(f.readline().split()[0])

    if full or Q <= probe:
        out = Path(workdir) / "base.full"
        r, _ = best_of(foundation, graph, queries, out, BASELINE_REPEATS, timeout,
                       enforce_limits=False)
        if not r.ok:
            return None, {"mode": "full", "status": r.status}
        return r.wall, {"mode": "full", "Q": Q, "status": "ok"}

    n1, n2 = max(50, probe // 4), probe
    p1 = Path(workdir) / "probe1.queries"
    p2 = Path(workdir) / "probe2.queries"
    n1 = prefix_query_file(queries, n1, p1)
    n2 = prefix_query_file(queries, n2, p2)
    out = Path(workdir) / "base.probe"

    r1, _ = best_of(foundation, graph, p1, out, BASELINE_REPEATS, timeout,
                    enforce_limits=False)
    r2, _ = best_of(foundation, graph, p2, out, BASELINE_REPEATS, timeout,
                    enforce_limits=False)
    if not (r1.ok and r2.ok):
        return None, {"mode": "probe", "status": (r1.status, r2.status)}

    b = (r2.wall - r1.wall) / max(1, (n2 - n1))     # seconds per query
    a = r2.wall - b * n2                            # parse + startup
    if b <= 0:
        # The probe is too short to separate the two terms; scale the larger
        # measurement instead of extrapolating a negative slope.
        return r2.wall * Q / n2, {"mode": "probe-degenerate", "Q": Q,
                                  "n2": n2, "t2": r2.wall}
    return a + b * Q, {"mode": "probe", "Q": Q, "n1": n1, "n2": n2,
                       "t1": r1.wall, "t2": r2.wall,
                       "parse_s": a, "per_query_s": b}


# --------------------------------------------------------------- scoring ----

def geomean(xs):
    xs = [x for x in xs if x > 0]
    if not xs:
        return 0.0
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solver", required=True, type=Path)
    ap.add_argument("--instances", required=True, type=Path)
    ap.add_argument("--foundation", type=Path,
                    help="prebuilt baseline; by default grade.py compiles "
                         f"{FOUNDATION_SRC} itself with pinned flags")
    ap.add_argument("--baseline-probe", type=int, default=2000,
                    help="queries used to estimate T_base (default 2000)")
    ap.add_argument("--baseline-full", action="store_true",
                    help="run the foundation over every query instead of "
                         "extrapolating; only feasible on the dev tier")
    ap.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT)
    ap.add_argument("--no-limits", action="store_true",
                    help="skip the memory and thread checks (debugging only)")
    ap.add_argument("--json", type=Path, help="write machine-readable results here")
    args = ap.parse_args()

    args.solver = args.solver.resolve()

    instances = []
    for line in args.instances.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 3 or parts[0] not in ALL_CATEGORIES:
            sys.exit(f"bad instances line (category must be one of "
                     f"{sorted(ALL_CATEGORIES)}): {line}")
        instances.append((parts[0], Path(parts[1]), Path(parts[2])))

    missing = [str(p) for _, g, q in instances for p in (g, q) if not p.exists()]
    if missing:
        hint = ("\nThe scored set is not in the repo -- fetch it first with:  "
                "scripts/download_large.sh"
                if any("_large" in m for m in missing) else "")
        sys.exit("missing instance files:\n  " + "\n  ".join(missing) + hint)

    speedups_by_cat = defaultdict(list)
    rows = []
    records = []

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        if args.foundation:
            foundation, src_digest = args.foundation.resolve(), None
        else:
            foundation, src_digest = build_foundation(td)

        for i, (cat, gpath, qpath) in enumerate(instances):
            stu_out = td / f"{i}.stu"
            rec = {"category": cat, "graph": str(gpath), "queries": str(qpath)}

            t_stu_res, digests = best_of(args.solver, gpath, qpath, stu_out,
                                         REPEATS, args.timeout,
                                         enforce_limits=not args.no_limits)

            note, speedup, t_stu, t_base = None, CLIP_LO, None, None

            if not t_stu_res.ok:
                note = t_stu_res.status.upper()
                if t_stu_res.stderr.strip():
                    rec["stderr"] = t_stu_res.stderr.strip()[:2000]
            elif len(digests) > 1:
                note = "NONDETERMINISTIC"
            else:
                t_stu = t_stu_res.wall
                key = answers_path(qpath)
                got = normalised_digest(stu_out)
                if key is not None:
                    if got != normalised_digest(key):
                        note = "WRONG"
                else:
                    ref = td / f"{i}.ref"
                    fr, _ = best_of(foundation, gpath, qpath, ref, 1,
                                    args.timeout, enforce_limits=False)
                    if not fr.ok:
                        note = "NOKEY"
                    elif got != normalised_digest(ref):
                        note = "WRONG"

                if note is None:
                    t_base, detail = estimate_baseline(
                        foundation, gpath, qpath, td, args.baseline_probe,
                        args.timeout, full=args.baseline_full)
                    rec["baseline"] = detail
                    if t_base is None:
                        note = "BASEFAIL"
                    else:
                        raw = t_base / max(t_stu, 1e-9)
                        speedup = max(CLIP_LO, min(CLIP_HI, raw))
                        note = "ok"
                        rec["raw_speedup"] = raw

            # Every instance contributes, including the failures: a failed
            # instance scores the floor and drags the mean down rather than
            # vanishing from it.
            speedups_by_cat[cat].append(speedup)

            rows.append((cat, gpath.name, t_base, t_stu, t_stu_res.maxrss,
                         speedup, note))
            rec.update(t_base=t_base, t_solver=t_stu, speedup=speedup,
                       note=note, max_rss_bytes=t_stu_res.maxrss,
                       cpu_s=t_stu_res.cpu)
            records.append(rec)

    print(f"{'category':<8} {'instance':<26} {'T_base':>11} {'T_solver':>10} "
          f"{'peakRSS':>9} {'speedup':>10}  note")
    print("-" * 92)
    for cat, name, tb, ts, rss, sp, note in rows:
        tbs = f"{tb:11.3f}" if tb is not None else f"{'-':>11}"
        tss = f"{ts:10.3f}" if ts is not None else f"{'-':>10}"
        print(f"{cat:<8} {name:<26} {tbs} {tss} {rss / 2**20:8.0f}M "
              f"{sp:>10.3f}  {note}")

    print()
    for cat in sorted(ALL_CATEGORIES):
        if speedups_by_cat[cat]:
            print(f"{cat:<8} geomean = {geomean(speedups_by_cat[cat]):.3f}  "
                  f"(n={len(speedups_by_cat[cat])})")

    glob = sum((speedups_by_cat[c] for c in GLOBAL_CATEGORIES), [])
    loc = sum((speedups_by_cat[c] for c in LOCAL_CATEGORIES), [])
    overall = glob + loc

    print()
    print(f"Global-throughput track geomean = {geomean(glob):.3f}  (n={len(glob)})")
    print(f"Local-latency track     geomean = {geomean(loc):.3f}  (n={len(loc)})")
    print(f"Overall                 geomean = {geomean(overall):.3f}  (n={len(overall)})")

    if args.json:
        payload = {
            "solver": str(args.solver),
            "solver_sha256": sha256_file(args.solver),
            "foundation_source_sha256": src_digest,
            "foundation_flags": FOUNDATION_FLAGS,
            "repeats": REPEATS,
            "baseline_repeats": BASELINE_REPEATS,
            "baseline_probe": args.baseline_probe,
            "limits_enforced": not args.no_limits,
            "mem_cap_bytes": MEM_CAP_BYTES,
            "instances": records,
            "global_track_geomean": geomean(glob),
            "local_track_geomean": geomean(loc),
            "overall_geomean": geomean(overall),
        }
        args.json.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"\nwrote {args.json}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
