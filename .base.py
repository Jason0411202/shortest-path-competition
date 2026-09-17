import sys, tempfile, json
sys.path.insert(0, '/root/spc')
import grade
from pathlib import Path
res = {}
with tempfile.TemporaryDirectory() as td:
    f, _ = grade.build_foundation(Path(td))
    for n in ['hugeq','local2d','road2d','scalefree','wide64','lattice3d']:
        g = Path(f'/root/spc/instances/{n}_large.graph'); q = Path(f'/root/spc/instances/{n}_large.queries')
        tb, det = grade.estimate_baseline(f, g, q, Path(td), 2000, 3600)
        res[n] = tb
        print(n, tb, det, flush=True)
json.dump(res, open('/root/spc/base_large.json','w'))
