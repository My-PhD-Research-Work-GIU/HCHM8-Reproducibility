"""Rerun part of the uncapped cycle measurement and compare it, run by run, with the committed data.

Usage: python reproduce_check.py <qbits>:<profiles>:<runs> [...] [--perturb[=<scale>]] [--expect-different]
Example: python reproduce_check.py 2:10:32 3:10:32 4:10:32 5:3:32

For every (profile, start) pair the tail length, cycle length and cycle identifier must equal the row in
DATA/profiles_Q<qbits>.csv. With --perturb the rerun uses the last-bit perturbation mode of cycle_scaling_profiles.py
(scale 1, or the given scale). The script exits non-zero if any run differs; --expect-different inverts this and is
used for the sanity check in which a large perturbation must change the orbits. If the environment variable
REPRODUCE_REPORT is set, a JSON report is written to that path.
"""
import csv, json, os, platform, subprocess, sys, tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent; DATA = HERE.parent / "DATA"
specs = [a for a in sys.argv[1:] if not a.startswith("--")]
pert = next((a[2:] for a in sys.argv[1:] if a.startswith("--perturb")), None)
report = {"environment": platform.platform(), "python": platform.python_version(), "perturbation": pert or "none", "sizes": []}; bad = 0; total = 0
for spec in specs:
    q, nprof, runs = (int(v) for v in spec.split(":"))
    tmp = Path(tempfile.mkdtemp()) / ("Q%d.csv" % q)
    cmd = [sys.executable, str(HERE / "cycle_scaling_profiles.py"), str(q), "0", str(nprof), str(runs), str(tmp), str(os.cpu_count() or 2)] + ([pert] if pert else [])
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    key = lambda r: (r["profile"], r["run"]); val = lambda r: (r["tail"], r["cycle"], r["rho"], r["cycle_id"])
    ref = {key(r): val(r) for r in csv.DictReader(open(DATA / ("profiles_Q%d.csv" % q)))}
    new = {key(r): val(r) for r in csv.DictReader(open(tmp))}
    diff = [k for k in new if ref.get(k) != new[k]]
    report["sizes"].append({"qbits": q, "profiles": nprof, "starts_per_profile": runs, "runs_compared": len(new), "runs_differing": len(diff)}); bad += len(diff); total += len(new)
    print("Q%d: %d of %d runs differ from the committed data" % (q, len(diff), len(new)))
report["runs_compared"] = total; report["runs_differing"] = bad
out = os.environ.get("REPRODUCE_REPORT")
if out:
    with open(out, "w") as f: json.dump(report, f, indent=1)
if "--expect-different" in sys.argv: sys.exit(0 if bad else 1)
sys.exit(1 if bad else 0)
