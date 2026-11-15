"""Distinct cycles reached from l starting states, against the prediction of Grebogi, Ott and Yorke (1988), and the total number of
cycles against the random-mapping value of Kruskal (1954) as quoted there (plan: ../PLAN.md, written before this was computed).

Data (stored, not rerun): the uncapped orbits of the rounded emulator, 08_FUTURE_WORK/2026-10-11_FINITE_STATE_CYCLE_SCALING/DATA/
profiles_Q{2..7}.csv (one row per profile and start; cycle_id identifies the cycle reached), and the exhaustive graphs of
08_FUTURE_WORK/2026-11-10_GRID_CONJUGACY/RESULTS_grid_conjugacy.json (number of cycles of F_Q at Q = 2 and 3).

Prediction for l starts on a mixing system with roundoff, valid for l <p> << 1: N(l) = sum_{j=1}^{l} 1/(2j - 1).
Total cycles of a random mapping on N points: (1/2) ln N + gamma asymptotically; the exact expectation sum_k (N)_k/(k N^k) is also given.
Usage: python distinct_cycles.py <out_dir>
"""
import csv, json, math, sys, time
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
FW = HERE.parents[2] / "08_FUTURE_WORK"
D = FW / "2026-10-11_FINITE_STATE_CYCLE_SCALING" / "DATA"
GRID = FW / "2026-11-10_GRID_CONJUGACY" / "RESULTS_grid_conjugacy.json"
GAMMA = 0.5772156649015329


def grebogi(l):
    return sum(1.0 / (2 * j - 1) for j in range(1, l + 1))


def exact_cycles(n):
    total, logp, k = 0.0, 0.0, 1
    while True:
        total += math.exp(logp) / k
        logp += math.log1p(-k / n); k += 1
        if logp < -60 or k > n:
            return total


out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
res = {"plan": "../PLAN.md", "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "distinct": [], "total": []}
for q in range(2, 8):
    by = defaultdict(set); starts = defaultdict(int)
    with open(D / ("profiles_Q%d.csv" % q)) as f:
        for r in csv.DictReader(f):
            by[int(r["profile"])].add(r["cycle_id"]); starts[int(r["profile"])] += 1
    ls = sorted(set(starts.values())); assert len(ls) == 1, (q, ls); l = ls[0]
    counts = [len(by[p]) for p in sorted(by)]
    mean = sum(counts) / len(counts); pred = grebogi(l)
    res["distinct"].append({"qbits": q, "state_bits": 8 * q, "profiles": len(counts), "starts_per_profile": l, "distinct_per_profile": counts,
                            "mean": mean, "min": min(counts), "max": max(counts), "grebogi_prediction": pred, "ratio": mean / pred,
                            "registered_band_met": 0.67 <= mean / pred <= 1.5 if 3 <= q <= 6 else None})
G = json.loads(GRID.read_text())
for q in (2, 3):
    c = [r["F"]["cycles"] for r in G["exhaustive"] if r.get("qbits") == q and "F" in r]
    n = 2 ** (8 * q)
    res["total"].append({"qbits": q, "cycles_per_profile": c, "mean": sum(c) / len(c), "random_exact": exact_cycles(n),
                         "kruskal_asymptotic": 0.5 * math.log(n) + GAMMA})
res["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
(out / "RESULTS_distinct_cycles.json").write_text(json.dumps(res, indent=1))
for r in res["distinct"]:
    print("Q%d: l=%d, mean distinct %.2f (range %d to %d) vs Grebogi %.2f, ratio %.2f, band %s" % (r["qbits"], r["starts_per_profile"], r["mean"], r["min"], r["max"], r["grebogi_prediction"], r["ratio"], r["registered_band_met"]))
for r in res["total"]:
    print("Q%d all cycles: mean %.1f vs random exact %.2f (Kruskal %.2f)" % (r["qbits"], r["mean"], r["random_exact"], r["kruskal_asymptotic"]))
