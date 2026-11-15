"""Summarise DATA/*.csv into RESULTS.json and FIGURES/cycle_scaling.png."""
import csv, glob, json, math, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
rows = []
for p in sorted(glob.glob(os.path.join(ROOT, "DATA", "*_Q*.csv"))):
    rows += list(csv.DictReader(open(p)))
for r in rows:
    for k in ("qbits", "run", "closed", "tail", "cycle", "rho"): r[k] = int(r[k])

def boot_ci(x, f=np.mean, n=4000, seed=1):
    g = np.random.default_rng(seed); x = np.asarray(x, float)
    b = [f(g.choice(x, len(x))) for _ in range(n)]
    return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]

out = {"hchm": [], "core": []}
for mode in ("hchm", "core"):
    for q in sorted({r["qbits"] for r in rows if r["mode"] == mode}):
        g = [r for r in rows if r["mode"] == mode and r["qbits"] == q]
        assert all(r["closed"] for r in g)
        N = 2.0 ** (8 * q); rm = math.sqrt(math.pi * N / 8)
        rho = [r["rho"] for r in g]; tail = [r["tail"] for r in g]; cyc = [r["cycle"] for r in g]
        e = {"qbits": q, "state_bits": 8 * q, "runs": len(g),
             "tail_mean": float(np.mean(tail)), "tail_median": float(np.median(tail)),
             "cycle_mean": float(np.mean(cyc)), "cycle_median": float(np.median(cyc)),
             "cycle_min": min(cyc), "cycle_max": max(cyc),
             "rho_mean": float(np.mean(rho)), "rho_mean_ci95": boot_ci(rho), "rho_median": float(np.median(rho)),
             "distinct_cycles_reached": len({r["cycle_id"] for r in g}),
             "distinct_cycle_lengths": sorted({r["cycle"] for r in g})[:12],
             "random_mapping_expected_tail_or_cycle": rm, "random_mapping_expected_rho": 2 * rm,
             "rho_mean_over_random_mapping": float(np.mean(rho)) / (2 * rm)}
        out[mode].append(e)

h = out["hchm"]
x = np.array([e["state_bits"] for e in h], float); y = np.log2([e["rho_mean"] for e in h])
fit = stats.linregress(x, y)
tcrit = stats.t.ppf(0.975, len(x) - 2)
out["hchm_scaling_fit"] = {"model": "log2(mean rho) = a + b * state_bits", "b": fit.slope, "a": fit.intercept,
    "b_ci95": [fit.slope - tcrit * fit.stderr, fit.slope + tcrit * fit.stderr], "r2": fit.rvalue ** 2,
    "points": len(x), "random_mapping_b": 0.5}
c = out["core"]
out["core_period_rule"] = {"observed": {str(e["qbits"]): e["cycle_max"] for e in c},
    "matches_3_times_2_pow_q_minus_1_for_q_ge_3": all(e["cycle_max"] == 3 * 2 ** (e["qbits"] - 1) for e in c if e["qbits"] >= 3)}
json.dump(out, open(os.path.join(ROOT, "RESULTS.json"), "w"), indent=2)

plt.figure(figsize=(7.4, 4.6))
xs = np.linspace(14, 66, 50)
plt.plot(xs, 0.5 * xs + math.log2(2 * math.sqrt(math.pi / 8)), "k--", lw=1, label=r"random-mapping expectation $2\sqrt{\pi N/8}$")
plt.errorbar(x, y, yerr=[y - np.log2([e["rho_mean_ci95"][0] for e in h]), np.log2([e["rho_mean_ci95"][1] for e in h]) - y],
             fmt="o", capsize=3, label="quantized HCHM-8: mean tail+cycle (95% bootstrap CI)")
plt.plot([e["state_bits"] for e in c], np.log2([e["cycle_max"] for e in c]), "s", label="hyperbolic core alone: maximum period")
plt.xlabel(r"state size $\log_2 N$ (8 coordinates $\times$ Q bits)"); plt.ylabel("$\log_2$ iterations until the orbit closes")
plt.legend(fontsize=8); plt.grid(alpha=.3); plt.tight_layout()
plt.savefig(os.path.join(ROOT, "FIGURES", "cycle_scaling.png"), dpi=220); plt.savefig(os.path.join(ROOT, "FIGURES", "cycle_scaling.pdf"))
for e in h: print("hchm Q%d N=2^%d runs=%d tail_med=%.0f cyc_med=%.0f rho_mean=%.4g ratio=%.2f distinct=%d lens=%s" % (e["qbits"], e["state_bits"], e["runs"], e["tail_median"], e["cycle_median"], e["rho_mean"], e["rho_mean_over_random_mapping"], e["distinct_cycles_reached"], e["distinct_cycle_lengths"]))
for e in c: print("core Q%d period min/max %d/%d distinct_lengths=%s" % (e["qbits"], e["cycle_min"], e["cycle_max"], e["distinct_cycle_lengths"]))
print(json.dumps(out["hchm_scaling_fit"], indent=1)); print(out["core_period_rule"])
