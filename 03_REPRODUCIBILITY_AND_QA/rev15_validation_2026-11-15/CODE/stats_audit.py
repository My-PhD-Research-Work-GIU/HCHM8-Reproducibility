"""Statistical audit of paper revision 14 from the raw records (revision 15 validation, 2026-11-15).

A. Ablation (08_FUTURE_WORK/2026-11-10_CHAOS_ONLY_ARM/DATA/chaos_only_metrics.csv, 4 arms x 4,096 tables = 256 keys x 16 rounds):
   recomputes the paper's tests (KS / Mann-Whitney, Holm within a comparison, Hedges' g with the two-sample normal interval);
   adds what the paper's intervals do not use: the dependence of tables of one key (intraclass correlation over keys), a
   cluster bootstrap over keys of the standardized paired difference, and the 90% interval bound that an equivalence claim at
   level 0.05 would need (reported as information only; no margin was fixed in advance, so no equivalence is claimed).
B. Orbit scaling (08_FUTURE_WORK/2026-10-11_FINITE_STATE_CYCLE_SCALING/DATA/profiles_Q2..Q7.csv): profile means of tail plus
   cycle, the fit of log2(mean) on state bits over 16 to 48 bits, a cluster bootstrap over profiles, residuals, the fitted range,
   other summaries (median, mean cycle, mean tail), and the 56-bit prediction (mean ratio against 2 sqrt(pi N / 8)).
C. Table 11 (08_FUTURE_WORK/2026-11-10_DERIVATION_COST/runs.csv, runs_per_table_ns.csv): every entry recomputed from the raw
   timings, with the definitions written out, and the relation between throughput and median time per table.

Usage: python stats_audit.py <08_FUTURE_WORK folder> <tab_cost.tex> <out.json>
"""
import csv, json, math, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy import stats

FW, TABCOST, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
rng = np.random.default_rng(20261115)
R = {}

# ------------------------------------------------------------------------------------------------ A. ablation
rows = list(csv.DictReader(open(FW / "2026-11-10_CHAOS_ONLY_ARM/DATA/chaos_only_metrics.csv")))
STATS = [("sac_mad", "ks"), ("sac_maxdev", "ks"), ("bic_max_abs_corr", "ks"), ("bic_mean_abs_corr", "ks"),
         ("fixed", "mw"), ("opposite", "mw"), ("cycles", "mw"), ("longest_cycle", "mw")]
arms = defaultdict(dict)
for r in rows:
    arms[r["arm"]][(int(r["key_index"]), int(r["round"]))] = r
keys = sorted(arms["S"])
R["ablation_counts"] = {a: len(v) for a, v in arms.items()}
R["ablation_distinct_tables"] = {a: len({r["sha256"] for r in v.values()}) for a, v in arms.items()}


def hedges(x, y):
    nx, ny = len(x), len(y); sp = math.sqrt(((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / (nx + ny - 2))
    d = (x.mean() - y.mean()) / sp; J = 1 - 3 / (4 * (nx + ny) - 9); g = J * d
    se = math.sqrt((nx + ny) / (nx * ny) + g * g / (2 * (nx + ny))); return g, se, sp


def icc_by_key(v):  # one-way ANOVA ICC(1) over 256 keys x 16 rounds
    a = np.array(v).reshape(256, 16); k = 16; gm = a.mean()
    msb = k * ((a.mean(axis=1) - gm) ** 2).sum() / 255; msw = ((a - a.mean(axis=1, keepdims=True)) ** 2).sum() / (256 * 15)
    return (msb - msw) / (msb + (k - 1) * msw)


comps = {}
for x_name, y_name in (("H", "S"), ("C", "S"), ("CL", "S"), ("C", "H")):
    tests = []
    for st, kind in STATS:
        x = np.array([float(arms[x_name][k][st]) for k in keys]); y = np.array([float(arms[y_name][k][st]) for k in keys])
        p = float(stats.ks_2samp(x, y).pvalue) if kind == "ks" else float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
        g, se, sp = hedges(x, y)
        diff = (x - y).reshape(256, 16)
        boot = []
        for _ in range(4000):
            idx = rng.integers(0, 256, 256); boot.append(diff[idx].mean() / sp)
        boot = np.array(boot)
        tests.append({"stat": st, "test": kind, "p": p, "g": g, "g_ci95_normal": [g - 1.96 * se, g + 1.96 * se],
                      "g_ci90_normal": [g - 1.645 * se, g + 1.645 * se],
                      "paired_std_diff": float((x - y).mean() / sp), "paired_std_diff_ci95_cluster_keys": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                      "paired_r": float(np.corrcoef(x, y)[0, 1]), "icc_keys_x": float(icc_by_key(x)), "icc_keys_y": float(icc_by_key(y))})
    order = sorted(range(len(tests)), key=lambda i: tests[i]["p"]); run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, tests[i]["p"] * (len(tests) - rank))); tests[i]["p_holm"] = run
    comps["%s_vs_%s" % (x_name, y_name)] = {
        "tests": tests, "smallest_p": min(t["p"] for t in tests), "smallest_holm": min(t["p_holm"] for t in tests),
        "max_abs_g": max(abs(t["g"]) for t in tests), "g_bound95_normal": max(max(abs(t["g_ci95_normal"][0]), abs(t["g_ci95_normal"][1])) for t in tests),
        "g_bound90_normal": max(max(abs(t["g_ci90_normal"][0]), abs(t["g_ci90_normal"][1])) for t in tests),
        "paired_bound95_cluster": max(max(abs(t["paired_std_diff_ci95_cluster_keys"][0]), abs(t["paired_std_diff_ci95_cluster_keys"][1])) for t in tests),
        "max_abs_icc": max(max(abs(t["icc_keys_x"]), abs(t["icc_keys_y"])) for t in tests), "max_abs_paired_r": max(abs(t["paired_r"]) for t in tests)}
R["ablation"] = comps

# ------------------------------------------------------------------------------------------------ B. orbit scaling
D = FW / "2026-10-11_FINITE_STATE_CYCLE_SCALING/DATA"
orb = defaultdict(lambda: defaultdict(list))
for q in range(2, 8):
    for r in csv.DictReader(open(D / ("profiles_Q%d.csv" % q))):
        orb[q][int(r["profile"])].append((int(r["tail"]), int(r["cycle"]), int(r["rho"])))
summ = {}
for q in orb:
    N = 2.0 ** (8 * q); rm = 2 * math.sqrt(math.pi * N / 8)
    pm = {p: np.mean([x[2] for x in v]) for p, v in orb[q].items()}
    summ[q] = {"profiles": len(pm), "starts_per_profile": sorted({len(v) for v in orb[q].values()}), "profile_mean_ratio": {p: m / rm for p, m in pm.items()},
               "mean_ratio": float(np.mean([m / rm for m in pm.values()])), "censored_runs": 0}
R["orbit_summary"] = {q: {k: v for k, v in s.items() if k != "profile_mean_ratio"} for q, s in summ.items()}


def fit(qs, how="mean", which=2):
    X, Y, P = [], [], []
    for q in qs:
        for p, v in orb[q].items():
            vals = [x[which] for x in v]
            m = np.mean(vals) if how == "mean" else np.median(vals)
            X.append(8 * q); Y.append(math.log2(m)); P.append(p)
    X, Y, P = np.array(X, float), np.array(Y), np.array(P)
    res = stats.linregress(X, Y)
    profs = sorted(set(P)); bs = []
    for _ in range(4000):
        pick = rng.choice(profs, len(profs)); xi = np.concatenate([X[P == p] for p in pick]); yi = np.concatenate([Y[P == p] for p in pick])
        bs.append(stats.linregress(xi, yi).slope)
    resid = Y - (res.intercept + res.slope * X)
    return {"points": int(len(X)), "slope": float(res.slope), "ols_ci95": [float(res.slope - 1.96 * res.stderr), float(res.slope + 1.96 * res.stderr)],
            "cluster_boot_ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], "r2": float(res.rvalue ** 2),
            "residual_sd": float(resid.std(ddof=2)), "max_abs_residual_log2": float(np.abs(resid).max()),
            "mean_residual_by_bits": {int(b): float(resid[X == b].mean()) for b in sorted(set(X))}}


R["scaling_fit_16_48_mean_rho"] = fit([2, 3, 4, 5, 6])
R["scaling_fit_sensitivity"] = {"16_40": fit([2, 3, 4, 5]), "24_48": fit([3, 4, 5, 6]), "32_48": fit([4, 5, 6]),
                                "median_rho_16_48": fit([2, 3, 4, 5, 6], "median"), "mean_cycle_16_48": fit([2, 3, 4, 5, 6], "mean", 1),
                                "mean_tail_16_48": fit([2, 3, 4, 5, 6], "mean", 0)}
f = R["scaling_fit_16_48_mean_rho"]
pred56 = {"mean_ratio": summ[7]["mean_ratio"], "profile_ratios": summ[7]["profile_mean_ratio"], "band": [0.7, 1.3],
          "met": 0.7 <= summ[7]["mean_ratio"] <= 1.3}
R["prediction_56_bits"] = pred56

# ------------------------------------------------------------------------------------------------ C. Table 11
runs = list(csv.DictReader(open(FW / "2026-11-10_DERIVATION_COST/runs.csv")))
per = defaultdict(list)
for r in csv.DictReader(open(FW / "2026-11-10_DERIVATION_COST/runs_per_table_ns.csv")):
    per[(int(r["sequence"]), r["arm"])].append((int(r["material_ns"]), int(r["expand_ns"]), int(r["table_ns"])))
tab = {}
for arm in ("S", "H", "L"):
    rr = [r for r in runs if r["arm"] == arm]
    thr, med_ms, mean_ms, share, overhead, peak, trace = [], [], [], [], [], [], []
    for r in rr:
        t = per[(int(r["sequence"]), arm)]; tot = np.array([a + b + c for a, b, c in t]); mat = np.array([a for a, _, _ in t])
        thr.append(int(r["tables"]) / (int(r["loop_ns"]) * 1e-9)); med_ms.append(np.median(tot) / 1e6); mean_ms.append(tot.mean() / 1e6)
        share.append(np.median(mat) / np.median(tot)); overhead.append(1 - tot.sum() / int(r["loop_ns"]))
        peak.append(int(r["peak_after_loop"]) / 2 ** 20); trace.append(int(r["tracemalloc_peak_16"]) / 1024)

    def bci(v, fn=np.median):
        v = np.array(v); b = [fn(v[rng.integers(0, len(v), len(v))]) for _ in range(10000)]
        return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]
    tab[arm] = {"runs": len(rr), "tables_per_run": sorted({int(r["tables"]) for r in rr}), "all_reference_checks": all(r["reference_check"] == "True" for r in rr),
                "throughput_median": float(np.median(thr)), "throughput_ci95": bci(thr), "time_per_table_ms_median_of_run_medians": float(np.median(med_ms)),
                "time_ci95": bci(med_ms), "time_per_table_ms_median_of_run_means": float(np.median(mean_ms)),
                "reciprocal_of_median_throughput_ms": 1000 / float(np.median(thr)), "loop_overhead_outside_stage_timers_median": float(np.median(overhead)),
                "material_share_median_ratio_of_medians": float(np.median(share)), "peak_working_set_MiB_median": float(np.median(peak)),
                "tracemalloc_KiB_median": float(np.median(trace))}
R["table11_recomputed"] = tab
R["table11_as_printed"] = [l for l in TABCOST.read_text(encoding="utf-8").splitlines() if "&" in l]
json.dump(R, open(OUT, "w"), indent=1, default=str)
for k in ("H_vs_S", "C_vs_S", "CL_vs_S", "C_vs_H"):
    c = comps[k]; print(k, {x: round(c[x], 4) for x in ("smallest_p", "smallest_holm", "max_abs_g", "g_bound95_normal", "g_bound90_normal", "paired_bound95_cluster", "max_abs_icc", "max_abs_paired_r")})
print("fit 16-48:", {k: v for k, v in R["scaling_fit_16_48_mean_rho"].items() if k != "mean_residual_by_bits"})
for k, v in R["scaling_fit_sensitivity"].items(): print(" ", k, round(v["slope"], 4), [round(x, 4) for x in v["cluster_boot_ci95"]])
print("residual by bits:", {k: round(v, 3) for k, v in R["scaling_fit_16_48_mean_rho"]["mean_residual_by_bits"].items()})
print("56-bit:", round(pred56["mean_ratio"], 3), pred56["met"])
for a, v in tab.items(): print("T11", a, {k: (round(x, 4) if isinstance(x, float) else x) for k, x in v.items() if not isinstance(x, list)}, [round(x, 3) for x in v["throughput_ci95"]], [round(x, 4) for x in v["time_ci95"]])
print("\n".join(R["table11_as_printed"]))
