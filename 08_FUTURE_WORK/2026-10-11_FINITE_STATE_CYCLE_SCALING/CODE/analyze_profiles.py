"""Summarise DATA/profiles_Q*.csv into RESULTS_PROFILES.json.

The parameter profile is the unit of analysis: starts within one profile share a functional graph,
so each profile contributes one value (its mean tail + cycle) to the scaling fit.

The fit uses 16 to 48 state bits only. The 56-bit runs were made afterwards as an out-of-sample test of the fit
(see PLAN.md) and are reported separately; they never enter the fit.
"""
import csv, glob, json, math, os
import numpy as np
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = []
for p in sorted(glob.glob(os.path.join(ROOT, "DATA", "profiles_Q*.csv"))):
    rows += [{k: (v if k == "cycle_id" else float(v) if k == "seconds" else int(v)) for k, v in r.items()} for r in csv.DictReader(open(p))]

FIT_MAX_Q = 6
out = {"unit_of_analysis": "parameter profile", "fit_uses_qbits_up_to": FIT_MAX_Q, "per_size": []}
fx, fy = [], []
for q in sorted({r["qbits"] for r in rows}):
    g = [r for r in rows if r["qbits"] == q]; profs = sorted({r["profile"] for r in g})
    N = 2.0 ** (8 * q); rm = 2 * math.sqrt(math.pi * N / 8)
    pm = [float(np.mean([r["rho"] for r in g if r["profile"] == k])) for k in profs]
    cyc = [r["cycle"] for r in g]; tails = [r["tail"] for r in g]
    out["per_size"].append({
        "qbits": q, "state_bits": 8 * q, "in_fit": q <= FIT_MAX_Q, "profiles": len(profs), "starts_per_profile": len(g) // len(profs), "runs": len(g),
        "profile_mean_rho": pm, "geomean_profile_mean_rho": float(np.exp(np.mean(np.log(pm)))),
        "random_mapping_expected_rho": rm,
        "ratio_to_random_mapping_min_median_max": [min(pm) / rm, float(np.median(pm)) / rm, max(pm) / rm],
        "tail_median": float(np.median(tails)), "cycle_median": float(np.median(cyc)),
        "shortest_cycle": min(cyc), "longest_cycle": max(cyc),
        "shortest_cycle_per_profile": [min(r["cycle"] for r in g if r["profile"] == k) for k in profs],
        "distinct_cycles_reached_per_profile": [len({r["cycle_id"] for r in g if r["profile"] == k}) for k in profs],
        "share_of_starts_on_cycles_shorter_than_sqrtN_over_100": float(np.mean([c < math.sqrt(N) / 100 for c in cyc])),
        "all_cycle_lengths": cyc})
    if q <= FIT_MAX_Q:
        fx += [8 * q] * len(pm); fy += list(np.log2(pm))
fit = stats.linregress(fx, fy); t = stats.t.ppf(0.975, len(fx) - 2)
out["scaling_fit"] = {"model": "log2(profile mean rho) = a + b * state_bits", "points": len(fx), "b": fit.slope, "a": fit.intercept,
                      "b_ci95": [fit.slope - t * fit.stderr, fit.slope + t * fit.stderr], "r2": fit.rvalue ** 2, "random_mapping_b": 0.5,
                      "random_mapping_a": math.log2(2 * math.sqrt(math.pi / 8))}
for e in out["per_size"]:
    if not e["in_fit"]:
        pred = fit.intercept + fit.slope * e["state_bits"]; ratios = [v / e["random_mapping_expected_rho"] for v in e["profile_mean_rho"]]
        out["out_of_sample_%d_bits" % e["state_bits"]] = {
            "prediction_fit_log2": pred, "prediction_random_mapping_log2": math.log2(e["random_mapping_expected_rho"]),
            "observed_log2_mean_of_profile_means": math.log2(float(np.mean(e["profile_mean_rho"]))),
            "observed_mean_log2_profile_mean": float(np.mean(np.log2(e["profile_mean_rho"]))),
            "ratio_to_random_mapping_mean_over_profiles": float(np.mean(ratios)), "ratio_to_random_mapping_median": float(np.median(ratios)),
            "ratio_min_max": [min(ratios), max(ratios)], "prespecified_interval_for_mean_ratio": [0.7, 1.3],
            "criterion_met": bool(0.7 <= float(np.mean(ratios)) <= 1.3), "profiles": e["profiles"], "starts_per_profile": e["starts_per_profile"]}
json.dump(out, open(os.path.join(ROOT, "RESULTS_PROFILES.json"), "w"), indent=1)
for e in out["per_size"]:
    print("Q%d  profiles=%d x %d  tail_med=%.0f  cyc_med=%.0f  shortest=%d  ratio(min/med/max)=%.2f/%.2f/%.2f  short-cycle share=%.2f  distinct/profile=%s" % (
        e["qbits"], e["profiles"], e["starts_per_profile"], e["tail_median"], e["cycle_median"], e["shortest_cycle"],
        *e["ratio_to_random_mapping_min_median_max"], e["share_of_starts_on_cycles_shorter_than_sqrtN_over_100"], e["distinct_cycles_reached_per_profile"]))
print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in out["scaling_fit"].items() if k != "b_ci95"}, [round(v, 4) for v in out["scaling_fit"]["b_ci95"]])
