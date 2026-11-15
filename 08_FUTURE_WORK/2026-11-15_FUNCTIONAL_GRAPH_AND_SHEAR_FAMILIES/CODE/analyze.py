"""Analysis of the registered runs of fg_study.py against the criteria of ../PLAN.md (written before the results were seen).

Reads DATA/RESULTS_fg_study_part1.json and DATA/RESULTS_fg_study_part2.json, writes DATA/SUMMARY.json and prints the verdicts.
Replication unit: the profile; intervals are percentile bootstraps over the ten profiles (10,000 resamples).

Usage: python analyze.py <DATA folder>
"""
import json, math, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

D = Path(sys.argv[1]).resolve(); rng = np.random.default_rng(20261115)
p1 = json.loads((D / "RESULTS_fg_study_part1.json").read_text())["part1"]
p2f = D / "RESULTS_fg_study_part2.json"
p2 = json.loads(p2f.read_text())["part2"] if p2f.exists() else []  # Part 1 can be analysed before Part 2 ends
POIS = [math.exp(-1) / math.factorial(k) for k in range(5)]


def boot(v):
    v = np.array(v, float); b = v[rng.integers(0, len(v), (10000, len(v)))].mean(axis=1)
    return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


def rows(q, cond, scale=None, fam=None):
    return [r for r in p1 if r.get("q") == q and r.get("condition") == cond and (scale is None or r["scale"] == scale) and (fam is None or r["family"] == fam)]


S = {}
# H1 implementation check
gch = [r["G_checks"] for r in p1 if "G_checks" in r]
S["H1"] = {"G_cases": len(gch), "all_pass": all(all(v is True or v == 0 for v in c.values()) for c in gch),
           "failures": [c for c in gch if not all(v is True or v == 0 for v in c.values())],
           "out_of_range_successors": sum(m.get("out_of_range", 0) for r in p1 if "maps" in r for m in r["maps"].values()),
           "C0_F_equals_A_every_state": all(m.get("equals_A_on_every_state") for r in p1 if r.get("condition") == "C0" for m in r["maps"].values())}


def indeg_freq(m, N):
    h = {int(k): v for k, v in m["indegree_hist"].items()}; return [h.get(k, 0) / N for k in range(5)]


# H2 in-degree of the emulator
for q in (2, 3):
    N = 2 ** (8 * q); rr = rows(q, "C1")
    fr = np.array([indeg_freq(r["maps"]["F_nearest"], N) for r in rr]); mean = fr.mean(axis=0)
    imf = [r["maps"]["F_nearest"]["image_fraction"] for r in rr]
    S["H2_q%d" % q] = {"indegree_freq_mean_k0_4": [float(x) for x in mean], "poisson": POIS, "max_abs_dev": float(np.max(np.abs(mean - POIS))),
                       "per_profile_max_abs_dev": [float(np.max(np.abs(f - POIS))) for f in fr], "image_fraction_mean": float(np.mean(imf)), "image_fraction_ci": boot(imf),
                       "max_indegree_range": [min(r["maps"]["F_nearest"]["max_indegree"] for r in rr), max(r["maps"]["F_nearest"]["max_indegree"] for r in rr)],
                       "states_2plus_pred_mean_fraction": float(np.mean([r["maps"]["F_nearest"]["states_with_2plus_predecessors"] / N for r in rr])),
                       "largest_basin_fraction_mean": float(np.mean([r["maps"]["F_nearest"]["largest_basin_fraction"] for r in rr])),
                       "cycles_mean": float(np.mean([r["maps"]["F_nearest"]["cycles"] for r in rr])),
                       "random_mapping_cycles": 0.5 * (math.log(2 * N) + 0.5772156649), "random_mapping_largest_component": 0.7578}
    if q == 3:
        S["H2_q3"]["criterion_met"] = bool(S["H2_q3"]["max_abs_dev"] <= 0.01 and 0.60 <= S["H2_q3"]["image_fraction_mean"] <= 0.66)


# H3 families at Q = 3 and slopes over Q = 2..5
def mean_rho_profile(q, cond, fam, scale, rule_name):
    if q <= 3:
        out = {}
        for r in p1:
            if r.get("q") == q and r.get("family") == fam and r.get("scale") == scale and rule_name in r.get("maps", {}) and r["condition"] in ("C1", "C2", "C3", "C4"):
                out[r["profile"]] = r["maps"][rule_name]["mean_rho"]
        return out
    acc = defaultdict(list)
    for r in p2:
        if r["q"] == q and r["condition"] == cond: acc[r["profile"]].append(r["rho"])
    return {p: float(np.mean(v)) for p, v in acc.items()}


fam_conds = {1: ("C1", "fam1_s1_nearest"), 2: ("C2", "fam2_s1_nearest"), 3: ("C3", "fam3_s1_nearest")}
S["H3"] = {}
for fam, (c1, c2) in fam_conds.items():
    rr = rows(3, c1, 1.0, fam)
    imf = [r["maps"]["F_nearest"]["image_fraction"] for r in rr]; rho = [r["maps"]["F_nearest"]["mean_rho_over_sqrt_piN_2"] for r in rr]
    per_q = {q: mean_rho_profile(q, c2, fam, 1.0, "F_nearest") for q in (2, 3, 4, 5)}
    slopes = []
    for p in range(10):
        xs = [8 * q for q in (2, 3, 4, 5) if p in per_q[q]]; ys = [math.log2(per_q[q][p]) for q in (2, 3, 4, 5) if p in per_q[q]]
        slopes.append(float(np.polyfit(xs, ys, 1)[0]) if len(xs) >= 4 else float("nan"))
    ratio = {q: float(np.mean([v / math.sqrt(math.pi * 2.0 ** (8 * q) / 2) for v in per_q[q].values()])) for q in per_q}
    S["H3"]["family_%d" % fam] = {"image_fraction_q3_mean": float(np.mean(imf)), "image_fraction_q3_ci": boot(imf), "rho_ratio_q3_mean": float(np.mean(rho)),
                                  "rho_ratio_q3_ci": boot(rho), "rho_ratio_by_q": ratio, "slope_mean": float(np.mean(slopes)), "slope_ci": boot(slopes),
                                  "slopes": slopes, "indegree_freq_q3_mean": [float(x) for x in np.mean([indeg_freq(r["maps"]["F_nearest"], 2 ** 24) for r in rr], axis=0)],
                                  "cycles_q3_mean": float(np.mean([r["maps"]["F_nearest"]["cycles"] for r in rr])),
                                  "criteria_met": bool(0.60 <= np.mean(imf) <= 0.66 and 0.7 <= np.mean(rho) <= 1.3 and 0.45 <= np.mean(slopes) <= 0.55)}
# H4 rounding rule
rr = rows(3, "C1")
a = [r["maps"]["F_nearest"]["image_fraction"] for r in rr]; b = [r["maps"]["F_truncation"]["image_fraction"] for r in rr]
ra = [r["maps"]["F_nearest"]["mean_rho_over_sqrt_piN_2"] for r in rr]; rb = [r["maps"]["F_truncation"]["mean_rho_over_sqrt_piN_2"] for r in rr]
trunc_q = {q: mean_rho_profile(q, "fam1_s1_truncation", 1, 1.0, "F_truncation") for q in (4, 5)}
S["H4"] = {"image_fraction_nearest": float(np.mean(a)), "image_fraction_truncation": float(np.mean(b)), "difference": float(np.mean(b) - np.mean(a)),
           "rho_ratio_nearest": float(np.mean(ra)), "rho_ratio_truncation": float(np.mean(rb)),
           "rho_ratio_truncation_q4_q5": {q: float(np.mean([v / math.sqrt(math.pi * 2.0 ** (8 * q) / 2) for v in d.values()])) for q, d in trunc_q.items()},
           "criterion_met": bool(abs(np.mean(b) - np.mean(a)) < 0.01 and 0.7 <= np.mean(ra) <= 1.3 and 0.7 <= np.mean(rb) <= 1.3)}
# exploratory: strength series, stepwise discretization
st = {}
for q in (2, 3):
    st[q] = {}
    c0 = rows(q, "C0"); st[q]["0"] = {"image_fraction": float(np.mean([r["maps"]["F_nearest"]["image_fraction"] for r in c0])),
                                      "rho_ratio": float(np.mean([r["maps"]["F_nearest"]["mean_rho_over_sqrt_piN_2"] for r in c0])),
                                      "cycles": float(np.mean([r["maps"]["F_nearest"]["cycles"] for r in c0]))}
    for sc, cond in ((0.25, "C4"), (0.5, "C4"), (1.0, "C1"), (2.0, "C4")):
        rr = rows(q, cond, sc, 1)
        st[q][str(sc)] = {"image_fraction": float(np.mean([r["maps"]["F_nearest"]["image_fraction"] for r in rr])),
                          "image_fraction_range": [min(r["maps"]["F_nearest"]["image_fraction"] for r in rr), max(r["maps"]["F_nearest"]["image_fraction"] for r in rr)],
                          "rho_ratio": float(np.mean([r["maps"]["F_nearest"]["mean_rho_over_sqrt_piN_2"] for r in rr])),
                          "cycles": float(np.mean([r["maps"]["F_nearest"]["cycles"] for r in rr]))}
S["strength_series"] = st
sp = {}
for q in (2, 3):
    N = 2 ** (8 * q); HN = math.log(N) + 0.5772156649
    rr = rows(q, "C1")
    sp[q] = {"S_nearest_equals_G_fraction_min": min(r["S_nearest_equals_G_fraction"] for r in rr), "F_equals_G_fraction_max": max(r["F_equals_G_fraction"] for r in rr)}
    for fam, cond in ((1, "C1"), (2, "C2"), (3, "C3")):
        r2 = rows(q, cond, 1.0, fam); m = [r["maps"]["S_truncation"] for r in r2]
        sp[q]["S_truncation_family_%d" % fam] = {"all_permutations": all(x["is_permutation"] for x in m), "cycles_mean": float(np.mean([x["cycles"] for x in m])),
                                                 "cycles_range": [min(x["cycles"] for x in m), max(x["cycles"] for x in m)], "random_permutation_cycles": HN,
                                                 "mean_cycle_length_of_random_state_over_N_plus_1_half": float(np.mean([x["mean_cycle_length_of_random_state"] / ((N + 1) / 2) for x in m])),
                                                 "longest_cycle_fraction_mean": float(np.mean([x["longest_cycle"] / N for x in m]))}
S["stepwise"] = sp
# reproduction of the registered revision 9 run (same F_Q and G_Q code path for family 1, s = 1)
old = json.loads((D.parents[1] / "2026-11-10_GRID_CONJUGACY" / "RESULTS_grid_conjugacy.json").read_text())
rep = {"cases": 0, "identical": 0, "differences": []}
for o in old["exhaustive"]:
    if "profile" not in o: continue
    new = [r for r in p1 if r.get("q") == o["qbits"] and r.get("profile") == o["profile"] and r.get("condition") == "C1"]
    if not new: continue
    f = new[0]["maps"]["F_nearest"]; g = new[0]["maps"]["G"]; rep["cases"] += 1
    same = (abs(f["image_fraction"] - o["F"]["image_fraction"]) < 1e-12 and f["cycles"] == o["F"]["cycles"] and abs(f["mean_rho"] - o["F"]["mean_rho"]) < 1e-9
            and f["max_tail"] == o["F"]["max_tail"] and f["longest_cycle"] == o["F"]["longest_cycle"] and g["cycles"] == o["G"]["cycles"]
            and abs(new[0]["F_equals_G_fraction"] - o["F_equals_G_fraction"]) < 1e-12)
    rep["identical"] += int(same)
    if not same: rep["differences"].append({"q": o["qbits"], "profile": o["profile"]})
S["reproduction_of_rev9"] = rep
(D / "SUMMARY.json").write_text(json.dumps(S, indent=1))
print("reproduction of rev9:", rep)
print(json.dumps({"H1": {k: v for k, v in S["H1"].items() if k != "failures"}, "H2_q3": {k: S["H2_q3"][k] for k in ("indegree_freq_mean_k0_4", "max_abs_dev", "image_fraction_mean", "criterion_met", "max_indegree_range", "largest_basin_fraction_mean", "cycles_mean")}}, indent=1))
for k, v in S["H3"].items(): print(k, {x: (round(y, 4) if isinstance(y, float) else y) for x, y in v.items() if x not in ("slopes", "indegree_freq_q3_mean", "rho_ratio_by_q")}, {q: round(y, 3) for q, y in v["rho_ratio_by_q"].items()})
print("H4", S["H4"])
print("strength", json.dumps(st))
print("stepwise", json.dumps(sp))
