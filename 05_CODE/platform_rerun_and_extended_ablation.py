"""Rerun the published experiments in the current environment, compare with the stored results, and extend the ablation.

Part A  map numerics: reruns run_map_tests() of run_validated_experiments.py and compares with 04_DATASETS/map_validated_summary.json.
Part B  portability: regenerates the 8,192-table corpus and compares every table's SHA-256 with 04_DATASETS/dynamic_family_direct_8192.csv,
        which was produced in the original (Linux) environment.
Part C  extended ablation, fixed before running: for all 4,096 tables per arm, eight statistics that are not affine-invariant
        (SAC mean absolute deviation, SAC maximum deviation, BIC maximum and mean absolute correlation, fixed points,
        opposite fixed points, number of cycles, longest cycle). Continuous statistics: two-sample Kolmogorov-Smirnov test;
        counts: Mann-Whitney U test; Holm correction over the eight tests.

The original module is imported from a temporary copy so that the directories it creates do not land in the repository.

Usage: python platform_rerun_and_extended_ablation.py <tag> <output_dir>
"""
import csv, hashlib, importlib, json, platform, shutil, sys, tempfile, time
from pathlib import Path
import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
tag, out = sys.argv[1], Path(sys.argv[2]); (out / "DATA").mkdir(parents=True, exist_ok=True)
tmp = Path(tempfile.mkdtemp(prefix="r01_")); (tmp / "code").mkdir()
shutil.copy(HERE / "run_validated_experiments.py", tmp / "code" / "run_validated_experiments.py")
sys.path.insert(0, str(tmp / "code"))
import matplotlib; matplotlib.use("Agg")
m = importlib.import_module("run_validated_experiments")

res = {"tag": tag, "environment": {"platform": platform.platform(), "machine": platform.machine(), "python": platform.python_version(),
                                   "numpy": np.__version__, "scipy": __import__("scipy").__version__, "numba": __import__("numba").__version__}}

# ---------------------------------------------------------------- Part A
t0 = time.time(); new = m.run_map_tests(); old = json.load(open(ROOT / "04_DATASETS" / "map_validated_summary.json"))
exact = np.array(old["exact_lyapunov_spectrum"])
res["map_numerics"] = {
    "numeric_le_mean_stored": old["numeric_le_mean"], "numeric_le_mean_rerun": [float(v) for v in new["numeric_le_mean"]],
    "max_abs_numeric_minus_exact_stored": float(np.max(np.abs(np.array(old["numeric_le_mean"]) - exact))),
    "max_abs_numeric_minus_exact_rerun": float(np.max(np.abs(np.array(new["numeric_le_mean"]) - exact))),
    "max_abs_rerun_minus_stored_mean": float(np.max(np.abs(np.array(new["numeric_le_mean"]) - np.array(old["numeric_le_mean"])))),
    "numeric_le_identical": bool(np.array_equal(np.array(new["numeric_le_mean"]), np.array(old["numeric_le_mean"]))),
    "max_roundtrip_error": [old["max_roundtrip_error_2000"], float(new["max_roundtrip_error_2000"])],
    "max_conjugacy_error": [old["max_conjugacy_error_2000"], float(new["max_conjugacy_error_2000"])],
    "detT_max_abs_error": [old["detT_max_abs_error_from_1"], float(new["detT_max_abs_error_from_1"])],
    "reference_sensitivity_thresholds_stored": [old["reference_sensitivity"][k] for k in ("first_gt_1e-12", "first_gt_1e-9", "first_gt_1e-6", "first_gt_1e-3")],
    "reference_sensitivity_thresholds_rerun": [int(new["reference_sensitivity"][k]) for k in ("first_gt_1e-12", "first_gt_1e-9", "first_gt_1e-6", "first_gt_1e-3")],
    "seconds": round(time.time() - t0, 1)}

# ---------------------------------------------------------------- Parts B and C
stored = {(r["mode"], int(r["key_index"]), int(r["round"])): r["sha256"] for r in csv.DictReader(open(ROOT / "04_DATASETS" / "dynamic_family_direct_8192.csv"))}
X = np.arange(256, dtype=np.uint8); rows = []; t0 = time.time()
for mode_flag, mode in ((True, "HCHM+SHAKE"), (False, "SHAKE-only")):
    for ki in range(256):
        key = m.key_for_index(40000 + ki)
        for r in range(16):
            base = key + r.to_bytes(4, "big") + b"|HAS1-FAMILY-V2|"
            if mode_flag: base += m.serialize_hchm(key, r)
            S = m.candidate(base, 0)[0]
            h = hashlib.sha256(S.tobytes()).hexdigest(); M = m.sac(S); bmax, bmean = m.bic_corr(S); cyc = m.cycle_structure(S)
            rows.append({"mode": mode, "key_index": ki, "round": r, "sha256": h, "same_as_stored": int(h == stored[(mode, ki, r)]),
                         "sac_mad": float(np.abs(M - .5).mean()), "sac_maxdev": float(np.abs(M - .5).max()), "bic_max_abs_corr": float(bmax), "bic_mean_abs_corr": float(bmean),
                         "fixed": int(np.sum(S == X)), "opposite": int(np.sum(S == (X ^ 255))), "cycles": len(cyc), "longest_cycle": int(cyc[0])})
with open(out / "DATA" / ("family_metrics_%s.csv" % tag), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
arm = lambda mode: [r for r in rows if r["mode"] == mode]
res["portability"] = {"stored_corpus_environment": "original package run (Linux); see 03_REPRODUCIBILITY_AND_QA",
                      "tables_identical_to_stored": {mode: [sum(r["same_as_stored"] for r in arm(mode)), len(arm(mode))] for mode in ("HCHM+SHAKE", "SHAKE-only")},
                      "distinct_tables_in_rerun": len({r["sha256"] for r in rows}), "seconds": round(time.time() - t0, 1)}
tests = []
for name, kind in (("sac_mad", "ks"), ("sac_maxdev", "ks"), ("bic_max_abs_corr", "ks"), ("bic_mean_abs_corr", "ks"), ("fixed", "mwu"), ("opposite", "mwu"), ("cycles", "mwu"), ("longest_cycle", "mwu")):
    a = np.array([r[name] for r in arm("HCHM+SHAKE")], float); b = np.array([r[name] for r in arm("SHAKE-only")], float)
    if kind == "ks": st = stats.ks_2samp(a, b); stat, p, tname = float(st.statistic), float(st.pvalue), "Kolmogorov-Smirnov"
    else: st = stats.mannwhitneyu(a, b, alternative="two-sided"); stat, p, tname = float(st.statistic), float(st.pvalue), "Mann-Whitney U"
    pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    tests.append({"statistic": name, "test": tname, "n_per_arm": len(a), "mean_hchm_shake": float(a.mean()), "mean_shake_only": float(b.mean()),
                  "standardized_mean_difference": float((a.mean() - b.mean()) / pooled), "test_statistic": stat, "p": p})
order = sorted(range(len(tests)), key=lambda i: tests[i]["p"]); run = 0.0
for rank, i in enumerate(order):
    run = max(run, min(1.0, tests[i]["p"] * (len(tests) - rank))); tests[i]["p_holm"] = run
res["extended_ablation"] = {"prespecified_statistics": 8, "tests": tests, "smallest_p": min(t["p"] for t in tests), "smallest_p_holm": min(t["p_holm"] for t in tests),
                            "any_significant_at_0.05_after_holm": any(t["p_holm"] < 0.05 for t in tests)}
# Added after the first run, when two statistics turned out to take very few distinct values in this family
# (SAC maximum deviation: two values; BIC maximum correlation: fifteen). The Kolmogorov-Smirnov test is weak for such data,
# so a test suited to each is reported as a supplement. These two tests were not part of the plan fixed in advance.
hm, sm = np.array([r["sac_maxdev"] for r in arm("HCHM+SHAKE")]), np.array([r["sac_maxdev"] for r in arm("SHAKE-only")])
top = max(hm.max(), sm.max()); tab = [[int((hm == top).sum()), int((hm != top).sum())], [int((sm == top).sum()), int((sm != top).sum())]]
res["extended_ablation"]["supplementary_post_hoc"] = {
    "distinct_values": {k: len({round(r[k], 9) for r in rows}) for k in ("sac_maxdev", "bic_max_abs_corr")},
    "sac_maxdev_counts_at_largest_value_vs_other": tab, "sac_maxdev_fisher_exact_p": float(stats.fisher_exact(tab)[1]),
    "bic_max_abs_corr_mann_whitney_p": float(stats.mannwhitneyu([r["bic_max_abs_corr"] for r in arm("HCHM+SHAKE")], [r["bic_max_abs_corr"] for r in arm("SHAKE-only")], alternative="two-sided").pvalue)}
# ---------------------------------------------------------------- Part D: lattice serialization, a tested mitigation
# The published generator serializes a Float64 orbit, so one differing last bit in any transcendental value changes the table.
# Here the same map is iterated on the lattice with LATTICE_Q bits per coordinate (the emulator of the finite-state section),
# where last-bit differences are absorbed by rounding except when a value falls within about 1e-10 lattice units of a boundary.
import math
LATTICE_Q = 16
BLK = (((2, 1), (1, 1)), ((3, 2), (1, 1)), ((5, 2), (2, 1)), ((7, 3), (2, 1)))
def lattice_serialize(key, rnd, burn=64, samples=16, pert=0):
    sysm = m.HCHM8(key, rnd); p = sysm.p; eps = [float(e) for e in p.epsilon]; i0 = float(sysm.i0_eta); tp = 2 * math.pi
    phi1, phi2, rho3, k4, eta6, k7, rho8 = float(p.phi1), float(p.phi2), float(p.rho3), float(p.kappa4), float(p.eta6), float(p.kappa7), float(p.rho8)
    def f(i, z):
        if i == 0: return math.sin(tp * z[2] + phi1)
        if i == 1:
            c = math.cos(tp * z[5] + phi2); return 4 * c * c * c - 3 * c
        if i == 2: return math.sin(tp * (z[7] - z[1])) / (1 + rho3 * math.cos(tp * z[1]) ** 2)
        if i == 3: return math.tanh(k4 * math.sin(tp * (z[0] - z[6])))
        if i == 4:
            q = .5 * (1 + math.sin(tp * (z[6] + z[7]))); return 4 * q * (1 - q) - .5
        if i == 5: return (math.exp(eta6 * math.cos(tp * (z[1] - z[0]))) - i0) / (math.exp(eta6) - math.exp(-eta6))
        if i == 6: return (2 / math.pi) * math.atan(k7 * math.sin(tp * (z[3] + z[5])))
        return math.log(1 + rho8 * math.sin(math.pi * (z[4] - z[2])) ** 2) / math.log(1 + rho8) - .5
    N = 1 << LATTICE_Q; s = [min(int(v * N), N - 1) for v in sysm.x0]; outb = bytearray()
    for n in range(burn + samples):
        h = 0
        if pert:
            h = 1469598103934665603
            for v in s: h = ((h ^ v) * 1099511628211) & 0xFFFFFFFFFFFFFFFF
        z = [v / N for v in s]
        for i in range(8):
            fv = f(i, z)
            if pert: fv *= 1.0 + pert * (((h >> (4 * i)) & 15) % 9 - 4) * 2.220446049250313e-16
            z[i] = (z[i] + eps[i] * fv) % 1.0
        y = []
        for k, ((a, b), (c, d)) in enumerate(BLK): y += [(a * z[2 * k] + b * z[2 * k + 1]) % 1.0, (c * z[2 * k] + d * z[2 * k + 1]) % 1.0]
        for i in range(7, -1, -1):
            fv = f(i, y)
            if pert: fv *= 1.0 + pert * (((h >> (4 * (i + 8))) & 15) % 9 - 4) * 2.220446049250313e-16
            y[i] = (y[i] - eps[i] * fv) % 1.0
        s = [int(math.floor(v * N + .5)) % N for v in y]
        if n >= burn:
            for v in s: outb += v.to_bytes(2, "big")
    return bytes(outb)

t0 = time.time(); lat_rows = []; same_pert = 0; changed_big = 0
for ki in range(256):
    key = m.key_for_index(40000 + ki)
    for r in range(16):
        base = key + r.to_bytes(4, "big") + b"|HAS1-FAMILY-V3-LATTICE16|"
        ser = lattice_serialize(key, r); S = m.candidate(base + ser, 0)[0]
        same_pert += int(lattice_serialize(key, r, pert=1) == ser)
        if r == 0: changed_big += int(lattice_serialize(key, r, pert=1e12) != ser)
        M = m.sac(S); bmax, bmean = m.bic_corr(S); cyc = m.cycle_structure(S)
        lat_rows.append({"mode": "HCHM-lattice+SHAKE", "key_index": ki, "round": r, "sha256": hashlib.sha256(S.tobytes()).hexdigest(),
                         "sac_mad": float(np.abs(M - .5).mean()), "sac_maxdev": float(np.abs(M - .5).max()), "bic_max_abs_corr": float(bmax), "bic_mean_abs_corr": float(bmean),
                         "fixed": int(np.sum(S == X)), "opposite": int(np.sum(S == (X ^ 255))), "cycles": len(cyc), "longest_cycle": int(cyc[0])})
with open(out / "DATA" / ("family_lattice_metrics_%s.csv" % tag), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(lat_rows[0])); w.writeheader(); w.writerows(lat_rows)
lat_tests = []
for name, kind in (("sac_mad", "ks"), ("sac_maxdev", "ks"), ("bic_max_abs_corr", "ks"), ("bic_mean_abs_corr", "ks"), ("fixed", "mwu"), ("opposite", "mwu"), ("cycles", "mwu"), ("longest_cycle", "mwu")):
    a = np.array([r[name] for r in lat_rows], float); b = np.array([r[name] for r in arm("SHAKE-only")], float)
    pv = float(stats.ks_2samp(a, b).pvalue) if kind == "ks" else float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    lat_tests.append({"statistic": name, "mean_lattice_arm": float(a.mean()), "mean_shake_only": float(b.mean()), "p": pv})
order = sorted(range(8), key=lambda i: lat_tests[i]["p"]); run = 0.0
for rank, i in enumerate(order):
    run = max(run, min(1.0, lat_tests[i]["p"] * (8 - rank))); lat_tests[i]["p_holm"] = run
res["lattice_serialization"] = {
    "bits_per_coordinate": LATTICE_Q, "serialized_bits_per_round": 16 * 8 * LATTICE_Q, "tables": len(lat_rows),
    "distinct_tables": len({r["sha256"] for r in lat_rows}),
    "serializations_unchanged_under_last_bit_perturbation": [same_pert, len(lat_rows)],
    "sanity_serializations_changed_under_perturbation_scaled_1e12": [changed_big, 256],
    "ablation_against_shake_only": lat_tests, "any_significant_at_0.05_after_holm": any(t["p_holm"] < 0.05 for t in lat_tests),
    "seconds": round(time.time() - t0, 1)}
for other in sorted((out / "DATA").glob("family_lattice_metrics_*.csv")):
    if other.name != "family_lattice_metrics_%s.csv" % tag:
        ref = {(r["key_index"], r["round"]): r["sha256"] for r in csv.DictReader(open(other))}
        res["lattice_serialization"]["tables_identical_to_" + other.stem.replace("family_lattice_metrics_", "")] = [sum(ref[(str(r["key_index"]), str(r["round"]))] == r["sha256"] for r in lat_rows), len(lat_rows)]
for other in sorted((out / "DATA").glob("family_metrics_*.csv")):
    if other.name != "family_metrics_%s.csv" % tag:
        ref = {(r["mode"], r["key_index"], r["round"]): r["sha256"] for r in csv.DictReader(open(other))}
        res["portability"]["tables_identical_to_" + other.stem.replace("family_metrics_", "")] = {mode: [sum(ref[(mode, str(r["key_index"]), str(r["round"]))] == r["sha256"] for r in arm(mode)), len(arm(mode))] for mode in ("HCHM+SHAKE", "SHAKE-only")}

json.dump(res, open(out / ("RESULTS_%s.json" % tag), "w"), indent=1)
print("lattice serialization:", json.dumps({k: v for k, v in res["lattice_serialization"].items() if k != "ablation_against_shake_only"}))
print("lattice arm against SHAKE-only, smallest p = %.4f, smallest Holm p = %.4f" % (min(t["p"] for t in lat_tests), min(t["p_holm"] for t in lat_tests)))
shutil.rmtree(tmp, ignore_errors=True)
print(json.dumps({"environment": res["environment"], "map_numerics": {k: res["map_numerics"][k] for k in ("max_abs_numeric_minus_exact_stored", "max_abs_numeric_minus_exact_rerun", "max_abs_rerun_minus_stored_mean", "numeric_le_identical")},
                  "portability": res["portability"]["tables_identical_to_stored"]}, indent=1))
for t in tests: print("%-18s %-19s mean %.6f vs %.6f  SMD %+.4f  p=%.4f  Holm p=%.4f" % (t["statistic"], t["test"], t["mean_hchm_shake"], t["mean_shake_only"], t["standardized_mean_difference"], t["p"], t["p_holm"]))
