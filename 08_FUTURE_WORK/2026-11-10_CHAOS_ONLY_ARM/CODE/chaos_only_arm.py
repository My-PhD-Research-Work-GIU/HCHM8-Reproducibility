"""Chaos-only arms of the S-box ablation (plan: ../PLAN.md, fixed before the runs).

Four arms over the same 256 keys x 16 rounds (key_for_index(40000 + k), rounds 0 to 15) as the extended ablation of the
paper; every table is candidate 0 of the derivation, so the arms differ only in the 258-byte stream that drives the affine
layers (make_rows, constants a and b of run_validated_experiments.candidate):

  S   SHAKE-only      stream = SHAKE256(key | round | label | counter | label)                     (paper, unchanged)
  H   HCHM+SHAKE      stream = SHAKE256(key | round | label | Float64 orbit | counter | label)     (paper, unchanged)
  C   chaos-only      stream = the first 258 bytes of the Float64 orbit of HCHM8(key, round) after 64 burn-in steps,
                      32 bits per coordinate as in serialize_hchm, with no hash
  CL  chaos-only, lattice  stream = the first 258 bytes of the 16-bit lattice orbit (the portable serialization of the
                      paper), 64 burn-in steps, with no hash

Statistics, tests and effect sizes are those of the extended ablation: SAC mean absolute deviation and maximum deviation, BIC
maximum and mean absolute correlation (two-sample Kolmogorov-Smirnov), fixed points, opposite fixed points, number of cycles,
longest cycle (Mann-Whitney U); Holm correction over the eight tests of each comparison; Hedges' g with a 95% interval and
Cliff's delta. Primary comparisons: C against S and CL against S. Secondary: C against H. Descriptive: byte frequencies of the
C and CL streams (chi-square against uniform), tables identical across arms.

Usage: python chaos_only_arm.py <out_dir> [--keys a:b] [--smoke]
"""
import argparse, csv, hashlib, importlib, json, math, platform, shutil, sys, tempfile, time
from pathlib import Path
import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[2]
ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--keys", default="0:256"); ap.add_argument("--smoke", action="store_true")
a = ap.parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True); k0, k1 = (int(v) for v in a.keys.split(":"))
tmp = Path(tempfile.mkdtemp(prefix="r01_")); (tmp / "code").mkdir()
shutil.copy(ROOT / "05_CODE" / "run_validated_experiments.py", tmp / "code" / "run_validated_experiments.py")
shutil.copy(ROOT / "05_CODE" / "platform_rerun_and_extended_ablation.py", tmp / "code" / "_platform_source.py")
sys.path.insert(0, str(tmp / "code"))
import matplotlib; matplotlib.use("Agg")
m = importlib.import_module("run_validated_experiments")

# The 16-bit lattice serialization of the paper, taken verbatim from platform_rerun_and_extended_ablation.py (Part D) so that
# CL uses exactly the portable orbit of the paper; only the number of samples differs (17 instead of 16, to give 258 bytes).
src = (tmp / "code" / "_platform_source.py").read_text(encoding="utf-8")
i0 = src.index("import math\nLATTICE_Q"); i1 = src.index("t0 = time.time(); lat_rows")
ns = {"m": m}; exec(compile(src[i0:i1], "lattice_serialize", "exec"), ns); lattice_serialize = ns["lattice_serialize"]

NBYTES = 258
def stream_S(key, r): return m.shake(key + r.to_bytes(4, "big") + b"|HAS1-FAMILY-V2|" + (0).to_bytes(2, "big") + b"|AFFINE-PROFILE-V2|", NBYTES)
def stream_H(key, r): return m.shake(key + r.to_bytes(4, "big") + b"|HAS1-FAMILY-V2|" + m.serialize_hchm(key, r) + (0).to_bytes(2, "big") + b"|AFFINE-PROFILE-V2|", NBYTES)
def stream_C(key, r): return m.serialize_hchm(key, r, burn=64, samples=9)[:NBYTES]
def stream_CL(key, r): return lattice_serialize(key, r, burn=64, samples=17)[:NBYTES]

def table(st):
    rin, idx = m.make_rows(st, 0); rout, idx = m.make_rows(st, idx); return m.build_sbox(rin, rout, st[idx % len(st)], st[(idx + 1) % len(st)])

# consistency with the paper's corpus: stream_S and stream_H must give candidate(base, 0) exactly
X = np.arange(256, dtype=np.uint8)
def metrics(S):
    M = m.sac(S); bmax, bmean = m.bic_corr(S); cyc = m.cycle_structure(S)
    return {"sac_mad": float(np.abs(M - .5).mean()), "sac_maxdev": float(np.abs(M - .5).max()), "bic_max_abs_corr": float(bmax), "bic_mean_abs_corr": float(bmean),
            "fixed": int(np.sum(S == X)), "opposite": int(np.sum(S == (X ^ 255))), "cycles": len(cyc), "longest_cycle": int(cyc[0])}

res = {"environment": {"platform": platform.platform(), "python": platform.python_version(), "numpy": np.__version__, "scipy": __import__("scipy").__version__},
       "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "keys": [k0, k1], "smoke_test": a.smoke,
       "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
rows = []; bytes_C = np.zeros(256, np.int64); bytes_CL = np.zeros(256, np.int64); check = 0; t0 = time.time()
for ki in range(k0, k1):
    key = m.key_for_index(40000 + ki)
    for r in range(16):
        streams = {"S": stream_S(key, r), "H": stream_H(key, r), "C": stream_C(key, r), "CL": stream_CL(key, r)}
        base = key + r.to_bytes(4, "big") + b"|HAS1-FAMILY-V2|"
        assert np.array_equal(table(streams["S"]), m.candidate(base, 0)[0]); assert np.array_equal(table(streams["H"]), m.candidate(base + m.serialize_hchm(key, r), 0)[0]); check += 1
        bytes_C += np.bincount(np.frombuffer(streams["C"], np.uint8), minlength=256); bytes_CL += np.bincount(np.frombuffer(streams["CL"], np.uint8), minlength=256)
        for arm, st in streams.items():
            S = table(st); rows.append({"arm": arm, "key_index": ki, "round": r, "sha256": hashlib.sha256(S.tobytes()).hexdigest(), **metrics(S)})
res["paper_arms_reproduced_exactly"] = check; res["seconds_generation"] = round(time.time() - t0, 1)
with open(out / ("chaos_only_metrics%s.csv" % ("_smoke" if a.smoke else "")), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

STATS = (("sac_mad", "ks"), ("sac_maxdev", "ks"), ("bic_max_abs_corr", "ks"), ("bic_mean_abs_corr", "ks"), ("fixed", "mwu"), ("opposite", "mwu"), ("cycles", "mwu"), ("longest_cycle", "mwu"))
def arm(name): return {(r["key_index"], r["round"]): r for r in rows if r["arm"] == name}
def compare(x_name, y_name):
    X_, Y_ = arm(x_name), arm(y_name); keys = sorted(X_); tests = []
    for s, kind in STATS:
        x = np.array([X_[k][s] for k in keys], float); y = np.array([Y_[k][s] for k in keys], float); n1, n2 = len(x), len(y)
        p = float(stats.ks_2samp(x, y).pvalue) if kind == "ks" else float(stats.mannwhitneyu(x, y, alternative="two-sided").pvalue)
        sp = math.sqrt(((n1 - 1) * x.var(ddof=1) + (n2 - 1) * y.var(ddof=1)) / (n1 + n2 - 2)); g = (1 - 3 / (4 * (n1 + n2) - 9)) * ((x.mean() - y.mean()) / sp if sp > 0 else 0.0)
        se = math.sqrt((n1 + n2) / (n1 * n2) + g * g / (2 * (n1 + n2)))
        ranks = stats.rankdata(np.concatenate([x, y])); U1 = ranks[:n1].sum() - n1 * (n1 + 1) / 2; delta = 2 * U1 / (n1 * n2) - 1
        tests.append({"statistic": s, "test": "Kolmogorov-Smirnov" if kind == "ks" else "Mann-Whitney U", "mean_x": float(x.mean()), "mean_y": float(y.mean()),
                      "p": p, "hedges_g": g, "g_ci95": [g - 1.96 * se, g + 1.96 * se], "cliffs_delta": float(delta)})
    order = sorted(range(len(tests)), key=lambda i: tests[i]["p"]); run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, tests[i]["p"] * (len(tests) - rank))); tests[i]["p_holm"] = run
    return {"x": x_name, "y": y_name, "n_per_arm": len(keys), "tests": tests, "smallest_p": min(t["p"] for t in tests), "smallest_p_holm": min(t["p_holm"] for t in tests),
            "any_significant_after_holm": any(t["p_holm"] < 0.05 for t in tests), "max_abs_g": max(abs(t["hedges_g"]) for t in tests),
            "g_interval_bound": max(max(abs(t["g_ci95"][0]), abs(t["g_ci95"][1])) for t in tests)}
res["primary"] = [compare("C", "S"), compare("CL", "S")]; res["secondary"] = [compare("C", "H")]; res["paper_comparison_recomputed"] = compare("H", "S")
def chi(c): e = c.sum() / 256; st = float(((c - e) ** 2 / e).sum()); return {"bytes": int(c.sum()), "chi2": st, "df": 255, "p": float(stats.chi2.sf(st, 255))}
res["stream_byte_uniformity"] = {"C": chi(bytes_C), "CL": chi(bytes_CL)}
res["distinct_tables"] = {name: len({r["sha256"] for r in rows if r["arm"] == name}) for name in ("S", "H", "C", "CL")}
res["tables_shared_between_arms"] = {f"{p}-{q}": sum(arm(p)[k]["sha256"] == arm(q)[k]["sha256"] for k in arm(p)) for p, q in (("C", "S"), ("CL", "S"), ("C", "H"), ("C", "CL"))}
res["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
(out / ("RESULTS_chaos_only%s.json" % ("_smoke" if a.smoke else ""))).write_text(json.dumps(res, indent=1))
shutil.rmtree(tmp, ignore_errors=True)
for c in res["primary"] + res["secondary"] + [res["paper_comparison_recomputed"]]:
    print("%s vs %s: smallest p %.4f, smallest Holm p %.4f, max |g| %.4f, g bound %.4f" % (c["x"], c["y"], c["smallest_p"], c["smallest_p_holm"], c["max_abs_g"], c["g_interval_bound"]))
print(json.dumps(res["stream_byte_uniformity"]), json.dumps(res["distinct_tables"]))
