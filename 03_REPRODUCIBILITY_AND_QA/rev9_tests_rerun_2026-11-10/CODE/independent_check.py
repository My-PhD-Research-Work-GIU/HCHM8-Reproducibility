"""Independent recomputation of the two tests of paper revision 9 (plan: ../PLAN.md, fixed before the runs).

This program shares no code with the scripts that produced the results, except that it reads their inputs as data: the
parameter vector of each profile (profile_params) and the starting states (initial) of the finite-state study, and, in Part C,
the byte streams of the four arms, which it rebuilds with the paper's derivation functions so that the tables can be checked
against the stored file. Everything that is measured is computed here from the equations of the paper:

Part A  Every lattice state at Q = 2 and 3, ten profiles. F_Q (rounded emulator) and G_Q (grid-exact conjugate) are evaluated
        with NumPy arrays from the shear equations of Section 3 (the Bessel value I0(eta6) from scipy.special, the blocks of A
        written out); their functional graphs are analysed by a path-following algorithm written for this check (the
        original peels the graph by in-degree). Compared with the stored RESULTS_grid_conjugacy.json, and, state by state,
        with the arrays of the original program.
Part B  The orbit periods of G_Q at Q = 4 to 7 from the starting states of the uncapped measurement, by iteration of the
        NumPy map; compared with the stored periods.
Part C  The 16,384 tables of the chaos-only test: SHA-256 and the eight statistics recomputed by functions written here and
        compared with DATA/chaos_only_metrics.csv; then every number of Table 10 recomputed from that file (Hedges' g and its
        interval, Cliff's delta, Holm), with the correlation within matched pairs and intervals adjusted for eight
        statistics (Bonferroni), which the paper did not report before.

Usage: python independent_check.py <out_dir> [--q 2,3] [--parts A,B,C]
"""
import argparse, csv, hashlib, importlib, importlib.util, json, math, os, platform, shutil, sys, tempfile, time
from pathlib import Path
import numpy as np
from scipy import special, stats

os.environ.setdefault("NUMBA_CACHE_DIR", tempfile.mkdtemp(prefix="r01_ind_"))
import numba
from numba import njit

HERE = Path(__file__).resolve().parent
B = HERE.parents[2]                                    # 05_PAPERS_AND_DATASETS/05_PAPERS_AND_DATASETS
FW = B / "08_FUTURE_WORK"
EMU = FW / "2026-10-11_FINITE_STATE_CYCLE_SCALING" / "CODE" / "cycle_scaling_profiles.py"
GRID = FW / "2026-11-10_GRID_CONJUGACY"
CHAOS = FW / "2026-11-10_CHAOS_ONLY_ARM"


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path); mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod


# ======================================================================================= the map, from the paper's equations
BLOCKS = ((2, 1, 1, 1), (3, 2, 1, 1), (5, 2, 2, 1), (7, 3, 2, 1))
TWO_PI = 2.0 * math.pi


class Profile:
    """Parameters read from the profile vector [eps1..eps8, phi1, phi2, rho3, kappa4, eta6, kappa7, rho8, ...]."""
    def __init__(self, p):
        self.eps = [float(v) for v in p[:8]]
        self.phi1, self.phi2, self.rho3, self.kappa4, self.eta6, self.kappa7, self.rho8 = (float(v) for v in p[8:15])
        self.i0 = float(special.i0(self.eta6))
        self.den6 = math.exp(self.eta6) - math.exp(-self.eta6)
        self.log8 = math.log(1.0 + self.rho8)

    def g(self, i, z):
        """Shear function g_{i+1} of Section 3 on an (n, 8) array of torus points; it does not read column i."""
        z1, z2, z3, z4, z5, z6, z7, z8 = (z[:, k] for k in range(8))
        if i == 0: return np.sin(TWO_PI * z3 + self.phi1)
        if i == 1:
            c = np.cos(TWO_PI * z6 + self.phi2); return 4.0 * c ** 3 - 3.0 * c
        if i == 2: return np.sin(TWO_PI * (z8 - z2)) / (1.0 + self.rho3 * np.cos(TWO_PI * z2) ** 2)
        if i == 3: return np.tanh(self.kappa4 * np.sin(TWO_PI * (z1 - z7)))
        if i == 4:
            q = 0.5 * (1.0 + np.sin(TWO_PI * (z7 + z8))); return 4.0 * q * (1.0 - q) - 0.5
        if i == 5: return (np.exp(self.eta6 * np.cos(TWO_PI * (z2 - z1))) - self.i0) / self.den6
        if i == 6: return (2.0 / math.pi) * np.arctan(self.kappa7 * np.sin(TWO_PI * (z4 + z6)))
        return np.log(1.0 + self.rho8 * np.sin(math.pi * (z5 - z3)) ** 2) / self.log8 - 0.5


def frac(t):
    return t - np.floor(t)


def core_float(z):
    y = np.empty_like(z)
    for k, (a, b, c, d) in enumerate(BLOCKS):
        y[:, 2 * k] = frac(a * z[:, 2 * k] + b * z[:, 2 * k + 1]); y[:, 2 * k + 1] = frac(c * z[:, 2 * k] + d * z[:, 2 * k + 1])
    return y


def core_int(s, Q):
    y = np.empty_like(s)
    for k, (a, b, c, d) in enumerate(BLOCKS):
        y[:, 2 * k] = (a * s[:, 2 * k] + b * s[:, 2 * k + 1]) % Q; y[:, 2 * k + 1] = (c * s[:, 2 * k] + d * s[:, 2 * k + 1]) % Q
    return y


def F_map(s, pr, qbits):
    """Rounded emulator: lattice point -> torus -> T = H^-1 A H in Float64 -> nearest lattice point."""
    Q = 1 << qbits; z = s / Q
    for i in range(8): z[:, i] = frac(z[:, i] + pr.eps[i] * pr.g(i, z))
    y = core_float(z)
    for i in range(7, -1, -1): y[:, i] = frac(y[:, i] - pr.eps[i] * pr.g(i, y))
    return np.floor(y * Q + 0.5).astype(np.int64) % Q


def HQ(s, pr, qbits, inverse=False):
    Q = 1 << qbits; z = s.copy()
    order = range(7, -1, -1) if inverse else range(8)
    for i in order:
        inc = np.floor(Q * pr.eps[i] * pr.g(i, z / Q) + 0.5).astype(np.int64)
        z[:, i] = (z[:, i] - inc) % Q if inverse else (z[:, i] + inc) % Q
    return z


def G_map(s, pr, qbits):
    return HQ(core_int(HQ(s, pr, qbits), 1 << qbits), pr, qbits, inverse=True)


def all_states(qbits, lo, hi):
    idx = np.arange(lo, hi, dtype=np.int64); m = (1 << qbits) - 1
    return np.stack([(idx >> (qbits * k)) & m for k in range(8)], axis=1)


def encode(s, qbits):
    v = np.zeros(len(s), np.int64)
    for k in range(8): v |= s[:, k] << (qbits * k)
    return v


def successor_array(fn, qbits, chunk=1 << 20):
    N = 1 << (8 * qbits); nxt = np.empty(N, np.int64)
    for lo in range(0, N, chunk):
        hi = min(N, lo + chunk); nxt[lo:hi] = encode(fn(all_states(qbits, lo, hi)), qbits)
    return nxt


# ======================================================================================= functional graphs, path following
@njit(cache=True)
def walk_graph(nxt):
    """tail and cycle length of every point, by following each unvisited path until it meets a known point."""
    N = nxt.shape[0]; state = np.zeros(N, np.int8)          # 0 new, 1 on the current path, 2 done
    tail = np.zeros(N, np.int64); cyc = np.zeros(N, np.int64); path = np.empty(N, np.int64)
    for s in range(N):
        if state[s] != 0: continue
        n = 0; v = s
        while state[v] == 0:
            state[v] = 1; path[n] = v; n += 1; v = nxt[v]
        if state[v] == 1:                                   # the path closed on itself: a new cycle starts at v
            k = n - 1
            while path[k] != v: k -= 1
            L = n - k
            for j in range(k, n):
                cyc[path[j]] = L; tail[path[j]] = 0; state[path[j]] = 2
            n = k
        for j in range(n - 1, -1, -1):                      # the rest of the path hangs on known points
            u = path[j]; w = nxt[u]; tail[u] = tail[w] + 1; cyc[u] = cyc[w]; state[u] = 2
    return tail, cyc


def graph_summary(nxt):
    N = len(nxt); tail, cyc = walk_graph(nxt)
    on_cycle = tail == 0; lens, counts = np.unique(cyc[on_cycle], return_counts=True)
    ctype = {int(L): int(c // L) for L, c in zip(lens, counts)}
    return {"image_fraction": float(np.count_nonzero(np.bincount(nxt, minlength=N)) / N), "cyclic_states": int(on_cycle.sum()),
            "cycles": int(sum(ctype.values())), "mean_rho": float((tail + cyc).mean()), "max_tail": int(tail.max()),
            "longest_cycle": int(lens.max())}, ctype


# ======================================================================================= Part A and Part B
def part_A(qs, emu, orig, stored):
    out = {"rows": [], "all_equal_to_stored": True}
    for q in qs:
        N = 1 << (8 * q); rm = math.sqrt(math.pi * N / 2)
        nA = successor_array(lambda s: core_int(s, 1 << q), q); sA, ctA = graph_summary(nA)
        for prof in range(10):
            P = emu.profile_params(prof); P[16] = 0.0; pr = Profile(P); t = time.time()
            nF = successor_array(lambda s: F_map(s, pr, q), q); nG = successor_array(lambda s: G_map(s, pr, q), q)
            sF, _ = graph_summary(nF); sG, ctG = graph_summary(nG)
            oF = orig.next_array(q, P, 0); oG = orig.next_array(q, P, 1)
            st = next(r for r in stored["exhaustive"] if r.get("qbits") == q and r.get("profile") == prof)
            sF["mean_rho_over_sqrt_piN_2"] = sF["mean_rho"] / rm; sF["cyclic_over_sqrt_piN_2"] = sF["cyclic_states"] / rm
            fields = ("image_fraction", "cyclic_states", "cycles", "mean_rho", "max_tail", "longest_cycle")
            row = {"qbits": q, "profile": prof,
                   "G_is_permutation": sG["image_fraction"] == 1.0, "G_cycle_type_equals_A": ctG == ctA,
                   "G_cycle_type_equals_stored": {str(k): v for k, v in ctG.items()} == {str(k): v for k, v in st["G"]["cycle_type"].items()},
                   "F_stats": sF, "F_stats_equal_stored": {f: (sF[f] == st["F"][f]) if isinstance(sF[f], int) else abs(sF[f] - st["F"][f]) < 1e-9 for f in fields},
                   "F_states_differing_from_original": int(np.count_nonzero(nF != oF)), "G_states_differing_from_original": int(np.count_nonzero(nG != oG)),
                   "F_equals_G_fraction": float(np.mean(nF == nG)), "F_equals_G_fraction_stored": st["F_equals_G_fraction"], "seconds": round(time.time() - t, 1)}
            row["all_equal"] = (row["G_is_permutation"] and row["G_cycle_type_equals_A"] and row["G_cycle_type_equals_stored"]
                                and all(row["F_stats_equal_stored"].values()) and abs(row["F_equals_G_fraction"] - row["F_equals_G_fraction_stored"]) < 1e-12)
            out["all_equal_to_stored"] &= row["all_equal"]; out["rows"].append(row)
            print("A Q%d profile %d: G perm %s, type=A %s, F stats equal %s, F differs on %d, G differs on %d states (%.0f s)" % (
                q, prof, row["G_is_permutation"], row["G_cycle_type_equals_A"], all(row["F_stats_equal_stored"].values()),
                row["F_states_differing_from_original"], row["G_states_differing_from_original"], row["seconds"]), flush=True)
        out["rows"].append({"qbits": q, "A_cycle_type": {str(k): v for k, v in ctA.items()}, "A_summary": sA})
    return out


def part_B(emu, stored):
    starts = {4: 32, 5: 32, 6: 13, 7: 2}; res = {"orbits": 0, "equal_to_stored": 0, "periods": {}}
    stored_p = {(r["qbits"], r["profile"], r["run"]): (r["period_G"], r["period_A_of_HQx0"]) for r in stored["orbits"]}
    for q, n in starts.items():
        core = 3 * 2 ** (q - 1); seen = set()
        for prof in range(10):
            P = emu.profile_params(prof); P[16] = 0.0; pr = Profile(P)
            x0 = np.stack([emu.initial(q, run) for run in range(n)]); x = x0.copy(); period = np.zeros(n, np.int64)
            y = HQ(x0, pr, q); ya = y.copy(); periodA = np.zeros(n, np.int64)
            for step in range(1, core + 2):
                x = G_map(x, pr, q); hit = (period == 0) & np.all(x == x0, axis=1); period[hit] = step
                ya = core_int(ya, 1 << q); hitA = (periodA == 0) & np.all(ya == y, axis=1); periodA[hitA] = step
            for run in range(n):
                res["orbits"] += 1; seen.add(int(period[run]))
                res["equal_to_stored"] += int(stored_p[(q, prof, run)] == (int(period[run]), int(periodA[run])) and period[run] == periodA[run] and core % period[run] == 0)
        res["periods"][str(q)] = sorted(seen)
        print("B Q%d: periods %s" % (q, sorted(seen)), flush=True)
    res["all_equal"] = res["orbits"] == res["equal_to_stored"]
    return res


# ======================================================================================= Part C
X = np.arange(256)
BITS = (X[:, None] >> np.arange(8)) & 1                     # bit j of x


def metrics(S):
    S = np.asarray(S, np.int64); out = {}
    av = np.stack([(S ^ S[X ^ (1 << i)])[:, None] >> np.arange(8) & 1 for i in range(8)])   # av[i, x, j]
    sac = av.mean(axis=1)
    out["sac_mad"] = float(np.abs(sac - 0.5).mean()); out["sac_maxdev"] = float(np.abs(sac - 0.5).max())
    corr = []
    for i in range(8):
        a = av[i].astype(float); a = a - a.mean(axis=0); nrm = np.sqrt((a * a).sum(axis=0)); c = (a.T @ a) / np.outer(nrm, nrm)
        corr.extend(abs(c[p, q]) for p in range(8) for q in range(p + 1, 8))
    out["bic_max_abs_corr"] = float(max(corr)); out["bic_mean_abs_corr"] = float(np.mean(corr))
    out["fixed"] = int(np.sum(S == X)); out["opposite"] = int(np.sum(S == (X ^ 255)))
    seen = np.zeros(256, bool); lens = []
    for x in range(256):
        L = 0; y = x
        while not seen[y]: seen[y] = True; y = int(S[y]); L += 1
        if L: lens.append(L)
    out["cycles"] = len(lens); out["longest_cycle"] = max(lens)
    return out


STATS = (("sac_mad", "ks"), ("sac_maxdev", "ks"), ("bic_max_abs_corr", "ks"), ("bic_mean_abs_corr", "ks"),
         ("fixed", "mwu"), ("opposite", "mwu"), ("cycles", "mwu"), ("longest_cycle", "mwu"))
LABEL = {("H", "S"): "HCHM+SHAKE vs SHAKE-only", ("C", "S"): "Float64 orbit, no hash, vs SHAKE-only",
         ("CL", "S"): "lattice orbit, no hash, vs SHAKE-only", ("C", "H"): "Float64 orbit, no hash, vs HCHM+SHAKE"}


def comparison(rows, xa, ya):
    X_ = {(r["key_index"], r["round"]): r for r in rows if r["arm"] == xa}; Y_ = {(r["key_index"], r["round"]): r for r in rows if r["arm"] == ya}
    keys = sorted(X_); tests = []
    # 1.96, as in the script of the test; the exact quantile 1.959964 moves the interval ends by less than 1e-5
    z95 = 1.96; z_bonf = stats.norm.ppf(1 - 0.05 / 16)
    for s, kind in STATS:
        x = np.array([float(X_[k][s]) for k in keys]); y = np.array([float(Y_[k][s]) for k in keys]); n1 = len(x); n2 = len(y)
        p = stats.ks_2samp(x, y).pvalue if kind == "ks" else stats.mannwhitneyu(x, y, alternative="two-sided").pvalue
        pooled = math.sqrt((np.sum((x - x.mean()) ** 2) + np.sum((y - y.mean()) ** 2)) / (n1 + n2 - 2))
        g = (x.mean() - y.mean()) / pooled * (1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0))
        se = math.sqrt((n1 + n2) / (n1 * n2) + g * g / (2.0 * (n1 + n2)))
        ys = np.sort(y); gt = np.searchsorted(ys, x, side="left").sum(); lt = n2 - np.searchsorted(ys, x, side="right")
        delta = (gt - lt.sum()) / (n1 * n2)
        r_pair = float(np.corrcoef(x, y)[0, 1]) if x.std() > 0 and y.std() > 0 else 0.0
        tests.append({"statistic": s, "p": float(p), "g": float(g), "ci95": [g - z95 * se, g + z95 * se],
                      "ci_bonferroni8": [g - z_bonf * se, g + z_bonf * se], "delta": float(delta), "r_within_pairs": r_pair})
    m = len(tests); order = sorted(range(m), key=lambda i: tests[i]["p"]); running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * tests[i]["p"])); tests[i]["p_holm"] = running
    return {"label": LABEL[(xa, ya)], "x": xa, "y": ya, "n": len(keys), "tests": tests,
            "smallest_p": min(t["p"] for t in tests), "smallest_p_holm": min(t["p_holm"] for t in tests),
            "max_abs_g": max(abs(t["g"]) for t in tests), "g_bound": max(max(abs(v) for v in t["ci95"]) for t in tests),
            "g_bound_bonferroni8": max(max(abs(v) for v in t["ci_bonferroni8"]) for t in tests),
            "max_abs_delta": max(abs(t["delta"]) for t in tests), "max_abs_r_within_pairs": max(abs(t["r_within_pairs"]) for t in tests)}


def part_C():
    rows = []
    with open(CHAOS / "DATA" / "chaos_only_metrics.csv") as f:
        for r in csv.DictReader(f):
            for k in ("key_index", "round", "fixed", "opposite", "cycles", "longest_cycle"): r[k] = int(r[k])
            for k in ("sac_mad", "sac_maxdev", "bic_max_abs_corr", "bic_mean_abs_corr"): r[k] = float(r[k])
            rows.append(r)
    # rebuild the streams with the paper's derivation (input data), compute the tables and the statistics here
    tmp = Path(tempfile.mkdtemp(prefix="r01_indc_")); shutil.copy(B / "05_CODE" / "run_validated_experiments.py", tmp / "run_validated_experiments.py")
    src = (B / "05_CODE" / "platform_rerun_and_extended_ablation.py").read_text(encoding="utf-8")
    sys.path.insert(0, str(tmp)); import matplotlib; matplotlib.use("Agg"); m = importlib.import_module("run_validated_experiments")
    ns = {"m": m}; exec(compile(src[src.index("import math\nLATTICE_Q"):src.index("t0 = time.time(); lat_rows")], "lattice_serialize", "exec"), ns)
    lat = ns["lattice_serialize"]; tail = (0).to_bytes(2, "big") + b"|AFFINE-PROFILE-V2|"
    stream = {"S": lambda k, r: m.shake(k + r.to_bytes(4, "big") + b"|HAS1-FAMILY-V2|" + tail, 258),
              "H": lambda k, r: m.shake(k + r.to_bytes(4, "big") + b"|HAS1-FAMILY-V2|" + m.serialize_hchm(k, r) + tail, 258),
              "C": lambda k, r: m.serialize_hchm(k, r, burn=64, samples=9)[:258],
              "CL": lambda k, r: lat(k, r, burn=64, samples=17)[:258]}
    by = {(r["arm"], r["key_index"], r["round"]): r for r in rows}; sha_ok = 0; met_ok = 0; worst = 0.0; t0 = time.time()
    for ki in range(256):
        key = m.key_for_index(40000 + ki)
        for rnd in range(16):
            for arm in ("S", "H", "C", "CL"):
                st = stream[arm](key, rnd); rin, idx = m.make_rows(st, 0); rout, idx = m.make_rows(st, idx)
                S = m.build_sbox(rin, rout, st[idx % 258], st[(idx + 1) % 258]); ref = by[(arm, ki, rnd)]
                sha_ok += hashlib.sha256(S.tobytes()).hexdigest() == ref["sha256"]
                mt = metrics(S); d = max(abs(mt[k] - ref[k]) for k in mt); worst = max(worst, d); met_ok += d < 1e-12
    shutil.rmtree(tmp, ignore_errors=True)
    out = {"tables": len(rows), "sha256_equal": int(sha_ok), "metrics_equal": int(met_ok), "largest_metric_difference": worst,
           "seconds": round(time.time() - t0, 1)}
    print("C tables %d: SHA-256 equal %d, statistics equal %d (largest difference %.3g)" % (len(rows), sha_ok, met_ok, worst), flush=True)
    stored = json.loads((CHAOS / "RESULTS_chaos_only.json").read_text())
    st_cmp = {(c["x"], c["y"]): c for c in stored["primary"] + stored["secondary"] + [stored["paper_comparison_recomputed"]]}
    out["comparisons"] = []
    for xa, ya in (("H", "S"), ("C", "S"), ("CL", "S"), ("C", "H")):
        c = comparison(rows, xa, ya); s = st_cmp[(xa, ya)]
        c["equal_to_stored"] = {"smallest_p": abs(c["smallest_p"] - s["smallest_p"]) < 1e-12, "smallest_p_holm": abs(c["smallest_p_holm"] - s["smallest_p_holm"]) < 1e-12,
                                "max_abs_g": abs(c["max_abs_g"] - s["max_abs_g"]) < 1e-12, "g_bound": abs(c["g_bound"] - s["g_interval_bound"]) < 1e-12,
                                "max_abs_delta": abs(c["max_abs_delta"] - max(abs(t["cliffs_delta"]) for t in s["tests"])) < 1e-12}
        out["comparisons"].append(c)
        print("C %s: p %.4f, Holm %.4f, |g| %.4f, bound %.4f, Bonferroni bound %.4f, |delta| %.4f, |r| within pairs %.4f, equal to stored %s" % (
            c["label"], c["smallest_p"], c["smallest_p_holm"], c["max_abs_g"], c["g_bound"], c["g_bound_bonferroni8"], c["max_abs_delta"],
            c["max_abs_r_within_pairs"], all(c["equal_to_stored"].values())), flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--q", default="2,3"); ap.add_argument("--parts", default="A,B,C")
    a = ap.parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True); parts = a.parts.split(",")
    res = {"environment": {"platform": platform.platform(), "python": platform.python_version(), "numpy": np.__version__,
                           "scipy": __import__("scipy").__version__, "numba": numba.__version__},
           "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "parts": parts,
           "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    emu = load(EMU, "cycle_scaling_profiles"); stored = json.loads((GRID / "RESULTS_grid_conjugacy.json").read_text())
    if "A" in parts:
        orig = load(GRID / "CODE" / "grid_conjugacy.py", "grid_conjugacy_original")
        res["A"] = part_A([int(v) for v in a.q.split(",")], emu, orig, stored)
    if "B" in parts: res["B"] = part_B(emu, stored)
    if "C" in parts: res["C"] = part_C()
    res["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    (out / "RESULTS_independent_check.json").write_text(json.dumps(res, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


if __name__ == "__main__":
    main()
