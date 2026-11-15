"""Functional graphs, in-degree distributions and shear-family robustness of three discretizations of HCHM-8 (plan: ../PLAN.md).

Part 1 enumerates every state of the lattice (Q = 2, 3); Part 2 measures sampled orbits (Q = 4, 5). The emulator of the paper
(family 1, s = 1, nearest rounding, end quantization) uses f_scalar of
08_FUTURE_WORK/2026-10-11_FINITE_STATE_CYCLE_SCALING/CODE/cycle_scaling_profiles.py, imported unchanged, with the same order of
floating-point operations as its step(), so that its successors are those of the paper.

Usage: python fg_study.py <out_dir> [--profiles 0:10] [--q1 2,3] [--q2 4,5] [--starts 8] [--smoke] [--part 1|2|12]
"""
import argparse, hashlib, importlib.util, json, math, os, platform, sys, tempfile, time
from pathlib import Path
import numpy as np

os.environ.setdefault("NUMBA_CACHE_DIR", tempfile.mkdtemp(prefix="r01_numba_"))
import numba
from numba import njit, prange

HERE = Path(__file__).resolve().parent
EMU = HERE.parents[1] / "2026-10-11_FINITE_STATE_CYCLE_SCALING" / "CODE" / "cycle_scaling_profiles.py"
spec = importlib.util.spec_from_file_location("cycle_scaling_profiles", EMU); emu = importlib.util.module_from_spec(spec); spec.loader.exec_module(emu)
f_scalar, A_INT, profile_params, initial = emu.f_scalar, emu.A_INT, emu.profile_params, emu.initial
IDX_A = np.array([2, 5, 7, 0, 6, 1, 3, 4], np.int64)
IDX_B = np.array([6, 3, 1, 6, 7, 0, 5, 2], np.int64)
TWO_PI = 2.0 * math.pi


def phases(family, profile):
    raw = hashlib.shake_256(b"R01-B3-PHASES-2026-11-15" + bytes([family]) + profile.to_bytes(4, "big")).digest(64)
    return TWO_PI * (np.frombuffer(raw, dtype=">u8").astype(np.float64) / 2.0**64)


@njit(cache=True)
def g(fam, i, z, P, PH):
    if fam == 1:
        return f_scalar(i, z, P)
    if fam == 2:
        return math.sin(TWO_PI * z[IDX_A[i]] + PH[i])
    if fam == 3:
        return math.sin(TWO_PI * (z[IDX_A[i]] - z[IDX_B[i]]) + PH[i])
    return 0.0


@njit(cache=True)
def mod1(x):
    return x - math.floor(x)


@njit(cache=True)
def qround(x, rule):
    if rule == 0:
        return int(math.floor(x + 0.5))
    return int(math.floor(x))


@njit(cache=True)
def step_F(s, o, qbits, fam, scale, rule, P, PH):
    Q = 1 << qbits; z = np.empty(8); y = np.empty(8)
    for i in range(8): z[i] = s[i] / Q
    for i in range(8):
        fv = g(fam, i, z, P, PH)
        z[i] = mod1(z[i] + (scale * P[i]) * fv)
    for i in range(8):
        t = 0.0
        for j in range(8): t += A_INT[i, j] * z[j]
        y[i] = mod1(t)
    for i in range(7, -1, -1):
        fv = g(fam, i, y, P, PH)
        y[i] = mod1(y[i] - (scale * P[i]) * fv)
    for i in range(8): o[i] = qround(y[i] * Q, rule) % Q


@njit(cache=True)
def shear_fwd(z, qbits, fam, scale, P, PH, rule):
    Q = 1 << qbits; zf = np.empty(8)
    for i in range(8): zf[i] = z[i] / Q
    for i in range(8):
        inc = qround(Q * (scale * P[i]) * g(fam, i, zf, P, PH), rule)
        z[i] = (z[i] + inc) % Q; zf[i] = z[i] / Q


@njit(cache=True)
def shear_bwd_exact(y, qbits, fam, scale, P, PH):
    Q = 1 << qbits; zf = np.empty(8)
    for i in range(8): zf[i] = y[i] / Q
    for i in range(7, -1, -1):
        inc = qround(Q * (scale * P[i]) * g(fam, i, zf, P, PH), 0)
        y[i] = (y[i] - inc) % Q; zf[i] = y[i] / Q


@njit(cache=True)
def shear_bwd_quantized(y, qbits, fam, scale, P, PH, rule):
    Q = 1 << qbits; zf = np.empty(8)
    for i in range(8): zf[i] = y[i] / Q
    for i in range(7, -1, -1):
        inc = qround(-Q * (scale * P[i]) * g(fam, i, zf, P, PH), rule)
        y[i] = (y[i] + inc) % Q; zf[i] = y[i] / Q


@njit(cache=True)
def core(z, y, qbits):
    Q = 1 << qbits
    for i in range(8):
        t = 0
        for j in range(8): t += A_INT[i, j] * z[j]
        y[i] = t % Q


@njit(cache=True)
def decode(idx, qbits, s):
    m = (1 << qbits) - 1
    for i in range(8): s[i] = (idx >> (qbits * i)) & m


@njit(cache=True)
def encode(s, qbits):
    v = 0
    for i in range(8): v |= s[i] << (qbits * i)
    return v


# kind: 0 F, 1 G, 2 S (stepwise), 3 A, 4 H_Q, 5 H_Q^-1
@njit(parallel=True, cache=True)
def successors(qbits, kind, fam, scale, rule, P, PH):
    N = 1 << (8 * qbits); nxt = np.empty(N, np.int64); bad = 0
    for idx in prange(N):
        s = np.empty(8, np.int64); o = np.empty(8, np.int64)
        decode(idx, qbits, s)
        if kind == 0:
            step_F(s, o, qbits, fam, scale, rule, P, PH)
        elif kind == 1:
            shear_fwd(s, qbits, fam, scale, P, PH, 0); core(s, o, qbits); shear_bwd_exact(o, qbits, fam, scale, P, PH)
        elif kind == 2:
            shear_fwd(s, qbits, fam, scale, P, PH, rule); core(s, o, qbits); shear_bwd_quantized(o, qbits, fam, scale, P, PH, rule)
        elif kind == 3:
            core(s, o, qbits)
        elif kind == 4:
            shear_fwd(s, qbits, fam, scale, P, PH, 0)
            for i in range(8): o[i] = s[i]
        else:
            shear_bwd_exact(s, qbits, fam, scale, P, PH)
            for i in range(8): o[i] = s[i]
        for i in range(8):
            if o[i] < 0 or o[i] >= (1 << qbits): bad += 1
        nxt[idx] = encode(o, qbits)
    return nxt, bad


@njit(cache=True)
def graph(nxt):
    N = nxt.shape[0]; indeg = np.zeros(N, np.int64)
    for v in range(N): indeg[nxt[v]] += 1
    maxd = 0
    for v in range(N):
        if indeg[v] > maxd: maxd = indeg[v]
    hist = np.zeros(maxd + 1, np.int64)
    for v in range(N): hist[indeg[v]] += 1
    order = np.empty(N, np.int64); head = 0; tail = 0; deg = indeg.copy()
    for v in range(N):
        if deg[v] == 0: order[tail] = v; tail += 1
    while head < tail:
        v = order[head]; head += 1; w = nxt[v]; deg[w] -= 1
        if deg[w] == 0: order[tail] = w; tail += 1
    cyclic = N - tail
    cyclen = np.zeros(N, np.int64); cid = np.full(N, -1, np.int64); seen = np.zeros(N, np.bool_)
    lens = np.empty(max(cyclic, 1), np.int64); nl = 0
    for v in range(N):
        if deg[v] > 0 and not seen[v]:
            L = 0; w = v
            while not seen[w]: seen[w] = True; w = nxt[w]; L += 1
            w = v
            for _ in range(L): cyclen[w] = L; cid[w] = nl; w = nxt[w]
            lens[nl] = L; nl += 1
    depth = np.zeros(N, np.int64)
    for k in range(tail - 1, -1, -1):
        v = order[k]; w = nxt[v]; depth[v] = depth[w] + 1; cyclen[v] = cyclen[w]; cid[v] = cid[w]
    maxtail = 0; rho = 0.0
    for v in range(N):
        rho += depth[v] + cyclen[v]
        if depth[v] > maxtail: maxtail = depth[v]
    th = np.zeros(maxtail + 1, np.int64)
    for v in range(N): th[depth[v]] += 1
    basin = np.zeros(nl, np.int64)
    for v in range(N): basin[cid[v]] += 1
    return hist, cyclic, lens[:nl], rho / N, maxtail, th, basin, cyclen


@njit(cache=True)
def brent_F(x0, qbits, fam, scale, rule, P, PH):
    tort = x0.copy(); hare = np.empty(8, np.int64); tmp = np.empty(8, np.int64)
    step_F(x0, hare, qbits, fam, scale, rule, P, PH)
    power = 1; lam = 1
    while True:
        same = True
        for i in range(8):
            if tort[i] != hare[i]: same = False; break
        if same: break
        if power == lam:
            tort[:] = hare; power *= 2; lam = 0
        step_F(hare, tmp, qbits, fam, scale, rule, P, PH); hare[:] = tmp; lam += 1
    tort[:] = x0; hare[:] = x0
    for _ in range(lam):
        step_F(hare, tmp, qbits, fam, scale, rule, P, PH); hare[:] = tmp
    mu = 0
    while True:
        same = True
        for i in range(8):
            if tort[i] != hare[i]: same = False; break
        if same: break
        step_F(tort, tmp, qbits, fam, scale, rule, P, PH); tort[:] = tmp
        step_F(hare, tmp, qbits, fam, scale, rule, P, PH); hare[:] = tmp; mu += 1
    return mu, lam


@njit(parallel=True, cache=True)
def brent_batch(X0, qbits, fam, scale, rule, PS, PHS):
    n = X0.shape[0]; MU = np.zeros(n, np.int64); LAM = np.zeros(n, np.int64)
    for k in prange(n):
        mu, lam = brent_F(X0[k].copy(), qbits, fam, scale, rule, PS[k], PHS[k])
        MU[k] = mu; LAM[k] = lam
    return MU, LAM


def summary(nxt, N):
    hist, cyclic, lens, rho, maxtail, th, basin, cyclen = graph(nxt)
    u, c = np.unique(lens, return_counts=True)
    rm = math.sqrt(math.pi * N / 2)
    return {"image_fraction": float(1 - hist[0] / N), "indegree_hist": {int(k): int(v) for k, v in enumerate(hist) if v},
            "states_with_2plus_predecessors": int(hist[2:].sum()), "max_indegree": int(len(hist) - 1),
            "is_permutation": bool(hist[0] == 0), "cyclic_states": int(cyclic), "cycles": int(len(lens)),
            "cycle_type": {int(a): int(b) for a, b in zip(u, c)} if len(u) <= 400 else {"distinct_lengths": int(len(u)), "longest": int(lens.max()), "shortest": int(lens.min())},
            "longest_cycle": int(lens.max()), "shortest_cycle": int(lens.min()), "mean_rho": float(rho), "mean_rho_over_sqrt_piN_2": float(rho / rm),
            "cyclic_over_sqrt_piN_2": float(cyclic / rm), "max_tail": int(maxtail), "mean_tail": float(np.dot(np.arange(len(th)), th) / N),
            "tail_hist_head": [int(x) for x in th[:12]], "largest_basin_fraction": float(basin.max() / N),
            "basin_quantiles": [int(x) for x in np.quantile(basin, [0, 0.5, 0.9, 1.0])],
            "mean_cycle_length_of_random_state": float(cyclen.mean()) if hist[0] == 0 else None}, cyclen


def part1(out, profiles, qs, log):
    rows = []
    for q in qs:
        N = 1 << (8 * q); Z = np.zeros(17)
        nA, badA = successors(q, 3, 0, 0.0, 0, Z, np.zeros(8)); sA, cycA = summary(nA, N)
        rows.append({"q": q, "map": "A", "summary": sA, "out_of_range": int(badA)})
        for prof in profiles:
            P = profile_params(prof); P[16] = 0.0
            conds = [("C0", 0, 1.0), ("C1", 1, 1.0), ("C2", 2, 1.0), ("C3", 3, 1.0), ("C4", 1, 0.25), ("C4", 1, 0.5), ("C4", 1, 2.0)]
            for cname, fam, sc in conds:
                PH = phases(fam, prof); t0 = time.time()
                rec = {"q": q, "profile": prof, "condition": cname, "family": fam, "scale": sc, "maps": {}}
                kinds = [("F_nearest", 0, 0)]
                if cname in ("C0", "C1"): kinds.append(("F_truncation", 0, 1))
                if cname != "C0": kinds.append(("G", 1, 0))
                if cname == "C1": kinds.append(("S_nearest", 2, 0))
                if cname in ("C1", "C2", "C3"): kinds.append(("S_truncation", 2, 1))
                keep = {}
                for name, kind, rule in kinds:
                    nx, bad = successors(q, kind, fam, sc, rule, P, PH)
                    sm, cyc = summary(nx, N); sm["out_of_range"] = int(bad)
                    if cname == "C0": sm["equals_A_on_every_state"] = bool(np.array_equal(nx, nA))
                    rec["maps"][name] = sm
                    if name in ("G", "S_nearest", "F_nearest"): keep[name] = (nx, cyc)
                if "G" in keep:
                    H, b1 = successors(q, 4, fam, sc, 0, P, PH); Hi, b2 = successors(q, 5, fam, sc, 0, P, PH)
                    nG, cycG = keep["G"]
                    ident = np.arange(N, dtype=np.int64)
                    rec["G_checks"] = {"H_is_bijection": bool(np.unique(H).size == N), "Hinv_of_H_is_identity": bool(np.array_equal(Hi[H], ident)),
                                       "G_equals_Hinv_A_H": bool(np.array_equal(nG, Hi[nA[H]])),
                                       "cycle_type_equals_A": rec["maps"]["G"]["cycle_type"] == sA["cycle_type"],
                                       "every_period_equals_A_period_of_H": bool(np.array_equal(cycG, cycA[H])),
                                       "out_of_range": int(b1 + b2)}
                    if "S_nearest" in keep:
                        rec["S_nearest_equals_G_fraction"] = float(np.mean(keep["S_nearest"][0] == nG))
                    rec["F_equals_G_fraction"] = float(np.mean(keep["F_nearest"][0] == nG))
                rec["seconds"] = round(time.time() - t0, 1); rows.append(rec)
                msg = "Q%d p%d %s fam%d s%.2f %s" % (q, prof, cname, fam, sc, " ".join("%s:if=%.4f,rho=%.3f,cyc=%d" % (k, v["image_fraction"], v["mean_rho_over_sqrt_piN_2"], v["cycles"]) for k, v in rec["maps"].items()))
                if "G_checks" in rec: msg += " Gchecks=" + str(all(v is True or v == 0 for v in rec["G_checks"].values()))
                print(msg, flush=True); log.write(msg + "\n"); log.flush()
            (out / "part1_partial.json").write_text(json.dumps(rows, indent=0))
    return rows


def part2(out, profiles, qs, nstarts, log):
    conds = [("fam1_s1_nearest", 1, 1.0, 0), ("fam2_s1_nearest", 2, 1.0, 0), ("fam3_s1_nearest", 3, 1.0, 0),
             ("fam1_s1_truncation", 1, 1.0, 1), ("fam1_s0.5_nearest", 1, 0.5, 0), ("fam1_s2_nearest", 1, 2.0, 0)]
    rows = []
    for q in qs:
        for cname, fam, sc, rule in conds:
            t0 = time.time()
            tasks = [(prof, run) for prof in profiles for run in range(nstarts)]
            X0 = np.array([initial(q, run) for _, run in tasks], np.int64)
            PS = np.array([profile_params(prof) for prof, _ in tasks]); PS[:, 16] = 0.0
            PHS = np.array([phases(fam, prof) for prof, _ in tasks])
            MU, LAM = brent_batch(X0, q, fam, sc, rule, PS, PHS)
            for (prof, run), mu, lam in zip(tasks, MU, LAM):
                rows.append({"q": q, "condition": cname, "profile": prof, "run": run, "tail": int(mu), "cycle": int(lam), "rho": int(mu + lam)})
            sub = [r for r in rows if r["q"] == q and r["condition"] == cname]
            msg = "part2 Q%d %s mean rho %.1f ratio %.3f (%.0f s)" % (q, cname, np.mean([r["rho"] for r in sub]), np.mean([r["rho"] for r in sub]) / math.sqrt(math.pi * 2.0 ** (8 * q) / 2), time.time() - t0)
            print(msg, flush=True); log.write(msg + "\n"); log.flush()
            (out / "part2_partial.json").write_text(json.dumps(rows, indent=0))
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--profiles", default="0:10")
    ap.add_argument("--q1", default="2,3"); ap.add_argument("--q2", default="4,5"); ap.add_argument("--starts", type=int, default=8)
    ap.add_argument("--smoke", action="store_true"); ap.add_argument("--part", default="12")
    a = ap.parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    p0, p1 = (int(v) for v in a.profiles.split(":")); profiles = list(range(p0, p1))
    log = open(out / "run_log.txt", "a", encoding="utf-8")
    plan = HERE.parent / "PLAN.md"
    head = {"started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "plan_sha256": hashlib.sha256(plan.read_bytes()).hexdigest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "emulator_sha256": hashlib.sha256(EMU.read_bytes()).hexdigest(),
            "argv": sys.argv[1:], "environment": {"platform": platform.platform(), "python": platform.python_version(), "numpy": np.__version__,
                                                 "numba": numba.__version__, "threads": numba.get_num_threads()}}
    log.write(json.dumps(head) + "\n"); log.flush()
    res = dict(head)
    if "1" in a.part:
        res["part1"] = part1(out, profiles, [int(v) for v in a.q1.split(",") if v], log)
    if "2" in a.part:
        res["part2"] = part2(out, profiles, [int(v) for v in a.q2.split(",") if v], a.starts, log)
    res["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    name = "RESULTS_fg_study%s_part%s.json" % ("_smoke" if a.smoke else "", a.part)
    (out / name).write_text(json.dumps(res, indent=1)); log.write("finished " + res["finished_utc"] + "\n"); log.close()


if __name__ == "__main__":
    main()
