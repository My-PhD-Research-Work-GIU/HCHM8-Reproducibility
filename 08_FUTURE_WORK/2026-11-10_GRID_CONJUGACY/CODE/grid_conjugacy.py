"""Grid-exact conjugate of HCHM-8 against the rounded emulator and the linear core (plan: ../PLAN.md, fixed before the runs).

Three maps on the lattice L_Q = (Z / 2^Q Z)^8 are compared for the parameter profiles of the finite-state section:

  F_Q  the emulator of the paper: lattice -> Float64 -> T = H^-1 A H -> round to the nearest lattice point
       (step() of 08_FUTURE_WORK/2026-10-11_FINITE_STATE_CYCLE_SCALING/CODE/cycle_scaling_profiles.py, imported unchanged);
  G_Q  H_Q^-1 A_Q H_Q, where H_Q applies the eight shears on the lattice with each increment rounded to a lattice step,
       z_i <- z_i + round(2^Q eps_i g_i(z / 2^Q)) mod 2^Q, and A_Q is the integer core modulo 2^Q. Since g_i does not depend
       on z_i, each shear is a bijection of L_Q and is undone exactly by subtracting the same increment, so G_Q is conjugate
       to A_Q by the bijection H_Q;
  A_Q  the integer core modulo 2^Q (all shear amplitudes zero).

Part 1 (Q = 2, 3; every lattice state): functional-graph statistics of F_Q, G_Q and A_Q: image fraction, number of cyclic
states, number of cycles, mean over all starting states of tail + cycle (rho), cycle type; F_Q = G_Q agreement.
Part 2 (Q = 4 to 7; the starting states of the uncapped measurement): orbit period of G_Q, compared with the period of the A_Q
orbit of H_Q(x0) and with 3 * 2^(Q-1).

Usage: python grid_conjugacy.py <out_dir> [--profiles a:b] [--exhaustive-q 2,3] [--orbit-q 4,5,6,7] [--smoke]
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
f_scalar, step_F, A_INT, profile_params, initial = emu.f_scalar, emu.step, emu.A_INT, emu.profile_params, emu.initial


@njit(cache=True)
def step_G(state, out, zf, qbits, P):
    Q = 1 << qbits
    z = state.copy()
    for i in range(8): zf[i] = z[i] / Q
    for i in range(8):
        inc = int(math.floor(Q * P[i] * f_scalar(i, zf, P) + 0.5))
        z[i] = (z[i] + inc) % Q; zf[i] = z[i] / Q
    y = np.empty(8, np.int64)
    for i in range(8):
        s = 0
        for j in range(8): s += A_INT[i, j] * z[j]
        y[i] = s % Q
    for i in range(8): zf[i] = y[i] / Q
    for i in range(7, -1, -1):
        inc = int(math.floor(Q * P[i] * f_scalar(i, zf, P) + 0.5))
        y[i] = (y[i] - inc) % Q; zf[i] = y[i] / Q
    for i in range(8): out[i] = y[i]


@njit(cache=True)
def H_Q(state, qbits, P):
    Q = 1 << qbits; z = state.copy(); zf = np.empty(8)
    for i in range(8): zf[i] = z[i] / Q
    for i in range(8):
        inc = int(math.floor(Q * P[i] * f_scalar(i, zf, P) + 0.5))
        z[i] = (z[i] + inc) % Q; zf[i] = z[i] / Q
    return z


@njit(cache=True)
def step_A(state, out, qbits):
    Q = 1 << qbits
    for i in range(8):
        s = 0
        for j in range(8): s += A_INT[i, j] * state[j]
        out[i] = s % Q


@njit(cache=True)
def decode(idx, qbits, s):
    m = (1 << qbits) - 1
    for i in range(8): s[i] = (idx >> (qbits * i)) & m


@njit(cache=True)
def encode(s, qbits):
    v = 0
    for i in range(8): v |= s[i] << (qbits * i)
    return v


@njit(parallel=True, cache=True)
def next_array(qbits, P, which):
    N = 1 << (8 * qbits); nxt = np.empty(N, np.int64)
    for idx in prange(N):
        s = np.empty(8, np.int64); o = np.empty(8, np.int64); z = np.empty(8); y = np.empty(8)
        decode(idx, qbits, s)
        if which == 0: step_F(s, o, z, y, qbits, P)
        elif which == 1: step_G(s, o, z, qbits, P)
        else: step_A(s, o, qbits)
        nxt[idx] = encode(o, qbits)
    return nxt


@njit(cache=True)
def graph_stats(nxt):
    N = nxt.shape[0]; indeg = np.zeros(N, np.int64)
    for v in range(N): indeg[nxt[v]] += 1
    image = 0
    for v in range(N):
        if indeg[v] > 0: image += 1
    order = np.empty(N, np.int64); head = 0; tail = 0; deg = indeg.copy()
    for v in range(N):
        if deg[v] == 0: order[tail] = v; tail += 1
    while head < tail:
        v = order[head]; head += 1; w = nxt[v]; deg[w] -= 1
        if deg[w] == 0: order[tail] = w; tail += 1
    cyclic = N - tail
    cyclen = np.zeros(N, np.int64); ncyc = 0; seen = np.zeros(N, np.bool_)
    lens = np.empty(cyclic if cyclic > 0 else 1, np.int64); nl = 0
    for v in range(N):
        if deg[v] > 0 and not seen[v]:
            L = 0; w = v
            while not seen[w]: seen[w] = True; w = nxt[w]; L += 1
            w = v
            for _ in range(L): cyclen[w] = L; w = nxt[w]
            lens[nl] = L; nl += 1; ncyc += 1
    depth = np.zeros(N, np.int64)
    for k in range(tail - 1, -1, -1):
        v = order[k]; w = nxt[v]; depth[v] = depth[w] + 1; cyclen[v] = cyclen[w]
    rho_sum = 0.0; maxtail = 0
    for v in range(N):
        rho_sum += depth[v] + cyclen[v]
        if depth[v] > maxtail: maxtail = depth[v]
    return image, cyclic, ncyc, rho_sum / N, maxtail, lens[:nl]


@njit(cache=True)
def orbit_period(x0, qbits, P, which, cap):
    s = x0.copy(); o = np.empty(8, np.int64); z = np.empty(8); y = np.empty(8)
    for n in range(1, cap + 1):
        if which == 1: step_G(s, o, z, qbits, P)
        else: step_A(s, o, qbits)
        s[:] = o
        same = True
        for i in range(8):
            if s[i] != x0[i]: same = False; break
        if same: return n
    return -1


def cycle_type(lens):
    u, c = np.unique(lens, return_counts=True)
    return {int(a): int(b) for a, b in zip(u, c)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--profiles", default="0:10")
    ap.add_argument("--exhaustive-q", default="2,3"); ap.add_argument("--orbit-q", default="4,5,6,7"); ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args(); out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    p0, p1 = (int(v) for v in a.profiles.split(":")); profiles = list(range(p0, p1))
    eq = [int(v) for v in a.exhaustive_q.split(",") if v]; oq = [int(v) for v in a.orbit_q.split(",") if v]
    starts = {2: 32, 3: 32, 4: 32, 5: 32, 6: 13, 7: 2}
    res = {"environment": {"platform": platform.platform(), "python": platform.python_version(), "numpy": np.__version__, "numba": numba.__version__},
           "emulator_source_sha256": hashlib.sha256(EMU.read_bytes()).hexdigest(), "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "profiles": profiles, "smoke_test": a.smoke, "exhaustive": [], "orbits": [], "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    for q in eq:
        N = 1 << (8 * q); rm = math.sqrt(math.pi * N / 2)
        tA = time.time(); nA = next_array(q, profile_params(0), 2); sA = graph_stats(nA); ctA = cycle_type(sA[5]); tA = time.time() - tA
        for prof in profiles:
            P = profile_params(prof); P[16] = 0.0; t = time.time()
            nF = next_array(q, P, 0); nG = next_array(q, P, 1)
            sF = graph_stats(nF); sG = graph_stats(nG); ctG = cycle_type(sG[5])
            row = {"qbits": q, "state_bits": 8 * q, "N": N, "profile": prof,
                   "F": {"image_fraction": sF[0] / N, "cyclic_states": int(sF[1]), "cycles": int(sF[2]), "mean_rho": sF[3], "max_tail": int(sF[4]),
                         "mean_rho_over_sqrt_piN_2": sF[3] / rm, "cyclic_over_sqrt_piN_2": sF[1] / rm, "longest_cycle": int(sF[5].max())},
                   "G": {"image_fraction": sG[0] / N, "is_permutation": bool(sG[0] == N), "cycles": int(sG[2]), "max_tail": int(sG[4]),
                         "cycle_type": ctG, "cycle_type_equals_A": ctG == ctA, "longest_cycle": int(sG[5].max())},
                   "F_equals_G_fraction": float(np.mean(nF == nG)), "seconds": round(time.time() - t, 1)}
            res["exhaustive"].append(row); print(json.dumps({k: row[k] for k in ("qbits", "profile", "F_equals_G_fraction")}), json.dumps(row["F"]), "G perm", row["G"]["is_permutation"], "type=A", row["G"]["cycle_type_equals_A"], flush=True)
        res["exhaustive"].append({"qbits": q, "A": {"is_permutation": bool(sA[0] == N), "cycle_type": ctA, "longest_cycle": int(sA[5].max()), "seconds": round(tA, 1)}})
    for q in oq:
        period = 12 if q == 2 else 3 * 2 ** (q - 1)
        for prof in profiles:
            P = profile_params(prof); P[16] = 0.0
            for run in range(starts[q]):
                x0 = initial(q, run); pg = orbit_period(x0, q, P, 1, period + 1); pa = orbit_period(H_Q(x0, q, P), q, P, 2, period + 1)
                res["orbits"].append({"qbits": q, "profile": prof, "run": run, "period_G": int(pg), "period_A_of_HQx0": int(pa),
                                      "equal": pg == pa, "divides_core_period": pg > 0 and period % pg == 0})
        sub = [r for r in res["orbits"] if r["qbits"] == q]
        print("Q%d orbits: %d, all equal %s, all divide %d: %s, periods %s" % (q, len(sub), all(r["equal"] for r in sub), period, all(r["divides_core_period"] for r in sub), sorted({r["period_G"] for r in sub})), flush=True)
    res["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    (out / ("RESULTS_grid_conjugacy%s.json" % ("_smoke" if a.smoke else ""))).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
