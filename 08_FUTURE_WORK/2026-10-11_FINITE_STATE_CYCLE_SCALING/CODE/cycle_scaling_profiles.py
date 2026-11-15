"""Uncapped cycle measurement of the quantized-state HCHM-8 emulator over several parameter profiles.

Extends cycle_scaling.py (reference profile only). Profile 0 is the reference profile;
profile k >= 1 uses the parameters of HCHM8(key_for_index(9000 + k)) exactly as derived in
05_CODE/run_validated_experiments.py. Initial states are the same lattice points as in
cycle_scaling.py, so profile 0 must reproduce the earlier CSV files row for row.

Usage: python cycle_scaling_profiles.py <qbits> <first_profile> <n_profiles> <runs> <out.csv> [workers] [perturb]

With a seventh argument "perturb", every transcendental shear value is multiplied by (1 + s*k*2^-52) with k in -4..4 chosen by a
hash of the lattice state and s = 1 ("perturb=<s>" sets another scale). This imitates a mathematical library whose results differ in the last bits, and is used to test
whether the measured orbits depend on the library.
"""
import csv, hashlib, math, os, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
import numpy as np
from numba import njit

REF_LABEL = b"GIU-PhD-HCHM8-heterogeneous-conjugacy-v1"
SEED_LABEL = b"R01-CYCLE-SCALING-2026-10-11"
A_INT = np.array([[2,1,0,0,0,0,0,0],[1,1,0,0,0,0,0,0],[0,0,3,2,0,0,0,0],[0,0,1,1,0,0,0,0],
                  [0,0,0,0,5,2,0,0],[0,0,0,0,2,1,0,0],[0,0,0,0,0,0,7,3],[0,0,0,0,0,0,2,1]], dtype=np.int64)


def _unit(label, count):
    raw = hashlib.shake_256(label).digest(8 * count)
    return np.frombuffer(raw, dtype=">u8").astype(np.float64) / 2.0**64


def profile_params(k):
    """Returns [eps1..eps8, phi1, phi2, rho3, kappa4, eta6, kappa7, rho8, I0(eta6), perturbation switch]."""
    if k == 0:
        r = _unit(REF_LABEL, 32)
    else:
        key = hashlib.shake_256(b"HCHM8-CONJUGATE-KEY" + (9000 + k).to_bytes(8, "big")).digest(64)
        r = _unit(key + b"|HCHM8-PARAM|" + (0).to_bytes(4, "big"), 32)
    p = np.zeros(17)
    p[0:8] = 0.07 + 0.08 * r[8:16]
    p[8] = 2 * np.pi * r[16]; p[9] = 2 * np.pi * r[17]
    p[10] = 0.75 + 2.25 * r[18]; p[11] = 1.25 + 2.75 * r[19]
    p[12] = 0.75 + 1.75 * r[20]; p[13] = 1.25 + 2.75 * r[21]; p[14] = 1.0 + 5.0 * r[22]
    p[15] = float(np.i0(p[12]))
    return p


@njit(cache=True)
def mod1s(x):
    return x - math.floor(x)


@njit(cache=True)
def f_scalar(i, z, P):
    if i == 0: return math.sin(2*math.pi*z[2] + P[8])
    if i == 1:
        c = math.cos(2*math.pi*z[5] + P[9]); return 4*c*c*c - 3*c
    if i == 2: return math.sin(2*math.pi*(z[7]-z[1])) / (1 + P[10]*math.cos(2*math.pi*z[1])**2)
    if i == 3: return math.tanh(P[11]*math.sin(2*math.pi*(z[0]-z[6])))
    if i == 4:
        q = .5*(1 + math.sin(2*math.pi*(z[6]+z[7]))); return 4*q*(1-q) - .5
    if i == 5: return (math.exp(P[12]*math.cos(2*math.pi*(z[1]-z[0]))) - P[15]) / (math.exp(P[12]) - math.exp(-P[12]))
    if i == 6: return (2/math.pi)*math.atan(P[13]*math.sin(2*math.pi*(z[3]+z[5])))
    return math.log(1 + P[14]*math.sin(math.pi*(z[4]-z[2]))**2) / math.log(1 + P[14]) - .5


@njit(cache=True)
def step(state, out, z, y, qbits, P):
    Q = 1 << qbits
    for i in range(8): z[i] = state[i] / Q
    h = 0
    if P[16] != 0.0:
        h = 1469598103934665603
        for j in range(8): h = (h ^ state[j]) * 1099511628211
    for i in range(8):
        fv = f_scalar(i, z, P)
        if P[16] != 0.0: fv *= 1.0 + P[16] * (((h >> (4 * i)) & 15) % 9 - 4) * 2.220446049250313e-16
        z[i] = mod1s(z[i] + P[i]*fv)
    for i in range(8):
        s = 0.0
        for j in range(8): s += A_INT[i, j]*z[j]
        y[i] = mod1s(s)
    for i in range(7, -1, -1):
        fv = f_scalar(i, y, P)
        if P[16] != 0.0: fv *= 1.0 + P[16] * (((h >> (4 * (i + 8))) & 15) % 9 - 4) * 2.220446049250313e-16
        y[i] = mod1s(y[i] - P[i]*fv)
    for i in range(8): out[i] = int(math.floor(y[i]*Q + .5)) % Q


@njit(cache=True)
def same(a, b):
    for i in range(8):
        if a[i] != b[i]: return False
    return True


@njit(cache=True)
def brent(x0, qbits, P):
    z = np.empty(8); y = np.empty(8); tmp = np.empty(8, np.int64)
    tort = x0.copy(); hare = np.empty(8, np.int64); step(x0, hare, z, y, qbits, P)
    power = 1; lam = 1
    while not same(tort, hare):
        if power == lam:
            tort[:] = hare; power *= 2; lam = 0
        step(hare, tmp, z, y, qbits, P); hare[:] = tmp; lam += 1
    tort[:] = x0; hare[:] = x0
    for _ in range(lam):
        step(hare, tmp, z, y, qbits, P); hare[:] = tmp
    mu = 0
    while not same(tort, hare):
        step(tort, tmp, z, y, qbits, P); tort[:] = tmp
        step(hare, tmp, z, y, qbits, P); hare[:] = tmp; mu += 1
    best = tort.copy(); cur = tort.copy()
    for _ in range(lam):
        step(cur, tmp, z, y, qbits, P); cur[:] = tmp
        for i in range(8):
            if cur[i] < best[i]:
                best[:] = cur; break
            if cur[i] > best[i]: break
    return mu, lam, best


def initial(qbits, run):
    raw = hashlib.shake_256(SEED_LABEL + bytes([qbits]) + run.to_bytes(4, "big")).digest(64)
    v = np.frombuffer(raw, dtype=">u8").astype(object)
    return np.array([int(x) % (1 << qbits) for x in v], dtype=np.int64)


def one(args):
    qbits, prof, run, pert = args
    t = time.time()
    P = profile_params(prof); P[16] = pert
    mu, lam, rep = brent(initial(qbits, run), qbits, P)
    return dict(qbits=qbits, profile=prof, run=run, tail=mu, cycle=lam, rho=mu + lam,
                cycle_id=hashlib.sha256(rep.tobytes()).hexdigest()[:16], seconds=round(time.time() - t, 2))


if __name__ == "__main__":
    qbits, first, nprof, runs, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    workers = int(sys.argv[6]) if len(sys.argv) > 6 else 12
    pert = 0.0
    if len(sys.argv) > 7 and sys.argv[7].startswith("perturb"):  # "perturb" = last-bit scale; "perturb=<scale>" for a sanity check
        pert = float(sys.argv[7].split("=")[1]) if "=" in sys.argv[7] else 1.0
    one((2, 0, 0, pert))  # compile before the workers start
    # submitted run-major so that every profile gets its first start before any gets a second one
    jobs = [(qbits, pr, r, pert) for r in range(runs) for pr in range(first, first + nprof)]
    for attempt in range(8):  # Windows process pools occasionally fail at start-up; results are deterministic
        try:
            with ProcessPoolExecutor(workers) as ex:
                rows = []
                with open(out + ".partial", "w", newline="") as pf:  # finished runs are kept if the job is interrupted
                    pw = None
                    for fut in as_completed([ex.submit(one, j) for j in jobs]):
                        r = fut.result(); rows.append(r)
                        if pw is None:
                            pw = csv.DictWriter(pf, fieldnames=list(r)); pw.writeheader()
                        pw.writerow(r); pf.flush()
            break
        except (BrokenProcessPool, OSError):
            time.sleep(2)
    else:
        raise SystemExit("process pool failed repeatedly")
    rows.sort(key=lambda r: (r["profile"], r["run"]))
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    os.remove(out + ".partial")
    print("Q", qbits, "profiles", first, "-", first + nprof - 1, "runs/profile", runs, "rows", len(rows),
          "median rho", float(np.median([r["rho"] for r in rows])), "min cycle", min(r["cycle"] for r in rows),
          "max s", max(r["seconds"] for r in rows))
