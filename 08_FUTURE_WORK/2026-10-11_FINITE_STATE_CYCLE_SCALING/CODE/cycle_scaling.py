"""Uncapped functional-graph measurement of the quantized-state HCHM-8 emulator.

Same update rule as quant_step() in 05_CODE/run_validated_experiments.py
(reference profile, Float64 transcendental evaluation, round-to-nearest to Q bits
per coordinate), but run until the orbit closes, so tail length (mu) and cycle
length (lam) are measured instead of right-censored.

Modes
  hchm : full T = H^-1 A H, quantized after each step
  core : shear amplitudes set to zero (pure toral automorphism A on the 2^Q lattice)

Usage: python cycle_scaling.py <qbits> <runs> <mode> <cap> <out.csv> [workers]
"""
import csv, hashlib, math, sys, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from numba import njit

REF_LABEL = b"GIU-PhD-HCHM8-heterogeneous-conjugacy-v1"
SEED_LABEL = b"R01-CYCLE-SCALING-2026-10-11"

def _unit(label, count):
    raw = hashlib.shake_256(label).digest(8 * count)
    return np.frombuffer(raw, dtype=">u8").astype(np.float64) / 2.0**64

_r = _unit(REF_LABEL, 32)
EPS_REF = (0.07 + 0.08 * _r[8:16]).astype(np.float64)
PHI1 = 2 * np.pi * _r[16]; PHI2 = 2 * np.pi * _r[17]
RHO3 = 0.75 + 2.25 * _r[18]; K4 = 1.25 + 2.75 * _r[19]
ETA6 = 0.75 + 1.75 * _r[20]; K7 = 1.25 + 2.75 * _r[21]; RHO8 = 1.0 + 5.0 * _r[22]
I0ETA = float(np.i0(ETA6))
A_INT = np.array([[2,1,0,0,0,0,0,0],[1,1,0,0,0,0,0,0],[0,0,3,2,0,0,0,0],[0,0,1,1,0,0,0,0],
                  [0,0,0,0,5,2,0,0],[0,0,0,0,2,1,0,0],[0,0,0,0,0,0,7,3],[0,0,0,0,0,0,2,1]], dtype=np.int64)

@njit(cache=True)
def mod1s(x):
    return x - math.floor(x)

@njit(cache=True)
def f_scalar(i, z):
    if i == 0: return math.sin(2*math.pi*z[2] + PHI1)
    if i == 1:
        c = math.cos(2*math.pi*z[5] + PHI2); return 4*c*c*c - 3*c
    if i == 2: return math.sin(2*math.pi*(z[7]-z[1])) / (1 + RHO3*math.cos(2*math.pi*z[1])**2)
    if i == 3: return math.tanh(K4*math.sin(2*math.pi*(z[0]-z[6])))
    if i == 4:
        q = .5*(1 + math.sin(2*math.pi*(z[6]+z[7]))); return 4*q*(1-q) - .5
    if i == 5: return (math.exp(ETA6*math.cos(2*math.pi*(z[1]-z[0]))) - I0ETA) / (math.exp(ETA6) - math.exp(-ETA6))
    if i == 6: return (2/math.pi)*math.atan(K7*math.sin(2*math.pi*(z[3]+z[5])))
    return math.log(1 + RHO8*math.sin(math.pi*(z[4]-z[2]))**2) / math.log(1 + RHO8) - .5

@njit(cache=True)
def step(state, out, z, y, qbits, eps):
    Q = 1 << qbits
    for i in range(8): z[i] = state[i] / Q
    for i in range(8): z[i] = mod1s(z[i] + eps[i]*f_scalar(i, z))
    for i in range(8):
        s = 0.0
        for j in range(8): s += A_INT[i, j]*z[j]
        y[i] = mod1s(s)
    for i in range(7, -1, -1): y[i] = mod1s(y[i] - eps[i]*f_scalar(i, y))
    for i in range(8): out[i] = int(math.floor(y[i]*Q + .5)) % Q

@njit(cache=True)
def same(a, b):
    for i in range(8):
        if a[i] != b[i]: return False
    return True

@njit(cache=True)
def brent(x0, qbits, eps, cap):
    z = np.empty(8); y = np.empty(8); tmp = np.empty(8, np.int64)
    tort = x0.copy(); hare = np.empty(8, np.int64); step(x0, hare, z, y, qbits, eps)
    power = 1; lam = 1; steps = 1
    while not same(tort, hare):
        if steps >= cap: return -1, -1, steps, tort
        if power == lam:
            tort[:] = hare; power *= 2; lam = 0
        step(hare, tmp, z, y, qbits, eps); hare[:] = tmp; lam += 1; steps += 1
    tort[:] = x0; hare[:] = x0
    for _ in range(lam):
        step(hare, tmp, z, y, qbits, eps); hare[:] = tmp
    mu = 0
    while not same(tort, hare):
        step(tort, tmp, z, y, qbits, eps); tort[:] = tmp
        step(hare, tmp, z, y, qbits, eps); hare[:] = tmp; mu += 1
    # canonical cycle representative = lexicographically smallest state on the cycle
    best = tort.copy(); cur = tort.copy()
    for _ in range(lam):
        step(cur, tmp, z, y, qbits, eps); cur[:] = tmp
        for i in range(8):
            if cur[i] < best[i]:
                best[:] = cur; break
            if cur[i] > best[i]: break
    return mu, lam, steps, best

def initial(qbits, run):
    raw = hashlib.shake_256(SEED_LABEL + bytes([qbits]) + run.to_bytes(4, "big")).digest(64)
    v = np.frombuffer(raw, dtype=">u8").astype(object)
    return np.array([int(x) % (1 << qbits) for x in v], dtype=np.int64)

def one(args):
    qbits, run, mode, cap = args
    eps = EPS_REF if mode == "hchm" else np.zeros(8)
    t = time.time()
    mu, lam, steps, rep = brent(initial(qbits, run), qbits, eps, cap)
    return dict(qbits=qbits, run=run, mode=mode, cap=cap, closed=int(lam >= 0), tail=mu, cycle=lam,
                rho=(mu + lam if lam >= 0 else -1), brent_steps=steps,
                cycle_id=hashlib.sha256(rep.tobytes()).hexdigest()[:16] if lam >= 0 else "", seconds=round(time.time()-t, 2))

if __name__ == "__main__":
    qbits, runs, mode, cap, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5]
    workers = int(sys.argv[6]) if len(sys.argv) > 6 else 12
    one((2, 0, mode, 1000))  # compile once before forking workers
    with ProcessPoolExecutor(workers) as ex:
        rows = list(ex.map(one, [(qbits, r, mode, cap) for r in range(runs)]))
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    closed = [r for r in rows if r["closed"]]
    print(mode, "Q", qbits, "closed", len(closed), "/", runs,
          "median tail", np.median([r["tail"] for r in closed]) if closed else None,
          "median cycle", np.median([r["cycle"] for r in closed]) if closed else None,
          "distinct cycles", len({r["cycle_id"] for r in closed}), "max s", max(r["seconds"] for r in rows))
