"""Independent checker for the Q = 2 results of fg_study.py (plan: ../PLAN.md).

Written separately: the maps are evaluated in vectorized numpy from the formulas of Section 3 of the manuscript (not from the
emulator's f_scalar), the functional graph is analysed by another algorithm (iterated colouring in pure Python), and only the
final successor arrays and statistics are compared with fg_study.py. Family 1 (HCHM-8), s = 1, profiles 0 and 1, maps F
(nearest), G and A on all 65,536 states of L_2.

Usage: python independent_check.py <out.json>
"""
import importlib.util, json, math, sys
from pathlib import Path
import numpy as np
from scipy.special import i0

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("fg", HERE / "fg_study.py"); fg = importlib.util.module_from_spec(spec); spec.loader.exec_module(fg)
q = 2; Q = 1 << q; N = Q ** 8
A = np.array([[2, 1, 0, 0, 0, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 3, 2, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0],
              [0, 0, 0, 0, 5, 2, 0, 0], [0, 0, 0, 0, 2, 1, 0, 0], [0, 0, 0, 0, 0, 0, 7, 3], [0, 0, 0, 0, 0, 0, 2, 1]])


def g(i, z, p):  # Section 3, vectorized, 0-based coordinates
    tp = 2 * np.pi
    if i == 0: return np.sin(tp * z[:, 2] + p[8])
    if i == 1: c = np.cos(tp * z[:, 5] + p[9]); return 4 * c ** 3 - 3 * c
    if i == 2: return np.sin(tp * (z[:, 7] - z[:, 1])) / (1 + p[10] * np.cos(tp * z[:, 1]) ** 2)
    if i == 3: return np.tanh(p[11] * np.sin(tp * (z[:, 0] - z[:, 6])))
    if i == 4: qq = 0.5 * (1 + np.sin(tp * (z[:, 6] + z[:, 7]))); return 4 * qq * (1 - qq) - 0.5
    if i == 5: return (np.exp(p[12] * np.cos(tp * (z[:, 1] - z[:, 0]))) - i0(p[12])) / (np.exp(p[12]) - np.exp(-p[12]))
    if i == 6: return (2 / np.pi) * np.arctan(p[13] * np.sin(tp * (z[:, 3] + z[:, 5])))
    return np.log(1 + p[14] * np.sin(np.pi * (z[:, 4] - z[:, 2])) ** 2) / np.log(1 + p[14]) - 0.5


idx = np.arange(N); S = np.stack([(idx >> (q * i)) & (Q - 1) for i in range(8)], axis=1)
enc = lambda M: sum(M[:, i].astype(np.int64) << (q * i) for i in range(8))


def F_map(p):
    z = S / Q
    for i in range(8): z[:, i] = np.mod(z[:, i] + p[i] * g(i, z, p), 1.0)
    y = np.mod(z @ A.T, 1.0)
    for i in range(7, -1, -1): y[:, i] = np.mod(y[:, i] - p[i] * g(i, y, p), 1.0)
    return enc(np.floor(y * Q + 0.5).astype(np.int64) % Q)


def G_map(p):
    z = S.copy()
    for i in range(8): z[:, i] = (z[:, i] + np.floor(Q * p[i] * g(i, z / Q, p) + 0.5).astype(np.int64)) % Q
    y = (z @ A.T) % Q
    for i in range(7, -1, -1): y[:, i] = (y[:, i] - np.floor(Q * p[i] * g(i, y / Q, p) + 0.5).astype(np.int64)) % Q
    return enc(y)


def A_map(): return enc((S @ A.T) % Q)


def graph_stats(nxt):  # iterated colouring: walk each unvisited state until a visited one
    nxt = [int(v) for v in nxt]; state = [0] * N; depth = [0] * N; cyc = [0] * N; lens = []
    for v in range(N):
        if state[v]: continue
        path = []; w = v
        while state[w] == 0: state[w] = 1; path.append(w); w = nxt[w]
        if state[w] == 1:  # new cycle found inside the current path
            k = path.index(w); L = len(path) - k; lens.append(L)
            for u in path[k:]: cyc[u] = L; depth[u] = 0; state[u] = 2
            path = path[:k]
        for u in reversed(path):
            depth[u] = depth[nxt[u]] + 1; cyc[u] = cyc[nxt[u]]; state[u] = 2
    indeg = np.bincount(np.array(nxt), minlength=N)
    u, c = np.unique(lens, return_counts=True)
    return {"indegree_hist": {int(k): int(v) for k, v in enumerate(np.bincount(indeg)) if v}, "cycles": len(lens),
            "cycle_type": {int(a): int(b) for a, b in zip(u, c)}, "mean_rho": float(np.mean(np.array(depth) + np.array(cyc))),
            "max_tail": int(max(depth))}


out = {}
for prof in (0, 1):
    p = fg.profile_params(prof); p[16] = 0.0; PH = np.zeros(8)
    rec = {}
    for name, mine, kind in (("F_nearest", F_map(p), 0), ("G", G_map(p), 1), ("A", A_map(), 3)):
        theirs, _ = fg.successors(q, kind, 1 if kind != 3 else 0, 1.0, 0, p, PH)
        mismatch = int(np.sum(mine != theirs))
        s_mine = graph_stats(mine); s_theirs, _ = fg.summary(theirs, N)
        rec[name] = {"successor_mismatches": mismatch, "independent": s_mine,
                     "agree_indegree": s_mine["indegree_hist"] == s_theirs["indegree_hist"], "agree_cycles": s_mine["cycles"] == s_theirs["cycles"],
                     "agree_cycle_type": s_mine["cycle_type"] == s_theirs["cycle_type"] if isinstance(s_theirs["cycle_type"], dict) and "distinct_lengths" not in s_theirs["cycle_type"] else None,
                     "agree_mean_rho": abs(s_mine["mean_rho"] - s_theirs["mean_rho"]) < 1e-9, "agree_max_tail": s_mine["max_tail"] == s_theirs["max_tail"]}
    out["profile_%d" % prof] = rec
    print(prof, {k: {kk: vv for kk, vv in v.items() if kk != "independent"} for k, v in rec.items()})
Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
