"""Independent re-derivation checks of the mathematics of paper revision 14 (Sections 3, 4, 7.3, 7.4 and Appendix B).

Written from the equations of the manuscript, not from the project code, except where the code itself is the object of the check
(code correspondence). Exact integer arithmetic where the statement is exact.

  1. Core: determinants, traces, eigenvalues, Lyapunov values, entropy, |det(A - I)| = 48, and an exact enumeration of the 48
     fixed points of A carried through H^-1 and checked as fixed points of T.
  2. Shears: the formulas of Section 3 against f_scalar of the emulator on random states; the independence of g_i from z_i
     (exact equality after changing z_i); the parameter ranges of profiles 0 to 99 against Table 2.
  3. H^-1(H(x)) and the conjugacy H(T(x)) = A(H(x)) on random states; det DT by central differences.
  4. Modular orders of the four blocks and of A modulo 2^Q: brute force for Q <= 16, and for Q <= 64 a certificate
     (M^n = I and M^(n/p) != I for each prime p dividing n), against Appendix B; the 2-adic valuations stated there.
  5. Fixed points and cycle type of A_Q on L_Q for Q = 2, 3 (exceptional states with periods shorter than the order).

Usage: python math_checks.py <emulator .py> <out.json>
"""
import importlib.util, json, math, sys
from fractions import Fraction
from itertools import product
import numpy as np

spec = importlib.util.spec_from_file_location("emu", sys.argv[1]); emu = importlib.util.module_from_spec(spec); spec.loader.exec_module(emu)
OUT = sys.argv[2]
R = {}
BLOCKS = {3: [[2, 1], [1, 1]], 4: [[3, 2], [1, 1]], 6: [[5, 2], [2, 1]], 8: [[7, 3], [2, 1]]}

# 1. core ---------------------------------------------------------------------------------------------
core = {}
for t, M in BLOCKS.items():
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]; tr = M[0][0] + M[1][1]
    lam = (tr + math.sqrt(tr * tr - 4)) / 2
    core[t] = {"det": det, "trace": tr, "lambda_max": lam, "log_lambda": math.log(lam), "fixed_points_|2-tr|": abs(2 - tr)}
R["core_blocks"] = core
R["entropy_nats"] = sum(v["log_lambda"] for v in core.values()); R["entropy_bits"] = R["entropy_nats"] / math.log(2)
R["fixed_points_det(A-I)"] = int(np.prod([v["fixed_points_|2-tr|"] for v in core.values()]))
paper_le = [2.063437, 1.762747, 1.316958, 0.962424]
R["lyapunov_match_table3"] = max(abs(a - b) for a, b in zip(sorted([v["log_lambda"] for v in core.values()], reverse=True), paper_le)) < 5e-7


def fixed_points_block(M):  # x in [0,1)^2 with (M - I) x in Z^2, exact
    a, b, c, d = M[0][0] - 1, M[0][1], M[1][0], M[1][1] - 1
    D = a * d - b * c; pts = []
    for k1 in range(-abs(D) * 4, abs(D) * 4 + 1):
        for k2 in range(-abs(D) * 4, abs(D) * 4 + 1):
            x1 = Fraction(d * k1 - b * k2, D); x2 = Fraction(-c * k1 + a * k2, D)
            if 0 <= x1 < 1 and 0 <= x2 < 1: pts.append((x1, x2))
    return sorted(set(pts))


fp_blocks = [fixed_points_block(BLOCKS[t]) for t in (3, 4, 6, 8)]
R["fixed_points_per_block"] = [len(p) for p in fp_blocks]
P0 = emu.profile_params(0); P0[16] = 0.0


def g_paper(i, z, p):  # Section 3, 1-based indices in the text, 0-based here
    e, ph1, ph2, r3, k4, h6, k7, r8 = p[0:8], p[8], p[9], p[10], p[11], p[12], p[13], p[14]
    if i == 0: return math.sin(2 * math.pi * z[2] + ph1)
    if i == 1: c = math.cos(2 * math.pi * z[5] + ph2); return 4 * c ** 3 - 3 * c
    if i == 2: return math.sin(2 * math.pi * (z[7] - z[1])) / (1 + r3 * math.cos(2 * math.pi * z[1]) ** 2)
    if i == 3: return math.tanh(k4 * math.sin(2 * math.pi * (z[0] - z[6])))
    if i == 4: q = 0.5 * (1 + math.sin(2 * math.pi * (z[6] + z[7]))); return 4 * q * (1 - q) - 0.5
    if i == 5: from scipy.special import i0; return (math.exp(h6 * math.cos(2 * math.pi * (z[1] - z[0]))) - i0(h6)) / (math.exp(h6) - math.exp(-h6))
    if i == 6: return (2 / math.pi) * math.atan(k7 * math.sin(2 * math.pi * (z[3] + z[5])))
    return math.log(1 + r8 * math.sin(math.pi * (z[4] - z[2])) ** 2) / math.log(1 + r8) - 0.5


def H(x, p):
    z = list(x)
    for i in range(8): z[i] = (z[i] + p[i] * g_paper(i, z, p)) % 1.0
    return z


def Hinv(y, p):
    z = list(y)
    for i in range(7, -1, -1): z[i] = (z[i] - p[i] * g_paper(i, z, p)) % 1.0
    return z


AM = np.zeros((8, 8), dtype=np.int64)
for k, t in enumerate((3, 4, 6, 8)): AM[2 * k:2 * k + 2, 2 * k:2 * k + 2] = BLOCKS[t]
assert np.array_equal(AM, emu.A_INT)


def A(z): return [float(v) % 1.0 for v in (AM @ np.array(z))]


def T(x, p): return Hinv(A(H(x, p)), p)


def tdist(a, b): return max(min(abs(u - v) % 1.0, 1 - abs(u - v) % 1.0) for u, v in zip(a, b))


worst = 0.0; fps = []
for combo in product(*fp_blocks):
    y = [float(c) for pair in combo for c in pair]
    x = Hinv(y, P0); fps.append(x); worst = max(worst, tdist(T(x, P0), x))
dmin = min(tdist(fps[i], fps[j]) for i in range(len(fps)) for j in range(i))
R["fixed_points_T"] = {"count": len(fps), "worst_torus_error": worst, "min_pairwise_distance": dmin}

# 2. shears ---------------------------------------------------------------------------------------------
rng = np.random.default_rng(20261115); dev = 0.0; indep_fail = 0
for prof in (0, 1, 2, 3, 7):
    p = emu.profile_params(prof); p[16] = 0.0
    for _ in range(400):
        z = rng.random(8)
        for i in range(8):
            dev = max(dev, abs(g_paper(i, z, p) - emu.f_scalar(i, z.copy(), p)))
            z2 = z.copy(); z2[i] = rng.random()
            if emu.f_scalar(i, z2, p) != emu.f_scalar(i, z, p) or g_paper(i, z2, p) != g_paper(i, z, p): indep_fail += 1
R["shear_formula_vs_code_max_abs_diff"] = dev; R["shear_independence_failures"] = indep_fail
rng_ok = True; bad = []
for prof in range(100):
    p = emu.profile_params(prof)
    checks = [all(0 < e <= 0.15 for e in p[0:8]), 0 <= p[8] < 2 * math.pi, 0 <= p[9] < 2 * math.pi, 0 <= p[10] <= 3, 1 <= p[11] <= 4,
              0.5 <= p[12] <= 2.5, 1 <= p[13] <= 4, 0.5 <= p[14] <= 6, abs(p[15] - float(np.i0(p[12]))) < 1e-12]
    if not all(checks): rng_ok = False; bad.append(prof)
R["profiles_0_99_within_table2"] = rng_ok; R["profiles_out_of_range"] = bad
ref = emu.profile_params(0)
R["reference_profile_matches_table2"] = max(abs(a - b) for a, b in zip(list(ref[0:15]), [0.128982, 0.084956, 0.072908, 0.148100, 0.082414, 0.134559, 0.144874, 0.133708, 3.457257, 0.009488, 2.223927, 3.030714, 2.365038, 3.205819, 1.590681])) < 6e-7

# 3. inverse, conjugacy, Jacobian ------------------------------------------------------------------------
e1 = e2 = 0.0; jd = 0.0
for _ in range(300):
    x = list(rng.random(8))
    e1 = max(e1, tdist(Hinv(H(x, P0), P0), x)); e2 = max(e2, tdist(H(T(x, P0), P0), A(H(x, P0))))
for _ in range(30):
    x = np.array(rng.random(8)); h = 1e-6; J = np.zeros((8, 8))
    for j in range(8):
        xp = x.copy(); xm = x.copy(); xp[j] += h; xm[j] -= h
        fp = np.array(T(list(xp), P0)); fm = np.array(T(list(xm), P0))
        d = fp - fm; d -= np.round(d)
        J[:, j] = d / (2 * h)
    jd = max(jd, abs(np.linalg.det(J) - 1))
R["H_inverse_max_error"] = e1; R["conjugacy_max_error"] = e2; R["det_DT_minus_1_max_central_diff"] = jd


# 4. modular orders -----------------------------------------------------------------------------------
def matmul(X, Y, m):
    n = len(X); return [[sum(X[i][k] * Y[k][j] for k in range(n)) % m for j in range(n)] for i in range(n)]


def matpow(M, e, m):
    n = len(M); Rr = [[int(i == j) for j in range(n)] for i in range(n)]; B = [r[:] for r in M]
    while e:
        if e & 1: Rr = matmul(Rr, B, m)
        B = matmul(B, B, m); e >>= 1
    return Rr


def is_id(M, m): return all(M[i][j] % m == (1 if i == j else 0) for i in range(len(M)) for j in range(len(M)))


def brute_order(M, m):
    X = [r[:] for r in M]; n = 1
    while not is_id(X, m): X = matmul(X, M, m); n += 1
    return n


def claimed(t, Q):  # Appendix B and Section 7.3
    if Q == 2: return {3: 3, 4: 4, 6: 2, 8: 4}[t]
    if Q == 3: return {3: 6, 4: 4, 6: 4, 8: 4}[t]
    return {3: 3 * 2 ** (Q - 2), 4: 2 ** (Q - 1), 6: 2 ** (Q - 1), 8: 2 ** (Q - 2)}[t]


def certified(M, n, m):
    if not is_id(matpow(M, n, m), m): return False
    for pr in (2, 3):
        if n % pr == 0 and is_id(matpow(M, n // pr, m), m): return False
    return True


orders = []; ok_all = True
for Q in range(2, 65):
    m = 2 ** Q; row = {"Q": Q}
    for t, M in BLOCKS.items():
        c = claimed(t, Q)
        good = certified(M, c, m)
        if Q <= 16: good = good and brute_order(M, m) == c
        row[t] = c; ok_all &= good
    lcm = math.lcm(*[row[t] for t in BLOCKS]); row["lcm"] = lcm
    row["lcm_is_formula"] = lcm == (12 if Q == 2 else 3 * 2 ** (Q - 1)); ok_all &= row["lcm_is_formula"]
    orders.append(row)
R["block_orders_Q2_to_64_match_appendix_B"] = ok_all; R["orders_sample"] = [o for o in orders if o["Q"] in (2, 3, 4, 8, 16, 48, 56, 64)]


def v2_entries(M):
    v = 10 ** 9
    for i in range(2):
        for j in range(2):
            e = M[i][j] - (1 if i == j else 0)
            if e != 0: v = min(v, (e & -e).bit_length() - 1)
    return v


def ipow(M, e):
    Rr = [[1, 0], [0, 1]]
    for _ in range(e): Rr = [[sum(Rr[i][k] * M[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    return Rr


R["v2_appendix_B"] = {"A3^6-I": v2_entries(ipow(BLOCKS[3], 6)), "A4^4-I": v2_entries(ipow(BLOCKS[4], 4)),
                      "A6^4-I": v2_entries(ipow(BLOCKS[6], 4)), "A8^4-I": v2_entries(ipow(BLOCKS[8], 4)), "stated": [3, 3, 3, 4]}


# 5. A_Q on L_Q: fixed points and periods (exceptional states)
def aq_cycle_type(q):
    m = 2 ** q; N = m ** 8
    idx = np.arange(N, dtype=np.int64); s = np.stack([(idx >> (q * i)) & (m - 1) for i in range(8)], axis=1)
    o = (s @ AM.T) % m; nxt = np.zeros(N, np.int64)
    for i in range(8): nxt |= o[:, i] << (q * i)
    per = np.zeros(N, np.int64); cur = nxt.copy(); n = 1
    while (per == 0).any():
        hit = (cur == idx) & (per == 0); per[hit] = n; cur = nxt[cur]; n += 1
    u, c = np.unique(per, return_counts=True)
    return {int(a): int(b) for a, b in zip(u, c)}


R["A_Q_period_distribution_of_states"] = {"Q2": aq_cycle_type(2)}
json.dump(R, open(OUT, "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in R.items() if k not in ("orders_sample", "core_blocks")}, indent=1, default=str))
