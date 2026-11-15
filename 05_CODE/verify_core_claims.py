"""Independent verification of the core analytical and S-box claims.

Written separately from run_validated_experiments.py and sharing no code with it. Uses only the
Python standard library and NumPy. Every check is exact or has an explicit tolerance.

Usage: python verify_core_claims.py [output.json]      exits non-zero if any check fails
"""
import hashlib, itertools, json, math, re, sys
from fractions import Fraction
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
checks = []
def check(name, ok, detail=""):
    checks.append({"name": name, "passed": bool(ok), "detail": str(detail)})

# ---------------------------------------------------------------- 1. linear core
BLOCKS = [((2, 1), (1, 1)), ((3, 2), (1, 1)), ((5, 2), (2, 1)), ((7, 3), (2, 1))]
mus = []
for (a, b), (c, d) in BLOCKS:
    tr, det = a + d, a * d - b * c
    check("block %s has determinant 1" % str(((a, b), (c, d))), det == 1, det)
    check("block %s is hyperbolic (trace > 2)" % str(((a, b), (c, d))), tr > 2, tr)
    mus.append((tr + math.sqrt(tr * tr - 4)) / 2)
closed = [(3 + math.sqrt(5)) / 2, 2 + math.sqrt(3), 3 + 2 * math.sqrt(2), 4 + math.sqrt(15)]
check("expanding eigenvalues equal the closed forms in Corollary 1", max(abs(x - y) for x, y in zip(mus, closed)) < 1e-12, mus)
spec = sorted((math.log(m) for m in mus), reverse=True)
check("positive Lyapunov exponents to 6 decimals", ["%.6f" % s for s in spec] == ["2.063437", "1.762747", "1.316958", "0.962424"], spec)
ent = sum(spec)
check("entropy 6.105566 nats = 8.808469 bits", "%.6f" % ent == "6.105566" and "%.6f" % (ent / math.log(2)) == "8.808469", (ent, ent / math.log(2)))
fixed_counts = [abs((a - 1) * (d - 1) - b * c) for (a, b), (c, d) in BLOCKS]
check("fixed-point counts per block are 1, 2, 4, 6 (product 48)", fixed_counts == [1, 2, 4, 6] and math.prod(fixed_counts) == 48, fixed_counts)

# exact order of A modulo 2^Q (the period of the linear core on the lattice with Q bits per coordinate)
def order_mod(block, m):
    (a, b), (c, d) = block; x = (a % m, b % m, c % m, d % m); k = 1
    while x != (1 % m, 0, 0, 1 % m):
        x = ((x[0] * a + x[1] * c) % m, (x[0] * b + x[1] * d) % m, (x[2] * a + x[3] * c) % m, (x[2] * b + x[3] * d) % m); k += 1
    return k
orders = {}
for q in range(2, 17):
    per = [order_mod(B, 1 << q) for B in BLOCKS]; orders[q] = math.lcm(*per)
check("order of A modulo 2^Q equals 3*2^(Q-1) for every Q from 3 to 16", all(orders[q] == 3 * 2 ** (q - 1) for q in range(3, 17)), orders)
check("order of A modulo 4 is 12", orders[2] == 12, orders[2])

# ---------------------------------------------------------------- 2. the map itself (reference profile)
def unit(label, n):
    raw = hashlib.shake_256(label).digest(8 * n)
    return [int.from_bytes(raw[8 * i:8 * i + 8], "big") / 2.0 ** 64 for i in range(n)]
r = unit(b"GIU-PhD-HCHM8-heterogeneous-conjugacy-v1", 32)
EPS = [0.07 + 0.08 * v for v in r[8:16]]
PHI1, PHI2 = 2 * math.pi * r[16], 2 * math.pi * r[17]
RHO3, K4, ETA6, K7, RHO8 = 0.75 + 2.25 * r[18], 1.25 + 2.75 * r[19], 0.75 + 1.75 * r[20], 1.25 + 2.75 * r[21], 1.0 + 5.0 * r[22]
I0 = float(np.i0(ETA6))
import current_versions
table1 = (current_versions.PAPER / "main.tex").read_text(encoding="utf-8")
printed = ", ".join("%.6f" % e for e in EPS)
check("Table 1 shear amplitudes reproduce from the SHAKE256 label", printed in table1, printed)
check("Table 1 remaining parameters reproduce", all(("%.6f" % v) in table1 for v in (PHI1, PHI2, RHO3, K4, ETA6, K7, RHO8)), (PHI1, PHI2, RHO3, K4, ETA6, K7, RHO8))

def f(i, z):
    tp = 2 * math.pi
    if i == 0: return math.sin(tp * z[2] + PHI1)
    if i == 1:
        c = math.cos(tp * z[5] + PHI2); return 4 * c ** 3 - 3 * c
    if i == 2: return math.sin(tp * (z[7] - z[1])) / (1 + RHO3 * math.cos(tp * z[1]) ** 2)
    if i == 3: return math.tanh(K4 * math.sin(tp * (z[0] - z[6])))
    if i == 4:
        q = 0.5 * (1 + math.sin(tp * (z[6] + z[7]))); return 4 * q * (1 - q) - 0.5
    if i == 5: return (math.exp(ETA6 * math.cos(tp * (z[1] - z[0]))) - I0) / (math.exp(ETA6) - math.exp(-ETA6))
    if i == 6: return (2 / math.pi) * math.atan(K7 * math.sin(tp * (z[3] + z[5])))
    return math.log(1 + RHO8 * math.sin(math.pi * (z[4] - z[2])) ** 2) / math.log(1 + RHO8) - 0.5
def H(x):
    z = list(x)
    for i in range(8): z[i] = (z[i] + EPS[i] * f(i, z)) % 1.0
    return z
def Hinv(y):
    z = list(y)
    for i in range(7, -1, -1): z[i] = (z[i] - EPS[i] * f(i, z)) % 1.0
    return z
def A(z):
    out = []
    for k, ((a, b), (c, d)) in enumerate(BLOCKS):
        u, v = z[2 * k], z[2 * k + 1]; out += [(a * u + b * v) % 1.0, (c * u + d * v) % 1.0]
    return out
def T(x): return Hinv(A(H(x)))
def tdist(u, v): return max(min(abs(p - q), 1 - abs(p - q)) for p, q in zip(u, v))

rng = np.random.default_rng(20261012); pts = rng.random((2000, 8)).tolist()
check("H^-1(H(x)) = x on 2,000 random states (tolerance 1e-14)", max(tdist(Hinv(H(p)), p) for p in pts) < 1e-14, max(tdist(Hinv(H(p)), p) for p in pts))
check("H(T(x)) = A(H(x)) on 2,000 random states (tolerance 1e-13)", max(tdist(H(T(p)), A(H(p))) for p in pts) < 1e-13, max(tdist(H(T(p)), A(H(p))) for p in pts))
# periodicity of every shear function in every argument (well-definedness on the torus)
worst = 0.0
for p in pts[:200]:
    for i in range(8):
        base = f(i, p)
        for j in range(8):
            q = list(p); q[j] += 1.0; worst = max(worst, abs(f(i, q) - base))
check("every shear function is 1-periodic in every coordinate (tolerance 1e-12)", worst < 1e-12, worst)
dep = max(abs(f(i, p) - f(i, [v if j != i else (v + 0.37) % 1 for j, v in enumerate(p)])) for p in pts[:200] for i in range(8))
check("shear i does not depend on coordinate i", dep == 0.0, dep)
# Jacobian determinant of H by central differences
def jac_det(fun, x, h=1e-6):
    J = np.zeros((8, 8))
    for j in range(8):
        a, b = list(x), list(x); a[j] += h; b[j] -= h
        fa, fb = fun(a), fun(b)
        J[:, j] = [((u - v + 0.5) % 1.0 - 0.5) / (2 * h) for u, v in zip(fa, fb)]
    return float(np.linalg.det(J))
dets = [jac_det(H, p) for p in pts[:50]]
check("det DH = 1 by finite differences at 50 states (tolerance 1e-6)", max(abs(d - 1) for d in dets) < 1e-6, max(abs(d - 1) for d in dets))
dets = [jac_det(T, p) for p in pts[:50]]
check("det DT = 1 by finite differences at 50 states (tolerance 1e-5)", max(abs(d - 1) for d in dets) < 1e-5, max(abs(d - 1) for d in dets))

# the 48 fixed points: enumerate those of A exactly, carry them through H^-1, test T(x) = x
def block_fixed(block):
    (a, b), (c, d) = block; det = (a - 1) * (d - 1) - b * c; out = set()
    for k1 in range(abs(det)):
        for k2 in range(abs(det)):
            y1 = Fraction((d - 1) * k1 - b * k2, det) % 1; y2 = Fraction(-c * k1 + (a - 1) * k2, det) % 1
            if ((a - 1) * y1 + b * y2).denominator == 1 and (c * y1 + (d - 1) * y2).denominator == 1: out.add((y1, y2))
    return sorted(out)
per_block = [block_fixed(B) for B in BLOCKS]
check("exact enumeration gives 1, 2, 4, 6 fixed points per block", [len(p) for p in per_block] == [1, 2, 4, 6], [len(p) for p in per_block])
fixed_A = [sum(c, ()) for c in itertools.product(*per_block)]
fixed_T = [Hinv([float(v) for v in y]) for y in fixed_A]
check("48 points x = H^-1(y) satisfy T(x) = x (tolerance 1e-13)", len(fixed_T) == 48 and max(tdist(T(x), x) for x in fixed_T) < 1e-13, max(tdist(T(x), x) for x in fixed_T))
sep = min(tdist(u, v) for u, v in itertools.combinations(fixed_T, 2))
check("the 48 fixed points are pairwise distinct", sep > 1e-3, sep)

# ---------------------------------------------------------------- 3. HAS-1 rebuilt from its formula
def gmul(a, b):
    out = 0
    for _ in range(8):
        if b & 1: out ^= a
        a <<= 1
        if a & 0x100: a ^= 0x11B
        b >>= 1
    return out
inv = [0] * 256
for a in range(1, 256):
    inv[a] = next(b for b in range(1, 256) if gmul(a, b) == 1)
par = [bin(v).count("1") & 1 for v in range(256)]
def lin(rows, x): return sum(par[rows[i] & x] << i for i in range(8))
LIN = [0xA1, 0x14, 0xD6, 0x71, 0xFB, 0x49, 0x21, 0xB6]; LOUT = [0x11, 0x56, 0x94, 0x43, 0x71, 0xEE, 0x19, 0x91]
S = [lin(LOUT, inv[lin(LIN, x) ^ 0xF0]) ^ 0x77 for x in range(256)]
rank = lambda rows: np.linalg.matrix_rank(np.array([[(r >> j) & 1 for j in range(8)] for r in rows]))  # over the reals; confirmed over GF(2) below
def gf2_rank(rows):
    rows = list(rows); rk = 0
    for bit in range(8):
        piv = next((i for i in range(rk, 8) if rows[i] >> bit & 1), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(8):
            if i != rk and rows[i] >> bit & 1: rows[i] ^= rows[rk]
        rk += 1
    return rk
check("L_in and L_out are invertible over GF(2)", gf2_rank(LIN) == 8 and gf2_rank(LOUT) == 8, (gf2_rank(LIN), gf2_rank(LOUT)))
sha = hashlib.sha256(bytes(S)).hexdigest()
check("HAS-1 SHA-256 checksum", sha == "450fbd045e90d588497162eb2badd692a68d16b03188b6682ee9f82e9949bc33", sha)
check("HAS-1 is a permutation; S(00)=09, S(01)=68", sorted(S) == list(range(256)) and S[0] == 0x09 and S[1] == 0x68, (S[0], S[1]))
app = re.search(r"\\begin\{verbatim\}(.*?)\\end\{verbatim\}", table1, re.S).group(1).split()
check("appendix table in main.tex equals the rebuilt table", [int(v, 16) for v in app] == S, len(app))
Sa = np.array(S); X = np.arange(256)
ddt = np.zeros((256, 256), int)
for a in range(256): ddt[a] = np.bincount(Sa ^ Sa[X ^ a], minlength=256)
check("differential uniformity 4; DDT has 32,895 zeros, 32,130 twos, 255 fours (rows a != 0)", ddt[1:].max() == 4 and [(ddt[1:] == v).sum() for v in (0, 2, 4)] == [32895, 32130, 255], [(int((ddt[1:] == v).sum())) for v in (0, 2, 4)])
stored = np.loadtxt(ROOT / "04_DATASETS" / "HAS1_DDT.csv", delimiter=",", dtype=int)
check("stored HAS1_DDT.csv equals the independently computed DDT", np.array_equal(stored, ddt))
parr = np.array(par)
walsh = np.array([[abs(int(np.sum(1 - 2 * (parr[a & X] ^ parr[b & Sa])))) for a in range(256)] for b in range(1, 256)])
check("vectorial nonlinearity 112 for all 255 components; max Walsh 32", (walsh.max(axis=1) == 32).all(), int(walsh.max()))
Si = np.zeros(256, int); Si[Sa] = X
bu = max(int(np.sum((Si[Sa ^ b] ^ Si[Sa[X ^ a] ^ b]) == a)) for a in range(1, 256) for b in range(1, 256))
check("boomerang uniformity 6", bu == 6, bu)
def degree(tt):
    c = list(tt)
    for i in range(8):
        for x in range(256):
            if x >> i & 1: c[x] ^= c[x ^ (1 << i)]
    return max(bin(x).count("1") for x in range(256) if c[x])
degs = {degree(parr[b & Sa]) for b in range(1, 256)}
check("algebraic degree 7 for all 255 components", degs == {7}, degs)
dl = max(abs(int(np.sum(1 - 2 * parr[(Sa ^ Sa[X ^ a]) & b]))) for a in range(1, 256) for b in range(1, 256))
check("DLCT absolute indicator 32", dl == 32, dl)
check("no fixed and no opposite fixed points", int((Sa == X).sum()) == 0 and int((Sa == (X ^ 255)).sum()) == 0)
seen, cyc = set(), []
for x in range(256):
    if x not in seen:
        n, y = 0, x
        while y not in seen: seen.add(y); y = S[y]; n += 1
        cyc.append(n)
check("permutation cycle lengths 116, 110, 30", sorted(cyc, reverse=True) == [116, 110, 30], sorted(cyc, reverse=True))
sacm = np.array([[np.mean(((Sa ^ Sa[X ^ (1 << i)]) >> j) & 1) for j in range(8)] for i in range(8)])
bvals = []
for i in range(8):
    dbits = np.array([[(v >> j) & 1 for j in range(8)] for v in (Sa ^ Sa[X ^ (1 << i)])], float); C = np.corrcoef(dbits, rowvar=False)
    bvals += [abs(C[p, q]) for p in range(8) for q in range(p + 1, 8)]
check("BIC maximum / mean absolute correlation", ("%.6f" % max(bvals), "%.6f" % (sum(bvals) / len(bvals))) == ("0.131696", "0.058254"), (max(bvals), sum(bvals) / len(bvals)))
check("SAC mean / mean absolute deviation / maximum deviation", ("%.6f" % sacm.mean(), "%.6f" % abs(sacm - .5).mean(), "%.6f" % abs(sacm - .5).max()) == ("0.497314", "0.021240", "0.046875"),
      (sacm.mean(), abs(sacm - .5).mean(), abs(sacm - .5).max()))

for c in checks: print("[%s] %s%s" % ("PASS" if c["passed"] else "FAIL", c["name"], "" if c["passed"] else "  -> " + c["detail"][:300]))
print("%d/%d checks passed" % (sum(c["passed"] for c in checks), len(checks)))
if len(sys.argv) > 1:
    json.dump({"checks": checks, "order_of_A_mod_2^Q": {str(k): v for k, v in orders.items()}}, open(sys.argv[1], "w"), indent=1)
sys.exit(0 if all(c["passed"] for c in checks) else 1)
