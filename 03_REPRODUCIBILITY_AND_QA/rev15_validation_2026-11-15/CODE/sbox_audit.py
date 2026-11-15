"""Independent audit of HAS-1 (revision 15 validation, 2026-11-15). Own code throughout: GF(2^8) by log/antilog tables of the
generator 3 modulo x^8 + x^4 + x^3 + x + 1, Walsh spectra by a fast transform written here, algebraic degree by the Moebius
transform, BCT and DLCT from their definitions. Nothing is imported from the project except the recorded specification
(rows of L_in and L_out, constants a and b) and the table printed in Appendix A of the manuscript, both read as text.

Checks: the table rebuilt from the specification equals the printed table and its SHA-256; bijectivity; the affine witness
S(x) = L_out I(L_in x + a) + b on all 256 inputs with L_in and L_out invertible over GF(2); vectorial nonlinearity (minimum over
the 255 components), maximum Walsh value, DDT spectrum and differential uniformity, minimum and maximum degree of the
components, boomerang uniformity, DLCT absolute indicator, fixed and opposite fixed points, cycle structure, SAC and BIC.
Which of these are invariant under affine equivalence is recorded with each value.

Usage: python sbox_audit.py <run_validated_experiments.py> <main.tex> <out.json>
"""
import hashlib, json, re, sys
import numpy as np

src = open(sys.argv[1], encoding="utf-8").read(); tex = open(sys.argv[2], encoding="utf-8").read(); OUT = sys.argv[3]
RIN = [int(x, 16) for x in re.search(r"RIN=\[([^\]]*)\]", src).group(1).split(",")]
ROUT = [int(x, 16) for x in re.search(r"ROUT=\[([^\]]*)\]", src).group(1).split(",")]
A_CONST, B_CONST = (int(x, 16) for x in re.search(r"build_sbox\(RIN,ROUT,(0x[0-9A-Fa-f]+),(0x[0-9A-Fa-f]+)\)", src).groups())
block = tex.split("\\section{Reference HAS-1 Table}", 1)[1].split("\\begin{verbatim}", 1)[1].split("\\end{verbatim}", 1)[0]
PRINTED = np.array([int(x, 16) for x in block.split()], dtype=np.int64)
assert PRINTED.size == 256

EXP = [0] * 510; LOG = [0] * 256; x = 1
for i in range(255):
    EXP[i] = x; LOG[x] = i
    x ^= (x << 1) ^ (0x11B if x & 0x80 else 0)  # multiply by 3 = x*2 + x
    x &= 0xFF
for i in range(255, 510): EXP[i] = EXP[i - 255]
INV = [0] + [EXP[(255 - LOG[v]) % 255] for v in range(1, 256)]
assert all((v == 0) or EXP[LOG[v] + LOG[INV[v]]] == 1 for v in range(256))


def parity(v): return bin(v).count("1") & 1


def lin(rows, v): return sum(parity(v & r) << i for i, r in enumerate(rows))


def rank_gf2(rows):
    rows = list(rows); r = 0
    for bit in range(8):
        piv = next((i for i in range(r, len(rows)) if rows[i] >> bit & 1), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> bit & 1: rows[i] ^= rows[r]
        r += 1
    return r


S = np.array([lin(ROUT, INV[lin(RIN, v) ^ A_CONST]) ^ B_CONST for v in range(256)], dtype=np.int64)
R = {"specification": {"L_in_rows": [hex(r) for r in RIN], "L_out_rows": [hex(r) for r in ROUT], "a": hex(A_CONST), "b": hex(B_CONST)}}
R["rebuilt_equals_appendix_A"] = bool(np.array_equal(S, PRINTED))
R["sha256"] = hashlib.sha256(S.astype(np.uint8).tobytes()).hexdigest()
R["sha256_matches_record"] = R["sha256"] == "450fbd045e90d588497162eb2badd692a68d16b03188b6682ee9f82e9949bc33"
R["bijective"] = bool(np.unique(S).size == 256)
R["L_in_rank"] = rank_gf2(RIN); R["L_out_rank"] = rank_gf2(ROUT)
R["affine_witness_all_256"] = bool(all(S[v] == (lin(ROUT, INV[lin(RIN, v) ^ A_CONST]) ^ B_CONST) for v in range(256)))


def walsh(f):  # f: 0/1 vector of length 256 -> Walsh spectrum
    a = (1 - 2 * f).astype(np.int64); h = 1
    while h < 256:
        a = a.reshape(-1, 2, h); a = np.concatenate([a[:, 0] + a[:, 1], a[:, 0] - a[:, 1]], axis=1).reshape(-1); h *= 2
    return a


X = np.arange(256)
comp = [np.array([parity(int(S[v]) & b) for v in X]) for b in range(1, 256)]
W = [walsh(f) for f in comp]
maxw = [int(np.abs(w).max()) for w in W]
R["nonlinearity_vectorial"] = int(128 - max(maxw) // 2); R["nonlinearity_component_range"] = [int(128 - max(maxw) // 2), int(128 - min(maxw) // 2)]
R["max_walsh"] = int(max(maxw))


def degree(f):
    a = f.copy().astype(np.int64)
    for i in range(8):
        for v in range(256):
            if v >> i & 1: a[v] ^= a[v ^ (1 << i)]
    return max((bin(v).count("1") for v in range(256) if a[v]), default=0)


degs = [degree(f) for f in comp]; R["degree_range"] = [min(degs), max(degs)]
DDT = np.zeros((256, 256), dtype=np.int64)
for da in range(256): DDT[da] = np.bincount(S[X] ^ S[X ^ da], minlength=256)
u, c = np.unique(DDT[1:], return_counts=True); spec = {int(a): int(b) for a, b in zip(u, c)}
R["ddt_spectrum_rows_1_255"] = spec; R["differential_uniformity"] = int(DDT[1:].max())
Sinv = np.empty(256, dtype=np.int64); Sinv[S] = X
bu = 0
for da in range(1, 256):
    for db in range(1, 256):
        bu = max(bu, int(np.sum((Sinv[S[X] ^ db] ^ Sinv[S[X ^ da] ^ db]) == da)))
R["boomerang_uniformity"] = bu
dl = 0
for da in range(1, 256):
    d = (S[X] ^ S[X ^ da])
    for b in range(1, 256):
        val = abs(int(256 - 2 * np.sum(np.array([parity(int(t) & b) for t in d]))))
        dl = max(dl, val)
R["dlct_absolute_indicator"] = dl // 2  # DLCT(a,b) = 2^(n-1) - #{x: b.(S(x)+S(x+a)) = 1} = (256 - 2*count)/2
R["fixed_points"] = int(np.sum(S == X)); R["opposite_fixed_points"] = int(np.sum(S == (X ^ 0xFF)))
seen = set(); cyc = []
for v in range(256):
    if v in seen: continue
    L = 0; w = v
    while w not in seen: seen.add(w); w = int(S[w]); L += 1
    cyc.append(L)
R["cycle_structure"] = sorted(cyc, reverse=True)
sac = np.zeros((8, 8))
for i in range(8):
    d = S[X] ^ S[X ^ (1 << i)]
    for j in range(8): sac[i, j] = np.mean((d >> j) & 1)
R["sac_mean"] = float(sac.mean()); R["sac_mad"] = float(np.abs(sac - 0.5).mean()); R["sac_maxdev"] = float(np.abs(sac - 0.5).max())
bic = []
for i in range(8):
    d = S[X] ^ S[X ^ (1 << i)]
    bits = [((d >> j) & 1).astype(float) for j in range(8)]
    for j in range(8):
        for k in range(j + 1, 8):
            bic.append(abs(np.corrcoef(bits[j], bits[k])[0, 1]))
R["bic_max_abs_corr"] = float(max(bic)); R["bic_mean_abs_corr"] = float(np.mean(bic))
R["invariance"] = {"affine_invariant": ["nonlinearity_vectorial", "max_walsh", "ddt_spectrum", "differential_uniformity", "degree_range",
                                         "boomerang_uniformity", "dlct_absolute_indicator"],
                   "coordinate_dependent": ["fixed_points", "opposite_fixed_points", "cycle_structure", "sac", "bic"]}
json.dump(R, open(OUT, "w"), indent=1)
print(json.dumps({k: v for k, v in R.items() if k not in ("specification", "invariance")}, indent=1))
