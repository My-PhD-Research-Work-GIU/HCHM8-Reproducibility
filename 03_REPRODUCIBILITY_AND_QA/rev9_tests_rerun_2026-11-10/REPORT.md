# Report: rerun and independent check of the two tests of paper revision 9 (2026-11-10)

Plan: `PLAN.md`, pushed in commit c436e59 before any of the runs below (workflow run 42764875635 created on GitHub at
2026-11-10T05:52:58Z). Tests checked: the grid-exact conjugate (revision 9 Section 7.5, revision 10 Section 7.4, Table 6;
`08_FUTURE_WORK/2026-11-10_GRID_CONJUGACY`) and the chaos-only arms of the ablation (Section 9.3, Table 10;
`08_FUTURE_WORK/2026-11-10_CHAOS_ONLY_ARM`).

**Verdict: both tests reproduce. No stored value was found wrong.** One defect of the reproduction instructions was found and
corrected; two defects of the checking program were corrected after its first run. Details below.

## 1. Reproduction instructions, new copy (Windows 11)

Fresh clone of commit 4b9507e into a new folder; a new environment created exactly as `REPRODUCE.md` says
(`uv venv --python 3.12`, then NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.1, Matplotlib 3.11.2, pandas 3.1.1; CPython 3.12.13); both
commands run as written. Outputs: `DATA/windows11_fresh_clone/`; comparison: `DATA/COMPARISON_windows11-fresh-clone.json`.

| Test | Result |
|---|---|
| Grid-exact conjugate | statistics of F_Q identical in 20 of 20 cases (Q = 2, 3; ten profiles); cycle type of G_Q identical in 20 of 20; 790 of 790 orbit periods identical |
| Chaos-only arms | 4,096 of 4,096 tables identical in each of the four arms; every number of Table 10 identical; criteria of the plan met |

## 2. Linux (GitHub-hosted runner, glibc 2.43, CPython 3.12.16, same pinned packages)

Workflow job `rev9-tests-rerun`. Runs: 42764875635 (commit c436e59) completed the grid test and stopped at `import pandas` in the
chaos-only test, because the first version of the instructions did not list pandas, which `05_CODE/run_validated_experiments.py`
imports; the instructions were corrected in commit 4b9507e without changing any command or code. Run 42765909602 (4b9507e)
completed both tests and the comparison and then failed in the independent check when writing its JSON (item 3). Run
42767491873 (commit 41a3c8b) completed every step; its artifact is `DATA/linux_ci/`, with `DATA/COMPARISON_linux-ci.json`.

| Test | Result |
|---|---|
| Grid-exact conjugate | statistics of F_Q identical in 20 of 20 cases; 790 of 790 periods identical; every check of Proposition 2 held |
| Chaos-only arms | tables identical to the stored ones: SHAKE-only 4,096, lattice orbit 4,096, HCHM+SHAKE 2,774 (67.7%), Float64 orbit without hash 2,882 (70.4%); the paper's derivation reproduced by the script for every key and round |

The comparisons on Linux, stored value then Linux value:

| Comparison | Smallest p | Smallest Holm p | Largest abs(g) | Bound of g |
|---|---|---|---|---|
| HCHM+SHAKE vs SHAKE-only | 0.4375 / 0.0714 | 1.00 / 0.571 | 0.0176 / 0.0423 | 0.0609 / 0.0856 |
| Float64 orbit vs SHAKE-only | 0.2276 / 0.1579 | 1.00 / 1.00 | 0.0281 / 0.0334 | 0.0714 / 0.0767 |
| lattice orbit vs SHAKE-only | 0.0162 / 0.0162 | 0.129 / 0.129 | 0.0534 / 0.0534 | 0.0967 / 0.0967 |
| Float64 orbit vs HCHM+SHAKE | 0.2447 / 0.1355 | 1.00 / 1.00 | 0.0280 / 0.0238 | 0.0713 / 0.0671 |

The rows that involve the Float64 orbit change, because about 30% of those tables differ between the two systems (the
portability defect of the paper's Section 9.4); the lattice row is unchanged, because its tables are identical. On Linux, as on
Windows, no comparison is significant after Holm correction and every interval of g lies within +/-0.10, so the registered
criteria are met in both environments.

## 3. Independent code

`CODE/independent_check.py` evaluates F_Q and G_Q with NumPy arrays from the equations of the paper (Bessel value from
`scipy.special.i0`, blocks of A written out), analyses functional graphs by following paths (the original removes states of
in-degree zero), and recomputes the chaos-only test from its stored raw file with separately written statistics.

| Part | Windows 11 (Q = 2 and 3) | Linux (Q = 2) |
|---|---|---|
| A: every state, ten profiles | all 20 cases (2^16 and 2^24 states, ten profiles): every stored statistic of Table 6 recomputed; G_Q a permutation with the cycle type of A_Q; F_Q and G_Q equal to the arrays of the original program on every state (0 states differ) | the same at Q = 2 |
| B: 790 orbit periods, Q = 4 to 7 | all equal to the stored periods | all equal |
| C: tables and Table 10 | all 16,384 tables regenerated: SHA-256 and the eight statistics equal to the stored file (largest difference 6e-17); every number of Table 10 recomputed equal | SHA-256 and statistics equal for 13,848 of 16,384 regenerated tables (the 2,536 others are Float64-orbit tables that differ on Linux, as in item 2); every number of Table 10 recomputed from the stored file equal |

Two defects of the checking program were corrected after its first complete run on Linux (commit 41a3c8b): it could not
serialize NumPy booleans to JSON, and it computed the 95% intervals of g with the exact normal quantile 1.959964, which moved
their ends by less than 1e-5 from those of the test (which uses 1.96) and was reported as a difference. Neither change alters
what is compared.

Reported for the first time (Part C, from the stored raw file): the correlation of each statistic within matched pairs is at
most 0.049 in absolute value (largest for HCHM+SHAKE against SHAKE-only), so the two-sample intervals of Table 10 are adequate;
widened by a Bonferroni factor of eight (two-sided 0.05 / 8), the intervals of g stay within +/-0.078 (HCHM+SHAKE vs SHAKE-only),
0.089, 0.114 (lattice orbit vs SHAKE-only) and 0.088.

## 4. Parameters

The scripts, stored results and the text agree: profiles 0 to 9 (`profile_params`), starting states 32, 32, 13 and 2 at Q = 4
to 7 (790 orbits), keys `key_for_index(40000 + k)` for k = 0 to 255, rounds 0 to 15, candidate 0, 4,096 tables per arm. Each
result file records the SHA-256 of the script that wrote it; the reruns used the same scripts (`scripts_identical` in both
comparison files).

## Not a finding

The number of cycles of F_Q at 2^24 states (mean 13.1 over the ten profiles, range 5 to 21) is above the exact expectation
for a random mapping of the same size (9.0). This was computed from the stored results after the run, is not part of either
plan and was not tested; revision 10 reports it as such.
