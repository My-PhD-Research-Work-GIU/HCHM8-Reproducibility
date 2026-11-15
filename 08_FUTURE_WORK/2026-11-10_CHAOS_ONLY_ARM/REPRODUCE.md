# Reproducing the chaos-only arms of the ablation (paper Section 9.3, Table 10)

Plan: `PLAN.md` (pushed before the registered run). Code: `CODE/chaos_only_arm.py`. Stored result: `RESULTS_chaos_only.json`,
the raw statistics of every table `DATA/chaos_only_metrics.csv` (one row per arm, key and round, with the SHA-256 of the
table) and the console output `DATA/run_log.txt`. The tables of the paper are written from the stored result by
`05_CODE/make_rev9_tables.py` (macros `\Co...` in `tables/rev9_numbers.tex`, table body `tables/tab_chaos_only.tex`).

## Environment of the registered run

Windows 11 (10.0.26300), CPython 3.12.13, NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.1, and Matplotlib 3.11.2 and pandas 3.1.1,
which the paper's module imports but this test does not use. (pandas was missing from the first version of these
instructions; the Linux rerun of the workflow stopped at that import and the line was added on 2026-11-10.) The script copies `05_CODE/run_validated_experiments.py` and `05_CODE/platform_rerun_and_extended_ablation.py` to a
temporary folder and takes `lattice_serialize` verbatim from Part D of the second file.

## Steps

From the folder `05_PAPERS_AND_DATASETS/05_PAPERS_AND_DATASETS` of a clone of the repository, with the environment of
`../2026-11-10_GRID_CONJUGACY/REPRODUCE.md`:

```
.venv-rerun/Scripts/python -X utf8 08_FUTURE_WORK/2026-11-10_CHAOS_ONLY_ARM/CODE/chaos_only_arm.py <out_dir>
```

The run took about six minutes (4.5 for the tables, the rest for the tests).

## Parameters

Keys `key_for_index(40000 + k)` for k = 0 to 255 (the keys of the paper's extended ablation), rounds 0 to 15, candidate 0:
4,096 tables per arm. Arms S (SHAKE-only), H (HCHM+SHAKE), C (first 258 bytes of the Float64 orbit, 64 burn-in steps, 32 bits
per coordinate, no hash) and CL (first 258 bytes of the 16-bit lattice orbit, 64 burn-in steps, no hash). Eight statistics;
Kolmogorov-Smirnov for the four continuous ones, Mann-Whitney U for the four counts; Holm over the eight tests of each
comparison; Hedges' g with a two-sample normal-approximation interval; Cliff's delta.

## Output and comparison

`<out_dir>/RESULTS_chaos_only.json` and `<out_dir>/chaos_only_metrics.csv`; compare with
`03_REPRODUCIBILITY_AND_QA/rev9_tests_rerun_2026-11-10/CODE/compare_rerun.py`. In the same environment every table and every
statistic is identical. In another environment the SHAKE-only tables must be identical; the lattice tables were identical
between the environments tested before; the H and C tables pass through a Float64 orbit and differ for some keys (paper
Section 9.4), so the numbers of the comparisons that involve them can change.
