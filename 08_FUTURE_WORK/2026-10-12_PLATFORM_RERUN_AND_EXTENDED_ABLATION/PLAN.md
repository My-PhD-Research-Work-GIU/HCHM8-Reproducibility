# Platform rerun and extended ablation — plan

Date: 2026-10-12. One script, `05_CODE/platform_rerun_and_extended_ablation.py <tag> <dir>`, run in each environment.

## Questions

1. Do the published floating-point results reproduce in a different computing environment?
2. Does the 8,192-table corpus reproduce, table by table?
3. Does HCHM-8 change any statistic of the generated tables that is not fixed by affine equivalence? (The published ablation tested SAC only, on 512 tables per arm, and did not measure BIC.)
4. If the corpus does not reproduce, can the dependence be removed?

## Fixed before running

- Part A reruns `run_map_tests()` unchanged and compares with `04_DATASETS/map_validated_summary.json`.
- Part B regenerates the corpus with the unchanged generator and compares each SHA-256 with `04_DATASETS/dynamic_family_direct_8192.csv`.
- Part C: eight statistics for all 4,096 tables per arm — SAC mean absolute deviation, SAC maximum deviation, BIC maximum and mean absolute correlation, fixed points, opposite fixed points, number of cycles, longest cycle. Kolmogorov–Smirnov for the first four, Mann–Whitney U for the counts, Holm correction over the eight.
- Part D (added when Part B showed a mismatch): serialize lattice states (16 bits per coordinate) in place of Float64 states; test under the last-bit perturbation; compare tables with the SHAKE-only arm by the same eight tests.

## Added after the first run

Two statistics take very few distinct values in this family (SAC maximum deviation: 2; BIC maximum correlation: 15). A Fisher exact test and a Mann–Whitney test were added for these as a supplement and are labelled post hoc.
