# Cost of three S-box derivations (2026-11-10)

Plan: `PLAN.md`, pushed in commit c436e59 before the registered run. Workflow run 42764875635 of that commit was created on
GitHub at 2026-11-10T05:52:58Z; the registered run started at 05:53:36Z and finished at 06:06:35Z (`RESULTS_derivation_cost.json`).

Command, as in the plan: `python -X utf8 CODE/derivation_cost.py run . --reps 20 --keys 0:16` from this folder.

## Files

| File | Content |
|---|---|
| `RESULTS_derivation_cost.json` | machine, packages, design, SHA-256 of the script and of the paper's modules, summary per arm and ratios |
| `runs.csv` | one row per run: arm, sequence, loop time, import time, the three memory peaks, traced peak, reference check, SHA-256 of the 256 tables |
| `runs_per_table_ns.csv` | 15,360 rows: the three stage times of every timed table |
| `DATA/run_log.txt` | console output of the registered run |
| `DATA/smoke/` | the smoke test disclosed in the plan (keys outside the registered set) |

## Result (laptop: Intel Core Ultra 7 255H, 32 GB, mains power, other programs open; Python 3.12.13, NumPy 2.5.4)

| Arm | Tables per second, median (95% bootstrap interval) | Time per table, ms | Orbit share | Peak working set, MiB |
|---|---|---|---|---|
| S, SHAKE-only | 3,260 (3,154 to 3,324) | 0.276 | 0.4% | 167 |
| H, Float64 orbit + SHAKE | 44.6 (44.3 to 44.9) | 22.2 | 97.9% | 167 |
| L, lattice orbit + SHAKE | 410 (404 to 415) | 2.34 | 84.8% | 167 |

Ratios of the median time per table: H/S 80.6 (79.7 to 81.8), L/S 8.50 (8.35 to 8.63), L/H 0.105 (0.104 to 0.107). The memory
peak is set by the imports (NumPy, SciPy, Matplotlib, pandas, Numba); it rose by at most 0.29 MiB during the derivations, and the
Python allocations of 16 derivations peaked at 8 to 9 KiB. Every timed table equalled the paper's reference function and the
tables were identical across the 20 runs of each arm. Spearman's correlation of throughput with the run order: -0.29 (S),
0.19 (H), -0.27 (L).

Compared with the expectation written in the plan: S was the fastest and H the slowest, as expected; the memory peaks differed
by less than 1 MiB. The smoke test had shown a higher throughput for H (73 tables per second, 16 tables per run) than the
registered run (44.6); the registered figure is the one reported.

These numbers describe research code in Python on one laptop; they bear on cost only, not on novelty, efficiency in general or
any device. Paper: revision 10, Section 9.5 and Table 11 (`05_CODE/make_rev10_tables.py` writes the macros and the table).
