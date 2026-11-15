# Plan: rerun and independent check of the two tests added in paper revision 9 (written 2026-11-10, before the runs)

## Purpose

Revision 9 added two tests whose plans were pushed before their runs: the grid-exact conjugate of HCHM-8 against the rounded
emulator (Proposition 2, Section 7.5, Table 6; `08_FUTURE_WORK/2026-11-10_GRID_CONJUGACY`) and the chaos-only arms of the
S-box ablation (Section 9.3, Table 10; `08_FUTURE_WORK/2026-11-10_CHAOS_ONLY_ARM`). Both were run once, on one machine, by
the code that defines them. This plan checks their code, parameters, raw outputs and reproduction instructions before the
paper relies on them further.

## Checks (fixed here)

1. **Reproduction instructions.** `REPRODUCE.md` in each test folder states the environment, the commands, the expected
   outputs and what may differ between environments. A fresh clone of the pushed commit that contains this plan is made in a
   new folder, a new Python environment is created exactly as `REPRODUCE.md` says (uv, pinned versions), and both commands are
   run as written. `CODE/compare_rerun.py` compares the outputs with the stored ones.
   Expectation: in the same environment every table and every statistic is identical (only time stamps and durations differ).
2. **Second operating system.** The job `rev9-tests-rerun` of the GitHub workflow runs the same commands on Linux with the same
   pinned packages and the same comparison, and then the independent check of item 3.
   Must hold (the job fails otherwise): G_Q is a permutation with the cycle type of A_Q for every profile; every G_Q period
   equals the A_Q period at H_Q(x0) and divides the core period; the S and H arms reproduce the paper's derivation; the
   SHAKE-only tables equal the stored ones.
   Expected and reported, not required: the registered criteria of both plans are met again; the lattice-arm tables are
   identical; the F_Q statistics are identical; the H and C tables, which pass through a Float64 orbit, differ for some keys
   (the portability defect of Section 9.4), so the numbers of the comparisons that involve them can change.
3. **Independent code.** `CODE/independent_check.py` recomputes the results with code written for this check from the
   equations of the paper: F_Q and G_Q with NumPy arrays, their functional graphs with a path-following algorithm (the original
   peels by in-degree), every lattice state at Q = 2 and 3 for the ten profiles (Part A), the 790 orbit periods at Q = 4 to 7
   (Part B), and for the 16,384 tables of the chaos-only test their SHA-256, their eight statistics and every number of
   Table 10 (Part C). Part C also reports what the paper did not: the correlation of each statistic within matched pairs (the
   arms share keys and rounds; the intervals of Table 10 treat them as independent samples) and the bound of the intervals of
   Hedges' g when they are widened for eight statistics (Bonferroni, two-sided 0.05 / 8).
   Expectation: every stored value is recomputed exactly; the within-pair correlations are small.
4. **Parameters.** The registered parameters (profiles 0 to 9, starting states 32, 32, 13, 2 at Q = 4 to 7, keys
   `key_for_index(40000 + k)` for k = 0 to 255, rounds 0 to 15, candidate 0) are checked against the scripts, the stored
   results and the text of Sections 7.5 and 9.3.

Whatever the outcome, it is reported in `REPORT.md` of this folder and carried into the next revision of the paper; a
difference that is not explained by the stated environment dependence is reported as an error of the stored result.

## Smoke test before this plan

`independent_check.py` was run once with Parts A (Q = 2 only) and B to check that it executes. Seen before this plan was
written: at Q = 2 every stored statistic was recomputed, G_Q was a permutation with the cycle type of A_Q for every profile,
F_Q and G_Q agreed with the arrays of the original program on every state, and the orbit periods at Q = 4 to 7 were those
stored. Q = 3 and Part C were not run.

## Not covered

An integer-only evaluation of the shear functions, other machines than the two named, and any change of the tests themselves.
