# Reproducing the grid-exact conjugate test (paper Section 7.5, Proposition 2, Table 6)

Plan: `PLAN.md` (pushed before the registered run). Code: `CODE/grid_conjugacy.py`. Stored result: `RESULTS_grid_conjugacy.json`
and the console output `DATA/run_log.txt`. The tables of the paper are written from the stored result by
`05_CODE/make_rev9_tables.py` (macros `\Gc...` in `tables/rev9_numbers.tex`, table body `tables/tab_grid_conjugacy.tex`).

## Environment of the registered run

Windows 11 (10.0.26300), CPython 3.12.13, NumPy 2.5.4, Numba 0.68.1 (llvmlite 0.50.0); 16 logical processors.
The script imports `step()`, `f_scalar()`, `profile_params()` and `initial()` unchanged from
`08_FUTURE_WORK/2026-10-11_FINITE_STATE_CYCLE_SCALING/CODE/cycle_scaling_profiles.py`; the SHA-256 of that file and of the
script are recorded in the result.

## Steps

From the folder `05_PAPERS_AND_DATASETS/05_PAPERS_AND_DATASETS` of a clone of the repository:

```
uv venv --python 3.12 .venv-rerun
uv pip install --python .venv-rerun numpy==2.5.4 scipy==1.18.1 numba==0.68.1 matplotlib==3.11.2 pandas==3.1.1
.venv-rerun/Scripts/python -X utf8 08_FUTURE_WORK/2026-11-10_GRID_CONJUGACY/CODE/grid_conjugacy.py <out_dir> --profiles 0:10 --exhaustive-q 2,3 --orbit-q 4,5,6,7
```

(on Linux the interpreter is `.venv-rerun/bin/python`). The run took about three minutes on the machine above; Numba compiles
the functions on first use and caches them in a temporary folder, not in the repository.

## Parameters

- Profiles 0 to 9: `profile_params(k)`; profile 0 is the reference profile, profile k >= 1 the parameters of
  `HCHM8(key_for_index(9000 + k))`. The last-bit perturbation switch is off.
- Part 1: every lattice state at Q = 2 (2^16 states) and Q = 3 (2^24 states).
- Part 2: Q = 4, 5, 6, 7 with 32, 32, 13 and 2 starting states per profile, `initial(Q, run)`: 790 orbits.

## Output and comparison

`<out_dir>/RESULTS_grid_conjugacy.json`. Compare with the stored file by
`03_REPRODUCIBILITY_AND_QA/rev9_tests_rerun_2026-11-10/CODE/compare_rerun.py` (which also needs a rerun of the chaos-only
test). In the same environment everything except the time stamps and the durations is identical. In another environment the
cycle type of G_Q and every orbit period must be identical (Proposition 2 does not depend on how the shear functions are
evaluated); the statistics of the rounded emulator F_Q can differ if a function value falls on a rounding boundary.
