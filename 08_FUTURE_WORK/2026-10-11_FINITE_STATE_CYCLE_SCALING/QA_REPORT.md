# Finite-state cycle scaling — QA report

Date: 2026-10-11. Environment: Python 3.12, NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.0, Windows 11, 16 logical cores.

## Checks

- [x] All 320 HCHM runs (5 state sizes × 64 starts) and all 384 core runs closed; none was censored.
- [x] Core control behaves as an exact permutation: tail length 0 in every run.
- [x] Core maximum period equals 3·2^(Q−1) for every tested Q ≥ 3 (12, 24, 48, 96, 384).
- [x] Consistency with the paper's Table 3: at Q6 the shortest measured tail is 129,862 iterations, so no orbit closes within the original 100,000 cap. The earlier "no repeat" observation is reproduced, and explained.
- [x] `analyze.py` regenerates `RESULTS.json` and the figure from `DATA/*.csv` alone.

## Reproduce

```bash
python CODE/cycle_scaling.py 6 64 hchm 400000000000 DATA/hchm_Q6.csv 13
python CODE/cycle_scaling.py 6 64 core 100000000 DATA/core_Q6.csv 8
python CODE/analyze.py
```

Q6 took about two minutes per start. On Windows the process pool occasionally fails at start-up (`WinError 6`); rerunning the same command succeeds and does not affect results, which are deterministic.

## Limitations

- One parameter profile (the reference profile). The 64 starts at each Q share one functional graph, so they are not independent samples of "a random mapping"; between-profile variance is not measured.
- Five state sizes; the slope interval is wide (see `RESULTS.json`).
- Transcendental functions are evaluated in Float64, as in the original emulator. This is not an integer-only implementation and results may differ across math libraries.
- Q ≥ 7 was not run (estimated at more than an hour per start).
