# Finite-state cycle scaling — manuscript notes

Evidence class: **observed experiment**, single parameter profile. Not yet in the paper or thesis.

## Results

| Q | State bits | Runs | Median tail | Median cycle | Mean tail+cycle | Ratio to random mapping | Distinct cycles reached |
|---|---|---|---|---|---|---|---|
| 2 | 16 | 64 | 190 | 84 | 261 | 0.81 | 2 |
| 3 | 24 | 64 | 5,482 | 76 | 5,378 | 1.05 | 4 |
| 4 | 32 | 64 | 38,832 | 53,404 | 91,560 | 1.11 | 3 |
| 5 | 40 | 64 | 438,098 | 1,015,648 | 1.42×10^6 | 1.08 | 2 |
| 6 | 48 | 64 | 7,147,454 | 8,523 | 9.04×10^6 | 0.43 | 5 |

Fit: log2(mean ρ) = 0.80 + 0.478 · (state bits), 95% CI for the slope 0.404–0.551, R² = 0.993, five points. The random-mapping slope 0.5 lies inside the interval.

Hyperbolic core alone (no shears): exact permutation, maximum period 12, 12, 24, 48, 96, 384 for Q = 2, 3, 4, 5, 6, 8.

At Q6, 41 of 64 starts end on one cycle of length 8,523; the others end on cycles of length 209,880 to 10,196,637.

## What may be written

- The quantized-state emulator behaves like a random mapping on its state space: orbits close after a number of iterations of order sqrt(N), not of order N.
- The earlier capped searches were below this scale and could not have observed a repeat; the "no repeat" rows carry no information about cycle length.
- On the same lattice the unconjugated toral automorphism has period at most 3·2^(Q−1). The nonlinear shears, together with rounding, lengthen the observed orbit by several orders of magnitude at every tested size (at Q6: about 9×10^6 against 96).
- Short cycles with large basins exist (8,523 at 48 state bits). A finite-state deployment would need an explicit cycle-avoidance or re-seeding rule.

## What may not be written

- Any statement for Q ≥ 7. The sqrt(N) behaviour at Q = 8–16 is an extrapolation.
- That the map "reaches" or "achieves" the random-mapping bound; one profile and five sizes do not establish that.
- Any cryptographic security conclusion.

## Suggested use in the paper

Replace Table 3 (or add beside it) with the uncapped measurements and the core control, and reword Section 5 accordingly. This gives the finite-precision section a measured result and gives a concrete, limited answer to "what do the shears contribute": in finite state they remove the short lattice period of the linear core, up to the sqrt(N) scale.

## Update 2026-10-12: several parameter profiles

`CODE/cycle_scaling_profiles.py`, `CODE/analyze_profiles.py`, `DATA/profiles_Q*.csv`, `RESULTS_PROFILES.json`. The parameter profile is the unit of analysis; each profile contributes its mean tail + cycle.

| Q | State bits | Profiles × starts | Median tail | Median cycle | Shortest cycle | Profile mean ÷ random mapping (min / median / max) |
|---|---|---|---|---|---|---|
| 2 | 16 | 10 × 32 | 137 | 99 | 2 | 0.68 / 0.84 / 1.06 |
| 3 | 24 | 10 × 32 | 2,763 | 1,189 | 76 | 0.67 / 0.86 / 1.53 |
| 4 | 32 | 10 × 32 | 38,708 | 52,531 | 254 | 0.54 / 1.06 / 1.56 |
| 5 | 40 | 10 × 32 | 548,432 | 770,929 | 939 | 0.74 / 1.06 / 1.52 |
| 6 | 48 | 8 × 13 | 6,669,904 | 11,725,216 | 8,523 | 0.44 / 1.02 / 1.55 |

Fit over 48 profile points: log2(mean tail + cycle) = −0.008 + 0.5096 · (state bits); 95% CI for the slope 0.4989–0.5203; R² = 0.995. Random mapping: slope 0.5, intercept 0.326.

What changes relative to the single-profile notes above:

- The sqrt(N) scaling is now supported across profiles, with a narrow interval that contains 0.5. "Consistent with a random mapping over 16–48 state bits" may be written. It remains an observed experiment, not a proof, and says nothing about Q ≥ 7.
- The single-profile ratio of 0.43 at Q6 was the low end of the between-profile spread (0.44–1.55), not a trend.
- The 8,523-cycle belongs to the reference profile. Short cycles are a property of individual profiles: 12% of pooled Q6 starts and 7% of pooled Q4 starts end on a cycle more than 100 times shorter than sqrt(N). Each profile reaches only one to five distinct cycles from its starts.
- Profile 0 reproduces the earlier reference-profile CSV files row for row (Q2–Q4 checked), which validates the re-parameterised code.

Still limited: eight profiles and 13 starts at Q6; Float64 transcendental evaluation; no integer-only implementation.

## Update 2026-10-12 (second): out-of-sample test, exact core period, library dependence

**Test at 56 state bits.** The prediction was written in `PLAN.md` and committed (commit `ae8fd81`, 21:43:54 +05:00) before the run started. Twelve profiles (0–11; 10 and 11 new), two starts each, 24 runs, none censored. Mean over profiles of (profile mean ÷ random-mapping expectation) = 0.97, median 0.97, range 0.38–1.91. The pre-specified interval for the mean was 0.7–1.3: met. Median tail 82,327,518; median cycle 120,137,022; shortest cycle 32,604,056. These runs are not in the fit (`analyze_profiles.py` restricts the fit to Q ≤ 6; the fit is unchanged: slope 0.5096, CI 0.4989–0.5203).

**Exact period of the linear core.** `05_CODE/verify_core_claims.py` computes the multiplicative order of A modulo 2^Q exactly: 12 for Q = 2 and 3·2^(Q−1) for every Q from 3 to 16. This replaces "largest observed period" by an exact statement.

**Dependence on the mathematical library.** `cycle_scaling_profiles.py ... perturb` multiplies every transcendental value by 1 + k·2^−52, k ∈ {−4..4} chosen by a hash of the lattice state.
- Last-bit scale: 0 of 1,056 runs changed (320 each at Q2, Q3, Q4; 96 at Q5). Report: `DATA/perturbation_last_bit_windows11.json`.
- Sanity check with the perturbation scaled by 10^10: 317 of 320 runs at Q4 changed, so the test detects a dependence when there is one. Report: `DATA/perturbation_sanity_scale_1e10_windows11.json`.
- `reproduce_check.py` reruns any part of the data and compares it row by row; `DATA/reproduction_windows11.json` records 0 of 1,056 differing on the machine that produced the data. The GitHub workflow runs the same check on Linux.

What may be written: the finite-state results at these lattice sizes did not depend on last-bit differences in the mathematical library in the tests made. What may not: that they are independent of the library in general, or that this holds for much finer lattices.
