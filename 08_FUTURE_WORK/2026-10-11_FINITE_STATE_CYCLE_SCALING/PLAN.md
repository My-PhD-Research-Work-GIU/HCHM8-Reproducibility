# Finite-state cycle scaling — plan

Date: 2026-10-11. Status: exploratory workstream; not merged into the paper or thesis.

## Question

Table 3 of the paper reports "no repeat" for quantized-state HCHM-8 within caps of 10^5–10^6 iterations. For Q6 the state space has 2^48 points, so even a random mapping would not be expected to close before about 2^24 ≈ 1.7×10^7 iterations. The capped test therefore cannot distinguish a good finite-state map from a poor one. This workstream measures what actually happens.

## Design (fixed before running)

- Same update rule as `quant_step()` in `05_CODE/run_validated_experiments.py`: reference profile, Float64 transcendental evaluation, round-to-nearest to Q bits per coordinate.
- No censoring: Brent cycle detection runs until the orbit closes; tail length μ and cycle length λ are both recorded.
- Q = 2, 3, 4, 5, 6 (state sizes 2^16 … 2^48), 64 initial states per Q, derived from SHAKE256 with label `R01-CYCLE-SCALING-2026-10-11`.
- Control: the hyperbolic core A alone on the same lattice (all shear amplitudes zero), Q = 2, 3, 4, 5, 6, 8.
- Reference: random-mapping expectation E[μ] = E[λ] = sqrt(πN/8) for a state space of N points.
- Outcomes: mean and median μ, λ, ρ = μ + λ; number of distinct cycles reached; slope of log2(mean ρ) against log2 N.

## Not covered

Other parameter profiles, Q ≥ 7, an integer-only implementation, and any cryptographic conclusion.

## Addition 2026-10-12: out-of-sample test at 56 state bits (written before the run)

The Q2–Q6 fit (`RESULTS_PROFILES.json`) was obtained on 16–48 state bits. Q7 (56 bits, N = 2^56) is run afterwards as a test of the extrapolation, and is **not** added to the fit.

- Design: profiles 0–11 (profiles 10 and 11 are new), 2 starts each, same derivation rules, `cycle_scaling_profiles.py 7 0 12 2`.
- Predictions fixed in advance:
  - random mapping: expected tail + cycle = sqrt(pi N / 2) = 3.36e8 iterations (log2 = 28.33);
  - Q2–Q6 fit: log2(profile mean) = -0.0076 + 0.5096 * 56 = 28.53, i.e. 3.88e8 iterations.
- Criterion: if random-mapping behaviour extends to 56 bits, the mean over profiles of (profile mean ÷ 3.36e8) should lie between 0.7 and 1.3. For Q2–Q6 the median of this ratio was 0.84–1.06 and individual profiles ranged from 0.44 to 1.56.
- Whatever the outcome, it is reported.
