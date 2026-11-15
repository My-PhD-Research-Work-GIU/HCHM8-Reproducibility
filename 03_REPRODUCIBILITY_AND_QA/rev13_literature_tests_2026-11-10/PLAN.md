# Plan: a test of Grebogi, Ott and Yorke (1988) on the stored orbits (written 2026-11-10, before the computation)

## Why

The full texts of Grebogi et al. (1988) and Binder and Jensen (1986), supplied by the researcher on 2026-11-10, give
quantitative predictions for rounded chaotic maps that the paper can test without new runs. Grebogi et al. predict, for a mixing
system with roundoff and l starting states (valid while l <p> << 1), an expected number of distinct periodic orbits
N(l) = sum_{j=1..l} 1/(2j - 1), independent of the roundoff; they also quote Kruskal's value (1/2) ln N + gamma for the total
number of cycles of a random mapping of N points. HCHM-8 is volume preserving and Anosov, so it is mixing with the uniform
measure, d = 8 and <p> = 1/N; their period scaling eps^(-d/2) then reads sqrt(N) for eps = 2^-Q.

## Computation (CODE/distinct_cycles.py; stored data only)

- Distinct cycles reached by the l starting states of each profile, from `cycle_id` in `profiles_Q2.csv` to `profiles_Q7.csv` of
  the uncapped measurement (l = 32 at Q = 2 to 5, 13 at Q = 6, 2 at Q = 7). Mean over profiles against N(l).
- Total number of cycles of F_Q on every state at Q = 2 and 3 (stored exhaustive graphs) against the exact random-mapping
  expectation and Kruskal's asymptotic value.

## Criterion (fixed here)

For Q = 3 to 6, the mean number of distinct cycles over profiles lies within 0.67 to 1.5 times N(l). Q = 2 (four points per
coordinate) and Q = 7 (l = 2) are reported only. The total-cycle comparison is descriptive.

## Disclosure

The data exist since 2026-10 and parts of them were seen: the paper already states that each profile reached one to five
distinct cycles, and that the total number of cycles at 2^24 states (mean 13.1) exceeds the random-mapping value (9.0). The ratio
to N(l) has not been computed before this plan. The test is therefore a registered analysis of existing data, not a new
experiment, and the paper will say so.
