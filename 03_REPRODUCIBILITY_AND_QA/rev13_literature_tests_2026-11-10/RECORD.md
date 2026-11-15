# Revision 13: four sources read in full, and a registered test of Grebogi et al. (2026-11-10)

The researcher supplied the full texts on 2026-11-10 (copies not stored in the repository; SHA-256 below). Before anything was
added, a critical judgment was given in the session: all four are worth citing precisely; none may be described as surpassed.

| Source | SHA-256 of the copy read | What it establishes (own words) | Use in revision 13 |
|---|---|---|---|
| Binder and Jensen (1986), Phys. Rev. A 34, 4460 | aaf04f0e...29c06 | logistic map truncated to m states: cycles and transients of order sqrt(m), about 90% of states on one cycle, at most three or four cycles; invariant density and Lyapunov exponent preserved on the cycle; reports Rannou (1974): a bijective, area-preserving discretization of the standard map behaves like a random permutation | Section 2.4; Section 7.4 contrast: the conjugacy-preserving bijection G_Q has the core's cycles instead (Rannou cited as reported) |
| Grebogi, Ott, Yorke (1988), Phys. Rev. A 38, 3688 | 9734aa24...1fea0 | mixing systems with roundoff eps: expected period eps^(-d/2), d the correlation dimension; period distribution; expected number of distinct periodic orbits from l starts, sum 1/(2j-1); quotes Kruskal's (1/2) ln N + gamma for all cycles of a random map; tested on the Ikeda map (d = 1.67 estimated) | Section 2.4 (d = 8 exactly for HCHM-8, so sqrt(N)); Section 7.4 and Table S7: the distinct-cycles prediction tested on the stored orbits |
| Hua et al. (2018), IEEE Trans. Cybern. 48, 463 | f0669e9d...c06fe | linear n-D cat maps; exponents = logs of eigenvalues of an integer matrix built by a similarity transformation of two basic cat matrices; randomness shown by entropy and phase plots on finite states | Section 2.4 and Table S7: linear exactness against exactness under nonlinear shears; in finite state such maps have the periods of their matrix |
| Ozkaynak (2020), Physica A 550, 124072 | 8c63905a...6b158 | S-boxes from the Lorenz system with periodic (b = 21) and chaotic (b = 28) parameters; periodic-source S-boxes no worse in average nonlinearity, SAC, BIC and XOR table; 54 of 100 better; no statistical test | Sections 2.3 and 9.2, Table S7, with the caveat that our tables have the inverse's nonlinearity and differential uniformity by construction |

## Registered test (PLAN.md, commit 5063cba; workflow run created 2026-11-10T10:38:20Z; computed afterwards)

Distinct cycles reached from l starts, mean over profiles, against N(l) = sum 1/(2j-1) (`DATA/RESULTS_distinct_cycles.json`):

| State bits | l | Mean (range) | N(l) | Ratio | Band 0.67 to 1.5 |
|---|---|---|---|---|---|
| 16 | 32 | 3.10 (2 to 5) | 2.71 | 1.14 | reported only |
| 24 | 32 | 2.80 (1 to 4) | 2.71 | 1.03 | met |
| 32 | 32 | 1.90 (1 to 3) | 2.71 | 0.70 | met |
| 40 | 32 | 2.50 (1 to 4) | 2.71 | 0.92 | met |
| 48 | 13 | 2.25 (1 to 3) | 2.26 | 0.99 | met |
| 56 | 2 | 1.33 (1 to 2) | 1.33 | 1.00 | reported only |

All cycles on every state (descriptive): 7.1 against 6.18 (exact random-mapping value) at 2^16 states, 13.1 against 8.95 at 2^24
states. The distinct-cycles prediction holds; the total number of cycles exceeds the random-mapping value, as already reported.

## Not claimed

That HCHM-8 or its discretizations improve on these works; that Lyapunov exponents are preserved on the cycles of F_Q (not
measured); Rannou's result beyond what Binder and Jensen report.
