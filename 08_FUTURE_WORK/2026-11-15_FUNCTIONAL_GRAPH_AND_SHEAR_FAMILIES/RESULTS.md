# Results: functional graphs, in-degree and shear families (2026-11-15)

Plan: `PLAN.md` (written before the registered runs; SHA-256 and start time in `DATA/run_log.txt`; not committed, because commits
were not allowed in this round). Raw results: `DATA/RESULTS_fg_study_part1.json`, `DATA/RESULTS_fg_study_part2.json`; analysis
`CODE/analyze.py` -> `DATA/SUMMARY.json`; independent checker `CODE/independent_check.py` -> `DATA/INDEPENDENT_CHECK_Q2.json`.

| Hypothesis (PLAN.md) | Result | Verdict |
|---|---|---|
| H1 implementation check of G_Q | 120 cases (Q = 2, 3; 10 profiles; families 1 to 3 and strengths 0.25, 0.5, 2): H_Q bijective, H_Q^-1(H_Q(x)) = x, G_Q = H_Q^-1 A_Q H_Q on every state, cycle type of A_Q, period of every x equal to the A_Q period of H_Q(x); no successor outside the lattice; with zero shears F_Q = A_Q on every state | met (consequence of conjugacy, not a new theorem) |
| H2 in-degree of F_Q, HCHM-8, 2^24 states | in-degree 0 to 4: 0.3674, 0.3683, 0.1843, 0.0612, 0.0152 against Poisson(1); largest difference 0.0005; image fraction 0.633; max in-degree 9 to 11; 26.4% of states have two or more predecessors; largest basin 0.73 of the states; 13.1 cycles (random mapping 8.95); mean tail 1.26 times sqrt(pi N / 8) | met |
| H3 other shear families | 2^24 states: image fraction 0.635 (single sine) and 0.637 (sine of a difference), rho ratio 0.89 and 0.79; slopes over 16 to 40 bits 0.493 (0.475 to 0.511) and 0.513 (0.497 to 0.528), HCHM-8 0.510 (0.499 to 0.520) in the same runs | met |
| H4 rounding rule | image fraction 0.633 (nearest) and 0.631 (truncation); rho ratio 0.92 and 0.87 at 2^24, 1.03 and 0.93 at 32 and 40 bits | met |
| Exploratory: strength | image fraction at 2^24 for amplitudes x0, 0.25, 0.5, 1, 2: 1.000, 0.722, 0.658, 0.633, 0.626 | no criterion |
| Exploratory: stepwise quantization | nearest: equal to G_Q on every state (no ties occurred); truncation: always a permutation (each step is a bijection), mean cycles 4,213 / 649,626 / 7,644 for the three families at 2^24, against 1,444,656 for A_Q and 17.2 for a random permutation | no criterion |
| Reproduction of the revision 9 registered run | 20 of 20 lattice and profile cases identical | reproduced |
| Independent checker (Q = 2, profiles 0 and 1) | 0 successor differences for F_Q, G_Q, A_Q; in-degree, cycles, cycle type, mean rho and longest tail equal | agreed |

Limits: exhaustive only at 2^16 and 2^24 states; Q = 4, 5 are sampled (8 starts per profile); families 2 and 3 not measured above
40 bits; increments and shears evaluated in Float64. Collision counts show non-injectivity directly; they do not explain or prove
the scaling exponent.
