# Plan: chaos-only arms of the S-box ablation (written 2026-11-10, before the registered runs)

## Question

The extended ablation of the paper compares HCHM+SHAKE with SHAKE-only and detects no difference in eight non-affine-invariant
statistics. A reviewer can ask whether SHAKE256 hides an effect of the chaotic layer. This test removes the hash instead of
the map: the affine layers are driven directly by the orbit.

## Arms (same 256 keys x 16 rounds as the paper: key_for_index(40000 + k), k = 0 to 255, rounds 0 to 15; candidate 0)

| Arm | 258-byte stream that drives `make_rows` and the constants a, b |
|---|---|
| S, SHAKE-only | SHAKE256 of key, round and labels (the paper's arm, reproduced exactly) |
| H, HCHM+SHAKE | SHAKE256 of key, round, labels and the Float64 orbit (the paper's arm, reproduced exactly) |
| C, chaos-only | the first 258 bytes of the Float64 orbit of HCHM8(key, round) after 64 burn-in steps, 32 bits per coordinate as in `serialize_hchm`, no hash |
| CL, chaos-only lattice | the first 258 bytes of the 16-bit lattice orbit of the paper (`lattice_serialize`, 64 burn-in steps), no hash |

The script asserts, for every key and round, that S and H give exactly `candidate(base, 0)` of the paper's derivation.

## Statistics and tests (as in the paper)

SAC mean absolute deviation and maximum deviation, BIC maximum and mean absolute correlation (two-sample Kolmogorov-Smirnov);
fixed points, opposite fixed points, number of cycles, longest cycle (Mann-Whitney U). Holm correction over the eight tests of
each comparison. Effect sizes: Hedges' g with a 95% interval, Cliff's delta.

- Primary comparisons: C against S, and CL against S.
- Secondary: C against H. The paper's comparison H against S is recomputed as a check.
- Descriptive: chi-square of the byte frequencies of the C and CL streams against uniform; distinct tables per arm; tables
  shared between arms.

## Predictions and criteria

The orbit of a volume-preserving map is equidistributed for almost every start, and the derivation uses 64 row operations per
affine layer, so the statistics are expected to be insensitive to the source of the bytes. Prediction, fixed here: in both
primary comparisons no test is significant after Holm correction (all Holm p >= 0.05) and every 95% interval of g lies within
+/-0.10. Whatever the outcome, it is reported; a significant difference would be reported as evidence that the raw orbit is a
weaker byte source for this derivation than SHAKE256.

Plan and script are committed and pushed to GitHub before the registered run, and the commit is opened as a pull request so
that the hosting service records the time of the first workflow run.

## Smoke test before this plan

The script was run once to check that it executes, on keys key_for_index(50000) and key_for_index(50001) (not among the
registered keys), 32 tables per arm. Its p-values and effect sizes are not informative at that size and were not used to set
the criteria. Seen before the plan was written: the assertions on S and H passed, the 32 tables of each arm were distinct, and
the byte chi-square p-values of C and CL were 0.27 and 0.23.

## Not covered

Other derivations (selection over several candidates, other affine-layer constructions), other maps, and any security
property. Every table of every arm is affine-equivalent to the field inverse by construction, as in the paper.
