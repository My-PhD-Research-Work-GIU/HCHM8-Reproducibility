# Plan: cost of three S-box derivations on one machine (written 2026-11-10, before the registered runs)

## Question

The paper finds that the HCHM-8 layer changes none of eight measured S-box statistics, that its floating-point evaluation is
not portable, and that a lattice serialization restored portability in the cases tested. What does each derivation cost?
This is a limited benchmark of the research implementations on one machine. It completes the cost side of the comparison;
it is not evidence of novelty, of efficiency in general, or of suitability for any device, and it is not compared with
published timings.

## Arms (one table = candidate 0 of the paper's derivation for one key and round)

| Arm | Key material in front of SHAKE256 | Source |
|---|---|---|
| S, SHAKE-only | key, round, label | `candidate()` of `05_CODE/run_validated_experiments.py` |
| H, HCHM+SHAKE | key, round, label, Float64 orbit (`serialize_hchm`, 64 burn-in, 16 states, 32 bits per coordinate) | same |
| L, lattice+SHAKE | key, round, label `V3-LATTICE16`, lattice orbit (`lattice_serialize`, 16 bits per coordinate) | Part D of `05_CODE/platform_rerun_and_extended_ablation.py`, verbatim |

Every timed table is compared, after the timing, with the reference function (`candidate(base, 0)`); a mismatch stops the
analysis and is reported.

## Measurement (fixed here)

Script: `CODE/derivation_cost.py run <out> --reps 20 --keys 0:16`.

- Keys `key_for_index(60000 + k)`, k = 0 to 15, rounds 0 to 15: 256 tables per run. These keys are used by no earlier test.
- 20 runs per arm, 60 runs in all. Each run is a fresh Python process that measures one arm; the order of the arms cycles
  through the six permutations of (S, H, L), one per repetition, so that drift of the machine is spread over the arms.
- Per run: 4 warm-up tables on keys `key_for_index(70000 + i)`, then `time.perf_counter_ns` around three stages of every
  table: material (parameters, orbit and serialization; empty for S), expansion by SHAKE256 to 258 bytes, and table
  construction (`make_rows` twice and `build_sbox`); and around the whole loop.
- Memory: peak working set of the process (Windows; peak resident set on Linux) after imports, after the warm-up and after the
  loop; and the peak of Python allocations traced by `tracemalloc` in a separate, untimed pass over 16 tables.
- The machine, its processor, memory, power source and battery state, and the package versions are recorded by the script.
  Nothing else is controlled: other programs may run, the processor clock is not fixed, and the run is made on the
  researcher's laptop as it is used.

## Reported quantities

- Throughput: tables per second per run (256 divided by the loop time); median over the 20 runs with a 95% bootstrap
  percentile interval (10,000 resamples of runs, seed 20261110), mean, standard deviation and range.
- Time per table: the median within each run, then the median over runs with its bootstrap interval; the pooled 5th and 95th
  percentiles; the median share of each stage.
- Ratios H/S, L/S and L/H of the median time per table, with bootstrap intervals in which the runs of each arm are
  resampled separately.
- Memory: median and range of the three peaks per arm, the increase of the peak during the loop, and the traced peak.
- A drift check: Spearman's correlation between the run sequence number and the throughput, per arm.

No hypothesis test is planned and no threshold is set: the benchmark describes these implementations on this machine.
Expectation, written down so that it can be compared, not a criterion: S is the fastest; H and L are dominated by the 80
iterations of the map, H more than L because `serialize_hchm` steps a NumPy array object per iteration; the memory peaks
differ little, because the imports dominate them.

## Smoke test before this plan

The script was run once to check that it executes, with keys `key_for_index(60100)` (outside the registered set), 2 runs
per arm, 16 tables per run. Seen before the plan was written: every reference check passed and the tables were identical
across runs; medians of 5,135 (S), 73 (H) and 569 (L) tables per second; peak working set 167 MiB in every arm. These
numbers come from too few runs to be reported and did not change the design.

## Not covered

Other machines, compiled or constant-time implementations, the 64-candidate selection of the family generator (only
candidate 0 is derived), the cost of the statistical tests, energy, and any comparison with published generators.
