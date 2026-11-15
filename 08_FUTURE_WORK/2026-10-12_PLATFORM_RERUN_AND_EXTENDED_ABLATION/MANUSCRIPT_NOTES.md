# Platform rerun and extended ablation — results and what may be written

Evidence class: observed experiment. Environment "windows11": Windows 11, Python 3.12.13, NumPy 2.5.4, SciPy 1.18.1, Numba 0.68.0. The stored results come from the original package run on Linux. Files: `RESULTS_windows11.json`, `DATA/family_metrics_windows11.csv`, `DATA/family_lattice_metrics_windows11.csv`. Results of the Linux workflow, once committed, carry the tag `linux-ci`.

## A. Map numerics

| Quantity | Stored (Linux) | Rerun (Windows) |
|---|---|---|
| Largest numerical − exact Lyapunov difference | 0.001391 | 0.001180 |
| Largest difference between the two sets of numerical means | — | 0.00038 |
| Iterations to exceed 1e-12, 1e-9, 1e-6, 1e-3 | 6, 9, 12, 16 | 6, 9, 12, 16 |
| Largest inverse / conjugacy error | 4.44e-16 / 8.88e-16 | same |

The digits of trajectory-dependent numbers differ between environments; the conclusions do not. HAS-1 matrices and metrics are byte-identical (integer arithmetic).

## B. Portability of the table corpus — a defect of the published generator

| Arm | Tables identical to the stored corpus |
|---|---|
| SHAKE-only | 4,096 of 4,096 |
| HCHM+SHAKE | 2,707 of 4,096 (66.1%) |

Cause: the generator serializes a Float64 orbit after 64–80 iterations of a map with largest exponent 2.06. One differing last bit in any transcendental value changes the serialized state completely, and SHAKE256 then gives an unrelated table. The floating-point derivation is deterministic only within one environment.

Must be written in the paper. Must not be written: that HCHM-8 can be used as a key-dependent layer in an interoperable system as it stands.

## C. Extended ablation (4,096 tables per arm, regenerated corpus)

No statistic differs: smallest p = 0.44, Holm-corrected 1.00. Values are in `RESULTS_windows11.json` and in the table of the paper, which is generated from that file. Supplementary tests for the two near-discrete statistics: Fisher exact p = 0.89 (SAC maximum deviation), Mann–Whitney p = 0.78 (BIC maximum correlation).

BIC has now been measured for the family; the earlier limitation "BIC was not measured" is closed.

## D. Lattice serialization (mitigation)

Lattice with 16 bits per coordinate, 2,048 serialized bits per round, domain label `HAS1-FAMILY-V3-LATTICE16`.

- Unchanged under the last-bit perturbation for 4,096 of 4,096 key and round pairs; a perturbation scaled by 10^12 changed 256 of 256 tested.
- 4,096 distinct tables. Against the SHAKE-only arm: smallest p = 0.065 (opposite fixed points), Holm-corrected 0.52.

May be written: lattice serialization removed the dependence in these tests. Must not be written: that it guarantees portability. A value within about 1e-10 lattice units of a rounding boundary can still flip; a guarantee needs integer-only or correctly rounded evaluation.

## Provenance note

`run_family_fast.py` as packaged writes a fixed-budget subset of 128 selections per arm (`fixed_budget_subset_256.csv`), whereas the stored `fixed_budget_subset_32.csv` and the numbers quoted in v1 and v2 (16 per arm, median counters 5.0 and 4.5) come from a different configuration that is not in the package. The v3 manuscript does not use those numbers; the selection behaviour is covered by the SBOX LAB study.
