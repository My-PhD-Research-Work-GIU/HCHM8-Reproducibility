# HCHM-8 reproducibility package

Code and data behind the results of the article *Two Discretizations of HCHM-8, an Eight-Dimensional Chaotic Map with Exactly
Known Invariants: Finite-State Behaviour and Consequences for Chaos-Based S-Boxes* (manuscript; authors and affiliations as listed
in the article). Release `v1.0.1` corresponds to paper revision 15 (2026-11-15); it corrects the checksum file of `v1.0.0`
without changing any script or result (see `CHANGELOG.md`).

The folder layout follows the project from which the package was taken, so that the scripts find each other by relative path.
Every result file was produced by the script next to it; plans written before the runs are kept beside their results.

## Contents and where each result appears in the article

| Folder | What it reproduces | Article |
|---|---|---|
| `05_CODE/run_validated_experiments.py` | HCHM-8, the HAS-1 reference S-box, the S-box family and its metrics | Sections 3, 4, 8, 9 |
| `05_CODE/verify_core_claims.py` | independent checks of the analytical statements (spectrum, entropy, fixed points, orders) | Sections 4, 7.3, App. B |
| `05_CODE/platform_rerun_and_extended_ablation.py` | extended ablation and portability across computing environments | Sections 9.2, 9.4 |
| `05_CODE/run_family_fast.py` | fast generation of the S-box family | Section 9 |
| `04_DATASETS/` | HAS-1 tables (DDT, LAT, BCT, DLCT, SAC, inverse), stored summaries | Sections 6, 8 |
| `08_FUTURE_WORK/2026-10-11_FINITE_STATE_CYCLE_SCALING` | uncapped orbits of the rounded emulator, 16 to 56 state bits; library-perturbation test | Sections 7.1, 7.2, 7.5 to 7.7 |
| `08_FUTURE_WORK/2026-11-10_GRID_CONJUGACY` | grid-exact conjugate G_Q against the rounded emulator F_Q and the core A_Q | Section 7.4, Proposition 2 |
| `08_FUTURE_WORK/2026-11-15_FUNCTIONAL_GRAPH_AND_SHEAR_FAMILIES` | functional graphs, in-degree distributions, other shear families, rounding rule and amplitude | Section 7.5, Online Resource 1 S11 |
| `08_FUTURE_WORK/2026-11-10_CHAOS_ONLY_ARM` | ablation in both directions (without the map, without the hash) | Section 9.3 |
| `08_FUTURE_WORK/2026-10-12_PLATFORM_RERUN_AND_EXTENDED_ABLATION` | stored metrics of the extended ablation in two environments | Sections 9.2, 9.4 |
| `08_FUTURE_WORK/2026-11-10_DERIVATION_COST` | cost benchmark, raw timings | Section 9.5, Table 11 |
| `08_FUTURE_WORK/2026-10-24_FAMILY_EQUIVALENCE` | plan and results of the complete affine-equivalence decision (sboxU, PEIGEN); the tool drivers are not included | Section 9.6 |
| `03_REPRODUCIBILITY_AND_QA/rev9_tests_rerun_2026-11-10` | rerun and independent check of the tests of Sections 7.4 and 9.3 | Online Resource 1 S9 |
| `03_REPRODUCIBILITY_AND_QA/rev13_literature_tests_2026-11-10` | distinct-cycles test of Grebogi et al. on the stored orbits | Section 7.4 |
| `03_REPRODUCIBILITY_AND_QA/rev15_validation_2026-11-15` | independent recomputation of the mathematics, the statistics and HAS-1 | revision 15 checks |

## Environment

Python 3.12 with the packages in `requirements.txt` (numba is needed for the finite-state experiments). The results were produced
on Windows 11 and checked on Linux; floating-point orbits of HCHM-8 can differ between mathematical libraries, as the article
reports, while every lattice and integer result is library independent.

```
python -m venv .venv && . .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Examples

```
python 03_REPRODUCIBILITY_AND_QA/rev15_validation_2026-11-15/CODE/math_checks.py 08_FUTURE_WORK/2026-10-11_FINITE_STATE_CYCLE_SCALING/CODE/cycle_scaling_profiles.py out/math.json
python 08_FUTURE_WORK/2026-11-10_GRID_CONJUGACY/CODE/grid_conjugacy.py out --profiles 0:10 --exhaustive-q 2,3 --orbit-q 4,5,6,7
python 08_FUTURE_WORK/2026-11-15_FUNCTIONAL_GRAPH_AND_SHEAR_FAMILIES/CODE/independent_check.py out/check_q2.json
python 03_REPRODUCIBILITY_AND_QA/rev15_validation_2026-11-15/CODE/sbox_audit.py 05_CODE/run_validated_experiments.py <main.tex of the article> out/sbox.json
```

`05_CODE/verify_core_claims.py` also compares its values with Table 1 of the manuscript source, which it locates through the
project module `current_versions.py`; neither is part of this release, so that comparison needs the article's LaTeX source.
Each experiment folder has its own `PLAN.md`, and most have `REPRODUCE.md` or a usage line at the top of their script. The
exhaustive runs at 2^24 states take minutes per profile on 16 threads; the uncapped orbits at 48 and 56 state bits take hours.

## Integrity

`SHA256SUMS.txt` lists the SHA-256 of every file of this release except itself. It was computed on the committed contents, and
`.gitattributes` (`* -text`) makes Git store and check out every file byte for byte, so the list verifies on a clone and on the
release zip:

```
sha256sum -c SHA256SUMS.txt
```

Erratum for `v1.0.0`: in that release Git converted the line endings of 53 result files (CSV, JSON and TXT) from CRLF to LF at
commit, after the checksums had been computed, so 53 entries of its checksum file do not verify. The files differ from the
checksummed originals only in line endings. `v1.0.1` restores the original bytes; `v1.0.0` is kept unchanged.

## Not included

The manuscript and its supplementary file, the affine-equivalence software used for auxiliary checks (EASAT) and the drivers of
sboxU and PEIGEN, and project material unrelated to the article.

## Licence

No licence is granted at present: the authors keep the copyright. The repository is public so that editors, reviewers and readers
can inspect the code and data and rerun them to verify the article. For any other use, contact the corresponding author named in
the article.
