# Platform rerun and extended ablation — QA report

- [x] The original module is imported unchanged from a temporary copy; nothing it writes lands in the repository.
- [x] SHAKE-only arm: 4,096 of 4,096 tables identical to the stored corpus. This confirms that the regeneration code path is the published one and that the mismatch in the other arm is not a scripting error.
- [x] 8,192 distinct tables in the regenerated corpus.
- [x] The last-bit perturbation test has power: the same perturbation scaled by 10^12 changes every tested serialization.
- [x] Results reproduce when the script is run twice in the same environment (the JSON was regenerated after adding the supplementary tests; values unchanged).
- [x] The paper table and number macros are generated from `RESULTS_windows11.json` by `05_CODE/make_paper_tables.py`; `check_manuscript.py` fails if they are stale.

- [x] The original 46-check QA script (`03_REPRODUCIBILITY_AND_QA/run_final_qa.py`), run inside the unpacked original package against the data regenerated on Windows 11, passes 46 of 46. Three recorded values differ from the original run, as expected for trajectory-dependent quantities: largest error of det DT (7.6e-14 against 5.4e-14), agreement of the numerical spectrum (0.001180 against 0.001391), and the p-value of the 512-table SAC comparison (0.384 against 0.524).
- [x] The capped quantized-state searches of the original package reproduce exactly on Windows 11 (`quantized_cycle_20runs.csv` byte-identical: no repeat within the caps).
- [x] Third environment (GitHub workflow, Linux): `RESULTS_linux-ci.json`. SHAKE-only 4,096 of 4,096 against both other environments; HCHM+SHAKE 3,962 of 4,096 against the stored corpus and 2,774 of 4,096 against Windows; lattice serialization 4,096 of 4,096 against Windows; extended ablation not significant (smallest Holm p 0.57).

Reproduce: `python 05_CODE/platform_rerun_and_extended_ablation.py <tag> 08_FUTURE_WORK/2026-10-12_PLATFORM_RERUN_AND_EXTENDED_ABLATION` (about one minute). The GitHub workflow runs it with tag `linux-ci`.

Limitations: one key set (indices 40000–40255, 16 rounds); the supplementary tests are post hoc; the lattice serialization was tested on 4,096 key and round pairs only.
