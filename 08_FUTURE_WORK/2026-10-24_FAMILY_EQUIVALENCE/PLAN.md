# E-F: complete affine-equivalence decision on the dynamic S-box family (pre-registration)

Written and committed on 2026-10-24 **before any run** on these tables; the commit that adds this file must precede the commit that
adds `DATA/`. Requested by the researcher on 2026-10-24, after E-A (`../2026-10-23_EQUIVALENCE_TOOL_PARITY/`).

## Question

The paper states that HAS-1 and every member of the key-dependent family are affine-equivalent to the GF(2^8) inverse function. This
follows from the construction `S(x) = L_out[(L_in(x) XOR a)^-1] XOR b`, but no complete decision procedure has checked the tables that
the implementation actually produces. Does one?

## Data (fixed)

- All 8,192 tables of the family, 4,096 per arm (HCHM+SHAKE and SHAKE-only; keys 40000 to 40255, rounds 0 to 15), regenerated with
  `05_CODE/run_validated_experiments.py` (`candidate(base, 0)`, exactly as in `05_CODE/platform_rerun_and_extended_ablation.py`) in
  the WSL Ubuntu environment of E-A.
- Reference: the GF(2^8) inverse with 0 -> 0 (AES polynomial).
- Observation only (no hypothesis): the share of regenerated tables whose SHA-256 equals `04_DATASETS/dynamic_family_direct_8192.csv`.
  The paper already reports that HCHM+SHAKE tables are not reproducible across environments, so the tables tested are the ones
  generated here, whatever that share is.

## Tools (pinned, as in E-A)

sboxU bedb6d1 (`affine_equivalence_permutations`) under SageMath 10.5; PEIGEN 6fe84ab (`is_AE`); EASAT 0.1.1 at 670dc63 with
`--reference inverse --witness`; the independent verifier of E-A (`run_parity.py`: `witness_holds`, `affine_parts`), importing none of
the tools.

## Hypotheses (fixed now)

- **H-F1.** sboxU decides every one of the 8,192 tables affine-equivalent to the inverse.
- **H-F2.** Every sboxU witness holds on all 256 inputs under the independent verifier and is accepted by EASAT as
  `PROVED_AFFINE_WITNESS`.
- **H-F3.** PEIGEN decides all 64 tables of a sample fixed here affine-equivalent to the inverse, and its witnesses pass the same two
  checks. Sample: 32 tables per arm, drawn with Python's `random.Random(20261024).sample` from the 4,096 (key index, round) pairs of
  each arm, HCHM+SHAKE first.

## Limits and acceptance

- Per case: 600 s for sboxU (enforced in the worker when control returns to Python) and 2,400 s for PEIGEN (process killed). A case
  that exceeds its limit is recorded as undecided and counts against the hypothesis.
- Accept each hypothesis only at 100 %. Every exception is reported by case; none is removed.
- If H-F1 holds, the paper may say that a complete decision procedure confirmed the affine equivalence of all 8,192 family tables to the
  inverse function in this environment. Run times are reported descriptively only.
