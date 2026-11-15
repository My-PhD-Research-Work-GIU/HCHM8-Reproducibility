"""Compares a rerun of the two tests of paper revision 9 with the stored results (plan: ../PLAN.md).

  python compare_rerun.py <grid_rerun_dir> <chaos_rerun_dir> <out.json> [--environment-label text]

<grid_rerun_dir> holds RESULTS_grid_conjugacy.json and <chaos_rerun_dir> holds RESULTS_chaos_only.json and
chaos_only_metrics.csv, written by the unchanged scripts with the commands of the plans. The stored files are those in
08_FUTURE_WORK/2026-11-10_GRID_CONJUGACY and 08_FUTURE_WORK/2026-11-10_CHAOS_ONLY_ARM.

Two kinds of result are separated:
  must hold in every environment (the exit code is 1 otherwise): G_Q is a permutation with the cycle type of A_Q for every
  profile (Proposition 2); every G_Q orbit period equals the A_Q period at H_Q(x0) and divides the core period; the S and H
  arms reproduce the paper's derivation; the SHAKE-only tables equal the stored ones (they involve no floating point);
  reported, because they may differ between environments: whether the registered criteria of both plans are met again, the
  statistics of F_Q, the tables of the lattice arm (portable in the tests made, not in principle) and of the arms that pass
  through a Float64 orbit (H and C), and the numbers of the four comparisons.
"""
import csv, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FW = HERE.parents[2] / "08_FUTURE_WORK"
GRID = FW / "2026-11-10_GRID_CONJUGACY"; CHAOS = FW / "2026-11-10_CHAOS_ONLY_ARM"


def rows_csv(p):
    with open(p) as f: return {(r["arm"], int(r["key_index"]), int(r["round"])): r for r in csv.DictReader(f)}


def main():
    gdir, cdir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    label = sys.argv[sys.argv.index("--environment-label") + 1] if "--environment-label" in sys.argv else ""
    gs = json.loads((GRID / "RESULTS_grid_conjugacy.json").read_text()); gr = json.loads((gdir / "RESULTS_grid_conjugacy.json").read_text())
    cs = json.loads((CHAOS / "RESULTS_chaos_only.json").read_text()); cr = json.loads((cdir / "RESULTS_chaos_only.json").read_text())
    res = {"environment_label": label, "rerun_environment": {"grid": gr["environment"], "chaos": cr["environment"]},
           "stored_environment": {"grid": gs["environment"], "chaos": cs["environment"]},
           "scripts_identical": {"grid": gr["script_sha256"] == gs["script_sha256"] and gr["emulator_source_sha256"] == gs["emulator_source_sha256"],
                                 "chaos": cr["script_sha256"] == cs["script_sha256"]}}
    must = {}

    # ---- grid-exact conjugate
    ex_s = {(r["qbits"], r["profile"]): r for r in gs["exhaustive"] if "profile" in r}
    ex_r = {(r["qbits"], r["profile"]): r for r in gr["exhaustive"] if "profile" in r}
    must["G_permutation_with_type_of_A"] = len(ex_r) == 20 and all(r["G"]["is_permutation"] and r["G"]["cycle_type_equals_A"] for r in ex_r.values())
    fields = ("image_fraction", "cyclic_states", "cycles", "mean_rho", "max_tail", "longest_cycle")
    f_same = sum(all(ex_r[k]["F"][f] == ex_s[k]["F"][f] for f in fields) and ex_r[k]["F_equals_G_fraction"] == ex_s[k]["F_equals_G_fraction"] for k in ex_s)
    g_same = sum(ex_r[k]["G"]["cycle_type"] == ex_s[k]["G"]["cycle_type"] for k in ex_s)
    orb_s = {(r["qbits"], r["profile"], r["run"]): r for r in gs["orbits"]}; orb_r = {(r["qbits"], r["profile"], r["run"]): r for r in gr["orbits"]}
    must["G_orbit_periods_equal_A_and_divide_core"] = len(orb_r) == 790 and all(r["equal"] and r["divides_core_period"] for r in orb_r.values())
    p_same = sum(orb_r[k]["period_G"] == orb_s[k]["period_G"] and orb_r[k]["period_A_of_HQx0"] == orb_s[k]["period_A_of_HQx0"] for k in orb_s)
    q3 = [ex_r[(3, p)]["F"] for p in range(10)]
    img = sum(r["image_fraction"] for r in q3) / 10; rho = sum(r["mean_rho_over_sqrt_piN_2"] for r in q3) / 10
    criteria = {"registered_F_criteria_at_Q3": 0.60 <= img <= 0.66 and 0.7 <= rho <= 1.3}
    res["grid"] = {"profiles_with_identical_F_statistics": f_same, "profiles_with_identical_G_cycle_type": g_same, "profiles": len(ex_s),
                   "orbits_with_identical_periods": p_same, "orbits": len(orb_s), "F_Q3_mean_image_fraction": img, "F_Q3_mean_rho_ratio": rho,
                   "periods_by_q": {str(q): sorted({r["period_G"] for r in orb_r.values() if r["qbits"] == q}) for q in (4, 5, 6, 7)}}

    # ---- chaos-only arms
    ts = rows_csv(CHAOS / "DATA" / "chaos_only_metrics.csv"); tr = rows_csv(cdir / "chaos_only_metrics.csv")
    same = {a: sum(tr[k]["sha256"] == ts[k]["sha256"] for k in ts if k[0] == a) for a in ("S", "H", "C", "CL")}
    must["paper_arms_reproduced_by_the_script"] = cr["paper_arms_reproduced_exactly"] == 4096
    must["SHAKE_only_tables_identical_to_stored"] = same["S"] == 4096
    prim = cr["primary"]
    criteria["registered_chaos_criteria"] = all(not c["any_significant_after_holm"] and c["g_interval_bound"] <= 0.10 for c in prim)
    cmp_s = {(c["x"], c["y"]): c for c in cs["primary"] + cs["secondary"] + [cs["paper_comparison_recomputed"]]}
    cmp_r = {(c["x"], c["y"]): c for c in cr["primary"] + cr["secondary"] + [cr["paper_comparison_recomputed"]]}
    res["chaos"] = {"tables_identical_to_stored": same, "tables_per_arm": 4096,
                    "comparisons": {"%s vs %s" % k: {q: {"stored": cmp_s[k][q], "rerun": cmp_r[k][q]}
                                                    for q in ("smallest_p", "smallest_p_holm", "max_abs_g", "g_interval_bound")} for k in cmp_s},
                    "byte_uniformity_rerun": cr["stream_byte_uniformity"], "distinct_tables_rerun": cr["distinct_tables"]}
    res["must_hold"] = must; res["all_must_hold"] = all(must.values()); res["registered_criteria_met_again"] = criteria
    out.parent.mkdir(parents=True, exist_ok=True); out.write_text(json.dumps(res, indent=1))
    print("grid: F statistics identical for %d of %d profiles, G cycle types identical %d, orbit periods identical %d of %d" % (
        f_same, len(ex_s), g_same, p_same, len(orb_s)))
    print("chaos: tables identical to stored per arm", same)
    for k, v in res["chaos"]["comparisons"].items():
        print("  %s: smallest Holm p %.4f -> %.4f, g bound %.4f -> %.4f" % (k, v["smallest_p_holm"]["stored"], v["smallest_p_holm"]["rerun"],
                                                                         v["g_interval_bound"]["stored"], v["g_interval_bound"]["rerun"]))
    for k, v in must.items(): print("must hold: %-55s %s" % (k, v))
    for k, v in criteria.items(): print("criterion met again: %-45s %s" % (k, v))
    sys.exit(0 if res["all_must_hold"] else 1)


if __name__ == "__main__":
    main()
