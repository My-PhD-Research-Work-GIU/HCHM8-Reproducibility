"""Cost of three S-box derivations on one machine (plan: ../PLAN.md, fixed before the registered runs).

Arms (one table = candidate 0 of the paper's derivation for one key and round; the arms differ only in the key material that
enters SHAKE256 before the affine layers are built):

  S  SHAKE-only          base = key | round | "|HAS1-FAMILY-V2|"
  H  HCHM+SHAKE          base = key | round | "|HAS1-FAMILY-V2|" | serialize_hchm(key, round)          (Float64 orbit)
  L  lattice+SHAKE       base = key | round | "|HAS1-FAMILY-V3-LATTICE16|" | lattice_serialize(key, round) (16-bit lattice)

then stream = SHAKE256(base | counter 0 | "|AFFINE-PROFILE-V2|", 258 bytes), make_rows twice, build_sbox. The functions are
those of the paper (05_CODE/run_validated_experiments.py; lattice_serialize taken verbatim from Part D of
05_CODE/platform_rerun_and_extended_ablation.py, as in the chaos-only test). Every timed table is checked, after the timing,
against the reference function of the paper.

Each run is a fresh Python process that measures one arm: imports, a warm-up of 4 tables on keys outside the measured set,
then the timed derivation of every table of the measured keys (16 rounds each) with time.perf_counter_ns around each stage:
  material  parameters, orbit and serialization (empty for S);
  expand    SHAKE256 expansion to 258 bytes;
  table     make_rows twice and build_sbox.
Memory: the process's peak working set (Windows) or peak resident set (Linux) after imports, after the warm-up and after the
timed loop, and, in a separate untimed pass over 16 tables, the peak of Python allocations traced by tracemalloc.

Usage:
  python derivation_cost.py run <out_dir> [--reps 20] [--keys 0:16] [--smoke]   driver: runs every worker, writes the summary
  python derivation_cost.py worker <code_dir> <arm> <k0> <k1> <out.json>          one measured process (called by the driver)
"""
import hashlib, importlib, itertools, json, math, os, platform, shutil, subprocess, sys, tempfile, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                      # 05_PAPERS_AND_DATASETS/05_PAPERS_AND_DATASETS
KEY_OFFSET = 60000                          # measured keys: key_for_index(60000 + k); not used by any earlier test
WARM_OFFSET = 70000                         # warm-up keys: key_for_index(70000 + i)
ROUNDS = 16
ARMS = ("S", "H", "L")
BOOT_B, BOOT_SEED = 10000, 20261110


# ------------------------------------------------------------------------------------------------ memory of this process
def memory_bytes():
    """(current, peak) resident memory of this process in bytes."""
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t), ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPagedPoolUsage", ctypes.c_size_t), ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaNonPagedPoolUsage", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t),
                        ("PeakPagefileUsage", ctypes.c_size_t)]
        c = PMC(); c.cb = ctypes.sizeof(PMC)
        k32 = ctypes.WinDLL("kernel32"); psapi = ctypes.WinDLL("psapi")
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
        assert psapi.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(c), c.cb)
        return int(c.WorkingSetSize), int(c.PeakWorkingSetSize)
    cur = peak = 0
    with open("/proc/self/status") as f:
        for line in f:
            if line.startswith("VmRSS:"): cur = int(line.split()[1]) * 1024
            if line.startswith("VmHWM:"): peak = int(line.split()[1]) * 1024
    return cur, peak


def machine():
    info = {"platform": platform.platform(), "machine": platform.machine(), "python": platform.python_version(),
            "logical_cpus": os.cpu_count()}
    if sys.platform == "win32":
        import ctypes, winreg
        try:
            k = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r"HARDWARE\DESCRIPTION\System\CentralProcessor\0")
            info["cpu"] = winreg.QueryValueEx(k, "ProcessorNameString")[0].strip()
            info["cpu_nominal_mhz"] = winreg.QueryValueEx(k, "~MHz")[0]
        except OSError:
            pass

        class MS(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong), ("ullTotalPhys", ctypes.c_ulonglong),
                        ("ullAvailPhys", ctypes.c_ulonglong), ("ullTotalPageFile", ctypes.c_ulonglong),
                        ("ullAvailPageFile", ctypes.c_ulonglong), ("ullTotalVirtual", ctypes.c_ulonglong),
                        ("ullAvailVirtual", ctypes.c_ulonglong), ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
        ms = MS(); ms.dwLength = ctypes.sizeof(MS); ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms))
        info["ram_bytes"] = int(ms.ullTotalPhys); info["memory_load_percent_at_start"] = int(ms.dwMemoryLoad)

        class SPS(ctypes.Structure):
            _fields_ = [("ACLineStatus", ctypes.c_ubyte), ("BatteryFlag", ctypes.c_ubyte), ("BatteryLifePercent", ctypes.c_ubyte),
                        ("SystemStatusFlag", ctypes.c_ubyte), ("BatteryLifeTime", ctypes.c_ulong), ("BatteryFullLifeTime", ctypes.c_ulong)]
        sps = SPS(); ctypes.windll.kernel32.GetSystemPowerStatus(ctypes.byref(sps))
        info["power_source"] = {0: "battery", 1: "mains", 255: "unknown"}.get(sps.ACLineStatus, "unknown")
        info["battery_percent"] = None if sps.BatteryLifePercent == 255 else int(sps.BatteryLifePercent)
        info["battery_saver_on"] = bool(sps.SystemStatusFlag & 1)
    else:
        try:
            for line in open("/proc/cpuinfo"):
                if line.startswith("model name"): info["cpu"] = line.split(":", 1)[1].strip(); break
            for line in open("/proc/meminfo"):
                if line.startswith("MemTotal:"): info["ram_bytes"] = int(line.split()[1]) * 1024; break
        except OSError:
            pass
    return info


# ------------------------------------------------------------------------------------------------ one measured process
def worker(code_dir, arm, k0, k1, out_json):
    t_start = time.perf_counter_ns()
    sys.path.insert(0, code_dir)
    import numpy as np
    import matplotlib; matplotlib.use("Agg")
    m = importlib.import_module("run_validated_experiments")
    src = (Path(code_dir) / "_platform_source.py").read_text(encoding="utf-8")
    i0 = src.index("import math\nLATTICE_Q"); i1 = src.index("t0 = time.time(); lat_rows")
    ns = {"m": m}; exec(compile(src[i0:i1], "lattice_serialize", "exec"), ns); lattice_serialize = ns["lattice_serialize"]
    import_ns = time.perf_counter_ns() - t_start
    mem_import = memory_bytes()

    tail = (0).to_bytes(2, "big") + b"|AFFINE-PROFILE-V2|"

    def derive(key, r, clock):
        t0 = clock()
        if arm == "S":
            base = key + r.to_bytes(4, "big") + b"|HAS1-FAMILY-V2|"
        elif arm == "H":
            base = key + r.to_bytes(4, "big") + b"|HAS1-FAMILY-V2|" + m.serialize_hchm(key, r)
        else:
            base = key + r.to_bytes(4, "big") + b"|HAS1-FAMILY-V3-LATTICE16|" + lattice_serialize(key, r)
        t1 = clock()
        st = m.shake(base + tail, 258)
        t2 = clock()
        rin, idx = m.make_rows(st, 0); rout, idx = m.make_rows(st, idx)
        S = m.build_sbox(rin, rout, st[idx % len(st)], st[(idx + 1) % len(st)])
        t3 = clock()
        return S, base, (t1 - t0, t2 - t1, t3 - t2)

    clock = time.perf_counter_ns
    for i in range(4):                                   # warm-up, keys outside the measured set
        derive(m.key_for_index(WARM_OFFSET + i), i, clock)
    mem_warm = memory_bytes()

    keys = [m.key_for_index(KEY_OFFSET + k) for k in range(k0, k1)]
    tables, bases, rows = [], [], []
    loop0 = clock()
    for key in keys:
        for r in range(ROUNDS):
            S, base, t = derive(key, r, clock)
            tables.append(S); bases.append(base); rows.append(t)
    loop_ns = clock() - loop0
    mem_loop = memory_bytes()

    # untimed check against the reference function of the paper
    ok = all(np.array_equal(S, m.candidate(base, 0)[0]) for S, base in zip(tables, bases))
    digest = hashlib.sha256(b"".join(S.tobytes() for S in tables)).hexdigest()

    import tracemalloc                                   # separate, untimed pass: Python allocations of 16 tables
    tracemalloc.start()
    for r in range(ROUNDS):
        derive(keys[0], r, clock)
    traced_peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()

    res = {"arm": arm, "keys": [k0, k1], "tables": len(tables), "reference_check_passed": bool(ok), "tables_sha256": digest,
           "import_ns": import_ns, "loop_ns": loop_ns, "per_table_ns": [list(t) for t in rows],
           "memory": {"after_import": mem_import, "after_warmup": mem_warm, "after_loop": mem_loop, "tracemalloc_peak_16_tables": traced_peak},
           "pid": os.getpid()}
    Path(out_json).write_text(json.dumps(res))


# ------------------------------------------------------------------------------------------------ statistics
def boot_ci(values, stat, rng):
    import numpy as np
    v = np.asarray(values, float); idx = rng.integers(0, len(v), size=(BOOT_B, len(v)))
    b = np.array([stat(v[i]) for i in idx]); return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


def boot_ratio(x, y, rng):
    """95% percentile interval of median(x) / median(y), runs resampled within each arm."""
    import numpy as np
    x = np.asarray(x, float); y = np.asarray(y, float)
    ix = rng.integers(0, len(x), size=(BOOT_B, len(x))); iy = rng.integers(0, len(y), size=(BOOT_B, len(y)))
    b = np.median(x[ix], axis=1) / np.median(y[iy], axis=1)
    return [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))]


def summarize(runs):
    import numpy as np
    from scipy import stats
    rng = np.random.default_rng(BOOT_SEED); out = {"arms": {}, "ratios": {}}
    per_arm_med = {}
    for arm in ARMS:
        rs = [r for r in runs if r["arm"] == arm]
        thr = [r["tables"] / (r["loop_ns"] / 1e9) for r in rs]
        med_ms = [float(np.median([sum(t) for t in r["per_table_ns"]])) / 1e6 for r in rs]
        stage = {name: [float(np.median([t[i] for t in r["per_table_ns"]])) / 1e6 for r in rs] for i, name in enumerate(("material", "expand", "table"))}
        pooled = np.array([sum(t) for r in rs for t in r["per_table_ns"]], float) / 1e6
        mem = {k: [r["memory"][k][1] / 2**20 for r in rs] for k in ("after_import", "after_warmup", "after_loop")}
        traced = [r["memory"]["tracemalloc_peak_16_tables"] / 2**10 for r in rs]
        order = [r["sequence"] for r in rs]
        out["arms"][arm] = {
            "runs": len(rs), "tables_per_run": rs[0]["tables"], "reference_check_passed_all": all(r["reference_check_passed"] for r in rs),
            "tables_identical_across_runs": len({r["tables_sha256"] for r in rs}) == 1,
            "throughput_tables_per_s": {"median": float(np.median(thr)), "ci95_median": boot_ci(thr, np.median, rng),
                                        "mean": float(np.mean(thr)), "sd": float(np.std(thr, ddof=1)), "min": float(min(thr)), "max": float(max(thr))},
            "per_table_ms": {"median_of_run_medians": float(np.median(med_ms)), "ci95": boot_ci(med_ms, np.median, rng),
                             "pooled_p05": float(np.percentile(pooled, 5)), "pooled_p95": float(np.percentile(pooled, 95)),
                             "pooled_max": float(pooled.max())},
            "stage_ms_median_of_run_medians": {k: float(np.median(v)) for k, v in stage.items()},
            "stage_share_of_total": {k: float(np.median(v)) / float(np.median(med_ms)) for k, v in stage.items()},
            "peak_memory_mib": {k: {"median": float(np.median(v)), "min": float(min(v)), "max": float(max(v))} for k, v in mem.items()},
            "peak_increase_during_loop_mib": {"median": float(np.median(np.array(mem["after_loop"]) - np.array(mem["after_warmup"]))),
                                              "max": float(max(np.array(mem["after_loop"]) - np.array(mem["after_warmup"])))},
            "tracemalloc_peak_kib_16_tables": {"median": float(np.median(traced)), "min": float(min(traced)), "max": float(max(traced))},
            "import_s_median": float(np.median([r["import_ns"] for r in rs])) / 1e9,
            "drift_spearman_throughput_vs_sequence": float(stats.spearmanr(order, thr).statistic) if len(rs) > 2 else None}
        per_arm_med[arm] = med_ms
    for x, y in (("H", "S"), ("L", "S"), ("L", "H")):
        out["ratios"]["%s/%s" % (x, y)] = {"median_per_table_time_ratio": float(np.median(per_arm_med[x]) / np.median(per_arm_med[y])),
                                          "ci95": boot_ratio(per_arm_med[x], per_arm_med[y], rng)}
    return out


# ------------------------------------------------------------------------------------------------ driver
def run(out_dir, reps, keys, smoke):
    import numpy as np, scipy
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True); k0, k1 = (int(v) for v in keys.split(":"))
    code = Path(tempfile.mkdtemp(prefix="r01_cost_"))
    shutil.copy(ROOT / "05_CODE" / "run_validated_experiments.py", code / "run_validated_experiments.py")
    shutil.copy(ROOT / "05_CODE" / "platform_rerun_and_extended_ablation.py", code / "_platform_source.py")
    env = dict(os.environ); env["NUMBA_CACHE_DIR"] = str(code / "numba_cache"); env["PYTHONIOENCODING"] = "utf-8"
    res = {"plan": "../PLAN.md", "smoke_test": smoke, "machine": machine(),
           "packages": {"numpy": np.__version__, "scipy": scipy.__version__},
           "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "reference_sources_sha256": {p: hashlib.sha256((ROOT / "05_CODE" / p).read_bytes()).hexdigest()
                                        for p in ("run_validated_experiments.py", "platform_rerun_and_extended_ablation.py")},
           "design": {"repetitions_per_arm": reps, "keys": [KEY_OFFSET + k0, KEY_OFFSET + k1], "rounds": ROUNDS,
                      "tables_per_run": (k1 - k0) * ROUNDS, "order": "the six permutations of (S, H, L) in turn, one per repetition",
                      "warmup_tables": 4, "bootstrap": {"resamples": BOOT_B, "seed": BOOT_SEED}},
           "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    perms = list(itertools.permutations(ARMS)); runs = []; seq = 0
    for rep in range(reps):
        for arm in perms[rep % len(perms)]:
            f = code / ("run_%03d_%s.json" % (seq, arm))
            subprocess.run([sys.executable, "-X", "utf8", str(Path(__file__).resolve()), "worker", str(code), arm, str(k0), str(k1), str(f)],
                           check=True, env=env)
            r = json.loads(f.read_text()); r["sequence"] = seq; r["repetition"] = rep; runs.append(r); seq += 1
            print("run %3d  %s  %.1f tables/s  peak %.0f MiB" % (seq, arm, r["tables"] / (r["loop_ns"] / 1e9), r["memory"]["after_loop"][1] / 2**20), flush=True)
    res["machine_at_end"] = {k: v for k, v in machine().items() if k in ("power_source", "battery_percent", "battery_saver_on", "memory_load_percent_at_start")}
    res["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    res["summary"] = summarize(runs)
    sfx = "_smoke" if smoke else ""
    (out / ("RESULTS_derivation_cost%s.json" % sfx)).write_text(json.dumps(res, indent=1))
    with open(out / ("runs_per_table_ns%s.csv" % sfx), "w", newline="") as fh:
        fh.write("sequence,repetition,arm,table,material_ns,expand_ns,table_ns\n")
        for r in runs:
            for i, t in enumerate(r["per_table_ns"]):
                fh.write("%d,%d,%s,%d,%d,%d,%d\n" % (r["sequence"], r["repetition"], r["arm"], i, *t))
    with open(out / ("runs%s.csv" % sfx), "w", newline="") as fh:
        fh.write("sequence,repetition,arm,tables,loop_ns,import_ns,peak_after_import,peak_after_warmup,peak_after_loop,tracemalloc_peak_16,reference_check,tables_sha256\n")
        for r in runs:
            mm = r["memory"]
            fh.write("%d,%d,%s,%d,%d,%d,%d,%d,%d,%d,%s,%s\n" % (r["sequence"], r["repetition"], r["arm"], r["tables"], r["loop_ns"], r["import_ns"],
                     mm["after_import"][1], mm["after_warmup"][1], mm["after_loop"][1], mm["tracemalloc_peak_16_tables"], r["reference_check_passed"], r["tables_sha256"]))
    shutil.rmtree(code, ignore_errors=True)
    s = res["summary"]
    for arm in ARMS:
        a = s["arms"][arm]
        print("%s: %.1f tables/s (95%% CI %.1f to %.1f), %.3f ms per table, peak %.0f MiB, traced %.0f KiB" % (
            arm, a["throughput_tables_per_s"]["median"], *a["throughput_tables_per_s"]["ci95_median"], a["per_table_ms"]["median_of_run_medians"],
            a["peak_memory_mib"]["after_loop"]["median"], a["tracemalloc_peak_kib_16_tables"]["median"]))
    for k, v in s["ratios"].items():
        print("%s time ratio %.2f (95%% CI %.2f to %.2f)" % (k, v["median_per_table_time_ratio"], *v["ci95"]))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        worker(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6])
    elif len(sys.argv) > 2 and sys.argv[1] == "run":
        import argparse
        ap = argparse.ArgumentParser(); ap.add_argument("mode"); ap.add_argument("out"); ap.add_argument("--reps", type=int, default=20)
        ap.add_argument("--keys", default="0:16"); ap.add_argument("--smoke", action="store_true")
        a = ap.parse_args(); run(a.out, a.reps, a.keys, a.smoke)
    else:
        print(__doc__); sys.exit(2)
