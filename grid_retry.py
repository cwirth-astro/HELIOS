#!/usr/bin/env python3
"""
Grid rainout retry: identifies T-P profiles with any convergence or
element-conservation failures in the original monitor.dat, then reruns
each failing profile in full with four alternative FastChem solver settings.

Passes (run in order; stops early if a pass fully fixes the profile):
  A  accuracyChem = 1e-7        (tighten chemistry convergence criterion)
  B  condUseSVD  = True         (SVD solver for condensate system)
  C  increase all iteration limits (~3× for gas-phase solvers, ~3× for cond)
  D  condIterChangeLimit = 50   (relax condensate step-size restriction)
  E  nbSwitchToJoint = 0        (activate joint Newton from the first combined
                                 iteration; bypasses hardcoded stagnation cutoff
                                 at 1500 combined iterations that fires before
                                 the default nb_switch_to_joint = 3000)

All passes use rainout condensation (cr).

Output per failing profile:
  output/grid_day/{T}K_{type}_day/chemistry/retry/
      monitor_passA.dat   chem_passA.dat   condensates_passA.dat
      monitor_passB.dat   chem_passB.dat   condensates_passB.dat
      monitor_passC.dat   chem_passC.dat   condensates_passC.dat
      monitor_passD.dat   chem_passD.dat   condensates_passD.dat
      summary.txt

Existing pass files are skipped (safe restart). If a pass fully resolves all
failures the remaining passes are skipped for that profile.

Usage:
    python grid_retry.py [--output-root PATH] [--fastchem-dir PATH]
                         [--pyfastchem-dir PATH] [--dry-run]
"""

import argparse
import signal
import sys
import numpy as np
from pathlib import Path

# ── import monitor parser from chemistry/ ─────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent / "input/chemistry"))
from monitor_utils import parse_monitor

# ── constants ─────────────────────────────────────────────────────────────────
TEMPERATURES = list(range(200, 1800, 200))   # 200, 400, ..., 1600 K

COMPOSITIONS = [
    "Bath(H4)",
    "EH3",
    "EL3",
    "Krymka(LL3)",
    "R3",
]

FASTCHEM_LOGK_GAS  = "input/logK/logK.dat"
FASTCHEM_LOGK_COND = "input/logK/logK_condensates.dat"

K_B_CGS = 1.380649e-16   # erg/K, CODATA 2018

# accuracyChem used in the original rainout runs (rainout_step_setup.py)
ORIGINAL_CHEM_ACCURACY = 1e-4

# Each pass: label → dict of setParameter() calls.
# Stops early if a pass produces zero failures.
PASSES = [
    ("A", {
        "accuracyChem": 1e-7,                   # tighten chemistry criterion
    }),
    ("B", {
        "condUseSVD": True,                     # SVD solver for condensate system
    }),
    ("C", {                                     # increase all iteration limits
        "nbIterationsChem":       10000,
        "nbIterationsBisection":  10000,
        "nbIterationsNewton":     10000,
        "nbIterationsNelderMead": 10000,
        "nbIterationsCond":       10000,
        "nbIterationsChemCond":   100000,
    }),
    ("D", {
        "condIterChangeLimit": 10.0,            # relax condensate step-size limit
    }),
    ("E", {
        # The equilibrium condensation loop has a hardcoded stagnation check that
        # breaks at combined iteration 1500 (checkpoint_interval=500, three strikes).
        # The joint Newton step (which couples gas+condensate simultaneously and
        # breaks the gas/cond oscillation) only activates when
        # nb_combined_iter >= nb_switch_to_joint (default 3000).
        # Because 3000 > 1500, joint Newton never fires in the default configuration.
        # Setting nbSwitchToJoint=0 forces joint Newton from the very first
        # combined iteration whenever active condensates are present.
        "nbSwitchToJoint": 0,
        "accuracyChem": 1e-7,                   # combine with Pass A's tighter chem criterion
    }),
]

# Default values to restore after each pass (FastChem's compiled-in defaults).
PASS_DEFAULTS = {
    "accuracyChem":           ORIGINAL_CHEM_ACCURACY,
    "condUseSVD":             False,
    "nbIterationsChem":       3000,
    "nbIterationsBisection":  3000,
    "nbIterationsNewton":     3000,
    "nbIterationsNelderMead": 3000,
    "nbIterationsCond":       3000,
    "nbIterationsChemCond":   30000,
    "condIterChangeLimit":    5.0,
    "nbSwitchToJoint":        3000,
}


# ── profile scanning ──────────────────────────────────────────────────────────

def profile_has_failures(monitor_path):
    """Return True if any grid point failed convergence or element conservation."""
    try:
        df = parse_monitor(monitor_path)
    except Exception as exc:
        print(f"    WARNING: could not parse {monitor_path}: {exc}")
        return False

    if "converged" in df.columns and (~df["converged"].astype(bool)).any():
        return True
    if "elem_conserved" in df.columns and (~df["elem_conserved"].astype(bool)).any():
        return True
    return False


def count_failures(monitor_path):
    """Return (n_conv_fail, n_elem_fail) from a monitor file."""
    df = parse_monitor(monitor_path)
    n_conv = int((~df["converged"].astype(bool)).sum())      if "converged"      in df.columns else 0
    n_elem = int((~df["elem_conserved"].astype(bool)).sum()) if "elem_conserved" in df.columns else 0
    return n_conv, n_elem


# ── config / tp parsing ───────────────────────────────────────────────────────

def parse_config(config_path):
    """Extract tp_file, run_type, and abundance_file from a FastChem config.dat.

    Non-comment, non-empty values appear in this order (see rainout_step_setup.py):
      0: p-T file path
      1: run type  (g / ce / cr)
      2: chem output paths (space-separated pair on one line)
      3: monitor output path
      4: verbose level
      5: output type  (MR / ND)
      6: element abundance file
      7: logK species data files
      8: chem_accuracy
      9: element_conserve_accuracy
     10: nb_max_fastchem_iter
     11: nb_max_internal_iter
    """
    values = []
    with open(config_path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "console" in stripped.lower() and "verbose" in stripped.lower():
                continue
            values.append(stripped)

    return {
        "tp_file":        values[0],
        "run_type":       values[1],
        "abundance_file": values[6],
    }


def resolve_path(path_str, fastchem_dir, subdir_hint):
    """Return a Path, falling back to searching under fastchem_dir/input/ by hint.

    subdir_hint: the directory name to look for in the original path (e.g.
    'element_abundances') so we can reconstruct the relative portion.
    """
    p = Path(path_str)
    if p.exists():
        return p

    parts = p.parts
    try:
        idx = next(i for i, part in enumerate(parts) if part == subdir_hint)
        fallback = fastchem_dir / "input" / Path(*parts[idx:])
        if fallback.exists():
            return fallback
    except StopIteration:
        pass

    raise FileNotFoundError(
        f"File not found at '{path_str}'\n"
        f"Also tried: {fastchem_dir}/input/{subdir_hint}/..."
    )


def read_tp(tp_path):
    """Read a tp.dat (one header row, then P[bar] T[K] columns)."""
    data = np.loadtxt(tp_path, skiprows=1)
    P = data[:, 0]  # bar
    T = data[:, 1]  # K
    return T, P


# ── FastChem helpers ──────────────────────────────────────────────────────────

def build_fastchem(abundance_file, fastchem_dir):
    import pyfastchem
    fc = pyfastchem.FastChem(
        str(abundance_file),
        str(fastchem_dir / FASTCHEM_LOGK_GAS),
        str(fastchem_dir / FASTCHEM_LOGK_COND),
        1,
    )
    # The compiled default for accuracyChem is 1e-5; the original rainout runs
    # used 1e-4, so set that here so Pass A's tightening is a fair comparison.
    fc.setParameter("accuracyChem", ORIGINAL_CHEM_ACCURACY)
    return fc


def run_profile(fastchem, T, P, run_type):
    import pyfastchem
    inp = pyfastchem.FastChemInput()
    inp.temperature           = T.tolist()
    inp.pressure              = P.tolist()
    inp.equilibrium_condensation = (run_type == "ce")
    inp.rainout_condensation     = (run_type == "cr")
    out = pyfastchem.FastChemOutput()
    fastchem.calcDensities(inp, out)
    return out


def output_failure_counts(output):
    """Return (n_conv_fail, n_elem_fail) from a pyfastchem output object."""
    import pyfastchem
    ec = np.array(output.element_conserved)
    n_conv = sum(1 for f in output.fastchem_flag if f != pyfastchem.FASTCHEM_SUCCESS)
    n_elem = int((~np.all(ec == 1, axis=1)).sum())
    return n_conv, n_elem


# ── output writers ────────────────────────────────────────────────────────────

def write_monitor(path, T, P, output, fastchem):
    import pyfastchem
    ec      = np.array(output.element_conserved)
    flags   = list(output.fastchem_flag)
    ok_fail = {True: "ok", False: "fail"}

    with open(path, "w") as f:
        header = [
            "#grid_point", "iterations", "chem_iter", "cond_iter",
            "converged", "elem_conserved", "p(bar)", "T(K)",
            "n_<tot>(cm-3)", "n_g(cm-3)", "m(u)",
        ]
        f.write("\t".join(header))
        for j in range(fastchem.getElementNumber()):
            f.write("\t" + fastchem.getElementSymbol(j))
        f.write("\n")

        for k in range(len(T)):
            conv   = ok_fail[flags[k] == pyfastchem.FASTCHEM_SUCCESS]
            all_ok = ok_fail[bool(np.all(ec[k]))]
            n_g    = P[k] * 1e6 / (K_B_CGS * T[k])
            row = [
                str(k),
                str(output.nb_iterations[k]),
                str(output.nb_chemistry_iterations[k]),
                str(output.nb_cond_iterations[k]),
                conv, all_ok,
                f"{P[k]:.10e}", f"{T[k]:.10e}",
                f"{output.total_element_density[k]:.10e}",
                f"{n_g:.10e}",
                f"{output.mean_molecular_weight[k]:.10e}",
            ]
            f.write("\t".join(row))
            for j in range(fastchem.getElementNumber()):
                f.write("\t" + ok_fail[bool(ec[k, j])])
            f.write("\n")


def write_chem(path, T, P, output, fastchem):
    nd      = np.array(output.number_densities)
    n_g     = P * 1e6 / (K_B_CGS * T)
    mr      = nd / n_g[:, None]
    nb_spec = fastchem.getGasSpeciesNumber()

    with open(path, "w") as f:
        header = ["#grid_point", "p(bar)", "T(K)", "n_<tot>(cm-3)", "n_g(cm-3)", "m(u)"]
        header += [fastchem.getGasSpeciesSymbol(j) for j in range(nb_spec)]
        f.write("\t".join(header) + "\n")

        for k in range(len(T)):
            row = [
                str(k),
                f"{P[k]:.10e}", f"{T[k]:.10e}",
                f"{output.total_element_density[k]:.10e}",
                f"{n_g[k]:.10e}",
                f"{output.mean_molecular_weight[k]:.10e}",
            ]
            row += [f"{mr[k, j]:.10e}" for j in range(nb_spec)]
            f.write("\t".join(row) + "\n")


def write_cond(path, T, P, output, fastchem):
    ec_deg  = np.array(output.element_cond_degree)
    nd_cond = np.array(output.number_densities_cond)
    nb_elem = fastchem.getElementNumber()
    nb_cond = fastchem.getCondSpeciesNumber()

    with open(path, "w") as f:
        header = ["#grid_point", "p(bar)", "T(K)"]
        header += [fastchem.getElementSymbol(j)     for j in range(nb_elem)]
        header += [fastchem.getCondSpeciesSymbol(j) for j in range(nb_cond)]
        f.write("\t".join(header) + "\n")

        for k in range(len(T)):
            row = [str(k), f"{P[k]:.10e}", f"{T[k]:.10e}"]
            row += [f"{ec_deg[k, j]:.10e}"  for j in range(nb_elem)]
            row += [f"{nd_cond[k, j]:.10e}" for j in range(nb_cond)]
            f.write("\t".join(row) + "\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    signal.signal(signal.SIGINT, signal.default_int_handler)

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output-root", default="output/grid_day",
        help="Root of the grid output tree (default: output/grid_day)",
    )
    parser.add_argument(
        "--fastchem-dir", default=None,
        help="Path to FastChem repo root (default: FastChem/ relative to this script)",
    )
    parser.add_argument(
        "--pyfastchem-dir", default=None,
        help="Directory containing pyfastchem .so (default: same as --fastchem-dir)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print failing profiles and exit without running anything",
    )
    args = parser.parse_args()

    output_root  = Path(args.output_root)
    script_dir   = Path(__file__).parent.resolve()
    fastchem_dir = (
        Path(args.fastchem_dir).resolve() if args.fastchem_dir
        else script_dir / "FastChem"
    )
    pyfastchem_dir = (
        Path(args.pyfastchem_dir).resolve() if args.pyfastchem_dir
        else fastchem_dir
    )
    sys.path.insert(0, str(pyfastchem_dir))

    # ── verify pyfastchem is importable before doing any real work ─────────
    import glob as _glob
    so_files = _glob.glob(str(pyfastchem_dir / "pyfastchem*.so"))
    print(f"Python:        {sys.executable}  ({sys.version.split()[0]})")
    print(f"pyfastchem dir: {pyfastchem_dir}")
    print(f"  .so files found: {so_files if so_files else 'NONE'}")
    try:
        import pyfastchem as _pfc  # noqa: F401
        print(f"  import OK: {_pfc.__file__}")
    except ModuleNotFoundError as exc:
        print(f"\nERROR: cannot import pyfastchem: {exc}", file=sys.stderr)
        print(
            "Hints:\n"
            "  1. Check that the .so ABI tag matches your Python version\n"
            "     (e.g. cpython-311 requires python3.11, not python3.10/3.12)\n"
            "  2. Check LD_LIBRARY_PATH includes the FastChem build dir:\n"
            f"     export LD_LIBRARY_PATH={fastchem_dir / 'build'}:$LD_LIBRARY_PATH\n"
            "  3. Confirm the .so file is in --pyfastchem-dir",
            file=sys.stderr,
        )
        sys.exit(1)
    print()

    # ── scan for failing profiles ──────────────────────────────────────────
    print("Scanning profiles for failures...")
    failing = []
    for comp in COMPOSITIONS:
        for T_eq in TEMPERATURES:
            chem_dir = output_root / f"{T_eq}K_{comp}_day" / "chemistry"
            monitor  = chem_dir / "monitor.dat"
            if not monitor.exists():
                print(f"  MISSING  {T_eq}K {comp}")
                continue
            if profile_has_failures(monitor):
                n_conv, n_elem = count_failures(monitor)
                failing.append((T_eq, comp, chem_dir))
                print(f"  FAILING  {T_eq}K {comp}  "
                      f"(conv_fail={n_conv}, elem_fail={n_elem})")

    print(f"\n{len(failing)} failing profile(s) found across "
          f"{len(COMPOSITIONS) * len(TEMPERATURES)} total.\n")

    if args.dry_run or not failing:
        return

    # ── retry each failing profile ─────────────────────────────────────────
    for T_eq, comp, chem_dir in failing:
        retry_dir = chem_dir / "retry"
        retry_dir.mkdir(exist_ok=True)

        print(f"{'='*62}")
        print(f"{T_eq}K  {comp}")

        # Parse config for this profile
        config_path = chem_dir / "config.dat"
        try:
            config = parse_config(config_path)
        except Exception as exc:
            print(f"  ERROR reading config: {exc} — skipping")
            continue

        run_type = config["run_type"]

        # Resolve T-P profile
        tp_path = Path(config["tp_file"])
        if not tp_path.exists():
            tp_path = chem_dir / "tp.dat"
        if not tp_path.exists():
            print(f"  ERROR: tp.dat not found at {tp_path} — skipping")
            continue
        T_arr, P_arr = read_tp(tp_path)
        print(f"  T-P profile: {len(T_arr)} grid points from {tp_path}")

        # Resolve element abundance file
        try:
            abundance_file = resolve_path(
                config["abundance_file"], fastchem_dir, "element_abundances"
            )
        except FileNotFoundError as exc:
            print(f"  ERROR: {exc} — skipping")
            continue
        print(f"  Abundances:  {abundance_file}")
        print(f"  Run type:    {run_type}")
        print(f"  Output:      {retry_dir}/\n")

        fastchem = build_fastchem(abundance_file, fastchem_dir)

        summary_lines = [
            f"{T_eq}K  {comp}",
            f"T-P profile: {len(T_arr)} grid points",
            f"",
        ]

        solved = False
        for pass_label, params in PASSES:
            if solved:
                break

            mon_out  = retry_dir / f"monitor_pass{pass_label}.dat"
            chem_out = retry_dir / f"chem_pass{pass_label}.dat"
            cond_out = retry_dir / f"condensates_pass{pass_label}.dat"

            param_str = ", ".join(f"{k}={v}" for k, v in params.items())

            if mon_out.exists() and chem_out.exists() and cond_out.exists():
                n_conv, n_elem = count_failures(mon_out)
                suffix = "  (existing)"
                if n_conv == 0 and n_elem == 0:
                    suffix += "  ✓ SOLVED"
                    solved = True
                line = (f"Pass {pass_label} [{param_str}]  "
                        f"conv_fail={n_conv}  elem_fail={n_elem}{suffix}")
                print(f"  {line}")
                summary_lines.append(line)
                continue

            print(f"  Pass {pass_label}: {param_str}")
            # The dev_log pybind11 wrapper uses .noconvert() on all overloads,
            # so Python bool routes to the bool overload and int to unsigned int.
            for k, v in params.items():
                fastchem.setParameter(k, v)

            output = run_profile(fastchem, T_arr, P_arr, run_type)

            # Reset to original-run defaults before next pass
            for k in params:
                fastchem.setParameter(k, PASS_DEFAULTS[k])

            n_conv, n_elem = output_failure_counts(output)
            suffix = ""
            if n_conv == 0 and n_elem == 0:
                suffix = "  ✓ SOLVED"
                solved = True
            line = (f"Pass {pass_label} [{param_str}]  "
                    f"conv_fail={n_conv}  elem_fail={n_elem}{suffix}")
            print(f"    → {line}")
            summary_lines.append(line)

            write_monitor(mon_out,  T_arr, P_arr, output, fastchem)
            write_chem(   chem_out, T_arr, P_arr, output, fastchem)
            write_cond(   cond_out, T_arr, P_arr, output, fastchem)

        if solved:
            summary_lines.append("(remaining passes skipped — profile fully converged)")
            print("  Profile fully converged — skipping remaining passes.")

        # summary.txt (always overwritten with up-to-date results)
        with open(retry_dir / "summary.txt", "w") as f:
            f.write("\n".join(summary_lines) + "\n")

        print()

    print("Done.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        sys.exit(130)
