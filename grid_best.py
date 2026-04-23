#!/usr/bin/env python3
"""
grid_best.py: for each T-P profile in the grid, select the best available
chemistry run (original or any retry pass) and write the winner's three
output files into the profile's chemistry/ directory as:

    monitor_best.dat
    chem_best.dat
    condensates_best.dat

"Best" means fewest grid points that failed on *either* convergence or element
conservation.  Ties are broken by candidate order: original → pass A → B → C
→ D → E (i.e. prefer the earliest / simplest solution).

If a profile has no retry directory (or an empty one), the original files are
the only candidates and are simply copied as *_best.dat.

Existing *_best.dat files are always overwritten.

Usage:
    python grid_best.py [--output-root PATH]
"""

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "input/chemistry"))
from monitor_utils import parse_monitor

# ── grid dimensions (must match grid_retry.py) ───────────────────────────────
TEMPERATURES = list(range(200, 1800, 200))   # 200, 400, ..., 1600 K

COMPOSITIONS = [
    "Bath(H4)",
    "EH3",
    "EL3",
    "Krymka(LL3)",
    "R3",
    "asplund_m3"
]

PASS_LABELS = ["A", "B", "C", "D", "E"]


# ── helpers ───────────────────────────────────────────────────────────────────

def count_any_failures(monitor_path):
    """Count grid points where converged=fail OR elem_conserved=fail.

    Returns None if the file cannot be parsed.
    """
    try:
        df = parse_monitor(monitor_path)
    except Exception as exc:
        print(f"    WARNING: could not parse {monitor_path}: {exc}")
        return None

    failed = 0
    if "converged" in df.columns:
        failed_mask = ~df["converged"].astype(bool)
    else:
        failed_mask = None

    if "elem_conserved" in df.columns:
        elem_mask = ~df["elem_conserved"].astype(bool)
        failed_mask = elem_mask if failed_mask is None else (failed_mask | elem_mask)

    if failed_mask is not None:
        failed = int(failed_mask.sum())

    return failed


def build_candidates(chem_dir):
    """Return an ordered list of (label, monitor, chem, condensates) tuples.

    Order determines tie-breaking: original first, then passes A–E.
    Only candidates where all three files exist are included.
    """
    candidates = []

    # Original run
    mon  = chem_dir / "monitor.dat"
    chem = chem_dir / "chem.dat"
    cond = chem_dir / "condensates.dat"
    if mon.exists() and chem.exists() and cond.exists():
        candidates.append(("original", mon, chem, cond))

    # Retry passes
    retry_dir = chem_dir / "retry"
    if retry_dir.is_dir():
        for label in PASS_LABELS:
            mon  = retry_dir / f"monitor_pass{label}.dat"
            chem = retry_dir / f"chem_pass{label}.dat"
            cond = retry_dir / f"condensates_pass{label}.dat"
            if mon.exists() and chem.exists() and cond.exists():
                candidates.append((f"pass{label}", mon, chem, cond))

    return candidates


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output-root", default="output/grid_day",
        help="Root of the grid output tree (default: output/grid_day)",
    )
    args = parser.parse_args()

    output_root = Path(args.output_root)

    n_total   = 0
    n_skipped = 0   # no candidates found at all
    n_written = 0

    for comp in COMPOSITIONS:
        for T_eq in TEMPERATURES:
            chem_dir = output_root / f"{T_eq}K_{comp}_day" / "chemistry"

            if not chem_dir.is_dir():
                print(f"  MISSING  {T_eq}K {comp}  (directory not found)")
                n_skipped += 1
                continue

            n_total += 1
            candidates = build_candidates(chem_dir)

            if not candidates:
                print(f"  SKIP     {T_eq}K {comp}  (no complete file sets found)")
                n_skipped += 1
                continue

            # Score each candidate
            scored = []
            for label, mon, chem, cond in candidates:
                n_fail = count_any_failures(mon)
                if n_fail is None:
                    continue   # unreadable — skip
                scored.append((n_fail, label, mon, chem, cond))

            if not scored:
                print(f"  SKIP     {T_eq}K {comp}  (all monitor files unreadable)")
                n_skipped += 1
                continue

            # Pick winner: lowest failure count; tie → earliest in list
            scored.sort(key=lambda x: x[0])   # stable sort preserves order on ties
            best_n, best_label, best_mon, best_chem, best_cond = scored[0]

            # Summarise all candidates for the user
            summary_parts = []
            for n_fail, label, *_ in scored:
                mark = " ✓" if label == best_label else ""
                summary_parts.append(f"{label}={n_fail}{mark}")
            print(f"  {T_eq}K {comp:<16s}  best={best_label}  [{', '.join(summary_parts)}]")

            # Write *_best.dat files into chem_dir (overwrite if present)
            shutil.copy2(best_mon,  chem_dir / "monitor_best.dat")
            shutil.copy2(best_chem, chem_dir / "chem_best.dat")
            shutil.copy2(best_cond, chem_dir / "condensates_best.dat")
            n_written += 1

    print(f"\nDone.  {n_written} profile(s) written, {n_skipped} skipped "
          f"(out of {n_total + n_skipped} total).")


if __name__ == "__main__":
    main()
