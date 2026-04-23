#!/usr/bin/env python3
"""
interpolate_fastchem_vmr.py

Reads a HELIOS _tp.dat temperature-pressure profile and a FastChem chem.dat
equilibrium chemistry grid, then interpolates each species VMR onto the
1D T-P profile using exactly the same bilinear scheme as HELIOS:

    scipy.interpolate.RectBivariateSpline(kx=1, ky=1)
    axes: T [K] (linear) and log10(P [bar]) (linear)
    out-of-bounds: clamp to grid edge (no extrapolation)

FastChem outputs its species values in the units of the file as-is (typically
number densities in cm^-3, or mixing ratios depending on FastChem config).
This script passes those values through without unit conversion; the output
column header marks the units as whatever FastChem wrote.

Output: tab-separated file, one row per atmospheric layer, columns:
    P[bar]  T[K]  <species1>  <species2>  ...
where species columns are named by HELIOS name (from species.dat).

Usage:
    python interpolate_fastchem_vmr.py \\
        --tp   path/to/run_tp.dat \\
        --chem path/to/chem.dat \\
        --species path/to/species.dat \\
        --output path/to/output.dat

Arguments:
    --tp        HELIOS _tp.dat output file (columns: layer T[K] P[1e-6 bar] ...)
    --chem      FastChem chem.dat grid file
    --species   HELIOS species.dat file listing species with source_for_vmr=FastChem
    --output    Output file path (default: vmr_interpolated.dat)
"""

import argparse
import sys
import numpy as np
from scipy.interpolate import RectBivariateSpline


# ---------------------------------------------------------------------------
# HELIOS species database
# Copied verbatim from HELIOS/source/species_database.py
# Maps HELIOS species name -> FastChem column name (fc_name)
# For CIA/pair species the fc_name is "A&B"; VMR = VMR_A * VMR_B
# ---------------------------------------------------------------------------
SPECIES_FC_NAME = {
    # neutral molecules
    "CO2":        "C1O2",
    "H2O":        "H2O1",
    "CO":         "C1O1",
    "O2":         "O2",
    "CH4":        "C1H4",
    "HCN":        "C1H1N1_1",
    "NH3":        "H3N1",
    "H2S":        "H2S1",
    "PH3":        "H3P1",
    "O3":         "O3",
    "O3_IR":      "O3",
    "O3_UV":      "O3",
    "NO":         "N1O1",
    "SO2":        "O2S1",
    "SH":         "H1S1",
    "H2":         "H2",
    "N2":         "N2",
    "SO":         "O1S1",
    "OH":         "H1O1",
    "COS":        "C1O1S1",
    "CS":         "C1S1",
    "HCHO":       "H2C1O1",
    "C2H4":       "C2H4",
    "C2H2":       "C2H2",
    "CH3":        "C1H3",
    "C3H":        "C3H1",
    "C2H":        "C2H1",
    "C2N2":       "C2N2",
    "C3O2":       "C3O2",
    "C4N2":       "C4N2",
    "C3":         "C3",
    "S2":         "S2",
    "S3":         "S3",
    "S2O":        "O1S2",
    "CS2":        "C1S2",
    "NO2":        "N1O2",
    "N2O":        "N2O1",
    "HNO3":       "H1N1O3",
    "SO3":        "O3S1",
    "H2SO4":      "H2O4S1",
    "TiO":        "O1Ti1",
    "VO":         "O1V1",
    "SiO":        "O1Si1",
    "AlO":        "Al1O1",
    "CaO":        "Ca1O1",
    "PO":         "O1P1",
    "SiH":        "H1Si1",
    "CaH":        "Ca1H1",
    "AlH":        "Al1H1",
    "MgH":        "H1Mg1",
    "CrH":        "Cr1H1",
    "NaH":        "H1Na1",
    "CaOH":       "Ca1H1O1",
    "HCl":        "Cl1H1",
    "NaCl":       "Cl1Na1",
    "NaOH":       "H1Na1O1",
    "SiS":        "S1Si1",
    "KCl":        "Cl1K1",
    "HF":         "F1H1",
    "KOH":        "H1K1O1",
    "NaF":        "F1Na1",
    "KF":         "F1K1",
    "PS":         "P1S1",
    "PN":         "N1P1",
    "PH":         "H1P1",
    "AlF":        "Al1F1",
    "MgF":        "F1Mg1",
    # neutral atoms
    "H":          "H",
    "He":         "He",
    "C":          "C",
    "N":          "N",
    "O":          "O",
    "F":          "F",
    "Na":         "Na",
    "Ne":         "Ne",
    "Ni":         "Ni",
    "Mg":         "Mg",
    "Mn":         "Mn",
    "Al":         "Al",
    "Ar":         "Ar",
    "Si":         "Si",
    "P":          "P",
    "S":          "S",
    "Cl":         "Cl",
    "K":          "K",
    "Ca":         "Ca",
    "Ti":         "Ti",
    "V":          "V",
    "Co":         "Co",
    "Cr":         "Cr",
    "Cu":         "Cu",
    "Fe":         "Fe",
    "Zn":         "Zn",
    # ions
    "H-_bf":      "H1-",
    "H-_ff":      "H&e-",    # product of H × e-
    "He-":        "He&e-",   # product of He × e-
    "Fe+":        "Fe1+",
    "Ti+":        "Ti1+",
    "e-":         "e-",
    # CIA pairs — VMR = VMR_A × VMR_B
    "CIA_H2H2":   "H2&H2",
    "CIA_H2He":   "H2&He",
    "CIA_CO2CO2": "C1O2&C1O2",
    "CIA_O2CO2":  "O2&C1O2",
    "CIA_O2O2":   "O2&O2",
    "CIA_O2N2":   "O2&N2",
    "CIA_N2N2":   "N2&N2",
    "CIA_N2H2":   "N2&H2",
}

# deletechars used by HELIOS when reading FastChem files
_DELETECHARS = " !#$%&'()*,./:;<=>?@[\\]^{|}~"


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def read_helios_tp(path):
    """
    Read a HELIOS _tp.dat file.

    Format (from HELIOS source/write.py):
      Line 0: free-text description
      Line 1: column header
      Line 2+: data rows — first column is layer index (or "BOA"),
               second column is T [K], third column is P [10^-6 bar]

    Returns
    -------
    temps  : np.ndarray, shape (N,), temperature [K]
    press  : np.ndarray, shape (N,), pressure [bar]
    """
    temps, press = [], []
    with open(path) as fh:
        next(fh)  # skip text header
        next(fh)  # skip column header
        for line in fh:
            cols = line.split()
            if not cols:
                continue
            try:
                temps.append(float(cols[1]))
                press.append(float(cols[2]) * 1e-6)  # convert microbar → bar
            except (IndexError, ValueError):
                continue
    return np.array(temps), np.array(press)


def read_species_dat(path):
    """
    Read a HELIOS species.dat file.

    Format:
      Line 0: header
      Line 1+: name  absorbing  scattering  mixing_ratio_source

    Returns list of (helios_name, source) tuples for all species entries.
    H- is expanded to H-_bf and H-_ff exactly as HELIOS does.
    """
    entries = []
    with open(path) as fh:
        next(fh)  # skip header
        for line in fh:
            cols = line.split()
            if not cols:
                continue
            name   = cols[0]
            source = cols[3] if len(cols) >= 4 else None
            if name == "H-":
                entries.append(("H-_bf", source))
                entries.append(("H-_ff", source))
            else:
                entries.append((name, source))
    return entries


def load_fastchem_data(path):
    """
    Read a FastChem chem.dat file using the same numpy.genfromtxt call as HELIOS.

    Returns
    -------
    data       : numpy structured array with one field per column
    press_col  : name of the pressure field  (bar)
    temp_col   : name of the temperature field (K)
    """
    data = np.genfromtxt(
        path,
        names=True,
        dtype=None,
        skip_header=0,
        deletechars=_DELETECHARS,
        encoding=None,
    )
    # Identify P and T column names robustly (FastChem versions differ in case)
    names_lower = {n.lower(): n for n in data.dtype.names}
    press_col = names_lower.get("pbar") or names_lower.get("p_bar")
    temp_col  = names_lower.get("tk")   or names_lower.get("t_k")
    if press_col is None or temp_col is None:
        raise KeyError(
            f"Could not find P/T columns in FastChem file.\n"
            f"  Available names (first 10): {data.dtype.names[:10]}\n"
            f"  Expected lowercase 'pbar' and 'tk'."
        )
    return data, press_col, temp_col


def build_fastchem_grid(data, press_col, temp_col):
    """
    Extract sorted unique T and P arrays from the FastChem structured array.
    Mirrors HELIOS load_fastchem_data().

    FastChem grids are ordered T-outer, P-inner (pressure varies fastest).

    Returns
    -------
    t_grid : np.ndarray  sorted unique temperatures [K]
    p_grid : np.ndarray  sorted unique pressures [bar]
    """
    t_grid = np.array(sorted(set(data[temp_col])))
    p_grid = np.array(sorted(set(data[press_col])))   # already in bar
    return t_grid, p_grid


def get_vmr_1d(data, fc_name):
    """
    Return the 1-D VMR array for a species from the FastChem structured array.

    For pair species (fc_name contains '&') the VMR is the product of the two
    component columns — exactly as HELIOS does for CIA / H-_ff / He-.

    The returned array is indexed as vmr[p_idx + n_p * t_idx].
    """
    if "&" in fc_name:
        a, b = fc_name.split("&", 1)
        try:
            vmr_a = data[a].astype(float)
        except ValueError:
            raise KeyError(f"FastChem column '{a}' not found (from fc_name '{fc_name}')")
        try:
            vmr_b = data[b].astype(float)
        except ValueError:
            raise KeyError(f"FastChem column '{b}' not found (from fc_name '{fc_name}')")
        return vmr_a * vmr_b
    else:
        try:
            return data[fc_name].astype(float)
        except ValueError:
            raise KeyError(f"FastChem column '{fc_name}' not found")


def interpolate_vmr_to_profile(vmr_1d, t_grid, p_grid, t_profile, p_profile_bar):
    """
    Interpolate a 2-D FastChem VMR grid onto a 1-D atmospheric T-P profile.

    Replicates the two-step HELIOS approach collapsed into one:
      1. Reshape 1-D array (T-outer, P-inner) to (nT, nP)
      2. Build RectBivariateSpline(kx=1, ky=1) in (T, log10 P)
      3. Clamp T and log10P to grid range (edge extrapolation = nearest value)
      4. Evaluate at each (T_i, P_i) point in the profile

    Parameters
    ----------
    vmr_1d       : 1-D array, length nT*nP, indexed vmr[p + nP*t]
    t_grid       : sorted unique T values from FastChem [K]
    p_grid       : sorted unique P values from FastChem [bar]
    t_profile    : T at each atmospheric layer [K]
    p_profile_bar: P at each atmospheric layer [bar]

    Returns
    -------
    vmr_profile  : 1-D array, length = len(t_profile)
    """
    n_t = len(t_grid)
    n_p = len(p_grid)

    vmr_2d = vmr_1d.reshape((n_t, n_p))           # shape: (nT, nP)
    log_p_grid = np.log10(p_grid)

    spline = RectBivariateSpline(t_grid, log_p_grid, vmr_2d, kx=1, ky=1)

    # Clamp to grid edges — matches HELIOS boundary behaviour
    t_clamped    = np.clip(t_profile, t_grid.min(), t_grid.max())
    logp_clamped = np.clip(np.log10(p_profile_bar),
                           log_p_grid.min(), log_p_grid.max())

    vmr_profile = np.array([
        spline(t_clamped[i], logp_clamped[i])[0, 0]
        for i in range(len(t_profile))
    ])
    return vmr_profile


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Interpolate FastChem VMRs onto a HELIOS T-P profile."
    )
    parser.add_argument("--tp",      required=True, help="HELIOS _tp.dat file")
    parser.add_argument("--chem",    required=True, help="FastChem chem.dat file")
    parser.add_argument("--species", required=True, help="HELIOS species.dat file")
    parser.add_argument("--output",  default="vmr_interpolated.dat",
                        help="Output file (default: vmr_interpolated.dat)")
    args = parser.parse_args()

    # ── 1. Read the T-P profile ────────────────────────────────────────────
    print(f"Reading T-P profile: {args.tp}")
    t_profile, p_profile = read_helios_tp(args.tp)
    print(f"  {len(t_profile)} layers, "
          f"T = [{t_profile.min():.1f}, {t_profile.max():.1f}] K, "
          f"P = [{p_profile.min():.3e}, {p_profile.max():.3e}] bar")

    # ── 2. Read the FastChem grid ──────────────────────────────────────────
    print(f"Reading FastChem grid: {args.chem}")
    fc_data, press_col, temp_col = load_fastchem_data(args.chem)
    t_grid, p_grid = build_fastchem_grid(fc_data, press_col, temp_col)
    print(f"  FastChem grid: {len(t_grid)} T-points × {len(p_grid)} P-points")
    print(f"  T = [{t_grid.min():.1f}, {t_grid.max():.1f}] K, "
          f"P = [{p_grid.min():.3e}, {p_grid.max():.3e}] bar")

    # ── 3. Read the species list ───────────────────────────────────────────
    print(f"Reading species list: {args.species}")
    all_species = read_species_dat(args.species)
    fc_species  = [(name, src) for name, src in all_species if src == "FastChem"]
    print(f"  {len(fc_species)} FastChem species: "
          f"{[n for n, _ in fc_species]}")

    # ── 4. Interpolate each species ────────────────────────────────────────
    results = {}       # helios_name -> vmr_profile array
    skipped = []

    for helios_name, _ in fc_species:
        fc_name = SPECIES_FC_NAME.get(helios_name)
        if fc_name is None:
            print(f"  WARNING: '{helios_name}' not found in built-in species "
                  f"database — skipping.")
            skipped.append(helios_name)
            continue

        # Check for "not included in FastChem" entries
        if "not included" in fc_name.lower():
            print(f"  WARNING: '{helios_name}' has no FastChem equivalent "
                  f"({fc_name}) — skipping.")
            skipped.append(helios_name)
            continue

        try:
            vmr_1d = get_vmr_1d(fc_data, fc_name)
        except KeyError as exc:
            print(f"  WARNING: {exc} — skipping '{helios_name}'.")
            skipped.append(helios_name)
            continue

        vmr_profile = interpolate_vmr_to_profile(
            vmr_1d, t_grid, p_grid, t_profile, p_profile
        )
        results[helios_name] = vmr_profile
        print(f"  {helios_name:20s}  (fc: {fc_name})")

    if skipped:
        print(f"\n  Skipped species: {skipped}")

    # ── 5. Write output ────────────────────────────────────────────────────
    print(f"\nWriting output: {args.output}")

    species_cols = list(results.keys())
    header_note = (
        "# Interpolated FastChem VMRs on HELIOS T-P profile.\n"
        "# Method: bilinear in (T [K], log10 P [bar]) — RectBivariateSpline(kx=1,ky=1),\n"
        "#         edge-clamped. Mirrors HELIOS on-the-fly FastChem interpolation.\n"
        "# VMR values are as output by FastChem (number densities [cm^-3] if FastChem\n"
        "#   was run with default settings, or mixing ratios if configured for that).\n"
        "# T-P profile source: " + args.tp + "\n"
        "# FastChem grid source: " + args.chem + "\n"
    )

    col_names = ["P[bar]", "T[K]"] + species_cols
    header_row = "\t".join(col_names)

    with open(args.output, "w") as fh:
        fh.write(header_note)
        fh.write(header_row + "\n")
        for i in range(len(t_profile)):
            row = [f"{p_profile[i]:.6e}", f"{t_profile[i]:.4f}"]
            for name in species_cols:
                row.append(f"{results[name][i]:.6e}")
            fh.write("\t".join(row) + "\n")

    print(f"Done. {len(t_profile)} rows × {len(species_cols)+2} columns "
          f"written to {args.output}.")
    if skipped:
        print(f"Note: {len(skipped)} species were skipped (see warnings above).")


if __name__ == "__main__":
    main()
