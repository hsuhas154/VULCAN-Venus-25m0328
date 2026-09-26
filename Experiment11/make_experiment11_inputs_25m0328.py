#!/usr/bin/env python3
"""
Generate the Experiment 11 absorber spectra from the ClSSCl and SSCl2
computed photoabsorption cross sections.
Project roll number: 25m0328

Experiment 11 substitutes the empirical UUV absorption spectrum with the
computed spectrum of a sulfur-chlorine isomer, normalised so that the
band-integrated optical depth over 300 to 698 nm is unchanged. Vertical
placement and A_UUV stay at nominal, so the only variable is spectral shape.
This is the same protocol as Experiment 8, which makes the five candidate
spectra directly comparable.

Unlike Experiment 8, whose candidate spectra were taken from cross sections
already inside the model, these come from an external source: the computed
cross sections of Trabelsi and Francisco (2026), supplied by the authors.
See data/PROVENANCE.md.

Writes into atm/. Run from the repository root:
    python3 Experiment11/make_experiment11_inputs_25m0328.py
"""

import os
import numpy as np
from scipy import interpolate
import pickle

A1_PATH = "output/Experiment2_A1_UUV-nominal-25m0328.vul"
NOMINAL_SPECTRUM = "atm/UV_absorber.txt"
DATA_DIR = "Experiment11/data"

GRID_LO_NM, GRID_HI_NM = 300.0, 698.0

CASES = {
    "11A_SSCl2":  "csec-sscl2.dat",
    "11B_ClSSCl": "csec-clsscl.dat",
}

# Published Figure 2 values, used as a sanity check on the supplied files.
EXPECTED = {
    "csec-sscl2.dat":  (3.37e-17, 264.0),
    "csec-clsscl.dat": (2.63e-17, 240.0),
}


def read_cross_section(path):
    """Two columns, wavelength nm and cross section cm2 per molecule,
    ordered long to short wavelength. Returned sorted ascending."""
    d = np.loadtxt(path)
    lam, sigma = d[:, 0], d[:, 1]
    order = np.argsort(lam)
    return lam[order], sigma[order]


def main():
    for name in CASES.values():
        p = os.path.join(DATA_DIR, name)
        if not os.path.exists(p):
            raise SystemExit(f"missing {p}. See {DATA_DIR}/PROVENANCE.md")

    with open(A1_PATH, "rb") as f:
        sol = pickle.load(f)
    bins = np.asarray(sol["variable"]["bins"])
    grid = bins[(bins >= GRID_LO_NM) & (bins <= GRID_HI_NM)]
    dlam = np.full_like(grid, 2.0)          # model resolution above 240 nm

    uv = np.loadtxt(NOMINAL_SPECTRUM, skiprows=1)
    q_nominal = interpolate.interp1d(uv[:, 0], uv[:, 1], kind="slinear")(grid)
    integral_nominal = float(np.sum(q_nominal * dlam))

    print(f"nominal spectrum: peak Q {q_nominal.max():.4f} at "
          f"{grid[np.argmax(q_nominal)]:.0f} nm, band integral "
          f"{integral_nominal:.4f} nm\n")

    for tag, fname in CASES.items():
        lam, sigma = read_cross_section(os.path.join(DATA_DIR, fname))

        peak, peak_lam = sigma.max(), lam[np.argmax(sigma)]
        exp_peak, exp_lam = EXPECTED[fname]
        ok = (abs(peak - exp_peak) / exp_peak < 0.05) and (abs(peak_lam - exp_lam) < 5.0)
        print(f"{tag}: global peak {peak:.4e} cm2 at {peak_lam:.1f} nm "
              f"(figure: {exp_peak:.2e} near {exp_lam:.0f} nm) "
              f"{'OK' if ok else 'MISMATCH, check the file'}")
        if not ok:
            raise SystemExit("supplied cross section does not match the published "
                             "figure; stopping rather than building a bad input")

        in_band = (lam >= 300) & (lam <= 620)
        print(f"        peak inside 300 to 620 nm: {sigma[in_band].max():.4e} at "
              f"{lam[in_band][np.argmax(sigma[in_band])]:.1f} nm")

        mapped = interpolate.interp1d(lam, sigma, kind="linear",
                                      bounds_error=False, fill_value=0.0)(grid)
        if mapped.max() <= 0:
            raise SystemExit(f"{tag}: no support on the model grid")

        scale = integral_nominal / float(np.sum(mapped * dlam))
        q = scale * mapped

        path = f"atm/UV_absorber_Experiment{tag}_25m0328.txt"
        with open(path, "w") as f:
            f.write("lambda(nm)    Q_abs\n")
            for wl, val in zip(grid, q):
                f.write(f"  {wl:9.4f}  {val:.10f}\n")

        support = grid[mapped > 0]
        print(f"        scaled to matched band opacity: peak Q {q.max():.4f} at "
              f"{grid[np.argmax(q)]:.0f} nm, scale {scale:.4e}")
        print(f"        support {support.min():.0f} to {support.max():.0f} nm, "
              f"integral {np.sum(q * dlam):.4f}")
        print(f"        wrote {path}\n")


if __name__ == "__main__":
    main()
