#!/usr/bin/env python3
"""
Generate the Experiment 8 and Experiment 9 input files.
Project roll number: 25m0328

Experiment 8 substitutes the empirical UUV absorption spectrum with the
absorption spectrum of a candidate species, normalised so that the
band-integrated optical depth over 300 to 698 nm is unchanged. Vertical
placement and A_UUV are left at nominal, so the only variable is spectral
shape.

Experiment 9 additionally replaces the mode-1 vertical profile with the
candidate's own modelled vertical distribution, with the column preserved,
so that both spectral shape and vertical placement are the candidate's.

The candidate spectra are taken from the cross sections already loaded inside
the model, so no external optical data is introduced and the comparison is
internally consistent.

Run from the repository root. Writes into atm/.
"""

import pickle
import numpy as np
from scipy import interpolate

A1_PATH = "output/Experiment2_A1_UUV-nominal-25m0328.vul"
NOMINAL_SPECTRUM = "atm/UV_absorber.txt"
NOMINAL_PROFILE = "atm/mode1+2.txt"

GRID_LO_NM, GRID_HI_NM = 300.0, 698.0

# tag -> species whose cross section supplies the spectral shape
SPECTRUM_CASES = {
    "8A_ClS2": "ClS2",
    "8B_OSSO_S2O2": "S2O2",
    "8C_S3": "S3",
}

# Experiment 9: species whose modelled profile supplies the vertical shape
PROFILE_CASE = "ClS2"


def main():
    with open(A1_PATH, "rb") as f:
        sol = pickle.load(f)
    var, atm = sol["variable"], sol["atm"]
    bins = np.asarray(var["bins"])
    z_km = np.asarray(atm["zmco"]) / 1e5
    species = list(var["species"])
    y = np.asarray(var["y"])

    grid = bins[(bins >= GRID_LO_NM) & (bins <= GRID_HI_NM)]
    dlam = np.full_like(grid, 2.0)  # model resolution above 240 nm

    uv = np.loadtxt(NOMINAL_SPECTRUM, skiprows=1)
    q_nominal = interpolate.interp1d(uv[:, 0], uv[:, 1], kind="slinear")(grid)
    integral_nominal = float(np.sum(q_nominal * dlam))

    print(f"nominal spectrum: peak Q {q_nominal.max():.4f} at "
          f"{grid[np.argmax(q_nominal)]:.0f} nm, band integral {integral_nominal:.4f} nm")
    print()

    for tag, name in SPECTRUM_CASES.items():
        sigma = np.asarray(var["cross"][name])[(bins >= GRID_LO_NM) & (bins <= GRID_HI_NM)]
        if sigma.max() <= 0:
            print(f"{tag}: no cross-section support, skipped")
            continue
        scale = integral_nominal / float(np.sum(sigma * dlam))
        q = scale * sigma
        path = f"atm/UV_absorber_Experiment{tag}_25m0328.txt"
        with open(path, "w") as f:
            f.write("lambda(nm)    Q_abs\n")
            for lam, val in zip(grid, q):
                f.write(f"  {lam:9.4f}  {val:.10f}\n")
        support = grid[sigma > 0]
        print(f"{tag:14s} peak Q {q.max():7.4f} at {grid[np.argmax(q)]:.0f} nm | "
              f"support {support.min():.0f} to {support.max():.0f} nm | "
              f"integral {np.sum(q * dlam):.4f} -> {path}")

    # Experiment 9 vertical profile
    src = np.loadtxt(NOMINAL_PROFILE, skiprows=1)
    alt, n1, n2 = src[:, 0], src[:, 1], src[:, 2]
    density = y[:, species.index(PROFILE_CASE)]
    shape = np.interp(alt, z_km, density)
    shape = shape * n1.sum() / shape.sum()

    path = f"atm/mode1+2_Experiment9_{PROFILE_CASE}profile_25m0328.txt"
    with open(path, "w") as f:
        f.write("alt   mode1(cm-3)  mode2(cm-3) # Experiment 9: mode-1 replaced by "
                f"modelled {PROFILE_CASE} shape, column preserved\n")
        for a, m1, m2 in zip(alt, shape, n2):
            f.write(f"{a:<5.0f} {m1:<18.10e} {m2:<12.6g}\n")

    centroid = lambda p: float(np.sum(alt * p) / np.sum(p))
    in_band = lambda p: 100.0 * float(np.sum(p[(alt >= 69) & (alt <= 93)]) / np.sum(p))
    print()
    print(f"nominal profile: column {n1.sum():.6e}, centroid {centroid(n1):.2f} km, "
          f"{in_band(n1):.2f}% in 69 to 93 km")
    print(f"{PROFILE_CASE} profile: column {shape.sum():.6e}, centroid {centroid(shape):.2f} km, "
          f"{in_band(shape):.2f}% in 69 to 93 km")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
