#!/usr/bin/env python3
"""
Candidate opacity screen for the VULCAN-Venus Unknown UV Absorber project.
Project roll number: 25m0328

Purpose
-------
Read-only test. For every sulfur- and chlorine-bearing species that the model
already tracks AND for which the model already stores an absorption cross
section, compute the column optical depth that the converged steady state
actually produces, and compare it with the optical depth the empirical UUV
requires.

This is the same test that has been applied to OSSO in the literature
(modelled abundance far too low to supply the observed absorption), extended
to every gas-phase candidate available inside the network.

No VULCAN rerun is required. Run from the repository root.
"""

import pickle
import numpy as np

A1_PATH = "output/Experiment2_A1_UUV-nominal-25m0328.vul"

# Total UUV optical depth at 320 nm implied by the nominal absorber loading,
# established in the UUV optical-depth reconstruction.
TAU_REQUIRED = 1.094968472081

# Cloud base for the column integral. The entire nominal UUV particle column
# lies above this altitude, so the required optical depth and the candidate
# column are compared over the same path.
Z_BASE_KM = 58.0

BAND_LO_NM, BAND_HI_NM = 300.0, 620.0

CANDIDATES = [
    "SCl2", "ClS2", "S2Cl2", "SCl",
    "S3", "S4", "S2", "S2O", "S2O2",
    "SO3", "Cl2", "ClO", "OSCl",
]


def main():
    with open(A1_PATH, "rb") as f:
        sol = pickle.load(f)

    var, atm = sol["variable"], sol["atm"]
    z_km = np.asarray(atm["zmco"]) / 1e5
    dz_cm = np.asarray(atm["dz"])
    bins = np.asarray(var["bins"])
    y = np.asarray(var["y"])
    ymix = np.asarray(var["ymix"])
    species = list(var["species"])

    above = z_km >= Z_BASE_KM
    band = (bins >= BAND_LO_NM) & (bins <= BAND_HI_NM)

    print(f"Required UUV optical depth at 320 nm : {TAU_REQUIRED:.6f}")
    print(f"Column integrated above              : {Z_BASE_KM:.0f} km")
    print(f"Wavelength window                    : {BAND_LO_NM:.0f} to {BAND_HI_NM:.0f} nm")
    print()
    header = (f"{'species':8s} {'peak lam':>9s} {'sigma peak':>12s} "
              f"{'column':>12s} {'tau max':>11s} {'tau/req':>10s} {'short by':>11s}")
    print(header)
    print("-" * len(header))

    rows = []
    for name in CANDIDATES:
        if name not in species or name not in var["cross"]:
            continue
        idx = species.index(name)
        sigma = np.asarray(var["cross"][name])
        column = float(np.sum(y[:, idx][above] * dz_cm[above]))
        tau = sigma * column
        tau_in_band = np.where(band, tau, 0.0)
        j = int(np.argmax(tau_in_band))
        if tau_in_band[j] <= 0.0:
            print(f"{name:8s}   no cross-section support in the window")
            continue
        rows.append((tau_in_band[j], name))
        print(f"{name:8s} {bins[j]:9.1f} {sigma[j]:12.4e} {column:12.4e} "
              f"{tau_in_band[j]:11.4e} {tau_in_band[j] / TAU_REQUIRED:10.2e} "
              f"{TAU_REQUIRED / tau_in_band[j]:10.4g}x")

    if not rows:
        return

    best_tau, best = max(rows)
    print()
    print(f"Closest gas-phase candidate: {best}")

    idx = species.index(best)
    sigma = np.asarray(var["cross"][best])
    column = float(np.sum(y[:, idx][above] * dz_cm[above]))
    j = int(np.argmax(np.where(band, sigma * column, 0.0)))
    needed = TAU_REQUIRED / sigma[j]

    print(f"  absorption support      : {bins[sigma > 0].min():.0f} to "
          f"{bins[sigma > 0].max():.0f} nm, peak at {bins[np.argmax(sigma)]:.0f} nm")
    print(f"  modelled column         : {column:.4e} cm-2")
    print(f"  column needed           : {needed:.4e} cm-2")
    print(f"  shortfall               : {needed / column:.4g}x")

    hcl = float(np.sum(y[:, species.index("HCl")][above] * dz_cm[above]))
    print(f"  as a fraction of the HCl column above {Z_BASE_KM:.0f} km: {needed / hcl:.3f}")

    layer = (z_km >= 58.0) & (z_km <= 70.0)
    path = float(np.sum(dz_cm[layer]))
    mean_M = float(np.average(np.asarray(atm["M"])[layer], weights=dz_cm[layer]))
    print(f"  required mean mixing ratio if held in 58 to 70 km: {needed / path / mean_M:.3e}")
    print(f"  modelled mixing ratio at 60 km                   : "
          f"{ymix[np.argmin(np.abs(z_km - 60)), idx]:.3e}")


if __name__ == "__main__":
    main()
