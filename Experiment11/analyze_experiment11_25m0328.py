#!/usr/bin/env python3
"""
Experiment 11 analysis: the full five-spectrum comparison, plus the corrected
sulfur-chlorine opacity arithmetic using the authors' computed cross sections.
Project roll number: 25m0328

Part 1 places the two sulfur-chlorine isomer spectra alongside the three
Experiment 8 candidates. All five carry the same band-integrated opacity at
the same vertical placement, so differences are due to spectral shape alone.

Part 2 repeats the candidate opacity screen for the isomers, using their real
computed cross sections rather than values read from a published figure, and
evaluating each inside the 300 to 620 nm range in which the absorber operates
rather than at its global peak.

Run from the repository root. Read-only, no VULCAN rerun.
"""

import os
import pickle
import numpy as np

REFERENCE = "output/Experiment8_baseline_A1-25m0328.vul"
NOMINAL = "output/Experiment2_A1_UUV-nominal-25m0328.vul"
DATA_DIR = "Experiment11/data"

CASES = {
    "8B  OSSO":    "output/Experiment8B_OSSO_S2O2-25m0328.vul",
    "8A  ClS2":    "output/Experiment8A_ClS2-25m0328.vul",
    "8C  S3":      "output/Experiment8C_S3-25m0328.vul",
    "11A SSCl2":   "output/Experiment11A_SSCl2-25m0328.vul",
    "11B ClSSCl":  "output/Experiment11B_ClSSCl-25m0328.vul",
}

PARENTS = ["S2O2", "SCl2", "ClS2", "S3", "S2Cl2", "SCl", "S2O", "S4", "SO2", "SO"]

TAU_REQUIRED = 1.094969
BAND_LO_NM, BAND_HI_NM = 300.0, 620.0
Z_BASE_KM = 58.0


def load(p):
    with open(p, "rb") as f:
        return pickle.load(f)


def main():
    ref = load(REFERENCE)
    z = np.asarray(ref["atm"]["zmco"]) / 1e5
    dz_cm = np.asarray(ref["atm"]["dz"])
    dz_km = dz_cm / 1e5
    bins = np.asarray(ref["variable"]["bins"])
    quad = np.where(bins < 240.0, 0.1, 2.0)
    f0 = np.asarray(ref["variable"]["aflux"])
    norm = np.sqrt(np.sum(f0 ** 2 * quad[None, :]))
    band_z = (z >= 69) & (z <= 93)

    print("=" * 76)
    print("PART 1: five candidate spectra, identical band opacity and placement")
    print("=" * 76)
    head = (f"{'case':12s} {'peak dF':>14s} {'z':>5s} {'nm':>6s} "
            f"{'global L2':>11s} {'mean alt':>9s}")
    print(head)
    print("-" * len(head))
    present = {}
    for label, path in CASES.items():
        if not os.path.exists(path):
            print(f"{label:12s} missing: {path}")
            continue
        case = load(path)
        present[label] = case
        d = np.asarray(case["variable"]["aflux"]) - f0
        ij = np.unravel_index(np.argmax(np.abs(d)), d.shape)
        l2 = np.sqrt(np.sum(d ** 2 * quad[None, :]))
        l1z = np.sum(np.abs(d) * quad[None, :], axis=1)
        p = (l1z * dz_km) / (l1z * dz_km).sum()
        print(f"{label:12s} {d[ij]:+14.5e} {z[ij[0]]:5.0f} {bins[ij[1]]:6.0f} "
              f"{100 * l2 / norm:10.4f}% {np.sum(z * p):8.2f}km")

    if present:
        print(f"\nFor scale, an 18 km vertical displacement of the same absorber "
              f"gives 2.2498%.")
        print("Spectral shape and vertical placement are separate controls; "
              "placement is the larger one.")

        print("\nSulfur photolysis response, max change in 69 to 93 km, dJ convention")
        cols = list(present)
        print(f"{'parent':8s} " + " ".join(f"{c.split()[1]:>10s}" for c in cols))
        for parent in PARENTS:
            key = (parent, 0)
            if key not in ref["variable"]["J_sp"]:
                continue
            a = np.asarray(ref["variable"]["J_sp"][key])
            row = f"{parent:8s} "
            for c in cols:
                d = np.asarray(present[c]["variable"]["J_sp"][key]) - a
                j = int(np.argmax(np.abs(d * band_z)))
                row += f" {100 * d[j] / a[j] if a[j] else np.nan:+9.3f}%"
            print(row)

    print("\n" + "=" * 76)
    print("PART 2: opacity screen for the isomers, computed cross sections")
    print("=" * 76)
    nom = load(NOMINAL)
    species = list(nom["variable"]["species"])
    y = np.asarray(nom["variable"]["y"])
    above = z >= Z_BASE_KM
    col_s2cl2 = float(np.sum(y[:, species.index("S2Cl2")][above] * dz_cm[above]))

    print(f"modelled S2Cl2 column above {Z_BASE_KM:.0f} km: {col_s2cl2:.4e} cm-2")
    print(f"required optical depth at 320 nm: {TAU_REQUIRED:.6f}\n")
    print(f"{'isomer':10s} {'global peak':>22s} {'peak in band':>24s} "
          f"{'tau':>11s} {'short by':>11s}")

    for tag, fname in [("SSCl2", "csec-sscl2.dat"), ("ClSSCl", "csec-clsscl.dat")]:
        p = os.path.join(DATA_DIR, fname)
        if not os.path.exists(p):
            print(f"{tag:10s} missing {p}")
            continue
        d = np.loadtxt(p)
        lam, sig = d[:, 0], d[:, 1]
        o = np.argsort(lam)
        lam, sig = lam[o], sig[o]
        m = (lam >= BAND_LO_NM) & (lam <= BAND_HI_NM)
        gpk, glam = sig.max(), lam[np.argmax(sig)]
        bpk, blam = sig[m].max(), lam[m][np.argmax(sig[m])]
        tau = bpk * col_s2cl2
        print(f"{tag:10s} {gpk:12.4e} at {glam:6.1f}nm {bpk:12.4e} at {blam:6.1f}nm "
              f"{tau:11.4e} {TAU_REQUIRED / tau:10.4g}x")

    print("\nNote: the global peak of each isomer lies below 300 nm, outside the range")
    print("in which the unknown absorber operates. Screening must use the in-band peak.")


if __name__ == "__main__":
    main()
