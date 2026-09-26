#!/usr/bin/env python3
"""
Analysis for Experiments 8 and 9.
Project roll number: 25m0328

Compares each candidate-spectrum run against the nominal A=1 solution using
the exact VULCAN wavelength quadrature, and reports the radiative response,
its vertical distribution, and the sulfur photolysis response.

IMPORTANT: the reference solution must have been produced on the same machine
and with the same library versions as the candidate runs. Cross-machine
comparison is not valid at this signal level; see the reproducibility note in
the Experiment 8 and 9 documentation.

Run from the repository root.
"""

import pickle
import os
import numpy as np

# The reference must come from the SAME machine and library versions as the
# candidate runs. run_experiment8_9_25m0328.sh produces that baseline first.
# The committed A=1 solution is used only as a fallback, with a warning,
# because a cross-machine comparison is not valid at this signal level.
BASELINE_LOCAL = "output/Experiment8_baseline_A1-25m0328.vul"
BASELINE_COMMITTED = "output/Experiment2_A1_UUV-nominal-25m0328.vul"

if os.path.exists(BASELINE_LOCAL):
    REFERENCE = BASELINE_LOCAL
else:
    REFERENCE = BASELINE_COMMITTED
    print("WARNING: the same-machine baseline was not found.")
    print(f"         Falling back to {BASELINE_COMMITTED}.")
    print("         Solver-path differences between machines can exceed the signal")
    print("         being measured. Run the baseline case before trusting these numbers.")
    print()

CASES = {
    "8A ClS2 spectrum":      "output/Experiment8A_ClS2-25m0328.vul",
    "8B OSSO spectrum":      "output/Experiment8B_OSSO_S2O2-25m0328.vul",
    "8C S3 spectrum":        "output/Experiment8C_S3-25m0328.vul",
    "9  ClS2 self-consistent": "output/Experiment9_ClS2-selfconsistent-25m0328.vul",
}

PARENTS = ["S2O2", "SCl2", "ClS2", "S3", "S2Cl2", "SCl", "S2O", "S4", "SO2", "SO", "SO3"]

BAND_LO_KM, BAND_HI_KM = 69.0, 93.0


def load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def main():
    print(f"reference solution: {REFERENCE}\n")
    ref = load(REFERENCE)
    z = np.asarray(ref["atm"]["zmco"]) / 1e5
    dz_km = np.asarray(ref["atm"]["dz"]) / 1e5
    bins = np.asarray(ref["variable"]["bins"])
    weights = np.where(bins < 240.0, 0.1, 2.0)   # exact VULCAN quadrature
    f0 = np.asarray(ref["variable"]["aflux"])

    norm_global = np.sqrt(np.sum(f0 ** 2 * weights[None, :]))
    norm_layer = np.sqrt(np.sum(f0 ** 2 * weights[None, :], axis=1))
    band = (z >= BAND_LO_KM) & (z <= BAND_HI_KM)

    results = {}
    print("Radiative response relative to the nominal A=1 solution\n")
    head = (f"{'case':26s} {'max|dF|':>13s} {'z':>5s} {'lam':>6s} "
            f"{'globalL2':>10s} {'mean alt':>9s} {'69-93km':>8s}")
    print(head)
    print("-" * len(head))

    for label, path in CASES.items():
        if not os.path.exists(path):
            print(f"{label:26s} missing: {path}")
            continue
        case = load(path)
        results[label] = case
        d = np.asarray(case["variable"]["aflux"]) - f0
        ij = np.unravel_index(np.argmax(np.abs(d)), d.shape)
        l2 = np.sqrt(np.sum(d ** 2 * weights[None, :]))
        l1z = np.sum(np.abs(d) * weights[None, :], axis=1)
        w = l1z * dz_km
        p = w / w.sum()
        print(f"{label:26s} {d[ij]:+13.4e} {z[ij[0]]:5.0f} {bins[ij[1]]:6.0f} "
              f"{100 * l2 / norm_global:9.4f}% {np.sum(z * p):8.2f}km "
              f"{100 * p[band].sum():7.2f}%")

    print()
    print("Maximum local fractional response")
    for label, case in results.items():
        d = np.asarray(case["variable"]["aflux"]) - f0
        l2z = np.sqrt(np.sum(d ** 2 * weights[None, :], axis=1))
        frac = l2z / np.maximum(norm_layer, 1e-300)
        print(f"  {label:26s} {100 * frac.max():7.2f}% at {z[np.argmax(frac)]:.0f} km")

    print()
    print(f"Sulfur photolysis response, maximum change in {BAND_LO_KM:.0f} to "
          f"{BAND_HI_KM:.0f} km, as a percentage of the nominal rate")
    print(f"{'parent':8s} " + " ".join(f"{k.split()[0]:>10s}" for k in results))
    for parent in PARENTS:
        key = (parent, 0)
        if key not in ref["variable"]["J_sp"]:
            continue
        a = np.asarray(ref["variable"]["J_sp"][key])
        row = f"{parent:8s} "
        for case in results.values():
            d = np.asarray(case["variable"]["J_sp"][key]) - a
            j = int(np.argmax(np.abs(d * band)))
            row += f" {100 * d[j] / a[j] if a[j] else float('nan'):+9.3f}%"
        print(row)


if __name__ == "__main__":
    main()
