#!/usr/bin/env python3
"""
Experiment 10: does the below-cloud response have a universal fingerprint?
Project roll number: 25m0328

Experiments 7B, 9 and 10 all place the absorber below the photochemically
active region, but differ in spectrum and in vertical distribution:

    7B   empirical UUV spectrum, nominal profile shifted down 18 km
    9    ClS2 spectrum,          ClS2 modelled vertical distribution
    10   S3 spectrum,            ClS2 modelled vertical distribution

If the sulfur photolysis response of each is the same pattern scaled by a
single number, then the below-cloud response is a property of the atmosphere
and the absorber sets only the amplitude. The test is whether the channel by
channel ratios between any two cases are constant.

Each case is compared against a reference from its own run batch:
    7B   against the committed A=1 solution
    9,10 against the Experiment 8 baseline

Cross-batch comparison of absolute converged states is not valid, but these
are normalised percentage responses, and the empirical agreement between 7B
and 9 (0.2% spread) demonstrates that the paired differences are stable
across batches. See the reproducibility note in the Experiment 8 and 9
documentation.

Run from the repository root. Read-only, no VULCAN rerun.
"""

import os
import pickle
import numpy as np

PAIRS = {
    "7B": ("output/Experiment2_A1_UUV-nominal-25m0328.vul",
           "output/Experiment7B_UUV-downward18-25m0328.vul"),
    "9":  ("output/Experiment8_baseline_A1-25m0328.vul",
           "output/Experiment9_ClS2-selfconsistent-25m0328.vul"),
    "10": ("output/Experiment8_baseline_A1-25m0328.vul",
           "output/Experiment10_S3-on-ClS2profile-25m0328.vul"),
}

# For contrast: an above-cloud case, which should NOT match the pattern.
CONTRAST = {
    "7A": ("output/Experiment2_A1_UUV-nominal-25m0328.vul",
           "output/Experiment7A_UUV-upward18-25m0328.vul"),
}

PARENTS = ["S2O2", "SCl2", "ClS2", "S3", "S2Cl2", "SCl", "S2O", "S4"]
BAND_LO_KM, BAND_HI_KM = 69.0, 93.0


def load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def fingerprint(ref_path, case_path):
    """Percentage change of each parent's total photolysis rate, at the
    altitude of maximum absolute change inside the response band."""
    ref, case = load(ref_path), load(case_path)
    z = np.asarray(ref["atm"]["zmco"]) / 1e5
    band = (z >= BAND_LO_KM) & (z <= BAND_HI_KM)
    out = {}
    for parent in PARENTS:
        key = (parent, 0)
        if key not in ref["variable"]["J_sp"]:
            continue
        a = np.asarray(ref["variable"]["J_sp"][key])
        d = np.asarray(case["variable"]["J_sp"][key]) - a
        j = int(np.argmax(np.abs(d * band)))
        out[parent] = 100.0 * d[j] / a[j] if a[j] else float("nan")
    return out


def radiative(ref_path, case_path):
    ref, case = load(ref_path), load(case_path)
    bins = np.asarray(ref["variable"]["bins"])
    w = np.where(bins < 240.0, 0.1, 2.0)
    f0 = np.asarray(ref["variable"]["aflux"])
    d = np.asarray(case["variable"]["aflux"]) - f0
    return 100.0 * np.sqrt(np.sum(d ** 2 * w[None, :])) / np.sqrt(np.sum(f0 ** 2 * w[None, :]))


def main():
    have, fps, rad = [], {}, {}
    for label, (ref, case) in {**PAIRS, **CONTRAST}.items():
        if not (os.path.exists(ref) and os.path.exists(case)):
            print(f"missing files for case {label}, skipping")
            continue
        have.append(label)
        fps[label] = fingerprint(ref, case)
        rad[label] = radiative(ref, case)

    below = [c for c in ("7B", "9", "10") if c in have]
    if len(below) < 2:
        print("\nNeed at least two below-cloud cases to test the scaling.")
        return

    print("\nSulfur photolysis fingerprint, percentage change in "
          f"{BAND_LO_KM:.0f} to {BAND_HI_KM:.0f} km\n")
    cols = [c for c in ("7B", "9", "10", "7A") if c in have]
    print(f"{'parent':8s} " + " ".join(f"{c:>10s}" for c in cols))
    print("-" * (8 + 11 * len(cols)))
    for parent in PARENTS:
        row = f"{parent:8s} "
        for c in cols:
            row += f" {fps[c].get(parent, float('nan')):+9.3f}%"
        print(row)

    print("\nGlobal radiative response, percentage of the reference norm")
    for c in cols:
        print(f"  {c:4s} {rad[c]:7.4f}%")

    base = below[0]
    print(f"\nChannel by channel ratios against case {base}")
    print("A constant column means the fingerprint shape is identical and only "
          "the amplitude differs.\n")
    others = below[1:] + [c for c in ("7A",) if c in have]
    print(f"{'parent':8s} " + " ".join(f"{c + '/' + base:>12s}" for c in others))
    print("-" * (8 + 13 * len(others)))
    ratios = {c: [] for c in others}
    for parent in PARENTS:
        row = f"{parent:8s} "
        for c in others:
            a, b = fps[c].get(parent), fps[base].get(parent)
            if a is None or b in (None, 0) or not np.isfinite(a) or not np.isfinite(b):
                row += f" {'n/a':>12s}"
                continue
            r = a / b
            ratios[c].append(r)
            row += f" {r:12.4f}"
        print(row)

    print()
    for c in others:
        vals = np.array(ratios[c])
        if vals.size < 2:
            continue
        spread = vals.max() - vals.min()
        verdict = ("CONSTANT, fingerprint shape preserved" if spread < 0.02
                   else "NOT constant, shapes differ")
        print(f"  {c}/{base}: mean {vals.mean():.4f}, stdev {vals.std():.5f}, "
              f"spread {spread:.5f}  ->  {verdict}")

    if "10" in have:
        print("\nInterpretation")
        print("  If 9/7B and 10/7B are both constant, three below-cloud cases with")
        print("  two different spectra and two different vertical distributions share")
        print("  one fingerprint. The below-cloud response is then a property of the")
        print("  atmosphere, with the absorber setting only its amplitude.")
        print("  The 7A column is the control: an above-cloud case should NOT be")
        print("  a constant multiple of the below-cloud pattern.")


if __name__ == "__main__":
    main()
