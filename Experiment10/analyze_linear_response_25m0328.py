#!/usr/bin/env python3
"""
Diagnosis of the below-cloud response scaling, and the limits of Experiment 10.
Project roll number: 25m0328

Experiments 7B, 9 and 10 place the absorber below the photochemically active
region and produce sulfur photolysis responses whose channel by channel ratios
are constant. This script establishes what that constancy does and does not
mean.

Three findings, each reproduced below.

1. Experiments 9 and 10 do not independently test spectral shape.
   Both use the ClS2 vertical profile, which retains essentially no absorber
   above 69 km. Both spectra are normalised to the same band-integrated Q.
   The retained optical depth above 69 km is therefore identical for the two
   cases, and the perturbation inside the response band is identical by
   construction. The agreement is a property of the experiment design, not a
   physical discovery. The two runs do differ substantially lower down, where
   the absorber actually sits.

2. The constant ratio is linear response, not a universal fingerprint.
   The channel ratio between cases equals the ratio of the local
   band-integrated actinic flux change at the altitude of maximum response.
   Every channel responds linearly to the same local flux change, so every
   channel ratio takes the same value. This is expected behaviour, and the
   spread across channels measures how well linearity holds.

3. Linearity fails when the absorber is concentrated in the band.
   Experiment 7A raises the retained optical depth above 69 km by almost an
   order of magnitude. There the response is nonlinear and channel specific,
   and the ratios scatter by more than a factor of ten. Results may be scaled
   within the linear regime and must not be scaled across it.

Run from the repository root. Read-only, no VULCAN rerun.
"""

import os
import pickle
import numpy as np
from scipy import interpolate

AREA_CM2 = np.pi * (0.98e-4) ** 2 / 4.0     # pi D^2 / 4 with D = D_1
BAND_LO_NM, BAND_HI_NM = 300.0, 620.0
Z_BAND_LO, Z_BAND_HI = 69.0, 93.0

REF_OLD = "output/Experiment2_A1_UUV-nominal-25m0328.vul"
REF_NEW = "output/Experiment8_baseline_A1-25m0328.vul"

CASES = {
    "7A": (REF_OLD, "output/Experiment7A_UUV-upward18-25m0328.vul",
           "atm/mode1+2_Experiment7A_shift18_25m0328.txt", "atm/UV_absorber.txt"),
    "7B": (REF_OLD, "output/Experiment7B_UUV-downward18-25m0328.vul",
           "atm/mode1+2_Experiment7B_shiftDown18_25m0328.txt", "atm/UV_absorber.txt"),
    "9":  (REF_NEW, "output/Experiment9_ClS2-selfconsistent-25m0328.vul",
           "atm/mode1+2_Experiment9_ClS2profile_25m0328.txt",
           "atm/UV_absorber_Experiment8A_ClS2_25m0328.txt"),
    "10": (REF_NEW, "output/Experiment10_S3-on-ClS2profile-25m0328.vul",
           "atm/mode1+2_Experiment9_ClS2profile_25m0328.txt",
           "atm/UV_absorber_Experiment8C_S3_25m0328.txt"),
}

PARENTS = ["S2O2", "SCl2", "ClS2", "S3", "S2Cl2", "SCl", "S2O", "S4"]


def load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def main():
    base = load(REF_NEW)
    z = np.clip(np.asarray(base["atm"]["zmco"]) / 1e5, 0.0, None)
    dz = np.asarray(base["atm"]["dz"])
    bins = np.asarray(base["variable"]["bins"])
    quad = np.where(bins < 240.0, 0.1, 2.0)
    band_w = (bins >= BAND_LO_NM) & (bins <= BAND_HI_NM)
    band_z = (z >= Z_BAND_LO) & (z <= Z_BAND_HI)
    aloft = z >= Z_BAND_LO

    def retained_tau(profile_path, spectrum_path):
        p = np.loadtxt(profile_path, skiprows=1)
        n1 = interpolate.interp1d(p[:, 0], p[:, 1], kind="slinear")(z)
        u = np.loadtxt(spectrum_path, skiprows=1)
        g = interpolate.interp1d(u[:, 0], u[:, 1], kind="slinear",
                                 bounds_error=False, fill_value=0.0)
        q = np.where((bins >= u[0, 0]) & (bins <= u[-1, 0]),
                     g(np.clip(bins, u[0, 0], u[-1, 0])), 0.0)
        per_lambda = np.sum((q[None, :] * AREA_CM2 * n1[:, None] * dz[:, None])[aloft], axis=0)
        return float(np.sum(per_lambda * quad * band_w))

    print("FINDING 1: retained UUV optical depth above "
          f"{Z_BAND_LO:.0f} km, band integrated\n")
    tau = {"A1": retained_tau("atm/mode1+2.txt", "atm/UV_absorber.txt")}
    for label, (_, _, prof, spec) in CASES.items():
        if os.path.exists(prof) and os.path.exists(spec):
            tau[label] = retained_tau(prof, spec)
    for label, value in tau.items():
        print(f"  {label:3s} {value:14.8e}")
    if "9" in tau and "10" in tau:
        same = "IDENTICAL" if tau["9"] == tau["10"] else f"differ by {abs(tau['10']-tau['9']):.2e}"
        print(f"\n  Cases 9 and 10: {same}. Both use the ClS2 profile, which retains")
        print("  essentially no absorber above the band, and both spectra carry the same")
        print("  band-integrated opacity. Their in-band perturbation is fixed by design.")

    print("\n\nFINDING 2: the channel ratio equals the local flux-change ratio\n")
    results = {}
    for label, (ref_p, case_p, _, _) in CASES.items():
        if not (os.path.exists(ref_p) and os.path.exists(case_p)):
            continue
        ref, case = load(ref_p), load(case_p)
        f0 = np.asarray(ref["variable"]["aflux"])
        d = np.asarray(case["variable"]["aflux"]) - f0
        fp = {}
        for parent in PARENTS:
            key = (parent, 0)
            if key not in ref["variable"]["J_sp"]:
                continue
            a = np.asarray(ref["variable"]["J_sp"][key])
            dj = np.asarray(case["variable"]["J_sp"][key]) - a
            j = int(np.argmax(np.abs(dj * band_z)))
            fp[parent] = (100.0 * dj[j] / a[j] if a[j] else np.nan, j)
        j81 = int(np.argmin(np.abs(z - 81.0)))
        results[label] = (fp, float(np.sum(d[j81] * quad * band_w)))

    if "7B" in results and "9" in results:
        fp_b, flux_b = results["7B"]
        for label in ("9", "10"):
            if label not in results:
                continue
            fp_c, flux_c = results[label]
            ratios = [fp_c[p][0] / fp_b[p][0] for p in PARENTS
                      if p in fp_c and p in fp_b and fp_b[p][0]]
            r = np.array(ratios)
            print(f"  {label} against 7B")
            print(f"    channel ratios : mean {r.mean():.4f}, spread {r.max()-r.min():.5f}")
            print(f"    flux ratio 81km: {flux_c / flux_b:.4f}")
            print(f"    agreement      : {abs(r.mean() - flux_c / flux_b):.2e}\n")

    print("\nFINDING 3: the same test on an in-band case\n")
    if "7A" in results and "7B" in results:
        fp_a, flux_a = results["7A"]
        fp_b, flux_b = results["7B"]
        r = np.array([fp_a[p][0] / fp_b[p][0] for p in PARENTS
                      if p in fp_a and p in fp_b and fp_b[p][0]])
        print(f"  7A against 7B: channel ratios mean {r.mean():.4f}, "
              f"spread {r.max()-r.min():.4f}")
        print(f"  retained optical depth aloft: A1 {tau['A1']:.3f} -> 7A {tau['7A']:.3f}")
        print("  The perturbation is large, the response is nonlinear, and the channel")
        print("  ratios scatter. Scaling is not valid across this regime.")


if __name__ == "__main__":
    main()
