#!/usr/bin/env python3
"""
Candidate A tested against the published bulk-absorbance requirement.
Project roll number: 25m0328

THE QUESTION
Jiang et al. (2024) proposed rhomboclase plus acid ferric sulfate as the Venus
unknown ultraviolet absorber, on the strength of laboratory spectra. They did
not test whether those minerals absorb strongly enough, because at the time
there was no published number to test against.

Spacek et al. (2026) supplied one. From Venus reflectance observations and a
multiple-scattering radiative transfer model, they derive the decadic
absorption coefficient the liquid in the cloud droplets must have:
a(375 nm) = 1278 cm-1. They apply it to ferric chloride and find it needs about
1.3 kg/L, which they call physically implausible. They cite Jiang but do not
apply their own test to those minerals.

This script does. It uses Jiang's own measured extinction coefficients and
Spacek's own linear extrapolation method.

THE ABSORBANCE CONVENTION
Spacek define a(lambda) = alpha(lambda) / ln(10), and state that the decadic
absorbance A = a*l is "the quantity measured in the laboratory by UV-visible
absorbance spectroscopy". Jiang measured in a 1 mm quartz cuvette on a
spectrophotometer, so their coefficients are decadic on the same footing.
The script reports the other convention as well, because the verdict should
not rest on a reading of someone else's notation. It does not: both readings
exceed the figure Spacek reject.

THE SELF-TEST
Before touching the minerals, the script reproduces Spacek's own ferric
chloride figure from their stated input. If that fails, the unit chain is
broken and nothing after it can be trusted, so the script stops.

WHY THE RESULT IS CONSERVATIVE
Jiang's coefficient is an EXTINCTION coefficient measured on a suspension, so
it contains scattering as well as absorption. Removing the scattering can only
lower it, which raises the required loading. The shortfall reported here is a
lower bound on the true shortfall.

Run from anywhere:
    python Experiment13/candidate_A_bulk_absorbance_25m0328.py
"""

import os
import sys

import numpy as np

# ---- Spacek et al. 2026, Astrobiology 26(9), 727-739 ------------------------
A_REQUIRED = 1278.0        # cm-1, decadic, bulk liquid of the cloud aerosol
LAMBDA_NM = 375.0          # wavelength at which the requirement peaks
MODEL_RANGE = (365.0, 455.0)
FECL3_A_PER_GPL = 1.0      # their stated a(375) for a 1 g/L aqueous solution
FECL3_STATED_KGPL = 1.3    # the loading they print
FECL3_DENSITY = 2.90       # g/cm3

# ---- Jiang et al. 2024, Science Advances 10, eadg8826, supplementary --------
RHO_SOLUTION = 1.8         # g/cm3, the sulfuric acid solution
MINERALS = {
    "rhomboclase":         dict(rho=2.23, loading=0.0100, label="rhomboclase"),
    "acid-ferric-sulfate": dict(rho=2.80, loading=0.0125, label="acid ferric sulfate"),
}
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def loading_needed(eps):
    """kg/L of material needed to reach the requirement. 1 g/cm3 is 1 kg/L, so
    the required mass per unit volume in g/cm3 IS the loading in kg/L."""
    return A_REQUIRED / eps


def main():
    print("=" * 78)
    print("CANDIDATE A AGAINST THE PUBLISHED BULK-ABSORBANCE REQUIREMENT")
    print("=" * 78)
    print(f"requirement : a({LAMBDA_NM:.0f} nm) = {A_REQUIRED:.0f} cm-1, decadic, bulk liquid")
    print(f"model range : {MODEL_RANGE[0]:.0f} to {MODEL_RANGE[1]:.0f} nm")
    print("source      : Spacek, Rimmer, Petkowski and Lee 2026, Astrobiology 26(9), 727")

    # ---------------------------------------------------------------- self-test
    eps_fecl3 = FECL3_A_PER_GPL / 0.001          # cm2/g, from 1 g/L = 0.001 g/cm3
    got = loading_needed(eps_fecl3)
    print("\nSELF-TEST: reproduce their ferric chloride figure")
    print(f"  stated   : 1 g/L gives a({LAMBDA_NM:.0f}) ~ {FECL3_A_PER_GPL:.0f} cm-1, "
          f"so eps = {eps_fecl3:.0f} cm2/g")
    print(f"  computed : {got:.3f} kg/L      printed in the paper: ~{FECL3_STATED_KGPL} kg/L")
    if abs(got - FECL3_STATED_KGPL) > 0.05:
        sys.exit(f"  SELF-TEST FAILED: {got:.3f} against {FECL3_STATED_KGPL}. "
                 f"The unit chain is wrong; nothing below can be trusted.")
    print("  passed\n")

    rows = []
    for key, spec in MINERALS.items():
        path = os.path.join(DATA, f"jiang_{key}_extracted.dat")
        if not os.path.exists(path):
            sys.exit(f"missing {path}. Run extract_jiang_spectra_25m0328.py first.")
        lam, eps = np.loadtxt(path).T
        e = float(np.interp(LAMBDA_NM, lam, eps))
        e_peak = float(eps.max())
        l_peak = float(lam[int(np.argmax(eps))])
        rows.append(dict(label=spec["label"], rho=spec["rho"], loading=spec["loading"],
                         eps=e, eps_peak=e_peak, lam_peak=l_peak))

    print("=" * 78)
    print(f"AT {LAMBDA_NM:.0f} nm, WHERE THE REQUIREMENT PEAKS")
    print("=" * 78)
    h = (f"{'material':22s} {'eps':>9s} {'vs FeCl3':>9s} {'needed':>10s} "
         f"{'% of droplet':>13s} {'pure solid':>11s} {'verdict':>16s}")
    print(h); print("-" * len(h))
    for r in rows:
        need = loading_needed(r["eps"])
        pure = r["eps"] * r["rho"]
        verdict = ("short by %.2fx" % (A_REQUIRED / pure)) if pure < A_REQUIRED else "reaches it"
        print(f"{r['label']:22s} {r['eps']:8.1f} {eps_fecl3/r['eps']:8.2f}x "
              f"{need:7.2f}kg/L {100*need/RHO_SOLUTION:12.0f}% {pure:10.0f} {verdict:>16s}")
    pure = eps_fecl3 * FECL3_DENSITY
    print(f"{'ferric chloride':22s} {eps_fecl3:8.0f} {1.0:8.2f}x "
          f"{loading_needed(eps_fecl3):7.2f}kg/L "
          f"{100*loading_needed(eps_fecl3)/RHO_SOLUTION:12.0f}% {pure:10.0f} "
          f"{'reaches it':>16s}")
    print("\n  'pure solid' is the coefficient a droplet of pure material would have,")
    print("  using the solid densities Jiang give. '% of droplet' uses their stated")
    print(f"  solution density of {RHO_SOLUTION} g/cm3. Above 100 per cent means no droplet")
    print("  composition reaches the requirement.")
    print("  Ferric chloride is shown for scale only. Spacek already reject its")
    print(f"  {FECL3_STATED_KGPL} kg/L as physically implausible.")

    print("\n" + "=" * 78)
    print("AT EACH MINERAL'S OWN PEAK, THE MOST FAVOURABLE CASE POSSIBLE")
    print("=" * 78)
    for r in rows:
        need = loading_needed(r["eps_peak"])
        note = "" if MODEL_RANGE[0] <= r["lam_peak"] <= MODEL_RANGE[1] \
            else "   (below the modelled range, so the requirement there is unpublished)"
        print(f"  {r['label']:22s} peak {r['eps_peak']:6.1f} cm2/g at {r['lam_peak']:5.1f} nm"
              f"  ->  {need:5.2f} kg/L, {100*need/RHO_SOLUTION:3.0f}% of the droplet{note}")
    print(f"\n  Both still exceed the {FECL3_STATED_KGPL} kg/L the same paper rejects.")

    print("\n" + "=" * 78)
    print("ROBUSTNESS TO THE ABSORBANCE CONVENTION")
    print("=" * 78)
    for r in rows:
        dec = loading_needed(r["eps"])
        nap = loading_needed(r["eps"] / np.log(10))
        print(f"  {r['label']:22s} decadic {dec:5.2f} kg/L    "
              f"if Jiang's tau were napierian {nap:6.2f} kg/L")
    print(f"\n  The verdict does not depend on the convention. Only the margin does.")

    print("\n" + "=" * 78)
    print("CONCLUSION")
    print("=" * 78)
    worst = max(rows, key=lambda r: loading_needed(r["eps"]))
    best = min(rows, key=lambda r: loading_needed(r["eps"]))
    print(f"  Neither mineral reaches {A_REQUIRED:.0f} cm-1 even as a pure solid. The better")
    print(f"  of the two, {best['label']}, would need {loading_needed(best['eps']):.2f} kg/L, which is")
    print(f"  {100*loading_needed(best['eps'])/RHO_SOLUTION:.0f} per cent of the droplet mass, and "
          f"{loading_needed(best['eps'])/FECL3_STATED_KGPL:.1f} times the loading")
    print(f"  that Spacek et al. already reject for ferric chloride.")
    print("\n  This constrains the minerals as bulk absorbers in the droplets. It says")
    print("  nothing about whether iron-sulfur chemistry occurs on Venus, and it rests")
    print(f"  entirely on the {A_REQUIRED:.0f} cm-1 figure, which is five weeks old and has not")
    print("  been independently checked.")


if __name__ == "__main__":
    main()
