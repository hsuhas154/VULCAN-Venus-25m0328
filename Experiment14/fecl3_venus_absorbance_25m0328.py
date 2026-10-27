#!/usr/bin/env python3
"""Experiment 14: does the published rejection of ferric chloride survive when
the absorbance is taken in the right solvent?

THE QUESTION
------------
Spacek et al. (2026, Astrobiology 26(9), 727) derive the decadic absorption
coefficient a(375) = 1278 cm-1 that the bulk liquid of a Venus cloud droplet must
have. They apply it to ferric chloride and reject it:

    "Our model predicts a375 = 1278 cm-1 for the venusian cloud liquid, while an
     aqueous solution of 1 g/L FeCl3 has a375 ~ 1 cm-1. Simple linear
     extrapolation would therefore require approximately 1.3 kg/L FeCl3."

1.3 kg/L is 71 per cent of the droplet mass at the droplet density of 1.8 g/cm3
used by Dai et al. (2023), a loading they call physically implausible.

The input to that extrapolation is an AQUEOUS measurement. Venus cloud droplets
are roughly 75 wt% sulfuric acid, and Egan et al. (2025) report that the spectrum
of FeCl3 in concentrated sulfuric acid "does not appear to have been measured
previously" and that "the choice of solvent has a significant effect on the
measured spectrum". Their Figure 5a gives the molar absorptivity of ferric
chloride in hydrochloric acid at five strengths, which they use as the ferric
chloride components when they fit their mixed sulfuric and hydrochloric acid
samples. This experiment substitutes those curves for the aqueous value and
recomputes. The curves are not a measurement in sulfuric acid, so the result is
conditional: it holds if ferric chloride absorbs in the droplets as it does in
hydrochloric acid.

WHAT IS AND IS NOT CONCLUDED
----------------------------
This recomputes one number in someone else's argument. It does not identify the
absorber, does not test FeCl3 against the altitude criterion of Experiment 9
(FeCl3 has no chemistry in the VULCAN network, so the model cannot place it), and
does not adjudicate Spacek's requirement itself, which Egan et al. (2026) put in
question by a separate route.

There is no single answer because there is no single FeCl3. The iron speciates
into chloro complexes set by chloride activity, so the result is reported as a
range over the five measured HCl concentrations.

Run from this directory, after the extraction script:
    python3 fecl3_venus_absorbance_25m0328.py
"""
import os, sys, glob, json
import numpy as np

A_REQUIRED = 1278.0        # cm-1, decadic, at 375 nm (Spacek et al. 2026)
LAMBDA_NM = 375.0
RHO_SOLUTION = 1.8         # g cm-3, droplet density used by Dai et al. (2023)
M_FECL3 = 162.204          # g mol-1
SPACEK_AQUEOUS_A_PER_GPL = 1.0     # their stated a(375) for 1 g/L aqueous
SPACEK_LOADING_KGPL = 1.3          # the figure they print and reject
EGAN2026_FIT_WT_PCT = (1.0, 1.2)   # FeCl3 in mode 1 droplets, Egan et al. (2026)

here = os.path.dirname(os.path.abspath(__file__))
files = sorted(glob.glob(os.path.join(here, "data", "fecl3_Cl*_eps.dat")))
if not files:
    sys.exit("No extracted curves found. Run extract_egan_fig5a_25m0328.py first.")

eps_mass_aqueous = SPACEK_AQUEOUS_A_PER_GPL / 0.001      # cm2/g = 1000
loading_aqueous = A_REQUIRED / eps_mass_aqueous           # g/cm3 = kg/L

# ---------------------------------------------------------------- self-test
# Reproduce the figure Spacek print, from their own stated input. If this does
# not come back at their 1.3 kg/L the arithmetic below cannot be trusted.
if abs(loading_aqueous - SPACEK_LOADING_KGPL) > 0.05:
    sys.exit(f"SELF-TEST FAILED: aqueous route gives {loading_aqueous:.3f} kg/L "
             f"against their {SPACEK_LOADING_KGPL} kg/L.")
print("SELF-TEST: reproducing the published aqueous figure")
print(f"  their input : 1 g/L aqueous FeCl3 gives a(375) ~ "
      f"{SPACEK_AQUEOUS_A_PER_GPL:.0f} cm-1, so eps = {eps_mass_aqueous:.0f} cm2/g")
print(f"  computed    : {loading_aqueous:.3f} kg/L     printed in the paper: "
      f"~{SPACEK_LOADING_KGPL} kg/L     PASS")

W = 104
print("\n" + "=" * W)
print("FeCl3 REQUIRED LOADING, USING THE ABSORBANCE MEASURED IN CONCENTRATED SULFURIC ACID")
print("=" * W)
print(f"requirement a({LAMBDA_NM:.0f}) = {A_REQUIRED:.0f} cm-1 decadic;  "
      f"droplet solution density {RHO_SOLUTION} g/cm3;  FeCl3 M = {M_FECL3} g/mol\n")
hdr = (f"{'curve':8s}{'HCl / M':>9s}{'eps_Fe':>10s}{'eps_mass':>11s}{'vs aqueous':>12s}"
       f"{'needed':>11s}{'% droplet':>11s}  verdict")
print(hdr); print("-" * len(hdr))

rows, out = [], {}
for f in files:
    k = os.path.basename(f).split("Cl")[1][0]
    hcl = float([l for l in open(f) if "HCl molarity" in l][0].split()[-1])
    d = np.loadtxt(f)
    eps = float(np.interp(LAMBDA_NM, d[:, 0], d[:, 1]))     # L mol-1 cm-1
    eps_mass = 1000.0 * eps / M_FECL3                       # cm2/g
    need = A_REQUIRED / eps_mass                            # kg/L
    pct = 100.0 * need / RHO_SOLUTION
    verdict = "exceeds droplet" if pct > 100 else ("major component" if pct > 20
              else "plausible minor component")
    rows.append((k, hcl, eps, eps_mass, need, pct))
    out[f"Cl{k}"] = dict(hcl_M=hcl, eps_molar=eps, eps_mass=eps_mass,
                         ratio_to_aqueous=eps_mass / eps_mass_aqueous,
                         loading_kg_per_L=need, pct_of_droplet=pct)
    print(f"Cl {k:5s}{hcl:9.2f}{eps:10.0f}{eps_mass:11.0f}{eps_mass/eps_mass_aqueous:11.1f}x"
          f"{need:11.3f}{pct:10.1f}%  {verdict}")

print(f"\n{'aqueous (Spacek)':17s}{'':9s}{162.2:10.0f}{eps_mass_aqueous:11.0f}"
      f"{1.0:11.1f}x{loading_aqueous:11.3f}{100*loading_aqueous/RHO_SOLUTION:10.1f}%  "
      f"rejected as implausible")

lo = min(r[4] for r in rows); hi = max(r[4] for r in rows)
plo, phi = 100*lo/RHO_SOLUTION, 100*hi/RHO_SOLUTION
print("\n" + "=" * W)
print("RESULT")
print("=" * W)
flo = min(r[3] for r in rows) / eps_mass_aqueous
fhi = max(r[3] for r in rows) / eps_mass_aqueous
print(f"  Measured in concentrated sulfuric acid, FeCl3 is {flo:.1f} to {fhi:.1f} times "
      f"more absorbing at {LAMBDA_NM:.0f} nm than the aqueous value used to reject it.")
print(f"  The required loading falls from {loading_aqueous:.2f} kg/L "
      f"({100*loading_aqueous/RHO_SOLUTION:.0f} per cent of the droplet) to "
      f"{lo:.3f} to {hi:.3f} kg/L ({plo:.1f} to {phi:.1f} per cent).")
print(f"  No measured curve requires FeCl3 to be more than {phi:.0f} per cent of the")
print(f"  droplet, and the most chloride-rich requires {plo:.1f} per cent.")
print()
print("  The rejection does not survive the change of solvent. It was an artefact")
print("  of extrapolating an aqueous absorbance to a sulfuric acid droplet.")
print()
fit_mid = sum(EGAN2026_FIT_WT_PCT) / 2.0
factor_after = plo / fit_mid
factor_before = 100 * loading_aqueous / RHO_SOLUTION / fit_mid
print(f"  For comparison, Egan et al. (2026) independently fit {EGAN2026_FIT_WT_PCT[0]:.1f} to "
      f"{EGAN2026_FIT_WT_PCT[1]:.1f} wt% FeCl3")
print(f"  in mode 1 droplets to MESSENGER reflectance. The most chloride-rich curve")
print(f"  here needs {plo:.1f} per cent, within a factor of {factor_after:.1f} of that fit,")
print(f"  against a factor of {factor_before:.0f} before the correction. The residue is the")
print("  difference between a bulk-average model and a multiple-scattering")
print("  calculation with the absorber confined to mode 1, which is precisely the")
print("  assumption Spacek et al. flag as needing further study.")

with open(os.path.join(here, "fecl3_result_25m0328.json"), "w") as fh:
    json.dump(dict(requirement_cm1=A_REQUIRED, lambda_nm=LAMBDA_NM,
                   rho_solution=RHO_SOLUTION, aqueous_loading_kg_per_L=loading_aqueous,
                   aqueous_eps_mass=eps_mass_aqueous,
                   aqueous_pct_of_droplet=100 * loading_aqueous / RHO_SOLUTION,
                   range_kg_per_L=[lo, hi], range_pct_droplet=[plo, phi],
                   ratio_lo=flo, ratio_hi=fhi,
                   egan2026_fit_mid_wt_pct=fit_mid,
                   factor_to_fit_before=factor_before, factor_to_fit_after=factor_after,
                   per_curve=out), fh, indent=1, sort_keys=True)
print(f"\nresult written to fecl3_result_25m0328.json")
