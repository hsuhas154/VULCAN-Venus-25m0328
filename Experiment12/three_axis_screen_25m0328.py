#!/usr/bin/env python3
"""
Three-axis admissibility screen over every species the network can represent.
Project roll number: 25m0328

Experiments 8 through 11 evaluated candidate absorbers on three criteria but
never all three over the same set of species:

  opacity   Experiment 8 screened 12 hand-picked species for the optical depth
            their converged column produces above 58 km.
  spectrum  Experiments 8 and 11 substituted 5 spectra at matched opacity.
  altitude  Experiment 9 tested exactly 1 species at its own distribution and
            showed that this is the axis that decides the outcome.

This screen evaluates all three for every species that appears in the network
and holds a stored cross section with support in the operating range. It is
read-only: no VULCAN rerun, no new inputs, no configuration change.

Two design decisions carry the result and are stated here rather than buried.

1. THE COMPARISON IS MADE ABOVE THE RESPONSE BAND, NOT ABOVE THE CLOUD BASE.
   Experiment 9 established that an absorber outside 69 to 93 km is radiatively
   invisible. A species is therefore credited only with the column it holds
   above 69 km. The column above 58 km is reported alongside so that the two
   failure modes, too little material and material in the wrong place, can be
   told apart.

2. THE COMPARISON IS MADE WAVELENGTH BY WAVELENGTH, NOT AT THE WINDOW PEAK.
   Ranking on the peak inside 300 to 620 nm promotes species whose absorption
   is collapsing from below 300 nm. SO2 is the case in point: its cross section
   falls from 5.44e-19 at 300 nm to 5.15e-24 at 400 nm, four orders of
   magnitude, which is precisely why an unknown absorber is needed above
   320 nm at all. Each species is therefore compared against the absorber's own
   optical depth at the same wavelength and over the same altitude range.

EXCLUSIONS. The species V and V2 are excluded and the reason is recorded in
the output. They are inert particle tracers, not molecules: their only network
entries are the null reactions V -> V and V2 -> V2, their columns track the
mode-1 and mode-2 cloud droplet profiles, and their cross sections of order
1e-8 to 1e-7 cm2 are geometric areas for micron-sized spheres. Those cross
sections vary by under 15 per cent across the whole window, which is grey
extinction. Including them would return the result that Venus has optically
thick clouds, which is not a finding about a near-ultraviolet absorber.

OBSERVATIONAL STANDING OF THE SURVIVING CANDIDATES. A species can pass all three
axes here and still be ruled out by measurement, so the screen's output is a
shortlist and not a result. NO2 is the case that needs stating explicitly,
because its screen position overstates its standing. Mahieux et al. (2024,
Icarus 409, 115862) place an upper limit on Venus NO2, but that limit is derived
at the terminator from solar occultation. NO2 photolyses rapidly in daylight, so
the dayside abundance in the cloud tops, which is the quantity this screen needs,
is lower than the terminator limit rather than bounded by it. The screen
therefore credits NO2 with a column that daytime photochemistry does not permit,
and no conclusion about NO2 should be drawn from its ranking here without a
dayside constraint.

Run from the repository root:
    python Experiment12/three_axis_screen_25m0328.py
"""

import os
import sys
import pickle
import numpy as np
from scipy import interpolate

NOMINAL = "output/Experiment2_A1_UUV-nominal-25m0328.vul"
PROFILE = "atm/mode1+2.txt"
SPECTRUM = "atm/UV_absorber.txt"

TAU_REQUIRED = 1.094968472081      # absorber optical depth at 320 nm, whole column
D_UV = 0.98e-4                     # cm, mode-1 absorber particle diameter
AREA = np.pi * D_UV ** 2 / 4.0     # cm2
Z_BASE_KM = 58.0                   # cloud base
BAND_Z = (69.0, 93.0)              # radiative response band, Experiment 6
SCREEN_NM = (300.0, 620.0)         # range in which the absorber operates
GRID_NM = (300.0, 698.0)           # full absorber wavelength grid
PEAK_NM = (320.0, 400.0)           # observed near-ultraviolet enhancement
PROBE_NM = [320, 340, 360, 380, 400, 440, 480]
RANK_NM = 360                      # middle of the enhancement
EXCLUDED = {"V": "inert mode-1 particle tracer", "V2": "inert mode-2 particle tracer"}


def main():
    if not (os.path.isdir("output") and os.path.isdir("atm")):
        sys.exit("ERROR: run this from the repository root, where output/ and atm/ live.")

    with open(NOMINAL, "rb") as f:
        sol = pickle.load(f)
    var, atm = sol["variable"], sol["atm"]

    z = np.asarray(atm["zmco"]) / 1e5
    dz = np.asarray(atm["dz"])
    zu = np.asarray(atm["zco"])[1:] / 1e5
    bins = np.asarray(var["bins"])
    y = np.asarray(var["y"])
    species = list(var["species"])
    cross = var["cross"]

    above = z >= Z_BASE_KM
    aloft = z >= BAND_Z[0]
    in_z = (z >= BAND_Z[0]) & (z <= BAND_Z[1])
    win = (bins >= SCREEN_NM[0]) & (bins <= SCREEN_NM[1])

    # the absorber's own optical depth, above the response band, wavelength by wavelength
    p = np.loadtxt(PROFILE, skiprows=1)
    n_uv = interpolate.interp1d(p[:, 0], p[:, 1], kind="slinear")(zu)
    col_uv = float(np.sum(n_uv[zu >= BAND_Z[0]] * dz[zu >= BAND_Z[0]]))
    uv = np.loadtxt(SPECTRUM, skiprows=1)
    q = interpolate.interp1d(uv[:, 0], uv[:, 1], kind="slinear",
                             bounds_error=False, fill_value=0.0)
    tau_req = {l: float(q(l)) * AREA * col_uv for l in PROBE_NM}

    grid = bins[(bins >= GRID_NM[0]) & (bins <= GRID_NM[1])]
    on_grid = (bins >= GRID_NM[0]) & (bins <= GRID_NM[1])
    q_emp = q(grid)
    peak_w = (grid >= PEAK_NM[0]) & (grid <= PEAK_NM[1])

    rows, nosupport, nocross = [], [], []
    for name in species:
        if name in EXCLUDED:
            continue
        if name not in cross:
            nocross.append(name)
            continue
        sigma = np.asarray(cross[name])
        if sigma.shape != bins.shape or sigma[win].max() <= 0:
            nosupport.append(name)
            continue
        n = y[:, species.index(name)]
        c58 = float(np.sum(n[above] * dz[above]))
        c69 = float(np.sum(n[aloft] * dz[aloft]))
        if c58 <= 0:
            nosupport.append(name)
            continue
        j = int(np.argmax(np.where(win, sigma * c58, 0.0)))
        sg = sigma[on_grid]
        tot = float(np.sum(sg * 2.0))
        rows.append(dict(
            name=name, lam=bins[j], sigma=sigma[j], c58=c58, c69=c69,
            centroid=float(np.sum(z[above] * n[above] * dz[above]) / c58),
            f_band=100.0 * float(np.sum(n[in_z] * dz[in_z])) / c58,
            tau58=sigma[j] * c58, tau69=sigma[j] * c69,
            corr=float(np.corrcoef(q_emp, sg)[0, 1]),
            share=100.0 * float(np.sum(sg[peak_w] * 2.0)) / tot if tot > 0 else np.nan,
            probe={l: (sigma[int(np.argmin(np.abs(bins - l)))] * c69) for l in PROBE_NM},
            cov400=bool(np.all(sigma[(bins >= 320) & (bins <= 400)] > 0)),
            cov480=bool(np.all(sigma[(bins >= 320) & (bins <= 480)] > 0)),
            sup_lo=float(bins[(sigma > 0) & win].min()),
            sup_hi=float(bins[(sigma > 0) & win].max())))

    W = 112
    print("=" * W)
    print("THREE-AXIS ADMISSIBILITY SCREEN")
    print("=" * W)
    print(f"response band          : {BAND_Z[0]:.0f} to {BAND_Z[1]:.0f} km")
    print(f"absorber column aloft  : {col_uv:.6e} cm-2")
    print(f"species screened       : {len(rows)} of {len(species)} tracked")
    print(f"excluded               : " + ", ".join(f"{k} ({v})" for k, v in EXCLUDED.items()))

    print("\n" + "=" * W)
    print("PART 1: the three axes, for every screened species")
    print("=" * W)
    print("Opacity is the peak inside the window from the column above the cloud base,")
    print(f"compared against the whole-column requirement of {TAU_REQUIRED:.6f}, as in")
    print("Experiment 8. Altitude is where that column sits. Spectrum is the shape")
    print("correlation against the empirical curve and its share inside 320 to 400 nm.")
    h = (f"{'species':8s} {'col>58km':>11s} {'centroid':>9s} {'in band':>8s} "
         f"{'peak':>6s} {'tau>58km':>11s} {'short':>10s} {'corr':>7s} {'320-400':>8s}")
    print(); print(h); print("-" * len(h))
    for r in sorted(rows, key=lambda r: -r["tau58"]):
        print(f"{r['name']:8s} {r['c58']:11.4e} {r['centroid']:8.2f}k {r['f_band']:7.3f}% "
              f"{r['lam']:5.0f}n {r['tau58']:11.4e} {TAU_REQUIRED/r['tau58']:10.4g} "
              f"{r['corr']:7.4f} {r['share']:7.2f}%")

    print("\n" + "=" * W)
    print("PART 2: the decisive test, shortfall at each wavelength from the column aloft")
    print("=" * W)
    print("Each species is credited only with the column it holds above the response band,")
    print("and compared against the absorber's own optical depth over the same range:")
    for l in PROBE_NM:
        print(f"    {l} nm   Q = {float(q(l)):.4f}   required tau above "
              f"{BAND_Z[0]:.0f} km = {tau_req[l]:.6f}")
    print()
    h = f"{'species':8s} {'support':>12s} {'cov':>4s} " + " ".join(f"{l:>11d}" for l in PROBE_NM)
    print(h); print("-" * len(h))
    ranked = sorted(rows, key=lambda r: (tau_req[RANK_NM] / r["probe"][RANK_NM])
                    if r["probe"][RANK_NM] > 0 else np.inf)
    for r in ranked:
        cells = []
        for l in PROBE_NM:
            t = r["probe"][l]
            cells.append(f"{tau_req[l]/t:11.4g}" if t > 0 else f"{'-':>11s}")
        flag = "full" if r["cov480"] else ("400" if r["cov400"] else "gap")
        print(f"{r['name']:8s} {r['sup_lo']:5.0f}-{r['sup_hi']:<6.0f} {flag:>4s} "
              + " ".join(cells))
    print(f"\nranked by shortfall at {RANK_NM} nm; a dash means no cross section there.")
    print("cov: full = continuous over 320 to 480 nm, 400 = continuous over 320 to 400 nm,")
    print("     gap = gapped or truncated inside the enhancement, so the species cannot be")
    print("     ranked across it and its absence from the top of this table means nothing.")
    full = [r["name"] for r in rows if r["cov480"]]
    part = [r["name"] for r in rows if r["cov400"] and not r["cov480"]]
    gap  = [r["name"] for r in rows if not r["cov400"]]
    print(f"\n  continuous over 320 to 480 nm : {len(full):2d}  {', '.join(sorted(full))}")
    print(f"  continuous over 320 to 400 nm : {len(part):2d}  {', '.join(sorted(part))}")
    print(f"  gapped or truncated           : {len(gap):2d}  {', '.join(sorted(gap))}")
    ranked_full = [r for r in ranked if r["cov480"]]
    if len(ranked_full) > 1:
        a, b = ranked_full[0], ranked_full[1]
        print(f"\n  Among the {len(full)} species evaluable across the whole enhancement, "
              f"{a['name']} is best,")
        print(f"  short by {tau_req[RANK_NM]/a['probe'][RANK_NM]:.4g} at {RANK_NM} nm against "
              f"{tau_req[RANK_NM]/b['probe'][RANK_NM]:.4g} for {b['name']}, "
              f"a separation of {(b['probe'][RANK_NM] and a['probe'][RANK_NM]/b['probe'][RANK_NM]):.4g}x.")

    print("\n" + "=" * W)
    print("PART 3: what the screen cannot reach")
    print("=" * W)
    print(f"{len(nosupport)} species hold a cross section with no support in "
          f"{SCREEN_NM[0]:.0f} to {SCREEN_NM[1]:.0f} nm:")
    print("  " + ", ".join(sorted(nosupport)))
    print(f"\n{len(nocross)} tracked species hold no cross section at all and cannot be "
          f"screened at any wavelength:")
    print("  " + ", ".join(sorted(nocross)))
    print("\n  This screen is exhaustive over what the network can represent, not over")
    print("  what may be present on Venus, and its reach is narrower than the species")
    print("  count suggests. Of the tracked species, most carry no cross section at all,")
    print("  and most of the remainder carry one that is gapped or truncated inside the")
    print("  320 to 480 nm enhancement. Only a handful can be evaluated across it, as")
    print("  the coverage summary in Part 2 reports. The sulfur allotropes S5 through S8")
    print("  and the chlorine-sulfur species OSCl and HSCl carry no cross section, so no")
    print("  statement about them follows from this screen at any wavelength.")

    best = ranked[0]
    print("\n" + "=" * W)
    print(f"Best at {RANK_NM} nm from the column aloft: {best['name']}, short by "
          f"{tau_req[RANK_NM]/best['probe'][RANK_NM]:.4g}")
    print("=" * W)


if __name__ == "__main__":
    main()
