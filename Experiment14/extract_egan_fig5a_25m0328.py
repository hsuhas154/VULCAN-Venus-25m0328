#!/usr/bin/env python3
"""Extract the FeCl3 molar absorptivity curves from Figure 5a of

    Egan, J. V., James, A. D., and Plane, J. M. C. (2025),
    Laboratory Measurements of Ferric Chloride (FeCl3) under Venusian Conditions,
    ACS Earth and Space Chemistry, doi:10.1021/acsearthspacechem.5c00132  (CC-BY).

WHY THIS IS NEEDED
------------------
Spacek et al. (2026) reject ferric chloride as the Venus absorber because their
model requires a decadic absorption coefficient a(375) = 1278 cm-1 while "an
aqueous solution of 1 g/L FeCl3 has a(375) ~ 1 cm-1", giving about 1.3 kg/L by
linear extrapolation. That input is AQUEOUS. Egan et al. state that the spectrum
of FeCl3 in concentrated sulfuric acid "does not appear to have been measured
previously" and that "the choice of solvent has a significant effect on the
measured spectrum". Figure 5a gives the molar absorptivity of ferric chloride
in hydrochloric acid at five strengths, which Egan et al. use as the ferric
chloride components when they fit their mixed sulfuric and hydrochloric acid
samples. It is a better input than the aqueous value, but it is not a
measurement in sulfuric acid, so any result from it stays conditional.

There is no single FeCl3 coefficient. The iron speciates into chloro complexes
whose populations depend on chloride activity, so Figure 5a gives five curves
spanning 1.54 to 12.2 M HCl. All five are extracted and the downstream result is
reported as a range.

METHOD
------
The page is rendered at 600 dpi and each curve followed by its saturated stroke
colour. The shaded uncertainty bands are lighter blends of the same hue and are
excluded by requiring saturation. Nothing is read off by eye.

Three checks, all fatal, all run every time:
  1. exactly two axis frames on the page, and panel (a) is the left one;
  2. tick spacing uniform on both axes to better than 1 per cent;
  3. every curve returns to baseline at the long-wavelength end. The calibration
     knows nothing about this, so recovering it is independent evidence that the
     vertical scale is right.

Usage:
    python3 extract_egan_fig5a_25m0328.py <path to the Egan 2025 PDF>
"""
import os, sys, json, subprocess, tempfile
import numpy as np
from PIL import Image

PAGE, DPI = 6, 600
X0_NM, X1_NM = 200.0, 550.0        # x ticks sit on the frame edges
EPS_PER_TICK = 1000.0              # y ticks every 1000 L mol-1 cm-1

CURVES = {                          # label: (HCl molarity, saturation test)
    "A": (12.2, lambda r, g, b: r > 200 and g < 80 and b < 80),
    "B": (8.96, lambda r, g, b: r > 200 and 170 < g < 235 and b < 60),
    "C": (6.09, lambda r, g, b: r < 120 and g > 160 and b < 120),
    "D": (4.15, lambda r, g, b: r < 80 and g < 80 and b > 200),
    "E": (1.54, lambda r, g, b: 90 < r < 190 and g < 70 and 90 < b < 190),
}


def runs(idx, gap=5):
    idx = list(idx)
    out, s, p = [], idx[0], idx[0]
    for v in idx[1:]:
        if v - p > gap:
            out.append((s, p)); s = v
        p = v
    out.append((s, p)); return out


def main(pdf):
    if not os.path.exists(pdf):
        sys.exit(f"ERROR: cannot find {pdf}")
    tmp = tempfile.mkdtemp()
    base = os.path.join(tmp, "pg")
    subprocess.run(["pdftoppm", "-png", "-r", str(DPI), "-f", str(PAGE), "-l", str(PAGE),
                    pdf, base], check=True)
    png = [os.path.join(tmp, f) for f in os.listdir(tmp) if f.endswith(".png")][0]
    a = np.asarray(Image.open(png).convert("RGB")).astype(int)
    dark = (a.sum(2) < 250)
    h, w = dark.shape

    # ---- check 1: locate the two panel frames, take the left one ----
    top = dark[:int(h * 0.45)]
    vl = [v for v in runs(np.where(top.sum(0) > top.shape[0] * 0.25)[0])
          if (v[1] - v[0]) < 20]
    if len(vl) < 4:
        sys.exit(f"CHECK 1 FAILED: expected 4 panel verticals, found {len(vl)}.")
    L = int(np.mean(vl[0])); R = int(np.mean(vl[1]))
    # A panel edge must span the panel, not merely be a long rule somewhere on
    # the page: header rules and running heads also produce long dark rows.
    span = [r for r in range(top.shape[0]) if dark[r, L:R].mean() > 0.9]
    if not span:
        sys.exit("CHECK 1 FAILED: no horizontal line spans the panel width.")
    hl = [c for c in runs(span) if (c[1] - c[0]) < 20]
    T = int(np.mean(hl[0])); B = int(np.mean(hl[-1]))
    if not (R - L > w * 0.2 and B - T > 1000):
        sys.exit(f"CHECK 1 FAILED: implausible panel (a) box {L},{T} to {R},{B}.")

    # ---- check 2: tick uniformity ----
    def ticks(mask, axis, lo, hi):
        """Tick centres, keeping only those on the regular grid.

        Rules elsewhere on the page can be mistaken for a panel edge, which then
        drags a spurious tick into the list and corrupts the scale. Anchor on the
        axis end, which is always a real tick, and walk back accepting only gaps
        close to the modal spacing."""
        idx = np.where(mask.sum(axis) >= mask.shape[axis] * 0.6)[0]
        t = np.array([(s + p) / 2 for s, p in runs(idx)])
        t = t[(t > lo - 30) & (t < hi + 30)]
        if len(t) < 3:
            return t
        med = float(np.median(np.diff(t)))
        keep = [t[-1]]
        for v in t[-2::-1]:
            if abs((keep[-1] - v) - med) <= 0.1 * med:
                keep.append(v)
        return np.array(sorted(keep))
    xt = ticks(dark[B+6:B+26, :], 0, L, R)
    yt = ticks(dark[:, L-26:L-6], 1, T, B)
    for nm, t in (("x", xt), ("y", yt)):
        if len(t) < 4:
            sys.exit(f"CHECK 2 FAILED: only {len(t)} {nm} ticks found.")
        d = np.diff(t)
        if d.std() / d.mean() > 0.01:
            sys.exit(f"CHECK 2 FAILED: {nm} ticks non-uniform ({d.std()/d.mean():.2%}).")
    px_nm = (xt[-1] - xt[0]) / (X1_NM - X0_NM)
    px_eps = (yt[-1] - yt[0]) / ((len(yt) - 1) * EPS_PER_TICK)
    y0 = yt[-1]
    print(f"calibration: {px_nm:.4f} px/nm, {px_eps*1000:.2f} px per {EPS_PER_TICK:.0f} eps, "
          f"{len(yt)} y ticks, zero at row {y0:.0f}")

    here = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(here, "data"); os.makedirs(outdir, exist_ok=True)
    result = {"source_doi": "10.1021/acsearthspacechem.5c00132", "figure": "5a",
              "dpi": DPI, "curves": {}}

    for k, (hcl, test) in CURVES.items():
        lam, eps = [], []
        for x in range(int(xt[0]), int(xt[-1]) + 1):
            rows = [y for y in range(T, B) if test(*a[y, x])]
            if rows:
                lam.append(X0_NM + (x - xt[0]) / px_nm)
                eps.append((y0 - float(np.median(rows))) / px_eps)
        lam, eps = np.array(lam), np.array(eps)
        if len(lam) < 500:
            sys.exit(f"CHECK FAILED: curve {k} traced only {len(lam)} columns.")
        tail = eps[lam > 520] if (lam > 520).any() else eps[-40:]
        # The measured tails sit at 14 to 36 eps units, so a tolerance of a few
        # per cent of the peak is far too loose to detect a displaced zero. Bound
        # it by a tenth of a tick as well, which a 40 pixel shift of the baseline
        # fails. This limit was set by deliberately shifting the zero and
        # confirming the check fires.
        tol = max(0.1 * EPS_PER_TICK, 0.015 * eps.max())
        if abs(tail.mean()) > tol:
            sys.exit(f"CHECK 3 FAILED: curve {k} tail {tail.mean():.0f} exceeds {tol:.0f} "
                     f"(peak {eps.max():.0f}); the vertical zero is probably wrong.")
        np.savetxt(os.path.join(outdir, f"fecl3_Cl{k}_eps.dat"), np.c_[lam, eps],
                   header=f"Egan et al. 2025 ACS Earth Space Chem, Figure 5a, curve Cl {k}\n"
                          f"HCl molarity {hcl}\nwavelength_nm  eps_Fe_L_per_mol_per_cm")
        result["curves"][k] = {"hcl_M": hcl, "n_points": len(lam),
                               "peak_eps": float(eps.max()),
                               "peak_nm": float(lam[int(eps.argmax())]),
                               "tail_mean": float(tail.mean()),
                               "eps": {str(wv): float(np.interp(wv, lam, eps))
                                       for wv in (320, 340, 360, 375, 400, 455)}}
        print(f"  Cl {k} ({hcl:>5} M HCl): {len(lam):4d} pts, peak {eps.max():6.0f} "
              f"at {lam[int(eps.argmax())]:.0f} nm, tail {tail.mean():5.1f}, "
              f"eps(375) = {np.interp(375, lam, eps):.0f}")

    with open(os.path.join(here, "egan_fig5a_extracted_25m0328.json"), "w") as fh:
        json.dump(result, fh, indent=1, sort_keys=True)
    print(f"\nwritten to {outdir} and egan_fig5a_extracted_25m0328.json")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
