#!/usr/bin/env python3
"""
Extract the rhomboclase and acid ferric sulfate extinction coefficients from the
supplementary figures of Jiang et al. (2024).
Project roll number: 25m0328

WHY THIS EXISTS
The numerical data behind Figures S9 and S10 of that paper were never published:
the data availability statement says everything is in the paper or supplement,
and in the supplement those coefficients appear only as plots. The plots are
vector graphics, not bitmaps, so the curve is stored in the PDF as an exact
ordered list of coordinates. This script reads that coordinate list. It does not
read pixels off a picture and it does not interpolate between eyeballed points.

The output files data/*.dat are committed, so the rest of the analysis runs
without this script. It exists so the numbers are reproducible and so the
validation below is part of the record rather than a claim in a document.

THREE VALIDATION CHECKS, ALL FATAL
1. Frame identification. The plot box is taken as the largest inner rectangle.
2. Wavelength scale. The extracted curve's endpoints must land on the stated
   axis limits. They are not fitted to them, so recovering 200 and 800 nm is an
   independent check that the horizontal calibration is right.
3. Extinction scale. The minor tick marks on the vertical axis give the scale a
   second time, independently of the frame. Five minor intervals make one major
   interval of 100 cm2/g, so the frame height divided by that must come out as
   the whole number of major divisions the axis shows. An earlier version of
   this extraction picked an inner clipping rectangle instead of the axis box
   and inflated the acid ferric sulfate peak by 6.4 per cent. Under this check
   that frame gives 7.52 divisions instead of 8 and the script stops.

REQUIREMENTS
pdftocairo, from poppler-utils. If it is missing, install it with
    sudo apt install poppler-utils
or simply use the committed data/*.dat files, which this script produced.

USAGE
    python Experiment13/extract_jiang_spectra_25m0328.py /path/to/jiang_sm.pdf
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

# page in the supplementary PDF, axis maximum in cm2/g, output name
FIGURES = {
    14: dict(ymax=800.0, name="acid-ferric-sulfate",
             desc="acid ferric sulfate (H3O)Fe(SO4)2, Figure S9"),
    15: dict(ymax=500.0, name="rhomboclase",
             desc="rhomboclase (H5O2)Fe(SO4)2.3H2O, Figure S10"),
}
LAM_LO, LAM_HI = 200.0, 800.0      # stated wavelength axis limits
CURVE_STROKE = "rgb(29.803467%, 44.7052%, 69.020081%)"   # the smoothed fit
MINORS_PER_MAJOR = 5
MAJOR_UNITS = 100.0
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def parse_path(d):
    """Return the anchor points of an SVG path. Curve segments contribute their
    endpoint, which is the vertex the plotting program emitted."""
    toks = d.replace(",", " ").split()
    pts, i, cmd = [], 0, None
    while i < len(toks):
        t = toks[i]
        if t in "MLCZmlcz":
            cmd, i = t, i + 1
            continue
        if cmd in "ML":
            pts.append((float(toks[i]), float(toks[i + 1]))); i += 2
        elif cmd == "C":
            pts.append((float(toks[i + 4]), float(toks[i + 5]))); i += 6
        else:
            i += 1
    return np.array(pts)


def transform(tr):
    m = re.search(r"matrix\(([^)]+)\)", tr or "")
    if not m:
        return np.eye(3)
    a, b, c, d, e, f = [float(x) for x in m.group(1).replace(",", " ").split()]
    return np.array([[a, c, e], [b, d, f], [0, 0, 1]])


def elements(svg):
    """Yield (element_text, points_in_page_coordinates) for every path."""
    for el in re.findall(r"<path[^>]*/>", svg):
        d = re.search(r'\sd="([^"]+)"', el)
        if not d:
            continue
        p = parse_path(d.group(1))
        if len(p) < 2:
            continue
        tr = re.search(r'transform="([^"]+)"', el)
        yield el, (np.c_[p, np.ones(len(p))] @ transform(tr.group(1) if tr else None).T)[:, :2]


def extract(svg, ymax, desc):
    segments, rects, ticks = [], [], []
    for el, q in elements(svg):
        w = float(q[:, 0].max() - q[:, 0].min())
        h = float(q[:, 1].max() - q[:, 1].min())
        if CURVE_STROKE in el:
            segments.append(q)
        elif len(q) <= 6 and 250 < w < 420 and h > 150:
            rects.append((h, w, q))
        elif len(q) == 2 and abs(q[1, 1] - q[0, 1]) < 0.5 and 1 < w < 15:
            ticks.append((float(q[:, 0].min()), float(q[:, 1].mean())))

    if not segments:
        raise SystemExit(f"{desc}: no curve found with the expected stroke colour")
    if not rects:
        raise SystemExit(f"{desc}: no axis rectangle found")

    # check 1: the axis box is the largest inner rectangle
    rects.sort(key=lambda r: -r[0])
    h, w, fr = rects[0]
    x0, x1 = float(fr[:, 0].min()), float(fr[:, 0].max())
    y0, y1 = float(fr[:, 1].min()), float(fr[:, 1].max())

    # check 3: the vertical scale, from tick spacing, independent of the frame
    left = [t for t in ticks if abs(t[0] - x0) < 3.0]
    ys = sorted({round(t[1], 3) for t in left})
    if len(ys) < 6:
        raise SystemExit(f"{desc}: found only {len(ys)} axis ticks, cannot verify the scale")
    gaps = np.diff(ys)
    minor = float(np.median(gaps[gaps < 1.5 * np.median(gaps)]))
    divisions = h / (MINORS_PER_MAJOR * minor)
    expected = ymax / MAJOR_UNITS
    if abs(divisions - expected) > 0.02 * expected:
        raise SystemExit(
            f"{desc}: vertical scale check FAILED. Tick spacing implies "
            f"{divisions:.2f} major divisions across the frame, the axis shows "
            f"{expected:.0f}. The wrong rectangle was probably selected.")

    a = np.vstack(segments)
    a = a[np.argsort(a[:, 0])]
    lam = LAM_LO + (a[:, 0] - x0) / (x1 - x0) * (LAM_HI - LAM_LO)
    eps = (y1 - a[:, 1]) / (y1 - y0) * ymax

    # check 2: the horizontal scale, from endpoints that were never fitted
    for got, want in ((lam.min(), LAM_LO), (lam.max(), LAM_HI)):
        if abs(got - want) > 1.0:
            raise SystemExit(f"{desc}: wavelength scale check FAILED, endpoint "
                             f"{got:.2f} nm against a stated axis limit of {want:.0f} nm")
    if eps.min() < -1.0:
        raise SystemExit(f"{desc}: extracted a negative extinction coefficient, {eps.min():.2f}")

    print(f"  {desc}")
    print(f"    frame {w:.2f} x {h:.2f} pt; tick spacing implies {divisions:.3f} major "
          f"divisions against {expected:.0f} shown")
    print(f"    endpoints recovered at {lam.min():.2f} and {lam.max():.2f} nm "
          f"against stated {LAM_LO:.0f} and {LAM_HI:.0f}")
    print(f"    {len(lam)} points; peak {eps.max():.1f} cm2/g at "
          f"{lam[int(np.argmax(eps))]:.1f} nm")
    return lam, eps


def main():
    if shutil.which("pdftocairo") is None:
        sys.exit("pdftocairo not found. Install poppler-utils, or use the committed "
                 "data/*.dat files that this script produced.")
    if len(sys.argv) < 2:
        sys.exit(__doc__.strip().splitlines()[-1])
    pdf = sys.argv[1]
    if not os.path.exists(pdf):
        sys.exit(f"not found: {pdf}")

    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"source: {pdf}\n")
    with tempfile.TemporaryDirectory() as tmp:
        for page, spec in sorted(FIGURES.items()):
            out = os.path.join(tmp, f"p{page}")
            subprocess.run(["pdftocairo", "-svg", "-f", str(page), "-l", str(page),
                            pdf, out + ".svg"], check=True)
            lam, eps = extract(open(out + ".svg").read(), spec["ymax"], spec["desc"])
            path = os.path.join(OUT_DIR, f"jiang_{spec['name']}_extracted.dat")
            np.savetxt(path, np.c_[lam, eps], fmt="%10.4f %12.5f",
                       header="wavelength_nm  extinction_cm2_per_g\n"
                              f"{spec['desc']}, Jiang et al. 2024 supplementary material\n"
                              "extracted from the vector paths of the published figure")
            print(f"    wrote {path}\n")
    print("all validation checks passed")


if __name__ == "__main__":
    main()
