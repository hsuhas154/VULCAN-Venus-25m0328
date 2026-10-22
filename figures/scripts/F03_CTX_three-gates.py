"""F3. The three tests of Experiment 12, as a funnel.

Every species the screen could test (23) passes through three tests in turn:
right altitude, right spectrum, enough material. The pass marks are stated on
the figure and defined once. Each is a mark the empirical absorber, the one
known to reproduce the observed albedo, itself meets. They are defined, in verify/recompute_25m0328.py (block "Experiment 12
gates"). The tests are drawn in the order that keeps the most species longest,
so that the last test visibly removes the last two.

Sources: registry e12.gate.n_in, e12.gate.n_placement, e12.gate.n_place_and_spec,
e12.gate.n_all, e12.gate.inband_min_pct (the empirical absorber's own share),, e12.gate.share_min_pct,
e12.gate.corr_min, e6.interval.90.lo, e12.S2O2.short.360, e12.HOCl.short.360.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

STAGES = [("e12.gate.n_in", "species screened", ""),
          ("e12.gate.n_placement", "right altitude",
           "at least " + S.val("e12.gate.inband_min_pct", ".2f")
           + " per cent of the\ncolumn above " + S.lit("58 km", "lower edge of the Experiment 12 screen")
           + " inside\nthe response band"),
          ("e12.gate.n_place_and_spec", "and right spectrum",
           "more than " + S.val("e12.gate.share_min_pct", ".0f")
           + " per cent\nof the absorption in\n" + S.lit("320 to 400 nm", "the observed band")
           + ", correlation\nabove " + S.val("e12.gate.corr_min", ".1f")),
          ("e12.gate.n_all", "and enough material",
           "column above " + S.val("e6.interval.90.lo", ".0f") + " km\nsupplies the required\noptical depth at " + S.lit("360 nm", "the test wavelength, as F15")),
          ]

fig, ax = plt.subplots(figsize=(13, 6.2))
n = [S.reg(k) for k, _, _ in STAGES]
N0 = n[0]
x0, w, gap = 0.0, 1.0, 0.55
for i, (k, name, rule) in enumerate(STAGES):
    h = max(n[i] / N0, 0.04)
    x = x0 + i * (w + gap)
    col = S.COLOR["nominal"] if n[i] > 0 else S.COLOR["neutral"]
    if n[i] > 0:
        ax.add_patch(plt.Rectangle((x, -h / 2), w, h, color=col, lw=0))
    else:  # nothing left: an empty outline, not a bar
        ax.add_patch(plt.Rectangle((x, -0.06), w, 0.12, fill=False, ec=col, lw=2, ls="--"))
    if i + 1 < len(STAGES):
        h2 = max(n[i + 1] / N0, 0.04)
        ax.add_patch(Polygon([(x + w, h / 2), (x + w + gap, h2 / 2),
                              (x + w + gap, -h2 / 2), (x + w, -h / 2)],
                             closed=True, color=S.BAND_COLOR, lw=0))
    ax.text(x + w / 2, max(h / 2, 0.08) + 0.04, S.val(k, ".0f"), ha="center", va="bottom",
            fontsize=22, fontweight="bold", color=col)
    ax.text(x + w / 2, -0.62, name, ha="center", va="top", fontsize=14)
    if rule:
        ax.text(x + w / 2, -0.72, rule, ha="center", va="top", fontsize=11.5, color="#444444")

# who the two near-misses are, and how far short they fall at the last test
x = x0 + 2 * (w + gap)
ax.text(x + w / 2, 0.36, S.chem("S2O2") + " (OSSO) and HOCl", ha="center",
        va="bottom", fontsize=12)
x = x0 + 3 * (w + gap)
ax.text(x + w / 2, 0.36, "short by factors of\n" + S.sci("e12.S2O2.short.360", 1)
        + " and " + S.sci("e12.HOCl.short.360", 1), ha="center", va="bottom", fontsize=12)
fig.text(0.02, 0.015, "Every pass mark is one the empirical absorber, which reproduces the observed "
         "brightness of Venus, itself meets.", fontsize=11, color="#444444")
ax.set_xlim(-0.2, x0 + 4 * w + 3 * gap + 0.2)
ax.set_ylim(-1.05, 0.75)
ax.axis("off")
ax.set_title("No species passes all three tests", loc="left")
S.save(fig, "F03_CTX_three-gates_ch1_backup")
