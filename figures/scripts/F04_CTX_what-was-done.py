"""F4. What was actually done: the fourteen experiments in four stages.

Switch the absorber on and check it, find where its effect lives (amount,
wavelength, reaction, altitude), find what decides the size of the effect, then
test real candidates against what was learned. Box wording follows the
experiments as described in report Chapters 5 and 6. No numbers.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

STAGES = [
    ("Switch it on", S.COLOR["control"], [
        ("Exp 1", "does the absorber\ndo anything?"),
        ("Exp 2", "what does it\nchange?")]),
    ("Find where it acts", S.COLOR["nominal"], [
        ("Exp 3", "is the response\nproportional?"),
        ("Exp 4", "which wavelengths?"),
        ("Exp 5", "which reactions?"),
        ("Exp 6", "which altitudes?")]),
    ("Find what decides it", S.COLOR["model"], [
        ("Exp 7", "moving the absorber\nup or down"),
        ("Exp 8", "changing its\nspectrum"),
        ("Exp 9", "letting chemistry\nplace it"),
        ("Exp 10", "is the result\nreproducible?")]),
    ("Test real candidates", S.COLOR["accent"], [
        ("Exp 11", "computed cross\nsections"),
        ("Exp 12", "screen of every\nspecies"),
        ("Exp 13", "iron sulfate\nminerals"),
        ("Exp 14", "ferric chloride\nin acid")]),
]
fig, ax = plt.subplots(figsize=(15, 7.2))
W, H, GX, GY = 2.6, 1.15, 1.5, 0.35
heads = []
for i, (title, col, boxes) in enumerate(STAGES):
    x = i * (W + GX)
    heads.append(ax.text(x + W / 2, 0.55, title, ha="center", va="bottom", fontsize=14.5,
                         fontweight="bold", color=col))
    for j, (exp, q) in enumerate(boxes):
        y = -j * (H + GY) - H
        ax.add_patch(FancyBboxPatch((x, y), W, H, boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc="white", ec=col, lw=2.2))
        ax.text(x + 0.15, y + H - 0.14, exp, ha="left", va="top", fontsize=12.5,
                fontweight="bold", color=col)
        ax.text(x + W / 2, y + H * 0.38, q, ha="center", va="center", fontsize=12.5)
ax.set_xlim(-0.2, 4 * W + 3 * GX + 0.2)
ax.set_ylim(-4 * (H + GY) - 0.1, 1.2)
ax.axis("off")
# arrows join the stage headings, gap to gap, measured from the drawn text
fig.canvas.draw()
inv = ax.transData.inverted()
ext = [inv.transform(h.get_window_extent()) for h in heads]
for e0, e1 in zip(ext, ext[1:]):
    ym = (e0[0][1] + e0[1][1]) / 2
    ax.annotate("", xy=(e1[0][0] - 0.08, ym), xytext=(e0[1][0] + 0.08, ym),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=2, mutation_scale=20))
ax.set_title("What was done, in four stages", loc="left")
S.save(fig, "F04_CTX_what-was-done_ch3_report")
