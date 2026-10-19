"""F13. Experiment 10: how small an effect the model can resolve.

Two runs whose absorbers are both effectively absent from the response band
(the ClS2 placement of Experiment 9 with two different spectra) should be
identical. Their radiation fields agree to 0.0017 per cent, but individual
species abundances differ by up to 2.53 per cent. So there are two floors: a
tiny one for radiation claims and a much larger one for chemistry claims.
Every radiative result of Phase 1 sits far above the radiation floor.

Sources: output/Experiment9_*.vul, output/Experiment10_*.vul; registry
e10.fielddiff_pct, res.e10.SO.inband, e7a.L2pct, e7b.L2pct, e8a.L2pct,
e8b.L2pct, e8c.L2pct, e9.L2pct.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

RESULTS = [("e7a.L2pct", "Exp 7, moved up"), ("e8c.L2pct", "Exp 8, " + S.chem("S3") + " spectrum"),
           ("e8a.L2pct", "Exp 8, " + S.chem("ClS2") + " spectrum"),
           ("e9.L2pct", "Exp 9, " + S.chem("ClS2") + " placed by chemistry"),
           ("e7b.L2pct", "Exp 7, moved down"), ("e8b.L2pct", "Exp 8, OSSO spectrum")]

fig, ax = plt.subplots(figsize=(11, 6.5))
rad = S.reg("e10.fielddiff_pct"); chem = S.reg("res.e10.SO.inband")
ax.axvspan(1e-5, rad, color="#BBBBBB", lw=0, zorder=0)
ax.axvline(chem, color=S.COLOR["control"], lw=2.2, ls="--")
y = np.arange(len(RESULTS))[::-1]
for yy, (k, name) in zip(y, RESULTS):
    v = S.reg(k)
    ax.plot([v], [yy], "o", color=S.COLOR["nominal"], ms=11)
    ax.text(v * 0.78, yy, name + ": " + S.val(k, ".4f") + " per cent", va="center", ha="right", fontsize=12.5)
ax.text(rad * 0.85, len(RESULTS) + 0.8, "radiation floor\n" + S.val("e10.fielddiff_pct", ".4f")
        + " per cent", ha="right", va="top", fontsize=12.5)
ax.text(chem * 1.15, len(RESULTS) + 0.8, "chemistry floor, " + S.val("res.e10.SO.inband", ".4f")
        + " per cent:\nabundance changes smaller than this are not resolved",
        ha="left", va="top", fontsize=12.5)
ax.set_xscale("log")
ax.set_xlim(1e-5, 1e3)
ax.set_ylim(-0.6, len(RESULTS) + 0.9)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlabel("Size of the change, per cent (log scale)")
ax.set_title("Every radiative result is far above the radiation floor", loc="left")
S.save(fig, "F13_E10_noise-floors_ch6_backup")
