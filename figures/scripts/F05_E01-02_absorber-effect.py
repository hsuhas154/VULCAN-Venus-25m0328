"""F5. Experiments 1 and 2: switching the absorber on has a measurable effect.

Left: optical depth from the top of the atmosphere down to each altitude, at
one wavelength inside the absorber's range, with the absorber off (A = 0) and
on (A = 1). Right: where in altitude and wavelength the radiation field changes.

Sources: output/Experiment1_A0_UUV-control-25m0328.vul,
output/Experiment2_A1_UUV-nominal-25m0328.vul; registry e6.interval.90.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

A0 = S.vul("Experiment1_A0_UUV-control-25m0328.vul")
A1 = S.vul("Experiment2_A1_UUV-nominal-25m0328.vul")
z, bins, quad = S.grid(A0)
zi = np.asarray(A0["atm"]["zco"]) / 1e5
LAM = 350.0
j = int(np.argmin(np.abs(bins - LAM)))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 6.4), gridspec_kw=dict(width_ratios=[1, 1.35]))
S.response_band(a1, label=True, where="bottom")
dt = np.asarray(A1["variable"]["tau"])[:, j] - np.asarray(A0["variable"]["tau"])[:, j]
a1.plot(dt, zi, color=S.COLOR["nominal"], lw=2.8)
a1.fill_betweenx(zi, 0, dt, color=S.COLOR["nominal"], alpha=0.2, lw=0)
a1.set_ylim(0, 112)
a1.set_xlim(left=0)
a1.set_xlabel("Extra optical depth from the top when the\nabsorber is on, at "
              + S.lit("350 nm", "chosen display wavelength"))
a1.set_ylabel("Altitude, km")
a1.set_title("The absorber blocks light high up", loc="left")

F0 = np.asarray(A0["variable"]["aflux"]); F1 = np.asarray(A1["variable"]["aflux"])
d = np.abs(F1 - F0)
sel = (bins >= 280) & (bins <= 700)
show = d[:, sel] / d[:, sel].max()
im = a2.pcolormesh(bins[sel], z, show, shading="nearest", cmap="cividis", vmin=0, vmax=1, rasterized=True, linewidth=0, antialiased=False)
a2.axhline(S.BAND_LO_KM, color="white", lw=1.2, ls="--")
a2.axhline(S.BAND_HI_KM, color="white", lw=1.2, ls="--")
a2.text(690, (S.BAND_LO_KM + S.BAND_HI_KM) / 2, "response band", color="white",
        ha="right", va="center", fontsize=12)
cb = fig.colorbar(im, ax=a2)
cb.set_label("Change in actinic flux, scaled to its largest value")
a2.set_ylim(0, 112)
a2.set_xlabel("Wavelength, nm")
a2.set_title("Where the radiation field changes", loc="left")
fig.tight_layout()
S.save(fig, "F05_E01-02_absorber-effect_ch5_s7")
