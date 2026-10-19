"""F7. Experiment 4: the response is spectrally local.

The absorber restricted to one wavelength band at a time (six runs). In each
run the radiation field changes only inside the band the absorber was given.

Sources: output/Experiment1_A0_*.vul, output/Experiment4B_B*.vul,
atm/UV_absorber_Experiment4B_B*_25m0328.txt; registry e4b.B1 to e4b.B6
support_lo and support_hi.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

A0 = S.vul("Experiment1_A0_UUV-control-25m0328.vul")
z, bins, quad = S.grid(A0)
F0 = np.asarray(A0["variable"]["aflux"])
BANDS = [("B1", "300-320"), ("B2", "320-340"), ("B3", "340-360"),
         ("B4", "360-380"), ("B5", "380-420"), ("B6", "420-698")]
sel = (bins >= 280) & (bins <= 700)

fig, axes = plt.subplots(len(BANDS), 1, figsize=(11, 8.5), sharex=True)
for ax, (b, r) in zip(axes, BANDS):
    d = np.abs(np.asarray(S.vul(f"Experiment4B_{b}_{r}-25m0328.vul")["variable"]["aflux"]) - F0)
    spec = d.sum(axis=0)            # summed over altitude, at each wavelength
    spec = spec / spec[sel].max()
    lo, hi = S.reg(f"e4b.{b}.support_lo"), S.reg(f"e4b.{b}.support_hi")
    ax.axvspan(lo, hi, color=S.BAND_COLOR, lw=0, zorder=0)
    ax.fill_between(bins[sel], 0, spec[sel], color=S.COLOR["nominal"], lw=0, step="mid")
    ax.text(1.01, 0.5, S.val(f"e4b.{b}.support_lo", ".0f") + " to "
            + S.val(f"e4b.{b}.support_hi", ".0f") + " nm", transform=ax.transAxes,
            va="center", fontsize=12.5)
    ax.set_ylim(0, 1.1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
axes[-1].set_xlim(280, 700)
axes[-1].set_xlabel("Wavelength, nm")
axes[0].set_title("Each run changes the radiation field only where its absorber absorbs\n"
                  "(blue: change in actinic flux, scaled to its peak; grey: the band given "
                  "to the absorber)", loc="left", fontsize=14)
S.save(fig, "F07_E04_spectrally-local_ch5_s9")
