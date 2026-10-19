"""F18. The two eddy diffusion profiles that ship with the model.

atm/atm_Venus_Bkzz.txt is the profile every Phase 1 run used. atm/atm_Venus_B+Dkzz.txt
is the same file with four levels raised to the Dai et al. (2023) values. Both are
drawn against the model's own altitude grid. The shaded range is what three
published studies derive for the convective cloud layer; the profile in use sits
far below it.

Sources: atm/atm_Venus_Bkzz.txt, atm/atm_Venus_B+Dkzz.txt, model altitudes from
output/Experiment2_A1_*.vul; registry p2.kzz.raised.ratio_min, p2.kzz.raised.ratio_max,
p2.kzz.raised.z_lo, p2.kzz.raised.z_hi; literature (figures/literature_values.json)
kzz_lef24_48_55km.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

z = np.asarray(S.vul("Experiment2_A1_UUV-nominal-25m0328.vul")["atm"]["zmco"]) / 1e5
KB = np.loadtxt(os.path.join(S.ROOT, "atm", "atm_Venus_Bkzz.txt"), skiprows=2)[:, 2]
KD = np.loadtxt(os.path.join(S.ROOT, "atm", "atm_Venus_B+Dkzz.txt"), skiprows=2)[:, 2]
import json
L = json.load(open(os.path.join(S.ROOT, S.LIT_PATH)))

fig, ax = plt.subplots(figsize=(10, 7.5))
S.response_band(ax, where="above")
lo, hi = L["kzz_lef24_48_55km"]["value"]
zl, zh = S.reg("p2.kzz.raised.z_lo"), S.reg("p2.kzz.raised.z_hi")
# drawn over the model layers that the raised levels occupy (layer centres +/- 1 km)
ax.fill_betweenx([zl - 1, zh + 1], lo, hi, color=S.COLOR["accent"], alpha=0.18, lw=0)
ax.text(np.sqrt(lo * hi), zh + 3, "convective cloud layer, "
        + S.lit("Lefèvre et al. (2024)", "citation year") + ":\n"
        + S.sci_lit("kzz_lef24_48_55km", 0, 0) + " to " + S.sci_lit("kzz_lef24_48_55km", 0, 1),
        ha="center", va="bottom", fontsize=12, color="#1b6b52")
ax.plot(KD, z, color=S.COLOR["model"], lw=2.4, ls="--")
ax.plot(KB, z, color=S.COLOR["nominal"], lw=3)
ax.text(KB[np.argmin(np.abs(z - 30))] * 1.6, 30, "used in every Phase 1 run\n(Bierson and Zhang profile)",
        color=S.COLOR["nominal"], fontsize=13, va="center")
ax.text(KD.max() * 0.6, zl - 4, "alternative file:\nraised at four levels,\nby " +
        S.val("p2.kzz.raised.ratio_min", ".0f") + " to " + S.val("p2.kzz.raised.ratio_max", ".0f")
        + " times", color=S.COLOR["model"], fontsize=12.5, va="top", ha="right")
ax.set_xscale("log")
ax.set_xlim(3e3, 3e9)
ax.set_ylim(0, 112)
ax.set_xlabel("Eddy diffusion coefficient, cm$^{2}$ per second (log scale)")
ax.set_ylabel("Altitude, km")
ax.set_title("The mixing used in Phase 1 is far below what is measured", loc="left")
S.save(fig, "F18_P2_two-mixing-profiles_ch9_report")
