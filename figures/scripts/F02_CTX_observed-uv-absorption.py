"""F2. The problem, in an observed spectrum.

Top: how much sunlight Venus reflects (I/F) at each wavelength, as measured by
MESSENGER/MASCS, digitised from Perez-Hoyos et al. (2018), Figure 8; the shaded
envelope is the drawn observation with its plotted error range. Bottom: how
strongly sulfur dioxide absorbs, from the cross section the model uses. Sulfur
dioxide absorbs strongly only below about 320 nm, yet the planet stays dark out
to about 400 nm. Something else absorbs in between: the unknown absorber.

Sources: figures/data/perez-hoyos-2018-fig8-IF.csv (provenance in the matching
.PROVENANCE.txt file); thermo/photo_cross/SO2/SO2_cross.csv.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

d = np.genfromtxt(os.path.join(S.ROOT, "figures", "data", "perez-hoyos-2018-fig8-IF.csv"),
                  delimiter=",", names=True)
x = np.loadtxt(os.path.join(S.ROOT, "thermo", "photo_cross", "SO2", "SO2_cross.csv"),
               delimiter=",", comments="#")
SO2_COL = "#882255"
lo, hi = 320, 400

fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 8.2), sharex=True,
                             gridspec_kw=dict(height_ratios=[1.3, 1], hspace=0.08))
for a in (a1, a2):
    a.axvspan(lo, hi, color=S.BAND_COLOR, lw=0, zorder=0)
a1.fill_between(d["lambda_nm"], d["IF_low"], d["IF_high"], color=S.COLOR["neutral"],
                alpha=0.35, lw=0)
a1.plot(d["lambda_nm"], d["IF_mid"], color="#333333", lw=2)
a1.text((lo + hi) / 2, 0.97, "unknown\nabsorber", transform=a1.get_xaxis_transform(),
        ha="center", va="top", fontsize=14, fontweight="bold")
a1.text(590, 0.42, "Venus reflects well\nin visible light", ha="center", fontsize=13)
a1.text(d["lambda_nm"][0] - 4, 0.12, "no data\nin this\nspectrum", ha="right",
        fontsize=12, color="#555555")
a1.set_ylim(0, 0.75)
a1.set_ylabel("Reflectivity, I/F")
a1.set_title("Venus is dark between " + S.lit("320 and 400 nm", "the observed band")
             + ", where sulfur dioxide absorbs only weakly", loc="left", fontsize=15)

m = (x[:, 0] >= 200) & (x[:, 0] <= 700) & (x[:, 1] > 0)
a2.plot(x[m, 0], x[m, 1], color=SO2_COL, lw=2.4)
a2.set_yscale("log")
a2.set_ylim(1e-25, 1e-16)
a2.text(215, 3e-18, "sulfur dioxide", color=SO2_COL, fontsize=14, va="bottom")
a2.set_ylabel("Sulfur dioxide\ncross section, cm$^{2}$\n(log scale)")
a2.set_xlabel("Wavelength, nm")
a2.set_xlim(200, 700)
fig.text(0.125, -0.02, "Observed: MESSENGER spectrum, digitised from "
         + S.lit("Pérez-Hoyos et al. (2018), Figure 8", "citation") + ".  Grey band: "
         + S.lit("320 to 400 nm", "the observed band") + ".", fontsize=11.5, color="#444444")
S.save(fig, "F02_CTX_observed-uv-absorption_ch1_s2")
