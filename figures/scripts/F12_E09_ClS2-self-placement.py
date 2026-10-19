"""F12. Experiment 9: the decisive result.

The empirical absorber profile against the profile the model's own chemistry
gives ClS2. Same molecule, same spectrum: placed by the chemistry, it settles
near the surface and its optical depth above the response band collapses.

Sources: atm/mode1+2.txt, atm/mode1+2_Experiment9_ClS2profile_25m0328.txt;
registry pnom.centroid, p9.centroid, spec_emp.tau69, spec_9.tau69,
tau.reduction, pnom.inband, p9.inband.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

def profile(path):
    d = np.loadtxt(os.path.join(S.ROOT, path), skiprows=1)
    return d[:, 0], d[:, 1]

z_e, n_e = profile("atm/mode1+2.txt")
z_m, n_m = profile("atm/mode1+2_Experiment9_ClS2profile_25m0328.txt")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 6.2), sharey=True, sharex=True)
panels = (
    (a1, z_e, n_e, "empirical", "Empirical profile, as prescribed in the model",
     "pnom.centroid", "spec_emp.tau69", None, "pnom.inband"),
    (a2, z_m, n_m, "model", S.chem("ClS2") + " placed by the model's own chemistry",
     "p9.centroid", "spec_9.tau69", "sci", "p9.inband"),
)
for ax, zz, nn, col, title, kc, kt, tfmt, kb in panels:
    S.response_band(ax, label=(ax is a2))
    ax.plot(nn, zz, color=S.COLOR[col], lw=2.8)
    ax.fill_betweenx(zz, 0, nn, color=S.COLOR[col], alpha=0.25, lw=0)
    c = S.reg(kc)
    ax.axhline(c, color=S.COLOR[col], ls="--", lw=1.6)
    ax.text(0.97, c + 2.0, "centroid " + S.val(kc, ".2f") + " km",
            transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            color="#222222", fontsize=13)
    tau = S.sci(kt, 4) if tfmt else S.val(kt, ".4f")
    ax.text(0.97, 0.985, S.lit("optical depth above 69 km: ", "lower edge of the response band, a definition") + tau + "\n"
            + "share of profile in the band: " + S.val(kb, ".4f") + " per cent",
            transform=ax.transAxes, ha="right", va="top", fontsize=12.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#BBBBBB"))
    ax.set_title(title, loc="left", fontsize=15)
    ax.set_xlabel("Number density, particles per cm$^{3}$")
    ax.set_xlim(left=0)
    ax.set_ylim(0, 112)
a1.set_ylabel("Altitude, km")
fig.suptitle("Same molecule, same spectrum: placed by its own chemistry, it falls "
             "out of the response band\n(optical depth above the band falls by a factor of "
             + S.sci("tau.reduction", 1) + ")", fontsize=15, y=1.02)
fig.tight_layout()
S.save(fig, "F12_E09_ClS2-self-placement_ch6_s13")
