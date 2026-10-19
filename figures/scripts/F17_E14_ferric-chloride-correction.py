"""F17. Experiment 14: the ferric chloride correction.

Left: the five measured ferric chloride curves in concentrated sulfuric acid,
as mass absorption coefficients, against the single aqueous value the published
rejection used. Right: the droplet mass fraction each would need.

Sources: Experiment14/data/fecl3_Cl*_eps.dat (Egan et al. 2025, Figure 5a);
Experiment14/fecl3_result_25m0328.json; registry e14.*.
"""
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

CURVES = ["ClA", "ClB", "ClC", "ClD", "ClE"]
# five shades of one hue, darkest for the most chloride-rich curve
SHADES = ["#08306B", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1"]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 6.2),
                             gridspec_kw=dict(width_ratios=[1.15, 1]))

# ---- left: spectra
lam375 = S.reg("e14.lambda_nm")
for c, sh in zip(CURVES, SHADES):
    d = np.loadtxt(os.path.join(S.ROOT, "Experiment14", "data", f"fecl3_{c}_eps.dat"))
    # molar to mass coefficient with the molar mass implied by the result file
    molar_mass = 1000.0 * S.reg(f"e14.{c}.eps_molar") / S.reg(f"e14.{c}.eps_mass")
    eps_mass = 1000.0 * d[:, 1] / molar_mass
    a1.plot(d[:, 0], eps_mass, color=sh, lw=2.4)
    # direct label beside the curve's own value at 375 nm, in the empty space
    # to the right where every curve has already dropped
    at375 = S.reg(f"e14.{c}.eps_mass")
    a1.plot([lam375], [at375], "o", color=sh, ms=7)
    a1.text(468, at375, "HCl " + S.val(f"e14.{c}.hcl_M", ".3g") + " M",
            color=sh, va="center", fontsize=12.5)
aq = S.reg("e14.aqueous_eps_mass")
a1.axhline(aq, color=S.COLOR["control"], lw=2.2, ls="--")
# below the dashed line on the left: every curve there is far above it, so
# the label needs no background box and covers nothing
a1.text(203, aq * 0.86, "aqueous value used in\nthe published rejection,\n"
        + S.val("e14.aqueous_eps_mass", ".0f") + " cm$^{2}$ g$^{-1}$",
        va="top", fontsize=10.5)
a1.axvline(lam375, color=S.COLOR["neutral"], lw=1.4, ls=":")
a1.text(lam375 + 3, 0.97, S.val("e14.lambda_nm", ".0f") + " nm",
        transform=a1.get_xaxis_transform(), va="top", fontsize=13)
a1.set_yscale("log")
a1.set_xlim(200, 540)
a1.set_xlabel("Wavelength, nm")
a1.set_ylabel("Mass absorption coefficient,\ncm$^{2}$ g$^{-1}$ (log scale)")
a1.set_title("Measured in concentrated sulfuric acid", loc="left", fontsize=15)

# ---- right: required share of the droplet
labels = ["aqueous"] + ["HCl " + S.val(f"e14.{c}.hcl_M", ".3g") + " M" for c in CURVES]
vals = [S.reg("e14.aqueous_pct_of_droplet")] + [S.reg(f"e14.{c}.pct_of_droplet") for c in CURVES]
texts = [S.val("e14.aqueous_pct_of_droplet", ".1f")] + \
        [S.val(f"e14.{c}.pct_of_droplet", ".1f") for c in CURVES]
colors = [S.COLOR["control"]] + SHADES
y = np.arange(len(vals))[::-1]
a2.barh(y, vals, color=colors, height=0.62)
for yy, v, t in zip(y, vals, texts):
    if v > 50:   # long bar: label inside, clear of the 100 per cent line
        a2.text(v - 1.5, yy, t + " per cent", va="center", ha="right", fontsize=13,
                color="white")
    else:
        a2.text(v + 1.5, yy, t + " per cent", va="center", fontsize=13)
a2.axvline(100, color="#222222", lw=2.6)
a2.text(98, y[0] + 0.55, S.lit("100 per cent: the whole droplet", "definition of the ceiling"),
        ha="right", va="bottom", fontsize=12.5)
a2.set_yticks(y)
a2.set_yticklabels(labels)
a2.set_xlim(0, 108)
a2.set_ylim(-0.7, len(vals) - 0.1)
a2.set_xlabel("Ferric chloride needed, per cent of droplet mass")
a2.set_title("Loading needed to meet the requirement", loc="left", fontsize=15)

fig.suptitle("The published rejection used an absorbance measured in the wrong solvent",
             fontsize=16, y=1.01)
fig.text(0.5, -0.03, "Most chloride-rich curve: a factor of "
         + S.val("e14.factor_to_fit_after", ".1f")
         + S.lit(" from the independent Egan et al. (2026) fit", "year of a citation")
         + ", against " + S.val("e14.factor_to_fit_before", ".1f") + " before the correction",
         ha="center", fontsize=13)
fig.tight_layout()
S.save(fig, "F17_E14_ferric-chloride-correction_ch6_s16")
