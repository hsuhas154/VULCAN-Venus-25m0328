"""F11. Experiment 8: spectral shape matters less than placement.

Three candidate spectra (ClS2, OSSO, S3) put in place of the empirical one, at
the nominal placement and column. The spread between best and worst is smaller
than the effect of moving the absorber (Experiment 7).

Sources: atm/UV_absorber.txt, atm/UV_absorber_Experiment8A_ClS2_25m0328.txt,
atm/UV_absorber_Experiment8B_OSSO_S2O2_25m0328.txt,
atm/UV_absorber_Experiment8C_S3_25m0328.txt; registry e8a.L2pct, e8b.L2pct,
e8c.L2pct, e7a.L2pct, e7b.L2pct (spread and ratio derived from these).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

def spec(path):
    d = np.loadtxt(os.path.join(S.ROOT, path), skiprows=1)
    return d[:, 0], d[:, 1]

CAND = [("e8a", "ClS2", "atm/UV_absorber_Experiment8A_ClS2_25m0328.txt", "#332288"),
        ("e8b", "OSSO (S2O2)", "atm/UV_absorber_Experiment8B_OSSO_S2O2_25m0328.txt", "#117733"),
        ("e8c", "S3",   "atm/UV_absorber_Experiment8C_S3_25m0328.txt", "#AA4499")]

fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(16, 6),
                                 gridspec_kw=dict(width_ratios=[1.5, 1, 0.8]))
lo, hi = 320, 400
a1.axvspan(lo, hi, color=S.BAND_COLOR, lw=0, zorder=0,
           label=S.lit("observed band, 320 to 400 nm",
                       "the observed absorption band, a definition"))
w, q = spec("atm/UV_absorber.txt")
a1.plot(w, q, color=S.COLOR["empirical"], lw=3.2, label="empirical")
for key, name, path, col in CAND:
    w, q = spec(path)
    a1.plot(w, q, color=col, lw=2.0, label=S.chem(name))
# a legend in the empty upper right, so no label sits on a curve
a1.legend(loc="upper right", fontsize=12, frameon=False, handlelength=1.6)
a1.set_xlim(280, 700)
a1.set_ylim(0, 0.6)
a1.set_xlabel("Wavelength, nm")
a1.set_ylabel("Absorption efficiency, as used by the model")
a1.set_title("Candidate spectra", loc="left")

x = np.arange(3)
vals = [S.reg(k + ".L2pct") for k, *_ in CAND]
a2.bar(x, vals, color=[c for *_, c in CAND], width=0.6)
for xi, (k, *_), v in zip(x, CAND, vals):
    a2.text(xi, v * 1.02, S.val(k + ".L2pct", ".4f"), ha="center", va="bottom", fontsize=12.5)
a2.set_xticks(x)
a2.set_xticklabels([S.chem(n.split(" ")[0]) for _, n, *_ in CAND])
a2.set_ylim(0, max(vals) * 1.25)
a2.set_ylabel("Change in the radiation field,\nper cent of the nominal field")
a2.set_title("Same placement, three spectra", loc="left")

spread = S.dval("max/min", ["e8a.L2pct", "e8b.L2pct", "e8c.L2pct"], ".2f")
ratio = S.dval("a/b", ["e7a.L2pct", "e7b.L2pct"], ".2f")
r_num = S.DERIVED["a/b"]([S.reg("e7a.L2pct"), S.reg("e7b.L2pct")])
s_num = S.DERIVED["max/min"](vals)
a3.bar([0, 1], [r_num, s_num], color=["#444444", "#AAAAAA"], width=0.6)
a3.text(0, r_num * 1.02, ratio, ha="center", va="bottom", fontsize=16, fontweight="bold")
a3.text(1, s_num * 1.02, spread, ha="center", va="bottom", fontsize=16, fontweight="bold")
a3.set_xticks([0, 1])
a3.set_xticklabels(["moving it\n(Exp 7)", "changing\nspectrum"])
a3.set_ylim(0, r_num * 1.25)
a3.set_ylabel("Largest over smallest effect")
a3.set_title("Placement wins", loc="left")
fig.tight_layout()
S.save(fig, "F11_E08_spectrum-vs-placement_ch6_s12")
