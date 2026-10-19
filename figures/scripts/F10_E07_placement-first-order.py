"""F10. Experiment 7: vertical placement is first-order.

The same absorber column and spectrum, moved up by about 18 km or down by
exactly 18 km. Moving it up changes the radiation field several times more than
moving it down.

Sources: atm/mode1+2.txt, atm/mode1+2_Experiment7A_shift18_25m0328.txt,
atm/mode1+2_Experiment7B_shiftDown18_25m0328.txt; registry
e7.7a.centroid_shift, e7.7b.centroid_shift, e7a.L2pct, e7b.L2pct,
e7.ratio_7A_over_7B.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

def profile(path):
    d = np.loadtxt(os.path.join(S.ROOT, path), skiprows=1)
    return d[:, 0], d[:, 1]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 6.4), gridspec_kw=dict(width_ratios=[1.25, 1]))
S.response_band(a1, where="above")
cases = [
    ("atm/mode1+2.txt", S.COLOR["nominal"], "-", "nominal", None),
    ("atm/mode1+2_Experiment7A_shift18_25m0328.txt", S.COLOR["up"], "-", "moved up ",
     "e7.7a.centroid_shift"),
    ("atm/mode1+2_Experiment7B_shiftDown18_25m0328.txt", S.COLOR["down"], "--", "moved down ",
     "e7.7b.centroid_shift"),
]
for path, col, ls, name, key in cases:
    zz, nn = profile(path)
    a1.plot(nn, zz, color=col, lw=2.8, ls=ls)
    k = int(np.argmax(nn))
    label = name if key is None else name + S.val(key, "+.4f") + " km"
    a1.text(nn[k] * 1.03, zz[k], label, color=col, va="center", fontsize=13)
a1.set_ylim(0, 112)
a1.set_xlim(left=0)
a1.set_xlabel("Absorber number density, per cm$^{3}$")
a1.set_ylabel("Altitude, km")
a1.set_title("Same column, same spectrum, moved", loc="left")

vals = [S.reg("e7a.L2pct"), S.reg("e7b.L2pct")]
txt = [S.val("e7a.L2pct", ".4f"), S.val("e7b.L2pct", ".4f")]
x = [0, 1]
a2.bar(x, vals, color=[S.COLOR["up"], S.COLOR["down"]], width=0.6)
for xi, v, t in zip(x, vals, txt):
    a2.text(xi, v * 1.02, t + "\nper cent", ha="center", va="bottom", fontsize=13)
a2.set_xticks(x)
a2.set_xticklabels(["moved up", "moved down"])
a2.set_ylabel("Change in the radiation field,\nper cent of the nominal field")
a2.set_ylim(0, vals[0] * 1.45)
a2.text(0.97, 0.90, "ratio " + S.val("e7.ratio_7A_over_7B", ".2f"), transform=a2.transAxes,
        ha="right", va="top", fontsize=24, fontweight="bold", color="#222222")
a2.set_title("Effect on the radiation field", loc="left")
fig.suptitle("Where the absorber sits matters more than anything else tested", y=1.01)
fig.tight_layout()
S.save(fig, "F10_E07_placement-first-order_ch5_s11")
