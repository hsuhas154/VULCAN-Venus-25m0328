"""F19. The transport timescale problem.

How long vertical mixing takes to carry a parcel across a layer, against how
long the atmosphere takes to go once round the planet and how long one Venus
solar day lasts. Crossing time is thickness squared over the geometric mean
eddy diffusion coefficient of the layer, from the model's own profiles; the
windows are defined in verify/recompute_25m0328.py (block "Phase 2").

Sources: registry p2.raised.Bkzz.t_days, p2.raised.BDkzz.t_minutes,
p2.gap.Bkzz.t_days, p2.band.Bkzz.t_days, p2.kzz.raised.z_lo, p2.kzz.raised.z_hi;
literature (figures/literature_values.json) zonal_circuit_period_days_limaye
(Limaye et al. 2018), solar_day_days (Bertaux et al. 2016).
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

L = json.load(open(os.path.join(S.ROOT, S.LIT_PATH)))
zr = S.val("p2.kzz.raised.z_lo", ".0f") + " to " + S.val("p2.kzz.raised.z_hi", ".0f") + " km"
c_lo, c_hi = L["zonal_circuit_period_days_limaye"]["value"]
ROWS = [  # label, value in days (or a range), printed value, colour
    ("vertical mixing across " + zr + ",\nwith the profile Phase 1 used",
     S.reg("p2.raised.Bkzz.t_days"), S.val("p2.raised.Bkzz.t_days", ".0f") + " days", S.COLOR["nominal"]),
    ("the same layer, with the\nraised (alternative) profile",
     S.reg("p2.raised.BDkzz.t_minutes") / 1440.0,
     S.val("p2.raised.BDkzz.t_minutes", ".0f") + " minutes", S.COLOR["model"]),
    ("vertical mixing across the\n" + S.lit("55 to 69 km", "the layer between the raised levels and the band")
     + " gap, either profile",
     S.reg("p2.gap.Bkzz.t_days"), S.val("p2.gap.Bkzz.t_days", ".0f") + " days", S.COLOR["nominal"]),
    ("vertical mixing across the\nresponse band, either profile",
     S.reg("p2.band.Bkzz.t_days"), S.val("p2.band.Bkzz.t_days", ".0f") + " days", S.COLOR["nominal"]),
    ("one trip round the planet\nin the cloud layer",
     (c_lo, c_hi), S.litval("zonal_circuit_period_days_limaye", ".0f", 0) + " to "
     + S.litval("zonal_circuit_period_days_limaye", ".0f", 1) + " days", S.COLOR["neutral"]),
    ("one Venus solar day",
     L["solar_day_days"]["value"], S.litval("solar_day_days", ".0f") + " days", S.COLOR["neutral"]),
]
fig, ax = plt.subplots(figsize=(12.5, 6.8))
y = np.arange(len(ROWS))[::-1]
for yi, (lab, v, txt, col) in zip(y, ROWS):
    if isinstance(v, tuple):
        ax.barh(yi, v[1], color=col, height=0.6, alpha=0.45)
        ax.barh(yi, v[0], color=col, height=0.6)
        end = v[1]
    else:
        ax.barh(yi, v, color=col, height=0.6)
        end = v
    ax.text(end * 1.25, yi, txt, va="center", fontsize=13)
# row labels as text, not tick labels, so the guard reads the numbers in them
ax.set_yticks([])
for yi, r in zip(y, ROWS):
    ax.text(-0.015, yi, r[0], transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=12.5)
ax.set_xscale("log")
ax.set_xlim(3e-3, 3e4)
ax.set_xlabel("Time, Earth days (log scale)")
note = ("While mixing carries a parcel across the gap, the winds take it round the planet\n"
        + S.dval("a/b", ["p2.gap.Bkzz.t_days", "lit:zonal_circuit_period_days_limaye:1"], ".0f")
        + " to " + S.dval("a/b", ["p2.gap.Bkzz.t_days", "lit:zonal_circuit_period_days_limaye:0"], ".0f")
        + " times, through " + S.dval("a/b", ["p2.gap.Bkzz.t_days", "lit:solar_day_days"], ".1f")
        + " Venus solar days")
fig.text(0.5, -0.04, note, ha="center", va="top", fontsize=13)
ax.set_title("Air goes round Venus far faster than it rises", loc="left")
S.save(fig, "F19_P2_transport-timescales_ch9_s19")
