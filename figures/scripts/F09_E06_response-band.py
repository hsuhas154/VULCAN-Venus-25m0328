"""F9. Experiment 6: where in altitude the absorber acts.

The absorber's effect on the radiation field, summed over wavelength, at each
altitude. It is concentrated in one band, 69 to 93 km, centred near 82 km, and
the centre does not move when the absorber is restricted to one wavelength band
at a time (the six runs of Experiment 4B, drawn faintly).

Sources: output/Experiment1_A0_*.vul, output/Experiment2_A1_*.vul,
output/Experiment4B_B*.vul; registry e6.A1.meanz, e6.A1.q25, e6.A1.q50,
e6.A1.q75, e6.B1.meanz to e6.B6.meanz, e6.interval.90.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

A0 = S.vul("Experiment1_A0_UUV-control-25m0328.vul")
A1 = S.vul("Experiment2_A1_UUV-nominal-25m0328.vul")
z, bins, quad = S.grid(A0)
F0 = np.asarray(A0["variable"]["aflux"])

def weight(case):
    """Actinic flux change summed over wavelength at each level, as Experiment 6."""
    d = np.abs(np.asarray(case["variable"]["aflux"]) - F0)
    w = (d * quad[None, :]).sum(axis=1)
    return w / w.sum()

BANDS = [("B1", "300-320"), ("B2", "320-340"), ("B3", "340-360"),
         ("B4", "360-380"), ("B5", "380-420"), ("B6", "420-698")]

fig, ax = plt.subplots(figsize=(8.5, 7))
S.response_band(ax, label=False)
# band label just under the shaded band, right side, clear of the curve and of
# the quartile column inside the band
ax.text(0.98, S.BAND_LO_KM - 1.0, "response band " + S.val("e6.interval.90.lo", ".0f")
        + " to " + S.val("e6.interval.90.hi", ".0f") + " km",
        transform=ax.get_yaxis_transform(), ha="right", va="top", fontsize=12,
        color="#555555")
for b, r in BANDS:
    w = weight(S.vul(f"Experiment4B_{b}_{r}-25m0328.vul"))
    ax.plot(100 * w, z, color=S.COLOR["neutral"], lw=7, alpha=0.35, solid_capstyle="round")
w = weight(A1)
ax.plot(100 * w, z, color=S.COLOR["nominal"], lw=3.0)

xmax = 100 * w.max()
c = S.reg("e6.A1.meanz")
# Annotations sit in a column to the right of the curve's peak. The centre line
# stops short of that column, and a bracket marks the middle half of the
# response (lower to upper quartile), so no text crosses a line or the curve.
xb = xmax * 1.18
ax.plot([0, xb], [c, c], color=S.COLOR["nominal"], lw=1.8, ls="--")
q25, q75 = S.reg("e6.A1.q25"), S.reg("e6.A1.q75")
ax.plot([xb, xb], [q25, q75], color="#222222", lw=1.6)
for zq in (q25, q75):
    ax.plot([xb - xmax * 0.03, xb], [zq, zq], color="#222222", lw=1.6)
xt = xb + xmax * 0.05
ax.text(xt, q75, "upper quartile " + S.val("e6.A1.q75", ".0f") + " km",
        va="center", fontsize=12)
ax.text(xt, c, "centre " + S.val("e6.A1.meanz", ".3f") + " km",
        va="center", fontsize=12, color=S.COLOR["nominal"])
ax.text(xt, q25, "lower quartile " + S.val("e6.A1.q25", ".0f") + " km",
        va="center", fontsize=12)
ax.text(xt, q25 - 4.5, "median " + S.val("e6.A1.q50", ".0f") + " km",
        va="center", fontsize=12, color="#444444")

mz = [S.reg(f"e6.{b}.meanz") for b, _ in BANDS]
lo_b = BANDS[int(np.argmin(mz))][0]; hi_b = BANDS[int(np.argmax(mz))][0]
ax.text(0.98, 0.03, "wide grey: absorber restricted to one wavelength band at a time;\n"
        "centre stays between " + S.val(f"e6.{lo_b}.meanz", ".3f") + " and "
        + S.val(f"e6.{hi_b}.meanz", ".3f") + " km",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=12, color="#444444")
ax.set_ylim(0, 112)
ax.set_xlim(0, xmax * 1.95)
ax.set_xlabel("Share of the radiation response at each level, per cent")
ax.set_ylabel("Altitude, km")
ax.set_title("The absorber acts in one altitude band", loc="left")
S.save(fig, "F09_E06_response-band_ch5_s10")
