"""F15. Experiment 12: the three-axis screen of every candidate the model carries.

Each species is placed by how far its column falls short of the opacity needed
(horizontal, log scale), by the share of its column inside the response band
(vertical) and, as colour, by how well its spectrum matches the empirical one.
A carrier would need to sit in the top left corner: enough column, in the band,
with the right spectrum. No species is there, and the three rankings disagree.

One definition of shortfall throughout: the required optical depth at 360 nm
(the screen's ranking wavelength) over the species' own optical depth at 360 nm
from its column above 69 km. Species with no absorption at 360 nm cannot be
placed on that axis and are listed beneath.

Sources: Experiment12/screen_results_25m0328.json; registry e12.<species>.short.360,
e12.<species>.frac_in_band_pct, e12.<species>.corr, e12.rank_nm via the JSON.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

with open(os.path.join(S.ROOT, "Experiment12", "screen_results_25m0328.json")) as fh:
    SC = json.load(fh)
NM = str(int(SC["rank_nm"]))
SP = sorted(SC["species"])
LABEL = {"NO2", "HNO2", "Cl2", "O3", "S2O2", "ClS2", "SO3", "S3", "SO2"}

def r(sp, f):
    k = f"e12.{sp}.{f}"
    return S.REG.get(k)

placed = [s for s in SP if r(s, f"short.{NM}") is not None]
missing = [s for s in SP if s not in placed]

fig = plt.figure(figsize=(16, 7.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1], wspace=0.15)
ax = fig.add_subplot(gs[:, 0])
x = np.array([r(s, f"short.{NM}") for s in placed])
y = np.array([r(s, "frac_in_band_pct") for s in placed])
c = np.array([r(s, "corr") for s in placed])
ax.axvspan(1e-1, 1, ymin=0.5, ymax=1, color="#E8F4E8", lw=0, zorder=0)
ax.text(0.15, 97, "a carrier\nwould be here", va="top", fontsize=12.5, color="#2E6B2E")
sc = ax.scatter(x, y, c=c, cmap="cividis", vmin=-1, vmax=1, s=140, edgecolor="black", zorder=3)
for s, xx, yy in zip(placed, x, y):
    if s in LABEL:
        off = {"O3": (8, -16), "S3": (8, -14), "ClS2": (8, 8)}.get(s, (8, 6))
        ax.annotate(S.chem(s), (xx, yy), xytext=off, textcoords="offset points", fontsize=13)
ax.set_xscale("log")
ax.set_xlim(1e-1, 1e15)
ax.set_ylim(-5, 108)
ax.set_xlabel("Column shortfall at " + S.jval("Experiment12/screen_results_25m0328.json",
              ["rank_nm"], ".0f") + " nm, times too little (log scale)")
ax.set_ylabel("Share of the column in the response band, per cent")
cb = fig.colorbar(sc, ax=ax, pad=0.01)
cb.set_label("Spectral match with the empirical absorber (correlation)")
ax.text(0.0, -0.17, "no absorption at this wavelength, so not placed: "
        + ", ".join(S.chem(s) for s in missing), transform=ax.transAxes, fontsize=11.5,
        va="top", wrap=True)
ax.set_title("No candidate passes all three tests", loc="left")

axr = fig.add_subplot(gs[:, 1])
axr.axis("off")
axr.set_xlim(0, 3); axr.set_ylim(-0.5, 9.6)
axr.set_title("Best eight on each test, best first", loc="left")
cols = [(f"short.{NM}", False, "smallest\nshortfall"),
        ("frac_in_band_pct", True, "most column\nin the band"),
        ("corr", True, "best spectral\nmatch")]
for ci, (key, reverse, head) in enumerate(cols):
    vals = sorted(((s_, r(s_, key)) for s_ in SP if r(s_, key) is not None),
                  key=lambda t: t[1], reverse=reverse)[:8]
    axr.text(ci + 0.05, 9.4, head, va="top", fontsize=13, fontweight="bold")
    for k, (s_, _) in enumerate(vals):
        axr.text(ci + 0.05, 7.4 - k, S.chem(s_), va="center", fontsize=13)
axr.text(0.05, -0.4, "Each test picks a different winner: the rankings disagree.",
         fontsize=12.5, style="italic")
S.save(fig, "F15_E12_three-axis-screen_ch6_s14")
