"""F1. Venus's atmosphere, from the surface to the top of the model.

The cloud layers and the surface temperature are from the Dai et al. (2025)
review; the response band is Experiment 6's; each marker is where a screened
species actually sits in the model (its column centroid, Experiment 12).

Sources: figures/literature_values.json (surface_temperature_K, lower_cloud_km,
middle_cloud_km, upper_cloud_km, cloud_top_range_km); registry e6.interval.90.lo,
e6.interval.90.hi, e12.<species>.centroid.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

L = json.load(open(os.path.join(S.ROOT, S.LIT_PATH)))
CLOUD = [("lower_cloud_km", "lower cloud", "#C6DBEF"),
         ("middle_cloud_km", "middle cloud", "#9ECAE1"),
         ("upper_cloud_km", "upper cloud", "#6BAED6")]
SPECIES = ["NO2", "HNO2", "Cl2", "O3", "S2O2", "ClS2", "SO3", "S3"]

fig, ax = plt.subplots(figsize=(9.5, 9.5))
S.response_band(ax, label=False)
for k, name, col in CLOUD:
    lo, hi = L[k]["value"]
    ax.axhspan(lo, hi, xmin=0.0, xmax=0.42, color=col, lw=0)
    yl, va = ((lo + hi) / 2, "center") if hi - lo > 3 else (lo - 0.8, "top")
    ax.text(0.21, yl, name + ", " + S.litval(k, ".0f", 0) + " to "
            + S.litval(k, ".0f", 1) + " km", transform=ax.get_yaxis_transform(),
            ha="center", va=va, fontsize=11.5)
t_lo, t_hi = L["cloud_top_range_km"]["value"]
ax.plot([0.44, 0.44], [t_lo, t_hi], color="#08519C", lw=3, transform=ax.get_yaxis_transform())
ax.text(0.43, t_hi + 0.6, "cloud top, " + S.litval("cloud_top_range_km", ".0f", 0) + " to "
        + S.litval("cloud_top_range_km", ".0f", 1) + " km", transform=ax.get_yaxis_transform(),
        ha="right", va="bottom", fontsize=11.5, color="#08519C")
ax.text(0.21, (S.BAND_LO_KM + S.BAND_HI_KM) / 2,
        "the absorber must be here:\nresponse band " + S.val("e6.interval.90.lo", ".0f") + " to "
        + S.val("e6.interval.90.hi", ".0f") + " km", transform=ax.get_yaxis_transform(),
        ha="center", va="center", fontsize=13, fontweight="bold")
ax.axhspan(-3, 0, color="#8C6D46", lw=0)
ax.text(0.21, 4, "surface, about " + S.litval("surface_temperature_K", ".0f") + " K",
        transform=ax.get_yaxis_transform(), ha="center", va="bottom", fontsize=12)

# where the screened candidates actually sit, spread across so labels never meet
zs = sorted(((S.reg(f"e12.{sp}.centroid"), sp) for sp in SPECIES))
for i, (zc, sp) in enumerate(zs):
    x = 0.53 + 0.058 * i
    ax.plot([x], [zc], "o", ms=10, color=S.COLOR["model"], transform=ax.get_yaxis_transform(),
            mec="black", mew=0.8)
    up = i % 2 == 1
    ax.text(x, zc + (2.0 if up else -2.0), S.chem(sp), transform=ax.get_yaxis_transform(),
            ha="center", va="bottom" if up else "top", fontsize=12)
ax.text(0.8, 30, "where screened species\nactually sit in the model\n(column centre)",
        transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=12,
        color=S.COLOR["model"])
ax.set_ylim(-3, 112)
ax.set_xticks([])
ax.spines["bottom"].set_visible(False)
ax.set_ylabel("Altitude, km")
ax.set_title("Venus's atmosphere, and where the absorber has to be", loc="left")
S.save(fig, "F01_CTX_venus-atmosphere_ch1_s3")
