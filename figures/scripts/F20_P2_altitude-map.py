"""F20. The altitude map for Phase 2.

The same picture as F1, now with the two convective layers and the gravity waves
between them. The cloud-top convective layer, about 70 km and present near noon
under the Sun, lies inside the response band: the layer that could carry the
absorber sits where the absorber acts. That overlap is the case for Phase 2.

Sources: figures/literature_values.json (deep_convective_layer_km,
cloud_top_convective_layer_km, cloud_deck_km; Lefevre, Marcq and Lefevre 2022;
Dai et al. 2025); registry e6.interval.90.lo, e6.interval.90.hi.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

L = json.load(open(os.path.join(S.ROOT, S.LIT_PATH)))
fig, ax = plt.subplots(figsize=(9.5, 9.5))
S.response_band(ax, label=False)
ax.text(0.97, S.BAND_HI_KM - 1.2, "response band " + S.val("e6.interval.90.lo", ".0f") + " to "
        + S.val("e6.interval.90.hi", ".0f") + " km:\nwhere the absorber acts",
        transform=ax.get_yaxis_transform(), ha="right", va="top", fontsize=13, fontweight="bold")
c_lo, c_hi = L["cloud_deck_km"]["value"]
ax.axhspan(c_lo, c_hi, xmax=0.12, color="#9ECAE1", lw=0)
ax.text(0.06, c_lo - 1.0, "clouds", transform=ax.get_yaxis_transform(), ha="center", va="top",
        fontsize=12)

d_lo, d_hi = L["deep_convective_layer_km"]["value"]
ax.axhspan(d_lo, d_hi, xmin=0.16, xmax=0.98, color=S.COLOR["up"], alpha=0.22, lw=0)
ax.text(0.57, (d_lo + d_hi) / 2, "deep convective layer, " + S.litval("deep_convective_layer_km", ".0f", 0)
        + " to " + S.litval("deep_convective_layer_km", ".0f", 1) + " km,\nthicker at night",
        transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=13)

zt = L["cloud_top_convective_layer_km"]["value"]
ax.add_patch(Ellipse((0.57, zt), 0.62, 4.5, transform=ax.get_yaxis_transform(),
                     fc=S.COLOR["up"], alpha=0.45, ec=S.COLOR["up"], lw=2))
ax.text(0.57, zt + 3.0, "cloud-top convective layer, about " + S.litval("cloud_top_convective_layer_km", ".0f")
        + " km,\nnear noon under the Sun", transform=ax.get_yaxis_transform(), ha="center",
        va="bottom", fontsize=12)

ax.text(0.57, (d_hi + zt) / 2 - 1.5, "gravity waves, trapped between\nthe two convective layers",
        transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=12,
        style="italic", color="#444444")
ax.annotate("", xy=(0.92, zt - 2.5), xytext=(0.92, d_hi + 0.5), xycoords=ax.get_yaxis_transform(),
            arrowprops=dict(arrowstyle="-|>", color="#444444", lw=2, mutation_scale=18))
ax.text(0.935, (d_hi + zt) / 2, "slow\nupward\nmixing", transform=ax.get_yaxis_transform(),
        ha="left", va="center", fontsize=11, color="#444444")
ax.axhspan(-3, 0, color="#8C6D46", lw=0)
ax.text(0.06, 2, "surface", transform=ax.get_yaxis_transform(), ha="center", va="bottom",
        fontsize=12)
ax.set_ylim(-3, 112)
ax.set_xlim(0, 1.05)
ax.set_xticks([])
ax.spines["bottom"].set_visible(False)
ax.set_ylabel("Altitude, km")
ax.set_title("Phase 2: the convective layer at the cloud top sits in the band", loc="left",
             fontsize=15)
S.save(fig, "F20_P2_altitude-map_ch9_s19")
