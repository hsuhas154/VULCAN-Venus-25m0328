"""F16. Experiment 13: the iron-sulfur minerals cannot reach the requirement.

The mass of each mineral a litre of droplet would need to hold to meet the
published bulk absorbance requirement at 375 nm. Both need more than the whole
droplet weighs. Reading the measured optical depth the other way (napierian)
makes it worse, not better. The ferric chloride line is a check that the method
reproduces the published figure, not a result.

Sources: Jiang et al. (2024) extracted spectra in Experiment13/data; registry
e13.*.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12, 6.2))
rows = [("rhomboclase", "rhomboclase"), ("acid_ferric_sulfate", "acid ferric sulfate")]
y0 = [4.2, 2.6]
for (k, name), yy in zip(rows, y0):
    v = S.reg(f"e13.{k}.loading_kg_per_L"); vn = S.reg(f"e13.{k}.loading_napierian")
    ax.barh(yy + 0.3, v, height=0.55, color=S.COLOR["nominal"])
    ax.barh(yy - 0.35, vn, height=0.55, color=S.COLOR["nominal"], alpha=0.35)
    ax.text(v + 0.15, yy + 0.3, S.val(f"e13.{k}.loading_kg_per_L", ".2f") + " kg per litre",
            va="center", fontsize=13)
    ax.text(vn + 0.15, yy - 0.35, S.val(f"e13.{k}.loading_napierian", ".2f")
            + " kg per litre if read as napierian", va="center", fontsize=12, color="#444444")
    ax.text(-0.2, yy, name, ha="right", va="center", fontsize=14)
rho = S.reg("e13.rho_solution")
ax.axvline(rho, color=S.COLOR["control"], lw=3.2)
ax.text(rho + 0.1, 5.25, "the whole droplet weighs " + S.val("e13.rho_solution", ".1f")
        + " kg per litre", fontsize=13, va="bottom")
# method calibration, visually separated
ax.axhspan(0.0, 1.4, color="#F2F2F2", lw=0, zorder=0)
ax.barh(0.7, S.reg("e13.fecl3.loading_kg_per_L"), height=0.5, color="none",
        edgecolor=S.COLOR["control"], hatch="//", lw=1.5)
ax.text(rho + 0.15, 0.7, "ferric chloride, aqueous: "
        + S.val("e13.fecl3.loading_kg_per_L", ".3f") + " against the published "
        + S.val("e13.fecl3.published_loading", ".1f") + " kg per litre", va="center", fontsize=12.5)
ax.text(-0.2, 0.7, "method\ncalibration", ha="right", va="center", fontsize=13, style="italic")
ax.set_ylim(-0.1, 5.9)
ax.set_xlim(0, 11.5)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlabel("Mineral needed per litre of droplet, kg")
ax.set_title("Both minerals would have to outweigh the droplet that carries them",
             loc="left")
S.save(fig, "F16_E13_minerals-cannot-reach_ch6_s15")
