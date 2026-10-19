"""F8. Experiment 5: which reactions carry the response.

Each species' chemical tendency in VULCAN is a sum of signed reaction terms.
Experiment 5 rebuilt that sum term by term for the nominal (A = 1) and control
(A = 0) runs. For every reaction, the bar is the largest change, over altitude,
in that reaction's contribution to the sulfur dioxide tendency. For a photolysis
term the contribution is J times the parent density, so this is the d(J y)
convention, not dJ. The twelve largest are shown, labelled with the reaction as
written in the network.

Sources: registry e5.SO2.R<id>.maxdelta and e5.SO2.R<id>.level; reaction text
from the network as stored in the output (variable "Rf").
"""
import os, sys, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

SP, TOP = "SO2", 12
sol = S.vul("Experiment2_A1_UUV-nominal-25m0328.vul")
Rf = sol["variable"]["Rf"]
z = np.asarray(sol["atm"]["zmco"]) / 1e5
keys = sorted([k for k in S.REG if k.startswith(f"e5.{SP}.R") and k.endswith(".maxdelta")],
              key=lambda k: -S.REG[k])[:TOP]

def label(rid):
    """Reaction text with species typeset, e.g. COS + SO3 + SO3 gives SO$_3$."""
    txt = Rf[rid].replace("->", "gives").replace("_l", " (liquid)")
    return " ".join(S.chem(t) if re.match(r"^[A-Z]", t) else t for t in txt.split())

fig, ax = plt.subplots(figsize=(12, 7.5))
y = np.arange(len(keys))[::-1]
for yi, k in zip(y, keys):
    rid = int(k.split(".")[2][1:])
    photo = len(Rf[rid].split("->")[0].split("+")) == 1   # one reactant: photolysis
    ax.barh(yi, S.REG[k], color=S.COLOR["accent"] if photo else S.COLOR["nominal"], height=0.65)
    ax.text(S.REG[k] * 1.15, yi, label(rid), va="center", fontsize=12)
ax.set_xscale("log")
ax.set_yticks(y)
ax.set_yticklabels([S.lit(k.split(".")[2], "reaction number in the VULCAN network") for k in keys],
                   fontsize=12, color="#444444")
ax.set_xlim(S.REG[keys[-1]] / 2, S.REG[keys[0]] * 1e3)
ax.set_xlabel("Largest change in the reaction's contribution to the " + S.chem(SP)
              + " tendency,\nmolecules per cm$^{3}$ per second (log scale; d(J y) for photolysis)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=S.COLOR["nominal"], label="chemical reaction"),
                   Patch(color=S.COLOR["accent"], label="photolysis")],
          loc="lower right", frameon=False)
ax.set_title("Reactions whose contribution to " + S.chem(SP) + " changes most when the absorber is on",
             loc="left")
S.save(fig, "F08_E05_reaction-pathways_ch5_report")
