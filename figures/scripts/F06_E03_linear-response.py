"""F6. Experiment 3: the response is proportional to the amount of absorber.

The absorber amplitude A scaled to 0, 0.5, 1 and 2. Over the absorber's own
wavelength range, 300 to 698 nm, the change in optical depth and in actinic
flux grows in proportion to A.

Sources: output/Experiment1_A0_*.vul, output/Experiment3_A0p5_*.vul,
output/Experiment2_A1_*.vul, output/Experiment3_A2_*.vul; registry
e3.lin.tau.inband.median, e3.lin.aflux.inband.median.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

RUNS = [(0.0, "Experiment1_A0_UUV-control-25m0328.vul"),
        (0.5, "Experiment3_A0p5_UUV-half-25m0328.vul"),
        (1.0, "Experiment2_A1_UUV-nominal-25m0328.vul"),
        (2.0, "Experiment3_A2_UUV-double-25m0328.vul")]
sols = {a: S.vul(f) for a, f in RUNS}
z, bins, quad = S.grid(sols[0.0])
inb = (bins >= 300.0) & (bins <= 698.0)

def size(a, field):
    """Median over every in-range element of the change at A, relative to the
    change at A = 1, which is the Experiment 3 definition of linearity."""
    f0 = np.asarray(sols[0.0]["variable"][field])
    num = np.abs(np.asarray(sols[a]["variable"][field]) - f0)
    den = np.abs(np.asarray(sols[1.0]["variable"][field]) - f0)
    m = (den > 0) & inb[None, :]
    return np.median(num[m] / den[m])

fig, ax = plt.subplots(figsize=(8, 6.5))
pts = {}
amps = [a for a, _ in RUNS]
ax.plot([0, 2.1], [0, 2.1], color=S.COLOR["neutral"], ls="--", lw=1.6)
ax.text(1.55, 1.42, "exactly proportional", rotation=38, color=S.COLOR["neutral"], fontsize=12,
        ha="center")
for field, col, mk, name in (("tau", S.COLOR["nominal"], "o", "optical depth"),
                              ("aflux", S.COLOR["control"], "s", "actinic flux")):
    y = [size(a, field) for a in amps]
    ax.plot(amps, y, mk, color=col, ms=11, mfc="none" if field == "aflux" else col, mew=2)
    pts[field] = y
lbl = (S.lit("each point: median over the 300 to 698 nm range",
             "the absorber's wavelength range, a definition") + "\n"
       "half amount against half the nominal effect:\n"
       "optical depth " + S.val("e3.lin.tau.inband.median", ".4f")
       + ",  actinic flux " + S.val("e3.lin.aflux.inband.median", ".4f"))
ax.text(0.03, 0.97, lbl, transform=ax.transAxes, va="top", fontsize=12.5)
ax.plot([], [], "o", color=S.COLOR["nominal"], ms=10, label="optical depth")
ax.plot([], [], "s", color=S.COLOR["control"], mfc="none", mew=2, ms=10, label="actinic flux")
ax.legend(loc="upper left", frameon=False, bbox_to_anchor=(0.0, 0.83))
ax.set_xlim(-0.1, 2.2); ax.set_ylim(-0.1, 2.3)
ax.set_xlabel("Absorber amount, multiple of the nominal")
ax.set_ylabel("Response, multiple of the nominal response")
ax.set_title("Twice the absorber, twice the effect", loc="left")

# Magnified view of the doubled case. Optical depth is exactly proportional to
# the amount by construction; actinic flux falls slightly short of double,
# because light is absorbed exponentially, so the open square sits below the
# filled circle here and only here. The inset's own tick labels give the size.
ins = fig.add_axes([0.66, 0.22, 0.22, 0.27])
ins.plot([1.9, 2.1], [1.9, 2.1], color=S.COLOR["neutral"], ls="--", lw=1.4)
ins.plot([2.0], [pts["tau"][-1]], "o", color=S.COLOR["nominal"], ms=10)
ins.plot([2.0], [pts["aflux"][-1]], "s", color=S.COLOR["control"], mfc="none", mew=2, ms=10)
ins.set_xlim(1.95, 2.05); ins.set_ylim(1.93, 2.02)
ins.set_xticks([2.0]); ins.tick_params(labelsize=11)
ins.set_title("at twice the amount, magnified", fontsize=11)

S.save(fig, "F06_E03_linear-response_ch5_s8")
