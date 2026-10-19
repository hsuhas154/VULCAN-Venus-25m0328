"""F14. Experiment 11: the two S2Cl2 isomer cross sections used in the model.

The computed photoabsorption cross sections of SSCl2 and ClSSCl supplied by
Trabelsi and Francisco, as used. Their peaks agree with the authors' published
figure to within the values annotated (a comparison of read-off points, not a
digitised curve).

Sources: Experiment11/data/csec-sscl2.dat, Experiment11/data/csec-clsscl.dat;
registry iso.SSCl2.global_sigma, iso.SSCl2.global_lam, iso.SSCl2.val_pct,
iso.SSCl2.val_nm, and the same for ClSSCl.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import style as S
import matplotlib.pyplot as plt

ISO = [("SSCl2", "csec-sscl2.dat", "#332288"), ("ClSSCl", "csec-clsscl.dat", "#AA4499")]
fig, ax = plt.subplots(figsize=(10, 6.2))
ax.axvspan(320, 400, color=S.BAND_COLOR, lw=0, zorder=0)
ax.text(323, 0.03, S.lit("observed band\n320 to 400 nm", "the observed absorption band"),
        transform=ax.get_xaxis_transform(), ha="left", va="bottom", fontsize=11, color="#555555")
for k, (tag, f, col) in enumerate(ISO):
    d = np.loadtxt(os.path.join(S.ROOT, "Experiment11", "data", f))
    o = np.argsort(d[:, 0]); lam, sig = d[o, 0], d[o, 1]
    m = (lam >= 200) & (lam <= 520)
    ax.plot(lam[m], sig[m], color=col, lw=2.6)
    pl, ps = S.reg(f"iso.{tag}.global_lam"), S.reg(f"iso.{tag}.global_sigma")
    ax.plot([pl], [ps], "o", color=col, ms=9)
    # top right, above every curve, where the plot is empty
    ax.text(0.985, [0.97, 0.87][k], S.chem(tag) + ": peak within "
            + S.val(f"iso.{tag}.val_pct", ".3f") + " per cent and "
            + S.val(f"iso.{tag}.val_nm", ".3f") + " nm\nof the published figure",
            transform=ax.transAxes, color=col, fontsize=11, va="top", ha="right")
ax.set_yscale("log")
ax.set_ylim(1e-22, 1e-15)
ax.set_xlim(200, 520)
ax.set_xlabel("Wavelength, nm")
ax.set_ylabel("Cross section, cm$^{2}$ per molecule (log scale)")
ax.set_title("The isomer spectra as used, checked against the published figure", loc="left")
S.save(fig, "F14_E11_published-cross-sections_ch6_report")
