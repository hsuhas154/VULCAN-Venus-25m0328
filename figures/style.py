"""Shared style, registry access and number bookkeeping for every Phase 1 figure.

Every number printed on a figure must come from verify/registry.json or from a
committed result file, through val() or jval() below. save() then walks every
piece of text in the figure and refuses to write the files if it finds a digit
that did not come through one of those routes. check_figures.py later re-reads
the registry and confirms each printed value still matches it.

Run every figure script from the repository root, after the registry has been
built with: python verify/recompute_25m0328.py
"""
import os, re, io, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG_DIR = os.path.join(ROOT, "figures")
MANIFEST_DIR = os.path.join(FIG_DIR, "manifest")
REG_PATH = os.path.join(ROOT, "verify", "registry.json")

if not os.path.exists(REG_PATH):
    raise SystemExit("verify/registry.json not found. Run: python verify/recompute_25m0328.py")
with open(REG_PATH) as fh:
    REG = json.load(fh)

# ---------------------------------------------------------------- appearance
# DejaVu ships with matplotlib, so the figures look identical on every machine.
plt.rcParams.update({
    "font.family": "DejaVu Serif",
    "mathtext.fontset": "dejavuserif",
    "font.size": 14,
    "axes.labelsize": 15,
    "axes.titlesize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 13,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "axes.unicode_minus": False,      # a hyphen for minus signs, never a dash
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
})

# One colour per meaning, fixed across the whole set of twenty figures.
# Okabe and Ito palette: colourblind safe, and the four stay distinct in greyscale.
COLOR = {
    "control":   "#000000",   # A = 0, no absorber
    "nominal":   "#0072B2",   # the nominal absorber, A = 1
    "empirical": "#E69F00",   # the empirical absorber profile or spectrum
    "model":     "#CC79A7",   # a profile the model's own chemistry produced
    "neutral":   "#7F7F7F",
    "accent":    "#009E73",
    "up":        "#D55E00",   # Experiment 7, absorber moved up
    "down":      "#56B4E9",   # Experiment 7, absorber moved down
}
BAND_COLOR = "#DDDDDD"

# The response band is the narrowest altitude interval holding 90 per cent of the
# Experiment 6 response (registry e6.interval.90.lo and .hi: 69 and 93 km).
BAND_LO_KM, BAND_HI_KM = REG["e6.interval.90.lo"], REG["e6.interval.90.hi"]

def response_band(ax, label=True, fontsize=12, where="middle"):
    ax.axhspan(BAND_LO_KM, BAND_HI_KM, color=BAND_COLOR, zorder=0, lw=0)
    if label:
        y, va = {"middle": ((BAND_LO_KM + BAND_HI_KM) / 2, "center"),
                 "bottom": (BAND_LO_KM + 0.8, "bottom"),
                 "above": (BAND_HI_KM + 0.8, "bottom")}[where]
        ax.text(0.98, y,
                "response band " + val("e6.interval.90.lo", ".0f") + " to "
                + val("e6.interval.90.hi", ".0f") + " km",
                transform=ax.get_yaxis_transform(), ha="right", va=va,
                fontsize=fontsize, color="#555555")

# ---------------------------------------------------------------- model output
import pickle
import numpy as np

def vul(name):
    """Load a committed VULCAN output from output/ (read only)."""
    with open(os.path.join(ROOT, "output", name), "rb") as fh:
        return pickle.load(fh)

def grid(sol):
    """Altitude of each layer in km and the quadrature weight of each bin, nm."""
    z = np.asarray(sol["atm"]["zmco"]) / 1e5
    bins = np.asarray(sol["variable"]["bins"])
    # VULCAN's two bin widths, 0.1 nm below the split at 240 nm and 2 nm above
    quad = np.where(bins < 240.0, REG["e8.base.v.dbin1"], REG["e8.base.v.dbin2"])
    return z, bins, quad

# ---------------------------------------------------------------- chemical names
_ATOM = re.compile(r"([A-Z][a-z]?)(\d+)")

def chem(name):
    """A species name with its counts as subscripts, for every label in the set:
    chem("ClS2") gives ClS$_2$, chem("H2SO4") gives H$_2$SO$_4$. The guard in
    save() reads subscript digits as part of a name, never as a value."""
    return _ATOM.sub(lambda m: f"{m.group(1)}$_{{{m.group(2)}}}$", name)

# ---------------------------------------------------------------- number bookkeeping
_VALUES = []     # (key, fmt, scale, printed string)
_LITERALS = []   # (printed string, reason)

def reg(key):
    if key not in REG:
        raise KeyError(f"registry key {key!r} not found; rebuild the registry")
    return REG[key]

def val(key, fmt, scale=1.0):
    """Format a registry value for printing and record it for the check."""
    s = format(reg(key) * scale, fmt)
    _VALUES.append(dict(source="registry", key=key, fmt=fmt, scale=scale, text=s))
    return s

def jval(path, keys, fmt, scale=1.0):
    """Format a value from a committed JSON result file and record it."""
    with open(os.path.join(ROOT, path)) as fh:
        d = json.load(fh)
    v = d
    for k in keys:
        v = v[k]
    s = format(v * scale, fmt)
    _VALUES.append(dict(source=path, key=list(keys), fmt=fmt, scale=scale, text=s))
    return s

def sci(key, digits, scale=1.0):
    """A registry value in scientific notation, as mathtext: 4.0153 x 10^-8."""
    v = reg(key) * scale
    m, e = format(v, f".{digits}e").split("e")
    s = f"${m}\\times10^{{{int(e)}}}$"
    _VALUES.append(dict(source="registry", key=key, fmt=f"sci{digits}", scale=scale, text=s))
    return s

# Quantities derived from registry keys rather than stored,
# computed here and re-computed the same way by check_figures.py.
DERIVED = {
    "max/min": lambda v: max(v) / min(v),
    "a/b":     lambda v: v[0] / v[1],
}

LIT_PATH = "figures/literature_values.json"

def _resolve(k):
    """A registry key, or "lit:<id>" / "lit:<id>:<index>" for a literature value."""
    if k.startswith("lit:"):
        parts = k.split(":")
        with open(os.path.join(ROOT, LIT_PATH)) as fh:
            v = json.load(fh)[parts[1]]["value"]
        return v[int(parts[2])] if len(parts) > 2 else v
    return reg(k)

def litval(lid, fmt, index=None, scale=1.0):
    """A literature value from figures/literature_values.json, recorded like jval."""
    keys = [lid, "value"] + ([index] if index is not None else [])
    return jval(LIT_PATH, keys, fmt, scale)

def sci_lit(lid, digits, index=None):
    """A literature value in scientific notation, as mathtext, recorded like sci()."""
    keys = [lid, "value"] + ([index] if index is not None else [])
    with open(os.path.join(ROOT, LIT_PATH)) as fh:
        v = json.load(fh)
    for k in keys:
        v = v[k]
    m, e = format(v, f".{digits}e").split("e")
    s = f"${m}\\times10^{{{int(e)}}}$"
    _VALUES.append(dict(source=LIT_PATH, key=keys, fmt=f"sci{digits}", scale=1.0, text=s))
    return s

def dval(op, keys, fmt):
    v = DERIVED[op]([_resolve(k) for k in keys])
    s = format(v, fmt)
    _VALUES.append(dict(source="derived", op=op, key=list(keys), fmt=fmt, scale=1.0, text=s))
    return s

def lit(text, reason):
    """Text with digits that are not data (a definition, a label), with the reason."""
    _LITERALS.append(dict(text=text, reason=reason))
    return text

# Chemical formulas carry digits that are not numbers: ClS2, FeCl3, H2SO4. So do
# unit exponents and subscripts in mathtext: cm$^{-1}$, SO$_2$.
_FORMULA = re.compile(r"\b(?:[A-Z][a-z]?\d*){1,6}\b")
# Experiment and phase numbers are names, not data: "Exp 7", "Experiments 8 to 11", "Phase 2".
_EXPNO = re.compile(r"\b(?:Exp|Experiments?|Phase)\s+\d+[AB]?(?:\s+(?:to|and)\s+\d+[AB]?)?")
_SUPSUB = re.compile(r"[\^_]\{?-?\d+\}?")

def _residual_digits(s):
    for v in sorted({x["text"] for x in _VALUES} | {x["text"] for x in _LITERALS},
                    key=len, reverse=True):
        s = s.replace(v, " ")
    s = _SUPSUB.sub(" ", s)
    s = _EXPNO.sub(" ", s)
    s = _FORMULA.sub(lambda m: " " if re.search(r"[A-Z]", m.group(0)) else m.group(0), s)
    return re.findall(r"\d[\d.,]*", s)

DASHES = (chr(0x2013), chr(0x2014), chr(0x2212))  # en dash, em dash, minus sign

def _texts(fig):
    ticks = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                ticks.update((id(tk.label1), id(tk.label2)))
        for axis in (ax.xaxis, ax.yaxis):
            ticks.add(id(axis.get_offset_text()))
    out = []
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if s and id(t) not in ticks:
            out.append(s)
    return out

def save(fig, name):
    """Write PDF and 300 dpi PNG, after refusing any unsourced number or dash."""
    fig.canvas.draw()
    texts = _texts(fig)
    problems = []
    for s in texts:
        if any(d in s for d in DASHES):
            problems.append(f"dash in text: {s!r}")
        bad = _residual_digits(s)
        if bad:
            problems.append(f"unsourced number(s) {bad} in text: {s!r}")
    # every recorded value must actually appear in the drawn text
    alltext = "\n".join(texts)
    for v in _VALUES:
        if v["text"] not in alltext:
            problems.append(f"value {v['text']} ({v['key']}) recorded but not drawn")
    if problems:
        raise SystemExit("FIGURE REFUSED, " + name + ":\n  " + "\n  ".join(problems))
    os.makedirs(MANIFEST_DIR, exist_ok=True)
    fig.savefig(os.path.join(FIG_DIR, name + ".pdf"))
    fig.savefig(os.path.join(FIG_DIR, name + ".png"), dpi=300)
    with open(os.path.join(MANIFEST_DIR, name + ".json"), "w") as fh:
        json.dump(dict(figure=name, values=_VALUES, literals=_LITERALS, texts=texts),
                  fh, indent=1)
    print(f"wrote figures/{name}.pdf and .png; {len(_VALUES)} sourced values, "
          f"{len(_LITERALS)} literal labels")
