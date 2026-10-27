#!/usr/bin/env python3
"""Figure-level check: every number printed on a figure still equals its source.

For each figure, save() recorded which registry key or result-file entry each
printed number came from, the format it was printed with, and all text drawn.
This script re-reads the CURRENT registry and result files, formats each value
again, and confirms the printed text is identical and still present. It also
fails on any en dash, em dash or minus-sign character. The figure equivalent
of verify/check_documents_25m0328.py (report Section 8.2).

Run from the repository root: python figures/check_figures.py
"""
import os, sys, json, glob

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REG = json.load(open(os.path.join(ROOT, "verify", "registry.json")))

def fmt(v, f, scale):
    v = v * scale
    if f.startswith("sci"):
        m, e = format(v, f".{int(f[3:])}e").split("e")
        return f"${m}\\times10^{{{int(e)}}}$"
    return format(v, f)

DERIVED = {
    "max/min": lambda v: max(v) / min(v),
    "a/b":     lambda v: v[0] / v[1],
}

def _resolve(k):
    if k.startswith("lit:"):
        parts = k.split(":")
        v = json.load(open(os.path.join(ROOT, "figures", "literature_values.json")))[parts[1]]["value"]
        return v[int(parts[2])] if len(parts) > 2 else v
    return REG[k]

def source_value(rec):
    if rec["source"] == "derived":
        return DERIVED[rec["op"]]([_resolve(k) for k in rec["key"]])
    if rec["source"] == "registry":
        return REG[rec["key"]]
    d = json.load(open(os.path.join(ROOT, rec["source"])))
    for k in rec["key"]:
        d = d[k]
    return d

total = bad = 0
for mf in sorted(glob.glob(os.path.join(ROOT, "figures", "manifest", "*.json"))):
    m = json.load(open(mf))
    name, alltext = m["figure"], "\n".join(m["texts"])
    errs = []
    if not os.path.exists(os.path.join(ROOT, "figures", name + ".pdf")):
        errs.append("PDF missing")
    for rec in m["values"]:
        total += 1
        try:
            now = fmt(source_value(rec), rec["fmt"], rec["scale"])
        except KeyError as e:
            errs.append(f"source key missing: {rec['key']}"); continue
        if now != rec["text"]:
            errs.append(f"{rec['key']}: printed {rec['text']}, source now gives {now}")
        if rec["text"] not in alltext:
            errs.append(f"{rec['key']}: {rec['text']} not found in the figure text")
    for t in m["texts"]:
        if any(c in t for c in (chr(0x2013), chr(0x2014), chr(0x2212))):
            errs.append(f"dash in text: {t!r}")
    bad += len(errs)
    print(f"{name:52s} {len(m['values']):3d} values  {'OK' if not errs else 'FAIL'}")
    for e in errs:
        print("     ", e)
print(f"\nTOTAL {total} printed values checked, {bad} problems")
sys.exit(1 if bad else 0)
