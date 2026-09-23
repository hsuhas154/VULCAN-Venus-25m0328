"""
Read-only four-point sensitivity analysis for Experiments 1, 2, and 3.

Important:
VULCAN .vul files are pickled dictionaries with the structure:
    data['variable']
    data['atm']
    data['parameter']

This script therefore accesses model outputs through data['variable'].
It does not modify any .vul file or the VULCAN configuration.
"""

import pickle
from pathlib import Path
import numpy as np

ROOT = Path(".")
CASES = {
    0.0: ROOT / "output/Experiment1_A0_UUV-control-25m0328.vul",
    0.5: ROOT / "output/Experiment3_A0p5_UUV-half-25m0328.vul",
    1.0: ROOT / "output/Experiment2_A1_UUV-nominal-25m0328.vul",
    2.0: ROOT / "output/Experiment3_A2_UUV-double-25m0328.vul",
}

def load_case(path: Path):
    with path.open("rb") as f:
        data = pickle.load(f)

    if not isinstance(data, dict):
        raise TypeError(f"{path}: top-level object is {type(data)}, expected dict")

    required_top = {"variable", "atm", "parameter"}
    missing = required_top - set(data)
    if missing:
        raise KeyError(f"{path}: missing top-level keys {sorted(missing)}")

    if not isinstance(data["variable"], dict):
        raise TypeError(f"{path}: data['variable'] is {type(data['variable'])}, expected dict")

    return data

print("===== EXPERIMENT 3 FOUR-POINT SENSITIVITY ANALYSIS =====")
print("Read-only analysis of saved .vul files.")
print()

print("===== FILE CHECK =====")
for amp, path in CASES.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing A={amp}: {path}")
    print(f"A={amp:g}: {path} ({path.stat().st_size:,} bytes)")

print()
print("===== LOADING CASES =====")
DATA = {amp: load_case(path) for amp, path in CASES.items()}

for amp in sorted(DATA):
    v = DATA[amp]["variable"]
    tau = v.get("tau")
    aflux = v.get("aflux")
    ymix = v.get("ymix")
    jsp = v.get("J_sp")
    print(
        f"A={amp:g}: "
        f"variable keys={len(v)}, "
        f"tau={np.shape(tau)}, "
        f"aflux={np.shape(aflux)}, "
        f"ymix={np.shape(ymix)}, "
        f"J_sp={len(jsp) if isinstance(jsp, dict) else 0}"
    )

print()
print("===== COMMON STRUCTURE CHECK =====")
VARIABLES = {amp: DATA[amp]["variable"] for amp in DATA}

required = ["tau", "aflux", "ymix", "J_sp", "species", "bins", "cross_J"]
for key in required:
    missing = [amp for amp in sorted(VARIABLES) if key not in VARIABLES[amp]]
    if missing:
        print(f"{key}: missing in A={missing}")
    else:
        shapes = {amp: np.shape(VARIABLES[amp][key]) for amp in sorted(VARIABLES)}
        print(f"{key}: {shapes}")

if not all(
    "species" in VARIABLES[amp]
    and "bins" in VARIABLES[amp]
    and "cross_J" in VARIABLES[amp]
    and "J_sp" in VARIABLES[amp]
    for amp in VARIABLES
):
    raise RuntimeError("One or more required common structures are missing.")

species_sets = [list(VARIABLES[amp]["species"]) for amp in sorted(VARIABLES)]
bins_arr = [np.asarray(VARIABLES[amp]["bins"]) for amp in sorted(VARIABLES)]

if any(s != species_sets[0] for s in species_sets[1:]):
    raise RuntimeError("Species ordering differs between cases.")

if any(not np.array_equal(b, bins_arr[0]) for b in bins_arr[1:]):
    raise RuntimeError("Wavelength grids differ between cases.")

print("Species ordering: identical.")
print("Wavelength grids: identical.")

print()
print("===== ARRAY FINITENESS CHECK =====")
for amp in sorted(VARIABLES):
    for key in ["tau", "aflux", "ymix"]:
        arr = np.asarray(VARIABLES[amp][key])
        print(
            f"A={amp:g} {key}: "
            f"finite={np.isfinite(arr).all()}, "
            f"nan={np.isnan(arr).sum()}, "
            f"inf={np.isinf(arr).sum()}"
        )

print()
print("===== SPECIES INDICES =====")
species = species_sets[0]
for sp in ["SO2", "SO3", "H2SO4", "H2SO4_l"]:
    if sp in species:
        print(f"{sp}: index {species.index(sp)}")
    else:
        print(f"{sp}: NOT PRESENT")

def bounded_relative(a, b, floor=0.0):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = np.maximum(np.abs(a), np.abs(b))
    out = np.full_like(denom, np.nan, dtype=float)
    mask = denom > floor
    out[mask] = np.abs(a[mask] - b[mask]) / denom[mask]
    return out

def describe_pair(name, a, b, label_a, label_b, meaningful_floor=0.0):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    d = a - b
    r = bounded_relative(a, b, floor=meaningful_floor)

    print(f"\n{name}: {label_a} vs {label_b}")
    print(f"  max |difference| = {np.nanmax(np.abs(d)):.6e}")
    print(f"  mean |difference| = {np.nanmean(np.abs(d)):.6e}")

    finite_r = r[np.isfinite(r)]
    if finite_r.size:
        print(f"  max bounded relative difference = {np.max(finite_r):.6e}")
        print(f"  mean bounded relative difference = {np.mean(finite_r):.6e}")

    flat_index = np.nanargmax(np.abs(d))
    idx = np.unravel_index(flat_index, d.shape)
    print(f"  location of max |difference| = {idx}")
    print(f"  {label_a} value = {a[idx]:.12e}")
    print(f"  {label_b} value = {b[idx]:.12e}")

# ----------------------------------------------------------------------
# A=0 control and amplitude cases
# ----------------------------------------------------------------------
BASE = VARIABLES[0.0]

print()
print("===== TOTAL OPTICAL DEPTH RESPONSE =====")
for amp in [0.5, 1.0, 2.0]:
    describe_pair(
        "tau",
        VARIABLES[amp]["tau"],
        BASE["tau"],
        f"A={amp:g}",
        "A=0",
        meaningful_floor=1e-12,
    )

print()
print("===== ACTINIC FLUX RESPONSE =====")
for amp in [0.5, 1.0, 2.0]:
    describe_pair(
        "aflux",
        VARIABLES[amp]["aflux"],
        BASE["aflux"],
        f"A={amp:g}",
        "A=0",
        meaningful_floor=1e-12,
    )

print()
print("===== SULFUR SPECIES RESPONSE =====")
for sp in ["SO2", "SO3", "H2SO4", "H2SO4_l"]:
    if sp not in species:
        continue
    i = species.index(sp)
    print(f"\n--- {sp} ---")
    # ymix dimensions are expected to be (nz, nspecies)
    for amp in [0.5, 1.0, 2.0]:
        describe_pair(
            sp,
            VARIABLES[amp]["ymix"][:, i],
            BASE["ymix"][:, i],
            f"A={amp:g}",
            "A=0",
            meaningful_floor=1e-30,
        )

print()
print("===== J_SP RESPONSE =====")
jsp_keys = list(BASE["J_sp"].keys())
common_jsp = [
    k for k in jsp_keys
    if all(k in VARIABLES[amp]["J_sp"] for amp in VARIABLES)
]
print(f"Common J_sp branches: {len(common_jsp)}")

def jsp_response(amp, key):
    return np.asarray(VARIABLES[amp]["J_sp"][key]) - np.asarray(BASE["J_sp"][key])

for amp in [0.5, 1.0, 2.0]:
    responses = []
    for key in common_jsp:
        d = np.asarray(jsp_response(amp, key), dtype=float)
        responses.append((float(np.nanmax(np.abs(d))), key))
    responses.sort(reverse=True)
    print(f"\nTop 10 J_sp branches by max |A={amp:g} - A=0|:")
    for value, key in responses[:10]:
        print(f"  {str(key):30s} {value:.6e}")

print()
print("===== J_SP BRANCH RELATIVE RESPONSE =====")
for amp in [0.5, 1.0, 2.0]:
    responses = []
    for key in common_jsp:
        a = np.asarray(VARIABLES[amp]["J_sp"][key], dtype=float)
        b = np.asarray(BASE["J_sp"][key], dtype=float)
        r = bounded_relative(a, b, floor=1e-30)
        finite = r[np.isfinite(r)]
        value = np.nanmax(finite) if finite.size else np.nan
        responses.append((value, key))
    responses.sort(reverse=True, key=lambda x: (-np.nan_to_num(x[0], nan=-1.0), str(x[1])))
    print(f"\nTop 10 J_sp branches by bounded relative response, A={amp:g}:")
    for value, key in responses[:10]:
        print(f"  {str(key):30s} {value:.6%}")

print()
print("===== AMPLITUDE RESPONSE / LINEARITY DIAGNOSTICS =====")
print("For each observable, compare Delta(A) with A * Delta(1).")
print("Large relative errors are not interpreted when the nominal response is numerically tiny.")

def summarize_linearity(name, arrs):
    x0 = np.asarray(arrs[0.0], dtype=float)
    d1 = np.asarray(arrs[1.0], dtype=float)
    base_response = d1 - x0
    scale = np.maximum(np.abs(d1), np.abs(x0))
    meaningful = scale > 1e-30

    print(f"\n{name}")
    for amp in [0.5, 2.0]:
        actual = np.asarray(arrs[amp], dtype=float) - x0
        expected = amp * base_response
        resid = actual - expected

        denom = np.maximum(np.abs(actual), np.abs(expected))
        valid = meaningful & (denom > 1e-30)

        if np.any(valid):
            med = np.median(np.abs(resid[valid]) / denom[valid])
            p95 = np.percentile(np.abs(resid[valid]) / denom[valid], 95)
            maxv = np.max(np.abs(resid[valid]) / denom[valid])
            print(
                f"  A={amp:g}: median normalized linearity residual = {med:.6e}, "
                f"p95 = {p95:.6e}, max = {maxv:.6e}"
            )
        else:
            print(f"  A={amp:g}: insufficient non-negligible response for a normalized test.")

summarize_linearity("tau", {amp: VARIABLES[amp]["tau"] for amp in VARIABLES})
summarize_linearity("aflux", {amp: VARIABLES[amp]["aflux"] for amp in VARIABLES})

for sp in ["SO2", "SO3", "H2SO4", "H2SO4_l"]:
    if sp in species:
        i = species.index(sp)
        summarize_linearity(
            sp,
            {amp: VARIABLES[amp]["ymix"][:, i] for amp in VARIABLES},
        )

print()
print("===== ABSOLUTE FINITE-DIFFERENCE SENSITIVITY =====")
print("S_X = [X(A_high) - X(A_low)] / [A_high - A_low]")

for name, arrs in [
    ("tau", {amp: VARIABLES[amp]["tau"] for amp in VARIABLES}),
    ("aflux", {amp: VARIABLES[amp]["aflux"] for amp in VARIABLES}),
]:
    s = (
        np.asarray(arrs[2.0], dtype=float)
        - np.asarray(arrs[0.5], dtype=float)
    ) / 1.5
    print(f"{name}: max |S_X| = {np.nanmax(np.abs(s)):.6e}, mean |S_X| = {np.nanmean(np.abs(s)):.6e}")

for sp in ["SO2", "SO3", "H2SO4", "H2SO4_l"]:
    if sp in species:
        i = species.index(sp)
        s = (
            VARIABLES[2.0]["ymix"][:, i]
            - VARIABLES[0.5]["ymix"][:, i]
        ) / 1.5
        print(f"{sp}: max |S_X| = {np.nanmax(np.abs(s)):.6e}, mean |S_X| = {np.nanmean(np.abs(s)):.6e}")

print()
print("===== ANALYSIS COMPLETE =====")
print("No source, configuration, or output file was modified.")
