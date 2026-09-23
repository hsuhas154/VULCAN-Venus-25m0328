"""
Final read-only four-point sensitivity diagnostics for VULCAN-Venus.

Uses:
    A = 0, 0.5, 1, 2

The script:
1. validates common saved structure;
2. reports absolute and bounded-relative radiative/chemical responses;
3. reports meaningful J_sp relative responses, excluding zero/undefined branches;
4. tests amplitude linearity with L2 response norms and residual norms.

No model source, configuration, or .vul output is modified.
"""

import pickle
from pathlib import Path
import numpy as np

CASES = {
    0.0: Path("output/Experiment1_A0_UUV-control-25m0328.vul"),
    0.5: Path("output/Experiment3_A0p5_UUV-half-25m0328.vul"),
    1.0: Path("output/Experiment2_A1_UUV-nominal-25m0328.vul"),
    2.0: Path("output/Experiment3_A2_UUV-double-25m0328.vul"),
}

ORDER = [0.0, 0.5, 1.0, 2.0]

def load(path):
    with path.open("rb") as f:
        d = pickle.load(f)
    if not isinstance(d, dict) or "variable" not in d:
        raise RuntimeError(f"Unexpected .vul structure: {path}")
    return d

def finite_abs_max(arr):
    arr = np.asarray(arr, dtype=float)
    good = np.isfinite(arr)
    return float(np.max(np.abs(arr[good]))) if np.any(good) else np.nan

def bounded_rel(a, b, floor):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    den = np.maximum(np.abs(a), np.abs(b))
    out = np.full_like(den, np.nan, dtype=float)
    m = np.isfinite(den) & (den > floor) & np.isfinite(a) & np.isfinite(b)
    out[m] = np.abs(a[m] - b[m]) / den[m]
    return out

def l2(arr):
    arr = np.asarray(arr, dtype=float)
    good = np.isfinite(arr)
    return float(np.linalg.norm(arr[good])) if np.any(good) else np.nan

print("===== EXPERIMENT 3 FINAL FOUR-POINT SENSITIVITY DIAGNOSTICS =====")
print("Read-only analysis.")
print()

D = {}
for A in ORDER:
    if not CASES[A].exists():
        raise FileNotFoundError(str(CASES[A]))
    D[A] = load(CASES[A])["variable"]
    print(f"A={A:g}: {CASES[A]}")

print()
print("===== COMMON STRUCTURE =====")
species = list(D[0.0]["species"])
bins = np.asarray(D[0.0]["bins"])
for A in ORDER:
    v = D[A]
    print(
        f"A={A:g}: tau={np.shape(v['tau'])}, aflux={np.shape(v['aflux'])}, "
        f"ymix={np.shape(v['ymix'])}, J_sp={len(v['J_sp'])}"
    )
    if list(v["species"]) != species:
        raise RuntimeError("Species ordering differs.")
    if not np.array_equal(np.asarray(v["bins"]), bins):
        raise RuntimeError("Wavelength grid differs.")
print("Species ordering: identical.")
print("Wavelength grid: identical.")

print()
print("===== CORE RESPONSE SUMMARY =====")
core = {
    "tau": lambda A: D[A]["tau"],
    "aflux": lambda A: D[A]["aflux"],
}
for name, getter in core.items():
    x0 = np.asarray(getter(0.0), dtype=float)
    print(f"\n{name}")
    for A in [0.5, 1.0, 2.0]:
        x = np.asarray(getter(A), dtype=float)
        delta = x - x0
        rel = bounded_rel(x, x0, floor=1e-20)
        good = np.isfinite(rel)
        print(
            f"  A={A:g}: max|d|={np.nanmax(np.abs(delta)):.6e}, "
            f"mean|d|={np.nanmean(np.abs(delta)):.6e}, "
            f"max bounded rel={np.nanmax(rel):.6e}, "
            f"mean bounded rel={np.nanmean(rel):.6e}"
        )

print()
print("===== SULFUR SPECIES RESPONSE =====")
for sp in ["SO2", "SO3", "H2SO4", "H2SO4_l"]:
    i = species.index(sp)
    print(f"\n{sp}")
    x0 = np.asarray(D[0.0]["ymix"][:, i], dtype=float)
    for A in [0.5, 1.0, 2.0]:
        x = np.asarray(D[A]["ymix"][:, i], dtype=float)
        delta = x - x0
        rel = bounded_rel(x, x0, floor=1e-30)
        print(
            f"  A={A:g}: max|d|={np.max(np.abs(delta)):.6e}, "
            f"mean|d|={np.mean(np.abs(delta)):.6e}, "
            f"max bounded rel={np.nanmax(rel):.6e}, "
            f"mean bounded rel={np.nanmean(rel):.6e}"
        )

print()
print("===== J_SP ABSOLUTE RESPONSE =====")
common = [k for k in D[0.0]["J_sp"] if all(k in D[A]["J_sp"] for A in ORDER)]
print(f"Common branches: {len(common)}")
for A in [0.5, 1.0, 2.0]:
    rows = []
    for k in common:
        a = np.asarray(D[A]["J_sp"][k], dtype=float)
        b = np.asarray(D[0.0]["J_sp"][k], dtype=float)
        if not np.isfinite(a).any() or not np.isfinite(b).any():
            continue
        rows.append((finite_abs_max(a-b), k))
    rows.sort(reverse=True, key=lambda x: x[0])
    print(f"\nTop 10 by max |J(A={A:g})-J(0)|:")
    for value, k in rows[:10]:
        print(f"  {str(k):30s} {value:.6e}")

print()
print("===== J_SP MEANINGFUL RELATIVE RESPONSE =====")
print("Branches that are zero/undefined over the useful range are excluded.")
for A in [0.5, 1.0, 2.0]:
    rows = []
    for k in common:
        a = np.asarray(D[A]["J_sp"][k], dtype=float)
        b = np.asarray(D[0.0]["J_sp"][k], dtype=float)

        scale = np.nanmax(
            np.concatenate([
                np.abs(a[np.isfinite(a)]),
                np.abs(b[np.isfinite(b)])
            ])
        ) if (np.isfinite(a).any() or np.isfinite(b).any()) else 0.0

        if not np.isfinite(scale) or scale <= 0.0:
            continue

        # Require a response denominator that is not merely machine-level.
        floor = max(1e-30, 1e-8 * scale)
        r = bounded_rel(a, b, floor=floor)
        good = np.isfinite(r)
        if not np.any(good):
            continue

        rows.append((
            float(np.max(r[good])),
            float(np.mean(r[good])),
            k,
            scale
        ))

    rows.sort(reverse=True, key=lambda x: x[0])
    print(f"\nTop 15 branches by max bounded relative response, A={A:g}:")
    for mx, mean, k, scale in rows[:15]:
        print(
            f"  {str(k):30s} max={mx:.6%} mean={mean:.6%} "
            f"branch-scale={scale:.6e}"
        )

print()
print("===== L2 AMPLITUDE LINEARITY TEST =====")
print("For each field, compare ||Delta(A)|| with A*||Delta(1)||.")
print("Residual ratio = ||Delta(A)-A*Delta(1)|| / ||A*Delta(1)||.")

fields = {
    "tau": lambda A: D[A]["tau"],
    "aflux": lambda A: D[A]["aflux"],
    "SO2": lambda A: D[A]["ymix"][:, species.index("SO2")],
    "SO3": lambda A: D[A]["ymix"][:, species.index("SO3")],
    "H2SO4": lambda A: D[A]["ymix"][:, species.index("H2SO4")],
    "H2SO4_l": lambda A: D[A]["ymix"][:, species.index("H2SO4_l")],
}

for name, getter in fields.items():
    x0 = np.asarray(getter(0.0), dtype=float)
    d1 = np.asarray(getter(1.0), dtype=float) - x0
    n1 = l2(d1)

    print(f"\n{name}: ||Delta(1)|| = {n1:.6e}")

    if not np.isfinite(n1) or n1 <= 0:
        print("  Insufficient nominal response for a linearity test.")
        continue

    for A in [0.5, 2.0]:
        dA = np.asarray(getter(A), dtype=float) - x0

        good = np.isfinite(dA) & np.isfinite(d1)
        if not np.any(good):
            print(f"  A={A:g}: no finite overlap.")
            continue

        actual_norm = np.linalg.norm(dA[good])
        expected = A * d1[good]
        expected_norm = np.linalg.norm(expected)

        ratio = actual_norm / expected_norm if expected_norm > 0 else np.nan
        residual = np.linalg.norm(dA[good] - expected)
        residual_ratio = residual / expected_norm if expected_norm > 0 else np.nan

        print(
            f"  A={A:g}: ||Delta(A)||/(A||Delta(1)||)={ratio:.6f}, "
            f"residual ratio={residual_ratio:.6f}"
        )

print()
print("===== ACTINIC FLUX: MOST RESPONSIVE LOCATION =====")
x0 = np.asarray(D[0.0]["aflux"], dtype=float)
d1 = np.asarray(D[1.0]["aflux"], dtype=float) - x0
idx = np.unravel_index(np.nanargmax(np.abs(d1)), d1.shape)
print(f"Location selected from max A=1 response: {idx}")
for A in ORDER:
    value = np.asarray(D[A]["aflux"])[idx]
    print(f"  A={A:g}: aflux={value:.12e}")

print()
print("===== CONCLUSION-READY DIAGNOSTICS =====")
print("1. The four cases have identical saved array dimensions and wavelength/species grids.")
print("2. The actinic-flux response can be tested with a global norm and a fixed sensitive location.")
print("3. Total optical depth includes the re-equilibrated gas state, so it is not a pure UUV scaling diagnostic.")
print("4. Sulfur mixing-ratio responses are reported both absolutely and with bounded relative changes.")
print("5. J_sp relative rankings exclude branches whose denominator is effectively zero.")
print()
print("No files were modified.")
