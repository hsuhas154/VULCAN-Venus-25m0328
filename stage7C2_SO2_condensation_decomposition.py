import os
import pickle
import numpy as np
import chem_funs


REF = "output/Nominal_Bkzz_SO2-25m0328-reference.vul"
BASE = "output/Nominal_Bkzz_SO2-25m0328-baseline.vul"

TARGET_LEVELS = [27, 28, 29, 30, 31]
REACTIONS = [585, 587, 589]


def load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def species_index():
    return {sp: i for i, sp in enumerate(chem_funs.spec_list)}


def abundance(data, sp, level):
    y = np.asarray(data["variable"]["y"])

    idx = species_index()[sp]

    # VULCAN saved y is expected to be (nz, nspecies).
    # Handle the transposed possibility defensively.
    if y.shape[1] == len(chem_funs.spec_list):
        return float(y[level, idx])

    if y.shape[0] == len(chem_funs.spec_list):
        return float(y[idx, level])

    raise RuntimeError(
        f"Unexpected y shape {y.shape}; "
        f"number of species = {len(chem_funs.spec_list)}"
    )


def get_M(data, level):
    M = np.asarray(data["atm"]["M"])
    return float(M[level])


def get_k(data, reaction, level):
    """
    VULCAN's saved k is a dictionary/array structure indexed by
    reaction number. Handle scalar and vertical-profile forms.
    """
    k = data["variable"]["k"][reaction]

    arr = np.asarray(k)

    if arr.ndim == 0:
        return float(arr)

    if arr.ndim == 1:
        if len(arr) == 57:
            return float(arr[level])

        # Some VULCAN rate coefficients are scalar-like 1D objects.
        if len(arr) == 1:
            return float(arr[0])

    raise RuntimeError(
        f"Unexpected k[{reaction}] shape: {arr.shape}"
    )


def get_reaction(reaction):
    """
    Return reactants/products from VULCAN's generated reaction dictionary.
    """
    reactants, products = chem_funs.re_dict[reaction]
    return list(reactants), list(products)


def stoich_counts(items):
    counts = {}
    for sp in items:
        counts[sp] = counts.get(sp, 0) + 1
    return counts


def format_reaction(reaction):
    reactants, products = get_reaction(reaction)

    def fmt(items):
        counts = stoich_counts(items)
        return " + ".join(
            f"{n}{sp}" if n != 1 else sp
            for sp, n in counts.items()
        )

    return f"{fmt(reactants)} -> {fmt(products)}"


def concentration_product(data, reaction, level):
    """
    Compute the abundance/Molecular-number-density product implied
    by the reaction stoichiometry.

    VULCAN's mass-action expressions use species number densities
    represented through y and M. For a reactant species X with
    mixing ratio y_X, its number density is y_X * M.

    M is therefore included once for every reactant molecule.
    """
    reactants, _ = get_reaction(reaction)
    counts = stoich_counts(reactants)

    M = get_M(data, level)

    product = 1.0

    for sp, count in counts.items():
        y = abundance(data, sp, level)
        product *= (y * M) ** count

    return product


def print_state(label, data):
    print()
    print("=" * 120)
    print(label)
    print("=" * 120)

    print()
    print("Atmospheric state around the SO2 discrepancy")
    print("-" * 120)

    header = (
        f"{'Level':>5} "
        f"{'T(K)':>10} "
        f"{'SO2':>16} "
        f"{'SO2_l':>16} "
        f"{'H2SO4':>16} "
        f"{'H2SO4_l':>16} "
        f"{'M':>16}"
    )

    print(header)
    print("-" * 120)

    for level in TARGET_LEVELS:
        T = float(data["atm"]["Tco"][level])

        print(
            f"{level:5d} "
            f"{T:10.3f} "
            f"{abundance(data,'SO2',level):16.8e} "
            f"{abundance(data,'SO2_l',level):16.8e} "
            f"{abundance(data,'H2SO4',level):16.8e} "
            f"{abundance(data,'H2SO4_l',level):16.8e} "
            f"{get_M(data,level):16.8e}"
        )


def print_rates(label, data):
    print()
    print("=" * 120)
    print(label)
    print("=" * 120)

    for reaction in REACTIONS:
        print()
        print(f"R{reaction}: {format_reaction(reaction)}")
        print("-" * 120)

        for level in TARGET_LEVELS:
            k = get_k(data, reaction, level)

            print(
                f"Level {level:2d} : "
                f"k = {k:.8e}"
            )


def compare_states(ref, base):
    print()
    print("=" * 120)
    print("REFERENCE vs BASELINE — REACTANT STATE COMPARISON")
    print("=" * 120)

    species = ["SO2", "SO2_l", "H2SO4", "H2SO4_l"]

    for level in TARGET_LEVELS:
        print()
        print(f"LEVEL {level}")
        print("-" * 120)

        for sp in species:
            r = abundance(ref, sp, level)
            b = abundance(base, sp, level)

            diff = b - r
            rel = abs(diff) / max(abs(r), 1e-300)

            print(
                f"{sp:8s} : "
                f"reference={r:.8e}  "
                f"baseline={b:.8e}  "
                f"diff={diff:.8e}  "
                f"rel={rel:.8e}"
            )

        rM = get_M(ref, level)
        bM = get_M(base, level)

        print(
            f"{'M':8s} : "
            f"reference={rM:.8e}  "
            f"baseline={bM:.8e}  "
            f"diff={bM-rM:.8e}  "
            f"rel={abs(bM-rM)/max(abs(rM),1e-300):.8e}"
        )


def compare_k(ref, base):
    print()
    print("=" * 120)
    print("REFERENCE vs BASELINE — RATE COEFFICIENT COMPARISON")
    print("=" * 120)

    for reaction in REACTIONS:
        print()
        print(f"R{reaction}: {format_reaction(reaction)}")
        print("-" * 120)

        for level in TARGET_LEVELS:
            r = get_k(ref, reaction, level)
            b = get_k(base, reaction, level)

            diff = b - r
            rel = abs(diff) / max(abs(r), 1e-300)

            print(
                f"Level {level:2d} : "
                f"reference={r:.8e}  "
                f"baseline={b:.8e}  "
                f"diff={diff:.8e}  "
                f"rel={rel:.8e}"
            )


def compare_reconstructed_rates(ref, base):
    print()
    print("=" * 120)
    print("REACTION-RATE RECONSTRUCTION")
    print("=" * 120)

    print(
        "\nThis section compares the saved rate coefficient and "
        "reactant state against the actual rate terms seen in Stage 7C.1."
    )

    for reaction in REACTIONS:
        print()
        print(f"R{reaction}: {format_reaction(reaction)}")
        print("-" * 120)

        for level in TARGET_LEVELS:
            rk = get_k(ref, reaction, level)
            bk = get_k(base, reaction, level)

            rprod = concentration_product(ref, reaction, level)
            bprod = concentration_product(base, reaction, level)

            rcalc = rk * rprod
            bcalc = bk * bprod

            print(
                f"Level {level:2d} : "
                f"ref k={rk:.6e}  "
                f"base k={bk:.6e}  "
                f"ref reactant-product={rprod:.6e}  "
                f"base reactant-product={bprod:.6e}  "
                f"ref reconstructed={rcalc:.6e}  "
                f"base reconstructed={bcalc:.6e}"
            )


def main():
    print("=" * 120)
    print("STAGE 7C.2 — SO2 REACTANT-STATE AND RATE-COEFFICIENT DECOMPOSITION")
    print("=" * 120)

    print()
    print("Reference:", REF)
    print("Baseline :", BASE)
    print("Target levels:", TARGET_LEVELS)
    print("Target reactions:", REACTIONS)

    if not os.path.exists(REF):
        raise FileNotFoundError(REF)

    if not os.path.exists(BASE):
        raise FileNotFoundError(BASE)

    ref = load(REF)
    base = load(BASE)

    print_state("REFERENCE STATE", ref)
    print_state("BASELINE STATE", base)

    print_rates("REFERENCE RATE COEFFICIENTS", ref)
    print_rates("BASELINE RATE COEFFICIENTS", base)

    compare_states(ref, base)
    compare_k(ref, base)
    compare_reconstructed_rates(ref, base)

    print()
    print("=" * 120)
    print("STAGE 7C.2 COMPLETE")
    print("=" * 120)
    print("No VULCAN integration was performed.")
    print("No model/configuration parameters were changed.")
    print("No Unknown UV absorber was added or activated.")
    print("=" * 120)


if __name__ == "__main__":
    main()
