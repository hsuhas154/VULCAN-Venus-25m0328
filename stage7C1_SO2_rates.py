import os
import re
import ast
import inspect
import pickle
import numpy as np

import chem_funs


REF = "output/Nominal_Bkzz_SO2-25m0328-reference.vul"
BASE = "output/Nominal_Bkzz_SO2-25m0328-baseline.vul"

TARGET_LEVELS = [27, 28, 29, 30, 31]


def load_vul(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def extract_so2_rate_terms():
    """
    Extract the k-index associated with each individual SO2
    production/loss term from chem_funs.rate_ans().
    """
    source = inspect.getsource(chem_funs.rate_ans)
    tree = ast.parse(source)

    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue

        for target in node.targets:
            if not isinstance(target, ast.Subscript):
                continue

            # Recover the species name from rate_str["SO2"]
            try:
                key = target.slice.value
            except AttributeError:
                continue

            if key != "SO2":
                continue

            if not isinstance(node.value, ast.List):
                raise RuntimeError("SO2 rate expression is not a list.")

            terms = []

            for term in node.value.elts:
                text = ast.unparse(term)

                matches = re.findall(r"k\[(\d+)\]", text)

                if len(matches) != 1:
                    raise RuntimeError(
                        f"Could not uniquely identify k index in: {text}"
                    )

                k_index = int(matches[0])

                # The generated expressions have +1* or -1*
                if re.search(r"^\s*-\s*1\s*\*", text):
                    sign = -1
                elif re.search(r"^\s*\+\s*1\s*\*", text):
                    sign = 1
                else:
                    # Handle expressions such as +1*k[...] without
                    # relying exclusively on whitespace.
                    if text.lstrip().startswith("-1*"):
                        sign = -1
                    else:
                        sign = 1

                terms.append((k_index, sign, text))

            return terms

    raise RuntimeError("Could not locate rate_str['SO2'].")


def reaction_equation(k_index):
    """
    chem_funs.re_dict contains the reactants/products for every
    forward and reverse reaction.
    """
    reac, prod = chem_funs.re_dict[k_index]

    def fmt(items):
        if not items:
            return "—"

        counts = {}
        for sp in items:
            counts[sp] = counts.get(sp, 0) + 1

        parts = []
        for sp, count in counts.items():
            if count == 1:
                parts.append(sp)
            else:
                parts.append(f"{count}{sp}")

        return " + ".join(parts)

    return f"{fmt(reac)} -> {fmt(prod)}"


def reaction_name(k_index, Rf):
    """
    Rf is indexed by odd/forward reaction number.
    Even k corresponds to the reverse of the preceding odd reaction.
    """
    odd = k_index if k_index % 2 == 1 else k_index - 1

    base = Rf.get(odd, f"<R{odd} unavailable>")

    if k_index % 2 == 1:
        return base
    else:
        return "reverse: " + base


def evaluate_terms(data, term_info):
    """
    Evaluate VULCAN's individual SO2 rate expressions.

    rate_ans() expects:
        y = species x altitude
        k = rate coefficient dictionary
        M = number density profile
    """
    variable = data["variable"]
    atm = data["atm"]

    y = np.asarray(variable["y"])
    if y.shape[0] == len(chem_funs.spec_list):
        y_species = y
    else:
        y_species = y.T

    chem_funs.y = y_species
    chem_funs.k = variable["k"]
    chem_funs.M = np.asarray(atm["M"])

    raw_terms = np.asarray(chem_funs.rate_ans("SO2"), dtype=float)

    if raw_terms.ndim == 1:
        raw_terms = raw_terms[np.newaxis, :]

    if raw_terms.shape[0] != len(term_info):
        raise RuntimeError(
            f"Term count mismatch: evaluated {raw_terms.shape[0]}, "
            f"but extracted {len(term_info)}."
        )

    return raw_terms


def summarize(label, data, term_info):
    Rf = data["variable"]["Rf"]

    values = evaluate_terms(data, term_info)

    print()
    print("=" * 120)
    print(f"{label} — SO2 REACTION-LEVEL RATES")
    print("=" * 120)

    print(f"Number of SO2 terms : {values.shape[0]}")
    print(f"Number of levels    : {values.shape[1]}")

    print()
    print("TARGET LEVELS")
    print("-" * 120)

    for level in TARGET_LEVELS:
        total = np.sum(values[:, level])

        production = np.sum(values[values[:, level] > 0, level])
        loss = np.sum(values[values[:, level] < 0, level])

        print(
            f"Level {level:2d} : "
            f"production = {production:.8e}   "
            f"loss = {loss:.8e}   "
            f"net = {total:.8e}"
        )

    print()
    print("TOP SO2 PRODUCTION TERMS")
    print("-" * 120)

    for level in TARGET_LEVELS:
        vals = values[:, level]
        order = np.argsort(vals)[::-1]

        print()
        print(f"LEVEL {level}")
        print("-" * 120)

        shown = 0
        for idx in order:
            if vals[idx] <= 0:
                continue

            k_index, sign, expression = term_info[idx]

            print(
                f"R{k_index:3d}  "
                f"rate={vals[idx]:.8e}  "
                f"{reaction_name(k_index, Rf)}"
            )

            shown += 1
            if shown == 10:
                break

    print()
    print("TOP SO2 LOSS TERMS")
    print("-" * 120)

    for level in TARGET_LEVELS:
        vals = values[:, level]
        order = np.argsort(vals)

        print()
        print(f"LEVEL {level}")
        print("-" * 120)

        shown = 0
        for idx in order:
            if vals[idx] >= 0:
                continue

            k_index, sign, expression = term_info[idx]

            print(
                f"R{k_index:3d}  "
                f"rate={vals[idx]:.8e}  "
                f"{reaction_name(k_index, Rf)}"
            )

            shown += 1
            if shown == 10:
                break

    return values


def compare(ref_values, base_values, term_info, ref_data, base_data):
    print()
    print("=" * 120)
    print("REFERENCE vs BASELINE — SO2 REACTION-RATE DIFFERENCES")
    print("=" * 120)

    ref_Rf = ref_data["variable"]["Rf"]
    base_Rf = base_data["variable"]["Rf"]

    if ref_Rf != base_Rf:
        print("WARNING: Rf dictionaries are not identical.")

    for level in TARGET_LEVELS:
        print()
        print(f"LEVEL {level}")
        print("-" * 120)

        rows = []

        for i, (k_index, sign, expression) in enumerate(term_info):
            ref = ref_values[i, level]
            base = base_values[i, level]

            diff = base - ref

            scale = max(abs(ref), 1e-300)
            rel = abs(diff) / scale

            rows.append(
                (
                    abs(diff),
                    rel,
                    k_index,
                    ref,
                    base,
                    diff,
                )
            )

        rows.sort(reverse=True)

        print(
            f"{'Reaction':>8}  "
            f"{'Reference':>16}  "
            f"{'Baseline':>16}  "
            f"{'Difference':>16}  "
            f"{'RelDiff':>12}  "
            f"Reaction"
        )
        print("-" * 120)

        for _, rel, k_index, ref, base, diff in rows[:15]:
            odd = k_index if k_index % 2 else k_index - 1
            name = ref_Rf.get(odd, "<unknown>")

            if k_index % 2 == 0:
                name = "reverse: " + name

            print(
                f"R{k_index:3d}  "
                f"{ref:16.8e}  "
                f"{base:16.8e}  "
                f"{diff:16.8e}  "
                f"{rel:12.4e}  "
                f"{name}"
            )

    print()
    print("=" * 120)
    print("SO2 NET CHEMICAL TENDENCY — REFERENCE vs BASELINE")
    print("=" * 120)

    for level in TARGET_LEVELS:
        ref_net = np.sum(ref_values[:, level])
        base_net = np.sum(base_values[:, level])

        print(
            f"Level {level:2d} : "
            f"reference={ref_net:.8e}   "
            f"baseline={base_net:.8e}   "
            f"difference={base_net-ref_net:.8e}"
        )


def main():
    print("=" * 120)
    print("STAGE 7C.1 — SO2 CHEMICAL SOURCE/SINK TRACE")
    print("=" * 120)

    print()
    print("Reference:", REF)
    print("Baseline :", BASE)
    print("Target levels:", TARGET_LEVELS)

    if not os.path.exists(REF):
        raise FileNotFoundError(REF)

    if not os.path.exists(BASE):
        raise FileNotFoundError(BASE)

    ref_data = load_vul(REF)
    base_data = load_vul(BASE)

    term_info = extract_so2_rate_terms()

    print()
    print(f"SO2 individual rate terms found: {len(term_info)}")

    ref_values = summarize(
        "REFERENCE",
        ref_data,
        term_info
    )

    base_values = summarize(
        "BASELINE",
        base_data,
        term_info
    )

    compare(
        ref_values,
        base_values,
        term_info,
        ref_data,
        base_data
    )

    print()
    print("=" * 120)
    print("STAGE 7C.1 COMPLETE")
    print("=" * 120)
    print("No VULCAN integration was performed.")
    print("No model/configuration parameters were changed.")
    print("No Unknown UV absorber was added or activated.")
    print("=" * 120)


if __name__ == "__main__":
    main()
