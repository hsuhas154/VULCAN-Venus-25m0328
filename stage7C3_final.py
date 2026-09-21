import pickle
import numpy as np
import chem_funs

REF = "output/Nominal_Bkzz_SO2-25m0328-reference.vul"
BASE = "output/Nominal_Bkzz_SO2-25m0328-baseline.vul"

LEVELS = [27, 28, 29, 30, 31]

SPEC = {sp: i for i, sp in enumerate(chem_funs.spec_list)}


def load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def state(data, sp, level):
    y = np.asarray(data["variable"]["y"])
    return float(y[level, SPEC[sp]])


def k(data, reaction, level):
    arr = np.asarray(data["variable"]["k"][reaction])

    if arr.ndim == 0:
        return float(arr)

    if arr.ndim == 1:
        if len(arr) == 1:
            return float(arr[0])
        return float(arr[level])

    raise RuntimeError(f"Unexpected k[{reaction}] shape {arr.shape}")


def M(data, level):
    return float(np.asarray(data["atm"]["M"])[level])


def components(data, fwd, rev, level):
    """
    Exact forward and reverse components corresponding to the
    generated VULCAN v_* definitions in chem_funs.py.
    """

    H2SO4 = state(data, "H2SO4", level)
    H2SO4_l = state(data, "H2SO4_l", level)
    SO3 = state(data, "SO3", level)
    H2O = state(data, "H2O", level)
    SO2 = state(data, "SO2", level)
    SO2_l = state(data, "SO2_l", level)

    if fwd == 169:
        forward = k(data,169,level) * H2SO4 * H2O
        reverse = k(data,170,level) * SO3 * H2O * H2O

    elif fwd == 723:
        forward = k(data,723,level) * SO3 * H2O * H2O
        reverse = k(data,724,level) * H2SO4 * H2O

    elif fwd == 727:
        forward = k(data,727,level) * H2SO4
        reverse = k(data,728,level) * H2SO4_l

    elif fwd == 791:
        forward = k(data,791,level) * H2SO4
        reverse = k(data,792,level) * SO3 * H2O

    elif fwd == 589:
        forward = k(data,589,level) * SO2_l * H2SO4
        reverse = k(data,590,level) * SO2 * H2SO4

    else:
        raise ValueError(fwd)

    return forward, reverse, forward - reverse


def compare(ref, base, fwd, rev, label):
    print()
    print("=" * 120)
    print(f"R{fwd}/R{rev}: {label}")
    print("=" * 120)

    for level in LEVELS:
        rf, rr, rn = components(ref, fwd, rev, level)
        bf, br, bn = components(base, fwd, rev, level)

        print(
            f"Level {level:2d} | "
            f"F ref={rf:.8e} base={bf:.8e} "
            f"| R ref={rr:.8e} base={br:.8e} "
            f"| NET ref={rn:.8e} base={bn:.8e}"
        )


def main():
    ref = load(REF)
    base = load(BASE)

    print("=" * 120)
    print("STAGE 7C.3 FINAL — H2SO4 SOURCE/SINK DECOMPOSITION")
    print("=" * 120)

    print("\nH2SO4 state comparison")
    print("-" * 120)

    for level in LEVELS:
        rh = state(ref, "H2SO4", level)
        bh = state(base, "H2SO4", level)

        rhl = state(ref, "H2SO4_l", level)
        bhl = state(base, "H2SO4_l", level)

        rs = state(ref, "SO3", level)
        bs = state(base, "SO3", level)

        print(
            f"Level {level:2d} | "
            f"H2SO4 ref={rh:.8e} base={bh:.8e} "
            f"| H2SO4_l ref={rhl:.8e} base={bhl:.8e} "
            f"| SO3 ref={rs:.8e} base={bs:.8e}"
        )

    compare(
        ref, base, 169, 170,
        "H2SO4 + H2O <-> SO3 + H2O + H2O"
    )

    compare(
        ref, base, 723, 724,
        "SO3 + H2O + H2O <-> H2SO4 + H2O"
    )

    compare(
        ref, base, 727, 728,
        "H2SO4 <-> H2SO4_l"
    )

    compare(
        ref, base, 791, 792,
        "H2SO4 <-> SO3 + H2O"
    )

    compare(
        ref, base, 589, 590,
        "SO2_l + H2SO4 <-> SO2 + H2SO4"
    )

    print()
    print("=" * 120)
    print("STAGE 7C.3 COMPLETE")
    print("=" * 120)
    print("No VULCAN integration performed.")
    print("No configuration changed.")
    print("No Unknown UV absorber activated.")
    print("=" * 120)


if __name__ == "__main__":
    main()
