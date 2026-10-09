"""Fit an X rate(T9) table to REACLIB's 7-parameter form, following the
physically-constrained fitting rules in the JINA REACLIB database paper
(Cyburt et al. 2010, ApJS 189, 240, Table 1). Most of a0..a6 are FIXED by
the reaction's known physics (charges, masses, orbital angular momentum,
resonance energy) rather than freely fit -- this is what keeps the result
from diverging outside the tabulated range the way an unconstrained
7-parameter fit can (see the exploratory fit that overflowed at t9=100).

Supported mechanism types (see unicorrn.reaction.MECHANISM_TYPES):
    charge_induced_nonresonant: a1=0, a2=fixed(Z1,Z2,A_reduced), a6=-2/3.
        a0 fixed from s0 if given, else floated. a3,a4,a5 always floated.
    neutron_induced_nonresonant: a1=0, a2=0, a6=l.
        a0 floated (our schema has no field for the theoretical
        normalization Table 1's formula needs). a3,a4,a5 always floated.
    narrow_resonant: a1=fixed(-11.6045*er), a2=a3=a4=a5=0, a6=-3/2.
        a0 floated (it encodes the resonance strength omega*gamma, which
        our schema doesn't take as an input either).

Only two-reactant reactions (REACLIB chapters 4, 5, 6, 7) are supported --
Table 1's formulas describe entrance-channel physics between two
particles and don't depend on the exit channel, but a fundamentally
different entrance channel (one reactant, or three) needs different
formulas this module doesn't implement.

A reaction's mechanisms are summed, not switched between:
    rate_total(T9) = sum_i exp(set_i(T9))
Fitting this jointly is genuinely nonlinear once there's more than one
component (log of a sum of exponentials isn't a sum), so this uses
scipy.optimize.least_squares rather than plain linear least-squares.
"""

import re

import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp

# --- physical constants (Cyburt et al. 2010, Table 1) ---
_B = 7.8318e9          # cm^3 s^-1 mol^-1 MeV^-1 barn^-1 (q-induced NR a0)
_D = 1.5394e11         # cm^3 s^-1 mol^-1 MeV^-1 (narrow resonant a0)

_ELEMENTS = [
    "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne",
    "na", "mg", "al", "si", "p", "s", "cl", "ar", "k", "ca",
    "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn",
    "ga", "ge", "as", "se", "br", "kr", "rb", "sr", "y", "zr",
    "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn",
    "sb", "te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd",
    "pm", "sm", "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb",
    "lu", "hf", "ta", "w", "re", "os", "ir", "pt", "au", "hg",
    "tl", "pb", "bi", "po", "at", "rn", "fr", "ra", "ac", "th",
    "pa", "u", "np", "pu", "am", "cm", "bk", "cf", "es", "fm",
    "md", "no", "lr", "rf", "db", "sg", "bh", "hs", "mt", "ds",
    "rg", "cn", "nh", "fl", "mc", "lv", "ts", "og",
]
_Z_BY_SYMBOL = {symbol: z for z, symbol in enumerate(_ELEMENTS, start=1)}

_SHORTHAND_ZA = {
    "n": (0, 1),
    "p": (1, 1),
    "d": (1, 2),
    "t": (1, 3),
    "a": (2, 4),
    "alpha": (2, 4),
}

_NUCLIDE_RE = re.compile(r"^([a-z]{1,2})(\d+)$")


def nuclide_z_a(name):
    """Parse a reactant/product name into (Z, A).

    Recognizes shorthand (p, d, t, a/alpha, n) and standard
    element-symbol + mass-number names (e.g. c12, ne20).
    """
    if name in _SHORTHAND_ZA:
        return _SHORTHAND_ZA[name]

    match = _NUCLIDE_RE.match(name)
    if not match:
        raise ValueError(f"cannot parse nuclide name {name!r}")
    symbol, mass_str = match.groups()
    if symbol not in _Z_BY_SYMBOL:
        raise ValueError(f"unknown element symbol {symbol!r} in nuclide name {name!r}")
    return _Z_BY_SYMBOL[symbol], int(mass_str)


PHOTON_TOKENS = frozenset({"gamma", "g", "photon"})


def strip_photon_products(products):
    """REACLIB doesn't track an emitted photon as a nuclide token at all
    (pynucastro's Nucleus class has no special case for "gamma"/"g" --
    unlike webnucleo, which does). c12(p,g)n13 is chapter 4 with tokens
    [c12, p, n13], not chapter 5 with an explicit gamma product.
    """
    return [p for p in products if p not in PHOTON_TOKENS]


CHAPTER_BY_COUNTS = {
    (1, 1): 1, (1, 2): 2, (1, 3): 3,
    (2, 1): 4, (2, 2): 5, (2, 3): 6, (2, 4): 7,
    (3, 1): 8, (3, 2): 9, (4, 2): 10, (1, 4): 11,
}


def reaclib_chapter(reactants, products):
    key = (len(reactants), len(products))
    if key not in CHAPTER_BY_COUNTS:
        raise ValueError(
            f"no REACLIB chapter for {len(reactants)} reactants -> "
            f"{len(products)} products"
        )
    return CHAPTER_BY_COUNTS[key]


def _mechanism_fixed_and_float(mechanism, reactants):
    """Return (fixed: {coeff_index: value}, float_idx: [coeff_index, ...])
    for one mechanism component, per Table 1 (Cyburt et al. 2010).
    """
    mech_type = mechanism["type"]

    if mech_type in ("charge_induced_nonresonant", "neutron_induced_nonresonant"):
        if len(reactants) != 2:
            raise ValueError(
                f"{mech_type} requires exactly 2 reactants, got {reactants}"
            )
        (z1, a1), (z2, a2) = nuclide_z_a(reactants[0]), nuclide_z_a(reactants[1])
        a_reduced = a1 * a2 / (a1 + a2)

        if mech_type == "charge_induced_nonresonant":
            fixed = {1: 0.0, 2: -4.2486 * (z1 * z1 * z2 * z2 * a_reduced) ** (1.0 / 3.0), 6: -2.0 / 3.0}
            float_idx = [3, 4, 5]
            if "s0" in mechanism:
                fixed[0] = np.log(_B * (z1 * z2 / a_reduced) ** (1.0 / 3.0) * mechanism["s0"])
            else:
                float_idx = [0, 3, 4, 5]
            return fixed, float_idx

        # neutron_induced_nonresonant
        ell = float(mechanism.get("l", 0))
        fixed = {1: 0.0, 2: 0.0, 6: ell}
        float_idx = [0, 3, 4, 5]
        return fixed, float_idx

    if mech_type == "narrow_resonant":
        er = mechanism["er"]
        fixed = {1: -11.6045 * er, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: -1.5}
        float_idx = [0]
        return fixed, float_idx

    raise ValueError(f"unsupported mechanism type: {mech_type!r}")


def _component_log_rate(t9, coeffs):
    return (
        coeffs[0]
        + coeffs[1] / t9
        + coeffs[2] * t9 ** (-1.0 / 3.0)
        + coeffs[3] * t9 ** (1.0 / 3.0)
        + coeffs[4] * t9
        + coeffs[5] * t9 ** (5.0 / 3.0)
        + coeffs[6] * np.log(t9)
    )


def fit_reaclib(rate_table, mechanisms, reactants):
    """Jointly fit all mechanism components against one measured total
    rate(T9) table.

    Args:
        rate_table: list of (t9, rate, sigma_low, sigma_high) rows (sigma
            entries are ignored -- REACLIB has no field for them).
        mechanisms: list of mechanism dicts (see module docstring).
        reactants: reactant name list, needed for Z/A in the fixed terms.

    Returns:
        dict with:
            "coeffs": list of length-7 numpy arrays, one per mechanism,
                in the same order as `mechanisms`.
            "rel_resid": array of (model/measured - 1) at each table row.
            "result": the scipy.optimize.OptimizeResult.
    """
    if len(rate_table) == 0:
        raise ValueError("rate_table is empty")

    t9 = np.array([row[0] for row in rate_table])
    rate = np.array([row[1] for row in rate_table])
    log_rate = np.log(rate)

    specs = [_mechanism_fixed_and_float(m, reactants) for m in mechanisms]

    free_params = [
        (ci, k) for ci, (_, float_idx) in enumerate(specs) for k in float_idx
    ]
    if not free_params:
        raise ValueError("no floated parameters -- nothing to fit")

    n_components = len(mechanisms)

    def build_coeffs(x):
        coeffs = [np.zeros(7) for _ in range(n_components)]
        for ci, (fixed, _) in enumerate(specs):
            for k, v in fixed.items():
                coeffs[ci][k] = v
        for (ci, k), val in zip(free_params, x):
            coeffs[ci][k] = val
        return coeffs

    def total_log_rate(t9_vals, coeffs):
        logs = np.stack([_component_log_rate(t9_vals, c) for c in coeffs], axis=0)
        return logsumexp(logs, axis=0)

    def residuals(x):
        return total_log_rate(t9, build_coeffs(x)) - log_rate

    x0 = np.array([
        np.mean(log_rate) / n_components if k == 0 else 0.0
        for _, k in free_params
    ])

    result = least_squares(residuals, x0)
    fitted_coeffs = build_coeffs(result.x)
    log_resid = residuals(result.x)
    rel_resid = np.expm1(log_resid)  # log(model)-log(meas) = log(model/meas)

    return {"coeffs": fitted_coeffs, "rel_resid": rel_resid, "result": result}
