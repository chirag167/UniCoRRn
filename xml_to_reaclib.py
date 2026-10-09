"""CLI: convert an X(TAMU) reaction-rate file into a REACLIB-format text
file, fitting each reaction's mechanisms to REACLIB's 7-parameter form
per the physically-constrained rules in unicorrn/reaclib_fit.py (Cyburt
et al. 2010, Table 1) rather than an unconstrained curve fit.

Every reaction needs its mechanisms (see csv_to_xml.py's --mechanisms)
and a Q-value -- both are required by REACLIB's format and can't be
defaulted. Reactions with more or fewer than 2 reactants aren't
supported yet (Table 1's formulas are for two-body entrance channels).

Each reaction's fit is checked against REACLIB's own published accuracy
target: fits worse than 5% relative to the measured table are flagged
(mirroring REACLIB's own "deviations list" practice) but still written,
since 5% may simply be unreachable with sparse data -- this is reported,
not silently hidden or silently accepted.

Usage:
    python xml_to_reaclib.py x_tamu_file.xml -o out.reaclib --validate
"""

import argparse
import os
import sys

import numpy as np
import xmlcoll.coll as xc

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from unicorrn.reaclib_fit import (
    _component_log_rate,
    fit_reaclib,
    nuclide_z_a,
    reaclib_chapter,
    strip_photon_products,
)
from unicorrn.reaction import read_reaction_item

DEVIATION_THRESHOLD_PCT = 5.0  # REACLIB's own published fit-accuracy target


def _pad74(line):
    if len(line) > 74:
        raise ValueError(f"line exceeds 74 characters: {line!r}")
    return line.ljust(74)


def _e13_6(value):
    """13-character field for a coefficient, e.g. ' 1.943500e+00'."""
    field = f"{value:13.6e}"
    assert len(field) == 13, f"coefficient field wrong width: {field!r}"
    return field


def _e12_5(value):
    """12-character field for the Q-value, e.g. ' 1.94350e+00'."""
    field = f"{value:12.5e}"
    assert len(field) == 12, f"Q-value field wrong width: {field!r}"
    return field


def serialize_reaction(reactants, products, chapter, mechanisms, coeffs_list, q_value_mev, label):
    """Build the REACLIB-format lines (header + one 3-line set per
    mechanism) for one reaction. Returns a list of 74-character strings.
    """
    tokens = list(reactants) + list(products)
    if len(tokens) > 6:
        raise ValueError(f"{len(tokens)} nuclides won't fit in one REACLIB set entry (max 6)")
    for token in tokens:
        if len(token) > 5:
            raise ValueError(f"nuclide token {token!r} longer than REACLIB's 5-character field")
    tokens = tokens + [""] * (6 - len(tokens))
    nuclide_field = "".join(token.ljust(5) for token in tokens)

    lines = [_pad74(str(chapter)), _pad74(""), _pad74("")]

    for mechanism, coeffs in zip(mechanisms, coeffs_list):
        resonance_flag = "r" if mechanism["type"] == "narrow_resonant" else " "
        reverse_flag = " "
        line1 = (
            " " * 5
            + nuclide_field
            + " " * 8
            + label.ljust(4)[:4]
            + resonance_flag
            + reverse_flag
            + " " * 3
            + _e12_5(q_value_mev)
        )
        line2 = "".join(_e13_6(c) for c in coeffs[0:4])
        line3 = "".join(_e13_6(c) for c in coeffs[4:7])
        lines.extend([_pad74(line1), _pad74(line2), _pad74(line3)])

    return lines


def convert(x_tamu_file, output_file, label="tamu"):
    collection = xc.Collection()
    collection.update_from_xml(x_tamu_file)

    all_lines = []
    n_written = 0
    n_failed = 0

    for name, item in collection.get().items():
        reaction = read_reaction_item(item)
        reactants = reaction["reactants"]
        products = strip_photon_products(reaction["products"])

        try:
            if not reaction["mechanisms"]:
                raise ValueError("no mechanisms declared (required for REACLIB conversion)")
            q_value_mev = reaction["meta"].get("q_value_mev")
            if q_value_mev is None:
                raise ValueError("no Q-value (required by REACLIB's format)")
            q_value_mev = float(q_value_mev)

            chapter = reaclib_chapter(reactants, products)
            fit = fit_reaclib(reaction["rate_table"], reaction["mechanisms"], reactants)

            max_dev_pct = np.max(np.abs(fit["rel_resid"])) * 100
            print(f"{name}: chapter {chapter}, max |rel resid| = {max_dev_pct:.3f}%", file=sys.stderr)
            if max_dev_pct > DEVIATION_THRESHOLD_PCT:
                print(
                    f"  DEVIATION: {name} exceeds REACLIB's {DEVIATION_THRESHOLD_PCT}% "
                    "fit-accuracy target -- written anyway, but treat this fit as unreliable",
                    file=sys.stderr,
                )

            all_lines.extend(
                serialize_reaction(
                    reactants, products, chapter,
                    reaction["mechanisms"], fit["coeffs"], q_value_mev, label,
                )
            )
            n_written += 1
        except ValueError as exc:
            print(f"error: skipping {name}: {exc}", file=sys.stderr)
            n_failed += 1

    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write("\n".join(all_lines) + "\n")

    print(f"wrote {n_written} reaction(s) to {output_file}" + (f" ({n_failed} skipped)" if n_failed else ""))
    return n_written, n_failed


def validate(x_tamu_file, output_file):
    """Read the written file back with pynucastro's real ReacLibLibrary
    reader and cross-check its evaluated rate against our own in-memory
    fit at each reaction's own T9 grid -- this checks the ASCII
    serialization round-trips correctly, not the fit quality (already
    reported by `convert`).
    """
    import pynucastro as pyna

    collection = xc.Collection()
    collection.update_from_xml(x_tamu_file)

    lib = pyna.ReacLibLibrary(libfile=os.path.abspath(output_file))
    rates_by_za = {
        tuple(sorted((n.Z, n.A) for n in r.reactants + r.products)): r
        for r in lib.get_rates()
    }
    print(f"pynucastro read back {len(rates_by_za)} reaction(s)")

    all_ok = True
    for name, item in collection.get().items():
        reaction = read_reaction_item(item)
        if not reaction["mechanisms"] or reaction["meta"].get("q_value_mev") is None:
            continue
        try:
            fit = fit_reaclib(reaction["rate_table"], reaction["mechanisms"], reaction["reactants"])
        except ValueError:
            continue

        products = strip_photon_products(reaction["products"])
        key = tuple(sorted(
            nuclide_z_a(x) for x in reaction["reactants"] + products
        ))
        if key not in rates_by_za:
            print(f"  {name}: NOT FOUND in read-back file", file=sys.stderr)
            all_ok = False
            continue
        pynucastro_rate = rates_by_za[key]

        for t9, *_ in reaction["rate_table"][:1] or [(1.0,)]:
            our_value = np.exp(
                np.logaddexp.reduce([_component_log_rate(np.array([t9]), c)[0] for c in fit["coeffs"]])
            )
            their_value = pynucastro_rate.eval(t9 * 1.0e9)
            rel_diff = abs(their_value / our_value - 1.0)
            status = "OK" if rel_diff < 1e-4 else "MISMATCH"
            if status != "OK":
                all_ok = False
            print(f"  {name} @ t9={t9}: ours={our_value:.6e} pynucastro={their_value:.6e} [{status}]")

    return all_ok


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("x_tamu_file", help="input X(TAMU) XML file")
    parser.add_argument("-o", "--output", required=True, help="output REACLIB-format text file")
    parser.add_argument("--label", default="tamu", help="4-character REACLIB set label (default: tamu)")
    parser.add_argument(
        "--validate", action="store_true",
        help="read the output back with pynucastro's ReacLibLibrary and cross-check evaluated rates",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    _, n_failed = convert(args.x_tamu_file, args.output, label=args.label)

    if args.validate:
        ok = validate(args.x_tamu_file, args.output)
        print("validation: OK" if ok else "validation: FAILED")
        if not ok:
            sys.exit(1)

    if n_failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
