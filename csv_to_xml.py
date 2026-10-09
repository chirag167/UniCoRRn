"""CLI: build one X(TAMU) reaction item from a delimited rate(T9) table and
write it into an X(TAMU) XML file (or add it alongside existing reactions
with --append).

Input rate-table format: a header row naming its columns, then one row per
T9 point. Required columns: t9 and one of value/rate. Optional: sigma_low,
sigma_high. Blank lines and lines starting with '#' are ignored.

    t9,value,sigma_low,sigma_high
    0.010,1.2e-22,1.0e-22,1.5e-22
    0.100,3.4e-11,3.0e-11,3.9e-11

Every reaction also needs a --mechanisms JSON file: a non-empty list of the
physical component(s) that sum to the total rate (a rate can be a
non-resonant background plus one or more narrow resonances contributing
simultaneously, not different mechanisms for different T9 ranges -- see
src/unicorrn/reaction.py's module docstring for the field schema per type).

    [{"type": "narrow_resonant", "er": 0.183},
     {"type": "charge_induced_nonresonant", "l": 0}]

Usage:
    python csv_to_xml.py rate.csv c12_pg_n13 \\
        --reactants c12,p --products n13,gamma \\
        --mechanisms mechanisms.json \\
        --q 1.9435 --source "TAMU 2026 run" --evaluator "C. Rathi" \\
        --date 2026-09-03 --status recommended \\
        -o out.xml

    # add a second reaction to the same file:
    python csv_to_xml.py rate2.csv o16_ag_ne20 \\
        --reactants o16,he4 --products ne20,gamma \\
        --mechanisms mechanisms2.json \\
        -o out.xml --append
"""

import argparse
import csv
import json
import os
import sys

from lxml import etree
import xmlcoll.coll as xc

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from unicorrn.reaction import (
    FORMAT_NAME,
    FORMAT_VERSION,
    MECHANISM_TYPES,
    build_reaction_item,
)


def read_mechanisms(path):
    """Read and validate a mechanisms JSON file: a non-empty list of dicts,
    one per physical component contributing to the reaction's total rate
    (see unicorrn/reaction.py's module docstring for field requirements
    per mechanism type).
    """
    with open(path, encoding="utf-8") as handle:
        mechanisms = json.load(handle)

    if not isinstance(mechanisms, list) or not mechanisms:
        raise ValueError(f"{path}: must contain a non-empty JSON list of mechanism objects")

    for i, mechanism in enumerate(mechanisms):
        mech_type = mechanism.get("type")
        if mech_type not in MECHANISM_TYPES:
            raise ValueError(
                f"{path}: mechanism {i}: type must be one of "
                f"{sorted(MECHANISM_TYPES)}, got {mech_type!r}"
            )
        if mech_type == "narrow_resonant" and "er" not in mechanism:
            raise ValueError(
                f"{path}: mechanism {i}: narrow_resonant requires 'er' "
                "(resonance energy in MeV) -- it has no default"
            )
        if mech_type in ("charge_induced_nonresonant", "neutron_induced_nonresonant"):
            mechanism.setdefault("l", 0)

    return mechanisms


def read_rate_table(path, delimiter=","):
    """Read a delimited t9/value[/sigma_low/sigma_high] table.

    Returns:
        list of (t9, rate, sigma_low, sigma_high) float tuples; sigma
        entries are None where the column is absent or the cell is blank.
    """
    with open(path, newline="", encoding="utf-8") as handle:
        lines = (
            line
            for line in handle
            if line.strip() and not line.lstrip().startswith("#")
        )
        reader = csv.reader(lines, delimiter=delimiter, skipinitialspace=True)

        header = [h.strip().lower() for h in next(reader)]
        value_col = next((c for c in ("value", "rate") if c in header), None)
        if "t9" not in header or value_col is None:
            raise ValueError(
                f"{path}: header must contain 't9' and one of 'value'/'rate' "
                f"(got {header})"
            )
        col = {name: i for i, name in enumerate(header)}

        def cell(row, name):
            i = col.get(name)
            if i is None or i >= len(row) or not row[i].strip():
                return None
            return float(row[i])

        rows = []
        for row in reader:
            if not row:
                continue
            t9 = cell(row, "t9")
            rate = cell(row, value_col)
            if t9 is None or rate is None:
                raise ValueError(f"{path}: row missing t9 or {value_col}: {row}")
            rows.append((t9, rate, cell(row, "sigma_low"), cell(row, "sigma_high")))

    if not rows:
        raise ValueError(f"{path}: no data rows found")
    return rows


def validate_xml(path):
    """Validate an X(TAMU) file against xmlcoll's real XSD (requires
    network access to resolve the schema's imported types namespace).
    """
    parser = etree.XMLParser(remove_blank_text=True)
    doc = etree.parse(path, parser)
    doc.xinclude()
    schema_path = os.path.join(os.path.dirname(xc.__file__), "xsd_pub", "xmlcoll.xsd")
    validator = etree.XMLSchema(etree.parse(schema_path))
    if validator.validate(doc):
        return True
    for error in validator.error_log:
        print(" ", error, file=sys.stderr)
    return False


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "rate_file", help="delimited rate table: header + t9,value[,sigma_low,sigma_high] rows"
    )
    parser.add_argument("name", help="unique item name, e.g. c12_pg_n13")
    parser.add_argument("--reactants", required=True, help="comma-separated, e.g. c12,p")
    parser.add_argument("--products", required=True, help="comma-separated, e.g. n13,gamma")
    parser.add_argument(
        "--mechanisms", required=True,
        help="JSON file listing the physical component(s) making up this reaction's "
        "total rate (charge_induced_nonresonant / neutron_induced_nonresonant / "
        "narrow_resonant); see README for the field schema",
    )
    parser.add_argument("--q", type=float, default=None, dest="q_value_mev", help="Q-value in MeV")
    parser.add_argument("--source", default=None)
    parser.add_argument("--evaluator", default=None)
    parser.add_argument("--date", default=None)
    parser.add_argument("--status", default=None, help="e.g. recommended, upper_limit, placeholder")
    parser.add_argument("--delimiter", default=",", help="rate file column delimiter (default: ',')")
    parser.add_argument("--lab", default="TAMU", help="lab tag written into new collection metadata")
    parser.add_argument("-o", "--output", required=True, help="output X(TAMU) XML file")
    parser.add_argument(
        "--append", action="store_true",
        help="add this reaction to the existing collection at --output instead of overwriting it",
    )
    parser.add_argument(
        "--validate", action="store_true",
        help="validate the written file against xmlcoll's XSD (needs network access)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    rate_table = read_rate_table(args.rate_file, delimiter=args.delimiter)
    mechanisms = read_mechanisms(args.mechanisms)

    collection = xc.Collection()
    if args.append and os.path.exists(args.output):
        collection.update_from_xml(args.output)
    else:
        collection.update_properties(
            {"format": FORMAT_NAME, "format_version": FORMAT_VERSION, "lab": args.lab}
        )

    if args.name in collection.get():
        sys.exit(f"error: item '{args.name}' already exists in {args.output}")

    item = build_reaction_item(
        name=args.name,
        reactants=[r.strip() for r in args.reactants.split(",")],
        products=[p.strip() for p in args.products.split(",")],
        rate_table=rate_table,
        mechanisms=mechanisms,
        q_value_mev=args.q_value_mev,
        source=args.source,
        evaluator=args.evaluator,
        date=args.date,
        status=args.status,
    )
    collection.add_item(item)
    collection.write_to_xml(args.output)
    print(f"wrote {len(rate_table)}-row rate for '{args.name}' to {args.output}")

    if args.validate:
        if validate_xml(args.output):
            print("schema validation: OK (xmlcoll.xsd)")
        else:
            sys.exit("schema validation: FAILED")


if __name__ == "__main__":
    main()
