"""Prototype: build one placeholder X(TAMU) reaction item, write it to XML,
validate it against the real xmlcoll schema, and round-trip it back to
plain Python data to confirm nothing is lost.

The rate numbers below are FABRICATED, for schema testing only -- not a
real measured 12C(p,g)13N rate.
"""

import os
import sys

from lxml import etree
import xmlcoll.coll as xc

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"),
)

from unicorrn.reaction import (
    FORMAT_NAME,
    FORMAT_VERSION,
    build_reaction_item,
    read_reaction_item,
)

OUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_x_tamu.xml")

RATE_TABLE = [
    # t9,    rate,     sigma_low, sigma_high
    (0.010, 1.20e-22, 1.00e-22, 1.50e-22),
    (0.100, 3.40e-11, 3.00e-11, 3.90e-11),
    (0.500, 2.10e-05, 1.90e-05, 2.30e-05),
    (1.000, 4.50e-03, 4.10e-03, 4.90e-03),
    (2.000, 1.80e-01, 1.60e-01, 2.00e-01),
    (5.000, 3.30e00, 3.00e00, 3.60e00),
]


def main():
    collection = xc.Collection()
    collection.update_properties(
        {
            "format": FORMAT_NAME,
            "format_version": FORMAT_VERSION,
            "lab": "TAMU",
        }
    )

    item = build_reaction_item(
        name="c12_pg_n13",
        reactants=["c12", "p"],
        products=["n13", "gamma"],
        rate_table=RATE_TABLE,
        mechanisms=[{"type": "charge_induced_nonresonant", "l": 1}],
        q_value_mev=1.9435,
        source="TAMU placeholder (schema test, not real data)",
        evaluator="C. Rathi",
        date="2026-09-03",
        status="placeholder",
    )
    collection.add_item(item)

    collection.write_to_xml(OUT_FILE)
    print(f"wrote {OUT_FILE}")

    # --- Real schema validation (xmlcoll's own validate() computes a
    # pass/fail bool internally but never returns or raises on it, so it
    # can't actually tell you whether the file is valid -- do it directly
    # with lxml instead so a violation is loud.) ---
    parser = etree.XMLParser(remove_blank_text=True)
    doc = etree.parse(OUT_FILE, parser)
    doc.xinclude()
    schema_path = os.path.join(os.path.dirname(xc.__file__), "xsd_pub", "xmlcoll.xsd")
    schema_doc = etree.parse(schema_path)
    validator = etree.XMLSchema(schema_doc)
    if validator.validate(doc):
        print("schema validation: OK (xmlcoll.xsd)")
    else:
        print("schema validation: FAILED")
        for error in validator.error_log:
            print(" ", error)

    # --- Round trip: read the XML back and confirm values match exactly ---
    check = xc.Collection()
    check.update_from_xml(OUT_FILE)
    restored = read_reaction_item(check.get()["c12_pg_n13"])

    assert restored["reactants"] == ["c12", "p"]
    assert restored["products"] == ["n13", "gamma"]
    assert restored["rate_table"] == RATE_TABLE
    assert restored["mechanisms"] == [{"type": "charge_induced_nonresonant", "l": 1}]
    print("round-trip check: OK, values match exactly")

    print("\nRestored data:")
    for key, value in restored.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
