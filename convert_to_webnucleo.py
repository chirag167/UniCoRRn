"""Convert an X(TAMU) reaction-rate file into a webnucleo reaction_data XML
file, using wnutils' own Reaction/New_Xml object model and its offline
schema validator -- no hand-written XML on either end.

Two lossy/mapping steps are unavoidable and are handled explicitly rather
than silently:

1. Naming convention. X(TAMU) reactant/product names are whatever the lab
   used (this example project uses shorthand like 'p', 'a', 'gamma').
   Webnucleo/libnucnet expects nuclide names as element symbol + mass
   number (a proton is 'h1', not 'p'; an alpha is 'he4', not 'a') plus a
   fixed set of non-nuclide names (gamma, electron, positron, the
   neutrinos). NAME_ALIASES below translates the common shorthand;
   anything not in it is passed through unchanged (assumed already in
   webnucleo form, e.g. 'c12', 'ne20').

2. Uncertainty. Webnucleo's native rate_table type is only (t9, rate, sef)
   triples -- there is no sigma_low/sigma_high field. If the X(TAMU) item
   carries uncertainties, they are summarized into the human-readable
   <source> text and a warning is printed; they are not machine-readable
   in the converted file. The original X(TAMU) file remains the record of
   record for uncertainties.
"""

import argparse
import os
import sys

import wnutils.xml as wx
import xmlcoll.coll as xc

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from unicorrn.reaction import read_reaction_item

NAME_ALIASES = {
    "p": "h1",
    "d": "h2",
    "t": "h3",
    "a": "he4",
    "alpha": "he4",
    "n": "n1",
    "g": "gamma",
    "e-": "electron",
    "e+": "positron",
}


def to_webnucleo_name(name):
    """Translate common nuclear-astrophysics shorthand to webnucleo's
    nuclide-name convention. Anything not recognized is passed through.
    """
    return NAME_ALIASES.get(name, name)


def build_webnucleo_reaction(reaction_dict):
    """Build a wnutils.xml.Reaction from an X(TAMU) reaction dict, as
    returned by x_tamu.reaction.read_reaction_item.
    """
    reaction = wx.Reaction()
    reaction.reactants = [to_webnucleo_name(r) for r in reaction_dict["reactants"]]
    reaction.products = [to_webnucleo_name(p) for p in reaction_dict["products"]]

    meta = reaction_dict["meta"]
    source_parts = []
    if "source" in meta:
        source_parts.append(meta["source"])
    if "evaluator" in meta:
        source_parts.append(f"evaluator: {meta['evaluator']}")
    if "date" in meta:
        source_parts.append(f"date: {meta['date']}")
    if "q_value_mev" in meta:
        source_parts.append(f"Q = {meta['q_value_mev']} MeV")

    rate_table = reaction_dict["rate_table"]
    has_sigma = any(row[2] is not None or row[3] is not None for row in rate_table)
    if has_sigma:
        source_parts.append(
            "NOTE: X(TAMU) source carries rate uncertainties "
            "(sigma_low/sigma_high) not representable in webnucleo's "
            "rate_table type; see the original X(TAMU) file for those."
        )
        print(
            f"warning: {reaction_dict['name']}: dropping sigma_low/sigma_high "
            "(webnucleo rate_table has no uncertainty field)",
            file=sys.stderr,
        )

    reaction.source = "; ".join(source_parts) if source_parts else "X(TAMU)"

    reaction.data = {
        "type": "rate_table",
        "t9": [row[0] for row in rate_table],
        "rate": [row[1] for row in rate_table],
        "sef": [1.0] * len(rate_table),
    }
    return reaction


def convert(x_tamu_file, output_file, do_validate=False):
    collection = xc.Collection()
    collection.update_from_xml(x_tamu_file)

    reactions = {}
    for name, item in collection.get().items():
        reaction_dict = read_reaction_item(item)
        reactions[name] = build_webnucleo_reaction(reaction_dict)

    if not reactions:
        sys.exit(f"error: no reaction items found in {x_tamu_file}")

    new_xml = wx.New_Xml(xml_type="reaction_data")
    new_xml.set_reaction_data(reactions)
    new_xml.write(output_file)
    print(f"wrote {len(reactions)} reaction(s) to {output_file}")

    if do_validate:
        wx.validate(output_file)
        print("schema validation: OK (libnucnet__reac.xsd)")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("x_tamu_file", help="input X(TAMU) XML file")
    parser.add_argument(
        "-o", "--output", required=True, help="output webnucleo reaction_data XML file"
    )
    parser.add_argument(
        "--validate", action="store_true", dest="do_validate",
        help="validate output against libnucnet's reaction schema (offline, uses wnutils' bundled XSDs)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    convert(args.x_tamu_file, args.output, do_validate=args.do_validate)


if __name__ == "__main__":
    main()
