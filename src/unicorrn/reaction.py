"""Build and parse X reaction-rate items on top of xmlcoll.

X is the generic reaction-rate format; a given file is one lab's or one
evaluation's instance of it (e.g. X(TAMU)), identified by its own `lab`
property rather than by this package's name.

One reaction (one lab, one version) is one xmlcoll Item. Because xmlcoll
properties are scalars keyed by (name, tag1..tag5), reactant/product lists
and the rate(t9) table are flattened into indexed/tagged properties rather
than stored as native XML lists:

    ('reactant', '1')            -> 'c12'
    ('reactant', '2')            -> 'p'
    ('product', '1')             -> 'n13'
    ('product', '2')             -> 'gamma'
    ('rate', 't9', '000')        -> 0.01
    ('rate', 'value', '000')     -> 1.2e-22
    ('rate', 'sigma_low', '000') -> 1.0e-22
    ('rate', 'sigma_high', '000')-> 1.5e-22
    ('mechanism', 'type', '000') -> 'narrow_resonant'
    ('mechanism', 'er', '000')   -> 0.183

Row index is zero-padded (3 digits) so table rows/mechanisms sort correctly
as strings.

A reaction's true rate can be the sum of several physically distinct
mechanisms (e.g. a non-resonant background plus one or more narrow
resonances) contributing simultaneously at every T9 -- not different
mechanisms for different T9 ranges. `mechanisms` is therefore a list, one
entry per physical component, each requiring different fields:

    charge_induced_nonresonant: l (optional, default 0), s0 (optional)
    neutron_induced_nonresonant: l (optional, default 0)
    narrow_resonant: er (required, MeV)
"""

import xmlcoll.coll as xc

FORMAT_NAME = "X"
FORMAT_VERSION = "0.1"

_ROW_WIDTH = 3

MECHANISM_TYPES = frozenset({
    "charge_induced_nonresonant",
    "neutron_induced_nonresonant",
    "narrow_resonant",
})

_MECHANISM_FIELD_CASTS = {
    "type": str,
    "l": int,
    "er": float,
    "s0": float,
}


def build_reaction_item(
    name,
    reactants,
    products,
    rate_table,
    mechanisms,
    q_value_mev=None,
    source=None,
    evaluator=None,
    date=None,
    status=None,
):
    """Build an xmlcoll Item representing one reaction rate.

    Args:
        name: unique item name, e.g. "c12_pg_n13".
        reactants: ordered list of reactant names, e.g. ["c12", "p"].
        products: ordered list of product names, e.g. ["n13", "gamma"].
        rate_table: list of (t9, rate, sigma_low, sigma_high) rows.
            sigma_low/sigma_high may be None if not available.
        mechanisms: list of dicts, one per physical component contributing
            to the total rate, each with a "type" key (one of
            MECHANISM_TYPES) plus that type's fields (see module docstring).
        q_value_mev, source, evaluator, date, status: optional metadata.

    Returns:
        xmlcoll.coll.Item
    """
    item = xc.Item(name)
    props = {}

    for i, reactant in enumerate(reactants, start=1):
        props[("reactant", str(i))] = reactant
    for i, product in enumerate(products, start=1):
        props[("product", str(i))] = product

    if q_value_mev is not None:
        props["q_value_mev"] = q_value_mev
    if source is not None:
        props["source"] = source
    if evaluator is not None:
        props["evaluator"] = evaluator
    if date is not None:
        props["date"] = date
    if status is not None:
        props["status"] = status

    for i, (t9, rate, sigma_low, sigma_high) in enumerate(rate_table):
        idx = str(i).zfill(_ROW_WIDTH)
        props[("rate", "t9", idx)] = t9
        props[("rate", "value", idx)] = rate
        if sigma_low is not None:
            props[("rate", "sigma_low", idx)] = sigma_low
        if sigma_high is not None:
            props[("rate", "sigma_high", idx)] = sigma_high

    for i, mechanism in enumerate(mechanisms):
        mech_type = mechanism.get("type")
        if mech_type not in MECHANISM_TYPES:
            raise ValueError(
                f"mechanism {i}: type must be one of {sorted(MECHANISM_TYPES)}, "
                f"got {mech_type!r}"
            )
        if mech_type == "narrow_resonant" and "er" not in mechanism:
            raise ValueError(
                f"mechanism {i}: narrow_resonant requires 'er' (resonance "
                "energy in MeV) -- it cannot be defaulted or derived from "
                "the rate table"
            )
        idx = str(i).zfill(_ROW_WIDTH)
        for field, value in mechanism.items():
            if field not in _MECHANISM_FIELD_CASTS:
                raise ValueError(f"mechanism {i}: unknown field {field!r}")
            props[("mechanism", field, idx)] = value

    item.update_properties(props)
    return item


def read_reaction_item(item):
    """Parse an Item built by build_reaction_item back into plain Python data.

    Returns:
        dict with keys: name, reactants, products, rate_table, mechanisms,
        meta. rate_table is a list of (t9, rate, sigma_low, sigma_high)
        float tuples, sorted by row index; sigma entries are None if
        absent. mechanisms is a list of dicts, one per physical component,
        sorted in the order they were written.

    Raises:
        ValueError: if a property has a tuple key with an unrecognized
        prefix -- this fails loudly rather than silently dropping data
        a future property type forgot to teach this function about.
    """
    props = item.get_properties()

    reactants, products, meta, rows, mechs = {}, {}, {}, {}, {}

    for key, value in props.items():
        if isinstance(key, tuple) and key[0] == "reactant":
            reactants[int(key[1])] = value
        elif isinstance(key, tuple) and key[0] == "product":
            products[int(key[1])] = value
        elif isinstance(key, tuple) and key[0] == "rate":
            field, idx = key[1], int(key[2])
            rows.setdefault(idx, {})[field] = value
        elif isinstance(key, tuple) and key[0] == "mechanism":
            field, idx = key[1], int(key[2])
            mechs.setdefault(idx, {})[field] = value
        elif not isinstance(key, tuple):
            meta[key] = value
        else:
            raise ValueError(f"read_reaction_item: unrecognized property key {key!r}")

    rate_table = [
        (
            float(rows[idx]["t9"]),
            float(rows[idx]["value"]),
            float(rows[idx]["sigma_low"]) if "sigma_low" in rows[idx] else None,
            float(rows[idx]["sigma_high"]) if "sigma_high" in rows[idx] else None,
        )
        for idx in sorted(rows)
    ]

    mechanisms = [
        {
            field: _MECHANISM_FIELD_CASTS[field](value)
            for field, value in mechs[idx].items()
        }
        for idx in sorted(mechs)
    ]

    return {
        "name": item.get_name(),
        "reactants": [reactants[i] for i in sorted(reactants)],
        "products": [products[i] for i in sorted(products)],
        "rate_table": rate_table,
        "mechanisms": mechanisms,
        "meta": meta,
    }
