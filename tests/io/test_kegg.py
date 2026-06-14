"""Test static functions for KEGG."""

import textwrap

from acore.io.kegg import (
    lookup_cid_to_kegg_id,
    parse_compound_pathway_mapping,
    parse_kegg_name_description,
)


def test_lookup_cid_to_kegg_id():
    res = lookup_cid_to_kegg_id([7311724, 7478, 5957, 5950, 5793])
    assert res.to_dict() == {
        5793: "C02855",
        5950: "C03047",
        5957: "C03057",
        7478: "C04936",
    }


def test_parse_compound_pathway_mapping():
    raw_mapping = (
        "cpd:C00031\tpath:map00010\n"
        "cpd:C00031\tpath:map00030\n"
        "cpd:C00031\tpath:map00052\n"
        "cpd:C00022\tpath:map00010\n"
        "cpd:C00022\tpath:map00030\n"
    )
    compound_to_pathways = parse_compound_pathway_mapping(raw_mapping)
    assert compound_to_pathways == {
        "cpd:C00031": ["path:map00010", "path:map00030", "path:map00052"],
        "cpd:C00022": ["path:map00010", "path:map00030"],
    }


def test_parse_kegg_name_description():
    # https://rest.kegg.jp/get/path:map00030+path:map00010
    raw_text = textwrap.dedent("""
    ENTRY       map00030                    Pathway
    NAME        Pentose phosphate pathway
    DESCRIPTION The pentose phosphate pathway is a process of glucose turnover that produces NADPH as reducing equivalents and pentoses as essential parts of nucleotides.

    (...)
    ///

    ENTRY       map00010                    Pathway
    NAME        Glycolysis / Gluconeogenesis
    DESCRIPTION Glycolysis is the process of converting glucose into pyruvate and generating small amounts of ATP (energy) and NADH (reducing power).

    (...)
    """)

    parsed_entries = parse_kegg_name_description(raw_text)
    print(parsed_entries)
    expected = {
        "map00030": {
            "NAME": "Pentose phosphate pathway",
            "DESCRIPTION": (
                "The pentose phosphate pathway is a process of glucose "
                "turnover that produces NADPH as reducing equivalents and pentoses as "
                "essential parts of nucleotides."
            ),
        },
        "map00010": {
            "NAME": "Glycolysis / Gluconeogenesis",
            "DESCRIPTION": (
                "Glycolysis is the process of converting glucose into "
                "pyruvate and generating small amounts of ATP (energy) "
                "and NADH (reducing power)."
            ),
        },
    }
    assert parsed_entries == expected
