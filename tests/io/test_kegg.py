"""Test static functions for KEGG."""

from acore.io.kegg import lookup_cid_to_kegg_id


def test_lookup_cid_to_kegg_id():
    res = lookup_cid_to_kegg_id([7311724, 7478, 5957, 5950, 5793])
    assert res.to_dict() == {
        5793: "C02855",
        5950: "C03047",
        5957: "C03057",
        7478: "C04936",
    }
