"""Tests which query the UniProt REST API. They are run in the separate
'Test API notebooks' workflow (.github/workflows/api_notebooks.yml)."""

from acore.io.uniprot import fetch_annotations, process_annotations

FIELDS = "go_p,go_c,go_f,go"


def test_fetch_and_process_annotations():
    ids = ["P62699", "Q9UID3"]
    annotations = fetch_annotations(ids, fields=FIELDS)
    assert set(annotations["From"]) == set(ids)

    annotations = process_annotations(annotations, fields=FIELDS)
    assert set(annotations["identifier"]) == set(ids)
    assert annotations["annotation"].str.contains("GO:").all()
