import pathlib

import pandas as pd

from acore.io.uniprot import filter_annotations, process_annotations

folder = pathlib.Path(__file__).parent


def test_process_annotations():
    annotations_fetched = pd.read_csv(folder / "annotations_50_fetched.csv")
    annotations_expected = pd.read_csv(folder / "annotations_50.csv")

    annotations_actual = process_annotations(
        annotations_fetched, fields="go_p,go_c,go_f,go"
    )

    sort_cols = ["identifier", "source", "annotation"]
    annotations_actual = annotations_actual.sort_values(sort_cols).reset_index(
        drop=True
    )
    annotations_expected = annotations_expected.sort_values(sort_cols).reset_index(
        drop=True
    )

    pd.testing.assert_frame_equal(annotations_actual, annotations_expected)


def test_filter_annotations_any_keyword():
    df = pd.DataFrame(
        {
            "identifier": ["P1", "P1", "P2"],
            "source": ["go_p", "go_p", "go_p"],
            "annotation": [
                "cell apoptosis [GO:0042981]",
                "cell division",
                "mitochondria",
            ],
        }
    )

    result = filter_annotations(df, keywords=["apoptosis", "mitochondr"])
    assert set(result["annotation"]) == {
        "cell apoptosis [GO:0042981]",
        "mitochondria",
    }
