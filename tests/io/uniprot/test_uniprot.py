import pathlib

import pandas as pd

from acore.io.uniprot import process_annotations

folder = pathlib.Path(__file__).parent


def test_process_annotations():
    annotations_fetched = pd.read_csv(folder / "annotations_50_fetched.csv")
    annotations_expected = pd.read_csv(folder / "annotations_50.csv")

    annotations_actual = process_annotations(
        annotations_fetched, fields="go_p,go_c,go_f,go"
    )

    pd.testing.assert_frame_equal(annotations_actual, annotations_expected)
