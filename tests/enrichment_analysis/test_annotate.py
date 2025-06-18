import numpy as np
import pandas as pd

import acore.enrichment_analysis as ea


def test_annotate_features():
    expected = pd.Series(
        [
            "foreground",
            "foreground",
            "background",
            "foreground",
            "background",
            "background",
            np.nan,
        ]
    )

    features = pd.Series(["G1", "G2", "G3", "G4", "G5", "G6", "G9"])
    in_foreground = ["G1", "G2", "G4"]
    in_background = ["G3", "G5", "G6"]
    actual = ea.annotate_features(features, in_foreground, in_background)
    pd.testing.assert_series_equal(expected, actual)


def test_annotate_features_with_duplicates():
    """for example if multiple peptides are associated with the same protein."""
    expected = pd.Series(
        [
            "foreground",
            "foreground",
            "background",
            "background",
            "foreground",
            "background",
            "background",
            np.nan,
        ]
    )

    features = pd.Series(["G1", "G2", "G3", "G3", "G4", "G5", "G6", "G9"])
    in_foreground = ["G1", "G2", "G4"]
    in_background = ["G3", "G5", "G6"]
    actual = ea.annotate_features(features, in_foreground, in_background)
    pd.testing.assert_series_equal(expected, actual)
